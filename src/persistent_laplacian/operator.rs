//! Fixed row-space projection and factorized persistent action.
use super::{LaplacianComponent, SpectralOptions, allocation, basis_ids, numerical, validate_pair};
use crate::complex::{Simplex, SimplexId};
use crate::execution::{Execution, WorkBudget};
use crate::filtration::{Coverage, FiltrationContext, SimplicialStage};
use crate::{Error, Result};

/// Storage and numerical limits for a factorized real persistent action.
#[derive(Clone, Copy, Debug)]
pub struct OperatorOptions {
    /// Numerical policy and per-buffer dense materialization cap.
    /// Row-space reduction additionally requires relative tolerance >= 128 eps.
    pub spectral: SpectralOptions,
    /// Maximum stored incidence entries, including temporary outside rows.
    pub max_incidence_entries: usize,
    /// Maximum cumulative f64 entries in the retained row-space basis.
    /// The basis can be dense; this is not a total process memory quota.
    pub max_projection_entries: usize,
    /// Maximum entries in any action, projection or coordinate vector.
    pub max_vector_entries: usize,
}
impl Default for OperatorOptions {
    fn default() -> Self {
        Self {
            spectral: SpectralOptions::default(),
            max_incidence_entries: 1_048_576,
            max_projection_entries: 65_536,
            max_vector_entries: 65_536,
        }
    }
}

/// Owned real pair operator applied through boundaries and a fixed projector.
///
/// Computes B1 (I-Q^T Q) B1^T + A^T A, where rows of Q approximate an
/// orthonormal basis of row(B2), all upper columns are retained, and A is
/// boundary_q(K). Two-pass modified Gram--Schmidt fixes Q during construction;
/// actions never use adaptive inner iterations. This is a numerical projection,
/// not an exact rank certificate. Discarded rows and orthogonality are checked.
/// There is no ridge, implicit dense fallback or spectral completeness claim.
///
/// Ordered simplex labels, degree, inclusive scales and source certificates are
/// owned. Standard inner products, increasing-vertex orientations and ordinary
/// nonaugmented chains match [`super::PersistentLaplacian`]. This operator does
/// not change that type's borrowed dense `matrix()` contract.
#[derive(Clone, Debug)]
pub struct PersistentLaplacianOperator {
    degree: usize,
    scales: (f64, f64),
    basis: Vec<Simplex>,
    target_basis: Vec<Simplex>,
    source_context: Option<FiltrationContext>,
    source_coverage: Coverage,
    construction_dimension: Option<usize>,
    upper: Vec<f64>,
    lower: Vec<f64>,
    upper_columns: usize,
    lower_rows: usize,
    row_basis: Vec<Vec<f64>>,
    projection_residual: f64,
    options: OperatorOptions,
}
impl PersistentLaplacianOperator {
    /// Snapshot a legal same-source inclusion and construct its fixed projector.
    ///
    /// Uses one budget for source traversal, incidence, orthogonalization and
    /// projector verification. Repeated action terminals start fresh budgets.
    /// All temporary factors have explicit entry caps; they may still be dense.
    /// # Errors
    /// Rejects invalid policy/inclusion/skeleton, size overflow, entry limits,
    /// allocation failure, cancellation, work exhaustion, ambiguous nonzero
    /// residuals between 64 eps and the rank cutoff, or failed projection checks.
    pub fn new(
        k: &SimplicialStage<'_>,
        l: &SimplicialStage<'_>,
        degree: usize,
        options: OperatorOptions,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        validate_pair(
            k,
            l,
            degree,
            LaplacianComponent::Full,
            &options.spectral,
            &mut budget,
        )?;
        if options.spectral.relative_tolerance < 128. * f64::EPSILON {
            return Err(Error::InvalidParameter {
                parameter: "relative_tolerance",
                reason: "row-space projection requires at least 128 machine epsilons",
            });
        }
        let target = basis_ids(l, degree, &mut budget)?;
        let n = target.iter().take_while(|id| id.index() < k.len()).count();
        let next = basis_ids(l, degree + 1, &mut budget)?;
        let previous = if degree == 0 {
            Vec::new()
        } else {
            basis_ids(k, degree - 1, &mut budget)?
        };
        for size in [target.len(), next.len(), previous.len()] {
            check_limit(size, options.max_vector_entries, "max_vector_entries")?;
        }
        let p = next.len();
        let entries = target
            .len()
            .checked_mul(p)
            .and_then(|v| previous.len().checked_mul(n).and_then(|w| v.checked_add(w)))
            .ok_or(Error::SizeOverflow {
                operation: "operator incidence entries",
            })?;
        check_limit(
            entries,
            options.max_incidence_entries,
            "max_incidence_entries",
        )?;
        let mut upper = super::dense(n, p, &options.spectral, &mut budget)?;
        let mut outside = super::dense(target.len() - n, p, &options.spectral, &mut budget)?;
        let mut lower = super::dense(previous.len(), n, &options.spectral, &mut budget)?;
        for (column, &id) in next.iter().enumerate() {
            budget.step()?;
            for term in l.boundary(id).ok_or(Error::InternalInvariant {
                reason: "upper boundary missing",
            })? {
                budget.step()?;
                let row =
                    target
                        .binary_search(&term.face)
                        .map_err(|_| Error::InternalInvariant {
                            reason: "upper face missing",
                        })?;
                if row < n {
                    upper[row * p + column] = f64::from(term.coefficient);
                } else {
                    outside[(row - n) * p + column] = f64::from(term.coefficient);
                }
            }
        }
        for (column, &id) in target[..n].iter().enumerate() {
            budget.step()?;
            for term in k.boundary(id).ok_or(Error::InternalInvariant {
                reason: "lower boundary missing",
            })? {
                budget.step()?;
                let row =
                    previous
                        .binary_search(&term.face)
                        .map_err(|_| Error::InternalInvariant {
                            reason: "lower face missing",
                        })?;
                lower[row * n + column] = f64::from(term.coefficient);
            }
        }
        let mut row_basis: Vec<Vec<f64>> = Vec::new();
        let mut projection_entries = 0usize;
        let mut projection_residual: f64 = 0.;
        for row in outside.chunks(p.max(1)) {
            let mut v = vector(p, &options, &mut budget)?;
            v.copy_from_slice(row);
            let original_norm = length(row, &mut budget)?;
            for _ in 0..2 {
                remove_rows(&mut v, &row_basis, &mut budget)?;
            }
            let size = length(&v, &mut budget)?;
            if size <= options.spectral.relative_tolerance * original_norm {
                if size > 64. * f64::EPSILON * original_norm.max(1.) {
                    return Err(numerical("ambiguous row-space rank"));
                }
                projection_residual = projection_residual.max(size / original_norm.max(1.));
                continue;
            }
            projection_entries = projection_entries
                .checked_add(p)
                .ok_or(Error::SizeOverflow {
                    operation: "projection basis entries",
                })?;
            check_limit(
                projection_entries,
                options.max_projection_entries,
                "max_projection_entries",
            )?;
            for value in &mut v {
                budget.step()?;
                *value /= size;
            }
            for q in &row_basis {
                projection_residual = projection_residual.max(dot(&v, q, &mut budget)?.abs());
            }
            projection_residual = projection_residual.max((dot(&v, &v, &mut budget)? - 1.).abs());
            row_basis.try_reserve(1).map_err(|_| allocation())?;
            row_basis.push(v);
        }
        if !projection_residual.is_finite()
            || projection_residual > options.spectral.residual_tolerance
        {
            return Err(numerical("row-space projection residual"));
        }
        let mut owned = |ids: &[SimplexId]| -> Result<Vec<Simplex>> {
            let mut basis = Vec::new();
            basis
                .try_reserve_exact(ids.len())
                .map_err(|_| allocation())?;
            for &id in ids {
                budget.step()?;
                basis.push(
                    l.simplex(id)
                        .ok_or(Error::InternalInvariant {
                            reason: "operator basis missing",
                        })?
                        .clone(),
                );
            }
            Ok(basis)
        };
        let basis = owned(&target[..n])?;
        let target_basis = owned(&target)?;
        budget.check()?;
        Ok(Self {
            degree,
            scales: (k.scale(), l.scale()),
            basis,
            target_basis,
            source_context: l.source_filtration().map(|s| s.context().clone()),
            source_coverage: l
                .source_filtration()
                .map_or(Coverage::Complete, |s| s.coverage()),
            construction_dimension: l.source_filtration().map(|s| s.max_simplex_dimension()),
            upper,
            lower,
            upper_columns: p,
            lower_rows: previous.len(),
            row_basis,
            projection_residual,
            options,
        })
    }

    /// Ordered q-basis of K, with original labels and births.
    pub fn basis(&self) -> &[Simplex] {
        &self.basis
    }
    /// Ordered q-basis of L before projection, with original labels and births.
    pub fn target_basis(&self) -> &[Simplex] {
        &self.target_basis
    }
    /// Real chain degree.
    pub fn degree(&self) -> usize {
        self.degree
    }
    /// Inclusive source/target scales.
    pub fn scales(&self) -> (f64, f64) {
        self.scales
    }
    /// Owned source context; None for bare supplied topology.
    pub fn source_context(&self) -> Option<&FiltrationContext> {
        self.source_context.as_ref()
    }
    /// Original source coverage, independent of the selected pair scales.
    pub fn source_coverage(&self) -> Coverage {
        self.source_coverage
    }
    /// Certified source construction dimension, if present.
    pub fn construction_dimension(&self) -> Option<usize> {
        self.construction_dimension
    }
    /// Applied limits and numerical policy.
    pub fn options(&self) -> OperatorOptions {
        self.options
    }
    /// Number of stored f64 incidence entries, excluding the projection basis.
    pub fn incidence_entries(&self) -> usize {
        self.upper.len() + self.lower.len()
    }
    /// Cumulative f64 entries retained in Q, which may be dense.
    pub fn projection_entries(&self) -> usize {
        self.row_basis.len() * self.upper_columns
    }
    /// Numerical row-space rank, not an exact real rank certificate.
    pub fn projection_rank(&self) -> usize {
        self.row_basis.len()
    }
    /// Maximum normalized discarded-row and basis-orthogonality diagnostic.
    pub fn projection_residual(&self) -> f64 {
        self.projection_residual
    }

    /// Apply in `basis()` coordinates without forming the persistent matrix.
    /// # Errors
    /// Rejects wrong shape, nonfinite input, overflow, allocation/cap failure,
    /// cancellation and work exhaustion. Each call shares one budget across
    /// input validation, both boundary actions and the fixed projection.
    pub fn action(&self, x: &[f64], execution: &Execution<'_>) -> Result<Vec<f64>> {
        let mut budget = WorkBudget::new(execution)?;
        self.action_with(x, &mut budget)
    }

    fn action_with(&self, x: &[f64], budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
        let n = self.basis.len();
        if x.len() != n {
            return Err(Error::ShapeMismatch {
                input: "Laplacian coordinates",
                expected: n,
                actual: x.len(),
            });
        }
        for (i, value) in x.iter().enumerate() {
            budget.step()?;
            if !value.is_finite() {
                return Err(Error::NonFiniteValue {
                    field: "Laplacian coordinates",
                    index: Some(i),
                });
            }
        }
        let p = self.upper_columns;
        let mut y = vector(p, &self.options, budget)?;
        for (i, &value) in x.iter().enumerate() {
            for (j, target) in y.iter_mut().enumerate() {
                budget.step()?;
                *target += self.upper[i * p + j] * value;
            }
        }
        // Compute all coefficients against the original y: I-Q^T Q is a fixed
        // symmetric matrix even with small orthogonality errors in Q.
        let mut coefficients = vector(self.row_basis.len(), &self.options, budget)?;
        for (coefficient, q) in coefficients.iter_mut().zip(&self.row_basis) {
            *coefficient = dot(q, &y, budget)?;
        }
        for (&coefficient, q) in coefficients.iter().zip(&self.row_basis) {
            for (value, &entry) in y.iter_mut().zip(q) {
                budget.step()?;
                *value -= coefficient * entry;
            }
        }
        let mut result = vector(n, &self.options, budget)?;
        for (i, value) in result.iter_mut().enumerate() {
            *value = dot(&self.upper[i * p..(i + 1) * p], &y, budget)?;
        }
        let mut lower_action = vector(self.lower_rows, &self.options, budget)?;
        for (row, value) in lower_action.iter_mut().enumerate() {
            *value = dot(&self.lower[row * n..(row + 1) * n], x, budget)?;
        }
        for (row, &value) in lower_action.iter().enumerate() {
            for (j, target) in result.iter_mut().enumerate() {
                budget.step()?;
                *target += self.lower[row * n + j] * value;
            }
        }
        for value in result.iter().chain(&y).chain(&lower_action) {
            budget.step()?;
            if !value.is_finite() {
                return Err(numerical("factorized Laplacian action"));
            }
        }
        budget.check()?;
        Ok(result)
    }

    /// Explicit bounded dense export in row-major `basis()` order.
    ///
    /// Both `max_entries` and the configured dense limit apply. One shared
    /// budget covers every basis action; no per-column budget is restarted.
    /// # Errors
    /// Includes shape overflow, cap/allocation failure, numerical errors,
    /// cancellation and cumulative work exhaustion. No partial matrix is returned.
    pub fn materialize(&self, max_entries: usize, execution: &Execution<'_>) -> Result<Vec<f64>> {
        let mut budget = WorkBudget::new(execution)?;
        let n = self.basis.len();
        let entries = n.checked_mul(n).ok_or(Error::SizeOverflow {
            operation: "operator materialization",
        })?;
        check_limit(entries, max_entries, "max_entries")?;
        let mut matrix = super::dense(n, n, &self.options.spectral, &mut budget)?;
        let mut x = vector(n, &self.options, &mut budget)?;
        for j in 0..n {
            x[j] = 1.;
            let result = self.action_with(&x, &mut budget)?;
            for i in 0..n {
                budget.step()?;
                matrix[i * n + j] = result[i];
            }
            x[j] = 0.;
        }
        budget.check()?;
        Ok(matrix)
    }
}

fn check_limit(entries: usize, limit: usize, parameter: &'static str) -> Result<()> {
    if entries > limit {
        return Err(Error::InvalidParameter {
            parameter,
            reason: "operator storage exceeds the declared limit",
        });
    }
    Ok(())
}
fn vector(n: usize, options: &OperatorOptions, budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
    check_limit(n, options.max_vector_entries, "max_vector_entries")?;
    budget.check()?;
    let mut v = Vec::new();
    v.try_reserve_exact(n).map_err(|_| allocation())?;
    v.resize(n, 0.);
    budget.check()?;
    Ok(v)
}
fn dot(a: &[f64], b: &[f64], budget: &mut WorkBudget<'_>) -> Result<f64> {
    let mut result = 0.;
    for (&x, &y) in a.iter().zip(b) {
        budget.step()?;
        result += x * y;
    }
    if !result.is_finite() {
        return Err(numerical("projection dot product"));
    }
    Ok(result)
}
fn length(v: &[f64], budget: &mut WorkBudget<'_>) -> Result<f64> {
    let mut norm: f64 = 0.;
    for &value in v {
        budget.step()?;
        norm = norm.hypot(value);
    }
    Ok(norm)
}
fn remove_rows(v: &mut [f64], rows: &[Vec<f64>], budget: &mut WorkBudget<'_>) -> Result<()> {
    for q in rows {
        let coefficient = dot(q, v, budget)?;
        for (value, &entry) in v.iter_mut().zip(q) {
            budget.step()?;
            *value -= coefficient * entry;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicBool;

    #[test]
    fn cancellation_during_actual_projection_work() {
        let flag = AtomicBool::new(false);
        let execution = Execution::new(Some(1000), Some(&flag));
        let mut budget = WorkBudget::new(&execution).unwrap();
        budget.cancel_at_work(20);
        let mut v = vec![1.; 100];
        let rows = vec![vec![0.1; 100]];
        assert!(matches!(
            remove_rows(&mut v, &rows, &mut budget),
            Err(Error::Cancelled)
        ));
    }
}
