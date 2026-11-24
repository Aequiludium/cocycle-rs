//! Fixed row-space projection and factorized persistent action.
use super::{LaplacianComponent, SpectralOptions, allocation, basis_ids, numerical, validate_pair};
use crate::complex::{Simplex, SimplexId};
use crate::execution::{Execution, WorkBudget};
use crate::filtration::{Coverage, FiltrationContext, SimplicialStage};
use crate::{Error, Result};

/// Which part of the fixed operator's spectrum to approximate.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum PersistentSpectrumTarget {
    /// Algebraically smallest modes, including harmonic zeros.
    Smallest,
    /// Algebraically largest modes, in descending order.
    Largest,
    /// Smallest modes of a numerically resolved boundary-image restriction.
    SmallestPositive,
}

/// Owned approximate eigenpairs in the original persistent-pair coordinates.
///
/// No full persistent matrix or complete spectrum is computed. Numerical range
/// preparation, when requested, can still require O(n*rank) storage.
pub struct PersistentEigenpairs {
    /// Original chain degree.
    pub degree: usize,
    /// Inclusive source and target scales.
    pub scales: (f64, f64),
    /// Owned ordered domain simplices, preserving labels and births.
    pub basis: Vec<Simplex>,
    /// Owned ordered target q-basis.
    pub target_basis: Vec<Simplex>,
    /// Original source context, when certified by a filtration.
    pub source_context: Option<FiltrationContext>,
    /// Original construction coverage.
    pub source_coverage: Coverage,
    /// Certified constructed skeleton dimension, when present.
    pub construction_dimension: Option<usize>,
    /// Requested spectral selection.
    pub target: PersistentSpectrumTarget,
    /// Ritz values, ascending except for Largest.
    pub values: Vec<f64>,
    /// Row-major basis.len() by values.len() vector columns.
    pub vectors: Vec<f64>,
    /// Fresh original-action residual norms, in eigenvalue units.
    pub residuals: Vec<f64>,
    /// Frobenius norm of V^T V-I in the original coordinates.
    pub orthogonality_error: f64,
    /// Completed block updates.
    pub iterations: usize,
    /// Original operator calls, including final verification.
    pub operator_applications: usize,
    /// Largest iterative Rayleigh-Ritz subspace.
    pub peak_subspace_dimension: usize,
    /// Resolved numerical image rank for positive requests; never an exact Betti certificate.
    pub range_rank: Option<usize>,
    /// Maximum normalized leakage of returned images outside the prepared range.
    pub range_residual: f64,
    /// Positive cutoff: absolute_tolerance + relative_tolerance * factor energy.
    /// Factor energy is a trace bound, not a computed full-spectrum norm.
    pub zero_threshold: Option<f64>,
    /// Iterative policy used for this request.
    pub options: crate::eigensolver::IterativeOptions,
}
impl PersistentEigenpairs {
    /// Always false: residual convergence does not certify extremal indexes.
    pub fn is_index_coverage_certified(&self) -> bool {
        false
    }
}

/// Storage and numerical limits for a factorized real persistent action.
#[derive(Clone, Copy, Debug)]
pub struct OperatorOptions {
    /// Numerical policy and per-buffer dense materialization cap.
    /// Row-space reduction additionally requires relative tolerance >= 128 eps.
    pub spectral: SpectralOptions,
    /// Maximum stored incidence entries, including temporary outside rows.
    pub max_incidence_entries: usize,
    /// Maximum cumulative scalar entries in the retained projection factor.
    /// Counts dense Q entries or compact signs plus component sizes; integer maps
    /// coexist. The factor can be dense; this is not a process memory quota.
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
/// Boundaries use sparse signed columns and their actual transposes. A single
/// action costs O(nnz + p*rank(B2) + coordinate counts); Q may fill in. Neither
/// construction nor action forms an N*N, n*n or p*p operator/projector buffer.
#[derive(Clone, Debug)]
pub struct PersistentLaplacianOperator {
    degree: usize,
    scales: (f64, f64),
    basis: Vec<Simplex>,
    target_basis: Vec<Simplex>,
    source_context: Option<FiltrationContext>,
    source_coverage: Coverage,
    construction_dimension: Option<usize>,
    upper: Vec<Vec<(usize, f64)>>,
    lower: Vec<Vec<(usize, f64)>>,
    upper_columns: usize,
    lower_rows: usize,
    row_basis: Vec<Vec<f64>>,
    compact: Option<super::nonbranching::KernelProjection>,
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
        Self::new_with(k, l, degree, options, &mut budget)
    }

    pub(super) fn new_with(
        k: &SimplicialStage<'_>,
        l: &SimplicialStage<'_>,
        degree: usize,
        options: OperatorOptions,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Self> {
        Self::prepare(k, l, degree, options, false, budget)
    }

    /// Prepare a compact signed-component projector for an eligible pair.
    ///
    /// Outside rows must have at most two signed unit entries. Disjoint kernel
    /// supports use O(p) storage and projection work; no dense R/E/V/Z or Q is
    /// formed. Other incidence, partial-spectrum and metadata contracts match
    /// `new`. The entry cap counts signs and component sizes; integer maps and
    /// bounded adjacency also coexist. This explicitly selected route never
    /// falls back. Incremental starts must remain non-branching.
    /// # Errors
    /// Includes `new` errors and ineligible outside constraints or adjacency caps.
    pub fn new_non_branching(
        k: &SimplicialStage<'_>,
        l: &SimplicialStage<'_>,
        degree: usize,
        options: OperatorOptions,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        Self::prepare(
            k,
            l,
            degree,
            options,
            true,
            &mut WorkBudget::new(execution)?,
        )
    }

    fn prepare(
        k: &SimplicialStage<'_>,
        l: &SimplicialStage<'_>,
        degree: usize,
        options: OperatorOptions,
        compact_projection: bool,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Self> {
        validate_pair(
            k,
            l,
            degree,
            LaplacianComponent::Full,
            &options.spectral,
            budget,
        )?;
        if options.spectral.relative_tolerance < 128. * f64::EPSILON {
            return Err(Error::InvalidParameter {
                parameter: "relative_tolerance",
                reason: "row-space projection requires at least 128 machine epsilons",
            });
        }
        let target = basis_ids(l, degree, budget)?;
        let n = target.iter().take_while(|id| id.index() < k.len()).count();
        let next = basis_ids(l, degree + 1, budget)?;
        let previous = if degree == 0 {
            Vec::new()
        } else {
            basis_ids(k, degree - 1, budget)?
        };
        for size in [target.len(), next.len(), previous.len()] {
            check_limit(size, options.max_vector_entries, "max_vector_entries")?;
        }
        let p = next.len();
        let mut upper = incidence(p, budget)?;
        let mut outside = incidence(target.len() - n, budget)?;
        let mut lower = incidence(n, budget)?;
        let mut entries = 0usize;
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
                    push_term(
                        &mut upper[column],
                        (row, f64::from(term.coefficient)),
                        &mut entries,
                        &options,
                    )?;
                } else {
                    push_term(
                        &mut outside[row - n],
                        (column, f64::from(term.coefficient)),
                        &mut entries,
                        &options,
                    )?;
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
                push_term(
                    &mut lower[column],
                    (row, f64::from(term.coefficient)),
                    &mut entries,
                    &options,
                )?;
            }
        }
        let mut row_basis: Vec<Vec<f64>> = Vec::new();
        let mut projection_entries = 0usize;
        let mut projection_residual: f64 = 0.;
        let compact = if compact_projection {
            Some(super::nonbranching::KernelProjection::new(
                &outside,
                p,
                options.max_projection_entries,
                options.max_incidence_entries,
                budget,
            )?)
        } else {
            None
        };
        for row in outside.iter().filter(|_| !compact_projection) {
            budget.step()?;
            if row.is_empty() {
                continue;
            }
            let mut v = vector(p, &options, budget)?;
            for &(column, value) in row {
                budget.step()?;
                v[column] = value;
            }
            add_projection_row(
                v,
                &mut row_basis,
                &mut projection_entries,
                &mut projection_residual,
                &options,
                budget,
            )?;
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
            compact,
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
    /// Number of stored nonzero incidence terms, excluding the projection basis.
    pub fn incidence_entries(&self) -> usize {
        self.upper.iter().chain(&self.lower).map(Vec::len).sum()
    }
    /// Retained dense Q entries or compact signs plus component sizes.
    pub fn projection_entries(&self) -> usize {
        self.compact
            .as_ref()
            .map_or(self.row_basis.len() * self.upper_columns, |kernel| {
                kernel.entries()
            })
    }
    /// Outside rank: numerical for Q, exact signed-component rank for the compact route.
    pub fn projection_rank(&self) -> usize {
        self.compact
            .as_ref()
            .map_or(self.row_basis.len(), |kernel| kernel.rank())
    }
    /// Maximum normalized discarded-row and basis-orthogonality diagnostic.
    pub fn projection_residual(&self) -> f64 {
        self.projection_residual
    }

    /// Consume this snapshot and move its start backward at the same endpoint.
    ///
    /// Newly excluded q-rows append constraints to the existing fixed row-space
    /// factor. Existing factors are retained, not recomputed. The domain and its
    /// lower boundary are truncated in the original filtration order; unused
    /// previous-degree rows remain zero. The immutable source is not needed.
    /// This is a fixed-endpoint workspace update, not online endpoint insertion.
    /// # Errors
    /// Rejects nonfinite/forward scales, factor/incidence limits, ambiguous
    /// projection updates, allocation failure, cancellation and work exhaustion.
    /// Failure consumes the input snapshot and returns no partially updated one.
    pub fn restrict_start(mut self, scale: f64, execution: &Execution<'_>) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        self.restrict_with(scale, &mut budget)?;
        Ok(self)
    }

    pub(super) fn reuse_upper(
        k: &SimplicialStage<'_>,
        l: &SimplicialStage<'_>,
        degree: usize,
        upper: Option<&Self>,
        options: OperatorOptions,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Self> {
        validate_pair(
            k,
            l,
            degree,
            LaplacianComponent::Full,
            &options.spectral,
            budget,
        )?;
        let mut op = Self::new_with(k, k, degree, options, budget)?;
        op.upper.clear();
        op.row_basis.clear();
        op.upper_columns = 0;
        op.projection_residual = 0.;
        if let Some(source) = upper {
            op.upper = incidence(source.upper.len(), budget)?;
            let mut entries: usize = op.lower.iter().map(Vec::len).sum();
            for (column, terms) in op.upper.iter_mut().zip(&source.upper) {
                for &term in terms {
                    budget.step()?;
                    push_term(column, term, &mut entries, &options)?;
                }
            }
            op.upper_columns = source.upper_columns;
            op.row_basis
                .try_reserve_exact(source.row_basis.len())
                .map_err(|_| allocation())?;
            for q in &source.row_basis {
                let mut row = vector(q.len(), &options, budget)?;
                budget.step_by(q.len())?;
                row.copy_from_slice(q);
                op.row_basis.push(row);
            }
            op.projection_residual = source.projection_residual;
            op.compact = source
                .compact
                .as_ref()
                .map(|kernel| kernel.copy(budget))
                .transpose()?;
        }
        let ids = basis_ids(l, degree, budget)?;
        check_limit(ids.len(), options.max_vector_entries, "max_vector_entries")?;
        op.target_basis.clear();
        op.target_basis
            .try_reserve_exact(ids.len())
            .map_err(|_| allocation())?;
        for id in ids {
            budget.step()?;
            op.target_basis.push(
                l.simplex(id)
                    .ok_or(Error::InternalInvariant {
                        reason: "reused upper target basis",
                    })?
                    .clone(),
            );
        }
        op.scales = (k.scale(), l.scale());
        budget.check()?;
        Ok(op)
    }

    fn restrict_with(&mut self, scale: f64, budget: &mut WorkBudget<'_>) -> Result<()> {
        budget.step()?;
        if !scale.is_finite() || scale > self.scales.0 {
            return Err(Error::InvalidParameter {
                parameter: "operator start scale",
                reason: "requires a finite nonincreasing scale",
            });
        }
        budget.step_by(self.basis.len())?;
        let n = self
            .basis
            .partition_point(|simplex| simplex.value() <= scale);
        if n == self.basis.len() {
            self.scales.0 = scale;
            return budget.check();
        }
        let mut outside = incidence(self.basis.len() - n, budget)?;
        let mut entries = self.incidence_entries();
        for (j, terms) in self.upper.iter().enumerate() {
            for &(row, value) in terms {
                budget.step()?;
                if row >= n {
                    push_term(
                        &mut outside[row - n],
                        (j, value),
                        &mut entries,
                        &self.options,
                    )?;
                }
            }
        }
        if let Some(kernel) = self.compact.take() {
            let mut constraints = kernel.constraints(budget)?;
            constraints
                .try_reserve_exact(outside.len())
                .map_err(|_| allocation())?;
            constraints.extend(outside);
            self.compact = Some(super::nonbranching::KernelProjection::new(
                &constraints,
                self.upper_columns,
                self.options.max_projection_entries,
                self.options.max_incidence_entries,
                budget,
            )?);
        } else {
            let mut projection_entries = self.projection_entries();
            for terms in outside {
                if terms.is_empty() {
                    continue;
                }
                let mut v = vector(self.upper_columns, &self.options, budget)?;
                for (column, value) in terms {
                    budget.step()?;
                    v[column] = value;
                }
                add_projection_row(
                    v,
                    &mut self.row_basis,
                    &mut projection_entries,
                    &mut self.projection_residual,
                    &self.options,
                    budget,
                )?;
            }
        }
        if self.projection_residual > self.options.spectral.residual_tolerance {
            return Err(numerical("incremental projection residual"));
        }
        for terms in &mut self.upper {
            budget.step_by(terms.len())?;
            terms.retain(|&(row, _)| row < n);
        }
        self.lower.truncate(n);
        self.basis.truncate(n);
        self.scales.0 = scale;
        budget.check()
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

    /// Approximate partial spectra directly through the fixed sparse action.
    ///
    /// `initial`, if supplied, is row-major n by min(block_size,n).
    /// For positive requests original-coordinate columns are projected into
    /// the range; at most min(block_size,range_rank) columns are used. Only vectors
    /// are reused: all Ritz values and residuals use this current operator.
    ///
    /// Smallest/Largest reuse the shared block LOBPCG without an n by n matrix.
    /// SmallestPositive streams columns of [B1 P, A^T] through two-pass QR,
    /// then applies the operator in this range without a dense reduced matrix
    /// or an initial complete eigensolve. The range is capped by
    /// max_basis_entries and can be O(n*rank); it is not an implicit sparse QR.
    /// Its numerical rank uses this operator's spectral policy relative to
    /// max(1, largest generator norm), for these standard signed-unit chains.
    /// This prevents normalization of roundoff-only annihilated columns. Positivity uses
    /// the iterative tolerances and the total squared generator norms, an upper
    /// bound for the exact PSD operator's spectral norm. Ambiguity is an error.
    ///
    /// A single budget covers range setup, starts, iteration, fresh original
    /// residuals and owned output metadata. There is no dense fallback. Small
    /// residuals do not certify smallest indexes or exact real nullity.
    /// # Errors
    /// Rejects invalid options/counts/starts, storage or work exhaustion,
    /// cancellation, ambiguous/insufficient positive range, nonconvergence,
    /// allocation/arithmetic failure or inaccurate original-coordinate output.
    pub fn eigenpairs(
        &self,
        count: usize,
        target: PersistentSpectrumTarget,
        initial: Option<&[f64]>,
        options: &crate::eigensolver::IterativeOptions,
        execution: &Execution<'_>,
    ) -> Result<PersistentEigenpairs> {
        let mut budget = WorkBudget::new(execution)?;
        self.eigenpairs_with(count, target, initial, options, &mut budget)
    }

    pub(super) fn eigenpairs_with(
        &self,
        count: usize,
        target: PersistentSpectrumTarget,
        initial: Option<&[f64]>,
        options: &crate::eigensolver::IterativeOptions,
        budget: &mut WorkBudget<'_>,
    ) -> Result<PersistentEigenpairs> {
        // Validate the shared policy before potentially expensive range setup.
        crate::eigensolver::lobpcg(self.basis.len(), 0, false, options, budget, |_, _| {
            Err(numerical("unexpected zero-count action"))
        })?;
        let (pairs, range_rank, range_residual, zero_threshold) =
            if target == PersistentSpectrumTarget::SmallestPositive {
                self.positive_pairs(count, initial, options, budget)?
            } else {
                (
                    crate::eigensolver::lobpcg_with_initial(
                        self.basis.len(),
                        count,
                        target == PersistentSpectrumTarget::Largest,
                        options,
                        initial,
                        budget,
                        |x, b| self.action_with(x, b),
                    )?,
                    None,
                    0.,
                    None,
                )
            };
        let mut owned = |input: &[Simplex]| -> Result<Vec<Simplex>> {
            let mut output = Vec::new();
            output
                .try_reserve_exact(input.len())
                .map_err(|_| allocation())?;
            for simplex in input {
                budget.step()?;
                output.push(simplex.clone());
            }
            Ok(output)
        };
        let basis = owned(&self.basis)?;
        let target_basis = owned(&self.target_basis)?;
        budget.check()?;
        Ok(PersistentEigenpairs {
            degree: self.degree,
            scales: self.scales,
            basis,
            target_basis,
            source_context: self.source_context.clone(),
            source_coverage: self.source_coverage,
            construction_dimension: self.construction_dimension,
            target,
            values: pairs.values,
            vectors: pairs.vectors,
            residuals: pairs.residuals,
            orthogonality_error: pairs.orthogonality_error,
            iterations: pairs.iterations,
            operator_applications: pairs.actions,
            peak_subspace_dimension: pairs.peak_subspace,
            range_rank,
            range_residual,
            zero_threshold,
            options: *options,
        })
    }

    fn project_upper(&self, y: &mut [f64], budget: &mut WorkBudget<'_>) -> Result<()> {
        if let Some(kernel) = &self.compact {
            return kernel.apply(y, budget);
        }
        let mut coefficients = vector(self.row_basis.len(), &self.options, budget)?;
        for (coefficient, q) in coefficients.iter_mut().zip(&self.row_basis) {
            *coefficient = dot(q, y, budget)?;
        }
        for (&coefficient, q) in coefficients.iter().zip(&self.row_basis) {
            for (value, &entry) in y.iter_mut().zip(q) {
                budget.step()?;
                *value -= coefficient * entry;
            }
        }
        Ok(())
    }

    #[allow(clippy::type_complexity)]
    fn positive_pairs(
        &self,
        count: usize,
        initial: Option<&[f64]>,
        options: &crate::eigensolver::IterativeOptions,
        budget: &mut WorkBudget<'_>,
    ) -> Result<(
        crate::eigensolver::Eigenpairs,
        Option<usize>,
        f64,
        Option<f64>,
    )> {
        let n = self.basis.len();
        if count > n || count > options.block_size {
            return Err(numerical("positive request exceeds dimension/block"));
        }
        let mut range = Vec::new();
        let mut pivots = Vec::new();
        let mut scale = 1_f64;
        let mut energy = 0_f64;
        // Stream the unsquared factor: Delta*g would suppress small modes twice.
        for j in 0..self.upper_columns {
            budget.step()?;
            let mut y = vector(self.upper_columns, &self.options, budget)?;
            y[j] = 1.;
            self.project_upper(&mut y, budget)?;
            let mut v = vector(n, &self.options, budget)?;
            for (terms, &value) in self.upper.iter().zip(&y) {
                for &(row, coefficient) in terms {
                    budget.step()?;
                    v[row] += coefficient * value;
                }
            }
            let size = length(&v, budget)?;
            energy += size * size;
            append_range(
                v,
                &mut range,
                &mut pivots,
                &mut scale,
                self,
                options,
                budget,
            )?;
        }
        for row in 0..self.lower_rows {
            budget.step()?;
            let mut v = vector(n, &self.options, budget)?;
            for (i, terms) in self.lower.iter().enumerate() {
                for &(index, value) in terms {
                    budget.step()?;
                    if index == row {
                        v[i] = value;
                    }
                }
            }
            let size = length(&v, budget)?;
            energy += size * size;
            append_range(
                v,
                &mut range,
                &mut pivots,
                &mut scale,
                self,
                options,
                budget,
            )?;
        }
        let rank = range.len();
        if count > rank
            || pivots
                .iter()
                .any(|v| *v <= self.options.spectral.relative_tolerance * scale)
        {
            return Err(numerical("insufficient or ambiguous positive range"));
        }
        let threshold = options.absolute_tolerance + options.relative_tolerance * energy;
        if !threshold.is_finite() {
            return Err(numerical("positive range energy"));
        }
        let block = options.block_size.min(rank);
        let mut projected = None;
        if let Some(start) = initial {
            let input_columns = options.block_size.min(n);
            let expected = n.checked_mul(input_columns).ok_or(Error::SizeOverflow {
                operation: "positive initial block",
            })?;
            check_limit(expected, options.max_basis_entries, "max_basis_entries")?;
            if start.len() != expected {
                return Err(Error::ShapeMismatch {
                    input: "positive initial block",
                    expected,
                    actual: start.len(),
                });
            }
            for value in start {
                budget.step()?;
                if !value.is_finite() {
                    return Err(numerical("positive initial block"));
                }
            }
            let mut values = vector(
                rank.checked_mul(block).ok_or(Error::SizeOverflow {
                    operation: "positive initial block",
                })?,
                &self.options,
                budget,
            )?;
            for (i, q) in range.iter().enumerate() {
                for j in 0..block {
                    budget.step_by(n)?;
                    values[i * block + j] =
                        (0..n).map(|r| q[r] * start[r * input_columns + j]).sum();
                }
            }
            projected = Some(values);
        }
        let mut pairs = crate::eigensolver::lobpcg_with_initial(
            rank,
            count,
            false,
            options,
            projected.as_deref(),
            budget,
            |x, b| {
                let v = lift_range(&range, x, n, &self.options, b)?;
                let image = self.action_with(&v, b)?;
                let mut reduced = vector(rank, &self.options, b)?;
                for (value, q) in reduced.iter_mut().zip(&range) {
                    *value = dot(q, &image, b)?;
                }
                Ok(reduced)
            },
        )?;
        if pairs
            .values
            .iter()
            .any(|v| !v.is_finite() || *v <= threshold)
        {
            return Err(numerical("ambiguous positive eigenvalue"));
        }
        let entries = n.checked_mul(count).ok_or(Error::SizeOverflow {
            operation: "positive output vectors",
        })?;
        check_limit(entries, options.max_basis_entries, "max_basis_entries")?;
        let mut vectors = Vec::new();
        vectors
            .try_reserve_exact(entries)
            .map_err(|_| allocation())?;
        vectors.resize(entries, 0.);
        let mut leakage = 0_f64;
        for j in 0..count {
            let mut x = vector(rank, &self.options, budget)?;
            for (i, value) in x.iter_mut().enumerate() {
                budget.step()?;
                *value = pairs.vectors[i * count + j];
            }
            let v = lift_range(&range, &x, n, &self.options, budget)?;
            let image = self.action_with(&v, budget)?;
            pairs.actions = pairs.actions.checked_add(1).ok_or(Error::SizeOverflow {
                operation: "positive action count",
            })?;
            let mut error = vector(n, &self.options, budget)?;
            for i in 0..n {
                budget.step()?;
                vectors[i * count + j] = v[i];
                error[i] = image[i] - pairs.values[j] * v[i];
            }
            pairs.residuals[j] = length(&error, budget)?;
            if pairs.residuals[j]
                > options.absolute_tolerance + options.relative_tolerance * pairs.values[j].abs()
            {
                return Err(numerical("positive original residual"));
            }
            let mut outside = vector(n, &self.options, budget)?;
            budget.step_by(n)?;
            outside.copy_from_slice(&image);
            remove_rows(&mut outside, &range, budget)?;
            leakage = leakage.max(length(&outside, budget)? / length(&image, budget)?.max(1.));
        }
        let mut orthogonality = 0_f64;
        for i in 0..count {
            for j in 0..count {
                budget.step_by(n)?;
                let product: f64 = (0..n)
                    .map(|r| vectors[r * count + i] * vectors[r * count + j])
                    .sum();
                orthogonality = orthogonality.hypot(product - f64::from(i == j));
            }
        }
        if !orthogonality.is_finite()
            || orthogonality > self.options.spectral.residual_tolerance
            || leakage > self.options.spectral.residual_tolerance
        {
            return Err(numerical("positive lifted eigensystem"));
        }
        pairs.vectors = vectors;
        pairs.orthogonality_error = orthogonality;
        budget.check()?;
        Ok((pairs, Some(rank), leakage, Some(threshold)))
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
        for (column, value) in self.upper.iter().zip(&mut y) {
            budget.step()?;
            for &(row, coefficient) in column {
                budget.step()?;
                *value += coefficient * x[row];
            }
        }
        // Compute all coefficients against the original y: I-Q^T Q is a fixed
        // symmetric matrix even with small orthogonality errors in Q.
        self.project_upper(&mut y, budget)?;
        let mut result = vector(n, &self.options, budget)?;
        for (column, &value) in self.upper.iter().zip(&y) {
            budget.step()?;
            for &(row, coefficient) in column {
                budget.step()?;
                result[row] += coefficient * value;
            }
        }
        let mut lower_action = vector(self.lower_rows, &self.options, budget)?;
        for (column, &value) in self.lower.iter().zip(x) {
            budget.step()?;
            for &(row, coefficient) in column {
                budget.step()?;
                lower_action[row] += coefficient * value;
            }
        }
        for (column, target) in self.lower.iter().zip(&mut result) {
            budget.step()?;
            for &(row, coefficient) in column {
                budget.step()?;
                *target += coefficient * lower_action[row];
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

/// Solve a nondecreasing start grid at one fixed endpoint using one workspace.
///
/// Results preserve input order, duplicates and pair metadata. Reverse traversal
/// appends outside-row constraints to one prepared projection. With reuse_vectors,
/// previous vectors are truncated in the same source basis and zero-filled to
/// the requested lookahead block; dependent starts are completed by the shared
/// solver. Every current value/residual is recomputed. This is a warm start,
/// not an extremal-index certificate or a promised iteration improvement.
/// All preparation, updates, solves and output copies share one work budget.
/// The sum of returned vector entries must fit max_basis_entries.
/// # Errors
/// Includes operator/iteration errors, invalid stage grids or per-pair counts,
/// cumulative output limits and failed incremental projection constraints.
#[allow(clippy::too_many_arguments)]
pub fn fixed_endpoint_eigenpairs(
    starts: &[SimplicialStage<'_>],
    endpoint: &SimplicialStage<'_>,
    degree: usize,
    count: usize,
    target: PersistentSpectrumTarget,
    operator_options: OperatorOptions,
    options: &crate::eigensolver::IterativeOptions,
    reuse_vectors: bool,
    execution: &Execution<'_>,
) -> Result<Vec<PersistentEigenpairs>> {
    let mut budget = WorkBudget::new(execution)?;
    validate_pair(
        endpoint,
        endpoint,
        degree,
        LaplacianComponent::Full,
        &operator_options.spectral,
        &mut budget,
    )?;
    crate::eigensolver::lobpcg(0, 0, false, options, &mut budget, |_, _| {
        Err(numerical("unexpected empty batch action"))
    })?;
    let mut total = 0usize;
    for (i, stage) in starts.iter().enumerate() {
        validate_pair(
            stage,
            endpoint,
            degree,
            LaplacianComponent::Full,
            &operator_options.spectral,
            &mut budget,
        )?;
        if i > 0 && starts[i - 1].scale() > stage.scale() {
            return Err(Error::InvalidParameter {
                parameter: "start stages",
                reason: "scales must be nondecreasing",
            });
        }
        let n = basis_ids(stage, degree, &mut budget)?.len();
        if count > n || count > options.block_size {
            return Err(numerical("batch request exceeds dimension/block"));
        }
        total = n
            .checked_mul(count)
            .and_then(|entries| total.checked_add(entries))
            .ok_or(Error::SizeOverflow {
                operation: "batch vector entries",
            })?;
    }
    check_limit(total, options.max_basis_entries, "max_basis_entries")?;
    let mut results: Vec<PersistentEigenpairs> = Vec::new();
    results
        .try_reserve_exact(starts.len())
        .map_err(|_| allocation())?;
    let Some(last) = starts.last() else {
        return Ok(results);
    };
    let mut op = PersistentLaplacianOperator::new_with(
        last,
        endpoint,
        degree,
        operator_options,
        &mut budget,
    )?;
    for stage in starts.iter().rev() {
        op.restrict_with(stage.scale(), &mut budget)?;
        let n = op.basis.len();
        let block = options.block_size.min(n);
        let mut initial = None;
        if reuse_vectors && let Some(previous) = results.last() {
            let entries = n.checked_mul(block).ok_or(Error::SizeOverflow {
                operation: "batch initial block",
            })?;
            check_limit(entries, options.max_basis_entries, "max_basis_entries")?;
            let mut values = Vec::new();
            values
                .try_reserve_exact(entries)
                .map_err(|_| allocation())?;
            values.resize(entries, 0.);
            for i in 0..n {
                for j in 0..count.min(block) {
                    budget.step()?;
                    values[i * block + j] = previous.vectors[i * count + j];
                }
            }
            initial = Some(values);
        }
        results.push(op.eigenpairs_with(
            count,
            target,
            initial.as_deref(),
            options,
            &mut budget,
        )?);
    }
    budget.step_by(results.len())?;
    results.reverse();
    budget.check()?;
    Ok(results)
}

fn add_projection_row(
    mut v: Vec<f64>,
    rows: &mut Vec<Vec<f64>>,
    entries: &mut usize,
    residual: &mut f64,
    options: &OperatorOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<()> {
    let original = length(&v, budget)?;
    for _ in 0..2 {
        remove_rows(&mut v, rows, budget)?;
    }
    let size = length(&v, budget)?;
    if size <= options.spectral.relative_tolerance * original {
        if size > 64. * f64::EPSILON * original.max(1.) {
            return Err(numerical("ambiguous row-space rank"));
        }
        *residual = residual.max(size / original.max(1.));
        return Ok(());
    }
    *entries = entries.checked_add(v.len()).ok_or(Error::SizeOverflow {
        operation: "projection basis entries",
    })?;
    check_limit(
        *entries,
        options.max_projection_entries,
        "max_projection_entries",
    )?;
    for value in &mut v {
        budget.step()?;
        *value /= size;
    }
    for q in rows.iter() {
        *residual = residual.max(dot(&v, q, budget)?.abs());
    }
    *residual = residual.max((dot(&v, &v, budget)? - 1.).abs());
    rows.try_reserve(1).map_err(|_| allocation())?;
    rows.push(v);
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn append_range(
    mut v: Vec<f64>,
    range: &mut Vec<Vec<f64>>,
    pivots: &mut Vec<f64>,
    scale: &mut f64,
    op: &PersistentLaplacianOperator,
    options: &crate::eigensolver::IterativeOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<()> {
    *scale = scale.max(length(&v, budget)?);
    for _ in 0..2 {
        remove_rows(&mut v, range, budget)?;
    }
    let size = length(&v, budget)?;
    if size <= 64. * f64::EPSILON * *scale {
        return Ok(());
    }
    if size <= op.options.spectral.relative_tolerance * *scale {
        return Err(numerical("ambiguous boundary-image rank"));
    }
    let entries = range
        .len()
        .checked_add(1)
        .and_then(|r| r.checked_mul(v.len()))
        .ok_or(Error::SizeOverflow {
            operation: "positive range entries",
        })?;
    check_limit(entries, options.max_basis_entries, "max_basis_entries")?;
    for value in &mut v {
        budget.step()?;
        *value /= size;
    }
    for q in range.iter() {
        if dot(q, &v, budget)?.abs() > op.options.spectral.residual_tolerance {
            return Err(numerical("positive range orthogonality"));
        }
    }
    pivots.try_reserve(1).map_err(|_| allocation())?;
    range.try_reserve(1).map_err(|_| allocation())?;
    pivots.push(size);
    range.push(v);
    Ok(())
}

fn lift_range(
    range: &[Vec<f64>],
    x: &[f64],
    n: usize,
    options: &OperatorOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let mut v = vector(n, options, budget)?;
    for (q, &coefficient) in range.iter().zip(x) {
        for (value, &entry) in v.iter_mut().zip(q) {
            budget.step()?;
            *value += coefficient * entry;
        }
    }
    Ok(v)
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
fn incidence(count: usize, budget: &mut WorkBudget<'_>) -> Result<Vec<Vec<(usize, f64)>>> {
    budget.check()?;
    let mut columns = Vec::new();
    columns.try_reserve_exact(count).map_err(|_| allocation())?;
    for _ in 0..count {
        budget.step()?;
        columns.push(Vec::new());
    }
    Ok(columns)
}
fn push_term(
    column: &mut Vec<(usize, f64)>,
    term: (usize, f64),
    entries: &mut usize,
    options: &OperatorOptions,
) -> Result<()> {
    *entries = entries.checked_add(1).ok_or(Error::SizeOverflow {
        operation: "operator incidence entries",
    })?;
    check_limit(
        *entries,
        options.max_incidence_entries,
        "max_incidence_entries",
    )?;
    column.try_reserve(1).map_err(|_| allocation())?;
    column.push(term);
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
    if norm.is_finite() {
        Ok(norm)
    } else {
        Err(numerical("projection norm"))
    }
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

    #[test]
    fn materialization_shares_one_action_budget() {
        use crate::complex::SimplicialComplex;
        let cells = [
            (vec![0], 0.),
            (vec![2], 0.),
            (vec![1], 1.),
            (vec![0, 1], 2.),
            (vec![1, 2], 2.),
        ];
        let source = SimplicialComplex::new(
            cells
                .into_iter()
                .map(|(v, t)| Simplex::new(v, t).unwrap())
                .collect(),
        )
        .unwrap();
        let op = PersistentLaplacianOperator::new(
            &source.stage(0.).unwrap(),
            &source.stage(2.).unwrap(),
            0,
            OperatorOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        let execution = Execution::default().max_work(10_000);
        let mut budget = WorkBudget::new(&execution).unwrap();
        op.action_with(&[1., 0.], &mut budget).unwrap();
        let one_action = budget.used();
        assert!(
            op.action(&[1., 0.], &Execution::default().max_work(one_action))
                .is_ok()
        );
        assert!(matches!(
            op.materialize(4, &Execution::default().max_work(one_action)),
            Err(Error::WorkLimitExceeded { .. })
        ));
    }
}
