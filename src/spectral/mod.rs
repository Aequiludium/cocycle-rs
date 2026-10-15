//! Real persistent Laplacians of same-source simplicial stage pairs.
//!
//! Standard inner products and ordinary, nonaugmented chains only. This bounded
//! dense reproduction exposes full, upper and lower operators, Algorithm 2 and
//! fixed-endpoint reverse
//! elimination of Mémoli--Wan--Wang, arXiv:2012.02808v3, Sections 4 and 5.2.
//! It does not implement weighted Grams, general maps or online endpoint updates.
//! Signed incidence is assembled directly; no dense boundary matrix is required.
//! Dense Gram/Schur buffers have explicit limits before allocation. A large
//! source is rejected rather than silently densified beyond those limits.

mod schur;

use crate::complex::{Simplex, SimplexId};
use crate::execution::{Execution, WorkBudget};
use crate::filtration::{Coverage, FiltrationContext, SimplicialStage};
use crate::{Error, Result};
use schur::{copy_matrix, dense, eigensystem, norm, schur_complement};

/// Mathematical component of an operator and its numerical spectrum.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum LaplacianComponent {
    /// Persistent upper term plus the source stage's lower term.
    Full,
    /// Persistent upper term A - E C+ E^T for K included in L.
    Up,
    /// Ordinary lower term Bq(K)^T Bq(K), independent of L.
    Down,
}

/// Numerical and storage policy for the bounded dense reproduction.
#[derive(Clone, Copy, Debug)]
pub struct SpectralOptions {
    /// Relative eigenvalue cutoff for the Schur pseudoinverse (default 1e-12).
    /// Also controls Jacobi convergence and the stepwise zero-pivot test.
    pub relative_tolerance: f64,
    /// Maximum relative eigensystem/range residual (default 1e-9).
    /// Must be at least `relative_tolerance` and less than one.
    pub residual_tolerance: f64,
    /// Maximum cyclic Jacobi sweeps, excluding the initial convergence check.
    pub max_sweeps: usize,
    /// Maximum f64 entries in any one dense buffer (default 65,536).
    /// The fixed-endpoint operation also caps the sum of returned matrix entries.
    /// This is not a process memory quota: temporary buffers coexist.
    pub max_dense_entries: usize,
}
impl Default for SpectralOptions {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-12,
            residual_tolerance: 1e-9,
            max_sweeps: 64,
            max_dense_entries: 65_536,
        }
    }
}
impl SpectralOptions {
    fn validate(&self) -> Result<()> {
        if !self.relative_tolerance.is_finite()
            || !self.residual_tolerance.is_finite()
            || self.relative_tolerance <= 0.0
            || self.residual_tolerance < self.relative_tolerance
            || self.residual_tolerance >= 1.0
        {
            return Err(Error::InvalidParameter {
                parameter: "spectral tolerances",
                reason: "require 0 < relative <= residual < 1, both finite",
            });
        }
        Ok(())
    }
}

/// Owned full or component operator on the ordered q-simplices of K.
///
/// Matrix entries and action use f64, real coefficients, standard inner products
/// and increasing-vertex orientation. Bases contain owned simplex labels and
/// birth values, not transferable owner-local handles. Scales and optional source
/// certificates describe this pair; they do not assert cross-source identity.
/// A bare supplied complex is the entire mathematical source. Certified sources
/// require the q+1 skeleton for full/upper results and the q-skeleton for lower
/// results, unless their construction exhausted the source.
/// All results are independent of the input lifetime.
#[derive(Clone, Debug)]
pub struct PersistentLaplacian {
    component: LaplacianComponent,
    degree: usize,
    scales: (f64, f64),
    basis: Vec<Simplex>,
    target_basis: Vec<Simplex>,
    source_context: Option<FiltrationContext>,
    source_coverage: Coverage,
    construction_dimension: Option<usize>,
    matrix: Vec<f64>,
    options: SpectralOptions,
    solver_residual: f64,
    reference_scale: f64,
}
impl PersistentLaplacian {
    /// Construct the actual pair-level Schur operator A - E C+ E^T + Bq^T Bq.
    ///
    /// Uses one budget across validation, incidence, Jacobi and Schur work.
    /// Work units charge traversals, dense row operations and Jacobi rotations;
    /// they are not timings and need not match the stepwise algorithm.
    /// # Errors
    /// Rejects illegal same-source inclusions, missing certified skeletons,
    /// invalid policy, size/allocation limits, cancellation, work exhaustion,
    /// nonconvergence, negative Schur eigenvalues or excessive range residuals.
    pub fn compute(
        k: &SimplicialStage<'_>,
        l: &SimplicialStage<'_>,
        degree: usize,
        options: SpectralOptions,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        Self::compute_component(k, l, degree, LaplacianComponent::Full, options, execution)
    }

    /// Construct only the persistent upper operator A - E C+ E^T.
    ///
    /// Its owned metadata and ordered coordinates describe the pair (K,L).
    /// It does not assemble the lower term or diagonalize the full operator.
    /// # Errors
    /// Has the same inclusion, skeleton, numeric and resource checks as `compute`.
    pub fn compute_up(
        k: &SimplicialStage<'_>,
        l: &SimplicialStage<'_>,
        degree: usize,
        options: SpectralOptions,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        Self::compute_component(k, l, degree, LaplacianComponent::Up, options, execution)
    }

    /// Construct only Bq(K)^T Bq(K), independently of any endpoint L.
    ///
    /// Both scales and both bases describe K: `scales()` is (s,s).
    /// Requires the certified q-skeleton, but no q+1 skeleton or Schur solve.
    /// At q=0 this is the ordinary nonaugmented zero operator.
    /// # Errors
    /// Rejects an insufficient certified q-skeleton, invalid policy, allocation
    /// or storage limits, degree overflow, cancellation and work exhaustion.
    pub fn compute_down(
        k: &SimplicialStage<'_>,
        degree: usize,
        options: SpectralOptions,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        Self::compute_component(k, k, degree, LaplacianComponent::Down, options, execution)
    }

    fn compute_component(
        k: &SimplicialStage<'_>,
        l: &SimplicialStage<'_>,
        degree: usize,
        component: LaplacianComponent,
        options: SpectralOptions,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        validate_pair(k, l, degree, component, &options, &mut budget)?;
        let target_ids = basis_ids(l, degree, &mut budget)?;
        let count = target_ids
            .iter()
            .take_while(|id| id.index() < k.len())
            .count();
        let (matrix, residual, reference_scale) = if component == LaplacianComponent::Down {
            (dense(count, count, &options, &mut budget)?, 0.0, 0.0)
        } else {
            let up = up_matrix(l, &target_ids, degree, &options, &mut budget)?;
            let reference_scale = norm(&up);
            let (matrix, residual) =
                schur_complement(&up, target_ids.len(), count, &options, &mut budget)?;
            (matrix, residual, reference_scale)
        };
        finish(
            k,
            l,
            degree,
            &target_ids[..count],
            &target_ids,
            matrix,
            residual,
            options,
            component,
            reference_scale,
            &mut budget,
        )
    }

    /// Full, persistent upper, or source-stage lower operator.
    pub fn component(&self) -> LaplacianComponent {
        self.component
    }
    /// Degree of the real chain space on which this operator acts.
    pub fn degree(&self) -> usize {
        self.degree
    }
    /// Inclusive source and target scales (s,t).
    pub fn scales(&self) -> (f64, f64) {
        self.scales
    }
    /// Complete ordered q-basis of K, with original vertex labels and births.
    pub fn basis(&self) -> &[Simplex] {
        &self.basis
    }
    /// Complete ordered q-basis of L before elimination.
    pub fn target_basis(&self) -> &[Simplex] {
        &self.target_basis
    }
    /// Owned original construction context, or None for a supplied complex.
    /// Source coverage is limited by its construction cutoff; the computation
    /// was checked for sufficient dimension at construction time.
    pub fn source_context(&self) -> Option<&FiltrationContext> {
        self.source_context.as_ref()
    }
    /// Original source coverage, independent of the selected pair scales.
    /// A bare supplied complex is complete as its own mathematical source.
    pub fn source_coverage(&self) -> Coverage {
        self.source_coverage
    }
    /// Certified construction dimension, or None for bare supplied topology.
    /// The operator only promises degree q and was checked for needed cofaces;
    /// this value alone is not a full-source spectral completeness certificate.
    pub fn construction_dimension(&self) -> Option<usize> {
        self.construction_dimension
    }
    /// Full owned matrix in row-major basis order, including an empty 0x0 matrix.
    pub fn matrix(&self) -> &[f64] {
        &self.matrix
    }
    /// Applied numerical/storage policy. No numerical rank is an exact certificate.
    pub fn options(&self) -> SpectralOptions {
        self.options
    }
    /// Maximum normalized Schur eigensystem/discarded-range, solve and symmetry
    /// residual. Solve residual uses the max entry of E^T as denominator;
    /// symmetry uses the original upper Gram scale. A zero denominator keeps
    /// the absolute error. For the fixed-endpoint algorithm this records the
    /// discarded-row residual.
    /// For a lower-only operator this is zero because no Schur solve is used;
    /// its eigensystem residual is checked separately by `spectrum`.
    pub fn solver_residual(&self) -> f64 {
        self.solver_residual
    }

    /// Max-entry scale used for numerical PSD acceptance and near-zero guidance.
    ///
    /// This is the maximum of the assembled matrix scale and the original target
    /// upper Gram scale (the latter is absent for lower-only construction).
    /// Retaining the construction scale prevents near-zero Schur cancellation
    /// errors from being judged relative only to their own tiny magnitude.
    /// It is neither a forward-error bound nor an exact rank certificate.
    pub fn reference_scale(&self) -> f64 {
        self.reference_scale
    }

    /// Apply this pair-level operator to finite coordinates in `basis()` order.
    /// # Errors
    /// Rejects wrong shape, nonfinite coordinates, numerical overflow,
    /// allocation failure, cancellation and work exhaustion.
    pub fn action(&self, x: &[f64], execution: &Execution<'_>) -> Result<Vec<f64>> {
        let mut budget = WorkBudget::new(execution)?;
        let n = self.basis.len();
        if x.len() != n {
            return Err(Error::ShapeMismatch {
                input: "Laplacian coordinates",
                expected: n,
                actual: x.len(),
            });
        }
        for (index, value) in x.iter().enumerate() {
            budget.step()?;
            if !value.is_finite() {
                return Err(Error::NonFiniteValue {
                    field: "Laplacian coordinates",
                    index: Some(index),
                });
            }
        }
        let mut result = dense(n, 1, &self.options, &mut budget)?;
        for (i, row) in self.matrix.chunks(n.max(1)).enumerate() {
            budget.step_by(n)?;
            result[i] = row.iter().zip(x).map(|(a, b)| a * b).sum();
            if !result[i].is_finite() {
                return Err(numerical("Laplacian action"));
            }
        }
        budget.check()?;
        Ok(result)
    }

    /// Compute this component's complete sorted spectrum, retaining multiplicity
    /// and small signed roundoff values. The residual is checked without clipping values.
    /// A near-zero count is not a real Betti certificate. Even exact upper/lower
    /// nullities individually are not the persistent Betti number of the pair.
    /// # Errors
    /// Includes resource, nonconvergence, residual and significant negativity errors.
    pub fn spectrum(&self, execution: &Execution<'_>) -> Result<LaplacianSpectrum> {
        let mut budget = WorkBudget::new(execution)?;
        let (mut eigenvalues, _, residual) =
            eigensystem(&self.matrix, self.basis.len(), &self.options, &mut budget)?;
        let scale = self.reference_scale;
        if eigenvalues
            .iter()
            .any(|v| *v < -self.options.residual_tolerance * scale)
        {
            return Err(numerical(
                "persistent Laplacian is not numerically positive semidefinite",
            ));
        }
        budget.check()?;
        eigenvalues.sort_by(f64::total_cmp);
        budget.check()?;
        Ok(LaplacianSpectrum {
            eigenvalues,
            residual,
            zero_threshold: self.options.relative_tolerance * scale,
            component: self.component,
        })
    }
}

/// Complete numerical spectrum; near-zero multiplicity is not exact homology.
#[derive(Clone, Debug)]
pub struct LaplacianSpectrum {
    /// Component whose spectrum was computed; this is not a pair/source identity.
    /// Keep the originating operator's degree, bases and source metadata when
    /// storing spectra from different requests. No thresholded spectral union
    /// or exact rank inference is performed implicitly.
    pub component: LaplacianComponent,
    /// All eigenvalues in increasing order, with multiplicity and no clipping.
    pub eigenvalues: Vec<f64>,
    /// Relative max-entry residual of the computed eigensystem.
    pub residual: f64,
    /// Suggested absolute near-zero threshold using the operator's reference scale.
    /// In particular, upper Schur cancellation is judged at the input Gram scale.
    pub zero_threshold: f64,
}

/// Fixed-endpoint reverse elimination (Section 5.2, equation 15).
///
/// `starts` must be in nondecreasing scale order and included in `target`.
/// Builds the target up-Laplacian once, eliminates q-coordinates in reverse
/// filtration order, then adds each start's own lower term. Returns owned
/// operators in the input order, including repeated and empty stages.
/// Tied births are eliminated together before a requested stage is returned.
/// This is an explicit alternative, not an automatically selected fast path.
/// It does not update the endpoint online. Spectrum computation is separate.
/// # Errors
/// Includes `PersistentLaplacian::compute` errors, unsorted starts, cumulative
/// result storage limits and incompatible near-zero pivot rows.
pub fn fixed_endpoint_laplacians(
    starts: &[SimplicialStage<'_>],
    target: &SimplicialStage<'_>,
    degree: usize,
    options: SpectralOptions,
    execution: &Execution<'_>,
) -> Result<Vec<PersistentLaplacian>> {
    fixed_endpoint_components(
        starts,
        target,
        degree,
        LaplacianComponent::Full,
        options,
        execution,
    )
}

/// Fixed-endpoint reverse elimination returning only persistent upper operators.
///
/// Like `fixed_endpoint_laplacians`, but without constructing any lower term.
/// Results retain pair metadata, input order and `LaplacianComponent::Up`.
/// # Errors
/// Has the same validation, storage, numeric and execution checks as the full route.
pub fn fixed_endpoint_up_laplacians(
    starts: &[SimplicialStage<'_>],
    target: &SimplicialStage<'_>,
    degree: usize,
    options: SpectralOptions,
    execution: &Execution<'_>,
) -> Result<Vec<PersistentLaplacian>> {
    fixed_endpoint_components(
        starts,
        target,
        degree,
        LaplacianComponent::Up,
        options,
        execution,
    )
}

fn fixed_endpoint_components(
    starts: &[SimplicialStage<'_>],
    target: &SimplicialStage<'_>,
    degree: usize,
    component: LaplacianComponent,
    options: SpectralOptions,
    execution: &Execution<'_>,
) -> Result<Vec<PersistentLaplacian>> {
    let mut budget = WorkBudget::new(execution)?;
    validate_pair(target, target, degree, component, &options, &mut budget)?;
    let mut total_entries = 0usize;
    for (index, k) in starts.iter().enumerate() {
        validate_pair(k, target, degree, component, &options, &mut budget)?;
        if index > 0 && starts[index - 1].scale() > k.scale() {
            return Err(Error::InvalidParameter {
                parameter: "start stages",
                reason: "scales must be nondecreasing",
            });
        }
        let mut n = 0usize;
        for simplex in k.simplices() {
            budget.step()?;
            if simplex.dimension() == degree {
                n += 1;
            }
        }
        total_entries = total_entries
            .checked_add(n.checked_mul(n).ok_or(Error::SizeOverflow {
                operation: "fixed-endpoint result size",
            })?)
            .ok_or(Error::SizeOverflow {
                operation: "fixed-endpoint result size",
            })?;
    }
    schur::check_entries(total_entries, &options)?;
    let mut result = Vec::new();
    result
        .try_reserve_exact(starts.len())
        .map_err(|_| allocation())?;
    if starts.is_empty() {
        return Ok(result);
    }
    let ids = basis_ids(target, degree, &mut budget)?;
    let order = ids.len();
    let mut up = up_matrix(target, &ids, degree, &options, &mut budget)?;
    let scale = norm(&up);
    let mut active = order;
    let mut residual: f64 = 0.0;
    for k in starts.iter().rev() {
        let count = ids.iter().take_while(|id| id.index() < k.len()).count();
        while active > count {
            budget.step()?;
            let last = active - 1;
            let pivot = up[last * order + last];
            if pivot < -options.residual_tolerance * scale {
                return Err(numerical("negative stepwise pivot"));
            }
            if pivot.abs() <= options.relative_tolerance * scale {
                let row = (0..last)
                    .map(|i| up[i * order + last].abs())
                    .fold(0.0, f64::max);
                residual = residual.max(if scale == 0.0 { row } else { row / scale });
                if residual > options.residual_tolerance {
                    return Err(numerical("incompatible zero stepwise pivot"));
                }
            } else {
                for i in 0..last {
                    budget.step_by(last)?;
                    for j in 0..=i {
                        let value =
                            up[i * order + j] - up[i * order + last] * up[last * order + j] / pivot;
                        if !value.is_finite() {
                            return Err(numerical("stepwise Schur subtraction"));
                        }
                        up[i * order + j] = value;
                        up[j * order + i] = value;
                    }
                }
            }
            active -= 1;
        }
        let matrix = copy_matrix(&up, order, count, &options, &mut budget)?;
        result.push(finish(
            k,
            target,
            degree,
            &ids[..count],
            &ids,
            matrix,
            residual,
            options,
            component,
            scale,
            &mut budget,
        )?);
    }
    result.reverse();
    budget.check()?;
    Ok(result)
}

fn validate_pair(
    k: &SimplicialStage<'_>,
    l: &SimplicialStage<'_>,
    degree: usize,
    component: LaplacianComponent,
    options: &SpectralOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<()> {
    budget.step()?;
    options.validate()?;
    let adjacent = degree.checked_add(1).ok_or(Error::SizeOverflow {
        operation: "adjacent spectral degree",
    })?;
    let _ = k.inclusion_into(l)?;
    if let Some(source) = l.source_filtration()
        && !source.is_dimension_complete()
        && source.max_simplex_dimension()
            < if component == LaplacianComponent::Down {
                degree
            } else {
                adjacent
            }
    {
        return Err(Error::InsufficientSkeleton {
            requested_homology_dimension: degree,
            constructed_simplex_dimension: source.max_simplex_dimension(),
        });
    }
    Ok(())
}
fn basis_ids(
    stage: &SimplicialStage<'_>,
    degree: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<SimplexId>> {
    let mut ids = Vec::new();
    for (index, simplex) in stage.simplices().iter().enumerate() {
        budget.step()?;
        if simplex.dimension() == degree {
            ids.try_reserve(1).map_err(|_| allocation())?;
            ids.push(SimplexId(index));
        }
    }
    Ok(ids)
}
fn up_matrix(
    stage: &SimplicialStage<'_>,
    ids: &[SimplexId],
    degree: usize,
    options: &SpectralOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let n = ids.len();
    let mut matrix = dense(n, n, options, budget)?;
    for (index, simplex) in stage.simplices().iter().enumerate() {
        budget.step()?;
        if simplex.dimension() != degree + 1 {
            continue;
        }
        let terms = stage
            .boundary(SimplexId(index))
            .ok_or(Error::InternalInvariant {
                reason: "stage boundary missing",
            })?;
        for a in terms {
            let i = ids
                .binary_search(&a.face)
                .map_err(|_| Error::InternalInvariant {
                    reason: "q-face absent",
                })?;
            for b in terms {
                budget.step()?;
                let j = ids
                    .binary_search(&b.face)
                    .map_err(|_| Error::InternalInvariant {
                        reason: "q-face absent",
                    })?;
                matrix[i * n + j] += f64::from(a.coefficient) * f64::from(b.coefficient);
            }
        }
    }
    Ok(matrix)
}
#[allow(clippy::too_many_arguments)]
fn finish(
    k: &SimplicialStage<'_>,
    l: &SimplicialStage<'_>,
    degree: usize,
    ids: &[SimplexId],
    target_ids: &[SimplexId],
    mut matrix: Vec<f64>,
    residual: f64,
    options: SpectralOptions,
    component: LaplacianComponent,
    reference_scale: f64,
    budget: &mut WorkBudget<'_>,
) -> Result<PersistentLaplacian> {
    let n = ids.len();
    if component != LaplacianComponent::Up {
        for (i, &a) in ids.iter().enumerate() {
            for (j, &b) in ids[..=i].iter().enumerate() {
                budget.step()?;
                let ba = k.boundary(a).ok_or(Error::InternalInvariant {
                    reason: "K boundary missing",
                })?;
                let bb = k.boundary(b).ok_or(Error::InternalInvariant {
                    reason: "K boundary missing",
                })?;
                let mut dot = 0.0;
                for x in ba {
                    for y in bb {
                        budget.step()?;
                        if x.face == y.face {
                            dot += f64::from(x.coefficient) * f64::from(y.coefficient);
                        }
                    }
                }
                matrix[i * n + j] += dot;
                matrix[j * n + i] = matrix[i * n + j];
            }
        }
    }
    if matrix.iter().any(|v| !v.is_finite()) {
        return Err(numerical("persistent Laplacian assembly"));
    }
    let owned_basis = |stage: &SimplicialStage<'_>,
                       handles: &[SimplexId],
                       budget: &mut WorkBudget<'_>|
     -> Result<Vec<Simplex>> {
        let mut basis = Vec::new();
        basis
            .try_reserve_exact(handles.len())
            .map_err(|_| allocation())?;
        for &id in handles {
            budget.step()?;
            basis.push(
                stage
                    .simplex(id)
                    .ok_or(Error::InternalInvariant {
                        reason: "basis simplex missing",
                    })?
                    .clone(),
            );
        }
        Ok(basis)
    };
    let basis = owned_basis(k, ids, budget)?;
    let target_basis = owned_basis(l, target_ids, budget)?;
    let reference_scale = reference_scale.max(norm(&matrix));
    budget.check()?;
    Ok(PersistentLaplacian {
        component,
        degree,
        scales: (k.scale(), l.scale()),
        basis,
        target_basis,
        source_context: l.source_filtration().map(|s| s.context().clone()),
        source_coverage: l
            .source_filtration()
            .map_or(Coverage::Complete, |s| s.coverage()),
        construction_dimension: l.source_filtration().map(|s| s.max_simplex_dimension()),
        matrix,
        options,
        solver_residual: residual,
        reference_scale,
    })
}
fn allocation() -> Error {
    Error::AllocationFailed {
        context: "persistent Laplacian",
    }
}
fn numerical(context: &'static str) -> Error {
    Error::NumericalFailure { context }
}
