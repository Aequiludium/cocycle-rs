//! Bounded, unoptimized real-coordinate rewrite of the #145 paper controls.
//! Both Grams and complete restricted decompositions are deliberately retained.
use super::schur::dense;
use super::{SpectralOptions, numerical};
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result};

/// Both complete Gram spectra and a numerical original-domain reconstruction.
///
/// Coordinates are those of the supplied real matrix A, not inferred simplex
/// labels or a source/pair certificate. Both square Grams are constructed even
/// when one is smaller. This is a bounded reference, not automatic dispatch.
#[derive(Clone, Debug)]
pub struct GramSpectrum {
    /// Number of rows of A.
    pub rows: usize,
    /// Original domain dimension, the number of columns of A.
    pub columns: usize,
    /// Whether rows < columns; advisory shape metadata, not automatic dispatch.
    pub dual_is_smaller: bool,
    /// Row-major A^T A, columns by columns.
    pub primal: Vec<f64>,
    /// Row-major A A^T, rows by rows.
    pub dual: Vec<f64>,
    /// Complete ascending primal spectrum, including signed roundoff.
    pub primal_spectrum: Vec<f64>,
    /// Complete ascending dual spectrum, including signed roundoff.
    pub dual_spectrum: Vec<f64>,
    /// Dual values strictly above the numerical rank threshold, with multiplicity.
    pub positive: Vec<f64>,
    /// Original-domain reconstruction, with columns-rank numerical zeros.
    /// This does not replace either unfiltered complete spectrum.
    pub reconstructed: Vec<f64>,
    /// Original columns by numerical_rank row-major lifted dual vectors.
    pub lifted_vectors: Vec<f64>,
    /// Threshold rank, never an exact rational rank certificate.
    pub numerical_rank: usize,
    /// relative_tolerance times the largest absolute eigenvalue of either Gram.
    pub threshold: f64,
    /// A positive raw primal value lies within ten times the rank cutoff.
    /// Roundoff zero signs can affect this advisory flag.
    pub rank_ambiguous: bool,
    /// Original-coordinate Frobenius residual divided by max(1,||V Lambda||F).
    pub lifted_residual: f64,
    /// Frobenius norm ||V^T V-I|| for the lifted positive vectors.
    pub orthogonality_error: f64,
}

/// Compute both bounded Grams and lift positive dual modes by A^T z/sqrt(lambda).
///
/// A is finite real f64, row-major, with standard inner products. Empty matrix
/// dimensions are explicit. Numerical zeros are used only in `reconstructed`;
/// both complete raw spectra retain signs and multiplicities. The rank cutoff
/// is numerical, and tiny lifted modes can fail the residual/orthogonality guard.
/// One budget covers validation, products, both eigensolves and fresh residuals.
/// # Errors
/// Rejects invalid options, shape/size overflow, nonfinite input/arithmetic,
/// per-buffer storage limits, cancellation, exhaustion, eigensolver failure,
/// negative Grams, inconsistent numerical ranks or inaccurate lifted vectors.
pub fn gram_spectrum(
    a: &[f64],
    rows: usize,
    columns: usize,
    options: SpectralOptions,
    execution: &Execution<'_>,
) -> Result<GramSpectrum> {
    let mut budget = WorkBudget::new(execution)?;
    options.validate()?;
    gram(a, rows, columns, &options, &mut budget)
}

fn gram(
    a: &[f64],
    rows: usize,
    columns: usize,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<GramSpectrum> {
    validate(a, rows, columns, o, b)?;
    let at = transpose(a, rows, columns, o, b)?;
    let primal = product(&at, a, columns, rows, columns, o, b)?;
    let dual = product(a, &at, rows, columns, rows, o, b)?;
    let (lp, _) = solve(&primal, columns, o, b)?;
    let (ld, zd) = solve(&dual, rows, o, b)?;
    let scale = lp.iter().chain(&ld).fold(0_f64, |x, v| x.max(v.abs()));
    let threshold = o.relative_tolerance * scale;
    if lp.iter().chain(&ld).any(|v| *v < -threshold) {
        return Err(numerical("negative Gram spectrum"));
    }
    let rank = ld.iter().filter(|v| **v > threshold).count();
    if rank > columns || lp.iter().filter(|v| **v > threshold).count() != rank {
        return Err(numerical("inconsistent Gram numerical ranks"));
    }
    let mut positive = dense(rank, 1, o, b)?;
    let mut lifted = dense(columns, rank, o, b)?;
    let mut k = 0;
    for (j, &value) in ld.iter().enumerate() {
        if value <= threshold {
            continue;
        }
        positive[k] = value;
        for i in 0..columns {
            b.step_by(rows)?;
            let v: f64 = (0..rows)
                .map(|r| a[r * columns + i] * zd[r * rows + j])
                .sum();
            lifted[i * rank + k] = finite(v / value.sqrt())?;
        }
        k += 1;
    }
    let residual = residual(&primal, columns, &lifted, &positive, o, b)?;
    let orthogonality = orthogonality(&lifted, columns, rank, o, b)?;
    if residual > o.residual_tolerance || orthogonality > o.residual_tolerance {
        return Err(numerical("Gram lifted eigensystem"));
    }
    let mut reconstructed = dense(columns, 1, o, b)?;
    reconstructed[columns - rank..].copy_from_slice(&positive);
    b.step_by(columns)?;
    let ambiguous = lp.iter().any(|v| *v > 0. && *v <= 10. * threshold);
    b.check()?;
    Ok(GramSpectrum {
        rows,
        columns,
        dual_is_smaller: rows < columns,
        primal,
        dual,
        primal_spectrum: lp,
        dual_spectrum: ld,
        positive,
        reconstructed,
        lifted_vectors: lifted,
        numerical_rank: rank,
        threshold,
        rank_ambiguous: ambiguous,
        lifted_residual: residual,
        orthogonality_error: orthogonality,
    })
}

/// Explicit project extension merging positive component multisets, not zeros.
#[derive(Clone, Debug)]
pub struct SplitSpectrum {
    /// Lower boundary Gram computation on the original n-column domain.
    pub lower: GramSpectrum,
    /// Gram computation of F^T, with original n-column domain.
    pub upper: GramSpectrum,
    /// Full n by n row-major F F^T + A^T A.
    pub matrix: Vec<f64>,
    /// Exactly n values: n-rank(A)-rank(F) numerical zeros and both positives.
    pub eigenvalues: Vec<f64>,
    /// Either component reports numerical rank ambiguity.
    pub rank_ambiguous: bool,
    /// ||A F||F/max(1,||A||F ||F||F), checked before recombination.
    pub chain_residual: f64,
}

/// Recombine the spectra of A^T A and F F^T when A F is numerically zero.
///
/// A is lower_rows by n and F is n by factor_columns, both row-major. F must
/// already include the orthonormal persistent-boundary metric. No kernel basis,
/// source inclusion or topology is inferred from these raw coordinates.
/// Both Gram computations remain unoptimized and share the operation budget.
/// # Errors
/// Includes `gram_spectrum` errors, incompatible chain products and rank sums.
pub fn split_spectrum(
    a: &[f64],
    lower_rows: usize,
    n: usize,
    f: &[f64],
    factor_columns: usize,
    options: SpectralOptions,
    execution: &Execution<'_>,
) -> Result<SplitSpectrum> {
    let mut b = WorkBudget::new(execution)?;
    options.validate()?;
    validate(a, lower_rows, n, &options, &mut b)?;
    validate(f, n, factor_columns, &options, &mut b)?;
    let af = product(a, f, lower_rows, n, factor_columns, &options, &mut b)?;
    let denominator = finite(length(a, &mut b)? * length(f, &mut b)?)?.max(1.);
    let chain_residual = finite(length(&af, &mut b)? / denominator)?;
    if chain_residual > options.residual_tolerance {
        return Err(numerical("split spectrum chain identity"));
    }
    let lower = gram(a, lower_rows, n, &options, &mut b)?;
    let ft = transpose(f, n, factor_columns, &options, &mut b)?;
    let upper = gram(&ft, factor_columns, n, &options, &mut b)?;
    let rank = lower
        .numerical_rank
        .checked_add(upper.numerical_rank)
        .ok_or(overflow())?;
    if rank > n {
        return Err(numerical("split spectrum rank sum"));
    }
    let mut eigenvalues = dense(n, 1, &options, &mut b)?;
    for (i, v) in lower.positive.iter().chain(&upper.positive).enumerate() {
        eigenvalues[n - rank + i] = *v;
    }
    b.step_by(n)?;
    eigenvalues.sort_by(f64::total_cmp);
    let mut matrix = dense(n, n, &options, &mut b)?;
    for (i, value) in matrix.iter_mut().enumerate() {
        b.step()?;
        *value = finite(lower.primal[i] + upper.primal[i])?;
    }
    b.check()?;
    Ok(SplitSpectrum {
        rank_ambiguous: lower.rank_ambiguous || upper.rank_ambiguous,
        lower,
        upper,
        matrix,
        eigenvalues,
        chain_residual,
    })
}

/// Owned bounded restriction to the complement of a supplied harmonic kernel.
///
/// Uses PETLS Section 4.5 with a verified numerical kernel, not arbitrary cycle
/// representatives. Retains the original real coordinate matrix; callers retain
/// their labels, degree and source/pair metadata. No source identity is inferred.
/// The complete reduced spectrum is retained, including tiny positive values.
#[derive(Clone, Debug)]
pub struct HarmonicRestriction {
    n: usize,
    k: usize,
    original: Vec<f64>,
    kernel: Vec<f64>,
    complement: Vec<f64>,
    matrix: Vec<f64>,
    values: Vec<f64>,
    vectors: Vec<f64>,
    scale: f64,
    threshold: f64,
    ambiguous: bool,
    leakage: f64,
    residual: f64,
    orthogonality: f64,
    options: SpectralOptions,
}

impl HarmonicRestriction {
    /// Normalize a complete kernel N, form X, and solve X^T H X completely.
    ///
    /// H is an exactly symmetric n by n row-major finite matrix. N is n by k,
    /// row-major and full rank; standard inner products only. Numerical kernel
    /// membership is checked relative to ||H||2. Independence uses one-sided
    /// Jacobi singular values, without squaring N's condition number. X is an
    /// orthonormal complement built from coordinate vectors. All returned data
    /// are owned; input buffers may be released. Each dense buffer obeys the cap.
    /// # Errors
    /// Rejects shape/options, nonfinite arithmetic, asymmetric H, dependent or
    /// nonharmonic N, a nonpositive restricted spectrum, inaccurate eigenvectors,
    /// storage limits, nonconvergence, cancellation and work exhaustion. Tiny
    /// positive reduced values remain present but set `rank_ambiguous()`.
    pub fn compute(
        h: &[f64],
        n: usize,
        kernel: &[f64],
        k: usize,
        options: SpectralOptions,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        let mut b = WorkBudget::new(execution)?;
        options.validate()?;
        validate(h, n, n, &options, &mut b)?;
        validate(kernel, n, k, &options, &mut b)?;
        if k > n {
            return Err(invalid("harmonic kernel dimension"));
        }
        let (_, _, scale) = symmetric_scale(h, n, &options, &mut b)?;
        let q = kernel_basis(kernel, n, k, &options, &mut b)?;
        let hq = product(h, &q, n, n, k, &options, &mut b)?;
        let leakage = finite(length(&hq, &mut b)? / scale.max(f64::MIN_POSITIVE))?;
        if leakage > options.residual_tolerance {
            return Err(numerical("supplied representatives are not harmonic"));
        }
        let r = n - k;
        let mut x = dense(n, r, &options, &mut b)?;
        let mut used = 0;
        for coordinate in 0..n {
            if used == r {
                break;
            }
            let mut v = dense(n, 1, &options, &mut b)?;
            v[coordinate] = 1.;
            for _ in 0..2 {
                subtract_columns(&mut v, &q, n, k, k, &mut b)?;
                subtract_columns(&mut v, &x, n, r, used, &mut b)?;
            }
            let norm = length(&v, &mut b)?;
            if norm <= options.relative_tolerance {
                continue;
            }
            for i in 0..n {
                b.step()?;
                x[i * r + used] = v[i] / norm;
            }
            used += 1;
        }
        if used != r {
            return Err(numerical("harmonic complement rank"));
        }
        let hx = product(h, &x, n, n, r, &options, &mut b)?;
        let xt = transpose(&x, n, r, &options, &mut b)?;
        let mut reduced = product(&xt, &hx, r, n, r, &options, &mut b)?;
        symmetrize(&mut reduced, r, &mut b)?;
        let (values, rotations) = solve(&reduced, r, &options, &mut b)?;
        if values.first().is_some_and(|v| *v <= 0.) {
            return Err(numerical("incomplete kernel or nonpositive restriction"));
        }
        let vectors = product(&x, &rotations, n, r, r, &options, &mut b)?;
        let residual = residual(h, n, &vectors, &values, &options, &mut b)?;
        let orthogonality = orthogonality(&vectors, n, r, &options, &mut b)?;
        if residual > options.residual_tolerance || orthogonality > options.residual_tolerance {
            return Err(numerical("harmonic restriction eigensystem"));
        }
        let original = copy(h, n, n, &options, &mut b)?;
        let threshold = options.relative_tolerance * scale;
        b.check()?;
        Ok(Self {
            n,
            k,
            original,
            kernel: q,
            complement: x,
            matrix: reduced,
            ambiguous: values.iter().any(|v| *v <= 10. * threshold),
            values,
            vectors,
            scale,
            threshold,
            leakage,
            residual,
            orthogonality,
            options,
        })
    }
    /// Original coordinate dimension n.
    pub fn dimension(&self) -> usize {
        self.n
    }
    /// Supplied kernel dimension k; not a certified Betti number.
    pub fn kernel_dimension(&self) -> usize {
        self.k
    }
    /// Original n by n matrix, preserving input coordinate order.
    pub fn original_matrix(&self) -> &[f64] {
        &self.original
    }
    /// Orthonormal n by k row-major kernel basis.
    pub fn kernel(&self) -> &[f64] {
        &self.kernel
    }
    /// Orthonormal n by (n-k) row-major complement X.
    pub fn complement(&self) -> &[f64] {
        &self.complement
    }
    /// Complete (n-k) by (n-k) row-major reduced matrix.
    pub fn matrix(&self) -> &[f64] {
        &self.matrix
    }
    /// Complete sorted reduced spectrum; tiny positive values are not clipped.
    pub fn eigenvalues(&self) -> &[f64] {
        &self.values
    }
    /// Complete n by (n-k) row-major lifted eigenvectors.
    pub fn vectors(&self) -> &[f64] {
        &self.vectors
    }
    /// Relative cutoff times ||H||2, in eigenvalue units.
    pub fn zero_threshold(&self) -> f64 {
        self.threshold
    }
    /// A reduced value is within ten times the cutoff; positive iteration rejects it.
    pub fn rank_ambiguous(&self) -> bool {
        self.ambiguous
    }
    /// ||H N||F/max(||H||2, MIN_POSITIVE) for the normalized supplied kernel.
    pub fn kernel_residual(&self) -> f64 {
        self.leakage
    }
    /// Fresh original-coordinate complete eigensystem residual.
    pub fn residual(&self) -> f64 {
        self.residual
    }
    /// Frobenius norm ||V^T V-I|| for the complete lifted eigenvectors.
    pub fn orthogonality_error(&self) -> f64 {
        self.orthogonality
    }

    /// Basic LOBPCG Algorithms 1/2 on the bounded positive restriction.
    ///
    /// Initial vectors are n by block_size row-major in ORIGINAL coordinates;
    /// they are projected with X^T. Require count <= block_size <= n-k. This is
    /// the checked paper-basic algorithm with B=I, identity preconditioning,
    /// diagonal metric scaling and Cholesky, without SVD search compression.
    /// Residual stopping is relative_tolerance * max(1,||X^T H X||2).
    /// Returned vectors are lifted and freshly checked against the original H
    /// under the SAME terminal budget. This does not certify smallest indexes.
    /// A zero count is rejected, matching the frozen paper-basic transcription;
    /// shape, finite input and execution controls are still validated.
    /// # Errors
    /// Includes invalid count/tolerance/shape, ambiguous positive rank, dependent
    /// initial/search blocks, metric condition > 1e12, failed Cholesky, iteration
    /// or work exhaustion, cancellation, storage limits and fresh residual errors.
    pub fn smallest_eigenpairs(
        &self,
        initial: &[f64],
        block_size: usize,
        count: usize,
        max_iterations: usize,
        relative_tolerance: f64,
        execution: &Execution<'_>,
    ) -> Result<PositiveEigenpairs> {
        let mut b = WorkBudget::new(execution)?;
        let r = self.n - self.k;
        validate(initial, self.n, block_size, &self.options, &mut b)?;
        if count == 0
            || count > block_size
            || block_size > r
            || (count > 0 && block_size == 0)
            || !relative_tolerance.is_finite()
            || relative_tolerance <= 0.
        {
            return Err(invalid("positive iteration count/block/tolerance"));
        }
        if count > 0 && self.ambiguous {
            return Err(numerical("ambiguous positive restriction"));
        }
        let xt = transpose(&self.complement, self.n, r, &self.options, &mut b)?;
        let start = product(&xt, initial, r, self.n, block_size, &self.options, &mut b)?;
        let (values, reduced_vectors, iterations) = basic_lobpcg(
            &self.matrix,
            r,
            &start,
            block_size,
            count,
            max_iterations,
            relative_tolerance,
            &self.options,
            &mut b,
        )?;
        let vectors = product(
            &self.complement,
            &reduced_vectors,
            self.n,
            r,
            count,
            &self.options,
            &mut b,
        )?;
        let mut residuals = dense(count, 1, &self.options, &mut b)?;
        let images = product(
            &self.original,
            &vectors,
            self.n,
            self.n,
            count,
            &self.options,
            &mut b,
        )?;
        for j in 0..count {
            let mut error = 0_f64;
            for i in 0..self.n {
                b.step()?;
                error = error.hypot(images[i * count + j] - vectors[i * count + j] * values[j]);
            }
            residuals[j] = finite(error / self.scale.max(1.))?;
        }
        let qt = transpose(&self.kernel, self.n, self.k, &self.options, &mut b)?;
        let leaks = product(&qt, &vectors, self.k, self.n, count, &self.options, &mut b)?;
        let leakage = length(&leaks, &mut b)?;
        let orthogonality = orthogonality(&vectors, self.n, count, &self.options, &mut b)?;
        if residuals.iter().any(|v| *v > relative_tolerance)
            || leakage > self.options.residual_tolerance
            || orthogonality > self.options.residual_tolerance
        {
            return Err(numerical("positive original-coordinate eigensystem"));
        }
        b.check()?;
        Ok(PositiveEigenpairs {
            original_dimension: self.n,
            requested_count: count,
            values,
            vectors,
            residuals,
            range_leakage: leakage,
            orthogonality_error: orthogonality,
            iterations,
            zero_threshold: self.threshold,
            rank_ambiguous: self.ambiguous,
        })
    }
}

/// Owned approximate positive modes in the original matrix coordinates.
#[derive(Clone, Debug)]
pub struct PositiveEigenpairs {
    /// Original coordinate dimension, before kernel removal.
    pub original_dimension: usize,
    /// Validated requested count, equal to values.len(); never padded.
    pub requested_count: usize,
    /// Ascending returned Ritz values with their retained multiplicities.
    pub values: Vec<f64>,
    /// Original_dimension by requested_count row-major vector columns.
    pub vectors: Vec<f64>,
    /// Fresh ||Hv-lambda v||2/max(1,||H||2) for each returned original vector.
    pub residuals: Vec<f64>,
    /// ||N^T V||F with N the orthonormal harmonic basis.
    pub range_leakage: f64,
    /// ||V^T V-I||F for returned vectors.
    pub orthogonality_error: f64,
    /// Updates after initial Rayleigh--Ritz; no heuristic iteration guarantee.
    pub iterations: usize,
    /// Numerical positivity threshold, not an exact rank certificate.
    pub zero_threshold: f64,
    /// False on successful positive requests; ambiguous prepared domains are errors.
    pub rank_ambiguous: bool,
}
impl PositiveEigenpairs {
    /// Always false: residuals do not establish extremal-index coverage.
    pub fn is_index_coverage_certified(&self) -> bool {
        false
    }
}

fn validate(
    a: &[f64],
    m: usize,
    n: usize,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<()> {
    let entries = m.checked_mul(n).ok_or(overflow())?;
    super::schur::check_entries(entries, o)?;
    if a.len() != entries {
        return Err(Error::ShapeMismatch {
            input: "split spectral matrix",
            expected: entries,
            actual: a.len(),
        });
    }
    for (index, v) in a.iter().enumerate() {
        b.step()?;
        if !v.is_finite() {
            return Err(Error::NonFiniteValue {
                field: "split spectral matrix",
                index: Some(index),
            });
        }
    }
    b.check()
}
fn copy(
    a: &[f64],
    m: usize,
    n: usize,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let mut out = dense(m, n, o, b)?;
    b.step_by(a.len())?;
    out.copy_from_slice(a);
    Ok(out)
}
fn transpose(
    a: &[f64],
    m: usize,
    n: usize,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let mut out = dense(n, m, o, b)?;
    if m == 0 || n == 0 {
        return Ok(out);
    }
    for i in 0..m {
        b.step_by(n)?;
        for j in 0..n {
            out[j * m + i] = a[i * n + j];
        }
    }
    Ok(out)
}
#[allow(clippy::too_many_arguments)]
fn product(
    a: &[f64],
    c: &[f64],
    m: usize,
    inner: usize,
    n: usize,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let mut out = dense(m, n, o, b)?;
    if m == 0 || n == 0 {
        return Ok(out);
    }
    for i in 0..m {
        for j in 0..n {
            b.step_by(inner)?;
            out[i * n + j] = finite((0..inner).map(|k| a[i * inner + k] * c[k * n + j]).sum())?;
        }
    }
    Ok(out)
}
fn solve(
    a: &[f64],
    n: usize,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<(Vec<f64>, Vec<f64>)> {
    super::schur::check_entries(n.checked_mul(n).ok_or(overflow())?, o)?;
    let (values, vectors, _) =
        crate::eigensolver::jacobi(a, n, 0., 4. * f64::EPSILON, o.max_sweeps, b)?;
    let mut order = Vec::new();
    order
        .try_reserve_exact(n)
        .map_err(|_| super::allocation())?;
    order.extend(0..n);
    b.step_by(n)?;
    order.sort_by(|&i, &j| values[i].total_cmp(&values[j]));
    let mut sorted = dense(n, 1, o, b)?;
    let mut q = dense(n, n, o, b)?;
    for (j, &index) in order.iter().enumerate() {
        sorted[j] = values[index];
        for i in 0..n {
            b.step()?;
            q[i * n + j] = vectors[i * n + index];
        }
    }
    Ok((sorted, q))
}
fn symmetric_scale(
    a: &[f64],
    n: usize,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<(Vec<f64>, Vec<f64>, f64)> {
    let (v, q) = solve(a, n, o, b)?;
    let scale = v.iter().fold(0_f64, |a, v| a.max(v.abs()));
    Ok((v, q, scale))
}
fn symmetrize(a: &mut [f64], n: usize, b: &mut WorkBudget<'_>) -> Result<()> {
    for i in 0..n {
        for j in 0..i {
            b.step()?;
            let value = finite(a[i * n + j] / 2. + a[j * n + i] / 2.)?;
            a[i * n + j] = value;
            a[j * n + i] = value;
        }
    }
    Ok(())
}
fn residual(
    h: &[f64],
    n: usize,
    v: &[f64],
    values: &[f64],
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<f64> {
    let r = values.len();
    let hv = product(h, v, n, n, r, o, b)?;
    let (mut error, mut scale) = (0_f64, 0_f64);
    for i in 0..n {
        for j in 0..r {
            b.step()?;
            let expected = finite(v[i * r + j] * values[j])?;
            error = error.hypot(hv[i * r + j] - expected);
            scale = scale.hypot(expected);
        }
    }
    finite(error / scale.max(1.))
}
fn orthogonality(
    v: &[f64],
    n: usize,
    r: usize,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<f64> {
    let vt = transpose(v, n, r, o, b)?;
    let mut metric = product(&vt, v, r, n, r, o, b)?;
    for i in 0..r {
        b.step()?;
        metric[i * r + i] -= 1.;
    }
    length(&metric, b)
}
fn length(a: &[f64], b: &mut WorkBudget<'_>) -> Result<f64> {
    b.step_by(a.len())?;
    finite(a.iter().fold(0_f64, |a, v| a.hypot(*v)))
}
fn finite(v: f64) -> Result<f64> {
    if v.is_finite() {
        Ok(v)
    } else {
        Err(numerical("split spectral arithmetic"))
    }
}
fn invalid(reason: &'static str) -> Error {
    Error::InvalidParameter {
        parameter: "split spectral input",
        reason,
    }
}
fn overflow() -> Error {
    Error::SizeOverflow {
        operation: "split spectral shape",
    }
}

fn subtract_columns(
    v: &mut [f64],
    q: &[f64],
    n: usize,
    stride: usize,
    count: usize,
    b: &mut WorkBudget<'_>,
) -> Result<()> {
    for j in 0..count {
        b.step_by(n)?;
        let dot: f64 = (0..n).map(|i| v[i] * q[i * stride + j]).sum();
        finite(dot)?;
        for i in 0..n {
            b.step()?;
            v[i] -= dot * q[i * stride + j];
        }
    }
    Ok(())
}

fn kernel_basis(
    a: &[f64],
    n: usize,
    k: usize,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let mut q = copy(a, n, k, o, b)?;
    if k == 0 {
        return Ok(q);
    }
    let scale = a.iter().fold(0_f64, |s, v| s.max(v.abs()));
    if scale == 0. {
        return Err(numerical("dependent harmonic basis"));
    }
    for v in &mut q {
        b.step()?;
        *v /= scale;
    }
    let mut converged = false;
    for sweep in 0..=o.max_sweeps {
        let mut changed = false;
        for p in 0..k {
            for j in p + 1..k {
                let (mut aa, mut bb, mut ab) = (0., 0., 0.);
                for i in 0..n {
                    b.step()?;
                    let (x, y) = (q[i * k + p], q[i * k + j]);
                    aa += x * x;
                    bb += y * y;
                    ab += x * y;
                }
                if ab.abs() <= 8. * f64::EPSILON * (aa * bb).sqrt() {
                    continue;
                }
                changed = true;
                if sweep == o.max_sweeps {
                    continue;
                }
                let tau = (bb - aa) / (2. * ab);
                let t = if tau >= 0. {
                    1. / (tau + tau.hypot(1.))
                } else {
                    -1. / (-tau + tau.hypot(1.))
                };
                let c = 1. / t.hypot(1.);
                let s = t * c;
                for i in 0..n {
                    b.step()?;
                    let (x, y) = (q[i * k + p], q[i * k + j]);
                    q[i * k + p] = c * x - s * y;
                    q[i * k + j] = s * x + c * y;
                }
            }
        }
        if !changed {
            converged = true;
            break;
        }
    }
    if !converged {
        return Err(numerical("harmonic basis Jacobi convergence"));
    }
    let mut norms = dense(k, 1, o, b)?;
    for j in 0..k {
        for i in 0..n {
            b.step()?;
            norms[j] = norms[j].hypot(q[i * k + j]);
        }
    }
    let maximum = norms.iter().copied().fold(0_f64, f64::max);
    if norms.iter().any(|v| *v <= o.relative_tolerance * maximum) {
        return Err(numerical("dependent harmonic basis"));
    }
    for i in 0..n {
        for j in 0..k {
            b.step()?;
            q[i * k + j] /= norms[j];
        }
    }
    if orthogonality(&q, n, k, o, b)? > o.residual_tolerance {
        return Err(numerical("harmonic basis orthogonality"));
    }
    Ok(q)
}

// Algorithms 1/2, B=I and K=I. No call to the existing stabilized Hodge/Dirac
// entry, no hidden search-space compression and no fresh nested WorkBudget.
#[allow(clippy::too_many_arguments)]
fn basic_lobpcg(
    h: &[f64],
    n: usize,
    initial: &[f64],
    block: usize,
    count: usize,
    max_iterations: usize,
    tolerance: f64,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<(Vec<f64>, Vec<f64>, usize)> {
    if count == 0 {
        return Err(invalid("basic LOBPCG requires positive count"));
    }
    let (_, _, scale) = symmetric_scale(h, n, o, b)?;
    let limit = finite(tolerance * scale.max(1.))?;
    let (c, mut values) = rayleigh_ritz(h, n, initial, block, block, o, b)?;
    let mut x = product(initial, &c, n, block, block, o, b)?;
    let mut p = Vec::new();
    for iteration in 0..=max_iterations {
        let mut residual = product(h, &x, n, n, block, o, b)?;
        let mut errors = dense(block, 1, o, b)?;
        for i in 0..n {
            for j in 0..block {
                b.step()?;
                residual[i * block + j] =
                    finite(residual[i * block + j] - x[i * block + j] * values[j])?;
                errors[j] = finite(errors[j].hypot(residual[i * block + j]))?;
            }
        }
        if errors[..count].iter().all(|v| *v <= limit) {
            let mut vectors = dense(n, count, o, b)?;
            let mut selected = dense(count, 1, o, b)?;
            selected.copy_from_slice(&values[..count]);
            for i in 0..n {
                b.step_by(count)?;
                vectors[i * count..(i + 1) * count]
                    .copy_from_slice(&x[i * block..i * block + count]);
            }
            return Ok((selected, vectors, iteration));
        }
        if iteration == max_iterations {
            return Err(numerical("basic LOBPCG iteration limit"));
        }
        let previous = if p.is_empty() { 0 } else { block };
        let ns = block
            .checked_mul(2)
            .and_then(|v| v.checked_add(previous))
            .ok_or(overflow())?;
        let mut s = dense(n, ns, o, b)?;
        for i in 0..n {
            b.step_by(ns)?;
            s[i * ns..i * ns + block].copy_from_slice(&x[i * block..(i + 1) * block]);
            s[i * ns + block..i * ns + 2 * block]
                .copy_from_slice(&residual[i * block..(i + 1) * block]);
            if previous > 0 {
                s[i * ns + 2 * block..(i + 1) * ns].copy_from_slice(&p[i * block..(i + 1) * block]);
            }
        }
        let (c, next_values) = rayleigh_ritz(h, n, &s, ns, block, o, b)?;
        x = product(&s, &c, n, ns, block, o, b)?;
        p = dense(n, block, o, b)?;
        for i in 0..n {
            for j in 0..block {
                b.step_by(ns - block)?;
                p[i * block + j] =
                    finite((block..ns).map(|k| s[i * ns + k] * c[k * block + j]).sum())?;
            }
        }
        values = next_values;
    }
    Err(numerical("basic LOBPCG iteration limit"))
}

fn rayleigh_ritz(
    h: &[f64],
    n: usize,
    s: &[f64],
    ns: usize,
    count: usize,
    o: &SpectralOptions,
    b: &mut WorkBudget<'_>,
) -> Result<(Vec<f64>, Vec<f64>)> {
    if ns > n || count > ns {
        return Err(numerical("dependent Rayleigh-Ritz search"));
    }
    let st = transpose(s, n, ns, o, b)?;
    let mut metric = product(&st, s, ns, n, ns, o, b)?;
    let mut d = dense(ns, 1, o, b)?;
    for i in 0..ns {
        b.step()?;
        if metric[i * ns + i] <= 0. {
            return Err(numerical("zero Rayleigh-Ritz search column"));
        }
        d[i] = finite(1. / metric[i * ns + i].sqrt())?;
    }
    for i in 0..ns {
        for j in 0..ns {
            b.step()?;
            metric[i * ns + j] = finite(d[i] * metric[i * ns + j] * d[j])?;
        }
    }
    symmetrize(&mut metric, ns, b)?;
    let (e, _) = solve(&metric, ns, o, b)?;
    if e.first().is_none_or(|v| *v <= 0.) || e[ns - 1] / e[0] > 1e12 {
        return Err(numerical("ill-conditioned Rayleigh-Ritz search"));
    }
    let mut r = dense(ns, ns, o, b)?;
    for i in 0..ns {
        for j in i..ns {
            b.step_by(i)?;
            let value = finite(
                metric[i * ns + j] - (0..i).map(|k| r[k * ns + i] * r[k * ns + j]).sum::<f64>(),
            )?;
            r[i * ns + j] = if i == j {
                if value <= 0. {
                    return Err(numerical("Rayleigh-Ritz Cholesky pivot"));
                }
                value.sqrt()
            } else {
                finite(value / r[i * ns + i])?
            };
        }
    }
    let mut inverse = dense(ns, ns, o, b)?;
    for j in 0..ns {
        for i in (0..ns).rev() {
            b.step_by(ns - i)?;
            let rhs = if i == j { 1. } else { 0. };
            inverse[i * ns + j] = finite(
                (rhs - (i + 1..ns)
                    .map(|k| r[i * ns + k] * inverse[k * ns + j])
                    .sum::<f64>())
                    / r[i * ns + i],
            )?;
        }
    }
    for i in 0..ns {
        for j in 0..ns {
            b.step()?;
            inverse[i * ns + j] = finite(d[i] * inverse[i * ns + j])?;
        }
    }
    let hs = product(h, s, n, n, ns, o, b)?;
    let shs = product(&st, &hs, ns, n, ns, o, b)?;
    let transformed = product(&shs, &inverse, ns, ns, ns, o, b)?;
    let it = transpose(&inverse, ns, ns, o, b)?;
    let mut projected = product(&it, &transformed, ns, ns, ns, o, b)?;
    symmetrize(&mut projected, ns, b)?;
    let (values, z) = solve(&projected, ns, o, b)?;
    let mut selected = dense(ns, count, o, b)?;
    for i in 0..ns {
        b.step_by(count)?;
        selected[i * count..(i + 1) * count].copy_from_slice(&z[i * ns..i * ns + count]);
    }
    let c = product(&inverse, &selected, ns, ns, count, o, b)?;
    let mut v = dense(count, 1, o, b)?;
    v.copy_from_slice(&values[..count]);
    Ok((c, v))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rayleigh_ritz_satisfies_an_independent_two_by_two_pencil() {
        let h = [
            1., 0., 0., 0., 0., 2., 0., 0., 0., 0., 5., 0., 0., 0., 0., 7.,
        ];
        let s = [1., 0., 1., 2., 0., 1., 2., 1.];
        let mut b = WorkBudget::new(&Execution::default()).unwrap();
        let (c, values) =
            rayleigh_ritz(&h, 4, &s, 2, 2, &SpectralOptions::default(), &mut b).unwrap();
        // det(Hp-lambda G)=20 lambda²-162 lambda+296, derived by hand.
        let expected = [(81. - 641_f64.sqrt()) / 20., (81. + 641_f64.sqrt()) / 20.];
        for (a, e) in values.iter().zip(expected) {
            assert!((a - e).abs() < 1e-12);
        }
        for i in 0..2 {
            for j in 0..2 {
                let (mut metric, mut projected) = (0., 0.);
                for row in 0..4 {
                    let x = s[2 * row] * c[i] + s[2 * row + 1] * c[2 + i];
                    let y = s[2 * row] * c[j] + s[2 * row + 1] * c[2 + j];
                    metric += x * y;
                    projected += x * h[5 * row] * y;
                }
                assert!((metric - if i == j { 1. } else { 0. }).abs() < 1e-12);
                assert!((projected - if i == j { values[i] } else { 0. }).abs() < 1e-11);
            }
        }
    }

    #[test]
    fn preparation_and_iteration_observe_cancellation_inside_work() {
        use std::sync::atomic::{AtomicBool, Ordering};
        let flag = AtomicBool::new(false);
        let e = Execution::default().max_work(u64::MAX).cancellation(&flag);
        let o = SpectralOptions::default();
        let mut b = WorkBudget::new(&e).unwrap();
        b.cancel_at_work(15);
        assert!(matches!(
            gram(&[1., 0., 0., 1.], 2, 2, &o, &mut b),
            Err(Error::Cancelled)
        ));
        flag.store(false, Ordering::Relaxed);
        let mut b = WorkBudget::new(&e).unwrap();
        b.cancel_at_work(15);
        assert!(matches!(
            kernel_basis(&[1., 1., 0., 1.], 2, 2, &o, &mut b),
            Err(Error::Cancelled)
        ));
        flag.store(false, Ordering::Relaxed);
        let mut b = WorkBudget::new(&e).unwrap();
        b.cancel_at_work(120);
        assert!(matches!(
            basic_lobpcg(&[1., 0., 0., 2.], 2, &[1., 1.], 1, 1, 10, 1e-9, &o, &mut b),
            Err(Error::Cancelled)
        ));
    }
}
