use super::*;
use std::ops::Range;

/// Explicit numeric and resource policy for the small dense real solver.
///
/// This reference path uses scaled cyclic Jacobi rotations, O(n²) memory and
/// O(max_sweeps * n³) work. It adds no runtime dependency and is not an iterative
/// sparse eigensolver. Defaults limit n to 128 and sweeps to 64. Tolerances bound
/// the off-diagonal Frobenius norm by `absolute + relative * ||L||_F`.
#[derive(Clone, Copy, Debug)]
pub struct SpectrumOptions {
    /// Maximum permitted domain dimension, checked before dense allocation.
    pub max_dimension: usize,
    /// Maximum complete cyclic sweeps; zero accepts already diagonal matrices.
    pub max_sweeps: usize,
    /// Finite nonnegative absolute convergence tolerance, in eigenvalue units.
    pub absolute_tolerance: f64,
    /// Finite nonnegative relative convergence tolerance.
    pub relative_tolerance: f64,
}
impl Default for SpectrumOptions {
    fn default() -> Self {
        Self {
            max_dimension: 128,
            max_sweeps: 64,
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
        }
    }
}

/// Owned selected eigenpairs of an ordinary real symmetric coordinate matrix.
///
/// Values are sorted algebraically with multiplicity; vector column j belongs
/// to value j and is Euclidean-normalized. Repeated eigenspaces may rotate and
/// signs are arbitrary. The solver computes the full dense decomposition before
/// selecting the requested index range; a partial result never certifies nullity.
/// This result retains the actual matrix, basis handles, degree and summands so
/// it remains a coordinate problem after source release. Handles are owner-local,
/// not cross-source names; retain source descriptions for geometric interpretation.
pub struct HodgeSpectrum<Id> {
    pub(super) degree: usize,
    pub(super) selection: HodgeSelection,
    pub(super) basis: Vec<Id>,
    matrix: Vec<f64>,
    range: Range<usize>,
    values: Vec<f64>,
    vectors: Vec<f64>,
    residuals: Vec<f64>,
    relative_residuals: Vec<Option<f64>>,
    orthogonality_error: f64,
    matrix_norm: f64,
    sweeps: usize,
    options: SpectrumOptions,
}
impl<Id> HodgeSpectrum<Id> {
    /// Chain degree.
    pub fn degree(&self) -> usize {
        self.degree
    }
    /// Summands used for the retained matrix and eigenpairs.
    pub fn selection(&self) -> HodgeSelection {
        self.selection
    }
    /// Owned ordered source-local handles for the entire domain.
    pub fn basis(&self) -> &[Id] {
        &self.basis
    }
    /// Actual n by n row-major standard-inner-product matrix.
    pub fn matrix(&self) -> &[f64] {
        &self.matrix
    }
    /// Requested and returned algebraic index range within the complete spectrum.
    pub fn index_range(&self) -> Range<usize> {
        self.range.clone()
    }
    /// Whether all n eigenpairs are retained, including an explicitly empty domain.
    pub fn is_complete(&self) -> bool {
        self.range == (0..self.basis.len())
    }
    /// Sorted eigenvalues, preserving multiplicity without clamping small negatives.
    pub fn values(&self) -> &[f64] {
        &self.values
    }
    /// n by r row-major buffer of paired vector columns; r = values().len().
    pub fn vectors(&self) -> &[f64] {
        &self.vectors
    }
    /// Absolute Euclidean residual norms ||L v - lambda v||, using sparse action.
    pub fn residuals(&self) -> &[f64] {
        &self.residuals
    }
    /// Residual divided by (||L||_F + |lambda|) ||v||; None if denominator is zero.
    pub fn relative_residuals(&self) -> &[Option<f64>] {
        &self.relative_residuals
    }
    /// Frobenius norm of V^T V - I for the returned columns.
    pub fn orthogonality_error(&self) -> f64 {
        self.orthogonality_error
    }
    /// Actual Frobenius norm of the original retained matrix, not a norm estimate.
    pub fn matrix_norm(&self) -> f64 {
        self.matrix_norm
    }
    /// Completed cyclic Jacobi sweeps, before algebraic selection.
    pub fn sweeps(&self) -> usize {
        self.sweeps
    }
    /// Policy used for this successful computation.
    pub fn options(&self) -> SpectrumOptions {
        self.options
    }
    /// Solver identity; successful construction means its off-diagonal test passed.
    pub fn method(&self) -> &'static str {
        "scaled cyclic Jacobi, real f64, ordinary symmetric eigenproblem"
    }

    /// Count |lambda| <= threshold only for a complete spectrum.
    ///
    /// For Full this is numerical harmonic nullity, not an exact Betti certificate.
    /// For Lower/Upper it counts the kernel of that summand only. Threshold units
    /// are eigenvalue units and are independent of convergence tolerance.
    /// # Errors
    /// Rejects partial spectra and nonfinite/negative thresholds.
    pub fn numerical_nullity(&self, threshold: f64) -> Result<usize> {
        Ok(self.numerical_kernel_columns(threshold)?.count())
    }

    /// Borrow the vector-column indexes selected by |lambda| <= threshold.
    ///
    /// Requires the complete spectrum. For Full these columns span the numerical
    /// harmonic space; for Lower/Upper they describe only that summand's kernel.
    /// Column j is paired with `values()[j]` and `residuals()[j]`; coordinate i is
    /// `vectors()[i * values().len() + j]`. No vector copy or solve is performed.
    /// Zero classification does not certify exact nullity, even with small residuals.
    /// # Errors
    /// Rejects partial spectra and nonfinite/negative thresholds.
    pub fn numerical_kernel_columns(
        &self,
        threshold: f64,
    ) -> Result<impl Iterator<Item = usize> + '_> {
        tolerance(threshold, "Hodge zero threshold")?;
        if !self.is_complete() {
            return Err(Error::InvalidParameter {
                parameter: "Hodge nullity",
                reason: "requires all eigenpairs",
            });
        }
        Ok(self
            .values
            .iter()
            .enumerate()
            .filter_map(move |(column, value)| (value.abs() <= threshold).then_some(column)))
    }
}

impl<Id: Copy + Eq + Hash> HodgeLaplacian<Id> {
    /// Compute selected eigenpairs with one matrix/solve/verification budget.
    ///
    /// Range indexes the full algebraically sorted spectrum and may be empty.
    /// An empty selection on a nonempty source is still partial. Eigenvectors
    /// always accompany values; no fake empty or zero spectrum represents failure.
    /// The result retains dense matrix storage even for a partial selection.
    /// # Errors
    /// Rejects invalid range/tolerances, dimensions exceeding max_dimension,
    /// allocation failure, arithmetic failure, failed convergence, cancellation
    /// and work exhaustion. All returned residuals/orthogonality are measured.
    pub fn spectrum(
        &self,
        range: Range<usize>,
        options: &SpectrumOptions,
        execution: &Execution<'_>,
    ) -> Result<HodgeSpectrum<Id>> {
        let mut budget = WorkBudget::new(execution)?;
        let n = self.basis.len();
        if range.start > range.end || range.end > n {
            return Err(Error::InvalidParameter {
                parameter: "Hodge eigenpair range",
                reason: "must be ordered within 0..domain size",
            });
        }
        tolerance(options.absolute_tolerance, "Hodge absolute tolerance")?;
        tolerance(options.relative_tolerance, "Hodge relative tolerance")?;
        if n > options.max_dimension {
            return Err(Error::InvalidParameter {
                parameter: "Hodge dense dimension",
                reason: "exceeds the explicit dense solver limit",
            });
        }
        let matrix = self.matrix_budget(&mut budget)?;
        budget.step_by(matrix.len())?;
        let matrix_norm = norm(&matrix);
        finite(matrix_norm)?;
        let (diagonal, eigenvectors, sweeps) = jacobi(&matrix, n, options, &mut budget)?;
        let mut order = Vec::new();
        order.try_reserve_exact(n).map_err(|_| allocation())?;
        budget.step_by(n)?;
        order.extend(0..n);
        order.sort_unstable_by(|&i, &j| diagonal[i].total_cmp(&diagonal[j]).then(i.cmp(&j)));
        budget.check()?;
        let r = range.len();
        let mut values = Vec::new();
        let mut vectors = zeros(n.checked_mul(r).ok_or(size_error())?, &mut budget)?;
        for (column, &index) in order[range.clone()].iter().enumerate() {
            budget.step()?;
            push(&mut values, diagonal[index])?;
            for row in 0..n {
                budget.step()?;
                vectors[row * r + column] = eigenvectors[row * n + index];
            }
        }
        let mut residuals = Vec::new();
        let mut relative_residuals = Vec::new();
        let mut v = zeros(n, &mut budget)?;
        for column in 0..r {
            for row in 0..n {
                budget.step()?;
                v[row] = vectors[row * r + column];
            }
            let mut residual = self.action_budget(&v, &mut budget)?;
            for row in 0..n {
                budget.step()?;
                residual[row] -= values[column] * v[row];
                finite(residual[row])?;
            }
            let residual = norm(&residual);
            finite(residual)?;
            let denominator = (matrix_norm + values[column].abs()) * norm(&v);
            finite(denominator)?;
            push(&mut residuals, residual)?;
            push(
                &mut relative_residuals,
                if denominator == 0.0 {
                    None
                } else {
                    Some(residual / denominator)
                },
            )?;
        }
        let mut orthogonality_error: f64 = 0.0;
        for i in 0..r {
            for j in 0..r {
                let mut dot = 0.0;
                for row in 0..n {
                    budget.step()?;
                    dot += vectors[row * r + i] * vectors[row * r + j];
                }
                orthogonality_error =
                    orthogonality_error.hypot(dot - if i == j { 1.0 } else { 0.0 });
            }
        }
        finite(orthogonality_error)?;
        let mut basis = Vec::new();
        basis.try_reserve_exact(n).map_err(|_| allocation())?;
        budget.step_by(n)?;
        basis.extend_from_slice(&self.basis);
        budget.check()?;
        Ok(HodgeSpectrum {
            degree: self.degree,
            selection: self.selection,
            basis,
            matrix,
            range,
            values,
            vectors,
            residuals,
            relative_residuals,
            orthogonality_error,
            matrix_norm,
            sweeps,
            options: *options,
        })
    }
}

fn tolerance(value: f64, parameter: &'static str) -> Result<()> {
    if !value.is_finite() || value < 0.0 {
        Err(Error::InvalidParameter {
            parameter,
            reason: "must be finite and nonnegative",
        })
    } else {
        Ok(())
    }
}
fn norm(values: &[f64]) -> f64 {
    values.iter().fold(0.0, |sum, &x| sum.hypot(x))
}

fn jacobi(
    matrix: &[f64],
    n: usize,
    options: &SpectrumOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<(Vec<f64>, Vec<f64>, usize)> {
    let mut a = zeros(matrix.len(), budget)?;
    let mut q = zeros(matrix.len(), budget)?;
    budget.step_by(matrix.len())?;
    let scale = matrix
        .iter()
        .fold(0.0_f64, |max, value| max.max(value.abs()));
    for row in 0..n {
        budget.step()?;
        q[row * n + row] = 1.0;
        for column in 0..n {
            budget.step()?;
            if matrix[row * n + column] != matrix[column * n + row] {
                return Err(Error::NumericalFailure {
                    context: "Hodge matrix symmetry",
                });
            }
            a[row * n + column] = if scale == 0.0 {
                0.0
            } else {
                matrix[row * n + column] / scale
            };
        }
    }
    let target = if scale == 0.0 {
        0.0
    } else {
        budget.step_by(a.len())?;
        options.absolute_tolerance / scale + options.relative_tolerance * norm(&a)
    };
    finite(target)?;
    let mut sweeps = 0;
    loop {
        let mut off: f64 = 0.0;
        for i in 0..n {
            for j in 0..n {
                budget.step()?;
                if i != j {
                    off = off.hypot(a[i * n + j]);
                }
            }
        }
        if off <= target {
            break;
        }
        if sweeps == options.max_sweeps {
            return Err(Error::NumericalFailure {
                context: "Hodge Jacobi convergence",
            });
        }
        for p in 0..n {
            for s in p + 1..n {
                budget.step()?;
                let aps = a[p * n + s];
                if aps == 0.0 {
                    continue;
                }
                let delta = (a[s * n + s] - a[p * n + p]) * 0.5;
                let t = if delta == 0.0 {
                    1.0
                } else {
                    aps / (delta + delta.signum() * delta.hypot(aps))
                };
                let c = 1.0 / 1.0_f64.hypot(t);
                let sn = t * c;
                a[p * n + p] -= t * aps;
                a[s * n + s] += t * aps;
                a[p * n + s] = 0.0;
                a[s * n + p] = 0.0;
                for i in 0..n {
                    budget.step()?;
                    if i != p && i != s {
                        let aip = a[i * n + p];
                        let ais = a[i * n + s];
                        a[i * n + p] = c * aip - sn * ais;
                        a[p * n + i] = a[i * n + p];
                        a[i * n + s] = sn * aip + c * ais;
                        a[s * n + i] = a[i * n + s];
                    }
                    let qip = q[i * n + p];
                    let qis = q[i * n + s];
                    q[i * n + p] = c * qip - sn * qis;
                    q[i * n + s] = sn * qip + c * qis;
                }
            }
        }
        sweeps += 1;
    }
    let mut diagonal = Vec::new();
    for i in 0..n {
        budget.step()?;
        let value = a[i * n + i] * scale;
        finite(value)?;
        push(&mut diagonal, crate::canonical_zero(value))?;
    }
    budget.check()?;
    Ok((diagonal, q, sweeps))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicBool, Ordering};

    #[test]
    fn cancellation_after_real_jacobi_work_returns_no_eigenpairs() {
        let flag = AtomicBool::new(false);
        let execution = Execution::default().max_work(1000).cancellation(&flag);
        let mut budget = WorkBudget::new(&execution).unwrap();
        budget.cancel_at_work(75);
        let result = jacobi(
            &[2., 1., -1., 1., 2., 1., -1., 1., 2.],
            3,
            &SpectrumOptions::default(),
            &mut budget,
        );
        assert!(matches!(result, Err(Error::Cancelled)));
        assert!(flag.load(Ordering::Relaxed));
        assert!(budget.used() >= 75);
    }

    fn eigenpair_errors(matrix: &[f64], values: &[f64], vectors: &[f64], n: usize) -> (f64, f64) {
        let mut residual: f64 = 0.;
        let mut orthogonality: f64 = 0.;
        for j in 0..n {
            let mut r: f64 = 0.;
            for i in 0..n {
                let lv: f64 = (0..n).map(|k| matrix[i * n + k] * vectors[k * n + j]).sum();
                r = r.hypot(lv - values[j] * vectors[i * n + j]);
            }
            residual = residual.max(r);
            for k in 0..n {
                let dot: f64 = (0..n)
                    .map(|i| vectors[i * n + j] * vectors[i * n + k])
                    .sum();
                orthogonality = orthogonality.hypot(dot - if j == k { 1. } else { 0. });
            }
        }
        (residual, orthogonality)
    }

    #[test]
    fn jacobi_cluster_small_positive_and_global_scale_cases() {
        // H is the fixed orthogonal Walsh-Hadamard transform. Eigenvalues are
        // the declared diagonal by similarity, independently of Jacobi.
        let h = [
            [1., 1., 1., 1.],
            [1., -1., 1., -1.],
            [1., 1., -1., -1.],
            [1., -1., -1., 1.],
        ];
        for (name, expected, scale, relative) in [
            ("cluster", [1., 1. + 1e-10, 2., 5.], 1., 1e-14),
            ("small_positive", [0., 1e-14, 1., 2.], 1., 1e-16),
            ("small_scale", [0., 1., 2., 5.], 1e-180, 1e-14),
            ("large_scale", [0., 1., 2., 5.], 1e180, 1e-14),
        ] {
            let mut matrix = vec![0.; 16];
            for i in 0..4 {
                for j in 0..4 {
                    matrix[i * 4 + j] = (0..4)
                        .map(|k| h[i][k] * expected[k] * h[j][k] / 4.)
                        .sum::<f64>()
                        * scale;
                }
            }
            let options = SpectrumOptions {
                absolute_tolerance: 0.,
                relative_tolerance: relative,
                ..Default::default()
            };
            let mut budget = WorkBudget::new(&Execution::default()).unwrap();
            let (values, vectors, sweeps) = jacobi(&matrix, 4, &options, &mut budget).unwrap();
            let normalized: Vec<_> = values.iter().map(|v| v / scale).collect();
            let mut sorted = normalized.clone();
            sorted.sort_by(f64::total_cmp);
            for (actual, expected) in sorted.iter().zip(expected) {
                assert!((actual - expected).abs() < 2e-15, "{name}: {sorted:?}");
            }
            if name == "small_positive" {
                assert_eq!(sorted.iter().filter(|v| v.abs() <= 1e-15).count(), 1);
            }
            let normalized_matrix: Vec<_> = matrix.iter().map(|v| v / scale).collect();
            let (residual, orthogonality) =
                eigenpair_errors(&normalized_matrix, &normalized, &vectors, 4);
            assert!(
                residual < 1e-12 && orthogonality < 1e-12,
                "{name}: r={residual},o={orthogonality}"
            );
            println!(
                "Hodge backend {name}: normalized_values={sorted:?}, normalized_residual={residual}, orthogonality={orthogonality}, sweeps={sweeps}"
            );
        }
    }

    #[test]
    fn convergence_and_small_residual_do_not_certify_nullity() {
        // Exact rank is three: the lower block has identical rows and the other
        // two diagonal entries are positive. Kernel = span(1,1,0,0).
        let matrix = [
            5e-13, -5e-13, 0., 0., -5e-13, 5e-13, 0., 0., 0., 0., 1e12, 0., 0., 0., 0., 2e12,
        ];
        let mut budget = WorkBudget::new(&Execution::default()).unwrap();
        let (values, vectors, sweeps) =
            jacobi(&matrix, 4, &SpectrumOptions::default(), &mut budget).unwrap();
        assert_eq!(sweeps, 0);
        let (residual, orthogonality) = eigenpair_errors(&matrix, &values, &vectors, 4);
        assert!(residual < 1e-12);
        assert_eq!(orthogonality, 0.);
        assert_eq!(values.iter().filter(|v| v.abs() <= 1e-14).count(), 0);
        // A global stopping tolerance swamps the small block. A stricter declared
        // policy resolves it; the default result is never an exact certificate.
        let strict = SpectrumOptions {
            absolute_tolerance: 0.,
            relative_tolerance: 1e-30,
            ..Default::default()
        };
        let (values, _, _) = jacobi(&matrix, 4, &strict, &mut budget).unwrap();
        assert_eq!(values.iter().filter(|v| v.abs() <= 1e-14).count(), 1);
        println!(
            "Hodge dynamic-range counterexample: default_converged=true, residual={residual}, threshold_nullity=0, exact_nullity=1; strict_nullity=1"
        );
    }
}
