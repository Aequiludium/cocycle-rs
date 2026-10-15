//! Bounded symmetric Jacobi and generalized Schur complement, independently
//! implemented from the mathematical formulas; no external source is copied.
use super::{SpectralOptions, allocation, numerical};
use crate::execution::WorkBudget;
use crate::{Error, Result};

pub(super) fn check_entries(entries: usize, options: &SpectralOptions) -> Result<()> {
    if entries > options.max_dense_entries {
        return Err(Error::InvalidParameter {
            parameter: "max_dense_entries",
            reason: "dense spectral storage exceeds the declared limit",
        });
    }
    Ok(())
}
pub(super) fn dense(
    rows: usize,
    columns: usize,
    options: &SpectralOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let count = rows.checked_mul(columns).ok_or(Error::SizeOverflow {
        operation: "spectral matrix shape",
    })?;
    check_entries(count, options)?;
    budget.check()?;
    let mut result = Vec::new();
    result.try_reserve_exact(count).map_err(|_| allocation())?;
    result.resize(count, 0.0);
    budget.check()?;
    Ok(result)
}
pub(super) fn norm(matrix: &[f64]) -> f64 {
    matrix.iter().map(|v| v.abs()).fold(0.0, f64::max)
}

pub(super) fn eliminate_last(
    matrix: &mut [f64],
    stride: usize,
    active: usize,
    scale: f64,
    options: &SpectralOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<()> {
    budget.step()?;
    let last = active - 1;
    let pivot = matrix[last * stride + last];
    if pivot < -options.residual_tolerance * scale {
        return Err(numerical("negative stepwise pivot"));
    }
    if pivot.abs() <= options.relative_tolerance * scale {
        // Small coupling is not a bound on c*c^T/pivot. Even a PSD rank-one
        // block [[1, eps], [eps, eps^2]] loses a unit correction if dropped.
        // Only an exactly uncoupled coordinate can be removed without a solve.
        for i in 0..last {
            budget.step()?;
            if matrix[i * stride + last] != 0.0 {
                return Err(numerical("coupled near-zero stepwise pivot"));
            }
        }
        return Ok(());
    }
    for i in 0..last {
        budget.step_by(last)?;
        for j in 0..=i {
            let correction = matrix[i * stride + last] * (matrix[last * stride + j] / pivot);
            let value = matrix[i * stride + j] - correction;
            if !value.is_finite() {
                return Err(numerical("stepwise Schur subtraction"));
            }
            matrix[i * stride + j] = value;
            matrix[j * stride + i] = value;
        }
    }
    Ok(())
}
pub(super) fn copy_matrix(
    matrix: &[f64],
    stride: usize,
    n: usize,
    options: &SpectralOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let mut result = dense(n, n, options, budget)?;
    for i in 0..n {
        budget.step_by(n)?;
        result[i * n..(i + 1) * n].copy_from_slice(&matrix[i * stride..i * stride + n]);
    }
    Ok(result)
}

pub(super) fn eigensystem(
    matrix: &[f64],
    n: usize,
    options: &SpectralOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<(Vec<f64>, Vec<f64>, f64)> {
    let mut a = copy_matrix(matrix, n, n, options, budget)?;
    let mut v = dense(n, n, options, budget)?;
    for i in 0..n {
        v[i * n + i] = 1.0;
    }
    let scale = norm(matrix);
    let threshold = options.relative_tolerance * scale;
    let mut converged = false;
    for sweep in 0..=options.max_sweeps {
        let mut off: f64 = 0.0;
        for p in 0..n {
            budget.step_by(n)?;
            for q in p + 1..n {
                off = off.max(a[p * n + q].abs());
            }
        }
        if off <= threshold {
            converged = true;
            break;
        }
        if sweep == options.max_sweeps {
            break;
        }
        for p in 0..n {
            for q in p + 1..n {
                budget.step()?;
                let apq = a[p * n + q];
                if apq == 0.0 {
                    continue;
                }
                // Stable Jacobi tangent: avoid squaring a potentially huge tau.
                let tau = (a[q * n + q] - a[p * n + p]) / (2.0 * apq);
                let t = if tau >= 0.0 {
                    1.0 / (tau + tau.hypot(1.0))
                } else {
                    -1.0 / (-tau + tau.hypot(1.0))
                };
                let c = 1.0 / (1.0 + t * t).sqrt();
                let s = t * c;
                a[p * n + p] -= t * apq;
                a[q * n + q] += t * apq;
                a[p * n + q] = 0.0;
                a[q * n + p] = 0.0;
                for i in 0..n {
                    budget.step()?;
                    if i != p && i != q {
                        let x = a[i * n + p];
                        let y = a[i * n + q];
                        a[i * n + p] = c * x - s * y;
                        a[p * n + i] = a[i * n + p];
                        a[i * n + q] = s * x + c * y;
                        a[q * n + i] = a[i * n + q];
                    }
                    let x = v[i * n + p];
                    let y = v[i * n + q];
                    v[i * n + p] = c * x - s * y;
                    v[i * n + q] = s * x + c * y;
                }
            }
        }
    }
    if !converged {
        return Err(numerical("spectral Jacobi did not converge"));
    }
    let mut values = dense(n, 1, options, budget)?;
    for i in 0..n {
        values[i] = a[i * n + i];
    }
    let mut error: f64 = 0.0;
    for i in 0..n {
        for j in 0..n {
            budget.step_by(n)?;
            let mut product = 0.0;
            for k in 0..n {
                product += matrix[i * n + k] * v[k * n + j];
            }
            error = error.max((product - v[i * n + j] * values[j]).abs());
        }
    }
    let residual = if scale == 0.0 { error } else { error / scale };
    if !residual.is_finite()
        || residual > options.residual_tolerance
        || values.iter().chain(&v).any(|x| !x.is_finite())
    {
        return Err(numerical("spectral eigensystem residual"));
    }
    budget.check()?;
    Ok((values, v, residual))
}

pub(super) fn schur_complement(
    up: &[f64],
    order: usize,
    keep: usize,
    options: &SpectralOptions,
    budget: &mut WorkBudget<'_>,
) -> Result<(Vec<f64>, f64)> {
    let mut result = copy_matrix(up, order, keep, options, budget)?;
    let m = order - keep;
    if m == 0 || keep == 0 {
        return Ok((result, 0.0));
    }
    let mut c = dense(m, m, options, budget)?;
    for i in 0..m {
        budget.step_by(m)?;
        c[i * m..(i + 1) * m]
            .copy_from_slice(&up[(i + keep) * order + keep..(i + keep) * order + order]);
    }
    let (values, v, mut residual) = eigensystem(&c, m, options, budget)?;
    let scale = norm(&c);
    let up_scale = norm(up);
    let cutoff = options.relative_tolerance * values.iter().map(|x| x.abs()).fold(0.0, f64::max);
    let mut coupling = dense(keep, 1, options, budget)?;
    // Match the reference SchurComp.m structure: solve C X = E^T in the
    // retained eigenspace, then form A - E X. Never materialize C+.
    let mut solution = dense(m, keep, options, budget)?;
    let mut discarded_correction = 0.0;
    for eigen in 0..m {
        let lambda = values[eigen];
        if lambda < -options.residual_tolerance * scale {
            return Err(numerical("negative Schur block eigenvalue"));
        }
        for i in 0..keep {
            budget.step_by(m)?;
            coupling[i] = (0..m)
                .map(|j| up[i * order + keep + j] * v[j * m + eigen])
                .sum();
        }
        if lambda <= cutoff {
            if lambda > 0.0 {
                budget.step_by(keep)?;
                // A small C X - E^T residual relative to another large RHS
                // does not bound a discarded c*c^T/lambda correction. Sum
                // the positive-direction entry bounds before accepting it.
                let scaled = norm(&coupling) / lambda.sqrt();
                let correction = scaled * scaled;
                discarded_correction += if up_scale == 0.0 {
                    correction
                } else {
                    correction / up_scale
                };
                if !discarded_correction.is_finite()
                    || discarded_correction > options.residual_tolerance
                {
                    return Err(numerical(
                        "discarded Schur direction has significant correction",
                    ));
                }
                residual = residual.max(discarded_correction);
            }
            let range = if up_scale == 0.0 {
                norm(&coupling)
            } else {
                norm(&coupling) / up_scale
            };
            residual = residual.max(range);
            if !range.is_finite() || range > options.residual_tolerance {
                return Err(numerical("discarded Schur direction has nonzero coupling"));
            }
            continue;
        }
        for row in 0..m {
            budget.step_by(keep)?;
            for i in 0..keep {
                // Apply the inverse eigenvalue before multiplying the two
                // couplings; squaring their original scale can lose the
                // correction even when X and E X are representable.
                solution[row * keep + i] += v[row * m + eigen] * (coupling[i] / lambda);
            }
        }
    }
    if solution.iter().any(|x| !x.is_finite()) {
        return Err(numerical("Schur minimum-norm solution"));
    }
    let mut solve_error: f64 = 0.0;
    let mut rhs_scale: f64 = 0.0;
    for row in 0..m {
        for i in 0..keep {
            budget.step_by(m)?;
            let rhs = up[(keep + row) * order + i];
            let product: f64 = (0..m)
                .map(|j| c[row * m + j] * solution[j * keep + i])
                .sum();
            let error = (product - rhs).abs();
            if !error.is_finite() {
                return Err(numerical("Schur solve residual"));
            }
            solve_error = solve_error.max(error);
            rhs_scale = rhs_scale.max(rhs.abs());
        }
    }
    let solve_residual = if rhs_scale == 0.0 {
        solve_error
    } else {
        solve_error / rhs_scale
    };
    if solve_residual > options.residual_tolerance {
        return Err(numerical("Schur solve residual exceeds tolerance"));
    }
    residual = residual.max(solve_residual);
    for i in 0..keep {
        for j in 0..keep {
            budget.step_by(m)?;
            let correction: f64 = (0..m)
                .map(|row| up[i * order + keep + row] * solution[row * keep + j])
                .sum();
            result[i * keep + j] -= correction;
        }
    }
    if result.iter().any(|x| !x.is_finite()) {
        return Err(numerical("Schur complement"));
    }
    for i in 0..keep {
        for j in 0..i {
            budget.step()?;
            let difference = (result[i * keep + j] - result[j * keep + i]).abs();
            let symmetry_residual = if up_scale == 0.0 {
                difference
            } else {
                difference / up_scale
            };
            if symmetry_residual > options.residual_tolerance {
                return Err(numerical(
                    "Schur reconstruction is not numerically symmetric",
                ));
            }
            residual = residual.max(symmetry_residual);
            let mean = 0.5 * result[i * keep + j] + 0.5 * result[j * keep + i];
            result[i * keep + j] = mean;
            result[j * keep + i] = mean;
        }
    }
    budget.check()?;
    Ok((result, residual))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::execution::Execution;
    use std::sync::atomic::AtomicBool;
    #[test]
    fn jacobi_can_cancel_after_actual_work() {
        let flag = AtomicBool::new(false);
        let execution = Execution::default().max_work(100).cancellation(&flag);
        let mut budget = WorkBudget::new(&execution).unwrap();
        budget.cancel_at_work(7);
        assert!(matches!(
            eigensystem(
                &[2.0, -1.0, -1.0, 2.0],
                2,
                &SpectralOptions::default(),
                &mut budget
            ),
            Err(Error::Cancelled)
        ));
    }
    #[test]
    fn scaled_singular_blocks_and_incompatible_cutoff() {
        for scale in [1e-120, 1.0, 1e120] {
            let u = [scale, -scale, 0.0, -scale, scale, 0.0, 0.0, 0.0, 0.0];
            let (s, r) = schur_complement(
                &u,
                3,
                1,
                &SpectralOptions::default(),
                &mut WorkBudget::new(&Execution::default()).unwrap(),
            )
            .unwrap();
            assert!(s[0].abs() <= scale * 1e-10);
            assert!(r < 1e-10);
        }
        let u = [1.0, 1e-4, 0.0, 1e-4, 1e-8, 0.0, 0.0, 0.0, 1.0];
        let options = SpectralOptions {
            relative_tolerance: 1e-7,
            residual_tolerance: 1e-6,
            ..Default::default()
        };
        assert!(
            schur_complement(
                &u,
                3,
                1,
                &options,
                &mut WorkBudget::new(&Execution::default()).unwrap()
            )
            .is_err()
        );
    }
    #[test]
    fn schur_correction_has_no_intermediate_square_overflow_or_underflow() {
        for scale in [1e-200, 1e200] {
            let matrix = [scale, -scale, -scale, scale];
            let (result, _) = schur_complement(
                &matrix,
                2,
                1,
                &SpectralOptions::default(),
                &mut WorkBudget::new(&Execution::default()).unwrap(),
            )
            .unwrap();
            // The exact Schur complement is zero at every positive scale.
            // Solving C X = E^T gives X=-1 without squaring the scale.
            assert!(
                result[0].abs() <= scale * 1e-12,
                "scale={scale}, result={result:?}"
            );
        }
    }
    #[test]
    fn discarded_direction_must_satisfy_the_actual_solve_residual() {
        // Coupling 1e-10 is small relative to the whole Gram, but dropping it
        // loses the entire unit Schur correction. The correction guard and
        // actual solve check reject what a whole-Gram coupling test misses.
        let matrix = [1.0, 1e-10, 0.0, 1e-10, 1e-20, 0.0, 0.0, 0.0, 1.0];
        let mut budget = WorkBudget::new(&Execution::default()).unwrap();
        assert!(matches!(
            schur_complement(&matrix, 3, 1, &SpectralOptions::default(), &mut budget),
            Err(Error::NumericalFailure {
                context: "discarded Schur direction has significant correction"
            })
        ));
        let options = SpectralOptions {
            relative_tolerance: 1e-22,
            ..Default::default()
        };
        let (result, _) = schur_complement(&matrix, 3, 1, &options, &mut budget).unwrap();
        assert!(result[0].abs() < 1e-12);
    }

    #[test]
    fn large_retained_rhs_cannot_hide_a_discarded_unit_correction() {
        // PSD Gram of [(.5,1), (1,0), (0,1e-10)]. Its exact Schur
        // complement on the first coordinate is zero. A 1e-10 residual
        // relative to the other RHS entry .5 can still lose a unit correction.
        let matrix = [1.25, 0.5, 1e-10, 0.5, 1.0, 0.0, 1e-10, 0.0, 1e-20];
        let mut budget = WorkBudget::new(&Execution::default()).unwrap();
        assert!(matches!(
            schur_complement(&matrix, 3, 1, &SpectralOptions::default(), &mut budget),
            Err(Error::NumericalFailure { .. })
        ));
        let options = SpectralOptions {
            relative_tolerance: 1e-22,
            ..Default::default()
        };
        let (result, _) = schur_complement(&matrix, 3, 1, &options, &mut budget).unwrap();
        assert!(result[0].abs() < 1e-12);
    }
    #[test]
    fn storage_overflow_and_nonconvergence() {
        let mut b = WorkBudget::new(&Execution::default()).unwrap();
        assert!(matches!(
            dense(usize::MAX, 2, &SpectralOptions::default(), &mut b),
            Err(Error::SizeOverflow { .. })
        ));
        let o = SpectralOptions {
            max_sweeps: 0,
            ..Default::default()
        };
        assert!(eigensystem(&[2.0, -1.0, -1.0, 2.0], 2, &o, &mut b).is_err());
    }

    #[test]
    fn small_coupling_cannot_justify_discarding_a_small_pivot() {
        let options = SpectralOptions::default();
        let execution = Execution::default();
        let mut budget = WorkBudget::new(&execution).unwrap();
        for epsilon in [1e-4_f64, 1e-8, 1e-10] {
            for delta in [0.0, 1e-4] {
                let original = [1.0, epsilon, epsilon, epsilon * epsilon * (1.0 + delta)];
                let (direct, _) = schur_complement(&original, 2, 1, &options, &mut budget).unwrap();
                assert!((direct[0] - delta / (1.0 + delta)).abs() < 1e-12);
                let mut stepwise = original;
                if original[3] <= options.relative_tolerance {
                    assert!(matches!(
                        eliminate_last(&mut stepwise, 2, 2, 1.0, &options, &mut budget),
                        Err(Error::NumericalFailure {
                            context: "coupled near-zero stepwise pivot"
                        })
                    ));
                } else {
                    eliminate_last(&mut stepwise, 2, 2, 1.0, &options, &mut budget).unwrap();
                    assert!((stepwise[0] - direct[0]).abs() < 1e-12);
                }
            }
        }
        // A genuinely uncoupled zero coordinate has exactly zero correction.
        let mut uncoupled = [1.0, 0.0, 0.0, 0.0];
        eliminate_last(&mut uncoupled, 2, 2, 1.0, &options, &mut budget).unwrap();
        assert_eq!(uncoupled[0], 1.0);
        for scale in [1e-200, 1.0, 1e200] {
            let mut matrix = [scale, -scale, -scale, scale];
            eliminate_last(&mut matrix, 2, 2, scale, &options, &mut budget).unwrap();
            assert_eq!(matrix[0], 0.0);
        }
    }
}
