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
        // loses the entire unit Schur correction. A relative C X - E^T check
        // catches what a whole-Gram coupling check alone misses.
        let matrix = [1.0, 1e-10, 0.0, 1e-10, 1e-20, 0.0, 0.0, 0.0, 1.0];
        let mut budget = WorkBudget::new(&Execution::default()).unwrap();
        assert!(matches!(
            schur_complement(&matrix, 3, 1, &SpectralOptions::default(), &mut budget),
            Err(Error::NumericalFailure {
                context: "Schur solve residual exceeds tolerance"
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
}
