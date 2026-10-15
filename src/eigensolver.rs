//! Real symmetric matrix-free block iteration and its small dense Ritz solver.
use crate::execution::WorkBudget;
use crate::{Error, Result};

/// Resource and residual policy for identity-preconditioned block LOBPCG.
///
/// Uses O(n*b + b²) storage, with b=min(block_size,n), and Ritz problems of at
/// most 3*b columns. Each retained block and the projected square matrix must
/// fit max_basis_entries; this is not an aggregate process-memory quota.
/// Residual stopping is absolute_tolerance + relative_tolerance * |lambda|.
/// Two-pass orthogonalization drops relative directions below 1e-13.
/// A deterministic random start cannot certify extremal spectral coverage.
#[derive(Clone, Copy, Debug)]
pub struct IterativeOptions {
    /// Maximum block updates after the initial Rayleigh-Ritz step.
    pub max_iterations: usize,
    /// Ritz block size, including lookahead; must be positive and >= requested count.
    pub block_size: usize,
    /// Maximum entries in a basis/image block or projected square buffer.
    pub max_basis_entries: usize,
    /// Finite nonnegative absolute residual tolerance, in eigenvalue units.
    pub absolute_tolerance: f64,
    /// Finite nonnegative relative residual tolerance.
    pub relative_tolerance: f64,
    /// Deterministic SplitMix64 starting-block seed; not a coverage certificate.
    pub seed: u64,
}
impl Default for IterativeOptions {
    fn default() -> Self {
        Self {
            max_iterations: 800,
            block_size: 8,
            max_basis_entries: 1_048_576,
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-9,
            seed: 73,
        }
    }
}

pub(crate) struct Eigenpairs {
    pub values: Vec<f64>,
    pub vectors: Vec<f64>,
    pub residuals: Vec<f64>,
    pub orthogonality_error: f64,
    pub iterations: usize,
    pub actions: usize,
    pub peak_subspace: usize,
}

/// The caller supplies a symmetric operator in one fixed real orthonormal basis.
pub(crate) fn lobpcg(
    n: usize,
    count: usize,
    largest: bool,
    options: &IterativeOptions,
    budget: &mut WorkBudget<'_>,
    mut action: impl FnMut(&[f64], &mut WorkBudget<'_>) -> Result<Vec<f64>>,
) -> Result<Eigenpairs> {
    if options.block_size == 0 || count > n || count > options.block_size {
        return Err(invalid(
            "iterative count/block size",
            "require count <= min(n, block_size), block_size > 0",
        ));
    }
    for value in [options.absolute_tolerance, options.relative_tolerance] {
        if !value.is_finite() || value < 0. {
            return Err(invalid(
                "iterative tolerance",
                "must be finite and nonnegative",
            ));
        }
    }
    if count == 0 {
        return Ok(Eigenpairs {
            values: Vec::new(),
            vectors: Vec::new(),
            residuals: Vec::new(),
            orthogonality_error: 0.,
            iterations: 0,
            actions: 0,
            peak_subspace: 0,
        });
    }
    let b = options.block_size.min(n);
    let max_columns = b.checked_mul(3).ok_or(overflow())?.min(n);
    for entries in [
        n.checked_mul(max_columns),
        max_columns.checked_mul(max_columns),
    ] {
        if entries.ok_or(overflow())? > options.max_basis_entries {
            return Err(invalid(
                "iterative basis entries",
                "exceeds the explicit workspace entry limit",
            ));
        }
    }
    let mut seed = options.seed;
    let mut x = Vec::new();
    for _ in 0..b {
        let mut v = zeros(n, budget)?;
        for value in &mut v {
            seed = seed.wrapping_add(0x9e3779b97f4a7c15);
            let mut z = seed;
            z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
            z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
            z ^= z >> 31;
            *value = (z >> 11) as f64 * (1. / ((1_u64 << 53) as f64)) - 0.5;
        }
        orthogonalize(v, &mut x, budget)?;
    }
    // Complete a rare deficient random block without changing the block size.
    for i in 0..n {
        if x.len() == b {
            break;
        }
        let mut v = zeros(n, budget)?;
        v[i] = 1.;
        orthogonalize(v, &mut x, budget)?;
    }
    if x.len() != b {
        return Err(failure("iterative starting rank"));
    }
    let mut residual_block = Vec::new();
    let mut directions = Vec::new();
    let mut actions = 0_usize;
    let mut peak_subspace = 0;
    for iteration in 0..=options.max_iterations {
        let mut basis = Vec::new();
        for column in x.iter().chain(&residual_block).chain(&directions) {
            if basis.len() == max_columns {
                break;
            }
            orthogonalize(copy(column, budget)?, &mut basis, budget)?;
        }
        if basis.len() < b {
            return Err(failure("iterative search rank"));
        }
        let m = basis.len();
        peak_subspace = peak_subspace.max(m);
        let mut images = Vec::new();
        for column in &basis {
            let image = action(column, budget)?;
            actions = actions.checked_add(1).ok_or(overflow())?;
            if image.len() != n {
                return Err(failure("iterative action shape"));
            }
            push(&mut images, image)?;
        }
        let mut projected = zeros(m.checked_mul(m).ok_or(overflow())?, budget)?;
        for i in 0..m {
            for j in i..m {
                let value = dot(&basis[i], &images[j], budget)?;
                projected[i * m + j] = value;
                projected[j * m + i] = value;
            }
        }
        let (diagonal, rotations, _) = jacobi(&projected, m, 0., 4. * f64::EPSILON, 64, budget)?;
        let mut order = Vec::new();
        for i in 0..m {
            push(&mut order, i)?;
        }
        budget.check()?;
        order.sort_unstable_by(|&i, &j| {
            let comparison = diagonal[i].total_cmp(&diagonal[j]);
            if largest {
                comparison.reverse().then(i.cmp(&j))
            } else {
                comparison.then(i.cmp(&j))
            }
        });
        budget.check()?;
        let old = x;
        x = Vec::new();
        residual_block = Vec::new();
        let mut values = Vec::new();
        let mut errors = Vec::new();
        for &index in &order[..b] {
            let mut vector = zeros(n, budget)?;
            let mut image = zeros(n, budget)?;
            for j in 0..m {
                budget.step_by(n)?;
                let c = rotations[j * m + index];
                for row in 0..n {
                    vector[row] += c * basis[j][row];
                    image[row] += c * images[j][row];
                }
            }
            let value = diagonal[index];
            budget.step_by(n)?;
            for row in 0..n {
                image[row] -= value * vector[row];
            }
            let error = norm(&image);
            finite(error)?;
            push(&mut values, value)?;
            push(&mut errors, error)?;
            push(&mut x, vector)?;
            push(&mut residual_block, image)?;
        }
        if errors[..count]
            .iter()
            .zip(&values)
            .all(|(&r, &v)| r <= options.absolute_tolerance + options.relative_tolerance * v.abs())
        {
            let mut actual_values = Vec::new();
            let mut actual_errors = Vec::new();
            // Fresh sparse applications verify the returned pairs independently
            // of cached Ritz image combinations and their rounding drift.
            for column in &x[..count] {
                let mut image = action(column, budget)?;
                actions = actions.checked_add(1).ok_or(overflow())?;
                let value = dot(column, &image, budget)?;
                budget.step_by(n)?;
                for i in 0..n {
                    image[i] -= value * column[i];
                }
                let error = norm(&image);
                finite(error)?;
                push(&mut actual_values, value)?;
                push(&mut actual_errors, error)?;
            }
            if actual_errors.iter().zip(&actual_values).all(|(&r, &v)| {
                r <= options.absolute_tolerance + options.relative_tolerance * v.abs()
            }) {
                let mut order: Vec<usize> = Vec::new();
                for i in 0..count {
                    push(&mut order, i)?;
                }
                order.sort_unstable_by(|&i, &j| {
                    let c = actual_values[i].total_cmp(&actual_values[j]);
                    if largest { c.reverse() } else { c }
                });
                let mut vectors = zeros(n.checked_mul(count).ok_or(overflow())?, budget)?;
                let mut sorted_values = Vec::new();
                let mut residuals = Vec::new();
                for (j, &index) in order.iter().enumerate() {
                    push(&mut sorted_values, actual_values[index])?;
                    push(&mut residuals, actual_errors[index])?;
                    budget.step_by(n)?;
                    for row in 0..n {
                        vectors[row * count + j] = x[index][row];
                    }
                }
                let mut orthogonality_error = 0_f64;
                for i in 0..count {
                    for j in 0..count {
                        orthogonality_error = orthogonality_error
                            .hypot(dot(&x[i], &x[j], budget)? - if i == j { 1. } else { 0. });
                    }
                }
                if orthogonality_error > 1e-8 {
                    return Err(failure("iterative orthogonality"));
                }
                budget.check()?;
                return Ok(Eigenpairs {
                    values: sorted_values,
                    vectors,
                    residuals,
                    orthogonality_error,
                    iterations: iteration,
                    actions,
                    peak_subspace,
                });
            }
        }
        if iteration == options.max_iterations {
            return Err(failure("LOBPCG convergence"));
        }
        directions = Vec::new();
        for column in &x {
            let mut direction = copy(column, budget)?;
            // P'=X' - X (X^T X') spans the non-X part of the local minimizer.
            for previous in &old {
                let c = dot(previous, column, budget)?;
                budget.step_by(n)?;
                for i in 0..n {
                    direction[i] -= c * previous[i];
                }
            }
            push(&mut directions, direction)?;
        }
    }
    Err(failure("LOBPCG convergence"))
}

fn orthogonalize(
    mut v: Vec<f64>,
    basis: &mut Vec<Vec<f64>>,
    budget: &mut WorkBudget<'_>,
) -> Result<()> {
    budget.step_by(v.len())?;
    let scale = norm(&v);
    finite(scale)?;
    if scale == 0. {
        return Ok(());
    }
    for x in &mut v {
        *x /= scale;
    }
    for _ in 0..2 {
        for q in basis.iter() {
            let c = dot(&v, q, budget)?;
            budget.step_by(v.len())?;
            for i in 0..v.len() {
                v[i] -= c * q[i];
            }
        }
    }
    budget.step_by(v.len())?;
    let length = norm(&v);
    finite(length)?;
    if length <= 1e-13 {
        return Ok(());
    }
    for x in &mut v {
        *x /= length;
    }
    push(basis, v)
}

fn dot(x: &[f64], y: &[f64], budget: &mut WorkBudget<'_>) -> Result<f64> {
    budget.step_by(x.len())?;
    if x.len() != y.len() {
        return Err(failure("iterative action shape"));
    }
    let value: f64 = x.iter().zip(y).map(|(a, b)| a * b).sum();
    finite(value)?;
    Ok(value)
}
fn norm(x: &[f64]) -> f64 {
    x.iter().fold(0., |a, &b| a.hypot(b))
}
fn zeros(n: usize, budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
    n.checked_mul(8)
        .filter(|&b| b <= isize::MAX as usize)
        .ok_or(overflow())?;
    budget.step_by(n)?;
    let mut v = Vec::new();
    v.try_reserve_exact(n).map_err(|_| allocation())?;
    v.resize(n, 0.);
    budget.check()?;
    Ok(v)
}
fn copy(x: &[f64], budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
    let mut v = zeros(x.len(), budget)?;
    v.copy_from_slice(x);
    Ok(v)
}
fn push<T>(v: &mut Vec<T>, x: T) -> Result<()> {
    v.try_reserve(1).map_err(|_| allocation())?;
    v.push(x);
    Ok(())
}
fn finite(x: f64) -> Result<()> {
    if x.is_finite() {
        Ok(())
    } else {
        Err(failure("symmetric arithmetic"))
    }
}
fn failure(context: &'static str) -> Error {
    Error::NumericalFailure { context }
}
fn overflow() -> Error {
    Error::SizeOverflow {
        operation: "symmetric eigensolver",
    }
}
fn allocation() -> Error {
    Error::AllocationFailed {
        context: "symmetric eigensolver",
    }
}
fn invalid(parameter: &'static str, reason: &'static str) -> Error {
    Error::InvalidParameter { parameter, reason }
}

pub(crate) fn jacobi(
    matrix: &[f64],
    n: usize,
    absolute_tolerance: f64,
    relative_tolerance: f64,
    max_sweeps: usize,
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
        absolute_tolerance / scale + relative_tolerance * norm(&a)
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
        if sweeps == max_sweeps {
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
    use crate::execution::Execution;
    use std::sync::atomic::{AtomicBool, Ordering};

    #[test]
    fn cancellation_after_an_operator_call_returns_no_pairs() {
        let flag = AtomicBool::new(false);
        let policy = Execution::default().cancellation(&flag);
        let mut budget = WorkBudget::new(&policy).unwrap();
        let mut calls = 0;
        let result = lobpcg(
            32,
            1,
            false,
            &IterativeOptions::default(),
            &mut budget,
            |x, _| {
                calls += 1;
                flag.store(true, Ordering::Relaxed);
                Ok(x.to_vec())
            },
        );
        assert!(matches!(result, Err(Error::Cancelled)));
        assert!(calls >= 1);
    }

    #[test]
    fn residuals_do_not_certify_extremal_index_coverage() {
        let n = 32;
        let options = IterativeOptions::default();
        let mut budget = WorkBudget::new(&Execution::default()).unwrap();
        let start = lobpcg(
            n,
            options.block_size,
            false,
            &options,
            &mut budget,
            |x, _| Ok(vec![0.; x.len()]),
        )
        .unwrap();
        let mut u = vec![0.; n];
        u[0] = 1.;
        for _ in 0..2 {
            for j in 0..options.block_size {
                let c: f64 = (0..n)
                    .map(|i| u[i] * start.vectors[i * options.block_size + j])
                    .sum();
                for (i, x) in u.iter_mut().enumerate() {
                    *x -= c * start.vectors[i * options.block_size + j];
                }
            }
        }
        let length = norm(&u);
        for x in &mut u {
            *x /= length;
        }
        // A=I-uu^T has an independently known zero mode u. The fixed start
        // misses it, so residual-small returned vectors can have eigenvalue 1.
        let result = lobpcg(n, 1, false, &options, &mut budget, |x, _| {
            let c: f64 = u.iter().zip(x).map(|(a, b)| a * b).sum();
            Ok(x.iter().zip(&u).map(|(x, u)| x - c * u).collect())
        })
        .unwrap();
        assert!((result.values[0] - 1.).abs() < 1e-10);
        assert!(result.residuals[0] < 1e-10);
        let c: f64 = u.iter().map(|x| x * x).sum();
        assert!((1. - c).abs() < 1e-14);
    }
}
