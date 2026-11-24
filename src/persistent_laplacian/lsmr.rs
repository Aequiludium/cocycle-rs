//! Undamped Section 2.8 LSMR, preserving the frozen Python reproduction scope.
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result};

/// Limits for the unoptimized dense-boundary LSMR reproduction.
#[derive(Clone, Copy, Debug)]
pub struct LsmrOptions {
    /// Stop when ||B2 z|| <= tolerance * ||B2||_F * ||y|| (default 1e-13).
    /// Must be finite and strictly between zero and one. Residuals are freshly
    /// evaluated; this is the project's rule, not the paper's norm estimator.
    pub relative_tolerance: f64,
    /// Iteration limit; None uses checked 4*min(rows,columns)+4.
    /// Zero accepts only an initially stationary least-squares problem.
    pub max_iterations: Option<usize>,
    /// Maximum entries in the borrowed row-major B2 (default 65,536).
    /// This bounds traversal; the caller already owns the input storage.
    pub max_matrix_entries: usize,
    /// Maximum entries in each temporary/result vector (default 65,536).
    /// Several vectors coexist; this is not an aggregate process-memory quota.
    pub max_vector_entries: usize,
}
impl Default for LsmrOptions {
    fn default() -> Self {
        Self {
            relative_tolerance: 1e-13,
            max_iterations: None,
            max_matrix_entries: 65_536,
            max_vector_entries: 65_536,
        }
    }
}

/// Numerical least-squares residual in caller-supplied upper-column coordinates.
///
/// No source identity, exact rank or fixed linear/self-adjoint projector is
/// certified. Preserve the originating boundary and ordered bases separately.
#[derive(Debug)]
pub struct LsmrProjection {
    /// Zero-start estimate w minimizing ||B2^T w-y||, in outside-row order.
    /// Minimum-norm behavior is numerical, not a rank or forward-error guarantee.
    pub solution: Vec<f64>,
    /// Fresh residual z=y-B2^T w, in upper-column order, without clipping.
    pub projected: Vec<f64>,
    /// Euclidean norm of z, including any unavoidable least-squares residual.
    pub residual_norm: f64,
    /// Euclidean norm of the freshly evaluated normal residual B2 z.
    /// Small normal residual alone does not bound projection error near singularity.
    pub normal_residual_norm: f64,
    /// Completed bidiagonalization/rotation updates; zero for initial stationarity.
    pub iterations: usize,
}

/// Estimate projection onto ker(B2) by undamped LSMR on B2^T w approximately y.
///
/// `boundary` is a borrowed, finite `outside_rows` by `upper_columns` row-major
/// matrix. It contains all outside-row upper-boundary columns in caller-declared
/// real coordinates. `y` has `upper_columns` entries. No source/chain validation,
/// weights, damping, sparse factorization or reorthogonalization is implicit.
///
/// Transcribes [Fong--Saunders (2011)](https://doi.org/10.1137/10079687X),
/// Section 2.8, steps 1-6, using the same fresh
/// normal-residual stopping rule as the frozen Python reproduction. Zero damping
/// and a zero initial solution target the minimum-norm least-squares solution.
/// Empty/zero matrices and a stationary initial residual return zero w and y.
/// Adaptive stopping generally depends on y: this diagnostic terminal does not
/// define the fixed action used by [`super::PersistentLaplacianOperator`].
///
/// One budget covers finite-input validation, all dense multiplies, rotations,
/// vector work and final residuals. Storage is O(rows+columns) beyond the borrowed
/// matrix; multiply work is O(rows*columns) per iteration. Ordinary library
/// builds retain no iteration trace. All resource limits are checked explicitly.
/// # Errors
/// Rejects invalid shape/policy/nonfinite data, size overflow, entry limits,
/// allocation failure, cancellation, work exhaustion, nonfinite arithmetic or
/// underflow of a required normal/rotation scale,
/// a zero rotation before convergence or an exhausted iteration limit.
/// No partial solution or successful unconverged projection is returned.
pub fn project_boundary_kernel_lsmr(
    boundary: &[f64],
    outside_rows: usize,
    upper_columns: usize,
    y: &[f64],
    options: LsmrOptions,
    execution: &Execution<'_>,
) -> Result<LsmrProjection> {
    let mut budget = WorkBudget::new(execution)?;
    solve(
        boundary,
        outside_rows,
        upper_columns,
        y,
        options,
        &mut budget,
        #[cfg(any(test, cocycle_lsmr_replay))]
        &mut |_| {},
    )
}

// Only the explicit standalone worker and private invariant tests compile this
// observation path. Production debug/release builds allocate no trace records.
#[cfg(any(test, cocycle_lsmr_replay))]
pub(crate) struct LsmrStep<'a> {
    pub iteration: usize,
    pub alpha: f64,
    pub beta: f64,
    pub rho: f64,
    pub rhobar: f64,
    pub zetabar: f64,
    pub u: &'a [f64],
    pub v: &'a [f64],
    pub solution: &'a [f64],
    pub residual: &'a [f64],
    pub normal_residual: &'a [f64],
}

#[cfg(cocycle_lsmr_replay)]
pub(crate) fn replay(
    boundary: &[f64],
    rows: usize,
    columns: usize,
    y: &[f64],
    options: LsmrOptions,
    execution: &Execution<'_>,
    observe: &mut impl FnMut(LsmrStep<'_>),
) -> Result<LsmrProjection> {
    solve(
        boundary,
        rows,
        columns,
        y,
        options,
        &mut WorkBudget::new(execution)?,
        observe,
    )
}

#[allow(clippy::too_many_arguments)]
fn solve(
    b: &[f64],
    rows: usize,
    columns: usize,
    y: &[f64],
    options: LsmrOptions,
    budget: &mut WorkBudget<'_>,
    #[cfg(any(test, cocycle_lsmr_replay))] observe: &mut impl FnMut(LsmrStep<'_>),
) -> Result<LsmrProjection> {
    let entries = rows.checked_mul(columns).ok_or(Error::SizeOverflow {
        operation: "LSMR boundary shape",
    })?;
    if b.len() != entries {
        return Err(Error::ShapeMismatch {
            input: "LSMR boundary",
            expected: entries,
            actual: b.len(),
        });
    }
    if y.len() != columns {
        return Err(Error::ShapeMismatch {
            input: "LSMR right side",
            expected: columns,
            actual: y.len(),
        });
    }
    if !options.relative_tolerance.is_finite()
        || options.relative_tolerance <= 0.
        || options.relative_tolerance >= 1.
    {
        return Err(invalid(
            "LSMR tolerance",
            "require a finite value strictly between zero and one",
        ));
    }
    if entries > options.max_matrix_entries {
        return Err(invalid(
            "max_matrix_entries",
            "LSMR input exceeds the declared entry limit",
        ));
    }
    if rows.max(columns) > options.max_vector_entries {
        return Err(invalid(
            "max_vector_entries",
            "LSMR coordinates exceed the declared entry limit",
        ));
    }
    let limit = match options.max_iterations {
        Some(value) => value,
        None => rows
            .min(columns)
            .checked_mul(4)
            .and_then(|v| v.checked_add(4))
            .ok_or(Error::SizeOverflow {
                operation: "LSMR iteration limit",
            })?,
    };
    for (field, values) in [("LSMR boundary", b), ("LSMR right side", y)] {
        for (index, value) in values.iter().enumerate() {
            budget.step()?;
            if !value.is_finite() {
                return Err(Error::NonFiniteValue {
                    field,
                    index: Some(index),
                });
            }
        }
    }
    let beta_initial = norm(y, budget)?;
    let boundary_norm = norm(b, budget)?;
    let normal_scale = boundary_norm * beta_initial;
    finite(normal_scale)?;
    if normal_scale == 0. && boundary_norm > 0. && beta_initial > 0. {
        return Err(failure("LSMR normal scale underflow"));
    }
    let mut u = copy(y, budget)?;
    normalize(&mut u, beta_initial, budget)?;
    let mut v = multiply(b, rows, columns, &u, false, budget)?;
    let mut alpha = norm(&v, budget)?;
    normalize(&mut v, alpha, budget)?;
    let mut x = zeros(rows, budget)?;
    let mut h = copy(&v, budget)?;
    let mut hbar = zeros(rows, budget)?;
    let mut alphabar = alpha;
    let mut zetabar = alpha * beta_initial;
    finite(zetabar)?;
    let (mut rho, mut rhobar, mut cbar, mut sbar) = (1., 1., 1., 0.);
    if alpha == 0. || beta_initial == 0. {
        let result = finish(b, rows, columns, y, x, 0, budget)?;
        if result.normal_residual_norm > options.relative_tolerance * normal_scale {
            return Err(failure("LSMR initial normal residual"));
        }
        return Ok(result);
    }
    for iteration in 0..limit {
        // M=B2^T: beta*u_next=M*v-alpha*u, alpha*v_next=M^T*u_next-beta*v.
        let mut raw_u = multiply(b, rows, columns, &v, true, budget)?;
        subtract_scaled(&mut raw_u, alpha, &u, budget)?;
        let beta = norm(&raw_u, budget)?;
        normalize(&mut raw_u, beta, budget)?;
        u = raw_u;
        let mut raw_v = multiply(b, rows, columns, &u, false, budget)?;
        subtract_scaled(&mut raw_v, beta, &v, budget)?;
        alpha = norm(&raw_v, budget)?;
        normalize(&mut raw_v, alpha, budget)?;
        v = raw_v;
        let (old_rho, old_rhobar) = (rho, rhobar);
        rho = alphabar.hypot(beta);
        if rho == 0. {
            return Err(failure("LSMR zero rotation before convergence"));
        }
        let (c, s) = (alphabar / rho, beta / rho);
        let theta = s * alpha;
        alphabar = c * alpha;
        let thetabar = sbar * rho;
        let rhotemp = cbar * rho;
        rhobar = rhotemp.hypot(theta);
        if rhobar == 0. {
            return Err(failure("LSMR zero secondary rotation"));
        }
        cbar = rhotemp / rhobar;
        sbar = theta / rhobar;
        let zeta = cbar * zetabar;
        zetabar *= -sbar;
        for value in [rho, rhobar, alphabar, zetabar, zeta] {
            finite(value)?;
        }
        let correction_numerator = thetabar * rho;
        let correction_denominator = old_rho * old_rhobar;
        let update_denominator = rho * rhobar;
        for value in [
            correction_numerator,
            correction_denominator,
            update_denominator,
        ] {
            finite(value)?;
        }
        if correction_denominator == 0. || update_denominator == 0. {
            return Err(failure("LSMR rotation scale underflow"));
        }
        let correction = correction_numerator / correction_denominator;
        let update = zeta / update_denominator;
        let next = theta / rho;
        for value in [correction, update, next] {
            finite(value)?;
        }
        for i in 0..rows {
            budget.step()?;
            hbar[i] = h[i] - correction * hbar[i];
            x[i] += update * hbar[i];
            h[i] = v[i] - next * h[i];
            for value in [hbar[i], x[i], h[i]] {
                finite(value)?;
            }
        }
        let residual = residual(b, rows, columns, y, &x, budget)?;
        let normal = multiply(b, rows, columns, &residual, false, budget)?;
        let normal_norm = norm(&normal, budget)?;
        #[cfg(any(test, cocycle_lsmr_replay))]
        observe(LsmrStep {
            iteration: iteration + 1,
            alpha,
            beta,
            rho,
            rhobar,
            zetabar,
            u: &u,
            v: &v,
            solution: &x,
            residual: &residual,
            normal_residual: &normal,
        });
        if normal_norm <= options.relative_tolerance * normal_scale {
            let residual_norm = norm(&residual, budget)?;
            budget.check()?;
            return Ok(LsmrProjection {
                solution: x,
                projected: residual,
                residual_norm,
                normal_residual_norm: normal_norm,
                iterations: iteration + 1,
            });
        }
    }
    Err(failure("LSMR iteration limit"))
}

fn finish(
    b: &[f64],
    rows: usize,
    columns: usize,
    y: &[f64],
    x: Vec<f64>,
    iterations: usize,
    budget: &mut WorkBudget<'_>,
) -> Result<LsmrProjection> {
    let projected = residual(b, rows, columns, y, &x, budget)?;
    let normal = multiply(b, rows, columns, &projected, false, budget)?;
    let residual_norm = norm(&projected, budget)?;
    let normal_residual_norm = norm(&normal, budget)?;
    budget.check()?;
    Ok(LsmrProjection {
        solution: x,
        projected,
        residual_norm,
        normal_residual_norm,
        iterations,
    })
}
fn residual(
    b: &[f64],
    rows: usize,
    columns: usize,
    y: &[f64],
    x: &[f64],
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let mut fitted = multiply(b, rows, columns, x, true, budget)?;
    for (value, &rhs) in fitted.iter_mut().zip(y) {
        budget.step()?;
        *value = rhs - *value;
        finite(*value)?;
    }
    Ok(fitted)
}
fn multiply(
    b: &[f64],
    rows: usize,
    columns: usize,
    x: &[f64],
    transpose: bool,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<f64>> {
    let mut result = zeros(if transpose { columns } else { rows }, budget)?;
    for i in 0..rows {
        for j in 0..columns {
            budget.step()?;
            if transpose {
                result[j] += b[i * columns + j] * x[i];
            } else {
                result[i] += b[i * columns + j] * x[j];
            }
        }
    }
    for &value in &result {
        budget.step()?;
        finite(value)?;
    }
    Ok(result)
}
fn subtract_scaled(
    target: &mut [f64],
    factor: f64,
    source: &[f64],
    budget: &mut WorkBudget<'_>,
) -> Result<()> {
    for (value, &entry) in target.iter_mut().zip(source) {
        budget.step()?;
        *value -= factor * entry;
        finite(*value)?;
    }
    Ok(())
}
fn normalize(v: &mut [f64], size: f64, budget: &mut WorkBudget<'_>) -> Result<()> {
    for value in v {
        budget.step()?;
        *value = if size == 0. { 0. } else { *value / size };
        finite(*value)?;
    }
    Ok(())
}
fn norm(v: &[f64], budget: &mut WorkBudget<'_>) -> Result<f64> {
    let mut value: f64 = 0.;
    for &entry in v {
        budget.step()?;
        value = value.hypot(entry);
    }
    finite(value)?;
    Ok(value)
}
fn zeros(count: usize, budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
    budget.check()?;
    let mut v = Vec::new();
    v.try_reserve_exact(count)
        .map_err(|_| Error::AllocationFailed {
            context: "LSMR vector",
        })?;
    v.resize(count, 0.);
    budget.check()?;
    Ok(v)
}
fn copy(v: &[f64], budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
    let mut result = zeros(v.len(), budget)?;
    for (value, &entry) in result.iter_mut().zip(v) {
        budget.step()?;
        *value = entry;
    }
    Ok(result)
}
fn finite(value: f64) -> Result<()> {
    if !value.is_finite() {
        return Err(failure("LSMR nonfinite arithmetic"));
    }
    Ok(())
}
fn failure(context: &'static str) -> Error {
    Error::NumericalFailure { context }
}
fn invalid(parameter: &'static str, reason: &'static str) -> Error {
    Error::InvalidParameter { parameter, reason }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicBool, Ordering};

    #[test]
    fn exact_breakdown_completes_the_rotation_and_fresh_residual() {
        let mut seen = 0;
        let mut observe = |step: LsmrStep<'_>| {
            seen += 1;
            assert_eq!(step.iteration, 1);
            assert_eq!(
                (step.alpha, step.beta, step.rho, step.rhobar, step.zetabar),
                (0., 0., 1., 1., 0.)
            );
            assert_eq!(step.u, &[0.]);
            assert_eq!(step.v, &[0.]);
            assert_eq!(step.solution, &[1.]);
            assert_eq!(step.residual, &[0.]);
            assert_eq!(step.normal_residual, &[0.]);
        };
        let result = solve(
            &[1.],
            1,
            1,
            &[1.],
            LsmrOptions::default(),
            &mut WorkBudget::new(&Execution::default()).unwrap(),
            &mut observe,
        )
        .unwrap();
        assert_eq!(seen, 1);
        assert_eq!(result.iterations, 1);
        assert_eq!(result.projected, [0.]);
    }

    #[test]
    fn cancellation_after_actual_update_discards_the_result() {
        let flag = AtomicBool::new(false);
        let execution = Execution::default().cancellation(&flag);
        let mut seen = 0;
        let mut observe = |_step: LsmrStep<'_>| {
            seen += 1;
            flag.store(true, Ordering::Relaxed);
        };
        let result = solve(
            &[1.],
            1,
            1,
            &[1.],
            LsmrOptions::default(),
            &mut WorkBudget::new(&execution).unwrap(),
            &mut observe,
        );
        assert!(matches!(result, Err(Error::Cancelled)));
        assert_eq!(seen, 1);
    }

    #[test]
    fn all_updates_and_final_residual_share_the_budget() {
        let matrix = [1., 0., 0., 0., 2., 0., 0., 0., 3.];
        let execution = Execution::default().max_work(100_000);
        let mut budget = WorkBudget::new(&execution).unwrap();
        let result = solve(
            &matrix,
            3,
            3,
            &[1.; 3],
            LsmrOptions::default(),
            &mut budget,
            &mut |_| {},
        )
        .unwrap();
        assert!(result.iterations >= 2);
        let used = budget.used();
        assert!(
            project_boundary_kernel_lsmr(
                &matrix,
                3,
                3,
                &[1.; 3],
                LsmrOptions::default(),
                &Execution::default().max_work(used)
            )
            .is_ok()
        );
        assert!(matches!(
            project_boundary_kernel_lsmr(
                &matrix,
                3,
                3,
                &[1.; 3],
                LsmrOptions::default(),
                &Execution::default().max_work(used - 1)
            ),
            Err(Error::WorkLimitExceeded { .. })
        ));
    }

    #[test]
    fn initial_stationarity_is_checked_in_original_coordinates() {
        // Adjacent y values can normalize to equal floats. Then B2*u is zero,
        // although the hand-derived original residual B2*y = a-b is nonzero.
        // Find the rounding collision with this platform's hypot implementation.
        for i in 0..1000 {
            let bits = 0x3ff8_0000_0000_0000 + i * 1_000_003;
            let a = f64::from_bits(bits);
            let b = f64::from_bits(bits + 1);
            let size = a.hypot(b);
            if a / size == b / size {
                assert_ne!(a - b, 0.);
                let options = LsmrOptions {
                    relative_tolerance: 1e-20,
                    ..Default::default()
                };
                assert!(matches!(
                    project_boundary_kernel_lsmr(
                        &[1., -1.],
                        1,
                        2,
                        &[a, b],
                        options,
                        &Execution::default()
                    ),
                    Err(Error::NumericalFailure {
                        context: "LSMR initial normal residual"
                    })
                ));
                return;
            }
        }
        panic!("failed to construct the normalized-coordinate rounding regression");
    }
}
