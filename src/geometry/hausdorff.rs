//! Direct comparison of finite Euclidean point sets.

use super::PointCloudView;
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result, canonical_zero};

/// Symmetric Euclidean Hausdorff distance between the supports of two clouds.
///
/// Computes `max(max_a min_b ||a-b||_2, max_b min_a ||b-a||_2)` without
/// constructing a filtration or computing persistence. Point counts may differ;
/// order and duplicate multiplicity do not affect the mathematical distance.
/// This is a metric on nonempty finite point sets, a pseudometric on clouds
/// retaining duplicates, and is not a rigid-motion-invariant shape distance.
///
/// Coordinate dimensions must agree. The caller must establish common coordinate
/// axes, origin, units and normalization: views contain no metadata certifying
/// them. No alignment or unit conversion is performed. The result has coordinate
/// units. Two empty clouds have distance zero; exactly one empty cloud gives
/// positive infinity, an extended-distance convention.
///
/// Uses binary64 stable norms with no algorithmic approximation. An overflowing
/// candidate is treated as infinity during nearest-point selection; a finite
/// nonempty-cloud answer that is unrepresentable is an error. Inputs are borrowed
/// and unchanged; scratch is O(1), work is O(n*m*d), with two directed scans.
///
/// # Errors
/// Returns [`Error::ShapeMismatch`] for unequal coordinate dimensions and
/// [`Error::NumericalFailure`] for an unrepresentable nonempty-cloud answer.
///
/// ```
/// use cocycle::geometry::{PointCloudView, hausdorff_distance};
/// let first = PointCloudView::new(&[0., 2.], 2, 1)?;
/// let second = PointCloudView::new(&[0., 1., 2.], 3, 1)?;
/// assert_eq!(hausdorff_distance(first, second)?, 1.);
/// # Ok::<(), cocycle::Error>(())
/// ```
pub fn hausdorff_distance(first: PointCloudView<'_>, second: PointCloudView<'_>) -> Result<f64> {
    hausdorff(first, second, &mut WorkBudget::unlimited())
}

/// Hausdorff distance with cooperative work and cancellation controls.
///
/// One fresh budget charges one compatibility check and one unit per coordinate
/// compared, including repeated visits in the two directed scans. A zero budget
/// rejects even empty clouds. Validated view construction is outside this budget.
/// Cancellation is polled at each coordinate; no partial scalar is returned and
/// inputs and the flag are unchanged. Counts are neither timings nor memory
/// quotas. Retrying starts a new budget.
///
/// # Errors
/// Includes [`hausdorff_distance`] errors, [`Error::Cancelled`] and
/// [`Error::WorkLimitExceeded`]. A pre-set cancellation flag is checked first.
pub fn hausdorff_distance_with(
    first: PointCloudView<'_>,
    second: PointCloudView<'_>,
    execution: &Execution<'_>,
) -> Result<f64> {
    if execution.is_unlimited() {
        hausdorff_distance(first, second)
    } else {
        hausdorff(first, second, &mut WorkBudget::new(execution)?)
    }
}

fn hausdorff<const CONTROLLED: bool>(
    first: PointCloudView<'_>,
    second: PointCloudView<'_>,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<f64> {
    budget.step()?;
    if first.dimension() != second.dimension() {
        return Err(Error::ShapeMismatch {
            input: "Hausdorff coordinate dimension",
            expected: first.dimension(),
            actual: second.dimension(),
        });
    }
    if first.is_empty() || second.is_empty() {
        return Ok(if first.is_empty() && second.is_empty() {
            0.0
        } else {
            f64::INFINITY
        });
    }
    let mut distance: f64 = 0.0;
    for (from, to) in [(first, second), (second, first)] {
        for a in from.coordinates().chunks_exact(from.dimension()) {
            let mut nearest = f64::INFINITY;
            for b in to.coordinates().chunks_exact(to.dimension()) {
                let mut norm: f64 = 0.0;
                for (&x, &y) in a.iter().zip(b) {
                    budget.step()?;
                    norm = norm.hypot(x - y);
                }
                nearest = nearest.min(norm);
            }
            distance = distance.max(nearest);
        }
    }
    if !distance.is_finite() {
        return Err(Error::NumericalFailure {
            context: "Hausdorff distance",
        });
    }
    Ok(canonical_zero(distance))
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicBool, Ordering};

    #[test]
    fn cancellation_after_actual_work_returns_no_scalar_and_retry_is_fresh() {
        let points = PointCloudView::new(&[0., 1.], 2, 1).unwrap();
        let flag = AtomicBool::new(false);
        let execution = Execution::default().max_work(100).cancellation(&flag);
        let mut budget = WorkBudget::new(&execution).unwrap();
        budget.cancel_at_work(3);
        assert_eq!(
            hausdorff(points, points, &mut budget),
            Err(Error::Cancelled)
        );
        assert!(flag.load(Ordering::Relaxed));
        flag.store(false, Ordering::Relaxed);
        assert_eq!(hausdorff_distance_with(points, points, &execution), Ok(0.));
    }
}
