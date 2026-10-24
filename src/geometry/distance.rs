//! Euclidean set comparison and checked scalars for distance consumers.

use super::{PointCloudView, pair_count};
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result, canonical_zero};

/// Owned symmetric Hausdorff distances with one ordered identity list for both axes.
///
/// Values use condensed lower-triangle order `(1,0), (2,0), (2,1), ...`.
/// The diagonal is implicitly positive zero. Distinct sample IDs remain distinct
/// even when their point sets have distance zero. Coordinates are not retained;
/// IDs are caller labels, not authenticated dataset provenance or units.
#[derive(Debug, PartialEq)]
pub struct HausdorffDistanceMatrix {
    sample_ids: Vec<String>,
    values: Vec<f64>,
}

impl HausdorffDistanceMatrix {
    /// Caller-supplied sample IDs, in input order, identifying both matrix axes.
    pub fn sample_ids(&self) -> &[String] {
        &self.sample_ids
    }

    /// Condensed lower triangle, containing `len() choose 2` finite values.
    ///
    /// `DissimilarityView::new(matrix.values(), matrix.len())` borrows these
    /// values without copying, with an explicit O(values length) validation.
    pub fn values(&self) -> &[f64] {
        &self.values
    }

    /// Number of samples, including distinct IDs with identical point sets.
    pub fn len(&self) -> usize {
        self.sample_ids.len()
    }

    /// Whether the matrix has no samples.
    pub fn is_empty(&self) -> bool {
        self.sample_ids.is_empty()
    }

    /// Symmetric distance or positive-zero diagonal; `None` for invalid indices.
    pub fn get(&self, i: usize, j: usize) -> Option<f64> {
        if i >= self.len() || j >= self.len() {
            return None;
        }
        if i == j {
            return Some(0.0);
        }
        let (row, column) = if i > j { (i, j) } else { (j, i) };
        self.values
            .get(pair_count(row)?.checked_add(column)?)
            .copied()
    }
}

/// Symmetric Hausdorff distance between two nonempty Euclidean point sets.
///
/// Returns `max(max_a min_b ||a-b||_2, max_b min_a ||a-b||_2)`, in coordinate
/// units. Both clouds must have the same positive coordinate dimension. The
/// caller establishes common coordinates, physical units and normalization;
/// this call performs no registration or rescaling. Point order and duplicate
/// rows do not affect the mathematical value. This is a metric on finite sets,
/// and only a pseudometric on records retaining duplicate rows. It is sensitive
/// to outliers and does not establish positive definiteness of a derived kernel.
///
/// No algorithmic approximation is used; binary64 arithmetic still rounds.
/// Stable norms avoid squaring overflow/underflow. An unrepresentable distant
/// point pair is ignored if a representable nearer neighbor exists; only an
/// unrepresentable final Hausdorff value is a numerical error.
///
/// Two directed scans take O(m*n*d) time and O(1) auxiliary space. Coordinates
/// are borrowed without copying or modifying either input.
///
/// # Errors
/// Rejects either empty cloud, unequal dimensions or an unrepresentable result.
///
/// ```
/// use cocycle::geometry::{PointCloudView, hausdorff_distance};
/// let a = PointCloudView::new(&[0., 2.], 2, 1)?;
/// let b = PointCloudView::new(&[1., 3.], 2, 1)?;
/// assert_eq!(hausdorff_distance(a, b)?, 1.);
/// # Ok::<(), cocycle::Error>(())
/// ```
pub fn hausdorff_distance(first: PointCloudView<'_>, second: PointCloudView<'_>) -> Result<f64> {
    hausdorff_distance_with(first, second, &Execution::default())
}

/// Hausdorff distance with one fresh cooperative work/cancellation budget.
///
/// Counts operation entry, operand checks, point/coordinate visits and output.
/// Previously validated coordinates are outside validation work. Cancellation is
/// checked before validation and during coordinate scans; zero work rejects the
/// operation. Counts are implementation-dependent, not time or memory limits.
/// Failure returns no partial result and never changes the cancellation flag.
///
/// # Errors
/// In addition to [`hausdorff_distance`], returns cancellation or work exhaustion.
pub fn hausdorff_distance_with(
    first: PointCloudView<'_>,
    second: PointCloudView<'_>,
    execution: &Execution<'_>,
) -> Result<f64> {
    if execution.is_unlimited() {
        compare(first, second, &mut WorkBudget::unlimited())
    } else {
        compare(first, second, &mut WorkBudget::new(execution)?)
    }
}

/// Pairwise Hausdorff matrix retaining owned sample IDs in input order.
///
/// Requires nonempty, unique IDs and nonempty clouds of one coordinate dimension.
/// IDs are compared exactly, without trimming or normalization. Every sample is
/// checked before comparison, including the sole sample in a one-sample call.
/// Zero samples produce an empty matrix; one valid sample has a zero diagonal.
/// All samples must share the caller-established frame, units and normalization.
///
/// Each unordered pair is computed once. Storage is O(s*s + total ID bytes);
/// coordinate buffers are borrowed, with O(1) per-comparison scratch. ID uniqueness
/// checks use O(s*s) string comparisons without a registry or sorting. No caching,
/// internal parallelism, approximate neighbors or partial matrices are provided.
///
/// # Errors
/// Rejects invalid IDs/clouds/dimensions, capacity overflow, failed allocation or
/// an unrepresentable distance, including failures after earlier pairs succeeded.
pub fn hausdorff_distance_matrix(
    samples: &[(&str, PointCloudView<'_>)],
) -> Result<HausdorffDistanceMatrix> {
    hausdorff_distance_matrix_with(samples, &Execution::default())
}

/// Pairwise matrix with one shared budget across validation, copies and all pairs.
///
/// Includes operation entry, sample/ID-byte checks, ID copying, coordinate work
/// and stored values. The numeric and ID buffers use checked, fallible reserves.
/// Allocations cannot be interrupted internally. There is no memory quota or
/// deadline; cancellation is cooperative. Failure returns no partial matrix.
/// Even a zero-sample call rejects zero work or a pre-set cancellation flag.
///
/// # Errors
/// In addition to [`hausdorff_distance_matrix`], returns cancellation or work exhaustion.
pub fn hausdorff_distance_matrix_with(
    samples: &[(&str, PointCloudView<'_>)],
    execution: &Execution<'_>,
) -> Result<HausdorffDistanceMatrix> {
    if execution.is_unlimited() {
        matrix(samples, &mut WorkBudget::unlimited())
    } else {
        matrix(samples, &mut WorkBudget::new(execution)?)
    }
}

fn validate_cloud<const C: bool>(
    cloud: PointCloudView<'_>,
    dimension: usize,
    budget: &mut WorkBudget<'_, C>,
) -> Result<()> {
    budget.step()?;
    if cloud.is_empty() {
        return Err(Error::InvalidParameter {
            parameter: "point cloud",
            reason: "Hausdorff distance requires nonempty point sets",
        });
    }
    if cloud.dimension() != dimension {
        return Err(Error::InvalidParameter {
            parameter: "dimension",
            reason: "Hausdorff point clouds must have equal coordinate dimensions",
        });
    }
    Ok(())
}

fn compare<const C: bool>(
    first: PointCloudView<'_>,
    second: PointCloudView<'_>,
    budget: &mut WorkBudget<'_, C>,
) -> Result<f64> {
    budget.step()?;
    validate_cloud(first, first.dimension(), budget)?;
    validate_cloud(second, first.dimension(), budget)?;
    compare_validated(first, second, budget)
}

fn compare_validated<const C: bool>(
    first: PointCloudView<'_>,
    second: PointCloudView<'_>,
    budget: &mut WorkBudget<'_, C>,
) -> Result<f64> {
    let value = directed(first, second, budget)?.max(directed(second, first, budget)?);
    budget.step()?;
    if !value.is_finite() {
        return Err(Error::NumericalFailure {
            context: "Hausdorff distance",
        });
    }
    budget.check()?;
    Ok(canonical_zero(value))
}

fn directed<const C: bool>(
    first: PointCloudView<'_>,
    second: PointCloudView<'_>,
    budget: &mut WorkBudget<'_, C>,
) -> Result<f64> {
    let mut maximum: f64 = 0.0;
    for a in first.coordinates().chunks_exact(first.dimension()) {
        budget.step()?;
        let mut nearest = f64::INFINITY;
        for b in second.coordinates().chunks_exact(second.dimension()) {
            budget.step()?;
            let mut norm: f64 = 0.0;
            for (&x, &y) in a.iter().zip(b) {
                budget.step()?;
                norm = norm.hypot(x - y);
            }
            // Finite coordinates cannot produce NaN here. Infinity is a valid
            // internal upper bound, not necessarily the nearest-neighbor result.
            nearest = nearest.min(norm);
        }
        maximum = maximum.max(nearest);
    }
    Ok(maximum)
}

fn matrix<const C: bool>(
    samples: &[(&str, PointCloudView<'_>)],
    budget: &mut WorkBudget<'_, C>,
) -> Result<HausdorffDistanceMatrix> {
    budget.step()?;
    let count = pair_count(samples.len()).ok_or(Error::SizeOverflow {
        operation: "sample count choose 2",
    })?;
    let dimension = samples.first().map_or(1, |(_, cloud)| cloud.dimension());
    for (i, &(id, cloud)) in samples.iter().enumerate() {
        validate_cloud(cloud, dimension, budget)?;
        if id.is_empty() {
            return Err(Error::InvalidParameter {
                parameter: "sample ID",
                reason: "sample IDs must be nonempty and unique",
            });
        }
        for &(previous, _) in &samples[..i] {
            budget.step()?;
            if id.len() == previous.len() {
                let mut equal = true;
                for (a, b) in id.bytes().zip(previous.bytes()) {
                    budget.step()?;
                    if a != b {
                        equal = false;
                        break;
                    }
                }
                if equal {
                    return Err(Error::InvalidParameter {
                        parameter: "sample ID",
                        reason: "sample IDs must be nonempty and unique",
                    });
                }
            }
        }
    }
    let mut sample_ids = Vec::new();
    budget.check()?;
    sample_ids
        .try_reserve_exact(samples.len())
        .map_err(|_| Error::AllocationFailed {
            context: "Hausdorff sample IDs",
        })?;
    for &(id, _) in samples {
        let mut owned = String::new();
        budget.check()?;
        owned
            .try_reserve_exact(id.len())
            .map_err(|_| Error::AllocationFailed {
                context: "Hausdorff sample ID",
            })?;
        // Copy by character so cancellation also covers unusually long IDs.
        for character in id.chars() {
            budget.step_by(character.len_utf8())?;
            owned.push(character);
        }
        sample_ids.push(owned);
    }
    let mut values = Vec::new();
    budget.check()?;
    values
        .try_reserve_exact(count)
        .map_err(|_| Error::AllocationFailed {
            context: "Hausdorff distances",
        })?;
    for i in 0..samples.len() {
        for j in 0..i {
            let value = compare_validated(samples[i].1, samples[j].1, budget)?;
            budget.step()?;
            values.push(value);
        }
    }
    budget.check()?;
    Ok(HausdorffDistanceMatrix { sample_ids, values })
}

pub(crate) fn nonnegative(value: f64, field: &'static str, index: Option<usize>) -> Result<f64> {
    if !value.is_finite() {
        return Err(Error::NonFiniteValue { field, index });
    }
    if value < 0.0 {
        return Err(Error::NegativeValue { field, index });
    }
    Ok(canonical_zero(value))
}

pub(crate) fn cutoff(value: Option<f64>) -> Result<Option<f64>> {
    value.map(|v| nonnegative(v, "max_edge", None)).transpose()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicBool, Ordering};

    #[test]
    fn every_budget_boundary_and_mid_operation_cancellation_are_atomic() {
        let a = PointCloudView::new(&[0., 2.], 2, 1).unwrap();
        let b = PointCloudView::new(&[1., 3.], 2, 1).unwrap();
        let samples = [("alpha", a), ("beta", b), ("gamma", a)];
        for is_matrix in [false, true] {
            let run = |budget: &mut WorkBudget<'_>| {
                if is_matrix {
                    matrix(&samples, budget).map(|m| m.values[0])
                } else {
                    compare(a, b, budget)
                }
            };
            let mut budget = WorkBudget::new(&Execution::default().max_work(u64::MAX)).unwrap();
            assert_eq!(run(&mut budget), Ok(1.));
            let total = budget.used();
            for limit in 0..total {
                let mut budget = WorkBudget::new(&Execution::default().max_work(limit)).unwrap();
                assert_eq!(run(&mut budget), Err(Error::WorkLimitExceeded { limit }));
            }
            let mut exact = WorkBudget::new(&Execution::default().max_work(total)).unwrap();
            assert_eq!(run(&mut exact), Ok(1.));
            for at in [1, total / 2, total] {
                let flag = AtomicBool::new(false);
                let execution = Execution::new(Some(u64::MAX), Some(&flag));
                let mut budget = WorkBudget::new(&execution).unwrap();
                budget.cancel_at_work(at);
                assert_eq!(run(&mut budget), Err(Error::Cancelled));
                assert!(flag.load(Ordering::Relaxed));
                assert_eq!(
                    run(&mut WorkBudget::new(&Execution::default()).unwrap()),
                    Ok(1.)
                );
            }
        }
    }
}
