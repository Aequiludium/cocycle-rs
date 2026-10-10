//! Checked scaling and reconstruction from original matching costs.
//! Adapted from Topp; the parent module retains the MIT attribution.

use super::{Metric, allocation, buffer, sum_size};
use crate::diagram::{DiagramDimension, IntervalEnd};
use crate::execution::WorkBudget;
use crate::{Error, Result};

pub(super) fn numerical() -> Error {
    Error::NumericalFailure {
        context: "Wasserstein distance",
    }
}

pub(super) fn finite(value: f64) -> Result<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(numerical())
    }
}

pub(super) fn restore_scale(value: f64, scale: f64) -> Result<f64> {
    let restored = finite(value * scale)?;
    if value > 0.0 && restored == 0.0 {
        return Err(numerical());
    }
    Ok(restored)
}

#[derive(Clone, Copy)]
pub(in crate::diagram_distances) struct Point {
    pub(super) coordinates: [f64; 2],
    pub(super) midpoint: f64,
    pub(super) half: f64,
}

pub(super) fn prepare_dimensions<const CONTROLLED: bool>(
    first: &DiagramDimension<'_>,
    second: &DiagramDimension<'_>,
    counts: (usize, usize),
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<(Vec<Point>, Vec<Point>, f64)> {
    // Find the common scale while filling the algorithm-owned vectors. Then
    // normalize those same vectors in place: no second logical-view scan and
    // no intermediate coordinate arrays, including on duplicate-heavy inputs.
    let (mut first, a, first_valid) =
        collect_points(first.iter().map(|i| (i.birth(), i.end())), counts.0, budget)?;
    let (mut second, b, second_valid) = collect_points(
        second.iter().map(|i| (i.birth(), i.end())),
        counts.1,
        budget,
    )?;
    let scale = scale_for(a.max(b));
    if scale == 1.0 {
        // Division by one leaves the fields computed during collection intact.
        // Wide exponent ranges still take the original checked scaling path.
        if !first_valid || !second_valid {
            return Err(numerical());
        }
    } else {
        normalize(&mut first, scale, budget)?;
        normalize(&mut second, scale, budget)?;
    }
    Ok((first, second, scale))
}

fn collect_points<const CONTROLLED: bool>(
    input: impl Iterator<Item = (f64, IntervalEnd)>,
    count: usize,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<(Vec<Point>, f64, bool)> {
    let mut result = Vec::new();
    if count == 0 {
        return Ok((result, 0.0, true));
    }
    result.try_reserve_exact(count).map_err(|_| allocation())?;
    let mut largest = 0.0_f64;
    let mut valid = true;
    for (birth, end) in input {
        budget.step()?;
        let IntervalEnd::Finite(death) = end else {
            continue;
        };
        largest = largest.max(birth.abs()).max(death.abs());
        let half = (death - birth) * 0.5;
        valid &= birth.is_finite()
            && death.is_finite()
            && birth < death
            && half.is_finite()
            && half > 0.0;
        result.push(Point {
            coordinates: [birth, death],
            midpoint: birth * 0.5 + death * 0.5,
            half,
        });
    }
    Ok((result, largest, valid))
}

fn normalize<const CONTROLLED: bool>(
    points: &mut [Point],
    scale: f64,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<()> {
    for chunk in points.chunks_mut(256) {
        budget.step_by(chunk.len())?;
        for point in chunk {
            let [birth, death] = point.coordinates;
            let b = birth / scale;
            let d = death / scale;
            // Scaling by a power of two must be reversible: no implicit quantization
            // or collapsed intervals are accepted when the exponent range is wide.
            if !b.is_finite()
                || !d.is_finite()
                || b * scale != birth
                || d * scale != death
                || b >= d
            {
                return Err(numerical());
            }
            let half = (d - b) * 0.5;
            if !half.is_finite() || half <= 0.0 {
                return Err(numerical());
            }
            *point = Point {
                coordinates: [b, d],
                midpoint: b * 0.5 + d * 0.5,
                half,
            };
        }
    }
    Ok(())
}

pub(super) fn prepare<const CONTROLLED: bool>(
    input: impl Iterator<Item = (f64, IntervalEnd)>,
    count: usize,
    scale: f64,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<Vec<Point>> {
    let (mut result, _, _) = collect_points(input, count, budget)?;
    normalize(&mut result, scale, budget)?;
    Ok(result)
}

pub(super) fn power_scale<const CONTROLLED: bool>(
    input: impl Iterator<Item = (f64, IntervalEnd)>,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<f64> {
    let mut largest = 0.0_f64;
    for (birth, end) in input {
        budget.step()?;
        if let IntervalEnd::Finite(death) = end {
            largest = largest.max(birth.abs()).max(death.abs());
        }
    }
    Ok(scale_for(largest))
}

fn scale_for(largest: f64) -> f64 {
    if largest == 0.0 || (2.0_f64.powi(-200)..=2.0_f64.powi(200)).contains(&largest) {
        return 1.0;
    }
    let exponent = ((largest.to_bits() >> 52) & 0x7ff) as i32 - 1023;
    // A normal scaling factor also scales the smallest subnormal up to 2^-52.
    2.0_f64.powi(exponent.max(-1022))
}

pub(super) fn square(value: f64) -> Result<f64> {
    let squared = finite(value * value)?;
    if value != 0.0 && squared == 0.0 {
        return Err(numerical());
    }
    Ok(squared)
}

pub(super) fn diagonal_power(point: Point, metric: Metric) -> Result<f64> {
    match metric {
        Metric::W1 => Ok(point.half),
        Metric::W2 => finite(2.0 * square(point.half)?),
    }
}

pub(super) fn cross_power(first: Point, second: Point, metric: Metric) -> Result<f64> {
    let b = finite(first.coordinates[0] - second.coordinates[0])?.abs();
    let d = finite(first.coordinates[1] - second.coordinates[1])?.abs();
    match metric {
        Metric::W1 => Ok(b.max(d)),
        Metric::W2 => finite(square(b)? + square(d)?),
    }
}

pub(super) fn saving(first: Point, second: Point, metric: Metric) -> Result<(f64, bool)> {
    // Rounded midpoints are only a search index, never the cost authority:
    // adjacent endpoint values need not have a representable midpoint.
    let baseline = finite(diagonal_power(first, metric)? + diagonal_power(second, metric)?)?;
    let cross = cross_power(first, second, metric)?;
    // In this regime, subtracting a small cross cost from the diagonal baseline
    // cannot retain enough significant bits to rank almost equal matchings.
    // Reconstructing the final cost is insufficient: choose the matching itself
    // using original costs. This is a numerical guard, not a changed tolerance
    // or a heuristic substitute for the normal saving-graph solver.
    let direct_cost_required = cross > 0.0 && cross <= 64.0 * f64::EPSILON * baseline;
    Ok((finite(baseline - cross)?, direct_cost_required))
}

pub(super) fn accumulate(sum: &mut f64, correction: &mut f64, cost: f64) -> Result<()> {
    // Kahan accumulation of original nonnegative costs avoids baseline-minus-
    // savings cancellation, especially for almost equal diagrams in W2.
    let adjusted = cost - *correction;
    let next = finite(*sum + adjusted)?;
    *correction = (next - *sum) - adjusted;
    *sum = next;
    Ok(())
}

pub(super) fn from_flows<const CONTROLLED: bool>(
    first: &[Point],
    second: &[Point],
    row_capacity: &[usize],
    column_capacity: &[usize],
    flows: &[(usize, usize, usize)],
    metric: Metric,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<f64> {
    let mut row_used = buffer(first.len(), 0)?;
    let mut col_used = buffer(second.len(), 0)?;
    let mut total = 0.0;
    let mut correction = 0.0;
    for (index, &(row, column, count)) in flows.iter().enumerate() {
        if CONTROLLED && index % 256 == 0 {
            budget.step_by((flows.len() - index).min(256))?;
        }
        row_used[row] = sum_size(row_used[row], count)?;
        col_used[column] = sum_size(col_used[column], count)?;
        accumulate(
            &mut total,
            &mut correction,
            finite(cross_power(first[row], second[column], metric)? * count as f64)?,
        )?;
    }
    for (row, &point) in first.iter().enumerate() {
        if CONTROLLED && row % 256 == 0 {
            budget.step_by((first.len() - row).min(256))?;
        }
        let remaining =
            row_capacity[row]
                .checked_sub(row_used[row])
                .ok_or(Error::InternalInvariant {
                    reason: "Wasserstein row capacity exceeded",
                })?;
        if remaining > 0 {
            accumulate(
                &mut total,
                &mut correction,
                finite(diagonal_power(point, metric)? * remaining as f64)?,
            )?;
        }
    }
    for (column, &point) in second.iter().enumerate() {
        if CONTROLLED && column % 256 == 0 {
            budget.step_by((second.len() - column).min(256))?;
        }
        let remaining = column_capacity[column]
            .checked_sub(col_used[column])
            .ok_or(Error::InternalInvariant {
                reason: "Wasserstein column capacity exceeded",
            })?;
        if remaining > 0 {
            accumulate(
                &mut total,
                &mut correction,
                finite(diagonal_power(point, metric)? * remaining as f64)?,
            )?;
        }
    }
    Ok(if metric == Metric::W2 {
        total.sqrt()
    } else {
        total
    })
}

pub(super) fn from_matching<const CONTROLLED: bool>(
    first: &[Point],
    second: &[Point],
    matching: Vec<Option<usize>>,
    metric: Metric,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<f64> {
    let mut flows = Vec::new();
    flows
        .try_reserve_exact(first.len().min(second.len()))
        .map_err(|_| allocation())?;
    for (row, column) in matching.into_iter().enumerate() {
        if CONTROLLED && row % 256 == 0 {
            budget.step_by((first.len() - row).min(256))?;
        }
        if let Some(column) = column {
            flows.push((row, column, 1));
        }
    }
    from_flows(
        first,
        second,
        &buffer(first.len(), 1)?,
        &buffer(second.len(), 1)?,
        &flows,
        metric,
        budget,
    )
}
