//! Standard real discrete Dirac operators on finite supplied chain complexes.
//!
//! This opt-in path reads signed integer incidence from [`FilteredComplex`],
//! independently of persistence and its coefficient field. It uses ordinary
//! nonaugmented chains and the Euclidean inner product: `K = d + d.transpose()`.
//! Degree truncation removes maps crossing either endpoint. No eigensolver,
//! weighted/complex adjoint, persistent Dirac or signal filter is provided.

use crate::complex::FilteredComplex;
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result};
use std::collections::HashMap;
use std::hash::Hash;
use std::ops::RangeInclusive;

/// An owned sparse Dirac operator with its ordered, source-local basis.
///
/// Coordinates are grouped by increasing degree; within a degree, source traversal
/// order is preserved. Empty degrees retain offsets. The operator freezes signed
/// incidence and IDs, but does not invent a global source identity or construction
/// coverage. Keep the original owner when interpreting IDs as cells/vertices.
/// The source may be released when only frozen coordinates are needed.
///
/// Only maps inside the selected degree interval are retained. Their boundary
/// square is checked exactly over the integers, before conversion to `f64`.
/// Consequently `K` is symmetric and `K squared` has the Hodge blocks of this
/// **same truncated chain interval**, not the full source's endpoint Laplacians.
/// Filtration values establish traversal validity; they are not weights/Grams.
///
/// Construction stores O(cells + incidence + degrees) data, with additional
/// temporary storage for duplicate aggregation and boundary-square validation.
/// Validation visits two-step incidence paths and may exceed linear work.
/// Action needs O(coordinates + incidence) work/storage; dense export is opt-in.
#[derive(Clone, Debug)]
pub struct DiracOperator<I> {
    degrees: RangeInclusive<usize>,
    basis: Vec<I>,
    offsets: Vec<usize>,
    boundaries: Vec<Vec<(usize, i64)>>,
    truncated: bool,
}

impl<I: Copy + Eq + Hash> DiracOperator<I> {
    /// Freeze a signed real operator on an inclusive interval of nonnegative degrees.
    ///
    /// No PH computation or implicit topology expansion is performed. Select a
    /// single scale explicitly by passing a stage as the source. Missing degrees
    /// are zero-dimensional spaces, including intervals entirely above the source.
    ///
    /// # Errors
    /// Rejects reversed/unrepresentable intervals, duplicate IDs, nonfinite or
    /// unordered filtration values, missing/later or wrong-degree boundary cells
    /// of selected cells, and nonzero boundary square within the retained interval.
    /// Zero incidence terms are ignored; duplicates are summed with checked `i64`
    /// arithmetic. Sums outside the exact `f64` integer range (+/-2^53), integer
    /// overflow in chain validation, and allocation failure are errors.
    /// Unselected cells' boundaries are not read or certified. Accessors are called
    /// only with this instance's IDs; a misbehaving source accessor may still panic.
    pub fn new<C: FilteredComplex<CellId = I>>(
        source: &C,
        degrees: RangeInclusive<usize>,
    ) -> Result<Self> {
        Self::new_with(source, degrees, &Execution::default())
    }

    /// Construct under one budget spanning reads, grouping and exact chain checks.
    ///
    /// Work units count cells, degree slots, incidence terms and two-step terms;
    /// sorting, allocation and source accessors cannot be interrupted internally.
    /// No partial operator is returned on cancellation or budget exhaustion.
    /// # Errors
    /// Includes [`Self::new`] errors, cancellation and work exhaustion.
    pub fn new_with<C: FilteredComplex<CellId = I>>(
        source: &C,
        degrees: RangeInclusive<usize>,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        let (&first, &last) = (degrees.start(), degrees.end());
        if degrees.is_empty() {
            return Err(Error::InvalidParameter {
                parameter: "degrees",
                reason: "degree interval must be nonempty and increasing",
            });
        }
        let offset_count = (last - first).checked_add(2).ok_or(Error::SizeOverflow {
            operation: "Dirac degree offsets",
        })?;
        budget.step_by(offset_count)?;
        let mut offsets = buffer(offset_count, 0_usize)?;
        // ID -> (source traversal position, degree, selected coordinate).
        let mut positions = HashMap::new();
        let mut selected = Vec::new();
        let mut previous = None;
        let mut truncated = false;
        for (position, cell) in source.cells().enumerate() {
            budget.step()?;
            let degree = source.dimension(cell);
            let value = source.filtration_value(cell);
            if !value.is_finite() || previous.is_some_and(|p| value < p) {
                return Err(invalid(position, "values must be finite and nondecreasing"));
            }
            previous = Some(value);
            positions.try_reserve(1).map_err(|_| allocation())?;
            if positions.insert(cell, (position, degree, None)).is_some() {
                return Err(invalid(position, "duplicate cell ID"));
            }
            if degrees.contains(&degree) {
                selected.try_reserve(1).map_err(|_| allocation())?;
                selected.push((degree, position, cell));
                offsets[degree - first + 1] += 1;
            } else {
                truncated = true;
            }
        }
        budget.check()?;
        selected.sort_unstable_by_key(|&(degree, position, _)| (degree, position));
        budget.check()?;
        for degree in 1..offsets.len() {
            offsets[degree] += offsets[degree - 1];
        }
        let mut basis = Vec::new();
        let mut boundaries: Vec<Vec<(usize, i64)>> = Vec::new();
        basis
            .try_reserve_exact(selected.len())
            .map_err(|_| allocation())?;
        boundaries
            .try_reserve_exact(selected.len())
            .map_err(|_| allocation())?;
        for (coordinate, &(_, _, cell)) in selected.iter().enumerate() {
            budget.step()?;
            // Every selected ID was inserted above and occurs only once.
            positions.get_mut(&cell).unwrap().2 = Some(coordinate);
            basis.push(cell);
        }
        for &(degree, position, cell) in &selected {
            budget.step()?;
            let mut terms = Vec::new();
            for (face, coefficient) in source.boundary(cell) {
                budget.step()?;
                if coefficient == 0 {
                    continue;
                }
                let &(face_position, face_degree, coordinate) = positions
                    .get(&face)
                    .ok_or(invalid(position, "boundary cell is missing"))?;
                if face_position >= position {
                    return Err(invalid(position, "boundary cell must precede its coface"));
                }
                if degree.checked_sub(1) != Some(face_degree) {
                    return Err(invalid(position, "boundary dimension must decrease by one"));
                }
                if let Some(row) = coordinate {
                    terms.try_reserve(1).map_err(|_| allocation())?;
                    terms.push((row, i64::from(coefficient)));
                }
            }
            budget.check()?;
            terms.sort_unstable_by_key(|&(row, _)| row);
            budget.check()?;
            let mut column: Vec<(usize, i64)> = Vec::new();
            column
                .try_reserve_exact(terms.len())
                .map_err(|_| allocation())?;
            for (row, coefficient) in terms {
                budget.step()?;
                if let Some((previous_row, sum)) = column.last_mut()
                    && *previous_row == row
                {
                    *sum = sum
                        .checked_add(coefficient)
                        .ok_or(Error::NumericalFailure {
                            context: "Dirac integer incidence sum",
                        })?;
                } else {
                    column.push((row, coefficient));
                }
            }
            column.retain(|&(_, coefficient)| coefficient != 0);
            if column
                .iter()
                .any(|&(_, c)| !(-(1_i64 << 53)..=(1_i64 << 53)).contains(&c))
            {
                return Err(Error::NumericalFailure {
                    context: "Dirac incidence is not exactly representable in f64",
                });
            }
            let mut twice: HashMap<usize, i128> = HashMap::new();
            for &(row, coefficient) in &column {
                for &(face, next) in &boundaries[row] {
                    budget.step()?;
                    twice.try_reserve(1).map_err(|_| allocation())?;
                    let sum = twice.entry(face).or_default();
                    // Each factor is bounded by 2^53, so its product fits i128.
                    *sum = sum
                        .checked_add(i128::from(coefficient) * i128::from(next))
                        .ok_or(Error::NumericalFailure {
                            context: "Dirac integer boundary-square sum",
                        })?;
                }
            }
            if twice.values().any(|&value| value != 0) {
                return Err(invalid(
                    position,
                    "boundary of boundary is nonzero over the integers",
                ));
            }
            boundaries.push(column);
        }
        budget.check()?;
        Ok(Self {
            degrees,
            basis,
            offsets,
            boundaries,
            truncated,
        })
    }

    /// Selected inclusive degree interval, including empty degrees.
    pub fn degrees(&self) -> RangeInclusive<usize> {
        self.degrees.clone()
    }

    /// Frozen source-local IDs in graded coordinate order.
    pub fn basis(&self) -> &[I] {
        &self.basis
    }

    /// Offsets of each selected degree, followed by the total coordinate count.
    /// Degree `j` occupies `offsets[j-a]..offsets[j-a+1]`, where `a` is the first
    /// selected degree. Empty degrees have equal adjacent offsets.
    pub fn degree_offsets(&self) -> &[usize] {
        &self.offsets
    }

    /// Whether at least one source cell was excluded by the selected degree range.
    /// This describes only the supplied topology, not a larger omitted skeleton.
    pub fn is_truncated(&self) -> bool {
        self.truncated
    }

    /// Nonzero `(row, column, value)` entries of `K`, without matrix allocation.
    /// Both each signed boundary entry and its transpose are emitted once.
    pub fn entries(&self) -> impl Iterator<Item = (usize, usize, f64)> + '_ {
        self.boundaries
            .iter()
            .enumerate()
            .flat_map(|(column, terms)| {
                terms.iter().flat_map(move |&(row, coefficient)| {
                    [
                        (row, column, coefficient as f64),
                        (column, row, coefficient as f64),
                    ]
                })
            })
    }

    /// Apply `K` to a finite vector in [`Self::basis`] order, returning owned coordinates.
    ///
    /// No whole Dirac matrix is allocated. Each boundary entry contributes to
    /// both adjacent degrees. The result has the same basis/offsets as `self`;
    /// retain that metadata when exporting the vector.
    /// # Errors
    /// Rejects wrong length, nonfinite input, nonfinite products/partial sums and
    /// allocation failure. Floating-point operations follow stored column order.
    pub fn action(&self, x: &[f64]) -> Result<Vec<f64>> {
        self.action_with(x, &Execution::default())
    }

    /// Apply under a fresh budget counting input/output coordinates and incidence.
    /// No partial vector is returned on failure; the operator remains reusable.
    /// # Errors
    /// Includes [`Self::action`] errors, cancellation and work exhaustion.
    pub fn action_with(&self, x: &[f64], execution: &Execution<'_>) -> Result<Vec<f64>> {
        let mut budget = WorkBudget::new(execution)?;
        if x.len() != self.basis.len() {
            return Err(Error::ShapeMismatch {
                input: "Dirac vector",
                expected: self.basis.len(),
                actual: x.len(),
            });
        }
        for (index, value) in x.iter().enumerate() {
            budget.step()?;
            if !value.is_finite() {
                return Err(Error::NonFiniteValue {
                    field: "Dirac vector",
                    index: Some(index),
                });
            }
        }
        budget.step_by(x.len())?;
        let mut result = buffer(x.len(), 0.0)?;
        for (column, terms) in self.boundaries.iter().enumerate() {
            for &(row, coefficient) in terms {
                budget.step()?;
                let c = coefficient as f64;
                result[row] = accumulate(result[row], c, x[column])?;
                result[column] = accumulate(result[column], c, x[row])?;
            }
        }
        budget.check()?;
        Ok(result)
    }

    /// Export owned row-major `f64` entries of the square matrix, bounded explicitly.
    /// Shape is `basis.len()` by `basis.len()`; metadata remains on `self`.
    /// `max_entries` caps the number of scalars **before** dense allocation.
    /// An empty operator exports an empty buffer even with a zero cap.
    /// # Errors
    /// Rejects shape overflow, a matrix exceeding the cap and allocation failure.
    pub fn matrix(&self, max_entries: usize) -> Result<Vec<f64>> {
        self.matrix_with(max_entries, &Execution::default())
    }

    /// Export under a fresh budget counting dense entries and signed incidence.
    /// # Errors
    /// Includes [`Self::matrix`] errors, cancellation and work exhaustion.
    pub fn matrix_with(&self, max_entries: usize, execution: &Execution<'_>) -> Result<Vec<f64>> {
        let mut budget = WorkBudget::new(execution)?;
        let n = self.basis.len();
        let count = n.checked_mul(n).ok_or(Error::SizeOverflow {
            operation: "Dirac dense matrix",
        })?;
        if count > max_entries {
            return Err(Error::InvalidParameter {
                parameter: "max_entries",
                reason: "Dirac dense matrix exceeds the entry limit",
            });
        }
        budget.step_by(count)?;
        let mut result = buffer(count, 0.0)?;
        for (row, column, value) in self.entries() {
            budget.step()?;
            result[row * n + column] = value;
        }
        budget.check()?;
        Ok(result)
    }
}

fn buffer<T: Clone>(count: usize, value: T) -> Result<Vec<T>> {
    let mut result = Vec::new();
    result.try_reserve_exact(count).map_err(|_| allocation())?;
    result.resize(count, value);
    Ok(result)
}

fn accumulate(sum: f64, coefficient: f64, value: f64) -> Result<f64> {
    let result = sum + coefficient * value;
    if !result.is_finite() {
        return Err(Error::NumericalFailure {
            context: "Dirac action",
        });
    }
    Ok(crate::canonical_zero(result))
}

fn invalid(cell: usize, reason: &'static str) -> Error {
    Error::InvalidComplex {
        cell: Some(cell),
        reason,
    }
}

fn allocation() -> Error {
    Error::AllocationFailed {
        context: "Dirac operator",
    }
}
