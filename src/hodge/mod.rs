//! Single-scale real Hodge Laplacians with the standard inner product.
//!
//! For the ordered window `C[k-1] <-A- C[k] <-D- C[k+1]`, the full operator
//! is `A^T A + D D^T`. Degree zero is ordinary, nonaugmented homology. This
//! module reads a finite supplied complex or stage without running persistence
//! or expanding an implicit Rips source. Filtration values are not weights.
//! General Gram matrices, complex coefficients, normalized and persistent
//! Laplacians are outside this API's explicitly real, standard-inner-product
//! definition. See Lim, *Hodge Laplacians on Graphs* (2020), Theorems 5.2–5.3,
//! <https://doi.org/10.1137/18M1223101>.

mod spectrum;
pub use spectrum::{HodgeSpectrum, SpectrumOptions};

use crate::complex::FilteredComplex;
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result};
use std::collections::HashMap;
use std::hash::Hash;

/// Default maximum entries in each dense Hodge square buffer: 1,048,576 (8 MiB).
///
/// Shared by [`HodgeLaplacian::matrix`] and [`SpectrumOptions`]. This bounds
/// each n by n f64 buffer, not aggregate allocations or process memory.
pub const DEFAULT_MAX_DENSE_ENTRIES: usize = 1_048_576;

/// Summands of the ordinary real Hodge Laplacian.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum HodgeSelection {
    /// `A^T A + D D^T` on the entire supplied finite complex or stage.
    Full,
    /// `A^T A` only; its nullity is generally not the Betti number.
    Lower,
    /// `D D^T` only; its nullity is generally not the Betti number.
    Upper,
}

/// An owned, validated sparse three-degree boundary snapshot.
///
/// Basis order is source traversal restricted to each degree; handles remain
/// owner-local. The owned coefficients and bases define the coordinate operator
/// after source release, but callers must retain their own source descriptions
/// to interpret handles as geometric cells or compare different owners.
/// Missing degrees in this supplied finite source are explicit zero-size maps.
/// No completeness claim is made about an omitted larger complex or skeleton.
/// Preparation scans all cell IDs and reads boundaries only in degrees k/k+1;
/// action uses `O(n[k-1]+n[k])` scratch and sparse multiply work. Dense storage
/// is allocated only by an explicit matrix or spectrum request.
pub struct HodgeLaplacian<Id> {
    degree: usize,
    selection: HodgeSelection,
    previous_basis: Vec<Id>,
    basis: Vec<Id>,
    next_basis: Vec<Id>,
    lower: Vec<Vec<(usize, i64)>>,
    upper: Vec<Vec<(usize, i64)>>,
}

impl<Id: Copy + Eq + Hash> HodgeLaplacian<Id> {
    /// Extract and check a real chain window, preserving signed incidence.
    ///
    /// Uses unlimited [`Execution::default`] controls for preparation. To limit
    /// preparation itself, call [`Self::new_with`]; controls passed to a later
    /// action/matrix/spectrum call do not retroactively cover this work.
    ///
    /// # Errors
    /// Rejects duplicate IDs, unknown/nonpreceding faces, wrong boundary degrees,
    /// nonzero vertex boundaries, integer boundary-square violations, coefficient
    /// or size overflow and allocation failure. Duplicate incidence terms are
    /// summed before testing `AD=0` exactly with checked i128 arithmetic.
    /// Source accessors must satisfy the valid-ID preconditions of FilteredComplex.
    pub fn new<C: FilteredComplex<CellId = Id>>(
        source: &C,
        degree: usize,
        selection: HodgeSelection,
    ) -> Result<Self> {
        Self::new_with(source, degree, selection, &Execution::default())
    }

    /// Extract under one cooperative traversal, incidence and chain-check budget.
    ///
    /// This constructor starts its own budget. Reusing the same Execution value
    /// for a later terminal starts a fresh budget, not a cumulative pipeline limit.
    ///
    /// # Errors
    /// Includes [`Self::new`] errors, cancellation and work exhaustion. Iterator
    /// bodies and sorting cannot be interrupted internally. Filtration values
    /// and unrelated chain windows are not revalidated by this operation.
    pub fn new_with<C: FilteredComplex<CellId = Id>>(
        source: &C,
        degree: usize,
        selection: HodgeSelection,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        let next_degree = degree.checked_add(1).ok_or(size_error())?;
        let mut cells = HashMap::new();
        let mut previous_basis = Vec::new();
        let mut basis = Vec::new();
        let mut next_basis = Vec::new();
        for (position, cell) in source.cells().enumerate() {
            budget.step()?;
            let dimension = source.dimension(cell);
            let coordinate = if degree.checked_sub(1) == Some(dimension) {
                let index = previous_basis.len();
                push(&mut previous_basis, cell)?;
                index
            } else if dimension == degree {
                let index = basis.len();
                push(&mut basis, cell)?;
                index
            } else if dimension == next_degree {
                let index = next_basis.len();
                push(&mut next_basis, cell)?;
                index
            } else {
                0
            };
            cells.try_reserve(1).map_err(|_| allocation())?;
            if cells
                .insert(cell, (position, dimension, coordinate))
                .is_some()
            {
                return Err(invalid(Some(position), "duplicate cell ID"));
            }
        }
        let mut lower = Vec::new();
        let mut upper = Vec::new();
        for &cell in &basis {
            let column = read_column(source, cell, &cells, &mut budget)?;
            push(&mut lower, column)?;
        }
        for &cell in &next_basis {
            let column = read_column(source, cell, &cells, &mut budget)?;
            push(&mut upper, column)?;
        }
        for column in &upper {
            let mut product: HashMap<usize, i128> = HashMap::new();
            for &(middle, right) in column {
                for &(row, left) in &lower[middle] {
                    budget.step()?;
                    product.try_reserve(1).map_err(|_| allocation())?;
                    let entry = product.entry(row).or_default();
                    *entry = entry
                        .checked_add(i128::from(left) * i128::from(right))
                        .ok_or(size_error())?;
                }
            }
            if product.values().any(|&value| value != 0) {
                return Err(invalid(
                    None,
                    "boundary square is nonzero over the integers",
                ));
            }
        }
        budget.check()?;
        Ok(Self {
            degree,
            selection,
            previous_basis,
            basis,
            next_basis,
            lower,
            upper,
        })
    }

    /// Requested chain degree.
    pub fn degree(&self) -> usize {
        self.degree
    }
    /// Selected full, lower or upper operator.
    pub fn selection(&self) -> HodgeSelection {
        self.selection
    }
    /// Owned source-local handles in coordinate order for `C[k]`.
    pub fn basis(&self) -> &[Id] {
        &self.basis
    }
    /// Owned source-local handles for `C[k-1]`; empty for ordinary degree zero.
    pub fn previous_basis(&self) -> &[Id] {
        &self.previous_basis
    }
    /// Owned source-local handles for `C[k+1]`.
    pub fn next_basis(&self) -> &[Id] {
        &self.next_basis
    }

    /// Apply the selected operator without constructing its matrix or running PH.
    ///
    /// x uses [`Self::basis`] coordinates and the standard real inner product.
    /// The result owns coordinates in that same basis. A computed empty degree
    /// accepts an empty vector and returns an empty vector.
    /// # Errors
    /// Rejects wrong length, nonfinite input/output, allocation failure, cancellation
    /// or work exhaustion. Each call starts a fresh budget; no partial output escapes.
    pub fn action(&self, x: &[f64], execution: &Execution<'_>) -> Result<Vec<f64>> {
        let mut budget = WorkBudget::new(execution)?;
        self.action_budget(x, &mut budget)
    }

    fn action_budget(&self, x: &[f64], budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
        if x.len() != self.basis.len() {
            return Err(Error::ShapeMismatch {
                input: "Hodge vector",
                expected: self.basis.len(),
                actual: x.len(),
            });
        }
        for (index, value) in x.iter().enumerate() {
            budget.step()?;
            if !value.is_finite() {
                return Err(Error::NonFiniteValue {
                    field: "Hodge vector",
                    index: Some(index),
                });
            }
        }
        let mut y = zeros(x.len(), budget)?;
        if self.selection != HodgeSelection::Upper {
            let mut ax = zeros(self.previous_basis.len(), budget)?;
            for (column, terms) in self.lower.iter().enumerate() {
                for &(row, coefficient) in terms {
                    budget.step()?;
                    ax[row] += coefficient as f64 * x[column];
                }
            }
            for (column, terms) in self.lower.iter().enumerate() {
                for &(row, coefficient) in terms {
                    budget.step()?;
                    y[column] += coefficient as f64 * ax[row];
                }
            }
        }
        if self.selection != HodgeSelection::Lower {
            for terms in &self.upper {
                let mut dot = 0.0;
                for &(row, coefficient) in terms {
                    budget.step()?;
                    dot += coefficient as f64 * x[row];
                }
                for &(row, coefficient) in terms {
                    budget.step()?;
                    y[row] += coefficient as f64 * dot;
                }
            }
        }
        for value in &mut y {
            budget.step()?;
            finite(*value)?;
            *value = crate::canonical_zero(*value);
        }
        budget.check()?;
        Ok(y)
    }

    /// Materialize the n by n operator in row-major [`Self::basis`] order.
    ///
    /// An empty degree returns length zero with n=0. This explicit conversion
    /// rejects n² above [`DEFAULT_MAX_DENSE_ENTRIES`] before dense allocation,
    /// including with unlimited Execution. Use [`Self::matrix_with_limit`] to
    /// explicitly adjust the entry limit. Work is O(n * (n + sparse incidence));
    /// Execution controls work/cancellation, not aggregate memory.
    /// # Errors
    /// Includes excessive dense entries, shape overflow, allocation/numerical
    /// failure and execution limits.
    pub fn matrix(&self, execution: &Execution<'_>) -> Result<Vec<f64>> {
        self.matrix_with_limit(DEFAULT_MAX_DENSE_ENTRIES, execution)
    }

    /// Materialize with an explicit maximum number of n by n output entries.
    ///
    /// Checks n² against max_dense_entries before allocating any dense buffer.
    /// Zero permits only an empty degree. Trusted callers can raise the limit;
    /// this per-buffer bound is not a total memory quota. Sparse action is
    /// unaffected, and allocation/byte-capacity checks still apply.
    /// # Errors
    /// Includes [`Self::matrix`] errors. Excess entries return
    /// [`Error::InvalidParameter`] with parameter `Hodge dense entries`.
    pub fn matrix_with_limit(
        &self,
        max_dense_entries: usize,
        execution: &Execution<'_>,
    ) -> Result<Vec<f64>> {
        let mut budget = WorkBudget::new(execution)?;
        self.matrix_budget(max_dense_entries, &mut budget)
    }

    fn matrix_budget(
        &self,
        max_dense_entries: usize,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Vec<f64>> {
        let n = self.basis.len();
        let entries = n.checked_mul(n).ok_or(size_error())?;
        if entries > max_dense_entries {
            return Err(Error::InvalidParameter {
                parameter: "Hodge dense entries",
                reason: "exceeds the explicit dense entry limit",
            });
        }
        let mut matrix = zeros(entries, budget)?;
        let mut x = zeros(n, budget)?;
        for column in 0..n {
            budget.step()?;
            x[column] = 1.0;
            let y = self.action_budget(&x, budget)?;
            for (row, value) in y.into_iter().enumerate() {
                budget.step()?;
                matrix[row * n + column] = value;
            }
            x[column] = 0.0;
        }
        budget.check()?;
        Ok(matrix)
    }
}

fn read_column<C: FilteredComplex>(
    source: &C,
    cell: C::CellId,
    cells: &HashMap<C::CellId, (usize, usize, usize)>,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<(usize, i64)>> {
    let &(position, degree, _) = &cells[&cell];
    let mut terms: HashMap<usize, i64> = HashMap::new();
    for (face, coefficient) in source.boundary(cell) {
        budget.step()?;
        if coefficient == 0 {
            continue;
        }
        let &(face_position, face_degree, row) = cells
            .get(&face)
            .ok_or(invalid(Some(position), "unknown boundary cell"))?;
        if degree.checked_sub(1) != Some(face_degree) || face_position >= position {
            return Err(invalid(
                Some(position),
                "boundary has wrong degree or traversal order",
            ));
        }
        terms.try_reserve(1).map_err(|_| allocation())?;
        let entry = terms.entry(row).or_default();
        *entry = entry
            .checked_add(i64::from(coefficient))
            .ok_or(size_error())?;
    }
    let mut column = Vec::new();
    column
        .try_reserve_exact(terms.len())
        .map_err(|_| allocation())?;
    for (row, coefficient) in terms {
        budget.step()?;
        // Exact integer-to-f64 conversion is part of this initial arithmetic contract.
        if coefficient.unsigned_abs() > (1_u64 << 53) {
            return Err(invalid(
                Some(position),
                "incidence exceeds exact f64 integer range",
            ));
        }
        if coefficient != 0 {
            column.push((row, coefficient));
        }
    }
    budget.check()?;
    column.sort_unstable_by_key(|&(row, _)| row);
    budget.check()?;
    Ok(column)
}

fn push<T>(values: &mut Vec<T>, value: T) -> Result<()> {
    values.try_reserve(1).map_err(|_| allocation())?;
    values.push(value);
    Ok(())
}
fn zeros(n: usize, budget: &mut WorkBudget<'_>) -> Result<Vec<f64>> {
    n.checked_mul(std::mem::size_of::<f64>())
        .filter(|&bytes| bytes <= isize::MAX as usize)
        .ok_or(size_error())?;
    budget.step_by(n)?;
    let mut values = Vec::new();
    values.try_reserve_exact(n).map_err(|_| allocation())?;
    values.resize(n, 0.0);
    budget.check()?;
    Ok(values)
}
fn finite(value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(Error::NumericalFailure {
            context: "real Hodge arithmetic",
        })
    }
}
fn invalid(cell: Option<usize>, reason: &'static str) -> Error {
    Error::InvalidComplex { cell, reason }
}
fn allocation() -> Error {
    Error::AllocationFailed {
        context: "Hodge window or result",
    }
}
fn size_error() -> Error {
    Error::SizeOverflow {
        operation: "Hodge shape or integer incidence",
    }
}
