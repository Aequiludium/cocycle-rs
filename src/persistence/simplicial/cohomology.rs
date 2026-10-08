//! Zero-born, dimension-generic implicit prime-field coboundary reduction with clearing.
//!
//! Only the current simplex dimension, pivot ownership, transformation columns
//! and one working coboundary are retained. Reduced coboundaries are regenerated
//! from transformations instead of stored. Higher-dimensional simplices are
//! visited on demand through the shared filtration access contract.
use crate::algebra::{PrimeField, column::Column};
use crate::complex::Simplex;
use crate::filtration::simplicial::{ZeroBornSimplicialAccess, next_dimension};
use crate::persistence::execution::WorkBudget;
use crate::persistence::{RawIntervals, union_find::UnionFind};
use crate::{Error, Result};
use std::collections::{HashMap, HashSet};

pub(in crate::persistence) fn compute(
    access: &impl ZeroBornSimplicialAccess,
    max_dimension: usize,
    field: PrimeField,
    budget: &mut WorkBudget<'_>,
) -> Result<RawIntervals> {
    budget.check()?;
    #[cfg(test)]
    record_dimension(0);
    let level = next_dimension(access, &access.vertices()?, &mut || budget.step())?;
    let mut forest = UnionFind::new(access.vertex_count())?;
    let mut cleared = HashSet::new();
    cleared.try_reserve(level.len()).map_err(|_| allocation())?;
    let mut raw = Vec::new();
    raw.try_reserve(access.vertex_count())
        .map_err(|_| allocation())?;
    for edge in &level {
        budget.step()?;
        if forest.merge(
            access.vertex_position(edge.vertices[0]),
            access.vertex_position(edge.vertices[1]),
        ) {
            raw.push((0, 0.0, Some(edge.value)));
            cleared.insert(edge.vertices.clone());
        }
    }
    raw.extend((0..forest.components()).map(|_| (0, 0.0, None)));
    reduce_dimensions(
        access,
        ReductionStart {
            dimension: 1,
            level,
            cleared,
            intervals: raw,
        },
        max_dimension,
        field,
        budget,
    )
}

/// Reduce only `min_dimension..=max_dimension`. The initial clearing set must
/// contain *all* pivots of the preceding dimension under this access and order,
/// including zero-lifetime pairs. An empty set is valid only if that rank is zero.
/// Enumerating lower levels constructs topology; it does not reduce H0/H1 again.
pub(in crate::persistence) fn compute_range(
    access: &impl ZeroBornSimplicialAccess,
    min_dimension: usize,
    max_dimension: usize,
    field: PrimeField,
    cleared: HashSet<Vec<usize>>,
    budget: &mut WorkBudget<'_>,
) -> Result<RawIntervals> {
    budget.check()?;
    if min_dimension == 0 {
        return Err(Error::InternalInvariant {
            reason: "cohomology continuation must start above H0",
        });
    }
    if min_dimension > max_dimension || min_dimension >= access.vertex_count() {
        return Ok(Vec::new());
    }
    let mut level = access.vertices()?;
    for _ in 0..min_dimension {
        level = next_dimension(access, &level, &mut || budget.step())?;
        if level.is_empty() {
            return Ok(Vec::new());
        }
    }
    reduce_dimensions(
        access,
        ReductionStart {
            dimension: min_dimension,
            level,
            cleared,
            intervals: Vec::new(),
        },
        max_dimension,
        field,
        budget,
    )
}

struct ReductionStart {
    dimension: usize,
    level: Vec<Simplex>,
    cleared: HashSet<Vec<usize>>,
    intervals: RawIntervals,
}

fn reduce_dimensions(
    access: &impl ZeroBornSimplicialAccess,
    start: ReductionStart,
    max_dimension: usize,
    field: PrimeField,
    budget: &mut WorkBudget<'_>,
) -> Result<RawIntervals> {
    let ReductionStart {
        dimension: first_dimension,
        mut level,
        mut cleared,
        intervals: mut raw,
    } = start;
    for dimension in first_dimension..=max_dimension.min(access.vertex_count().saturating_sub(1)) {
        if level.is_empty() {
            break;
        }
        #[cfg(test)]
        record_dimension(dimension);
        let mut owners: HashMap<Vec<usize>, usize> = HashMap::new();
        let mut columns: Vec<Column<usize>> = Vec::new();
        for position in (0..level.len()).rev() {
            budget.step()?;
            let simplex = &level[position];
            if cleared.contains(&simplex.vertices) {
                continue;
            }
            let mut working = Column::new();
            append(access, simplex, 1, field, &mut working, budget)?;
            let mut transform = Column::unit(position);
            loop {
                budget.step()?;
                let Some((pivot, &coefficient)) = working.first() else {
                    push_interval(&mut raw, (dimension, simplex.value, None))?;
                    break;
                };
                if let Some(&owner) = owners.get(&pivot.vertices) {
                    let factor = field.negate(coefficient);
                    for (&source, &value) in columns[owner].entries() {
                        budget.step()?;
                        append(
                            access,
                            &level[source],
                            field.multiply(factor, value),
                            field,
                            &mut working,
                            budget,
                        )?;
                    }
                    transform.add_scaled(&columns[owner], factor, field, &mut || budget.step())?;
                } else {
                    push_interval(&mut raw, (dimension, simplex.value, Some(pivot.value)))?;
                    owners.try_reserve(1).map_err(|_| allocation())?;
                    columns.try_reserve(1).map_err(|_| allocation())?;
                    transform.scale(field.inverse(coefficient)?, field, &mut || budget.step())?;
                    owners.insert(pivot.vertices.clone(), columns.len());
                    columns.push(transform);
                    break;
                }
            }
        }
        budget.check()?;
        if dimension == max_dimension {
            break;
        }
        cleared.clear();
        cleared
            .try_reserve(owners.len())
            .map_err(|_| allocation())?;
        cleared.extend(owners.into_keys());
        // Cleared simplices still participate in clique enumeration: clearing
        // removes reduction columns, never topology needed by the next level.
        level = next_dimension(access, &level, &mut || budget.step())?;
    }
    budget.check()?;
    Ok(raw)
}
fn append(
    access: &impl ZeroBornSimplicialAccess,
    simplex: &Simplex,
    factor: u32,
    field: PrimeField,
    working: &mut Column<Simplex>,
    budget: &mut WorkBudget<'_>,
) -> Result<()> {
    access.visit_cofacets(simplex, false, &mut || budget.step(), |row| {
        // The added vertex is omitted to recover this oriented facet.
        let omitted = row
            .vertices
            .iter()
            .position(|v| simplex.vertices.binary_search(v).is_err())
            .unwrap();
        let coefficient = field.multiply(factor, field.orientation(omitted));
        working.add_term(row, coefficient, field);
        Ok(())
    })
}
fn push_interval(raw: &mut RawIntervals, interval: (usize, f64, Option<f64>)) -> Result<()> {
    raw.try_reserve(1).map_err(|_| allocation())?;
    raw.push(interval);
    Ok(())
}
fn allocation() -> Error {
    Error::AllocationFailed {
        context: "dimension-generic cohomology",
    }
}

// Per-thread diagnostics verify backend selection without production counters.
#[cfg(test)]
std::thread_local! {
    static REDUCED_DIMENSIONS: std::cell::RefCell<Vec<usize>> = const { std::cell::RefCell::new(Vec::new()) };
}
#[cfg(test)]
fn record_dimension(dimension: usize) {
    REDUCED_DIMENSIONS.with(|trace| trace.borrow_mut().push(dimension));
}
#[cfg(test)]
pub(in crate::persistence) fn take_reduced_dimensions() -> Vec<usize> {
    REDUCED_DIMENSIONS.with(|trace| std::mem::take(&mut *trace.borrow_mut()))
}
