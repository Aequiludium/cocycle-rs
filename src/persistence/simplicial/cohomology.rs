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
    let level = next_dimension(access, &access.vertices()?, &mut || budget.step())?;
    let mut forest = UnionFind::new(access.vertex_count())?;
    let mut cleared = HashSet::new();
    cleared.try_reserve(level.len()).map_err(|_| allocation())?;
    let mut raw = Vec::new();
    raw.try_reserve(access.vertex_count())
        .map_err(|_| allocation())?;
    for edge in &level {
        budget.step()?;
        #[cfg(test)]
        profiling::record(0, |s| s.processed += 1);
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
        1..=max_dimension,
        field,
        raw,
        level,
        cleared,
        budget,
    )
}

/// Continue exact F2 flag persistence without redoing H0 or H1 reduction.
/// Death keys own original vertices; topology enumeration still includes them.
pub(in crate::persistence) fn continue_from_h1(
    access: &impl ZeroBornSimplicialAccess,
    max_dimension: usize,
    raw: RawIntervals,
    deaths: Vec<[usize; 3]>,
    budget: &mut WorkBudget<'_>,
) -> Result<RawIntervals> {
    budget.check()?;
    let mut cleared = HashSet::new();
    cleared
        .try_reserve(deaths.len())
        .map_err(|_| allocation())?;
    #[cfg(test)]
    let handoff_capacity_bytes = deaths.capacity() * std::mem::size_of::<[usize; 3]>();
    #[cfg(test)]
    let handoff_len = deaths.len();
    for triangle in deaths {
        budget.step()?;
        let mut key = Vec::new();
        key.try_reserve_exact(3).map_err(|_| allocation())?;
        key.extend(triangle);
        cleared.insert(key);
        #[cfg(test)]
        if cleared.len() == handoff_len {
            profiling::workspace_event(
                "handoff_conversion_overlap",
                2,
                &[
                    ("handoff", handoff_capacity_bytes),
                    (
                        "clearing_keys",
                        cleared
                            .iter()
                            .map(|key| key.capacity() * std::mem::size_of::<usize>())
                            .sum(),
                    ),
                ],
                &[
                    ("cleared", cleared.len()),
                    ("clearing_slots", cleared.capacity()),
                ],
            )?;
        }
    }
    #[cfg(test)]
    profiling::workspace_level(
        "handoff_converted",
        2,
        &[],
        &cleared,
        &Vec::new(),
        &HashMap::new(),
        &raw,
    )?;
    let edges = next_dimension(access, &access.vertices()?, &mut || budget.step())?;
    let triangles = next_dimension(access, &edges, &mut || budget.step())?;
    #[cfg(test)]
    profiling::workspace_level(
        "triangle_assembly",
        2,
        &[("edges", &edges), ("triangles", &triangles)],
        &cleared,
        &Vec::new(),
        &HashMap::new(),
        &raw,
    )?;
    drop(edges);
    reduce_dimensions(
        access,
        2..=max_dimension,
        PrimeField::new(2)?,
        raw,
        triangles,
        cleared,
        budget,
    )
}

fn reduce_dimensions(
    access: &impl ZeroBornSimplicialAccess,
    dimensions: std::ops::RangeInclusive<usize>,
    field: PrimeField,
    mut raw: RawIntervals,
    mut level: Vec<Simplex>,
    mut cleared: HashSet<Vec<usize>>,
    budget: &mut WorkBudget<'_>,
) -> Result<RawIntervals> {
    let max_dimension = *dimensions.end();
    for dimension in dimensions.take_while(|&d| d < access.vertex_count()) {
        if level.is_empty() {
            break;
        }
        let mut owners: HashMap<Vec<usize>, usize> = HashMap::new();
        let mut columns: Vec<Column<usize>> = Vec::new();
        #[cfg(test)]
        profiling::workspace_level(
            "reduction_start",
            dimension,
            &[("level", &level)],
            &cleared,
            &columns,
            &owners,
            &raw,
        )?;
        #[cfg(test)]
        profiling::record(dimension, |s| s.simplices = level.len());
        for position in (0..level.len()).rev() {
            budget.step()?;
            let simplex = &level[position];
            if cleared.contains(&simplex.vertices) {
                #[cfg(test)]
                profiling::record(dimension, |s| s.cleared += 1);
                continue;
            }
            #[cfg(test)]
            profiling::record(dimension, |s| s.processed += 1);
            let mut working = Column::new();
            append(access, simplex, 1, field, &mut working, budget)?;
            let mut transform = Column::unit(position);
            loop {
                budget.step()?;
                let Some((pivot, &coefficient)) = working.first() else {
                    push_interval(&mut raw, (dimension, simplex.value, None))?;
                    break;
                };
                #[cfg(test)]
                profiling::record(dimension, |s| s.pivot_lookups += 1);
                if let Some(&owner) = owners.get(&pivot.vertices) {
                    #[cfg(test)]
                    profiling::record(dimension, |s| {
                        s.pivot_hits += 1;
                        s.column_additions += 1;
                        s.reconstructions += 1;
                    });
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
                    #[cfg(test)]
                    profiling::record(dimension, |s| {
                        s.stored_transform_entries += transform.entries().count();
                        s.largest_transform = s.largest_transform.max(transform.entries().count());
                    });
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
        #[cfg(test)]
        profiling::workspace_level(
            "reduction_end",
            dimension,
            &[("level", &level)],
            &cleared,
            &columns,
            &owners,
            &raw,
        )?;
        budget.check()?;
        if dimension == max_dimension {
            break;
        }
        cleared.clear();
        cleared
            .try_reserve(owners.len())
            .map_err(|_| allocation())?;
        cleared.extend(owners.into_keys());
        #[cfg(test)]
        profiling::workspace_level(
            "clear_next_extracted",
            dimension,
            &[("level", &level)],
            &cleared,
            &columns,
            &HashMap::new(),
            &raw,
        )?;
        // Cleared simplices still participate in clique enumeration: clearing
        // removes reduction columns, never topology needed by the next level.
        let next_level = next_dimension(access, &level, &mut || budget.step())?;
        #[cfg(test)]
        profiling::workspace_level(
            "next_level_assembled",
            dimension,
            &[("old_level", &level), ("next_level", &next_level)],
            &cleared,
            &columns,
            &HashMap::new(),
            &raw,
        )?;
        level = next_level;
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
        #[cfg(test)]
        profiling::record(simplex.dimension(), |s| {
            s.peak_working_column = s.peak_working_column.max(working.entries().count())
        });
        Ok(())
    })
}

#[cfg(test)]
mod profiling;
#[cfg(test)]
pub(in crate::persistence) use profiling::take_counts;
#[cfg(test)]
pub(in crate::persistence) use profiling::{
    workspace_event, workspace_force_failure, workspace_heap_growth,
};
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
