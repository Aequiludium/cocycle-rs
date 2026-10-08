//! F2 Rips H1 via implicit coboundaries and stored change-of-basis columns.
//! Independently implemented from the invariants in docs/reference/mathematics.md section 9;
//! the explicit boundary reducer remains an independent test oracle.

use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

#[cfg(test)]
use super::clearing::IgnoreTriangles;
use super::clearing::TriangleClearing;
use super::edges::{ColumnPosition, EdgePosition};
#[cfg(test)]
use crate::filtration::flag::DenseFlag;
use crate::filtration::flag::{FlagAccess, SimplexEntry};
#[cfg(test)]
use crate::filtration::rips::cone_radius;
#[cfg(test)]
use crate::geometry::DissimilarityView;
use crate::persistence::execution::WorkBudget;
use crate::{Error, Result};

use crate::persistence::RawIntervals;
type Coboundary = BinaryHeap<Reverse<WorkingEntry>>;

// Working rows carry canonical nonnegative finite values, stored as their exact
// f64 bits. Integer equality preserves F2 cancellation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
struct WorkingEntry {
    id: usize,
    value_bits: u64,
}
impl From<SimplexEntry> for WorkingEntry {
    fn from(row: SimplexEntry) -> Self {
        debug_assert!(row.value.is_finite() && !row.value.is_sign_negative());
        Self {
            id: row.id,
            value_bits: row.value.to_bits(),
        }
    }
}
impl WorkingEntry {
    fn simplex(self) -> SimplexEntry {
        SimplexEntry {
            id: self.id,
            value: f64::from_bits(self.value_bits),
        }
    }
}
impl Ord for WorkingEntry {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        self.simplex().cmp(&other.simplex())
    }
}
impl PartialOrd for WorkingEntry {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
    fn lt(&self, other: &Self) -> bool {
        crate::complex::finite_filtration_le(
            f64::from_bits(self.value_bits),
            f64::from_bits(other.value_bits),
            self.id > other.id,
        )
    }
    fn le(&self, other: &Self) -> bool {
        crate::complex::finite_filtration_le(
            f64::from_bits(self.value_bits),
            f64::from_bits(other.value_bits),
            self.id >= other.id,
        )
    }
    fn gt(&self, other: &Self) -> bool {
        other.lt(self)
    }
    fn ge(&self, other: &Self) -> bool {
        other.le(self)
    }
}

const NO_SHORTCUTS: u8 = 0;
const APPARENT: u8 = 1;
const EMERGENT: u8 = 2;
// Omit zero apparent columns and reconstruct their pivots during reduction.
const VIRTUAL_APPARENT: u8 = 4;
const APPARENT_EMERGENT: u8 = APPARENT | EMERGENT;
const ALL_SHORTCUTS: u8 = APPARENT_EMERGENT | VIRTUAL_APPARENT;
const PRODUCTION_SHORTCUTS: u8 = ALL_SHORTCUTS;

struct TransformColumn {
    // The diagonal entry of V is implicit. Other entries index later edges in
    // the forward filtration (earlier columns in the reverse computation).
    edge: EdgePosition,
    additions: Vec<EdgePosition>,
    #[cfg(test)]
    reduced: Vec<SimplexEntry>,
}

/// Instrumentation is compiled out of production builds.
#[derive(Default, Debug)]
struct Stats {
    #[cfg(test)]
    two_pass_initialization: bool,
    #[cfg(test)]
    edges: usize,
    #[cfg(test)]
    cofacets: usize,
    #[cfg(test)]
    initial_candidates: usize,
    #[cfg(test)]
    reconstruction_candidates: usize,
    #[cfg(test)]
    apparent_candidates: usize,
    #[cfg(test)]
    skipped_apparent: usize,
    #[cfg(test)]
    virtual_additions: usize,
    #[cfg(test)]
    stored_columns: usize,
    #[cfg(test)]
    verify_transforms: bool,
    #[cfg(test)]
    checked_transforms: usize,
    #[cfg(test)]
    column_additions: usize,
    #[cfg(test)]
    shortcuts: usize,
    #[cfg(test)]
    peak_heap: usize,
    #[cfg(test)]
    stored_entries: usize,
    #[cfg(test)]
    largest_transform: usize,
    #[cfg(test)]
    peak_transform_heap: usize,
}

/// Retain every H1 death triangle for a subsequent H2 reduction.
pub(super) fn compute_with_clearing<const CONTROLLED: bool>(
    rips: &impl FlagAccess,
    edges: Vec<SimplexEntry>,
    budget: &mut WorkBudget<'_, CONTROLLED>,
    clearing: &mut impl TriangleClearing,
) -> Result<RawIntervals> {
    reduce_edges::<true, true, PRODUCTION_SHORTCUTS, CONTROLLED>(
        rips,
        edges,
        &mut Stats::default(),
        budget,
        clearing,
    )
}

// Retain independent optimization configurations for the dense test oracle.
#[cfg(test)]
fn run<const IMPLICIT: bool, const CLEAR: bool, const CONE: bool, const SHORTCUTS: u8>(
    input: DissimilarityView<'_>,
    cutoff: f64,
    stats: &mut Stats,
) -> Result<RawIntervals> {
    let mut budget = WorkBudget::new(&crate::persistence::ExecutionLimits::default())?;
    let stop = if CONE {
        cutoff.min(cone_radius(input.into(), &mut || budget.step())?)
    } else {
        cutoff
    };
    run_access::<IMPLICIT, CLEAR, SHORTCUTS, true>(
        &DenseFlag::new(input.into(), stop)?,
        stats,
        &mut budget,
    )
}

#[cfg(test)]
fn run_access<
    const IMPLICIT: bool,
    const CLEAR: bool,
    const SHORTCUTS: u8,
    const CONTROLLED: bool,
>(
    rips: &impl FlagAccess,
    stats: &mut Stats,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<RawIntervals> {
    let edges = rips.edges(&mut || budget.step())?;
    reduce_edges::<IMPLICIT, CLEAR, SHORTCUTS, CONTROLLED>(
        rips,
        edges,
        stats,
        budget,
        &mut IgnoreTriangles,
    )
}

fn reduce_edges<
    const IMPLICIT: bool,
    const CLEAR: bool,
    const SHORTCUTS: u8,
    const CONTROLLED: bool,
>(
    rips: &impl FlagAccess,
    edges: Vec<SimplexEntry>,
    stats: &mut Stats,
    budget: &mut WorkBudget<'_, CONTROLLED>,
    clearing: &mut impl TriangleClearing,
) -> Result<RawIntervals> {
    budget.check()?;
    #[cfg(test)]
    {
        stats.edges = edges.len();
    }
    let super::edges::EdgePreparation {
        cycles: cycle_edges,
        intervals: mut raw,
    } = super::edges::classify(rips, &edges, budget)?;

    // Map triangle ids to stored transformation columns; edge and
    // additions in those columns are positions in edges, not simplex ids.
    let mut pivot_owners: HashMap<usize, ColumnPosition> = HashMap::new();
    let mut columns: Vec<TransformColumn> = Vec::new();
    let mut working = Coboundary::new();
    let mut scratch = Vec::new();
    let mut transform = BinaryHeap::new();
    for j in (0..edges.len()).rev() {
        budget.step()?;
        if CLEAR && !cycle_edges[j] {
            continue;
        }
        working.clear();
        transform.clear();
        let edge = edges[j];
        let (shortcut, apparent_pair) = initialize_coboundary::<SHORTCUTS, CONTROLLED>(
            rips,
            edge,
            &pivot_owners,
            &mut working,
            stats,
            budget,
        )?;
        if apparent_pair {
            if !cycle_edges[j] {
                return Err(Error::InternalInvariant {
                    reason: "H0 death edge in zero apparent pair",
                });
            }
            #[cfg(test)]
            {
                stats.skipped_apparent += 1;
            }
            clearing.record(shortcut.ok_or(Error::InternalInvariant {
                reason: "apparent pair missing its triangle",
            })?)?;
            continue;
        }
        let pivot = if let Some(pivot) = shortcut {
            Some(pivot)
        } else {
            loop {
                budget.step()?;
                let Some(Reverse(entry)) = pop_parity(&mut working, budget)? else {
                    break None;
                };
                let pivot = entry.simplex();
                let Some(&owner) = pivot_owners.get(&pivot.id) else {
                    if SHORTCUTS & VIRTUAL_APPARENT != 0
                        && let Some(edge) = zero_apparent_facet(rips, pivot, stats, budget)?
                    {
                        let position =
                            edges
                                .binary_search(&edge)
                                .map_err(|_| Error::InternalInvariant {
                                    reason: "zero apparent facet missing from edges",
                                })?;
                        // Its original coboundary has this pivot as its first
                        // row and its column precedes ours in reverse order.
                        if position <= j || !cycle_edges[position] {
                            return Err(Error::InternalInvariant {
                                reason: "zero apparent facet violates reduction order",
                            });
                        }
                        push_heap(&mut working, Reverse(pivot.into()))?;
                        append_coboundary(rips, edge, &mut working, &mut scratch, stats, budget)?;
                        push_heap(&mut transform, EdgePosition(position))?;
                        #[cfg(test)]
                        {
                            stats.virtual_additions += 1;
                            stats.peak_transform_heap =
                                stats.peak_transform_heap.max(transform.len());
                        }
                        continue;
                    }
                    break Some(pivot);
                };
                // Put the pivot back: adding the owner's column cancels it.
                push_heap(&mut working, Reverse(pivot.into()))?;
                let column = &columns[owner.0];
                #[cfg(test)]
                {
                    stats.column_additions += 1;
                }
                if IMPLICIT {
                    append_coboundary(
                        rips,
                        edges[column.edge.0],
                        &mut working,
                        &mut scratch,
                        stats,
                        budget,
                    )?;
                    push_heap(&mut transform, column.edge)?;
                    for &k in &column.additions {
                        append_coboundary(
                            rips,
                            edges[k.0],
                            &mut working,
                            &mut scratch,
                            stats,
                            budget,
                        )?;
                        push_heap(&mut transform, k)?;
                    }
                } else {
                    #[cfg(test)]
                    for &row in &column.reduced {
                        push_heap(&mut working, Reverse(row.into()))?;
                    }
                }
                #[cfg(test)]
                {
                    stats.peak_heap = stats.peak_heap.max(working.len());
                    stats.peak_transform_heap = stats.peak_transform_heap.max(transform.len());
                }
            }
        };
        if let Some(pivot) = pivot {
            clearing.record(pivot)?;
            if !cycle_edges[j] {
                return Err(Error::InternalInvariant {
                    reason: "H0 death edge paired in H1",
                });
            }
            let mut additions = Vec::new();
            while let Some(k) = pop_parity(&mut transform, budget)? {
                additions
                    .try_reserve(1)
                    .map_err(|_| allocation("Rips transform column"))?;
                additions.push(k);
            }
            #[cfg(test)]
            if IMPLICIT && stats.verify_transforms && shortcut.is_none() {
                tests::check_transform(rips, &edges, edge, &additions, &working, pivot);
                stats.checked_transforms += 1;
            }
            #[cfg(test)]
            let mut reduced = Vec::new();
            #[cfg(test)]
            if !IMPLICIT {
                if shortcut.is_some() {
                    append_coboundary(rips, edge, &mut working, &mut scratch, stats, budget)?;
                } else {
                    push_heap(&mut working, Reverse(pivot.into()))?;
                }
                while let Some(Reverse(entry)) = pop_parity(&mut working, budget)? {
                    let row = entry.simplex();
                    reduced.push(row);
                }
            }
            #[cfg(test)]
            {
                stats.stored_entries += additions.len() + reduced.len();
                stats.largest_transform = stats.largest_transform.max(additions.len());
            }
            pivot_owners
                .try_reserve(1)
                .map_err(|_| allocation("Rips pivot owners"))?;
            columns
                .try_reserve(1)
                .map_err(|_| allocation("Rips transform columns"))?;
            pivot_owners.insert(pivot.id, ColumnPosition(columns.len()));
            columns.push(TransformColumn {
                edge: EdgePosition(j),
                additions,
                #[cfg(test)]
                reduced,
            });
            #[cfg(test)]
            {
                stats.stored_columns += 1;
            }
            // Keep the pivot and transformation even when this public bar is
            // empty: later columns can still need them for cancellation.
            if edge.value != pivot.value {
                raw.try_reserve(1)
                    .map_err(|_| allocation("Rips intervals"))?;
                raw.push((1, edge.value, Some(pivot.value)));
            }
        } else if cycle_edges[j] {
            raw.try_reserve(1)
                .map_err(|_| allocation("Rips intervals"))?;
            raw.push((1, edge.value, None));
        }
    }
    Ok(raw)
}

/// Initialize the original column, retaining independently testable strategies.
fn initialize_coboundary<const SHORTCUTS: u8, const CONTROLLED: bool>(
    rips: &impl FlagAccess,
    edge: SimplexEntry,
    pivot_owners: &HashMap<usize, ColumnPosition>,
    working: &mut Coboundary,
    stats: &mut Stats,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<(Option<SimplexEntry>, bool)> {
    #[cfg(test)]
    if stats.two_pass_initialization {
        return initialize_two_pass::<SHORTCUTS, CONTROLLED>(
            rips,
            edge,
            pivot_owners,
            working,
            stats,
            budget,
        );
    }
    // Recover the existing heap allocation as an unsorted scratch buffer.
    let mut rows = std::mem::take(working).into_vec();
    rows.clear();
    let mut found = None;
    let mut apparent_pair = false;
    let mut check_shortcut = SHORTCUTS != NO_SHORTCUTS;
    #[cfg(test)]
    let mut candidates = 0;
    // The access callbacks execute sequentially. An eligible cofacet can
    // require a nested apparent-facet scan using the same execution budget.
    let budget = std::cell::RefCell::new(budget);
    // Both access implementations visit triangles in decreasing combinatorial
    // ID order. The first equal-valued cofacet is the original column pivot.
    rips.visit_cofacets(
        edge,
        &mut || {
            #[cfg(test)]
            {
                candidates += 1;
            }
            budget.borrow_mut().step()
        },
        |triangle| {
            count_cofacet(stats);
            if check_shortcut && triangle.value == edge.value {
                // An occupied first equal-valued cofacet requires ordinary
                // reduction; never try a later equal-valued cofacet instead.
                check_shortcut = false;
                if SHORTCUTS & VIRTUAL_APPARENT != 0 && rips.latest_facet(triangle) == edge {
                    if pivot_owners.contains_key(&triangle.id) {
                        return Err(Error::InternalInvariant {
                            reason: "zero apparent pivot has an ordinary owner",
                        });
                    }
                    found = Some(triangle);
                    apparent_pair = true;
                    return Ok(false);
                }
                let eligible = SHORTCUTS & EMERGENT != 0
                    || (SHORTCUTS & APPARENT != 0 && rips.latest_facet(triangle) == edge);
                if eligible
                    && !pivot_owners.contains_key(&triangle.id)
                    && (SHORTCUTS & VIRTUAL_APPARENT == 0
                        || zero_apparent_facet(rips, triangle, stats, &mut budget.borrow_mut())?
                            .is_none())
                {
                    #[cfg(test)]
                    {
                        stats.shortcuts += 1;
                    }
                    found = Some(triangle);
                    return Ok(false);
                }
            }
            rows.try_reserve(1)
                .map_err(|_| allocation("Rips working heap"))?;
            rows.push(Reverse(triangle.into()));
            Ok(true)
        },
    )?;
    let budget = budget.into_inner();
    #[cfg(test)]
    {
        stats.initial_candidates += candidates;
    }
    if found.is_some() {
        rows.clear();
    }
    // Heap construction, like sorting, is checked at its phase boundaries.
    budget.check()?;
    *working = BinaryHeap::from(rows);
    budget.check()?;
    #[cfg(test)]
    {
        stats.peak_heap = stats.peak_heap.max(working.len());
    }
    Ok((found, apparent_pair))
}

/// Probe without caching a prefix; enumerate the full column on failure.
/// Omission stays independent so tests can compare all four combinations.
#[cfg(test)]
fn initialize_two_pass<const SHORTCUTS: u8, const CONTROLLED: bool>(
    rips: &impl FlagAccess,
    edge: SimplexEntry,
    owners: &HashMap<usize, ColumnPosition>,
    working: &mut Coboundary,
    stats: &mut Stats,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<(Option<SimplexEntry>, bool)> {
    working.clear();
    let mut found = None;
    let mut omitted = false;
    let mut candidates = 0;
    let budget = std::cell::RefCell::new(budget);
    if SHORTCUTS != NO_SHORTCUTS {
        rips.visit_cofacets(
            edge,
            &mut || {
                candidates += 1;
                budget.borrow_mut().step()
            },
            |triangle| {
                count_cofacet(stats);
                if triangle.value != edge.value {
                    return Ok(true);
                }
                if SHORTCUTS & VIRTUAL_APPARENT != 0 && rips.latest_facet(triangle) == edge {
                    if owners.contains_key(&triangle.id) {
                        return Err(Error::InternalInvariant {
                            reason: "zero apparent pivot has an ordinary owner",
                        });
                    }
                    found = Some(triangle);
                    omitted = true;
                } else if !owners.contains_key(&triangle.id)
                    && (SHORTCUTS & EMERGENT != 0
                        || (SHORTCUTS & APPARENT != 0 && rips.latest_facet(triangle) == edge))
                    && (SHORTCUTS & VIRTUAL_APPARENT == 0
                        || zero_apparent_facet(rips, triangle, stats, &mut budget.borrow_mut())?
                            .is_none())
                {
                    found = Some(triangle);
                    stats.shortcuts += 1;
                }
                // Stop the probe even when the first equal candidate is owned.
                Ok(false)
            },
        )?;
    }
    let budget = budget.into_inner();
    if found.is_none() {
        rips.visit_cofacets(
            edge,
            &mut || {
                candidates += 1;
                budget.step()
            },
            |triangle| {
                count_cofacet(stats);
                push_heap(working, Reverse(triangle.into()))?;
                Ok(true)
            },
        )?;
    }
    stats.initial_candidates += candidates;
    stats.peak_heap = stats.peak_heap.max(working.len());
    budget.check()?;
    Ok((found, omitted))
}

/// A zero apparent pair is determined by mutual earliest-cofacet/latest-facet
/// tests, not by equal filtration values alone. No ownership is stored for it.
fn zero_apparent_facet<const CONTROLLED: bool>(
    rips: &impl FlagAccess,
    triangle: SimplexEntry,
    stats: &mut Stats,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<Option<SimplexEntry>> {
    let edge = rips.latest_facet(triangle);
    if edge.value != triangle.value {
        return Ok(None);
    }
    let mut found = None;
    #[cfg(test)]
    let mut candidates = 0;
    rips.visit_cofacets(
        edge,
        &mut || {
            #[cfg(test)]
            {
                candidates += 1;
            }
            budget.step()
        },
        |row| {
            count_cofacet(stats);
            if row.value == edge.value {
                if row == triangle {
                    found = Some(edge);
                }
                return Ok(false);
            }
            Ok(true)
        },
    )?;
    #[cfg(test)]
    {
        stats.apparent_candidates += candidates;
    }
    Ok(found)
}

fn allocation(context: &'static str) -> Error {
    Error::AllocationFailed { context }
}

fn push_heap<T: Ord>(heap: &mut BinaryHeap<T>, value: T) -> Result<()> {
    heap.try_reserve(1)
        .map_err(|_| allocation("Rips working heap"))?;
    heap.push(value);
    Ok(())
}

/// F2 cancellation, including repeated entries introduced by multiple additions.
fn pop_parity<T: Ord + Copy, const CONTROLLED: bool>(
    heap: &mut BinaryHeap<T>,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<Option<T>> {
    while !heap.is_empty() {
        budget.step()?;
        // Nonempty heap established above; count before removing the entry.
        let value = heap.pop().unwrap();
        let mut odd = true;
        while heap.peek() == Some(&value) {
            budget.step()?;
            heap.pop();
            odd = !odd;
        }
        if odd {
            return Ok(Some(value));
        }
    }
    Ok(None)
}

fn count_cofacet(_stats: &mut Stats) {
    #[cfg(test)]
    {
        _stats.cofacets += 1;
    }
}

fn append_coboundary<const CONTROLLED: bool>(
    rips: &impl FlagAccess,
    edge: SimplexEntry,
    heap: &mut Coboundary,
    scratch: &mut Vec<Reverse<WorkingEntry>>,
    stats: &mut Stats,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<()> {
    scratch.clear();
    #[cfg(test)]
    let mut candidates = 0;
    rips.visit_cofacets(
        edge,
        &mut || {
            #[cfg(test)]
            {
                candidates += 1;
            }
            budget.step()
        },
        |row| {
            count_cofacet(stats);
            scratch
                .try_reserve(1)
                .map_err(|_| allocation("Rips cofacet batch"))?;
            scratch.push(Reverse(row.into()));
            Ok(true)
        },
    )?;
    budget.check()?;
    heap.try_reserve(scratch.len())
        .map_err(|_| allocation("Rips working heap"))?;
    heap.extend(scratch.drain(..));
    budget.check()?;
    #[cfg(test)]
    {
        stats.reconstruction_candidates += candidates;
        stats.peak_heap = stats.peak_heap.max(heap.len());
    }
    Ok(())
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod profiling;
