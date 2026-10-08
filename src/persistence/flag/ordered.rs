//! F2 H1 via ordered coboundary cursors and an implicit working column.
#[cfg(test)]
use super::clearing::IgnoreTriangles;
use super::clearing::TriangleClearing;
use super::edges::{ColumnPosition, EdgePosition};
use crate::execution::WorkBudget;
use crate::filtration::flag::{DenseFlag, FlagAccess, SimplexEntry, SimplexIndex};
use crate::geometry::{DissimilarityMatrixView, MatrixLayout};
use crate::persistence::RawIntervals;
use crate::{Error, Result};
use std::cmp::Reverse;
use std::collections::{BinaryHeap, HashMap};

#[derive(Clone, Copy)]
struct Neighbor {
    vertex: usize,
    value: f64,
}
struct OrderedCoboundaryAccess<'a> {
    input: DissimilarityMatrixView<'a>,
    rows: Vec<Vec<Neighbor>>,
    index: SimplexIndex,
}
impl<'a> OrderedCoboundaryAccess<'a> {
    fn new<const CONTROLLED: bool>(
        input: DissimilarityMatrixView<'a>,
        stop: f64,
        budget: &mut WorkBudget<'_, CONTROLLED>,
    ) -> Result<Self> {
        budget.check()?;
        let n = input.len();
        let index = SimplexIndex::new(n)?;
        let mut rows = Vec::new();
        rows.try_reserve_exact(n).map_err(|_| allocation())?;
        for a in 0..n {
            let mut row = Vec::new();
            row.try_reserve_exact(n).map_err(|_| allocation())?;
            for v in 0..n {
                budget.step()?;
                if a != v {
                    let value = input.get(a, v).unwrap();
                    if value <= stop {
                        row.push(Neighbor { vertex: v, value });
                    }
                }
            }
            budget.check()?;
            row.sort_unstable_by(|x, y| {
                x.value
                    .total_cmp(&y.value)
                    .then_with(|| y.vertex.cmp(&x.vertex))
            });
            budget.check()?;
            rows.push(row);
        }
        Ok(Self { input, rows, index })
    }
    // CofacetCursor endpoints and neighbor IDs are validated, distinct input vertices.
    // Reuse checked index prefixes rather than recomputing fallible pair counts.
    #[inline]
    fn distance(&self, i: usize, j: usize) -> f64 {
        let values = self.input.values();
        let n = self.input.len();
        let index = match self.input.layout() {
            MatrixLayout::LowerTriangle => self.index.edge(i, j),
            MatrixLayout::UpperTriangle => {
                let (a, b) = if i < j { (i, j) } else { (j, i) };
                self.index.edge(0, n) - self.index.edge(0, n - a) + b - a - 1
            }
            MatrixLayout::Square => i * n + j,
        };
        crate::canonical_zero(values[index])
    }
    fn triangle(&self, a: usize, b: usize, v: usize) -> usize {
        self.index.triangle(a, b, v)
    }
}
// For one edge, the initial equal-valued group is visited in decreasing colex
// order. Later cofacets arrive when the later of their two other edges appears;
// equal-weight events are merged by decreasing third-vertex ID, with one emitter.
// Thus each cursor follows the original full filtration order without sorting
// or storing its entire triangle coboundary.
#[derive(Clone)]
struct CofacetCursor {
    edge: SimplexEntry,
    a: usize,
    b: usize,
    base: usize,
    left: usize,
    right: usize,
}
impl CofacetCursor {
    fn new(
        data: &OrderedCoboundaryAccess<'_>,
        rips: &impl FlagAccess,
        edge: SimplexEntry,
        bound: Option<SimplexEntry>,
    ) -> Self {
        let [a, b] = rips.edge_vertices(edge.id);
        let lower = bound.map_or(edge.value, |x| x.value.max(edge.value));
        let start =
            |row: &[Neighbor]| row.partition_point(|x| x.value <= edge.value || x.value < lower);
        Self {
            edge,
            a,
            b,
            base: if lower > edge.value {
                0
            } else {
                data.input.len()
            },
            left: start(&data.rows[a]),
            right: start(&data.rows[b]),
        }
    }
    fn next<const CONTROLLED: bool>(
        &mut self,
        data: &OrderedCoboundaryAccess<'_>,
        bound: Option<SimplexEntry>,
        budget: &mut WorkBudget<'_, CONTROLLED>,
    ) -> Result<Option<SimplexEntry>> {
        while self.base > 0 {
            budget.step()?;
            self.base -= 1;
            let v = self.base;
            if v == self.a || v == self.b {
                continue;
            }
            if data.distance(self.a, v) <= self.edge.value
                && data.distance(self.b, v) <= self.edge.value
            {
                let row = SimplexEntry {
                    id: data.triangle(self.a, self.b, v),
                    value: self.edge.value,
                };
                if bound.is_none_or(|x| row >= x) {
                    return Ok(Some(row));
                }
            }
        }
        loop {
            budget.step()?;
            let left = data.rows[self.a].get(self.left).copied();
            let right = data.rows[self.b].get(self.right).copied();
            let from_left = match (left, right) {
                (None, None) => return Ok(None),
                (Some(_), None) => true,
                (None, Some(_)) => false,
                (Some(x), Some(y)) => {
                    x.value < y.value || (x.value == y.value && x.vertex >= y.vertex)
                }
            };
            let candidate = if from_left {
                self.left += 1;
                left.unwrap()
            } else {
                self.right += 1;
                right.unwrap()
            };
            let v = candidate.vertex;
            if v == self.a || v == self.b {
                continue;
            }
            let other = data.distance(if from_left { self.b } else { self.a }, v);
            if if from_left {
                other > candidate.value
            } else {
                other >= candidate.value
            } {
                continue;
            }
            let row = SimplexEntry {
                id: data.triangle(self.a, self.b, v),
                value: candidate.value,
            };
            if bound.is_none_or(|x| row >= x) {
                return Ok(Some(row));
            }
        }
    }
}
/// Position in the current column's cursor array, distinct from edge positions.
#[derive(Clone, Copy, PartialEq, Eq)]
struct CursorPosition(usize);
impl CursorPosition {
    // A reinserted pivot is consumed once and has no cursor to advance. Retain
    // the compact heap representation while keeping its sentinel meaning local.
    const SINGLE_ROW: Self = Self(usize::MAX);
}
#[derive(Clone, Copy)]
struct CursorHead {
    row: SimplexEntry,
    cursor: CursorPosition,
}
impl PartialEq for CursorHead {
    fn eq(&self, o: &Self) -> bool {
        self.row == o.row
    }
}
impl Eq for CursorHead {}
impl Ord for CursorHead {
    fn cmp(&self, o: &Self) -> std::cmp::Ordering {
        self.row.cmp(&o.row)
    }
}
impl PartialOrd for CursorHead {
    fn partial_cmp(&self, o: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(o))
    }
}
struct TransformColumn {
    edge: EdgePosition,
    additions: Vec<EdgePosition>,
}
fn allocation() -> Error {
    Error::AllocationFailed {
        context: "ordered Rips coboundaries",
    }
}
#[derive(Default, Clone)]
struct WorkingColumn {
    cursors: Vec<CofacetCursor>,
    heap: BinaryHeap<Reverse<CursorHead>>,
}
impl WorkingColumn {
    fn clear(&mut self) {
        self.cursors.clear();
        self.heap.clear();
    }
    fn add<const CONTROLLED: bool>(
        &mut self,
        data: &OrderedCoboundaryAccess<'_>,
        rips: &impl FlagAccess,
        edge: SimplexEntry,
        bound: Option<SimplexEntry>,
        budget: &mut WorkBudget<'_, CONTROLLED>,
    ) -> Result<()> {
        let mut cursor = CofacetCursor::new(data, rips, edge, bound);
        if let Some(row) = cursor.next(data, bound, budget)? {
            self.cursors.try_reserve(1).map_err(|_| allocation())?;
            self.heap.try_reserve(1).map_err(|_| allocation())?;
            let id = self.cursors.len();
            self.cursors.push(cursor);
            self.heap.push(Reverse(CursorHead {
                row,
                cursor: CursorPosition(id),
            }));
        }
        Ok(())
    }
    fn pop<const CONTROLLED: bool>(
        &mut self,
        data: &OrderedCoboundaryAccess<'_>,
        budget: &mut WorkBudget<'_, CONTROLLED>,
    ) -> Result<Option<SimplexEntry>> {
        while let Some(Reverse(first)) = self.heap.pop() {
            let row = first.row;
            let mut odd = false;
            let mut head = first;
            loop {
                budget.step()?;
                odd = !odd;
                if head.cursor != CursorPosition::SINGLE_ROW
                    && let Some(next) = self.cursors[head.cursor.0].next(data, None, budget)?
                {
                    self.heap.push(Reverse(CursorHead {
                        row: next,
                        cursor: head.cursor,
                    }));
                }
                if self.heap.peek().is_some_and(|x| x.0.row == row) {
                    head = self.heap.pop().unwrap().0;
                } else {
                    break;
                }
            }
            if odd {
                return Ok(Some(row));
            }
        }
        Ok(None)
    }
}

#[cfg(test)]
fn compute<const CONTROLLED: bool>(
    input: DissimilarityMatrixView<'_>,
    stop: f64,
    edges: Vec<SimplexEntry>,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<RawIntervals> {
    compute_with_clearing(input, stop, edges, budget, &mut IgnoreTriangles)
}

/// Retain every H1 death triangle, including omitted zero apparent pairs.
pub(super) fn compute_with_clearing<const CONTROLLED: bool>(
    input: DissimilarityMatrixView<'_>,
    stop: f64,
    edges: Vec<SimplexEntry>,
    budget: &mut WorkBudget<'_, CONTROLLED>,
    clearing: &mut impl TriangleClearing,
) -> Result<RawIntervals> {
    reduce(input, stop, edges, budget, clearing, &mut Stats::default())
}

/// Test diagnostics disappear entirely from ordinary library builds.
#[derive(Default)]
struct Stats {
    #[cfg(test)]
    verify_transforms: bool,
    #[cfg(test)]
    checked_transforms: usize,
    #[cfg(test)]
    stored_additions: usize,
    #[cfg(test)]
    virtual_additions: usize,
    #[cfg(test)]
    omitted_pairs: usize,
}

// An ordinary owner stores V; an omitted apparent owner has implicit V = e_j.
enum CancellationSource {
    Stored(ColumnPosition),
    Apparent(EdgePosition),
}

fn cancellation_source<const CONTROLLED: bool>(
    data: &OrderedCoboundaryAccess<'_>,
    rips: &impl FlagAccess,
    edges: &[SimplexEntry],
    pivot: SimplexEntry,
    owners: &HashMap<usize, ColumnPosition>,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<Option<CancellationSource>> {
    if let Some(&owner) = owners.get(&pivot.id) {
        return Ok(Some(CancellationSource::Stored(owner)));
    }
    let facet = rips.latest_facet(pivot);
    if facet.value != pivot.value {
        return Ok(None);
    }
    let mut apparent = CofacetCursor::new(data, rips, facet, None);
    if apparent.next(data, None, budget)? != Some(pivot) {
        return Ok(None);
    }
    let position = edges
        .binary_search(&facet)
        .map_err(|_| Error::InternalInvariant {
            reason: "ordered virtual facet missing",
        })?;
    Ok(Some(CancellationSource::Apparent(EdgePosition(position))))
}

fn reduce<const CONTROLLED: bool>(
    input: DissimilarityMatrixView<'_>,
    stop: f64,
    edges: Vec<SimplexEntry>,
    budget: &mut WorkBudget<'_, CONTROLLED>,
    clearing: &mut impl TriangleClearing,
    _stats: &mut Stats,
) -> Result<RawIntervals> {
    budget.check()?;
    let rips = DenseFlag::new(input, stop)?;
    let data = OrderedCoboundaryAccess::new(input, stop, budget)?;
    let super::edges::EdgePreparation {
        cycles,
        intervals: mut raw,
    } = super::edges::classify(&rips, &edges, budget)?;
    let mut owners: HashMap<usize, ColumnPosition> = HashMap::new();
    let mut columns: Vec<TransformColumn> = Vec::new();
    let mut working = WorkingColumn::default();
    let mut transform = BinaryHeap::new();
    for j in (0..edges.len()).rev() {
        budget.step()?;
        if !cycles[j] {
            continue;
        }
        working.clear();
        transform.clear();
        let edge = edges[j];
        let mut cursor = CofacetCursor::new(&data, &rips, edge, None);
        let first = cursor.next(&data, None, budget)?;
        if let Some(row) = first
            && row.value == edge.value
            && rips.latest_facet(row) == edge
        {
            clearing.record(row)?;
            #[cfg(test)]
            {
                _stats.omitted_pairs += 1;
            }
            continue;
        }
        if let Some(row) = first {
            working.cursors.try_reserve(1).map_err(|_| allocation())?;
            working.heap.try_reserve(1).map_err(|_| allocation())?;
            working.cursors.push(cursor);
            working.heap.push(Reverse(CursorHead {
                row,
                cursor: CursorPosition(0),
            }));
        }
        let pivot = loop {
            let Some(pivot) = working.pop(&data, budget)? else {
                break None;
            };
            let Some(owner) = cancellation_source(&data, &rips, &edges, pivot, &owners, budget)?
            else {
                break Some(pivot);
            };
            let source = match owner {
                CancellationSource::Stored(position) => columns[position.0].edge,
                CancellationSource::Apparent(position) => position,
            };
            working.heap.try_reserve(1).map_err(|_| allocation())?;
            working.heap.push(Reverse(CursorHead {
                row: pivot,
                cursor: CursorPosition::SINGLE_ROW,
            }));
            // The owner's reduced column has no nonzero rows before this pivot.
            // Skipping that prefix in every source cursor preserves their XOR.
            working.add(&data, &rips, edges[source.0], Some(pivot), budget)?;
            transform.try_reserve(1).map_err(|_| allocation())?;
            transform.push(source);
            match owner {
                CancellationSource::Stored(owner) => {
                    #[cfg(test)]
                    {
                        _stats.stored_additions += 1;
                    }
                    for &k in &columns[owner.0].additions {
                        working.add(&data, &rips, edges[k.0], Some(pivot), budget)?;
                        transform.try_reserve(1).map_err(|_| allocation())?;
                        transform.push(k);
                    }
                }
                CancellationSource::Apparent(_) => {
                    if source.0 <= j || !cycles[source.0] {
                        return Err(Error::InternalInvariant {
                            reason: "ordered apparent order",
                        });
                    }
                    #[cfg(test)]
                    {
                        _stats.virtual_additions += 1;
                    }
                }
            }
        };
        if let Some(pivot) = pivot {
            clearing.record(pivot)?;
            let mut additions = Vec::new();
            while !transform.is_empty() {
                budget.step()?;
                let k = transform.pop().unwrap();
                let mut odd = true;
                while transform.peek() == Some(&k) {
                    budget.step()?;
                    transform.pop();
                    odd = !odd;
                }
                if odd {
                    additions.try_reserve(1).map_err(|_| allocation())?;
                    additions.push(k);
                }
            }
            #[cfg(test)]
            if _stats.verify_transforms {
                tests::check_transform(&data, stop, &edges, edge, &additions, &working, pivot);
                _stats.checked_transforms += 1;
            }
            owners.try_reserve(1).map_err(|_| allocation())?;
            columns.try_reserve(1).map_err(|_| allocation())?;
            owners.insert(pivot.id, ColumnPosition(columns.len()));
            columns.push(TransformColumn {
                edge: EdgePosition(j),
                additions,
            });
            if edge.value != pivot.value {
                raw.try_reserve(1).map_err(|_| allocation())?;
                raw.push((1, edge.value, Some(pivot.value)));
            }
        } else {
            raw.try_reserve(1).map_err(|_| allocation())?;
            raw.push((1, edge.value, None));
        }
    }
    Ok(raw)
}

#[cfg(test)]
mod tests;
