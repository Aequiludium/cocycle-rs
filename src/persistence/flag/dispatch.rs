//! Default selection for exact flag inputs; algorithms do not call this module.
use super::clearing::{IgnoreTriangles, TriangleClearing, TrianglePivots};
use super::{cohomology, h0};
use crate::Result;
use crate::algebra::PrimeField;
use crate::complex::{WeightedEdge, WeightedGraph};
use crate::execution::WorkBudget;
use crate::filtration::flag::{
    BitsetFlag, CliqueAccess, DenseFlag, FlagAccess, SimplexEntry, SparseFlag,
};
use crate::geometry::DissimilarityMatrixView;
use crate::persistence::{RawIntervals, simplicial};

pub(in crate::persistence) fn compute_dense(
    input: DissimilarityMatrixView<'_>,
    dimension: usize,
    cutoff: f64,
    field: PrimeField,
    budget: &mut WorkBudget<'_>,
) -> Result<RawIntervals> {
    budget.check()?;
    if dimension == 0 {
        return h0::compute(
            input.len(),
            (0..input.len()).flat_map(|b| (0..b).map(move |a| (a, b, input.get(a, b).unwrap()))),
            cutoff,
            budget,
        );
    }
    let stop = if field.characteristic() == 2 && budget.is_unlimited() {
        stopping_scale(input, cutoff, &mut WorkBudget::unlimited())?
    } else {
        stopping_scale(input, cutoff, budget)?
    };
    let access = CliqueAccess::Dense(input, stop);
    if field.characteristic() != 2
        || (dimension > 1 && !super::selection::supports_h1_indices(input.len()))
    {
        return simplicial::cohomology::compute(&access, dimension, field, budget);
    }
    #[cfg(cocycle_h2_bench)]
    if dimension == 2 {
        return cohomology::compute_h2(&DenseFlag::new(input, stop)?, &access, budget);
    }
    if dimension == 1 {
        return dense_h1(input, stop, budget, &mut IgnoreTriangles);
    }
    let mut clearing = TrianglePivots::default();
    let intervals = dense_h1(input, stop, budget, &mut clearing)?;
    finish_higher(
        &access,
        input.len(),
        dimension,
        field,
        intervals,
        clearing,
        budget,
    )
}

fn stopping_scale<const CONTROLLED: bool>(
    input: DissimilarityMatrixView<'_>,
    cutoff: f64,
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<f64> {
    Ok(
        cutoff.min(crate::filtration::rips::cone_radius(input, &mut || {
            budget.step()
        })?),
    )
}

fn dense_h1(
    input: DissimilarityMatrixView<'_>,
    stop: f64,
    budget: &mut WorkBudget<'_>,
    clearing: &mut impl TriangleClearing,
) -> Result<RawIntervals> {
    // Real limits/cancellation share the original budget; only unlimited calls
    // specialize the F2 H1 hot loop to compile-time no-op checkpoints.
    if budget.is_unlimited() {
        compute_dense_h1(input, stop, &mut WorkBudget::unlimited(), clearing)
    } else {
        compute_dense_h1(input, stop, budget, clearing)
    }
}

fn compute_dense_h1<const CONTROLLED: bool>(
    input: DissimilarityMatrixView<'_>,
    stop: f64,
    budget: &mut WorkBudget<'_, CONTROLLED>,
    clearing: &mut impl TriangleClearing,
) -> Result<RawIntervals> {
    let access = DenseFlag::new(input, stop)?;
    let edges = access.edges(&mut || budget.step())?;
    budget.check()?;
    let pairs = crate::geometry::pair_count(input.len()).ok_or(crate::Error::SizeOverflow {
        operation: "matrix flag edge count",
    })?;
    if super::selection::prefer_sparse(input.len(), edges.len(), pairs) {
        let mut graph_edges = Vec::new();
        graph_edges
            .try_reserve_exact(edges.len())
            .map_err(|_| crate::Error::AllocationFailed {
                context: "matrix flag graph edges",
            })?;
        for edge in &edges {
            budget.step()?;
            graph_edges.push(WeightedEdge {
                vertices: access.edge_vertices(edge.id),
                value: edge.value,
            });
        }
        budget.check()?;
        let graph = WeightedGraph::new(input.len(), graph_edges)?;
        budget.check()?;
        return reduce_graph_edges(&graph, stop, Some(edges), budget, clearing);
    }
    if super::selection::prefer_ordered(input.len(), &edges, pairs) {
        return super::ordered::compute_with_clearing(input, stop, edges, budget, clearing);
    }
    cohomology::compute_prepared_edges(&access, edges, budget, clearing)
}

pub(in crate::persistence) fn compute_graph(
    graph: &WeightedGraph,
    dimension: usize,
    cutoff: f64,
    field: PrimeField,
    budget: &mut WorkBudget<'_>,
) -> Result<RawIntervals> {
    budget.check()?;
    if dimension == 0 {
        return h0::compute(
            graph.vertex_count(),
            graph
                .edges()
                .iter()
                .map(|e| (e.vertices[0], e.vertices[1], e.value)),
            cutoff,
            budget,
        );
    }
    let access = CliqueAccess::Sparse(graph, cutoff);
    if field.characteristic() != 2
        || (dimension > 1 && !super::selection::supports_h1_indices(graph.vertex_count()))
    {
        return simplicial::cohomology::compute(&access, dimension, field, budget);
    }
    #[cfg(cocycle_h2_bench)]
    if dimension == 2 {
        return cohomology::compute_h2(&SparseFlag::new(graph, cutoff)?, &access, budget);
    }
    if dimension == 1 {
        return graph_h1(graph, cutoff, budget, &mut IgnoreTriangles);
    }
    let mut clearing = TrianglePivots::default();
    let intervals = graph_h1(graph, cutoff, budget, &mut clearing)?;
    finish_higher(
        &access,
        graph.vertex_count(),
        dimension,
        field,
        intervals,
        clearing,
        budget,
    )
}

fn graph_h1(
    graph: &WeightedGraph,
    cutoff: f64,
    budget: &mut WorkBudget<'_>,
    clearing: &mut impl TriangleClearing,
) -> Result<RawIntervals> {
    if budget.is_unlimited() {
        reduce_graph_edges(graph, cutoff, None, &mut WorkBudget::unlimited(), clearing)
    } else {
        reduce_graph_edges(graph, cutoff, None, budget, clearing)
    }
}
/// Matrix adapters and supplied/threshold graphs share the exact cache policy.
/// The memory gate belongs to BitsetFlag; missing graph edges remain absent.
/// Prepared edges, when supplied by the matrix adapter, must describe exactly
/// this access's retained set in forward order. Otherwise prepare them once
/// through the chosen access, avoiding a second sparse index allocation.
fn reduce_graph_edges<const CONTROLLED: bool>(
    graph: &WeightedGraph,
    cutoff: f64,
    prepared_edges: Option<Vec<SimplexEntry>>,
    budget: &mut WorkBudget<'_, CONTROLLED>,
    clearing: &mut impl TriangleClearing,
) -> Result<RawIntervals> {
    budget.check()?;
    if let Some(bitset) = BitsetFlag::new_if_useful(graph, cutoff, &mut || budget.step())? {
        #[cfg(test)]
        GRAPH_ACCESS_SELECTIONS.with(|trace| trace.borrow_mut().push(true));
        let edges = match prepared_edges {
            Some(edges) => edges,
            None => bitset.edges(&mut || budget.step())?,
        };
        return cohomology::compute_prepared_edges(&bitset, edges, budget, clearing);
    }
    #[cfg(test)]
    GRAPH_ACCESS_SELECTIONS.with(|trace| trace.borrow_mut().push(false));
    let sparse = SparseFlag::new(graph, cutoff)?;
    let edges = match prepared_edges {
        Some(edges) => edges,
        None => sparse.edges(&mut || budget.step())?,
    };
    cohomology::compute_prepared_edges(&sparse, edges, budget, clearing)
}

fn finish_higher(
    access: &CliqueAccess<'_>,
    vertices: usize,
    dimension: usize,
    field: PrimeField,
    intervals: RawIntervals,
    clearing: TrianglePivots,
    budget: &mut WorkBudget<'_>,
) -> Result<RawIntervals> {
    debug_assert_eq!(field.characteristic(), 2);
    budget.check()?;
    // Every triangle has a nonzero F2 boundary, so any retained triangle gives
    // rank(delta_1) > 0. No H1 pivots therefore certify no triangles, hence no
    // higher cliques. This avoids repeating O(n) preparation on isolated graphs.
    if clearing.is_empty() {
        return Ok(intervals);
    }
    let deaths = clearing.into_triangles(vertices, budget)?;
    #[cfg(test)]
    simplicial::cohomology::workspace_event(
        "h1_released",
        1,
        &[
            (
                "intervals",
                intervals.capacity() * std::mem::size_of::<(usize, f64, Option<f64>)>(),
            ),
            (
                "handoff",
                deaths.capacity() * std::mem::size_of::<[usize; 3]>(),
            ),
        ],
        &[("deaths", deaths.len())],
    )?;
    simplicial::cohomology::continue_from_h1(access, dimension, intervals, deaths, budget)
}

#[cfg(test)]
std::thread_local! {
    static GRAPH_ACCESS_SELECTIONS: std::cell::RefCell<Vec<bool>> = const { std::cell::RefCell::new(Vec::new()) };
}

#[cfg(test)]
mod tests;
