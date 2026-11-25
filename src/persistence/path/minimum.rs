//! Bounded exhaustive F2 H1 optimization with exact integer support costs.
use super::*;
use crate::complex::DirectedGraph;
use std::collections::BTreeSet;

/// A globally minimum F2 H1 basis, returned only after exhaustive candidate search.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MinimumPathBasis {
    mode: PathBoundaryMode,
    edges: Vec<(usize, usize)>,
    representatives: Vec<SparseVector>,
    costs: Vec<u128>,
    total_cost: u128,
}
impl MinimumPathBasis {
    /// Ambient boundary convention; coefficients are always F2.
    pub fn mode(&self) -> PathBoundaryMode {
        self.mode
    }
    /// Lexicographic directed edge coordinates, retaining both reverse edges.
    pub fn edges(&self) -> &[(usize, usize)] {
        &self.edges
    }
    /// Selected basis cycles in these directed edge coordinates.
    pub fn representatives(&self) -> impl ExactSizeIterator<Item = PathVector<'_>> + Clone {
        self.representatives.iter().map(PathVector)
    }
    /// Exact support costs, one per selected cycle; zero costs are permitted.
    pub fn costs(&self) -> &[u128] {
        &self.costs
    }
    /// Exact total cost, accumulated without floating-point rounding.
    pub fn total_cost(&self) -> u128 {
        self.total_cost
    }
}

/// Compute a global minimum F2 H1 basis with nonnegative u64 directed-edge costs.
///
/// Enumerates every elementary cycle in the underlying undirected multigraph,
/// including two-edge cycles from reverse edges. Nonnegative costs make these
/// candidates sufficient for greedy independent homology-class selection. This
/// is exponential: candidate/search caps must be positive. `max_entries` bounds
/// candidate support and frontier entries separately, plus exact algebra caps.
///
/// # Errors
/// Returns no partial optimum on exhausted caps, failed reservations, cost/size
/// overflow, wrong cost count or a failed final cycle/basis certificate.
pub fn compute_minimum_path_basis(
    graph: &DirectedGraph,
    costs: &[u64],
    mode: PathBoundaryMode,
    max_paths: usize,
    max_entries: usize,
    max_candidates: usize,
    max_search_steps: usize,
) -> Result<MinimumPathBasis> {
    compute_minimum_path_basis_with(
        graph,
        costs,
        mode,
        max_paths,
        max_entries,
        max_candidates,
        max_search_steps,
        &Execution::default(),
    )
}

/// Share one work budget across path enumeration, projection and exhaustive search.
///
/// # Errors
/// Returns the errors of [`compute_minimum_path_basis`], cancellation or work exhaustion.
#[allow(clippy::too_many_arguments)]
pub fn compute_minimum_path_basis_with(
    graph: &DirectedGraph,
    costs: &[u64],
    mode: PathBoundaryMode,
    max_paths: usize,
    limit: usize,
    max_candidates: usize,
    max_search_steps: usize,
    execution: &Execution<'_>,
) -> Result<MinimumPathBasis> {
    let mut budget = WorkBudget::new(execution)?;
    if costs.len() != graph.edges().len() || max_candidates == 0 || max_search_steps == 0 {
        return Err(Error::InvalidParameter {
            parameter: "minimum path basis",
            reason: "one cost per edge and positive search/candidate caps are required",
        });
    }
    let source = graph.allowed_paths_budget(2, max_paths, &mut budget)?;
    let field = PrimeField::default();
    let h = homology(&source, 1, field, mode, limit, &mut budget)?;
    let projection = h.projection_budget(limit, &mut budget)?;
    let beta = h.representatives[1].len();
    let mut edges = Vec::new();
    for &e in graph.edges() {
        budget.step()?;
        push(&mut edges, e)?;
    }
    let mut adjacency = Vec::new();
    for _ in 0..graph.vertex_count() {
        budget.step()?;
        push(&mut adjacency, Vec::new())?;
    }
    for (i, &(a, b)) in edges.iter().enumerate() {
        budget.step()?;
        push(&mut adjacency[a], (b, i))?;
        push(&mut adjacency[b], (a, i))?;
    }
    let (mut candidates, mut candidate_entries, mut search_steps) = (BTreeSet::new(), 0, 0);
    if beta > 0 {
        for root in 0..graph.vertex_count() {
            budget.step()?;
            let mut vertices = Vec::new();
            push(&mut vertices, root)?;
            let mut stack = Vec::new();
            push(&mut stack, (root, vertices, Vec::<usize>::new()))?;
            let mut frontier = 1;
            while let Some((vertex, vertices, route)) = stack.pop() {
                budget.step()?;
                frontier -= vertices.len() + route.len();
                search_steps = sum(search_steps, 1)?;
                cap(search_steps, max_search_steps, "max_search_steps")?;
                for &(next, edge) in &adjacency[vertex] {
                    budget.step_by(sum(vertices.len(), route.len())?)?;
                    if next < root || route.contains(&edge) {
                        continue;
                    }
                    if next == root && !route.is_empty() {
                        let mut cycle = copy_indices(&route, &mut budget)?;
                        push(&mut cycle, edge)?;
                        budget.check()?;
                        cycle.sort_unstable();
                        budget.check()?;
                        if !candidates.contains(&cycle) {
                            cap(sum(candidates.len(), 1)?, max_candidates, "max_candidates")?;
                            candidate_entries = sum(candidate_entries, cycle.len())?;
                            count(candidate_entries, limit)?;
                            candidates.insert(cycle);
                        }
                    } else if !vertices.contains(&next) {
                        frontier = sum(sum(sum(frontier, vertices.len())?, route.len())?, 2)?;
                        count(frontier, limit)?;
                        let mut v = copy_indices(&vertices, &mut budget)?;
                        push(&mut v, next)?;
                        let mut r = copy_indices(&route, &mut budget)?;
                        push(&mut r, edge)?;
                        push(&mut stack, (next, v, r))?;
                    }
                }
            }
        }
    }
    let mut ordered = Vec::new();
    for cycle in candidates {
        let mut cost = 0u128;
        for &i in &cycle {
            budget.step()?;
            cost = cost
                .checked_add(u128::from(costs[i]))
                .ok_or(Error::SizeOverflow {
                    operation: "path cycle cost",
                })?;
        }
        push(&mut ordered, (cost, cycle))?;
    }
    budget.check()?;
    ordered.sort();
    budget.check()?;
    let mut solver = Echelon::new(field, limit);
    let omega = Echelon::from_columns(&h.omega_bases[1], field, limit, &mut budget)?;
    let (mut reps, mut selected_costs, mut total_cost) = (Vec::new(), Vec::new(), 0u128);
    for (cost, cycle) in ordered {
        budget.step()?;
        let mut vector = SparseVector::new();
        for i in cycle {
            budget.step()?;
            vector.insert(i, 1);
        }
        let coordinates = omega.coordinates(&vector, &mut budget)?;
        let mut annotations = apply(
            &projection.annotations[1],
            std::slice::from_ref(&coordinates),
            field,
            limit,
            &mut budget,
        )?;
        let annotation = annotations.pop().ok_or(invariant("one cycle annotation"))?;
        if solver
            .insert(&annotation, reps.len(), 0, &mut budget)?
            .is_none()
        {
            total_cost = total_cost.checked_add(cost).ok_or(Error::SizeOverflow {
                operation: "minimum path basis cost",
            })?;
            push(&mut reps, vector)?;
            matrix(&reps, limit)?;
            push(&mut selected_costs, cost)?;
        }
        if reps.len() == beta {
            break;
        }
    }
    let mut coordinates = Vec::new();
    for v in &reps {
        push(&mut coordinates, omega.coordinates(v, &mut budget)?)?;
        matrix(&coordinates, limit)?;
    }
    if reps.len() != beta
        || apply(&h.boundaries[1], &coordinates, field, limit, &mut budget)?
            .iter()
            .any(|v| !v.is_empty())
    {
        return Err(invariant(
            "minimum path basis is complete and consists of cycles",
        ));
    }
    budget.check()?;
    Ok(MinimumPathBasis {
        mode,
        edges,
        representatives: reps,
        costs: selected_costs,
        total_cost,
    })
}
fn copy_indices(input: &[usize], budget: &mut WorkBudget<'_>) -> Result<Vec<usize>> {
    budget.step_by(input.len())?;
    let mut out = Vec::new();
    out.try_reserve_exact(input.len())
        .map_err(|_| Error::AllocationFailed {
            context: "path cycle search",
        })?;
    out.extend_from_slice(input);
    Ok(out)
}
fn cap(count: usize, limit: usize, parameter: &'static str) -> Result<()> {
    if count > limit {
        Err(Error::InvalidParameter {
            parameter,
            reason: "exhausted before a complete optimum was found",
        })
    } else {
        Ok(())
    }
}
