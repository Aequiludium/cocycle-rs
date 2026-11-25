//! Filtered H0/H1 boundary presentations without constructing an Omega2 basis.
use super::*;
use crate::algebra::column::Column;
use crate::complex::DirectedGraph;
use crate::diagram::{ComputationContext, PersistenceData, PersistenceDiagram};
use crate::filtration::{Coverage, FiltrationKind};
use crate::persistence::boundary::{self, BoundaryInput};

/// Canonical H0/H1 results plus dependent generators of the filtered one-boundaries.
#[derive(Clone, Debug, PartialEq)]
pub struct LowDimensionalPathPersistence {
    data: PersistenceData,
    times: Vec<f64>,
    edges: Vec<(usize, usize)>,
    generators: Vec<PathVector>,
    births: Vec<usize>,
    two_path_count: usize,
}
impl LowDimensionalPathPersistence {
    /// Ordinary H0/H1 diagram with complete coverage and preserved multiplicity.
    pub fn diagram(&self) -> &PersistenceDiagram {
        self.data.diagram()
    }
    /// Selected field and regular/full directed-path model.
    pub fn context(&self) -> &ComputationContext {
        self.data.context()
    }
    /// Exact represented events, always including zero.
    pub fn times(&self) -> &[f64] {
        &self.times
    }
    /// Lexicographic directed-edge coordinates for every generator.
    pub fn edges(&self) -> &[(usize, usize)] {
        &self.edges
    }
    /// Generators of partial Omega2 in directed-edge coordinates, not an Omega basis.
    pub fn boundary_generators(&self) -> &[PathVector] {
        &self.generators
    }
    /// Birth event index of every boundary generator.
    pub fn births(&self) -> &[usize] {
        &self.births
    }
    /// Number of traversed allowed two-paths, including returning paths.
    pub fn two_path_count(&self) -> usize {
        self.two_path_count
    }
    /// Retain common data while discarding the presentation certificate.
    pub fn into_data(self) -> PersistenceData {
        self.data
    }
}
impl AsRef<PersistenceData> for LowDimensionalPathPersistence {
    fn as_ref(&self) -> &PersistenceData {
        &self.data
    }
}

/// Compute graph H0/H1 from bigon, triangle and anchored quadrangle boundaries.
///
/// Input births/models follow [`compute_path_persistence`]. The positive path cap
/// counts vertices, edges and every traversed two-path. The positive entry cap
/// bounds presentation matrices and column counts; ordinary reduction fill-in is
/// governed by Execution instead. No arboricity bound or minimum basis is claimed.
///
/// # Errors
/// Rejects invalid directed input/births, exhausted limits, failed reservations
/// or a broken presentation chain identity, without a partial result.
pub fn compute_low_dimensional_path_persistence(
    vertex_count: usize,
    edges: &[(usize, usize, f64)],
    field: PrimeField,
    mode: PathBoundaryMode,
    max_paths: usize,
    max_entries: usize,
) -> Result<LowDimensionalPathPersistence> {
    compute_low_dimensional_path_persistence_with(
        vertex_count,
        edges,
        field,
        mode,
        max_paths,
        max_entries,
        &Execution::default(),
    )
}

/// Share execution controls through input validation, presentation and ordinary reduction.
///
/// # Errors
/// Returns the errors of [`compute_low_dimensional_path_persistence`], cancellation
/// or work exhaustion. B-tree node allocations are not fallible reservations.
pub fn compute_low_dimensional_path_persistence_with(
    vertex_count: usize,
    weighted: &[(usize, usize, f64)],
    field: PrimeField,
    mode: PathBoundaryMode,
    max_paths: usize,
    limit: usize,
    execution: &Execution<'_>,
) -> Result<LowDimensionalPathPersistence> {
    let mut budget = WorkBudget::new(execution)?;
    if max_paths == 0 || limit == 0 {
        return Err(Error::InvalidParameter {
            parameter: "path limits",
            reason: "must be positive",
        });
    }
    let mut path_count = sum(vertex_count, weighted.len())?;
    count(path_count, max_paths)?;
    count(vertex_count, limit)?;
    count(weighted.len(), limit)?;
    let (mut pairs, mut weights, mut times) = (Vec::new(), BTreeMap::new(), Vec::new());
    push(&mut times, 0.0)?;
    for (i, &(a, b, w)) in weighted.iter().enumerate() {
        budget.step()?;
        if !w.is_finite() || w < 0.0 {
            return Err(Error::InvalidGraph {
                edge: i,
                reason: "birth must be finite and nonnegative",
            });
        }
        push(&mut pairs, (a, b))?;
        weights.insert((a, b), crate::canonical_zero(w));
        push(&mut times, crate::canonical_zero(w))?;
    }
    let graph = DirectedGraph::new_budget(vertex_count, pairs, &mut budget)?;
    budget.check()?;
    times.sort_by(f64::total_cmp);
    budget.check()?;
    times.dedup();
    count(times.len(), limit)?;
    let mut edges = Vec::new();
    let mut edge_births = Vec::new();
    let mut index = BTreeMap::new();
    for (i, &edge) in graph.edges().iter().enumerate() {
        budget.step()?;
        push(&mut edges, edge)?;
        let birth = times
            .binary_search_by(|t| t.total_cmp(&weights[&edge]))
            .map_err(|_| invariant("edge birth is an event"))?;
        push(&mut edge_births, birth)?;
        index.insert(edge, i);
    }
    let mut outgoing = Vec::new();
    for _ in 0..vertex_count {
        budget.step()?;
        push(&mut outgoing, Vec::new())?;
    }
    for (i, &(a, b)) in edges.iter().enumerate() {
        budget.step()?;
        push(&mut outgoing[a], (b, i))?;
    }
    let mut anchors = BTreeMap::new();
    for (ab, &(a, b)) in edges.iter().enumerate() {
        for &(c, bc) in &outgoing[b] {
            budget.step()?;
            path_count = sum(path_count, 1)?;
            count(path_count, max_paths)?;
            if a != c || mode == PathBoundaryMode::NonRegular {
                let route = (edge_births[ab].max(edge_births[bc]), b);
                let anchor = anchors.entry((a, c)).or_insert(route);
                *anchor = (*anchor).min(route);
                count(anchors.len(), limit)?;
            }
        }
    }
    let two_path_count = path_count - vertex_count - edges.len();
    let (mut generators, mut births) = (Vec::new(), Vec::new());
    let mut entries = 0;
    let minus = field.negate(1);
    for (ab, &(a, b)) in edges.iter().enumerate() {
        for &(c, bc) in &outgoing[b] {
            budget.step()?;
            let birth = edge_births[ab].max(edge_births[bc]);
            let mut candidates = Vec::new();
            if a == c {
                if mode == PathBoundaryMode::Regular && a < b {
                    push(&mut candidates, (birth, BTreeMap::from([(ab, 1), (bc, 1)])))?;
                } else if mode == PathBoundaryMode::NonRegular {
                    let &(_, middle) = anchors
                        .get(&(a, a))
                        .ok_or(invariant("return route anchor"))?;
                    if b != middle {
                        push(
                            &mut candidates,
                            (
                                birth,
                                BTreeMap::from([
                                    (ab, 1),
                                    (bc, 1),
                                    (index[&(a, middle)], minus),
                                    (index[&(middle, a)], minus),
                                ]),
                            ),
                        )?;
                    }
                }
            } else {
                let shortcut = index.get(&(a, c)).copied();
                if let Some(ac) = shortcut {
                    push(
                        &mut candidates,
                        (
                            birth.max(edge_births[ac]),
                            BTreeMap::from([(ab, 1), (bc, 1), (ac, minus)]),
                        ),
                    )?;
                }
                let &(_, middle) = anchors.get(&(a, c)).ok_or(invariant("two-path anchor"))?;
                if b != middle && shortcut.is_none_or(|ac| edge_births[ac] > birth) {
                    push(
                        &mut candidates,
                        (
                            birth,
                            BTreeMap::from([
                                (ab, 1),
                                (bc, 1),
                                (index[&(a, middle)], minus),
                                (index[&(middle, c)], minus),
                            ]),
                        ),
                    )?;
                }
            }
            for (birth, column) in candidates {
                budget.step()?;
                entries = sum(entries, column.len())?;
                count(entries, limit)?;
                count(sum(generators.len(), 1)?, limit)?;
                push(&mut generators, column)?;
                push(&mut births, birth)?;
            }
        }
    }
    // Dependent two-boundaries form a presentation with the same im(partial2).
    // Its artificial H2 is not requested. Feed only H0/H1 to the existing engine.
    count(
        sum(sum(vertex_count, edges.len())?, generators.len())?,
        limit,
    )?;
    count(sum(sum(edges.len(), edges.len())?, entries)?, limit)?;
    let mut cells = Vec::new();
    for i in 0..vertex_count {
        budget.step()?;
        push(&mut cells, (0, 0, i))?;
    }
    for (i, &t) in edge_births.iter().enumerate() {
        budget.step()?;
        push(&mut cells, (t, 1, i))?;
    }
    for (i, &t) in births.iter().enumerate() {
        budget.step()?;
        push(&mut cells, (t, 2, i))?;
    }
    count(cells.len(), limit)?;
    budget.check()?;
    cells.sort_unstable();
    budget.check()?;
    let mut edge_rows = Vec::new();
    for _ in &edges {
        push(&mut edge_rows, 0usize)?;
    }
    let mut vertex_rows = Vec::new();
    for _ in 0..vertex_count {
        push(&mut vertex_rows, 0usize)?;
    }
    for (row, &(_, q, j)) in cells.iter().enumerate() {
        budget.step()?;
        if q == 0 {
            vertex_rows[j] = row;
        } else if q == 1 {
            edge_rows[j] = row;
        }
    }
    let mut input = BoundaryInput::new();
    for (position, &(t, q, j)) in cells.iter().enumerate() {
        budget.step()?;
        let mut column = Column::new();
        if q == 1 {
            let (a, b) = edges[j];
            column.add_term(vertex_rows[a], minus, field);
            column.add_term(vertex_rows[b], 1, field);
        }
        if q == 2 {
            let mut square = Column::new();
            for (&edge, &value) in &generators[j] {
                budget.step()?;
                let row = edge_rows[edge];
                if row >= position || edge_births[edge] > t {
                    return Err(invariant("relation respects filtration"));
                }
                column.add_term(row, value, field);
                let (a, b) = edges[edge];
                square.add_term(a, field.negate(value), field);
                square.add_term(b, value, field);
            }
            if !square.is_empty() {
                return Err(invariant("one-boundary generator is a cycle"));
            }
        }
        input.push((q, j), q, times[t], column)?;
    }
    let diagram = boundary::diagram(input, 1, field, Coverage::Complete, &mut budget)?;
    let kind = match mode {
        PathBoundaryMode::Regular => FiltrationKind::RegularDirectedPaths,
        PathBoundaryMode::NonRegular => FiltrationKind::NonRegularDirectedPaths,
    };
    let context = ComputationContext::new(field, kind, vertex_count, None, None, None);
    budget.check()?;
    Ok(LowDimensionalPathPersistence {
        data: PersistenceData::new(diagram, context),
        times,
        edges,
        generators,
        births,
        two_path_count,
    })
}
