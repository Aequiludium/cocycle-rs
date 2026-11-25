//! Compatible Omega bases feeding the existing private ordinary boundary engine.
use super::*;
use crate::algebra::column::Column;
use crate::complex::DirectedGraph;
use crate::diagram::{ComputationContext, PersistenceData, PersistenceDiagram};
use crate::filtration::{Coverage, FiltrationKind};
use crate::persistence::boundary::{self, BoundaryInput};
use std::collections::BTreeSet;

/// Ordinary diagram/context composed with an owned filtration-compatible path chain.
///
/// A basis prefix at an event is that event's complete Omega space. This is a
/// chain certificate, not a collection of persistent simplicial representatives.
#[derive(Clone, Debug, PartialEq)]
pub struct PathPersistence {
    data: PersistenceData,
    times: Vec<f64>,
    paths: Vec<Vec<Path>>,
    omega_bases: Vec<Vec<PathVector>>,
    births: Vec<Vec<usize>>,
    boundaries: Vec<Vec<PathVector>>,
}

impl PathPersistence {
    /// Canonical ordinary intervals, computed dimensions and complete coverage.
    pub fn diagram(&self) -> &PersistenceDiagram {
        self.data.diagram()
    }
    /// Bound field, boundary-model source, vertices and scale interpretation.
    pub fn context(&self) -> &ComputationContext {
        self.data.context()
    }
    /// Sorted exact represented f64 events, always including canonical zero.
    pub fn times(&self) -> &[f64] {
        &self.times
    }
    /// Final allowed paths in one degree, through the retained upper degree d+1.
    pub fn paths(&self, dimension: usize) -> Option<&[Path]> {
        self.paths.get(dimension).map(Vec::as_slice)
    }
    /// Compatible Omega columns in final allowed-path coordinates.
    pub fn omega_bases(&self, dimension: usize) -> Option<&[PathVector]> {
        self.omega_bases.get(dimension).map(Vec::as_slice)
    }
    /// Birth event index for each compatible Omega column.
    pub fn births(&self, dimension: usize) -> Option<&[usize]> {
        self.births.get(dimension).map(Vec::as_slice)
    }
    /// Coordinate boundary columns in the preceding compatible Omega basis.
    pub fn boundaries(&self, dimension: usize) -> Option<&[PathVector]> {
        self.boundaries.get(dimension).map(Vec::as_slice)
    }
    /// Retain common result data, discarding the path-chain certificate without cloning.
    pub fn into_data(self) -> PersistenceData {
        self.data
    }
}
impl AsRef<PersistenceData> for PathPersistence {
    fn as_ref(&self) -> &PersistenceData {
        &self.data
    }
}

/// Compute regular/full non-reduced persistence on weighted directed edges.
///
/// Edge triples are (source, target, birth). Vertices are born at zero; weights
/// are finite/nonnegative and affect appearance only. Paths are born at their
/// maximum edge birth. Exact f64 equality defines ties; signed zeros coincide.
/// No cutoff is requested: all input graph events are included, yielding complete
/// coverage in the contiguous homology domain 0..=dimension.
///
/// The positive path cap bounds every stage. `max_entries` bounds construction
/// matrices, columns/events and Omega elimination workspaces, not ordinary
/// reduction fill-in or total memory. Use the controlled entry for bounded work.
///
/// # Errors
/// Rejects invalid directed edges, births, limits, dimension overflow, allocation
/// failure or an internal filtration/chain violation, without a partial result.
#[allow(clippy::too_many_arguments)]
pub fn compute_path_persistence(
    vertex_count: usize,
    edges: &[(usize, usize, f64)],
    dimension: usize,
    field: PrimeField,
    mode: PathBoundaryMode,
    max_paths: usize,
    max_entries: usize,
) -> Result<PathPersistence> {
    compute_path_persistence_with(
        vertex_count,
        edges,
        dimension,
        field,
        mode,
        max_paths,
        max_entries,
        &Execution::default(),
    )
}

/// Share one Execution budget across input validation, every stage and final reduction.
///
/// Stages recompute exact Omega then extend earlier bases without replacing them.
/// The ordinary engine and interval assembly are the existing private upstream
/// route; no modular chain is exposed as an integer FilteredComplex.
///
/// # Errors
/// Returns the errors of [`compute_path_persistence`], cancellation or work exhaustion.
#[allow(clippy::too_many_arguments)]
pub fn compute_path_persistence_with(
    vertex_count: usize,
    edges: &[(usize, usize, f64)],
    dimension: usize,
    field: PrimeField,
    mode: PathBoundaryMode,
    max_paths: usize,
    limit: usize,
    execution: &Execution<'_>,
) -> Result<PathPersistence> {
    let mut budget = WorkBudget::new(execution)?;
    if max_paths == 0 || limit == 0 {
        return Err(Error::InvalidParameter {
            parameter: "path limits",
            reason: "must be positive",
        });
    }
    let upper = dimension.checked_add(1).ok_or(Error::SizeOverflow {
        operation: "path chain degrees",
    })?;
    let layers = upper.checked_add(1).ok_or(Error::SizeOverflow {
        operation: "path chain degrees",
    })?;
    count(layers, limit)?;
    crate::complex::path_count_cap(sum(vertex_count, edges.len())?, max_paths)?;
    count(edges.len(), limit)?;
    let mut times = Vec::new();
    push(&mut times, 0.0)?;
    let mut seen = BTreeSet::new();
    let mut pairs = Vec::new();
    for (i, &(a, b, birth)) in edges.iter().enumerate() {
        budget.step()?;
        if a >= vertex_count || b >= vertex_count || a == b || !seen.insert((a, b)) {
            return Err(Error::InvalidGraph {
                edge: i,
                reason: "invalid or duplicate directed edge",
            });
        }
        if !birth.is_finite() || birth < 0.0 {
            return Err(Error::InvalidGraph {
                edge: i,
                reason: "birth must be finite and nonnegative",
            });
        }
        push(&mut pairs, (a, b))?;
        push(&mut times, crate::canonical_zero(birth))?;
    }
    budget.check()?;
    times.sort_by(f64::total_cmp);
    budget.check()?;
    times.dedup();
    count(times.len(), limit)?;
    let final_complex = DirectedGraph::new_budget(vertex_count, pairs, &mut budget)?
        .allowed_paths_budget(upper, max_paths, &mut budget)?;
    let paths = copy_paths(&final_complex, layers, &mut budget)?;
    let mut indices = Vec::new();
    for layer in &paths {
        let mut index = BTreeMap::new();
        for (i, path) in layer.iter().enumerate() {
            budget.step_by(path.len())?;
            index.insert(path.as_slice(), i);
        }
        push(&mut indices, index)?;
    }
    let (mut bases, mut births, mut boundaries) = (Vec::new(), Vec::new(), Vec::new());
    for _ in 0..layers {
        budget.step()?;
        push(&mut bases, Vec::<PathVector>::new())?;
        push(&mut births, Vec::<usize>::new())?;
        push(&mut boundaries, Vec::<PathVector>::new())?;
    }
    let mut cells = Vec::new();
    for (stage, &time) in times.iter().enumerate() {
        budget.step()?;
        let mut pairs = Vec::new();
        for &(a, b, birth) in edges {
            budget.step()?;
            if birth <= time {
                push(&mut pairs, (a, b))?;
            }
        }
        let complex = DirectedGraph::new_budget(vertex_count, pairs, &mut budget)?
            .allowed_paths_budget(upper, max_paths, &mut budget)?;
        let local = homology(&complex, dimension, field, mode, limit, &mut budget)?;
        let mut embedding = Vec::new();
        for (q, index) in indices.iter().enumerate() {
            let mut layer = Vec::new();
            for path in complex.paths(q) {
                budget.step_by(path.len())?;
                push(
                    &mut layer,
                    *index
                        .get(path.as_slice())
                        .ok_or(invariant("stage embeds in final paths"))?,
                )?;
            }
            push(&mut embedding, layer)?;
        }
        for q in 0..layers {
            budget.step()?;
            let mut solver = Echelon::from_columns(&bases[q], field, limit, &mut budget)?;
            let lower = if q == 0 {
                None
            } else {
                Some(Echelon::from_columns(
                    &bases[q - 1],
                    field,
                    limit,
                    &mut budget,
                )?)
            };
            let images = if q == 0 {
                Vec::new()
            } else {
                apply(
                    &local.omega_bases[q - 1],
                    &local.boundaries[q],
                    field,
                    limit,
                    &mut budget,
                )?
            };
            for (i, vector) in local.omega_bases[q].iter().enumerate() {
                let mapped = embed(vector, &embedding[q], &mut budget)?;
                let j = bases[q].len();
                if solver.insert(&mapped, j, 0, &mut budget)?.is_some() {
                    continue;
                }
                count(sum(cells.len(), 1)?, limit)?;
                let column = if let Some(lower) = &lower {
                    let image = embed(&images[i], &embedding[q - 1], &mut budget)?;
                    lower.coordinates(&image, &mut budget)?
                } else {
                    PathVector::new()
                };
                push(&mut bases[q], mapped)?;
                push(&mut boundaries[q], column)?;
                matrix(&bases[q], limit)?;
                matrix(&boundaries[q], limit)?;
                push(&mut births[q], stage)?;
                push(&mut cells, (stage, q, j))?;
            }
            if bases[q].len() != local.omega_bases[q].len() {
                return Err(invariant("compatible prefix spans each stage"));
            }
        }
    }
    let mut global = Vec::new();
    for layer in &births {
        let mut indices = Vec::new();
        for _ in layer {
            budget.step()?;
            push(&mut indices, 0usize)?;
        }
        push(&mut global, indices)?;
    }
    for (i, &(_, q, j)) in cells.iter().enumerate() {
        budget.step()?;
        global[q][j] = i;
    }
    let mut input = BoundaryInput::new();
    for (i, &(stage, q, j)) in cells.iter().enumerate() {
        budget.step()?;
        let mut column = Column::new();
        for (&row, &value) in &boundaries[q][j] {
            budget.step()?;
            if q == 0 || row >= global[q - 1].len() {
                return Err(invariant("graded Omega boundary range"));
            }
            let index = global[q - 1][row];
            if index >= i || cells[index].0 > stage || cells[index].1 + 1 != q {
                return Err(invariant("Omega face precedes coface"));
            }
            column.add_term(index, value, field);
        }
        if q > 0
            && apply(
                &boundaries[q - 1],
                std::slice::from_ref(&boundaries[q][j]),
                field,
                limit,
                &mut budget,
            )?
            .iter()
            .any(|v| !v.is_empty())
        {
            return Err(invariant("compatible boundary square"));
        }
        input.push((q, j), q, times[stage], column)?;
    }
    let diagram = boundary::diagram(input, dimension, field, Coverage::Complete, &mut budget)?;
    let kind = match mode {
        PathBoundaryMode::Regular => FiltrationKind::RegularDirectedPaths,
        PathBoundaryMode::NonRegular => FiltrationKind::NonRegularDirectedPaths,
    };
    let context = ComputationContext::new(field, kind, vertex_count, None, None, None);
    drop(indices);
    budget.check()?;
    Ok(PathPersistence {
        data: PersistenceData::new(diagram, context),
        times,
        paths,
        omega_bases: bases,
        births,
        boundaries,
    })
}

fn embed(
    vector: &PathVector,
    indices: &[usize],
    budget: &mut WorkBudget<'_>,
) -> Result<PathVector> {
    let mut mapped = PathVector::new();
    for (&row, &value) in vector {
        budget.step()?;
        mapped.insert(
            *indices
                .get(row)
                .ok_or(invariant("stage path coordinate range"))?,
            value,
        );
    }
    Ok(mapped)
}
