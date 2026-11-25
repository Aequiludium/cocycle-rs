//! Directed regular allowed paths and explicit endpoint-closed complexes.
//!
//! This storage is independent of undirected graphs and integer filtered-cell
//! complexes. A boundary may contain forbidden faces; no Omega or persistence
//! result is claimed by this module.
use crate::execution::{Execution, WorkBudget};
use crate::{Error, Result};
use std::collections::{BTreeMap, BTreeSet};

/// An elementary vertex word; dimension is length minus one.
pub type Path = Vec<usize>;

/// Ambient convention for a path boundary.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PathBoundaryMode {
    /// Quotient out words with equal adjacent vertices.
    Regular,
    /// Keep the full ambient word space, including nonregular faces.
    NonRegular,
}

/// A validated finite directed graph retaining reverse edges and isolates.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct DirectedGraph {
    vertex_count: usize,
    edges: Vec<(usize, usize)>,
}

impl DirectedGraph {
    /// Validate ordered edges and sort them lexicographically without reversal.
    ///
    /// # Errors
    /// Rejects self-loops, duplicate ordered edges and endpoints outside 0..n.
    pub fn new(vertex_count: usize, edges: Vec<(usize, usize)>) -> Result<Self> {
        Self::new_with(vertex_count, edges, &Execution::default())
    }

    /// Construct with cooperative work/cancellation controls around edge visits and sorting.
    ///
    /// # Errors
    /// Returns input, work-limit or cancellation errors without partial storage.
    pub fn new_with(
        vertex_count: usize,
        edges: Vec<(usize, usize)>,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        let mut budget = WorkBudget::new(execution)?;
        Self::new_budget(vertex_count, edges, &mut budget)
    }

    pub(crate) fn new_budget(
        vertex_count: usize,
        mut edges: Vec<(usize, usize)>,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Self> {
        for (edge, &(a, b)) in edges.iter().enumerate() {
            budget.step()?;
            if a >= vertex_count || b >= vertex_count || a == b {
                return Err(Error::InvalidGraph {
                    edge,
                    reason: "invalid directed endpoint or self-loop",
                });
            }
        }
        budget.check()?;
        edges.sort_unstable();
        budget.check()?;
        for (edge, pair) in edges.windows(2).enumerate() {
            budget.step()?;
            if pair[0] == pair[1] {
                return Err(Error::InvalidGraph {
                    edge: edge + 1,
                    reason: "duplicate ordered edge",
                });
            }
        }
        Ok(Self {
            vertex_count,
            edges,
        })
    }

    /// Number of declared vertices, including isolates.
    pub fn vertex_count(&self) -> usize {
        self.vertex_count
    }
    /// Ordered edges in lexicographic order; reverse edges remain distinct.
    pub fn edges(&self) -> &[(usize, usize)] {
        &self.edges
    }

    /// Enumerate regular allowed paths to a finite dimension with a positive total count cap.
    ///
    /// Nonadjacent repetitions such as 010 are retained. Singleton paths count
    /// toward the cap. Exhaustion returns an error; the cap is not a byte limit.
    ///
    /// # Errors
    /// Rejects zero or exhausted path caps and failed fallible reservations.
    pub fn allowed_paths(&self, max_dimension: usize, max_paths: usize) -> Result<PathComplex> {
        self.allowed_paths_with(max_dimension, max_paths, &Execution::default())
    }

    /// Enumerate with one shared work budget for vertex/path and edge-candidate visits.
    ///
    /// # Errors
    /// Returns count, allocation, work-limit or cancellation errors without partial output.
    pub fn allowed_paths_with(
        &self,
        max_dimension: usize,
        max_paths: usize,
        execution: &Execution<'_>,
    ) -> Result<PathComplex> {
        count_cap(self.vertex_count, max_paths)?;
        let mut budget = WorkBudget::new(execution)?;
        self.allowed_paths_budget(max_dimension, max_paths, &mut budget)
    }

    pub(crate) fn allowed_paths_budget(
        &self,
        max_dimension: usize,
        max_paths: usize,
        budget: &mut WorkBudget<'_>,
    ) -> Result<PathComplex> {
        count_cap(self.vertex_count, max_paths)?;
        let mut layers = Vec::new();
        let mut total = self.vertex_count;
        if total > 0 {
            let mut first = Vec::new();
            reserve(&mut first, total)?;
            for v in 0..total {
                budget.step()?;
                let mut singleton = Vec::new();
                reserve(&mut singleton, 1)?;
                singleton.push(v);
                first.push(singleton);
            }
            reserve(&mut layers, 1)?;
            layers.push(first);
            for _ in 0..max_dimension {
                let mut next = Vec::new();
                for path in layers
                    .last()
                    .expect("a nonempty initial layer was established")
                {
                    budget.step()?;
                    let last = path[path.len() - 1];
                    let start = self.edges.partition_point(|&(a, _)| a < last);
                    for &(a, v) in &self.edges[start..] {
                        if a != last {
                            break;
                        }
                        budget.step()?;
                        count_cap(
                            total.checked_add(1).ok_or(Error::SizeOverflow {
                                operation: "path count",
                            })?,
                            max_paths,
                        )?;
                        let mut extended = Vec::new();
                        reserve(
                            &mut extended,
                            path.len().checked_add(1).ok_or(Error::SizeOverflow {
                                operation: "path length",
                            })?,
                        )?;
                        extended.extend_from_slice(path);
                        extended.push(v);
                        reserve(&mut next, 1)?;
                        next.push(extended);
                        total += 1;
                    }
                }
                if next.is_empty() {
                    break;
                }
                reserve(&mut layers, 1)?;
                layers.push(next);
            }
        }
        let mut truncated_at = None;
        if let Some(last_layer) = layers.last() {
            for path in last_layer {
                budget.step()?;
                let last = path[path.len() - 1];
                let start = self.edges.partition_point(|&(a, _)| a < last);
                if self.edges.get(start).is_some_and(|&(a, _)| a == last) {
                    truncated_at = Some(max_dimension);
                    break;
                }
            }
        }
        Ok(PathComplex {
            vertex_count: self.vertex_count,
            layers,
            truncated_at,
        })
    }
}

/// Complete finite regular allowed paths, or a declared graph truncation.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PathComplex {
    vertex_count: usize,
    layers: Vec<Vec<Path>>,
    truncated_at: Option<usize>,
}

impl PathComplex {
    /// Validate complete explicit paths, including every singleton and both endpoint truncations.
    ///
    /// Internal faces may be absent; high paths are not inferred from edges.
    /// This input declares a complete finite complex, not a graph truncation.
    ///
    /// # Errors
    /// Rejects empty/nonregular/duplicate paths, invalid vertices, absent singleton
    /// or endpoint truncations, zero/exhausted count caps or failed reservations.
    pub fn new(vertex_count: usize, paths: &[Path], max_paths: usize) -> Result<Self> {
        Self::new_with(vertex_count, paths, max_paths, &Execution::default())
    }

    /// Validate and copy explicit paths with cooperative controls over word visits.
    ///
    /// # Errors
    /// Returns the errors of [`Self::new`], cancellation or work exhaustion.
    pub fn new_with(
        vertex_count: usize,
        paths: &[Path],
        max_paths: usize,
        execution: &Execution<'_>,
    ) -> Result<Self> {
        count_cap(vertex_count.max(paths.len()), max_paths)?;
        let mut budget = WorkBudget::new(execution)?;
        let mut seen = BTreeSet::new();
        for path in paths {
            budget.step()?;
            budget.step_by(path.len())?;
            regular(path)?;
            if path.iter().any(|&v| v >= vertex_count) || !seen.insert(path.as_slice()) {
                return Err(invalid("invalid vertex or duplicate path"));
            }
        }
        for v in 0..vertex_count {
            budget.step()?;
            if !seen.contains([v].as_slice()) {
                return Err(invalid("missing singleton"));
            }
        }
        for path in &seen {
            budget.step_by(path.len())?;
            if path.len() > 1
                && (!seen.contains(&path[1..]) || !seen.contains(&path[..path.len() - 1]))
            {
                return Err(invalid("missing endpoint truncation"));
            }
        }
        let mut layers = Vec::new();
        for path in seen {
            budget.step_by(path.len())?;
            while layers.len() < path.len() {
                reserve(&mut layers, 1)?;
                layers.push(Vec::new());
            }
            let layer = &mut layers[path.len() - 1];
            reserve(layer, 1)?;
            let mut copied = Vec::new();
            reserve(&mut copied, path.len())?;
            copied.extend_from_slice(path);
            layer.push(copied);
        }
        Ok(Self {
            vertex_count,
            layers,
            truncated_at: None,
        })
    }
    /// Declared vertices, including isolated singletons.
    pub fn vertex_count(&self) -> usize {
        self.vertex_count
    }
    /// Paths in one dimension, or an empty slice for an absent layer.
    pub fn paths(&self, dimension: usize) -> &[Path] {
        self.layers.get(dimension).map_or(&[], Vec::as_slice)
    }
    /// Total retained elementary paths across all layers.
    pub fn path_count(&self) -> usize {
        self.layers.iter().map(Vec::len).sum()
    }
    /// Largest retained dimension, or None for the empty source.
    pub fn max_dimension(&self) -> Option<usize> {
        self.layers.len().checked_sub(1)
    }
    /// Last graph cutoff when unseen allowed extensions exist; None for a complete finite source.
    pub fn truncated_at(&self) -> Option<usize> {
        self.truncated_at
    }
}

/// Non-reduced integer ambient boundary, with equal faces combined and zeros omitted.
///
/// Singleton boundary is zero. Regular mode requires a regular input and drops
/// nonregular faces. NonRegular mode accepts any nonempty word for iterated
/// ambient calculations. Rows need not be allowed paths of a PathComplex.
///
/// # Errors
/// Rejects empty words, nonregular inputs in Regular mode or coefficient overflow.
pub fn path_boundary(path: &[usize], mode: PathBoundaryMode) -> Result<Vec<(Path, i64)>> {
    path_boundary_with(path, mode, &Execution::default())
}

/// Compute an ambient boundary with cooperative controls over face/coordinate visits.
///
/// # Errors
/// Returns the errors of [`path_boundary`], failed reservations, cancellation or
/// work exhaustion. Collection nodes may still abort on allocation failure.
pub fn path_boundary_with(
    path: &[usize],
    mode: PathBoundaryMode,
    execution: &Execution<'_>,
) -> Result<Vec<(Path, i64)>> {
    let mut budget = WorkBudget::new(execution)?;
    boundary_budget(path, mode, &mut budget)
}

pub(crate) fn boundary_budget(
    path: &[usize],
    mode: PathBoundaryMode,
    budget: &mut WorkBudget<'_>,
) -> Result<Vec<(Path, i64)>> {
    budget.step_by(path.len())?;
    if mode == PathBoundaryMode::Regular {
        regular(path)?;
    } else if path.is_empty() {
        return Err(invalid("empty path"));
    }
    let mut boundary = BTreeMap::<Path, i64>::new();
    if path.len() > 1 {
        for j in 0..path.len() {
            budget.step()?;
            if mode == PathBoundaryMode::Regular
                && j > 0
                && j + 1 < path.len()
                && path[j - 1] == path[j + 1]
            {
                continue;
            }
            budget.step_by(path.len() - 1)?;
            let mut face = Vec::new();
            reserve(&mut face, path.len() - 1)?;
            face.extend_from_slice(&path[..j]);
            face.extend_from_slice(&path[j + 1..]);
            let value = boundary.entry(face).or_insert(0);
            *value =
                value
                    .checked_add(if j % 2 == 0 { 1 } else { -1 })
                    .ok_or(Error::SizeOverflow {
                        operation: "path boundary coefficient",
                    })?;
        }
    }
    let mut result = Vec::new();
    reserve(&mut result, boundary.len())?;
    for term in boundary {
        budget.step()?;
        if term.1 != 0 {
            result.push(term);
        }
    }
    Ok(result)
}

fn reserve<T>(v: &mut Vec<T>, additional: usize) -> Result<()> {
    v.try_reserve(additional)
        .map_err(|_| Error::AllocationFailed {
            context: "directed paths",
        })
}
pub(crate) fn count_cap(count: usize, limit: usize) -> Result<()> {
    if limit == 0 || count > limit {
        Err(Error::InvalidParameter {
            parameter: "max_paths",
            reason: "must be positive and cover the complete retained path count",
        })
    } else {
        Ok(())
    }
}
fn invalid(reason: &'static str) -> Error {
    Error::InvalidComplex { cell: None, reason }
}
fn regular(path: &[usize]) -> Result<()> {
    if path.is_empty() || path.windows(2).any(|v| v[0] == v[1]) {
        Err(invalid("expected a nonempty regular path"))
    } else {
        Ok(())
    }
}
