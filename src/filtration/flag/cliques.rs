//! Clique enumeration from distances or a supplied graph.
use crate::complex::{Simplex, WeightedGraph};
use crate::filtration::simplicial::ZeroBornSimplicialAccess;
use crate::geometry::DissimilarityMatrixView;
use crate::{Error, Result};
pub(crate) enum CliqueAccess<'a> {
    Dense(DissimilarityMatrixView<'a>, f64),
    Sparse(&'a WeightedGraph, f64),
}

#[cfg(any(test, cocycle_h2_bench))]
impl CliqueAccess<'_> {
    /// Full triangle coboundary: dense candidates or a three-way reverse
    /// adjacency intersection. Check and count rejected candidates as well.
    pub(crate) fn visit_tetrahedra(
        &self,
        triangle: super::TupleEntry<3>,
        checkpoint: &mut impl FnMut() -> Result<()>,
        mut visitor: impl FnMut(super::TupleEntry<4>) -> Result<()>,
    ) -> Result<()> {
        let [a, b, c] = triangle.vertices;
        let mut emit = |v: usize, value: f64| {
            let mut vertices = [a, b, c, v];
            vertices.sort_unstable();
            #[cfg(test)]
            count(2, 2);
            visitor(super::TupleEntry { vertices, value })
        };
        match self {
            Self::Dense(matrix, cutoff) => {
                for v in (0..matrix.len()).rev() {
                    checkpoint()?;
                    #[cfg(test)]
                    count(2, 0);
                    if triangle.vertices.contains(&v) {
                        continue;
                    }
                    let mut value = triangle.value;
                    for u in triangle.vertices {
                        checkpoint()?;
                        #[cfg(test)]
                        count(2, 1);
                        value = value.max(matrix.get(u, v).unwrap());
                        if value > *cutoff {
                            break;
                        }
                    }
                    if value <= *cutoff {
                        emit(v, value)?;
                    }
                }
            }
            Self::Sparse(graph, cutoff) => {
                // All three vertices belong to a validated triangle. These
                // lists are sorted and omit self-neighbors, so a three-way
                // match cannot be a vertex of the triangle itself.
                let lists = [a, b, c].map(|v| graph.neighbors(v).unwrap());
                let mut positions = lists.map(<[_]>::len);
                while positions.iter().all(|&p| p > 0) {
                    checkpoint()?;
                    #[cfg(test)]
                    count(2, 0);
                    let entries = std::array::from_fn::<_, 3, _>(|i| lists[i][positions[i] - 1]);
                    let minimum = entries.iter().map(|e| e.vertex).min().unwrap();
                    let mut common = true;
                    for i in 0..3 {
                        checkpoint()?;
                        #[cfg(test)]
                        count(2, 1);
                        if entries[i].vertex > minimum {
                            positions[i] -= 1;
                            common = false;
                        }
                    }
                    if common {
                        for p in &mut positions {
                            *p -= 1;
                        }
                        let value = entries.iter().fold(triangle.value, |v, e| v.max(e.value));
                        if value <= *cutoff {
                            emit(minimum, value)?;
                        }
                    }
                }
            }
        }
        Ok(())
    }
}
impl ZeroBornSimplicialAccess for CliqueAccess<'_> {
    fn vertex_count(&self) -> usize {
        match self {
            Self::Dense(matrix, _) => matrix.len(),
            Self::Sparse(graph, _) => graph.vertex_count(),
        }
    }
    fn visit_cofacets(
        &self,
        simplex: &Simplex,
        unique: bool,
        checkpoint: &mut impl FnMut() -> Result<()>,
        mut visitor: impl FnMut(Simplex) -> Result<()>,
    ) -> Result<()> {
        let mut candidate = |vertex: usize| -> Result<()> {
            checkpoint()?;
            #[cfg(test)]
            count(simplex.dimension(), 0);
            if (unique && vertex <= *simplex.vertices.last().unwrap())
                || simplex.vertices.binary_search(&vertex).is_ok()
            {
                return Ok(());
            }
            let mut value = simplex.value;
            for &v in &simplex.vertices {
                checkpoint()?;
                #[cfg(test)]
                count(simplex.dimension(), 1);
                let (edge, cutoff) = match self {
                    Self::Dense(matrix, cutoff) => (matrix.get(v, vertex), cutoff),
                    Self::Sparse(graph, cutoff) => (graph.edge_value(v, vertex), cutoff),
                };
                let Some(edge) = edge.filter(|w| w <= cutoff) else {
                    return Ok(());
                };
                value = value.max(edge);
            }
            let mut vertices = simplex.vertices.clone();
            vertices.try_reserve(1).map_err(|_| allocation())?;
            let position = vertices.partition_point(|&v| v < vertex);
            vertices.insert(position, vertex);
            #[cfg(test)]
            count(simplex.dimension(), 2);
            visitor(Simplex { vertices, value })
        };
        match self {
            Self::Dense(matrix, _) => {
                for vertex in (0..matrix.len()).rev() {
                    candidate(vertex)?;
                }
            }
            Self::Sparse(graph, _) => {
                // Every common neighbor must occur in the shortest adjacency list.
                let neighbors = simplex
                    .vertices
                    .iter()
                    .map(|&v| graph.neighbors(v).unwrap())
                    .min_by_key(|neighbors| neighbors.len())
                    .unwrap();
                for neighbor in neighbors.iter().rev() {
                    candidate(neighbor.vertex)?;
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
thread_local! {
    static COUNTERS: std::cell::RefCell<Vec<[usize; 3]>> = const { std::cell::RefCell::new(Vec::new()) };
}
#[cfg(test)]
fn count(dimension: usize, metric: usize) {
    COUNTERS.with_borrow_mut(|rows| {
        rows.resize(rows.len().max(dimension + 1), [0; 3]);
        rows[dimension][metric] += 1;
    });
}
/// Candidates include rejection and unique-extension scans, not just cofacets.
#[cfg(test)]
pub(crate) fn cofacet_counts(reset: bool) -> Vec<[usize; 3]> {
    COUNTERS.with_borrow_mut(|rows| {
        if reset {
            std::mem::take(rows)
        } else {
            rows.clone()
        }
    })
}

fn allocation() -> Error {
    Error::AllocationFailed {
        context: "clique enumeration",
    }
}
