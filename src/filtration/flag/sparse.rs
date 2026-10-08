//! Sparse cofacets by reverse sorted neighbor intersection.
use super::index::SimplexIndex;
use super::{FlagAccess, SimplexEntry};
use crate::complex::WeightedGraph;
use crate::{Error, Result};

pub(crate) struct SparseFlag<'a> {
    graph: &'a WeightedGraph,
    index: SimplexIndex,
    cutoff: f64,
}

impl<'a> SparseFlag<'a> {
    pub(crate) fn new(graph: &'a WeightedGraph, cutoff: f64) -> Result<Self> {
        Ok(Self {
            graph,
            index: SimplexIndex::new(graph.vertex_count())?,
            cutoff,
        })
    }
    fn edge(&self, a: usize, b: usize) -> SimplexEntry {
        // Called only for facets of triangles produced by neighbor intersection.
        SimplexEntry {
            id: self.index.edge(a, b),
            value: self.graph.edge_value(a, b).unwrap(),
        }
    }
}
impl FlagAccess for SparseFlag<'_> {
    fn vertex_count(&self) -> usize {
        self.graph.vertex_count()
    }
    fn edges(&self, checkpoint: &mut impl FnMut() -> Result<()>) -> Result<Vec<SimplexEntry>> {
        let mut edges = Vec::new();
        for edge in self.graph.edges() {
            checkpoint()?;
            if edge.value <= self.cutoff {
                edges.try_reserve(1).map_err(|_| Error::AllocationFailed {
                    context: "sparse flag edges",
                })?;
                edges.push(SimplexEntry {
                    id: self.index.edge(edge.vertices[0], edge.vertices[1]),
                    value: edge.value,
                });
            }
        }
        edges.sort_unstable();
        Ok(edges)
    }
    fn edge_vertices(&self, id: usize) -> [usize; 2] {
        self.index.edge_vertices(id)
    }
    #[cfg(any(test, cocycle_h2_bench))]
    fn triangle_vertices(&self, id: usize) -> [usize; 3] {
        self.index.triangle_vertices(id)
    }
    fn latest_facet(&self, triangle: SimplexEntry) -> SimplexEntry {
        let [a, b, c] = self.index.triangle_vertices(triangle.id);
        self.edge(a, b).max(self.edge(a, c)).max(self.edge(b, c))
    }
    #[inline]
    fn visit_cofacets(
        &self,
        edge: SimplexEntry,
        checkpoint: &mut impl FnMut() -> Result<()>,
        mut visitor: impl FnMut(SimplexEntry) -> Result<bool>,
    ) -> Result<()> {
        let [a, b] = self.edge_vertices(edge.id);
        let left = self.graph.neighbors(a).unwrap();
        let right = self.graph.neighbors(b).unwrap();
        let (Some(first_left), Some(first_right)) = (left.first(), right.first()) else {
            return Ok(());
        };
        // Sorted vertex ranges that do not overlap have no common neighbor.
        // Reject the whole intersection instead of walking either adjacency.
        if first_left.vertex > right.last().unwrap().vertex
            || first_right.vertex > left.last().unwrap().vertex
        {
            checkpoint()?;
            return Ok(());
        }
        let mut left = left.iter().rev();
        let mut right = right.iter().rev();
        let (mut x, mut y) = (left.next(), right.next());
        while let (Some(a_neighbor), Some(b_neighbor)) = (x, y) {
            checkpoint()?;
            match a_neighbor.vertex.cmp(&b_neighbor.vertex) {
                std::cmp::Ordering::Greater => x = left.next(),
                std::cmp::Ordering::Less => y = right.next(),
                std::cmp::Ordering::Equal => {
                    let value = edge.value.max(a_neighbor.value).max(b_neighbor.value);
                    if value <= self.cutoff
                        && !visitor(SimplexEntry {
                            id: self.index.triangle(a, b, a_neighbor.vertex),
                            value,
                        })?
                    {
                        break;
                    }
                    x = left.next();
                    y = right.next();
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::complex::WeightedEdge;

    #[test]
    fn sparse_cofacets_match_independent_triples_with_disjoint_and_overlapping_ranges() {
        for n in 0..=12 {
            for pattern in 0..4 {
                let mut edges = Vec::new();
                for b in 0..n {
                    for a in 0..b {
                        let retained = match pattern {
                            0 => a < n / 2 && b >= n / 2,
                            1 => a % 2 != b % 2,
                            2 => (a + b) % 3 != 0,
                            _ => true,
                        };
                        if retained {
                            edges.push(WeightedEdge {
                                vertices: [a, b],
                                value: ((a * 7 + b * 3) % 5) as f64 / 2.,
                            });
                        }
                    }
                }
                let graph = WeightedGraph::new(n, edges).unwrap();
                for cutoff in [0., 0.5, 1., 2.] {
                    let access = SparseFlag::new(&graph, cutoff).unwrap();
                    let mut triangles = Vec::new();
                    let mut id = 0;
                    for c in 2..n {
                        for b in 1..c {
                            for a in 0..b {
                                if let (Some(ab), Some(ac), Some(bc)) = (
                                    graph.edge_value(a, b),
                                    graph.edge_value(a, c),
                                    graph.edge_value(b, c),
                                ) {
                                    let value = ab.max(ac).max(bc);
                                    if value <= cutoff {
                                        triangles.push(([a, b, c], SimplexEntry { id, value }));
                                    }
                                }
                                id += 1;
                            }
                        }
                    }
                    for edge in access.edges(&mut || Ok(())).unwrap() {
                        let [a, b] = access.edge_vertices(edge.id);
                        let expected: Vec<_> = triangles
                            .iter()
                            .rev()
                            .filter(|(vertices, _)| vertices.contains(&a) && vertices.contains(&b))
                            .map(|(_, row)| *row)
                            .collect();
                        let mut actual = Vec::new();
                        access
                            .visit_cofacets(edge, &mut || Ok(()), |row| {
                                actual.push(row);
                                Ok(true)
                            })
                            .unwrap();
                        assert_eq!(actual, expected, "n={n} pattern={pattern} cutoff={cutoff}");
                    }
                }
            }
        }
    }
}
