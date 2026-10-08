//! Reject empty weighted coboundaries with bounded neighbor bitsets.
use super::{FlagAccess, SimplexEntry, SparseFlag};
use crate::complex::WeightedGraph;
use crate::{Error, Result};

// Private empirical cost gates from the run-004 adaptive-cache experiment.
// The word-count bound separately guarantees O(n+m) auxiliary storage.
const MIN_CACHE_VERTICES: usize = 64;
const MAX_CACHE_VERTICES: usize = 4096;
const PROBE_EDGES: usize = 32;
const SCAN_TO_WORD_RATIO: usize = 2;
const EMPTY_PROBE_DENOMINATOR: usize = 2;

pub(crate) struct BitsetFlag<'a> {
    base: SparseFlag<'a>,
    graph: &'a WeightedGraph,
    rows: Vec<u64>,
    words: usize,
}
impl<'a> BitsetFlag<'a> {
    pub(crate) fn new_if_useful(
        graph: &'a WeightedGraph,
        cutoff: f64,
        checkpoint: &mut impl FnMut() -> Result<()>,
    ) -> Result<Option<Self>> {
        let n = graph.vertex_count();
        if !(MIN_CACHE_VERTICES..=MAX_CACHE_VERTICES).contains(&n)
            || n * n.div_ceil(64) > graph.edge_count()
        {
            return Ok(None);
        }
        let words = n.div_ceil(64);
        let mut sampled = 0;
        let mut costly_empty = 0;
        // Sampling selects an exact access strategy; it never omits topology.
        // At most 32 bounded scans avoid paying for a cache on easy inputs.
        for edge in graph
            .edges()
            .iter()
            .step_by((graph.edge_count() / PROBE_EDGES).max(1))
            .take(PROBE_EDGES)
        {
            checkpoint()?;
            sampled += 1;
            let [a, b] = edge.vertices;
            let left = graph.neighbors(a).unwrap();
            let right = graph.neighbors(b).unwrap();
            if edge.value > cutoff
                || left.is_empty()
                || right.is_empty()
                || left.len().min(right.len()) < SCAN_TO_WORD_RATIO * words
                || left.first().unwrap().vertex > right.last().unwrap().vertex
                || right.first().unwrap().vertex > left.last().unwrap().vertex
            {
                continue;
            }
            let (mut x, mut y) = (0, 0);
            let mut common = false;
            while x < left.len() && y < right.len() {
                checkpoint()?;
                match left[x].vertex.cmp(&right[y].vertex) {
                    std::cmp::Ordering::Less => x += 1,
                    std::cmp::Ordering::Greater => y += 1,
                    std::cmp::Ordering::Equal => {
                        if left[x].value <= cutoff && right[y].value <= cutoff {
                            common = true;
                            break;
                        }
                        x += 1;
                        y += 1;
                    }
                }
            }
            if !common {
                costly_empty += 1;
            }
        }
        if EMPTY_PROBE_DENOMINATOR * costly_empty < sampled {
            return Ok(None);
        }
        Self::new(graph, cutoff, checkpoint)
    }

    fn new(
        graph: &'a WeightedGraph,
        cutoff: f64,
        checkpoint: &mut impl FnMut() -> Result<()>,
    ) -> Result<Option<Self>> {
        let n = graph.vertex_count();
        if !(MIN_CACHE_VERTICES..=MAX_CACHE_VERTICES).contains(&n) {
            return Ok(None);
        }
        let words = n.div_ceil(64);
        let entries = n.checked_mul(words).ok_or(Error::SizeOverflow {
            operation: "neighbor bitsets",
        })?;
        // Keep the auxiliary word count bounded by the number of graph edges.
        // Very sparse graphs continue to use O(n+m) sorted adjacency alone.
        if entries > graph.edge_count() {
            return Ok(None);
        }
        checkpoint()?;
        let mut rows = Vec::new();
        rows.try_reserve_exact(entries)
            .map_err(|_| Error::AllocationFailed {
                context: "neighbor bitsets",
            })?;
        rows.resize(entries, 0_u64);
        for edge in graph.edges() {
            checkpoint()?;
            if edge.value <= cutoff {
                let [a, b] = edge.vertices;
                rows[a * words + b / 64] |= 1 << (b % 64);
                rows[b * words + a / 64] |= 1 << (a % 64);
            }
        }
        checkpoint()?;
        Ok(Some(Self {
            base: SparseFlag::new(graph, cutoff)?,
            graph,
            rows,
            words,
        }))
    }
}
impl FlagAccess for BitsetFlag<'_> {
    fn vertex_count(&self) -> usize {
        self.base.vertex_count()
    }
    fn edges(&self, checkpoint: &mut impl FnMut() -> Result<()>) -> Result<Vec<SimplexEntry>> {
        self.base.edges(checkpoint)
    }
    fn edge_vertices(&self, id: usize) -> [usize; 2] {
        self.base.edge_vertices(id)
    }
    fn latest_facet(&self, triangle: SimplexEntry) -> SimplexEntry {
        self.base.latest_facet(triangle)
    }
    #[inline]
    fn visit_cofacets(
        &self,
        edge: SimplexEntry,
        checkpoint: &mut impl FnMut() -> Result<()>,
        visitor: impl FnMut(SimplexEntry) -> Result<bool>,
    ) -> Result<()> {
        let [a, b] = self.base.edge_vertices(edge.id);
        let left = self.graph.neighbors(a).unwrap();
        let right = self.graph.neighbors(b).unwrap();
        if left.is_empty() || right.is_empty() {
            return Ok(());
        }
        // Preserve the O(1) sorted-range rejection and cheap low-degree scans.
        if left.first().unwrap().vertex > right.last().unwrap().vertex
            || right.first().unwrap().vertex > left.last().unwrap().vertex
            || left.len().min(right.len()) < SCAN_TO_WORD_RATIO * self.words
        {
            return self.base.visit_cofacets(edge, checkpoint, visitor);
        }
        let a = &self.rows[a * self.words..(a + 1) * self.words];
        let b = &self.rows[b * self.words..(b + 1) * self.words];
        for (&x, &y) in a.iter().zip(b) {
            checkpoint()?;
            if x & y != 0 {
                return self.base.visit_cofacets(edge, checkpoint, visitor);
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::complex::WeightedEdge;

    fn graph(n: usize, pattern: usize) -> (WeightedGraph, Vec<Option<f64>>) {
        let mut values = Vec::new();
        let mut edges = Vec::new();
        for b in 0..n {
            for a in 0..b {
                let present = match pattern {
                    0 => (a < n / 2) != (b < n / 2),
                    1 => a % 2 != b % 2,
                    2 => ((a * 37 + 11) % n < n / 2) != ((b * 37 + 11) % n < n / 2),
                    _ => (a * 17 + b * 13) % 7 != 0,
                };
                let value = if pattern < 3 {
                    1.
                } else {
                    ((a * 7 + b * 11) % 5) as f64 / 2.
                };
                values.push(present.then_some(value));
                if present {
                    edges.push(WeightedEdge {
                        vertices: [a, b],
                        value,
                    });
                }
            }
        }
        (WeightedGraph::new(n, edges).unwrap(), values)
    }

    #[test]
    fn bitset_rejections_preserve_independent_cofacets_and_early_stopping() {
        for n in [64, 65, 129] {
            for pattern in 0..4 {
                let (graph, values) = graph(n, pattern);
                let distance = |a: usize, b: usize| {
                    let (a, b) = if a < b { (a, b) } else { (b, a) };
                    values[b * (b - 1) / 2 + a]
                };
                for cutoff in [0., 0.5, 1., 2.] {
                    let access = BitsetFlag::new(&graph, cutoff, &mut || Ok(()))
                        .unwrap()
                        .unwrap();
                    for edge in access.edges(&mut || Ok(())).unwrap() {
                        let [a, b] = access.edge_vertices(edge.id);
                        let mut expected = Vec::new();
                        for v in (0..n).rev() {
                            if v == a || v == b {
                                continue;
                            }
                            if let (Some(av), Some(bv)) = (distance(a, v), distance(b, v)) {
                                let value = edge.value.max(av).max(bv);
                                if value <= cutoff {
                                    let mut vertices = [a, b, v];
                                    vertices.sort_unstable();
                                    let [a, b, c] = vertices;
                                    // Independent enumeration of the preceding triangle blocks.
                                    let id =
                                        (0..c).map(|k| k * k.saturating_sub(1) / 2).sum::<usize>()
                                            + b * (b - 1) / 2
                                            + a;
                                    expected.push(SimplexEntry { id, value });
                                }
                            }
                        }
                        let mut actual = Vec::new();
                        access
                            .visit_cofacets(edge, &mut || Ok(()), |row| {
                                actual.push(row);
                                Ok(true)
                            })
                            .unwrap();
                        assert_eq!(actual, expected, "n={n} pattern={pattern} cutoff={cutoff}");
                        let mut first = Vec::new();
                        access
                            .visit_cofacets(edge, &mut || Ok(()), |row| {
                                first.push(row);
                                Ok(false)
                            })
                            .unwrap();
                        assert_eq!(first, expected.into_iter().take(1).collect::<Vec<_>>());
                    }
                }
            }
        }
    }

    #[test]
    fn bitset_preparation_and_empty_scans_propagate_failures_and_recover() {
        let (graph, _) = graph(65, 1);
        for fail_at in [0, 1, 100, graph.edge_count() + 1] {
            let mut work = 0;
            let result = BitsetFlag::new(&graph, 1., &mut || {
                if work == fail_at {
                    return Err(Error::Cancelled);
                }
                work += 1;
                Ok(())
            });
            assert!(matches!(result, Err(Error::Cancelled)));
        }
        let access = BitsetFlag::new(&graph, 1., &mut || Ok(()))
            .unwrap()
            .unwrap();
        let edge = access
            .edges(&mut || Ok(()))
            .unwrap()
            .into_iter()
            .find(|edge| access.edge_vertices(edge.id) == [0, 1])
            .unwrap();
        assert_eq!(
            access.visit_cofacets(edge, &mut || Err(Error::Cancelled), |_| Ok(true)),
            Err(Error::Cancelled)
        );
        let mut visited = 0;
        access
            .visit_cofacets(edge, &mut || Ok(()), |_| {
                visited += 1;
                Ok(true)
            })
            .unwrap();
        assert_eq!(visited, 0);
        let sparse = WeightedGraph::new(
            4096,
            vec![WeightedEdge {
                vertices: [0, 4095],
                value: 1.,
            }],
        )
        .unwrap();
        assert!(
            BitsetFlag::new(&sparse, 1., &mut || panic!(
                "sparse fallback must not prepare bitsets"
            ))
            .unwrap()
            .is_none()
        );
    }
    #[test]
    fn adaptive_cache_avoids_easy_ranges_and_retains_costly_empty_scans() {
        for n in [64, 65, 129] {
            for pattern in 0..4 {
                let (graph, _) = graph(n, pattern);
                let cache = BitsetFlag::new_if_useful(&graph, 2., &mut || Ok(())).unwrap();
                assert_eq!(cache.is_some(), pattern == 1 || pattern == 2);
            }
        }
        let (graph, _) = graph(65, 1);
        assert!(
            BitsetFlag::new_if_useful(&graph, 0.5, &mut || Ok(()))
                .unwrap()
                .is_none()
        );
        assert!(matches!(
            BitsetFlag::new_if_useful(&graph, 1., &mut || Err(Error::Cancelled)),
            Err(Error::Cancelled)
        ));
    }
}
