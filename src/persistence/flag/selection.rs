//! Empirical cost policy for exact F2 H1 matrix access, not topology rules.
//! Thresholds retain the run-003/run-004 native experiments' selection policy;
//! see docs/development/persistence-reduction.md for evidence and limitations.
use crate::filtration::flag::SimplexEntry;

const MIN_ADAPTIVE_VERTICES: usize = 64;
const SPARSE_MISSING_EDGE_DIVISOR: usize = 3;
const ORDERED_MISSING_EDGE_DIVISOR: usize = 10;

/// Higher-dimensional tuple access remains available when H1's compact IDs
/// cannot represent all triangles, preserving the former high-dimensional domain.
pub(super) fn supports_h1_indices(vertices: usize) -> bool {
    crate::filtration::flag::choose(vertices, 3).is_ok()
}

pub(super) fn prefer_sparse(vertices: usize, retained: usize, pairs: usize) -> bool {
    vertices >= MIN_ADAPTIVE_VERTICES && retained <= pairs - pairs / SPARSE_MISSING_EDGE_DIVISOR
}

pub(super) fn prefer_ordered(vertices: usize, edges: &[SimplexEntry], pairs: usize) -> bool {
    if vertices < MIN_ADAPTIVE_VERTICES
        || edges.len() < pairs - pairs / ORDERED_MISSING_EDGE_DIVISOR
    {
        return false;
    }
    // Sorted edges must contain at least three distinct retained weights.
    // One/two-level filtrations favor the existing cheap heap traversal.
    let Some(first) = edges.first() else {
        return false;
    };
    let next = edges.partition_point(|edge| edge.value == first.value);
    next < edges.len() && edges[next].value != edges.last().unwrap().value
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cost_boundaries_preserve_small_sparse_and_degenerate_routes() {
        assert!(supports_h1_indices(64));
        assert!(!supports_h1_indices(usize::MAX));
        assert!(!prefer_sparse(63, 0, 1953));
        assert!(prefer_sparse(64, 1344, 2016));
        assert!(!prefer_sparse(64, 1345, 2016));
        let mut edges: Vec<_> = (0..1815).map(|id| SimplexEntry { id, value: 1. }).collect();
        assert!(!prefer_ordered(64, &edges, 2016));
        edges[1814].value = 3.;
        assert!(!prefer_ordered(64, &edges, 2016));
        edges[1813].value = 2.;
        assert!(prefer_ordered(64, &edges, 2016));
        assert!(!prefer_ordered(63, &edges, 2016));
        assert!(!prefer_ordered(64, &edges[..1814], 2016));
    }
}
