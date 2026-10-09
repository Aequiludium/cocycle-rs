//! Deterministic order within a simplex dimension.

use std::cmp::Ordering;

/// Within one dimension: increasing value, then decreasing combinatorial id.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SimplexEntry {
    pub(crate) id: usize,
    pub(crate) value: f64,
}

// Values come only from validated, canonicalized finite distances.
impl Eq for SimplexEntry {}
impl Ord for SimplexEntry {
    fn cmp(&self, other: &Self) -> Ordering {
        crate::complex::compare_filtration(
            self.value,
            other.value,
            Ordering::Equal,
            other.id.cmp(&self.id),
        )
    }
}
impl PartialOrd for SimplexEntry {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Fixed vertex tuples for the private H2 experiment, with the shared order.
/// No binomial index or narrowed vertex width is needed at this level.
#[cfg(any(test, cocycle_h2_bench))]
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct TupleEntry<const N: usize> {
    pub(crate) vertices: [usize; N],
    pub(crate) value: f64,
}
#[cfg(any(test, cocycle_h2_bench))]
impl<const N: usize> Eq for TupleEntry<N> {}
#[cfg(any(test, cocycle_h2_bench))]
impl<const N: usize> Ord for TupleEntry<N> {
    fn cmp(&self, other: &Self) -> Ordering {
        crate::complex::compare_filtration(
            self.value,
            other.value,
            Ordering::Equal,
            other.vertices.iter().rev().cmp(self.vertices.iter().rev()),
        )
    }
}
#[cfg(any(test, cocycle_h2_bench))]
impl<const N: usize> PartialOrd for TupleEntry<N> {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}
