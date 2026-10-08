//! Forward edge positions and H0 classification shared by the F2 H1 engines.
//! Input entries obey `FlagAccess` ordering and describe its retained edge set.
//! Clearing marks H0 death edges; it never removes edges from clique topology.
use crate::execution::WorkBudget;
use crate::filtration::flag::{FlagAccess, SimplexEntry};
use crate::persistence::{RawIntervals, union_find::UnionFind};
use crate::{Error, Result};

/// Position in the forward edge array, distinct from a simplex ID.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub(super) struct EdgePosition(pub(super) usize);

/// Position in an engine's stored transformation-column array.
#[derive(Clone, Copy)]
pub(super) struct ColumnPosition(pub(super) usize);

pub(super) struct EdgePreparation {
    pub(super) cycles: Vec<bool>,
    pub(super) intervals: RawIntervals,
}

#[inline]
pub(super) fn classify<const CONTROLLED: bool>(
    access: &impl FlagAccess,
    edges: &[SimplexEntry],
    budget: &mut WorkBudget<'_, CONTROLLED>,
) -> Result<EdgePreparation> {
    let mut forest = UnionFind::new(access.vertex_count())?;
    let mut cycles = Vec::new();
    cycles
        .try_reserve_exact(edges.len())
        .map_err(|_| Error::AllocationFailed {
            context: "flag cycle edges",
        })?;
    // H0 contributes n intervals before zero bars are removed. Grow for actual
    // H1 output later instead of reserving space for zero-lifetime pairs.
    let mut intervals = Vec::new();
    intervals
        .try_reserve(access.vertex_count())
        .map_err(|_| Error::AllocationFailed {
            context: "flag intervals",
        })?;
    for edge in edges {
        budget.step()?;
        let [a, b] = access.edge_vertices(edge.id);
        let merged = forest.merge(a, b);
        cycles.push(!merged);
        if merged {
            intervals.push((0, 0.0, Some(edge.value)));
        }
    }
    intervals.extend((0..forest.components()).map(|_| (0, 0.0, None)));
    Ok(EdgePreparation { cycles, intervals })
}
