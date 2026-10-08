//! H1 pivot handoff to H2, including zero-lifetime and omitted apparent pairs.
//! A pivot is a death simplex, not an H2 birth. Diagram filtering must not erase
//! this information. H1-only callers use a zero-sized sink with no allocations.
use crate::execution::WorkBudget;
use crate::filtration::flag::{SimplexEntry, SimplexIndex};
use crate::{Error, Result};
use std::collections::HashSet;

pub(super) trait TriangleClearing {
    fn record(&mut self, triangle: SimplexEntry) -> Result<()>;
}

pub(super) struct IgnoreTriangles;
impl TriangleClearing for IgnoreTriangles {
    #[inline]
    fn record(&mut self, _triangle: SimplexEntry) -> Result<()> {
        Ok(())
    }
}

#[derive(Default)]
pub(super) struct TrianglePivots {
    ids: Vec<usize>,
}
impl TriangleClearing for TrianglePivots {
    fn record(&mut self, triangle: SimplexEntry) -> Result<()> {
        self.ids.try_reserve(1).map_err(|_| allocation())?;
        self.ids.push(triangle.id);
        Ok(())
    }
}
impl TrianglePivots {
    pub(super) fn is_empty(&self) -> bool {
        self.ids.is_empty()
    }

    pub(super) fn into_vertices(
        self,
        vertices: usize,
        budget: &mut WorkBudget<'_>,
    ) -> Result<HashSet<Vec<usize>>> {
        budget.check()?;
        let index = SimplexIndex::new(vertices)?;
        let mut cleared = HashSet::new();
        cleared
            .try_reserve(self.ids.len())
            .map_err(|_| allocation())?;
        for id in self.ids {
            budget.step()?;
            let mut simplex = Vec::new();
            simplex.try_reserve_exact(3).map_err(|_| allocation())?;
            simplex.extend(index.triangle_vertices(id));
            if !cleared.insert(simplex) {
                return Err(Error::InternalInvariant {
                    reason: "H1 clearing pivots must be unique",
                });
            }
        }
        budget.check()?;
        Ok(cleared)
    }
}
fn allocation() -> Error {
    Error::AllocationFailed {
        context: "H1 triangle clearing",
    }
}
