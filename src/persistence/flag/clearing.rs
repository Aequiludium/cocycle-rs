//! H1 pivot handoff to H2, including zero-lifetime and omitted apparent pairs.
//! A pivot is a death simplex, not an H2 birth. Diagram filtering must not erase
//! this information. H1-only callers use a zero-sized sink with no allocations.
use crate::execution::WorkBudget;
use crate::filtration::flag::{SimplexEntry, SimplexIndex};
use crate::{Error, Result};

pub(super) trait TriangleClearing {
    fn record(&mut self, triangle: SimplexEntry) -> Result<()>;
    #[cfg(test)]
    fn storage(&self) -> (usize, usize) {
        (0, 0)
    }
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
    #[cfg(test)]
    fn storage(&self) -> (usize, usize) {
        (
            self.ids.len(),
            self.ids.capacity() * std::mem::size_of::<usize>(),
        )
    }
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

    pub(super) fn into_triangles(
        self,
        vertices: usize,
        budget: &mut WorkBudget<'_>,
    ) -> Result<Vec<[usize; 3]>> {
        budget.check()?;
        let index = SimplexIndex::new(vertices)?;
        let mut triangles = Vec::new();
        triangles
            .try_reserve_exact(self.ids.len())
            .map_err(|_| allocation())?;
        for id in self.ids {
            budget.step()?;
            triangles.push(index.triangle_vertices(id));
        }
        budget.check()?;
        Ok(triangles)
    }
}
fn allocation() -> Error {
    Error::AllocationFailed {
        context: "H1 triangle clearing",
    }
}
