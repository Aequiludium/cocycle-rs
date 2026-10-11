//! Explicit simplicial filtration with certified source coverage.
use crate::complex::{BoundaryTerm, Simplex, SimplexId, SimplicialComplex};
use crate::filtration::{Coverage, FiltrationContext};
use crate::{Error, Result};
/// Owned, frozen filtered topology with certified scale and dimension coverage.
/// Only an explicit construction creates this object; arbitrary topology cannot
/// be assigned an unchecked original-Rips completeness claim.
#[derive(Clone, Debug)]
pub struct SimplicialFiltration {
    pub(crate) complex: SimplicialComplex,
    pub(crate) context: FiltrationContext,
    pub(crate) coverage: Coverage,
    pub(crate) dimension: usize,
    pub(crate) complete: bool,
}
impl SimplicialFiltration {
    /// Borrow stored simplices, values, lookup and oriented incidence.
    pub fn complex(&self) -> &SimplicialComplex {
        &self.complex
    }
    /// Original construction facts, independent of any later analysis.
    pub fn context(&self) -> &FiltrationContext {
        &self.context
    }
    /// Certified range of this mathematical source.
    pub fn coverage(&self) -> Coverage {
        self.coverage
    }
    /// Requested maximum simplex dimension, not the homology dimension.
    pub fn max_simplex_dimension(&self) -> usize {
        self.dimension
    }
    /// Whether expansion exhausted the source's simplices within its scale range.
    pub fn is_dimension_complete(&self) -> bool {
        self.complete
    }
    /// Borrow the closed sublevel set at a finite, possibly negative scale.
    ///
    /// Retains this source's coverage and dimension certificates. Selection
    /// neither allocates nor computes persistence; all tied cells enter together.
    /// # Errors
    /// Rejects nonfinite scales and scales beyond certified construction coverage.
    pub fn stage(&self, scale: f64) -> Result<SimplicialStage<'_>> {
        SimplicialStage::new(&self.complex, Some(self), scale)
    }
}

/// A borrowed sublevel set of one frozen simplicial source.
///
/// Contains exactly the stored simplices with value <= [`Self::scale`], in the
/// original filtration order. IDs, vertex orientation and signed incidence are
/// unchanged. Selection takes O(log n) time and constant space; iteration borrows
/// a prefix and does not copy simplices or construct missing dimensions.
///
/// [`crate::persistence::PersistenceExt::persistence`] analyzes the original
/// source through this scale, retaining its coverage, context and skeleton
/// checks. In contrast, [`crate::persistence::PersistenceBuilder::from_complex`]
/// explicitly treats the stage as a supplied complex in its own right, discarding
/// the relationship to later events or a larger construction.
///
/// Only same-source inclusions are supported, not arbitrary chain maps,
/// relabeling, multi-parameter persistence or a persistent spectral operator.
#[derive(Clone, Copy, Debug)]
pub struct SimplicialStage<'a> {
    pub(crate) complex: &'a SimplicialComplex,
    pub(crate) filtration: Option<&'a SimplicialFiltration>,
    scale: f64,
    end: usize,
}
impl<'a> SimplicialStage<'a> {
    pub(crate) fn new(
        complex: &'a SimplicialComplex,
        filtration: Option<&'a SimplicialFiltration>,
        scale: f64,
    ) -> Result<Self> {
        if !scale.is_finite() {
            return Err(Error::NonFiniteValue {
                field: "stage scale",
                index: None,
            });
        }
        let scale = crate::canonical_zero(scale);
        if let Some(source) = filtration
            && let Coverage::Through(through) = source.coverage()
            && scale > through
        {
            return Err(Error::IncompleteFiltration {
                requested: scale,
                through,
            });
        }
        let end = complex.simplices().partition_point(|s| s.value() <= scale);
        Ok(Self {
            complex,
            filtration,
            scale,
            end,
        })
    }
    /// Inclusive finite sublevel scale, in the original source's units.
    pub fn scale(&self) -> f64 {
        self.scale
    }
    /// Number of simplices in this stage, including its vertices.
    pub fn len(&self) -> usize {
        self.end
    }
    /// Whether no stored simplex has entered at this scale.
    pub fn is_empty(&self) -> bool {
        self.end == 0
    }
    /// Borrow the stage's simplices in original filtration order.
    /// Slice positions remain source-local [`SimplexId`] indices.
    pub fn simplices(&self) -> &'a [Simplex] {
        &self.complex.simplices()[..self.end]
    }
    /// Look up a simplex that has entered this stage, without renumbering it.
    /// Invalid vertex order or absent simplices return `None`.
    pub fn find(&self, vertices: &[usize]) -> Option<SimplexId> {
        self.complex
            .find(vertices)
            .filter(|id| id.index() < self.end)
    }
    /// Borrow a simplex at an ID from this source, or `None` if it is not in the stage.
    /// IDs do not encode ownership; never pass handles from another complex.
    pub fn simplex(&self, id: SimplexId) -> Option<&'a Simplex> {
        if id.index() < self.end {
            self.complex.simplex(id)
        } else {
            None
        }
    }
    /// Borrow signed boundary terms, or `None` for a cell outside this stage.
    /// Face monotonicity guarantees that every returned face is also in the stage.
    pub fn boundary(&self, id: SimplexId) -> Option<&'a [BoundaryTerm]> {
        if id.index() < self.end {
            self.complex.boundary(id)
        } else {
            None
        }
    }
    /// Original full storage, including cells beyond this stage.
    /// Using it as a supplied source explicitly discards construction certificates.
    pub fn source_complex(&self) -> &'a SimplicialComplex {
        self.complex
    }
    /// Original source certificates, when selected from a [`SimplicialFiltration`].
    /// `None` means the supplied complex itself is the complete mathematical source;
    /// it makes no claim about a larger Rips, Alpha or other construction.
    pub fn source_filtration(&self) -> Option<&'a SimplicialFiltration> {
        self.filtration
    }
    /// Read the identity inclusion into a later stage of this same source.
    ///
    /// Each pair is `(source_id, target_id)`; the two IDs are identical. The
    /// map preserves dimension, orientation and boundary over the integers and
    /// every prime field. Checking the stages takes constant time; traversing
    /// the map takes O(self.len()) time with no allocation. Equal scales are valid.
    ///
    /// # Errors
    /// Rejects decreasing scales, distinct source owners (even equal clones),
    /// and mixing certified and bare storage interpretations of the same object.
    pub fn inclusion_into(
        &self,
        target: &Self,
    ) -> Result<impl ExactSizeIterator<Item = (SimplexId, SimplexId)> + use<>> {
        let same_certificate = match (self.filtration, target.filtration) {
            (None, None) => true,
            (Some(a), Some(b)) => std::ptr::eq(a, b),
            _ => false,
        };
        if !std::ptr::eq(self.complex, target.complex) || !same_certificate {
            return Err(Error::InvalidParameter {
                parameter: "stage inclusion",
                reason: "stages must borrow the same source and certificates",
            });
        }
        if self.scale > target.scale {
            return Err(Error::InvalidParameter {
                parameter: "stage inclusion",
                reason: "target scale must not precede source scale",
            });
        }
        Ok((0..self.end).map(|index| (SimplexId(index), SimplexId(index))))
    }
}

mod access;
pub(crate) use access::{ZeroBornExplicitAccess, ZeroBornSimplicialAccess, next_dimension};
