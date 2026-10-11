//! Consumer-owned source support; engine dispatch stays out of public traits.
use super::{PersistenceBuilder, PersistenceOptions, flag, rips, simplicial};
use crate::complex::SimplicialComplex;
use crate::diagram::PersistenceResult;
use crate::execution::WorkBudget;
use crate::filtration::{
    ApproximateRipsBuilder, FlagFiltration, RipsBuilder, RipsExpansion, SimplicialFiltration,
    SimplicialStage, SparseRips, SparseRipsExpansion, ThresholdRips,
};
use crate::persistence::RepresentativeRequest;
use crate::{Error, Result};

pub trait Sealed {}

/// Add a configurable persistence operation to supported library sources.
///
/// Import this trait to use `.persistence()`. The method neither computes nor
/// mutates a cache. Implementations are sealed; implicit simplex access is private.
///
/// Unexpanded sources deliberately provide no simplex queries:
/// ```compile_fail
/// use cocycle::{filtration::RipsBuilder, geometry::PointCloudView};
/// let data = [0., 1.];
/// let points = PointCloudView::new(&data, 2, 1).unwrap();
/// RipsBuilder::from_points(points).find(&[0, 1]);
/// ```
/// Stateful callbacks must be consumed into a prepared source or explicit complex:
/// ```compile_fail
/// use cocycle::{filtration::RipsBuilder, persistence::PersistenceExt};
/// let items = [0., 1.];
/// let callback = RipsBuilder::from_distance_fn(&items, |_, _| Ok(1.));
/// callback.persistence();
/// ```
/// Analysis requests cannot outlive their source:
/// ```compile_fail
/// use cocycle::{filtration::RipsBuilder, geometry::PointCloudView, persistence::PersistenceExt};
/// let values = [0., 1.];
/// let points = PointCloudView::new(&values, 2, 1).unwrap();
/// let request;
/// { let rips = RipsBuilder::from_points(points); request = rips.persistence(); }
/// request.compute();
/// ```
pub trait PersistenceExt: Sealed + Sized {
    /// Borrow this source for analysis; expensive work begins only at `compute`.
    fn persistence(&self) -> PersistenceBuilder<'_, 'static, Self>;
}
macro_rules! source {
    ($ty:ty, $compute:path) => {
        impl Sealed for $ty {}
        impl PersistenceExt for $ty {
            fn persistence(&self) -> PersistenceBuilder<'_, 'static, Self> {
                PersistenceBuilder::new(self, $compute)
            }
        }
    };
}
source!(RipsBuilder<'_>, rips::compute_builder_budget);
source!(
    ApproximateRipsBuilder<'_>,
    rips::approximation::compute_builder_budget
);
source!(ThresholdRips, rips::compute_threshold_rips_budget);
source!(FlagFiltration, flag::compute_flag_budget);
source!(SparseRips, rips::approximation::compute_sparse_rips_budget);
source!(RipsExpansion, rips::expanded::compute_expanded_rips_budget);
source!(
    SparseRipsExpansion,
    rips::approximation::compute_expanded_sparse_rips_budget
);
impl Sealed for SimplicialFiltration {}
impl PersistenceExt for SimplicialFiltration {
    fn persistence(&self) -> PersistenceBuilder<'_, 'static, Self> {
        PersistenceBuilder::new_filtered(self, simplicial::compute_filtration_budget)
    }
}
impl Sealed for SimplicialComplex {}
impl PersistenceExt for SimplicialComplex {
    fn persistence(&self) -> PersistenceBuilder<'_, 'static, Self> {
        PersistenceBuilder::new_filtered(self, simplicial::compute_supplied_budget)
    }
}
impl Sealed for SimplicialStage<'_> {}
impl PersistenceExt for SimplicialStage<'_> {
    fn persistence(&self) -> PersistenceBuilder<'_, 'static, Self> {
        PersistenceBuilder::new_filtered(self, stage)
    }
}

fn stage(
    input: &SimplicialStage<'_>,
    options: &PersistenceOptions,
    requests: &[RepresentativeRequest],
    budget: &mut WorkBudget<'_>,
) -> Result<PersistenceResult> {
    let cutoff = options.max_edge().unwrap_or(input.scale());
    if cutoff > input.scale() {
        return Err(Error::IncompleteFiltration {
            requested: cutoff,
            through: input.scale(),
        });
    }
    let effective =
        PersistenceOptions::for_filtration(options.max_homology_dimension(), Some(cutoff))?
            .with_field(options.field());
    match input.filtration {
        Some(source) => simplicial::compute_filtration_budget(source, &effective, requests, budget),
        None => simplicial::compute_supplied_budget(input.complex, &effective, requests, budget),
    }
}
