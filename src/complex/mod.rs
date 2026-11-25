//! Validated topology storage and a read-only filtered-cell contract.
//! A simplex is one cell; a simplicial complex is a face-closed collection.

mod graph;
pub use graph::{Neighbor, WeightedEdge, WeightedGraph};

mod simplicial;
pub use simplicial::{BoundaryTerm, Simplex, SimplexId, SimplicialComplex};

pub(crate) use simplicial::compare_filtration;

mod filtered;
pub use filtered::FilteredComplex;
mod path;
pub(crate) use path::boundary_budget as path_boundary_budget;
pub use path::{
    DirectedGraph, Path, PathBoundaryMode, PathComplex, path_boundary, path_boundary_with,
};
/// Compatibility name for [`SimplicialComplex`].
pub type FilteredSimplicialComplex = SimplicialComplex;
