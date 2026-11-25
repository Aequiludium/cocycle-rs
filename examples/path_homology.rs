//! Owned exact path coordinates, full prime fields and boundary-model selection.
use cocycle::algebra::PrimeField;
use cocycle::complex::{DirectedGraph, PathBoundaryMode};
use cocycle::persistence::compute_path_homology;

fn main() -> cocycle::Result<()> {
    let source = DirectedGraph::new(3, vec![(0, 1), (1, 2), (2, 0)])?.allowed_paths(2, 100)?;
    let result = compute_path_homology(
        &source,
        1,
        PrimeField::new(4_294_967_291)?,
        PathBoundaryMode::Regular,
        10_000,
    )?;
    drop(source);
    assert_eq!(result.betti(), [1, 1]);
    let dual = result.cohomology(10_000)?;
    let projection = result.projection(10_000)?;
    assert_eq!(dual.representatives(1).unwrap().len(), 1);
    assert_eq!(
        projection.homotopies(1).unwrap().len(),
        result.omega_bases(1).unwrap().len()
    );
    assert_eq!(
        result.paths(1).unwrap(),
        [vec![0, 1], vec![1, 2], vec![2, 0]]
    );
    println!(
        "F{} Betti {:?}",
        result.field().characteristic(),
        result.betti()
    );
    Ok(())
}
