//! Weighted directed path persistence using shared diagrams and descriptors.
use cocycle::algebra::PrimeField;
use cocycle::complex::PathBoundaryMode;
use cocycle::persistence::compute_path_persistence;

fn main() -> cocycle::Result<()> {
    let edges = [(0, 1, 1.0), (1, 2, 1.0), (2, 0, 1.0), (0, 2, 2.0)];
    let result = compute_path_persistence(
        3,
        &edges,
        1,
        PrimeField::new(3)?,
        PathBoundaryMode::Regular,
        100,
        10_000,
    )?;
    let curve = cocycle::descriptors::betti_curve(result.diagram(), 1, &[0.0, 1.0, 2.0])?;
    assert_eq!(curve, [0, 1, 0]);
    println!("H1 Betti curve: {curve:?}");
    Ok(())
}
