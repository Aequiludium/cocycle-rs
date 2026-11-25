//! Real path chains and factorized action without modular or integer-cell adapters.
use cocycle::complex::{DirectedGraph, PathBoundaryMode};
use cocycle::hodge::{compute_persistent_path_laplacian, compute_real_path_chain};

fn main() -> cocycle::Result<()> {
    let source = DirectedGraph::new(2, vec![(0, 1)])?.allowed_paths(2, 10)?;
    let chain = compute_real_path_chain(&source, 1, PathBoundaryMode::Regular, 1000, 1e-10)?;
    let operator = chain.laplacian(0)?;
    drop(chain);
    drop(source);
    assert_eq!(operator.apply(&[1.0, 0.0])?, [1.0, -1.0]);
    let weighted = [(0, 1, 1.0), (1, 2, 1.0), (2, 0, 1.0), (0, 2, 2.0)];
    let persistent = compute_persistent_path_laplacian(
        3,
        &weighted,
        1,
        (1, 2),
        PathBoundaryMode::Regular,
        (100, 10000),
        1e-10,
    )?;
    println!(
        "Events {:?}, Omega dimension {}, boundary residual {}",
        persistent.times(),
        persistent.basis().len(),
        persistent.boundary_residual()
    );
    Ok(())
}
