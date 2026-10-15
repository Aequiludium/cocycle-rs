//! Reproducible standard real Hodge fixtures; stdout is one JSON record per window.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::hodge::{HodgeLaplacian, HodgeSelection, IterativeOptions, SpectrumOptions};
use std::collections::BTreeSet;

fn from_facets(facets: &[Vec<usize>]) -> cocycle::Result<SimplicialComplex> {
    let mut simplices = BTreeSet::new();
    for facet in facets {
        for mask in 1..(1_usize << facet.len()) {
            simplices.insert(
                facet
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| mask & (1 << i) != 0)
                    .map(|(_, &v)| v)
                    .collect::<Vec<_>>(),
            );
        }
    }
    SimplicialComplex::new(
        simplices
            .into_iter()
            .map(|v| Simplex::new(v, 0.))
            .collect::<cocycle::Result<_>>()?,
    )
}

fn emit(name: &str, source: &SimplicialComplex) -> cocycle::Result<()> {
    let simplices: Vec<_> = source.simplices().iter().map(|s| s.vertices()).collect();
    for degree in 0..=source.dimension().unwrap_or(0) + 1 {
        let laplacian = HodgeLaplacian::new(source, degree, HodgeSelection::Full)?;
        let basis: Vec<_> = laplacian
            .basis()
            .iter()
            .map(|&id| source.simplex(id).unwrap().vertices())
            .collect();
        let n = basis.len();
        let spectrum =
            laplacian.spectrum(0..n, &SpectrumOptions::default(), &Execution::default())?;
        let x: Vec<_> = (0..n)
            .map(|i| (i as f64 + 1.) * (-1.0_f64).powi((i % 2) as i32))
            .collect();
        let action = laplacian.action(&x, &Execution::default())?;
        let kernel_columns: Vec<_> = spectrum.numerical_kernel_columns(1e-10)?.collect();
        let iterative = laplacian.smallest_eigenpairs(
            n.min(3),
            &IterativeOptions::default(),
            &Execution::default(),
        )?;
        println!(
            "{{\"name\":{name:?},\"degree\":{degree},\"simplices\":{simplices:?},\"basis\":{basis:?},\"matrix\":{:?},\"values\":{:?},\"vectors\":{:?},\"residuals\":{:?},\"orthogonality_error\":{},\"x\":{x:?},\"action\":{action:?},\"kernel_columns\":{kernel_columns:?},\"sweeps\":{},\"iterative_values\":{:?},\"iterative_vectors\":{:?},\"iterative_residuals\":{:?},\"iterative_iterations\":{},\"iterative_actions\":{},\"iterative_peak_subspace\":{}}}",
            spectrum.matrix(),
            spectrum.values(),
            spectrum.vectors(),
            spectrum.residuals(),
            spectrum.orthogonality_error(),
            spectrum.sweeps(),
            iterative.values(),
            iterative.vectors(),
            iterative.residuals(),
            iterative.iterations(),
            iterative.operator_applications(),
            iterative.peak_subspace_dimension()
        );
    }
    Ok(())
}

fn main() -> cocycle::Result<()> {
    for (name, facets) in [
        ("edge", vec![vec![0, 1]]),
        (
            "triangle_boundary",
            vec![vec![0, 1], vec![0, 2], vec![1, 2]],
        ),
        ("triangle_filled", vec![vec![0, 1, 2]]),
        ("path_relabelled", vec![vec![2, 19], vec![19, 103]]),
        ("isolated_vertices", vec![vec![1], vec![20]]),
        (
            "tetrahedron_surface",
            vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3]],
        ),
        ("tetrahedron_filled", vec![vec![0, 1, 2, 3]]),
        ("disconnected", vec![vec![0, 1, 2], vec![7, 8], vec![10]]),
        // Six-vertex triangulation already cited in the repository bibliography.
        (
            "rp2",
            vec![
                vec![0, 1, 2],
                vec![0, 1, 3],
                vec![0, 2, 4],
                vec![0, 3, 5],
                vec![0, 4, 5],
                vec![1, 2, 5],
                vec![1, 3, 4],
                vec![1, 4, 5],
                vec![2, 3, 4],
                vec![2, 3, 5],
            ],
        ),
        ("empty", vec![]),
        ("path_128", (0..127).map(|v| vec![v, v + 1]).collect()),
    ] {
        emit(name, &from_facets(&facets)?)?;
    }
    // A fixed small corpus with nontrivial and repeated spectra. The oracle
    // reads the emitted topology and builds its own boundaries independently.
    for seed in 0..12_u64 {
        let mut state = seed + 73;
        let mut facets = (0..6).map(|v| vec![v]).collect::<Vec<_>>();
        for i in 0..6 {
            for j in i + 1..6 {
                for k in j + 1..6 {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    if state >> 62 != 0 {
                        facets.push(vec![i, j, k]);
                    }
                }
            }
        }
        emit(&format!("seed_{seed}"), &from_facets(&facets)?)?;
    }
    Ok(())
}
