//! Cross-domain identities in common owned simplex coordinates.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::dirac::DiracOperator;
use cocycle::execution::Execution;
use cocycle::hodge::{HodgeLaplacian, HodgeSelection};
use cocycle::persistent_laplacian::{PersistentLaplacian, SpectralOptions};

#[test]
fn dirac_square_hodge_and_same_stage_persistent_agree() -> cocycle::Result<()> {
    let source = SimplicialComplex::new(
        [
            (vec![9], 0.),
            (vec![2], 0.),
            (vec![17], 0.),
            (vec![2, 9], 1.),
            (vec![9, 17], 1.),
            (vec![2, 17], 1.),
            (vec![2, 9, 17], 2.),
        ]
        .into_iter()
        .map(|(vertices, birth)| Simplex::new(vertices, birth))
        .collect::<cocycle::Result<Vec<_>>>()?,
    )?;
    let execution = Execution::default();
    for scale in [-1., 0., 1., 2.] {
        let stage = source.stage(scale)?;
        let dirac = DiracOperator::new(&stage, 0..=2)?;
        let n = dirac.basis().len();
        let matrix = dirac.matrix(n * n, &execution)?;
        let square: Vec<f64> = (0..n * n)
            .map(|entry| {
                let (i, j) = (entry / n, entry % n);
                (0..n).map(|k| matrix[i * n + k] * matrix[k * n + j]).sum()
            })
            .collect();
        for degree in 0..=2 {
            let hodge = HodgeLaplacian::new(&stage, degree, HodgeSelection::Full)?;
            let ordinary = hodge.matrix(&execution)?;
            let persistent = PersistentLaplacian::compute(
                &stage,
                &stage,
                degree,
                SpectralOptions::default(),
                &execution,
            )?;
            let offsets = dirac.degree_offsets();
            let start = offsets[degree];
            let end = offsets[degree + 1];
            assert_eq!(&dirac.basis()[start..end], hodge.basis());
            for (id, simplex) in hodge.basis().iter().zip(persistent.basis()) {
                assert_eq!(stage.simplex(*id).unwrap(), simplex);
            }
            for i in start..end {
                for j in 0..n {
                    let expected = if (start..end).contains(&j) {
                        ordinary[(i - start) * (end - start) + j - start]
                    } else {
                        0.
                    };
                    assert!((square[i * n + j] - expected).abs() < 1e-12);
                }
            }
            assert_eq!(ordinary.len(), persistent.matrix().len());
            for (a, b) in ordinary.iter().zip(persistent.matrix()) {
                assert!((a - b).abs() < 1e-12);
            }
        }
    }
    Ok(())
}
