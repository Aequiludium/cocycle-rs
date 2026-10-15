//! Fixed public-API Dirac fixtures for the independent Python checker.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::spectral::DiracOperator;
fn main() -> cocycle::Result<()> {
    let cases: Vec<(&str, Vec<Vec<usize>>, usize, usize)> = vec![
        ("empty", vec![], 0, 2),
        ("isolates", vec![vec![5], vec![42]], 0, 2),
        ("edge", vec![vec![5, 42]], 0, 1),
        ("edge_only_c1", vec![vec![5, 42]], 1, 1),
        ("edge_missing_degrees", vec![vec![5, 42]], 3, 5),
        ("path", vec![vec![2, 9], vec![9, 30]], 0, 1),
        (
            "triangle_ring",
            vec![vec![0, 1], vec![0, 2], vec![1, 2]],
            0,
            2,
        ),
        ("filled_triangle", vec![vec![0, 1, 2]], 0, 2),
        ("filled_triangle_c1_c2", vec![vec![0, 1, 2]], 1, 2),
        ("tetrahedron", vec![vec![0, 1, 2, 3]], 0, 3),
        (
            "tetrahedron_boundary",
            vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3]],
            0,
            2,
        ),
    ];
    for (name, facets, a, b) in cases {
        let mut all = std::collections::BTreeSet::new();
        for f in facets {
            for mask in 1..(1 << f.len()) {
                all.insert(
                    f.iter()
                        .enumerate()
                        .filter_map(|(i, &v)| (mask & (1 << i) != 0).then_some(v))
                        .collect::<Vec<_>>(),
                );
            }
        }
        let source = SimplicialComplex::new(
            all.into_iter()
                .map(|v| Simplex::new(v, 0.))
                .collect::<cocycle::Result<_>>()?,
        )?;
        let op = DiracOperator::new(&source, a..=b)?;
        let basis: Vec<_> = op
            .basis()
            .iter()
            .map(|&id| source.simplex(id).unwrap().vertices().to_vec())
            .collect();
        let x: Vec<_> = (0..basis.len()).map(|i| (i as f64 - 3.) / 4.).collect();
        let y = op.action(&x)?;
        println!(
            "{{\"name\":\"{name}\",\"basis\":{basis:?},\"offsets\":{:?},\"truncated\":{},\"matrix\":{:?},\"x\":{x:?},\"action\":{y:?},\"square_action\":{:?}}}",
            op.degree_offsets(),
            op.is_truncated(),
            op.matrix(4096)?,
            op.action(&y)?
        );
    }
    Ok(())
}
