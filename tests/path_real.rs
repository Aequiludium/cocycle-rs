//! Public numerical path-chain contracts and independent ambient identities.
use cocycle::Error;
use cocycle::complex::{DirectedGraph, PathBoundaryMode as Mode};
use cocycle::execution::Execution;
use cocycle::hodge::{
    PathRealChain, compute_persistent_path_laplacian, compute_persistent_path_laplacian_with,
    compute_real_path_chain, compute_real_path_chain_with,
};
use std::collections::BTreeMap;
use std::sync::atomic::AtomicBool;

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn apply(columns: &[Vec<f64>], v: &[f64], rows: usize) -> Vec<f64> {
    (0..rows)
        .map(|i| columns.iter().zip(v).map(|(c, s)| c[i] * s).sum())
        .collect()
}
fn check(chain: &PathRealChain, d: usize) {
    assert_eq!(chain.tolerance(), 1e-10);
    assert!(chain.boundary_residual().is_finite() && chain.boundary_residual() < 1e-8);
    for q in 0..=d + 1 {
        let paths = chain.paths(q).unwrap();
        let basis = chain.bases(q).unwrap();
        for (j, v) in basis.iter().enumerate() {
            assert_eq!(v.len(), paths.len());
            assert!((dot(v, v) - 1.0).abs() < 1e-8);
            for u in &basis[..j] {
                assert!(dot(u, v).abs() < 1e-8);
            }
            let mut ambient = BTreeMap::new();
            for (word, s) in paths.iter().zip(v) {
                if word.len() > 1 {
                    for i in 0..word.len() {
                        let mut face = word.clone();
                        face.remove(i);
                        if chain.mode() == Mode::Regular && face.windows(2).any(|v| v[0] == v[1]) {
                            continue;
                        }
                        *ambient.entry(face).or_insert(0.0) +=
                            s * if i % 2 == 0 { 1.0 } else { -1.0 };
                    }
                }
            }
            if q == 0 {
                assert!(chain.boundaries(q).unwrap()[j].is_empty());
                continue;
            }
            let lower = chain.paths(q - 1).unwrap();
            for (face, value) in &ambient {
                if !lower.contains(face) {
                    assert!(value.abs() < 1e-8);
                }
            }
            let actual = apply(
                chain.bases(q - 1).unwrap(),
                &chain.boundaries(q).unwrap()[j],
                lower.len(),
            );
            for (face, x) in lower.iter().zip(&actual) {
                assert!((ambient.get(face).copied().unwrap_or(0.0) - x).abs() < 1e-8);
            }
            let rows = if q == 1 {
                0
            } else {
                chain.bases(q - 2).unwrap().len()
            };
            assert!(
                apply(
                    chain.boundaries(q - 1).unwrap(),
                    &chain.boundaries(q).unwrap()[j],
                    rows
                )
                .iter()
                .all(|x| x.abs() < 1e-8)
            );
        }
    }
}

#[test]
fn both_real_models_satisfy_independent_ambient_chain_identities() {
    let pairs = [(0, 1), (0, 2), (1, 0), (1, 2), (2, 0), (2, 1)];
    for mask in 0..64 {
        let edges: Vec<_> = pairs
            .iter()
            .enumerate()
            .filter_map(|(i, &e)| (mask & (1 << i) != 0).then_some(e))
            .collect();
        let source = DirectedGraph::new(3, edges)
            .unwrap()
            .allowed_paths(3, 1000)
            .unwrap();
        for mode in [Mode::Regular, Mode::NonRegular] {
            let chain = compute_real_path_chain(&source, 2, mode, 100000, 1e-10).unwrap();
            check(&chain, 2);
            for q in 0..=2 {
                let operator = chain.laplacian(q).unwrap();
                let n = operator.basis().len();
                let mut matrix = Vec::new();
                for i in 0..n {
                    let mut v = vec![0.0; n];
                    v[i] = 1.0;
                    matrix.push(operator.apply(&v).unwrap());
                }
                for (i, column) in matrix.iter().enumerate() {
                    for (j, value) in column.iter().enumerate() {
                        assert!((value - matrix[j][i]).abs() < 1e-8);
                    }
                }
                let vector: Vec<_> = (0..n).map(|i| i as f64 - 1.5).collect();
                assert!(dot(&vector, &operator.apply(&vector).unwrap()) >= -1e-8);
            }
        }
    }
}

#[test]
fn hand_derived_edge_and_reciprocal_operators_keep_the_real_model() {
    for mode in [Mode::Regular, Mode::NonRegular] {
        let source = DirectedGraph::new(2, vec![(0, 1)])
            .unwrap()
            .allowed_paths(2, 10)
            .unwrap();
        let chain = compute_real_path_chain(&source, 1, mode, 1000, 1e-10).unwrap();
        assert_eq!(
            chain.laplacian(0).unwrap().apply(&[1.0, 0.0]).unwrap(),
            [1.0, -1.0]
        );
        assert_eq!(chain.laplacian(1).unwrap().apply(&[3.0]).unwrap(), [6.0]);
        let source = DirectedGraph::new(2, vec![(0, 1), (1, 0)])
            .unwrap()
            .allowed_paths(2, 10)
            .unwrap();
        let operator = compute_real_path_chain(&source, 1, mode, 1000, 1e-10)
            .unwrap()
            .laplacian(1)
            .unwrap();
        let image = operator.apply(&[1.0, 1.0]).unwrap();
        let expected = if mode == Mode::Regular { 4.0 } else { 0.0 };
        assert!(image.iter().all(|x| (x - expected).abs() < 1e-8));
    }
}

#[test]
fn persistent_real_domains_agree_on_the_diagonal_and_preserve_owned_coordinates() {
    let edges = [(0, 1, 1.0), (1, 2, 1.0), (2, 0, 1.0), (0, 2, 2.0)];
    for mode in [Mode::Regular, Mode::NonRegular] {
        for event in 0..=2 {
            for q in 0..=1 {
                let pairs: Vec<_> = edges
                    .iter()
                    .filter(|e| e.2 <= event as f64)
                    .map(|e| (e.0, e.1))
                    .collect();
                let source = DirectedGraph::new(3, pairs)
                    .unwrap()
                    .allowed_paths(q + 1, 100)
                    .unwrap();
                let diagonal = compute_real_path_chain(&source, q, mode, 10000, 1e-10)
                    .unwrap()
                    .laplacian(q)
                    .unwrap();
                let persistent = compute_persistent_path_laplacian(
                    3,
                    &edges,
                    q,
                    (event, event),
                    mode,
                    (100, 10000),
                    1e-10,
                )
                .unwrap();
                assert_eq!(persistent.basis().len(), diagonal.basis().len());
                let n = diagonal.basis().len();
                for j in 0..n {
                    let mut v = vec![0.0; n];
                    v[j] = 1.0;
                    let a = diagonal.apply(&v).unwrap();
                    let b = persistent.apply(&v).unwrap();
                    assert!(a.iter().zip(b).all(|(a, b)| (a - b).abs() < 1e-8));
                }
            }
        }
    }
    let operator =
        compute_persistent_path_laplacian(3, &edges, 1, (1, 2), Mode::Regular, (100, 10000), 1e-10)
            .unwrap();
    assert_eq!(operator.times(), (1.0, 2.0));
    assert_eq!(operator.tolerance(), 1e-10);
    assert!(operator.boundary_residual().is_finite() && operator.boundary_residual() < 1e-8);
    assert_eq!(operator.paths().len(), 4);
    let v = vec![1.0; operator.basis().len()];
    assert!(operator.apply(&v).unwrap().iter().all(|x| x.is_finite()));
}

#[test]
fn real_omega_is_not_lifted_from_the_characteristic_two_basis() {
    let edges = vec![
        (0, 1),
        (0, 2),
        (0, 3),
        (1, 4),
        (2, 5),
        (3, 6),
        (1, 5),
        (2, 6),
        (3, 4),
        (4, 7),
        (5, 8),
        (6, 9),
        (4, 8),
        (5, 9),
        (6, 7),
        (7, 10),
        (8, 10),
        (9, 10),
        (1, 7),
        (2, 8),
        (3, 9),
        (1, 9),
        (2, 7),
        (3, 8),
    ];
    let source = DirectedGraph::new(11, edges)
        .unwrap()
        .allowed_paths(4, 10000)
        .unwrap();
    let chain = compute_real_path_chain(&source, 3, Mode::Regular, 1000000, 1e-10).unwrap();
    assert!(chain.bases(4).unwrap().is_empty());
    check(&chain, 3);
}

#[test]
fn real_limits_nonfinite_inputs_truncation_and_cancellation_fail_without_partial_results() {
    let graph = DirectedGraph::new(2, vec![(0, 1), (1, 0)]).unwrap();
    let short = graph.allowed_paths(1, 10).unwrap();
    assert!(compute_real_path_chain(&short, 1, Mode::Regular, 1000, 1e-10).is_err());
    let source = graph.allowed_paths(2, 10).unwrap();
    for tolerance in [0.0, -1.0, 1e-3, f64::NAN, f64::INFINITY] {
        assert!(compute_real_path_chain(&source, 1, Mode::Regular, 1000, tolerance).is_err());
    }
    assert!(compute_real_path_chain(&source, 1, Mode::Regular, 1, 1e-10).is_err());
    let flag = AtomicBool::new(true);
    assert!(matches!(
        compute_real_path_chain_with(
            &source,
            1,
            Mode::Regular,
            1000,
            1e-10,
            &Execution::default().cancellation(&flag)
        ),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        compute_real_path_chain_with(
            &source,
            1,
            Mode::Regular,
            1000,
            1e-10,
            &Execution::default().max_work(100)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let operator = compute_real_path_chain(&source, 1, Mode::Regular, 1000, 1e-10)
        .unwrap()
        .laplacian(1)
        .unwrap();
    assert!(operator.apply(&[1.0]).is_err());
    assert!(operator.apply(&[f64::NAN, 0.0]).is_err());
    assert!(operator.apply(&[f64::MAX, f64::MAX]).is_err());
    assert!(matches!(
        operator.apply_with(&[1.0, 0.0], &Execution::default().max_work(1)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(
        compute_persistent_path_laplacian(
            2,
            &[(0, 1, 1.0)],
            1,
            (1, 0),
            Mode::Regular,
            (100, 1000),
            1e-10
        )
        .is_err()
    );
    assert!(matches!(
        compute_persistent_path_laplacian_with(
            2,
            &[(0, 1, 1.0)],
            1,
            (0, 1),
            Mode::Regular,
            (100, 1000),
            1e-10,
            &Execution::default().max_work(1)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
}
