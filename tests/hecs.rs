//! Public ordinary HeCS contracts with hand-derived simplicial controls.
use cocycle::Error;
use cocycle::complex::FilteredComplex;
use cocycle::execution::Execution;
use cocycle::hodge::{HecsOptions, HecsSolution, HodgeLaplacian, HodgeSelection};
use std::sync::atomic::AtomicBool;

struct Cells(Vec<(usize, Vec<(usize, i32)>)>);
impl FilteredComplex for Cells {
    type CellId = usize;
    fn cells(&self) -> impl Iterator<Item = usize> + '_ {
        0..self.0.len()
    }
    fn dimension(&self, id: usize) -> usize {
        self.0[id].0
    }
    fn filtration_value(&self, _: usize) -> f64 {
        17.
    } // never a face weight
    fn boundary(&self, id: usize) -> impl Iterator<Item = (usize, i32)> + '_ {
        self.0[id].1.iter().copied()
    }
}
fn cells(n: usize, edges: &[[usize; 2]], faces: &[[usize; 3]]) -> Cells {
    let mut result = Cells((0..n).map(|_| (0, vec![])).collect());
    for e in edges {
        result.0.push((1, vec![(e[0], -1), (e[1], 1)]));
    }
    for f in faces {
        let mut terms = Vec::new();
        for k in 0..3 {
            let edge: Vec<_> = f
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != k)
                .map(|(_, &v)| v)
                .collect();
            let e = edges.iter().position(|e| e.as_slice() == edge).unwrap();
            terms.push((n + e, if k % 2 == 0 { 1 } else { -1 }));
        }
        result.0.push((2, terms));
    }
    result
}
const EDGES: [[usize; 2]; 6] = [[0, 1], [0, 2], [0, 3], [1, 2], [1, 3], [2, 3]];
const FACES: [[usize; 3]; 4] = [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]];
fn solve(
    source: &Cells,
    weights: &[f64],
    rhs: &[f64],
    options: &HecsOptions,
) -> cocycle::Result<HecsSolution<usize>> {
    HodgeLaplacian::new(source, 1, HodgeSelection::Upper)?.solve_hecs(
        weights,
        rhs,
        options,
        &Execution::default(),
    )
}
fn close(a: &[f64], b: &[f64], tol: f64) {
    assert_eq!(a.len(), b.len());
    for (&a, &b) in a.iter().zip(b) {
        assert!((a - b).abs() < tol, "{a} != {b}");
    }
}

#[test]
fn weighted_shell_preserves_factor_skeleton_and_owned_coordinates() {
    let result = {
        let source = cells(4, &EDGES, &FACES);
        solve(
            &source,
            &[1., 2., 4., 8.],
            &[2., -22., 20., -22., 24., -44.],
            &HecsOptions::default(),
        )
        .unwrap()
    };
    assert_eq!(result.basis(), &[4, 5, 6, 7, 8, 9]);
    assert_eq!(result.triangle_basis(), &[10, 11, 12, 13]);
    assert_eq!(result.selected_triangles(), &[3, 2, 1]);
    assert_eq!(result.rank(), 3);
    assert!(result.residual_history().last().unwrap() <= &1e-9);
    let b = [
        [1., 1., 0., 0.],
        [-1., 0., 1., 0.],
        [0., -1., -1., 0.],
        [1., 0., 0., 1.],
        [0., 1., 0., -1.],
        [0., 0., 1., 1.],
    ];
    for i in 0..6 {
        for j in 0..6 {
            let expected: f64 = (1..4)
                .map(|t| {
                    b[result.edge_permutation()[i]][t]
                        * b[result.edge_permutation()[j]][t]
                        * result.weights()[t]
                })
                .sum();
            let actual: f64 = (0..3)
                .map(|k| result.factor()[i * 3 + k] * result.factor()[j * 3 + k])
                .sum();
            assert!((actual - expected).abs() < 1e-12);
        }
    }
    let z: Vec<_> = result
        .collapse_pairs()
        .iter()
        .map(|p| (if p[1] == 2 { -1. } else { 1. }) / result.weights()[p[1]].sqrt())
        .collect();
    for i in 0..3 {
        for j in 0..3 {
            assert!(
                (result.active_matrix()[i * 3 + j] - (f64::from(i == j) + z[i] * z[j])).abs()
                    < 1e-12
            );
        }
    }
    // Minimum-norm solution is orthogonal to all vertex gradients.
    let mut divergence = [0.; 4];
    for (e, &x) in EDGES.iter().zip(result.values()) {
        divergence[e[0]] -= x;
        divergence[e[1]] += x;
    }
    close(&divergence, &[0.; 4], 1e-11);
}

#[test]
fn orientations_and_zero_rhs_are_explicit() {
    let mut source = cells(
        3,
        &EDGES[..1]
            .iter()
            .copied()
            .chain([[0, 2], [1, 2]])
            .collect::<Vec<_>>(),
        &[[0, 1, 2]],
    );
    source.0[3].1.iter_mut().for_each(|(_, v)| *v = -*v);
    source.0[6].1[2].1 = -source.0[6].1[2].1; // e01 is the last boundary term
    let result = solve(&source, &[2.], &[-3., -3., 3.], &HecsOptions::default()).unwrap();
    close(result.values(), &[-0.5, -0.5, 0.5], 1e-12);
    let zero = solve(&source, &[2.], &[0.; 3], &HecsOptions::default()).unwrap();
    assert_eq!(zero.iterations(), 0);
    close(zero.values(), &[0.; 3], 1e-15);
    assert_eq!(zero.residual_history(), &[0.]);
}

#[test]
fn rank_loss_is_rejected_over_the_reals() {
    let edges: Vec<_> = (0..6)
        .flat_map(|i| (i + 1..6).map(move |j| [i, j]))
        .collect();
    let faces = [
        [0, 1, 2],
        [0, 1, 3],
        [0, 2, 4],
        [0, 3, 5],
        [0, 4, 5],
        [1, 2, 5],
        [1, 3, 4],
        [1, 4, 5],
        [2, 3, 4],
        [2, 3, 5],
    ];
    for peel in [false, true] {
        let error = solve(
            &cells(6, &edges, &faces),
            &[1.; 10],
            &[0.; 15],
            &HecsOptions {
                peel_core: peel,
                ..Default::default()
            },
        )
        .err()
        .unwrap();
        assert!(matches!(
            error,
            Error::InvalidParameter {
                parameter: "HeCS rank",
                ..
            }
        ));
    }
}

#[test]
fn empty_skeleton_incompatible_rhs_and_limits() {
    let options = HecsOptions::default();
    assert!(
        solve(&Cells(vec![]), &[], &[], &options)
            .unwrap()
            .values()
            .is_empty()
    );
    let skeleton = cells(3, &[[0, 1], [1, 2]], &[]);
    close(
        solve(&skeleton, &[], &[0.; 2], &options).unwrap().values(),
        &[0.; 2],
        1e-15,
    );
    assert!(matches!(
        solve(&skeleton, &[], &[1., 0.], &options),
        Err(Error::InvalidParameter {
            parameter: "HeCS right-hand side",
            ..
        })
    ));
    let source = cells(3, &[[0, 1], [0, 2], [1, 2]], &[[0, 1, 2]]);
    let op = HodgeLaplacian::new(&source, 1, HodgeSelection::Upper).unwrap();
    assert!(matches!(
        op.solve_hecs(
            &[1.],
            &[1., -1., 1.],
            &HecsOptions {
                max_dense_entries: 8,
                ..options
            },
            &Execution::default().max_work(0)
        ),
        Err(Error::InvalidParameter {
            parameter: "HeCS dense entries",
            ..
        })
    ));
    assert!(matches!(
        op.solve_hecs(
            &[1.],
            &[1., -1., 1.],
            &options,
            &Execution::default().max_work(20)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let flag = AtomicBool::new(true);
    assert!(matches!(
        op.solve_hecs(
            &[1.],
            &[1., -1., 1.],
            &options,
            &Execution::default().cancellation(&flag)
        ),
        Err(Error::Cancelled)
    ));
    assert!(
        solve(
            &source,
            &[1.],
            &[1., -1., 1.],
            &HecsOptions {
                max_iterations: 0,
                ..options
            }
        )
        .is_err()
    );
    assert!(solve(&source, &[1.], &[1., -1., 1.], &options).is_ok());
}

#[test]
fn rejects_invalid_weights_windows_and_cellular_topology() {
    let source = cells(3, &[[0, 1], [0, 2], [1, 2]], &[[0, 1, 2]]);
    let options = HecsOptions::default();
    for w in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(solve(&source, &[w], &[0.; 3], &options).is_err());
    }
    assert!(solve(&source, &[], &[0.; 3], &options).is_err());
    assert!(solve(&source, &[1.], &[0.; 2], &options).is_err());
    assert!(solve(&source, &[1.], &[f64::NAN; 3], &options).is_err());
    for tol in [0., -1., f64::NAN, f64::INFINITY] {
        assert!(
            solve(
                &source,
                &[1.],
                &[0.; 3],
                &HecsOptions {
                    residual_tolerance: tol,
                    ..options
                }
            )
            .is_err()
        );
    }
    let full = HodgeLaplacian::new(&source, 1, HodgeSelection::Full).unwrap();
    assert!(
        full.solve_hecs(&[1.], &[0.; 3], &options, &Execution::default())
            .is_err()
    );
    let parallel = cells(2, &[[0, 1], [0, 1]], &[]);
    assert!(matches!(
        solve(&parallel, &[], &[0.; 2], &options),
        Err(Error::InvalidParameter {
            parameter: "HeCS topology",
            ..
        })
    ));
    let mut multiple = source;
    multiple.0[6].1.iter_mut().for_each(|(_, v)| *v *= 2);
    assert!(matches!(
        solve(&multiple, &[1.], &[0.; 3], &options),
        Err(Error::InvalidParameter {
            parameter: "HeCS topology",
            ..
        })
    ));
}
