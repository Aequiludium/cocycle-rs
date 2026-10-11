use super::*;
use crate::complex::Simplex;
use std::collections::BTreeSet;
use std::sync::atomic::{AtomicBool, Ordering};

type Matrix = Vec<Vec<u32>>;

fn fixture(seed: usize, ties: bool) -> SimplicialComplex {
    // Complete 3-skeleton on five vertices: finite H0/H1/H2 and essential H3.
    let cells = (1usize..32)
        .filter(|mask| mask.count_ones() <= 4)
        .map(|mask| {
            let vertices: Vec<_> = (0..5).filter(|v| mask & (1 << v) != 0).collect();
            let jitter = if ties {
                (mask + seed) % 3
            } else {
                (mask * 17 + seed * 11) % 37
            };
            Simplex::new(vertices.clone(), (vertices.len() * 40 + jitter) as f64).unwrap()
        })
        .collect();
    SimplicialComplex::new(cells).unwrap()
}

// Independently enumerate incidence by vertex containment, without boundary().
fn boundary(complex: &SimplicialComplex, dual: bool) -> Matrix {
    let n = complex.len();
    let mut d = vec![vec![0; n]; n];
    for (j, coface) in complex.simplices().iter().enumerate() {
        for (i, face) in complex.simplices().iter().enumerate() {
            if coface.dimension() == face.dimension() + 1
                && face
                    .vertices()
                    .iter()
                    .all(|v| coface.vertices().contains(v))
            {
                let (row, col) = if dual { (n - 1 - j, n - 1 - i) } else { (i, j) };
                d[row][col] = 1;
            }
        }
    }
    d
}

fn dense_columns(columns: &[Column<usize>]) -> Matrix {
    let n = columns.len();
    let mut m = vec![vec![0; n]; n];
    for (j, column) in columns.iter().enumerate() {
        for (&i, &coefficient) in column.entries() {
            m[i][j] = coefficient;
        }
    }
    m
}

// Full Gauss-Jordan inverse: no triangular solve, sparse columns or reducer.
fn inverse(v: &Matrix) -> Matrix {
    let n = v.len();
    let mut a = v.clone();
    let mut u = vec![vec![0; n]; n];
    for (i, row) in u.iter_mut().enumerate() {
        row[i] = 1;
    }
    for i in 0..n {
        let pivot = (i..n).find(|&k| a[k][i] == 1).unwrap();
        a.swap(i, pivot);
        u.swap(i, pivot);
        for k in 0..n {
            if k != i && a[k][i] != 0 {
                for j in 0..n {
                    a[k][j] ^= a[i][j];
                    u[k][j] ^= u[i][j];
                }
            }
        }
    }
    u
}

// Pairing-only set reducer used after each Algorithm-2 block swap. No V or U.
fn lows(d: &Matrix, order: &[usize]) -> Vec<Option<usize>> {
    let mut owners: BTreeMap<usize, BTreeSet<usize>> = BTreeMap::new();
    order
        .iter()
        .map(|&j| {
            let mut c: BTreeSet<_> = (0..d.len()).filter(|&i| d[i][j] != 0).collect();
            while let Some(&pivot) = c.last() {
                if let Some(other) = owners.get(&pivot) {
                    c = c.symmetric_difference(other).copied().collect();
                } else {
                    owners.insert(pivot, c);
                    return Some(pivot);
                }
            }
            None
        })
        .collect()
}

fn swap_oracle(
    d: &Matrix,
    focus: usize,
    values: &[f64],
    target: f64,
    block: &[usize],
) -> Vec<usize> {
    let mut order = block.to_vec();
    let pivot = lows(d, &order)[order.iter().position(|&j| j == focus).unwrap()].unwrap();
    let mut critical = vec![focus];
    let candidates: Vec<_> = if target > values[focus] {
        block
            .iter()
            .copied()
            .filter(|&j| j > focus && values[j] <= target)
            .collect()
    } else {
        block
            .iter()
            .rev()
            .copied()
            .filter(|&j| j < focus && values[j] >= target)
            .collect()
    };
    for candidate in candidates {
        let original = order.clone();
        let position = order.iter().position(|&j| j == candidate).unwrap();
        order.remove(position);
        let positions: Vec<_> = critical
            .iter()
            .map(|id| order.iter().position(|j| j == id).unwrap())
            .collect();
        let insert = if target > values[focus] {
            *positions.iter().min().unwrap()
        } else {
            *positions.iter().max().unwrap() + 1
        };
        order.insert(insert, candidate);
        let paired = lows(d, &order)
            .iter()
            .position(|&low| low == Some(pivot))
            .unwrap();
        if order[paired] == candidate {
            critical.push(candidate);
            order = original;
        }
    }
    critical.sort_unstable();
    critical
}

#[test]
fn sparse_decompositions_and_bounded_u_match_full_inverse() {
    for seed in 0..8 {
        let complex = fixture(seed, seed % 2 == 0);
        let n = complex.len();
        for dual in [false, true] {
            let mut budget = WorkBudget::new(&Execution::default()).unwrap();
            let mut state = decompose(&complex, dual, &mut budget).unwrap();
            let d = boundary(&complex, dual);
            let r = dense_columns(&state.reduction.reduced);
            let v = dense_columns(&state.reduction.transforms);
            let u = inverse(&v);
            for i in 0..n {
                for j in 0..n {
                    assert_eq!((0..n).fold(0, |sum, k| sum ^ (d[i][k] & v[k][j])), r[i][j]);
                    assert_eq!(
                        (0..n).fold(0, |sum, k| sum ^ (u[i][k] & v[k][j])),
                        u32::from(i == j)
                    );
                    if i > j {
                        assert_eq!(v[i][j], 0);
                    }
                    if i < j && v[i][j] != 0 {
                        // Lazy reduction adds only negative columns whose low is
                        // greater than the resulting low (zero uses low=-1).
                        let low_i = (0..n).rev().find(|&k| r[k][i] != 0).unwrap();
                        assert!(
                            (0..n)
                                .rev()
                                .find(|&k| r[k][j] != 0)
                                .is_none_or(|low_j| low_i > low_j)
                        );
                    }
                }
                let source = if dual { n - 1 - i } else { i };
                let dimension = complex.simplices()[source].dimension();
                for offset in [0., 1., 10., 200.] {
                    let target =
                        complex.simplices()[source].value() + if dual { -offset } else { offset };
                    let row = state
                        .u_row(&complex, dual, i, dimension, target, &mut budget)
                        .unwrap();
                    for (j, &expected) in u[i].iter().enumerate() {
                        let source_j = if dual { n - 1 - j } else { j };
                        let within = if dual {
                            complex.simplices()[source_j].value() >= target
                        } else {
                            complex.simplices()[source_j].value() <= target
                        };
                        assert_eq!(row.get(&j), if within { expected } else { 0 });
                    }
                }
            }
        }
    }
}

#[test]
fn all_four_finite_directions_match_independent_block_swaps() {
    for seed in 0..10 {
        let complex = fixture(seed, seed % 2 == 0);
        let n = complex.len();
        let mut workspace = CriticalSetWorkspace::new(&complex).unwrap();
        let pairs: Vec<_> = workspace.finite_pairs().collect();
        for (birth, death) in pairs {
            for (endpoint, dual) in [(birth, true), (death, false)] {
                let focus = if dual { n - 1 - endpoint.0 } else { endpoint.0 };
                let dimension = complex.simplex(endpoint).unwrap().dimension();
                let block: Vec<_> = (0..n)
                    .filter(|&j| {
                        complex.simplices()[if dual { n - 1 - j } else { j }].dimension()
                            == dimension
                    })
                    .collect();
                let values: Vec<_> = (0..n)
                    .map(|j| {
                        let value = complex.simplices()[if dual { n - 1 - j } else { j }].value();
                        if dual { -value } else { value }
                    })
                    .collect();
                let d = boundary(&complex, dual);
                for offset in [-200., -7., -0.5, 0.5, 7., 200.] {
                    let target = values[focus] + offset;
                    let expected = swap_oracle(&d, focus, &values, target, &block);
                    let mut expected: Vec<_> = expected
                        .into_iter()
                        .map(|j| SimplexId(if dual { n - 1 - j } else { j }))
                        .collect();
                    expected.sort_unstable();
                    assert_eq!(
                        workspace
                            .critical_set(endpoint, if dual { -target } else { target })
                            .unwrap(),
                        expected,
                        "seed={seed}, dual={dual}, endpoint={endpoint:?}, offset={offset}"
                    );
                }
            }
        }
        let dual = workspace.dual.as_ref().unwrap();
        for (b, d) in workspace.finite_pairs() {
            assert_eq!(dual.reduction.deaths[n - 1 - d.0], Some(n - 1 - b.0));
        }
    }
}

#[test]
fn lazy_preparation_and_interrupted_caches_are_retryable() {
    let complex = fixture(3, true);
    let mut workspace = CriticalSetWorkspace::new(&complex).unwrap();
    let (birth, death) = workspace.finite_pairs().next().unwrap();
    let value = complex.simplex(birth).unwrap().value();
    assert!(workspace.critical_set(birth, value).unwrap().is_empty());
    assert!(workspace.dual.is_none());
    assert!(workspace.primal.rows.is_empty());
    let death_value = complex.simplex(death).unwrap().value();
    workspace.critical_set(death, death_value - 1.).unwrap();
    assert!(workspace.dual.is_none());
    assert!(workspace.primal.rows.is_empty());
    assert!(matches!(
        workspace.critical_set_with(death, death_value + 100., &Execution::default().max_work(5)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(workspace.primal.rows.is_empty());
    let expected = workspace.critical_set(death, death_value + 100.).unwrap();
    assert_eq!(
        workspace.primal.rows.keys().copied().collect::<Vec<_>>(),
        [complex.simplex(death).unwrap().dimension()]
    );
    assert_eq!(
        workspace.critical_set(death, death_value + 100.).unwrap(),
        expected
    );
    assert!(matches!(
        workspace.critical_set_with(birth, value - 100., &Execution::default().max_work(5)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(workspace.dual.is_none());
    workspace.critical_set(birth, value - 100.).unwrap();
    assert!(workspace.dual.is_some());
    let flag = AtomicBool::new(false);
    let execution = Execution::default().max_work(10_000).cancellation(&flag);
    let mut budget = WorkBudget::new(&execution).unwrap();
    budget.cancel_at_work(4);
    assert!(matches!(
        decompose(&complex, true, &mut budget),
        Err(Error::Cancelled)
    ));
    assert!(flag.load(Ordering::Relaxed));
    flag.store(false, Ordering::Relaxed);
    assert!(CriticalSetWorkspace::new_with(&complex, &execution).is_ok());
}

#[test]
fn max_displacement_mixed_directions_and_ties_keep_first() {
    let complex = SimplicialComplex::new(vec![Simplex::new(vec![0], 0.).unwrap()]).unwrap();
    let id = complex.find(&[0]).unwrap();
    let mut workspace = CriticalSetWorkspace::new(&complex).unwrap();
    assert_eq!(
        workspace.targets(&[(id, -4.), (id, 2.)]).unwrap(),
        [(id, -4.)]
    );
    assert_eq!(
        workspace.targets(&[(id, 2.), (id, -4.)]).unwrap(),
        [(id, -4.)]
    );
    assert_eq!(
        workspace.targets(&[(id, -2.), (id, 2.)]).unwrap(),
        [(id, -2.)]
    );
    assert_eq!(
        workspace.targets(&[(id, 2.), (id, -2.)]).unwrap(),
        [(id, 2.)]
    );
    assert_eq!(workspace.essential_births().collect::<Vec<_>>(), [id]);
    assert_eq!(workspace.critical_set(id, -1.).unwrap(), [id]);
    assert_eq!(workspace.critical_set(id, 1.).unwrap(), [id]);
    assert!(workspace.primal.rows.is_empty());
    assert!(workspace.dual.as_ref().unwrap().rows.is_empty());
    // Warm caches and measure a singleton through private shared-budget composition.
    let execution = Execution::default().max_work(u64::MAX);
    let mut budget = WorkBudget::new(&execution).unwrap();
    workspace.critical(id, 1., &mut budget).unwrap();
    let limit = budget.used() + 3;
    assert!(
        workspace
            .targets_with(&[(id, 1.)], &Execution::default().max_work(limit))
            .is_ok()
    );
    assert!(matches!(
        workspace.targets_with(&[(id, 1.), (id, 1.)], &Execution::default().max_work(limit)),
        Err(Error::WorkLimitExceeded { .. })
    ));
}

#[test]
fn empty_invalid_nonfinite_and_overflow_contracts() {
    let empty = SimplicialComplex::new(vec![]).unwrap();
    let mut workspace = CriticalSetWorkspace::new(&empty).unwrap();
    assert_eq!(workspace.finite_pairs().count(), 0);
    assert_eq!(workspace.essential_births().count(), 0);
    assert!(workspace.targets(&[]).unwrap().is_empty());
    assert!(matches!(
        workspace.critical_set(SimplexId(0), 1.),
        Err(Error::InvalidComplex { .. })
    ));
    let complex = SimplicialComplex::new(vec![Simplex::new(vec![0], -f64::MAX).unwrap()]).unwrap();
    let id = complex.find(&[0]).unwrap();
    let mut workspace = CriticalSetWorkspace::new(&complex).unwrap();
    for target in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            workspace.critical_set(id, target),
            Err(Error::NonFiniteValue { .. })
        ));
    }
    assert!(matches!(
        workspace.targets(&[(id, f64::MAX)]),
        Err(Error::NumericalFailure { .. })
    ));
    assert_eq!(workspace.critical_set(id, f64::MAX).unwrap(), [id]);
    assert!(matches!(
        CriticalSetWorkspace::new_with(&complex, &Execution::default().max_work(0)),
        Err(Error::WorkLimitExceeded { .. })
    ));
}
