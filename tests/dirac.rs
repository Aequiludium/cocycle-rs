//! Public real Dirac contracts, independent hand matrices and invalid chain inputs.
use cocycle::complex::{FilteredComplex, Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::spectral::{DiracOperator, MAX_DIRAC_DEGREES};
use cocycle::{Error, Result};
use std::sync::atomic::AtomicBool;

struct Cell {
    id: u32,
    degree: usize,
    value: f64,
    boundary: Vec<(u32, i32)>,
}

struct Chain(Vec<Cell>);

impl FilteredComplex for Chain {
    type CellId = u32;
    fn cells(&self) -> impl Iterator<Item = u32> + '_ {
        self.0.iter().map(|c| c.id)
    }
    fn dimension(&self, id: u32) -> usize {
        self.0.iter().find(|c| c.id == id).unwrap().degree
    }
    fn filtration_value(&self, id: u32) -> f64 {
        self.0.iter().find(|c| c.id == id).unwrap().value
    }
    fn boundary(&self, id: u32) -> impl Iterator<Item = (u32, i32)> + '_ {
        self.0
            .iter()
            .find(|c| c.id == id)
            .unwrap()
            .boundary
            .iter()
            .copied()
    }
}

fn chain(degrees: &[usize], boundaries: &[&[(u32, i32)]]) -> Chain {
    Chain(
        degrees
            .iter()
            .enumerate()
            .map(|(i, &degree)| Cell {
                id: (i as u32 + 1) * 10,
                degree,
                value: 0.0,
                boundary: boundaries[i].to_vec(),
            })
            .collect(),
    )
}

fn edge() -> Chain {
    chain(&[0, 0, 1], &[&[], &[], &[(10, -1), (20, 1)]])
}

fn triangle(filled: bool) -> Chain {
    let mut result = chain(
        &[0, 0, 0, 1, 1, 1],
        &[
            &[],
            &[],
            &[],
            &[(10, -1), (20, 1)],
            &[(10, -1), (30, 1)],
            &[(20, -1), (30, 1)],
        ],
    );
    if filled {
        result.0.push(Cell {
            id: 70,
            degree: 2,
            value: 0.0,
            boundary: vec![(40, 1), (50, -1), (60, 1)],
        });
    }
    result
}

fn multiply(a: &[f64], b: &[f64], n: usize) -> Vec<f64> {
    (0..n * n)
        .map(|i| (0..n).map(|j| a[i / n * n + j] * b[j * n + i % n]).sum())
        .collect()
}

#[test]
fn one_edge_matrix_action_and_analytic_eigenvectors() -> Result<()> {
    let operator = DiracOperator::new(&edge(), 0..=1)?;
    assert_eq!(operator.basis(), &[10, 20, 30]);
    assert_eq!(operator.degree_offsets(), &[0, 2, 3]);
    assert!(!operator.is_truncated());
    let k = operator.matrix(9)?;
    assert_eq!(k, [0., 0., -1., 0., 0., 1., -1., 1., 0.]);
    assert_eq!(operator.action(&[2., 5., 7.])?, [-7., 7., 3.]);
    assert_eq!(multiply(&k, &k, 3), [1., -1., 0., -1., 1., 0., 0., 0., 2.]);
    assert_eq!(operator.action(&[1., 1., 0.])?, [0., 0., 0.]);
    for lambda in [-2.0_f64.sqrt(), 2.0_f64.sqrt()] {
        let vector = [-1., 1., lambda];
        for (actual, component) in operator.action(&vector)?.iter().zip(vector) {
            assert!((actual - lambda * component).abs() < 1e-14);
        }
    }
    Ok(())
}

#[test]
fn triangle_square_has_independent_hodge_blocks_and_parity() -> Result<()> {
    for filled in [false, true] {
        let operator = DiracOperator::new(&triangle(filled), 0..=2)?;
        let n = operator.basis().len();
        let k = operator.matrix(n * n)?;
        let mut expected = vec![0.; n * n];
        let l0 = [2., -1., -1., -1., 2., -1., -1., -1., 2.];
        let l1 = if filled {
            [3., 0., 0., 0., 3., 0., 0., 0., 3.]
        } else {
            [2., 1., -1., 1., 2., 1., -1., 1., 2.]
        };
        for i in 0..3 {
            for j in 0..3 {
                expected[i * n + j] = l0[i * 3 + j];
                expected[(i + 3) * n + j + 3] = l1[i * 3 + j];
            }
        }
        if filled {
            expected[6 * n + 6] = 3.;
        }
        assert_eq!(multiply(&k, &k, n), expected);
        for i in 0..n {
            for j in 0..n {
                assert_eq!(k[i * n + j], k[j * n + i]);
            }
        }
        let x: Vec<_> = (0..n).map(|i| i as f64 - 2.).collect();
        let parity: Vec<_> = (0..n)
            .map(|i| if (3..6).contains(&i) { -1. } else { 1. })
            .collect();
        let graded: Vec<_> = x.iter().zip(&parity).map(|(a, p)| a * p).collect();
        let action = operator.action(&x)?;
        for ((actual, expected), p) in operator.action(&graded)?.iter().zip(action).zip(parity) {
            assert_eq!(*actual, -p * expected);
        }
        assert_eq!(
            operator.action(&operator.action(&x)?)?,
            (0..n)
                .map(|i| (0..n).map(|j| expected[i * n + j] * x[j]).sum())
                .collect::<Vec<f64>>()
        );
    }
    Ok(())
}

#[test]
fn truncation_removes_both_endpoint_maps() -> Result<()> {
    let edges = DiracOperator::new(&edge(), 1..=1)?;
    assert!(edges.is_truncated());
    assert_eq!(edges.degree_offsets(), &[0, 1]);
    assert_eq!(edges.matrix(1)?, [0.]);
    assert_eq!(edges.action(&[7.])?, [0.]); // Full source L1 is [2].
    let top = DiracOperator::new(&triangle(true), 1..=2)?;
    assert_eq!(top.degree_offsets(), &[0, 3, 4]);
    let k = top.matrix(16)?;
    assert_eq!(
        multiply(&k, &k, 4),
        [
            1., -1., 1., 0., -1., 1., -1., 0., 1., -1., 1., 0., 0., 0., 0., 3.,
        ]
    );
    Ok(())
}

#[test]
fn empty_degrees_and_owned_coordinates() -> Result<()> {
    let operator = {
        let source = chain(&[0, 3, 0], &[&[], &[], &[]]);
        DiracOperator::new(&source, 0..=4)?
    };
    assert_eq!(operator.basis(), &[10, 30, 20]);
    assert_eq!(operator.degree_offsets(), &[0, 2, 2, 2, 3, 3]);
    assert_eq!(operator.action(&[1., 2., 3.])?, [0.; 3]);
    for source in [chain(&[], &[]), edge()] {
        let empty = DiracOperator::new(&source, 4..=6)?;
        assert_eq!(empty.degree_offsets(), &[0, 0, 0, 0]);
        assert!(empty.basis().is_empty());
        assert_eq!(empty.matrix(0)?, []);
        assert_eq!(
            empty.action_with(&[], &Execution::default().max_work(0))?,
            []
        );
    }
    Ok(())
}

#[test]
fn degree_span_is_bounded_before_work_or_allocation() -> Result<()> {
    let source = chain(&[], &[]);
    for first in [0, usize::MAX - (MAX_DIRAC_DEGREES - 1)] {
        let degrees = first..=first + (MAX_DIRAC_DEGREES - 1);
        let operator = DiracOperator::new(&source, degrees.clone())?;
        assert_eq!(operator.degrees(), degrees);
        assert_eq!(operator.degree_offsets(), vec![0; MAX_DIRAC_DEGREES + 1]);
        assert!(operator.basis().is_empty());
    }
    assert!(matches!(
        DiracOperator::new(&source, 0..=MAX_DIRAC_DEGREES),
        Err(Error::InvalidParameter {
            parameter: "degrees",
            ..
        })
    ));
    for execution in [Execution::default(), Execution::default().max_work(0)] {
        assert!(matches!(
            DiracOperator::new_with(&source, 0..=MAX_DIRAC_DEGREES, &execution),
            Err(Error::InvalidParameter {
                parameter: "degrees",
                ..
            })
        ));
    }
    // A zero work budget keeps this huge-input regression safe if the guard
    // regresses: it must reject the parameter, not merely exhaust work.
    assert!(matches!(
        DiracOperator::new_with(
            &source,
            0..=usize::MAX - 2,
            &Execution::default().max_work(0),
        ),
        Err(Error::InvalidParameter {
            parameter: "degrees",
            ..
        })
    ));
    Ok(())
}

#[test]
fn duplicate_signed_terms_and_reorientation() -> Result<()> {
    let mut source = edge();
    source.0[2].boundary = vec![(20, 3), (10, -2), (20, -2), (10, 1), (999, 0)];
    let operator = DiracOperator::new(&source, 0..=1)?;
    assert_eq!(operator.entries().count(), 4);
    assert_eq!(
        operator.matrix(9)?,
        DiracOperator::new(&edge(), 0..=1)?.matrix(9)?
    );
    source.0[2].boundary = vec![(20, -1), (10, 1)];
    let reversed = DiracOperator::new(&source, 0..=1)?;
    for (a, b) in operator.matrix(9)?.iter().zip(reversed.matrix(9)?) {
        assert_eq!(*a, -b);
    }
    Ok(())
}

#[test]
fn rejects_integer_chain_failure_even_when_f2_square_is_zero() {
    let source = chain(&[0, 1, 2], &[&[], &[(10, 1)], &[(20, 2)]]);
    assert!(matches!(
        DiracOperator::new(&source, 0..=2),
        Err(Error::InvalidComplex { .. })
    ));
    // Removing C0 makes this a valid truncated interval: only B2 remains.
    assert!(DiracOperator::new(&source, 1..=2).is_ok());
}

#[test]
fn rejects_bad_source_shapes_order_and_values() {
    let invalid = [
        chain(&[0, 1], &[&[], &[(999, 1)]]),
        chain(&[0, 2], &[&[], &[(10, 1)]]),
        chain(&[1, 0], &[&[(20, 1)], &[]]),
        chain(&[0, 0], &[&[], &[(10, 1)]]),
    ];
    for source in invalid {
        assert!(matches!(
            DiracOperator::new(&source, 0..=2),
            Err(Error::InvalidComplex { .. })
        ));
    }
    let mut source = edge();
    source.0[1].id = 10;
    assert!(DiracOperator::new(&source, 0..=1).is_err());
    for value in [f64::NAN, f64::INFINITY, -1.] {
        let mut source = edge();
        source.0[1].value = value;
        assert!(DiracOperator::new(&source, 0..=1).is_err());
    }
    assert!(DiracOperator::new(&edge(), std::ops::RangeInclusive::new(2, 1)).is_err());
    let mut exhausted = 0..=0;
    exhausted.next();
    assert!(DiracOperator::new(&edge(), exhausted).is_err());
    assert!(matches!(
        DiracOperator::new(&edge(), 0..=usize::MAX),
        Err(Error::SizeOverflow { .. })
    ));
}

#[test]
fn controls_dense_cap_nonfinite_and_overflow_errors_are_recoverable() -> Result<()> {
    let operator = DiracOperator::new(&edge(), 0..=1)?;
    assert!(matches!(
        operator.matrix(8),
        Err(Error::InvalidParameter {
            parameter: "max_entries",
            ..
        })
    ));
    assert!(matches!(
        operator.action(&[1.]),
        Err(Error::ShapeMismatch { .. })
    ));
    for bad in [f64::NAN, f64::INFINITY] {
        assert!(matches!(
            operator.action(&[bad, 0., 0.]),
            Err(Error::NonFiniteValue { index: Some(0), .. })
        ));
    }
    assert!(matches!(
        operator.action(&[-f64::MAX, f64::MAX, 0.]),
        Err(Error::NumericalFailure { .. })
    ));
    let limited = Execution::default().max_work(1);
    assert!(matches!(
        DiracOperator::new_with(&edge(), 0..=1, &limited),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(matches!(
        operator.action_with(&[0.; 3], &limited),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(matches!(
        operator.matrix_with(9, &limited),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let flag = AtomicBool::new(true);
    let cancelled = Execution::default().cancellation(&flag);
    assert!(matches!(
        DiracOperator::new_with(&edge(), 0..=1, &cancelled),
        Err(Error::Cancelled)
    ));
    assert_eq!(
        operator.action_with(&[0.; 3], &cancelled),
        Err(Error::Cancelled)
    );
    assert_eq!(operator.matrix_with(9, &cancelled), Err(Error::Cancelled));
    assert_eq!(operator.action(&[2., 5., 7.])?, [-7., 7., 3.]);
    Ok(())
}

#[test]
fn supplied_stage_uses_original_owner_ids_and_no_ph() -> Result<()> {
    let complex = SimplicialComplex::new(vec![
        Simplex::new(vec![9], -2.)?,
        Simplex::new(vec![42], -1.)?,
        Simplex::new(vec![9, 42], 3.)?,
    ])?;
    let early = complex.stage(0.)?;
    let late = complex.stage(3.)?;
    let isolated = DiracOperator::new(&early, 0..=1)?;
    assert_eq!(isolated.degree_offsets(), &[0, 2, 2]);
    assert_eq!(isolated.matrix(4)?, [0.; 4]);
    let connected = DiracOperator::new(&late, 0..=1)?;
    for &id in connected.basis() {
        assert!(late.simplex(id).is_some());
    }
    assert_eq!(connected.matrix(9)?, [0., 0., -1., 0., 0., 1., -1., 1., 0.]);
    Ok(())
}

#[test]
fn sparse_action_does_not_require_dense_export() -> Result<()> {
    struct Path(usize);
    impl FilteredComplex for Path {
        type CellId = usize;
        fn cells(&self) -> impl Iterator<Item = usize> + '_ {
            0..2 * self.0 - 1
        }
        fn dimension(&self, cell: usize) -> usize {
            usize::from(cell >= self.0)
        }
        fn filtration_value(&self, _: usize) -> f64 {
            0.0
        }
        fn boundary(&self, cell: usize) -> impl Iterator<Item = (usize, i32)> + '_ {
            let vertex = cell.saturating_sub(self.0);
            [(vertex, -1), (vertex + 1, 1)]
                .into_iter()
                .take(if cell < self.0 { 0 } else { 2 })
        }
    }
    let n = 20_000;
    let operator = DiracOperator::new(&Path(n), 0..=1)?;
    assert_eq!(operator.degree_offsets(), &[0, n, 2 * n - 1]);
    assert_eq!(operator.entries().count(), 4 * (n - 1));
    assert!(operator.matrix(1_000_000).is_err());
    let result = operator.action(&vec![1.; 2 * n - 1])?;
    assert_eq!(result[0], -1.);
    assert_eq!(result[n - 1], 1.);
    assert!(
        result[1..n - 1]
            .iter()
            .chain(&result[n..])
            .all(|&v| v == 0.)
    );
    Ok(())
}

#[test]
fn all_small_work_limits_are_retryable() -> Result<()> {
    let source = triangle(true);
    let expected = DiracOperator::new(&source, 0..=2)?;
    let matrix = expected.matrix(49)?;
    let vector = [1., 2., 3., 4., 5., 6., 7.];
    let action = expected.action(&vector)?;
    for limit in 0..128 {
        let controls = Execution::default().max_work(limit);
        match DiracOperator::new_with(&source, 0..=2, &controls) {
            Ok(operator) => assert_eq!(operator.matrix(49)?, matrix),
            Err(Error::WorkLimitExceeded { .. }) => {}
            Err(error) => panic!("unexpected construction error: {error}"),
        }
        for result in [
            expected.action_with(&vector, &controls),
            expected.matrix_with(49, &controls),
        ] {
            match result {
                Ok(values) => assert!(values == action || values == matrix),
                Err(Error::WorkLimitExceeded { .. }) => {}
                Err(error) => panic!("unexpected query error: {error}"),
            }
        }
        assert_eq!(expected.action(&vector)?, action);
    }
    Ok(())
}
