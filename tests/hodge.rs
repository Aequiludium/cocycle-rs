//! Independent finite expectations and public real Hodge contracts.
use cocycle::Error;
use cocycle::complex::{FilteredComplex, Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::hodge::{HodgeLaplacian, HodgeSelection, SpectrumOptions};
use std::sync::atomic::AtomicBool;

type Cell = (usize, usize, Vec<(usize, i32)>);
struct Cells(Vec<Cell>);
impl FilteredComplex for Cells {
    type CellId = usize;
    fn cells(&self) -> impl Iterator<Item = usize> + '_ {
        self.0.iter().map(|c| c.0)
    }
    fn dimension(&self, id: usize) -> usize {
        self.0.iter().find(|c| c.0 == id).unwrap().1
    }
    fn filtration_value(&self, _: usize) -> f64 {
        panic!("Hodge must not read filtration weights")
    }
    fn boundary(&self, id: usize) -> impl Iterator<Item = (usize, i32)> + '_ {
        self.0.iter().find(|c| c.0 == id).unwrap().2.iter().copied()
    }
}
fn triangle(filled: bool) -> Cells {
    let mut cells = vec![
        (0, 0, vec![]),
        (1, 0, vec![]),
        (2, 0, vec![]),
        (3, 1, vec![(0, -1), (1, 1)]),
        (4, 1, vec![(0, -1), (2, 1)]),
        (5, 1, vec![(1, -1), (2, 1)]),
    ];
    if filled {
        cells.push((6, 2, vec![(3, 1), (4, -1), (5, 1)]));
    }
    Cells(cells)
}
fn close(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (&x, &y) in a.iter().zip(b) {
        assert!((x - y).abs() < 1e-10, "{a:?} != {b:?}");
    }
}
fn op(cells: &Cells, degree: usize, selection: HodgeSelection) -> HodgeLaplacian<usize> {
    HodgeLaplacian::new(cells, degree, selection).unwrap()
}

#[test]
fn edge_ordinary_h0_and_top_degree() {
    let edge = Cells(vec![
        (7, 0, vec![]),
        (19, 0, vec![]),
        (103, 1, vec![(7, -1), (19, 1)]),
    ]);
    let l = op(&edge, 0, HodgeSelection::Full);
    assert_eq!(l.basis(), &[7, 19]);
    assert!(l.previous_basis().is_empty());
    assert_eq!(l.matrix(&Execution::default()).unwrap(), [1., -1., -1., 1.]);
    assert_eq!(
        l.action(&[2., 5.], &Execution::default()).unwrap(),
        [-3., 3.]
    );
    let s = l
        .spectrum(0..2, &SpectrumOptions::default(), &Execution::default())
        .unwrap();
    close(s.values(), &[0., 2.]);
    assert_eq!(s.numerical_nullity(1e-10).unwrap(), 1);
    assert!(s.residuals().iter().all(|&r| r < 1e-10));
    assert!(s.orthogonality_error() < 1e-10);
    assert_eq!(
        op(&edge, 1, HodgeSelection::Full)
            .matrix(&Execution::default())
            .unwrap(),
        [2.]
    );
    assert_eq!(
        op(&edge, 0, HodgeSelection::Lower)
            .matrix(&Execution::default())
            .unwrap(),
        [0.; 4]
    );
}

#[test]
fn triangle_harmonic_space_filling_and_summands() {
    let l = op(&triangle(false), 1, HodgeSelection::Full);
    assert_eq!(
        l.matrix(&Execution::default()).unwrap(),
        [2., 1., -1., 1., 2., 1., -1., 1., 2.]
    );
    assert_eq!(
        l.action(&[1., -1., 1.], &Execution::default()).unwrap(),
        [0.; 3]
    );
    let s = l
        .spectrum(0..3, &SpectrumOptions::default(), &Execution::default())
        .unwrap();
    close(s.values(), &[0., 3., 3.]);
    assert_eq!(s.numerical_nullity(1e-10).unwrap(), 1);
    assert!(s.residuals().iter().all(|&r| r < 1e-10));
    assert!(s.orthogonality_error() < 1e-10);
    // Compare the harmonic projector, independent of vector sign.
    let z = [1., -1., 1.];
    for i in 0..3 {
        for j in 0..3 {
            assert!((s.vectors()[i * 3] * s.vectors()[j * 3] - z[i] * z[j] / 3.).abs() < 1e-10);
        }
    }
    let full = op(&triangle(true), 1, HodgeSelection::Full);
    let upper = op(&triangle(true), 1, HodgeSelection::Upper);
    assert_eq!(
        full.matrix(&Execution::default()).unwrap(),
        [3., 0., 0., 0., 3., 0., 0., 0., 3.]
    );
    assert_eq!(
        upper.matrix(&Execution::default()).unwrap(),
        [1., -1., 1., -1., 1., -1., 1., -1., 1.]
    );
    close(
        full.spectrum(0..3, &SpectrumOptions::default(), &Execution::default())
            .unwrap()
            .values(),
        &[3.; 3],
    );
    let x = [2., -4., 7.];
    let a = l.action(&x, &Execution::default()).unwrap();
    let b = upper.action(&x, &Execution::default()).unwrap();
    close(
        &full.action(&x, &Execution::default()).unwrap(),
        &a.iter().zip(b).map(|(a, b)| a + b).collect::<Vec<_>>(),
    );
    let matrix = full.matrix(&Execution::default()).unwrap();
    for row in 0..3 {
        assert_eq!(
            matrix[row * 3..row * 3 + 3]
                .iter()
                .zip(x)
                .map(|(a, b)| a * b)
                .sum::<f64>(),
            3. * x[row]
        );
    }
}

#[test]
fn stage_and_tetrahedron_sphere_top_homology() {
    let simplices = (1_usize..16)
        .map(|mask| {
            Simplex::new(
                (0..4).filter(|v| mask & (1 << v) != 0).collect(),
                if mask == 15 { 1. } else { 0. },
            )
            .unwrap()
        })
        .collect();
    let source = SimplicialComplex::new(simplices).unwrap();
    let stage = source.stage(0.).unwrap();
    let l = HodgeLaplacian::new(&stage, 2, HodgeSelection::Full).unwrap();
    let s = l
        .spectrum(0..4, &SpectrumOptions::default(), &Execution::default())
        .unwrap();
    close(s.values(), &[0., 4., 4., 4.]);
    assert_eq!(s.numerical_nullity(1e-10).unwrap(), 1);
    let filled = HodgeLaplacian::new(&source, 2, HodgeSelection::Full).unwrap();
    close(
        filled
            .spectrum(0..4, &SpectrumOptions::default(), &Execution::default())
            .unwrap()
            .values(),
        &[4.; 4],
    );
    let edges = HodgeLaplacian::new(&source, 1, HodgeSelection::Full).unwrap();
    close(
        edges
            .spectrum(0..6, &SpectrumOptions::default(), &Execution::default())
            .unwrap()
            .values(),
        &[4.; 6],
    );
}

#[test]
fn orientation_basis_permutation_and_duplicate_incidence() {
    let original = triangle(true);
    let mut changed = triangle(true);
    // Reverse edge 01, consistently in both adjacent maps.
    for term in &mut changed.0[3].2 {
        term.1 = -term.1;
    }
    changed.0[6].2[0].1 = -changed.0[6].2[0].1;
    changed.0.swap(3, 5);
    changed.0[4].2.extend([(0, 17), (0, -17)]);
    let l = op(&changed, 1, HodgeSelection::Full);
    assert_eq!(l.basis(), &[5, 4, 3]);
    assert_eq!(
        l.matrix(&Execution::default()).unwrap(),
        op(&original, 1, HodgeSelection::Full)
            .matrix(&Execution::default())
            .unwrap()
    );
    let changed = op(&Cells(changed.0[..6].to_vec()), 1, HodgeSelection::Full);
    let original = op(&triangle(false), 1, HodgeSelection::Full);
    let x = [2., -3., 7.];
    let y = original.action(&x, &Execution::default()).unwrap();
    close(
        &changed
            .action(&[x[2], x[1], -x[0]], &Execution::default())
            .unwrap(),
        &[y[2], y[1], -y[0]],
    );
}

#[test]
fn empty_is_explicit_and_partial_spectrum_cannot_certify_nullity() {
    let l = op(&triangle(false), 9, HodgeSelection::Full);
    assert!(l.basis().is_empty());
    assert!(l.matrix(&Execution::default()).unwrap().is_empty());
    assert!(l.action(&[], &Execution::default()).unwrap().is_empty());
    let s = l
        .spectrum(0..0, &SpectrumOptions::default(), &Execution::default())
        .unwrap();
    assert!(s.is_complete());
    assert_eq!(s.numerical_nullity(0.).unwrap(), 0);
    let l = op(&triangle(false), 1, HodgeSelection::Full);
    for range in [0..0, 1..3] {
        let s = l
            .spectrum(
                range.clone(),
                &SpectrumOptions::default(),
                &Execution::default(),
            )
            .unwrap();
        assert_eq!(s.index_range(), range);
        assert!(!s.is_complete());
        assert!(s.numerical_nullity(1e-10).is_err());
    }
    let zero = op(&Cells(vec![(3, 0, vec![])]), 0, HodgeSelection::Full);
    let s = zero
        .spectrum(0..1, &SpectrumOptions::default(), &Execution::default())
        .unwrap();
    assert_eq!(s.relative_residuals(), &[None]);
    assert_eq!(s.residuals(), &[0.]);
}

#[test]
fn rejects_invalid_windows_in_integer_domain_and_invalid_vectors() {
    let cases = [
        Cells(vec![(0, 0, vec![]), (0, 0, vec![])]),
        Cells(vec![(0, 0, vec![]), (1, 1, vec![(99, 1)])]),
        Cells(vec![(0, 0, vec![]), (1, 2, vec![(0, 1)])]),
        Cells(vec![(1, 1, vec![(0, 1)]), (0, 0, vec![])]),
        Cells(vec![(0, 0, vec![(0, 1)])]),
        // This product is zero over F2 but nonzero over R/Z.
        Cells(vec![
            (0, 0, vec![]),
            (1, 1, vec![(0, 2)]),
            (2, 2, vec![(1, 1)]),
        ]),
    ];
    for (i, cells) in cases.iter().enumerate() {
        let degree = if i == 4 { 0 } else { 1 };
        assert!(
            matches!(
                HodgeLaplacian::new(cells, degree, HodgeSelection::Full),
                Err(Error::InvalidComplex { .. })
            ),
            "case {i}"
        );
    }
    let l = op(&triangle(false), 1, HodgeSelection::Full);
    assert!(matches!(
        l.action(&[], &Execution::default()),
        Err(Error::ShapeMismatch { .. })
    ));
    assert!(matches!(
        l.action(&[0., f64::NAN, 0.], &Execution::default()),
        Err(Error::NonFiniteValue { .. })
    ));
    assert!(matches!(
        l.action(&[f64::MAX; 3], &Execution::default()),
        Err(Error::NumericalFailure { .. })
    ));
    assert!(matches!(
        HodgeLaplacian::new(&triangle(false), usize::MAX, HodgeSelection::Full),
        Err(Error::SizeOverflow { .. })
    ));
    // A torsion cell attachment has different real and F2 nullities.
    let torsion = Cells(vec![(0, 0, vec![]), (1, 1, vec![]), (2, 2, vec![(1, 2)])]);
    assert_eq!(
        op(&torsion, 1, HodgeSelection::Full)
            .matrix(&Execution::default())
            .unwrap(),
        [4.]
    );
}

#[test]
fn policy_failure_budget_and_cancellation_are_errors() {
    let cells = triangle(false);
    let l = op(&cells, 1, HodgeSelection::Full);
    let execution = Execution::default().max_work(0);
    assert!(matches!(
        HodgeLaplacian::new_with(&cells, 1, HodgeSelection::Full, &execution),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(matches!(
        l.action(&[1.; 3], &execution),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(matches!(
        l.matrix(&execution),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(matches!(
        l.spectrum(0..3, &SpectrumOptions::default(), &execution),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let cancelled = AtomicBool::new(true);
    let execution = Execution::default().cancellation(&cancelled);
    assert!(matches!(
        HodgeLaplacian::new_with(&cells, 1, HodgeSelection::Full, &execution),
        Err(Error::Cancelled)
    ));
    assert_eq!(l.action(&[1.; 3], &execution), Err(Error::Cancelled));
    assert_eq!(l.matrix(&execution), Err(Error::Cancelled));
    assert!(matches!(
        l.spectrum(0..3, &SpectrumOptions::default(), &execution),
        Err(Error::Cancelled)
    ));
    for options in [
        SpectrumOptions {
            max_dimension: 2,
            ..Default::default()
        },
        SpectrumOptions {
            absolute_tolerance: f64::NAN,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            l.spectrum(0..3, &options, &Execution::default()),
            Err(Error::InvalidParameter { .. })
        ));
    }
    assert!(matches!(
        l.spectrum(
            0..3,
            &SpectrumOptions {
                max_sweeps: 0,
                ..Default::default()
            },
            &Execution::default()
        ),
        Err(Error::NumericalFailure { .. })
    ));
    assert!(
        l.spectrum(2..4, &SpectrumOptions::default(), &Execution::default())
            .is_err()
    );
    // Failure leaves the prepared operator usable.
    close(
        l.spectrum(0..3, &SpectrumOptions::default(), &Execution::default())
            .unwrap()
            .values(),
        &[0., 3., 3.],
    );
}

#[test]
fn preparation_and_terminal_have_separate_counters() {
    let cells = Cells(vec![
        (7, 0, vec![]),
        (19, 0, vec![]),
        (103, 1, vec![(7, -1), (19, 1)]),
    ]);
    let controls = Execution::default().max_work(10);
    // Preparation consumes seven units and action ten; sharing the immutable
    // policy does not create a seventeen-unit cumulative operation.
    let l = HodgeLaplacian::new_with(&cells, 0, HodgeSelection::Full, &controls).unwrap();
    assert_eq!(l.action(&[2., 5.], &controls).unwrap(), [-3., 3.]);
    assert!(matches!(
        HodgeLaplacian::new_with(
            &cells,
            0,
            HodgeSelection::Full,
            &Execution::default().max_work(6)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
}

#[test]
fn rp2_real_harmonics_differ_from_f2_persistence() {
    use cocycle::algebra::PrimeField;
    use cocycle::descriptors::betti_curve;
    use cocycle::persistence::PersistenceBuilder;
    use std::collections::BTreeSet;
    let facets = [
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
    let mut faces = BTreeSet::new();
    for facet in facets {
        for mask in 1..8 {
            faces.insert(
                facet
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| mask & (1 << i) != 0)
                    .map(|(_, &v)| v)
                    .collect::<Vec<_>>(),
            );
        }
    }
    let source = SimplicialComplex::new(
        faces
            .into_iter()
            .map(|v| Simplex::new(v, 0.).unwrap())
            .collect(),
    )
    .unwrap();
    let l = HodgeLaplacian::new(&source, 1, HodgeSelection::Full).unwrap();
    assert_eq!(l.basis().len(), 15);
    let s = l
        .spectrum(0..15, &SpectrumOptions::default(), &Execution::default())
        .unwrap();
    assert_eq!(s.numerical_nullity(1e-10).unwrap(), 0);
    assert!(s.numerical_kernel_columns(1e-10).unwrap().next().is_none());
    assert!(s.values()[0] > 0.7);
    for (prime, expected) in [(2, 1), (3, 0)] {
        let result = PersistenceBuilder::from_complex(&source)
            .max_homology_dimension(2)
            .field(PrimeField::new(prime).unwrap())
            .compute()
            .unwrap();
        assert_eq!(betti_curve(result.diagram(), 1, &[0.]).unwrap(), [expected]);
    }
}

#[test]
fn harmonic_columns_keep_vector_value_and_residual_pairing() {
    let s = op(&triangle(false), 1, HodgeSelection::Full)
        .spectrum(0..3, &SpectrumOptions::default(), &Execution::default())
        .unwrap();
    let columns: Vec<_> = s.numerical_kernel_columns(1e-10).unwrap().collect();
    assert_eq!(columns, [0]);
    for column in columns {
        let v: Vec<_> = (0..3).map(|row| s.vectors()[row * 3 + column]).collect();
        assert!(s.values()[column].abs() <= 1e-10);
        assert!(s.residuals()[column] < 1e-10);
        assert!(
            op(&triangle(false), 1, HodgeSelection::Full)
                .action(&v, &Execution::default())
                .unwrap()
                .iter()
                .all(|x| x.abs() < 1e-10)
        );
    }
    let partial = op(&triangle(false), 1, HodgeSelection::Full)
        .spectrum(1..3, &SpectrumOptions::default(), &Execution::default())
        .unwrap();
    assert!(partial.numerical_kernel_columns(1e-10).is_err());
    assert!(s.numerical_kernel_columns(f64::NAN).is_err());
}

#[test]
fn default_maximum_dimension_path_converges_and_next_size_is_rejected() {
    for n in [128, 129] {
        let mut cells: Vec<Cell> = (0..n).map(|i| (i, 0, vec![])).collect();
        cells.extend((0..n - 1).map(|i| (n + i, 1, vec![(i, -1), (i + 1, 1)])));
        let l = op(&Cells(cells), 0, HodgeSelection::Full);
        let result = l.spectrum(0..n, &SpectrumOptions::default(), &Execution::default());
        if n == 129 {
            assert!(matches!(result, Err(Error::InvalidParameter { .. })));
            continue;
        }
        let s = result.unwrap();
        for (j, &lambda) in s.values().iter().enumerate() {
            let expected = 2. - 2. * (std::f64::consts::PI * j as f64 / n as f64).cos();
            assert!((lambda - expected).abs() < 1e-9);
        }
        assert_eq!(s.numerical_nullity(1e-10).unwrap(), 1);
        assert!(s.residuals().iter().all(|&r| r < 1e-9));
        assert!(s.orthogonality_error() < 1e-10);
        println!(
            "Hodge n=128 path: sweeps={}, max_residual={}, orthogonality={}",
            s.sweeps(),
            s.residuals().iter().copied().fold(0., f64::max),
            s.orthogonality_error()
        );
    }
}
