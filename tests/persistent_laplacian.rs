//! Public contracts and hand-derived persistent Laplacian examples.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::filtration::RipsBuilder;
use cocycle::geometry::PointCloudView;
use cocycle::spectral::{PersistentLaplacian, SpectralOptions, fixed_endpoint_laplacians};
use cocycle::{Error, Result};
use std::sync::atomic::AtomicBool;

fn complex(cells: &[(&[usize], f64)]) -> SimplicialComplex {
    SimplicialComplex::new(
        cells
            .iter()
            .map(|(v, t)| Simplex::new(v.to_vec(), *t).unwrap())
            .collect(),
    )
    .unwrap()
}
fn close(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() < 1e-9, "{a:?} != {b:?}");
    }
}
fn pair(source: &SimplicialComplex, s: f64, t: f64, q: usize) -> PersistentLaplacian {
    PersistentLaplacian::compute(
        &source.stage(s).unwrap(),
        &source.stage(t).unwrap(),
        q,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap()
}

#[test]
fn path_kron_and_owned_action() {
    let result = {
        let source = complex(&[
            (&[0], 0.),
            (&[2], 0.),
            (&[1], 1.),
            (&[0, 1], 2.),
            (&[1, 2], 2.),
        ]);
        pair(&source, 0., 2., 0)
    };
    close(result.matrix(), &[0.5, -0.5, -0.5, 0.5]);
    close(
        &result.spectrum(&Execution::default()).unwrap().eigenvalues,
        &[0., 1.],
    );
    close(
        &result.action(&[1., -1.], &Execution::default()).unwrap(),
        &[1., -1.],
    );
    assert_eq!(result.scales(), (0., 2.));
    assert_eq!(result.basis()[0].vertices(), &[2]);
    assert_eq!(result.target_basis().len(), 3);
    assert_eq!(result.degree(), 0);
    assert!(result.source_context().is_none());
    assert_eq!(
        result.source_coverage(),
        cocycle::filtration::Coverage::Complete
    );
    assert_eq!(result.construction_dimension(), None);
}

#[test]
fn ring_to_filled_triangle_and_hodge_degeneracy() {
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 0.),
        (&[0, 1], 1.),
        (&[0, 2], 1.),
        (&[1, 2], 1.),
        (&[0, 1, 2], 2.),
    ]);
    let ring = pair(&source, 1., 1., 1);
    close(
        &ring.spectrum(&Execution::default()).unwrap().eigenvalues,
        &[0., 3., 3.],
    );
    let filled = pair(&source, 1., 2., 1);
    close(filled.matrix(), &[3., 0., 0., 0., 3., 0., 0., 0., 3.]);
    close(
        &filled.spectrum(&Execution::default()).unwrap().eigenvalues,
        &[3., 3., 3.],
    );
    close(pair(&source, 2., 2., 1).matrix(), filled.matrix());
    close(pair(&source, 2., 2., 2).matrix(), &[3.]);
    assert!(pair(&source, 2., 2., 3).matrix().is_empty());
}

#[test]
fn singular_removed_component_and_entirely_zero_block() {
    let source = complex(&[(&[0], 0.), (&[1], 1.), (&[2], 1.), (&[1, 2], 2.)]);
    close(pair(&source, 0., 2., 0).matrix(), &[0.]);
    let isolated = complex(&[(&[0], 0.), (&[1], 1.)]);
    close(pair(&isolated, 0., 1., 0).matrix(), &[0.]);
}

#[test]
fn fixed_endpoint_uses_only_upper_term_with_ties_and_empty_stages() {
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 1.),
        (&[3], 1.),
        (&[0, 1], 2.),
        (&[0, 2], 2.),
        (&[1, 2], 2.),
        (&[2, 3], 2.),
        (&[0, 1, 2], 3.),
    ]);
    let stages = [-1., 0., 1., 2., 2., 3.].map(|s| source.stage(s).unwrap());
    for q in 0..=3 {
        let results = fixed_endpoint_laplacians(
            &stages,
            &stages[5],
            q,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        for (k, result) in stages.iter().zip(results) {
            close(result.matrix(), pair(&source, k.scale(), 3., q).matrix());
            assert_eq!(result.scales(), (k.scale(), 3.));
            let eig = result.spectrum(&Execution::default()).unwrap();
            assert!(eig.eigenvalues.iter().all(|v| *v >= -1e-9));
        }
    }
    assert!(
        fixed_endpoint_laplacians(
            &[],
            &stages[5],
            0,
            SpectralOptions::default(),
            &Execution::default()
        )
        .unwrap()
        .is_empty()
    );
}

#[test]
fn basis_permutation_and_orientation_under_relabeling() {
    let a = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 0.),
        (&[0, 1], 1.),
        (&[0, 2], 1.),
        (&[1, 2], 2.),
        (&[0, 1, 2], 3.),
    ]);
    // 0->20, 1->7, 2->40 reverses the orientation of edge 01.
    let b = complex(&[
        (&[20], 0.),
        (&[7], 0.),
        (&[40], 0.),
        (&[7, 20], 1.),
        (&[20, 40], 1.),
        (&[7, 40], 2.),
        (&[7, 20, 40], 3.),
    ]);
    let x = pair(&a, 1., 3., 1);
    let y = pair(&b, 1., 3., 1);
    close(
        &x.spectrum(&Execution::default()).unwrap().eigenvalues,
        &y.spectrum(&Execution::default()).unwrap().eigenvalues,
    );
    let mapping = [20, 7, 40];
    for (i, s) in x.basis().iter().enumerate() {
        let raw: Vec<_> = s.vertices().iter().map(|v| mapping[*v]).collect();
        let sign_i = if raw[0] < raw[1] { 1. } else { -1. };
        let mut mapped = raw;
        mapped.sort();
        let ii = y
            .basis()
            .iter()
            .position(|v| v.vertices() == mapped)
            .unwrap();
        for (j, t) in x.basis().iter().enumerate() {
            let raw: Vec<_> = t.vertices().iter().map(|v| mapping[*v]).collect();
            let sign_j = if raw[0] < raw[1] { 1. } else { -1. };
            let mut mapped = raw;
            mapped.sort();
            let jj = y
                .basis()
                .iter()
                .position(|v| v.vertices() == mapped)
                .unwrap();
            close(
                &[x.matrix()[i * 2 + j]],
                &[sign_i * sign_j * y.matrix()[ii * 2 + jj]],
            );
        }
    }
}

#[test]
fn invalid_inclusion_policy_and_coordinates() {
    let a = complex(&[(&[0], 0.), (&[1], 1.)]);
    let b = a.clone();
    let k = a.stage(0.).unwrap();
    let l = a.stage(1.).unwrap();
    fn compute(
        k: &cocycle::filtration::SimplicialStage<'_>,
        l: &cocycle::filtration::SimplicialStage<'_>,
        o: SpectralOptions,
        e: &Execution<'_>,
    ) -> Result<PersistentLaplacian> {
        PersistentLaplacian::compute(k, l, 0, o, e)
    }
    assert!(
        compute(
            &k,
            &b.stage(1.).unwrap(),
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    assert!(compute(&l, &k, SpectralOptions::default(), &Execution::default()).is_err());
    assert!(
        fixed_endpoint_laplacians(
            &[l, k],
            &l,
            0,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    for tolerance in [0., -1., f64::NAN, f64::INFINITY, 1.] {
        assert!(
            compute(
                &k,
                &l,
                SpectralOptions {
                    relative_tolerance: tolerance,
                    ..Default::default()
                },
                &Execution::default()
            )
            .is_err()
        );
    }
    let r = pair(&a, 0., 1., 0);
    assert!(matches!(
        r.action(&[], &Execution::default()),
        Err(Error::ShapeMismatch { .. })
    ));
    assert!(matches!(
        r.action(&[f64::NAN], &Execution::default()),
        Err(Error::NonFiniteValue { .. })
    ));
    let source = complex(&[(&[0], 0.), (&[1], 0.), (&[0, 1], 1.)]);
    assert!(matches!(
        pair(&source, 0., 1., 0).action(&[f64::MAX, -f64::MAX], &Execution::default()),
        Err(Error::NumericalFailure { .. })
    ));
}

#[test]
fn cancellation_work_and_storage_limits() {
    let a = complex(&[(&[0], 0.), (&[1], 1.), (&[0, 1], 2.)]);
    let k = a.stage(0.).unwrap();
    let l = a.stage(2.).unwrap();
    let flag = AtomicBool::new(true);
    for execution in [
        Execution::default().cancellation(&flag),
        Execution::default().max_work(0),
        Execution::default().max_work(8),
    ] {
        assert!(
            PersistentLaplacian::compute(&k, &l, 0, SpectralOptions::default(), &execution)
                .is_err()
        );
        assert!(
            fixed_endpoint_laplacians(&[k, l], &l, 0, SpectralOptions::default(), &execution)
                .is_err()
        );
    }
    let small = SpectralOptions {
        max_dense_entries: 3,
        ..Default::default()
    };
    assert!(PersistentLaplacian::compute(&k, &l, 0, small, &Execution::default()).is_err());
    let o = SpectralOptions {
        max_dense_entries: 4,
        ..Default::default()
    };
    assert!(fixed_endpoint_laplacians(&[l, l], &l, 0, o, &Execution::default()).is_err());
    let r = pair(&a, 0., 2., 0);
    assert!(matches!(
        r.spectrum(&Execution::default().cancellation(&flag)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        r.action(&[1.], &Execution::default().max_work(0)),
        Err(Error::WorkLimitExceeded { .. })
    ));
}

#[test]
fn certified_missing_skeleton_is_not_an_empty_upper_boundary() -> Result<()> {
    let p = PointCloudView::new(&[0., 0., 1., 0., 0., 1.], 3, 2)?;
    let source = RipsBuilder::from_points(p).build_complex(1)?;
    let k = source.stage(2.)?;
    assert!(matches!(
        PersistentLaplacian::compute(&k, &k, 1, SpectralOptions::default(), &Execution::default()),
        Err(Error::InsufficientSkeleton { .. })
    ));
    let bare = source.complex().stage(2.)?;
    assert!(
        PersistentLaplacian::compute(
            &bare,
            &k,
            0,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    assert!(
        PersistentLaplacian::compute(
            &bare,
            &bare,
            1,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_ok()
    );
    let result =
        PersistentLaplacian::compute(&k, &k, 0, SpectralOptions::default(), &Execution::default())?;
    assert!(result.source_context().is_some());
    assert_eq!(result.construction_dimension(), Some(1));
    Ok(())
}

#[test]
fn default_dense_limit_rejects_large_input_before_gram_allocation() {
    let source = SimplicialComplex::new(
        (0..257)
            .map(|v| Simplex::new(vec![v], 0.).unwrap())
            .collect(),
    )
    .unwrap();
    let k = source.stage(0.).unwrap();
    assert!(matches!(
        PersistentLaplacian::compute(&k, &k, 0, SpectralOptions::default(), &Execution::default()),
        Err(Error::InvalidParameter {
            parameter: "max_dense_entries",
            ..
        })
    ));
}

#[test]
fn real_rp2_homology_does_not_inherit_f2_torsion() {
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
    let mut faces = std::collections::BTreeSet::new();
    for facet in facets {
        for mask in 1..8 {
            faces.insert(
                facet
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| mask & (1 << i) != 0)
                    .map(|(_, v)| *v)
                    .collect::<Vec<_>>(),
            );
        }
    }
    let source = SimplicialComplex::new(
        faces
            .into_iter()
            .map(|v| {
                let value = (v.len() - 1) as f64;
                Simplex::new(v, value).unwrap()
            })
            .collect(),
    )
    .unwrap();
    // H1(RP2;R)=0, whereas H1(RP2;F2)=F2. The oracle tool verifies
    // integer/Q ranks and separately exhibits GUDHI's two field counts.
    assert!(
        pair(&source, 2., 2., 1)
            .spectrum(&Execution::default())
            .unwrap()
            .eigenvalues
            .iter()
            .all(|v| *v > 1e-8)
    );
}

#[test]
fn long_reverse_elimination_matches_hand_derived_path_and_cone() {
    // Eliminating 254 interior vertices gives conductance 1/255 between
    // the retained endpoints; this stresses many successive positive pivots.
    let mut cells: Vec<_> = (0..256)
        .map(|i| Simplex::new(vec![i], if i == 0 || i == 255 { 0.0 } else { 1.0 }).unwrap())
        .collect();
    cells.extend((0..255).map(|i| Simplex::new(vec![i, i + 1], 2.0).unwrap()));
    let path = SimplicialComplex::new(cells).unwrap();
    let k = path.stage(0.0).unwrap();
    let l = path.stage(2.0).unwrap();
    let expected = [1.0 / 255.0, -1.0 / 255.0, -1.0 / 255.0, 1.0 / 255.0];
    close(pair(&path, 0.0, 2.0, 0).matrix(), &expected);
    let reverse = fixed_endpoint_laplacians(
        &[k],
        &l,
        0,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    close(reverse[0].matrix(), &expected);

    // Cone over a path: 255 edges, 127 faces. The removed upper block is
    // singular. K contains two spokes forming a tree, so its persistent upper
    // term is exactly zero and the full operator is its lower Gram.
    let mut cells: Vec<_> = (0..129)
        .map(|v| Simplex::new(vec![v], 0.0).unwrap())
        .collect();
    cells
        .extend((1..129).map(|v| {
            Simplex::new(vec![0, v], if v == 1 || v == 128 { 0.0 } else { 1.0 }).unwrap()
        }));
    cells.extend((1..128).map(|v| Simplex::new(vec![v, v + 1], 1.0).unwrap()));
    cells.extend((1..128).map(|v| Simplex::new(vec![0, v, v + 1], 2.0).unwrap()));
    let cone = SimplicialComplex::new(cells).unwrap();
    let k = cone.stage(0.0).unwrap();
    let l = cone.stage(2.0).unwrap();
    let expected = [2.0, 1.0, 1.0, 2.0];
    close(pair(&cone, 0.0, 2.0, 1).matrix(), &expected);
    let reverse = fixed_endpoint_laplacians(
        &[k],
        &l,
        1,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    close(reverse[0].matrix(), &expected);
}
