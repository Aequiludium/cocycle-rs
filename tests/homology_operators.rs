//! Public F2 operator contracts and an independent exhaustive tiny-span oracle.
use cocycle::Error;
use cocycle::execution::Execution;
use cocycle::homology_operators::{
    ChainWindow, F2Matrix, HomologyOperator, OperatorFamily, StageInterval, WeightPolicy,
};
use std::collections::HashSet;
use std::sync::atomic::AtomicBool;

fn matrix(rows: &[&[u8]], n: usize) -> F2Matrix {
    F2Matrix::from_rows(&rows.iter().map(|r| r.to_vec()).collect::<Vec<_>>(), n).unwrap()
}

#[test]
fn minimum_class_matches_independent_coset_enumeration_and_lexicographic_ties() {
    let controls = Execution::default();
    for mask in 0..64_u32 {
        let rows: Vec<Vec<u8>> = (0..3)
            .map(|i| (0..2).map(|j| ((mask >> (2 * i + j)) & 1) as u8).collect())
            .collect();
        let op = HomologyOperator::feasible(
            ChainWindow::new(
                0,
                [
                    F2Matrix::zero(0, 3).unwrap(),
                    F2Matrix::from_rows(&rows, 2).unwrap(),
                ],
                [vec![], vec![0, 1, 2], vec![3, 4]],
                vec![5, 1, 3],
                "tiny",
                "cost",
            )
            .unwrap(),
        )
        .unwrap();
        for bits in 0..8 {
            let z: Vec<u8> = (0..3).map(|i| ((bits >> i) & 1) as u8).collect();
            let expected = (0..4)
                .map(|c| {
                    let chain: Vec<u8> = (0..3)
                        .map(|i| z[i] ^ (rows[i][0] * (c & 1)) ^ (rows[i][1] * ((c >> 1) & 1)))
                        .collect();
                    let mass = chain
                        .iter()
                        .zip([5_u64, 1, 3])
                        .map(|(&b, w)| u64::from(b) * w)
                        .sum::<u64>();
                    (mass, chain)
                })
                .min()
                .unwrap();
            let minimum = op.minimum_class_mass(&z, &controls).unwrap();
            assert_eq!(
                (minimum.mass(), minimum.chain()),
                (expected.0, expected.1.as_slice())
            );
            let residual = op.window().boundaries()[1]
                .apply(minimum.boundary_witness().unwrap())
                .unwrap();
            assert_eq!(
                residual,
                z.iter()
                    .zip(minimum.chain())
                    .map(|(a, b)| a ^ b)
                    .collect::<Vec<_>>()
            );
        }
    }
    let op = HomologyOperator::feasible(
        ChainWindow::new(
            0,
            [F2Matrix::zero(0, 2).unwrap(), matrix(&[&[1], &[1]], 1)],
            [vec![], vec![0, 1], vec![2]],
            vec![10, 1],
            "edge",
            "cost",
        )
        .unwrap(),
    )
    .unwrap();
    assert_eq!(op.class_query(&[1, 0]).unwrap().selected_mass(), 10);
    assert_eq!(op.minimum_class_mass(&[1, 0], &controls).unwrap().mass(), 1);
    for controls in [
        Execution::default().max_work(0),
        Execution::default().cancellation(&AtomicBool::new(true)),
    ] {
        assert!(op.minimum_class_mass(&[1, 0], &controls).is_err());
    }
    assert!(
        op.minimum_class_mass(&[2, 0], &Execution::default())
            .is_err()
    );
    assert!(
        HomologyOperator::feasible(triangle(false))
            .unwrap()
            .minimum_class_mass(&[1, 0, 0], &Execution::default())
            .is_err()
    );
}

#[test]
fn minimum_filling_uses_upper_costs_and_distinguishes_zero_infeasible_and_failure() {
    let d = matrix(&[&[1, 1, 1], &[1, 1, 1]], 3);
    let lower = ChainWindow::new(
        0,
        [F2Matrix::zero(0, 2).unwrap(), d.clone()],
        [vec![], vec![0, 1], vec![2, 3, 4]],
        vec![100, 100],
        "edge",
        "lower cost",
    )
    .unwrap();
    let upper = ChainWindow::new(
        1,
        [d, F2Matrix::zero(3, 0).unwrap()],
        [vec![0, 1], vec![2, 3, 4], vec![]],
        vec![7, 2, 2],
        "edge",
        "upper cost",
    )
    .unwrap();
    let op = HomologyOperator::feasible(lower).unwrap();
    let controls = Execution::default();
    let minimum = op
        .minimum_filling(&[1, 1], &upper, &controls)
        .unwrap()
        .unwrap();
    assert_eq!((minimum.chain(), minimum.mass()), (&[0, 1, 0][..], 2));
    assert_eq!(minimum.boundary(), &[1, 1]);
    assert_eq!(minimum.window().unit(), "upper cost");
    assert_eq!(
        op.minimum_filling(&[0, 0], &upper, &controls)
            .unwrap()
            .unwrap()
            .mass(),
        0
    );
    assert!(
        op.minimum_filling(&[1, 0], &upper, &controls)
            .unwrap()
            .is_none()
    );
    assert!(
        op.minimum_filling(&[1, 1], &upper, &controls.max_work(0))
            .is_err()
    );
    let foreign = ChainWindow::new(
        1,
        upper.boundaries().clone(),
        upper.bases().clone(),
        upper.weights().to_vec(),
        "foreign",
        "upper cost",
    )
    .unwrap();
    assert!(
        op.minimum_filling(&[1, 1], &foreign, &Execution::default())
            .is_err()
    );
    let huge = HomologyOperator::feasible(
        ChainWindow::new(
            0,
            [F2Matrix::zero(0, 2).unwrap(), F2Matrix::zero(2, 0).unwrap()],
            [vec![], vec![0, 1], vec![]],
            vec![u64::MAX, 1],
            "huge",
            "cost",
        )
        .unwrap(),
    )
    .unwrap();
    assert!(matches!(
        huge.minimum_class_mass(&[1, 1], &Execution::default()),
        Err(Error::NumericalFailure { .. })
    ));
}

#[test]
fn feature_histories_correct_earlier_generators_and_bind_family_runs() {
    let window = |filled| {
        ChainWindow::new(
            0,
            [
                F2Matrix::zero(0, 2).unwrap(),
                if filled {
                    matrix(&[&[1], &[1]], 1)
                } else {
                    F2Matrix::zero(2, 0).unwrap()
                },
            ],
            [vec![], vec![0, 1], if filled { vec![2] } else { vec![] }],
            vec![1, 1],
            "merge",
            "cost",
        )
        .unwrap()
    };
    let ops = vec![
        HomologyOperator::feasible(window(false)).unwrap(),
        HomologyOperator::feasible(window(true)).unwrap(),
    ];
    let family = OperatorFamily::new(vec![0., 1.], ops.clone()).unwrap();
    let features = family.features(&Execution::default()).unwrap();
    assert_eq!(features.len(), 2);
    let dying = features.iter().find(|f| f.death() == Some(1)).unwrap();
    // The second original basis vector alone does not die: the full corrected
    // history is e0+e1, whose raw boundary preimage is the new edge.
    assert_eq!(
        family
            .resolve_feature(dying, 0, &Execution::default())
            .unwrap()
            .unwrap()
            .representative(),
        &[1, 1]
    );
    assert!(
        family
            .resolve_feature(dying, 1, &Execution::default())
            .unwrap()
            .is_none()
    );
    let foreign = OperatorFamily::new(vec![0., 1.], ops).unwrap();
    assert!(
        foreign
            .resolve_feature(dying, 0, &Execution::default())
            .is_err()
    );
    assert!(
        family
            .clone()
            .resolve_feature(dying, 0, &Execution::default())
            .is_ok()
    );
    assert!(family.features(&Execution::default().max_work(0)).is_err());
    let same = OperatorFamily::new(
        vec![0.],
        vec![HomologyOperator::feasible(window(false)).unwrap()],
    )
    .unwrap();
    let repeated = same.features(&Execution::default()).unwrap();
    assert_eq!(
        repeated.iter().map(|f| f.occurrence()).collect::<Vec<_>>(),
        vec![0, 1]
    );
}

#[test]
fn simplex_features_pair_cycles_and_volumes_preserve_original_coordinates() {
    use cocycle::complex::{Simplex, SimplicialComplex};
    let source = SimplicialComplex::new(vec![
        Simplex::new(vec![10], 0.).unwrap(),
        Simplex::new(vec![20], 0.).unwrap(),
        Simplex::new(vec![30], 0.).unwrap(),
        Simplex::new(vec![10, 20], 1.).unwrap(),
        Simplex::new(vec![10, 30], 1.).unwrap(),
        Simplex::new(vec![20, 30], 1.).unwrap(),
        Simplex::new(vec![10, 20, 30], 2.).unwrap(),
    ])
    .unwrap();
    let controls = Execution::default();
    let family =
        OperatorFamily::from_simplices(&source, 1, vec![1; 3], "triangle", "edges", &controls)
            .unwrap();
    let features = family.features(&controls).unwrap();
    assert_eq!(features.len(), 1);
    let f = &features[0];
    assert_eq!((f.birth(), f.death()), (5, Some(6)));
    let (birth, death) = family.feature_pairing(f).unwrap();
    assert_eq!(source.simplex(birth).unwrap().dimension(), 1);
    assert_eq!(
        source.simplex(death.unwrap()).unwrap().vertices(),
        &[10, 20, 30]
    );
    let cycle = family.optimal_1_cycle(f, &controls).unwrap();
    assert_eq!(cycle.mass(), 3);
    let location = family.locate_chain(f.birth(), cycle.chain()).unwrap();
    assert!(
        location.contains(&vec![10, 20])
            && location.contains(&vec![10, 30])
            && location.contains(&vec![20, 30])
    );
    let upper = family
        .filling_window(f, vec![7], "face cost", &controls)
        .unwrap();
    let volume = family
        .optimal_volume(f, &upper, &controls)
        .unwrap()
        .unwrap();
    assert_eq!(
        (volume.chain(), volume.boundary(), volume.mass()),
        (&[1][..], &[1, 1, 1][..], 7)
    );
    let stable = family
        .stable_volume(f, &upper, 1, "filtration time", &controls)
        .unwrap()
        .unwrap();
    assert_eq!(stable.chain(), &[1]);
    assert_eq!(stable.examined(), 1); // Positive-cost death-anchor certificate.
    assert!(
        family
            .stable_volume(f, &upper, 0, "time", &controls)
            .is_err()
    );
    assert!(family.stable_volume(f, &upper, 1, "", &controls).is_err());
    assert!(
        family
            .optimal_volume(f, &upper, &controls.max_work(0))
            .is_err()
    );
    let fractional = SimplicialComplex::new(
        source
            .simplices()
            .iter()
            .map(|s| Simplex::new(s.vertices().to_vec(), s.value() / 2.).unwrap())
            .collect(),
    )
    .unwrap();
    let fractional =
        OperatorFamily::from_simplices(&fractional, 1, vec![1; 3], "triangle", "edge", &controls)
            .unwrap();
    let fs = fractional.features(&controls).unwrap();
    let upper = fractional
        .filling_window(&fs[0], vec![1], "faces", &controls)
        .unwrap();
    assert!(
        fractional
            .stable_volume(&fs[0], &upper, 1, "time", &controls)
            .is_err()
    );
}
fn triangle(filled: bool) -> ChainWindow<usize> {
    let a = matrix(&[&[1, 1, 0], &[1, 0, 1], &[0, 1, 1]], 3);
    let d = if filled {
        matrix(&[&[1], &[1], &[1]], 1)
    } else {
        F2Matrix::zero(3, 0).unwrap()
    };
    ChainWindow::new(
        1,
        [a, d],
        [
            vec![10, 20, 30],
            vec![40, 50, 60],
            if filled { vec![70] } else { vec![] },
        ],
        vec![2, 3, 5],
        "triangle-snapshot",
        "cost",
    )
    .unwrap()
}

#[test]
fn ring_and_filling_have_hand_derived_class_queries_and_transport() {
    let ring = HomologyOperator::feasible(triangle(false)).unwrap();
    let filled = HomologyOperator::feasible(triangle(true)).unwrap();
    assert_eq!(ring.betti(), 1);
    assert_eq!(filled.betti(), 0);
    assert_eq!(ring.class_query(&[1, 1, 1]).unwrap().selected_mass(), 10);
    assert_eq!(filled.class_query(&[1, 1, 1]).unwrap().selected_mass(), 0);
    assert!(ring.class_query(&[1, 0, 0]).is_err());
    assert!(!ring.is_boundary(&[1, 1, 1]).unwrap());
    assert!(filled.is_boundary(&[1, 1, 1]).unwrap());
    assert_eq!(
        ring.stretch_with(&Execution::default()).unwrap(),
        Some((10, 10))
    );
    let family = OperatorFamily::new(vec![1., 1., 2.], vec![ring.clone(), ring, filled]).unwrap();
    assert_eq!(family.transport_rank(0, 1).unwrap(), 1);
    assert_eq!(family.transport_rank(0, 2).unwrap(), 0);
    assert_eq!(
        family
            .track_class(&[1, 1, 1], 0, 2)
            .unwrap()
            .selected_mass(),
        0
    );
    assert_eq!(
        family.barcode().unwrap(),
        vec![StageInterval {
            birth: 0,
            death: Some(2),
            multiplicity: 1
        }]
    );
    assert!(matches!(
        family.transport_with(0, 1, &Execution::default().max_work(0)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(matches!(
        family.barcode_with(&Execution::default().max_work(1)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(family.transport(2, 0).is_err());
}

#[test]
fn empty_spaces_h0_and_terminal_unknown_are_distinct() {
    let empty = ChainWindow::<usize>::new(
        0,
        [F2Matrix::zero(0, 0).unwrap(), F2Matrix::zero(0, 0).unwrap()],
        [vec![], vec![], vec![]],
        vec![],
        "empty",
        "cost",
    )
    .unwrap();
    let empty = HomologyOperator::feasible(empty).unwrap();
    assert_eq!(empty.betti(), 0);
    assert_eq!(empty.class_query(&[]).unwrap().selected_mass(), 0);
    assert_eq!(empty.stretch_with(&Execution::default()).unwrap(), None);
    let h0 = ChainWindow::new(
        0,
        [F2Matrix::zero(0, 2).unwrap(), matrix(&[&[1], &[1]], 1)],
        [vec![], vec![0, 1], vec![2]],
        vec![1, 2],
        "path",
        "cost",
    )
    .unwrap();
    let op = HomologyOperator::feasible(h0).unwrap();
    assert_eq!(op.betti(), 1);
    let family = OperatorFamily::new(vec![0.], vec![op]).unwrap();
    assert_eq!(
        family.barcode().unwrap(),
        vec![StageInterval {
            birth: 0,
            death: None,
            multiplicity: 1
        }]
    );
    assert!(F2Matrix::zero(0, 3).unwrap().solve(&[]).unwrap().is_some());
    assert_eq!(
        F2Matrix::zero(2, 0).unwrap().solve(&[0, 0]).unwrap(),
        Some(vec![])
    );
    assert_eq!(F2Matrix::zero(2, 0).unwrap().solve(&[1, 0]).unwrap(), None);
}

#[test]
fn imported_bad_projection_and_identity_mixing_are_rejected() {
    let ring = triangle(false);
    assert!(
        HomologyOperator::from_projection(
            ring.clone(),
            F2Matrix::zero(3, 3).unwrap(),
            &Execution::default()
        )
        .is_err()
    );
    let first = HomologyOperator::feasible(ring.clone()).unwrap();
    let second = HomologyOperator::from_projection(
        ring,
        first.projection_matrix().unwrap(),
        &Execution::default(),
    )
    .unwrap();
    let q = first.class_query(&[1, 1, 1]).unwrap();
    assert!(
        q.same_class(&second.class_query(&[1, 1, 1]).unwrap())
            .is_err()
    );
    assert!(
        q.same_class(&first.clone().class_query(&[1, 1, 1]).unwrap())
            .unwrap()
    );
    drop(first);
    assert_eq!(q.window().bases()[1], vec![40, 50, 60]);
    assert_eq!(q.selected_mass(), 10);
    assert!(F2Matrix::from_rows(&[vec![2]], 1).is_err());
    assert!(
        ChainWindow::new(
            1,
            [matrix(&[&[1]], 1), matrix(&[&[1]], 1)],
            [vec![0], vec![1], vec![2]],
            vec![1],
            "bad",
            "cost"
        )
        .is_err()
    );
    assert!(
        ChainWindow::new(
            1,
            [F2Matrix::zero(0, 2).unwrap(), F2Matrix::zero(2, 0).unwrap()],
            [vec![], vec![1, 1], vec![]],
            vec![1, 1],
            "bad",
            "cost"
        )
        .is_err()
    );
}

#[test]
fn cancellation_limits_recovery_and_integer_overflow_are_explicit() {
    let flag = AtomicBool::new(true);
    assert!(matches!(
        HomologyOperator::feasible_with(triangle(false), &Execution::default().cancellation(&flag)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        HomologyOperator::feasible_with(triangle(false), &Execution::default().max_work(0)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let op = HomologyOperator::feasible(triangle(false)).unwrap();
    assert!(matches!(
        op.class_queries_with(&[&[1, 1, 1]], &Execution::default().max_work(0)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(op.class_query(&[1, 1, 1]).is_ok());
    assert!(matches!(
        op.stretch_with(&Execution::default().max_work(1)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let window = ChainWindow::new(
        1,
        [F2Matrix::zero(0, 2).unwrap(), F2Matrix::zero(2, 0).unwrap()],
        [vec![], vec![0, 1], vec![]],
        vec![u64::MAX, 1],
        "overflow",
        "cost",
    )
    .unwrap();
    assert!(matches!(
        HomologyOperator::feasible(window)
            .unwrap()
            .class_query(&[1, 1]),
        Err(Error::NumericalFailure { .. })
    ));
}

#[test]
fn multiword_coordinates_cross_64_and_128_without_padding_leaks() {
    for n in [63, 64, 65, 127, 128, 129] {
        let mut rows = vec![vec![0; n]; 2];
        rows[0][0] = 1;
        rows[0][n - 1] = 1;
        rows[1][n / 2] = 1;
        let a = F2Matrix::from_rows(&rows, n).unwrap();
        assert_eq!(a.rank().unwrap(), 2);
        assert_eq!(a.kernel_basis().unwrap().len(), n - 2);
        for z in a.kernel_basis().unwrap() {
            assert_eq!(a.apply(&z).unwrap(), vec![0, 0]);
        }
        let w = ChainWindow::new(
            1,
            [a, F2Matrix::zero(n, 0).unwrap()],
            [vec![0, 1], (0..n).collect(), vec![]],
            vec![1; n],
            "multiword",
            "cost",
        )
        .unwrap();
        let op = HomologyOperator::feasible(w).unwrap();
        assert_eq!(op.betti(), n - 2);
        let mut x = vec![0; n];
        x[n - 1] = 1;
        let p = op.project(&x).unwrap();
        assert_eq!(p[0], 1);
        assert_eq!(p[n - 1], 1);
        assert_eq!(op.project(&p).unwrap(), p);
    }
}

#[test]
fn tiny_exhaustive_span_oracle_checks_quotients_without_elimination() {
    // This expectation enumerates sets directly; it uses no production rank,
    // inverse, kernel-basis or projection construction as an expected generator.
    for n in 0..=5 {
        for seed in 0..18_u64 {
            let bits = |x: usize| (0..n).map(|j| ((x >> j) & 1) as u8).collect::<Vec<_>>();
            let a_rows = vec![
                bits((seed as usize * 7 + 3) % (1 << n)),
                bits((seed as usize * 11 + 1) % (1 << n)),
            ];
            let cycles: Vec<_> = (0..1_usize << n)
                .filter(|&x| {
                    a_rows.iter().all(|row| {
                        row.iter()
                            .enumerate()
                            .fold(0_u8, |p, (j, &b)| p ^ (b * ((x >> j) & 1) as u8))
                            == 0
                    })
                })
                .collect();
            let ds = [
                cycles[(seed as usize) % cycles.len()],
                cycles[(seed as usize * 3 + 1) % cycles.len()],
            ];
            let boundary_set: HashSet<_> = [0, ds[0], ds[1], ds[0] ^ ds[1]].into_iter().collect();
            let a = F2Matrix::from_rows(&a_rows, n).unwrap();
            let d = F2Matrix::from_columns(&ds.iter().map(|&x| bits(x)).collect::<Vec<_>>(), n)
                .unwrap();
            let w = ChainWindow::new(
                1,
                [a, d],
                [vec![0, 1], (0..n).collect(), vec![0, 1]],
                vec![1; n],
                "tiny",
                "cost",
            )
            .unwrap();
            let op = HomologyOperator::feasible(w).unwrap();
            assert_eq!(1_usize << op.betti(), cycles.len() / boundary_set.len());
            let projected: HashSet<_> = cycles
                .iter()
                .map(|&x| op.project(&bits(x)).unwrap())
                .collect();
            assert_eq!(projected.len(), cycles.len() / boundary_set.len());
            for &x in &cycles {
                let p = op.project(&bits(x)).unwrap();
                let p_mask = p
                    .iter()
                    .enumerate()
                    .fold(0, |v, (j, &b)| v | ((b as usize) << j));
                assert!(boundary_set.contains(&(x ^ p_mask)));
                assert_eq!(op.is_boundary(&bits(x)).unwrap(), boundary_set.contains(&x));
            }
        }
    }
}

#[test]
fn incompatible_stages_costs_and_basis_permutations_are_checked() {
    let ring = triangle(false);
    let a = ring.boundaries()[0].clone();
    let d = ring.boundaries()[1].clone();
    let changed = ChainWindow::new(
        1,
        [a, d],
        ring.bases().clone(),
        vec![2, 3, 6],
        ring.source(),
        ring.unit(),
    )
    .unwrap();
    let ops = vec![
        HomologyOperator::feasible(ring.clone()).unwrap(),
        HomologyOperator::feasible(changed).unwrap(),
    ];
    assert!(OperatorFamily::new(vec![1., 2.], ops.clone()).is_err());
    assert!(OperatorFamily::with_weight_policy(vec![1., 2.], ops, WeightPolicy::Variable).is_ok());
    let a = matrix(&[&[0, 1, 1], &[1, 0, 1], &[1, 1, 0]], 3);
    let permuted = ChainWindow::new(
        1,
        [a, F2Matrix::zero(3, 0).unwrap()],
        [vec![10, 20, 30], vec![60, 50, 40], vec![]],
        vec![5, 3, 2],
        ring.source(),
        ring.unit(),
    )
    .unwrap();
    let family = OperatorFamily::new(
        vec![1., 2.],
        vec![
            HomologyOperator::feasible(ring).unwrap(),
            HomologyOperator::feasible(permuted).unwrap(),
        ],
    )
    .unwrap();
    assert_eq!(family.transport_rank(0, 1).unwrap(), 1);
    assert_eq!(
        family
            .track_class(&[1, 1, 1], 0, 1)
            .unwrap()
            .selected_mass(),
        10
    );
}

#[test]
fn adjacent_transport_composition_matches_direct_maps_with_reordering() {
    // Three H0 vertices, then one merge, then two merges. Every chain is a cycle.
    let bases = [vec![0, 1, 2], vec![2, 0, 1], vec![1, 2, 0]];
    let boundaries = [
        F2Matrix::zero(3, 0).unwrap(),
        matrix(&[&[0], &[1], &[1]], 1),
        matrix(&[&[1, 1], &[0, 1], &[1, 0]], 2),
    ];
    let operators = (0..3)
        .map(|i| {
            HomologyOperator::feasible(
                ChainWindow::new(
                    0,
                    [F2Matrix::zero(0, 3).unwrap(), boundaries[i].clone()],
                    [vec![], bases[i].clone(), (0..i).map(|j| 10 + j).collect()],
                    vec![1; 3],
                    "merge-tree",
                    "cost",
                )
                .unwrap(),
            )
            .unwrap()
        })
        .collect();
    let family = OperatorFamily::new(vec![0., 0., 1.], operators).unwrap();
    let first = family.transport(0, 1).unwrap();
    let second = family.transport(1, 2).unwrap();
    let direct = family.transport(0, 2).unwrap();
    for j in 0..direct.ncols() {
        assert_eq!(
            second.apply(&first.column(j).unwrap()).unwrap(),
            direct.column(j).unwrap()
        );
    }
    assert_eq!(family.transport_rank(0, 1).unwrap(), 2);
    assert_eq!(family.transport_rank(0, 2).unwrap(), 1);
    assert_eq!(
        family.barcode().unwrap(),
        vec![
            StageInterval {
                birth: 0,
                death: Some(1),
                multiplicity: 1
            },
            StageInterval {
                birth: 0,
                death: Some(2),
                multiplicity: 1
            },
            StageInterval {
                birth: 0,
                death: None,
                multiplicity: 1
            },
        ]
    );
}

#[test]
fn incremental_stretch_readout_keeps_binary_order_and_exact_objective() {
    let window = ChainWindow::new(
        0,
        [
            F2Matrix::zero(0, 3).unwrap(),
            matrix(&[&[1], &[1], &[0]], 1),
        ],
        [vec![], vec![0, 1, 2], vec![3]],
        vec![9, 2, 4],
        "stretch",
        "cost",
    )
    .unwrap();
    let operator = HomologyOperator::feasible(window).unwrap();
    let weights = operator.window().weights();
    let mut expected = (0_u64, 1_u64);
    for code in 1..8 {
        let x: Vec<_> = (0..3).map(|j| ((code >> j) & 1) as u8).collect();
        let projected = operator.project(&x).unwrap();
        let mass = |x: &[u8]| {
            x.iter()
                .zip(weights)
                .filter(|(bit, _)| **bit != 0)
                .map(|(_, w)| *w)
                .sum::<u64>()
        };
        let ratio = (mass(&projected), mass(&x));
        if u128::from(ratio.0) * u128::from(expected.1)
            > u128::from(expected.0) * u128::from(ratio.1)
        {
            expected = ratio;
        }
    }
    assert_eq!(
        operator.stretch_with(&Execution::default()).unwrap(),
        Some(expected)
    );
}

fn simplicial_triangle() -> cocycle::complex::SimplicialComplex {
    use cocycle::complex::{Simplex, SimplicialComplex};
    SimplicialComplex::new(
        [
            (vec![10], 0.),
            (vec![20], 0.),
            (vec![30], 0.),
            (vec![10, 20], 1.),
            (vec![10, 30], 1.),
            (vec![20, 30], 1.),
            (vec![10, 20, 30], 2.),
        ]
        .into_iter()
        .map(|(vertices, value)| Simplex::new(vertices, value).unwrap())
        .collect(),
    )
    .unwrap()
}

#[test]
fn borrowed_stage_family_owns_results_and_replays_boundary_witness() {
    let family = {
        let source = simplicial_triangle();
        let stages = [-1., 0., 1., 1., 2.].map(|scale| source.stage(scale).unwrap());
        OperatorFamily::from_stages(
            &stages,
            1,
            vec![vec![], vec![], vec![2, 3, 5], vec![2, 3, 5], vec![2, 3, 5]],
            "triangle",
            "cost",
            &Execution::default(),
        )
        .unwrap()
    }; // No stage/source borrow survives in the family or query.
    assert_eq!(family.scales(), &[-1., 0., 1., 1., 2.]);
    assert_eq!(
        (0..5)
            .map(|i| family.stage(i).unwrap().betti())
            .collect::<Vec<_>>(),
        vec![0, 0, 1, 1, 0]
    );
    assert_eq!(
        family.barcode().unwrap(),
        vec![StageInterval {
            birth: 2,
            death: Some(4),
            multiplicity: 1
        }]
    );
    let query = family.stage(2).unwrap().class_query(&[1, 1, 1]).unwrap();
    assert_eq!(
        family.track_query(&query, 2, 3).unwrap().selected_mass(),
        10
    );
    let target = family.stage(4).unwrap();
    let tracked = family.track_query(&query, 2, 4).unwrap();
    assert_eq!(tracked.representative(), &[0, 0, 0]);
    let filling = target
        .boundary_preimage(query.representative())
        .unwrap()
        .unwrap();
    // A hand-derived unique face fills the ring; replay actual D independently.
    assert_eq!(filling, vec![1]);
    assert_eq!(
        target.window().boundaries()[1].apply(&filling).unwrap(),
        vec![1, 1, 1]
    );
    assert_eq!(
        family
            .stage(2)
            .unwrap()
            .boundary_preimage(&[1, 1, 1])
            .unwrap(),
        None
    );
    assert_eq!(
        family
            .stage(2)
            .unwrap()
            .boundary_preimage(&[0, 0, 0])
            .unwrap(),
        Some(vec![])
    );
    assert!(target.boundary_preimage(&[1, 0, 0]).unwrap().is_none());
    assert!(target.boundary_preimage(&[2, 0, 0]).is_err());
    assert!(target.boundary_preimage(&[0]).is_err());
}

#[test]
fn stage_factory_rejects_foreign_owners_and_query_runs() {
    use cocycle::filtration::RipsBuilder;
    use cocycle::geometry::PointCloudView;
    let coordinates = [0., 1., 2.];
    let certified = RipsBuilder::from_points(PointCloudView::new(&coordinates, 3, 1).unwrap())
        .max_edge_length(1.)
        .build_complex(1)
        .unwrap();
    let certified_stages = [certified.stage(0.).unwrap(), certified.stage(1.).unwrap()];
    let h0 = OperatorFamily::from_stages(
        &certified_stages,
        0,
        vec![vec![1; 3]; 2],
        "path",
        "cost",
        &Execution::default(),
    )
    .unwrap();
    assert_eq!(
        (h0.stage(0).unwrap().betti(), h0.stage(1).unwrap().betti()),
        (3, 1)
    );
    // Same storage, but mixing bare and certified stage interpretations is invalid.
    assert!(
        OperatorFamily::from_stages(
            &[certified_stages[0], certified.complex().stage(1.).unwrap()],
            0,
            vec![vec![1; 3]; 2],
            "path",
            "cost",
            &Execution::default(),
        )
        .is_err()
    );
    let source = simplicial_triangle();
    let copy = source.clone();
    let stages = [source.stage(1.).unwrap(), source.stage(2.).unwrap()];
    let costs = vec![vec![2, 3, 5]; 2];
    assert!(
        OperatorFamily::from_stages(
            &[stages[0], copy.stage(2.).unwrap()],
            1,
            costs.clone(),
            "same-label",
            "cost",
            &Execution::default(),
        )
        .is_err()
    );
    assert!(
        OperatorFamily::from_stages(
            &[stages[1], stages[0]],
            1,
            costs.clone(),
            "same-label",
            "cost",
            &Execution::default()
        )
        .is_err()
    );
    assert!(
        OperatorFamily::from_stages(
            &stages,
            1,
            vec![vec![2, 3, 5]],
            "triangle",
            "cost",
            &Execution::default()
        )
        .is_err()
    );
    assert!(
        OperatorFamily::from_stages(
            &stages,
            1,
            vec![vec![2, 3, 5], vec![2, 4, 5]],
            "triangle",
            "cost",
            &Execution::default()
        )
        .is_err()
    );
    let family =
        OperatorFamily::from_stages(&stages, 1, costs, "triangle", "cost", &Execution::default())
            .unwrap();
    let source_run = family.stage(0).unwrap().clone();
    let query = source_run.class_query(&[1, 1, 1]).unwrap();
    assert_eq!(family.track_query(&query, 0, 1).unwrap().selected_mass(), 0);
    let foreign_run = HomologyOperator::feasible(source_run.window().clone()).unwrap();
    let foreign_query = foreign_run.class_query(&[1, 1, 1]).unwrap();
    assert!(family.track_query(&foreign_query, 0, 1).is_err());
    assert!(family.track_query(&query, 1, 1).is_err());
    assert!(family.track_query(&query, 0, 2).is_err());
    assert_eq!(
        family
            .clone()
            .track_query(&query, 0, 1)
            .unwrap()
            .selected_mass(),
        0
    );
}

#[test]
fn composed_stage_and_query_operations_share_controls_and_recover() {
    let source = simplicial_triangle();
    let stages = [source.stage(1.).unwrap(), source.stage(2.).unwrap()];
    let costs = vec![vec![2, 3, 5]; 2];
    // Find a quota that completes each extraction and solve independently.
    // The composed operation must not reset that quota between terminals.
    let separate_limit = (0..5000)
        .find(|&limit| {
            stages.iter().all(|stage| {
                let controls = Execution::default().max_work(limit);
                ChainWindow::from_complex_with(
                    stage,
                    1,
                    None,
                    vec![2, 3, 5],
                    "triangle",
                    "cost",
                    &controls,
                )
                .is_ok_and(|window| HomologyOperator::feasible_with(window, &controls).is_ok())
            })
        })
        .expect("small hand fixture fits a finite cooperative quota");
    assert!(matches!(
        OperatorFamily::from_stages(
            &stages,
            1,
            costs.clone(),
            "triangle",
            "cost",
            &Execution::default().max_work(separate_limit)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let limits = Execution::default().max_work(0);
    assert!(matches!(
        ChainWindow::from_complex_with(
            &stages[0],
            1,
            None,
            costs[0].clone(),
            "triangle",
            "cost",
            &limits
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(matches!(
        OperatorFamily::from_stages(&stages, 1, costs.clone(), "triangle", "cost", &limits),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let cancelled = AtomicBool::new(true);
    let cancellation = Execution::default().cancellation(&cancelled);
    assert!(matches!(
        OperatorFamily::from_stages(&stages, 1, costs.clone(), "triangle", "cost", &cancellation),
        Err(Error::Cancelled)
    ));
    let family =
        OperatorFamily::from_stages(&stages, 1, costs, "triangle", "cost", &Execution::default())
            .unwrap();
    let ring = family.stage(0).unwrap();
    let query = ring.class_query(&[1, 1, 1]).unwrap();
    for controls in [limits, cancellation] {
        assert!(ring.project_with(&[1, 1, 1], &controls).is_err());
        assert!(ring.class_query_with(&[1, 1, 1], &controls).is_err());
        assert!(
            family
                .track_class_with(&[1, 1, 1], 0, 1, &controls)
                .is_err()
        );
        assert!(family.track_query_with(&query, 0, 1, &controls).is_err());
        assert!(
            family
                .stage(1)
                .unwrap()
                .boundary_preimage_with(&[1, 1, 1], &controls)
                .is_err()
        );
    }
    // A failed budget never corrupts the accepted family or its existing query.
    assert_eq!(family.track_query(&query, 0, 1).unwrap().selected_mass(), 0);
}
