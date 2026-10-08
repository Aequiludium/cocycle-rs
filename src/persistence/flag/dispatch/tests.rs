use super::*;
use crate::execution::Execution;
use crate::geometry::MatrixLayout;
use std::sync::atomic::AtomicBool;

#[test]
fn capped_points_matrix_and_prepared_graphs_share_bitset_selection() {
    use crate::filtration::{FlagFiltration, RipsBuilder};
    use crate::geometry::PointCloudView;
    use crate::persistence::PersistenceExt;
    // Two orthogonal centered regular simplices realize K_32,32 below
    // 1.4. Interleaved labels prevent cheap disjoint-range rejection.
    let n = 64;
    let mut coordinates = vec![0.; n * n];
    for v in 0..n {
        for k in 0..n / 2 {
            coordinates[v * n + (v % 2) * (n / 2) + k] = -1. / 32.;
        }
        coordinates[v * n + (v % 2) * (n / 2) + v / 2] += 1.;
    }
    let points = PointCloudView::new(&coordinates, n, n).unwrap();
    let distances = crate::geometry::euclidean_distances(points).unwrap();
    let matrix = DissimilarityMatrixView::new(&distances, n, MatrixLayout::LowerTriangle).unwrap();
    let points_source = RipsBuilder::from_points(points).max_edge_length(1.4);
    let matrix_source = RipsBuilder::from_distance_matrix(matrix).max_edge_length(1.4);
    let prepared = points_source.prepare().unwrap();
    let flag = FlagFiltration::new(prepared.graph().clone());
    let take_selection =
        || GRAPH_ACCESS_SELECTIONS.with(|trace| std::mem::take(&mut *trace.borrow_mut()));
    for dimension in [1, 2] {
        take_selection();
        let expected = points_source
            .persistence()
            .max_homology_dimension(dimension)
            .compute()
            .unwrap();
        assert_eq!(take_selection(), [true]);
        let dense = matrix_source
            .persistence()
            .max_homology_dimension(dimension)
            .compute()
            .unwrap();
        assert_eq!(take_selection(), [true]);
        assert_eq!(dense.diagram(), expected.diagram());
        let graph = prepared
            .persistence()
            .max_homology_dimension(dimension)
            .compute()
            .unwrap();
        assert_eq!(take_selection(), [true]);
        assert_eq!(graph.diagram(), expected.diagram());
        let supplied = flag
            .persistence()
            .max_homology_dimension(dimension)
            .compute()
            .unwrap();
        assert_eq!(take_selection(), [true]);
        assert_eq!(expected.diagram().dimension(1).unwrap().len(), 31 * 31);
        assert_eq!(supplied.diagram().dimension(1).unwrap().len(), 31 * 31);
        if dimension == 2 {
            assert!(expected.diagram().dimension(2).unwrap().is_empty());
            assert!(supplied.diagram().dimension(2).unwrap().is_empty());
        }
    }
}

#[test]
fn higher_requests_reduce_only_h2_and_above_in_generic_engine() {
    let n = 6;
    // The octahedral sphere has a genuine H2 class, whereas a filled
    // triangle's zero H1 pair must never become an H2 birth.
    for values in [
        vec![1.; 15],
        (0..n)
            .flat_map(|b| (0..b).map(move |a| if a / 2 == b / 2 { 2. } else { 1. }))
            .collect(),
    ] {
        let input = DissimilarityMatrixView::new(&values, n, MatrixLayout::LowerTriangle).unwrap();
        for cutoff in [0., 1., 2.] {
            let graph =
                crate::filtration::threshold_rips_from_distances(input, Some(cutoff)).unwrap();
            for field in [PrimeField::default(), PrimeField::new(3).unwrap()] {
                for dimension in [0, 1, 2, 3] {
                    for sparse in [false, true] {
                        simplicial::cohomology::take_reduced_dimensions();
                        let mut budget = WorkBudget::new(&Execution::default()).unwrap();
                        let actual = if sparse {
                            compute_graph(graph.graph(), dimension, cutoff, field, &mut budget)
                        } else {
                            compute_dense(input, dimension, cutoff, field, &mut budget)
                        }
                        .unwrap();
                        let trace = simplicial::cohomology::take_reduced_dimensions();
                        if field.characteristic() == 2 || dimension == 0 {
                            assert!(
                                trace.iter().all(|&q| q >= 2),
                                "generic lower dimensions ran: {trace:?}"
                            );
                            if dimension < 2 {
                                assert!(trace.is_empty());
                            }
                            if cutoff >= 1. && dimension >= 2 {
                                assert_eq!(trace.first(), Some(&2));
                            }
                        } else {
                            assert_eq!(trace.first(), Some(&0));
                            if cutoff >= 1. {
                                assert!(trace.contains(&1));
                            }
                        }
                        let access = if sparse {
                            CliqueAccess::Sparse(graph.graph(), cutoff)
                        } else {
                            CliqueAccess::Dense(input, cutoff)
                        };
                        let expected = simplicial::cohomology::compute(
                            &access,
                            dimension,
                            field,
                            &mut WorkBudget::new(&Execution::default()).unwrap(),
                        )
                        .unwrap();
                        assert_eq!(
                            crate::persistence::assemble_diagram(
                                dimension,
                                crate::diagram::Coverage::Complete,
                                actual
                            )
                            .unwrap(),
                            crate::persistence::assemble_diagram(
                                dimension,
                                crate::diagram::Coverage::Complete,
                                expected
                            )
                            .unwrap(),
                            "sparse={sparse} dimension={dimension} cutoff={cutoff} field={field:?}",
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn combined_budget_covers_low_dimensions_handoff_and_higher_reduction() {
    let values = [1.; 15];
    let input = DissimilarityMatrixView::new(&values, 6, MatrixLayout::LowerTriangle).unwrap();
    let graph = crate::filtration::threshold_rips_from_distances(input, None).unwrap();
    for sparse in [false, true] {
        let run = |budget: &mut WorkBudget<'_>| {
            if sparse {
                compute_graph(graph.graph(), 3, 1., PrimeField::default(), budget)
            } else {
                compute_dense(input, 3, 1., PrimeField::default(), budget)
            }
        };
        let execution = Execution::default().max_work(u64::MAX);
        let mut complete = WorkBudget::new(&execution).unwrap();
        let expected = run(&mut complete).unwrap();
        let work = complete.used();
        // Every checkpoint, including handoff and later dimensions, fails
        // atomically and permits retrying the same read-only source.
        for limit in 0..work {
            let execution = Execution::default().max_work(limit);
            assert_eq!(
                run(&mut WorkBudget::new(&execution).unwrap()),
                Err(crate::Error::WorkLimitExceeded { limit })
            );
        }
        for interrupt in 0..work {
            let flag = AtomicBool::new(false);
            let execution = Execution::default().max_work(u64::MAX).cancellation(&flag);
            let mut budget = WorkBudget::new(&execution).unwrap();
            budget.cancel_at_work(interrupt);
            assert_eq!(run(&mut budget), Err(crate::Error::Cancelled));
        }
        assert_eq!(
            run(&mut WorkBudget::new(&Execution::default()).unwrap()).unwrap(),
            expected
        );
        let execution = Execution::default().max_work(work);
        assert_eq!(
            run(&mut WorkBudget::new(&execution).unwrap()).unwrap(),
            expected
        );
    }
}

#[test]
fn unlimited_specialization_preserves_results_and_mid_operation_cancellation() {
    let n = 64;
    for pattern in 0..3 {
        let value = |a: usize, b: usize| {
            if a == b {
                0.
            } else {
                let (a, b) = if a < b { (a, b) } else { (b, a) };
                match pattern {
                    0 => {
                        if a % 2 != b % 2 {
                            1.
                        } else {
                            2.
                        }
                    }
                    1 => ((a * 17 + b * 11) % 7) as f64 / 4. + 1.,
                    _ => ((a * 17 + b * 11) % 13) as f64 / 4. + 1.,
                }
            }
        };
        let cutoff = if pattern == 0 {
            1.
        } else if pattern == 1 {
            3.
        } else {
            2.5
        };
        for layout in [
            MatrixLayout::LowerTriangle,
            MatrixLayout::UpperTriangle,
            MatrixLayout::Square,
        ] {
            let values: Vec<_> = match layout {
                MatrixLayout::LowerTriangle => (0..n)
                    .flat_map(|b| (0..b).map(move |a| value(a, b)))
                    .collect(),
                MatrixLayout::UpperTriangle => (0..n)
                    .flat_map(|a| (a + 1..n).map(move |b| value(a, b)))
                    .collect(),
                MatrixLayout::Square => (0..n)
                    .flat_map(|a| (0..n).map(move |b| value(a, b)))
                    .collect(),
            };
            let input = DissimilarityMatrixView::new(&values, n, layout).unwrap();
            let unlimited = || {
                compute_dense(
                    input,
                    1,
                    cutoff,
                    PrimeField::default(),
                    &mut WorkBudget::new(&Execution::default()).unwrap(),
                )
                .unwrap()
            };
            let expected = unlimited();
            let execution = Execution::default().max_work(u64::MAX);
            let mut controlled = WorkBudget::new(&execution).unwrap();
            let actual =
                compute_dense(input, 1, cutoff, PrimeField::default(), &mut controlled).unwrap();
            assert_eq!(actual, expected, "pattern={pattern} layout={layout:?}");
            let work = controlled.used();
            assert!(work > 0);
            for interrupt in [0, work / 3, 2 * work / 3, work - 1] {
                let flag = AtomicBool::new(false);
                let execution = Execution::default().max_work(u64::MAX).cancellation(&flag);
                let mut budget = WorkBudget::new(&execution).unwrap();
                budget.cancel_at_work(interrupt);
                assert_eq!(
                    compute_dense(input, 1, cutoff, PrimeField::default(), &mut budget),
                    Err(crate::Error::Cancelled)
                );
                assert_eq!(unlimited(), expected);
            }
        }
    }
}
