//! End-to-end public workflows, source contracts and owned outputs.
use cocycle::{
    Error, Result,
    algebra::PrimeField,
    complex::{WeightedEdge, WeightedGraph},
    diagram::IntervalEnd,
    filtration::{ApproximateRipsBuilder, Coverage, FiltrationKind, FlagFiltration, RipsBuilder},
    geometry::{DissimilarityMatrixView, MatrixLayout, MetricPolicy, PointCloudView},
    persistence::{PersistenceExt, RepresentativeRequest, RepresentativeSelection},
};
use std::cell::Cell;

fn matrix(values: &[f64], n: usize) -> DissimilarityMatrixView<'_> {
    DissimilarityMatrixView::new(values, n, MatrixLayout::LowerTriangle).unwrap()
}

#[test]
fn higher_requests_agree_across_points_layouts_prepared_graphs_and_explicit_sources() -> Result<()>
{
    // Six octahedral vertices: a genuine H2 sphere appears at sqrt(2) and
    // dies at 2. Compare all public exact sources, including capped points.
    let coordinates = [
        1., 0., 0., -1., 0., 0., 0., 1., 0., 0., -1., 0., 0., 0., 1., 0., 0., -1.,
    ];
    let points = PointCloudView::new(&coordinates, 6, 3)?;
    let weight = |a: usize, b: usize| {
        if a == b {
            0.
        } else if a / 2 == b / 2 {
            2.
        } else {
            2.0_f64.sqrt()
        }
    };
    for construction in [None, Some(1.5), Some(2.)] {
        let point_source = if let Some(cap) = construction {
            RipsBuilder::from_points(points).max_edge_length(cap)
        } else {
            RipsBuilder::from_points(points)
        };
        let prepared_points = point_source.prepare()?;
        let explicit_points = point_source.build_complex(4)?;
        for cutoff in [None, Some(0.), Some(1.5), Some(2.)] {
            if construction
                .zip(cutoff)
                .is_some_and(|(cap, requested)| requested > cap)
            {
                continue;
            }
            for field in [PrimeField::default(), PrimeField::new(3)?] {
                let request = |q| {
                    let mut request = point_source
                        .persistence()
                        .max_homology_dimension(q)
                        .field(field);
                    if let Some(t) = cutoff {
                        request = request.max_filtration_value(t);
                    }
                    request.compute()
                };
                let low = request(1)?;
                let expected = request(3)?;
                for q in [0, 1] {
                    assert_eq!(
                        low.diagram().dimension(q)?.iter().collect::<Vec<_>>(),
                        expected.diagram().dimension(q)?.iter().collect::<Vec<_>>()
                    );
                }
                let mut prepared_request = prepared_points
                    .persistence()
                    .max_homology_dimension(3)
                    .field(field);
                let mut explicit_request = explicit_points
                    .persistence()
                    .max_homology_dimension(3)
                    .field(field);
                if let Some(t) = cutoff {
                    prepared_request = prepared_request.max_filtration_value(t);
                    explicit_request = explicit_request.max_filtration_value(t);
                }
                assert_eq!(prepared_request.compute()?.diagram(), expected.diagram());
                assert_eq!(explicit_request.compute()?.diagram(), expected.diagram());
                let h2: Vec<_> = expected.diagram().dimension(2)?.iter().collect();
                let stop = construction.unwrap_or(2.).min(cutoff.unwrap_or(2.));
                if stop < 2.0_f64.sqrt() {
                    assert!(h2.is_empty());
                } else {
                    assert_eq!(h2.len(), 1);
                    assert_eq!(h2[0].birth(), 2.0_f64.sqrt());
                    assert_eq!(
                        h2[0].end(),
                        if stop < 2. {
                            IntervalEnd::RightCensored { through: stop }
                        } else {
                            IntervalEnd::Finite(2.)
                        }
                    );
                }
                for layout in [
                    MatrixLayout::LowerTriangle,
                    MatrixLayout::UpperTriangle,
                    MatrixLayout::Square,
                ] {
                    let values: Vec<_> = match layout {
                        MatrixLayout::LowerTriangle => (0..6)
                            .flat_map(|b| (0..b).map(move |a| weight(a, b)))
                            .collect(),
                        MatrixLayout::UpperTriangle => (0..6)
                            .flat_map(|a| (a + 1..6).map(move |b| weight(a, b)))
                            .collect(),
                        MatrixLayout::Square => (0..6)
                            .flat_map(|a| (0..6).map(move |b| weight(a, b)))
                            .collect(),
                    };
                    let mut source = RipsBuilder::from_distance_matrix(
                        DissimilarityMatrixView::new(&values, 6, layout)?,
                    );
                    if let Some(cap) = construction {
                        source = source.max_edge_length(cap);
                    }
                    let mut request = source.persistence().max_homology_dimension(3).field(field);
                    if let Some(t) = cutoff {
                        request = request.max_filtration_value(t);
                    }
                    assert_eq!(request.compute()?.diagram(), expected.diagram());
                }
            }
        }
        // The supplied graph certifies its own completion: absent opposite
        // edges make H2 essential rather than right-censored as in capped Rips.
        let flag = FlagFiltration::new(prepared_points.graph().clone());
        let supplied = flag.persistence().max_homology_dimension(3).compute()?;
        let h2 = supplied.diagram().dimension(2)?.iter().next().unwrap();
        assert_eq!(
            h2.end(),
            if construction == Some(1.5) {
                IntervalEnd::Essential
            } else {
                IntervalEnd::Finite(2.)
            }
        );
    }
    Ok(())
}

#[test]
fn truncated_matrix_preserves_zero_edges_isolates_and_closed_coverage_in_every_layout() -> Result<()>
{
    // A square with one zero edge and sixty isolated vertices below scale two.
    let n = 64;
    let weight = |a: usize, b: usize| {
        let (a, b) = (a.min(b), a.max(b));
        match (a, b) {
            _ if a == b => 0.,
            (0, 1) => -0.,
            (1, 2) | (2, 3) | (0, 3) => 1.,
            _ => 2.,
        }
    };
    let lower: Vec<_> = (0..n)
        .flat_map(|b| (0..b).map(move |a| weight(a, b)))
        .collect();
    let upper: Vec<_> = (0..n)
        .flat_map(|a| (a + 1..n).map(move |b| weight(a, b)))
        .collect();
    let square: Vec<_> = (0..n)
        .flat_map(|a| (0..n).map(move |b| weight(a, b)))
        .collect();
    for (values, layout) in [
        (&lower[..], MatrixLayout::LowerTriangle),
        (&upper[..], MatrixLayout::UpperTriangle),
        (&square[..], MatrixLayout::Square),
    ] {
        let source =
            RipsBuilder::from_distance_matrix(DissimilarityMatrixView::new(values, n, layout)?);
        for cutoff in [0., 1., f64::from_bits(1_f64.to_bits() + 1), 2.] {
            let result = source
                .persistence()
                .max_filtration_value(cutoff)
                .compute()?;
            let diagram = result.diagram();
            let complete = cutoff == 2.;
            assert_eq!(
                diagram.coverage(),
                if complete {
                    Coverage::Complete
                } else {
                    Coverage::Through(cutoff)
                }
            );
            let survivor = if complete {
                IntervalEnd::Essential
            } else {
                IntervalEnd::RightCensored { through: cutoff }
            };
            let h0: Vec<_> = diagram.dimension(0)?.iter().collect();
            assert_eq!(
                h0.iter().filter(|bar| bar.end() == survivor).count(),
                if cutoff == 0. {
                    63
                } else if complete {
                    1
                } else {
                    61
                }
            );
            assert_eq!(
                h0.iter()
                    .filter(|bar| bar.end() == IntervalEnd::Finite(1.))
                    .count(),
                if cutoff == 0. { 0 } else { 2 }
            );
            assert_eq!(
                h0.iter()
                    .filter(|bar| bar.end() == IntervalEnd::Finite(2.))
                    .count(),
                if complete { 60 } else { 0 }
            );
            assert!(h0.iter().all(|bar| bar.birth() == 0.));
            assert_eq!(h0.len(), 63);
            let h1: Vec<_> = diagram.dimension(1)?.iter().collect();
            if cutoff == 0. {
                assert!(h1.is_empty());
            } else {
                assert_eq!(h1.len(), 1);
                assert_eq!(h1[0].birth(), 1.);
                assert_eq!(
                    h1[0].end(),
                    if complete {
                        IntervalEnd::Finite(2.)
                    } else {
                        survivor
                    }
                );
            }
            assert_eq!(result.context().construction_cutoff(), None);
            assert_eq!(result.context().requested_cutoff(), Some(cutoff));
        }
    }
    Ok(())
}

#[test]
fn square_workflows_preserve_fields_bases_context_and_source_ownership() -> Result<()> {
    let result = {
        let values = [1., 2., 1., 1., 2., 1.];
        let rips = RipsBuilder::from_distance_matrix(matrix(&values, 4));
        let filtration = rips.build_complex(2)?;
        assert_eq!(filtration.complex().len(), 14);
        let triangle = filtration.complex().find(&[0, 1, 2]).unwrap();
        assert_eq!(filtration.complex().boundary(triangle).unwrap().len(), 3);
        for characteristic in [2, 3, 251, 65537] {
            let field = PrimeField::new(characteristic)?;
            let requests = [RepresentativeRequest::new(
                1,
                1.,
                RepresentativeSelection::Both,
            )?];
            let direct = rips
                .persistence()
                .field(field)
                .representatives(&requests)
                .compute()?;
            let explicit = filtration
                .persistence()
                .field(field)
                .representatives(&requests)
                .compute()?;
            assert_eq!(direct.diagram(), explicit.diagram());
            assert_eq!(direct.context(), explicit.context());
            assert_eq!(direct.representatives().unwrap().len(), 2);
            assert_eq!(explicit.representatives().unwrap().len(), 2);
            let interval = direct.diagram().dimension(1)?.iter().next().unwrap();
            assert_eq!(
                (interval.birth(), interval.end()),
                (1., IntervalEnd::Finite(2.))
            );
        }
        rips.persistence().compute()?
    };
    assert_eq!(result.diagram().coverage(), Coverage::Complete);
    assert_eq!(result.diagram().dimension(1)?.iter().count(), 1);
    Ok(())
}
#[test]
fn matrix_layouts_and_points_agree_without_losing_construction_caps() -> Result<()> {
    let coordinates = [0., 0., 1., 0., 1., 1., 0., 1.];
    let points = PointCloudView::new(&coordinates, 4, 2)?;
    let diagonal = 2.0_f64.sqrt();
    let lower = [1., diagonal, 1., 1., diagonal, 1.];
    let upper = [1., diagonal, 1., 1., diagonal, 1.];
    let square = [
        0., 1., diagonal, 1., 1., 0., 1., diagonal, diagonal, 1., 0., 1., 1., diagonal, 1., 0.,
    ];
    for cap in [None, Some(1.), Some(2.)] {
        let mut points_builder = RipsBuilder::from_points(points);
        if let Some(t) = cap {
            points_builder = points_builder.max_edge_length(t);
        }
        let expected = points_builder.persistence().compute()?;
        assert_eq!(expected.context().construction_cutoff(), cap);
        for (values, layout) in [
            (&lower[..], MatrixLayout::LowerTriangle),
            (&upper[..], MatrixLayout::UpperTriangle),
            (&square[..], MatrixLayout::Square),
        ] {
            let mut rips =
                RipsBuilder::from_distance_matrix(DissimilarityMatrixView::new(values, 4, layout)?);
            if let Some(t) = cap {
                rips = rips.max_edge_length(t);
            }
            let direct = rips.persistence().compute()?;
            let prepared = rips.prepare()?.persistence().compute()?;
            let explicit = rips.build_complex(2)?.persistence().compute()?;
            assert_eq!(direct.diagram(), expected.diagram());
            assert_eq!(direct.diagram(), prepared.diagram());
            assert_eq!(direct.diagram(), explicit.diagram());
            assert_eq!(direct.context(), explicit.context());
        }
    }
    let rips = RipsBuilder::from_points(points);
    let analysis_only = rips.persistence().max_filtration_value(1.).compute()?;
    assert_eq!(analysis_only.context().construction_cutoff(), None);
    assert_eq!(analysis_only.context().requested_cutoff(), Some(1.));
    Ok(())
}

#[test]
fn independent_point_caps_preserve_source_coverage_and_representatives() -> Result<()> {
    let coordinates = [0., 1., 2., 3.];
    let points = PointCloudView::new(&coordinates, 4, 1)?;
    let distances = matrix(&[1., 2., 1., 3., 2., 1.], 4);
    for construction in [None, Some(0.5), Some(1.), Some(3.), Some(4.)] {
        let mut rips = RipsBuilder::from_points(points);
        let mut dense = RipsBuilder::from_distance_matrix(distances);
        if let Some(cap) = construction {
            rips = rips.max_edge_length(cap);
            dense = dense.max_edge_length(cap);
        }
        let prepared = rips.prepare()?;
        let expanded = rips.build_complex(2)?;
        for analysis in [None, Some(0.), Some(0.5), Some(1.5), Some(3.), Some(5.)] {
            for prime in [2, 3, 65537] {
                let field = PrimeField::new(prime)?;
                let scale = analysis.or(construction).unwrap_or(5.);
                let requests = [RepresentativeRequest::new(
                    0,
                    scale,
                    RepresentativeSelection::Both,
                )?];
                for representatives in [&[][..], &requests[..]] {
                    let mut direct = rips
                        .persistence()
                        .field(field)
                        .representatives(representatives);
                    let mut stored = prepared
                        .persistence()
                        .field(field)
                        .representatives(representatives);
                    let mut explicit = expanded
                        .persistence()
                        .field(field)
                        .representatives(representatives);
                    let mut matrix = dense
                        .persistence()
                        .field(field)
                        .representatives(representatives);
                    if let Some(cap) = analysis {
                        direct = direct.max_filtration_value(cap);
                        stored = stored.max_filtration_value(cap);
                        explicit = explicit.max_filtration_value(cap);
                        matrix = matrix.max_filtration_value(cap);
                    }
                    let results = [
                        direct.compute(),
                        stored.compute(),
                        explicit.compute(),
                        matrix.compute(),
                    ];
                    if construction.is_some_and(|cap| cap < 3. && analysis.is_some_and(|t| t > cap))
                    {
                        for result in results {
                            assert!(matches!(result, Err(Error::IncompleteFiltration { .. })));
                        }
                        continue;
                    }
                    let [direct, stored, explicit, matrix] = results;
                    let direct = direct?;
                    for other in [stored?, explicit?, matrix?] {
                        assert_eq!(direct.diagram(), other.diagram());
                        assert_eq!(
                            direct.context().construction_cutoff(),
                            other.context().construction_cutoff()
                        );
                        assert_eq!(
                            direct.context().requested_cutoff(),
                            other.context().requested_cutoff()
                        );
                        assert_eq!(
                            direct.representatives().map(<[_]>::len),
                            other.representatives().map(<[_]>::len)
                        );
                    }
                    assert_eq!(direct.context().construction_cutoff(), construction);
                    assert_eq!(direct.context().requested_cutoff(), analysis);
                    let effective = match (construction, analysis) {
                        (Some(a), Some(b)) => Some(a.min(b)),
                        (a, b) => a.or(b),
                    };
                    assert_eq!(
                        direct.diagram().coverage(),
                        match effective {
                            Some(t) if t < 3. => Coverage::Through(t),
                            _ => Coverage::Complete,
                        }
                    );
                }
            }
        }
    }
    Ok(())
}

#[test]
fn direct_point_analysis_does_not_process_edges_above_its_cutoff() -> Result<()> {
    let coordinates: Vec<_> = (0..32).map(f64::from).collect();
    let rips =
        RipsBuilder::from_points(PointCloudView::new(&coordinates, 32, 1)?).max_edge_length(100.);
    // Enough for every pair's distance checks and an edgeless H0 computation,
    // but not another pass over all 496 edges allowed by the construction cap.
    let execution = cocycle::execution::Execution::default().max_work(1200);
    for prime in [2, 3, 65537] {
        let result = rips
            .persistence()
            .max_homology_dimension(0)
            .max_filtration_value(0.5)
            .field(PrimeField::new(prime)?)
            .compute_with(&execution)?;
        assert_eq!(result.diagram().len(), 32);
        assert_eq!(result.diagram().coverage(), Coverage::Through(0.5));
        assert_eq!(result.context().construction_cutoff(), Some(100.));
    }
    Ok(())
}
#[test]
fn scale_and_dimension_truncation_are_independent() -> Result<()> {
    let values: Vec<_> = (0..6)
        .flat_map(|b| (0..b).map(move |a| if a / 2 == b / 2 { 2. } else { 1. }))
        .collect();
    let rips = RipsBuilder::from_distance_matrix(matrix(&values, 6));
    let skeleton = rips.build_complex(2)?;
    assert!(matches!(
        skeleton.persistence().max_homology_dimension(2).compute(),
        Err(Error::InsufficientSkeleton { .. })
    ));
    let explicit = rips
        .build_complex(3)?
        .persistence()
        .max_homology_dimension(2)
        .compute()?;
    let direct = rips.persistence().max_homology_dimension(2).compute()?;
    assert_eq!(explicit.diagram(), direct.diagram());
    let h2 = direct.diagram().dimension(2)?.iter().next().unwrap();
    assert_eq!((h2.birth(), h2.end()), (1., IntervalEnd::Finite(2.)));
    let truncated = rips.max_edge_length(1.);
    assert!(matches!(
        truncated.persistence().max_filtration_value(2.).compute(),
        Err(Error::IncompleteFiltration { .. })
    ));
    assert!(matches!(
        truncated
            .prepare()?
            .persistence()
            .max_filtration_value(2.)
            .compute(),
        Err(Error::IncompleteFiltration { .. })
    ));
    assert!(matches!(
        truncated
            .build_complex(3)?
            .persistence()
            .max_filtration_value(2.)
            .compute(),
        Err(Error::IncompleteFiltration { .. })
    ));
    let zero = rips.build_complex(0)?;
    assert_eq!(zero.complex().len(), 6);
    assert!(matches!(
        zero.persistence().max_homology_dimension(0).compute(),
        Err(Error::InsufficientSkeleton { .. })
    ));
    Ok(())
}
#[test]
fn supplied_flag_expansion_preserves_essentiality_and_isolates() -> Result<()> {
    let graph = WeightedGraph::new(
        5,
        [[0, 1], [1, 2], [2, 3], [0, 3]]
            .into_iter()
            .map(|vertices| WeightedEdge {
                vertices,
                value: 1.,
            })
            .collect(),
    )?;
    let flag = FlagFiltration::new(graph);
    let filtration = flag.build_complex(1)?;
    assert!(filtration.is_dimension_complete());
    assert_eq!(
        filtration.context().filtration_kind(),
        FiltrationKind::SuppliedFlag
    );
    assert_eq!(filtration.context().vertex_count(), 5);
    let direct = flag.persistence().compute()?;
    let explicit = filtration.persistence().compute()?;
    assert_eq!(direct, explicit);
    assert_eq!(
        explicit
            .diagram()
            .dimension(0)?
            .iter()
            .filter(|i| i.end() == IntervalEnd::Essential)
            .count(),
        2
    );
    assert_eq!(
        explicit
            .diagram()
            .dimension(1)?
            .iter()
            .next()
            .unwrap()
            .end(),
        IntervalEnd::Essential
    );
    Ok(())
}
#[test]
fn approximation_keeps_blockers_mapping_hypotheses_and_callback_samples() -> Result<()> {
    let points = [0., 0., 1., 2., 4., 7.];
    let input = PointCloudView::new(&points, 6, 1)?;
    for policy in [
        MetricPolicy::Check,
        MetricPolicy::Assume,
        MetricPolicy::Unchecked,
    ] {
        for epsilon in [0.5, 1.5] {
            let approximate = ApproximateRipsBuilder::from_points(input, epsilon, policy)
                .min_insertion_radius(1.)
                .max_filtration_value(4.);
            let direct = approximate
                .persistence()
                .max_homology_dimension(2)
                .compute()?;
            let filtration = approximate.build_complex(3)?;
            let explicit = filtration
                .persistence()
                .max_homology_dimension(2)
                .compute()?;
            assert_eq!(direct, explicit);
            let calls = Cell::new(0);
            let prepared = ApproximateRipsBuilder::from_distance_fn(
                &points,
                |a, b| {
                    calls.set(calls.get() + 1);
                    Ok((a - b).abs())
                },
                epsilon,
                policy,
            )
            .min_insertion_radius(1.)
            .max_filtration_value(4.)
            .prepare()?;
            assert_eq!(calls.get(), 15);
            assert_eq!(
                prepared
                    .persistence()
                    .max_homology_dimension(2)
                    .compute()?
                    .diagram(),
                direct.diagram()
            );
            assert_eq!(
                prepared.approximation(),
                direct.context().approximation().unwrap()
            );
            assert!(prepared.approximation().retained_vertices().len() < points.len());
            assert_eq!(calls.get(), 15);
        }
    }
    Ok(())
}
#[test]
fn empty_duplicate_and_singleton_inputs_are_owned_and_reusable() -> Result<()> {
    for n in [0, 1, 4] {
        let points = vec![0.; n];
        let rips = RipsBuilder::from_points(PointCloudView::new(&points, n, 1)?);
        let explicit = rips.build_complex(usize::MAX)?;
        let expected = usize::from(n > 0);
        assert_eq!(explicit.context().vertex_count(), n);
        assert!(explicit.is_dimension_complete());
        let result = rips.persistence().max_homology_dimension(3).compute()?;
        assert_eq!(result.diagram().len(), expected);
        assert_eq!(
            result,
            explicit.persistence().max_homology_dimension(3).compute()?
        );
        std::thread::scope(|scope| {
            let a = scope.spawn(|| rips.persistence().compute().unwrap());
            let b = scope.spawn(|| rips.persistence().compute().unwrap());
            assert_eq!(a.join().unwrap(), b.join().unwrap());
        });
    }
    Ok(())
}
#[test]
fn setters_replace_values_and_invalid_requests_fail_before_callbacks() -> Result<()> {
    let values = [1.];
    let rips = RipsBuilder::from_distance_matrix(matrix(&values, 2))
        .max_edge_length(f64::NAN)
        .full_range();
    let requests = [RepresentativeRequest::new(
        0,
        0.,
        RepresentativeSelection::Both,
    )?];
    assert!(
        rips.persistence()
            .max_filtration_value(f64::NAN)
            .available_range()
            .representatives(&requests)
            .representatives(&[])
            .compute()?
            .representatives()
            .is_none()
    );
    let calls = Cell::new(0);
    let failure = ApproximateRipsBuilder::from_distance_fn(
        &values,
        |_, _| {
            calls.set(calls.get() + 1);
            Ok(0.)
        },
        0.5,
        MetricPolicy::Check,
    )
    .start_vertex(4)
    .prepare();
    assert!(failure.is_err());
    assert_eq!(calls.get(), 0);
    let failure = RipsBuilder::from_distance_fn(&[0, 1], |_, _| Err(Error::Cancelled)).prepare();
    assert!(matches!(failure, Err(Error::Cancelled)));
    Ok(())
}

#[test]
fn representative_scales_follow_certified_coverage_not_the_requested_cap() -> Result<()> {
    let values = [1.];
    let rips = RipsBuilder::from_distance_matrix(matrix(&values, 2));
    let requests = [RepresentativeRequest::new(
        0,
        3.,
        RepresentativeSelection::Both,
    )?];
    // All changes occurred by 1: a requested cap of 2 still certifies the entire
    // filtration, so the surviving component has a valid representative at 3.
    let result = rips
        .persistence()
        .max_filtration_value(2.)
        .representatives(&requests)
        .compute()?;
    assert_eq!(result.diagram().coverage(), Coverage::Complete);
    assert_eq!(result.representatives().unwrap().len(), 2);
    let explicit = rips
        .build_complex(1)?
        .persistence()
        .max_filtration_value(2.)
        .representatives(&requests)
        .compute()?;
    assert_eq!(result, explicit);
    let points = [0., 1.];
    let input = PointCloudView::new(&points, 2, 1)?;
    let point_result = RipsBuilder::from_points(input)
        .persistence()
        .max_filtration_value(2.)
        .representatives(&requests)
        .compute()?;
    assert_eq!(point_result.diagram(), result.diagram());
    let approximate = ApproximateRipsBuilder::from_points(input, 0.5, MetricPolicy::Check);
    assert_eq!(
        approximate
            .persistence()
            .max_filtration_value(2.)
            .representatives(&requests)
            .compute()?
            .representatives()
            .unwrap()
            .len(),
        2
    );
    assert!(matches!(
        rips.persistence()
            .max_filtration_value(0.5)
            .representatives(&requests)
            .compute(),
        Err(Error::QueryOutsideCoverage { .. })
    ));
    Ok(())
}

#[test]
fn expanded_negative_cutoff_does_not_create_zero_born_components() -> Result<()> {
    let values = [1.];
    let expanded = RipsBuilder::from_distance_matrix(matrix(&values, 2)).build_complex(1)?;
    let requests = [RepresentativeRequest::new(
        0,
        -1.,
        RepresentativeSelection::Both,
    )?];
    let early = expanded.persistence().max_filtration_value(-1.).compute()?;
    assert!(early.diagram().is_empty());
    assert_eq!(early.diagram().coverage(), Coverage::Through(-1.));
    let bases = expanded
        .persistence()
        .max_filtration_value(-1.)
        .representatives(&requests)
        .compute()?;
    assert_eq!(early.diagram(), bases.diagram());
    assert!(bases.representatives().unwrap().is_empty());
    Ok(())
}
