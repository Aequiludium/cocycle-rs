use super::*;
use crate::geometry::{DissimilarityView, MatrixLayout};
use crate::persistence::rips::resolve_rips_range;
use crate::persistence::{RipsOptions, assemble_diagram, reference};

/// Replay C V using independent vertex triples and set XOR. Cursor suffixes
/// must reconstruct the *full* reduced column, including cancelled prefixes.
pub(super) fn check_transform(
    data: &OrderedCoboundaryAccess<'_>,
    stop: f64,
    edges: &[SimplexEntry],
    edge: SimplexEntry,
    additions: &[EdgePosition],
    working: &WorkingColumn,
    pivot: SimplexEntry,
) {
    let input = data.input;
    let n = input.len();
    let mut expected = std::collections::BTreeSet::new();
    for source in std::iter::once(edge).chain(additions.iter().map(|position| edges[position.0])) {
        // Enumerate IDs independently of SimplexIndex and production access.
        let mut id = 0;
        for c in 2..n {
            for b in 1..c {
                for a in 0..b {
                    let triangle_edges = [
                        b * (b - 1) / 2 + a,
                        c * (c - 1) / 2 + a,
                        c * (c - 1) / 2 + b,
                    ];
                    let value = input
                        .get(a, b)
                        .unwrap()
                        .max(input.get(a, c).unwrap())
                        .max(input.get(b, c).unwrap());
                    if triangle_edges.contains(&source.id) && value <= stop {
                        let row = SimplexEntry { id, value };
                        if !expected.insert(row) {
                            expected.remove(&row);
                        }
                    }
                    id += 1;
                }
            }
        }
    }
    let mut actual = std::collections::BTreeSet::from([pivot]);
    let mut remaining = working.clone();
    while let Some(row) = remaining.pop(data, &mut WorkBudget::unlimited()).unwrap() {
        assert!(actual.insert(row), "parity pop must emit distinct rows");
    }
    assert_eq!(actual, expected, "ordered suffixes must preserve R = C V");
    assert_eq!(actual.first(), Some(&pivot));
}

#[test]
fn transformations_replay_full_columns_through_stored_and_virtual_owners() {
    let mut total = Stats::default();
    for n in [8, 12] {
        for pattern in 0..3 {
            let lower: Vec<_> = (0..n)
                .flat_map(|b| {
                    (0..b).map(move |a| match pattern {
                        0 => (b - a).min(n - b + a) as f64,
                        1 => ((a * 17 + b * 11) % 7) as f64 / 4.,
                        _ => ((a * 37 + b * 13) % 23) as f64 / 8.,
                    })
                })
                .collect();
            for stop in [1., 2., f64::MAX] {
                for layout in [
                    MatrixLayout::LowerTriangle,
                    MatrixLayout::UpperTriangle,
                    MatrixLayout::Square,
                ] {
                    let values = matrix_values(&lower, n, layout);
                    let input = DissimilarityMatrixView::new(&values, n, layout).unwrap();
                    let edges = DenseFlag::new(input, stop)
                        .unwrap()
                        .edges(&mut || Ok(()))
                        .unwrap();
                    let mut stats = Stats {
                        verify_transforms: true,
                        ..Stats::default()
                    };
                    let raw = reduce(
                        input,
                        stop,
                        edges,
                        &mut WorkBudget::unlimited(),
                        &mut IgnoreTriangles,
                        &mut stats,
                    )
                    .unwrap();
                    let original = DissimilarityView::new(&lower, n).unwrap();
                    let options = RipsOptions::new(1, Some(stop)).unwrap();
                    let (_, coverage) = resolve_rips_range(original, &options);
                    assert_eq!(
                        assemble_diagram(1, coverage, raw).unwrap(),
                        reference::compute(original, &options).unwrap()
                    );
                    total.checked_transforms += stats.checked_transforms;
                    total.stored_additions += stats.stored_additions;
                    total.virtual_additions += stats.virtual_additions;
                    total.omitted_pairs += stats.omitted_pairs;
                }
            }
        }
    }
    assert!(total.checked_transforms > 0);
    assert!(total.stored_additions > 0);
    assert!(total.virtual_additions > 0);
    assert!(total.omitted_pairs > 0);
}

#[test]
fn ordered_triangle_handoff_matches_full_generic_through_h3() {
    use crate::filtration::flag::CliqueAccess;
    use crate::persistence::flag::clearing::TrianglePivots;
    use crate::persistence::simplicial::cohomology;
    for n in 0_usize..=10 {
        for pattern in 0..3 {
            let lower: Vec<_> = (0..n)
                .flat_map(|b| {
                    (0..b).map(move |a| match pattern {
                        0 => 1.,
                        1 => {
                            if a / 2 == b / 2 {
                                2.
                            } else {
                                1.
                            }
                        }
                        _ => ((a * 17 + b * 11) % 7) as f64 / 4.,
                    })
                })
                .collect();
            for stop in [0., 1., 2., f64::MAX] {
                for layout in [
                    MatrixLayout::LowerTriangle,
                    MatrixLayout::UpperTriangle,
                    MatrixLayout::Square,
                ] {
                    let values = matrix_values(&lower, n, layout);
                    let input = DissimilarityMatrixView::new(&values, n, layout).unwrap();
                    let rips = DenseFlag::new(input, stop).unwrap();
                    let edges = rips.edges(&mut || Ok(())).unwrap();
                    let mut pivots = TrianglePivots::default();
                    let mut budget =
                        WorkBudget::new(&crate::execution::Execution::default()).unwrap();
                    let mut actual =
                        compute_with_clearing(input, stop, edges, &mut budget, &mut pivots)
                            .unwrap();
                    let cleared = pivots.into_vertices(n, &mut budget).unwrap();
                    let access = CliqueAccess::Dense(input, stop);
                    actual.extend(
                        cohomology::compute_range(
                            &access,
                            2,
                            3,
                            crate::algebra::PrimeField::default(),
                            cleared,
                            &mut budget,
                        )
                        .unwrap(),
                    );
                    let expected = cohomology::compute(
                        &access,
                        3,
                        crate::algebra::PrimeField::default(),
                        &mut budget,
                    )
                    .unwrap();
                    assert_eq!(
                        assemble_diagram(3, crate::diagram::Coverage::Complete, actual).unwrap(),
                        assemble_diagram(3, crate::diagram::Coverage::Complete, expected).unwrap(),
                        "n={n} pattern={pattern} stop={stop} layout={layout:?}"
                    );
                }
            }
        }
    }
}

fn matrix_values(lower: &[f64], n: usize, layout: MatrixLayout) -> Vec<f64> {
    let get = |a: usize, b: usize| {
        if a == b {
            0.
        } else {
            let (a, b) = if a < b { (a, b) } else { (b, a) };
            lower[b * (b - 1) / 2 + a]
        }
    };
    match layout {
        MatrixLayout::LowerTriangle => lower.to_vec(),
        MatrixLayout::UpperTriangle => (0..n)
            .flat_map(|a| (a + 1..n).map(move |b| get(a, b)))
            .collect(),
        MatrixLayout::Square => (0..n)
            .flat_map(|a| (0..n).map(move |b| get(a, b)))
            .collect(),
    }
}

#[test]
fn ordered_cursors_and_lower_bound_seeks_match_independent_triangle_enumeration() {
    let extreme = [
        0.,
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        1.,
        f64::from_bits(1.0_f64.to_bits() + 1),
        1e300,
        f64::MAX,
    ];
    for n in 0_usize..=10 {
        for pattern in 0..4 {
            let lower: Vec<_> = (0..n * (n.saturating_sub(1)) / 2)
                .map(|i| match pattern {
                    0 => ((i * 7 + i / 3) % 5) as f64 / 2.,
                    1 => 1.,
                    2 => {
                        if i % 3 == 0 {
                            -0.
                        } else {
                            1.
                        }
                    }
                    _ => extreme[(i * 3) % extreme.len()],
                })
                .collect();
            for cutoff in [0., 0.5, 1., 2., f64::MAX] {
                for layout in [
                    MatrixLayout::LowerTriangle,
                    MatrixLayout::UpperTriangle,
                    MatrixLayout::Square,
                ] {
                    let values = matrix_values(&lower, n, layout);
                    let input = DissimilarityMatrixView::new(&values, n, layout).unwrap();
                    let mut budget =
                        WorkBudget::new(&crate::execution::Execution::default()).unwrap();
                    let rips = DenseFlag::new(input, cutoff).unwrap();
                    let data = OrderedCoboundaryAccess::new(input, cutoff, &mut budget).unwrap();
                    for edge in rips.edges(&mut || Ok(())).unwrap() {
                        let [x, y] = rips.edge_vertices(edge.id);
                        let mut expected = Vec::new();
                        let mut id = 0;
                        for c in 2..n {
                            for b in 1..c {
                                for a in 0..b {
                                    let value = input
                                        .get(a, b)
                                        .unwrap()
                                        .max(input.get(a, c).unwrap())
                                        .max(input.get(b, c).unwrap());
                                    if [a, b, c].contains(&x)
                                        && [a, b, c].contains(&y)
                                        && value <= cutoff
                                    {
                                        expected.push(SimplexEntry { id, value });
                                    }
                                    id += 1;
                                }
                            }
                        }
                        expected.sort_unstable();
                        let mut bounds = vec![None];
                        bounds.extend(expected.iter().copied().map(Some));
                        for bound in bounds {
                            let mut cursor = CofacetCursor::new(&data, &rips, edge, bound);
                            let mut actual = Vec::new();
                            while let Some(row) = cursor.next(&data, bound, &mut budget).unwrap() {
                                actual.push(row);
                            }
                            let expected: Vec<_> = expected
                                .iter()
                                .copied()
                                .filter(|&row| bound.is_none_or(|x| row >= x))
                                .collect();
                            assert_eq!(
                                actual, expected,
                                "n={n} pattern={pattern} layout={layout:?} bound={bound:?}"
                            );
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn cursor_reduction_matches_explicit_reference_for_ties_zeros_extremes_and_layouts() {
    for n in 0_usize..=10 {
        for pattern in 0..4 {
            let extreme = [
                0.,
                f64::from_bits(1),
                f64::MIN_POSITIVE,
                1.,
                f64::from_bits(1.0_f64.to_bits() + 1),
                1e300,
                f64::MAX,
            ];
            let lower: Vec<_> = (0..n * (n.saturating_sub(1)) / 2)
                .map(|i| match pattern {
                    0 => ((i * 7 + i / 3) % 5) as f64 / 2.,
                    1 => 1.,
                    2 => {
                        if i % 3 == 0 {
                            -0.
                        } else {
                            1.
                        }
                    }
                    _ => extreme[(i * 3) % extreme.len()],
                })
                .collect();
            let original = DissimilarityView::new(&lower, n).unwrap();
            for cutoff in [0., 0.5, 1., 2., f64::MAX] {
                let options = RipsOptions::new(1, Some(cutoff)).unwrap();
                let (_, coverage) = resolve_rips_range(original, &options);
                let expected = reference::compute(original, &options).unwrap();
                for layout in [
                    MatrixLayout::LowerTriangle,
                    MatrixLayout::UpperTriangle,
                    MatrixLayout::Square,
                ] {
                    let values = matrix_values(&lower, n, layout);
                    let input = DissimilarityMatrixView::new(&values, n, layout).unwrap();
                    let access = DenseFlag::new(input, cutoff).unwrap();
                    let edges = access.edges(&mut || Ok(())).unwrap();
                    for raw in [
                        compute(
                            input,
                            cutoff,
                            edges.clone(),
                            &mut WorkBudget::new(&crate::execution::Execution::default()).unwrap(),
                        )
                        .unwrap(),
                        compute(input, cutoff, edges, &mut WorkBudget::unlimited()).unwrap(),
                    ] {
                        assert_eq!(
                            assemble_diagram(1, coverage, raw).unwrap(),
                            expected,
                            "n={n} pattern={pattern} layout={layout:?}"
                        );
                    }
                }
            }
        }
    }
}

#[test]
fn cursor_reduction_fails_cleanly_on_work_exhaustion_and_cancellation() {
    use std::sync::atomic::{AtomicBool, Ordering};
    let lower: Vec<_> = (0..28).map(|i| ((i * 7) % 13) as f64 / 4.).collect();
    let input = DissimilarityMatrixView::new(&lower, 8, MatrixLayout::LowerTriangle).unwrap();
    let access = DenseFlag::new(input, 4.).unwrap();
    let edges = access.edges(&mut || Ok(())).unwrap();
    for limit in [0, 1, 8, 64] {
        let execution = crate::execution::Execution::default().max_work(limit);
        assert_eq!(
            compute(
                input,
                4.,
                edges.clone(),
                &mut WorkBudget::new(&execution).unwrap()
            ),
            Err(Error::WorkLimitExceeded { limit })
        );
    }
    let cancelled = AtomicBool::new(false);
    let execution = crate::execution::Execution::default().cancellation(&cancelled);
    let mut budget = WorkBudget::new(&execution).unwrap();
    cancelled.store(true, Ordering::Relaxed);
    assert_eq!(
        compute(input, 4., edges.clone(), &mut budget),
        Err(Error::Cancelled)
    );
    assert!(
        compute(
            input,
            4.,
            edges,
            &mut WorkBudget::new(&crate::execution::Execution::default()).unwrap()
        )
        .is_ok()
    );
}
