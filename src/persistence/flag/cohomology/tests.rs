use super::*;
use crate::diagram::{Coverage, IntervalEnd};
use crate::persistence::reference::FilteredBoundary;
use crate::persistence::rips::resolve_rips_range;
use crate::persistence::{RipsOptions, assemble_diagram, reference};

// H2 expectations use explicit vertex subsets and a forward F2 boundary
// matrix. No production clique visitor, tuple ordering or parity heap is used.
fn explicit_h2(
    n: usize,
    edge: impl Fn(usize, usize) -> Option<f64>,
    cutoff: f64,
    coverage: Coverage,
) -> crate::diagram::PersistenceDiagram {
    use std::collections::{BTreeSet, HashMap};
    let mut cells = Vec::new();
    for mask in 1_usize..1 << n {
        if mask.count_ones() > 4 {
            continue;
        }
        let mut value = 0_f64;
        let mut present = true;
        for b in 0..n {
            for a in 0..b {
                if mask & (1 << a) != 0 && mask & (1 << b) != 0 {
                    if let Some(w) = edge(a, b) {
                        value = value.max(w);
                    } else {
                        present = false;
                    }
                }
            }
        }
        if present && value <= cutoff {
            cells.push((mask, value));
        }
    }
    cells.sort_by(|a, b| {
        a.1.total_cmp(&b.1)
            .then(a.0.count_ones().cmp(&b.0.count_ones()))
            .then(b.0.cmp(&a.0))
    });
    let positions: HashMap<_, _> = cells.iter().enumerate().map(|(i, c)| (c.0, i)).collect();
    let mut reduced: Vec<BTreeSet<usize>> = Vec::new();
    let mut owners = HashMap::new();
    let mut births = BTreeSet::new();
    let mut raw = Vec::new();
    for (j, &(mask, value)) in cells.iter().enumerate() {
        let mut column = BTreeSet::new();
        if mask.count_ones() > 1 {
            for v in 0..n {
                if mask & (1 << v) != 0 {
                    column.insert(positions[&(mask ^ (1 << v))]);
                }
            }
        }
        while let Some(&pivot) = column.last() {
            if let Some(&owner) = owners.get(&pivot) {
                column = column
                    .symmetric_difference(&reduced[owner])
                    .copied()
                    .collect();
            } else {
                owners.insert(pivot, j);
                births.remove(&pivot);
                raw.push((
                    cells[pivot].0.count_ones() as usize - 1,
                    cells[pivot].1,
                    Some(value),
                ));
                break;
            }
        }
        if column.is_empty() {
            births.insert(j);
        }
        reduced.push(column);
    }
    raw.extend(
        births
            .into_iter()
            .map(|i| (cells[i].0.count_ones() as usize - 1, cells[i].1, None)),
    );
    assemble_diagram(2, coverage, raw).unwrap()
}

fn independent_tetrahedra(
    access: &crate::filtration::flag::CliqueAccess<'_>,
    triangle: crate::filtration::flag::TupleEntry<3>,
) -> Vec<crate::filtration::flag::TupleEntry<4>> {
    use crate::filtration::flag::{CliqueAccess, TupleEntry};
    let (n, cutoff) = match access {
        CliqueAccess::Dense(m, c) => (m.len(), *c),
        CliqueAccess::Sparse(g, c) => (g.vertex_count(), *c),
    };
    let mut result = Vec::new();
    for v in 0..n {
        if triangle.vertices.contains(&v) {
            continue;
        }
        let mut vertices = Vec::from(triangle.vertices);
        vertices.push(v);
        vertices.sort_unstable();
        let mut value = 0_f64;
        let mut present = true;
        for b in 1..4 {
            for a in 0..b {
                let edge = match access {
                    CliqueAccess::Dense(m, _) => m.get(vertices[a], vertices[b]),
                    CliqueAccess::Sparse(g, _) => g.edge_value(vertices[a], vertices[b]),
                };
                if let Some(w) = edge {
                    value = value.max(w);
                } else {
                    present = false;
                }
            }
        }
        if present && value <= cutoff {
            result.push(TupleEntry {
                vertices: vertices.try_into().unwrap(),
                value,
            });
        }
    }
    result
}

pub(super) fn check_h2_transform(
    access: &crate::filtration::flag::CliqueAccess<'_>,
    level: &[crate::filtration::flag::TupleEntry<3>],
    j: usize,
    column: &[usize],
    working: &BinaryHeap<Reverse<crate::filtration::flag::TupleEntry<4>>>,
    pivot: crate::filtration::flag::TupleEntry<4>,
) {
    use std::collections::BTreeSet;
    assert_eq!(column.iter().filter(|&&k| k == j).count(), 1);
    assert!(column.iter().all(|&k| k >= j));
    assert!(column.windows(2).all(|w| w[0] > w[1]));
    let mut expected = BTreeSet::new();
    for &k in column {
        for row in independent_tetrahedra(access, level[k]) {
            if !expected.insert(row) {
                expected.remove(&row);
            }
        }
    }
    let mut actual = BTreeSet::new();
    for row in working.iter().map(|r| r.0).chain([pivot]) {
        if !actual.insert(row) {
            actual.remove(&row);
        }
    }
    assert_eq!(
        actual, expected,
        "R = C V from independent all-edge enumeration"
    );
    assert_eq!(expected.first(), Some(&pivot));
}

fn check_h2_graph(
    graph: &crate::complex::WeightedGraph,
    cutoff: f64,
    coverage: Coverage,
) -> H2Stats {
    use crate::filtration::flag::{CliqueAccess, SparseFlag};
    let access = CliqueAccess::Sparse(graph, cutoff);
    let mut stats = H2Stats {
        verify_transforms: true,
        ..H2Stats::default()
    };
    let execution = crate::execution::Execution::default();
    let raw = run_h2(
        &SparseFlag::new(graph, cutoff).unwrap(),
        &access,
        &mut WorkBudget::new(&execution).unwrap(),
        &mut stats,
    )
    .unwrap();
    let diagram = assemble_diagram(2, coverage, raw).unwrap();
    assert_eq!(
        diagram,
        explicit_h2(
            graph.vertex_count(),
            |a, b| graph.edge_value(a, b),
            cutoff,
            coverage
        )
    );
    let generic = crate::persistence::simplicial::cohomology::compute(
        &access,
        2,
        crate::algebra::PrimeField::new(2).unwrap(),
        &mut WorkBudget::new(&execution).unwrap(),
    )
    .unwrap();
    assert_eq!(diagram, assemble_diagram(2, coverage, generic).unwrap());
    stats
}

#[test]
fn h2_fixed_tuple_order_and_all_cofacets_match_explicit_vertices() {
    use crate::filtration::flag::{CliqueAccess, TupleEntry};
    let n = 9;
    let values: Vec<_> = (0..36).map(|i| (i % 5) as f64).collect();
    let input = DissimilarityView::new(&values, n).unwrap();
    let graph = crate::complex::WeightedGraph::new(
        n,
        (0..n)
            .flat_map(|b| (0..b).map(move |a| (a, b)))
            .filter(|&(a, b)| (a + b) % 3 != 0)
            .map(|(a, b)| crate::complex::WeightedEdge {
                vertices: [a, b],
                value: input.get(a, b).unwrap(),
            })
            .collect(),
    )
    .unwrap();
    for cutoff in [0., 2., 4.] {
        for access in [
            CliqueAccess::Dense(input.into(), cutoff),
            CliqueAccess::Sparse(&graph, cutoff),
        ] {
            for c in 2..n {
                for b in 1..c {
                    for a in 0..b {
                        let value = match &access {
                            CliqueAccess::Dense(m, _) => Some(
                                m.get(a, b)
                                    .unwrap()
                                    .max(m.get(a, c).unwrap())
                                    .max(m.get(b, c).unwrap()),
                            ),
                            CliqueAccess::Sparse(g, _) => g
                                .edge_value(a, b)
                                .zip(g.edge_value(a, c))
                                .zip(g.edge_value(b, c))
                                .map(|((x, y), z)| x.max(y).max(z)),
                        };
                        let Some(value) = value.filter(|&v| v <= cutoff) else {
                            continue;
                        };
                        let triangle = TupleEntry {
                            vertices: [a, b, c],
                            value,
                        };
                        let mut actual = Vec::new();
                        access
                            .visit_tetrahedra(triangle, &mut || Ok(()), |row| {
                                actual.push(row);
                                Ok(())
                            })
                            .unwrap();
                        actual.sort_unstable();
                        let mut expected = independent_tetrahedra(&access, triangle);
                        expected.sort_unstable();
                        assert_eq!(actual, expected);
                        for row in &actual {
                            for other in &actual {
                                assert_eq!(
                                    row.cmp(other),
                                    crate::complex::Simplex::new(row.vertices.to_vec(), row.value)
                                        .unwrap()
                                        .cmp(
                                            &crate::complex::Simplex::new(
                                                other.vertices.to_vec(),
                                                other.value
                                            )
                                            .unwrap()
                                        )
                                );
                            }
                        }
                    }
                }
            }
        }
    }
}

#[test]
fn h2_prototype_matches_explicit_and_generic_dense_sparse_ties_cutoffs() {
    use crate::filtration::flag::CliqueAccess;
    let mut seed = 20261004_u64;
    let mut additions = 0;
    let mut transforms = 0;
    for n in 0_usize..=10 {
        for sample in 0..20 {
            let values: Vec<_> = (0..n * n.saturating_sub(1) / 2)
                .map(|_| {
                    seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                    if sample % 2 == 0 {
                        ((seed >> 32) % 5) as f64
                    } else {
                        (seed >> 11) as f64 / (1_u64 << 50) as f64
                    }
                })
                .collect();
            let input = DissimilarityView::new(&values, n).unwrap();
            for cutoff in [0., 1., 2., 4.] {
                let coverage = Coverage::Through(cutoff);
                let execution = crate::execution::Execution::default();
                let access = CliqueAccess::Dense(input.into(), cutoff);
                let mut stats = H2Stats {
                    verify_transforms: true,
                    ..H2Stats::default()
                };
                let raw = run_h2(
                    &DenseFlag::new(input.into(), cutoff).unwrap(),
                    &access,
                    &mut WorkBudget::new(&execution).unwrap(),
                    &mut stats,
                )
                .unwrap();
                let diagram = assemble_diagram(2, coverage, raw).unwrap();
                assert_eq!(
                    diagram,
                    explicit_h2(n, |a, b| input.get(a, b), cutoff, coverage)
                );
                let generic = crate::persistence::simplicial::cohomology::compute(
                    &access,
                    2,
                    crate::algebra::PrimeField::new(2).unwrap(),
                    &mut WorkBudget::new(&execution).unwrap(),
                )
                .unwrap();
                assert_eq!(diagram, assemble_diagram(2, coverage, generic).unwrap());
                additions += stats.additions;
                transforms += stats.checked_transforms;
                let graph = crate::complex::WeightedGraph::new(
                    n,
                    (0..n)
                        .flat_map(|b| (0..b).map(move |a| (a, b)))
                        .filter(|&(a, b)| sample % 3 == 0 || (a + b + sample) % 4 != 0)
                        .map(|(a, b)| crate::complex::WeightedEdge {
                            vertices: [a, b],
                            value: input.get(a, b).unwrap(),
                        })
                        .collect(),
                )
                .unwrap();
                let sparse = check_h2_graph(&graph, cutoff, coverage);
                additions += sparse.additions;
                transforms += sparse.checked_transforms;
            }
        }
    }
    assert!(additions > 0 && transforms > 0);
}

#[test]
fn h2_octahedra_preserve_finite_essential_censored_and_repeated_intervals() {
    use crate::complex::{WeightedEdge, WeightedGraph};
    use crate::filtration::flag::CliqueAccess;
    let n = 12;
    let values: Vec<_> = (0..n)
        .flat_map(|b| {
            (0..b).map(move |a| {
                if a / 6 != b / 6 {
                    3.
                } else if a / 2 == b / 2 {
                    2.
                } else {
                    1.
                }
            })
        })
        .collect();
    let input = DissimilarityView::new(&values, n).unwrap();
    for cutoff in [1., 1.5, 2., 3.] {
        let coverage = if cutoff == 3. {
            Coverage::Complete
        } else {
            Coverage::Through(cutoff)
        };
        let raw = compute_h2(
            &DenseFlag::new(input.into(), cutoff).unwrap(),
            &CliqueAccess::Dense(input.into(), cutoff),
            &mut WorkBudget::new(&crate::execution::Execution::default()).unwrap(),
        )
        .unwrap();
        let diagram = assemble_diagram(2, coverage, raw).unwrap();
        assert_eq!(
            diagram,
            explicit_h2(n, |a, b| input.get(a, b), cutoff, coverage)
        );
        let bars: Vec<_> = diagram.dimension(2).unwrap().iter().collect();
        assert_eq!(bars.len(), 2);
        let end = if cutoff < 2. {
            IntervalEnd::RightCensored { through: cutoff }
        } else {
            IntervalEnd::Finite(2.)
        };
        assert!(bars.iter().all(|bar| bar.birth() == 1. && bar.end() == end));
    }
    let graph = WeightedGraph::new(
        n,
        (0..n)
            .flat_map(|b| (0..b).map(move |a| (a, b)))
            .filter(|&(a, b)| a / 6 == b / 6 && a / 2 != b / 2)
            .map(|(a, b)| WeightedEdge {
                vertices: [a, b],
                value: 1.,
            })
            .collect(),
    )
    .unwrap();
    check_h2_graph(&graph, f64::INFINITY, Coverage::Complete);
    let expected = explicit_h2(
        n,
        |a, b| graph.edge_value(a, b),
        f64::INFINITY,
        Coverage::Complete,
    );
    assert_eq!(expected.dimension(2).unwrap().len(), 2);
    assert!(
        expected
            .dimension(2)
            .unwrap()
            .iter()
            .all(|b| b.end() == IntervalEnd::Essential)
    );
}

#[test]
fn h2_parity_heap_cancels_multiplicity_and_shared_budget_interrupts_every_boundary() {
    let execution = crate::execution::Execution::default().max_work(u64::MAX);
    let mut budget = WorkBudget::new(&execution).unwrap();
    let mut heap = BinaryHeap::from([3, 3, 3, 2, 2, 1, 1, 1, 1, 0]);
    assert_eq!(pop_parity(&mut heap, &mut budget).unwrap(), Some(3));
    assert_eq!(pop_parity(&mut heap, &mut budget).unwrap(), Some(0));
    assert_eq!(pop_parity(&mut heap, &mut budget).unwrap(), None);
    let values = [1.; 15];
    let input = DissimilarityView::new(&values, 6).unwrap();
    let rips = DenseFlag::new(input.into(), 1.).unwrap();
    let access = crate::filtration::flag::CliqueAccess::Dense(input.into(), 1.);
    compute_h2(&rips, &access, &mut budget).unwrap();
    let total = budget.used();
    for limit in 0..total - 10 {
        // includes handoff, level construction and reduction
        let flag = std::sync::atomic::AtomicBool::new(false);
        let limits = execution.max_work(limit);
        assert!(matches!(
            compute_h2(&rips, &access, &mut WorkBudget::new(&limits).unwrap()),
            Err(Error::WorkLimitExceeded { .. })
        ));
        let limits = execution.cancellation(&flag);
        let mut budget = WorkBudget::new(&limits).unwrap();
        budget.cancel_at_work(limit);
        assert_eq!(
            compute_h2(&rips, &access, &mut budget),
            Err(Error::Cancelled)
        );
    }
    assert!(compute_h2(&rips, &access, &mut WorkBudget::new(&execution).unwrap()).is_ok());
}

#[test]
fn h2_extreme_and_adjacent_f64_endpoints_remain_exact() {
    use crate::filtration::flag::CliqueAccess;
    for birth in [f64::from_bits(1), 1., 1.0e300] {
        let death = f64::from_bits(birth.to_bits() + 1);
        let values: Vec<_> = (0..6)
            .flat_map(|b| (0..b).map(move |a| if a / 2 == b / 2 { death } else { birth }))
            .collect();
        let input = DissimilarityView::new(&values, 6).unwrap();
        for cutoff in [birth, death] {
            let coverage = Coverage::Through(cutoff);
            let raw = compute_h2(
                &DenseFlag::new(input.into(), cutoff).unwrap(),
                &CliqueAccess::Dense(input.into(), cutoff),
                &mut WorkBudget::new(&crate::execution::Execution::default()).unwrap(),
            )
            .unwrap();
            let diagram = assemble_diagram(2, coverage, raw).unwrap();
            assert_eq!(
                diagram,
                explicit_h2(6, |a, b| input.get(a, b), cutoff, coverage)
            );
            let end = if cutoff == birth {
                IntervalEnd::RightCensored { through: birth }
            } else {
                IntervalEnd::Finite(death)
            };
            assert!(
                diagram
                    .dimension(2)
                    .unwrap()
                    .iter()
                    .any(|bar| bar.birth() == birth && bar.end() == end)
            );
        }
    }
}

fn compare_all(values: &[f64], n: usize, cutoff: Option<f64>) {
    let input = DissimilarityView::new(values, n).unwrap();
    let options = RipsOptions::new(1, cutoff).unwrap();
    let expected = reference::compute(input, &options).unwrap();
    let (cutoff, coverage) = resolve_rips_range(input, &options);
    macro_rules! check {
        ($implicit:literal, $clear:literal, $cone:literal, $short:ident) => {
            check!($implicit, $clear, $cone, $short, false)
        };
        ($implicit:literal, $clear:literal, $cone:literal, $short:ident, $two_pass:literal) => {{
            let mut stats = Stats {
                two_pass_initialization: $two_pass,
                ..Stats::default()
            };
            let raw = run::<$implicit, $clear, $cone, $short>(input, cutoff, &mut stats).unwrap();
            assert_eq!(
                assemble_diagram(1, coverage, raw).unwrap(),
                expected,
                "n={n} cutoff={cutoff} implicit={} clear={} cone={} shortcuts={} two_pass={} values={values:?}",
                $implicit,
                $clear,
                $cone,
                $short,
                $two_pass
            );
        }};
    }
    check!(false, false, false, NO_SHORTCUTS);
    check!(false, true, false, NO_SHORTCUTS);
    check!(true, true, false, NO_SHORTCUTS);
    check!(true, true, true, NO_SHORTCUTS);
    check!(true, true, true, APPARENT);
    check!(true, true, true, EMERGENT);
    check!(true, true, true, APPARENT_EMERGENT);
    check!(true, true, true, APPARENT_EMERGENT, true);
    check!(true, true, true, VIRTUAL_APPARENT);
    check!(true, true, true, ALL_SHORTCUTS);
    check!(true, true, true, ALL_SHORTCUTS, true);
    check!(false, true, false, ALL_SHORTCUTS);
}

#[test]
fn every_three_level_four_vertex_filtration_matches_reference_with_each_optimization() {
    for code in 0..3_usize.pow(6) {
        let mut remaining = code;
        let values: Vec<_> = (0..6)
            .map(|_| {
                let value = (remaining % 3) as f64;
                remaining /= 3;
                value
            })
            .collect();
        for cutoff in [None, Some(0.0), Some(1.0), Some(1.5)] {
            compare_all(&values, 4, cutoff);
        }
    }
}

#[test]
fn randomized_f64_nonmetric_filtrations_and_ties_match_each_optimization() {
    let mut state = 173_u64;
    for n in 0_usize..=12 {
        for sample in 0..24 {
            let values: Vec<_> = (0..n * n.saturating_sub(1) / 2)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    if sample % 2 == 0 {
                        ((state >> 32) % 9) as f64 / 4.0
                    } else {
                        ((state >> 11) as f64) / ((1_u64 << 52) as f64)
                    }
                })
                .collect();
            for cutoff in [None, Some(0.0), Some(0.75), Some(1.0), Some(1.8)] {
                compare_all(&values, n, cutoff);
            }
        }
    }
}

#[test]
fn extreme_scales_and_adjacent_f64_values_do_not_quantize_or_overflow() {
    for base in [f64::from_bits(1), 1.0, 1.0e300] {
        let next = f64::from_bits(base.to_bits() + 1);
        let values = [base, next, base, base, next, base];
        for cutoff in [None, Some(base), Some(next)] {
            compare_all(&values, 4, cutoff);
        }
    }
    compare_all(&[-0.0, f64::MAX, 0.0, f64::MAX, f64::MAX, -0.0], 4, None);
}

#[test]
fn cone_stopping_preserves_user_coverage_and_closed_boundary() {
    // Vertex zero is a cone point at 1, while the full input diameter is 5.
    let values = [1., 1., 5., 1., 5., 5.];
    let input = DissimilarityView::new(&values, 4).unwrap();
    assert_eq!(cone_radius(input.into(), &mut || Ok(())).unwrap(), 1.0);
    for cutoff in [0.5, 1.0, 2.0, 5.0] {
        compare_all(&values, 4, Some(cutoff));
        let options = RipsOptions::new(1, Some(cutoff)).unwrap();
        let result = crate::persistence::rips_from_dissimilarities(input, &options).unwrap();
        if (1.0..5.0).contains(&cutoff) {
            assert_eq!(result.coverage(), Coverage::Through(cutoff));
            assert_eq!(result.len(), 4);
            assert!(
                result
                    .intervals()
                    .any(|bar| bar.end() == IntervalEnd::RightCensored { through: cutoff })
            );
        }
    }
}

#[test]
fn large_cycles_and_disconnected_components_survive_truncation() {
    for n in [16, 24] {
        let values: Vec<_> = (0..n)
            .flat_map(|b| {
                (0..b).map(move |a| {
                    // Two disjoint cycle graphs until scale 3; a complete graph at 4.
                    let half = n / 2;
                    if a / half != b / half {
                        4.0
                    } else if b - a == 1 || b - a == half - 1 {
                        1.0
                    } else {
                        3.0
                    }
                })
            })
            .collect();
        for cutoff in [None, Some(0.0), Some(1.0), Some(2.0), Some(3.0)] {
            compare_all(&values, n, cutoff);
        }
    }
}

#[test]
fn heap_entries_cancel_by_parity_instead_of_set_deduplication() {
    let mut heap = BinaryHeap::from(vec![5, 5, 4, 4, 4, 2, 2]);
    assert_eq!(
        pop_parity(
            &mut heap,
            &mut WorkBudget::new(&crate::persistence::ExecutionLimits::default()).unwrap()
        )
        .unwrap(),
        Some(4)
    );
    assert_eq!(
        pop_parity(
            &mut heap,
            &mut WorkBudget::new(&crate::persistence::ExecutionLimits::default()).unwrap()
        )
        .unwrap(),
        None
    );
}

#[test]
fn original_column_initialization_reuses_storage_and_checks_only_the_first_equal_cofacet() {
    // Edge 01 has cofacets 014 at value 2, then 013 and 012 at value 1.
    let values = [1., 1., 1., 1., 1., 2., 2., 2., 2., 2.];
    let input = DissimilarityView::new(&values, 5).unwrap();
    let access = DenseFlag::new(input.into(), 2.).unwrap();
    let edge = SimplexEntry { id: 0, value: 1. };
    let first_equal = SimplexEntry { id: 1, value: 1. };
    let mut working = Coboundary::with_capacity(16);
    working.push(Reverse(SimplexEntry { id: 99, value: 9. }.into()));
    let capacity = working.capacity();
    let mut stats = Stats::default();
    let mut budget = WorkBudget::new(&crate::persistence::ExecutionLimits::default()).unwrap();
    assert_eq!(
        initialize_coboundary::<APPARENT_EMERGENT, true>(
            &access,
            edge,
            &HashMap::new(),
            &mut working,
            &mut stats,
            &mut budget,
        )
        .unwrap(),
        (Some(first_equal), false)
    );
    assert!(working.is_empty());
    assert_eq!(working.capacity(), capacity);
    assert_eq!(stats.initial_candidates, 2);
    assert_eq!(stats.cofacets, 2);

    // The next equal cofacet is unowned, but cannot replace an occupied pivot.
    let owners = HashMap::from([(first_equal.id, ColumnPosition(0))]);
    let mut stats = Stats::default();
    assert_eq!(
        initialize_coboundary::<APPARENT_EMERGENT, true>(
            &access,
            edge,
            &owners,
            &mut working,
            &mut stats,
            &mut budget,
        )
        .unwrap(),
        (None, false)
    );
    assert_eq!(stats.initial_candidates, 5);
    assert_eq!(stats.cofacets, 3);
    assert_eq!(working.pop(), Some(Reverse(first_equal.into())));
    assert_eq!(
        working.pop(),
        Some(Reverse(SimplexEntry { id: 0, value: 1. }.into()))
    );
    assert_eq!(
        working.pop(),
        Some(Reverse(SimplexEntry { id: 4, value: 2. }.into()))
    );
    assert!(working.is_empty());
}

#[test]
fn original_column_fallback_handles_empty_no_equal_and_apparent_only_rejection() {
    for (values, cutoff, edge, expected_rows) in [
        (vec![1., 2., 2.], 1., SimplexEntry { id: 0, value: 1. }, 0),
        (vec![1., 2., 2.], 2., SimplexEntry { id: 0, value: 1. }, 1),
        (vec![1., 1., 1.], 1., SimplexEntry { id: 1, value: 1. }, 1),
    ] {
        let input = DissimilarityView::new(&values, 3).unwrap();
        let access = DenseFlag::new(input.into(), cutoff).unwrap();
        let mut working = Coboundary::new();
        let mut stats = Stats::default();
        let mut budget = WorkBudget::new(&crate::persistence::ExecutionLimits::default()).unwrap();
        assert_eq!(
            initialize_coboundary::<APPARENT, true>(
                &access,
                edge,
                &HashMap::new(),
                &mut working,
                &mut stats,
                &mut budget,
            )
            .unwrap(),
            (None, false)
        );
        assert_eq!(working.len(), expected_rows);
        assert_eq!(stats.initial_candidates, 3);
    }
}

#[test]
fn single_pass_visits_failed_prefixes_once_and_reuses_alternating_buffers() {
    let mut working = Coboundary::with_capacity(128);
    let capacity = working.capacity();
    let mut large = vec![2.; 64 * 63 / 2];
    large[0] = 1.;
    // Alternate a large column, a rejected and a successful shortcut after a
    // high prefix, an empty column, and a column without an equal cofacet.
    for (n, values, cutoff, occupied, visits) in [
        (64, large, 2., false, [128, 64]),
        (
            5,
            vec![1., 1., 1., 1., 1., 2., 2., 2., 2., 2.],
            2.,
            true,
            [7, 5],
        ),
        (
            5,
            vec![1., 1., 1., 1., 1., 2., 2., 2., 2., 2.],
            2.,
            false,
            [2, 2],
        ),
        (3, vec![1., 2., 2.], 1., false, [6, 3]),
        (3, vec![1., 2., 2.], 2., false, [6, 3]),
    ] {
        let input = DissimilarityView::new(&values, n).unwrap();
        let access = DenseFlag::new(input.into(), cutoff).unwrap();
        let mut owners = HashMap::new();
        if occupied {
            owners.insert(1, ColumnPosition(0));
        }
        let mut outcomes = Vec::new();
        for (two_pass, expected_visits) in [true, false].into_iter().zip(visits) {
            let mut stats = Stats {
                two_pass_initialization: two_pass,
                ..Stats::default()
            };
            let mut budget =
                WorkBudget::new(&crate::persistence::ExecutionLimits::default()).unwrap();
            let result = initialize_coboundary::<APPARENT_EMERGENT, true>(
                &access,
                SimplexEntry { id: 0, value: 1. },
                &owners,
                &mut working,
                &mut stats,
                &mut budget,
            )
            .unwrap();
            assert_eq!(stats.initial_candidates, expected_visits);
            assert_eq!(working.capacity(), capacity);
            outcomes.push((result, working.clone().into_sorted_vec()));
        }
        assert_eq!(outcomes[0], outcomes[1]);
    }
}

#[test]
fn initial_scan_work_failure_does_not_poison_reused_input_or_heap() {
    let values = [1., 2., 2.];
    let input = DissimilarityView::new(&values, 3).unwrap();
    let access = DenseFlag::new(input.into(), 2.).unwrap();
    let edge = SimplexEntry { id: 0, value: 1. };
    let mut working = Coboundary::new();
    let mut budget =
        WorkBudget::new(&crate::persistence::ExecutionLimits::new(Some(1), None)).unwrap();
    assert_eq!(
        initialize_coboundary::<APPARENT_EMERGENT, true>(
            &access,
            edge,
            &HashMap::new(),
            &mut working,
            &mut Stats::default(),
            &mut budget,
        )
        .unwrap_err(),
        Error::WorkLimitExceeded { limit: 1 }
    );
    let mut budget = WorkBudget::new(&crate::persistence::ExecutionLimits::default()).unwrap();
    assert_eq!(
        initialize_coboundary::<APPARENT_EMERGENT, true>(
            &access,
            edge,
            &HashMap::new(),
            &mut working,
            &mut Stats::default(),
            &mut budget,
        )
        .unwrap(),
        (None, false)
    );
    assert_eq!(
        working.pop(),
        Some(Reverse(SimplexEntry { id: 0, value: 2. }.into()))
    );
}

#[test]
fn raw_h1_output_omits_zero_bars_without_omitting_repeated_positive_bars() {
    let equal = [1.; 15];
    let input = DissimilarityView::new(&equal, 6).unwrap();
    let mut stats = Stats::default();
    let raw = run::<true, true, false, APPARENT_EMERGENT>(input, 1., &mut stats).unwrap();
    assert!(stats.shortcuts > 0);
    assert!(raw.iter().all(|&(dimension, _, _)| dimension == 0));

    // K(3,3) has four independent cycles at 1, all dying at 2.
    let values: Vec<_> = (0..6)
        .flat_map(|b| (0..b).map(move |a| if a / 3 == b / 3 { 2. } else { 1. }))
        .collect();
    let input = DissimilarityView::new(&values, 6).unwrap();
    let raw =
        run::<true, true, false, APPARENT_EMERGENT>(input, 2., &mut Stats::default()).unwrap();
    assert_eq!(
        raw.iter().filter(|&&bar| bar == (1, 1., Some(2.))).count(),
        4
    );
    assert!(
        raw.iter()
            .all(|&(dimension, birth, death)| dimension == 0 || death != Some(birth))
    );
}

fn check_virtual_access(
    access: &impl FlagAccess,
    coverage: Coverage,
    expected: &crate::diagram::PersistenceDiagram,
) -> Stats {
    let mut single_pass = Stats::default();
    for two_pass in [true, false] {
        let mut ordinary = Stats {
            two_pass_initialization: two_pass,
            verify_transforms: true,
            ..Stats::default()
        };
        let mut optimized = Stats {
            two_pass_initialization: two_pass,
            verify_transforms: true,
            ..Stats::default()
        };
        for (shortcuts, stats) in [
            (APPARENT_EMERGENT, &mut ordinary),
            (ALL_SHORTCUTS, &mut optimized),
        ] {
            let mut budget =
                WorkBudget::new(&crate::persistence::ExecutionLimits::default()).unwrap();
            let raw = if shortcuts == APPARENT_EMERGENT {
                run_access::<true, true, APPARENT_EMERGENT, true>(access, stats, &mut budget)
            } else {
                run_access::<true, true, ALL_SHORTCUTS, true>(access, stats, &mut budget)
            }
            .unwrap();
            assert_eq!(&assemble_diagram(1, coverage, raw).unwrap(), expected);
        }
        assert!(optimized.stored_columns <= ordinary.stored_columns);
        if optimized.skipped_apparent > 0 {
            assert!(optimized.stored_columns < ordinary.stored_columns);
        }
        single_pass = optimized;
    }
    single_pass
}

/// Check R = C V with independent set XOR rather than the production heap.
/// This catches missing virtual terms even when re-eliminating their pivots
/// happens to leave the final interval multiset unchanged.
pub(super) fn check_transform(
    access: &impl FlagAccess,
    edges: &[SimplexEntry],
    edge: SimplexEntry,
    additions: &[EdgePosition],
    working: &Coboundary,
    pivot: SimplexEntry,
) {
    let mut expected = std::collections::BTreeSet::new();
    for source in std::iter::once(edge).chain(additions.iter().map(|position| edges[position.0])) {
        access
            .visit_cofacets(source, &mut || Ok(()), |row| {
                if !expected.insert(row) {
                    expected.remove(&row);
                }
                Ok(true)
            })
            .unwrap();
    }
    let mut actual = std::collections::BTreeSet::new();
    for row in working
        .iter()
        .map(|row| row.0.simplex())
        .chain(std::iter::once(pivot))
    {
        if !actual.insert(row) {
            actual.remove(&row);
        }
    }
    assert_eq!(
        actual, expected,
        "stored transformation must reconstruct its full column"
    );
}

#[test]
fn virtual_pairs_cancel_and_reconstruct_on_dense_and_sparse_tied_cycles() {
    use crate::filtration::flag::SparseFlag;
    use crate::filtration::threshold_rips_from_distances;

    let n = 12;
    // Integer geodesic distances on a cycle provide ties without rounding.
    let values: Vec<_> = (0..n)
        .flat_map(|b| (0..b).map(move |a| (b - a).min(n - b + a) as f64))
        .collect();
    let input = DissimilarityView::new(&values, n).unwrap();
    let mut virtual_additions = 0;
    let mut reconstructions = 0;
    let mut checked_transforms = 0;
    for cutoff in [1., 2., 4., 6.] {
        let options = RipsOptions::new(1, Some(cutoff)).unwrap();
        let expected = reference::compute(input, &options).unwrap();
        let (_, coverage) = resolve_rips_range(input, &options);
        let dense = DenseFlag::new(input.into(), cutoff).unwrap();
        let graph = threshold_rips_from_distances(input.into(), Some(cutoff)).unwrap();
        let sparse = SparseFlag::new(graph.graph(), cutoff).unwrap();
        for stats in [
            check_virtual_access(&dense, coverage, &expected),
            check_virtual_access(&sparse, coverage, &expected),
        ] {
            virtual_additions += stats.virtual_additions;
            reconstructions += stats.column_additions;
            checked_transforms += stats.checked_transforms;
        }
    }
    assert!(
        virtual_additions > 1,
        "virtual additions={virtual_additions}"
    );
    assert!(
        reconstructions > 0,
        "ordinary reconstructions={reconstructions}"
    );
    assert!(checked_transforms > 0);
}

struct CancelAccess<'a, A> {
    inner: &'a A,
    flag: &'a std::sync::atomic::AtomicBool,
    cancel_at: usize,
    candidates: std::cell::Cell<usize>,
}

impl<A: FlagAccess> FlagAccess for CancelAccess<'_, A> {
    fn vertex_count(&self) -> usize {
        self.inner.vertex_count()
    }
    fn edges(&self, checkpoint: &mut impl FnMut() -> Result<()>) -> Result<Vec<SimplexEntry>> {
        self.inner.edges(checkpoint)
    }
    fn edge_vertices(&self, id: usize) -> [usize; 2] {
        self.inner.edge_vertices(id)
    }
    fn triangle_vertices(&self, id: usize) -> [usize; 3] {
        self.inner.triangle_vertices(id)
    }
    fn latest_facet(&self, triangle: SimplexEntry) -> SimplexEntry {
        self.inner.latest_facet(triangle)
    }
    fn visit_cofacets(
        &self,
        edge: SimplexEntry,
        checkpoint: &mut impl FnMut() -> Result<()>,
        visitor: impl FnMut(SimplexEntry) -> Result<bool>,
    ) -> Result<()> {
        self.inner.visit_cofacets(
            edge,
            &mut || {
                let count = self.candidates.get() + 1;
                self.candidates.set(count);
                if count == self.cancel_at {
                    self.flag.store(true, std::sync::atomic::Ordering::Relaxed);
                }
                checkpoint()
            },
            visitor,
        )
    }
}

#[test]
fn cancellation_at_each_cofacet_checkpoint_preserves_input_and_reentrancy() {
    use std::sync::atomic::{AtomicBool, Ordering};

    let n = 8;
    let values: Vec<_> = (0..n)
        .flat_map(|b| (0..b).map(move |a| (b - a).min(n - b + a) as f64))
        .collect();
    let input = DissimilarityView::new(&values, n).unwrap();
    let dense = DenseFlag::new(input.into(), input.diameter()).unwrap();
    let expected = reference::compute(input, &RipsOptions::default()).unwrap();
    let flag = AtomicBool::new(false);
    let limits = crate::persistence::ExecutionLimits::new(None, Some(&flag));
    let mut access = CancelAccess {
        inner: &dense,
        flag: &flag,
        cancel_at: usize::MAX,
        candidates: std::cell::Cell::new(0),
    };
    let mut stats = Stats::default();
    run_access::<true, true, ALL_SHORTCUTS, true>(
        &access,
        &mut stats,
        &mut WorkBudget::new(&limits).unwrap(),
    )
    .unwrap();
    assert!(stats.apparent_candidates > 0 && stats.virtual_additions > 0);
    let total = access.candidates.get();
    for cancel_at in 1..=total {
        access.cancel_at = cancel_at;
        access.candidates.set(0);
        let error = run_access::<true, true, ALL_SHORTCUTS, true>(
            &access,
            &mut Stats::default(),
            &mut WorkBudget::new(&limits).unwrap(),
        )
        .unwrap_err();
        assert_eq!(error, Error::Cancelled, "checkpoint {cancel_at}");
        assert!(flag.load(Ordering::Relaxed));
        flag.store(false, Ordering::Relaxed);
        let raw = run_access::<true, true, ALL_SHORTCUTS, true>(
            &dense,
            &mut Stats::default(),
            &mut WorkBudget::new(&limits).unwrap(),
        )
        .unwrap();
        assert_eq!(
            assemble_diagram(1, Coverage::Complete, raw).unwrap(),
            expected
        );
    }
}

// Independent forward boundary matrix, including zero-length index pairs.
// Neither production cofacet enumeration nor reduction is used for expectations.
fn boundary_death_triangles(
    n: usize,
    values: &[f64],
    cutoff: f64,
) -> std::collections::BTreeSet<[usize; 3]> {
    use std::collections::{BTreeSet, HashMap};
    let mut cells = Vec::new();
    for mask in 1_usize..1 << n {
        if mask.count_ones() > 3 {
            continue;
        }
        let mut value = 0_f64;
        for b in 0..n {
            for a in 0..b {
                if mask & (1 << a) != 0 && mask & (1 << b) != 0 {
                    value = value.max(values[b * (b - 1) / 2 + a]);
                }
            }
        }
        if value <= cutoff {
            cells.push((mask, value));
        }
    }
    cells.sort_by(|a, b| {
        a.1.total_cmp(&b.1)
            .then(a.0.count_ones().cmp(&b.0.count_ones()))
            .then(b.0.cmp(&a.0))
    });
    let positions: HashMap<_, _> = cells.iter().enumerate().map(|(i, c)| (c.0, i)).collect();
    let mut reduced: Vec<BTreeSet<usize>> = Vec::new();
    let mut owners: HashMap<usize, usize> = HashMap::new();
    let mut deaths = BTreeSet::new();
    for (position, &(mask, _)) in cells.iter().enumerate() {
        let mut column = BTreeSet::new();
        if mask.count_ones() > 1 {
            for v in 0..n {
                if mask & (1 << v) != 0 {
                    column.insert(positions[&(mask ^ (1 << v))]);
                }
            }
        }
        while let Some(&pivot) = column.last() {
            if let Some(&owner) = owners.get(&pivot) {
                column = column
                    .symmetric_difference(&reduced[owner])
                    .copied()
                    .collect();
            } else {
                owners.insert(pivot, position);
                if mask.count_ones() == 3 {
                    let vertices: Vec<_> = (0..n).filter(|&v| mask & (1 << v) != 0).collect();
                    deaths.insert(vertices.try_into().unwrap());
                }
                break;
            }
        }
        reduced.push(column);
    }
    deaths
}

#[test]
fn complete_handoff_matches_independent_death_indices_with_all_shortcuts() {
    let mut seed = 20261004_u64;
    let mut total_virtual = 0;
    let mut total_stored = 0;
    let mut total_shortcuts = 0;
    for n in 3_usize..=6 {
        let count = n * (n - 1) / 2;
        let samples = if n <= 4 {
            3_usize.pow(count as u32)
        } else {
            64
        };
        for sample in 0..samples {
            let mut code = sample;
            let values: Vec<_> = (0..count)
                .map(|_| {
                    if n <= 4 {
                        let v = code % 3;
                        code /= 3;
                        v as f64
                    } else {
                        seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
                        ((seed >> 32) % 4) as f64
                    }
                })
                .collect();
            let input = DissimilarityView::new(&values, n).unwrap();
            for cutoff in [0., 1., 2.] {
                let access = DenseFlag::new(input.into(), cutoff).unwrap();
                let expected = boundary_death_triangles(n, &values, cutoff);
                let execution = crate::execution::Execution::default();
                macro_rules! check {
                    ($shortcuts:ident) => {{
                        let mut stats = Stats::default();
                        let mut keys = super::super::clearing::TrianglePivots::default();
                        reduce_edges::<true, true, $shortcuts, true>(
                            &access,
                            access.edges(&mut || Ok(())).unwrap(),
                            &mut stats,
                            &mut WorkBudget::new(&execution).unwrap(),
                            &mut keys,
                        )
                        .unwrap();
                        let keys = keys
                            .into_triangles(n, &mut WorkBudget::new(&execution).unwrap())
                            .unwrap();
                        assert_eq!(keys.len(), stats.stored_columns + stats.skipped_apparent);
                        let actual: std::collections::BTreeSet<_> = keys.iter().copied().collect();
                        assert_eq!(keys.len(), actual.len(), "duplicate death key");
                        assert_eq!(
                            actual, expected,
                            "n={n} sample={sample} cutoff={cutoff} shortcuts={}",
                            $shortcuts
                        );
                        total_virtual += stats.skipped_apparent;
                        total_stored += stats.stored_columns;
                        total_shortcuts += stats.shortcuts;
                    }};
                }
                check!(NO_SHORTCUTS);
                check!(APPARENT_EMERGENT);
                check!(ALL_SHORTCUTS);
                // IgnoreTriangles is zero-sized and cannot retain clearing keys.
                assert_eq!(std::mem::size_of::<IgnoreTriangles>(), 0);
                run_access::<true, true, ALL_SHORTCUTS, true>(
                    &access,
                    &mut Stats::default(),
                    &mut WorkBudget::new(&execution).unwrap(),
                )
                .unwrap();
            }
        }
    }
    assert!(total_virtual > 0 && total_stored > 0 && total_shortcuts > 0);
}

struct Matrix(Vec<Vec<usize>>);
impl FilteredBoundary for Matrix {
    fn len(&self) -> usize {
        self.0.len()
    }
    fn dimension(&self, _: usize) -> usize {
        unreachable!()
    }
    fn value(&self, _: usize) -> f64 {
        unreachable!()
    }
    fn write_boundary(&self, j: usize, out: &mut Vec<usize>) -> Result<()> {
        out.clear();
        out.extend_from_slice(&self.0[j]);
        Ok(())
    }
}

#[test]
fn reversed_transpose_pairs_and_unpaired_indices_map_back_to_boundary_reduction() {
    // Enumerate subsets independently of production combinatorial indexing.
    for n in 1_usize..=6 {
        let values: Vec<_> = (0..n * (n - 1) / 2)
            .map(|i| ((i * 13) % 5) as f64)
            .collect();
        let input = DissimilarityView::new(&values, n).unwrap();
        for cutoff in [0.0, 2.0, 4.0] {
            let mut cells: Vec<_> = (1_usize..(1 << n))
                .filter(|mask| mask.count_ones() <= 3)
                .map(|mask| {
                    let v: Vec<_> = (0..n).filter(|i| mask & (1 << i) != 0).collect();
                    let value = v
                        .iter()
                        .flat_map(|&a| v.iter().map(move |&b| input.get(a, b).unwrap()))
                        .fold(0.0_f64, f64::max);
                    (mask, v.len(), value)
                })
                .filter(|x| x.2 <= cutoff)
                .collect();
            // Within a fixed dimension, mask order equals combinatorial order.
            cells.sort_by(|a, b| a.2.total_cmp(&b.2).then(a.1.cmp(&b.1)).then(b.0.cmp(&a.0)));
            let size = cells.len();
            let d = Matrix(
                cells
                    .iter()
                    .map(|&(mask, dim, _)| {
                        cells
                            .iter()
                            .enumerate()
                            .filter_map(|(i, &(face, fdim, _))| {
                                (fdim + 1 == dim && face & mask == face).then_some(i)
                            })
                            .collect()
                    })
                    .collect(),
            );
            let mut c = Matrix(vec![vec![]; size]);
            for (j, rows) in d.0.iter().enumerate() {
                for &i in rows {
                    c.0[size - 1 - i].push(size - 1 - j);
                }
            }
            for rows in &mut c.0 {
                rows.sort_unstable();
            }
            let forward = crate::persistence::reference::reduction::reduce(&d).unwrap();
            let dual = crate::persistence::reference::reduction::reduce(&c).unwrap();
            let mut expected = forward.pairs;
            let mut actual: Vec<_> = dual
                .pairs
                .into_iter()
                .map(|(i, j)| (size - 1 - j, size - 1 - i))
                .collect();
            expected.sort_unstable();
            actual.sort_unstable();
            assert_eq!(actual, expected);
            let mut actual: Vec<_> = dual.unpaired.into_iter().map(|i| size - 1 - i).collect();
            actual.sort_unstable();
            assert_eq!(actual, forward.unpaired);
        }
    }
}

#[test]
fn shortcut_and_reconstruction_paths_are_exercised() {
    let n = 24;
    let values: Vec<_> = (0..n)
        .flat_map(|b| {
            (0..b).map(move |a| {
                let angle = std::f64::consts::PI * (b - a) as f64 / n as f64;
                2.0 * angle.sin()
            })
        })
        .collect();
    let input = DissimilarityView::new(&values, n).unwrap();
    let mut stats = Stats::default();
    let raw =
        run::<true, true, true, APPARENT_EMERGENT>(input, input.diameter(), &mut stats).unwrap();
    let expected = reference::compute(input, &RipsOptions::default()).unwrap();
    assert_eq!(
        assemble_diagram(1, Coverage::Complete, raw).unwrap(),
        expected
    );
    assert!(
        stats.shortcuts > 0 && stats.column_additions > 0 && stats.stored_entries > 0,
        "{stats:?}"
    );
}

#[test]
fn working_rows_preserve_f64_bits_order_and_f2_multiplicity() {
    // Values and IDs are listed in the required forward order independently.
    let values = [
        0.0,
        f64::from_bits(1),
        f64::MIN_POSITIVE,
        1.0,
        f64::from_bits(1.0_f64.to_bits() + 1),
        1.0e300,
        f64::MAX,
    ];
    let rows: Vec<_> = values
        .into_iter()
        .flat_map(|value| [usize::MAX, 257, 1, 0].map(|id| SimplexEntry { id, value }))
        .collect();
    for (i, &left) in rows.iter().enumerate() {
        let encoded = WorkingEntry::from(left);
        assert_eq!(encoded.simplex().value.to_bits(), left.value.to_bits());
        assert_eq!(encoded.simplex().id, left.id);
        for (j, &right) in rows.iter().enumerate() {
            let right = WorkingEntry::from(right);
            assert_eq!(encoded.cmp(&right), i.cmp(&j));
            assert_eq!(encoded < right, i < j);
            assert_eq!(encoded <= right, i <= j);
            assert_eq!(encoded > right, i > j);
            assert_eq!(encoded >= right, i >= j);
            assert_eq!(encoded == right, i == j);
        }
    }
    let mut heap = Coboundary::new();
    for (i, &row) in rows.iter().enumerate().rev() {
        // Two copies vanish, while three copies retain their common row.
        for _ in 0..2 + usize::from(i % 3 == 0) {
            heap.push(Reverse(row.into()));
        }
    }
    let mut budget = WorkBudget::new(&crate::execution::Execution::default()).unwrap();
    for &expected in rows.iter().step_by(3) {
        assert_eq!(
            pop_parity(&mut heap, &mut budget).unwrap(),
            Some(Reverse(expected.into()))
        );
    }
    assert!(heap.is_empty());
}

#[test]
fn cofacet_batches_cancel_shared_triangles_and_recover_after_interruption() {
    // 01 and 02 share triangle 012. The remaining triangles 013 and 023 survive.
    let input = DissimilarityView::new(&[1.; 6], 4).unwrap();
    let access = DenseFlag::new(input.into(), 1.).unwrap();
    let mut heap = Coboundary::new();
    let mut scratch = Vec::with_capacity(8);
    let capacity = scratch.capacity();
    let mut stats = Stats::default();
    let mut budget = WorkBudget::new(&crate::execution::Execution::default()).unwrap();
    for id in [0, 1] {
        append_coboundary(
            &access,
            SimplexEntry { id, value: 1. },
            &mut heap,
            &mut scratch,
            &mut stats,
            &mut budget,
        )
        .unwrap();
        assert!(scratch.is_empty());
        assert_eq!(scratch.capacity(), capacity);
    }
    for id in [2, 1] {
        assert_eq!(
            pop_parity(&mut heap, &mut budget).unwrap(),
            Some(Reverse(SimplexEntry { id, value: 1. }.into()))
        );
    }
    assert_eq!(pop_parity(&mut heap, &mut budget).unwrap(), None);
    let edge = SimplexEntry { id: 0, value: 1. };
    let mut limited = WorkBudget::new(&crate::execution::Execution::default().max_work(1)).unwrap();
    assert_eq!(
        append_coboundary(
            &access,
            edge,
            &mut heap,
            &mut scratch,
            &mut stats,
            &mut limited
        ),
        Err(Error::WorkLimitExceeded { limit: 1 })
    );
    assert!(heap.is_empty());
    // Reuse both buffers after a partial enumeration; its prefix must be cleared.
    append_coboundary(
        &access,
        edge,
        &mut heap,
        &mut scratch,
        &mut stats,
        &mut budget,
    )
    .unwrap();
    for id in [1, 0] {
        assert_eq!(
            pop_parity(&mut heap, &mut budget).unwrap(),
            Some(Reverse(SimplexEntry { id, value: 1. }.into()))
        );
    }
    assert!(heap.is_empty());
    let flag = std::sync::atomic::AtomicBool::new(false);
    let controls = crate::execution::Execution::default()
        .max_work(8)
        .cancellation(&flag);
    let mut cancelled = WorkBudget::new(&controls).unwrap();
    // Cancel after all four candidates, before committing the prepared batch.
    cancelled.cancel_at_work(4);
    assert_eq!(
        append_coboundary(
            &access,
            edge,
            &mut heap,
            &mut scratch,
            &mut stats,
            &mut cancelled
        ),
        Err(Error::Cancelled)
    );
    assert!(heap.is_empty());
}
