//! Test-only work accounting; never used for ordinary latency measurements.
use super::*;
use std::cell::RefCell;

// Safe capacity observations, not allocator accounting. Hash-table buckets and
// BTreeMap nodes have no stable allocation-size API and remain unmeasured.
thread_local! {
    static WORKSPACE_TRACE: RefCell<Option<std::time::Instant>> = const { RefCell::new(None) };
    static WORKSPACE_FAILURE: RefCell<Option<(&'static str, Error)>> = const { RefCell::new(None) };
    static HEAP_GROWTH: RefCell<[usize; 2]> = const { RefCell::new([0; 2]) };
}
pub(in crate::persistence) fn workspace_heap_growth(before: usize, after: usize) {
    WORKSPACE_TRACE.with_borrow(|trace| {
        if trace.is_some() && after > before {
            HEAP_GROWTH.with_borrow_mut(|counts| counts[usize::from(before != 0)] += 1);
        }
    });
}
pub(in crate::persistence) fn workspace_event(
    event: &str,
    dimension: usize,
    capacities: &[(&str, usize)],
    counts: &[(&str, usize)],
) -> Result<()> {
    WORKSPACE_TRACE.with_borrow(|started| {
        let Some(started) = started else { return };
        let status = std::fs::read_to_string("/proc/self/status").unwrap_or_default();
        let memory = |field: &str| {
            status.lines().find_map(|line| {
                line.strip_prefix(field)?.split_whitespace().next()?.parse::<usize>().ok()
            }).map_or_else(|| "null".to_owned(), |value| value.to_string())
        };
        let fields = |values: &[(&str, usize)]| {
            values.iter().map(|(name, value)| format!("\"{name}\":{value}"))
                .collect::<Vec<_>>().join(",")
        };
        let growth = HEAP_GROWTH.with_borrow(|counts| *counts);
        println!(
            "workspace_event={{\"event\":\"{event}\",\"dimension\":{dimension},\"elapsed_ns\":{},\"rss_kib\":{},\"hwm_kib\":{},\"vec_capacity_bytes\":{{{}}},\"counts\":{{{}}},\"heap_capacity_new\":{},\"heap_capacity_grow\":{}}}",
            started.elapsed().as_nanos(), memory("VmRSS:"), memory("VmHWM:"),
            fields(capacities), fields(counts), growth[0], growth[1],
        );
    });
    WORKSPACE_FAILURE.with_borrow(|failure| match failure {
        Some((at, error)) if *at == event => Err(error.clone()),
        _ => Ok(()),
    })
}

pub(super) fn workspace_level(
    event: &str,
    dimension: usize,
    levels: &[(&str, &Vec<Simplex>)],
    cleared: &HashSet<Vec<usize>>,
    columns: &Vec<Column<usize>>,
    owners: &HashMap<Vec<usize>, usize>,
    raw: &RawIntervals,
) -> Result<()> {
    let mut capacities = vec![
        (
            "intervals",
            raw.capacity() * std::mem::size_of::<(usize, f64, Option<f64>)>(),
        ),
        (
            "columns",
            columns.capacity() * std::mem::size_of::<Column<usize>>(),
        ),
        (
            "clearing_keys",
            cleared
                .iter()
                .map(|key| key.capacity() * std::mem::size_of::<usize>())
                .sum(),
        ),
        (
            "owner_keys",
            owners
                .keys()
                .map(|key| key.capacity() * std::mem::size_of::<usize>())
                .sum(),
        ),
    ];
    for &(name, level) in levels {
        capacities.push((
            name,
            level.capacity() * std::mem::size_of::<Simplex>()
                + level
                    .iter()
                    .map(|simplex| simplex.vertices.capacity() * std::mem::size_of::<usize>())
                    .sum::<usize>(),
        ));
    }
    workspace_event(
        event,
        dimension,
        &capacities,
        &[
            ("cleared", cleared.len()),
            ("clearing_slots", cleared.capacity()),
            ("owners", owners.len()),
            ("owner_slots", owners.capacity()),
            ("stored_columns", columns.len()),
            (
                "tree_entries",
                columns.iter().map(|column| column.entries().count()).sum(),
            ),
        ],
    )
}

#[test]
fn workspace_boundary_failures_return_no_partial_result_and_retry_cleanly() {
    use crate::complex::{WeightedEdge, WeightedGraph};
    use crate::diagram::Coverage;
    use crate::persistence::assemble_diagram;
    // Octahedral sphere plus an isolate: two H0 classes and one essential H2.
    let graph = WeightedGraph::new(
        7,
        (0..6)
            .flat_map(|b| {
                (0..b).filter_map(move |a| {
                    (a / 2 != b / 2).then_some(WeightedEdge {
                        vertices: [a, b],
                        value: 1.,
                    })
                })
            })
            .collect(),
    )
    .unwrap();
    let compute = || {
        let execution = crate::execution::Execution::default();
        crate::persistence::flag::compute_graph(
            &graph,
            3,
            1.,
            PrimeField::new(2).unwrap(),
            &mut WorkBudget::new(&execution).unwrap(),
        )
    };
    let expected = assemble_diagram(3, Coverage::Complete, compute().unwrap()).unwrap();
    assert_eq!(expected.dimension(0).unwrap().len(), 7);
    assert_eq!(
        expected
            .intervals()
            .filter(
                |bar| bar.dimension() == 0 && bar.end() == crate::diagram::IntervalEnd::Essential
            )
            .count(),
        2
    );
    assert_eq!(expected.dimension(1).unwrap().len(), 0);
    assert_eq!(expected.dimension(2).unwrap().len(), 1);
    assert_eq!(expected.dimension(3).unwrap().len(), 0);
    for at in [
        "h1_handoff_complete",
        "h1_released",
        "handoff_converted",
        "triangle_assembly",
        "reduction_start",
        "reduction_end",
        "clear_next_extracted",
        "next_level_assembled",
    ] {
        for error in [
            Error::Cancelled,
            Error::AllocationFailed {
                context: "workspace boundary test",
            },
        ] {
            WORKSPACE_FAILURE.with_borrow_mut(|failure| *failure = Some((at, error.clone())));
            let interrupted = compute();
            WORKSPACE_FAILURE.with_borrow_mut(|failure| *failure = None);
            assert_eq!(interrupted, Err(error), "{at}");
            assert_eq!(
                assemble_diagram(3, Coverage::Complete, compute().unwrap()).unwrap(),
                expected,
                "{at}"
            );
        }
    }
}

#[derive(Default, Debug)]
pub(super) struct Stats {
    pub(super) simplices: usize,
    pub(super) processed: usize,
    pub(super) cleared: usize,
    pub(super) pivot_lookups: usize,
    pub(super) pivot_hits: usize,
    pub(super) column_additions: usize,
    pub(super) reconstructions: usize,
    pub(super) stored_transform_entries: usize,
    pub(super) largest_transform: usize,
    pub(super) peak_working_column: usize,
}
thread_local! {
    static COUNTERS: RefCell<Vec<Stats>> = const { RefCell::new(Vec::new()) };
}
pub(super) fn record(dimension: usize, update: impl FnOnce(&mut Stats)) {
    COUNTERS.with_borrow_mut(|rows| {
        rows.resize_with(rows.len().max(dimension + 1), Stats::default);
        update(&mut rows[dimension]);
    });
}
pub(in crate::persistence) fn take_counts() -> Vec<[usize; 3]> {
    COUNTERS.with_borrow_mut(|rows| {
        std::mem::take(rows)
            .into_iter()
            .map(|s| [s.simplices, s.processed, s.cleared])
            .collect()
    })
}

#[test]
#[ignore = "explicit counter process; ordinary latency uses pipeline release workers"]
fn profile_h2() {
    use crate::complex::{WeightedEdge, WeightedGraph};
    use crate::filtration::flag::CliqueAccess;
    use crate::geometry::{DissimilarityMatrixView, MatrixLayout};
    let text = std::fs::read_to_string(std::env::var("COCYCLE_H2_FIXTURE").unwrap()).unwrap();
    let mut words = text.split_whitespace();
    let mode = words.next().unwrap();
    let n = words.next().unwrap().parse().unwrap();
    let q = words.next().unwrap().parse().unwrap();
    let cutoff = words.next().unwrap();
    let cutoff = if cutoff == "none" {
        f64::INFINITY
    } else {
        cutoff.parse().unwrap()
    };
    let count: usize = words.next().unwrap().parse().unwrap();
    let field = PrimeField::new(words.next().unwrap().parse().unwrap()).unwrap();
    let dispatch = std::env::var("COCYCLE_H2_ROUTE").is_ok_and(|v| v == "dispatch");
    if std::env::var_os("COCYCLE_WORKSPACE_TRACE").is_some() {
        WORKSPACE_TRACE.with_borrow_mut(|trace| *trace = Some(std::time::Instant::now()));
        HEAP_GROWTH.with_borrow_mut(|counts| *counts = [0; 2]);
    }
    let execution = crate::execution::Execution::default();
    let mut budget = WorkBudget::new(&execution).unwrap();
    COUNTERS.with_borrow_mut(Vec::clear);
    crate::filtration::flag::cofacet_counts(true);
    if mode == "dense" {
        let values: Vec<f64> = words.map(|v| v.parse().unwrap()).collect();
        assert_eq!(values.len(), count);
        let matrix = DissimilarityMatrixView::new(&values, n, MatrixLayout::LowerTriangle).unwrap();
        let stop = cutoff
            .min(crate::filtration::rips::cone_radius(matrix, &mut || budget.step()).unwrap());
        if dispatch {
            crate::persistence::flag::compute_dense(matrix, q, cutoff, field, &mut budget).unwrap();
        } else {
            compute(&CliqueAccess::Dense(matrix, stop), q, field, &mut budget).unwrap();
        }
    } else {
        assert_eq!(mode, "flag");
        let edges = (0..count)
            .map(|_| WeightedEdge {
                vertices: [
                    words.next().unwrap().parse().unwrap(),
                    words.next().unwrap().parse().unwrap(),
                ],
                value: words.next().unwrap().parse().unwrap(),
            })
            .collect();
        let graph = WeightedGraph::new(n, edges).unwrap();
        if dispatch {
            crate::persistence::flag::compute_graph(&graph, q, cutoff, field, &mut budget).unwrap();
        } else {
            compute(&CliqueAccess::Sparse(&graph, cutoff), q, field, &mut budget).unwrap();
        }
    }
    COUNTERS.with_borrow(|rows| {
        for (dimension, stats) in rows.iter().enumerate().skip(1) {
            println!("dimension={dimension} {stats:?}");
            assert_eq!(stats.processed + stats.cleared, stats.simplices);
        }
    });
    println!(
        "cofacet_counts[candidates,edge_queries,emitted]={:?}",
        crate::filtration::flag::cofacet_counts(false)
    );
}
