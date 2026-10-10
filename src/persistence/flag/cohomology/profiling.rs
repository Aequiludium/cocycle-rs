//! Explicit developer profiling; excluded from normal test runs.

use super::*;
use crate::diagram::{Coverage, IntervalEnd};
use crate::persistence::rips::resolve_rips_range;
use crate::persistence::{RipsOptions, assemble_diagram, reference};

#[test]
fn h2_workspace_failures_discard_scratch_and_retry_repeated_intervals() {
    use crate::complex::{WeightedEdge, WeightedGraph};
    use crate::filtration::flag::{CliqueAccess, SparseFlag};
    use crate::persistence::simplicial::cohomology::workspace_force_failure;
    // Two disjoint octahedral spheres and one isolate: H2 has multiplicity two.
    let graph = WeightedGraph::new(
        13,
        (0..12)
            .flat_map(|b| {
                (0..b).filter_map(move |a| {
                    (a / 6 == b / 6 && a / 2 != b / 2).then_some(WeightedEdge {
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
        compute_h2(
            &SparseFlag::new(&graph, 1.).unwrap(),
            &CliqueAccess::Sparse(&graph, 1.),
            &mut WorkBudget::new(&execution).unwrap(),
        )
    };
    let expected = assemble_diagram(2, Coverage::Complete, compute().unwrap()).unwrap();
    assert_eq!(expected.dimension(2).unwrap().len(), 2);
    for at in [
        "h2_handoff_conversion_overlap",
        "h2_handoff_converted",
        "h2_triangle_assembly",
        "h2_reduction_start",
        "h2_column_complete",
        "h2_reduction_end",
    ] {
        for error in [
            Error::Cancelled,
            Error::AllocationFailed {
                context: "H2 workspace test",
            },
        ] {
            workspace_force_failure(Some((at, error.clone())));
            let interrupted = compute();
            workspace_force_failure(None);
            assert_eq!(interrupted, Err(error), "{at}");
            assert_eq!(
                assemble_diagram(2, Coverage::Complete, compute().unwrap()).unwrap(),
                expected,
                "{at}"
            );
        }
    }
}

pub(super) fn h2_workspace_event(
    event: &str,
    raw: &RawIntervals,
    triangles: &Vec<crate::filtration::flag::TupleEntry<3>>,
    cleared: &std::collections::HashSet<[usize; 3]>,
    owners: &HashMap<[usize; 4], usize>,
    columns: &Vec<Vec<usize>>,
    scratch: (usize, usize),
) -> Result<()> {
    crate::persistence::simplicial::cohomology::workspace_event(
        event,
        2,
        &[
            (
                "intervals",
                raw.capacity() * std::mem::size_of::<(usize, f64, Option<f64>)>(),
            ),
            (
                "triangles",
                triangles.capacity()
                    * std::mem::size_of::<crate::filtration::flag::TupleEntry<3>>(),
            ),
            (
                "columns",
                columns.capacity() * std::mem::size_of::<Vec<usize>>(),
            ),
            (
                "transform_payload",
                columns
                    .iter()
                    .map(|column| column.capacity() * std::mem::size_of::<usize>())
                    .sum(),
            ),
            ("working", scratch.0),
            ("transform_scratch", scratch.1),
        ],
        &[
            ("cleared", cleared.len()),
            ("clearing_slots", cleared.capacity()),
            ("owners", owners.len()),
            ("owner_slots", owners.capacity()),
            ("stored_columns", columns.len()),
        ],
    )
}

/// Dedicated counters and ownership landmarks. Ordinary latency uses the
/// pipeline worker with --cfg cocycle_h2_bench, without test instrumentation.
#[test]
#[ignore = "explicit H2 prototype counters; one fixture per fresh process"]
fn profile_h2_prototype() {
    use crate::complex::{WeightedEdge, WeightedGraph};
    use crate::filtration::flag::{CliqueAccess, SparseFlag};
    use crate::geometry::{DissimilarityMatrixView, MatrixLayout};
    let text = std::fs::read_to_string(std::env::var("COCYCLE_H2_FIXTURE").unwrap()).unwrap();
    let mut words = text.split_whitespace();
    let mode = words.next().unwrap();
    let n = words.next().unwrap().parse().unwrap();
    assert_eq!(words.next().unwrap(), "2");
    let cutoff = words.next().unwrap();
    let cutoff = if cutoff == "none" {
        f64::INFINITY
    } else {
        cutoff.parse().unwrap()
    };
    let count: usize = words.next().unwrap().parse().unwrap();
    assert_eq!(words.next().unwrap(), "2");
    let execution = crate::execution::Execution::default();
    let mut budget = WorkBudget::new(&execution).unwrap();
    let mut stats = H2Stats::default();
    crate::filtration::flag::cofacet_counts(true);
    let raw = if mode == "dense" {
        let values: Vec<f64> = words.map(|v| v.parse().unwrap()).collect();
        assert_eq!(values.len(), count);
        let matrix = DissimilarityMatrixView::new(&values, n, MatrixLayout::LowerTriangle).unwrap();
        let stop = cutoff.min(cone_radius(matrix, &mut || budget.step()).unwrap());
        run_h2(
            &DenseFlag::new(matrix, stop).unwrap(),
            &CliqueAccess::Dense(matrix, stop),
            &mut budget,
            &mut stats,
        )
        .unwrap()
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
        run_h2(
            &SparseFlag::new(&graph, cutoff).unwrap(),
            &CliqueAccess::Sparse(&graph, cutoff),
            &mut budget,
            &mut stats,
        )
        .unwrap()
    };
    println!("h2-prototype {stats:?}");
    println!(
        "cofacet_counts[candidates,edge_queries_or_intersection_comparisons,emitted]={:?}",
        crate::filtration::flag::cofacet_counts(false)
    );
    // Raw intervals permit independent comparison with the frozen pipeline
    // result; do not accept diagnostics merely because the process exited zero.
    println!("raw_intervals={raw:?}");
}

/// One fresh test process per stage/fixture; tools/profile_rips.py orchestrates.
#[test]
#[ignore = "explicit development profiling, no timing assertions"]
fn profile_stage() {
    use std::hint::black_box;
    use std::time::Instant;
    let bytes = std::fs::read(std::env::var("COCYCLE_ABLATION_FIXTURE").unwrap()).unwrap();
    assert!(bytes.len() >= 48 && &bytes[..8] == b"COCYCLE1");
    let integer = |i| u64::from_le_bytes(bytes[i..i + 8].try_into().unwrap()) as usize;
    assert_eq!(integer(8), 0, "distance fixtures only");
    assert_eq!(integer(32), 1, "H1 fixtures only");
    let cutoff = f64::from_le_bytes(bytes[40..48].try_into().unwrap());
    let values: Vec<_> = bytes[48..]
        .as_chunks::<8>()
        .0
        .iter()
        .map(|x| f64::from_le_bytes(*x))
        .collect();
    let input = DissimilarityView::new(&values, integer(16)).unwrap();
    let options = RipsOptions::new(1, (!cutoff.is_nan()).then_some(cutoff)).unwrap();
    let (cutoff, coverage) = resolve_rips_range(input, &options);
    let stage = std::env::var("COCYCLE_ABLATION_STAGE").unwrap();
    let execute = |stats: &mut Stats| {
        stats.two_pass_initialization = matches!(stage.as_str(), "two-pass" | "virtual-two-pass");
        match stage.as_str() {
            "explicit" => run::<false, false, false, NO_SHORTCUTS>(input, cutoff, stats),
            "clearing" => run::<false, true, false, NO_SHORTCUTS>(input, cutoff, stats),
            "implicit" => run::<true, true, false, NO_SHORTCUTS>(input, cutoff, stats),
            "cone" => run::<true, true, true, NO_SHORTCUTS>(input, cutoff, stats),
            "apparent" => run::<true, true, true, APPARENT>(input, cutoff, stats),
            "emergent" | "two-pass" => {
                run::<true, true, true, APPARENT_EMERGENT>(input, cutoff, stats)
            }
            "virtual" | "virtual-two-pass" => {
                run::<true, true, true, ALL_SHORTCUTS>(input, cutoff, stats)
            }
            _ => panic!("unknown stage"),
        }
    };
    let rss = |field: &str| -> Option<usize> {
        std::fs::read_to_string("/proc/self/status")
            .ok()?
            .lines()
            .find_map(|line| {
                line.strip_prefix(field)?
                    .split_whitespace()
                    .next()?
                    .parse()
                    .ok()
            })
    };
    let before = rss("VmHWM:");
    let mut stats = Stats::default();
    let warmup = execute(&mut stats).unwrap();
    let mut samples = Vec::new();
    for _ in 0..5 {
        let start = Instant::now();
        drop(black_box(execute(&mut Stats::default()).unwrap()));
        samples.push(start.elapsed().as_secs_f64() * 1000.0);
    }
    let peak = rss("VmHWM:");
    // Reference execution is after measurement so it cannot pollute stage RSS.
    assert_eq!(
        assemble_diagram(1, coverage, warmup).unwrap(),
        reference::compute(input, &options).unwrap()
    );
    let json_number = |x: Option<usize>| x.map_or_else(|| "null".into(), |x| x.to_string());
    println!(
        "ABLATION {{\"stage\":\"{stage}\",\"samples_ms\":{samples:?},\"hwm_before_kib\":{},\"peak_rss_kib\":{},\"cofacets\":{},\"initial_candidates\":{},\"reconstruction_candidates\":{},\"apparent_candidates\":{},\"skipped_apparent\":{},\"virtual_additions\":{},\"stored_columns\":{},\"column_additions\":{},\"shortcuts\":{},\"peak_heap_entries\":{},\"stored_entries\":{}}}",
        json_number(before),
        json_number(peak),
        stats.cofacets,
        stats.initial_candidates,
        stats.reconstruction_candidates,
        stats.apparent_candidates,
        stats.skipped_apparent,
        stats.virtual_additions,
        stats.stored_columns,
        stats.column_additions,
        stats.shortcuts,
        stats.peak_heap,
        stats.stored_entries
    );
}

/// Large-input counters without running the cubic reference algorithm. The
/// Python controller MUST check the emitted diagram against a validated public
/// benchmark result before accepting these diagnostic counters.
#[test]
#[ignore = "development workload counters; validated by tools/profile_scaling.py"]
fn profile_workload() {
    let bytes = std::fs::read(std::env::var("COCYCLE_ABLATION_FIXTURE").unwrap()).unwrap();
    assert!(bytes.len() >= 48 && &bytes[..8] == b"COCYCLE1");
    let integer = |i| u64::from_le_bytes(bytes[i..i + 8].try_into().unwrap()) as usize;
    assert_eq!(integer(8), 0);
    assert_eq!(integer(32), 1);
    let cutoff = f64::from_le_bytes(bytes[40..48].try_into().unwrap());
    let values: Vec<_> = bytes[48..]
        .as_chunks::<8>()
        .0
        .iter()
        .map(|x| f64::from_le_bytes(*x))
        .collect();
    let input = DissimilarityView::new(&values, integer(16)).unwrap();
    let options = RipsOptions::new(1, (!cutoff.is_nan()).then_some(cutoff)).unwrap();
    let (cutoff, coverage) = resolve_rips_range(input, &options);
    let mut stats = Stats::default();
    let raw = run::<true, true, true, PRODUCTION_SHORTCUTS>(input, cutoff, &mut stats).unwrap();
    let diagram = assemble_diagram(1, coverage, raw).unwrap();
    print!(
        "WORKLOAD {{\"statistics\":{{\"edges\":{},\"cofacets\":{},\"initial_candidates\":{},\"reconstruction_candidates\":{},\"apparent_candidates\":{},\"skipped_apparent\":{},\"virtual_additions\":{},\"stored_columns\":{},\"column_additions\":{},\"shortcuts\":{},\"peak_coboundary_heap\":{},\"stored_transform_entries\":{},\"largest_transform\":{},\"peak_transform_heap\":{}}},",
        stats.edges,
        stats.cofacets,
        stats.initial_candidates,
        stats.reconstruction_candidates,
        stats.apparent_candidates,
        stats.skipped_apparent,
        stats.virtual_additions,
        stats.stored_columns,
        stats.column_additions,
        stats.shortcuts,
        stats.peak_heap,
        stats.stored_entries,
        stats.largest_transform,
        stats.peak_transform_heap
    );
    match coverage {
        Coverage::Complete => print!("\"coverage\":[\"complete\",null],"),
        Coverage::Through(t) => print!("\"coverage\":[\"through\",{t}],"),
    }
    print!("\"intervals\":[");
    for (i, bar) in diagram.intervals().enumerate() {
        if i > 0 {
            print!(",");
        }
        let (kind, endpoint) = match bar.end() {
            IntervalEnd::Finite(d) => ("F", d),
            IntervalEnd::Essential => ("E", 0.0),
            IntervalEnd::RightCensored { through } => ("C", through),
        };
        print!(
            "[{},{},\"{kind}\",{endpoint}]",
            bar.dimension(),
            bar.birth()
        );
    }
    println!("]}}");
}
