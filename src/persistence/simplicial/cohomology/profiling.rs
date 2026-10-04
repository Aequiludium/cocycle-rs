//! Test-only work accounting; never used for ordinary latency measurements.
use super::*;
use std::cell::RefCell;

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
