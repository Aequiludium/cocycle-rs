//! Analyze hand-built interval multisets without constructing a complex.

use cocycle::Result;
use cocycle::descriptors::{betti_curve, finite_lifetime_summary};
use cocycle::diagram::{
    ComputedDimensions, Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval,
};

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args_os().skip(1).collect();
    match args.as_slice() {
        [] => analyze_diagrams()?,
        [flag, path] if flag == "--export" => {
            let mut output = std::io::BufWriter::new(std::fs::File::create(path)?);
            export_graph_features(&mut output)?;
            std::io::Write::flush(&mut output)?;
        }
        _ => return Err("usage: diagram_analysis [--export PATH]".into()),
    }
    Ok(())
}

fn analyze_diagrams() -> Result<()> {
    // Two repeated H1 intervals of lifetime two, plus an essential class.
    // H0 is a computed empty dimension. Input order need not be sorted.
    let finite = PersistenceInterval::new(1, -1.0, IntervalEnd::Finite(1.0))?;
    let complete = PersistenceDiagram::new(
        1,
        Coverage::Complete,
        vec![
            PersistenceInterval::new(1, 0.0, IntervalEnd::Essential)?,
            finite,
            finite,
        ],
    )?;
    let grid = [-2.0, -1.0, 0.0, 1.0, 2.0];
    let counts = betti_curve(&complete, 1, &grid)?;
    let summary = finite_lifetime_summary(&complete, 1)?;
    assert_eq!(counts, [0, 2, 3, 1, 1]);
    assert_eq!(summary.total_persistence(), 4.0);
    assert!((summary.entropy().unwrap() - 2.0_f64.ln()).abs() < 1e-14);
    assert_eq!(summary.excluded_essential_count(), 1);
    println!("Complete H1 Betti curve at {grid:?}: {counts:?}");
    println!("Observed finite H1 lifetimes: {summary:?}");

    // This is a separate supplied diagram, observed only through scale one.
    // The class born at zero is censored: its eventual death is unknown.
    let truncated = PersistenceDiagram::new(
        1,
        Coverage::Through(1.0),
        vec![
            finite,
            finite,
            PersistenceInterval::new(1, 0.0, IntervalEnd::RightCensored { through: 1.0 })?,
        ],
    )?;
    assert_eq!(betti_curve(&truncated, 1, &grid[..4])?, [0, 2, 3, 1]);
    let summary = finite_lifetime_summary(&truncated, 1)?;
    assert_eq!(summary.excluded_censored_count(), 1);
    println!("Truncated H1 finite lifetimes (censored class excluded): {summary:?}");
    println!(
        "Computed empty H0: {:?}",
        betti_curve(&truncated, 0, &[0.0])?
    );
    // An H1-only algorithm must not declare H0 computed merely by omitting bars.
    let h1_only = PersistenceDiagram::with_dimensions(
        ComputedDimensions::new(vec![1])?,
        Coverage::Complete,
        vec![finite],
    )?;
    assert_eq!(betti_curve(&h1_only, 1, &[0.])?, [1]);
    assert!(betti_curve(&h1_only, 0, &[0.]).is_err());
    println!(
        "H1-only computed dimensions: {:?}",
        h1_only.computed_dimensions().iter().collect::<Vec<_>>()
    );
    Ok(())
}

// This branch reads graph storage directly: no diagram, PH or distance call.
fn graph_features(graph: &cocycle::complex::WeightedGraph) -> Option<[f64; 4]> {
    let n = graph.vertex_count();
    if n == 0 {
        return None;
    }
    let leaves = (0..n)
        .filter(|&vertex| graph.neighbors(vertex).unwrap().len() == 1)
        .count(); // All queried vertices are in the validated owner's range.
    let m = graph.edge_count();
    Some([
        n as f64,
        m as f64,
        2.0 * m as f64 / n as f64,
        leaves as f64 / n as f64,
    ])
}

fn teaching_graph(n: usize, cycle: bool) -> Result<cocycle::complex::WeightedGraph> {
    use cocycle::complex::{WeightedEdge, WeightedGraph};
    let mut edges: Vec<_> = (1..n)
        .map(|vertex| WeightedEdge {
            vertices: [vertex - 1, vertex],
            value: 1.0,
        })
        .collect();
    if cycle {
        edges.push(WeightedEdge {
            vertices: [n - 1, 0],
            value: 1.0,
        });
    }
    WeightedGraph::new(n, edges)
}

fn export_graph_features(
    output: &mut impl std::io::Write,
) -> std::result::Result<(), Box<dyn std::error::Error>> {
    writeln!(
        output,
        "sample_id,split,label,source,vertex_count,edge_count,mean_degree,leaf_fraction"
    )?;
    // Split is fixed before any data is fitted; labels describe construction.
    for (n, split) in [
        (3, "train"),
        (4, "train"),
        (6, "train"),
        (5, "test"),
        (7, "test"),
    ] {
        for cycle in [false, true] {
            let label = if cycle { "cycle" } else { "path" };
            let graph = teaching_graph(n, cycle)?;
            let [vertices, edges, mean, leaves] = graph_features(&graph).unwrap();
            // These teaching fixtures all have at least three vertices.
            writeln!(
                output,
                "{label}-{n},{split},{label},unweighted-graph-v1,{vertices},{edges},{mean},{leaves}"
            )?;
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use cocycle::complex::WeightedGraph;

    #[test]
    fn direct_graph_features_match_hand_counts_without_mutating_storage() -> Result<()> {
        let path = teaching_graph(4, false)?;
        assert_eq!(graph_features(&path), Some([4.0, 3.0, 1.5, 0.5]));
        assert_eq!(path.edge_count(), 3);
        assert_eq!(path.neighbors(0).unwrap().len(), 1);
        assert_eq!(
            graph_features(&teaching_graph(4, true)?),
            Some([4.0, 4.0, 2.0, 0.0])
        );
        assert_eq!(graph_features(&WeightedGraph::new(0, vec![])?), None);
        assert_eq!(
            graph_features(&WeightedGraph::new(3, vec![])?),
            Some([3.0, 0.0, 0.0, 0.0])
        );
        Ok(())
    }

    #[test]
    fn csv_preserves_sample_label_and_split_order() {
        let mut output = Vec::new();
        export_graph_features(&mut output).unwrap();
        let csv = String::from_utf8(output).unwrap();
        let rows: Vec<_> = csv.lines().skip(1).collect();
        assert_eq!(rows.len(), 10);
        assert_eq!(
            rows[0],
            "path-3,train,path,unweighted-graph-v1,3,2,1.3333333333333333,0.6666666666666666"
        );
        assert_eq!(rows[7], "cycle-5,test,cycle,unweighted-graph-v1,5,5,2,0");
        assert!(rows[9].starts_with("cycle-7,test,cycle,"));
    }
}
