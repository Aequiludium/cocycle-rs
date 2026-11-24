//! Ordinary q=1 HeCS example and frozen-input adapter; no persistent semantics.
use cocycle::Error;
use cocycle::complex::FilteredComplex;
use cocycle::execution::Execution;
use cocycle::hodge::{HecsOptions, HodgeLaplacian, HodgeSelection};
use std::io::{self, Read};

struct Cells(Vec<(usize, Vec<(usize, i32)>)>);
impl FilteredComplex for Cells {
    type CellId = usize;
    fn cells(&self) -> impl Iterator<Item = usize> + '_ {
        0..self.0.len()
    }
    fn dimension(&self, id: usize) -> usize {
        self.0[id].0
    }
    fn filtration_value(&self, _: usize) -> f64 {
        0.
    }
    fn boundary(&self, id: usize) -> impl Iterator<Item = (usize, i32)> + '_ {
        self.0[id].1.iter().copied()
    }
}

fn rows(a: &[f64], count: usize, columns: usize) -> Vec<&[f64]> {
    (0..count)
        .map(|i| &a[i * columns..(i + 1) * columns])
        .collect()
}
fn emit(input: &str, solution_only: bool) -> std::result::Result<(), Box<dyn std::error::Error>> {
    let mut tokens = input.split_whitespace();
    let mut take = || tokens.next().ok_or("missing fixture token");
    let n: usize = take()?.parse()?;
    let m: usize = take()?.parse()?;
    let t: usize = take()?.parse()?;
    let peel: usize = take()?.parse()?;
    if n > 256 || m > 512 || t > 1024 || peel > 1 {
        return Err("fixture limits exceeded".into());
    }
    let mut cells = Cells((0..n).map(|_| (0, vec![])).collect());
    let mut edges = Vec::new();
    for _ in 0..m {
        let edge = [take()?.parse::<usize>()?, take()?.parse::<usize>()?];
        if edge[0] >= edge[1] || edge[1] >= n {
            return Err("invalid fixture edge".into());
        }
        edges.push(edge);
        cells.0.push((1, vec![(edge[0], -1), (edge[1], 1)]));
    }
    let mut weights = Vec::new();
    for _ in 0..t {
        let face = [
            take()?.parse::<usize>()?,
            take()?.parse::<usize>()?,
            take()?.parse::<usize>()?,
        ];
        if face[0] >= face[1] || face[1] >= face[2] || face[2] >= n {
            return Err("invalid fixture face".into());
        }
        let mut boundary = Vec::new();
        for removed in 0..3 {
            let edge: Vec<_> = face
                .iter()
                .enumerate()
                .filter(|(i, _)| *i != removed)
                .map(|(_, &v)| v)
                .collect();
            let e = edges
                .iter()
                .position(|x| x.as_slice() == edge)
                .ok_or("missing lower edge")?;
            boundary.push((n + e, if removed % 2 == 0 { 1 } else { -1 }));
        }
        cells.0.push((2, boundary));
        weights.push(take()?.parse::<f64>()?);
    }
    let rhs = (0..m)
        .map(|_| Ok(take()?.parse::<f64>()?))
        .collect::<std::result::Result<Vec<_>, Box<dyn std::error::Error>>>()?;
    if tokens.next().is_some() {
        return Err("trailing fixture input".into());
    }
    let execution = Execution::default();
    let hodge = HodgeLaplacian::new_with(&cells, 1, HodgeSelection::Upper, &execution)?;
    let options = HecsOptions {
        peel_core: peel == 1,
        ..Default::default()
    };
    let solution = match hodge.solve_hecs(&weights, &rhs, &options, &execution) {
        Err(
            error @ Error::InvalidParameter {
                parameter: "HeCS rank",
                ..
            },
        ) => {
            println!(
                "{{\"status\":\"rank_loss\",\"error\":{:?}}}",
                error.to_string()
            );
            return Ok(());
        }
        other => other?,
    };
    if solution_only {
        println!(
            "{{\"values\":{:?},\"iterations\":{}}}",
            solution.values(),
            solution.residual_history().len().saturating_sub(1)
        );
        return Ok(());
    }
    let r = solution.rank();
    let trace: Vec<_> = solution
        .residual_history()
        .iter()
        .enumerate()
        .map(|(i, value)| format!("{{\"iteration\":{i},\"original_residual\":{value}}}"))
        .collect();
    println!(
        "{{\"selected\":{:?},\"pairs\":{:?},\"edge_permutation\":{:?},\"triangle_permutation\":{:?},\"factor\":{:?},\"active\":{:?},\"solves\":{{\"hecs\":{{\"status\":\"converged\",\"original_x\":{:?},\"trace\":[{}]}}}}}}",
        solution.selected_triangles(),
        solution.collapse_pairs(),
        solution.edge_permutation(),
        solution.triangle_permutation(),
        rows(solution.factor(), m, r),
        rows(solution.active_matrix(), r, r),
        solution.values(),
        trace.join(",")
    );
    Ok(())
}
fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args == ["--stdin"] || args == ["--stdin", "--solution-only"] {
        let mut input = String::new();
        io::stdin().read_to_string(&mut input)?;
        emit(&input, args.len() == 2)
    } else if args.is_empty() {
        // Four lexicographic faces, weights [1,2,4,8], compatible hand-derived b.
        emit(
            "4 6 4 0  0 1 0 2 0 3 1 2 1 3 2 3  0 1 2 1 0 1 3 2 0 2 3 4 1 2 3 8  2 -22 20 -22 24 -44",
            false,
        )
    } else {
        Err("expected --stdin [--solution-only] or no arguments".into())
    }
}
