//! Native resource worker: supplied unweighted graph, preparation, then solve.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::hodge::{HodgeLaplacian, HodgeSelection, IterativeOptions};
use std::time::Instant;

fn main() -> cocycle::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let family = &args[1];
    let n: usize = args[2].parse().unwrap();
    let block: usize = args[3].parse().unwrap();
    let iterations: usize = args[4].parse().unwrap();
    let prepare = Instant::now();
    let mut cells: Vec<_> = (0..n).map(|i| Simplex::new(vec![i], 0.).unwrap()).collect();
    match family.as_str() {
        "path" => cells.extend((0..n - 1).map(|i| Simplex::new(vec![i, i + 1], 0.).unwrap())),
        "star" => cells.extend((1..n).map(|i| Simplex::new(vec![0, i], 0.).unwrap())),
        "isolates" => (),
        _ => panic!("unsupported resource fixture"),
    }
    let source = SimplicialComplex::new(cells)?;
    let operator = HodgeLaplacian::new(&source, 0, HodgeSelection::Full)?;
    let prepare_ns = prepare.elapsed().as_nanos();
    let policy = IterativeOptions {
        block_size: block,
        max_iterations: iterations,
        ..Default::default()
    };
    let begin = Instant::now();
    let result = operator.smallest_eigenpairs(3, &policy, &Execution::default());
    let solve_ns = begin.elapsed().as_nanos();
    match result {
        Ok(pairs) => println!(
            "{{\"status\":\"success\",\"prepare_ns\":{prepare_ns},\"solve_ns\":{solve_ns},\"values\":{:?},\"residuals\":{:?},\"iterations\":{},\"actions\":{},\"peak_subspace\":{},\"index_coverage_certified\":{}}}",
            pairs.values(),
            pairs.residuals(),
            pairs.iterations(),
            pairs.operator_applications(),
            pairs.peak_subspace_dimension(),
            pairs.is_index_coverage_certified()
        ),
        Err(error) => println!(
            "{{\"status\":\"numerical_failure\",\"prepare_ns\":{prepare_ns},\"solve_ns\":{solve_ns},\"error\":{:?}}}",
            error.to_string()
        ),
    }
    Ok(())
}
