//! Native resource worker for sparse paths, clustered roots and ill-conditioned chains.
use cocycle::complex::FilteredComplex;
use cocycle::dirac::{DiracOperator, IterativeOptions};
use cocycle::execution::Execution;
use std::time::Instant;

struct Chain {
    n: usize,
    family: String,
}
impl FilteredComplex for Chain {
    type CellId = usize;
    fn cells(&self) -> impl Iterator<Item = usize> + '_ {
        0..self.n
            + if self.family == "jordan" {
                self.n
            } else {
                self.n - 1
            }
    }
    fn dimension(&self, cell: usize) -> usize {
        usize::from(cell >= self.n)
    }
    fn filtration_value(&self, _: usize) -> f64 {
        0.
    }
    fn boundary(&self, cell: usize) -> impl Iterator<Item = (usize, i32)> + '_ {
        let i = cell.saturating_sub(self.n);
        let entries = if self.family == "jordan" {
            [(i, 1), (i.saturating_sub(1), -2)]
        } else if self.family == "star" {
            [(0, -1), (i + 1, 1)]
        } else {
            [(i, -1), (i + 1, 1)]
        };
        entries.into_iter().take(if cell < self.n {
            0
        } else if self.family == "jordan" && i == 0 {
            1
        } else {
            2
        })
    }
}
fn main() -> cocycle::Result<()> {
    let args: Vec<_> = std::env::args().collect();
    let family = args[1].clone();
    let n: usize = args[2].parse().unwrap();
    let count: usize = args[3].parse().unwrap();
    let block: usize = args[4].parse().unwrap();
    let source = Chain { n, family };
    let prepare = Instant::now();
    let operator = DiracOperator::new(&source, 0..=1)?;
    let prepare_ns = prepare.elapsed().as_nanos();
    let policy = IterativeOptions {
        block_size: block,
        ..Default::default()
    };
    let begin = Instant::now();
    let result = operator.boundary_spectrum(1, count, 0., &policy, &Execution::default());
    let solve_ns = begin.elapsed().as_nanos();
    match result {
        Ok(pairs) => println!(
            "{{\"status\":\"success\",\"prepare_ns\":{prepare_ns},\"solve_ns\":{solve_ns},\"values\":{:?},\"residuals\":{:?},\"iterations\":{},\"actions\":{},\"peak_subspace\":{},\"index_coverage_certified\":{}}}",
            pairs.singular_values(),
            pairs.residuals(),
            pairs.iterations(),
            pairs.gram_applications(),
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
