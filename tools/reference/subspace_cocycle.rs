//! Standalone Issue 148 filter parity and private supplied-block LOBPCG worker.
//! Build with rustc against the current release library; no product FFT API.
pub use cocycle::{Error, Result};
#[allow(dead_code)]
#[path = "../../src/eigensolver.rs"]
mod eigensolver;
#[allow(dead_code)]
#[path = "../../src/execution/mod.rs"]
mod execution;

fn canonical_zero(value: f64) -> f64 {
    if value == 0. { 0. } else { value }
}

fn read<T: std::str::FromStr>(words: &mut std::str::SplitWhitespace<'_>) -> T
where
    T::Err: std::fmt::Debug,
{
    words
        .next()
        .expect("missing fixture field")
        .parse()
        .expect("fixture field")
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let text = std::fs::read_to_string(std::env::args().nth(1).ok_or("fixture path")?)?;
    let mut words = text.split_whitespace();
    let mode = words.next().ok_or("mode")?;
    let n: usize = read(&mut words);
    let b: usize = read(&mut words);
    if n == 0 || n > 256 || b == 0 || b > n {
        return Err("control shape: require 0 < b <= n <= 256".into());
    }
    if mode == "filter" {
        let lower: f64 = read(&mut words);
        let upper: f64 = read(&mut words);
        let anchor: f64 = read(&mut words);
        let degrees: Vec<usize> = (0..b).map(|_| read(&mut words)).collect();
        if !anchor.is_finite()
            || !lower.is_finite()
            || !upper.is_finite()
            || !(anchor < lower && lower < upper)
            || degrees.iter().any(|&m| m > 64)
        {
            return Err("filter bounds/degree".into());
        }
        let a: Vec<f64> = (0..n * n).map(|_| read(&mut words)).collect();
        let y: Vec<f64> = (0..n * b).map(|_| read(&mut words)).collect();
        if a.iter().chain(&y).any(|v| !v.is_finite()) || words.next().is_some() {
            return Err("filter values/trailing fields".into());
        }
        let center = (lower + upper) / 2.;
        let radius = (upper - lower) / 2.;
        let sigma1 = radius / (anchor - center);
        let mut output = y.clone();
        for column in 0..b {
            let mut previous: Vec<f64> = (0..n).map(|i| y[i * b + column]).collect();
            let mut current = previous.clone();
            let mut sigma = sigma1;
            for degree in 0..degrees[column] {
                let mut next = vec![0.; n];
                let next_sigma = if degree == 0 {
                    sigma1
                } else {
                    1. / (2. / sigma1 - sigma)
                };
                for row in 0..n {
                    let image: f64 = (0..n).map(|j| a[row * n + j] * current[j]).sum();
                    let image = (image - center * current[row]) / radius;
                    next[row] = if degree == 0 {
                        sigma1 * image
                    } else {
                        2. * next_sigma * image - next_sigma * sigma * previous[row]
                    };
                }
                previous = current;
                current = next;
                sigma = next_sigma;
            }
            for row in 0..n {
                output[row * b + column] = current[row];
            }
        }
        if output.iter().any(|v| !v.is_finite()) {
            return Err("nonfinite filter output".into());
        }
        println!("{{\"filtered\":{output:?}}}");
    } else if mode == "lobpcg" {
        let count: usize = read(&mut words);
        let warm: usize = read(&mut words);
        let a: Vec<f64> = (0..n * n).map(|_| read(&mut words)).collect();
        let initial: Vec<f64> = if warm == 1 {
            (0..n * b).map(|_| read(&mut words)).collect()
        } else {
            Vec::new()
        };
        if warm > 1 || words.next().is_some() || a.iter().any(|v| !v.is_finite()) {
            return Err("LOBPCG fixture".into());
        }
        let options = eigensolver::IterativeOptions {
            block_size: b,
            absolute_tolerance: 1e-11,
            relative_tolerance: 1e-11,
            ..Default::default()
        };
        let result = eigensolver::lobpcg_with_initial(
            n,
            count,
            false,
            &options,
            (warm == 1).then_some(initial.as_slice()),
            &mut execution::WorkBudget::new(&execution::Execution::default())?,
            |x, budget| {
                budget.step_by(n * n)?;
                Ok((0..n)
                    .map(|i| (0..n).map(|j| a[i * n + j] * x[j]).sum())
                    .collect())
            },
        )?;
        println!(
            "{{\"values\":{:?},\"vectors\":{:?},\"residuals\":{:?},\"actions\":{},\"iterations\":{}}}",
            result.values, result.vectors, result.residuals, result.actions, result.iterations
        );
    } else {
        return Err("unknown control mode".into());
    }
    Ok(())
}
