//! Explicit opt-in replay of the production LSMR source with iteration tracing.
#![allow(dead_code)]
#[path = "../../src/execution/mod.rs"]
mod execution;
#[path = "../../src/persistent_laplacian/lsmr.rs"]
mod lsmr;
use cocycle::{Error, Result};
use execution::Execution;
use lsmr::{LsmrOptions, LsmrStep};
use std::io::{self, Read};

fn values(name: &str, entries: &[f64]) {
    print!("{name}");
    for value in entries {
        print!(" {value:.17e}");
    }
    println!();
}
fn step(record: LsmrStep<'_>) {
    println!(
        "STEP {} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}",
        record.iteration, record.alpha, record.beta, record.rho, record.rhobar, record.zetabar
    );
    values("U", record.u);
    values("V", record.v);
    values("X", record.solution);
    values("R", record.residual);
    values("NORMAL", record.normal_residual);
}
fn next<T: std::str::FromStr>(
    words: &mut std::str::SplitWhitespace<'_>,
) -> std::result::Result<T, Box<dyn std::error::Error>> {
    words
        .next()
        .ok_or("missing LSMR input")?
        .parse()
        .map_err(|_| "invalid LSMR input".into())
}
fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let mut data = String::new();
    io::stdin().read_to_string(&mut data)?;
    let mut words = data.split_whitespace();
    while let Some(label) = words.next() {
        let rows: usize = next(&mut words)?;
        let columns: usize = next(&mut words)?;
        let tolerance: f64 = next(&mut words)?;
        let limit: isize = next(&mut words)?;
        let entries = rows.checked_mul(columns).ok_or("input shape overflow")?;
        if entries > 1_048_576 || rows.max(columns) > 65_536 {
            return Err("worker input cap".into());
        }
        let boundary: Vec<f64> = (0..entries)
            .map(|_| next(&mut words))
            .collect::<std::result::Result<_, _>>()?;
        let y: Vec<f64> = (0..columns)
            .map(|_| next(&mut words))
            .collect::<std::result::Result<_, _>>()?;
        let options = LsmrOptions {
            relative_tolerance: tolerance,
            max_iterations: if limit == -1 {
                None
            } else {
                Some(usize::try_from(limit)?)
            },
            max_matrix_entries: 1_048_576,
            ..Default::default()
        };
        println!("CASE {label}");
        let result = lsmr::replay(
            &boundary,
            rows,
            columns,
            &y,
            options,
            &Execution::default(),
            &mut step,
        );
        match result {
            Ok(result) => {
                println!(
                    "RESULT {} {:.17e} {:.17e}",
                    result.iterations, result.residual_norm, result.normal_residual_norm
                );
                values("SOLUTION", &result.solution);
                values("PROJECTED", &result.projected);
            }
            Err(error) => println!("ERROR {error:?}"),
        }
        println!("END");
    }
    Ok(())
}
