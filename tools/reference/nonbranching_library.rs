//! Frozen-matrix replay of the production weak-reduction source with its real
//! budget and dense allocator. Public stage APIs are checked by the ordinary
//! persistent_laplacian example and integration tests, not this matrix protocol.
#![forbid(unsafe_code)]

use cocycle::persistent_laplacian::SpectralOptions;
use cocycle::{Error, Result};
use std::io::{self, Read};

// These source owners also contain operations unused by this standalone worker.
#[allow(dead_code)]
#[path = "../../src/execution/mod.rs"]
mod execution;
#[path = "../../src/persistent_laplacian/nonbranching.rs"]
mod nonbranching;
#[allow(dead_code)]
#[path = "../../src/persistent_laplacian/schur.rs"]
mod schur;

fn allocation() -> Error {
    Error::AllocationFailed {
        context: "persistent Laplacian",
    }
}
fn numerical(context: &'static str) -> Error {
    Error::NumericalFailure { context }
}

fn main() -> std::result::Result<(), Box<dyn std::error::Error>> {
    let mut options = SpectralOptions::default();
    let mut execution = execution::Execution::default();
    let mut args = std::env::args().skip(1);
    while let Some(arg) = args.next() {
        let value = args.next().ok_or("missing option value")?;
        match arg.as_str() {
            "--max-dense-entries" => options.max_dense_entries = value.parse()?,
            "--max-work" => execution = execution.max_work(value.parse()?),
            _ => return Err("unknown option".into()),
        }
    }
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let mut words = input.split_whitespace();
    let mut dimension = || -> std::result::Result<usize, Box<dyn std::error::Error>> {
        Ok(words.next().ok_or("missing dimension")?.parse()?)
    };
    let (n, m, p) = (dimension()?, dimension()?, dimension()?);
    // Match the frozen standalone protocol, including zero-product shapes.
    if n > 65_536 || m > 65_536 || p > 65_536 {
        return Err("dimension limit exceeded".into());
    }
    let mut budget = execution::WorkBudget::new(&execution)?;
    let mut b1 = schur::dense(n, p, &options, &mut budget)?;
    let mut d = schur::dense(m, p, &options, &mut budget)?;
    for value in b1.iter_mut().chain(&mut d) {
        *value = words.next().ok_or("missing matrix entry")?.parse()?;
    }
    if words.next().is_some() {
        return Err("trailing input".into());
    }
    let reduction = nonbranching::weak_reduction(&d, m, p, &options, &mut budget)?;
    let (z, metric) = reduction.kernel(p, &options, &mut budget)?;
    let upper = nonbranching::upper(&b1, &d, n, m, p, &options, &mut budget)?;
    let c = reduction.nonzero.len();
    let mut e = schur::dense(c, c, &options, &mut budget)?;
    for (i, &value) in reduction.e.iter().enumerate() {
        e[i * c + i] = value;
    }
    println!(
        "{{\"R\":{:?},\"E\":{:?},\"V\":{:?},\"Z\":{:?},\"nonzero\":{:?},\"original_zero\":{:?},\"transformed_zero\":{:?},\"metric\":{:?},\"upper\":{:?}}}",
        reduction.r,
        e,
        reduction.v,
        z,
        reduction.nonzero,
        reduction.original_zero,
        reduction.zero_columns,
        metric,
        upper
    );
    Ok(())
}
