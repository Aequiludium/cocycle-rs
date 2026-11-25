//! Standalone numerical path certificate worker; no runtime dependencies.
//! Input: n d max_paths max_dense_entries tolerance source_event target_event,
//! followed by directed weighted edge triples. Optional --non-regular.
#![forbid(unsafe_code)]
use cocycle::complex::PathBoundaryMode;
use cocycle::hodge::compute_persistent_path_laplacian;
use std::error::Error;
use std::io::{self, Read};

fn main() -> Result<(), Box<dyn Error>> {
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args.iter().any(|a| a != "--non-regular") {
        return Err("unknown option".into());
    }
    let mode = if args.is_empty() {
        PathBoundaryMode::Regular
    } else {
        PathBoundaryMode::NonRegular
    };
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let t: Vec<_> = input.split_whitespace().collect();
    if t.len() < 7 || !(t.len() - 7).is_multiple_of(3) {
        return Err("expected n d paths dense tolerance source target then edge triples".into());
    }
    let edges = t[7..]
        .as_chunks::<3>()
        .0
        .iter()
        .map(|e| Ok((e[0].parse()?, e[1].parse()?, e[2].parse()?)))
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    let operator = compute_persistent_path_laplacian(
        t[0].parse()?,
        &edges,
        t[1].parse()?,
        (t[5].parse()?, t[6].parse()?),
        mode,
        (t[2].parse()?, t[3].parse()?),
        t[4].parse()?,
    )?;
    let mut columns = Vec::new();
    for i in 0..operator.basis().len() {
        let mut v = vec![0.0; operator.basis().len()];
        v[i] = 1.0;
        columns.push(operator.apply(&v)?);
    }
    println!(
        "{{\"protocol\":\"cocycle-path-real-certificate-v1\",\"field\":\"real\",\"regular\":{},\"dimension\":{},\"source_time\":{},\"target_time\":{},\"paths\":{:?},\"basis\":{:?},\"laplacian_columns\":{:?}}}",
        mode == PathBoundaryMode::Regular,
        operator.dimension(),
        operator.times().0,
        operator.times().1,
        operator.paths(),
        operator.basis().collect::<Vec<_>>(),
        columns
    );
    Ok(())
}
