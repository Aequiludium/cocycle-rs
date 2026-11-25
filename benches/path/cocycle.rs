//! Native worker for matched ordinary H0/H1 path-diagram terminals.
#![forbid(unsafe_code)]
use cocycle::algebra::PrimeField;
use cocycle::complex::PathBoundaryMode;
use cocycle::diagram::IntervalEnd;
use cocycle::execution::Execution;
use cocycle::persistence::{
    compute_low_dimensional_path_persistence_with, compute_path_persistence_with,
};
use std::error::Error;
use std::io::{self, Read};

fn main() -> Result<(), Box<dyn Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let t: Vec<_> = input.split_whitespace().collect();
    if t.len() < 7 || !(t.len() - 7).is_multiple_of(3) {
        return Err("expected engine n field mode paths entries work then edge triples".into());
    }
    let n = t[1].parse()?;
    let field = PrimeField::new(t[2].parse()?)?;
    let mode = match t[3] {
        "regular" => PathBoundaryMode::Regular,
        "full" => PathBoundaryMode::NonRegular,
        _ => return Err("invalid model".into()),
    };
    let paths = t[4].parse()?;
    let entries = t[5].parse()?;
    let execution = Execution::default().max_work(t[6].parse()?);
    let edges = t[7..]
        .as_chunks::<3>()
        .0
        .iter()
        .map(|e| Ok((e[0].parse()?, e[1].parse()?, e[2].parse()?)))
        .collect::<Result<Vec<_>, Box<dyn Error>>>()?;
    // Parsing/field validation are shared preparation; each public terminal
    // includes graph validation, its actual algorithm and canonical assembly.
    let started = std::time::Instant::now();
    let result = match t[0] {
        "generic" => {
            compute_path_persistence_with(n, &edges, 1, field, mode, paths, entries, &execution)?
        }
        "direct" => compute_low_dimensional_path_persistence_with(
            n, &edges, field, mode, paths, entries, &execution,
        )?,
        _ => return Err("invalid engine".into()),
    };
    let mut bars = Vec::new();
    for bar in result.diagram().intervals() {
        let death = match bar.end() {
            IntervalEnd::Finite(x) => Some(x),
            IntervalEnd::Essential => None,
            _ => return Err("unexpected censored interval".into()),
        };
        bars.push((bar.dimension(), bar.birth(), death));
    }
    drop(result); // Equal payload export and result/workspace destruction are timed.
    let elapsed = started.elapsed().as_nanos();
    print!(
        "{{\"protocol\":\"cocycle-path-diagram-v1\",\"engine\":\"{}\",\"field\":{},\"mode\":\"{}\",\"dimension\":1,\"coverage\":\"complete\",\"elapsed_ns\":{},\"intervals\":[",
        t[0],
        field.characteristic(),
        t[3],
        elapsed
    );
    for (i, (q, b, d)) in bars.iter().enumerate() {
        if i > 0 {
            print!(",");
        }
        print!("[{q},{b},");
        if let Some(d) = d {
            print!("{d}");
        } else {
            print!("null");
        }
        print!("]");
    }
    println!("]}}");
    Ok(())
}
