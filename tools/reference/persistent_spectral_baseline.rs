//! Optional #111 reproduction worker using existing Betti/lifetime descriptors.
//! Input: one four-scale grid line, then `degree birth death` (inf = essential).
//! Output: a JSON feature array; no spectral or training dependency is introduced.

use std::io::{self, Read};

use cocycle::descriptors::{betti_curve, finite_lifetime_summary};
use cocycle::diagram::{Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let mut lines = input.lines();
    let grid = lines
        .next()
        .ok_or("missing grid")?
        .split_whitespace()
        .map(str::parse)
        .collect::<Result<Vec<f64>, _>>()?;
    if grid.len() != 4 {
        return Err("expected four grid points".into());
    }
    let mut intervals = Vec::new();
    for line in lines.filter(|line| !line.trim().is_empty()) {
        let fields: Vec<_> = line.split_whitespace().collect();
        if fields.len() != 3 {
            return Err("expected degree birth death".into());
        }
        let degree = fields[0].parse()?;
        let birth = fields[1].parse()?;
        let end = if fields[2] == "inf" {
            IntervalEnd::Essential
        } else {
            IntervalEnd::Finite(fields[2].parse()?)
        };
        intervals.push(PersistenceInterval::new(degree, birth, end)?);
    }
    let diagram = PersistenceDiagram::new(1, Coverage::Complete, intervals)?;
    let mut features = Vec::new();
    for degree in 0..=1 {
        features.extend(
            betti_curve(&diagram, degree, &grid)?
                .iter()
                .map(|&n| n as f64),
        );
        let summary = finite_lifetime_summary(&diagram, degree)?;
        features.extend([
            summary.finite_count() as f64,
            summary.total_persistence(),
            summary.max_persistence().unwrap_or(0.),
            summary.entropy().unwrap_or(0.),
            summary.excluded_essential_count() as f64,
            f64::from(summary.max_persistence().is_some()),
            f64::from(summary.entropy().is_some()),
        ]);
    }
    println!("{features:?}");
    Ok(())
}
