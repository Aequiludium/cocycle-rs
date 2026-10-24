//! Optional #111 worker using existing descriptors and real persistent spectra.
//! Input: one four-scale grid line, then `degree birth death` (inf = essential).
//! `--spectral`: sample/source IDs, four scales, cell count, then birth/arity/vertices.
//! Output: JSON; no external spectral or training dependency is introduced.

use std::fmt::Write as _;
use std::io::{self, Read};

use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::descriptors::{betti_curve, finite_lifetime_summary};
use cocycle::diagram::{Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval};
use cocycle::execution::Execution;
use cocycle::persistent_laplacian::{PersistentLaplacian, SpectralOptions};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin().take(65_537).read_to_string(&mut input)?;
    if input.len() > 65_536 {
        return Err("input exceeds 64 KiB".into());
    }
    let args: Vec<_> = std::env::args().skip(1).collect();
    if args == ["--spectral"] {
        println!("{}", spectral_features(&input)?);
        return Ok(());
    }
    if !args.is_empty() {
        return Err("usage: persistent_spectral_baseline [--spectral]".into());
    }
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

// Application-owned, bounded export of #124/#109 owned results. Raw spectra and
// pair bases stay alongside the truncated features; there is no new core API.
fn spectral_features(input: &str) -> Result<String, Box<dyn std::error::Error>> {
    let mut words = input.split_whitespace();
    let sample = words.next().ok_or("missing sample ID")?;
    let source_id = words.next().ok_or("missing source ID")?;
    for id in [sample, source_id] {
        if id.len() > 128
            || !id
                .bytes()
                .all(|c| c.is_ascii_alphanumeric() || b"-_.".contains(&c))
        {
            return Err(
                "IDs require 1..128 ASCII letters, digits, hyphens, dots or underscores".into(),
            );
        }
    }
    let mut grid = [0_f64; 4];
    for scale in &mut grid {
        *scale = next(&mut words)?;
    }
    if grid.iter().any(|x| !x.is_finite()) || grid.windows(2).any(|x| x[0] >= x[1]) {
        return Err("require four finite strictly increasing scales".into());
    }
    let count: usize = next(&mut words)?;
    if count > 256 {
        return Err("at most 256 supplied simplices".into());
    }
    let mut cells = Vec::new();
    for _ in 0..count {
        let birth = next(&mut words)?;
        let length: usize = next(&mut words)?;
        if !(1..=3).contains(&length) {
            return Err("only supplied vertices, edges and triangles are supported".into());
        }
        let vertices = (0..length)
            .map(|_| next(&mut words))
            .collect::<Result<Vec<usize>, _>>()?;
        cells.push(Simplex::new(vertices, birth)?);
    }
    if words.next().is_some() {
        return Err("unexpected trailing input".into());
    }
    let source = SimplicialComplex::new(cells)?;
    let stages = grid
        .iter()
        .map(|&scale| source.stage(scale))
        .collect::<cocycle::Result<Vec<_>>>()?;
    let options = SpectralOptions {
        max_dense_entries: 4096,
        ..Default::default()
    };
    let execution = Execution::default().max_work(20_000_000);
    let mut records = Vec::new();
    let mut values = Vec::new();
    let mut mask = Vec::new();
    for degree in 0..=1 {
        for i in 0..4 {
            for j in i..4 {
                let operator = PersistentLaplacian::compute(
                    &stages[i], &stages[j], degree, options, &execution,
                )?;
                let spectrum = operator.spectrum(&execution)?;
                // An absolute application policy, separate from the producer's
                // suggested threshold. Keep all signed raw eigenvalues below.
                if spectrum.eigenvalues.iter().any(|&x| x < -1e-10) {
                    return Err("negative eigenvalue exceeds feature zero tolerance".into());
                }
                for rank in 0..8 {
                    let value = spectrum.eigenvalues.get(rank).copied();
                    values.push(value.map_or(0., |x| if x.abs() <= 1e-10 { 0. } else { x }));
                    mask.push(value.is_some());
                }
                let basis: Vec<_> = operator.basis().iter().map(Simplex::vertices).collect();
                let births: Vec<_> = operator.basis().iter().map(Simplex::value).collect();
                let target_basis: Vec<_> = operator
                    .target_basis()
                    .iter()
                    .map(Simplex::vertices)
                    .collect();
                let mut record = String::new();
                write!(
                    record,
                    "{{\"degree\":{degree},\"i\":{i},\"j\":{j},\"sample_id\":{sample:?},\"source_id\":{source_id:?},\"status\":\"computed\",\"coverage\":\"full\",\"domain_size\":{},\"basis\":{basis:?},\"basis_births\":{births:?},\"target_basis\":{target_basis:?},\"values\":{:?},\"matrix\":{:?},\"solver_residual\":{},\"eigen_residual\":{},\"reference_scale\":{},\"suggested_zero_threshold\":{}}}",
                    basis.len(),
                    spectrum.eigenvalues,
                    operator.matrix(),
                    operator.solver_residual(),
                    spectrum.residual,
                    operator.reference_scale(),
                    spectrum.zero_threshold
                )?;
                records.push(record);
            }
        }
    }
    Ok(format!(
        "{{\"sample_id\":{sample:?},\"source_id\":{source_id:?},\"grid\":{grid:?},\"operator\":\"persistent_laplacian\",\"field\":\"R\",\"inner_product\":\"identity\",\"complex_type\":\"simplicial\",\"scale_unit\":\"synthetic_scale\",\"eigenvalue_unit\":\"dimensionless\",\"source_coverage\":\"complete_supplied_complex\",\"numeric_origin\":{{\"producer\":\"cocycle.persistent_laplacian\",\"solver\":\"cyclic_jacobi\",\"dtype\":\"float64\",\"schur_rcond\":1e-12,\"zero_tolerance\":1e-10}},\"features\":{values:?},\"validity_mask\":{mask:?},\"spectra\":[{}]}}",
        records.join(",")
    ))
}

fn next<T: std::str::FromStr>(
    words: &mut std::str::SplitWhitespace<'_>,
) -> Result<T, Box<dyn std::error::Error>>
where
    T::Err: std::error::Error + 'static,
{
    Ok(words.next().ok_or("incomplete input")?.parse()?)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_source_exports_computed_empty_spectra() {
        let output = spectral_features("empty supplied 0 1 2 3 0").unwrap();
        assert_eq!(output.matches("\"domain_size\":0").count(), 20);
        assert_eq!(output.matches("false").count(), 160);
    }

    #[test]
    fn rejects_invalid_identity_grid_cells_and_excess_input() {
        for input in [
            "bad\"id supplied 0 1 2 3 0",
            "s source 0 0 1 2 0",
            "s source 0 1 2 NaN 0",
            "s source 0 1 2 3 257",
            "s source 0 1 2 3 1 0 4",
            "s source 0 1 2 3 1 0 2 0 1", // missing faces
            "s source 0 1 2 3 0 trailing",
        ] {
            assert!(spectral_features(input).is_err(), "{input}");
        }
    }
}
