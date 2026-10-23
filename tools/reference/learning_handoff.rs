//! Optional #92 supplied-graph/stage exporter; no PH or learning runtime.
//! Whitespace protocol: key snapshot order n m g, g scales, n (id birth x y),
//! then m (a b birth). Multiple samples may follow. Keys use ASCII word chars.
//! Trusted fixture tool; rejects input above 1 MiB or 128 samples in total.

use std::collections::HashMap;
use std::io::{self, Read};
use std::str::{FromStr, SplitWhitespace};

use cocycle::complex::{FilteredComplex, Simplex, SimplicialComplex};
use cocycle::filtration::SimplicialStage;

const MAX_INPUT_BYTES: u64 = 1024 * 1024;
const MAX_SAMPLES: usize = 128;

fn read<T: FromStr>(words: &mut SplitWhitespace<'_>) -> Result<T, Box<dyn std::error::Error>> {
    words
        .next()
        .ok_or_else(|| "incomplete sample".into())
        .and_then(|word| word.parse().map_err(|_| "invalid value".into()))
}

fn key(word: &str) -> Result<&str, Box<dyn std::error::Error>> {
    if word.is_empty()
        || !word
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'-' || b == b'_')
    {
        return Err("keys must use ASCII letters, digits, '-' or '_'".into());
    }
    Ok(word)
}

fn export(
    source: &SimplicialStage<'_>,
    features: &HashMap<usize, [f64; 2]>,
) -> Result<String, Box<dyn std::error::Error>> {
    let vertices: Vec<_> = source
        .cells()
        .filter(|&id| source.dimension(id) == 0)
        .collect();
    let edges: Vec<_> = source
        .cells()
        .filter(|&id| source.dimension(id) == 1)
        .collect();
    // A stage carries its own owner; callers cannot supply a different owner.
    let vertex_keys = vertices
        .iter()
        .map(|&id| {
            source
                .simplex(id)
                .map(Simplex::vertices)
                .ok_or("missing stage vertex")
        })
        .collect::<Result<Vec<_>, _>>()?;
    let edge_keys = edges
        .iter()
        .map(|&id| {
            source
                .simplex(id)
                .map(Simplex::vertices)
                .ok_or("missing stage edge")
        })
        .collect::<Result<Vec<_>, _>>()?;
    let values = vertex_keys
        .iter()
        .map(|v| {
            let vertex = v.first().ok_or("empty vertex key")?;
            features
                .get(vertex)
                .copied()
                .ok_or("missing supplied vertex features")
        })
        .collect::<Result<Vec<_>, _>>()?;
    if values.iter().flatten().any(|x| !x.is_finite()) {
        return Err("nonfinite supplied vertex features".into());
    }
    let mut terms = Vec::new();
    for (column, &edge) in edges.iter().enumerate() {
        for term in source.boundary(edge).ok_or("missing stage boundary")? {
            let row = vertices
                .iter()
                .position(|&id| id == term.face)
                .ok_or("boundary face absent from vertex basis")?;
            let sign = term.coefficient;
            terms.push(format!("[{row},{column},{sign}]"));
        }
    }
    Ok(format!(
        "\"bases\":{{\"0\":{vertex_keys:?},\"1\":{edge_keys:?},\"2\":[]}},\
         \"boundary\":{{\"shape\":[{},{}],\"terms\":[{}]}},\
         \"features\":{{\"shape\":[{},2],\"row_keys\":{vertex_keys:?},\
         \"values\":{values:?},\"meaning\":\"supplied-coordinate\",\"units\":\"dimensionless\"}}",
        vertices.len(),
        edges.len(),
        terms.join(","),
        vertices.len()
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin()
        .take(MAX_INPUT_BYTES + 1)
        .read_to_string(&mut input)?;
    if input.len() as u64 > MAX_INPUT_BYTES {
        return Err("fixture input exceeds 1 MiB".into());
    }
    let mut words = input.split_whitespace();
    let mut sample_count = 0;
    while let Some(sample) = words.next() {
        if sample_count == MAX_SAMPLES {
            return Err("fixture input exceeds 128 samples".into());
        }
        sample_count += 1;
        let sample = key(sample)?;
        let snapshot = key(words.next().ok_or("missing snapshot")?)?;
        let order: usize = read(&mut words)?;
        let n: usize = read(&mut words)?;
        let m: usize = read(&mut words)?;
        let g: usize = read(&mut words)?;
        // This is a bounded fixture tool, not a general-purpose file parser.
        if n > 64 || m > 2016 || g == 0 || g > 16 {
            return Err("fixture limits: 64 vertices, 2016 edges, 1..16 scales".into());
        }
        let grid = (0..g)
            .map(|_| read(&mut words))
            .collect::<Result<Vec<f64>, _>>()?;
        if grid.iter().any(|x| !x.is_finite()) || grid.windows(2).any(|x| x[0] >= x[1]) {
            return Err("scales must be finite and strictly increasing".into());
        }
        let mut cells = Vec::new();
        let mut features = HashMap::new();
        for _ in 0..n {
            let id: usize = read(&mut words)?;
            let birth = read(&mut words)?;
            let x: f64 = read(&mut words)?;
            let y: f64 = read(&mut words)?;
            if !x.is_finite() || !y.is_finite() || features.insert(id, [x, y]).is_some() {
                return Err("duplicate vertex or nonfinite feature".into());
            }
            cells.push(Simplex::new(vec![id], birth)?);
        }
        for _ in 0..m {
            let a = read(&mut words)?;
            let b = read(&mut words)?;
            let birth = read(&mut words)?;
            cells.push(Simplex::new(vec![a, b], birth)?);
        }
        let owner = SimplicialComplex::new(cells)?;
        let full = owner.stage(owner.max_filtration_value().unwrap_or(0.))?;
        let common = format!(
            "\"format\":\"cocycle-graph-learning-v1\",\"sample_key\":\"{sample}\",\
             \"snapshot\":\"{snapshot}\",\"sample_order\":{order},\
             \"domain\":\"signed-integer\",\"scale_unit\":\"fixture-time\",\
             \"topology\":\"supplied-graph-no-clique-expansion\""
        );
        println!(
            "{{{common},\"selection\":{{\"kind\":\"topology\",\"index\":null,\"scale\":null}},\
             \"coverage\":\"complete-supplied-source\",{}}}",
            export(&full, &features)?
        );
        for (index, &scale) in grid.iter().enumerate() {
            let stage = owner.stage(scale)?;
            println!(
                "{{{common},\"selection\":{{\"kind\":\"stage\",\"index\":{index},\"scale\":{scale}}},\
                 \"coverage\":\"inclusive-supplied-sublevel\",{}}}",
                export(&stage, &features)?
            );
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn missing_or_nonfinite_features_return_errors() -> Result<(), Box<dyn std::error::Error>> {
        let owner = SimplicialComplex::new(vec![Simplex::new(vec![10], 0.)?])?;
        let stage = owner.stage(0.)?;
        assert!(export(&stage, &HashMap::new()).is_err());
        assert!(export(&stage, &HashMap::from([(10, [f64::NAN, 0.])])).is_err());
        Ok(())
    }

    #[test]
    fn empty_stage_needs_no_vertex_features() -> Result<(), Box<dyn std::error::Error>> {
        let owner = SimplicialComplex::new(vec![Simplex::new(vec![10], 0.)?])?;
        let output = export(&owner.stage(-1.)?, &HashMap::new())?;
        assert!(output.contains("\"shape\":[0,2]"));
        Ok(())
    }
}
