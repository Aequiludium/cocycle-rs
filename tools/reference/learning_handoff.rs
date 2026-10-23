//! Optional #92 supplied-graph/stage exporter; no PH or learning runtime.
//! Whitespace protocol: key snapshot order n m g, g scales, n (id birth x y),
//! then m (a b birth). Multiple samples may follow. Keys use ASCII word chars.

use std::collections::HashMap;
use std::io::{self, Read};
use std::str::{FromStr, SplitWhitespace};

use cocycle::complex::{FilteredComplex, Simplex, SimplexId, SimplicialComplex};

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
    source: &impl FilteredComplex<CellId = SimplexId>,
    owner: &SimplicialComplex,
    features: &HashMap<usize, [f64; 2]>,
) -> String {
    let vertices: Vec<_> = source
        .cells()
        .filter(|&id| source.dimension(id) == 0)
        .collect();
    let edges: Vec<_> = source
        .cells()
        .filter(|&id| source.dimension(id) == 1)
        .collect();
    // These handles all belong to owner, including when source is its stage.
    let vertex_keys: Vec<_> = vertices
        .iter()
        .map(|&id| owner.simplex(id).unwrap().vertices())
        .collect();
    let edge_keys: Vec<_> = edges
        .iter()
        .map(|&id| owner.simplex(id).unwrap().vertices())
        .collect();
    let values: Vec<_> = vertex_keys.iter().map(|v| features[&v[0]]).collect();
    let mut terms = Vec::new();
    for (column, &edge) in edges.iter().enumerate() {
        for (face, sign) in source.boundary(edge) {
            let row = vertices.iter().position(|&id| id == face).unwrap();
            terms.push(format!("[{row},{column},{sign}]"));
        }
    }
    format!(
        "\"bases\":{{\"0\":{vertex_keys:?},\"1\":{edge_keys:?},\"2\":[]}},\
         \"boundary\":{{\"shape\":[{},{}],\"terms\":[{}]}},\
         \"features\":{{\"shape\":[{},2],\"row_keys\":{vertex_keys:?},\
         \"values\":{values:?},\"meaning\":\"supplied-coordinate\",\"units\":\"dimensionless\"}}",
        vertices.len(),
        edges.len(),
        terms.join(","),
        vertices.len()
    )
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mut input = String::new();
    io::stdin().read_to_string(&mut input)?;
    let mut words = input.split_whitespace();
    while let Some(sample) = words.next() {
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
        let common = format!(
            "\"format\":\"cocycle-graph-learning-v1\",\"sample_key\":\"{sample}\",\
             \"snapshot\":\"{snapshot}\",\"sample_order\":{order},\
             \"domain\":\"signed-integer\",\"scale_unit\":\"fixture-time\",\
             \"topology\":\"supplied-graph-no-clique-expansion\""
        );
        println!(
            "{{{common},\"selection\":{{\"kind\":\"topology\",\"index\":null,\"scale\":null}},\
             \"coverage\":\"complete-supplied-source\",{}}}",
            export(&owner, &owner, &features)
        );
        for (index, &scale) in grid.iter().enumerate() {
            let stage = owner.stage(scale)?;
            println!(
                "{{{common},\"selection\":{{\"kind\":\"stage\",\"index\":{index},\"scale\":{scale}}},\
                 \"coverage\":\"inclusive-supplied-sublevel\",{}}}",
                export(&stage, &owner, &features)
            );
        }
    }
    Ok(())
}
