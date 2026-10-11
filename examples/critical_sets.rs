//! One frozen-target lower-star simplification step; no optimizer dependency.
use cocycle::Result;
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::optimization::CriticalSetWorkspace;

const EDGES: [[usize; 2]; 4] = [[0, 1], [1, 2], [2, 3], [0, 3]];

fn lower_star(values: &[f64; 4]) -> Result<SimplicialComplex> {
    let mut cells = Vec::new();
    for (i, &value) in values.iter().enumerate() {
        cells.push(Simplex::new(vec![i], value)?);
    }
    for vertices in EDGES {
        cells.push(Simplex::new(
            vertices.to_vec(),
            values[vertices[0]].max(values[vertices[1]]),
        )?);
    }
    SimplicialComplex::new(cells)
}

fn step(values: &[f64; 4]) -> Result<[f64; 4]> {
    let source = lower_star(values)?;
    let mut workspace = CriticalSetWorkspace::new(&source)?;
    let mut endpoints = Vec::new();
    // Collapse short finite H0 points toward the diagonal. This fixture has
    // bounded values; a general caller must also check midpoint arithmetic.
    for (birth, death) in workspace.finite_pairs() {
        let b = source.simplex(birth).unwrap();
        let d = source.simplex(death).unwrap();
        if b.dimension() == 0 && d.value() > b.value() && d.value() - b.value() <= 1.5 {
            let midpoint = (b.value() + d.value()) / 2.;
            endpoints.extend([(birth, midpoint), (death, midpoint)]);
        }
    }
    let mut gradient = [0.; 4];
    for (id, target) in workspace.targets(&endpoints)? {
        let simplex = source.simplex(id).unwrap(); // IDs belong to this source.
        // Lower-star chain rule: select the first maximum vertex at a tie.
        // This is a subgradient convention, not a smooth derivative at ties.
        let vertex = *simplex
            .vertices()
            .iter()
            .find(|&&i| values[i] == simplex.value())
            .unwrap();
        gradient[vertex] += 2. * (simplex.value() - target);
    }
    let mut next = *values;
    for i in 0..4 {
        next[i] -= 0.25 * gradient[i];
    }
    // Rebuild from data: directly changing just the returned simplex values
    // could violate face <= coface. All IDs must be looked up again afterwards.
    lower_star(&next)?;
    Ok(next)
}

fn main() -> Result<()> {
    let values = [-2., 1., -1., 0.];
    let next = step(&values)?;
    println!("Vertex values: {values:?} -> {next:?}");
    for (label, values) in [("before", values), ("after", next)] {
        let source = lower_star(&values)?;
        let workspace = CriticalSetWorkspace::new(&source)?;
        let finite: Vec<_> = workspace
            .finite_pairs()
            .filter_map(|(b, d)| {
                let b = source.simplex(b).unwrap();
                let d = source.simplex(d).unwrap();
                (d.value() > b.value()).then_some((b.dimension(), b.value(), d.value()))
            })
            .collect();
        println!("{label}: positive-length finite intervals {finite:?}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn frozen_step_reduces_the_hand_derived_short_h0_interval() -> Result<()> {
        let next = step(&[-2., 1., -1., 0.])?;
        assert_eq!(next, [-2., 1., -0.75, -0.5]);
        let source = lower_star(&next)?;
        let workspace = CriticalSetWorkspace::new(&source)?;
        let lengths: Vec<_> = workspace
            .finite_pairs()
            .filter_map(|(b, d)| {
                let length =
                    source.simplex(d).unwrap().value() - source.simplex(b).unwrap().value();
                (length > 0.).then_some(length)
            })
            .collect();
        assert_eq!(lengths, [0.25]);
        Ok(())
    }
}
