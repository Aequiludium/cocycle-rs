//! A single-scale class query and finite transport from a ring to its filling.
use cocycle::Result;
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::homology_operators::OperatorFamily;

fn main() -> Result<()> {
    let source = SimplicialComplex::new(vec![
        Simplex::new(vec![10], 0.)?,
        Simplex::new(vec![20], 0.)?,
        Simplex::new(vec![30], 0.)?,
        Simplex::new(vec![10, 20], 1.)?,
        Simplex::new(vec![10, 30], 1.)?,
        Simplex::new(vec![20, 30], 1.)?,
        Simplex::new(vec![10, 20, 30], 2.)?,
    ])?;
    let family = OperatorFamily::from_stages(
        &[source.stage(1.)?, source.stage(2.)?],
        1,
        vec![vec![2, 3, 5]; 2],
        "triangle-v1",
        "abstract cost",
        &Execution::default(),
    )?;
    let ring = family.stage(0)?;
    let query = ring.class_query(&[1, 1, 1])?;
    assert_eq!(query.selected_mass(), 10);
    assert_eq!(ring.betti(), 1);
    let tracked = family.track_query(&query, 0, 1)?;
    assert_eq!(tracked.selected_mass(), 0);
    let residual: Vec<_> = query
        .representative()
        .iter()
        .zip(tracked.representative())
        .map(|(a, b)| a ^ b)
        .collect();
    let target = family.stage(1)?;
    let filling = target
        .boundary_preimage(&residual)?
        .expect("class dies through the filling face");
    assert_eq!(target.window().boundaries()[1].apply(&filling)?, residual);
    assert_eq!(family.transport_rank(0, 1)?, 0);
    println!(
        "selected ring mass: {}; finite stage intervals: {:?}",
        query.selected_mass(),
        family.barcode()?
    );
    Ok(())
}
