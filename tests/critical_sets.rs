//! Public source, endpoint and execution contracts for critical sets.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::optimization::CriticalSetWorkspace;
use cocycle::{Error, Result};

#[test]
fn finite_h0_and_essential_circle_use_source_ids_and_keep_source_immutable() -> Result<()> {
    let cells: &[(&[usize], f64)] = &[
        (&[0], -2.),
        (&[1], 1.),
        (&[2], -1.),
        (&[3], 0.),
        (&[0, 1], 1.),
        (&[1, 2], 1.),
        (&[2, 3], 0.),
        (&[0, 3], 0.),
    ];
    let source = SimplicialComplex::new(
        cells
            .iter()
            .map(|&(v, f)| Simplex::new(v.to_vec(), f))
            .collect::<Result<_>>()?,
    )?;
    let original: Vec<_> = source
        .simplices()
        .iter()
        .map(|s| (s.vertices().to_vec(), s.value()))
        .collect();
    let mut workspace = CriticalSetWorkspace::new(&source)?;
    let pairs: Vec<_> = workspace.finite_pairs().collect();
    assert_eq!(pairs.len(), 3); // Zero bars remain available to endpoint callers.
    let (birth, death) = pairs
        .into_iter()
        .find(|&(b, d)| source.simplex(b).unwrap().value() < source.simplex(d).unwrap().value())
        .unwrap();
    assert_eq!(source.simplex(birth).unwrap().vertices(), [2]);
    assert_eq!(source.simplex(death).unwrap().vertices(), [0, 3]);
    assert_eq!(workspace.critical_set(birth, -0.5)?, [birth]);
    let decreased = workspace.critical_set(death, -0.5)?;
    assert_eq!(decreased, [source.find(&[2, 3]).unwrap(), death]);
    let essential: Vec<_> = workspace.essential_births().collect();
    assert_eq!(essential.len(), 2); // H0 minimum and H1 circle, relative to source.
    let circle = essential[1];
    assert_eq!(workspace.critical_set(circle, -3.)?.len(), 4);
    assert_eq!(workspace.critical_set(circle, 2.)?, [circle]);
    assert!(matches!(
        workspace.critical_set_with(birth, -0.5, &Execution::default().max_work(0)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert_eq!(workspace.critical_set(birth, -0.5)?, [birth]);
    assert_eq!(
        source
            .simplices()
            .iter()
            .map(|s| (s.vertices().to_vec(), s.value()))
            .collect::<Vec<_>>(),
        original
    );
    Ok(())
}
