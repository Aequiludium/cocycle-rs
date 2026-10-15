//! Standard real Dirac action on a supplied filled triangle, without computing PH.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::dirac::DiracOperator;

fn main() -> cocycle::Result<()> {
    let cells = [
        vec![0],
        vec![1],
        vec![2],
        vec![0, 1],
        vec![0, 2],
        vec![1, 2],
        vec![0, 1, 2],
    ]
    .into_iter()
    .map(|v| Simplex::new(v, 0.))
    .collect::<cocycle::Result<Vec<_>>>()?;
    let source = SimplicialComplex::new(cells)?;
    let dirac = DiracOperator::new(&source, 0..=2)?;
    println!("degree offsets: {:?}", dirac.degree_offsets());
    for (coordinate, &id) in dirac.basis().iter().enumerate() {
        println!(
            "coordinate {coordinate}: {:?}",
            source.simplex(id).unwrap().vertices()
        );
    }
    let x = [1., 2., 3., 4., 5., 6., 7.];
    println!("K x: {:?}", dirac.action(&x)?);
    println!("K squared x: {:?}", dirac.action(&dirac.action(&x)?)?);
    println!("row-major K: {:?}", dirac.matrix(49)?);
    // Truncating to C1 removes both adjacent maps, so K is zero on three edges.
    let edges = DiracOperator::new(&source, 1..=1)?;
    assert_eq!(edges.action(&[1., 2., 3.])?, [0.; 3]);
    Ok(())
}
