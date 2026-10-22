//! A bounded graph-sheaf investigation using the existing real Hodge terminal.
//! Stalk coordinates below are algebraic coordinates, not geometric simplex IDs.
use cocycle::complex::FilteredComplex;
use cocycle::execution::Execution;
use cocycle::hodge::{HodgeLaplacian, HodgeSelection, SpectrumOptions};
use cocycle::{Error, Result};

// One oriented edge, standard stalk inner products, integer restriction entries.
// The dual chain boundary is delta^T; this is not the graph's ordinary boundary.
struct EdgeSheaf {
    vertices: [usize; 2],
    edge: usize,
    restrictions: [Vec<i32>; 2],
}

impl EdgeSheaf {
    fn new(vertices: [usize; 2], edge: usize, restrictions: [Vec<i32>; 2]) -> Result<Self> {
        let invalid = || Error::InvalidParameter {
            parameter: "sheaf probe",
            reason: "restriction shape, coordinate size or signed coefficient is invalid",
        };
        vertices[0]
            .checked_add(vertices[1])
            .and_then(|n| n.checked_add(edge))
            .ok_or_else(invalid)?;
        for (matrix, width) in restrictions.iter().zip(vertices) {
            if edge.checked_mul(width) != Some(matrix.len()) {
                return Err(invalid());
            }
        }
        if restrictions[0].contains(&i32::MIN) {
            return Err(invalid());
        }
        Ok(Self {
            vertices,
            edge,
            restrictions,
        })
    }
}

impl FilteredComplex for EdgeSheaf {
    type CellId = usize;
    fn cells(&self) -> impl Iterator<Item = usize> + '_ {
        0..self.vertices.iter().sum::<usize>() + self.edge
    }
    fn dimension(&self, cell: usize) -> usize {
        usize::from(cell >= self.vertices.iter().sum())
    }
    fn filtration_value(&self, _: usize) -> f64 {
        0.
    }
    fn boundary(&self, cell: usize) -> impl Iterator<Item = (usize, i32)> + '_ {
        let n = self.vertices.iter().sum();
        (0..n).filter_map(move |i| {
            if cell < n {
                return None;
            }
            let row = cell - n;
            let value = if i < self.vertices[0] {
                -self.restrictions[0][row * self.vertices[0] + i]
            } else {
                self.restrictions[1][row * self.vertices[1] + i - self.vertices[0]]
            };
            (value != 0).then_some((i, value))
        })
    }
}

fn main() -> Result<()> {
    for (name, vertices, edge, left, right) in [
        ("nonconstant", [1, 1], 1, vec![1], vec![2]),
        ("constant", [1, 1], 1, vec![1], vec![1]),
        ("rectangular", [2, 1], 1, vec![1, 0], vec![2]),
        ("zero_edge", [1, 1], 0, vec![], vec![]),
        ("zero_vertices", [0, 0], 1, vec![], vec![]),
        ("empty", [0, 0], 0, vec![], vec![]),
    ] {
        let source = EdgeSheaf::new(vertices, edge, [left, right])?;
        let execution = Execution::default().max_work(100_000);
        let operator = HodgeLaplacian::new_with(&source, 0, HodgeSelection::Full, &execution)?;
        let n = operator.basis().len();
        let spectrum = operator.spectrum(0..n, &SpectrumOptions::default(), &execution)?;
        let x: Vec<_> = (1..=n).map(|i| i as f64).collect();
        println!(
            "{{\"kind\":\"sheaf\",\"name\":{name:?},\"vertices\":{vertices:?},\"edge\":{edge},\"restrictions\":{:?},\"matrix\":{:?},\"values\":{:?},\"x\":{x:?},\"action\":{:?},\"residuals\":{:?}}}",
            source.restrictions,
            spectrum.matrix(),
            spectrum.values(),
            operator.action(&x, &execution)?,
            spectrum.residuals(),
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use cocycle::complex::{Simplex, SimplicialComplex};
    use std::sync::atomic::AtomicBool;

    #[test]
    fn nonconstant_sections_and_constant_specialization_have_independent_expectations() -> Result<()>
    {
        let source = EdgeSheaf::new([1, 1], 1, [vec![1], vec![2]])?;
        let op = HodgeLaplacian::new(&source, 0, HodgeSelection::Full)?;
        drop(source);
        assert_eq!(op.matrix(&Execution::default())?, [1., -2., -2., 4.]);
        assert_eq!(op.action(&[2., 1.], &Execution::default())?, [0., 0.]);
        assert_ne!(op.action(&[1., 1.], &Execution::default())?, [0., 0.]);
        let ordinary = SimplicialComplex::new(vec![
            Simplex::new(vec![0], 0.)?,
            Simplex::new(vec![1], 0.)?,
            Simplex::new(vec![0, 1], 0.)?,
        ])?;
        let constant = EdgeSheaf::new([1, 1], 1, [vec![1], vec![1]])?;
        assert_eq!(
            HodgeLaplacian::new(&constant, 0, HodgeSelection::Full)?
                .matrix(&Execution::default())?,
            HodgeLaplacian::new(&ordinary, 0, HodgeSelection::Full)?
                .matrix(&Execution::default())?,
        );
        Ok(())
    }

    #[test]
    fn rectangular_and_zero_stalks_preserve_coordinate_shapes() -> Result<()> {
        for (vertices, edge, maps, expected) in [
            (
                [2, 1],
                1,
                [vec![1, 0], vec![2]],
                vec![1., 0., -2., 0., 0., 0., -2., 0., 4.],
            ),
            ([1, 1], 0, [vec![], vec![]], vec![0.; 4]),
            ([0, 0], 1, [vec![], vec![]], vec![]),
            ([0, 0], 0, [vec![], vec![]], vec![]),
        ] {
            let source = EdgeSheaf::new(vertices, edge, maps)?;
            let op = HodgeLaplacian::new(&source, 0, HodgeSelection::Full)?;
            assert_eq!(op.matrix(&Execution::default())?, expected);
            assert_eq!(op.next_basis().len(), edge);
        }
        Ok(())
    }

    #[test]
    fn malformed_maps_and_execution_fail_without_a_successful_result() -> Result<()> {
        assert!(EdgeSheaf::new([1, 1], 1, [vec![], vec![2]]).is_err());
        assert!(EdgeSheaf::new([1, 1], 1, [vec![i32::MIN], vec![2]]).is_err());
        assert!(EdgeSheaf::new([usize::MAX, 1], 0, [vec![], vec![]]).is_err());
        let source = EdgeSheaf::new([1, 1], 1, [vec![1], vec![2]])?;
        assert!(matches!(
            HodgeLaplacian::new_with(
                &source,
                0,
                HodgeSelection::Full,
                &Execution::default().max_work(0)
            ),
            Err(Error::WorkLimitExceeded { .. })
        ));
        let cancelled = AtomicBool::new(true);
        assert!(matches!(
            HodgeLaplacian::new_with(
                &source,
                0,
                HodgeSelection::Full,
                &Execution::default().cancellation(&cancelled)
            ),
            Err(Error::Cancelled)
        ));
        HodgeLaplacian::new(&source, 0, HodgeSelection::Full)?;
        Ok(())
    }
}
