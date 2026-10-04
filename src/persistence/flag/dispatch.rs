//! Default selection for exact flag inputs; algorithms do not call this module.
use super::{cohomology, h0};
use crate::algebra::PrimeField;
use crate::complex::WeightedGraph;
use crate::execution::WorkBudget;
use crate::filtration::flag::{CliqueAccess, DenseFlag, SparseFlag};
use crate::geometry::DissimilarityMatrixView;
use crate::persistence::{RawIntervals, simplicial};
use crate::{Error, Result};

pub(in crate::persistence) fn compute_dense(
    input: DissimilarityMatrixView<'_>,
    dimension: usize,
    cutoff: f64,
    field: PrimeField,
    budget: &mut WorkBudget<'_>,
) -> Result<RawIntervals> {
    budget.check()?;
    if dimension == 0 {
        h0::compute(
            input.len(),
            (0..input.len()).flat_map(|b| (0..b).map(move |a| (a, b, input.get(a, b).unwrap()))),
            cutoff,
            budget,
        )
    } else {
        let stop = cutoff.min(crate::filtration::rips::cone_radius(input, &mut || {
            budget.step()
        })?);
        if field.characteristic() != 2 {
            return simplicial::cohomology::compute(
                &CliqueAccess::Dense(input, stop),
                dimension,
                field,
                budget,
            );
        }
        let access = match DenseFlag::new(input, stop) {
            Ok(access) => access,
            Err(Error::SizeOverflow { .. }) if dimension > 1 => {
                return simplicial::cohomology::compute(
                    &CliqueAccess::Dense(input, stop),
                    dimension,
                    field,
                    budget,
                );
            }
            Err(error) => return Err(error),
        };
        budget.check()?;
        if dimension == 1 {
            cohomology::compute(&access, budget)
        } else {
            let (raw, cleared) = cohomology::compute_with_clearing(&access, budget)?;
            simplicial::cohomology::continue_from_h1(
                &CliqueAccess::Dense(input, stop),
                dimension,
                raw,
                cleared,
                budget,
            )
        }
    }
}

pub(in crate::persistence) fn compute_graph(
    graph: &WeightedGraph,
    dimension: usize,
    cutoff: f64,
    field: PrimeField,
    budget: &mut WorkBudget<'_>,
) -> Result<RawIntervals> {
    budget.check()?;
    if dimension == 0 {
        h0::compute(
            graph.vertex_count(),
            graph
                .edges()
                .iter()
                .map(|e| (e.vertices[0], e.vertices[1], e.value)),
            cutoff,
            budget,
        )
    } else {
        if field.characteristic() != 2 {
            return simplicial::cohomology::compute(
                &CliqueAccess::Sparse(graph, cutoff),
                dimension,
                field,
                budget,
            );
        }
        let access = match SparseFlag::new(graph, cutoff) {
            Ok(access) => access,
            Err(Error::SizeOverflow { .. }) if dimension > 1 => {
                return simplicial::cohomology::compute(
                    &CliqueAccess::Sparse(graph, cutoff),
                    dimension,
                    field,
                    budget,
                );
            }
            Err(error) => return Err(error),
        };
        budget.check()?;
        if dimension == 1 {
            cohomology::compute(&access, budget)
        } else {
            let (raw, cleared) = cohomology::compute_with_clearing(&access, budget)?;
            simplicial::cohomology::continue_from_h1(
                &CliqueAccess::Sparse(graph, cutoff),
                dimension,
                raw,
                cleared,
                budget,
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::diagram::Coverage;
    use crate::geometry::DissimilarityView;
    use crate::persistence::assemble_diagram;

    #[test]
    fn sparse_continuation_and_odd_prime_fallback_keep_their_routes() {
        let graph = WeightedGraph::new(
            3,
            vec![
                crate::complex::WeightedEdge {
                    vertices: [0, 1],
                    value: 1.,
                },
                crate::complex::WeightedEdge {
                    vertices: [0, 2],
                    value: 1.,
                },
                crate::complex::WeightedEdge {
                    vertices: [1, 2],
                    value: 1.,
                },
            ],
        )
        .unwrap();
        let execution = crate::execution::Execution::default();
        for characteristic in [2, 3] {
            simplicial::cohomology::take_counts();
            let field = PrimeField::new(characteristic).unwrap();
            let raw = compute_graph(
                &graph,
                2,
                1.,
                field,
                &mut WorkBudget::new(&execution).unwrap(),
            )
            .unwrap();
            let counts = simplicial::cohomology::take_counts();
            if characteristic == 2 {
                assert_eq!(counts[0], [0; 3]);
                assert_eq!(counts[1], [0; 3]);
            } else {
                assert_eq!(counts[0][1], 3);
                assert_eq!(counts[1], [3, 1, 2]);
            }
            assert_eq!(counts[2], [1, 0, 1]);
            let reference = simplicial::cohomology::compute(
                &CliqueAccess::Sparse(&graph, 1.),
                2,
                field,
                &mut WorkBudget::new(&execution).unwrap(),
            )
            .unwrap();
            assert_eq!(
                assemble_diagram(2, Coverage::Complete, raw).unwrap(),
                assemble_diagram(2, Coverage::Complete, reference).unwrap()
            );
            simplicial::cohomology::take_counts();
        }
    }

    #[test]
    fn continuation_clears_virtual_deaths_without_generic_h0_h1_work() {
        // K3's only triangle dies in H1 at zero lifetime; H2 must be empty.
        let input = DissimilarityView::new(&[1.; 3], 3).unwrap();
        let execution = crate::execution::Execution::default();
        let field = PrimeField::new(2).unwrap();
        simplicial::cohomology::take_counts();
        let raw = compute_dense(
            input.into(),
            2,
            1.,
            field,
            &mut WorkBudget::new(&execution).unwrap(),
        )
        .unwrap();
        let counts = simplicial::cohomology::take_counts();
        assert_eq!(counts[0], [0; 3]);
        assert_eq!(counts[1], [0; 3]);
        assert_eq!(counts[2], [1, 0, 1]);
        let diagram = assemble_diagram(2, Coverage::Complete, raw.clone()).unwrap();
        assert!(diagram.dimension(2).unwrap().is_empty());

        let reference = simplicial::cohomology::compute(
            &CliqueAccess::Dense(input.into(), 1.),
            2,
            field,
            &mut WorkBudget::new(&execution).unwrap(),
        )
        .unwrap();
        let counts = simplicial::cohomology::take_counts();
        assert_eq!(counts[0][1], 3);
        assert_eq!(counts[1], [3, 1, 2]);
        assert_eq!(
            assemble_diagram(2, Coverage::Complete, reference).unwrap(),
            diagram
        );

        // The documented empty-clearing proposal would create a false H2 birth.
        let wrong = simplicial::cohomology::continue_from_h1(
            &CliqueAccess::Dense(input.into(), 1.),
            2,
            raw,
            Vec::new(),
            &mut WorkBudget::new(&execution).unwrap(),
        )
        .unwrap();
        assert!(wrong.contains(&(2, 1., None)));
        simplicial::cohomology::take_counts();
    }
}
