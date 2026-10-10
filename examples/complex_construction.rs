//! Construct one complex for persistence and direct signed-boundary reads.

use cocycle::complex::{FilteredComplex, Simplex, SimplexId, SimplicialComplex};
use cocycle::persistence::PersistenceExt;
use cocycle::{Error, Result};

/// Assign each simplex the maximum of its vertex values.
///
/// Vertex IDs are indices into `vertex_values`; all vertices, including isolates,
/// are inserted. Supply each nonvertex simplex exactly once, including every
/// nonvertex face, with strictly increasing IDs. No additional faces are inferred.
/// This example does not expose execution controls or a new library constructor.
fn lower_star_complex(
    vertex_values: &[f64],
    nonvertex_simplices: &[&[usize]],
) -> Result<SimplicialComplex> {
    let capacity = vertex_values
        .len()
        .checked_add(nonvertex_simplices.len())
        .ok_or(Error::SizeOverflow {
            operation: "lower-star simplex count",
        })?;
    let mut simplices = Vec::new();
    simplices
        .try_reserve_exact(capacity)
        .map_err(|_| Error::AllocationFailed {
            context: "lower-star simplices",
        })?;
    // Simplex::new also rejects nonfinite values, including isolated vertices.
    for (vertex, &value) in vertex_values.iter().enumerate() {
        simplices.push(Simplex::new(vec![vertex], value)?);
    }
    for &vertices in nonvertex_simplices {
        if vertices.len() < 2 {
            return Err(Error::InvalidParameter {
                parameter: "nonvertex_simplices",
                reason: "supply only edges and higher-dimensional simplices",
            });
        }
        let mut value = f64::NEG_INFINITY;
        for &vertex in vertices {
            let &birth = vertex_values.get(vertex).ok_or(Error::InvalidParameter {
                parameter: "vertex ID",
                reason: "no value supplied for this vertex",
            })?;
            value = value.max(birth);
        }
        simplices.push(Simplex::new(vertices.to_vec(), value)?);
    }
    SimplicialComplex::new(simplices)
}

/// Borrow one adjacent-degree boundary, exporting IDs with a small dense matrix.
///
/// This tutorial scans stored cells and allocates only the two bases and matrix.
/// Matrix rows/columns follow source traversal restricted to their degree; their
/// counts retain zero-row and zero-column shapes. No PH computation, topology
/// expansion or filtration-value interpretation is needed. Large consumers
/// should read sparse boundary terms instead of exporting a dense matrix.
fn signed_boundary(
    source: &SimplicialComplex,
    degree: usize,
) -> (Vec<SimplexId>, Vec<SimplexId>, Vec<Vec<i32>>) {
    let rows: Vec<_> = source
        .cells()
        .filter(|&id| Some(source.simplex(id).unwrap().dimension()) == degree.checked_sub(1))
        .collect();
    let columns: Vec<_> = source
        .cells()
        .filter(|&id| source.simplex(id).unwrap().dimension() == degree)
        .collect();
    let mut matrix = vec![vec![0; columns.len()]; rows.len()];
    for (column, &id) in columns.iter().enumerate() {
        // IDs come from this frozen owner; face closure proves each row exists.
        for term in source.boundary(id).unwrap() {
            let row = rows.iter().position(|&id| id == term.face).unwrap();
            matrix[row][column] = i32::from(term.coefficient);
        }
    }
    (rows, columns, matrix)
}

fn main() -> Result<()> {
    // A circle with two local minima. These four edges are the entire topology.
    let complex = lower_star_complex(
        &[-2.0, 1.0, -1.0, 0.0],
        &[&[0, 1], &[1, 2], &[2, 3], &[0, 3]],
    )?;
    for simplex in complex.simplices() {
        println!("{:?} enters at {}", simplex.vertices(), simplex.value());
    }
    let edge = complex.find(&[0, 3]).expect("the supplied edge exists");
    for term in complex.boundary(edge).expect("valid simplex ID") {
        let face = complex.simplex(term.face).expect("a boundary face exists");
        println!(
            "Boundary term: {} * {:?}",
            term.coefficient,
            face.vertices()
        );
    }
    // Reader 1 uses all stored topology at one scale, with standard real inner
    // products. This circle has no faces, so its edge Hodge matrix is B1^T B1.
    // Stored lower-star values are not used as conductances or Gram entries.
    let (vertices, edges, boundary) = signed_boundary(&complex, 1);
    let edge_laplacian: Vec<Vec<i32>> = (0..edges.len())
        .map(|a| {
            (0..edges.len())
                .map(|b| boundary.iter().map(|row| row[a] * row[b]).sum())
                .collect()
        })
        .collect();
    println!(
        "Direct boundary shape: {} by {}",
        vertices.len(),
        edges.len()
    );
    println!("Edge basis (owner-local IDs): {edges:?}");
    println!("Single-scale edge Hodge matrix: {edge_laplacian:?}");

    // Reader 2 borrows the very same complex, IDs, orientation and entry values.
    // The engine reads actual vertex births; no Rips-specific assumptions apply.
    let result = complex.persistence().max_homology_dimension(1).compute()?;
    println!("Intervals:");
    for interval in result.diagram().intervals() {
        println!("{interval:?}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use cocycle::algebra::PrimeField;
    use cocycle::diagram::{Coverage, IntervalEnd};

    #[test]
    fn circle_has_hand_derived_pairs_over_several_fields() -> Result<()> {
        let complex = lower_star_complex(
            &[-2.0, 1.0, -1.0, 0.0],
            &[&[0, 1], &[1, 2], &[2, 3], &[0, 3]],
        )?;
        assert_eq!(complex.len(), 8);
        assert_eq!(complex.vertex_count(), 4);
        // Check the construction independently of any persistence diagram.
        let expected: &[(&[usize], f64)] = &[
            (&[0], -2.0),
            (&[1], 1.0),
            (&[2], -1.0),
            (&[3], 0.0),
            (&[0, 1], 1.0),
            (&[1, 2], 1.0),
            (&[2, 3], 0.0),
            (&[0, 3], 0.0),
        ];
        for &(vertices, value) in expected {
            let id = complex.find(vertices).expect("expected simplex exists");
            assert_eq!(complex.simplex(id).unwrap().value(), value);
        }
        for prime in [2, 3, 251] {
            let result = complex
                .persistence()
                .field(PrimeField::new(prime)?)
                .compute()?;
            let pairs: Vec<_> = result
                .diagram()
                .intervals()
                .map(|i| (i.dimension(), i.birth(), i.end()))
                .collect();
            // Minima -2 and -1 merge at 0. At 1 the final two edges close a loop.
            assert_eq!(
                pairs,
                [
                    (0, -2.0, IntervalEnd::Essential),
                    (0, -1.0, IntervalEnd::Finite(0.0)),
                    (1, 1.0, IntervalEnd::Essential),
                ]
            );
        }
        let result = complex.persistence().max_filtration_value(-0.5).compute()?;
        assert_eq!(result.diagram().coverage(), Coverage::Through(-0.5));
        assert_eq!(result.diagram().len(), 2);
        for interval in result.diagram().intervals() {
            assert_eq!(interval.dimension(), 0);
            assert_eq!(interval.end(), IntervalEnd::RightCensored { through: -0.5 });
        }
        Ok(())
    }

    #[test]
    fn tied_triangle_has_expected_oriented_boundary_and_no_persistent_loop() -> Result<()> {
        // Deliberately supply the coface first. The constructor sorts by value
        // and dimension, so all faces still precede it when every value is tied.
        let complex = lower_star_complex(&[2.0; 3], &[&[0, 1, 2], &[1, 2], &[0, 2], &[0, 1]])?;
        assert_eq!(complex.len(), 7);
        let triangle = complex.find(&[0, 1, 2]).unwrap();
        let terms: Vec<_> = complex
            .boundary(triangle)
            .unwrap()
            .iter()
            .map(|t| {
                (
                    complex.simplex(t.face).unwrap().vertices().to_vec(),
                    t.coefficient,
                )
            })
            .collect();
        assert_eq!(terms, [(vec![1, 2], 1), (vec![0, 2], -1), (vec![0, 1], 1)]);
        for edge in [&[0, 1][..], &[0, 2], &[1, 2]] {
            assert_eq!(
                complex.cofacets(complex.find(edge).unwrap()).unwrap(),
                &[triangle]
            );
        }
        // Sum boundary-of-boundary over the integers, using stored incidences.
        let mut coefficients = [0_i32; 3];
        for edge in complex.boundary(triangle).unwrap() {
            for vertex in complex.boundary(edge.face).unwrap() {
                let label = complex.simplex(vertex.face).unwrap().vertices()[0];
                coefficients[label] += i32::from(edge.coefficient) * i32::from(vertex.coefficient);
            }
        }
        assert_eq!(coefficients, [0; 3]);
        let result = complex.persistence().field(PrimeField::new(3)?).compute()?;
        assert_eq!(result.diagram().len(), 1);
        let interval = &result.diagram().interval(0).unwrap();
        assert_eq!(
            (interval.dimension(), interval.birth(), interval.end()),
            (0, 2.0, IntervalEnd::Essential)
        );
        Ok(())
    }

    #[test]
    fn empty_inputs_and_isolates_keep_their_meaning() -> Result<()> {
        assert!(lower_star_complex(&[], &[])?.is_empty());
        let complex = lower_star_complex(&[-3.0, 4.0], &[])?;
        assert_eq!(complex.dimension(), Some(0));
        assert_eq!(complex.vertex_count(), 2);
        let result = complex.persistence().compute()?;
        let births: Vec<_> = result.diagram().intervals().map(|i| i.birth()).collect();
        assert_eq!(births, [-3.0, 4.0]);
        assert!(
            result
                .diagram()
                .intervals()
                .all(|i| i.end() == IntervalEnd::Essential)
        );
        Ok(())
    }

    #[test]
    fn malformed_topology_and_nonfinite_values_are_rejected() {
        for value in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
            assert!(matches!(
                lower_star_complex(&[value], &[]),
                Err(Error::NonFiniteValue { .. })
            ));
        }
        for cells in [&[&[][..]][..], &[&[0]], &[&[0, 2]]] {
            assert!(matches!(
                lower_star_complex(&[0.0; 2], cells),
                Err(Error::InvalidParameter { .. })
            ));
        }
        for cells in [
            &[&[1, 0][..]][..],
            &[&[0, 0]],
            &[&[0, 1], &[0, 1]],
            &[&[0, 1, 2], &[0, 1], &[1, 2]], // Missing edge [0, 2].
        ] {
            assert!(matches!(
                lower_star_complex(&[0.0; 3], cells),
                Err(Error::InvalidComplex { .. })
            ));
        }
    }

    #[test]
    fn direct_reader_retains_labels_orientation_and_filtration_identity() -> Result<()> {
        let complex = SimplicialComplex::new(vec![
            Simplex::new(vec![10, 20, 30], 2.)?,
            Simplex::new(vec![10, 20], 1.)?,
            Simplex::new(vec![10, 30], 1.)?,
            Simplex::new(vec![20, 30], 1.)?,
            Simplex::new(vec![10], 0.)?,
            Simplex::new(vec![20], 0.)?,
            Simplex::new(vec![30], 0.)?,
        ])?;
        let ids_before: Vec<_> = complex.cells().collect();
        let storage = complex.simplices().as_ptr();
        let (vertices, edges, a) = signed_boundary(&complex, 1);
        let (face_rows, faces, d) = signed_boundary(&complex, 2);
        assert_eq!(edges, face_rows);
        let labels = |basis: &[SimplexId]| -> Vec<Vec<usize>> {
            basis
                .iter()
                .map(|&id| complex.simplex(id).unwrap().vertices().to_vec())
                .collect()
        };
        assert_eq!(labels(&vertices), [vec![30], vec![20], vec![10]]);
        assert_eq!(labels(&edges), [vec![20, 30], vec![10, 30], vec![10, 20]]);
        assert_eq!(labels(&faces), [vec![10, 20, 30]]);
        assert_eq!(a, [[1, 1, 0], [-1, 0, 1], [0, -1, -1]]);
        assert_eq!(d, [[1], [-1], [1]]);
        for row in &a {
            assert_eq!(row.iter().zip(&d).map(|(x, y)| x * y[0]).sum::<i32>(), 0);
        }
        // A reversed orientation would change these integer signs even over F2.
        assert_eq!(complex.simplex(faces[0]).unwrap().value(), 2.);
        for prime in [2, 3] {
            let result = complex
                .persistence()
                .field(PrimeField::new(prime)?)
                .compute()?;
            let h1 = result.diagram().dimension(1)?.iter().next().unwrap();
            assert_eq!((h1.birth(), h1.end()), (1., IntervalEnd::Finite(2.)));
        }
        assert_eq!(complex.cells().collect::<Vec<_>>(), ids_before);
        assert_eq!(complex.simplices().as_ptr(), storage);
        Ok(())
    }

    #[test]
    fn direct_reader_keeps_empty_shapes_and_invalid_ids_distinct() -> Result<()> {
        let empty = lower_star_complex(&[], &[])?;
        for degree in [0, 1, 2, usize::MAX] {
            let (rows, columns, matrix) = signed_boundary(&empty, degree);
            assert!(rows.is_empty() && columns.is_empty() && matrix.is_empty());
        }
        let isolates = lower_star_complex(&[-3., 4.], &[])?;
        let (rows, columns, matrix) = signed_boundary(&isolates, 0);
        assert!(rows.is_empty() && matrix.is_empty());
        assert_eq!(columns.len(), 2); // Ordinary H0 boundary: 0 by 2.
        let (rows, columns, matrix) = signed_boundary(&isolates, 1);
        assert_eq!(rows.len(), 2);
        assert!(columns.is_empty());
        assert_eq!(matrix, [Vec::<i32>::new(), Vec::new()]); // 2 by 0.
        let other = lower_star_complex(&[0.; 3], &[&[0, 1], &[1, 2]])?;
        let out_of_range = other.cells().last().unwrap();
        assert_eq!(isolates.simplex(out_of_range), None);
        assert_eq!(isolates.boundary(out_of_range), None);
        // An in-range ID from another owner cannot be detected: never reuse it.
        Ok(())
    }
}
