//! Application-owned composition of geometry, computed results and aligned filtrations.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::diagram::{IntervalEnd, PersistenceDiagram, PersistenceResult};
use cocycle::diagram_distances::{
    bottleneck_distance, bottleneck_distance_results_with, wasserstein_1_infinity,
    wasserstein_1_infinity_results_with, wasserstein_2_euclidean,
    wasserstein_2_euclidean_results_with,
};
use cocycle::execution::Execution;
use cocycle::filtration::RipsBuilder;
use cocycle::geometry::{DissimilarityView, PointCloudView, hausdorff_distance_matrix_with};
use cocycle::persistence::PersistenceExt;
use cocycle::{Error, Result};

// Caller declarations stay with both operands. They are not authenticated by
// the crate, and a source snapshot need not be the same across compared samples.
#[derive(Clone, Copy)]
struct Sample<'a> {
    id: &'a str,
    snapshot: &'a str,
    units: Option<&'a str>,
    alignment: Option<&'a str>,
    normalization: Option<&'a str>,
}

fn compatible(first: Sample<'_>, second: Sample<'_>) -> Result<()> {
    for sample in [first, second] {
        if sample.id.is_empty() || sample.snapshot.is_empty() {
            return Err(invalid(
                "sample",
                "sample ID and source snapshot must be declared",
            ));
        }
    }
    for (parameter, a, b) in [
        ("units", first.units, second.units),
        ("alignment", first.alignment, second.alignment),
        ("normalization", first.normalization, second.normalization),
    ] {
        match (a, b) {
            (Some(a), Some(b)) if !a.is_empty() && a == b => {}
            _ => {
                return Err(invalid(
                    parameter,
                    "comparison requires equal known declarations",
                ));
            }
        }
    }
    Ok(())
}

fn invalid(parameter: &'static str, reason: &'static str) -> Error {
    Error::InvalidParameter { parameter, reason }
}

fn result_distances(
    first: (Sample<'_>, &PersistenceResult),
    second: (Sample<'_>, &PersistenceResult),
    degree: usize,
    execution: &Execution<'_>,
) -> Result<[f64; 3]> {
    compatible(first.0, second.0)?;
    // Each existing terminal starts its own budget. This application composition
    // does not claim a cumulative transaction or a pipeline-wide memory limit.
    Ok([
        bottleneck_distance_results_with(first.1, second.1, degree, execution)?,
        wasserstein_1_infinity_results_with(first.1, second.1, degree, execution)?,
        wasserstein_2_euclidean_results_with(first.1, second.1, degree, execution)?,
    ])
}

fn raw_distances(
    first: &PersistenceDiagram,
    second: &PersistenceDiagram,
    degree: usize,
) -> Result<[f64; 3]> {
    Ok([
        bottleneck_distance(first, second, degree)?,
        wasserstein_1_infinity(first, second, degree)?,
        wasserstein_2_euclidean(first, second, degree)?,
    ])
}

// A bounded application recipe, not a public comparison API. Both owners are
// complete assignments on exactly the supplied simplex-key domain. Source-local
// IDs cannot align owners: changing values can change their filtration ordering.
fn aligned_filtration_distance(
    first: (Sample<'_>, &SimplicialComplex),
    second: (Sample<'_>, &SimplicialComplex),
) -> Result<f64> {
    compatible(first.0, second.0)?;
    if first.1.len() != second.1.len() {
        return Err(invalid(
            "cell domain",
            "filtrations must have identical simplex keys",
        ));
    }
    let mut maximum: f64 = 0.;
    for cell in first.1.simplices() {
        let image = second
            .1
            .find(cell.vertices())
            .and_then(|id| second.1.simplex(id))
            .ok_or_else(|| {
                invalid(
                    "cell domain",
                    "filtrations must have identical simplex keys",
                )
            })?;
        maximum = maximum.max((cell.value() - image.value()).abs());
    }
    if !maximum.is_finite() {
        return Err(Error::NumericalFailure {
            context: "aligned filtration sup norm",
        });
    }
    Ok(maximum)
}

fn triangle(edges: [f64; 3], face: f64) -> Result<SimplicialComplex> {
    let cells = [
        (vec![10], 0.),
        (vec![20], 0.),
        (vec![30], 0.),
        (vec![10, 20], edges[0]),
        (vec![10, 30], edges[1]),
        (vec![20, 30], edges[2]),
        (vec![10, 20, 30], face),
    ];
    SimplicialComplex::new(
        cells
            .into_iter()
            .map(|(v, f)| Simplex::new(v, f))
            .collect::<Result<_>>()?,
    )
}

fn diagram_json(diagram: &PersistenceDiagram, degree: usize) -> Result<String> {
    let mut entries = Vec::new();
    for interval in diagram.dimension(degree)?.iter() {
        let death = match interval.end() {
            IntervalEnd::Finite(value) => format!("{value:?}"),
            IntervalEnd::Essential => String::from("null"),
            IntervalEnd::RightCensored { .. } => {
                return Err(invalid("coverage", "fixture must be complete"));
            }
        };
        entries.push(format!("[{:?},{death}]", interval.birth()));
    }
    Ok(format!("[{}]", entries.join(",")))
}

fn cells_json(source: &SimplicialComplex) -> String {
    let cells: Vec<_> = source
        .simplices()
        .iter()
        .map(|cell| format!("[{:?},{:?}]", cell.vertices(), cell.value()))
        .collect();
    format!("[{}]", cells.join(","))
}

fn routes(emit: bool) -> Result<()> {
    let geometry = Sample {
        id: "a",
        snapshot: "measurements/a/v1",
        units: Some("metres"),
        alignment: Some("lab-cartesian-v1"),
        normalization: Some("none"),
    };
    let second = Sample {
        id: "b",
        snapshot: "measurements/b/v1",
        ..geometry
    };
    let third = Sample {
        id: "c",
        snapshot: "measurements/c/v1",
        ..geometry
    };
    let prepare = std::time::Instant::now();
    let coordinates = [[0., 0., 2., 0.], [0., 1., 2., 1.], [0., 3., 2., 3.]];
    let views = coordinates
        .iter()
        .map(|v| PointCloudView::new(v, 2, 2))
        .collect::<Result<Vec<_>>>()?;
    let prepare_ns = prepare.elapsed().as_nanos();
    let compare = std::time::Instant::now();
    compatible(geometry, second)?;
    compatible(geometry, third)?;
    let matrix = hausdorff_distance_matrix_with(
        &[
            (geometry.id, views[0]),
            (second.id, views[1]),
            (third.id, views[2]),
        ],
        &Execution::default(),
    )?;
    let compare_ns = compare.elapsed().as_nanos();
    let analysis = std::time::Instant::now();
    let nearest = (1..matrix.len())
        .min_by(|&i, &j| {
            matrix
                .get(0, i)
                .unwrap()
                .total_cmp(&matrix.get(0, j).unwrap())
        })
        .unwrap();
    let analysis_ns = analysis.elapsed().as_nanos();
    let convert = std::time::Instant::now();
    let borrowed = DissimilarityView::new(matrix.values(), matrix.len())?;
    let alias = borrowed.values().as_ptr() == matrix.values().as_ptr();
    let conversion_ns = convert.elapsed().as_nanos();
    assert_eq!(matrix.values(), [1., 3., 2.]);
    assert_eq!(matrix.sample_ids()[nearest], "b");
    assert!(alias);
    if emit {
        println!(
            "{{\"route\":\"geometry\",\"ids\":{:?},\"snapshots\":{:?},\"units\":\"metres\",\"alignment\":\"lab-cartesian-v1\",\"normalization\":\"none\",\"points\":{coordinates:?},\"dimension\":2,\"values\":{:?},\"nearest\":{:?},\"converted_alias\":{alias},\"numeric_bytes\":{},\"prepare_ns\":{prepare_ns},\"compare_ns\":{compare_ns},\"analysis_ns\":{analysis_ns},\"conversion_ns\":{conversion_ns}}}",
            matrix.sample_ids(),
            [geometry.snapshot, second.snapshot, third.snapshot],
            matrix.values(),
            matrix.sample_ids()[nearest],
            std::mem::size_of_val(matrix.values())
        );
    }

    let prepare = std::time::Instant::now();
    let left_coordinates = [0., 2.];
    let right_coordinates = [0., 3.];
    let left = RipsBuilder::from_points(PointCloudView::new(&left_coordinates, 2, 1)?);
    let right = RipsBuilder::from_points(PointCloudView::new(&right_coordinates, 2, 1)?);
    let prepare_ns = prepare.elapsed().as_nanos();
    let compute = std::time::Instant::now();
    let first = left.persistence().max_homology_dimension(0).compute()?;
    let other = right.persistence().max_homology_dimension(0).compute()?;
    let compute_ns = compute.elapsed().as_nanos();
    let first_meta = Sample {
        id: "rips-a",
        snapshot: "rips/a/v1",
        ..geometry
    };
    let other_meta = Sample {
        id: "rips-b",
        snapshot: "rips/b/v1",
        ..geometry
    };
    let compare = std::time::Instant::now();
    let distances = result_distances(
        (first_meta, &first),
        (other_meta, &other),
        0,
        &Execution::default(),
    )?;
    let compare_ns = compare.elapsed().as_nanos();
    let convert = std::time::Instant::now();
    // Moving the stored diagrams discards context explicitly, without recomputing
    // PH or cloning their buffers. The operand metadata remains in this scope.
    let first = first.into_diagram();
    let other = other.into_diagram();
    let conversion_ns = convert.elapsed().as_nanos();
    let converted_compare = std::time::Instant::now();
    let converted = raw_distances(&first, &other, 0)?;
    let converted_compare_ns = converted_compare.elapsed().as_nanos();
    let analysis = std::time::Instant::now();
    let within_one_metre = distances[0] <= 1.;
    let analysis_ns = analysis.elapsed().as_nanos();
    assert_eq!(distances, converted);
    assert_eq!(distances, [1., 1., 1.]);
    if emit {
        println!(
            "{{\"route\":\"results\",\"ids\":{:?},\"snapshots\":{:?},\"units\":\"metres\",\"alignment\":\"lab-cartesian-v1\",\"normalization\":\"none\",\"field\":2,\"degree\":0,\"first_points\":{left_coordinates:?},\"second_points\":{right_coordinates:?},\"first_diagram\":{},\"second_diagram\":{},\"distances\":{distances:?},\"converted\":{converted:?},\"within_one_metre\":{within_one_metre},\"prepare_ns\":{prepare_ns},\"compute_ns\":{compute_ns},\"compare_ns\":{compare_ns},\"conversion_ns\":{conversion_ns},\"converted_compare_ns\":{converted_compare_ns},\"analysis_ns\":{analysis_ns}}}",
            [first_meta.id, other_meta.id],
            [first_meta.snapshot, other_meta.snapshot],
            diagram_json(&first, 0)?,
            diagram_json(&other, 0)?
        );
    }

    let meta = Sample {
        id: "filtration-a",
        snapshot: "triangle/a/v1",
        units: Some("seconds"),
        alignment: Some("triangle-labels-v1"),
        normalization: Some("none"),
    };
    let other_meta = Sample {
        id: "filtration-b",
        snapshot: "triangle/b/v1",
        ..meta
    };
    let prepare = std::time::Instant::now();
    let first = triangle([1., 2., 3.], 6.)?;
    let other = triangle([3., 1., 2.], 7.)?;
    let prepare_ns = prepare.elapsed().as_nanos();
    let compare = std::time::Instant::now();
    let sup = aligned_filtration_distance((meta, &first), (other_meta, &other))?;
    let compare_ns = compare.elapsed().as_nanos();
    let compute = std::time::Instant::now();
    let a = first.persistence().max_homology_dimension(1).compute()?;
    let b = other.persistence().max_homology_dimension(1).compute()?;
    let compute_ns = compute.elapsed().as_nanos();
    // Supplied filtrations have an unspecified parameter convention. Even two
    // such results must not bypass the context-aware adapter's rejection.
    assert!(matches!(
        bottleneck_distance_results_with(&a, &b, 1, &Execution::default()),
        Err(Error::IncompatibleDiagramContext { .. })
    ));
    // Explicit application-established field, seconds and common cell domain
    // permit the existing raw diagram call, not a forged edge-length context.
    assert_eq!(a.context().characteristic(), b.context().characteristic());
    let result_compare = std::time::Instant::now();
    let result_distances = raw_distances(a.diagram(), b.diagram(), 1)?;
    let result_compare_ns = result_compare.elapsed().as_nanos();
    assert_eq!(sup, 2.);
    assert_eq!(result_distances, [1., 1., 1.]);
    assert!(result_distances[0] <= sup);
    assert_ne!(first.find(&[10, 20]), other.find(&[10, 20]));
    if emit {
        println!(
            "{{\"route\":\"filtrations\",\"ids\":{:?},\"snapshots\":{:?},\"units\":\"seconds\",\"alignment\":\"triangle-labels-v1\",\"normalization\":\"none\",\"field\":2,\"degree\":1,\"first_cells\":{},\"second_cells\":{},\"first_diagram\":{},\"second_diagram\":{},\"sup\":{sup},\"distances\":{result_distances:?},\"bound_holds\":true,\"prepare_ns\":{prepare_ns},\"compute_ns\":{compute_ns},\"compare_ns\":{compare_ns},\"result_compare_ns\":{result_compare_ns}}}",
            [meta.id, other_meta.id],
            [meta.snapshot, other_meta.snapshot],
            cells_json(&first),
            cells_json(&other),
            diagram_json(a.diagram(), 1)?,
            diagram_json(b.diagram(), 1)?
        );
    }
    Ok(())
}

fn main() -> Result<()> {
    match std::env::args().nth(1).as_deref() {
        None => routes(true),
        Some("--measure") => {
            let start = std::time::Instant::now();
            for _ in 0..100 {
                routes(false)?;
            }
            println!(
                "{{\"iterations\":100,\"seconds_per_route_set\":{}}}",
                start.elapsed().as_secs_f64() / 100.
            );
            Ok(())
        }
        _ => Err(invalid("argument", "expected no argument or --measure")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::{AtomicBool, Ordering};

    fn meta() -> Sample<'static> {
        Sample {
            id: "first",
            snapshot: "test/v1",
            units: Some("metres"),
            alignment: Some("labels/v1"),
            normalization: Some("none"),
        }
    }

    #[test]
    fn three_real_routes_match_hand_values() {
        routes(false).unwrap();
    }

    #[test]
    fn unknown_mismatched_declarations_and_missing_identity_fail() {
        let first = meta();
        for second in [
            Sample {
                units: None,
                ..first
            },
            Sample {
                units: Some("centimetres"),
                ..first
            },
            Sample {
                alignment: Some("other"),
                ..first
            },
            Sample {
                alignment: None,
                ..first
            },
            Sample {
                normalization: None,
                ..first
            },
            Sample {
                normalization: Some("unit diameter"),
                ..first
            },
            Sample { id: "", ..first },
            Sample {
                snapshot: "",
                ..first
            },
        ] {
            assert!(compatible(first, second).is_err());
            assert!(compatible(second, first).is_err());
        }
        compatible(
            first,
            Sample {
                id: "other",
                snapshot: "other/v2",
                ..first
            },
        )
        .unwrap();
    }

    #[test]
    fn actual_result_context_coverage_and_degree_errors_propagate() -> Result<()> {
        let points = PointCloudView::new(&[0., 2.], 2, 1)?;
        let source = RipsBuilder::from_points(points);
        let complete = source.persistence().max_homology_dimension(0).compute()?;
        let odd = source
            .persistence()
            .max_homology_dimension(0)
            .field(cocycle::algebra::PrimeField::new(3)?)
            .compute()?;
        assert!(matches!(
            result_distances(
                (meta(), &complete),
                (meta(), &odd),
                0,
                &Execution::default()
            ),
            Err(Error::IncompatibleDiagramContext { .. })
        ));
        let truncated = RipsBuilder::from_points(points)
            .max_edge_length(1.)
            .persistence()
            .max_homology_dimension(0)
            .compute()?;
        assert!(matches!(
            result_distances(
                (meta(), &complete),
                (meta(), &truncated),
                0,
                &Execution::default()
            ),
            Err(Error::IncompleteDiagram { .. })
        ));
        assert!(matches!(
            result_distances(
                (meta(), &complete),
                (meta(), &complete),
                1,
                &Execution::default()
            ),
            Err(Error::DimensionNotComputed { .. })
        ));
        Ok(())
    }

    #[test]
    fn comparison_failure_does_not_recompute_or_discard_owned_results() -> Result<()> {
        let source = RipsBuilder::from_points(PointCloudView::new(&[0., 2.], 2, 1)?);
        let result = source.persistence().max_homology_dimension(0).compute()?;
        let diagram_ptr = result.diagram() as *const _;
        let flag = AtomicBool::new(true);
        for execution in [
            Execution::default().max_work(0),
            Execution::default().cancellation(&flag),
        ] {
            assert!(result_distances((meta(), &result), (meta(), &result), 0, &execution).is_err());
            assert_eq!(diagram_ptr, result.diagram() as *const _);
        }
        assert!(flag.load(Ordering::Relaxed));
        assert_eq!(
            result_distances(
                (meta(), &result),
                (meta(), &result),
                0,
                &Execution::default()
            )?,
            [0.; 3]
        );
        Ok(())
    }

    #[test]
    fn filtration_alignment_uses_keys_and_rejects_equal_shape_different_domains() -> Result<()> {
        let first = triangle([1., 2., 3.], 6.)?;
        let other = triangle([3., 1., 2.], 7.)?;
        assert_ne!(first.find(&[10, 20]), other.find(&[10, 20]));
        assert_eq!(
            aligned_filtration_distance((meta(), &first), (meta(), &other))?,
            2.
        );
        assert_eq!(
            aligned_filtration_distance((meta(), &other), (meta(), &first))?,
            2.
        );
        let cells = first
            .simplices()
            .iter()
            .map(|s| Simplex::new(s.vertices().iter().map(|x| x + 1).collect(), s.value()))
            .collect::<Result<_>>()?;
        let relabeled = SimplicialComplex::new(cells)?;
        assert!(aligned_filtration_distance((meta(), &first), (meta(), &relabeled)).is_err());
        let missing = SimplicialComplex::new(vec![Simplex::new(vec![10], 0.)?])?;
        assert!(aligned_filtration_distance((meta(), &first), (meta(), &missing)).is_err());
        Ok(())
    }

    #[test]
    fn filtration_self_shift_empty_and_unrepresentable_values() -> Result<()> {
        let first = triangle([1., 2., 3.], 6.)?;
        let shifted = SimplicialComplex::new(
            first
                .simplices()
                .iter()
                .rev()
                .map(|s| Simplex::new(s.vertices().to_vec(), s.value() + 2.))
                .collect::<Result<_>>()?,
        )?;
        assert_eq!(
            aligned_filtration_distance((meta(), &first), (meta(), &first))?,
            0.
        );
        assert_eq!(
            aligned_filtration_distance((meta(), &first), (meta(), &shifted))?,
            2.
        );
        let empty = SimplicialComplex::new(vec![])?;
        assert_eq!(
            aligned_filtration_distance((meta(), &empty), (meta(), &empty))?,
            0.
        );
        let low = SimplicialComplex::new(vec![Simplex::new(vec![0], -f64::MAX)?])?;
        let high = SimplicialComplex::new(vec![Simplex::new(vec![0], f64::MAX)?])?;
        assert!(matches!(
            aligned_filtration_distance((meta(), &low), (meta(), &high)),
            Err(Error::NumericalFailure { .. })
        ));
        Ok(())
    }
}
