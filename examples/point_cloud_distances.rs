//! Compare points directly and print a distance matrix with sample identities.
use cocycle::execution::Execution;
use cocycle::geometry::{PointCloudView, hausdorff_distance, hausdorff_distance_with};
use cocycle::{Error, Result};

struct Sample<'a> {
    id: &'a str,
    unit: &'a str,
    coordinate_frame: &'a str,
    points: PointCloudView<'a>,
}

// This example retains the original sample order and metadata alongside values.
// It is a local composition, not a universal batch or kernel-matrix API.
struct DistanceMatrix<'a> {
    samples: &'a [Sample<'a>],
    values: Vec<f64>,
}

fn distance_matrix<'a>(
    samples: &'a [Sample<'a>],
    execution: &Execution<'_>,
) -> Result<DistanceMatrix<'a>> {
    for (i, sample) in samples.iter().enumerate() {
        if sample.id.is_empty() || samples[..i].iter().any(|other| other.id == sample.id) {
            return Err(Error::InvalidParameter {
                parameter: "sample ID",
                reason: "sample IDs must be nonempty and unique",
            });
        }
        if sample.unit.is_empty() || sample.coordinate_frame.is_empty() {
            return Err(Error::InvalidParameter {
                parameter: "sample coordinates",
                reason: "unit and coordinate frame must be declared",
            });
        }
        let reference = &samples[0];
        if sample.unit != reference.unit || sample.coordinate_frame != reference.coordinate_frame {
            return Err(Error::InvalidParameter {
                parameter: "sample coordinates",
                reason: "units and coordinate frames must agree",
            });
        }
        if sample.points.dimension() != reference.points.dimension() {
            return Err(Error::ShapeMismatch {
                input: "sample coordinate dimension",
                expected: reference.points.dimension(),
                actual: sample.points.dimension(),
            });
        }
    }
    let count = samples
        .len()
        .checked_mul(samples.len())
        .ok_or(Error::SizeOverflow {
            operation: "sample count squared",
        })?;
    let mut values = Vec::new();
    values
        .try_reserve_exact(count)
        .map_err(|_| Error::AllocationFailed {
            context: "sample distance matrix",
        })?;
    values.resize(count, 0.);
    for (i, first) in samples.iter().enumerate() {
        for (j, second) in samples.iter().enumerate().skip(i) {
            // Each scalar call has its own fresh budget, including diagonal calls.
            // Errors drop the private buffer; a partially filled matrix never escapes.
            let value = hausdorff_distance_with(first.points, second.points, execution)?;
            values[i * samples.len() + j] = value;
            values[j * samples.len() + i] = value;
        }
    }
    Ok(DistanceMatrix { samples, values })
}

fn main() -> Result<()> {
    let a = PointCloudView::new(&[0., 2.], 2, 1)?;
    let b = PointCloudView::new(&[0., 1., 2.], 3, 1)?;
    let c = PointCloudView::new(&[3., 5.], 2, 1)?;
    assert_eq!(hausdorff_distance(a, b)?, 1.);
    let samples = [
        Sample {
            id: "endpoints",
            unit: "m",
            coordinate_frame: "survey-axis",
            points: a,
        },
        Sample {
            id: "with-midpoint",
            unit: "m",
            coordinate_frame: "survey-axis",
            points: b,
        },
        Sample {
            id: "translated",
            unit: "m",
            coordinate_frame: "survey-axis",
            points: c,
        },
    ];
    let matrix = distance_matrix(&samples, &Execution::default())?;
    assert_eq!(matrix.values, [0., 1., 3., 1., 0., 3., 3., 3., 0.]);
    println!(
        "Hausdorff distance; unit={}; frame={}",
        matrix.samples[0].unit, matrix.samples[0].coordinate_frame
    );
    println!(
        "sample,{}",
        matrix
            .samples
            .iter()
            .map(|sample| sample.id)
            .collect::<Vec<_>>()
            .join(",")
    );
    for (sample, row) in matrix
        .samples
        .iter()
        .zip(matrix.values.chunks_exact(samples.len()))
    {
        println!(
            "{},{}",
            sample.id,
            row.iter().map(f64::to_string).collect::<Vec<_>>().join(",")
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::sync::atomic::AtomicBool;

    #[test]
    fn matrix_retains_sample_order_and_metadata() {
        let a = PointCloudView::new(&[0.], 1, 1).unwrap();
        let b = PointCloudView::new(&[3.], 1, 1).unwrap();
        let samples = [
            Sample {
                id: "second",
                unit: "m",
                coordinate_frame: "axis",
                points: b,
            },
            Sample {
                id: "first",
                unit: "m",
                coordinate_frame: "axis",
                points: a,
            },
        ];
        let matrix = distance_matrix(&samples, &Execution::default()).unwrap();
        assert_eq!(matrix.samples[0].id, "second");
        assert_eq!(matrix.samples[1].id, "first");
        assert_eq!(matrix.samples[0].unit, "m");
        assert_eq!(matrix.values, [0., 3., 3., 0.]);
        assert!(
            distance_matrix(&[], &Execution::default())
                .unwrap()
                .values
                .is_empty()
        );
    }

    #[test]
    fn matrix_checks_every_sample_before_computing_and_returns_no_partial_values() {
        let point = PointCloudView::new(&[0.], 1, 1).unwrap();
        let mut samples = [
            Sample {
                id: "a",
                unit: "m",
                coordinate_frame: "axis",
                points: point,
            },
            Sample {
                id: "b",
                unit: "cm",
                coordinate_frame: "axis",
                points: point,
            },
        ];
        assert!(distance_matrix(&samples, &Execution::default()).is_err());
        samples[1].unit = "m";
        samples[1].coordinate_frame = "different-axis";
        assert!(distance_matrix(&samples, &Execution::default()).is_err());
        samples[1].coordinate_frame = "axis";
        samples[1].id = "a";
        assert!(distance_matrix(&samples, &Execution::default()).is_err());
        samples[1].id = "b";
        samples[1].points = PointCloudView::new(&[0., 0.], 1, 2).unwrap();
        assert!(matches!(
            distance_matrix(&samples, &Execution::default()),
            Err(Error::ShapeMismatch { .. })
        ));
        samples[1].points = point;
        assert!(matches!(
            distance_matrix(&samples, &Execution::default().max_work(0)),
            Err(Error::WorkLimitExceeded { .. })
        ));
        let flag = AtomicBool::new(true);
        assert!(matches!(
            distance_matrix(&samples, &Execution::default().cancellation(&flag)),
            Err(Error::Cancelled)
        ));
        // Diagonal call succeeds (three units); the cross pair fails later.
        samples[1].points = PointCloudView::new(&[0., 1.], 2, 1).unwrap();
        assert!(matches!(
            distance_matrix(&samples, &Execution::default().max_work(3)),
            Err(Error::WorkLimitExceeded { .. })
        ));
        assert!(distance_matrix(&samples, &Execution::default()).is_ok());
    }
}
