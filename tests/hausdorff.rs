//! Public Euclidean set-comparison contracts, independent of persistence.

use cocycle::execution::Execution;
use cocycle::geometry::{
    DissimilarityView, PointCloudView, hausdorff_distance, hausdorff_distance_matrix,
    hausdorff_distance_matrix_with, hausdorff_distance_with,
};
use cocycle::{Error, Result};
use std::sync::atomic::{AtomicBool, Ordering};

fn points(coordinates: &[f64], dimension: usize) -> PointCloudView<'_> {
    PointCloudView::new(coordinates, coordinates.len() / dimension, dimension).unwrap()
}

#[test]
fn hand_distances_include_both_directions() -> Result<()> {
    assert_eq!(
        hausdorff_distance(points(&[0., 0.], 2), points(&[3., 4.], 2))?,
        5.
    );
    let a = points(&[0., 2.], 1);
    let b = points(&[0.], 1);
    assert_eq!(hausdorff_distance(a, b)?, 2.);
    assert_eq!(hausdorff_distance(b, a)?, 2.);
    assert_eq!(hausdorff_distance(a, points(&[1., 3.], 1))?, 1.);
    assert_eq!(
        hausdorff_distance(points(&[0., 1.], 1), points(&[0., 1., 100.], 1))?,
        99.
    );
    Ok(())
}

#[test]
fn ordering_duplicates_and_signed_zero_preserve_sets_and_buffers() -> Result<()> {
    let coordinates = [-0., 2., 2., 0.];
    let original: Vec<_> = coordinates.iter().map(|x: &f64| x.to_bits()).collect();
    let a = points(&coordinates, 1);
    let b = points(&[2., 0.], 1);
    assert_eq!(hausdorff_distance(a, b)?.to_bits(), 0.0_f64.to_bits());
    assert_eq!(a.coordinates().as_ptr(), coordinates.as_ptr());
    assert_eq!(
        coordinates.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        original
    );
    assert_eq!(hausdorff_distance(a, a)?, 0.);
    Ok(())
}

#[test]
fn finite_lattice_sets_satisfy_metric_properties() -> Result<()> {
    // Exhaust all 31 nonempty subsets of {-2,-1,0,1,2}; all distances here
    // are exactly representable integers, so no floating tolerance is needed.
    let clouds: Vec<Vec<f64>> = (1_u32..32)
        .map(|mask| {
            (0..5)
                .filter(|i| mask & (1 << i) != 0)
                .map(|i| f64::from(i - 2))
                .collect()
        })
        .collect();
    let samples: Vec<_> = clouds.iter().map(|v| points(v, 1)).collect();
    let mut distances = vec![vec![0.; clouds.len()]; clouds.len()];
    for (i, &a) in samples.iter().enumerate() {
        for (j, &b) in samples.iter().enumerate() {
            distances[i][j] = hausdorff_distance(a, b)?;
            assert!(distances[i][j] >= 0.);
            assert_eq!(distances[i][j] == 0., i == j);
        }
    }
    for i in 0..clouds.len() {
        for j in 0..clouds.len() {
            assert_eq!(distances[i][j], distances[j][i]);
            for k in 0..clouds.len() {
                assert!(distances[i][k] <= distances[i][j] + distances[j][k]);
            }
        }
    }
    Ok(())
}

#[test]
fn common_translation_and_scaling_have_declared_units() -> Result<()> {
    let a = points(&[0., 0., 1., 0.], 2);
    let b = points(&[0., 2., 1., 2.], 2);
    let shifted_a = points(&[10., -4., 11., -4.], 2);
    let shifted_b = points(&[10., -2., 11., -2.], 2);
    assert_eq!(
        hausdorff_distance(a, b)?,
        hausdorff_distance(shifted_a, shifted_b)?
    );
    assert_eq!(
        hausdorff_distance(points(&[0., 0., 3., 0.], 2), points(&[0., 6., 3., 6.], 2))?,
        6.
    );
    Ok(())
}

#[test]
fn stable_norm_handles_large_tiny_and_subnormal_distances() -> Result<()> {
    for value in [1e200, 1e-200, f64::from_bits(1), f64::MAX] {
        assert_eq!(
            hausdorff_distance(points(&[0.], 1), points(&[value], 1))?,
            value
        );
    }
    let value = hausdorff_distance(points(&[0., 0.], 2), points(&[3e200, 4e200], 2))?;
    assert!((value / 5e200 - 1.).abs() < 1e-15);
    Ok(())
}

#[test]
fn irrelevant_overflow_is_ignored_in_either_order() -> Result<()> {
    let a = points(&[-f64::MAX, f64::MAX], 1);
    let b = points(&[f64::MAX, -f64::MAX], 1);
    assert_eq!(hausdorff_distance(a, b)?, 0.);
    assert_eq!(hausdorff_distance(b, a)?, 0.);
    assert_eq!(hausdorff_distance(a, a)?, 0.);
    let left = points(&[-f64::MAX, 0., f64::MAX, 0.], 2);
    let right = points(&[f64::MAX, 1., -f64::MAX, 1.], 2);
    assert_eq!(hausdorff_distance(left, right)?, 1.);
    Ok(())
}

#[test]
fn final_overflow_is_an_error_for_scalar_and_matrix() {
    let a = points(&[-f64::MAX], 1);
    let b = points(&[f64::MAX], 1);
    assert!(matches!(
        hausdorff_distance(a, b),
        Err(Error::NumericalFailure { .. })
    ));
    let samples = [("a", a), ("same", a), ("b", b)];
    assert!(matches!(
        hausdorff_distance_matrix(&samples),
        Err(Error::NumericalFailure { .. })
    ));
    assert!(matches!(
        hausdorff_distance(points(&[0., 0.], 2), points(&[f64::MAX, f64::MAX], 2)),
        Err(Error::NumericalFailure { .. })
    ));
}

#[test]
fn invalid_operands_are_checked_in_both_positions() {
    let a = points(&[0.], 1);
    let empty = points(&[], 1);
    let different_dimension = points(&[0., 0.], 2);
    for (left, right) in [
        (a, empty),
        (empty, a),
        (empty, empty),
        (a, different_dimension),
        (different_dimension, a),
    ] {
        assert!(matches!(
            hausdorff_distance(left, right),
            Err(Error::InvalidParameter { .. })
        ));
    }
}

#[test]
fn matrix_preserves_axes_condensed_order_and_zero_distinct_samples() -> Result<()> {
    let a = points(&[0.], 1);
    let b = points(&[1.], 1);
    let c = points(&[3.], 1);
    let matrix = hausdorff_distance_matrix(&[("c", c), ("a", a), ("b", b), ("a-copy", a)])?;
    assert_eq!(matrix.sample_ids(), ["c", "a", "b", "a-copy"]);
    assert_eq!(matrix.values(), [3., 2., 1., 3., 0., 1.]);
    assert_eq!(matrix.len(), 4);
    assert!(!matrix.is_empty());
    let views = [c, a, b, a];
    let borrowed = DissimilarityView::new(matrix.values(), matrix.len())?;
    assert_eq!(borrowed.values().as_ptr(), matrix.values().as_ptr());
    for i in 0..4 {
        for j in 0..4 {
            assert_eq!(
                matrix.get(i, j),
                Some(hausdorff_distance(views[i], views[j])?)
            );
            assert_eq!(matrix.get(i, j), borrowed.get(i, j));
        }
        assert_eq!(matrix.get(i, i).unwrap().to_bits(), 0.0_f64.to_bits());
    }
    assert_eq!(matrix.get(4, 0), None);
    assert_eq!(matrix.get(0, usize::MAX), None);
    Ok(())
}

#[test]
fn result_owns_ids_and_outlives_coordinates() -> Result<()> {
    let matrix = {
        let coordinates = [0., 2.];
        let id = String::from("owned ID");
        hausdorff_distance_matrix(&[(&id, points(&coordinates, 1))])?
    };
    assert_eq!(matrix.sample_ids(), ["owned ID"]);
    assert_eq!(matrix.get(0, 0), Some(0.));
    Ok(())
}

#[test]
fn empty_and_single_sample_matrices_are_distinct() -> Result<()> {
    let matrix = hausdorff_distance_matrix(&[])?;
    assert_eq!(matrix.len(), 0);
    assert!(matrix.is_empty());
    assert!(matrix.values().is_empty());
    assert_eq!(matrix.get(0, 0), None);
    let matrix = hausdorff_distance_matrix(&[("one", points(&[0.], 1))])?;
    assert_eq!(matrix.len(), 1);
    assert!(!matrix.is_empty());
    assert!(matrix.values().is_empty());
    assert_eq!(matrix.get(0, 0), Some(0.));
    Ok(())
}

#[test]
fn matrix_validates_even_when_no_pair_would_be_evaluated() {
    let a = points(&[0.], 1);
    for samples in [
        vec![("", a)],
        vec![("empty", points(&[], 1))],
        vec![("a", a), ("a", a)],
        vec![("a", a), ("b", points(&[0., 0.], 2))],
    ] {
        assert!(matches!(
            hausdorff_distance_matrix(&samples),
            Err(Error::InvalidParameter { .. })
        ));
    }
    // IDs are exact UTF-8 labels, not normalized identifiers.
    assert!(hausdorff_distance_matrix(&[("é", a), ("e\u{301}", a), (" ", a)]).is_ok());
}

#[test]
fn controlled_calls_match_default_and_whole_matrix_budget_is_shared() -> Result<()> {
    let a = points(&[0.], 1);
    let b = points(&[2.], 1);
    let samples = [("a", a), ("b", b), ("c", a)];
    let execution = Execution::default().max_work(1_000);
    assert_eq!(
        hausdorff_distance_with(a, b, &execution)?,
        hausdorff_distance(a, b)?
    );
    assert_eq!(
        hausdorff_distance_matrix_with(&samples, &execution)?,
        hausdorff_distance_matrix(&samples)?
    );
    // This limit permits any one scalar call, but not validation and three pairs.
    let small = Execution::default().max_work(16);
    assert_eq!(hausdorff_distance_with(a, b, &small)?, 2.);
    assert!(matches!(
        hausdorff_distance_matrix_with(&samples, &small),
        Err(Error::WorkLimitExceeded { limit: 16 })
    ));
    assert_eq!(
        hausdorff_distance_matrix_with(&samples, &execution)?.get(1, 0),
        Some(2.)
    );
    Ok(())
}

#[test]
fn cancellation_precedes_validation_and_zero_budget_rejects_empty_matrix() {
    let flag = AtomicBool::new(true);
    let execution = Execution::default().cancellation(&flag);
    let empty = points(&[], 1);
    assert_eq!(
        hausdorff_distance_with(empty, empty, &execution),
        Err(Error::Cancelled)
    );
    assert_eq!(
        hausdorff_distance_matrix_with(&[("", empty)], &execution),
        Err(Error::Cancelled)
    );
    assert_eq!(
        hausdorff_distance_matrix_with(&[], &execution),
        Err(Error::Cancelled)
    );
    assert!(flag.load(Ordering::Relaxed));
    let zero = Execution::default().max_work(0);
    assert_eq!(
        hausdorff_distance_with(empty, empty, &zero),
        Err(Error::WorkLimitExceeded { limit: 0 })
    );
    assert_eq!(
        hausdorff_distance_matrix_with(&[], &zero),
        Err(Error::WorkLimitExceeded { limit: 0 })
    );
}
