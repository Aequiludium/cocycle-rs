//! Euclidean distance stability, geometric invariance and independent calls.

use cocycle::Error;
use cocycle::descriptors::betti_curve;
use cocycle::diagram::{Coverage, IntervalEnd};
use cocycle::execution::Execution;
use cocycle::geometry::{DissimilarityView, PointCloudView};
use cocycle::geometry::{hausdorff_distance, hausdorff_distance_with};
use cocycle::persistence::{RipsOptions, rips_from_dissimilarities, rips_from_points};

#[test]
fn square_matches_distances_and_survives_release_of_coordinate_storage() {
    let diagram = {
        let points = vec![0., 0., 1., 0., 1., 1., 0., 1.];
        rips_from_points(
            PointCloudView::new(&points, 4, 2).unwrap(),
            &RipsOptions::default(),
        )
        .unwrap()
    };
    let diagonal = 2.0_f64.sqrt();
    let values = [1., diagonal, 1., 1., diagonal, 1.];
    let from_distances = rips_from_dissimilarities(
        DissimilarityView::new(&values, 4).unwrap(),
        &RipsOptions::default(),
    )
    .unwrap();
    assert_eq!(diagram, from_distances);
    assert_eq!(
        betti_curve(&diagram, 0, &[0., 1., diagonal, 2.]).unwrap(),
        [4, 1, 1, 1]
    );
    assert_eq!(
        betti_curve(&diagram, 1, &[0., 1., diagonal, 2.]).unwrap(),
        [0, 1, 0, 0]
    );
    let complete_with_cap = rips_from_dissimilarities(
        DissimilarityView::new(&values, 4).unwrap(),
        &RipsOptions::new(1, Some(2.)).unwrap(),
    )
    .unwrap();
    assert_eq!(complete_with_cap.coverage(), Coverage::Complete);
    assert_eq!(betti_curve(&complete_with_cap, 1, &[3.]).unwrap(), [0]);
}

#[test]
fn stable_norm_handles_large_and_tiny_representable_distances() {
    for scale in [1e154, 1e-200] {
        let coordinates = [0., 0., 3. * scale, 4. * scale];
        let diagram = rips_from_points(
            PointCloudView::new(&coordinates, 2, 2).unwrap(),
            &RipsOptions::new(0, None).unwrap(),
        )
        .unwrap();
        let IntervalEnd::Finite(death) = diagram.interval(0).unwrap().end() else {
            panic!()
        };
        assert!((death / (5. * scale) - 1.).abs() < 1e-14);
    }
    for coordinates in [vec![-f64::MAX, f64::MAX], vec![0., 0., f64::MAX, f64::MAX]] {
        assert!(matches!(
            rips_from_points(
                PointCloudView::new(&coordinates, 2, coordinates.len() / 2).unwrap(),
                &RipsOptions::default()
            ),
            Err(Error::NumericalFailure { .. })
        ));
    }
}

#[test]
fn empty_singleton_and_duplicate_clouds_remain_distinct_inputs() {
    for (coordinates, n, expected) in [(vec![], 0, 0), (vec![0., 0.], 1, 1), (vec![0.; 6], 3, 1)] {
        let diagram = rips_from_points(
            PointCloudView::new(&coordinates, n, 2).unwrap(),
            &RipsOptions::default(),
        )
        .unwrap();
        assert_eq!(diagram.len(), expected);
        assert_eq!(diagram.coverage(), Coverage::Complete);
    }
}

#[test]
fn translation_rotation_and_concurrent_calls_preserve_diagrams() {
    let original = [0., 0., 3., 0., 3., 4., 0., 4.];
    let transformed: Vec<_> = original
        .as_chunks::<2>()
        .0
        .iter()
        .flat_map(|p| [10. - p[1], 7. + p[0]])
        .collect();
    let options = RipsOptions::default();
    let input = PointCloudView::new(&original, 4, 2).unwrap();
    let expected = rips_from_points(input, &options).unwrap();
    assert_eq!(
        rips_from_points(PointCloudView::new(&transformed, 4, 2).unwrap(), &options).unwrap(),
        expected
    );
    std::thread::scope(|scope| {
        let handles: Vec<_> = (0..8)
            .map(|_| scope.spawn(|| rips_from_points(input, &options).unwrap()))
            .collect();
        for handle in handles {
            assert_eq!(handle.join().unwrap(), expected);
        }
    });
}

#[test]
fn hausdorff_hand_cases_keep_support_but_ignore_order_and_multiplicity() {
    let first = PointCloudView::new(&[0., 2.], 2, 1).unwrap();
    let denser = PointCloudView::new(&[0., 1., 2.], 3, 1).unwrap();
    let duplicate = PointCloudView::new(&[2., 0., 2.], 3, 1).unwrap();
    assert_eq!(hausdorff_distance(first, denser), Ok(1.));
    assert_eq!(hausdorff_distance(denser, first), Ok(1.));
    assert_eq!(hausdorff_distance(first, duplicate), Ok(0.));
    let origin = PointCloudView::new(&[0., 0.], 1, 2).unwrap();
    let shifted = PointCloudView::new(&[3., 4.], 1, 2).unwrap();
    assert_eq!(hausdorff_distance(origin, shifted), Ok(5.));
    assert_eq!(origin.coordinates(), &[0., 0.]);
}

#[test]
fn hausdorff_empty_and_incompatible_inputs_have_distinct_meanings() {
    let empty = PointCloudView::new(&[], 0, 1).unwrap();
    let point = PointCloudView::new(&[-0.], 1, 1).unwrap();
    assert_eq!(hausdorff_distance(empty, empty), Ok(0.));
    assert_eq!(hausdorff_distance(empty, point), Ok(f64::INFINITY));
    assert_eq!(hausdorff_distance(point, empty), Ok(f64::INFINITY));
    assert_eq!(hausdorff_distance(point, point).unwrap().to_bits(), 0);
    assert!(matches!(
        hausdorff_distance(empty, PointCloudView::new(&[], 0, 2).unwrap()),
        Err(Error::ShapeMismatch { .. })
    ));
    assert!(PointCloudView::new(&[f64::NAN], 1, 1).is_err());
    assert!(PointCloudView::new(&[f64::INFINITY], 1, 1).is_err());
}

#[test]
fn hausdorff_stable_norm_and_overflowing_candidates() {
    for scale in [1e154, 1e-200] {
        let origin = PointCloudView::new(&[0., 0.], 1, 2).unwrap();
        let coordinates = [3. * scale, 4. * scale];
        let other = PointCloudView::new(&coordinates, 1, 2).unwrap();
        let distance = hausdorff_distance(origin, other).unwrap();
        assert!((distance / (5. * scale) - 1.).abs() < 1e-14);
    }
    let negative = PointCloudView::new(&[-f64::MAX], 1, 1).unwrap();
    let positive = PointCloudView::new(&[f64::MAX], 1, 1).unwrap();
    assert!(matches!(
        hausdorff_distance(negative, positive),
        Err(Error::NumericalFailure { .. })
    ));
    // Distant candidates overflow, but every point has a representable nearest neighbor.
    let both = PointCloudView::new(&[-f64::MAX, f64::MAX], 2, 1).unwrap();
    assert_eq!(hausdorff_distance(both, both), Ok(0.));
}

#[test]
fn hausdorff_matches_independent_integer_squared_distance_table() {
    // A fixed 2-D grid, every nonempty subset on each side. Expectations use
    // integer squared distances and take one square root after selecting extrema.
    let grid = [[0_i64, 0], [1, 0], [0, 2], [3, 4]];
    for left_mask in 1..16 {
        for right_mask in 1..16 {
            let left: Vec<_> = (0..4).filter(|i| left_mask & (1 << i) != 0).collect();
            let right: Vec<_> = (0..4).filter(|i| right_mask & (1 << i) != 0).collect();
            let squared: Vec<Vec<i64>> = left
                .iter()
                .map(|&i| {
                    right
                        .iter()
                        .map(|&j| {
                            let dx = grid[i][0] - grid[j][0];
                            let dy = grid[i][1] - grid[j][1];
                            dx * dx + dy * dy
                        })
                        .collect()
                })
                .collect();
            let row_max = squared
                .iter()
                .map(|row| *row.iter().min().unwrap())
                .max()
                .unwrap();
            let col_max = (0..right.len())
                .map(|j| squared.iter().map(|row| row[j]).min().unwrap())
                .max()
                .unwrap();
            let expected = (row_max.max(col_max) as f64).sqrt();
            let a: Vec<f64> = left
                .iter()
                .flat_map(|&i| grid[i].map(|x| x as f64))
                .collect();
            let b: Vec<f64> = right
                .iter()
                .flat_map(|&i| grid[i].map(|x| x as f64))
                .collect();
            let actual = hausdorff_distance(
                PointCloudView::new(&a, left.len(), 2).unwrap(),
                PointCloudView::new(&b, right.len(), 2).unwrap(),
            )
            .unwrap();
            assert!(
                (actual - expected).abs() <= 1e-14,
                "masks {left_mask}/{right_mask}"
            );
        }
    }
}

#[test]
fn hausdorff_controls_charge_both_scans_and_retry_without_mutation() {
    use std::sync::atomic::{AtomicBool, Ordering};
    let point = PointCloudView::new(&[0., 0.], 1, 2).unwrap();
    let other = PointCloudView::new(&[3., 4.], 1, 2).unwrap();
    // One compatibility check plus two directed pairs of two coordinates.
    let insufficient = Execution::default().max_work(4);
    assert_eq!(
        hausdorff_distance_with(point, other, &insufficient),
        Err(Error::WorkLimitExceeded { limit: 4 })
    );
    let enough = Execution::default().max_work(5);
    assert_eq!(hausdorff_distance_with(point, other, &enough), Ok(5.));
    assert_eq!(hausdorff_distance_with(point, other, &enough), Ok(5.));
    let empty = PointCloudView::new(&[], 0, 2).unwrap();
    assert_eq!(
        hausdorff_distance_with(empty, empty, &Execution::default().max_work(0)),
        Err(Error::WorkLimitExceeded { limit: 0 })
    );
    let flag = AtomicBool::new(true);
    assert_eq!(
        hausdorff_distance_with(point, other, &Execution::default().cancellation(&flag)),
        Err(Error::Cancelled)
    );
    assert!(flag.load(Ordering::Relaxed));
    assert_eq!(point.coordinates(), &[0., 0.]);
}
