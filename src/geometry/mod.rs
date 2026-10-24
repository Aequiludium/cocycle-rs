//! Borrowed, validated numerical inputs.
//!
//! Point clouds use row-major coordinates. `DissimilarityView` uses a condensed
//! lower triangle: `[d(1,0), d(2,0), d(2,1), ...]`. `DissimilarityMatrixView` adds
//! explicit upper-triangle and full-square layouts. No view input is copied,
//! deduplicated or modified. Finite coordinates and nonnegative dissimilarities
//! do not, on their own, certify a metric.
//!
//! [`hausdorff_distance`] compares nonempty Euclidean point sets directly, without
//! persistence. [`hausdorff_distance_matrix`] retains ordered sample identities.

mod dissimilarity;
mod euclidean;
mod point_cloud;

pub use dissimilarity::DissimilarityView;
pub(crate) use dissimilarity::pair_count;
pub use distance::{
    HausdorffDistanceMatrix, hausdorff_distance, hausdorff_distance_matrix,
    hausdorff_distance_matrix_with, hausdorff_distance_with,
};
pub(crate) use euclidean::{euclidean_distance, euclidean_distances, euclidean_distances_with};
pub use point_cloud::PointCloudView;

pub(crate) mod distance;
mod matrix;
mod metric;
pub use matrix::{DissimilarityMatrixView, MatrixLayout};
pub(crate) use metric::validate_metric_with_checkpoints;
pub use metric::{MetricPolicy, MetricValidation, validate_metric};
