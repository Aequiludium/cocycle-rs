//! Independent analytic and execution contracts for the bounded #145 rewrite.
use cocycle::Error;
use cocycle::execution::Execution;
use cocycle::persistent_laplacian::{
    HarmonicRestriction, SpectralOptions, gram_spectrum, split_spectrum,
};
use std::sync::atomic::{AtomicBool, Ordering};

fn close(actual: &[f64], expected: &[f64]) {
    assert_eq!(actual.len(), expected.len());
    for (&a, &e) in actual.iter().zip(expected) {
        assert!((a - e).abs() < 1e-10, "{a} != {e}");
    }
}
fn identity(n: usize) -> Vec<f64> {
    (0..n * n)
        .map(|i| if i / n == i % n { 1. } else { 0. })
        .collect()
}

#[test]
fn rectangular_grams_preserve_original_zeros_and_repeated_values() {
    let o = SpectralOptions::default();
    let e = Execution::default();
    let a = [-1., -1., 0., 1., 0., -1., 0., 1., 1.];
    let g = gram_spectrum(&a, 3, 3, o, &e).unwrap();
    close(&g.reconstructed, &[0., 3., 3.]);
    assert_eq!(g.numerical_rank, 2);
    assert!(g.lifted_residual < 1e-12);
    assert!(g.orthogonality_error < 1e-12);
    let tall = gram_spectrum(&[-1., 0., 1., -1., 0., 1., 0., 0.], 4, 2, o, &e).unwrap();
    close(&tall.reconstructed, &[1., 3.]);
    assert_eq!(tall.dual_spectrum.len(), 4);
    let wide = gram_spectrum(&[1., 0., 0., 1., 0., 0.], 2, 3, o, &e).unwrap();
    close(&wide.reconstructed, &[0., 0., 2.]);
}

#[test]
fn one_sided_gram_lifts_modes_without_the_large_square() {
    use cocycle::persistent_laplacian::gram_eigenpairs;
    let e = Execution::default();
    // A A^T = diag(1,4), while A^T A has 298 additional zero modes.
    let mut wide = vec![0.; 600];
    wide[0] = 1.;
    wide[301] = 2.;
    let options = SpectralOptions {
        max_dense_entries: 600,
        ..Default::default()
    };
    assert!(gram_spectrum(&wide, 2, 300, options, &e).is_err());
    let g = gram_eigenpairs(&wide, 2, 300, options, &e).unwrap();
    assert!(g.used_dual);
    assert_eq!(g.numerical_rank, 2);
    close(&g.eigenvalues[298..], &[1., 4.]);
    assert!(g.eigenvalues[..298].iter().all(|v| *v == 0.));
    assert_eq!(g.raw_spectrum, vec![1., 4.]);
    assert!(g.residual < 1e-12 && g.orthogonality_error < 1e-12);
    for (m, n) in [(0, 3), (3, 0), (0, 0), (2, 5)] {
        let g = gram_eigenpairs(&vec![0.; m * n], m, n, options, &e).unwrap();
        assert_eq!(g.eigenvalues, vec![0.; n]);
        assert!(g.vectors.is_empty());
    }
    let tall = gram_eigenpairs(&[-1., 0., 1., -1., 0., 1., 0., 0.], 4, 2, options, &e).unwrap();
    assert!(!tall.used_dual);
    close(&tall.eigenvalues, &[1., 3.]);
    assert!(gram_eigenpairs(&[f64::NAN], 1, 1, options, &e).is_err());
    assert!(gram_eigenpairs(&[], usize::MAX, 2, options, &e).is_err());
    assert!(gram_eigenpairs(&wide, 2, 300, options, &Execution::default().max_work(1)).is_err());
}

#[test]
fn empty_shapes_and_zero_operators_are_explicit() {
    for (m, n) in [(0, 3), (3, 0), (0, 0), (2, 5)] {
        let g = gram_spectrum(
            &vec![0.; m * n],
            m,
            n,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        assert_eq!(g.reconstructed, vec![0.; n]);
        assert_eq!(g.primal.len(), n * n);
        assert_eq!(g.dual.len(), m * m);
    }
    let r = HarmonicRestriction::compute(
        &[0.; 9],
        3,
        &identity(3),
        3,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    assert!(r.eigenvalues().is_empty());
    assert!(r.complement().is_empty());
    assert!(
        r.smallest_eigenpairs(&[1., 0., 0.], 1, 1, 10, 1e-9, &Execution::default())
            .is_err()
    );
    assert!(
        r.smallest_eigenpairs(&[], 0, 0, 0, 1e-9, &Execution::default())
            .is_err()
    );
    let empty = HarmonicRestriction::compute(
        &[],
        0,
        &[],
        0,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    assert_eq!(empty.dimension(), 0);
    assert!(empty.eigenvalues().is_empty());
}

#[test]
fn triangle_union_does_not_duplicate_component_zeros() {
    let a = [-1., -1., 0., 1., 0., -1., 0., 1., 1.];
    let f = [1., -1., 1.];
    let s = split_spectrum(
        &a,
        3,
        3,
        &f,
        1,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    close(&s.lower.primal_spectrum, &[0., 3., 3.]);
    close(&s.upper.primal_spectrum, &[0., 0., 3.]);
    close(&s.eigenvalues, &[3., 3., 3.]);
    close(&s.matrix, &[3., 0., 0., 0., 3., 0., 0., 0., 3.]);
    assert_eq!(s.chain_residual, 0.);
    assert!(
        split_spectrum(
            &a,
            3,
            3,
            &[1., 1., 1.],
            1,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
}

#[test]
fn homology_representative_is_rejected_but_actual_harmonic_kernel_works() {
    let h = [1., -1., -1., 1.];
    assert!(
        HarmonicRestriction::compute(
            &h,
            2,
            &[1., 0.],
            1,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    let r = HarmonicRestriction::compute(
        &h,
        2,
        &[2., 2.],
        1,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    close(r.eigenvalues(), &[2.]);
    let p = r
        .smallest_eigenpairs(&[1., -1.], 1, 1, 120, 1e-9, &Execution::default())
        .unwrap();
    close(&p.values, &[2.]);
    assert!(p.residuals[0] < 1e-12);
    assert!(p.range_leakage < 1e-12);
    assert!(!p.is_index_coverage_certified());
}

#[test]
fn nonorthogonal_kernel_is_normalized_without_using_its_raw_gram_as_projection() {
    let h = [0., 0., 0., 0., 0., 0., 0., 0., 2.];
    let n = [2., 1., 0., 3., 0., 0.];
    let r = HarmonicRestriction::compute(
        &h,
        3,
        &n,
        2,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    close(r.eigenvalues(), &[2.]);
    assert!(r.kernel_residual() < 1e-12);
    assert!(r.orthogonality_error() < 1e-12);
}

#[test]
fn incomplete_dependent_and_negative_domains_return_errors() {
    let o = SpectralOptions::default();
    let e = Execution::default();
    assert!(HarmonicRestriction::compute(&[0., 0., 0., 1.], 2, &[], 0, o, &e).is_err());
    assert!(
        HarmonicRestriction::compute(&[0., 0., 0., 1.], 2, &[1., 1., 0., 0.], 2, o, &e).is_err()
    );
    assert!(HarmonicRestriction::compute(&[-1., 0., 0., 1.], 2, &[], 0, o, &e).is_err());
    assert!(HarmonicRestriction::compute(&[1., 1., 0., 1.], 2, &[], 0, o, &e).is_err());
}

#[test]
fn tiny_positive_remains_in_complete_spectrum_and_partial_request_is_ambiguous() {
    let h = [0., 0., 0., 0., 1e-14, 0., 0., 0., 1.];
    let r = HarmonicRestriction::compute(
        &h,
        3,
        &[1., 0., 0.],
        1,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    assert_eq!(r.eigenvalues()[0], 1e-14);
    assert!(r.rank_ambiguous());
    assert!(
        r.smallest_eigenpairs(&[0., 1., 1.], 1, 1, 120, 1e-9, &Execution::default())
            .is_err()
    );
    let g = gram_spectrum(
        &[1e-7, 0., 0., 1.],
        2,
        2,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    assert!(g.primal_spectrum[0] > 0.);
    assert_eq!(g.numerical_rank, 1);
    assert!(g.rank_ambiguous);
}

#[test]
fn requested_count_exhaustion_and_search_rank_loss_are_errors() {
    let h = [1., 0., 0., 0., 2., 0., 0., 0., 3.];
    let r = HarmonicRestriction::compute(
        &h,
        3,
        &[],
        0,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    assert!(
        r.smallest_eigenpairs(&[1.; 6], 2, 3, 120, 1e-9, &Execution::default())
            .is_err()
    );
    assert!(
        r.smallest_eigenpairs(&[1.; 6], 2, 2, 120, 1e-9, &Execution::default())
            .is_err()
    );
    assert!(
        r.smallest_eigenpairs(&[1.; 3], 1, 1, 0, 1e-9, &Execution::default())
            .is_err()
    );
    let p = r
        .smallest_eigenpairs(&[1., 0., 0.], 1, 1, 0, 1e-9, &Execution::default())
        .unwrap();
    assert_eq!(p.iterations, 0);
    close(&p.values, &[1.]);
}

#[test]
fn input_ownership_is_independent_of_callers_and_resource_failure_is_recoverable() {
    let mut h = vec![1., -1., -1., 1.];
    let mut kernel = vec![1., 1.];
    let r = HarmonicRestriction::compute(
        &h,
        2,
        &kernel,
        1,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    h.fill(9.);
    kernel.fill(0.);
    close(r.original_matrix(), &[1., -1., -1., 1.]);
    let flag = AtomicBool::new(true);
    assert!(matches!(
        r.smallest_eigenpairs(
            &[1., -1.],
            1,
            1,
            10,
            1e-9,
            &Execution::default().cancellation(&flag)
        ),
        Err(Error::Cancelled)
    ));
    flag.store(false, Ordering::Relaxed);
    assert!(matches!(
        r.smallest_eigenpairs(
            &[1., -1.],
            1,
            1,
            10,
            1e-9,
            &Execution::default().max_work(3)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
    assert!(
        r.smallest_eigenpairs(&[1., -1.], 1, 1, 10, 1e-9, &Execution::default())
            .is_ok()
    );
}

#[test]
fn shapes_nonfinite_options_overflow_and_storage_limits_are_checked() {
    let e = Execution::default();
    let o = SpectralOptions::default();
    for a in [&[f64::NAN][..], &[f64::INFINITY][..]] {
        assert!(matches!(
            gram_spectrum(a, 1, 1, o, &e),
            Err(Error::NonFiniteValue { .. })
        ));
    }
    assert!(gram_spectrum(&[f64::MAX], 1, 1, o, &e).is_err());
    assert!(matches!(
        gram_spectrum(&[], 1, 1, o, &e),
        Err(Error::ShapeMismatch { .. })
    ));
    assert!(matches!(
        gram_spectrum(&[], usize::MAX, 0, o, &e),
        Err(Error::SizeOverflow { .. })
    ));
    assert!(
        gram_spectrum(
            &[1.],
            1,
            1,
            SpectralOptions {
                max_dense_entries: 0,
                ..o
            },
            &e
        )
        .is_err()
    );
    assert!(
        gram_spectrum(
            &[1.],
            1,
            1,
            SpectralOptions {
                relative_tolerance: f64::NAN,
                ..o
            },
            &e
        )
        .is_err()
    );
    assert!(
        HarmonicRestriction::compute(
            &identity(2),
            2,
            &[],
            0,
            SpectralOptions {
                max_dense_entries: 3,
                ..o
            },
            &e
        )
        .is_err()
    );
}
