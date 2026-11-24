//! Public contracts and hand-derived persistent Laplacian examples.
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::filtration::RipsBuilder;
use cocycle::geometry::PointCloudView;
use cocycle::persistent_laplacian::{
    LaplacianComponent, LsmrOptions, OperatorOptions, PersistentLaplacian,
    PersistentLaplacianOperator, SpectralOptions, fixed_endpoint_laplacians,
    fixed_endpoint_up_laplacians, project_boundary_kernel_lsmr, relevant_pair_laplacians,
};
use cocycle::{Error, Result};
use std::sync::atomic::AtomicBool;

fn complex(cells: &[(&[usize], f64)]) -> SimplicialComplex {
    SimplicialComplex::new(
        cells
            .iter()
            .map(|(v, t)| Simplex::new(v.to_vec(), *t).unwrap())
            .collect(),
    )
    .unwrap()
}
fn close(a: &[f64], b: &[f64]) {
    assert_eq!(a.len(), b.len());
    for (x, y) in a.iter().zip(b) {
        assert!((x - y).abs() < 1e-9, "{a:?} != {b:?}");
    }
}

#[test]
fn partial_spectra_use_fixed_actions_and_current_operator_for_reused_vectors() {
    use cocycle::persistent_laplacian::{IterativeOptions, PersistentSpectrumTarget as Target};
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 0.),
        (&[0, 1], 1.),
        (&[1, 2], 1.),
    ]);
    let e = Execution::default();
    let op = PersistentLaplacianOperator::new(
        &source.stage(1.).unwrap(),
        &source.stage(1.).unwrap(),
        0,
        OperatorOptions::default(),
        &e,
    )
    .unwrap();
    let options = IterativeOptions {
        block_size: 3,
        ..Default::default()
    };
    let small = op
        .eigenpairs(3, Target::Smallest, None, &options, &e)
        .unwrap();
    close(&small.values, &[0., 1., 3.]);
    let large = op
        .eigenpairs(2, Target::Largest, None, &options, &e)
        .unwrap();
    close(&large.values, &[3., 1.]);
    let positive = op
        .eigenpairs(2, Target::SmallestPositive, None, &options, &e)
        .unwrap();
    close(&positive.values, &[1., 3.]);
    assert_eq!(positive.range_rank, Some(2));
    assert!(positive.range_residual < 1e-12);
    assert!(positive.residuals.iter().all(|v| *v < 1e-10));
    assert!(!positive.is_index_coverage_certified());
    // Same source basis, changed endpoint: only vectors, not old Ritz values, transfer.
    let old = PersistentLaplacianOperator::new(
        &source.stage(0.).unwrap(),
        &source.stage(0.).unwrap(),
        0,
        OperatorOptions::default(),
        &e,
    )
    .unwrap();
    let seed = old
        .eigenpairs(3, Target::Smallest, None, &options, &e)
        .unwrap();
    let warm = op
        .eigenpairs(3, Target::Smallest, Some(&seed.vectors), &options, &e)
        .unwrap();
    close(&warm.values, &[0., 1., 3.]);
    let warm_positive = op
        .eigenpairs(
            2,
            Target::SmallestPositive,
            Some(&positive.vectors),
            &IterativeOptions {
                block_size: 2,
                ..options
            },
            &e,
        )
        .unwrap();
    close(&warm_positive.values, &[1., 3.]);
}

#[test]
fn positive_spectrum_exceeds_dense_cap_without_first_computing_zero_modes() {
    use cocycle::persistent_laplacian::{IterativeOptions, PersistentSpectrumTarget as Target};
    let e = Execution::default();
    let op = {
        let mut cells: Vec<_> = (0..300)
            .map(|v| Simplex::new(vec![v], 0.).unwrap())
            .collect();
        cells.push(Simplex::new(vec![0, 1], 1.).unwrap());
        let source = SimplicialComplex::new(cells).unwrap();
        PersistentLaplacianOperator::new(
            &source.stage(1.).unwrap(),
            &source.stage(1.).unwrap(),
            0,
            OperatorOptions {
                spectral: SpectralOptions {
                    max_dense_entries: 0,
                    ..Default::default()
                },
                ..Default::default()
            },
            &e,
        )
        .unwrap()
    };
    assert!(op.materialize(90_000, &e).is_err());
    let options = IterativeOptions {
        block_size: 2,
        max_basis_entries: 600,
        ..Default::default()
    };
    // 299 exact harmonic dimensions, one positive eigenvalue 2 from a single edge.
    let result = op
        .eigenpairs(1, Target::SmallestPositive, None, &options, &e)
        .unwrap();
    close(&result.values, &[2.]);
    assert_eq!(result.basis.len(), 300);
    assert_eq!(result.range_rank, Some(1));
    assert_eq!(result.vectors.len(), 300);
    assert!(
        op.eigenpairs(2, Target::SmallestPositive, None, &options, &e)
            .is_err()
    );
    assert!(
        op.eigenpairs(
            1,
            Target::SmallestPositive,
            None,
            &IterativeOptions {
                max_basis_entries: 299,
                ..options
            },
            &e
        )
        .is_err()
    );
    assert!(
        op.eigenpairs(
            1,
            Target::SmallestPositive,
            None,
            &IterativeOptions {
                absolute_tolerance: 3.,
                ..options
            },
            &e
        )
        .is_err()
    );
    for input in [vec![0.; 299], vec![f64::NAN; 300], vec![f64::INFINITY; 300]] {
        assert!(
            op.eigenpairs(1, Target::SmallestPositive, Some(&input), &options, &e)
                .is_err()
        );
    }
    let flag = AtomicBool::new(true);
    assert!(matches!(
        op.eigenpairs(
            1,
            Target::SmallestPositive,
            None,
            &options,
            &Execution::default().cancellation(&flag)
        ),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        op.eigenpairs(
            1,
            Target::SmallestPositive,
            None,
            &options,
            &Execution::default().max_work(1)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
}

#[test]
fn positive_range_does_not_normalize_annihilated_upper_roundoff() {
    use cocycle::persistent_laplacian::{IterativeOptions, PersistentSpectrumTarget as Target};
    let source = complex(&[
        (&[0], 1.),
        (&[1], 0.),
        (&[2], 0.),
        (&[3], 0.),
        (&[0, 1], 2.),
        (&[0, 2], 1.),
        (&[0, 3], 2.),
        (&[1, 2], 0.),
        (&[1, 3], 2.),
        (&[2, 3], 0.),
        (&[0, 1, 2], 3.),
        (&[0, 1, 3], 2.),
        (&[0, 2, 3], 4.),
        (&[1, 2, 3], 3.),
        (&[0, 1, 2, 3], 4.),
    ]);
    let e = Execution::default();
    let op = PersistentLaplacianOperator::new(
        &source.stage(0.).unwrap(),
        &source.stage(3.).unwrap(),
        1,
        OperatorOptions::default(),
        &e,
    )
    .unwrap();
    // B2 is injective. Only the two tree-edge lower boundary images remain.
    let result = op
        .eigenpairs(
            2,
            Target::SmallestPositive,
            None,
            &IterativeOptions::default(),
            &e,
        )
        .unwrap();
    close(&result.values, &[1., 3.]);
    assert_eq!(result.range_rank, Some(2));
}

#[test]
fn incremental_start_projection_and_batch_spectra_preserve_every_pair() {
    use cocycle::persistent_laplacian::{
        IterativeOptions, PersistentSpectrumTarget as Target, fixed_endpoint_eigenpairs,
    };
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 1.),
        (&[3], 2.),
        (&[0, 1], 2.),
        (&[1, 2], 2.),
        (&[2, 3], 2.),
    ]);
    let e = Execution::default();
    let stages: Vec<_> = [0., 1., 1., 2.]
        .into_iter()
        .map(|s| source.stage(s).unwrap())
        .collect();
    let endpoint = source.stage(2.).unwrap();
    let options = IterativeOptions {
        block_size: 4,
        ..Default::default()
    };
    let mut op =
        PersistentLaplacianOperator::new(&endpoint, &endpoint, 0, OperatorOptions::default(), &e)
            .unwrap();
    for stage in stages.iter().rev() {
        op = op.restrict_start(stage.scale(), &e).unwrap();
        let direct =
            PersistentLaplacian::compute(stage, &endpoint, 0, SpectralOptions::default(), &e)
                .unwrap();
        close(&op.materialize(16, &e).unwrap(), direct.matrix());
        assert_eq!(op.basis(), direct.basis());
    }
    assert!(op.clone().restrict_start(1., &e).is_err());
    assert!(op.restrict_start(f64::NAN, &e).is_err());
    for reuse in [false, true] {
        let batch = fixed_endpoint_eigenpairs(
            &stages,
            &endpoint,
            0,
            1,
            Target::SmallestPositive,
            OperatorOptions::default(),
            &options,
            reuse,
            &e,
        )
        .unwrap();
        for (stage, pairs) in stages.iter().zip(&batch) {
            let direct =
                PersistentLaplacian::compute(stage, &endpoint, 0, SpectralOptions::default(), &e)
                    .unwrap();
            let spectrum = direct.spectrum(&e).unwrap();
            let expected = spectrum
                .eigenvalues
                .iter()
                .copied()
                .find(|v| *v > 1e-8)
                .unwrap();
            close(&pairs.values, &[expected]);
            assert_eq!(pairs.scales, (stage.scale(), 2.));
        }
    }
    assert!(
        fixed_endpoint_eigenpairs(
            &[stages[3], stages[0]],
            &endpoint,
            0,
            1,
            Target::Largest,
            OperatorOptions::default(),
            &options,
            true,
            &e
        )
        .is_err()
    );
    assert!(matches!(
        fixed_endpoint_eigenpairs(
            &stages,
            &endpoint,
            0,
            1,
            Target::Largest,
            OperatorOptions::default(),
            &options,
            true,
            &Execution::default().max_work(10)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
}

#[test]
fn compact_nonbranching_projection_uses_linear_storage_and_the_same_spectral_path() {
    use cocycle::persistent_laplacian::{IterativeOptions, PersistentSpectrumTarget as Target};
    let mut cells: Vec<_> = (0..300)
        .map(|v| Simplex::new(vec![v], v as f64).unwrap())
        .collect();
    for v in 0..299 {
        cells.push(Simplex::new(vec![v, v + 1], 300.).unwrap());
    }
    let source = SimplicialComplex::new(cells).unwrap();
    let e = Execution::default();
    let options = OperatorOptions {
        max_projection_entries: 600,
        spectral: SpectralOptions {
            max_dense_entries: 0,
            ..Default::default()
        },
        ..Default::default()
    };
    let start = source.stage(1.).unwrap();
    let end = source.stage(300.).unwrap();
    assert!(PersistentLaplacianOperator::new(&start, &end, 0, options, &e).is_err());
    let op = PersistentLaplacianOperator::new_non_branching(&start, &end, 0, options, &e).unwrap();
    assert!(op.projection_entries() <= 600);
    // Every outside path vertex forces its incident edge flow to zero.
    close(&op.action(&[1., -1.], &e).unwrap(), &[2., -2.]);
    let pairs = op
        .eigenpairs(
            1,
            Target::SmallestPositive,
            None,
            &IterativeOptions::default(),
            &e,
        )
        .unwrap();
    close(&pairs.values, &[2.]);
    let smaller = op.restrict_start(0., &e).unwrap();
    close(&smaller.action(&[1.], &e).unwrap(), &[0.]);
    assert_eq!(smaller.projection_rank(), 299);
    let triangle = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 0.),
        (&[0, 1], 1.),
        (&[0, 2], 1.),
        (&[1, 2], 1.),
    ]);
    let all = triangle.stage(1.).unwrap();
    let compact = PersistentLaplacianOperator::new_non_branching(
        &all,
        &all,
        0,
        OperatorOptions::default(),
        &e,
    )
    .unwrap()
    .restrict_start(0., &e)
    .unwrap();
    close(&compact.action(&[1., -1., 0.], &e).unwrap(), &[3., -3., 0.]);
    assert!(
        PersistentLaplacianOperator::new_non_branching(
            &start,
            &end,
            0,
            OperatorOptions {
                max_projection_entries: 1,
                ..options
            },
            &e
        )
        .is_err()
    );
}

#[test]
fn relevant_partial_spectra_reuse_only_upper_terms_and_keep_request_order() {
    use cocycle::persistent_laplacian::{
        IterativeOptions, PersistentSpectrumTarget as Target, relevant_pair_eigenpairs,
    };
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 0.),
        (&[3], 1.),
        (&[0, 1], 0.),
        (&[1, 2], 0.),
        (&[0, 2], 1.),
        (&[2, 3], 2.),
        (&[0, 1, 2], 3.),
    ]);
    let stages: Vec<_> = [0., 1., 1., 2., 3.]
        .into_iter()
        .map(|s| source.stage(s).unwrap())
        .collect();
    let mut requests: Vec<_> = (0..stages.len())
        .flat_map(|t| (0..=t).map(move |s| (s, t)))
        .collect();
    requests.reverse();
    requests.push((0, 0));
    let options = IterativeOptions {
        block_size: 8,
        ..Default::default()
    };
    let e = Execution::default();
    for degree in [0, 1] {
        for reuse in [false, true] {
            let batch = relevant_pair_eigenpairs(
                &stages,
                &requests,
                degree,
                1,
                Target::Largest,
                OperatorOptions::default(),
                &options,
                reuse,
                &e,
            )
            .unwrap();
            for (&(s, t), result) in requests.iter().zip(&batch) {
                let direct = PersistentLaplacian::compute(
                    &stages[s],
                    &stages[t],
                    degree,
                    SpectralOptions::default(),
                    &e,
                )
                .unwrap();
                close(
                    &result.values,
                    &[*direct.spectrum(&e).unwrap().eigenvalues.last().unwrap()],
                );
                assert_eq!(result.basis, direct.basis());
                assert_eq!(result.target_basis, direct.target_basis());
                assert_eq!(result.scales, (stages[s].scale(), stages[t].scale()));
            }
        }
    }
    assert!(
        relevant_pair_eigenpairs(
            &stages,
            &[(4, 0)],
            0,
            1,
            Target::Largest,
            OperatorOptions::default(),
            &options,
            false,
            &e
        )
        .is_err()
    );
    assert!(
        relevant_pair_eigenpairs(
            &stages,
            &requests,
            0,
            1,
            Target::Largest,
            OperatorOptions::default(),
            &IterativeOptions {
                max_basis_entries: 4,
                ..options
            },
            false,
            &e
        )
        .is_err()
    );
}

#[test]
fn lsmr_returns_minimum_norm_fit_and_inconsistent_residual() {
    // M=B2^T has dependent columns [1,2,0] and [2,4,0].
    // Its minimum-norm fit of [3,6,7] is [3/5,6/5], residual [0,0,7].
    let b2 = [1., 2., 0., 2., 4., 0.];
    let result = project_boundary_kernel_lsmr(
        &b2,
        2,
        3,
        &[3., 6., 7.],
        LsmrOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    close(&result.solution, &[0.6, 1.2]);
    close(&result.projected, &[0., 0., 7.]);
    assert!((result.residual_norm - 7.).abs() < 1e-10);
    assert!(result.normal_residual_norm < 1e-10);
    assert!(result.iterations > 0);
}

#[test]
fn lsmr_empty_zero_and_stationary_inputs_have_explicit_shapes() {
    for (b2, rows, columns, y, expected_x, expected_z) in [
        (vec![], 0, 2, vec![1., 2.], vec![], vec![1., 2.]),
        (vec![], 2, 0, vec![], vec![0., 0.], vec![]),
        (
            vec![0.; 6],
            2,
            3,
            vec![1., 2., 3.],
            vec![0., 0.],
            vec![1., 2., 3.],
        ),
        (vec![1., 0.], 1, 2, vec![0., 3.], vec![0.], vec![0., 3.]),
    ] {
        let result = project_boundary_kernel_lsmr(
            &b2,
            rows,
            columns,
            &y,
            LsmrOptions {
                max_iterations: Some(0),
                ..Default::default()
            },
            &Execution::default(),
        )
        .unwrap();
        close(&result.solution, &expected_x);
        close(&result.projected, &expected_z);
        assert_eq!(result.iterations, 0);
        assert_eq!(result.normal_residual_norm, 0.);
    }
}

#[test]
fn lsmr_rejects_invalid_inputs_caps_and_unconverged_output() {
    let matrix = [1., 0., 0., 0., 2., 0., 0., 0., 4.];
    for limit in [0, 1] {
        assert!(matches!(
            project_boundary_kernel_lsmr(
                &matrix,
                3,
                3,
                &[1.; 3],
                LsmrOptions {
                    max_iterations: Some(limit),
                    ..Default::default()
                },
                &Execution::default()
            ),
            Err(Error::NumericalFailure {
                context: "LSMR iteration limit"
            })
        ));
    }
    for options in [
        LsmrOptions {
            max_matrix_entries: 8,
            ..Default::default()
        },
        LsmrOptions {
            max_vector_entries: 2,
            ..Default::default()
        },
        LsmrOptions {
            relative_tolerance: f64::NAN,
            ..Default::default()
        },
        LsmrOptions {
            relative_tolerance: 0.,
            ..Default::default()
        },
    ] {
        assert!(matches!(
            project_boundary_kernel_lsmr(&matrix, 3, 3, &[1.; 3], options, &Execution::default()),
            Err(Error::InvalidParameter { .. })
        ));
    }
    assert!(matches!(
        project_boundary_kernel_lsmr(
            &[],
            usize::MAX,
            2,
            &[],
            LsmrOptions::default(),
            &Execution::default()
        ),
        Err(Error::SizeOverflow { .. })
    ));
    assert!(matches!(
        project_boundary_kernel_lsmr(
            &[],
            1,
            1,
            &[1.],
            LsmrOptions::default(),
            &Execution::default()
        ),
        Err(Error::ShapeMismatch { .. })
    ));
    assert!(matches!(
        project_boundary_kernel_lsmr(
            &[1.],
            1,
            1,
            &[],
            LsmrOptions::default(),
            &Execution::default()
        ),
        Err(Error::ShapeMismatch { .. })
    ));
    for (matrix, y) in [([f64::INFINITY], [1.]), ([1.], [f64::NAN])] {
        assert!(matches!(
            project_boundary_kernel_lsmr(
                &matrix,
                1,
                1,
                &y,
                LsmrOptions::default(),
                &Execution::default()
            ),
            Err(Error::NonFiniteValue { .. })
        ));
    }
    assert!(matches!(
        project_boundary_kernel_lsmr(
            &[f64::MAX],
            1,
            1,
            &[f64::MAX],
            LsmrOptions::default(),
            &Execution::default()
        ),
        Err(Error::NumericalFailure { .. })
    ));
    // The raw paper recurrence forms products of rotation scales. This simple
    // rewrite rejects overflow/underflow, rather than silently losing an update.
    for (boundary, rhs, context) in [
        (1e200, 1e-100, "LSMR nonfinite arithmetic"),
        (1e-200, 1e100, "LSMR rotation scale underflow"),
        (1e-200, 1e-200, "LSMR normal scale underflow"),
    ] {
        assert!(
            matches!(project_boundary_kernel_lsmr(&[boundary], 1, 1, &[rhs],
            LsmrOptions::default(), &Execution::default()),
            Err(Error::NumericalFailure { context: actual }) if actual == context)
        );
    }
    assert!(matches!(
        project_boundary_kernel_lsmr(
            &matrix,
            3,
            3,
            &[1.; 3],
            LsmrOptions::default(),
            &Execution::default().max_work(0)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let flag = AtomicBool::new(true);
    assert!(matches!(
        project_boundary_kernel_lsmr(
            &matrix,
            3,
            3,
            &[1.; 3],
            LsmrOptions::default(),
            &Execution::default().cancellation(&flag)
        ),
        Err(Error::Cancelled)
    ));
}
#[test]
fn long_thin_filtrations_keep_direct_and_reverse_operators_consistent() {
    for n in [16, 64, 128] {
        let mut cells: Vec<_> = (0..n)
            .map(|v| Simplex::new(vec![v], if v == 0 || v == n - 1 { 0. } else { 1. }).unwrap())
            .collect();
        cells.extend((0..n - 1).map(|v| Simplex::new(vec![v, v + 1], 2.).unwrap()));
        let source = SimplicialComplex::new(cells).unwrap();
        let starts: Vec<_> = [-1., 0., 1., 2.]
            .into_iter()
            .map(|s| source.stage(s).unwrap())
            .collect();
        let target = source.stage(2.).unwrap();
        let steps = fixed_endpoint_laplacians(
            &starts,
            &target,
            0,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        for (k, step) in starts.iter().zip(&steps) {
            let direct = PersistentLaplacian::compute(
                k,
                &target,
                0,
                SpectralOptions::default(),
                &Execution::default(),
            )
            .unwrap();
            close(step.matrix(), direct.matrix());
        }
        // The exact two-terminal Schur matrix has conductance 1/(n-1).
        let c = 1. / (n - 1) as f64;
        close(steps[1].matrix(), &[c, -c, -c, c]);
    }
}

fn pair(source: &SimplicialComplex, s: f64, t: f64, q: usize) -> PersistentLaplacian {
    PersistentLaplacian::compute(
        &source.stage(s).unwrap(),
        &source.stage(t).unwrap(),
        q,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap()
}

#[test]
fn non_branching_path_retains_owned_pair_contracts() {
    let (specialized, direct) = {
        let source = complex(&[
            (&[0], 0.),
            (&[2], 0.),
            (&[1], 1.),
            (&[0, 1], 2.),
            (&[1, 2], 2.),
        ]);
        let k = source.stage(0.).unwrap();
        let l = source.stage(2.).unwrap();
        let result = PersistentLaplacian::compute_non_branching(
            &k,
            &l,
            0,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        let up = PersistentLaplacian::compute_up_non_branching(
            &k,
            &l,
            0,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        close(up.matrix(), result.matrix());
        assert_eq!(up.component(), LaplacianComponent::Up);
        (result, pair(&source, 0., 2., 0))
    };
    close(specialized.matrix(), &[0.5, -0.5, -0.5, 0.5]);
    close(
        &specialized
            .action(&[1., -1.], &Execution::default())
            .unwrap(),
        &[1., -1.],
    );
    close(
        &specialized
            .spectrum(&Execution::default())
            .unwrap()
            .eigenvalues,
        &[0., 1.],
    );
    assert_eq!(specialized.component(), LaplacianComponent::Full);
    assert_eq!(specialized.basis(), direct.basis());
    assert_eq!(specialized.target_basis(), direct.target_basis());
    assert_eq!(specialized.scales(), direct.scales());
    assert_eq!(specialized.source_coverage(), direct.source_coverage());
    assert_eq!(specialized.source_context(), direct.source_context());
    assert_eq!(
        specialized.construction_dimension(),
        direct.construction_dimension()
    );
    assert_eq!(specialized.reference_scale(), direct.reference_scale());
    assert_eq!(specialized.solver_residual(), 0.);
    assert!(!specialized.used_direct_fallback());
}

#[test]
fn non_branching_surface_controls_match_schur_and_hand_spectra() {
    let small = (5. - 5_f64.sqrt()) / 2.;
    let large = (5. + 5_f64.sqrt()) / 2.;
    let examples = [
        (
            vec![vec![0, 1, 2]],
            vec![vec![0, 1], vec![0, 2], vec![1, 2]],
            vec![3.; 3],
        ),
        (
            vec![vec![0, 1, 2], vec![0, 2, 3]],
            vec![vec![0, 1], vec![0, 3], vec![1, 2], vec![2, 3]],
            vec![2., 2., 2., 4.],
        ),
        (
            vec![vec![0, 1, 2], vec![0, 1, 3], vec![0, 2, 3], vec![1, 2, 3]],
            vec![
                vec![0, 1],
                vec![0, 2],
                vec![0, 3],
                vec![1, 2],
                vec![1, 3],
                vec![2, 3],
            ],
            vec![4.; 6],
        ),
        (
            vec![
                vec![1, 3, 4],
                vec![1, 2, 4],
                vec![2, 4, 5],
                vec![2, 3, 5],
                vec![1, 3, 5],
            ],
            vec![vec![3, 4], vec![1, 2], vec![4, 5], vec![2, 3], vec![1, 5]],
            vec![0., small, small, large, large],
        ),
    ];
    for (facets, retained_edges, spectrum) in examples {
        let mut cells = std::collections::BTreeMap::new();
        for facet in facets {
            for mask in 1..(1 << facet.len()) {
                let simplex: Vec<_> = facet
                    .iter()
                    .enumerate()
                    .filter_map(|(i, &v)| (mask & (1 << i) != 0).then_some(v))
                    .collect();
                let birth = if simplex.len() == 1 {
                    0.
                } else if retained_edges.contains(&simplex) {
                    1.
                } else {
                    2.
                };
                cells.insert(simplex, birth);
            }
        }
        let source = SimplicialComplex::new(
            cells
                .into_iter()
                .map(|(v, t)| Simplex::new(v, t).unwrap())
                .collect(),
        )
        .unwrap();
        let k = source.stage(1.).unwrap();
        let l = source.stage(2.).unwrap();
        let full = PersistentLaplacian::compute_non_branching(
            &k,
            &l,
            1,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        let up = PersistentLaplacian::compute_up_non_branching(
            &k,
            &l,
            1,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        let direct_up = PersistentLaplacian::compute_up(
            &k,
            &l,
            1,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        close(full.matrix(), pair(&source, 1., 2., 1).matrix());
        close(up.matrix(), direct_up.matrix());
        close(
            &full.spectrum(&Execution::default()).unwrap().eigenvalues,
            &spectrum,
        );
        let reversed = PersistentLaplacian::compute_non_branching(
            &l,
            &l,
            1,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        close(reversed.matrix(), pair(&source, 2., 2., 1).matrix());
    }
}

#[test]
fn non_branching_eligibility_is_pair_local_and_never_falls_back() {
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 0.),
        (&[3], 1.),
        (&[0, 3], 2.),
        (&[1, 3], 2.),
        (&[2, 3], 2.),
    ]);
    let k = source.stage(0.).unwrap();
    let l = source.stage(2.).unwrap();
    for result in [
        PersistentLaplacian::compute_non_branching(
            &k,
            &l,
            0,
            SpectralOptions::default(),
            &Execution::default(),
        ),
        PersistentLaplacian::compute_up_non_branching(
            &k,
            &l,
            0,
            SpectralOptions::default(),
            &Execution::default(),
        ),
    ] {
        assert!(matches!(
            result,
            Err(Error::InvalidParameter {
                parameter: "non-branching boundary",
                ..
            })
        ));
    }
    assert!(
        PersistentLaplacian::compute_up(
            &k,
            &l,
            0,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_ok()
    );
    // The same full boundary branches in a retained row but D is empty.
    let equal_domain = PersistentLaplacian::compute_up_non_branching(
        &l,
        &l,
        0,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    close(equal_domain.matrix(), pair(&source, 2., 2., 0).matrix());
}

#[test]
fn non_branching_limits_and_stage_certificates_are_preserved() {
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 0.),
        (&[3], 0.),
        (&[0, 1], 1.),
        (&[0, 2], 1.),
        (&[0, 3], 1.),
        (&[1, 2], 1.),
        (&[1, 3], 1.),
        (&[2, 3], 1.),
    ]);
    let stage = source.stage(1.).unwrap();
    // N*N=16 and N*p=24 fit; the p*p identity kernel requires 36 entries.
    let bounded = SpectralOptions {
        max_dense_entries: 24,
        ..Default::default()
    };
    assert!(
        PersistentLaplacian::compute_up(&stage, &stage, 0, bounded, &Execution::default()).is_ok()
    );
    assert!(matches!(
        PersistentLaplacian::compute_up_non_branching(
            &stage,
            &stage,
            0,
            bounded,
            &Execution::default()
        ),
        Err(Error::InvalidParameter {
            parameter: "max_dense_entries",
            ..
        })
    ));
    let flag = AtomicBool::new(true);
    for execution in [
        Execution::default().max_work(0),
        Execution::default().cancellation(&flag),
    ] {
        assert!(
            PersistentLaplacian::compute_non_branching(
                &stage,
                &stage,
                0,
                SpectralOptions::default(),
                &execution
            )
            .is_err()
        );
    }
    let cloned = source.clone();
    assert!(
        PersistentLaplacian::compute_non_branching(
            &stage,
            &cloned.stage(1.).unwrap(),
            0,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    assert!(matches!(
        PersistentLaplacian::compute_non_branching(
            &stage,
            &stage,
            usize::MAX,
            SpectralOptions::default(),
            &Execution::default()
        ),
        Err(Error::SizeOverflow { .. })
    ));
    let points = PointCloudView::new(&[0., 0., 1., 0., 0., 1.], 3, 2).unwrap();
    let certified = RipsBuilder::from_points(points).build_complex(1).unwrap();
    let stage = certified.stage(2.).unwrap();
    assert!(matches!(
        PersistentLaplacian::compute_up_non_branching(
            &stage,
            &stage,
            1,
            SpectralOptions::default(),
            &Execution::default()
        ),
        Err(Error::InsufficientSkeleton { .. })
    ));
    let result = PersistentLaplacian::compute_up_non_branching(
        &stage,
        &stage,
        0,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    assert_eq!(result.construction_dimension(), Some(1));
    assert_eq!(result.source_coverage(), certified.coverage());
    assert!(result.source_context().is_some());
    let bare = certified.complex().stage(2.).unwrap();
    assert!(
        PersistentLaplacian::compute_up_non_branching(
            &bare,
            &stage,
            0,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    for tolerance in [0., f64::NAN, 1.] {
        assert!(
            PersistentLaplacian::compute_up_non_branching(
                &stage,
                &stage,
                0,
                SpectralOptions {
                    relative_tolerance: tolerance,
                    ..Default::default()
                },
                &Execution::default()
            )
            .is_err()
        );
    }
}

#[test]
fn non_branching_zero_and_absent_degrees_are_not_missing_skeletons() {
    let source = complex(&[(&[10], 0.), (&[40], 1.)]);
    for (s, t, q) in [(-1., -1., 0), (-1., 1., 0), (0., 1., 0), (1., 1., 5)] {
        let k = source.stage(s).unwrap();
        let l = source.stage(t).unwrap();
        let result = PersistentLaplacian::compute_non_branching(
            &k,
            &l,
            q,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        close(result.matrix(), pair(&source, s, t, q).matrix());
    }
}

#[test]
fn path_kron_and_owned_action() {
    let result = {
        let source = complex(&[
            (&[0], 0.),
            (&[2], 0.),
            (&[1], 1.),
            (&[0, 1], 2.),
            (&[1, 2], 2.),
        ]);
        pair(&source, 0., 2., 0)
    };
    close(result.matrix(), &[0.5, -0.5, -0.5, 0.5]);
    close(
        &result.spectrum(&Execution::default()).unwrap().eigenvalues,
        &[0., 1.],
    );
    close(
        &result.action(&[1., -1.], &Execution::default()).unwrap(),
        &[1., -1.],
    );
    assert_eq!(result.scales(), (0., 2.));
    assert_eq!(result.basis()[0].vertices(), &[2]);
    assert_eq!(result.target_basis().len(), 3);
    assert_eq!(result.degree(), 0);
    assert!(result.source_context().is_none());
    assert_eq!(
        result.source_coverage(),
        cocycle::filtration::Coverage::Complete
    );
    assert_eq!(result.construction_dimension(), None);
}

#[test]
fn relevant_pairs_combine_buckets_and_extend_the_requested_domain() {
    let results = {
        let source = complex(&[
            (&[0], 0.),
            (&[1], 1.),
            (&[2], 0.),
            (&[3], 0.),
            (&[4], 4.),
            (&[0, 1], 2.),
            (&[2, 3], 3.),
        ]);
        let stages = [0., 1., 2., 3., 4.].map(|s| source.stage(s).unwrap());
        let requests = [(4, 4), (0, 0), (4, 4), (1, 3)];
        let results = relevant_pair_laplacians(
            &stages,
            &requests,
            0,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        for (&(s, t), result) in requests.iter().zip(&results) {
            close(
                result.operator().matrix(),
                pair(&source, stages[s].scale(), stages[t].scale(), 0).matrix(),
            );
        }
        results
    };
    // No single bucket is at this pair: two independent edge events trigger it.
    assert_eq!(results[0].representative(), Some((1, 3)));
    assert_eq!(results[0].upper_embedding(), 0..4);
    assert_eq!(results[0].operator().basis().len(), 5);
    assert_eq!(results[0].operator().target_basis().len(), 5);
    close(&results[0].operator().matrix()[20..], &[0.; 5]);
    assert_eq!(results[0].operator().scales(), (4., 4.));
    assert_eq!(results[1].representative(), None);
    assert!(results[1].upper_embedding().is_empty());
    close(results[1].operator().matrix(), &[0.; 9]);
    close(
        results[0].operator().matrix(),
        results[2].operator().matrix(),
    );
}

#[test]
fn relevant_pairs_match_every_direct_pair_with_ties_cycles_and_lower_changes() {
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 1.),
        (&[3], 1.),
        (&[0, 1], 2.),
        (&[0, 2], 2.),
        (&[1, 2], 2.),
        (&[2, 3], 3.),
        (&[0, 1, 2], 4.),
    ]);
    let stages = [-1., 0., 1., 2., 2., 3., 4.].map(|s| source.stage(s).unwrap());
    let requests: Vec<_> = (0..stages.len())
        .flat_map(|t| (0..=t).map(move |s| (s, t)))
        .collect();
    for q in 0..=3 {
        let results = relevant_pair_laplacians(
            &stages,
            &requests,
            q,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        for (&(s, t), result) in requests.iter().zip(&results) {
            let operator = result.operator();
            close(
                operator.matrix(),
                pair(&source, stages[s].scale(), stages[t].scale(), q).matrix(),
            );
            assert_eq!(operator.scales(), (stages[s].scale(), stages[t].scale()));
            assert_eq!(operator.component(), LaplacianComponent::Full);
            close(
                &operator
                    .spectrum(&Execution::default())
                    .unwrap()
                    .eigenvalues,
                &pair(&source, stages[s].scale(), stages[t].scale(), q)
                    .spectrum(&Execution::default())
                    .unwrap()
                    .eigenvalues,
            );
        }
    }
    let coarse = [source.stage(0.).unwrap(), source.stage(4.).unwrap()];
    for q in 0..=3 {
        let results = relevant_pair_laplacians(
            &coarse,
            &[(0, 0), (0, 1), (1, 1)],
            q,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        for (&(s, t), result) in [(0, 0), (0, 1), (1, 1)].iter().zip(&results) {
            close(
                result.operator().matrix(),
                pair(&source, coarse[s].scale(), coarse[t].scale(), q).matrix(),
            );
        }
    }
}

#[test]
fn relevant_pairs_validate_the_whole_grid_and_bound_storage_and_work() {
    let source = complex(&[(&[0], 0.), (&[1], 1.), (&[0, 1], 2.)]);
    let other = source.clone();
    let stages = [0., 1., 2.].map(|s| source.stage(s).unwrap());
    let options = SpectralOptions::default();
    let run = |grid: &[_], requests: &[_], options, execution: &_| {
        relevant_pair_laplacians(grid, requests, 0, options, execution)
    };
    assert!(run(&stages, &[(2, 1)], options, &Execution::default()).is_err());
    assert!(run(&stages, &[(0, 3)], options, &Execution::default()).is_err());
    assert!(
        run(
            &[stages[1], stages[0], stages[2]],
            &[(0, 2)],
            options,
            &Execution::default()
        )
        .is_err()
    );
    assert!(
        run(
            &[stages[0], other.stage(2.).unwrap()],
            &[(0, 1)],
            options,
            &Execution::default()
        )
        .is_err()
    );
    assert!(
        run(
            &stages,
            &[(0, 2)],
            SpectralOptions {
                max_dense_entries: 5,
                ..options
            },
            &Execution::default()
        )
        .is_err()
    );
    assert!(
        run(
            &stages,
            &[(2, 2); 2],
            SpectralOptions {
                max_dense_entries: 6,
                ..options
            },
            &Execution::default()
        )
        .is_err()
    );
    for limit in [0, 10, 50] {
        assert!(matches!(
            run(
                &stages,
                &[(2, 2)],
                options,
                &Execution::default().max_work(limit)
            ),
            Err(Error::WorkLimitExceeded { .. })
        ));
    }
    let flag = AtomicBool::new(true);
    assert!(matches!(
        run(
            &stages,
            &[(0, 2)],
            options,
            &Execution::default().cancellation(&flag)
        ),
        Err(Error::Cancelled)
    ));
    assert!(
        run(&[], &[], options, &Execution::default())
            .unwrap()
            .is_empty()
    );
    assert!(run(&[], &[(0, 0)], options, &Execution::default()).is_err());
}

#[test]
fn fixed_projection_action_and_bounded_materialization() {
    let operator = {
        let source = complex(&[
            (&[0], 0.),
            (&[2], 0.),
            (&[1], 1.),
            (&[0, 1], 2.),
            (&[1, 2], 2.),
        ]);
        PersistentLaplacianOperator::new(
            &source.stage(0.).unwrap(),
            &source.stage(2.).unwrap(),
            0,
            OperatorOptions::default(),
            &Execution::default(),
        )
        .unwrap()
    };
    close(
        &operator.materialize(4, &Execution::default()).unwrap(),
        &[0.5, -0.5, -0.5, 0.5],
    );
    close(
        &operator.action(&[1., -1.], &Execution::default()).unwrap(),
        &[1., -1.],
    );
    assert_eq!(operator.projection_rank(), 1);
    assert_eq!(operator.projection_entries(), 2);
    assert_eq!(operator.basis()[0].vertices(), &[2]);
    assert_eq!(operator.degree(), 0);
    assert_eq!(operator.scales(), (0., 2.));
    assert!(operator.source_context().is_none());
    assert_eq!(
        operator.source_coverage(),
        cocycle::filtration::Coverage::Complete
    );
    assert_eq!(operator.construction_dimension(), None);
    assert!(matches!(
        operator.materialize(3, &Execution::default()),
        Err(Error::InvalidParameter { .. })
    ));
    assert!(matches!(
        operator.action(&[], &Execution::default()),
        Err(Error::ShapeMismatch { .. })
    ));
    assert!(matches!(
        operator.action(&[f64::NAN, 0.], &Execution::default()),
        Err(Error::NonFiniteValue { .. })
    ));
    assert!(matches!(
        operator.action(&[f64::MAX, f64::MAX], &Execution::default()),
        Err(Error::NumericalFailure { .. })
    ));
    assert!(
        operator
            .action(&[1., -1.], &Execution::new(Some(0), None))
            .is_err()
    );
    let flag = AtomicBool::new(true);
    assert!(matches!(
        operator.action(&[1., -1.], &Execution::new(None, Some(&flag))),
        Err(Error::Cancelled)
    ));
}

#[test]
fn factorized_degenerate_limits_and_inclusion() {
    let empty = complex(&[]);
    let stage = empty.stage(0.).unwrap();
    let op = PersistentLaplacianOperator::new(
        &stage,
        &stage,
        2,
        OperatorOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    assert!(op.action(&[], &Execution::default()).unwrap().is_empty());
    assert!(op.materialize(0, &Execution::default()).unwrap().is_empty());
    let source = complex(&[
        (&[0], 0.),
        (&[2], 0.),
        (&[1], 1.),
        (&[0, 1], 2.),
        (&[1, 2], 2.),
    ]);
    let k = source.stage(0.).unwrap();
    let l = source.stage(2.).unwrap();
    let other = source.clone();
    let make =
        |options| PersistentLaplacianOperator::new(&k, &l, 0, options, &Execution::default());
    for options in [
        OperatorOptions {
            max_incidence_entries: 0,
            ..Default::default()
        },
        OperatorOptions {
            max_projection_entries: 0,
            ..Default::default()
        },
        OperatorOptions {
            max_vector_entries: 1,
            ..Default::default()
        },
    ] {
        assert!(matches!(make(options), Err(Error::InvalidParameter { .. })));
    }
    assert!(
        PersistentLaplacianOperator::new(
            &k,
            &other.stage(2.).unwrap(),
            0,
            OperatorOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    assert!(
        PersistentLaplacianOperator::new(
            &l,
            &k,
            0,
            OperatorOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    assert!(
        PersistentLaplacianOperator::new(
            &k,
            &l,
            usize::MAX,
            OperatorOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    assert!(matches!(
        PersistentLaplacianOperator::new(
            &k,
            &l,
            0,
            OperatorOptions::default(),
            &Execution::new(Some(0), None)
        ),
        Err(Error::WorkLimitExceeded { .. })
    ));
}

#[test]
fn sparse_persistent_action_exceeds_dense_cap_with_nontrivial_projection() {
    let n = 257;
    let mut cells: Vec<_> = (0..n).map(|v| Simplex::new(vec![v], 0.).unwrap()).collect();
    cells.push(Simplex::new(vec![n], 1.).unwrap());
    cells.extend((0..n).map(|v| Simplex::new(vec![v, n], 2.).unwrap()));
    let source = SimplicialComplex::new(cells).unwrap();
    let k = source.stage(0.).unwrap();
    let l = source.stage(2.).unwrap();
    let op = PersistentLaplacianOperator::new(
        &k,
        &l,
        0,
        OperatorOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    assert_eq!(op.options().spectral.max_dense_entries, 65_536);
    assert_eq!(op.projection_rank(), 1);
    assert_eq!(op.projection_entries(), n);
    assert_eq!(op.incidence_entries(), n);
    let x: Vec<_> = (0..n).map(|i| (i % 17) as f64 / 17.).collect();
    let mean = x.iter().sum::<f64>() / n as f64;
    let expected: Vec<_> = x.iter().map(|x| x - mean).collect();
    close(&op.action(&x, &Execution::default()).unwrap(), &expected);
    assert!(matches!(
        op.materialize(usize::MAX, &Execution::default()),
        Err(Error::InvalidParameter { .. })
    ));
    assert!(matches!(
        PersistentLaplacian::compute(&k, &l, 0, SpectralOptions::default(), &Execution::default()),
        Err(Error::InvalidParameter { .. })
    ));
    // Projection storage is an independent cumulative cap, not the dense cap.
    let small = OperatorOptions {
        max_projection_entries: n - 1,
        ..Default::default()
    };
    assert!(matches!(
        PersistentLaplacianOperator::new(&k, &l, 0, small, &Execution::default()),
        Err(Error::InvalidParameter {
            parameter: "max_projection_entries",
            ..
        })
    ));
}

#[test]
fn ring_to_filled_triangle_and_hodge_degeneracy() {
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 0.),
        (&[0, 1], 1.),
        (&[0, 2], 1.),
        (&[1, 2], 1.),
        (&[0, 1, 2], 2.),
    ]);
    let ring = pair(&source, 1., 1., 1);
    close(
        &ring.spectrum(&Execution::default()).unwrap().eigenvalues,
        &[0., 3., 3.],
    );
    let filled = pair(&source, 1., 2., 1);
    close(filled.matrix(), &[3., 0., 0., 0., 3., 0., 0., 0., 3.]);
    close(
        &filled.spectrum(&Execution::default()).unwrap().eigenvalues,
        &[3., 3., 3.],
    );
    close(pair(&source, 2., 2., 1).matrix(), filled.matrix());
    close(pair(&source, 2., 2., 2).matrix(), &[3.]);
    assert!(pair(&source, 2., 2., 3).matrix().is_empty());
}

#[test]
fn singular_removed_component_and_entirely_zero_block() {
    let source = complex(&[(&[0], 0.), (&[1], 1.), (&[2], 1.), (&[1, 2], 2.)]);
    close(pair(&source, 0., 2., 0).matrix(), &[0.]);
    let isolated = complex(&[(&[0], 0.), (&[1], 1.)]);
    close(pair(&isolated, 0., 1., 0).matrix(), &[0.]);
}

#[test]
fn fixed_endpoint_uses_only_upper_term_with_ties_and_empty_stages() {
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 1.),
        (&[3], 1.),
        (&[0, 1], 2.),
        (&[0, 2], 2.),
        (&[1, 2], 2.),
        (&[2, 3], 2.),
        (&[0, 1, 2], 3.),
    ]);
    let stages = [-1., 0., 1., 2., 2., 3.].map(|s| source.stage(s).unwrap());
    for q in 0..=3 {
        let results = fixed_endpoint_laplacians(
            &stages,
            &stages[5],
            q,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        let upper = fixed_endpoint_up_laplacians(
            &stages,
            &stages[5],
            q,
            SpectralOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        for ((k, result), up) in stages.iter().zip(results).zip(upper) {
            close(result.matrix(), pair(&source, k.scale(), 3., q).matrix());
            let direct_up = PersistentLaplacian::compute_up(
                k,
                &stages[5],
                q,
                SpectralOptions::default(),
                &Execution::default(),
            )
            .unwrap();
            close(up.matrix(), direct_up.matrix());
            assert_eq!(up.component(), LaplacianComponent::Up);
            assert_eq!(result.component(), LaplacianComponent::Full);
            assert_eq!(result.scales(), (k.scale(), 3.));
            let eig = result.spectrum(&Execution::default()).unwrap();
            assert!(eig.eigenvalues.iter().all(|v| *v >= -1e-9));
        }
    }
    assert!(
        fixed_endpoint_laplacians(
            &[],
            &stages[5],
            0,
            SpectralOptions::default(),
            &Execution::default()
        )
        .unwrap()
        .is_empty()
    );
}

#[test]
fn basis_permutation_and_orientation_under_relabeling() {
    let a = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 0.),
        (&[0, 1], 1.),
        (&[0, 2], 1.),
        (&[1, 2], 2.),
        (&[0, 1, 2], 3.),
    ]);
    // 0->20, 1->7, 2->40 reverses the orientation of edge 01.
    let b = complex(&[
        (&[20], 0.),
        (&[7], 0.),
        (&[40], 0.),
        (&[7, 20], 1.),
        (&[20, 40], 1.),
        (&[7, 40], 2.),
        (&[7, 20, 40], 3.),
    ]);
    let x = pair(&a, 1., 3., 1);
    let y = pair(&b, 1., 3., 1);
    for (source, expected) in [(&a, &x), (&b, &y)] {
        let operator = PersistentLaplacianOperator::new(
            &source.stage(1.).unwrap(),
            &source.stage(3.).unwrap(),
            1,
            OperatorOptions::default(),
            &Execution::default(),
        )
        .unwrap();
        close(
            &operator.materialize(4, &Execution::default()).unwrap(),
            expected.matrix(),
        );
    }
    close(
        &x.spectrum(&Execution::default()).unwrap().eigenvalues,
        &y.spectrum(&Execution::default()).unwrap().eigenvalues,
    );
    let mapping = [20, 7, 40];
    for (i, s) in x.basis().iter().enumerate() {
        let raw: Vec<_> = s.vertices().iter().map(|v| mapping[*v]).collect();
        let sign_i = if raw[0] < raw[1] { 1. } else { -1. };
        let mut mapped = raw;
        mapped.sort();
        let ii = y
            .basis()
            .iter()
            .position(|v| v.vertices() == mapped)
            .unwrap();
        for (j, t) in x.basis().iter().enumerate() {
            let raw: Vec<_> = t.vertices().iter().map(|v| mapping[*v]).collect();
            let sign_j = if raw[0] < raw[1] { 1. } else { -1. };
            let mut mapped = raw;
            mapped.sort();
            let jj = y
                .basis()
                .iter()
                .position(|v| v.vertices() == mapped)
                .unwrap();
            close(
                &[x.matrix()[i * 2 + j]],
                &[sign_i * sign_j * y.matrix()[ii * 2 + jj]],
            );
        }
    }
}

#[test]
fn invalid_inclusion_policy_and_coordinates() {
    let a = complex(&[(&[0], 0.), (&[1], 1.)]);
    let b = a.clone();
    let k = a.stage(0.).unwrap();
    let l = a.stage(1.).unwrap();
    fn compute(
        k: &cocycle::filtration::SimplicialStage<'_>,
        l: &cocycle::filtration::SimplicialStage<'_>,
        o: SpectralOptions,
        e: &Execution<'_>,
    ) -> Result<PersistentLaplacian> {
        PersistentLaplacian::compute(k, l, 0, o, e)
    }
    assert!(
        compute(
            &k,
            &b.stage(1.).unwrap(),
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    assert!(compute(&l, &k, SpectralOptions::default(), &Execution::default()).is_err());
    assert!(
        fixed_endpoint_laplacians(
            &[l, k],
            &l,
            0,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    for tolerance in [0., -1., f64::NAN, f64::INFINITY, 1.] {
        assert!(
            compute(
                &k,
                &l,
                SpectralOptions {
                    relative_tolerance: tolerance,
                    ..Default::default()
                },
                &Execution::default()
            )
            .is_err()
        );
    }
    let r = pair(&a, 0., 1., 0);
    assert!(matches!(
        r.action(&[], &Execution::default()),
        Err(Error::ShapeMismatch { .. })
    ));
    assert!(matches!(
        r.action(&[f64::NAN], &Execution::default()),
        Err(Error::NonFiniteValue { .. })
    ));
    let source = complex(&[(&[0], 0.), (&[1], 0.), (&[0, 1], 1.)]);
    assert!(matches!(
        pair(&source, 0., 1., 0).action(&[f64::MAX, -f64::MAX], &Execution::default()),
        Err(Error::NumericalFailure { .. })
    ));
}

#[test]
fn cancellation_work_and_storage_limits() {
    let a = complex(&[(&[0], 0.), (&[1], 1.), (&[0, 1], 2.)]);
    let k = a.stage(0.).unwrap();
    let l = a.stage(2.).unwrap();
    let flag = AtomicBool::new(true);
    for execution in [
        Execution::default().cancellation(&flag),
        Execution::default().max_work(0),
        Execution::default().max_work(8),
    ] {
        assert!(
            PersistentLaplacian::compute(&k, &l, 0, SpectralOptions::default(), &execution)
                .is_err()
        );
        assert!(
            fixed_endpoint_laplacians(&[k, l], &l, 0, SpectralOptions::default(), &execution)
                .is_err()
        );
    }
    let small = SpectralOptions {
        max_dense_entries: 3,
        ..Default::default()
    };
    assert!(PersistentLaplacian::compute(&k, &l, 0, small, &Execution::default()).is_err());
    let o = SpectralOptions {
        max_dense_entries: 4,
        ..Default::default()
    };
    assert!(fixed_endpoint_laplacians(&[l, l], &l, 0, o, &Execution::default()).is_err());
    let r = pair(&a, 0., 2., 0);
    assert!(matches!(
        r.spectrum(&Execution::default().cancellation(&flag)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        r.action(&[1.], &Execution::default().max_work(0)),
        Err(Error::WorkLimitExceeded { .. })
    ));
}

#[test]
fn certified_missing_skeleton_is_not_an_empty_upper_boundary() -> Result<()> {
    let p = PointCloudView::new(&[0., 0., 1., 0., 0., 1.], 3, 2)?;
    let source = RipsBuilder::from_points(p).build_complex(1)?;
    let k = source.stage(2.)?;
    assert!(matches!(
        PersistentLaplacianOperator::new(
            &k,
            &k,
            1,
            OperatorOptions::default(),
            &Execution::default()
        ),
        Err(Error::InsufficientSkeleton { .. })
    ));
    assert!(matches!(
        PersistentLaplacian::compute(&k, &k, 1, SpectralOptions::default(), &Execution::default()),
        Err(Error::InsufficientSkeleton { .. })
    ));
    assert!(matches!(
        PersistentLaplacian::compute_up(
            &k,
            &k,
            1,
            SpectralOptions::default(),
            &Execution::default()
        ),
        Err(Error::InsufficientSkeleton { .. })
    ));
    let lower = PersistentLaplacian::compute_down(
        &k,
        1,
        SpectralOptions::default(),
        &Execution::default(),
    )?;
    assert_eq!(lower.component(), LaplacianComponent::Down);
    assert_eq!(lower.construction_dimension(), Some(1));
    assert!(lower.source_context().is_some());
    assert!(matches!(
        PersistentLaplacian::compute_down(&k, 2, SpectralOptions::default(), &Execution::default()),
        Err(Error::InsufficientSkeleton { .. })
    ));
    let bare = source.complex().stage(2.)?;
    assert!(
        PersistentLaplacian::compute(
            &bare,
            &k,
            0,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_err()
    );
    assert!(
        PersistentLaplacian::compute(
            &bare,
            &bare,
            1,
            SpectralOptions::default(),
            &Execution::default()
        )
        .is_ok()
    );
    let result =
        PersistentLaplacian::compute(&k, &k, 0, SpectralOptions::default(), &Execution::default())?;
    assert!(result.source_context().is_some());
    assert_eq!(result.construction_dimension(), Some(1));
    Ok(())
}

#[test]
fn component_matrices_and_nonzero_spectral_multisets() {
    // A filled triangle and a separate unfilled ring. In C1(K), up has rank 1,
    // down has rank 4, and the full operator has one harmonic zero direction.
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 0.),
        (&[3], 0.),
        (&[4], 0.),
        (&[5], 0.),
        (&[0, 1], 1.),
        (&[0, 2], 1.),
        (&[1, 2], 1.),
        (&[3, 4], 1.),
        (&[3, 5], 1.),
        (&[4, 5], 1.),
        (&[0, 1, 2], 2.),
    ]);
    let k = source.stage(1.).unwrap();
    let l = source.stage(2.).unwrap();
    let options = SpectralOptions::default();
    let execution = Execution::default();
    let up = PersistentLaplacian::compute_up(&k, &l, 1, options, &execution).unwrap();
    let down = PersistentLaplacian::compute_down(&k, 1, options, &execution).unwrap();
    let full = PersistentLaplacian::compute(&k, &l, 1, options, &execution).unwrap();
    assert_eq!(up.component(), LaplacianComponent::Up);
    assert_eq!(down.component(), LaplacianComponent::Down);
    assert_eq!(full.component(), LaplacianComponent::Full);
    assert_eq!(up.scales(), (1., 2.));
    assert_eq!(down.scales(), (1., 1.));
    assert_eq!(up.basis(), down.basis());
    assert_eq!(up.basis(), full.basis());
    let n = full.basis().len();
    for i in 0..n {
        for j in 0..n {
            close(
                &[up.matrix()[i * n + j] + down.matrix()[i * n + j]],
                &[full.matrix()[i * n + j]],
            );
            let product = (0..n)
                .map(|a| up.matrix()[i * n + a] * down.matrix()[a * n + j])
                .sum();
            close(&[product], &[0.]);
        }
    }
    let up_spec = up.spectrum(&execution).unwrap();
    let down_spec = down.spectrum(&execution).unwrap();
    assert_eq!(up_spec.component, LaplacianComponent::Up);
    assert_eq!(down_spec.component, LaplacianComponent::Down);
    close(&up_spec.eigenvalues, &[0., 0., 0., 0., 0., 3.]);
    close(&down_spec.eigenvalues, &[0., 0., 3., 3., 3., 3.]);
    close(
        &full.spectrum(&execution).unwrap().eigenvalues,
        &[0., 3., 3., 3., 3., 3.],
    );
    let x = [1., -2., 3., 4., -5., 6.];
    let u = up.action(&x, &execution).unwrap();
    let d = down.action(&x, &execution).unwrap();
    let sum: Vec<_> = u.iter().zip(&d).map(|(a, b)| a + b).collect();
    close(&sum, &full.action(&x, &execution).unwrap());
    let empty = source.stage(-1.).unwrap();
    assert!(
        PersistentLaplacian::compute_up(&empty, &l, 1, options, &execution)
            .unwrap()
            .matrix()
            .is_empty()
    );
    assert!(
        PersistentLaplacian::compute_down(&empty, 1, options, &execution)
            .unwrap()
            .spectrum(&execution)
            .unwrap()
            .eigenvalues
            .is_empty()
    );
    close(
        PersistentLaplacian::compute_down(&k, 0, options, &execution)
            .unwrap()
            .matrix(),
        &[0.; 36],
    );
}

#[test]
fn lower_component_skips_target_storage_and_schur_work() {
    let source = complex(&[
        (&[0], 0.),
        (&[1], 0.),
        (&[2], 2.),
        (&[0, 1], 1.),
        (&[0, 2], 3.),
        (&[1, 2], 3.),
        (&[0, 1, 2], 4.),
    ]);
    let k = source.stage(1.).unwrap();
    let l = source.stage(4.).unwrap();
    let options = SpectralOptions {
        max_sweeps: 0,
        max_dense_entries: 1,
        ..Default::default()
    };
    let execution = Execution::default();
    let down = PersistentLaplacian::compute_down(&k, 1, options, &execution).unwrap();
    close(down.matrix(), &[2.]);
    assert_eq!(down.target_basis(), down.basis());
    assert_eq!(down.solver_residual(), 0.);
    assert!(PersistentLaplacian::compute_up(&k, &l, 1, options, &execution).is_err());
    let zero_sweeps = SpectralOptions {
        max_sweeps: 0,
        ..Default::default()
    };
    assert!(matches!(
        PersistentLaplacian::compute_up(&k, &l, 1, zero_sweeps, &execution),
        Err(Error::NumericalFailure { .. })
    ));
    let flag = AtomicBool::new(true);
    assert!(matches!(
        PersistentLaplacian::compute_down(&k, 1, options, &execution.cancellation(&flag)),
        Err(Error::Cancelled)
    ));
    assert!(matches!(
        PersistentLaplacian::compute_down(&k, 1, options, &execution.max_work(0)),
        Err(Error::WorkLimitExceeded { .. })
    ));
    let insufficient = SpectralOptions {
        max_dense_entries: 0,
        ..Default::default()
    };
    assert!(PersistentLaplacian::compute_down(&k, 1, insufficient, &execution).is_err());
}

#[test]
fn upper_zero_cancellation_keeps_the_construction_scale() {
    // K is the two-edge path 1--2--3, so its cycle space is zero and U=0.
    // Three later triangles give a nontrivial Schur solve. Its tiny signed
    // cancellation error must be judged against the original upper Gram.
    let source = complex(&[
        (&[0], 1.),
        (&[1], 0.),
        (&[2], 0.),
        (&[3], 0.),
        (&[0, 1], 2.),
        (&[0, 2], 1.),
        (&[0, 3], 2.),
        (&[1, 2], 0.),
        (&[1, 3], 2.),
        (&[2, 3], 0.),
        (&[0, 1, 2], 3.),
        (&[0, 1, 3], 2.),
        (&[1, 2, 3], 3.),
    ]);
    let k = source.stage(0.).unwrap();
    let l = source.stage(3.).unwrap();
    let up = PersistentLaplacian::compute_up(
        &k,
        &l,
        1,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    close(up.matrix(), &[0.; 4]);
    assert_eq!(up.reference_scale(), 2.);
    let spectrum = up.spectrum(&Execution::default()).unwrap();
    assert_eq!(
        spectrum.zero_threshold,
        up.options().relative_tolerance * up.reference_scale()
    );
    assert!(
        spectrum
            .eigenvalues
            .iter()
            .all(|value| value.abs() <= spectrum.zero_threshold)
    );
}

#[test]
fn default_dense_limit_rejects_large_input_before_gram_allocation() {
    let source = SimplicialComplex::new(
        (0..257)
            .map(|v| Simplex::new(vec![v], 0.).unwrap())
            .collect(),
    )
    .unwrap();
    let k = source.stage(0.).unwrap();
    assert!(matches!(
        PersistentLaplacian::compute(&k, &k, 0, SpectralOptions::default(), &Execution::default()),
        Err(Error::InvalidParameter {
            parameter: "max_dense_entries",
            ..
        })
    ));
}

#[test]
fn real_rp2_homology_does_not_inherit_f2_torsion() {
    let facets = [
        [0, 1, 2],
        [0, 1, 3],
        [0, 2, 4],
        [0, 3, 5],
        [0, 4, 5],
        [1, 2, 5],
        [1, 3, 4],
        [1, 4, 5],
        [2, 3, 4],
        [2, 3, 5],
    ];
    let mut faces = std::collections::BTreeSet::new();
    for facet in facets {
        for mask in 1..8 {
            faces.insert(
                facet
                    .iter()
                    .enumerate()
                    .filter(|(i, _)| mask & (1 << i) != 0)
                    .map(|(_, v)| *v)
                    .collect::<Vec<_>>(),
            );
        }
    }
    let source = SimplicialComplex::new(
        faces
            .into_iter()
            .map(|v| {
                let value = (v.len() - 1) as f64;
                Simplex::new(v, value).unwrap()
            })
            .collect(),
    )
    .unwrap();
    // H1(RP2;R)=0, whereas H1(RP2;F2)=F2. The oracle tool verifies
    // integer/Q ranks and separately exhibits GUDHI's two field counts.
    assert!(
        pair(&source, 2., 2., 1)
            .spectrum(&Execution::default())
            .unwrap()
            .eigenvalues
            .iter()
            .all(|v| *v > 1e-8)
    );
}

#[test]
fn long_reverse_elimination_matches_hand_derived_path_and_cone() {
    // Eliminating 254 interior vertices gives conductance 1/255 between
    // the retained endpoints; this stresses many successive positive pivots.
    let mut cells: Vec<_> = (0..256)
        .map(|i| Simplex::new(vec![i], if i == 0 || i == 255 { 0.0 } else { 1.0 }).unwrap())
        .collect();
    cells.extend((0..255).map(|i| Simplex::new(vec![i, i + 1], 2.0).unwrap()));
    let path = SimplicialComplex::new(cells).unwrap();
    let k = path.stage(0.0).unwrap();
    let l = path.stage(2.0).unwrap();
    let expected = [1.0 / 255.0, -1.0 / 255.0, -1.0 / 255.0, 1.0 / 255.0];
    close(pair(&path, 0.0, 2.0, 0).matrix(), &expected);
    let reverse = fixed_endpoint_laplacians(
        &[k],
        &l,
        0,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    close(reverse[0].matrix(), &expected);

    // Cone over a path: 255 edges, 127 faces. The removed upper block is
    // singular. K contains two spokes forming a tree, so its persistent upper
    // term is exactly zero and the full operator is its lower Gram.
    let mut cells: Vec<_> = (0..129)
        .map(|v| Simplex::new(vec![v], 0.0).unwrap())
        .collect();
    cells
        .extend((1..129).map(|v| {
            Simplex::new(vec![0, v], if v == 1 || v == 128 { 0.0 } else { 1.0 }).unwrap()
        }));
    cells.extend((1..128).map(|v| Simplex::new(vec![v, v + 1], 1.0).unwrap()));
    cells.extend((1..128).map(|v| Simplex::new(vec![0, v, v + 1], 2.0).unwrap()));
    let cone = SimplicialComplex::new(cells).unwrap();
    let k = cone.stage(0.0).unwrap();
    let l = cone.stage(2.0).unwrap();
    let expected = [2.0, 1.0, 1.0, 2.0];
    close(pair(&cone, 0.0, 2.0, 1).matrix(), &expected);
    let reverse = fixed_endpoint_laplacians(
        &[k],
        &l,
        1,
        SpectralOptions::default(),
        &Execution::default(),
    )
    .unwrap();
    close(reverse[0].matrix(), &expected);
}
