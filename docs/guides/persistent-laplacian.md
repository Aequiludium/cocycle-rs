# Real persistent Laplacians

[Documentation](../README.md) / Usage

This **Unreleased** bounded dense implementation computes an actual stage-pair
operator, its action and its complete numerical spectrum. It uses real
coefficients, standard inner products, increasing-vertex orientation and ordinary
nonaugmented chains. It requires a checked same-source inclusion; cloned or
unrelated owners, reversed scales and mixed certificates are rejected.

```rust
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::spectral::{
    LaplacianComponent, PersistentLaplacian, SpectralOptions, fixed_endpoint_laplacians,
};

let source = SimplicialComplex::new(vec![
    Simplex::new(vec![0], 0.)?, Simplex::new(vec![2], 0.)?,
    Simplex::new(vec![1], 1.)?, Simplex::new(vec![0, 1], 2.)?,
    Simplex::new(vec![1, 2], 2.)?,
])?;
let k = source.stage(0.)?;
let l = source.stage(2.)?;
let options = SpectralOptions::default();
let execution = Execution::default();
let operator = PersistentLaplacian::compute(&k, &l, 0, options, &execution)?;
assert_eq!(operator.matrix(), &[0.5, -0.5, -0.5, 0.5]);
assert_eq!(operator.action(&[1., -1.], &execution)?, [1., -1.]);
let spectrum = operator.spectrum(&execution)?;
assert!(spectrum.eigenvalues[0].abs() < 1e-10);
assert!((spectrum.eigenvalues[1] - 1.).abs() < 1e-10);
let up = PersistentLaplacian::compute_up(&k, &l, 0, options, &execution)?;
let down = PersistentLaplacian::compute_down(&k, 0, options, &execution)?;
assert_eq!(up.matrix(), operator.matrix()); // ordinary B0 is zero
assert_eq!(down.matrix(), &[0.; 4]);
assert_eq!(down.scales(), (0., 0.)); // no endpoint L is required
assert_eq!(up.component(), LaplacianComponent::Up);
assert_eq!(down.spectrum(&execution)?.component, LaplacianComponent::Down);
let starts = [source.stage(-1.)?, k, l];
let reverse = fixed_endpoint_laplacians(&starts, &l, 0, options, &execution)?;
assert_eq!(reverse[1].matrix(), operator.matrix());
# Ok::<(), cocycle::Error>(())
```

Run `cargo run --locked --example persistent_laplacian`. Matrix rows and action
coordinates follow `basis()`, whose simplices retain original labels and births.
The owned result also carries the target q-basis, degree, pair scales, source
context, source coverage and construction dimension. It survives source release.
Metadata describes the source and requested pair, not cross-source equality.
No persistence diagram conversion is implicit.

`compute_up` and `compute_down` return the same owned result type with an
explicit `LaplacianComponent` tag. Each has its own matrix, action and spectrum.
The lower constructor accepts K alone, records (s,s) and K as both bases, and
requires only the certified q-skeleton. It does not form a target upper Gram or
run Schur elimination. This allows the lower result/spectrum to be reused across
different endpoints. Full and upper constructors still require the certified
q+1 skeleton. Spectrum component tags distinguish Full/Up/Down, but do not
authenticate a pair: consumers must retain the operator's degree, ordered bases
and source metadata. Spectra of equally sized or numerically equal matrices
from different requests are not implicitly identified.

In this finite-dimensional standard-inner-product setting, write U for the
persistent upper term and D for the lower term. The chain identity gives
UD=DU=0. The nonzero spectrum of U+D is the multiset union of the two nonzero
component spectra; the full zero multiplicity is z(U)+z(D)-dim(C_q(K)).
Concatenating both complete spectra would incorrectly double the dimension.
See Wolf--Fan--Monod, [arXiv:2509.20220v1](https://arxiv.org/html/2509.20220v1#S5.SS3),
Theorem 5.11. The current `spectrum` implementation diagonalizes the selected
matrix and never infers exact ranks or merges near-zero eigenvalues implicitly.
Only the full operator's exact nullity is the persistent Betti number; component
nullities need not be. This API does not assert arbitrary eigenvalue arrays are
monotone or stable: the paper's spectral-function theorems have their stated
filtration, inclusion and comparison conditions. In particular, fixing K and
enlarging L keeps D fixed, so it does not lose the upper term's monotonicity.

The direct route follows the author's [SchurComp.m](https://github.com/ndag/Persistent-Laplacian/blob/ebff6e3ad6f5eb90eb8058e4f5567755b62f0c3e/SchurComp.m): solve C X = E^T in minimum
norm, then form A - E X without materializing C+. A cyclic symmetric Jacobi
eigensystem supplies the retained solve directions. This order avoids squaring
couplings before division, which can lose a representable correction at extreme
scales. Relative eigenvalue cutoff, eigensystem, discarded-range, actual solve
and reconstructed symmetry residuals are checked. The complete spectrum retains multiplicity and small
signed roundoff values. Neither a thresholded pseudoinverse nor a near-zero count
is an exact Betti certificate. `solver_residual()` describes the Schur solve;
`LaplacianSpectrum::residual` describes the final eigensystem. These are max-entry
residuals normalized by their respective scales, not forward-error or
condition-number certificates. Solve error uses max(abs(E)); a zero denominator
keeps the absolute error. Reconstruction symmetry uses the original upper Gram
scale. See the exact policy in the source rustdoc.

`reference_scale()` is the larger max-entry scale of the assembled matrix and
its original target upper Gram; lower-only results use their own Gram scale.
PSD acceptance uses `residual_tolerance` times this reference, and the suggested
near-zero threshold uses `relative_tolerance` times it. Thus an exactly zero
upper term with tiny signed Schur cancellation errors is not judged solely at
the scale of that rounding noise. Eigenvalues remain unclipped. Eigensystem
residuals are still normalized by the selected matrix, and neither these
thresholds nor the reference scale certify forward error or exact ranks.

`fixed_endpoint_laplacians` builds the target upper term once and removes
q-coordinates in reverse filtration order. It adds each start's own lower term,
including ties, repeated scales and empty stages, and returns input order.
Near-zero pivots require a small eliminated row; incompatible rows fail explicitly.
The algorithm is caller-selected. Performance evidence does not automatically
adopt it as the default or establish online endpoint-update support.
`fixed_endpoint_up_laplacians` performs the same elimination while returning
only the tagged upper results, without constructing any source lower terms.

The default cap is 65,536 f64 entries per dense buffer, limiting the target
q-basis to 256 entries. Fixed-endpoint results additionally cap their cumulative
matrix entries. Buffers coexist; the cap is not a total RAM quota. Storage is
quadratic. Each Jacobi sweep on an m-by-m matrix costs O(m^3); the configured
`max_sweeps` bounds the number of sweeps, and residual evaluation also requires
matrix products. Schur construction, incidence traversal and owned output have
additional costs. A paper's construction bound is not an end-to-end bound for
this iterative implementation or its complete spectrum.
Large inputs are rejected before Gram allocation. This is a core reproduction,
not a scalable sparse or matrix-free eigensolver. `Execution` checks cooperative
cancellation and work limits throughout incidence, dense rows and rotations;
sorting and allocation cannot be interrupted internally.

For certified constructions, a missing q+1 skeleton is an error unless expansion
exhausted the source. Missing cofaces are not silently read as zeros. Bare supplied
complexes instead define the entire mathematical source, with genuinely absent
adjacent degrees treated as zero. Weighted Grams, complex coefficients, general
chain maps, multiparameter spectra and dynamic insertion/deletion have no entry
point in this API. Multiscale consumers and spectral feature learning remain
separate integration work.

See the [mathematical specification](../reference/mathematics.md#22-real-persistent-laplacians)
and [independent reproduction tool](https://github.com/Aequiludium/cocycle-rs/blob/main/tools/compare_persistent_laplacians.py).
The tool uses NumPy Schur and SVD nullspace routes, exact rational restricted-domain
matrices and real persistent ranks. Its optional GUDHI checks use explicitly
torsion-free fixtures, since a finite-field diagram is not general real truth.
