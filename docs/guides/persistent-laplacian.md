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
use cocycle::spectral::{PersistentLaplacian, SpectralOptions, fixed_endpoint_laplacians};

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

`fixed_endpoint_laplacians` builds the target upper term once and removes
q-coordinates in reverse filtration order. It adds each start's own lower term,
including ties, repeated scales and empty stages, and returns input order.
Near-zero pivots require a small eliminated row; incompatible rows fail explicitly.
The algorithm is caller-selected. Performance evidence does not automatically
adopt it as the default or establish online endpoint-update support.

The default cap is 65,536 f64 entries per dense buffer, limiting the target
q-basis to 256 entries. Fixed-endpoint results additionally cap their cumulative
matrix entries. Buffers coexist; the cap is not a total RAM quota. Storage is
quadratic and eigensystem/Schur work is cubic in the relevant basis sizes.
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
