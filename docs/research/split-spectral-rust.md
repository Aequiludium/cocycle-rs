# Bounded Rust parity for split spectra and harmonic restriction

[Documentation](../index.md) / Research

This is Round 2 of [Issue #145](https://github.com/Aequiludium/cocycle-rs/issues/145).
The independent Python reproduction remains in unmerged
[PR #151](https://github.com/Aequiludium/cocycle-rs/pull/151). Its source is pinned
to commit `8ee8f49b1aed55816456a0ef46eeec2501151235` and SHA-256
`9d33727de04b19db2d595c75f420e69a4a0c0457d7e0574fd86ea4d02c9f3b62`.
The second-round branch starts separately from main
`a3c3e1132800624743ad0f159c9303090e0b3b3f`; it does not merge, replace or
change the first-round PR. No PETLS author source is copied.

## Implemented scope

The callable Rust path is in `persistent_laplacian`, with its
[usage and numeric contract](../guides/persistent-laplacian.md#bounded-gram-and-harmonic-restriction-operations).
It constructs both Grams, restores numerical zero multiplicity in the original
domain, lifts dual positive vectors, checks the separate up/down spectral-union
extension, and restricts using a complete supplied harmonic kernel. Kernel
normalization uses one-sided Jacobi singular values; comparison is with the
independent first-round SciPy SVD/QR output. Bases may rotate, so tests compare
lifted matrices and subspace projectors rather than complement coordinates.

The basic LOBPCG transcription retains Algorithms 1/2, including Cholesky of
the scaled search metric and the coefficient update for P. Initial vectors
are projected from original coordinates. The existing Hodge/Dirac solver entry
is unchanged. The separate first-round SVD search correction and robust
Algorithm 8 are **not** included in this Rust scope. The clustered paper-basic
case must retain its ill-conditioned-search failure; a passing corrected Python
run is not evidence that the Rust basic run succeeded.

The implementation is deliberately unoptimized: both Grams and full dense
decompositions remain. There is no cache or automatic dispatch. It does not
remove the existing target-domain Schur allocation or establish a scalable
positive solver. Existing complete `spectrum()` output is preserved.

## Freeze and replay

`tools/compare_split_spectra.py` verifies the first-round source hash before
importing it. `--frozen DIRECTORY` consumes an existing first-round directory,
checking its source identity and all recorded artifact hashes before copying it
into the new run. The inputs are passed unchanged to the Rust worker. This is
the local historical replay route; no expectation comes from Rust.

Without `--frozen`, the fixed first-round program first generates and freezes a
fresh platform-specific replay, including its GUDHI checks. Both sides then
consume the same saved inputs. Windows and Linux BLAS can produce different
last bits in regenerated seeded rotations; historical replay and fresh-platform
replay remain distinct evidence with recorded hashes.

```sh
python3 -m pip install numpy==2.2.6 scipy==1.15.3 sympy==1.14.0 gudhi==3.12.0
git clone https://github.com/Aequiludium/cocycle-rs.git target/split-round1-reference
git -C target/split-round1-reference checkout --detach 8ee8f49b1aed55816456a0ef46eeec2501151235
cargo build --locked --release --example split_spectra
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 tools/compare_split_spectra.py --worker target/release/examples/split_spectra --reference target/split-round1-reference/tools/reproduce_split_spectra.py --output target/split-rust-parity-001
```

Use `split_spectra.exe` on Windows. Append `--frozen` with the preserved first-round
directory for historical replay. Every output directory must be new. The CLI
retains the exact worker requests, stdout/stderr, nonzero failures, inputs,
full spectra, lifted matrices, reference traces, hashes and source identities.
It records stage timers as diagnostics, not comparable speed measurements.

## Validation scope and limits

The 29 comparison cases cover rectangular, empty and zero boundaries, triangle
and tetrahedron-shell spectra, both frozen pair controls, real harmonic
restriction, the nonharmonic-representative counterexample, large nullity, tiny
positive ambiguity, basic separated/repeated spectra, clustered-search failure,
rank/count/iteration errors and execution/storage controls. Raw complete values
are compared before interpreting zeros. Numerical zero signs can differ between
LAPACK and cyclic Jacobi, including advisory ambiguity flags on roundoff zeros;
positive cutoff controls, ranks, multiplicities and unfiltered spectra remain
explicit. No exact rank certificate is inferred from either backend.

Matrix/projector errors use Frobenius normalization `max(1,||reference||F)`
with tolerance 1e-10. Iterative value comparisons use 1e-8, plus independently
checked original-coordinate normalized residuals <=1e-9, kernel leakage and
orthogonality. Iteration counts are recorded, not asserted equal. GUDHI only
corroborates the declared torsion-free topology; it cannot validate real spectral
values, arbitrary real ranks or eigenvectors. SciPy dense/iterative and exact
rational/analytic controls provide those bounded independent checks.

Public contract tests cover ownership, shapes, nonfinite values, negative or
incomplete domains, kernel normalization, zero-count validation and resource
recovery. Private tests independently check a hand-derived generalized
two-by-two characteristic polynomial and cancellation during preparation and
iteration. The CI job pins the first-round source in a separate checkout and
uploads raw results for 14 days. Inspect its actual head and outcome; configured
checks are not completed verification. Rust parity does not close the separate
optimization round or resolve the original missing-Boost PETLS author build.
