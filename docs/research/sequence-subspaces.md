# Controlled sequence-subspace reproduction

[Documentation](../README.md) / Research

This is the bounded first implementation for
[Issue 148](https://github.com/Aequiludium/cocycle-rs/issues/148). The optional
repository tools reproduce a **corrected, fixed-degree, real-symmetric**
SCSF/ChFSI control and replay its filter inputs in Rust. The private
[LOBPCG owner](../../src/eigensolver.rs) accepts a supplied starting block;
existing public callers still use its unchanged cold-start wrapper. There is
no public FFT, ChFSI, persistent warm-start API or change to spectral coverage.

## Frozen sources and discrepancies

| Source | Inspected scope | PDF SHA-256 |
| --- | --- | --- |
| Wang et al., [arXiv:2510.23215v2](https://arxiv.org/abs/2510.23215v2), 2026-04-27 | Algorithms 1-3, physical PDF pages 3-6; equations (2), (5)-(6) | `eee5de671463c0a98ef7887696cc8ae5051fe4c4798707c88d50aedb22d44617` |
| Berljafa, Wortmann and Di Napoli, [arXiv:1404.4161v2](https://arxiv.org/abs/1404.4161v2), 2014-07-06 | Section 3.1, Algorithms 1-2; physical PDF pages 9-12 | `464a08f4f812986d3eb40ad99028a4e083de39ceef559199883be781a7e845e2` |

No author implementation is pinned or claimed equivalent. These are independent
mathematical implementations, not copied author source or original datasets.
The companion paper has the same printed current-iterate subtraction as SCSF;
it cannot independently resolve that discrepancy by transcription alone.

For a declared suppression interval `[alpha,beta]`, let
`c=(alpha+beta)/2`, `e=(beta-alpha)/2`, `H=(A-cI)/e` and
`t0=(anchor-c)/e`, with `anchor < alpha < beta`. The defining control is
`Y_j=T_j(H)Y_0/T_j(t0)`. Therefore:

```text
sigma_1 = 1/t0
Y_1 = sigma_1 H Y_0
sigma_(j+1) = 1/(2/sigma_1 - sigma_j)
Y_(j+1) = 2 sigma_(j+1) H Y_j - sigma_(j+1) sigma_j Y_(j-1)
```

Degree zero returns `Y_0`. The corrected filter supports a separate degree for
each column without confusing degree and column indexes. The retained printed
control subtracts `Y_j` instead of `Y_(j-1)`; it is explicitly a polynomial
diagnostic with the ambiguous column split removed, **not an executable faithful
replay of the complete printed algorithm**. Independent NumPy scalar Chebyshev
evaluation checks degrees 0, 1, 2, 3 and 8, plus mixed column degrees.

Other explicit interpretations/corrections are:

- SCSF Algorithm 1 overloads degree and column-split notation; the companion
  Algorithm 2 also uses a row-size symbol in a column slice. The control tracks
  each column's requested degree explicitly.
- SCSF Algorithm 2's finite distance sentinel can leave its next index undefined
  when every distance exceeds 1000. The implementation minimizes over all unused
  descriptors. Start at original index 0; equal distances choose the lowest
  original index. Use the top-left `p0`-by-`p0` block of an unshifted, unnormalized
  FFT2. This fixes a truncation convention rather than assuming it is unique.
- With `l` locked and `a` active vectors, concatenate an `n`-by-`(l+a)` block for
  QR, retain its last `a` columns as `Q_active`, form an `a`-by-`a` quotient
  `Q_active.T A Q_active`, and rotate **Q_active** by the `a`-by-`a` Ritz vectors.
  The printed SCSF Algorithm 3 mixes `Q`, `Q_0` and the pre-QR vectors. Every
  control iteration saves these shapes and actual blocks.
- A rank-deficient filtered block fails explicitly. Converged leading pairs are
  locked only after fresh current-operator residuals pass. The retained block
  includes lookahead; only the previous block is transferred to the next sample.

## Controlled numerical contract

The sequence controls have 24 coordinates, three requested eigenpairs and six
block columns. A known spectrum contains a repeated lowest eigenvalue; an
independently constructed orthogonal rotation provides its exact projector.
Eight parameter fields give a shuffled smooth sequence. A second family uses
independent abrupt rotations with the same descriptors as an adverse control.
All four routes use the same matrices and accuracy: dense `eigh`, cold ChFSI,
original-order warm ChFSI and FFT-ordered SCSF. Outputs are restored to original
sample order. Setup, ordering, solve and verification diagnostics remain separate.

The filter degree is fixed at 8; `[2,12.1]` and anchor `0.3` are declared from
these known controls. The fixed-degree iteration does **not** implement the
companion paper's Lanczos bound estimator or residual-based degree optimizer.
Nor does it resolve general complex Hermitian problems, interior or
indefinite smallest-magnitude targeting, MPI/Elemental, or a full paper dataset.
The selected Rust rewrite covers the verified filter and the separate existing
LOBPCG initialization, not the complete FFT/QR/ChFSI pipeline.

Acceptance uses eigenvalue error at most `1e-10`, normalized original residual
`||AV-VLambda||F/max(1,||A||F)` at most `1e-9`, projector error at most `1e-8`,
and orthogonality error at most `1e-10`. The iteration's individual residual
stop is `1e-11*max(1,abs(lambda))`. These are numerical controls, not exact
rank/index certificates; repeated eigenvectors are compared as projectors.

The equation (2) metric control uses the printed 2-by-2 interior Dirichlet stencil,
unit mesh spacing and positive field `[1,2,3,4]`. For `A=D L`, use
`S=sqrt(D) L sqrt(D)` and restore `u=sqrt(D)v`. Verify residuals on **A** and
orthogonality in `D^-1`; `A` itself is not Euclidean symmetric. This separate
model check preserves ordering by eigenvalue magnitude. It does not turn the
positive algebraically smallest sequence controls into a full PDE dataset replay.

## Supplied blocks and changing bases

The private helper takes finite row-major `n`-by-`min(block_size,n)` input under
the existing entry cap and shared `WorkBudget`. It reuses two-pass
orthogonalization, drops relative directions below `1e-13`, completes rank with
the deterministic cold sequence and coordinate directions, and recomputes Ritz
values and final residuals using the **new** action. An all-zero input replays
the cold span exactly. No old values, residuals, global cache or history are used.
An invariant supplied block can miss smaller eigenvalues even at zero residual;
the existing uncertified-index semantics still apply.

The separate tiny basis-transfer experiment uses the existing independent
exact-rational persistent-matrix oracle for path, triangle and tetrahedron
filtrations. It checks fixed-start and fixed-target pairs in degrees 0-2,
including birth, filling, empty stages and rank changes. Transfer uses exact
source-local canonical simplex tuples and orientation signs, rejects different
sources/degrees/metrics, reorthogonalizes and fills missing directions under an
entry cap. Tests also check reordered/reversed bases, rank loss and invalid
metadata. A newly computed dense positive-range projector checks removal of new
harmonic directions separately; insufficient projected rank fails explicitly.

These are full-spectrum tiny Rust LOBPCG checks with a separate dense
positive-projection control. They are **not** a scalable positive-spectrum solver,
a certified source-identity API, a warm persistent product route, or a matrix
reuse performance comparison. Finite-field GUDHI corroborates nullities only
on these declared torsion-free fixtures; it does not validate real spectra.

## Reproduce

Use Python 3.10+ and NumPy. Optional `--gudhi` requires GUDHI and fails when it is
missing. The retained reference environment uses NumPy 2.2.6 and GUDHI 3.12.0.
The Rust library adds no dependency. Download the two pinned PDFs into a separate
ignored directory as `scsf-v2.pdf` and `chfsi-v2.pdf` and check the hashes above.

From the repository checkout:

```sh
cargo build --locked --release
rustc --edition 2024 -D warnings -O tools/reference/subspace_cocycle.rs --extern cocycle=target/release/libcocycle.rlib -L dependency=target/release/deps -o target/subspace_cocycle
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 tools/reproduce_subspaces.py --worker target/subspace_cocycle --papers target/subspace-papers --gudhi --output target/subspace-run-001
python3 -m unittest discover -s tools -p test_reproduce_subspaces.py
```

On Windows append `.exe` to the worker and set thread variables through the
shell. Output directories must be new. The controller preserves source SHA and
dirty hashes, package versions, worker/PDF hashes, exact input hashes, all
iteration blocks, printed failures, fresh residuals and raw worker fixtures and
stdout/stderr. Default tool tests skip these numerical checks when NumPy is
unavailable; an explicitly invoked reproduction never silently skips NumPy.

The `Controlled sequence subspace comparisons` CI job reruns the numerical
controls, Rust replay and GUDHI checks with NumPy 2.2.6 and GUDHI 3.12.0. It
uploads the raw `subspace-reference-<sha>-<attempt>` artifact even on failure,
with 14-day retention. CI uses the documented mathematical source pins without
downloading the PDFs; `--papers` additionally enforces their exact hashes locally.

## Retained numerical run

The local Windows run used clean implementation commit
`271c7ef2b4b0779946220358bde1c7a9f1886056`, Python 3.10, NumPy 2.2.6 and
GUDHI 3.12.0, with the three declared thread variables set to 1. Its manifest
retains exact versions, code/PDF/worker hashes and fixture hashes. The selected
rewrite replayed **334** fixed/mixed-degree and actual iteration filter inputs;
maximum normalized filter disagreement was `9.91e-16`. All **32** cold/warm
Rust sequence solves passed the known spectrum and repeated-cluster projector
checks. All **29** tiny persistent pair/degree controls passed exact-rational
expectations and GUDHI topology; the largest mapped current-operator residual
was `2.16e-15`.

Diagnostic counts sum vector actions in filtering, Rayleigh-Ritz and fresh
residuals across eight samples, with the same accuracy and block size:

| Family | Cold ChFSI | Original-order warm | FFT-ordered SCSF |
| --- | --- | --- | --- |
| Smooth rotations | 3194 | 3004 | 2764 |
| Abrupt independent rotations | 3104 | 3124 | 3104 |

These counts do not measure sorting, mapping, verification or total runtime;
they do not select a faster default. The adverse family retains the lack of an
FFT advantage. The defining-polynomial check rejects the printed recurrence
at degree 2 (maximum entry error about 6.32), while the corrected recurrence
passes the independent scalar controls.

Local validation passed Clippy, complete debug/release and Rust 1.91 tests/checks,
the CONTRIBUTING examples/doctests, strict rustdoc/Sphinx and source/link checks.
The 157-test Python suite passed on Ubuntu WSL with 39 optional skips; all six
new numerical tests separately passed in the pinned Windows reference
environment. Its Windows full-suite run retains an existing pipeline test's
unsupported `resource` import rather than modifying unrelated code. Hosted CI
results and artifact identities belong to the associated PR, separately from
this local run.

Timing fields include diagnostic trace preparation and are not formal benchmark
results. No optimization or changed default is selected. The remaining Issue
148 work is adaptive companion ChFSI, the broader paper model/data replay,
cross-track positive-projection integration, and committed, serial, matched
four-way matrix/vector reuse measurements under the repository's reporting
protocol. Keep the issue open for those gates.
