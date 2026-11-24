# Gram flips and harmonic restriction: bounded reproduction

[Documentation](../index.md) / Research

This is the first-round paper reproduction for
[Issue #145](https://github.com/Aequiludium/cocycle-rs/issues/145), part of
[Issue #143](https://github.com/Aequiludium/cocycle-rs/issues/143). The callable
entry is `tools/reproduce_split_spectra.py`. It is optional verification tooling,
outside the Rust runtime and crate package. There is no new Rust spectral API,
default dispatch, cache or scalable solver in this change.

## Sources and scope

The independent mathematical transcriptions cover:

- Jones and Wei, *PETLS: PErsistent Topological Laplacian Software*,
  [arXiv:2508.11560v2](https://arxiv.org/abs/2508.11560v2), Sections 4.4 and 4.5,
  PDF pages 28 and 29--31. Neither section has a separately numbered algorithm.
- Duersch, Shao, Yang and Gu, *A Robust and Efficient Implementation of LOBPCG*,
  [arXiv:1704.07458v3](https://arxiv.org/abs/1704.07458v3), Algorithms 1 and 2,
  PDF page 4: basic LOBPCG and Rayleigh--Ritz for B=I and identity preconditioning.
  This does not reproduce robust Algorithm 8.

The pinned PDF SHA-256 values are respectively
`5e14d20b0230d98550be0b802d0dfe9635d232bf794a33672b85fa151706fb86` and
`01d447acaea3d0ec49a427e554e499791a5e23c0f38fcdb8c50e159ac21bbc07`.
Passing `--references` verifies the actual local files against these values.
Without it, the output explicitly records that paper identity was not verified
in that run.

Author source inspection uses
[PETLS commit 5b6564e9e05af9d24928313b92a014696247bdd7](https://github.com/bdjones13/PETLS/tree/5b6564e9e05af9d24928313b92a014696247bdd7),
particularly `cpp/src/core/Complex.cpp`, `cpp/include/petls_headers/core/Complex.hpp`,
`cpp/include/petls_headers/core/up_algorithms.hpp` and the eigenvalue routines.
The root contains no license granting reuse of PETLS's original implementation
in this inspected snapshot; dependency licenses do not grant that permission.
No author code is copied into this contribution.

The original Windows CMake configuration was actually attempted, with Alpha and
LAPACK disabled. It failed while finding Boost `program_options` and `thread`.
The retained `petls-configure.log` documents the failure; neither the author's
`get_L_top_dim_flipped` nor `nonzero_spectra` runtime is established. The verdict
is **reproduced within a narrower scope using independent transcriptions**,
not verified execution of the PETLS software or reproduction of its datasets.

## Mathematical checks

For an integer boundary A, both complete Grams A^T A and A A^T are retained.
Strictly positive eigenvalues are compared as multisets. Original-domain zeros
are restored using numerical rank, checked independently against exact rational
rank. Positive dual modes are lifted by `v=A^T z/sqrt(lambda)` and checked in
the original coordinates. The comparison of repeated modes uses the whole
positive-space projector, rather than vector signs. Empty matrix shapes are
preserved explicitly in JSON.

The separate project extension constructs `F=B1 Q` with Q orthonormal in
`ker(B2)` and checks `A F=0`, `UD=DU=0`, `U=F F^T`, `D=A^T A`. An independent
rational kernel Z with its induced metric validates U. Only positive component
multisets are merged, with `n-rank(F)-rank(A)` zeros. The filled triangle has
U spectrum `[0,0,3]`, D spectrum `[0,3,3]`, and exactly three full values
`[3,3,3]`. The persistent control has source vertices 0--4, a filled triangle
on 234, a later path 0--5--1 and a disconnected later triangle on 678. All eight
upper-boundary columns are retained. Its degree-zero full spectrum is
`[0,0,1,3,3]`.

Harmonic restriction starts from a separately known complete real kernel N,
normalizes it and constructs an orthonormal complement X. It compares the full
positive spectrum with `X^T Delta X`, then lifts vectors to the original domain.
This is the orthonormal special case of the paper's change of basis. The
two-vertex example explicitly rejects `z=[1,0]` as a harmonic constraint for
`Delta=[[1,-1],[-1,1]]`: `Delta z=[1,-1]`. It accepts the actual kernel
`span([1,1])`. This is an input-contract counterexample, not a software bug
claim about an unexecuted PETLS route.

Basic LOBPCG independently performs initial Rayleigh--Ritz, fresh
`R=H X-X Lambda`, `W=R`, a generalized projected solve on `S=[X,W,P]`, and the
paper's coefficient updates for X and P. Algorithm 2 includes diagonal metric
scaling and checked Cholesky. Dependent or ill-conditioned search spaces raise
an error. An already converged initial block returns before constructing zero
search columns. A separately selectable project correction applies rank-revealing
SVD to normalized search columns and maps coefficients back to raw S. The
`compress` branch in `rayleigh_ritz` is the complete correction; it is not
silently substituted for the paper-basic run or labeled Algorithm 8. Both runs
retain X, P, R, W, S, the projected metrics, Cholesky factors and coefficients.

## Replay

Use an isolated Python environment; dependencies are not Rust dependencies:

```sh
python3 -m pip install numpy==2.2.6 scipy==1.15.3 sympy==1.14.0 gudhi==3.12.0
python3 -m unittest discover -s tools -p 'test_reproduce_split_spectra.py' -v
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 python3 tools/reproduce_split_spectra.py --gudhi --output target/split-spectral-run-001
```

For a paper/source-verified replay, prepare a separate ignored references
directory containing the two PDFs under the names `petls-v2.pdf` and
`lobpcg-v3.pdf`, downloaded from the versioned arXiv PDF URLs. Pass its path
with `--references target/split-spectral-references`. Optional author-source
checking additionally expects `petls-source.zip` from the fixed-commit
[GitHub archive](https://github.com/bdjones13/PETLS/archive/5b6564e9e05af9d24928313b92a014696247bdd7.zip)
and its extracted directory
`petls-source/PETLS-5b6564e9e05af9d24928313b92a014696247bdd7`.
Selected extracted author files are checked byte-for-byte against the archive;
their hashes and the archive hash are recorded. Extraction on Windows must skip
the archive's `:Zone.Identifier` metadata entries. Author configuration logs
are copied when present. This does not run or certify the author software.

Every output directory must be new. Inputs, shapes, hashes, environment, Git
revision/dirty status, source hash, complete spectra, lifted vectors, iteration
traces, warnings and failures are retained. Floating outputs use double
precision. The rank cutoff is `1e-12 * spectral_scale`; comparison errors use
Frobenius norm divided by `max(1, reference Frobenius norm)`. Iterative residuals
use `||H v-lambda v||2/max(1,||H||2)` with tolerance 1e-9. No tolerance is changed
after seeing a failure. Tiny positive values are retained and reported as
ambiguous near the rank threshold; they are not replaced by zeros in the
restriction output.

## Observed results and limits

The 2026-11-24 Windows run used Python 3.10.11, NumPy 2.2.6, SciPy 1.15.3,
SymPy 1.14.0 and GUDHI 3.12.0. The base was main commit
`a3c3e1132800624743ad0f159c9303090e0b3b3f`, plus the script fingerprint recorded
in each run. Sixteen focused tests and all fifteen reproduction groups passed,
including expected failure behavior.

| Solver control | Paper-basic result | Separate rank correction |
| --- | --- | --- |
| Separated spectrum, 48 coordinates, block 3, request 2 | Converged, 44 updates | Converged, 44 updates |
| Clustered positives, 48 coordinates, block/request 3 | Stopped on ill-conditioned search | Converged, 76 updates |
| Lowest root multiplicity 5, block/request 2 | Converged, 76 updates | Converged, 76 updates |

Successful returned modes had normalized original-operator residuals below
1e-9 and agreed with dense `eigh` and SciPy LOBPCG. SciPy did not take its small
dense fallback on these 48-coordinate controls. The large-zero-block control
also runs a restricted iteration and validates its lifted vectors against the
original 35-coordinate operator. Zero rank, excessive requests, incomplete or
dependent kernels, tiny positives, failed Cholesky and exhausted iterations
are covered separately. GUDHI corroborates only torsion-free triangle, shell
and edge topology; finite-field PH cannot establish real spectral values.

Construction, exact-kernel work, factor generation, eigsolve, recombination and
restriction projection have separate timers. Group totals also include oracle
checks and serialization; the summary excludes final manifest export. These
single-run diagnostic times are not comparable optimization measurements or
evidence of a speedup. Numerical cutoff ranks and small residuals are not exact
rank or smallest-index certificates. A request smaller than a repeated cluster
does not establish its complete multiplicity.

The `Gram and positive-spectrum reproduction` CI job replays the mathematical
controls and uploads raw results with 14-day retention. A configured job is not
a passing hosted run; inspect its exact revision. Historical local raw runs
are not a public archive. Rust parity and subsequent optimization remain
separate rounds: do not adopt the uncorrected basic solver on the basis of
the corrected run, change `spectrum()` completeness, or claim that Gram
flipping removes the existing target-domain dense Schur allocation.
