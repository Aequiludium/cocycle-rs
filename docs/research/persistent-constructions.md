# Persistent constructions and rectangular least squares

[Documentation](../README.md) / Research

Issue [#144](https://github.com/Aequiludium/cocycle-rs/issues/144), SPSC-A.
This reproduction specializes to ordinary, nonaugmented signed simplicial
chains over the reals with standard inner products and an actual same-source
inclusion. Filtration values select cells; they are not weights.

## Algorithm mapping and independent expectations

Mémoli, Wan and Wang, *Persistent Laplacians: Properties, Algorithms and
Implications*, [arXiv:2012.02808v3](https://arxiv.org/pdf/2012.02808v3),
Algorithm 1 (physical PDF page 10) column-reduces the outside-row boundary B2,
retaining Y and its zero-column kernel Z. With T=B1 Z its upper term is
T (Z^T Z)^-1 T^T. Algorithm 2 (page 13) instead takes the generalized Schur
complement of the outside block of B B^T. All upper columns are retained.
The lower term is A^T A for A=boundary_q(K) in both algorithms.
The induced Gram is necessary: a nonorthogonal Z does not give projector Z Z^T.

The independent expectation uses Fraction row RREF/free-variable kernels,
exact Gram inversion and rational matrices, independently of the paper's
column reduction and floating Schur solve. A separate NumPy SVD projector and
the existing Rust dense Schur implementation are additional comparators.
Input/oracle generation never calls the candidate Rust operator.

Fong and Saunders, *LSMR: An Iterative Algorithm for Sparse Least-Squares
Problems*, [author PDF](https://web.stanford.edu/group/SOL/software/lsmr/LSMR-SISC-2011.pdf),
Section 2.8 (physical page 5), steps 1-6 are independently transcribed for
M=B2^T and zero damping. Fresh residuals provide a project stopping rule;
Section 3's inexpensive norm estimators are outside this reproduction.
Each solver-control iterate is compared with SciPy LSMR at the same iteration.
Golub--Kahan identities, iterates, true residuals and normal residuals are saved.
Iteration exhaustion raises an error rather than returning a projection.

## Replay and handoff

Use Python 3.10+ with NumPy 2.2.6, SciPy 1.15.3 and GUDHI 3.12.0. These are
reference-tool dependencies, not Rust runtime dependencies. Fetch the two PDFs
above and the unchanged `persistLap.m` / `SchurComp.m` from
[author commit ebff6e3](https://github.com/ndag/Persistent-Laplacian/tree/ebff6e3ad6f5eb90eb8058e4f5567755b62f0c3e)
into `target/round1-sources`. The controller rejects a PDF hash or Git blob
mismatch. Frozen PDF SHA-256 values:

- MMW v3: `a4fff3f81ca41d78a43d3ac88e718540b21337e97f026a324f2942bf9751aebd`.
- LSMR: `d456214a23c51a8506a5e4f4cbc935a1f05ea0f2559755eba9847927e9974b82`.

```sh
cargo build --locked --release --example persistent_laplacian
python3 -m unittest discover -s tools -p 'test_reproduce_persistent_constructions.py'
python3 tools/reproduce_persistent_constructions.py --sources target/round1-sources --output target/spsc-a-round1 --worker target/release/examples/persistent_laplacian --gudhi
```

Use a fresh output directory. On Windows the worker has an `.exe` suffix.
`inputs.json` freezes labeled face-closed filtrations before execution;
`policy.json` freezes tolerances. `index.json` identifies each pair artifact,
which contains ordered bases, retained/outside row maps, integer boundaries,
exact kernel transformation, induced Gram, Schur blocks, exact ranks,
U/D/Delta, full unclipped spectra and eigenvectors, frozen action vectors,
LSMR traces and separately labeled composition diagnostics.
`manifest.json` hashes the files; `environment.json` records versions, sources,
command and checkout status. Failures remain in `failures.json` and make the
command fail. Other work packages can read this format directly.

The corpus includes paths, filled/unfilled triangles, a tetrahedron and its
shell, empty/singular blocks, K=L, tied/delayed births, fixed seed 124
face-closed filtrations and RP2. Separate seed 144 signed/permuted coordinate
controls preserve the coordinate conjugation. GUDHI corroborates persistent
ranks only on declared torsion-free controls. RP2 records the F2/F11 distinction;
finite-field persistence is not a real-spectrum oracle.

## Bounded first-round result

On 2026-11-24, Windows Python 3.10.11 with the versions above checked 801 pairs
in q=0/1/2 from 16 filtrations, 12 rectangular solver controls and 48 signed/
permuted controls. Two constructions, exact rational expectations, SVD and
existing Rust Schur agree. Maximum matrix Frobenius difference was
1.9793003835235033e-14; maximum normalized fresh eigenpair residual was
1.178123921054284e-15. The predeclared matrix gate is
1e-10*max(1,norm(reference,Fro)); fresh eigenpair residual is
norm(Delta*v-lambda*v)/max(1,norm(Delta,2)+abs(lambda)), bounded by 1e-9.
These gates cover unit-scale project controls, not ill-conditioned ranks.

This is a Python formula reproduction; MATLAB/Octave execution was unavailable.
The unchanged author source was inspected and its exact blobs verified. PETLS
was not run. This is not author-dataset or author-runtime reproduction.
The older 1,114-pair suite remains a separate regression, not a new result here.

The LSMR-based outer composition is a project experiment. Repeatability,
linearity, bilinear symmetry, nonnegative quadratic forms, kernel residuals,
idempotence and tighter-inner-action diagnostics passed the selected controls.
They do not certify a fixed linear/self-adjoint operator in general. Final
eigenpair residuals use the independent rational-oracle matrix. No performance
or scalable complete-spectrum claim follows from this reproduction.

## Rust rewrite and representation admission

The unoptimized Rust rewrite is commit
`b70e984` (`fixed-projection action reference`): dense rectangular incidence,
two-pass row-space reduction and a fixed I-Q^T Q action. Its 801-pair comparison
uses the same frozen project inputs and exact-oracle matrices. A separate
representation change replaces only dense incidence with sparse signed columns;
the projection policy, coordinate metadata and mathematical output stay fixed.
At that fixed-action increment LSMR remained in reference tooling. The separate
Round-2 diagnostic below adds a native terminal without changing the fixed action.

The nontrivial scale gate uses K containing 257 leaves and L adding one late
center and all leaf-center edges, q=0, scales (0,2). Its hand-derived action is
x-mean(x), i.e. I-11^T/257. The sparse path retains 257 incidence terms and 257
projection entries; the default dense cap remains 65,536. Dense materialization
and the original dense constructor reject this pair under that cap. This is
operator-only evidence; full spectra and general sparse projector scalability
remain outside the acceptance claim. Q still requires O(p*rank(B2)) storage.

The proposed bounded adoption criterion is: preserve all 801 oracle comparisons
and failures, pass the nontrivial 257 gate under unchanged dense limits, and
reduce incidence storage with no action-time regression on matched star controls
of 64, 128 and 256 leaves. Measure one and 1000 actions separately, with a
retained warmup and ten fresh serial processes per cell. Source preparation,
operator setup, actions, process wall time and absolute peak working set are
separate. These project controls do not support a native cross-library ranking
or a universal speed claim. Formal measurements follow the existing
`spectral-resource-v1` controller; retain every sample/failure and compare paired
samples, rather than dividing unrelated medians.

## Round-2 Rust LSMR rewrite

`project_boundary_kernel_lsmr` in the existing persistent-Laplacian owner is
the unoptimized dense-boundary rewrite of Section 2.8. Its options preserve the
frozen zero-start, zero-damping recurrence and fresh normal-residual stopping
rule. The result retains numerical solution/residual coordinates and diagnostics;
it does not certify exact rank, projection forward error or self-adjointness.
All validation, multiplies, rotations and final residuals share one work budget.
No partial successful result is returned on nonconvergence or interruption.

The opted-in worker compiles the same production LSMR and budget source with
`cocycle_lsmr_replay`. It records each update without adding trace allocation or
observation callbacks to ordinary library debug/release builds. This is an
equivalence round, with no sparse backend, changed recurrence or timing claim.

After the earlier paper command has frozen its output, run:

```sh
cargo build --locked --release
rustc --edition 2024 --cfg cocycle_lsmr_replay -C opt-level=3 tools/reference/persistent_lsmr.rs --extern cocycle=target/release/libcocycle.rlib -L dependency=target/release/deps -o target/persistent-lsmr-replay
python3 tools/reproduce_persistent_constructions.py --replay-lsmr-worker target/persistent-lsmr-replay --frozen target/spsc-a-round1 --output target/spsc-a-lsmr-rust
```

The original manifest is verified before reading frozen expected values. The
controller never calls Rust to generate them. It retains input/worker/source
hashes, complete stdout/stderr, update comparisons and failures in a fresh output
directory. Windows adds `.exe` to executable paths. Dedicated Rust tests cover
minimum-norm rank-deficient fits, inconsistent residuals, empty/stationary inputs,
caps, input errors, numerical overflow, iteration exhaustion, cancellation after
an actual update and cumulative final-residual charging.
The simple unscaled recurrence explicitly rejects overflowing rotation products
and underflowed required scales; no rescaling optimization is included.

On 2026-11-24, the original 801 frozen pairs and 12 solver controls preserve
update counts and solutions. Five explicit failure controls reject iteration
exhaustion or invalid/nonfinite input. All 208 updates satisfy independently
evaluated Golub--Kahan relations (maximum error 1.021e-15), fresh x/residual/
normal-residual checks and dense least-squares comparisons. Maximum accepted
iterate/rotation comparison error is 1.302e-13; the separately labeled outer
composition's maximum action error is 3.471e-15 against the original rational
matrices. The original 1e-10 comparison tolerance is unchanged.

### Retained terminal-breakdown discrepancy

An initial strict comparison failed on 25 cases' terminal alpha or alpha*v.
Their terminal beta is at rounding scale. Normalizing an almost-zero
M*v-alpha*u can produce unrelated u directions; alpha and v subsequently
depend on that unresolved direction. The mathematical recurrence is identical,
but these raw continuation variables are not reproduced across the NumPy BLAS
and Rust scalar evaluation orders. The failed attempt remains retained.

The accepted scope is explicitly narrower: normal steps compare raw variables;
terminal steps with both beta values below sqrt(machine epsilon)*max(1,norm(M,Fro))
compare beta*u, beta*alpha*v, rotations, solutions and fresh residuals, and check
each implementation's own Golub--Kahan relations. Such a step must be the final
update under the unchanged stopping rule. No nonterminal breakdown is exempted.
There are 46 such terminal steps in the local replay, of which the 25 raw-variable
discrepancies are retained in `terminal-breakdown-differences.json`. This is a
scoped numerical reproduction, not bitwise recurrence-state equivalence or an
exact breakdown/rank certificate. A caller needing a fixed action continues to
use the prepared row-space operator.
