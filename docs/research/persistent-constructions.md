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
