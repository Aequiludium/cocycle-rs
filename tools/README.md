# Developer tools

The optional [sequence-subspace reproduction](../docs/research/sequence-subspaces.md)
checks corrected fixed-degree SCSF/ChFSI controls, frozen Rust Fourier/QR/Ritz/sequence parity,
private supplied-block LOBPCG and tiny changing-basis experiments. It preserves
printed-recursion discrepancies and adverse cases; it is not a public spectral
API, full author-data reproduction or performance adoption result.

The maintained resource snapshots for [persistent Laplacians](../benches/reports/persistent-laplacian-resources.md),
[Dirac](../benches/reports/dirac-resources.md), [Hodge](../benches/reports/hodge-resources.md) and
[spectral features](../benches/reports/spectral-feature-resources.md) record
measured revisions, limits and adverse cases.

The optional [spectral resource controller](benchmark_spectral_resources.py)
runs at least ten fresh serial Windows processes per declared cell, retains one
discarded warmup and a seeded cell order, and records solve/preparation times,
whole-process wall time, absolute peak working set, failures and immutable source
identities. Commit measured code first. Supply a JSON list of explicit fixture
descriptions, commands, independent expected values and tolerances through
`--specification PATH --output NEW_DIRECTORY`. Specifications and raw outputs stay
under ignored `target/benchmarks/commit-<sha12>/spectral-resources/run-<NNN>/`.
See the [performance reporting rules](../benches/reporting.md): these resource
snapshots do not establish native cross-library rankings or general speedups.

These tools are outside the Rust library and crate package. See
[contribution rules](../CONTRIBUTING.md) and [benchmark navigation](../benches/README.md).

## Ordinary HeCS reproduction

`reproduce_hecs.py` implements the bounded ordinary q=1 upper-Hodge experiment
from [Issue #150](https://github.com/Aequiludium/cocycle-rs/issues/150).
It is experiment infrastructure; the bounded library terminal is documented in
the [Hodge guide](../docs/guides/hodge-laplacian.md#bounded-ordinary-triangular-hecs-solve).
Persistent/eigensolver integration is separate. The mathematical source is Savostianov,
Tudisco and Guglielmi, [arXiv:2401.15492v1](https://arxiv.org/abs/2401.15492v1),
PDF Algorithms 5.1/6.1 (physical pages 10/15), Sections 6.2, 7.1 and 7.3.
The PDF SHA-256 is
`20f7fa2ba91d308a97b1ad98f9d688a4a6ae00ce630bd457dd9e71f1e6c6d537`.

Keep downloaded sources and generated evidence under ignored `target/`:

```sh
python3 -m venv target/hecs-venv
target/hecs-venv/bin/python -m pip install numpy==2.2.6 scipy==1.15.3 gudhi==3.11.0
git clone https://github.com/COMPiLELab/HeCS target/hecs-author
git -C target/hecs-author checkout cd4d23e5d53900ece219fc0cdafc5c5e1eb81bf3
curl -L https://arxiv.org/pdf/2401.15492v1 -o target/hecs-paper.pdf
target/hecs-venv/bin/python tools/reproduce_hecs.py --output target/hecs-round1 --author-source target/hecs-author --paper target/hecs-paper.pdf --julia /path/to/julia
```

On Windows the venv interpreter is `Scripts/python.exe`; the tested author
runtime is Julia 1.10.10. The adapter uses only its LinearAlgebra, SparseArrays
and SHA standard libraries with one BLAS thread. It checks source hashes and
includes the unchanged pinned files; it does not run the plotting/large loops
or install their unrelated packages. No author source is redistributed; a
top-level license was not verified. Rust and Python implement the mathematics
independently. `--fixtures PATH` replays saved JSON and `--case NAME` restricts a
diagnostic run. Every output directory must be new.

Both direct heavy selection and Section 6.2's preliminary 2-core peeling are
exercised. The complete original vertex/edge skeleton stays in each fixture.
Weight-descending order uses original triangle ID to resolve ties. Python/Rust
choose the smallest free edge, while Julia's Set may choose another valid one;
each certificate is replayed independently from vertex containment. Save every
trial, rejection, selected triangle, collapse, permutation and factor. With
`B=boundary_2`, weights mean `U=B diag(w) B^T`, never `B diag(w^2) B^T`.

The compact factor is `C=(B[:,T] diag(sqrt(w[T])))[edge_order,:]`, shape m by r,
with a nonsingular lower triangular r by r leading minor. Its weighted Gram
identity and selected rank are checked independently. For the weighted shell
`w=[1,2,4,8]`, faces 1,2,3 (zero-based IDs) are retained; `B0=B1-B2+B3` also gives
an independent analytic active-matrix check. The ambient m by m upper Hodge is
singular: only the r by r matrix `C+ U_permuted C+^T` is solved. A dense SVD oracle
checks both the original equation and the lifted minimum-norm HeCS solution.
GUDHI corroborates supported topology and simplex counts, never real spectra.

The real-projective-plane control has rational rank 10, but heavy selection has
rank 9. It must return `rank_loss`, rather than solve an incomplete active space.
The Rust worker certifies image coverage by checked i128 unit-triangular solves
and exact reconstruction of every original boundary column in every edge row.
Overflow rejects the certificate. Its dense active Gram solve rejects ambiguous
pivots; it adds no ridge and performs no fallback to the ambient inverse.

CGLS starts at zero. Identity, diagonal and HeCS use the same compatible original
RHS and fresh `||Ux-b||2 / max(1,||b||2) <= 1e-9` stopping target, with at most 4000
iterations. The diagonal scale is `1/sqrt(Uii)` for positive entries, identity
for zero entries. It may produce a different kernel component; compare its
solution modulo ker(U). HeCS uses the Moore-Penrose lift. Traces retain iterates,
alpha/beta/gamma and original/normal/recursive residuals. Zero RHS, breakdown,
iteration limit and nonfinite arithmetic have explicit outcomes. The author's
zero-RHS nonfinite result and empty-collapse exception are retained as source
limitations. The adapter separately records that the README's `Perm` orientation
and acceptance-order triangle lookup need not yield a triangular factor; the
paper's explicit row/column indexing is the parity target.

After Round 1 passes, compile the independent Rust worker and replay those bytes:

```sh
rustc --edition 2024 -D warnings --test tools/reference/hecs.rs -o target/hecs-tests
target/hecs-tests
rustc --edition 2024 -D warnings -C opt-level=3 tools/reference/hecs.rs -o target/hecs
target/hecs-venv/bin/python tools/reproduce_hecs.py --parity-from target/hecs-round1 --rust-worker target/hecs --output target/hecs-round2
```

For Round-2 library parity, build the public example and replay the SAME saved
inputs with `--native-worker` (the standalone worker can also be supplied):

```sh
cargo build --locked --release --example hecs
target/hecs-venv/bin/python tools/reproduce_hecs.py --parity-from target/hecs-round1 --native-worker target/release/examples/hecs --output target/hecs-native-parity
target/hecs-venv/bin/python tools/reproduce_hecs.py --check-native --native-worker target/release/examples/hecs --output target/hecs-native-controls
```

`--check-native` regenerates the independent rational/analytic/SVD/GUDHI
controls for CI and saves them with actual dependency/source identities. It
does not execute Julia or claim a new author-program replay. Native rank loss
is an explicit library error, while successful factors, active matrices,
solutions and original residuals are compared with the independent oracle.

Worker stdin is whitespace-separated: `n m t peel`, m increasing vertex pairs,
t increasing vertex triples followed by positive finite weights, and m finite
RHS entries. Vertex IDs are in `0..n`; duplicate simplices, missing lower edges,
trailing input and limits n>256, m>512, t>1024 are rejected. `peel` is 0 or 1.
The standalone reference does not expose library execution controls. The native
terminal adds checked shared execution budgets and per-buffer entry caps;
neither route provides an aggregate memory quota. Existing Hodge/Dirac spectral
wrappers retain their behavior.

The only Round-3 candidate is `--adjacency` on the worker: dynamic adjacency sets
replace rescanning for free edges. Round-2 parity checks both paths against the
same certificates/factors/active matrices/statuses and original residuals. For
a measured comparison, commit the worker/controller first, then run:

```sh
target/hecs-venv/bin/python tools/reproduce_hecs.py --compare-from target/hecs-round1 --rust-worker target/hecs --samples 10 --output target/hecs-collapse-comparison
```

The paired serial experiment retains one warmup and ten fresh-process samples
per route on nine regenerated enriched triangulations and both peeling modes.
Seed 150 fixes points, deleted/added edges, triangle weights and sample order;
all arrays are exported. These are smaller Section 7.3-style inputs, not the
authors' original data. Incomplete Cholesky is not run, so the comparison is
narrower than the paper. Phase timings, factor nnz and failures are retained;
process startup, I/O and final JSON export are outside worker timings. Dense
reference costs and diagnostic trace construction are included. Peak process
memory, amortized reuse and persistent inner/outer adaptation are unmeasured.
No general speedup, sparse-scaling, production adoption or persistent guarantee
follows from this experiment. Run source/document checks, the Python tool tests,
standalone rustfmt/rustc tests and Rust MSRV checks for changes to this path.

The existing Hodge checker also accepts `--path-convergence`: it records
independent SciPy identity, inverse-diagonal and shifted-tridiagonal-inverse
controls on the 1024-node path against its complete analytic spectrum, using the
same SplitMix64 starting span. Warnings and incomplete lookahead convergence
remain explicit; these controls are numerical diagnosis, not native timing.

## Ordinary-Hodge sparsification reproduction

The optional [sampling/KID controller](reproduce_sparsification.py) and
[bounded Rust worker](reference/sparsification.rs) implement the reproducible
subroutines in [Issue #149](https://github.com/Aequiludium/cocycle-rs/issues/149).
They leave the library's exact, unweighted operators unchanged. Sampling retains
the complete lower skeleton and original coordinates; filtration values do not
become weights. These tools are correctness experiments, outside the crate.

The sources are Osting, Palande and Wang,
[arXiv:1708.08436v3, Section 3/Theorem 3.1](https://arxiv.org/abs/1708.08436v3),
and Savostianov, Schaub, Guglielmi and Tudisco,
[arXiv:2502.07558v2](https://arxiv.org/abs/2502.07558v2), PDF Algorithm 3.1
(physical page 8), Section 4.3 (pages 11-13) and Algorithm 4.1 (page 14).
PDF versions and SHA-256 are enforced before a run. No author implementation
was verified in these paper links; this is an independent mathematical rewrite.

For an integer boundary B, an independent rational row-space projector supplies
the exact rank/kernel and column leverage scores. A numerical pseudoinverse of
G = B B^T supplies a separate comparison. Positive-rank probabilities are
p_j = tau_j/rank. Sampling with replacement retains the seed, every draw,
repeated counts and Gram weights count_j/(k p_j); factor columns use the square
roots of these weights. Zero rank returns zero without probability normalization.
Uniform sampling is a separate control. Every trial retains the complete
generalized spectrum on the positive space of G and explicit kernel checks,
including rank loss and failed PSD sandwiches. Epsilon targets 0.25/0.5 and delta
target 0.05 are experiment metadata: no theorem-backed sample count or failure
probability is asserted, individually or jointly across inputs.

The four-face tetrahedron shell has B^T B = 4I - ss^T, s = [1,-1,1,-1], rank 3,
scores 3/4 and probabilities 1/4. For H = B^T B/4, all 16 Rademacher probes
recover odd moments 2 diag(T_m(H)) = 3/2 at degrees 1,3,5,7. A one-face draw
loses rank and is retained as a failure. The printed Algorithm 4.1 is replayed
separately: its upper-Gram/probe dimensions disagree on this rectangular input;
even after substituting the dual Gram, its indexing, conditional recurrence
update and factor differ from Section 4.3. These are observed discrepancies,
not author-confirmed errata. The corrected branch updates every degree and
collects only odd moments, checking scalar polynomials before sampling.

Full paper KID probability reconstruction remains unresolved. Its
`histogramIntegral` does not specify the executable bin/integration policy.
The separate **project diagnostic** fixes 64 positive bins of width 1/64,
undamped Chebyshev expansion, analytic bin integration and normalization by
the raw score sum. Its bin averages are not asserted to reproduce the paper's
convolution kernel. Raw scores, sums relative to exact rank, negative bins and
negative probabilities are retained. Invalid sampling probabilities are rejected
without clipping. This branch's probability errors and empirical sandwiches
cannot establish the paper's probability theorem or full KID reproduction.

The frozen project family contains eight points in two seeded Gaussian clusters,
four inclusive Rips thresholds and q=0/1. It is not author data. The full sweep
uses maximum polynomial degrees 8/16/32 (4/8/16 odd moments), probes 32/128/512
and seeds 0 through 19; `--quick` uses seed 0 only. The same nested probes and
draw seeds are retained across settings. GUDHI checks input simplices, while
GUDHI and Ripser field-2 Betti0/1 corroborate exact rational ranks on these inputs;
neither library is a real-spectrum oracle. Rust parity covers numerical scores,
fixed-draw Gram assembly and corrected moments at all nine degree/probe settings
with seed 0 on every fixture, plus shell, empty and zero controls. It does not
claim native RNG parity, a Rust histogram implementation or a persistent extension.
The worker reuses the existing checked Hodge eigensolver through a formal
two-degree chain window; its numerical rank threshold is 1e-10 max(1,lambda_max).

One tested environment uses Python 3.12.13, NumPy 2.4.6, SciPy 1.17.1,
GUDHI 3.12.0 and Ripser 0.6.12. Create an isolated environment, save its complete
dependency lock, and download both pinned PDFs to `target/sparsification-sources/`
using their arXiv PDF URLs. Then run from the repository root:

```sh
python3 -m venv target/sparsification-venv
target/sparsification-venv/bin/python -m pip install numpy==2.4.6 scipy==1.17.1 gudhi==3.12.0 ripser==0.6.12
target/sparsification-venv/bin/python -m pip freeze > target/sparsification-lock.txt
cargo build --locked --release
rustfmt --edition 2024 --check tools/reference/sparsification.rs
rustc --edition 2024 -O tools/reference/sparsification.rs --extern cocycle=target/release/libcocycle.rlib -L dependency=target/release/deps -o target/sparsification-worker
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 target/sparsification-venv/bin/python tools/reproduce_sparsification.py --sources target/sparsification-sources --worker target/sparsification-worker --output target/sparsification-run-001
target/sparsification-venv/bin/python -m unittest discover -s tools -p test_reproduce_sparsification.py
```

On Windows use `Scripts/python.exe`, append `.exe` to the worker and set the
three thread environment variables in PowerShell before starting Python.
The PDFs must hash to `373640de5a96c8fce16784aecc072ef240cdacb96476fffde7ef3263426fa8f8`
(2502.07558v2) and `1470ef2d8785dc9872f525479253b00e6d7d84b6b5b08ef8804b981e0983305e`
(1708.08436v3). An existing output directory is rejected. Sources, inputs, probes,
moments, all trials, Rust requests/results/rejections and environment identities
are retained under the new output. Phase durations and Python-allocation peaks
are diagnostics; they are not formal benchmarks, process RSS or speedup evidence.
Core integration, persistent Schur sampling and optimization remain deferred
until the full probability/error contract and matched cost evidence are resolved.

### Round 2: complete the scoped Rust rewrite

The same worker now also accepts
`native n p degree probes draws seed leverage|uniform`, followed by the
row-major integer B and p-by-probes Rademacher Z. Rust computes the scores,
probabilities and draws itself. Leverage probabilities divide the numerical
scores by their measured sum; the raw sum and numerical rank are returned so
the normalization and trace defect remain visible. Uniform sampling is a
separate control. Zero rank returns zero probabilities and no draws.

Native randomness uses SplitMix64 with explicit wrapping 64-bit arithmetic,
its high 53 bits divided by 2^53, and the ordered half-open inverse CDF. The last
positive bin absorbs floating-point endpoint roundoff. This is a declared finite
random-number implementation, with independent golden-word and CDF controls;
it does not claim the same draws as NumPy PCG64 or a probability certificate.
Every result retains the route, seed, actual draws, scores, probabilities,
weights and matrix. The original supplied-probability/draw protocol remains
available for exact Round-1 replay. Negative input probabilities are rejected.

After generating the complete Round-1 output above, use its saved input arrays,
probes, moments, weights and draws for the Round-2 command:

```sh
rustc --edition 2024 -D warnings --test tools/reference/sparsification.rs --extern cocycle=target/release/libcocycle.rlib -L dependency=target/release/deps -o target/sparsification-tests
target/sparsification-tests
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 target/sparsification-venv/bin/python tools/reproduce_sparsification.py --replay-round1 target/sparsification-run-001 --worker target/sparsification-worker --output target/sparsification-round2-001
```

The frozen eight-fixture input SHA-256 is enforced. The full command checks all
1,440 stored moment configurations against both frozen values and independent
dense scalar-polynomial evaluation on the same probes. It replays all 7,036
sampling trials, including the externally supplied project-diagnostic
probabilities, checking the whole generalized spectrum, rank-loss and sandwich
pass/fail states. It separately verifies Rust-generated leverage/uniform draws
for 1,280 seed/count/route cases with an independent Python integer-RNG/CDF
control and rational leverage expectations. The frozen invalid probability
estimate is replayed as a rejection. Shell/empty/zero analytic controls and
protocol rejections remain separate from these counts. `--quick` restricts the
frozen and native sweeps to seed 0.

`--batch` frames independent worker requests, one per line, with the original
2 MB per-request cap and dimension/probe/degree/draw caps. Each request still
reconstructs and solves the same simple dense reference; no factor cache or
algorithm optimization is added. A bad request terminates with a nonzero status;
earlier complete request outputs remain visible. The controller rejects partial
batches and retains exact requests, stdout/stderr, comparisons, fingerprints,
native trials and all failure diagnostics in a new output directory.

This completes only the reproduced mathematical scope. Full paper KID
histogram/probability reconstruction, persistent transfer, public weighted APIs
and Round-3 optimization remain separate admission gates. Native draw outcomes
are empirical evidence, not a replacement for a justified sample-count theorem.
CI compiles and runs the standalone analytic/protocol tests; the optional full
numerical replay additionally needs the declared Python reference environment
and an existing complete Round-1 output.

## Spectral reference CI

The [Gram and harmonic-restriction reproduction](../docs/research/split-spectral-reproduction.md)
uses `reproduce_split_spectra.py` for PETLS Sections 4.4/4.5 and basic LOBPCG
Algorithms 1/2. The `Gram and positive-spectrum reproduction` CI job runs its
focused numerical tests and controls with pinned NumPy/SciPy/SymPy/GUDHI.
It retains raw spectra, iteration traces and the paper-basic clustered-search
failure separately from the explicit project rank correction. It does not
certify author-software execution, Rust parity or production optimization.
The [bounded split-spectral Rust parity tool](../docs/research/split-spectral-rust.md)
compares `examples/split_spectra.rs` with the hash-pinned first-round #145
transcription. It supports immutable historical inputs with `--frozen` and a
freshly frozen platform replay. The `Split-spectral Rust parity` CI job retains
both ordinary results and expected basic-solver failures; it does not use the
Python SVD correction as a Rust basic-solver expectation.

The `Real Hodge reference comparisons` CI job requires `--iterative
--toponetx PATH --perslap PATH --gudhi` on fresh Rust example output. It checks
matrices, actions, dense eigenpairs and matrix-free Ritz pairs against independent
NumPy/rational expectations, unchanged pinned TopoNetX and perslap functions on
their supported subsets, and GUDHI field-101 topology. Unsupported reference
degrees and the RP2 real/F2 distinction remain explicit.

The job uses Python 3.12 with NumPy 2.2.6, SciPy 1.15.3 and GUDHI 3.12.0.
The `Persistent spectral feature reproduction` CI job additionally runs the
focused numerical tests and complete pinned-author reproduction: 160 worked pair
comparisons, independent SVD/exact ranks, GUDHI topology, the original feature
loop, Rust descriptor workers and both group-held-out random-forest experiments.
It also uses scikit-learn 1.7.2 and threadpoolctl 3.6.0. Full datasets and
new-family generalization remain outside this reproduction.

The `Real persistent Laplacian reference comparisons` CI job requires
`--components --gudhi --reference-python PATH` on a freshly built release worker.
It checks all 1,114 stage pairs against independent rational/NumPy Schur and
restricted-boundary SVD expectations, including upper/lower actions, spectra,
orthogonality and fixed-endpoint results. The pinned author's Schur function
runs on its supported subset; GUDHI supplies the declared finite-field topology
checks, including the explicit RP2 field counterexample. The larger optional
MATLAB author-data replay and performance experiments remain separate.

The `Real Dirac reference comparisons` CI job requires `--pairs --toponetx
PATH --gudhi` on every push and pull request. It builds a fresh release worker,
checks all fourteen fixtures with NumPy, executes the fixed TopoNetX methods on
their supported subset and corroborates complete torsion-free kernels with
GUDHI F3. Empty/truncated reference exclusions remain explicit in the report.

The job uses Python 3.12 with NumPy 2.2.6, SciPy 1.15.3 and GUDHI 3.12.0.
External sources are checked out at the exact revisions linked in the relevant
oracle documentation; no runtime dependency is added. Explicitly requested
missing packages or comparisons fail the job.

The artifacts named `dirac-reference-<sha>-<attempt>`, `hodge-reference-<sha>-<attempt>`,
`persistent-laplacian-reference-<sha>-<attempt>` and
`persistent-spectral-reference-<sha>-<attempt>` retain the actual
checkout SHA, event/run identity, reference commits and source hashes, dependency
versions, Rust outputs and comparison results/logs, including failure evidence.
Artifacts are uploaded even when a comparison fails and expire after 14 days.
These are numerical correctness checks; finite-field corroboration, supported
external function subsets and bounded learning controls do not establish exact
real topology, whole-package compatibility or performance rankings.

## Real persistent Laplacian reproduction

The bounded [relevant-pair reproduction](reproduce_relevant_pairs.py) transcribes
Kirchner--Vittorietti arXiv:2611.09610v1 Algorithms 1-2 over Q. It freezes the PDF
SHA-256, input bytes, R/V columns, event buckets, selected pairs, predecessor
embeddings and every exact upper/full matrix. Its independent row-RREF restricted
kernel oracle comes from `compare_persistent_laplacians.py`, not the new column
reducer. GUDHI 3.12.0 corroborates topology over F11; this is not a real-spectrum
oracle. NumPy 2.2.6 supplies numerical cross-checks. These remain external tools.

```sh
cargo build --locked --release --example persistent_laplacian
python3 tools/reproduce_relevant_pairs.py --paper target/2611.09610v1.pdf \
  --paper-style --rust-worker target/release/examples/persistent_laplacian \
  --output target/relevant-pair-reproduction
```

Obtain the exact PDF from [arXiv v1](https://arxiv.org/pdf/2611.09610v1); its
required SHA-256 is `81731c460bbeb12a210fd5745fa862919d4f9992e30e23c07e08419c483009cc`.
Output directories must be fresh. The tool fails on a different PDF or missing
dependencies. No author implementation has been verified; this is a mathematical
transcription. Controls include joint-bucket events, delayed vertices, ties,
zero columns, RP2 and tetrahedral shells. Every tiny-grid pair, including internal
empty stages, is checked. The optional paper-style inputs use 12 planar alpha
points and 8 spatial Rips points at radius 0.3, project-selected seeds 0/1/2 and
degrees 0-3. These are scaled instances, not author data or Figure 9 timings.
The Rust worker checks every real-grid pair's full matrix, spectrum, ordered
labels, scales, degree and evaluated representative against frozen rational
results. Its `--relevant` mode accepts the existing whitespace fixture protocol,
treating the start list as the entire grid and evaluating all legal pairs.

Round 2 can replay an existing Round-1 directory directly, without regenerating
expected R/V or matrices from the current Python algorithm:

```sh
python3 tools/reproduce_relevant_pairs.py --paper target/2611.09610v1.pdf \
  --replay-frozen target/relevant-pair-reproduction \
  --output target/relevant-pair-rust-replay
```

The controller builds and freezes the library test executable, then explicitly
runs its ignored `replay_frozen_cases` test. Test-only capture hooks in the actual
public constructor export exact rational R/V, events and selected pairs. They
are absent from ordinary library builds. Frozen C intermediates, predecessors,
embeddings, full matrices/spectra and owned coordinate metadata are compared;
empty grids retain their empty output. Cancellation, work/storage exhaustion,
invalid indices and degree overflow must return errors without partial output.
Private arithmetic tests separately exercise checked i128 overflow.

Add `--shared-baseline PATH` for package A's published R0 directory containing
`manifest.json`, `inputs.json`, `index.json` and the indexed `pair-*.json` files.
The tool verifies their original SHA-256 manifest and exact `exact_delta`
matrices, then compares the same Rust batch outputs in matched coordinates.
A's restricted-boundary Y/R are not the final-boundary C transforms, so only
C's frozen transforms are compared directly; R=B V and upper-triangular V are
checked independently for all exported traces. No A code or runtime API is
imported. All read artifact hashes, the test binary, source fingerprints and
raw input/output/logs are retained, and inputs are checked for concurrent edits.

New Round-1 runs include raw-byte input hashes and an artifact manifest. Legacy
Windows runs hashed UTF-8 LF logical input text before CRLF file writing; replay
verifies precisely that legacy hash while also freezing the original raw bytes.
Expected outputs from legacy runs without manifests receive current read-time
hashes, not a claim of a historical signed manifest. Replay is correctness
evidence only and is kept separate from reproduction/resource timing modes.

Serial single-run Python timings compare exhaustive Schur, selective Schur and
new-chain upper construction with identical returned matrices/spectra. They
exclude input construction, oracle validation and export, and remain exploratory
phase diagnostics. They do not establish an adopted optimization, formal timing
acceptance, total process memory or complete-spectrum scalability.

On Windows, add `--benchmark-baseline PATH --baseline-revision FULL_SHA` to
compare a frozen unoptimized Rust worker with `--rust-worker`. Commit the
implementation/harness first and run without concurrent builds/tests. The
controller retains one warmup and ten fresh-process paired samples per cell,
seeded alternating order, full-output equality, construction time, process wall
time and absolute peak working set. All-pair alpha/Rips, shell/filling controls
and an alpha single-pair case retain adverse/no-benefit outcomes. It is a bounded
resource snapshot; process memory includes input/runtime/output and is not an
allocation peak. The raw manifest records both commits, binary/source/input
hashes and the uncontrolled host-load limitation.
Add `--benchmark-direct` with the same current worker/revision as both inputs
to compare exhaustive native Rust Schur against the unoptimized new-chain batch.
That comparison checks full matrices/spectra numerically and matches ordered
coordinates/scales exactly; representative indices and solve counts differ by
design. It measures whether saved solves repay preparation for these inputs,
without selecting an automatic default.
The separate [SPSC-A construction and LSMR reproduction](../docs/research/persistent-constructions.md)
freezes shared inputs, exact oracles and intermediate paper quantities for
Issue #144. Run `reproduce_persistent_constructions.py` using that page's pinned
sources and a new output directory. It does not replace the existing core suite.
Its `--replay-lsmr-worker PATH --frozen DIRECTORY` mode performs Round-2 native
LSMR parity using hash-verified original Python outputs. The opted-in worker
also checks terminal-breakdown discrepancies under the explicitly scoped policy.

[compare_persistent_laplacians.py](compare_persistent_laplacians.py) compares the
release [example worker](../examples/persistent_laplacian.rs) with independently
assembled signed boundaries, NumPy Schur and SVD nullspace computations, exact
rational restricted-domain matrices and real persistent ranks. NumPy is required;
`--gudhi` adds finite-field corroboration on declared torsion-free examples only.
These dependencies remain outside the Rust runtime.

The local reproduction used Python 3.10.11, NumPy 2.2.6, SciPy 1.15.3 and
GUDHI 3.12.0. Create an isolated oracle environment for the same package versions:

```sh
python3 -m venv target/persistent-oracle-venv
target/persistent-oracle-venv/bin/python -m pip install numpy==2.2.6 scipy==1.15.3 gudhi==3.12.0
cargo build --locked --release --example persistent_laplacian
target/persistent-oracle-venv/bin/python tools/compare_persistent_laplacians.py --worker target/release/examples/persistent_laplacian --output target/persistent-run-001 --gudhi --benchmark
```

On Windows use `target/persistent-oracle-venv/Scripts/python.exe` and append
`.exe` to the worker path. Outputs must be new. Fixtures, raw worker results,
matrix/spectral errors and adverse timings remain under the output directory.
`--reference-python PATH` evaluates only the fixed external Python
`schur_complement` function, without importing its application dependencies.
The optional benchmark measures matrix construction for identical start/target
pairs, excluding source preparation and spectrum computation. Process elapsed
time and Windows peak working set are recorded separately; the latter includes
the process baseline and is not allocation peak. Linux peak measurement needs
an isolated resource harness. No timing threshold is an ordinary unit test.

Add `--components` to verify the independent upper/lower matrices, actions and
spectra, as well as fixed-endpoint upper results. The oracle checks U+D against
the full operator, UD=DU=0, and the nonzero spectral union with zero multiplicity
determined by independent exact rational ranks. It does not substitute numerical
threshold ranks for these exact fixture ranks. The recorded `component_pairs`
and per-pair component errors distinguish this coverage from full-only replay.

`--reference-matlab DIR` verifies the Git blobs of the author's fixed
`ebff6e3ad6f5eb90eb8058e4f5567755b62f0c3e` sources, 35-point input and
precomputed q=2 boundaries. It recovers the original 580 oriented K-triangle
labels, matches Rust coordinates, replays the least-norm Schur formula with
SciPy, compares both routes and eigenvalues, and checks the rebuilt skeleton
with GUDHI. It needs SciPy/GUDHI and an explicit worker storage override already
present in the example. This is author-data/formula replay, not execution of
MATLAB or Dionysus. Original `.m`, `Point.txt` and `.mat` files stay under ignored
artifact storage; no dependency or reference code enters the Rust runtime.

To prepare the optional author-data replay, keep a pinned reference checkout
under `target/`. Its root must contain `SchurComp.m`, `persistLap.m`,
`boundary.py`, `VR_persistLap.m`, `Point.txt`, `q1Boundary.mat` and `qBoundary.mat`;
the controller rejects any of these files whose Git blob identity differs.

```sh
git clone https://github.com/ndag/Persistent-Laplacian.git target/persistent-matlab-reference
git -C target/persistent-matlab-reference checkout --detach ebff6e3ad6f5eb90eb8058e4f5567755b62f0c3e
target/persistent-oracle-venv/bin/python tools/compare_persistent_laplacians.py --worker target/release/examples/persistent_laplacian --output target/persistent-author-run-001 --gudhi --reference-matlab target/persistent-matlab-reference
```

The optional Python subset can be downloaded from the
[fixed MIT reference](https://github.com/tomogwen/persistentlaplaciandatascience/blob/b36b080f739653d84284f3f7ca9ccea5fad286ae/src/perslap.py)
to `target/persistent-python-reference.py`, then supplied with
`--reference-python target/persistent-python-reference.py`. This runs only
`schur_complement`; it does not validate the reference application's pair,
feature or learning pipeline. Outputs also freeze the actual worker and controller,
preserving the worker's POSIX executable permissions, and record the worker hash;
preserve these together with the result JSON. The reference CI job runs a POSIX
regression that launches the frozen worker and checks retained failure output.

## Native cross-library benchmarks

Use [benchmark_native.py](benchmark_native.py) for Cocycle, GUDHI C++ and upstream
Ripser C++. Python's standard library orchestrates separate native executables;
no Python TDA binding is imported into the measured computation. Source setup,
commands and artifacts are documented in the [native guide](../benches/native/README.md),
with execution rules in the [H0/H1 protocol](../benches/protocol.md). The separate
[Rips pipeline suite](#rips-workflow-resources) measures broader public workflows.
These suites and the [distance suite](#diagram-distance-correctness-and-ablations)
follow the shared [reporting rules](../benches/reporting.md).

[build_native.py](build_native.py) verifies pinned sources and compiles workers.
[benchmark_inputs.py](benchmark_inputs.py) owns deterministic fixture generation,
shared with historical controllers without sharing their Python worker path.

## Critical-set performance and persistence context

[benchmark_critical_sets.py](benchmark_critical_sets.py) runs native sparse,
pinned dense and Oineus critical-set queries on Linux/WSL. It also measures
concrete/generic persistence, requested cycles, the internal reference reducer,
and the implicit flag path on compatible inputs. Different requested outputs
remain separate comparison scopes. The [critical-set guide](../benches/optimization/README.md)
owns baseline availability, setup, deterministic fixtures, timing and validation.

## Diagram-distance correctness and ablations

[compare_distances.py](compare_distances.py) checks bottleneck/L-infinity,
W1/L-infinity and W2/Euclidean against independent rational small-diagram
matching, pinned Topp, repaired native GUDHI bottleneck and isolated GUDHI/POT
weighted workers. The default supported suite is full; `--quick` selects a
smaller CI check, and `--suite stress` retains numerical reference disagreements.

[benchmark_distances.py](benchmark_distances.py) compares native Rust/Topp
single-pair calls and controlled memory/search variants on Linux.
[distance_common.py](distance_common.py) owns their f64 transport, pins/builds,
independent oracle and process protocol. The Rust worker compiles private
instrumentation with `--cfg cocycle_distance_bench`; its public-library reference
uses the separately built ordinary crate. No experiment feature is added to the
public API. The [distance guide](../benches/distances/README.md) owns setup,
commands, measurement boundaries, exclusions and artifact interpretation.

## Native Rips correctness checks

[compare_rips.py](compare_rips.py) compiles the workers under `tools/reference/`
and compares exact threshold edges and prime-field interval multisets through H4. GUDHI
also checks complete expanded simplex/value sets against Rust and independent
subset enumeration; Ripser contributes implicit persistence, not explicit complex
export. Dense Rust fixtures compare matrix, sparse and explicit paths. All Rust
fixtures also check that representative-enabled computation preserves the diagram.
Inputs cover
matrices, supplied graphs, isolates, zero edges, ties, cutoffs and random nonmetric
weights. Rust dense/threshold paths are also compared with each other. No Python
TDA wrapper is used.

```sh
python3 tools/compare_rips.py --output target/rips-reference
```

Supply `--boost-include target/native-sources/boost/usr/include` for the local Boost
setup. Output must be new. This shares [source pins](../benches/native/sources.json)
and verification/output instrumentation helpers with native benchmarks, but is a
correctness suite, not timing evidence. Artifacts include fixtures, source and
consumed-header hashes, binary hashes, build logs, outputs and a completion summary.
Failed comparisons leave their outputs in `results.json` and return nonzero.
Ripser's tiny-input adapter exclusions and float64-only fixtures are explicit.
Ripser is built with `USE_COEFFICIENTS`; its pinned 8-bit coefficient storage
supports primes through 251. GUDHI's pinned `Field_Zp` caps primes at 46337.
Larger-field cases record per-reference exclusions and are not counted as native
agreement. Rust still validates those fixtures against their structural contracts.
Field-sensitive fixtures include the flag subdivision of RP2. Representative
closure/basis checks live in Rust tests, not in these C++ output adapters.

The text header is `mode n max_homology_dimension cutoff count characteristic`;
readers accept older five-field headers as F2. Field values are never narrowed
silently. Actual executable protocol checks reject invalid/composite/overflowing
moduli, verify five-field compatibility and exercise reference limits before
comparing mathematical fixtures. Outputs and fixture records include the actual
characteristic; `environment.json` retains the protocol-check outcomes.
Its raw unpaired endpoints are compared along with separately checked Rust coverage.

## Real Dirac correctness checks

[check_dirac.py](check_dirac.py) compares the fixed public-API Rust fixtures
with independently assembled signed matrices, analytic spectra and boundary
SVD using NumPy. Optional fixed TopoNetX methods and GUDHI F3 Betti checks
retain their source identity and shape/field exclusions. See the
[Dirac guide](../docs/guides/dirac.md#independent-reproduction-and-observed-results)
for the standalone worker and commands. Keep output and external sources under
ignored `target/`; this is a correctness check, not a timing benchmark.

## Documentation and tool verification

[check_algorithm.py](check_algorithm.py) provides focused `diagram-analysis`,
`diagram-distances`, `complex-construction` and `persistence-reduction` checks for
algorithm contributors. See the
[command and scope](../CONTRIBUTING.md#focused-algorithm-checks) and
[walkthroughs](../docs/development/algorithm-contributions.md). It uses the local Rust
toolchain and Python standard library, stops on failure, and does not invoke
native comparisons. These developer tools are available in source checkouts;
crate packages include the tutorial and runnable Rust example, not `tools/`.

```sh
python3 tools/check_source.py
python3 tools/check_artifacts.py
python3 tools/check_docs.py
python3 -m unittest discover -s tools -p 'test_*.py'
rustfmt --edition 2024 --check tools/diagram_dump.rs tools/benchmark_driver.rs benches/native/cocycle.rs tools/reference/rips_cocycle.rs tools/reference/sparse_cocycle.rs benches/pipeline/cocycle.rs benches/distances/cocycle.rs tools/reference/dirac_cocycle.rs
```

Run the artifact check after staging. It rejects generated paths and log/archive/
binary files in the Git index, including historical or forcibly added outputs.

Source checks cover maintained text encoding/whitespace and parse Python without
importing tools. They exclude raw evidence, downloaded sources, and build output.
Documentation checks recurse through `docs/`, benchmark reports, and native setup
pages while excluding raw results. They validate inline and reference local links
and headings, flag duplicate/missing explicit reference labels, and scan for CJK
text. External URLs and English prose quality need manual review. Supported
Markdown conventions are in the [code conventions](../docs/development/conventions.md#tests-and-documentation).
Tool unit tests need only the Python standard
library. Native compilation and diagram comparisons require the C++ prerequisites
listed in the native guide. No tests assert machine-dependent speed thresholds.

## Algorithm diagnostics

[profile_rips.py](profile_rips.py) invokes private H1 optimization stages and
checks each against the independent explicit reducer. It varies single/two-pass
initialization independently of virtual apparent-pair reconstruction; the
[testing guide](../docs/development/testing.md#independent-invariants-and-properties)
lists the stages. Use it to inspect algorithm work, not to rank public API
performance.

[profile_scaling.py](profile_scaling.py) reads the historical scaling artifact
schema and requires matching source hashes and validated diagrams. It is not an
adapter for `cocycle-native-v1` results. Its counts are entries and operations,
not allocated bytes. Detailed commands and counter definitions are retained in
[diagnostic instructions](legacy-benchmarks.md#private-h1-optimization-profiling).

## Optional wrapper tools

The optional [persistent spectral reproduction](reproduce_persistent_spectral_features.py)
implements ICML 2023's smallest-eigenvalue truncation, zero padding and flattening
on real, unweighted pair spectra. It retains sample/source identities, bases,
scale pairs, numerical origin, full/partial coverage and a padding mask. Its tiny
Python Schur oracle, independent nullspace/rational-rank checks and GUDHI cases
are research tools; production Persistent Laplacian and result APIs remain
separate work. The [Rust baseline worker](reference/persistent_spectral_baseline.rs)
calls existing Betti/lifetime descriptors; external scikit-learn training stays
outside the crate. No dependency is added to Cargo.

Use Python 3.10+ and an optional environment with the tested versions
`numpy==2.2.6 scipy==1.15.3 gudhi==3.12.0 scikit-learn==1.7.2 threadpoolctl==3.6.0`.
Save the unmodified [pinned author source](https://github.com/tomogwen/persistentlaplaciandatascience/blob/b36b080f739653d84284f3f7ca9ccea5fad286ae/src/perslap.py)
as `target/persistent-spectral-reference/perslap.py` (create that directory first).
The command checks its Git blob before executing only the original matrix-level
functions with explicit boundary/prefix adapters, plus the unchanged feature
loop with trace callbacks. It does not run the author's image constructor or
full MNIST pipeline. Rust 1.91+ is needed for the existing-descriptor worker.
The small vertex-cubical fixture checks matrix conventions; the complete cubical
image-processing pipeline and the paper's MNIST/MoleculeNet experiments have
not been reproduced here.

```sh
python3 -m unittest discover -s tools -p 'test_reproduce_persistent_spectral_features.py'
python3 tools/reproduce_persistent_spectral_features.py --author-source target/persistent-spectral-reference/perslap.py --output target/persistent-spectral-run-001
```

The default inclusive grid is `[0, 1/3, 2/3, 1]`, with degrees zero and one and
eight retained eigenvalues per pair. Columns follow
`(degree, start_index, end_index, eigenvalue_rank)`, including diagonal pairs;
four points give 160 columns. Each spectrum is sorted ascending and truncated;
only a known exhausted domain is padded. The separate mask distinguishes actual
zero eigenvalues from padding. This uniform grid is an explicit convention:
the pinned author's four-point loop uses `[0, 0.25, 0.5, 1]`.

Outputs include 160-column spectral features, padding masks, raw spectra and
matrices, column definitions, group-held-out predictions, dependency/source
identities and file hashes. The output directory must be new. Tree controls add
BFS distance, grounded-solve effective resistance and degree-zero/one ablations.
The `higher-order/` experiment uses graph-isomorphic filtrations with one-face
versus three-face cycle fillings and matched disjoint nuisance trees. Its
degree-one persistent spectra distinguish pairs with equal PH, ordinary spectra
and degree-zero persistent spectra. These constructed examples establish
additional information; they do not reproduce MNIST accuracy or establish
general predictive improvement. The author feature loop is audited with trace
callbacks, with its sampling and aggregate failure retained separately.
Hash-verified `SOURCE_SNAPSHOT.json` also supports offline review-bundle replay
without a Git checkout. See the [paper](https://proceedings.mlr.press/v202/davies23c.html)
and [Issue #111](https://github.com/Aequiludium/cocycle-rs/issues/111) for scope.
Default standard-library tool tests skip the optional numerical checks if their
packages are absent; the reproduction command requires them and fails visibly.

The fixed seed is `20261012`. Each experiment has 64 matched pairs and splits
by whole pair groups into 96 training and 32 test rows. The forest uses 100
estimators, maximum depth four and `n_jobs=1`; IDs, labels and masks are not model
inputs. The tree labels are defined by early-vertex distance, so shortest path
and effective resistance supply the same separating information as the
degree-zero persistent spectrum. In the disk control, the one-face and
three-face fillings give degree-one spectra `{3,3,3}` and `{1,3,3}`, while
degree-zero persistent spectra and positive-length PH agree. The shared disk
motif occurs in both splits, so the experiment does not test new motif families.

Representation limits of this reproduction remain explicit:

- Sorting retains eigenvalue multiplicities and removes dependence on input
  spectrum order. Repeated eigenvalues are not deduplicated. For symmetric
  matrices of the same size, ordered eigenvalues satisfy Weyl's perturbation
  bound; sorting alone is not a reason to claim instability. Truncation still
  discards the remaining spectrum, including all positive values if the nullity
  reaches the retained length.
- Fixed feature width does not make different complex sizes equivalent.
  Domain sizes, coverage and bases are exported, but the model does not consume
  them. Neither the spectra nor their traces are automatically size-normalized.
- The mask separates padding from computed entries after numerical thresholding.
  It does not distinguish exact zero eigenvalues from small values mapped to zero
  by `zero_tolerance`. Raw floating-point spectra and the threshold are retained;
  exact-rank checks certify nullity only for the worked oracle cases.
- Only value features enter the current forests. Values-only, values-plus-mask,
  mask-only and domain-size controls on the same grouped splits remain future
  ablations, especially for varying complex sizes. No robustness guarantee under
  filtration perturbations has been established; threshold crossings, changes in
  sampled stages, domain dimension and retained coverage require separate checks.

For subsequent representation work, see Jung, Kang and Park,
[Stable Vectorization of Persistent Laplacians via Spectral Descriptors](https://arxiv.org/abs/2512.05463v2)
(v2, 11 August 2026). It constructs scalar-signature persistent Laplacian
diagrams (PLDs) and smoothed persistent Laplacian images (PLIs), and reports
trace-based improvements in its MNIST/QM7 experiments. Its stability statement
is at the PLD Wasserstein-distance level under signature admissibility conditions.
Before adoption, check the data-class/size bounds, how the matching metric
controls descriptor-amplitude changes, and the link from input/filtration
perturbations to that metric. This is not a stability guarantee for this tool.
Trace requires the full operator diagonal or full spectrum; the sum of a
truncated smallest-eigenvalue prefix is not the trace. Size normalization, empty
domains and partial-coverage policies need explicit definitions and validation.
These descriptor/PLI experiments are follow-up work; the present algorithm and
reported model inputs remain the ICML 2023 reproduction baseline.

Validation combines hand-derived matrices, restricted-boundary SVD nullspaces,
exact rational incidence ranks and GUDHI F2 corroboration on torsion-free
examples. The main run checks 160 pair matrices; the focused tests add 400
seeded pair comparisons. Schur `rcond=1e-12`, zero tolerance `1e-10` and absolute
matrix tolerance `1e-9` apply to these small unweighted examples. The oracle
computes spectra only in degrees zero and one; the tetrahedron boundary checks
H2 topology without computing its degree-two spectrum. Production spectral
integration remains separate under issues #124 and #109; this tool does not
complete all #111 engineering acceptance criteria.

Recorded timings describe a single traced workflow. Spectrum-generation timers
exclude vectorization; Rust baseline time includes compilation, GUDHI interval
export and worker processes. Parent CPU excludes child CPU, and Python
`tracemalloc` peak is not process-tree peak RSS. The pre-import BLAS inventory
may be empty and does not establish numerical single-thread execution; only
the forest's `n_jobs=1` setting is certain. Treat these measurements as bounded
workflow observations, not kernel speedups or production scaling evidence.

The [legacy tools guide](legacy-benchmarks.md) documents `compare_ripser.py`,
`benchmark_gudhi.py`, `benchmark_scaling.py`, their pinned Python environments,
and optional correctness checks. Generate fresh inputs under `target/`; historical
outputs are not shipped. Current native tools include `benchmark_native.py`, `compare_rips.py`,
`compare_sparse_rips.py` and `benchmark_rips_pipeline.py`, with separate
correctness and timing responsibilities. Do not mix old wrapper samples into
native rankings.

## Non-branching weak reduction

The [non-branching reproduction](../docs/research/nonbranching.md) separately
checks Botnan-Dong Algorithms 1-3, induced kernel metrics, pinned author cubical
code and standalone Rust baseline/compact parity. Its real signed algebraic
controls and GUDHI topology checks do not add a general Rips solver. Use the
documented fixed paper/source inputs and a new output directory; literal-paper
failures and corrected implementations remain separately labeled.

`reproduce_nonbranching.py --library --frozen-round1 PRIOR_REPLAY` additionally
replays the actual production kernel source and budget against the prior
R/E/V/Z, index maps, metric and upper matrices. The ordinary reference checker
accepts `--non-branching --gudhi` to verify both public constructors, action,
complete spectra and metadata against rational projection and generic reverse
elimination. Its explicit exclusion manifest records branching outside rows;
those pairs are not counted as successful comparisons. The example's
`--non-branching` mode keeps the default 65,536-entry cap.

## Native sparse Rips checks

```sh
python3 tools/compare_sparse_rips.py --output target/sparse-rips-reference
```

Uses the same pinned GUDHI source checkout and optional `--boost-include` argument
as exact comparisons. The output directory must be new. No Python TDA package is
used: Python builds processes, writes metric fixtures and checks their outputs.

The adapter generates a copy of the pinned sparse header under the output build
directory, replacing exactly one constructor sampling call. A separate exhaustive
C++ sampler fixes the start and smallest-ID tie policy without reading Rust
outputs. On unique positive-radius choices it also checks the unmodified upstream
metric sampler. GUDHI sparse edges, higher-simplex blockers and persistent
cohomology remain unchanged. For explicit dimension zero the adapter prunes the
initial graph to vertices to match Cocycle's documented expansion contract.

The suite compares full simplex/value sets and interval multisets over shared
fields, as well as deterministic sampling and an independent subset enumeration.
It retains fixture, source, dependency/header, compiler and binary hashes, plus
instrumentation details. `environment.json`, `results.json`, `summary.json` and
`build/build.log` are correctness artifacts, not timing/memory measurements.
Upstream Ripser is explicitly excluded because it lacks this sparse approximation
constructor; feeding only the modified graph to its flag engine would omit the
higher-simplex blockers. Continue running the exact native suite separately.

## Rips workflow resources

`benchmark_rips_pipeline.py` measures validation/conversion, construction,
expansion, public computation and export separately, plus end-to-end time and
per-process peak RSS. See the [pipeline protocol](../benches/pipeline/README.md)
for the full scopes, native comparability limits, failure retention and commands.
The existing H0/H1 timing protocol remains separate. Generated results and logs
are stored outside Git under the shared reporting policy.
