# Developer tools

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

The existing Hodge checker also accepts `--path-convergence`: it records
independent SciPy identity, inverse-diagonal and shifted-tridiagonal-inverse
controls on the 1024-node path against its complete analytic spectrum, using the
same SplitMix64 starting span. Warnings and incomplete lookahead convergence
remain explicit; these controls are numerical diagnosis, not native timing.

## Spectral reference CI

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
