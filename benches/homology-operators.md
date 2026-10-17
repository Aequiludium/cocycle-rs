# F2 operator verification and measurement protocol

[Benchmarks](README.md) / `cocycle-f2-operator-v1` / [Reporting rules](reporting.md)

This protocol concerns feasible weighted F2 projections, separately from native
persistence/diagram-distance rankings. The reference is homology-operator
`a11bf50b8c8320651e9ac015566778d957af9c17`, with MIT attribution in the
[third-party notices](../THIRD_PARTY_NOTICES.md). Use a clean reference checkout
at that commit. It is external infrastructure, never a Rust runtime dependency.

## Build and check

```sh
cargo build --locked --release
rustc --edition 2024 -C opt-level=3 tools/reference/homology_operator.rs \
  --extern cocycle=target/release/libcocycle.rlib \
  -L dependency=target/release/deps -o target/homology-worker
python3 tools/check_homology_operators.py --reference /path/to/pinned-reference \
  --worker target/homology-worker --oracle-python /path/to/gudhi-python \
  --geometry --output target/operator-correctness
```

On Windows the worker filename ends in `.exe`. `--native-dir` may identify a
directory containing the pinned upstream's built `_homology_native` extension,
built with its own locked manifest and matching Python ABI. The pinned reference
requires Rust 1.98.1, PyO3 0.29.3 and maturin 1.15.0 for that optional extension;
these are comparison infrastructure, separate from Cocycle's Rust 1.91 MSRV. This adds an actual
`NativeFactorizedSolver` comparison; missing native support fails visibly and
never falls back to Python. Retain its binary hash and toolchain versions.

The controller verifies the reference SHA/clean tracked files, fingerprints
source/worker/reference/native binary, records commands/environment, and
preserves complete worker stdout/stderr/status and shared fixture bytes. Each
run requires a fresh output directory. Generated data must never enter Git.

## Correctness domain

50 reproducible fixtures cover exhaustive binary inputs in dimensions 0-5,
general and homology-only windows, singular matrices, empty shapes, arbitrary
noncycle actions, positive integer weights and sizes 63/64/65/127/128/129 plus sparse 511/512/513/1025 cases.
Compare Betti numbers, every explicit P entry, all requested P actions and all
cycle selected masses exactly. A state other than `Computed` cannot pass.
The independent tiny oracle enumerates the actual cycle/boundary sets and
quotient cardinality, without production matrix elimination or expected P.

GUDHI F2 verifies explicit triangle stages and H1 birth/death; supplied simplices
and filtration are checked after insertion to reject silent repair. A finite
family additionally has hand-derived transport ranks and half-open stage bars.
Public Rust tests check imported projections, identity mixing, cancellation,
work exhaustion, overflow, basis permutations and variable/inherited costs.
The example exercises the actual `FilteredComplex` extraction path.

`--geometry` additionally reads the pinned `geometry_queries.json` T1--T5
corpus, verifies its semantic digest, and gives every implementation Cocycle's
same frozen simplex order. It compares full historical generators, endpoints,
occurrences, source vertices, exact class minima, H1 unit birth cycles, fixed
fillings and H1/H2 optimal/stable volumes. Upper costs include both unit counts
and deterministic positive integer weights. Independent raw-simplex subset
enumeration checks every filling and volume feasible set and minimum with
epsilon=1,2,3. Raw worker output and corpus byte/semantic hashes are retained.
The weighted queries execute the Python reference, not a new native geometry
kernel. The optional native extension still compares feasible P/actions.

The geometry comparison imports the pinned source's test builders. They are
verification infrastructure only. The source provides exact integer/Fraction
scales; Rust accepts the shared exact integer subset for stable volume. Neither
this comparison nor GUDHI establishes a convex-fan perturbation guarantee.

GUDHI has no comparable weighted projection/action oracle for these general
chain windows. There are explicit exclusions for globally minimum-stretch projections, arbitrary
precision/rational/float costs, Python snapshots, general linear maps and the
upstream's constant final extension. Ripser/Topp are not operator solvers and
are inapplicable to this new domain. Existing persistence and distance kernels
are unchanged.

## Owned family and history resources

`tools/benchmark_operator_family.py` implements `cocycle-f2-family-resources-v1`
on Linux/WSL. It compares two clean, committed Rust source snapshots with the
same standalone worker source. Build the library at the candidate and baseline
roots with `cargo build --locked --release` and compile
`tools/reference/homology_operator.rs` twice using `rustc --edition 2024
-C opt-level=3`, linking each root's release rlib/dependency directory. Keep the
executables separate and pass their actual paths:

```sh
python3 tools/benchmark_operator_family.py \
  --worker /path/to/candidate-worker \
  --baseline-worker /path/to/baseline-worker \
  --baseline-root /path/to/clean-baseline \
  --samples 10 --output target/benchmarks/commit-SHA/operator-family/run-001
```

The default supplied complexes are triangulated disk fans with 4/8/16/32/64
vertices (11--251 simplices) and complete graphs with 4/6/8/12/16 vertices
(10--136 simplices). Values are 0/1/2 for vertices/edges/faces, and same-value
insertions use Cocycle's face-first order. F2 H1 and unit costs are fixed.
These two families exercise finite deaths and many simultaneous surviving
generators; they do not represent every complex or establish large-TDA admission.

The `construct` mode times `from_simplices` from an already validated source,
including extraction, every stage's solve/validation, family checks and owned
result retention. The `features` mode prepares that owned family before timing,
then measures complete feature construction, backward correction, owned history
output and independent validation. Neither timer includes fixture parsing,
source validation, process startup or JSON export. Both keep the result alive
at the last memory reading. Features' prepared-family cost is not subtracted
from construction medians or described as free.

Each cell uses one retained/discarded warmup per engine and ten fresh serial
processes per engine, with seeded position-balanced order. Do not run other
builds/tests concurrently with the measurement. Default process limits are 30
seconds and 512 MiB address space; these are distinct from cooperative library
controls. CPU affinity, environment/compiler, source/worker/harness hashes,
fixture bytes, warmups and every measured stdout/stderr are retained. Host load
and CPU frequency are uncontrolled. Report median and min/max, retaining
regressions, timeouts and failures. Fewer than ten samples are resource snapshots.

Linux `/proc/self/status` supplies pre-call RSS/HWM, post-call retained RSS and
absolute HWM before export. Peaks include source, process runtime, prepared
families, workspaces and retained results; they are not live algorithm bytes.
In particular, prepared-family HWM can dominate a feature call, so zero HWM
growth does not mean zero allocations. Windows RSS fields are unavailable and
must not be reported as zero.

Every measured stage Betti sequence and feature interval list is checked against
an independent raw-simplex PH reducer. Historical chains are checked for cycles,
adjacent homology preservation and death in raw boundary spans. Candidate and
baseline export all histories for exact comparison. The oracle never supplies
library results. The [report](reports/homology-operators.md) identifies the
measured commits and their actual envelope, separately from historical solve
measurements and new correctness tests.

## Timing scope

Add `--baseline-worker target/homology-before --samples 12 --sizes 8,32,65` for serial independent processes with one
retained/discarded warmup per engine/cell and seeded shuffled engine order.
Generate one fixture per size/family and share its exact bytes. Timing measures
one fresh public feasible result from a prevalidated window, including legality
checks and result construction; fixture parsing, input validation, startup,
explicit P export, serialization and topology construction are outside time.
Rust's owned window clone is inside its measured solve. The optional baseline
worker executes the same solve boundary and every measured output is checked
against both versions. General, homology-only and sparse windows share fixtures.

Python reference solves through `FeasibleSolver`, or actual upstream
`NativeFactorizedSolver` with `matrix_free_output=True`. Both perform the
reference dispatch/independent validation and construct `HomologyOperator` with
Python identities and provenance. The Rust result has its documented smaller
native context/run contract. These measurements compare complete public
workflows delivering the same feasible P, not isolated equivalent-language
matrix kernels or complete feature parity. Never attribute a ratio solely to
language, packing, projection layout or the removal of Python bindings.

The worker's `query` mode is available for focused diagnosis, but formal
controller measurements here are solve-only. No batch/streaming/default
optimization admission follows from those solve numbers.

Retain all samples, warmups, output mismatches, timeouts and failures. Report
per-cell median and range, with source/build/environment identities. Linux
`VmRSS`/`VmHWM` record baseline RSS/HWM and absolute post-call peak before export;
Python and Rust process baselines differ. Windows memory fields are unavailable,
not zero. Do not pool Windows/Rust 1.98 and Linux/Rust 1.91 measurements or call
absolute process RSS an algorithm allocation peak. Follow the
[reporting rules](reporting.md). For a local uncommitted run, retain exact source
bytes, hashes and dirty/integration state and label the results accordingly; a
base commit alone does not identify modified code. The 2026-10-09 report remains
a historical comparison against its original reference pin.

The standalone worker also accepts family-bench n stages iterations for a
fixed H0 path-merge sequence with repeated scales. It measures only barcode()
on a previously accepted family; operator/input construction, startup and
serialization are outside the timer. Compare before/after workers compiled from
the same worker source, keep 12 fresh pairs plus discarded warmups, and check
the hand-derived interval multiplicities before interpreting timing. This is a
separate boundary from feasible solve and is not a full PH workflow comparison.

## Hosted correctness evidence

The existing CI workflow's `F2 operator reference comparisons` job builds the
pinned actual native extension and compares 50 cases against Python and native,
including the independent tiny-span and GUDHI F2 checks. It uses `--samples 0`,
so its outputs are correctness evidence, not a performance study. Failed worker
startup or solver availability is retained and cannot be fallback success.
Artifacts are named `operator-correctness-<sha>-<attempt>` and retained for 14 days.
Use the exact run/head linked from the PR; an expired artifact or unrun job is
not a passed check or durable raw-performance archive.
