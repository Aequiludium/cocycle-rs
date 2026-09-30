# Rips pipeline measurements

[Benchmarks](../README.md) / Native pipeline

This suite complements the existing [H0/H1 native protocol](../protocol.md) and
[construction correctness workers](../../tools/README.md#native-rips-correctness-checks).
The shared [reporting rules](../reporting.md) govern evidence classification and
claims; use the [report template](../report-template.md) for retained experiments.
The current machine-emitted protocol is `cocycle-rips-pipeline-v2`. Earlier
unversioned, fixed-order measurements must not be pooled with this protocol.

It measures the public APIs added by the complete Rips subsystem. GUDHI and
Ripser execute as native C++ processes; Python only prepares fixtures, builds,
limits processes, checks outputs and summarizes samples.

```sh
python3 tools/benchmark_rips_pipeline.py --quick --samples 1 --output target/rips-pipeline-smoke
python3 tools/benchmark_rips_pipeline.py --samples 12 --output target/rips-pipeline
```

Commit the harness and measured sources before running; a dirty working tree is
rejected. `--kernel-revision <full-sha>` can identify an earlier kernel commit
when its `src/` and Cargo manifests are identical to the current tree. The harness
commit is recorded separately. Each output directory must be new. The pinned source setup and optional
`--boost-include` match [the native setup](../native/README.md). Linux is required
for `/proc/self/status` memory counters and per-child address-space limits.
Library correctness tests and public examples remain portable.

## Workflows and boundaries

One fresh process performs one workflow. One discarded warmup and the requested
measured samples use separate processes. All warmups precede measured rounds.
Each round executes each supported backend once. `--order-seed` (default zero)
selects shuffled cyclic-order blocks; positions balance exactly over each block
of two or three rounds. The default 12 measured rounds balance both group sizes.
Each case adds its index to the seed and records its schedule and sample positions.
`--cpu <id>` pins workers to an allowed Linux CPU; without it they inherit the
controller's affinity. Frequency and competing host load remain uncontrolled.
Measurements run serially; do not run other compilation or benchmark workloads
concurrently. Repetition and balanced order alone do not establish stable rankings.
There are five phase bins:

| Bin | Rust | GUDHI C++ | Upstream Ripser C++ |
| --- | --- | --- | --- |
| Input | Borrowed input validation; upper/square fixture conversion or graph validation/adjacency | Finite-distance scan | Finite-distance scan and float32-exact check |
| Construction | Threshold or approximate graph; streamed Euclidean distances for points; metric checking when selected | Exact graph into tree, or sparse graph with original metric sampler | Dense float32 buffer conversion, or sparse adjacency |
| Expansion | Explicit frozen simplices and bidirectional incidence, when requested | Tree clique/blocker expansion; sparse graph insertion is included here | Not applicable |
| Compute | Public persistence call including owned diagram normalization and requested bases | Coefficient initialization and persistent cohomology | Implicit barcode computation and numeric pair capture |
| Export | Interval JSON payload construction | Pair extraction, sorting and interval JSON payload | Sorting and interval JSON payload |

End-to-end time includes these phases, intervening bookkeeping and destruction
of intermediates released during the workflow. It excludes text fixture I/O,
parsing, the final metrics envelope/transport and destruction of retained results.
Phase sums need not equal the whole workflow. Rust cannot expose a separate
reducer-only timing through its public API; its compute bin includes result
assembly. A GUDHI tree is not the same representation as a Rust frozen complex;
compare complete workflows before interpreting individual bins.

Matrix fixture conversion is explicitly measured, not hidden. Normal Rust users
can borrow an existing upper/square matrix without this conversion. Raw fixture
buffers remain part of process memory; they are not counted as reducer workspace.
Point workflows discard their unused precomputed values before measuring the
public point API, but parsing can already have raised the process high-water mark.

The input geometries include circle, uniform, nonmetric, cross-polytope sphere,
bipartite graphs, isolated vertices, exact point distances and dyadic metric
approximation. Cutoffs, dimensions, fields and representative requests are retained
in each result. Quantization to float32 occurs once when constructing shared
exact fixtures; all three workers receive those exact numbers. Sparse approximation
uses exact dyadic Manhattan metrics without additional rounding.

## Native comparability

GUDHI's sparse constructor uses its original metric sampler, with only the initial
vertex fixed to zero. A read-only accessor exposes the resulting permutation.
Fixtures have unique greedy choices, which the controller independently checks;
it rejects unequal sampling or diagrams. This performance adapter does **not**
use the exhaustive reference sampler from the construction correctness suite.
That separate suite still covers ties, duplicate points and full simplex sets.

Ripser has no sparse approximation constructor with higher-simplex blockers, so
those workers are explicitly excluded. Coefficient limits also produce explicit
exclusions. Native representative computation is not provided by these adapters;
GUDHI/Ripser runs for a representative case check diagram correctness only and
must not be used to claim a comparable basis-computation speedup. The same applies
to native precomputed-distance references for Rust point workflows, native lower
matrices compared with Rust upper/square conversion, Ripser diagram-only runs
compared with requested explicit complexes, and GUDHI comparisons against Rust's
exhaustive metric-checking workflow. Representative payloads are retained and
counted but only intervals are serialized in the measured export phase.

Full interval multisets are compared on every warmup and measured run, preserving
multiplicities and unpaired endpoints. Rust coverage is retained; reference infinity
is interpreted using the fixture's range, never as an unconditional essential class.
Approximate permutations and explicit simplex counts are additionally checked.
This timing harness does not replace independent mathematical or native topology
correctness suites.

## Memory, failures and artifacts

Each worker reads its own `VmRSS` and `VmHWM` before the measured workflow and its
`VmHWM` after computation/export. Peak RSS includes runtime, parser/input storage,
allocator retention, temporaries and outputs. The difference between high-water
marks is **not** allocated bytes or an algorithm-only peak. Do not use the Python
controller's RSS or cumulative child high-water mark as an individual sample.

`--address-space-mib` is a process virtual-address-space cap, not a library RSS
budget. `--timeout` kills that sample process. Failures and timeouts remain in
`results.json`, fail the suite and cannot be silently converted to exclusions.
After a backend fails, its remaining scheduled attempts are retained as `not_run`;
other backends finish their schedules. All statistics for that case are withheld.
The library's cooperative work limits are tested separately and do not claim
hard memory or wall-time enforcement.

`environment.json` records source/pin/header/binary hashes, compiler commands,
full kernel/harness commit identities, protocol ID, invocation, CPU/platform
details, affinity, order seed and limits. `fixtures/` retains exact inputs, `results.json`
all samples and validation outcomes, `measurements.json` phase medians, ranges,
peak RSS and comparison scopes, and `summary.json` the overall outcome.
`build/build.log` preserves compiler diagnostics. A comparison mismatch makes
all measurements of that case unvalidated; retained timings are not successful
performance evidence. Source or commit changes during a run also fail the overall summary; a report
must check that summary before using per-case statistics.

Generated outputs belong under ignored `target/` or in external artifact storage.
They must not be committed as loose files or archives; see the
[storage policy](../reporting.md#storage-and-evidence-lifecycle).

## Private assembly diagnostics

`tools/profile_assembly.py` implements the separate `cocycle-rips-assembly-v1`
diagnostic protocol. It copies a fixed kernel into the output's `build/profile-source`,
adds private clocks and counters there, and compares fresh ordinary and instrumented
Rust workers on identical inputs. Production source and public APIs are unchanged.
Use a clean committed checkout matching `--kernel-revision` for formal evidence:

```sh
python3 tools/profile_assembly.py --kernel-revision <R0-full-SHA> --samples 12 --cpu 0 --output target/assembly-coarse-001
```

The full matrix covers existing dense/point H1 and dense H0, large forest/isolated
H0, triangle-free bipartite H1, generic cohomology and nonzero-born supplied
complexes that require boundary reduction. Independent analytic interval multisets
validate the output-heavy controls; every sample also agrees with an ordinary
worker. One discarded warmup and balanced worker order precede measured samples.
`--quick --samples 1 --exploratory` is only a development smoke run.

The kernel clock covers the actual raw producer, including raw generation and
producer-local cleanup. The assembly clock covers filtering, endpoint conversion,
interval validation, fallible Vec growth, raw-buffer release, diagram validation
and its single canonical sort. Public analysis also includes access construction,
dispatch and cleanup; those costs remain a separate remainder. Both assembly/public
analysis and assembly/(kernel + assembly) sample fractions are retained. Export,
JSON transport, descriptors, distances and retained-result destruction are excluded.
Median phase values need not sum to a median total.

Only if coarse evidence warrants it, `--decompose --cases <case-substrings>` adds
materialization, diagram-validation and sorting clocks and counts actual Vec
capacity changes. Materialization groups filtering, conversion, interval-local
validation, growth and raw release; per-interval clocks would distort short loops.
Capacity changes are not allocator-call counts. Stage A and B are separate timing
pools; compare instrumented public analysis with ordinary controls for observer effects.

Runtime `size_of` and Vec capacities quantify overlapping owned raw/final output
storage during conversion. Raw storage is released before diagram validation and
sorting. This excludes reducer state, allocator metadata and transient reallocation
storage; it is neither whole live-memory peak nor measured RSS savings. Process RSS
is supplementary and includes fixture parsing, workspace and allocator retention.

Each attempt uses a fresh output directory. `environment.json`, source and binary
hashes, generated sources, build commands/logs, fixtures and every raw process record
are retained. `results.json` indexes raw files by hash; `measurements.json` and
`matrix.csv` contain validated statistics. Missing hooks, unexpected route/call
counts, mismatched diagrams, inconsistent capacities or impossible timing bounds
fail the run. Check `summary.json` before using any statistics.
