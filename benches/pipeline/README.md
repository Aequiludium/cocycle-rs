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

## Dimension-boundary workspace ablations

After freezing the private H2 prototype, M1 uses the same pipeline v2 worker and
fixture convention with a separate controller:

```sh
python3 tools/benchmark_workspace.py --native-environment target/rips-pipeline/environment.json --samples 15 --cpu 0 --output target/workspace-run-001
```

The native environment identifies hash-checked, pinned GUDHI/Ripser workers from
an earlier completed pipeline build. The controller generates three source-bound
Rust variants under its fresh output directory: joint baseline, explicit release
of generic transformation columns before next-level assembly, and reuse of the
two identically typed H2 heaps within a dimension. It refuses source-hook drift.
These are controlled experiments, not ordinary-library selectors or production
admission. The H1-to-H2 release route remains identical where existing scopes
already release the reducer. H3 controls audit the existing generic continuation;
they do not implement T9 or M2.

All builds finish before serial measurement begins. Executables reside on Linux
local storage; fixtures, source/build/binary hashes and every warmup/measured
sample remain in the output. Complete interval multisets are checked against both
native references; a failure withholds rankings for the entire case.

Test-only `workspace_event` traces run in separate processes after each latency
schedule. They observe boundary timestamps, Vec/heap capacity bytes including
nested vertex/transform payloads, hash-table logical entries/slots, BTreeMap
entries, sampled VmRSS and VmHWM. They do not measure allocator metadata, exact
live allocator bytes or allocator calls; those unavailable fields remain null.
Trace timings and RSS include the diagnostic work and must not be pooled with
ordinary worker latency or interpreted as process-memory equations. Boundary
failure tests simulate recoverable allocation errors and cancellation; they do
not claim universal recovery from real allocator exhaustion.

The [M1 report](../reports/h2-workspace.md) records the completed ownership audit,
two independent runs, capacity retention and failed control screens. The ablations
remain generated experiments; ordinary library dispatch is unchanged. Executables
are retained in each run's `binaries/` directory before measurement; repeat with
`--reuse-build` pointing to the first build-bearing run, using a new output and
order seed. Temporary execution paths are not the retained evidence store.
Capacity bytes describe test-build objects, including H1's empty test-only
reduced-column Vec headers; they do not measure ordinary-worker live allocation.

## Private H2 baseline

For the T3 fixed-tuple experiment, build the same worker/controller with
`RUSTFLAGS="--cfg cocycle_h2_bench"`. This private configuration changes only
exact F2 diagram-only maximum-dimension-2 dispatch. Normal builds retain generic
H2 continuation; H1-only, H3+, odd primes, representatives and approximation
blockers keep their established paths. The configuration is not a Cargo feature
or public engine selector. No H2 shortcut is enabled in this baseline.

Commit the clean implementation before running `--phase3` at both sizes; retain
the actual RUSTFLAGS, worker hashes and source/harness identity. Counters and
ownership landmarks use a separate release test executable:

```sh
COCYCLE_H2_FIXTURE=path/to/fixture.txt cargo test --release --lib \
  persistence::flag::cohomology::profiling::profile_h2_prototype \
  -- --ignored --exact --nocapture
```

Run one fixture per fresh counter process. Validate emitted raw intervals against
the frozen pipeline result before interpreting counts. These test builds contain
instrumentation and do not supply ordinary latency. H1 extraction is interleaved
with pairing and finishes before the H1-return event; capacity/event snapshots
do not measure total live allocations or establish RSS release. The mathematical
[baseline contract](../../docs/reference/mathematics.md#experimental-fixed-tuple-h2-baseline)
and [H2 report](../reports/h2-foundation.md) distinguish validation from admission.

The existing profiling controller can also compare two frozen, validated
`--phase3` pipeline runs without compiling during measurement:

```sh
python3 tools/profile_rips.py --h2-baseline path/to/t2-run \
  --h2-candidate path/to/t3-run --samples 96 --cpu 0 --quick \
  --output target/benchmarks/commit-<sha12>/localfs-paired/run-001
```

Omit `--quick` for the n32 runs. It verifies fixture/executable hashes, copies
byte-identical workers to Linux `/tmp`, retains every sample, checks diagrams
and resamples matched rounds for a ratio-of-medians 95% interval (2000 replicates).
The controller commit/hash and both original build environments are recorded.
Keep regressions visible; a passing screen alone does not admit default dispatch.
