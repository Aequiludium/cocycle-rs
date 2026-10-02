# Native diagram-distance experiments

[Benchmarks](../README.md) / [Reporting rules](../reporting.md)

The separate [preparation investigation](../reports/distance-preparation.md)
provides a deterministic public-API counterexample and runtime allocation
controls. Its warmed-call protocol differs from this suite; do not pool samples.

## Concurrent resource study

`tools/benchmark_distance_resources.py` is the separate Linux protocol
`cocycle-distance-resources-v2`. It generates a worker from the existing native
adapter plus `resources.rs`, on immutable Phase-2 R0. Production code is unchanged.
Each warm worker parses inputs and completes one discarded in-process warmup before
the ready barrier. Raw validation/copying, ordinary preparation, solve and cleanup
stay inside each native job clock. Candidate options retain existing diagnostics.
Finite inputs, endpoint precision and the distance definitions match the native
suite. Every job is checked against separately linked ordinary public R0; tiny
inputs additionally use the independent exact oracle.

Process groups use independently pinned workers; thread groups share one raw
input in one process with affinity over the same allowed CPU set. Jobs keep fresh
pair-local state. These are different concurrency contracts and stay separate.
One parent monotonic dispatch-to-completion envelope measures group throughput,
including go/result transport and group scheduling, excluding startup/input and
ready warmups. Native per-job monotonic clocks exclude result transport. Keep
both; this protocol does not pool samples with single-shot or lifecycle studies.
Each cell has a discarded fresh group and 12 measured balanced groups per route.
Per-group p95 is nearest-rank over its native jobs; summaries use group medians.

`--trace` enables a separate 2-ms requested sampling cadence. `/proc` CPU time is
scheduled core-seconds (including stalls while scheduled), not retired-instruction
work or a direct stall measurement. Sampled PSS sums proportionally account shared
pages and are a proxy, not cgroup aggregate memory; RSS sums remain separate.
Sequential procfs reads and observer scheduling limit temporal resolution. Trace
integrals use trapezoids over actual monotonic sample times, including the final
observer boundary. Warm/input/allocator residence is included. Individual process
VmHWM includes premeasurement work and is not a simultaneous group peak.

The current WSL environment has no usable perf/IMC counters and no writable cgroup
root. DRAM bytes and LLC misses are recorded as unavailable (`null`), never zero or
invented cache proxies. Trace wall times do not select routes. The
[completed available-host study](../reports/distance-resources.md) records finite
N95 grids, mixed interference, phase replay, held-out model errors and independent
repeats. Its adoption gate retains R0; the sampled two-resource model does not
identify physical saturation or justify a production resource-admission rule.

V2 adds cold single-call traces (`--cold --jobs 1`), per-child dispatch/completion
and CPU observations, trace node PSS, explicit staggered dispatch, and JSON
heterogeneous process plans (`--plan`). Plans list named groups with `nodes`, each
specifying generator `family`, `size`, `split`, `metric`, `variant` and optionally
`jobs`; a group can specify `stagger_ms`. They use the same numerical contract.
Completed children remain resident until group ack. Phase replay retains that
post-completion residence and subtracts single-process initial PSS before adding
increments to the measured group-ready PSS, limiting shared-page double counting.
Replay remains a sampled proxy, not a hard memory-feasibility certificate.

`--target-ms` chooses jobs from a four-call baseline calibration before any
measured group. It never chooses a route; counts are clamped to 1..4096 and retained.
Counts remain fixed across concurrency within a run. `--order-seed` controls fresh
balanced schedules. Processes exceeding the available CPU count use explicit
round-robin pinning; threads share the selected CPU mask. Oversubscription and
uncontrolled frequency/host load remain part of the measured contract.

`tools/analyze_distance_resources.py` fits tuning-only M0 peak, M1 static CPU/PSS,
M2 linear load-dependent and M3 pair-interference contrasts. It rejects holdout
records at fitting boundaries. Demand and peak budgets are separate. The current
two-resource empirical model cannot identify DRAM or a physical bottleneck.
Frozen model parameters, held-out errors and any no-go decision stay in artifacts.
The earlier first tranche uses V1; its measurements are never pooled with V2.

```sh
python3 tools/benchmark_distance_resources.py --families uniform sparse --sizes 128 --workers 1 --samples 3 --jobs 64 --trace --output target/resources-trace
python3 tools/benchmark_distance_resources.py --worker-dir target/resources-trace/build --families uniform sparse --sizes 128 --workers 1 2 4 --samples 12 --jobs 64 --output target/resources-process
python3 tools/benchmark_distance_resources.py --worker-dir target/resources-trace/build --families uniform --sizes 128 --workers 2 --modes threads --samples 12 --jobs 64 --output target/resources-threads
```

Formal runs require a committed clean harness, exact source/binary fingerprints
and the unchanged R0 kernel. Preserve failures in distinct output directories.

## Persistent lifecycle study

`tools/profile_distance_lifecycle.py` uses protocol
`cocycle-distance-lifecycle-v1`. It generates private copies of immutable R0
`b2c3bd5eebf0c1193f30ea651012b8ae61b274c1`; ordinary R0 is separately linked.
The production source and public API do not change. See the lifecycle report
in the [report index](../reports/README.md) for conclusions after measurement.

Diagrams are constructed before timing. Two initial operations warm generated
code, then all candidate preparation and retained buffers are dropped. The
whole measured sequence includes frequency counting, cache construction,
per-operation semantic validation, finite solve, output allocation and scalar
equivalence checks. Inner clocks are disabled for amortization curves; use
`--detailed` only for separate diagnostic decomposition. Phase timers overlap
with solver time and must not be summed as exclusive components.

R1 prepares the fixed reused operand; R2 prepares every operand occurring more
than once. R3 adds batch-local cleared scratch capacity, R4 retains it between
batches, R5 also retains the batch output vector, and `bounded` drops any scratch
slot exceeding 4096 bytes. This is a per-slot cutoff, not a total memory quota.
Matching, dual state, candidate graphs, residual networks and heaps remain
pair-local. Generation checks prevent newly dropped sparse scratch from being
reused by a later augmentation in the same operation. Thus ordinary sparse
intra-operation policy stays intact.

K is the requested operation count for same-pair, one-to-many, many-to-one and
repeated-batch workloads. Repeated batches have eight operations and distinct
diagrams per batch. All-pairs uses the full directed matrix (including diagonal)
of `max(2, ceil(sqrt(K)))` diagrams; inspect actual operation and operand counts.
Single has K=1 only. Every formal cell has a separately discarded process and
12 balanced measured processes per variant. Existing fresh-process distance
results remain separate.

```sh
python3 tools/profile_distance_lifecycle.py --output target/lifecycle-curves-001 --cpu 0
python3 tools/profile_distance_lifecycle.py --worker-dir target/lifecycle-curves-001/build --samples 1 --detailed --families uniform duplicates --sizes 8 32 128 512 --patterns same_pair --ks 4 --output target/lifecycle-phases-001 --cpu 0
python3 tools/profile_distance_lifecycle.py --worker-dir target/lifecycle-curves-001/build --samples 12 --poison --huge-size 4096 --families dense --sizes 32 --patterns same_pair --ks 1 --metrics bottleneck --variants r2 r4 bounded --output target/lifecycle-poison-001 --cpu 0
```

Poison traces solve reference operations in a separate process, avoiding an
unmeasured huge solve in the measured allocator. They record three small pairs,
one huge pair, then eight small pairs, with RSS snapshots outside micro-clocks.
RSS is read from `/proc` on Linux and `GetProcessMemoryInfo` on Windows while
the worker waits at an output handshake. Retained capacities count owned Vec
payloads, excluding allocator/header/HashMap overhead and temporary in-use
buffers. They are not RSS. `--workers 2` and `--workers 4` release independent
processes at a ready barrier and record throughput plus aggregate endpoint RSS
and sum of individual peak RSS (the latter is not a simultaneous RSS peak).
With `--cpu 0`, Linux uses consecutive CPU IDs for concurrent workers.

Correctness includes ordinary-R0 equality on every operation, an independent
exact tiny oracle for at most eight points, and worker selftests for cached
orientation, scale, coverage, computed-empty/uncomputed dimensions, essential
counts, contexts and valid-failed-valid recovery. Formal runs reject dirty
harnesses, changed worker hashes and kernels differing from R0. Smoke runs may
use `--exploratory`; their timings do not support selection.

Use this suite to check bottleneck, W1 and W2 against independent references,
then compare memory and search choices within Rust and C++. Protocol
`cocycle-distance-v1` preserves f64 inputs and has its own timing boundary;
do not pool its samples with either Rips suite. For library usage, start with
the [distance example](../../examples/diagram_distances.rs) and
[matching contract](../../docs/reference/mathematics.md#16-diagram-matching-distances).

## Prepare the reference environment

Requirements are Rust 1.91+, Python 3.10+, a GCC-compatible C++20 compiler, Git,
CGAL and Boost headers for the independent fixed GUDHI bottleneck oracle.
The pinned NumPy oracle environment requires Python 3.11 or later; the standard-
library controller itself supports Python 3.10.
The benchmark controller requires Linux, including for resource smoke checks.
Windows can run correctness checks with a GCC-compatible toolchain; this builder
does not support MSVC.
The library gains no dependencies. [sources.json](sources.json) pins Topp and the
exact GUDHI/POT/NumPy versions. Use a dedicated clean external source checkout:

```sh
git clone https://github.com/proffitteoy/Topp.git target/native-sources/topp
git -C target/native-sources/topp checkout --detach ffa1da051ca7ac5e313c74cc9fb92a2bcb20c234
git clone https://github.com/GUDHI/gudhi-devel.git target/native-sources/gudhi-distance
git -C target/native-sources/gudhi-distance checkout --detach 4ec34ac55e6d2e8cfd7c322e84c1b0a56d516d51
python3 -m venv target/distance-oracle-venv
target/distance-oracle-venv/bin/python -m pip install gudhi==3.11.0 numpy==2.4.6 POT==0.9.6.post1
```

These commands use Linux virtual-environment paths; on Windows the interpreter
is `target/distance-oracle-venv/Scripts/python.exe`.
`--topp-source`, `--cargo`, `--rustc`, `--cxx`
and `--gudhi-python` accept explicit paths. The builder rejects wrong/dirty Topp
and GUDHI sources. `--gudhi-source`, `--cgal-include` and `--boost-include` locate
the fixed native oracle and non-system headers. The builder reads external
sources and never installs packages, resets, or modifies an external checkout.
It fingerprints source, native binaries, commands and toolchains.

## Check correctness and the harness

Run from the repository root. Every output directory must be new; change the
attempt suffix when repeating a command.

```sh
python3 tools/compare_distances.py --quick --gudhi-python target/distance-oracle-venv/bin/python --output target/distance-correctness-001
```

This checks the small supported suite. Success prints `Distance correctness:
passed (...)` and writes `summary.json` with `status: "passed"`. Omit `--quick`
for the full supported suite. Missing references and numerical disagreements
return a nonzero exit status. Once workers run, their failures and comparison
results are retained in `results.json`; build failures remain in `build/build.log`.

On Linux, check benchmark workers with a small resource snapshot:

```sh
python3 tools/benchmark_distances.py --quick --samples 1 --exploratory --groups baseline arena --families uniform sparse --gudhi-python target/distance-oracle-venv/bin/python --output target/distance-smoke-001
```

Success prints `Distance benchmark: passed (...)`. Both `--quick` and
`--exploratory` exclude the run from algorithm selection; `--exploratory` also
allows uncommitted sources. Numerical stress is a separate diagnostic:

```sh
python3 tools/compare_distances.py --suite stress --gudhi-python target/distance-oracle-venv/bin/python --output target/distance-stress-001
```

The pinned references have [known stress disagreements](#reference-backends-and-numerical-limits),
so a nonzero stress exit must be inspected separately from supported acceptance.

## Native build and routing controls

The builder compiles the ordinary public crate without instrumentation, then
compiles the standalone Rust worker with `--cfg cocycle_distance_bench`. This
private cfg enables the counter schemas and forced policies also used by kernel
tests. Conditional compilation removes counter statements, capacity traversals,
binary-search ablation and experimental arena storage from ordinary library
builds, including debug. The worker's `public` variant calls the separately
linked ordinary crate; native ablations call the instrumented private kernels.
Build metadata records this boundary. There is no public Cargo experiment feature.

The diagnostic `logical` Rust worker variant constructs diagrams and calls the
instrumented logical-view facade. It reports `algorithm_preparation_bytes` from
retained Vec capacities and `algorithm_preparation_buffers` from nonempty owned
preparation buffers. Bottleneck includes coordinates and its four index/cost
arrays per operand; Wasserstein includes the two prepared Point vectors. These
are not allocator-call counts. Grouping and solver scratch remain separate.
The facade has no owned point buffers, so `boundary_materialization_bytes` and
`boundary_materialization_buffers` are zero for this variant (null for others).
This boundary accounting excludes the worker's raw-to-diagram conversion.
The existing Bottleneck workspace counter overlaps preparation indexes; do not
sum it with preparation bytes as a process-memory estimate. `logical` timing
includes diagram construction and is diagnostic, not part of native rankings.
For public API before/after timings, construct both diagrams outside the clock,
use the ordinary crate, and record that distinct timing boundary and revisions.

Topp uses a portable scalar build with MSVC-only AVX2 dispatch disabled. Weighted
Topp matching uses compiler-dependent `long double`; Rust uses f64. Preserve
these representation differences when interpreting results.
The correctness `default` variant keeps Topp's unrestricted adaptive defaults.
Its `dense_parallel` candidate route can spawn threads at 262144 pairs. Measured
C++ variants instead use a serial configuration: the same pinned four-row density
test replaces only that route with `dense_blocked`. Prepared diagrams and the
extra routing check are counted inside the clock; inputs are prepared once.
Original Topp source remains unchanged. Records identify `cpp_threads: 1`, whether
the serial override applied, and a post-call observed thread count. The latter
alone is not a peak-thread measurement or the reason to assert serial behavior.

## Inputs, outputs and independent correctness

All workers read the same binary fixture: `COCDST1\0` (eight bytes), little-endian
u64 left/right counts, then left and right interleaved little-endian f64 endpoint
pairs. Exact length is checked; no f32 quantization occurs. Original bits, order
and multiplicity remain in the fixture.

Correctness compares one complete computed dimension with finite births and
ordered finite deaths or positive infinity. The Rust public adapter explicitly
removes finite diagonal points when making typed intervals to match Topp's raw
convention; this is not acceptance of diagonal points by `PersistenceInterval`.
The library tests separately exercise censored coverage, uncomputed dimensions,
invalid endpoints and provenance rejection. Performance fixtures are strictly
finite and off-diagonal.

The pinned GUDHI weighted wrapper fixes POT's network-simplex limit at 2,000,000
iterations. An iteration-limit warning is retained as an unavailable reference,
even if POT returned a scalar. Such a cell remains incomplete and is never
ranked; no automatic budget increase or tolerance relaxation is performed.

Each process receives `fixture metric variant`. Metrics are `bottleneck`
(L-infinity), `w1` (order 1, L-infinity), and `w2` (order 2, Euclidean, final square
root included). JSON identifies protocol/backend/metric/variant and retains the
scalar, time, process memory and counters. Infinity is `value_kind: "infinite"`
with `value: null`. NaN, negative results, crashes and malformed output are errors.

The tiny independent oracle enumerates every partial injection, including
unmatched points' diagonal costs, using exact rational arithmetic before the W2
square root. It shares no production search, graph or matcher code. Hand-derived,
empty, repeated, unequal-size, near-diagonal, negative-scale, tied and essential
examples supplement deterministic random inputs. Every supported backend is
checked against independent expectations, not merely Rust/Topp mutual agreement.

## Reference backends and numerical limits

The primary GUDHI bottleneck reference is the native `e=0` implementation at the
fixed head of [merged PR #1367](https://github.com/GUDHI/gudhi-devel/pull/1367),
which corrects the premature matching shortcut. It uses double-coordinate KD
trees and builds with `CGAL_DISABLE_GMP`; it is a correctness oracle only. Full
source identity and header hashes are retained. The correctness suite additionally
checks GUDHI 3.11.0 Hera bottleneck `delta=0`. Wasserstein uses explicit order/internal norm with
`keep_essential_parts=True`, no autodiff and POT's exact transport solver. The
isolated worker uses NumPy and loads only the requested metric's native backend;
Wasserstein disables optional POT GPU/autodiff imports. Hera requires removal of
diagonal points, which does not change the raw diagram distance. The pinned Linux
wheel's default `gudhi.bottleneck_distance(e=0)` returned 2 instead of the independent
oracle's 1.375 on an ordinary tiny fixture, even before POT was loaded. Its failed
attempt is retained separately; it is not counted as agreement. Hera's zero-delta
mode is independently checked on supported inputs and retains its own extreme-
value limitations in stress results. GUDHI/POT Python workers run in isolated
processes for correctness and never qualify as native timing references.
Backend identities are fixed in `sources.json`.
GUDHI 3.11.0 predates PR #1367 (merged 2026-08-27 and labeled 3.14.0); its default
backend is not the acceptance reference. Do not treat that known old-version
failure as a new Rust defect or silently call Hera the repaired implementation.
Versions/settings are retained; a missing reference makes validation fail.
Small dyadic bottleneck/W1 cases use zero tolerance. Other cases use the fixed
bound `64 * f64_epsilon * (n+m+1) * max(cost_scale, |expected|)`, with scale from
endpoint differences and diagonal costs, not the absolute coordinate origin.
Never increase tolerances after observing a mismatch.

`--suite supported` is the default. `--suite stress` checks very large/small and
adjacent-float cases; `--suite all` includes both. Pinned Topp's rounded-midpoint
diagonal projection disagrees with `(death-birth)/2` at adjacent floats. Stress
retains the real disagreement and returns nonzero, rather than labeling it as
agreement or broadening tolerance. Rust correctness and reference disagreement
remain separate evidence.

## Timing and sampling

One fresh native process computes one pair. Parsing and initial raw-pair layout
finish before timing. Validation, diagonal projection, preparation, solve and
temporary cleanup are counted; startup, fixture I/O and JSON formatting/transport
are excluded. All Rust native variants share the same finite validation/copying
boundary; the separate `public` adapter is correctness-only. Statistics counters
run inside the clock for both languages: these are instrumented time-to-result
measurements, not uninstrumented kernels. Raw inputs survive computation; solver
temporaries are destroyed before the clock stops. Copies/parser high-water marks
remain part of process memory.

Retain/discard one fresh warmup per cell; formal comparisons use at least twelve
measured processes. Sample counts round upward to a multiple of active variants.
Seeded shuffled cyclic blocks balance execution positions. Workers run serially,
single-threaded, without concurrent compilation or testing. Optional `--cpu` pins
to an allowed Linux CPU. Selected/inherited affinity and uncontrolled frequency
and host load are recorded.

Development seed is `20260922`, holdout seed is `20260923`. Families include
uniform, clustered, near-diagonal, duplicates, imbalance, separated, threshold
shell and sparse/dense adversarial cases. Default sizes are 8/32/128/512; request
2048/4096 explicitly where limits permit. `--families`, `--sizes`, `--metrics` and
`--groups` preregister a smaller study; its conclusions remain scoped to that set.

Each family varies its actual geometry across the two seeds, including repeated
templates and regular sparse/threshold cases. The controller rejects identical
tuning and holdout fixture hashes for the same family and size.

## Controlled optimization contrasts

| Group | Controlled contrast |
| --- | --- |
| `baseline` | Topp serial adaptive configuration and Rust adaptation, always retained; unrestricted default remains a correctness reference |
| `search` | Forced quickselect versus binary with other settings held fixed within each language; adaptive baseline remains a separate row |
| `scratch` | Rust forced refinement with/without scratch reuse |
| `matching` | Rust forced refinement with/without matching reuse |
| `clipping` | Forced quickselect with/without candidate clipping in both languages |
| `arena` | Forced sparse vectors versus the arena plus scratch/heap-reuse candidate in both languages, plus Rust `adaptive_arena` under its ordinary default routing |

Topp has no corresponding exposed scratch/matching-reuse controls or adaptive
arena layout switch; these missing counterparts are limitations, never invented
equivalent experiments. Local sparse variants disable duplicate/component/greedy
shortcuts and keep those settings fixed within their pair. They cannot select an
adaptive default alone: `adaptive_arena` tests the actual candidate integration.
Route/counter records show whether the intended path ran. Rust
`direct_cost_fallbacks > 0`, `sparse_solves == 0`, no positive edges, or zero
residual-storage capacity does not test arena layout. Empty sparse calls can
return before constructing a network. The weighted arena contrast includes
scratch/heap reuse across augmentations; it does not isolate that reuse from
contiguous storage. Bottleneck scratch reuse has its own separate control.
Per-language retained/disabled ratios answer whether an optimization survives
migration; cross-language absolute times answer a different question.

## Optional kernel diagnostics

The benchmark-only `--profile-rust` flag copies the private Rust distance sources
and worker into the fresh artifact build directory, then inserts safe RAII
timing scopes at uniquely checked function markers. Production files and the
ordinary worker build stay unchanged. Generated sources retain `forbid(unsafe_code)`,
and build metadata records every original/generated source hash and hook.

Each Rust sample retains calls and integer `elapsed_ns` by phase alongside raw
stderr. Scopes include their nested calls and hook overhead; they must not be
summed as disjoint phases. Early returns close their scopes through `Drop`.
Preparation, candidate generation, KD construction/decisions, weighted graph
generation/components and dense/sparse/direct-cost solves can be inspected.
Public-library correctness references remain uninstrumented. A profiling run
always has evidence class `kernel_diagnostics` and cannot select an algorithm,
even with 12 samples. C++ is a numerical control here, not a timing comparator.

## Memory metrics

Workers read pre-call `VmRSS`/`VmHWM` and final `VmHWM` after cleanup, before JSON.
Report absolute peak and high-water growth separately; zero growth does not mean
no allocations. RSS includes runtime, input, allocator retention and outputs.
Explicit container capacities are not allocator peaks.
Rust Wasserstein records separate maximum capacities for retained graph storage,
residual storage (including vector headers or arena offsets), and search scratch
(vectors and heap). These maxima exclude construction temporaries and must not
be summed into a simultaneous memory peak. Compare them only within matching
categories and representation sizes; C++ long-double edges can be larger.
`--address-space-mib`
(default 2048 MiB) caps virtual address space, and `--timeout` (default 60 seconds)
limits whole-process wall time. Neither is a library budget; timeout is censored
evidence, not a measured runtime equal to the limit.

## Selection gates and repeat procedure

`results.json` retains every sample, warmup, order, exit and mismatch. After a
worker fails, its remaining attempts are `not_run` while others continue. No
mismatched/incomplete case is ranked. Valid summaries include median/min/max and
RSS. Changed sources invalidate a run. `--quick` and `--exploratory` are resource
snapshots and never select algorithms; formal measurements reject dirty sources.

Selection compares candidate/R0 time and peak-RSS ratios, with equal size weight
within each family and equal family weight. R0 is the Rust adaptive baseline
(`cocycle:baseline`). Ratios use per-case medians, then geometric means across
sizes and families. The held-out composite
`sqrt(time_ratio * peak_rss_ratio)` must be at most 0.95; neither metric may exceed
1.10 in any default-applicable family. A first qualifying result is only
`eligible_pending_independent_repeat`: repeat the same gates in a second fresh
run on identical sources/fixtures before deciding. Noise and fixed process
overhead are inconclusive; retain the simpler safe baseline on ties. Forced
local ablations are not automatically global default candidates.

Commit implementation and harness before formal measurement, then record the
full commit SHA. On Linux, the following runs all default families, metrics and
ablation groups at sizes 8/32/128/512 with development and held-out inputs:

```sh
distance_revision=$(git rev-parse --short=12 HEAD)
python3 tools/benchmark_distances.py --samples 12 --order-seed 2401 --gudhi-python target/distance-oracle-venv/bin/python --output "target/benchmarks/commit-${distance_revision}/distance/run-001"
```

For selection, repeat on the same committed sources and fixture hashes in a new
output directory with a different order seed. Apply the gates to both rounds;
the controller reports only eligibility and never changes the library default.
Specify `--families`, `--sizes`, `--metrics` and `--groups` before execution if
limiting the study. Formal commands are reproduction instructions, not evidence
that a run has completed. Shared-host load and fixed RSS overhead can leave
differences inconclusive.

## Inspect and retain artifacts

| Artifact | Interpretation |
| --- | --- |
| `summary.json` | Overall validation, source consistency, and benchmark selection eligibility |
| `results.json` | Raw reference results, samples, warmups, execution order, failures and route counters |
| `measurements.json` | Benchmark summaries in milliseconds and KiB; empty if overall validation fails |
| `environment.json` | Commit, source fingerprint, toolchains, protocol and evidence class |
| `fixtures/` | Shared binary inputs; SHA-256 hashes are recorded in results |
| `build/build.json`, `build/build.log` | Build commands, source/binary hashes and compiler output |

Keep formal artifacts under fresh
`target/benchmarks/commit-<sha12>/distance/run-<NNN>/` directories. Preserve raw
fixtures, JSON, logs and metadata outside Git. A report identifies its measured
commit separately from its documentation commit and links any published evidence.
A local artifact path does not provide public access.

For changes to this suite, run the applicable [contribution checks](../../CONTRIBUTING.md#verification),
including the focused controller tests and worker formatting:

```sh
python3 -m unittest discover -s tools -p 'test_*distances.py'
rustfmt --edition 2024 --check benches/distances/cocycle.rs
```
