# Diagram-distance concurrent resources

[Reports](README.md) / [Protocol](../distances/README.md#concurrent-resource-study)

## Question and conclusion

This first tranche of [Issue #19](https://github.com/Aequiludium/cocycle-rs/issues/19)
asks how the current adaptive route and existing forced sparse layouts behave
under homogeneous concurrency. Separate comparative performance groups and
resource snapshots cover finite W1/W2 pairs of size 128 on one Linux/WSL host.
All planned groups passed; production Rust and Cargo files remain identical to
the [integrated R0](phase2-baseline.md).

On the held-out sparse input, adaptive component decomposition avoids a full
sparse solve and is more than thirty times faster than either forced route.
The arena contrast has a smaller observed benefit over nested vectors, with
overlapping sample ranges and no independent formal repeat. It cannot select a
production default. Four workers still increase throughput; these three points
do not locate a saturation knee or identify the limiting physical resource.

This is an initial concurrency study, not completion of Issue #19. No model,
router or production optimization is delivered by this tranche. The
[lifecycle conclusion](distance-lifecycle.md) remains a separate lifetime axis:
prepared-state opportunities are scoped, and the tested scratch/output retention
does not pass its joint incremental gate.

## Measured identities and environment

| Identity | Recorded value |
| --- | --- |
| PR | No associated experiment PR as last inspected; local branch `codex/phase2-resource-study` |
| Production baseline and private candidate source | `b2c3bd5eebf0c1193f30ea651012b8ae61b274c1` |
| Frozen worker/controller commit | `0156860b614fed9f13d60e254e53435858e97d3f` |
| Protocol | `cocycle-distance-resources-v1` |
| Source fingerprint | `5127c769399a67e3758852e197faceb3617572302b5e4d2d0679f67fa4def679` |
| Generated worker SHA-256 | `1ff3ec10002c8dde13fce6d118dc503db8725fa2380ac2a4f59c3d173b3adc92` |
| Worker binary SHA-256 | `f004033d51442c23d5f1120924722eabb7552237e8b9cf0333cf88f17d0c5354` |
| Compiler/build | Rust 1.91.0, release library, `rustc -O --edition=2024 --cfg cocycle_distance_bench` |
| Host | Intel Core Ultra 7 155H; WSL2 Linux 6.6.87.2; 22 allowed logical CPUs |
| Placement | Processes pinned to CPUs 0, 1, 2, 3 as needed; thread process affinity spans CPUs 0 and 1 |
| Limits and controls | 120-second result wait; no explicit resource-study address-space cap; frequency and external host load uncontrolled |

The source fingerprint hashes canonical sorted JSON containing build metadata
`originals` and `owners` maps. The audit checks every original against R0, every
owner against the frozen harness commit, and generated source/binary hashes.
The report revision is distinct from both measured commits. WSL logical CPU
numbers do not establish physical-core isolation, NUMA placement or cache topology.
No compilation, tests, compression or other formal run overlapped these groups.

Run IDs are `linux-trace-uniform`, `linux-process-uniform`, `linux-trace-sparse`,
`linux-process-sparse`, and `linux-thread-uniform`, executed in that order under
one orchestration script. Manifests retain UTC creation metadata; completion UTC
is not recorded. Elapsed clocks are monotonic and do not depend on realtime UTC.

## Workload and coverage

Each family has one 128-by-128 finite f64 diagram pair per split, four distinct
pairs in total, using the maintained generator and seed 20260922 for tuning and
20260923 for holdout. W1 and W2 both use the internal
L-infinity cost convention. Raw parsing and an in-process warmup precede readiness;
every measured call includes validation/copying, preparation, solve and cleanup.
Calls retain no prepared representation or explicit solver workspace across jobs.
The allocator, input and code caches can remain warm within a group.

| Suite | Cells | Measured groups | Discarded groups | Jobs per worker/group |
| --- | ---: | ---: | ---: | ---: |
| Single-process traces, both families | 8 | 72 | 24 | uniform 64; sparse 2048 |
| Process performance, N = 1/2/4 | 24 | 864 | 72 | uniform 64; sparse 2048 |
| Uniform native thread performance, N = 2 | 4 | 144 | 12 | 64 |

The independent audit finds 36 complete cells, 1080 measured groups, 108 discarded
groups, 2232 measured child processes, 204 discarded-group child processes, and
2223360 measured calls. No formal group failed or was excluded. The 1080 groups
include the 72 trace snapshots, which never enter performance rankings.
Source identity, execution schedule, input hashes, job counts, throughput,
native quantiles, diagnostic consistency and trace integrals passed rechecking.

Every call is asserted against separately linked public R0. Two earlier 16-cell
tiny-input smokes exercise the independent exact oracle, handshake, processes,
threads and traces. The four-cell exploratory pilot sets observation duration
only. Neither smokes nor pilot enter these tables. No new external Topp/GUDHI
comparison is claimed; this tranche changes instrumentation and private options,
with the unchanged kernel's external comparisons owned by the R0 report.

| Fixture | SHA-256 |
| --- | --- |
| sparse-holdout | `f3cb4eea268e2675431860e024fc4cd74c2f34e0a618f532d8e89c128a89c155` |
| sparse-tuning | `3ecb8b53d33c935956f229fd0190b6e9dfa1a50b01bb35d6072f6ec3c085c9ba` |
| uniform-holdout | `f722e6b7d741e634b61840364f594fa1cd12db82ec250716b789eeb1260a8a1c` |
| uniform-tuning | `9a6427be2e36f858c481b00edc7f08ca8a120e6d63fffae5c26f3435a1e20439` |

## Timing and resource boundaries

Each cell/route has one fresh discarded group and twelve measured groups, except
the three-group trace snapshots. A seeded shuffled cyclic schedule balances
route positions. Tables use medians across groups and retain every sample.
Native p95 is nearest-rank within each group, followed by its median across
groups; it is not a pooled-call p95 or confidence bound. The selected tables below
show all held-out family/metric/route combinations. Tuning and complete resource
columns remain in `analysis/matrix.json` with the raw records.

Parent monotonic dispatch-to-result time includes go/result transport and group
scheduling; it excludes process creation, parsing and readiness warmups. Native
per-job clocks exclude result transport. The native group clock includes its
barrier release and thread joins. Parent time divided by the slowest child's
native group time has cell medians 1.0002-1.0142 in the process performance runs.
These boundaries are close for these job counts, not interchangeable.

CPU demand is scheduled core-seconds from procfs ticks, including stalls while
scheduled. The observed endpoint includes result production/transport before ack;
CPU time is neither retired work nor a stall counter. Tick quantization limits
small-demand precision. PSS apportions shared pages across processes and is a
proxy for aggregate memory, not cgroup accounting. RSS is separately retained.
Endpoint memory is after result transport and may include serialization buffers,
parser/input buffers, runtime and allocator retention. It is not live solver memory.

Trace sampling requests 2 ms. Actual intervals across measured samples have
minimum 0.174 ms, median 2.224 ms and maximum 11.652 ms, including final boundary
intervals. Sequential procfs reads and observer scheduling limit resolution.
The final observation follows sampler shutdown, so its window can exceed the
dispatch-to-result window. Integrals and peak rectangles use that same observed
window. They must not substitute the throughput clock.

DRAM bytes and LLC misses are unavailable (`null`) on this host; cgroup root is
not writable. No DRAM capacity, bandwidth or three-resource fit is inferred.

## Initial process curves

Throughput is jobs/s, displayed as median [minimum, maximum] over twelve groups.
Each worker runs the same pair serially; concurrency is between callers, not
parallel execution within one numerical solve.

| Held-out family / metric | Route | N1 | N2 | N4 | N4/N1 |
| --- | --- | ---: | ---: | ---: | ---: |
| uniform / W1 | baseline | 273.8 [257.3, 298.1] | 493.1 [460.1, 522.1] | 795.2 [612.9, 890.2] | 2.904 |
| uniform / W1 | vectors | 273.1 [251.1, 286.5] | 490.5 [458.7, 517.4] | 781.8 [661.0, 828.2] | 2.863 |
| uniform / W1 | arena | 284.5 [259.6, 292.6] | 513.4 [440.5, 529.2] | 813.6 [747.1, 902.4] | 2.860 |
| uniform / W2 | baseline | 280.0 [273.9, 291.8] | 508.9 [477.9, 525.1] | 801.0 [746.3, 839.6] | 2.860 |
| uniform / W2 | vectors | 266.1 [256.0, 277.1] | 463.6 [438.9, 493.9] | 744.6 [706.9, 791.9] | 2.799 |
| uniform / W2 | arena | 275.6 [261.5, 293.1] | 488.2 [462.4, 509.7] | 756.7 [660.3, 812.8] | 2.745 |
| sparse / W1 | baseline | 26656.0 [23452.6, 29282.3] | 49137.2 [33773.9, 52669.1] | 80426.8 [73326.4, 84212.6] | 3.017 |
| sparse / W1 | vectors | 763.2 [732.3, 802.9] | 1413.8 [1316.0, 1483.1] | 2341.9 [2228.2, 2466.6] | 3.069 |
| sparse / W1 | arena | 818.3 [787.5, 845.4] | 1535.3 [1432.6, 1570.2] | 2522.4 [2169.9, 2584.0] | 3.082 |
| sparse / W2 | baseline | 25244.8 [22702.5, 27902.5] | 47439.3 [45225.0, 51790.7] | 78218.7 [75545.1, 80552.5] | 3.098 |
| sparse / W2 | vectors | 728.6 [697.5, 743.0] | 1404.1 [1297.2, 1456.7] | 2255.6 [2171.3, 2384.1] | 3.096 |
| sparse / W2 | arena | 768.2 [749.3, 798.5] | 1496.0 [1418.5, 1534.2] | 2430.0 [2263.7, 2511.4] | 3.163 |

Held-out N4/N1 gains range from 2.745 to 3.163, below ideal fourfold scaling.
Tuning gains range from 2.864 to 3.326. These observations motivate a wider,
densified sweep; they do not prove a bandwidth, cache or CPU saturation mechanism.

The next table shows native median/p95 at N4 in microseconds, scheduled CPU
microseconds per completed call at N1/N4, and N4 endpoint PSS. CPU demand divides
each group's core-seconds by its fixed completed-call count before taking medians.

| Held-out family / metric | Route | N4 job median / p95 (us) | CPU us/call N1 / N4 | N4 PSS (MiB) |
| --- | --- | ---: | ---: | ---: |
| uniform / W1 | baseline | 4656.2 / 6243.5 | 3593.8 / 4746.1 | 2.403 |
| uniform / W1 | vectors | 4766.6 / 6383.5 | 3671.9 / 4863.3 | 2.414 |
| uniform / W1 | arena | 4611.2 / 5913.9 | 3437.5 / 4648.4 | 2.400 |
| uniform / W2 | baseline | 4655.1 / 6184.6 | 3593.8 / 4785.2 | 2.407 |
| uniform / W2 | vectors | 5057.0 / 6679.4 | 3750.0 / 5195.3 | 2.406 |
| uniform / W2 | arena | 4961.6 / 6425.1 | 3671.9 / 5058.6 | 2.405 |
| sparse / W1 | baseline | 46.0 / 58.3 | 34.2 / 44.6 | 2.329 |
| sparse / W1 | vectors | 1662.6 / 2259.5 | 1296.4 / 1689.5 | 2.720 |
| sparse / W1 | arena | 1537.8 / 1990.0 | 1210.9 / 1560.7 | 2.496 |
| sparse / W2 | baseline | 47.4 / 59.1 | 36.6 / 47.0 | 2.321 |
| sparse / W2 | vectors | 1716.6 / 2286.6 | 1374.5 / 1741.3 | 2.728 |
| sparse / W2 | arena | 1599.6 / 2116.4 | 1306.2 / 1622.9 | 2.490 |

Scheduled demand rises with concurrency on these cells. CPU placement,
frequency, scheduling and contention can all contribute; this is evidence
against silently assuming constant observed demand, not an attribution of cause.

## Same-process thread contrast

The thread process shares raw input and allocator state, with native worker
threads under a two-CPU process affinity. Processes have separate inputs and
allocators with individual pinning. Separate calls still use fresh pair-local
solver state. These contracts must not be pooled.

| Held-out uniform metric | Route | Two processes jobs/s | Two threads jobs/s [min, max] | Thread median / p95 (us) | Process / thread endpoint PSS (MiB) |
| --- | --- | ---: | ---: | ---: | ---: |
| W1 | baseline | 493.1 | 517.8 [496.7, 542.2] | 3780.3 / 4448.5 | 1.585 / 1.274 |
| W1 | vectors | 490.5 | 507.7 [485.3, 529.8] | 3809.8 / 4749.0 | 1.592 / 1.272 |
| W1 | arena | 513.4 | 522.7 [499.7, 539.6] | 3710.8 / 4461.1 | 1.583 / 1.282 |
| W2 | baseline | 508.9 | 517.1 [474.0, 537.3] | 3794.5 / 4586.5 | 1.589 / 1.261 |
| W2 | vectors | 463.6 | 486.0 [406.0, 493.8] | 3980.7 / 4857.4 | 1.589 / 1.297 |
| W2 | arena | 488.2 | 495.3 [473.1, 514.5] | 3899.6 / 4735.1 | 1.586 / 1.265 |

Threads have lower observed endpoint PSS here, while performance ranges overlap
several process counterparts. This single N2 family contrast cannot select a
general execution model or quantify high-concurrency allocator interference.

## Separate trace snapshots

The three measured groups per cell estimate repeated-call group trajectories,
not resolved single-operation phases. Warm input/code/allocator residence is
included. Peak is sampled simultaneous PSS at N1, not individual-process VmHWM.
For each group, I/R = trapezoid PSS integral / (sampled peak PSS times observed
window). Columns below are medians across three groups.

| Held-out family / metric | Route | CPU core-seconds/group | Sampled peak PSS (MiB) | Integral (MiB s) | I/R |
| --- | --- | ---: | ---: | ---: | ---: |
| uniform / W1 | baseline | 0.250 | 1.408 | 0.3523 | 0.9784 |
| uniform / W1 | vectors | 0.280 | 1.376 | 0.3768 | 0.9870 |
| uniform / W1 | arena | 0.250 | 1.240 | 0.3234 | 0.9911 |
| uniform / W2 | baseline | 0.270 | 1.338 | 0.3566 | 0.9797 |
| uniform / W2 | vectors | 0.290 | 1.333 | 0.3837 | 0.9734 |
| uniform / W2 | arena | 0.280 | 1.213 | 0.3380 | 0.9932 |
| sparse / W1 | baseline | 0.080 | 1.114 | 0.0925 | 0.9774 |
| sparse / W1 | vectors | 2.860 | 1.196 | 3.1719 | 0.9775 |
| sparse / W1 | arena | 2.500 | 1.165 | 2.9414 | 0.9785 |
| sparse / W2 | baseline | 0.080 | 1.100 | 0.0889 | 0.9534 |
| sparse / W2 | vectors | 2.730 | 1.172 | 3.0766 | 0.9549 |
| sparse / W2 | arena | 2.600 | 1.143 | 2.7971 | 0.9550 |

Across tuning and holdout cell medians, I/R is 0.9534-0.9932. The peak rectangle
therefore exceeds the observed integral modestly in these repeated-call groups.
This does not validate the static peak model on single-operation phases, large
inputs or staggered concurrent peaks. Sub-cadence allocation peaks remain unseen.

## Actual routes and capacity diagnostics

Representative last-job counters on held-out N1 process groups are consistent
across checked groups/children. Each worker reports its first thread's final job;
this is not a separately sampled counter trace for every job. All reported dense
solver counts are zero. Neither family exercises a forced dense comparison.

| Held-out family / metric | Route | Candidate pairs / positive edges | Components | Sparse solves | Residual capacity (bytes) | Search scratch capacity (bytes) |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| uniform / W1 | baseline | 1933 / 1344 | 2 | 1 | 143304 | 14024 |
| uniform / W1 | vectors | 1933 / 1344 | 0 | 1 | 143608 | 14072 |
| uniform / W1 | arena | 1933 / 1344 | 0 | 1 | 103536 | 37240 |
| uniform / W2 | baseline | 1860 / 1268 | 2 | 1 | 140088 | 14264 |
| uniform / W2 | vectors | 1860 / 1268 | 0 | 1 | 140392 | 14312 |
| uniform / W2 | arena | 1860 / 1268 | 0 | 1 | 99392 | 38760 |
| sparse / W1 | baseline | 128 / 128 | 128 | 0 | 0 | 0 |
| sparse / W1 | vectors | 128 / 128 | 0 | 1 | 47152 | 8240 |
| sparse / W1 | arena | 128 / 128 | 0 | 1 | 26648 | 18480 |
| sparse / W2 | baseline | 128 / 128 | 128 | 0 | 0 | 0 |
| sparse / W2 | vectors | 128 / 128 | 0 | 1 | 47152 | 8240 |
| sparse / W2 | arena | 128 / 128 | 0 | 1 | 26648 | 18480 |

Adaptive sparse input has 128 components, each resolved without a full sparse
network. Forced routes disable reductions and construct one full sparse solve.
Their comparison against baseline changes policy and reductions as well as layout.
The controlled storage comparison is vectors versus arena, including arena's
existing scratch/heap reuse. Arena reduces residual capacity while increasing
reported scratch capacity; category maxima are not summed into a simultaneous
memory peak. Endpoint PSS need not track these capacities proportionally.

## Evidence, reproduction and remaining work

Local-only raw evidence is
`target/benchmarks/commit-b2c3bd5eebf0/phase2-resources/`. It contains manifests,
fixtures, all warmups/samples, child stderr, run logs, host metadata, preregistration
and audited `analysis/matrix.json` / `analysis/audit.json`. No durable external
artifact URL is configured. This report does not provide public raw-data access.

The initial-tranche archive is
`target/deliverables/phase2-resources-initial-b2c3bd5eebf0/phase2-resources-initial-evidence.zip`:
16746202 bytes, SHA-256
`c41354e30b78187029ddd1c7f1a89a3d34e0c2d2944017e2712c8269e492a27a`.
All 4603 inventoried members passed length/SHA-256 verification. The package
contains formal raw records, preregistration, orchestration, original smokes and
pilot, generated source/binary/build metadata, plots, analyses and Git snapshots
of R0 and the frozen harness. Cargo caches are excluded.

The frozen protocol can rebuild its worker in a fresh trace output directory.
The following shell commands reproduce the planned matrix with fresh directories;
run them serially on the unchanged kernel and frozen worker/controller sources:

```sh
for family in uniform sparse; do
    jobs=64
    if [ "$family" = sparse ]; then jobs=2048; fi
    python3 tools/benchmark_distance_resources.py --families "$family" --sizes 128 --workers 1 --jobs "$jobs" --samples 3 --trace --output "target/resources-repeat-trace-$family"
    python3 tools/benchmark_distance_resources.py --worker-dir "target/resources-repeat-trace-$family/build" --families "$family" --sizes 128 --workers 1 2 4 --jobs "$jobs" --samples 12 --output "target/resources-repeat-process-$family"
done
python3 tools/benchmark_distance_resources.py --worker-dir target/resources-repeat-trace-uniform/build --families uniform --sizes 128 --workers 2 --modes threads --jobs 64 --samples 12 --output target/resources-repeat-thread-uniform
```

These commands are a repeat recipe, not a claim of an independent formal repeat.
The original worker was reused from the second smoke build after verifying all
source and binary fingerprints. Earlier validation passed 96 Linux tool tests,
two Windows focused tests, both tiny smokes, standalone formatting and project
source/documentation checks. Hosted CI for this harness is not claimed.

Issue #19 still requires resolved single-operation/phase traces, broader inputs
and thread curves, a wider sweep densified around an observed knee, mixed-route
interference, phase replay and held-out M0 peak / M1 static integral / M2
load-dependent / M3 interference comparisons. DRAM evidence requires a suitable
measurement host or an explicitly limited model. Final routing must incorporate
Issue #18 lifetime costs and pass an independent repeat and regret comparison.
Only stable, predictive gains can justify a cheap offline production change;
negative or inconclusive results retain the current simpler routing.
