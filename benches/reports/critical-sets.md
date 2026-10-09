# Cocycle, dense reproduction and Oineus: critical-set comparison

[Benchmarks](../README.md) / [Protocol](../optimization/README.md) /
[Reporting rules](../reporting.md)

## Question and conclusion

This local comparative study asks what the public sparse critical-set workspace
costs relative to the original dense reproduction, serial Oineus Partial/full U,
and the repository's existing persistence paths. The sparse workspace improves
cold cost on the larger grid, delayed fan and complete 2-skeleton relative to
dense reproduction, but regresses on mixed queries of the complete 3-skeleton.
Its wide repeated queries are slower than Oineus on all eight median rows.
Partial U reduces Oineus cold cost; full U wins repeated-query cost. There is
no universal sparse or Partial U speedup on these workloads.

Existing persistence methods are measured in a separate table because they
produce different outputs from critical-set queries. These are synthetic
endpoint batches, not an optimization-convergence comparison. All evidence is
local-only; no associated PR, remote CI or public artifact archive exists.

## Measured revision and environment

| Item | Recorded value |
| --- | --- |
| Measured candidate and harness | `ceea0f4e47e4c296902549e459432eaa3bbb2684`, clean local feature commit |
| Production implementation | `a380cba76d417c62be46d50495e4d7cce5fef2b3`; all production Rust bytes unchanged in measured commit |
| Production source fingerprint | SHA-256 `6a88f068d0a6eff552717ea414e65d17a00a6c9a3274f89aeb775cf57c2c64f3`, over the sorted compact JSON map of relative `src/**/*.rs` paths to SHA-256 |
| Dense baseline | `167026e7ee8dbaafe8abaf5e325149af6e93f660`, blob `53bd137d16e122c2ad00a9859188ff317a362472`; generated tie-order and birth-mask adaptations |
| Oineus | 0.9.39, `e52814a1ffb5b8a81e71ff1e93b4c14194673f0f`; all 47 consumed upstream files verified against that commit's archive |
| Suite / attempt | `critical-sets-native-v1` / `commit-ceea0f4e47e4/critical-sets/run-001` |
| UTC timing interval | 2026-10-09 11:01:52.809649 to 11:03:09.486949; build precedes this interval |
| Platform | Intel Core Ultra 7 155H; Ubuntu 24.04 under WSL2, Linux 6.6.87.2, x86_64, 22 exposed vCPUs |
| Rust | 1.91.1 (`ed61e7d7e242494fb7057f2657300d9e77bb4fcb`), LLVM 21.1.2; edition 2024, `-C opt-level=3` |
| C++ | g++ 13.3.0; C++20, `-O3 -pthread -DOINEUS_DISABLE_ICECREAM`; consumed dependencies recorded with `-MMD` |
| Controls | Every worker pinned to CPU 0; serial execution, 60-second timeout, 2-GiB address-space cap; frequency and host load uncontrolled |

Rust uses 64-bit indices, explicit field coefficients and its production sparse
column representation. Oineus uses int32 indices and implicit F2 sparse columns.
Oineus reduction is serial, retains V and disables clearing. The Partial/full U
paths share one binary. The C++ adapter applies the same absolute-displacement
merge as Rust; it does not time Oineus `TopologyOptimizer::combine_loss`,
autograd or parallel ELZ restoration. No Oineus implementation is copied into
production Rust. The dense implementation eagerly computes both primal and
dual R/V/U by dimension, even for a single primal query.

## Workloads and validation

All workers receive identical f64 supplied filtrations over F2, with full
coverage and preserved multiplicity. Critical queries include zero-length
pairs; exported diagrams exclude zero bars and retain essential intervals.
There is no approximation, censoring, distance evaluation or metric assumption.
Fixture seed is 220316748. Native preparation differs by backend and precedes
timing; this comparison starts from prepared source representations.

| Fixture | Vertices | Edges | Triangles | Tetrahedra | Requested homology |
| --- | ---: | ---: | ---: | ---: | --- |
| grid-8 | 64 | 161 | 98 | 0 | H0/H1 |
| grid-16 | 256 | 705 | 450 | 0 | H0/H1 |
| delayed-fan-64 | 64 | 126 | 63 | 0 | H0/H1 |
| delayed-fan-256 | 256 | 510 | 255 | 0 | H0/H1 |
| clique-12 | 12 | 66 | 220 | 0 | H0/H1 |
| clique-20 | 20 | 190 | 1140 | 0 | H0/H1 |
| tetra-skeleton-8 | 8 | 28 | 56 | 70 | H0/H1/H2 |
| tetra-skeleton-14 | 14 | 91 | 364 | 1001 | H0/H1/H2 |

`primal-U1` moves one finite death upward by 0.037 of the filtration span.
`mixed-Q8` and `mixed-Q64` rotate death up/down and birth up/down, with 8/64
requests and displacements of 0.037/0.371 of the span. Pair selection is
`(i * 17) % pair_count` and can repeat. Maximum absolute displacement wins
conflicts; equal displacement keeps the first proposal. These are fixed
synthetic requests, not distributions sampled from an application optimizer.

GUDHI 3.12.0 independently matched the interval multisets on all eight full
fixtures before timing. The saved preflight fixtures were checked against the
actual formal-run fixtures after timing. Every measured and warmup output
matched an independent Python vertex-face/set-XOR pairing oracle. Critical
methods also matched exact pair identities, essential births and dense merged
targets. The native smoke covered 132/132 quick cells, and all 86 Linux Python
tool tests passed. Production tests belong to the earlier implementation commit;
this run changes benchmark code only.

Representatives request cycles in every homology dimension at the middle
simplex's value. Independent face-XOR checks validate closure, term dimension,
term filtration and the active-interval count. The resulting cycle counts are
20, 62, 21, 81, 1, 1, 26 and 30 in the table's fixture order. These certificates
do not independently prove homology-basis independence.

## Measurement and sampling

The [protocol](../optimization/README.md#timing-and-validation) specifies all
boundaries. Cold time is same-process build plus first targets plus workspace
destruction. Pairing transport is excluded; owned targets remain alive. Warm
time rebuilds a separate workspace, prewarms the same batch, then averages
1024/128/16 repeated batches for U1/Q8/Q64, including each result's destruction.
Warm time excludes workspace build and destruction. Context time includes
public computation, normalized intervals and result destruction; representative
certificate validation is excluded. Input parsing, prepared-source construction,
startup, JSON formatting and controller checks are outside all timers.

The saved plan precedes worker execution. One discarded warmup and 12 fresh
measured processes per cell ran serially, with shuffled/rotated backend order
(seed 220316749). Positions differ by at most one for the five-backend context
groups. All 132 cells completed: 1584/1584 measured workers and 132/132 warmups,
zero mismatches, crashes or timeouts. An additional 24 untallied dense calls
calibrated expected targets. Four incompatible `flag` context cells were
excluded by design, representing 48 unrequested measured calls.

Tables show median [minimum, maximum] with all 12 observations retained. The
selected critical tables cover the wide mixed workload; the full matrix also
contains U1, Q8, separate phases, memory and per-round ratios in local artifacts.
Paired ratios match process rounds; a ratio `critical_ms / reference_ms` above
one means the critical workspace is slower. There is no pooled speed score.

## Wide mixed queries

Cold workflow, milliseconds:

| Fixture | Sparse critical | Dense | Oineus Partial U | Oineus full U |
| --- | ---: | ---: | ---: | ---: |
| grid-8 | 0.941 [0.861, 1.461] | 0.792 [0.733, 1.052] | 1.339 [1.081, 2.343] | 1.603 [1.511, 1.938] |
| grid-16 | 3.824 [3.574, 4.352] | 14.174 [13.328, 17.911] | 2.098 [1.909, 2.450] | 2.931 [2.695, 3.176] |
| delayed-fan-64 | 0.732 [0.680, 0.863] | 0.517 [0.476, 1.097] | 1.160 [1.055, 1.641] | 1.510 [1.371, 2.585] |
| delayed-fan-256 | 1.629 [1.509, 2.008] | 5.338 [4.940, 6.047] | 1.484 [1.314, 2.207] | 2.043 [1.724, 2.451] |
| clique-12 | 1.481 [1.182, 2.164] | 1.121 [0.909, 1.234] | 1.485 [1.246, 2.560] | 1.954 [1.548, 2.554] |
| clique-20 | 5.246 [4.992, 5.516] | 27.904 [26.459, 29.881] | 2.329 [2.102, 3.234] | 2.879 [2.758, 3.637] |
| tetra-skeleton-8 | 1.170 [1.016, 1.728] | 0.532 [0.464, 0.709] | 1.577 [1.121, 2.564] | 2.031 [1.510, 3.934] |
| tetra-skeleton-14 | 120.760 [117.379, 127.485] | 94.440 [89.458, 102.085] | 6.806 [6.339, 9.049] | 9.573 [8.961, 10.617] |

Warm repeated batch, milliseconds:

| Fixture | Sparse critical | Dense | Oineus Partial U | Oineus full U |
| --- | ---: | ---: | ---: | ---: |
| grid-8 | 0.020 [0.019, 0.024] | 0.033 [0.031, 0.040] | 0.012 [0.011, 0.019] | 0.007 [0.007, 0.008] |
| grid-16 | 0.053 [0.050, 0.120] | 0.106 [0.099, 0.149] | 0.024 [0.022, 0.029] | 0.015 [0.014, 0.017] |
| delayed-fan-64 | 0.013 [0.013, 0.015] | 0.025 [0.024, 0.036] | 0.011 [0.010, 0.022] | 0.007 [0.006, 0.007] |
| delayed-fan-256 | 0.029 [0.027, 0.037] | 0.081 [0.077, 0.169] | 0.019 [0.017, 0.039] | 0.015 [0.014, 0.019] |
| clique-12 | 0.111 [0.090, 0.123] | 0.053 [0.047, 0.069] | 0.048 [0.040, 0.062] | 0.024 [0.020, 0.029] |
| clique-20 | 0.161 [0.155, 0.260] | 0.147 [0.137, 0.234] | 0.076 [0.073, 0.125] | 0.047 [0.043, 0.063] |
| tetra-skeleton-8 | 0.065 [0.061, 0.123] | 0.029 [0.026, 0.040] | 0.023 [0.020, 0.041] | 0.013 [0.012, 0.021] |
| tetra-skeleton-14 | 1.594 [1.514, 1.998] | 0.289 [0.262, 0.419] | 0.325 [0.287, 0.379] | 0.191 [0.180, 0.297] |

The larger grid/fan/clique all beat dense cold time in every paired round.
On the 3-skeleton-14, sparse/dense median paired ratios are 1.297 cold and
5.401 warm, with zero sparse wins in either group. Sparse/Oineus Partial U is
17.788 cold and 4.846 warm, also with zero sparse wins. Small-input and
clique-20 warm results retain overlap/variation rather than imply a stable
blanket ranking.

The 3-skeleton-14 sparse phase medians are 67.413 ms build, 51.797 ms cold
targets and 1.182 ms destruction. Oineus Partial U is 4.190, 2.461 and 0.096 ms.
Phase medians do not necessarily sum to workflow medians. This directs further
investigation toward primal/dual construction and reduction as well as warm
support solving; no profiler attribution to a specific data structure is made.

For U1, Oineus full/Partial U ratios of cold medians are 1.15-1.54 across the
eight fixtures. For Q64 they are 1.20-1.41. Full U has lower warm medians in
every corresponding row: Q64 Partial/full ratios are 1.27-1.97. A cached full
inverse trades additional construction for cheaper repeated support lookup.

## Existing repository persistence paths

Milliseconds, median [minimum, maximum]. The reference reducer is test-only.
Concrete, filtered and flag columns request the same diagram on their admitted
inputs; representatives additionally request cycle payloads. None is a
critical-set speedup denominator.

| Fixture | Reference F2 | Concrete diagram | Generic filtered | With cycles | Implicit flag |
| --- | ---: | ---: | ---: | ---: | ---: |
| grid-8 | 0.041 [0.039, 0.048] | 0.359 [0.336, 0.696] | 0.267 [0.251, 0.322] | 0.499 [0.463, 0.602] | 0.192 [0.178, 0.319] |
| grid-16 | 0.187 [0.176, 0.329] | 0.859 [0.768, 1.259] | 0.855 [0.797, 1.160] | 1.777 [1.602, 2.396] | 0.443 [0.394, 0.722] |
| delayed-fan-64 | 0.026 [0.025, 0.027] | 0.333 [0.279, 0.512] | 0.258 [0.220, 0.476] | 0.385 [0.358, 0.443] | excluded |
| delayed-fan-256 | 0.077 [0.075, 0.081] | 0.542 [0.523, 0.695] | 0.496 [0.469, 0.953] | 0.807 [0.713, 0.939] | excluded |
| clique-12 | 0.064 [0.060, 0.149] | 0.390 [0.360, 0.663] | 0.355 [0.324, 0.979] | 0.521 [0.456, 0.673] | 0.163 [0.127, 0.766] |
| clique-20 | 0.425 [0.395, 0.824] | 0.977 [0.865, 1.257] | 1.648 [1.513, 1.778] | 2.424 [2.180, 2.689] | 0.246 [0.206, 0.489] |
| tetra-skeleton-8 | 0.057 [0.048, 0.061] | 0.496 [0.415, 0.704] | 0.406 [0.319, 0.527] | 0.645 [0.479, 0.875] | excluded |
| tetra-skeleton-14 | 5.195 [4.976, 7.326] | 61.024 [58.896, 65.448] | 25.814 [24.459, 27.768] | 68.346 [66.577, 75.371] | excluded |

The implicit flag path avoids explicit boundary work and has lower medians
than the concrete diagram path on all four admitted cases. It cannot represent
delayed triangle values or the requested H2 workload. The generic filtered path
beats concrete diagram time on the large 3-skeleton, but is slower on clique-20.
The internal reference's narrower pairing/result work is useful reducer cost
context; it does not establish equivalent public API, policy or context costs.

## Process memory

Absolute peak process RSS for Q64 on the larger fixtures, MiB, median [min, max]:

| Fixture | Sparse critical | Dense | Oineus Partial U | Oineus full U |
| --- | ---: | ---: | ---: | ---: |
| grid-16 | 4.64 [4.64, 4.64] | 7.50 [7.49, 7.52] | 4.81 [4.64, 4.81] | 4.98 [4.98, 5.16] |
| delayed-fan-256 | 3.78 [3.78, 3.78] | 5.34 [5.29, 5.35] | 4.47 [4.30, 4.47] | 4.64 [4.64, 4.81] |
| clique-20 | 4.64 [4.64, 4.64] | 8.94 [8.94, 8.94] | 4.81 [4.64, 4.81] | 4.98 [4.98, 5.16] |
| tetra-skeleton-14 | 7.91 [7.91, 7.91] | 8.77 [8.77, 8.77] | 5.67 [5.50, 5.67] | 5.67 [5.67, 5.84] |

Prepared-input median RSS spans 2.58-2.92 MiB for these Rust rows and is
3.61 MiB for Oineus; prepared-input HWM is recorded separately. Absolute peaks
cover cold/warm work and validation, including allocator retention and outputs.
These observations support lower process peaks than dense on these four
fixtures, but not a live-allocation bound. The different index/storage widths
and prepared inputs prevent interpreting the table as equal algorithm ownership.

## Evidence and reproduction

All original samples, warmups, plans, fixture hashes, generated original/adapted
sources, binary/header hashes and build logs remain at the ignored local path:

`target/benchmarks/commit-ceea0f4e47e4/critical-sets/run-001/`

`summary.json` contains all 132 cells; `analysis.json` contains paired ratios
and RSS spreads; `oracle/gudhi.json` retains independent diagram results.
`evidence-manifest.json` has SHA-256
`6cda5846f85d489f8e0e8361345d1dc461fffa5be6a25838ea6b0466ed1e66ea`.
The separately retained fixed Oineus archive has SHA-256
`614cc6972ad8b9c3a8edf78bd38560aa08ee37d3c4b4e253d6a07be1bb3ce23a`;
`upstream-verification.json` records its consumed-file checks. These files are
local disposable evidence, without a public retrieval URL or durable retention.

Smoke-001 failed before measurements because WSL Git could not resolve the
Windows worktree path. Smoke-002/003 were adapter exploration; smoke-004
validated the final replayable adapter. All remain under the earlier
`target/benchmarks/commit-a380cba76d41/critical-sets/` directory and are not used
for rankings. The full fixture GUDHI preflight also remains there. No sample
or failed attempt was removed to obtain favorable results.

Run the measured commit in its own checkout. From the original worktree's WSL
directory, this exact replay uses a new output directory and the saved pinned
dense source (other paths refer to the existing local environments):

```sh
LD_LIBRARY_PATH=target/linux-rust/root/usr/lib/x86_64-linux-gnu \
python3 tools/benchmark_critical_sets.py \
  --rustc target/linux-rust/root/usr/bin/rustc-1.91 \
  --oineus-source /mnt/f/cocycle-rs/target/phase4-oineus-source/pinned \
  --boost-include target/oineus-native/boost/usr/include \
  --dense-source target/benchmarks/commit-ceea0f4e47e4/critical-sets/run-001/workers/dense-original.rs \
  --samples 12 --cpu 0 --timeout 60 \
  --output target/benchmarks/commit-ceea0f4e47e4/critical-sets/run-002
```

The [protocol guide](../optimization/README.md#sources-and-build) describes setup
when these local environments or the unpublished dense object are unavailable.
New hardware, compiler, ownership boundaries or fixtures require a new
measurement. This report does not measure a later main/merge revision.
