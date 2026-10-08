# H1-to-H2 lifetime and scratch evidence

[Benchmarks](../README.md) / [Reporting rules](../reporting.md)

## Question and conclusion

[M1 #47](https://github.com/Aequiludium/cocycle-rs/issues/47) audits the lifetime
contract of the [T3 baseline](h2-foundation.md), before H2 shortcut experiments.
Two complete local comparative runs validate all 49 workloads against native
GUDHI and Ripser. The H1 reducer already drops before H2 allocation: this is a
verified no-change lifetime result. Generic old transformation columns can be
dropped before assembling the next level, but the measured runtime benefit is
unstable. Within-H2 heap reuse reduces capacity-growth events and often time,
but retains capacity and fails four H1 control screens in the first run.
Both changes remain generated experimental ablations; neither is adopted into
normal library dispatch. This resolves M1's lifetime prerequisite for T4;
M2's later generic reuse audit remains dependent on T9 and is not completed here.

The result establishes ownership on the inspected source and correctness on the
tested finite domain. It does not establish allocator live-byte totals, a process
memory bound, universal speedup or admission of T3's private kernel.

## Measured revision and environment

| Item | Recorded value |
| --- | --- |
| Stack | `main <-` [#50](https://github.com/Aequiludium/cocycle-rs/pull/50) `<-` [#51](https://github.com/Aequiludium/cocycle-rs/pull/51) `<- M1` |
| Preceding T3 head | `c5441935250145d0a37874689a0dbe3f4bb13d0e`; algorithm introduced at `59f14f3`, controller fix at `d6042f1895257d99258c61996eb0999f5be0792a` |
| Measured kernel and controller | `16ad170facf711437cca497ffeef71b3dfd6aedc`, clean before/after both runs |
| Joint source SHA-256 | `1162e75ea1d09275cd1ca6dd86e30d885b644ead442a6ea99ca5462d8dccdb7b` |
| Lifetime source SHA-256 | `3b45617fb1f041937e97df8b352829e02b1d5040ddc40bcdc01fccd85503e698` |
| Scratch source SHA-256 | `37a9778004eee44044df030cb02530ffa042d749dfd8b22578c1afaee5383bbe` |
| Controller SHA-256 | `538b9674712c12d618f6163c71d70f787b6c206e97e72168648df3be7a39740e` |
| Native environment SHA-256 | `a01c83c14ecc0688b9e912c427548d5abba0dfa4e7a2ff2b6e8b69de192c3a56` |
| Native pins | GUDHI `cba915e3ab8e1f5b1fe26eb44b407285f7af4e78`; Ripser `01add51ff64aaf40889483260cc5c3b7d0f2a1e7` |
| First run UTC | 2026-10-04 11:47:57.833736 to 11:49:19.281466 |
| Repeat | `run-003`, seed 20261005; first seed 20261004; exact same six Rust executable hashes |
| Repeat UTC | 2026-10-04 11:51:33.711235 to 11:52:54.966419 |
| Protocol | `cocycle-workspace-v1`, unchanged `cocycle-rips-pipeline-v2` workflow latency |
| Machine/build | Core Ultra 7 155H, WSL2 Linux 6.6.87.2/glibc 2.39, Rust 1.91.0; Cargo release, worker `rustc -O`, private `--cfg cocycle_h2_bench` |
| Controls | CPU 0; fresh serial processes, Linux-local executables; frequency and competing host load uncontrolled; 30 s timeout, 2048 MiB address-space limit |

Source fingerprints hash relative paths and bytes of Cargo.toml, Cargo.lock and
all src/**/*.rs. The generated source trees, commands, executable identities,
fixture hashes and inherited native build metadata are retained in the evidence.
Only the library and Rust worker are rebuilt; native workers are hash-verified
copies from the completed T3 pipeline environment. Post-measurement controller
changes retain executables outside temporary storage; they do not change these
kernel sources, sample order, worker protocol or recorded measurements.

## Ownership and last use

| State | Owner, nested storage and references | Last use / release and minimal survivor |
| --- | --- | --- |
| Shared input/access | Borrowed matrix or graph, graph adjacency and original IDs; DenseFlag/SparseFlag own an O(n) index prefix; context/options/Execution remain in caller frames | Access remains available for H2 clique enumeration; it is separate from reducer state and is not counted as released H1 workspace |
| H1 topology and forest | `run_access` owns edge Vec, cycle flags, UnionFind parents/sizes; EdgePosition keys index this live edge level | Pairing and death decoding finish inside `run_access`; all these owners drop on its return |
| H1 reduction | Pivot-owner HashMap, TransformColumn Vec with nested additions Vec, working and transformation heaps; virtual apparent owners are reconstructed rather than stored | The complete death-key extraction includes ordinary, emergent and virtual pairs before reducer return; no position key crosses the boundary |
| Handoff | Owned raw interval Vec and owned `[usize; 3]` original-vertex death tuples | `compute_with_clearing` returns only these two reduction products; conversion temporarily overlaps the death Vec with an owned clearing HashSet, then drops the Vec |
| H2 topology | Temporary edge Vec and triangle Vec; generic Simplex levels additionally own vertex Vec payloads | Edge/triangle construction necessarily overlaps until the triangle level is built; edges drop before reduction starts |
| Private H2 reduction | Owned clearing keys, tuple-owner HashMap, Vec<Vec<usize>> transformations and two tuple/position heaps | Owners and nested transformations remain throughout H2; joint drops the heaps each column, scratch clears logical contents and retains capacity until H2 returns |
| Generic next level | `owners.into_keys()` transfers owned vertex keys to clearing; old Simplex level is the enumeration source; columns own BTreeMap nodes | Old/next topology must coexist during assembly. Previous transformation columns have no use after clear-next extraction; lifetime alone drops them before assembly |
| Generic full fallback | Its outer `compute` frame also owns UnionFind; odd-prime/overflow fallback has different surrounding state | Not silently attributed to released specialized H1 state, and not changed by either M1 ablation |

The high_fill32 trace has 353 owned death keys at `h1_released`, with only 2048
interval-capacity bytes and 12288 handoff-capacity bytes reported there. H1 edge,
cycle, forest, owner and transformation storage has ended before the H2 conversion,
edge/triangle assembly and reduction-start events. At assembly the temporary edge
Vec reports 8192 bytes and triangles 131072 bytes; H2 start reports empty owners,
columns and heaps. The source scopes establish release; unchanged RSS across
these landmarks does not negate release and does not measure allocator retention.

Capacity observations describe the separate **test executable**. In particular,
H1 TransformColumn includes an extra empty test-only reduced-column Vec header;
its reported outer capacity is not the ordinary worker's allocation size. The
implicit route leaves that nested reduced payload empty. Reported production
payloads include the nested additions/vertex/transformation Vec capacities.
Hash-table buckets, BTreeMap node bytes, allocator metadata and tracing temporaries
are not included in reported Vec totals; their entries/slots are recorded when
available. Direct allocator live bytes and allocation/reallocation/free totals
remain unavailable (`null`). A capacity-new/growth event records an observed heap
capacity increase, not a malloc/realloc/free call.

## Workloads, sampling and correctness

The frozen phase3 generator supplies all n16/n32 families: uniform, circle,
clusters, duplicate, sphere/noisy sphere, nonmetric ties, equal clique, exact sparse
low-degree/geometric/clique-heavy/high-fill, plus one shared six-vertex octahedron.
There are 25 H2 cases and 20 H1 controls. Four n16 q3 cases (uniform, equal clique,
low-degree, high-fill) isolate existing generic next-level assembly, without a T9
or cross-type reuse claim. Every workflow is exact F2 and diagram-only. Dense
lower-triangle distances are quantized once to float32-exact values and stored as
f64 by Rust/GUDHI; supplied sparse graphs preserve missing edges and original IDs.
No approximation or representative payload is requested in the timing matrix.

Each case runs five backends: joint, lifetime, scratch, native GUDHI and native
Ripser. One discarded warmup plus 15 measured rounds use fresh processes, with
positions balanced over blocks of five rounds. Both runs pass all 3920 workflow
processes and all 147 separate diagnostic processes: 7840 and 294 respectively,
with no measured failure, timeout, mismatch or exclusion. Full interval multisets
retain multiplicity; Rust coverage and other shared metadata are compared across
all backend anchors. Every trace's emitted interval multiset must also agree.
No statistics are admitted for a failed workflow or trace group.

Latency includes the unchanged pipeline input/construction/expansion/compute/export
boundary and intermediate destruction; fixture I/O/parsing, final transport and
retained-result destruction are excluded as described in the
[pipeline contract](../pipeline/README.md#workflows-and-boundaries). VmHWM below
is the ordinary worker's absolute process peak. Separate event RSS/HWM includes
test instrumentation and must not be pooled with ordinary latency. Native
workers establish correctness here; no cross-library ranking is selected.

Each generated variant first passes all 64 release private unit tests (four
profiling tests ignored), including explicit boundary oracles and complete
ordinary/emergent/virtual death-index tests inherited from T3. New boundary tests
use an octahedral sphere plus an isolate, and two disjoint spheres plus an isolate
with H2 multiplicity two. Cancellation and synthetic recoverable allocation errors
at handoff/assembly/reduction boundaries return Err, then a clean retry restores
the full expected result. Inherited every-work-budget tests cover work exhaustion.
Synthetic failures do not establish universal real-OOM recovery. Normal
odd-prime, representative, cutoff/coverage and fallback tests pass in the complete
suite; these contracts are not performance-tested by this exact-F2 matrix.

Ratios are candidate median / joint median; below one favors the candidate.
Paired 95% intervals use the existing ratio-of-medians helper with 2000 resamples
of matched rounds, seed 20261004 plus case index, with no outlier deletion.
Predeclared screens are H1 median <=1.05 and upper95 <=1.10, H2 median <=1.10.
They are regression screens, not sufficient production-admission criteria.
The tables retain every H2 family, every failed H1 screen, selected large/small
runtime spreads and all n32 H2 heap observations. The complete 49-case matrices,
including both native references and all lifetime contrasts, are in results.json
and analysis.json in the local evidence.
## H2 scratch runtime ratios

| Family | n16 first / repeat | n32 first / repeat |
| --- | --- | --- |
| uniform | 0.952 / 0.899 | 0.918 / 0.894 |
| circle | 0.855 / 0.866 | 0.929 / 0.865 |
| clusters | 0.994 / 0.906 | 0.866 / 0.896 |
| duplicates | 0.873 / 0.789 | 0.908 / 0.897 |
| sphere | 0.904 / 0.864 | 0.894 / 0.877 |
| noisy_sphere | 0.905 / 0.867 | 0.861 / 0.903 |
| nonmetric_ties | 0.909 / 0.914 | 0.928 / 0.920 |
| equal_clique | 0.869 / 0.882 | 0.870 / 0.902 |
| low_degree | 1.025 / 0.827 | 1.011 / 0.845 |
| geometric | 0.937 / 0.884 | 0.858 / 0.882 |
| clique_heavy | 0.949 / 0.867 | 0.877 / 0.857 |
| high_fill | 0.938 / 0.931 | 0.838 / 0.898 |
| octahedron (n6 once) | 1.094 / 0.879 | same case |

## Every failed H1 screen in the first run

| Variant / case | Ratio [paired 95% interval] first | Repeat |
| --- | --- | --- |
| scratch / uniform16_h1_dense | 1.078 [1.025, 1.162] | 0.898 [0.866, 0.969] |
| scratch / clusters16_h1_dense | 1.048 [1.013, 1.101] | 0.917 [0.843, 0.977] |
| scratch / sphere16_h1_dense | 1.051 [1.019, 1.117] | 0.898 [0.837, 1.066] |
| lifetime / low_degree16_h1_flag | 1.027 [0.966, 1.202] | 0.927 [0.871, 1.043] |
| scratch / clusters32_h1_dense | 1.012 [0.983, 1.103] | 0.936 [0.896, 1.009] |

## Runtime spread and ordinary process peaks

| Case / run | Joint median [min, max] ms | Lifetime median [min, max] ms | Scratch median [min, max] ms | Max VmHWM KiB joint / lifetime / scratch |
| --- | --- | --- | --- | --- |
| high_fill32_flag / run-001 | 4.3791 [3.1734, 7.2323] | 4.3453 [3.3795, 5.5518] | 3.6700 [2.8781, 5.0408] | 3344 / 3344 / 3168 |
| high_fill32_flag / run-003 | 3.5686 [3.2812, 4.5173] | 3.5676 [3.3078, 5.3100] | 3.2060 [2.7384, 3.6982] | 3344 / 3344 / 3168 |
| circle32_dense / run-001 | 7.7443 [7.2979, 10.6075] | 7.9525 [7.1758, 11.4631] | 7.1951 [6.4648, 9.3573] | 3696 / 3696 / 3696 |
| circle32_dense / run-003 | 7.4873 [6.8644, 7.7292] | 7.4447 [6.9264, 8.6777] | 6.4797 [5.9960, 7.3811] | 3696 / 3696 / 3696 |
| low_degree16_flag / run-001 | 0.0440 [0.0405, 0.0550] | 0.0447 [0.0431, 0.0797] | 0.0451 [0.0313, 0.0492] | 2640 / 2640 / 2640 |
| low_degree16_flag / run-003 | 0.0446 [0.0361, 0.0760] | 0.0413 [0.0298, 0.0529] | 0.0369 [0.0272, 0.0442] | 2640 / 2640 / 2640 |
| octahedron_dense / run-001 | 0.0502 [0.0412, 0.0553] | 0.0496 [0.0466, 0.0562] | 0.0549 [0.0476, 0.0599] | 2640 / 2640 / 2640 |
| octahedron_dense / run-003 | 0.0548 [0.0470, 0.0812] | 0.0511 [0.0438, 0.0672] | 0.0482 [0.0409, 0.0517] | 2640 / 2640 / 2640 |
| uniform16_h3_dense / run-001 | 0.6366 [0.5883, 0.8391] | 0.6531 [0.6025, 1.5327] | 0.6344 [0.5993, 1.3303] | 2816 / 2816 / 2816 |
| uniform16_h3_dense / run-003 | 0.6351 [0.6109, 0.8348] | 0.6605 [0.6250, 1.1839] | 0.6151 [0.5894, 0.7455] | 2816 / 2816 / 2816 |
| high_fill16_h3_flag / run-001 | 0.5251 [0.5061, 0.6274] | 0.5273 [0.5001, 1.0646] | 0.5468 [0.4829, 0.8761] | 2816 / 2816 / 2816 |
| high_fill16_h3_flag / run-003 | 0.5386 [0.4967, 0.9476] | 0.5341 [0.5208, 0.6778] | 0.5492 [0.5015, 0.9871] | 2816 / 2816 / 2816 |

## H2 heap capacity events and retention, n32

| Family | New-capacity events joint / scratch | Growth events joint / scratch | Scratch retained heap capacity bytes at H2 end | Sampled diagnostic max RSS KiB joint / scratch |
| --- | --- | --- | --- | --- |
| uniform | 2984 / 2 | 3078 / 6 | 5184 | 4400 / 4224 |
| circle | 8990 / 2 | 13681 / 9 | 20608 | 5104 / 4928 |
| clusters | 3630 / 2 | 4207 / 4 | 2592 | 4752 / 4400 |
| duplicates | 5022 / 2 | 8221 / 5 | 5152 | 4576 / 4576 |
| sphere | 6098 / 2 | 8950 / 14 | 82944 | 4752 / 4752 |
| noisy_sphere | 5942 / 2 | 8686 / 14 | 82944 | 4752 / 4576 |
| nonmetric_ties | 8990 / 2 | 14026 / 19 | 659456 | 5280 / 5456 |
| equal_clique | 8990 / 2 | 13485 / 3 | 1312 | 4928 / 5104 |
| low_degree | 0 / 0 | 0 / 0 | 0 | 4048 / 3872 |
| geometric | 1196 / 2 | 857 / 7 | 5248 | 4400 / 4048 |
| clique_heavy | 280 / 2 | 155 / 4 | 1344 | 4048 / 4224 |
| high_fill | 3868 / 2 | 4298 / 13 | 41984 | 4752 / 4576 |

## Generic assembly lifetime ablation (q3, n16)

| Family | Old / next level vector capacity bytes | Removed outer-column capacity bytes / tree entries | Lifetime runtime ratio first / repeat | Joint / lifetime diagnostic max RSS KiB first |
| --- | --- | --- | --- | --- |
| uniform | 15232 / 31168 | 6144 / 149 | 1.026 / 1.040 | 4224 / 4224 |
| equal_clique | 50688 / 152896 | 12288 / 455 | 0.998 / 1.018 | 4928 / 4752 |
| low_degree | 960 / 0 | 0 / 0 | 1.009 / 0.881 | 3872 / 4048 |
| high_fill | 14304 / 18176 | 3072 / 175 | 1.004 / 0.992 | 4224 / 4224 |

## Interpretation and limits

H1-to-H2 release needs no lifetime fix. The generic assembly table shows precisely
which old columns can stop overlapping, while both topology levels and clearing
keys remain necessary. Dropping columns removes their outer Vec capacity and
BTreeMap entries at assembly; node bytes are unknown. It does not necessarily
reduce the maximum reported vector capacity or process peak. Uniform q3 runtime
regresses in both runs; equal-clique changes sign. The low-degree control has no
columns to release, yet changes time substantially, so its favorable repeat is
not evidence of a release benefit.

H2 scratch reuse preserves all owners, clearing and transformations, and resets
both heaps before each noncleared column. New/growth counts in its table subtract
the reduction-start counters, excluding preceding H1 and preparation events.
The low-degree fixture has no H2 scratch work. Other n32 families reduce initial
heap-capacity events to two, but retaining the largest heaps raises end capacity:
nonmetric_ties retains 659456 bytes and increases diagnostic max RSS from 5280
to 5456 KiB in the first run. At high_fill32, reported vector-capacity maximum
rises from 288736 to 291968 bytes despite a lower observed process RSS.
Owners and nested stored-column payloads are identical across scratch/joint.
`clear()` is logical reset, not freeing or cross-type capacity transfer.

The first run fails four scratch H1 controls and one lifetime control; the
independent repeat passes them. H1 code is untouched by these ablations, and the
repeat also shifts unchanged control routes. Binary layout and uncontrolled
frequency/load are not isolated. Preserve this contradictory evidence instead
of claiming H1 improvement or clearing the failed first-run screen. Decision:
**experimental, no default adoption** for both ablations. The verified existing
H1 lifetime is the M1 no-change result; M2 remains a later task.

## Evidence and reproduction

All evidence is **local-only** under
`target/benchmarks/commit-16ad170facf7/workspace/`. The successful first/repeat
runs are run-001 and run-003. The retained archive `m1-evidence.tar.gz` has 774
members, 23538179 bytes and SHA-256
`dc72abd28be206c351918bcf5acaa99b0da6256b5cf15c56240723263682c975`.
It contains a member checksum manifest, frozen controller/worker sources,
generated variants, fixtures, every raw sample/trace, summaries/bootstrap analysis,
build logs, reference metadata and all measured executable copies. No durable
public archive is configured; this report does not make the original data
publicly retrievable.

Two earlier build attempts are retained under commit-72dd1333eed1 and
commit-7d5668ce5855: a missing declared bench file and DrvFS executable detection,
respectively. Neither reached sampling. The measured revision's run-002 stopped
before sampling because temporary binaries were absent. Before run-003, the
original rustc command relinked the retained unchanged rlib and the stored test
executables were copied; **all six hashes matched run-001**. The recovery script
and binary-recovery.json are archived. No sample was discarded or replaced.
The final controller automatically persists executables before measurement and
reads those copies for repetition; this retention fix follows the measured run.

With the pinned [native setup](../native/README.md), a clean source checkout and
its completed pipeline environment, reproduce the experiment using the controller:

```sh
python3 tools/benchmark_workspace.py --native-environment path/to/environment.json \
  --samples 15 --cpu 0 \
  --output target/benchmarks/commit-<sha12>/workspace/run-001
python3 tools/benchmark_workspace.py --native-environment path/to/environment.json \
  --samples 15 --cpu 0 --order-seed 20261005 \
  --reuse-build target/benchmarks/commit-<sha12>/workspace/run-001 \
  --output target/benchmarks/commit-<sha12>/workspace/run-002
```

Reuse the first build-bearing directory, not a repeat directory. Fresh runs
record their actual SHA and controller hash. Repetition also requires identical
current source and full harness fingerprints, in addition to the clean checkout,
commit, reference and retained-executable checks. A failed diagnostic keeps its
logs and process samples in `results.json`, writes a failed `summary.json` and
withholds every backend's ranking for that case. Timeout, missing-event,
interval-mismatch and build-identity regression tests exercise these paths
without running native binaries; they do not add performance measurements.
The archived measured controller
and recovery script reproduce the original artifact lifecycle; the final
retention fix does not relabel those original measurements. Local quality/MSRV,
examples, guide doctests and 91 Python tests passed. The local collaboration
checker passes its source checks but still reports missing historical local
links; tracked documentation checks pass. Hosted CI and merging remain separate
PR evidence, and no PR was merged by this experiment.
