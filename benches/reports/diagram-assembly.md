# Persistence diagram assembly overhead

[Benchmarks](../README.md) / [Reports](README.md)

## Question and conclusion

This fixed-R0 diagnostic study answers [Issue #17](https://github.com/Aequiludium/cocycle-rs/issues/17):
how much persistence time and output storage belong to raw-to-canonical diagram
assembly? **Recommend a focused H0 materialization/memory follow-up.** Large forest
and isolated-vertex graphs repeatedly exceed the output-heavy engineering gate;
normal Rips H1 assembly remains small. These are synthetic, plausible sparse graph
controls, not application traces or evidence of a project-wide bottleneck.

The recommendation is to investigate reservation and overlapping output storage
first. This study adds diagnostic tools and a report; production Rust, public APIs,
dispatch, checked constructors and canonical ordering are unchanged. No
DiagramBuilder or trusted-constructor design is justified here. A strong
coarse/decomposition disagreement for dense H0 remains unresolved and is excluded
as independent justification for optimization. There is no measured candidate
speedup or RSS saving.

## Measured revision and environment

| Item | Recorded value |
| --- | --- |
| Associated PR | None for this study |
| Kernel | Phase-2 R0 `b2c3bd5eebf0c1193f30ea651012b8ae61b274c1`; identity and mathematical acceptance in [the baseline report](phase2-baseline.md) |
| Harness | Local commit `8c23aacd78a0637437471c4bf2e39d29e75a0324` |
| Source fingerprint | `94d44e0f6773fb4bcad2d93b34f05a2d228e89dc815843b69beecb7917947db9`; original kernel, pipeline worker/controller/helpers and assembly controller |
| Source state | All formal starts/ends clean; source fingerprints unchanged; generated Rust and worker hashes retained separately |
| Protocol | `cocycle-rips-assembly-v1`, separate from pipeline/native baseline timing pools; [frozen diagnostic contract](https://github.com/Aequiludium/cocycle-rs/blob/8c23aacd78a0637437471c4bf2e39d29e75a0324/benches/pipeline/README.md#private-assembly-diagnostics) |
| Runs | `stage-a-001`, `stage-a-002`, conditional `stage-b-001`, adaptive temporal control `stage-a-003` |
| UTC window | 2026-09-30 14:27:57 to 14:36:08, including builds and gaps between runs |
| Machine | Intel Core Ultra 7 155H host; Ubuntu 24.04.4 on WSL2, x86_64, kernel 6.6.87.2, glibc 2.39; guest reports 11 cores / 22 threads |
| Compilers | Rust 1.91.0 `f8297e351`, Cargo 1.91.0 `ea2d97820`; Python 3.12.3 |
| Build | `cargo build --release --locked --offline --lib`; standalone workers `rustc --edition=2024 -O -D warnings`; no features/custom RUSTFLAGS |
| Execution | Serial fresh single-threaded workers on guest CPU 0; OMP/OpenBLAS/MKL thread settings 1 |
| Limits | 60 s per process, 2048 MiB virtual-address-space cap; these are harness limits |
| Uncontrolled | Host contention, frequency, thermals, NUMA and WSL scheduling; CPU 0 physical-core class is unrecorded |

Ordinary and profiled workers share the same public calls and worker extension
for nonzero-born supplied complexes. Instrumentation exists only in generated
source copies under each run's `build/profile-source`. The ordinary worker uses
the unmodified R0 crate. No native reference library is timed by this diagnostic.
R0's existing oracle acceptance remains distinct from per-sample equality here.

## Workloads and comparison contract

The full matrix has 25 workloads: ten existing dense/point H1 cases, three dense
H0 sizes, forests and isolates at 10,000/100,000 vertices, three triangle-free
complete bipartite graphs, three disjoint filled-triangle complexes, and two
expanded generic cases. Dense input uses the existing lower-matrix fixture layout;
the point case calls the public Euclidean point workflow. Existing shared dense
fixtures are quantized once to float32-exact values. Point coordinates and all
Rust computation use f64. Graph weights and supplied birth values are exact dyadic
or integer values. There is no approximation or representative request.

Normal cases use F2 and H0 or H1 as named; the filled triangles use F3/H1, and
expanded nonmetric/sphere cases use F3/H2. The existing H1 cutoff is retained.
The expanded sphere has cutoff 1; other controls have complete coverage. Routes
are recorded from actual execution: H0 union-find, specialized F2/H1, generic
cohomology and boundary reduction. Expanded zero-born complexes use generic
cohomology; vertices born at 1 force the filled-triangle boundary route.

Every warmup and measured process agrees on the full canonical interval sequence,
multiplicity, coverage and available simplex metadata. Independent expectations
also check every output-heavy control: forest edges give one finite H0 bar per
edge and one essential bar; isolates give one essential H0 bar per vertex;
K(s,s) gives 2s-1 finite H0 bars and (s-1)^2 essential H1 bars; each filled triangle
gives two H0 [1,2), one essential H0 born 1, and one H1 [2,3).

## Measurement and sampling

The kernel clock encloses the selected raw producer, raw-bar generation and
cleanup completed before raw output returns. For boundary reduction it stops
immediately before assembly; later reducer-state destruction is outside that clock.
Input/access construction, public dispatch, context assembly and remaining caller
cleanup can form a public-analysis remainder. The assembly clock covers filtering,
endpoint conversion, interval-local validation, fallible final Vec growth, raw
buffer release, diagram validation and canonical sorting. Export/serialization,
descriptors, distance preparation and retained-result destruction are excluded.
Diagnostic stderr emission occurs after all worker clocks and memory readings.

Two distinct sample ratios are retained: assembly/public analysis and
assembly/(kernel + assembly). The first uses the public worker's compute bin;
the second follows the Issue's coarse-phase formula. Fractions are medians of
same-call ratios, not ratios of independent phase medians. Phase medians need
not sum. The parser checks that private clocks fit their enclosing public call.

Each worker/case has one discarded warmup and 12 measured fresh processes,
with balanced ordinary/profile positions. A1/A2 seeds are 1701/2701. Both cover
all 25 cases: 1,200 measured processes. The engineering gate selected 14 cases
for Stage B (seed 3701), with 336 measured processes. The dense-H0 disagreement
then prompted a separately recorded seven-case coarse temporal control (A3,
seed 4701), with 168 measured processes. **Total: 1,704 measured processes plus
142 discarded warmups; all 71 run/case cells passed with no failure or exclusion.**
Timing pools are separate; no mixed Stage A/B aggregate is reported.

## Coarse results

A1/A2 designate the first two complete rounds. All times are milliseconds and
medians. Counts and capacities agree across runs. Workload labels omit the path
suffix; `nonmetric24` and `sphere_h2` are expanded, filled triangles are supplied
boundary complexes, and named forest/isolate/bipartite controls use flag graphs.

| Workload | Raw/final bars | Kernel A2 ms | Assembly A1/A2 ms | Analysis A2 ms | Public % A1/A2 | Kernel+assembly % A1/A2 |
| --- | --- | --- | --- | --- | --- | --- |
| uniform_h1_32 | 39/39 | 0.4475 | 0.0026/0.0032 | 1.2908 | 0.23/0.26 | 0.63/0.73 |
| uniform_h1_64 | 78/78 | 0.8011 | 0.0040/0.0053 | 1.5523 | 0.28/0.33 | 0.56/0.64 |
| uniform_h1_128 | 163/163 | 2.6211 | 0.0130/0.0090 | 3.5396 | 0.26/0.26 | 0.33/0.33 |
| circle_h1_128 | 129/129 | 14.0682 | 0.0040/0.0038 | 15.0390 | 0.03/0.02 | 0.03/0.03 |
| clusters_h1_128 | 147/147 | 3.3262 | 0.0073/0.0122 | 4.5299 | 0.29/0.27 | 0.40/0.38 |
| duplicates_h1_128 | 142/78 | 2.7105 | 0.0060/0.0072 | 3.9203 | 0.13/0.16 | 0.17/0.24 |
| equal_h1_128 | 128/128 | 1.2932 | 0.0021/0.0023 | 2.3361 | 0.11/0.11 | 0.21/0.21 |
| nonmetric_h1_64 | 190/190 | 6.7744 | 0.0116/0.0122 | 7.9008 | 0.14/0.15 | 0.16/0.18 |
| uniform_h1_128_cutoff | 162/162 | 1.0562 | 0.0076/0.0099 | 1.9409 | 0.43/0.49 | 0.85/0.91 |
| uniform_h1_128_points | 163/163 | 2.7018 | 0.0148/0.0142 | 2.9264 | 0.49/0.52 | 0.52/0.57 |
| uniform_h0_128 | 128/128 | 0.5433 | 0.1433/0.1200 | 1.2581 | 9.66/9.30 | 18.18/17.34 |
| uniform_h0_512 | 512/512 | 7.5817 | 0.2564/0.2972 | 8.5159 | 2.72/3.29 | 2.93/3.53 |
| uniform_h0_1024 | 1024/1024 | 33.5673 | 0.3444/0.4408 | 34.7161 | 1.06/1.21 | 1.08/1.23 |
| forest10000_h0 | 10000/10000 | 0.7197 | 0.1170/0.1318 | 1.2273 | 10.35/10.51 | 15.22/15.49 |
| isolates10000_h0 | 10000/10000 | 0.1468 | 0.3009/0.2923 | 1.0027 | 28.53/27.92 | 67.25/67.55 |
| forest100000_h0 | 100000/100000 | 7.4783 | 2.2474/2.4058 | 10.4592 | 23.72/23.18 | 25.05/24.65 |
| isolates100000_h0 | 100000/100000 | 1.7486 | 1.8093/2.4152 | 4.9185 | 48.42/46.30 | 57.87/56.18 |
| bipartite32_h1 | 1025/1025 | 0.6176 | 0.0227/0.0207 | 1.0370 | 1.96/1.92 | 3.13/3.18 |
| bipartite64_h1 | 4097/4097 | 1.7228 | 0.0685/0.0933 | 2.2020 | 3.87/4.24 | 4.98/5.08 |
| bipartite128_h1 | 16385/16385 | 6.9606 | 0.2388/0.2644 | 7.6525 | 3.43/3.58 | 3.61/3.76 |
| filled_triangles16 | 64/64 | 0.1629 | 0.0024/0.0029 | 0.1911 | 1.43/1.55 | 1.66/1.78 |
| filled_triangles128 | 512/512 | 0.2043 | 0.0108/0.0120 | 0.3587 | 3.43/3.50 | 5.38/5.71 |
| filled_triangles1024 | 4096/4096 | 0.4383 | 0.1136/0.0885 | 1.4635 | 6.32/5.93 | 18.62/16.42 |
| nonmetric24 | 2048/61 | 8.1254 | 0.0150/0.0129 | 8.1428 | 0.16/0.17 | 0.16/0.17 |
| sphere_h2 | 14/7 | 0.1396 | 0.0032/0.0041 | 0.1453 | 2.77/2.71 | 2.80/2.73 |

The following complete A2 dispersion/control table retains the observer effect.
Profile/ordinary is the ratio of public-analysis medians; values below one do not
demonstrate a speedup. RSS is each variant's largest individual process high-water
mark in KiB, not allocated bytes or a construction-only peak. All sample ranges,
other phases and both variants remain in the full matrices.

| Workload | A2 assembly min-max ms | A2 public % min-max | A2 profile/ordinary analysis | A2 profile/ordinary max RSS KiB |
| --- | --- | --- | --- | --- |
| uniform_h1_32 | 0.0024-0.0083 | 0.23-0.61 | 1.312 | 2464/2464 |
| uniform_h1_64 | 0.0039-0.0098 | 0.26-0.43 | 1.172 | 2640/2640 |
| uniform_h1_128 | 0.0073-0.0123 | 0.18-0.37 | 1.046 | 2816/2816 |
| circle_h1_128 | 0.0028-0.0066 | 0.02-0.03 | 1.015 | 4048/4048 |
| clusters_h1_128 | 0.0074-0.0162 | 0.22-0.38 | 1.249 | 2816/2816 |
| duplicates_h1_128 | 0.0045-0.0141 | 0.11-0.37 | 1.126 | 2816/2816 |
| equal_h1_128 | 0.0018-0.0038 | 0.08-0.14 | 1.037 | 2816/2816 |
| nonmetric_h1_64 | 0.0086-0.0203 | 0.11-0.27 | 1.030 | 3168/3168 |
| uniform_h1_128_cutoff | 0.0073-0.0402 | 0.47-1.97 | 1.161 | 2816/2816 |
| uniform_h1_128_points | 0.0098-0.0275 | 0.29-1.08 | 1.082 | 3168/3168 |
| uniform_h0_128 | 0.0952-0.2718 | 7.71-20.58 | 0.935 | 2816/2816 |
| uniform_h0_512 | 0.1297-0.4340 | 1.82-5.43 | 0.980 | 8164/8164 |
| uniform_h0_1024 | 0.3402-0.7195 | 0.92-2.19 | 0.969 | 25140/25140 |
| forest10000_h0 | 0.0893-0.2182 | 8.48-14.31 | 0.913 | 3872/3872 |
| isolates10000_h0 | 0.2490-0.4232 | 19.46-30.13 | 1.100 | 3168/3168 |
| forest100000_h0 | 1.7596-4.0737 | 17.40-31.11 | 1.099 | 17108/17112 |
| isolates100000_h0 | 1.6102-3.1460 | 41.20-55.27 | 1.079 | 9952/9948 |
| bipartite32_h1 | 0.0149-0.0225 | 1.57-2.20 | 0.806 | 2640/2640 |
| bipartite64_h1 | 0.0671-0.1450 | 3.22-6.27 | 0.956 | 3168/3168 |
| bipartite128_h1 | 0.2338-0.3346 | 3.07-3.95 | 0.934 | 4572/4524 |
| filled_triangles16 | 0.0023-0.0041 | 1.05-2.07 | 0.521 | 2640/2640 |
| filled_triangles128 | 0.0096-0.0170 | 2.77-4.18 | 0.641 | 2992/2992 |
| filled_triangles1024 | 0.0750-0.0961 | 5.33-6.95 | 0.873 | 5632/5632 |
| nonmetric24 | 0.0118-0.0234 | 0.11-0.26 | 1.035 | 7920/7920 |
| sphere_h2 | 0.0028-0.0052 | 1.75-3.63 | 1.119 | 2464/2464 |

## Conditional decomposition and temporal control

Stage B groups filtering, endpoint conversion, interval-local validation, Vec
growth and raw release into materialization. Per-interval clocks or allocator
hooks were not introduced because they would change the short conversion loop.
Capacity growths count actual successful final-Vec capacity changes, not allocator
calls. Diagram-level validation and the single canonical sort have separate
clocks. Stage B is an explanation of its own generated copy, not a subtraction
from Stage A. Columns below are phase medians in milliseconds.

| Workload | Assembly ms | Materialization ms | Diagram validation ms | Sort ms | Capacity growths | Public % |
| --- | --- | --- | --- | --- | --- | --- |
| uniform_h1_128 | 0.0081 | 0.0016 | 0.0002 | 0.0062 | 7 | 0.26 |
| uniform_h0_128 | 0.0019 | 0.0011 | 0.0001 | 0.0005 | 6 | 0.16 |
| uniform_h0_512 | 0.0076 | 0.0050 | 0.0004 | 0.0021 | 8 | 0.08 |
| uniform_h0_1024 | 0.0142 | 0.0088 | 0.0008 | 0.0041 | 9 | 0.04 |
| forest10000_h0 | 0.1244 | 0.0878 | 0.0058 | 0.0303 | 13 | 10.95 |
| isolates10000_h0 | 0.2048 | 0.1803 | 0.0054 | 0.0162 | 13 | 27.09 |
| forest100000_h0 | 2.5815 | 2.1697 | 0.0963 | 0.3029 | 16 | 25.88 |
| isolates100000_h0 | 1.6111 | 1.3643 | 0.0921 | 0.1767 | 16 | 43.59 |
| bipartite64_h1 | 0.0760 | 0.0674 | 0.0021 | 0.0062 | 12 | 4.74 |
| filled_triangles16 | 0.0036 | 0.0009 | 0.0001 | 0.0025 | 5 | 0.94 |
| filled_triangles128 | 0.0111 | 0.0058 | 0.0003 | 0.0048 | 8 | 2.07 |
| filled_triangles1024 | 0.1412 | 0.0847 | 0.0037 | 0.0525 | 11 | 5.22 |
| nonmetric24 | 0.0139 | 0.0093 | 0.0002 | 0.0042 | 5 | 0.17 |
| sphere_h2 | 0.0013 | 0.0005 | 0.0002 | 0.0003 | 2 | 3.98 |

Materialization dominates the large H0 controls, while diagram validation is
small. Sorting matters for some producers but its absolute normal-H1 cost is
small. The source has one final canonical sort, and no evidence of a duplicate
canonical sort. Interval-local and diagram-level checks enforce different
invariants; these timings do not justify removing either layer.

Dense H0 n=128 shows 0.1200 ms in A2 but 0.0019 ms in B; the discrepancy extends
to n=512/1024. The additional coarse A3 run reproduces the earlier coarse pattern,
while its large-H0 effect also persists:

| Workload | A3 assembly ms | A3 public % | A3 kernel+assembly % | A3 profile/ordinary analysis |
| --- | --- | --- | --- | --- |
| uniform_h1_128 | 0.0127 | 0.26 | 0.33 | 1.112 |
| uniform_h0_128 | 0.1118 | 9.10 | 17.77 | 1.087 |
| uniform_h0_512 | 0.1873 | 2.30 | 2.54 | 0.873 |
| uniform_h0_1024 | 0.3531 | 0.96 | 0.98 | 0.916 |
| forest100000_h0 | 2.5190 | 24.42 | 26.23 | 1.014 |
| isolates100000_h0 | 1.6847 | 44.86 | 55.78 | 0.918 |
| sphere_h2 | 0.0034 | 2.54 | 2.55 | 1.060 |

This repeat supports a measurement-mode-sensitive discrepancy but does not
establish its mechanism. Code generation, allocator behavior and host scheduling
have not been discriminated. An ordinary public-call control cannot localize
this private-phase discrepancy. Dense H0's apparent normal-workload gate is
therefore insufficient on its own. The small sphere is also observer-sensitive:
Stage B profile/ordinary public-analysis medians have ratio 0.231. No precise
production assembly fraction is inferred for these short paths.

## Output storage and construction-path findings

Both raw tuple and final interval are 32 bytes on this measured target. The
following capacity accounting describes overlapping owned output allocations
during conversion. Raw storage is released at the end of the consuming loop,
before final diagram validation/sorting. It excludes allocator metadata,
transient reallocation storage and other live solver/input state. The boundary
reducer's state can remain live through assembly; output overlap is consequently
not whole live-memory peak. RSS cannot establish the saving from a future change.

| Workload | Raw/final capacity | Raw/final KiB | Overlap KiB | Final payload KiB |
| --- | --- | --- | --- | --- |
| uniform_h0_1024 | 1024/1024 | 32.0/32.0 | 64.0 | 32.0 |
| forest10000_h0 | 10000/16384 | 312.5/512.0 | 824.5 | 312.5 |
| isolates10000_h0 | 10000/16384 | 312.5/512.0 | 824.5 | 312.5 |
| forest100000_h0 | 100000/131072 | 3125.0/4096.0 | 7221.0 | 3125.0 |
| isolates100000_h0 | 100000/131072 | 3125.0/4096.0 | 7221.0 | 3125.0 |
| bipartite128_h1 | 32768/32768 | 1024.0/1024.0 | 2048.0 | 512.0 |
| filled_triangles1024 | 4096/4096 | 128.0/128.0 | 256.0 | 128.0 |
| nonmetric24 | 3072/64 | 96.0/2.0 | 98.0 | 1.9 |

For 100,000 output bars, raw capacity contributes 3,125 KiB and final capacity
4,096 KiB: 7,221 KiB overlap, versus 3,125 KiB final payload. Stage B records 16
final capacity growths. The bipartite case overallocates both buffers; generic
nonmetric24 filters 2,048 raw bars to 61, so reserving raw length everywhere could
waste memory. These are concrete opportunities to test, not measured savings.

Raw-to-final materialization and temporary dual storage are confirmed. The final
canonical sort is singular. Interval-local validation covers finite endpoints,
positive lifetime and endpoint rules; diagram validation covers computed dimension
membership and coverage compatibility. Contiguous computed dimensions remain
compact through `ComputedDimensions::through`, without an expanded dimension list.
Producer ordering differs; specialized H1 output cannot be assumed canonical.

## Decision and scoped follow-up

Campaign disposition: the B1-B4 summary adopts no production change. The
focused H0 recommendation below is a diagnostic lead, not a measured optimization
candidate or an implementation delivered by this closeout. The original
dense-H0 uncertainty remains unresolved.

**Outcome B: recommend investigating H0 materialization and construction memory.**
The repeatable large graph controls cross the Issue's output-heavy gate under
both denominators and show millisecond-scale absolute cost. They also expose
substantial temporary output storage and repeated capacity growth. Normal H1
medians remain below 1% under both denominators; retain its current path. Bipartite
H1 and generic/boundary results do not warrant a broad architecture change. The
filled-triangle example is synthetic and below 10% of public analysis, despite a
larger kernel-only fraction. The dense-H0 discrepancy remains a follow-up caveat.

Concrete local implementation proposal: **Reduce H0 raw-to-diagram capacity growth
and overlapping output storage.** First test bounded reservation from a justified
producer count/size hint, taking zero-bar filtering into account. Preserve fallible
allocation, coverage, multiplicity and canonical order. Do not blindly reserve raw
length for high-filtering generic producers. Buffer fusion, trusted validation and
new construction abstractions require additional evidence before consideration.

Acceptance for that follow-up: compare an ordinary candidate against immutable R0
on these exact large-H0 controls, a representative downstream graph and normal
H0/H1 regressions; retain time and memory trade-offs plus all failures; distinguish
capacity savings from RSS; resolve the dense-H0 timing discrepancy before using
it as an acceptance gate. Select a production change only from that comparison.
This report records the local recommendation; remote Issue closure/creation and
submission are pending the user's overall Direction-B review.

## Verification and evidence

Linux Python tool tests: 90 passed. Both 16-case development diagnostic modes
passed, including analytic controls and generated Rust 1.91 compilation with
strict standalone-worker warnings. Source, documentation, local collaboration and
staged-artifact checks passed. Windows tool-suite invocation encountered the existing
Linux `resource` requirement and a platform-path assertion; Linux is the supported
measurement/test platform and passed the full suite. Production Rust did not change;
the expensive fixed-R0 oracle suites were not repeated. R0's hosted CI belongs to
the baseline report; there is no new hosted CI claim for this local harness.

All evidence is **local-only**, with no retention expiry assigned or public upload.
The ignored study root is
`target/benchmarks/commit-b2c3bd5eebf0/phase2-assembly/`. It retains exact fixtures,
every raw process record, hashes, environment and build commands, transformed
source, binaries, validation, matrices and summaries. Development failures and
successful smoke attempts remain separate from formal timing evidence.

Portable archive:
`target/deliverables/phase2-assembly-b2c3bd5eebf0/phase2-assembly-evidence.zip`.
Size: **51,613,290 bytes**. SHA-256:
`62dbc0676ae780f512c013591a2748cfbd083de9239a941a3129251f831cceb5`.
All **3,077 archive members**, **1,846 formal process records** and **444 original
source-file checks** were verified; original Rust files match the pinned R0 Git
archive with checkout line-ending normalization. Kernel/harness Git snapshots
and generated diagnostic sources are included; intermediate Cargo caches are
excluded. The archive inventory verifies each member's length and SHA-256.

From a clean checkout of the recorded local harness, with R0 kernel files and
the recorded toolchain on Linux, reproduce into fresh directories:

```sh
python3 tools/profile_assembly.py --samples 12 --cpu 0 --order-seed 1701 --output target/assembly-replay/a1
python3 tools/profile_assembly.py --samples 12 --cpu 0 --order-seed 2701 --output target/assembly-replay/a2
python3 tools/profile_assembly.py --samples 12 --cpu 0 --order-seed 3701 --decompose --cases uniform_h0 forest isolates bipartite64 filled_triangles uniform_h1_128_dense_lower nonmetric24_expanded sphere_h2_expanded --output target/assembly-replay/b1
python3 tools/profile_assembly.py --samples 12 --cpu 0 --order-seed 4701 --cases uniform_h0 sphere_h2_expanded uniform_h1_128_dense_lower forest100000 isolates100000 --output target/assembly-replay/a3
```

These commands reproduce the contract and fixtures; uncontrolled host conditions
mean another run need not reproduce the observed timings or the unresolved discrepancy.
