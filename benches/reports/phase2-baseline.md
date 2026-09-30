# Integrated Phase-2 performance baseline

[Benchmarks](../README.md) / [Reports](README.md)

## Question and conclusion

This study establishes one immutable resource baseline for the prepared Rips,
complete Rips pipeline and diagram-distance research in
[Issue #16](https://github.com/Aequiludium/cocycle-rs/issues/16). All three suites
passed on the same clean integrated source. They retain **549 performance rows,
6,588 measured processes and 216 separate logical-facade diagnostic calls**.
Their timing contracts remain separate. This is a baseline resource study with
native cross-library controls, not a before/after optimization or route selection.
Substantial sample variation and uncontrolled WSL/host conditions limit rankings.

**Phase-2 R0 = `b2c3bd5eebf0c1193f30ea651012b8ae61b274c1`.**

The user authorized freezing the completed local Issue #14 implementation into
this Direction-B research snapshot while its production disposition remains under
discussion. R0 is an integrated research revision, not a claim that #14 or this
revision has been merged or released. Later boundary changes must be measured as
candidates against this R0; they must not overwrite or relabel its evidence.
Existing historical Rips and preparation reports retain their original identities.
No new production algorithm, policy or default was selected by this study.

## Measured revision and environment

| Item | Recorded value |
| --- | --- |
| Associated PR | None for this study |
| Kernel and harness | Both [`b2c3bd5eebf0c1193f30ea651012b8ae61b274c1`](https://github.com/Aequiludium/cocycle-rs/commit/b2c3bd5eebf0c1193f30ea651012b8ae61b274c1); clean integrated merge snapshot |
| Integrated main parent | `9e6715f4c2e0118ab738486bd747ec3e361397de` |
| Local Issue #14 parent | `5c4f7b87995141d1a95e876bc894cd89f616f0ac` |
| Before/after candidate | None; this run establishes R0 itself |
| Run identity | `commit-b2c3bd5eebf0/phase2-baseline`; native, pipeline and distance each use `run-001` |
| UTC execution window | 2026-09-30 07:25:46 to 09:46:17, including builds, reference checks and final diagnostics |
| Host | Intel Core Ultra 7 155H, 16 physical cores / 22 logical CPUs; 33,779,150,848 bytes RAM; Windows 11 build 26100 |
| Measurement guest | Ubuntu 24.04.4 LTS, x86_64, WSL2 kernel `6.6.87.2-microsoft-standard-WSL2`, glibc 2.39; 16,093,844 KiB guest RAM |
| Guest CPU topology | WSL reports 11 cores / 22 threads; this is distinct from the host's physical topology |
| Toolchains | rustc 1.91.0 `f8297e351`, Cargo 1.91.0 `ea2d97820`, g++ 13.3.0, Python controller 3.12.3 |
| Rust build | `cargo build --release --locked --offline --lib`; optimized standalone workers; no Cargo features or custom `RUSTFLAGS` |
| Distance instrumentation | Private `--cfg cocycle_distance_bench`, `-C opt-level=3`; separately linked ordinary public crate supplies correctness results |
| C++ build | C++17 `-O3 -DNDEBUG` for Rips and corrected GUDHI bottleneck; Topp C++20 `-O3 -DNDEBUG -pthread`, portable scalar build; no custom `CXXFLAGS` |
| Headers and adapters | Boost 1.83; consumed header, transformed adapter and executable SHA-256 inventories preserved in build metadata |
| Execution | Serial fresh worker processes on guest CPU 0; OMP/OpenBLAS/MKL thread settings all 1; measured Topp variants use its serial override |
| Limits | 2048 MiB `RLIMIT_AS`; native/distance 60 s and pipeline 30 s process timeout; these are harness limits, not library quotas |
| Uncontrolled factors | Host contention, CPU frequency, thermals, NUMA placement and WSL scheduling; guest CPU 0 is not a recorded physical-core class |

The integrated source inventory covers tracked `src/`, Cargo manifests, `benches/`
and `tools/`, with per-file SHA-256 values in `identity.json`. Its aggregate hash is
`65414fc37b6132143518366e0033c9489e05c2a19f4a6c52986c6e72f4db452f`.
The local orchestration hash is
`7231cda80317b0a831e7be6874a329e11cfba214c949fcaf84351f4d002b50a8`.
Every suite confirmed unchanged measured sources. The report is written after
measurement; its containing revision is not a new measured kernel revision.

| Suite | Worker/controller source SHA-256 |
| --- | --- |
| Prepared Rips | `a2511c455a0225fb6dbf0631b094a7d6d74a125638db9c5a8c89cec2645d2dca` |
| Rips pipeline | `39501e6b39d03713d2c2e726e1d98a8ad62e6f72c3eb13c8f00abace21ff34fc` |
| Diagram distances | `f384a82e9ad0a0538dad14be56da82060881fb589062d5ca05a7fc88f5f64d02` |

Reference pins are the actual sources consumed by the adapters:

| Reference | Revision or package version |
| --- | --- |
| GUDHI Rips/CAM/collapse | `cba915e3ab8e1f5b1fe26eb44b407285f7af4e78` |
| Upstream Ripser | `01add51ff64aaf40889483260cc5c3b7d0f2a1e7` |
| Topp | `ffa1da051ca7ac5e313c74cc9fb92a2bcb20c234`, version 1.0.1 |
| Corrected native GUDHI bottleneck | `4ec34ac55e6d2e8cfd7c322e84c1b0a56d516d51`, `e=0`, `CGAL_DISABLE_GMP=1` |
| Python correctness references | GUDHI 3.11.0, NumPy 2.4.6, POT 0.9.6.post1; Hera bottleneck `delta=0` in the supported correctness suite |

Formal distance references used an ext4 copy of the existing pinned Python venv
to reduce NTFS import overhead. The identity records its source, snapshot and
matching SHA-256 values for every copied native extension. No package was upgraded.
Python reference latency is not a native performance comparison. The known old
GUDHI default bottleneck and extreme-float disagreements remain documented in the
[distance protocol](../distances/README.md#reference-backends-and-numerical-limits);
this study uses its supported acceptance domain and repaired reference.

## Workloads and measurement contracts

The existing corpus and [reporting rules](../reporting.md) were reused. Fixtures,
full hashes, input layouts, cutoffs, dimensions, fields, requested outputs and
validation records are preserved, rather than inferred from case names.

| Suite | Corpus and timing boundary | Sampling |
| --- | --- | --- |
| [Prepared Rips](../protocol.md) `cocycle-native-v1` | Ordinary F2 H0/H1, shared float32-exact matrices; construction, computation, normalization and temporary cleanup after native input preparation | 12 measured processes per performance cell, no warmup; backend order seed 20260920, input seed 1729; semantic controls once |
| [Complete pipeline](../pipeline/README.md) `cocycle-rips-pipeline-v2` | Public validation/conversion, construction, optional expansion, compute including assembly, and interval export | One discarded fresh-process warmup and 12 measured processes; balanced rounds with order seed 2401 plus case index |
| [Distances](../distances/README.md) `cocycle-distance-v1` | Raw finite f64 endpoints; validation, algorithm preparation, candidates/matching, solve and temporary cleanup in instrumented native workers | One discarded warmup and 12 measured processes; balanced order seed 2401 plus fixture index; only adaptive `baseline` group |

All internal clocks exclude process startup, fixture I/O and final metrics
transport. Results/raw inputs remain alive as their protocols specify. Final
result destruction is not included in the Rips clocks. Distance public-facade
correctness calls and later `logical` diagnostics include raw-to-diagram
construction; they are not ordinary public-API latency measurements or timing rows.

Prepared Rips covers uniform H1 at 32/64/128 vertices, circle/clustered/duplicate/
equal-weight H1 at 128, nonmetric H1 at 64, a cutoff control, and uniform H0 at
128/512/1024. Native controls are direct GUDHI CAM, one-pass GUDHI collapse plus
CAM, and upstream Ripser. Shared values preserve exact interval comparisons;
Ripser's f32 representation and the f64 Rust/GUDHI representations remain distinct.

The 25 pipeline cases cover lower/upper/square dense matrices, threshold and
explicit complexes, F3, nonmetric H2, a cross-polytope H2 sphere, representatives,
bipartite graphs including 10,000 vertices with isolates, H0 controls, sparse Rips
approximation at 32/64/128, exhaustive metric checking and public point construction.
Approximation uses epsilon 0.5, dyadic Manhattan metrics and fixed initial vertex
zero; native greedy choices and permutations are checked independently.
Representative output is retained, but only intervals are exported by this timer.

Distance covers nine families (uniform, clustered, near-diagonal, duplicates,
imbalanced, separated, threshold shell, dense and sparse), sizes 8/32/128/512,
tuning seed 20260922 and holdout seed 20260923. There are 72 distinct pair fixtures
and 216 fixture/metric cells. Imbalance uses a smaller second operand; actual
operand counts are retained. Metrics are Bottleneck with L-infinity cost, W1 with
L-infinity cost, and W2 with Euclidean cost and its final square root. All measured
pairs are finite and off-diagonal. Public Rust, unrestricted Topp and repaired
GUDHI/POT references validate each cell; tiny cases also use the independent
partial-injection oracle. Topp uses compiler-dependent `long double` weighted
arithmetic; Rust uses f64. No tolerance was widened after observing results.

## Coverage and results

Every performance row has exactly 12 retained measured samples. Medians use
`statistics.median`; spread is the observed minimum/maximum, with no trimming.
RSS values below are the maximum absolute process `VmHWM` across measured samples,
in KiB. Raw `VmRSS`, prepared/pre-workflow `VmHWM`, individual peaks and HWM growth
are retained separately; the complete matrix includes these summaries and pipeline
phase statistics. They are process measurements, not live algorithm allocations.

| Suite | Validated coverage | Performance rows | Measured processes | Exclusions/failures |
| --- | --- | ---: | ---: | --- |
| Prepared Rips | 32 cases: 20 semantic, 12 performance; 112 validated comparisons | 48 | 576 | Four Ripser empty/singleton adapter exclusions; no failed case |
| Rips pipeline | 25 cases, 69 executed backend cells | 69 | 828 | Six Ripser sparse-approximation exclusions; no failed case |
| Distances | 216 fixture/metric cells, two timed native workers each | 432 | 5184 | No exclusion, crash, timeout, incomplete sample group or mismatch |
| Logical facade diagnostics | 216/216 public scalar matches | Not timing rows | Not performance samples | No mismatch |

Of the 69 pipeline rows, **15 are `correctness_reference_only`** because native
adapters do not perform equivalent representative extraction, point evaluation,
matrix-layout conversion, explicit-complex retention or exhaustive metric checking.
Their observations remain in the full matrix with that label; they are not eligible
for workflow speedup comparisons. Capability exclusions have no synthetic timings.
The distance controller's unrequested ablation `selection` entries are `incomplete`
by design: baseline-only coverage does not select quickselect, binary or arena policies.

The following tables display Rust observations. Native controls and every measured
row are in the local artifact's `phase2-baseline/matrix.csv`. This presentation makes
no cross-library ranking or aggregate speedup claim.

### Prepared Rips

| Case | Median ms | Min-max ms | Peak KiB |
| --- | ---: | --- | ---: |
| uniform H1 32 | 2.031 | 1.167-3.759 | 2288 |
| uniform H1 64 | 2.982 | 1.724-5.154 | 2464 |
| uniform H1 128 | 6.185 | 4.946-9.418 | 2640 |
| circle H1 128 | 14.500 | 11.655-42.287 | 3696 |
| clustered H1 128 | 2.324 | 1.930-6.736 | 2640 |
| duplicates H1 128 | 5.509 | 4.995-12.598 | 2640 |
| equal-weight H1 128 | 2.452 | 1.654-4.270 | 2640 |
| nonmetric H1 64 | 11.221 | 9.007-15.843 | 2992 |
| uniform H1 128 cutoff | 2.828 | 1.937-7.463 | 2464 |
| uniform H0 128 | 2.039 | 1.381-4.744 | 2640 |
| uniform H0 512 | 15.522 | 11.599-21.134 | 7152 |
| uniform H0 1024 | 57.272 | 26.877-78.483 | 21840 |

### Complete Rips pipeline

| Case | Median ms | Min-max ms | Peak KiB |
| --- | ---: | --- | ---: |
| circle64 dense lower | 8.934 | 6.386-19.986 | 2816 |
| circle64 dense upper | 9.034 | 7.172-12.227 | 2816 |
| circle64 dense square | 7.643 | 5.159-29.300 | 2816 |
| circle64 cutoff threshold | 5.005 | 4.372-7.638 | 2640 |
| circle64 cutoff expanded | 7.302 | 5.024-10.315 | 2816 |
| circle64 cutoff F3 threshold | 3.876 | 2.964-40.944 | 2640 |
| nonmetric24 dense | 16.195 | 12.765-22.191 | 3344 |
| nonmetric24 expanded | 30.028 | 18.824-61.494 | 7920 |
| sphere H2 threshold | 1.276 | 0.980-1.750 | 2464 |
| sphere H2 expanded | 1.370 | 1.085-1.752 | 2464 |
| sphere H2 threshold bases | 1.302 | 1.071-2.236 | 2464 |
| circle64 cutoff threshold bases | 21.230 | 12.747-107.407 | 3168 |
| bipartite16 flag | 7.082 | 4.449-49.493 | 2640 |
| bipartite16 F3 flag | 7.582 | 5.434-11.484 | 2640 |
| bipartite with isolates | 6.387 | 3.933-10.076 | 3520 |
| uniform32 H0 dense | 2.266 | 1.588-6.681 | 2464 |
| metric32 approximate | 5.260 | 3.387-8.162 | 2816 |
| uniform64 H0 dense | 2.751 | 1.876-4.423 | 2640 |
| metric64 approximate | 13.531 | 8.696-46.615 | 3168 |
| metric64 approximate checked | 12.996 | 8.758-19.251 | 3168 |
| metric64 approximate expanded | 14.207 | 11.008-22.948 | 4752 |
| metric64 approximate bases | 32.395 | 23.946-120.077 | 8096 |
| uniform128 H0 dense | 6.072 | 3.081-9.187 | 2816 |
| metric128 approximate | 21.332 | 14.430-27.246 | 3872 |
| line points | 2.178 | 1.705-7.033 | 2640 |

Representative phase medians, in ms, illustrate existing instrumentation:

| Workflow | Input | Construction | Expansion | Compute | Export |
| --- | ---: | ---: | ---: | ---: | ---: |
| circle64 dense lower | 0.922 | 0.000 | 0.000 | 6.755 | 0.018 |
| nonmetric24 expanded | 0.757 | 1.839 | 13.858 | 11.354 | 0.645 |
| metric64 approximate checked | 0.646 | 3.441 | 0.000 | 7.500 | 0.042 |
| metric64 approximate bases | 1.249 | 2.541 | 0.000 | 27.293 | 0.083 |

Compute includes diagram assembly and requested bases. Phase medians need not
sum to the total median. These rows do not establish assembly's independent cost;
that is the separate question in Issue #17.

### Diagram distances

This selected table contains **all holdout families at nominal size 512**, for
all three metrics. The full matrix also contains tuning and sizes 8/32/128.
Each cell is `median [min-max] ms; maximum peak KiB`, from 12 measured processes.

| Family | Bottleneck | W1 | W2 |
| --- | --- | --- | --- |
| uniform | 29.027 [25.982-36.643]; 2640 | 139.997 [118.203-219.360]; 5280 | 199.550 [149.532-248.098]; 5104 |
| clustered | 104.873 [58.989-115.700]; 2464 | 36.566 [31.970-43.174]; 5456 | 78.222 [56.073-112.448]; 5456 |
| near-diagonal | 2.898 [2.435-3.532]; 2640 | 5.364 [3.551-7.069]; 2464 | 3.087 [2.673-4.071]; 2464 |
| duplicates | 1.644 [0.845-3.773]; 2464 | 1.386 [0.791-7.066]; 2464 | 1.002 [0.681-2.667]; 2464 |
| imbalanced | 11.365 [9.450-24.983]; 2640 | 9.143 [7.108-11.074]; 2816 | 8.677 [7.591-11.672]; 2816 |
| separated | 0.841 [0.597-1.520]; 2464 | 1.787 [1.234-2.609]; 2464 | 2.243 [1.594-6.571]; 2464 |
| threshold shell | 5.473 [4.295-8.684]; 4576 | 338.415 [267.836-453.760]; 8624 | 379.255 [309.012-467.470]; 8624 |
| dense | 25.850 [21.519-38.378]; 2640 | 104.421 [97.681-111.553]; 8624 | 591.549 [369.326-658.771]; 8624 |
| sparse | 2.451 [1.544-4.033]; 2640 | 2.853 [1.951-3.269]; 2640 | 2.281 [1.503-4.200]; 2640 |

## Capacity and route diagnostics

Existing distance counters are retained in every raw sample and in the matrix's
`stats_ranges`. Bottleneck exercised `MandatorySparse`, `Multiplicity`, `NoCross`,
`Quickselect` and `Refinement`; 18/72 cells observed matching reuse and 51/72
scratch reuse. W1 observed dense solves in 27/72 cells and sparse solves in 13/72;
W2 observed them in 30/72 and 12/72 respectively. These are counts of cells with a
positive counter, not mutually exclusive route partitions. No direct-cost fallback
was observed in these W1/W2 cells. The Bottleneck `route` field is not a weighted
solver classifier; weighted observations use their own counters.

After all timed suites, each of the 216 fixtures/metrics was called once through
the instrumented `logical` facade and its scalar checked exactly against the
ordinary public-facade reference. Capacity accounting was:

| Metric | Diagnostic cells | Boundary owned bytes/buffers | Algorithm preparation capacity bytes | Nonempty preparation buffers |
| --- | ---: | --- | --- | ---: |
| Bottleneck | 72 | 0 / 0 throughout | 480-49152 | 10 |
| W1 | 72 | 0 / 0 throughout | 288-32768 | 2 |
| W2 | 72 | 0 / 0 throughout | 288-32768 | 2 |

These counters exclude the worker's raw-to-diagram construction. They count
retained Vec capacity and nonempty owned buffers, not allocation calls, RSS or
simultaneous peaks. Bottleneck workspace accounting overlaps preparation indexes;
different container maxima cannot be added into a peak-memory estimate. Zero
facade materialization does not establish a zero-copy pipeline or a speedup.
Rips records existing execution-path metadata, phases and process memory; this
baseline does not add a new Rips workspace or retained-state instrumentation project.

## Verification and interpretation

The result/context ownership review used the integrated source: `PersistenceData`
owns diagram and mathematical context, `PersistenceResult` adds owned optional
representatives, and coverage/computed dimensions belong to the diagram. Distance
consumers borrow logical views and create a fresh work budget for a controlled
terminal. This study did not change these boundaries.

Local Windows verification passed 30 commands, including formatting, Clippy,
debug/release all-features, example checks, strict rustdoc, Markdown doctests,
source/document checks and packaging. Linux Rust 1.91 debug/release all-features,
all-targets and construction-example checks passed; Python tool tests passed 86/86.
Independent full acceptance comparisons passed:

| Check | Result |
| --- | --- |
| Exact Rips against pinned GUDHI/Ripser | 588 cases, 1108 comparisons, 39 protocol checks; 68 explicit exclusions |
| Sparse Rips against pinned GUDHI | 156 topology cases, 153 persistence comparisons, 145 original metric-sampling checks; Ripser capability excluded |
| Supported diagram distances | 108/108 oracle cells, no failure; separate stress cases outside this acceptance run |

[Hosted CI run 36679172098](https://github.com/Aequiludium/cocycle-rs/actions/runs/36679172098)
passed all nine jobs for the exact R0 SHA, including Windows/macOS/Linux, MSRV and
package checks. Hosted native smoke evidence is separate from the local formal
baseline; it is not the source of the performance numbers above. The hosted run
does not establish merge or release status.

The preserved spread is material: circle64 F3 threshold ranges from 2.964 to
40.944 ms, and approximate representative extraction from 23.946 to 120.077 ms.
The checked/unchecked approximation medians also overlap in their distributions;
their ordering does not establish that checking is cheaper. There is one formal
attempt per suite and no independent second round or controlled physical-core
mapping. Future candidates need matched fixtures/protocols and contemporaneous
R0 replays, with repeated ordering where a ranking is claimed. Historical results
on another host cannot be used as a same-machine architectural regression audit.

R0 now supplies a replayable reference for Issue #17 assembly measurement, #18
lifecycle/reuse and #19 concurrency/resource-routing research. Those studies must
retain this identity and report their own candidate revisions and changed timing
boundaries. No reuse, concurrency or routing experiment was performed here.

## Evidence and reproduction

Evidence is **local-only**; no durable external upload or public retrieval is
configured for this formal run. The ignored archive is
`target/deliverables/phase2-r0-b2c3bd5eebf0/phase2-r0-evidence.zip` (14,968,263 bytes).
Its SHA-256 is:

```text
c2ccafa7ba9ed1d85849af83a8504620982414a4fc70b9a6322dd2aaf1d0332e
```

All 1043 inventoried archive members were read back and matched their hashes.
The archive contains the R0 source tar, fixtures, all raw samples/warmups,
validation outcomes/logs, environment and source/header/binary fingerprints,
full matrix, logical diagnostics and the five local orchestration/analysis scripts.
Cargo build caches and reference installations are excluded. `SHA256SUMS.txt` and
`archive-inventory.json` accompany the ZIP. `phase2-baseline/evidence-manifest.json`
provides the run inventory; `coverage.json` records coverage/counter summaries.

The original formal output remains at
`target/benchmarks/commit-b2c3bd5eebf0/phase2-baseline/`, and independent verification
at `target/issue16-validation/`. These are disposable ignored local paths, not
repository links or permanent public evidence. Keep the archive separately before
cleaning `target/`. No raw data, binaries or archive is a tracked contribution.

To reproduce, check out the full R0 SHA with a clean tree, use Rust 1.91.0 and the
recorded native pins/header setup, and prepare the pinned Python reference env
using the [distance setup](../distances/README.md#prepare-the-reference-environment).
The following are the executed controller commands; use a fresh attempt directory
for each rerun. The ext4 reference path is machine-local and may be replaced by an
equivalent pinned environment whose identity is recorded separately.

```sh
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
gudhi=target/vr-h1-safe/target/native-sources/gudhi
ripser=target/vr-h1-safe/target/native-sources/ripser
boost=target/vr-h1-safe/target/native-sources/boost/usr/include
out=target/benchmarks/commit-b2c3bd5eebf0/phase2-baseline

python3 tools/benchmark_native.py --samples 12 --cpu 0 \
  --gudhi-source "$gudhi" --ripser-source "$ripser" --boost-include "$boost" \
  --output "$out/native/run-001"
python3 tools/benchmark_rips_pipeline.py --samples 12 --order-seed 2401 --cpu 0 \
  --kernel-revision b2c3bd5eebf0c1193f30ea651012b8ae61b274c1 \
  --gudhi-source "$gudhi" --ripser-source "$ripser" --boost-include "$boost" \
  --output "$out/pipeline/run-001"
python3 tools/benchmark_distances.py --samples 12 --order-seed 2401 --cpu 0 \
  --groups baseline --topp-source target/native-sources/topp-linux \
  --cgal-include target/native-sources/distance-deps/root/usr/include \
  --boost-include "$boost" \
  --gudhi-python /tmp/cocycle-phase2-oracle-20260930-b2c3bd5/bin/python \
  --output "$out/distance/run-001"
```

The archived `orchestration/phase2-baseline-linux.sh` additionally records the
integrated identity and replays the logical diagnostics; the summary script
asserts the successful suite statuses and all 12-sample groups before producing
the matrix. Reproductions must preserve failures and obtain a new evidence identity
instead of replacing this archive.
