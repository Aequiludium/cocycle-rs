# Adaptive F2 H1 kernels

[Benchmark reports](README.md)

Evidence class: comparative performance study on one host. Question: does the
adaptive H1 implementation improve the existing F2 continuation branch, and how
does its H0/H1 time to result compare with upstream Ripser C++?

The historical matrix H1 studies below showed lower candidate medians on
all 11 declared workloads in both repeats. Their corrected direction is retained.
The current integrated graph/H2/H3 pipeline study is now complete and exposes
repeated workflow regressions, so this PR remains a draft with a no-go admission
decision for the measured integrated candidate. No threshold is retuned here.

## Current graph and high-dimensional continuation decision

Evidence class: comparative before/after workflow study on one host. This new
section measures the integrated PR #56 candidate; the older matrix H1 sections
below retain their original revisions and native-v1 timing boundary.

**No-go for adoption of the current integrated default candidate.** All outputs
pass, but every measured direct-graph H1 cell is slower than current main in
both repeats: candidate/baseline median ratios 1.269-1.504. H2 ratios span
0.973-1.718 and H3 ratios 0.957-1.639, with repeated adverse rows retained below.
The odd-field H3 control has ratios below one, while several F2 graph and handoff
controls regress. This is an observed pipeline cost, not a demonstrated kernel
cause, a universal slowdown or a new threshold decision. Keep the PR draft;
correctness and historical large-matrix gains do not waive this result.
No policy threshold was retuned against these measured workloads.

### Current measured identities and scope

| Identity | Value |
| --- | --- |
| PR | [Aequiludium/cocycle-rs #56](https://github.com/Aequiludium/cocycle-rs/pull/56) |
| Candidate and committed harness | `87821d5288c56fc1f8d2b30a690b722d8f717c35` (clean independent Linux clone) |
| Baseline | `a7ee9e911dcdeaea137b682bf383f7c2e1d090d5` (current main at sampling) |
| Candidate kernel fingerprint | `8fa0fbd1f3ee2209d60512028f705374a358ab1dee91ebec8ed2a0cd04bc2532` |
| Baseline kernel fingerprint | `fdc1ae64d0df670b4ad0868687c48a3d736985fe2efd5836077ea2bd57ab8a64` |
| Candidate plus measurement sources | `f9d63acf643c2aeca5f000fdcb437608f3b64342ea4a84c138286c10a8c97729` |
| Shared Rust pipeline worker | `454a00446355db07b23dac69f456ed18ba13140835f940b3594d4322dc743e1f` |
| References | GUDHI `cba915e3ab8e1f5b1fe26eb44b407285f7af4e78`; Ripser `01add51ff64aaf40889483260cc5c3b7d0f2a1e7` |
| Run | `target/benchmarks/commit-87821d5288c5/f2-h1-continuation/run-001/` |
| UTC start/end | 2026-10-11T03:57:13.322455+00:00 / 2026-10-11T04:00:27.701892+00:00 |
| Raw results SHA-256 | `bfc0bca22f678ddabde8dd2bdacc468241f61ca819a58733cf3c513a0f08380b` |
| Environment SHA-256 | `f957635f50c6943dbd8593db4107742e5f16e14b94deaf4543992cc657cd7c42` |
| Supplemental run contract SHA-256 | `9ec487b0ebb986c9e999a73b10025b49c4756a6555d062ad69da941a53222baa` |

Rust 1.91.0 and GCC 13.3.0 compile ordinary release libraries without experimental
H2 cfg or diagnostics. Both Rust revisions use the identical candidate pipeline
worker and compiler flags. The host is Ubuntu 24.04 / Linux WSL2 on Intel Core
Ultra 7 155H; workers are serial and pinned to CPU 0, with OMP/OpenBLAS threads
set to one. Host load/frequency are uncontrolled. The native dependency manifest,
source inventories, commands and binary hashes remain in the environment file.

Use [cocycle-rips-pipeline-v2](../pipeline/README.md), not native-v1. The timer
includes public validation/conversion, graph or dense-source preparation,
construction, computation, interval export and intervening workflow costs;
fixture I/O and final JSON transport are excluded. Graph preparation and the H1
handoff are therefore charged. These timings cannot be pooled with the older
matrix H1 tables or described as isolated reduction-kernel timings.

The committed workload generator fixes 52 cases before sampling: 45 distinct
Phase-3 dense/flag H1/H2 cases at 16/32 vertices, four H3 continuation controls,
and three eight-vertex cross-polytope cases (complete F2, cutoff F2 and odd F3).
The Phase-3 seed is 20261004. Values are common float32-exact f64 inputs;
coefficients are F2 except the explicit F3 control. All requests are diagram-only,
without representatives or approximation. Full interval multiplicity and
finite/essential/censored coverage are retained. The cross-polytope independently
requires seven H0 [0,1) bars, one essential H0 bar and H3 [1,2), or right-censored
H3 at cutoff 1.

Each of two repeats has 12 fresh measured processes per backend and one separate
discarded warmup, with position-balanced order seed 20261011 + 100*case + repeat.
All 104 case/repeat records and all 5,408 processes completed and passed output
validation (4,992 measured, 416 warmups). Each sample was compared with the
baseline and both pinned native references; there were no failures, exclusions,
reruns or discarded outliers. Reference timings remain in raw data; this section
makes no cross-library ranking claim. Builds finished before sampling.

### Complete current workflow results

B = current main, C = integrated candidate. Times are **median [min,max] ms**
from 12 processes. Ratios are C/B, independently for each repeat. RSS is the
maximum process peak KiB across both repeats for B/C; it is neither retained
workspace bytes nor allocator-live memory. Short timings and wide spreads limit
interpretation, but the repeated graph/H1 regressions do not support admission.

| Workload / max homology | B1 | C1 | B2 | C2 | C/B 1/2 | Peak RSS B/C KiB |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| uniform16_dense / H2 | 1.742 [1.396,2.066] | 2.139 [1.826,2.518] | 2.049 [1.638,2.539] | 2.413 [2.229,2.837] | 1.228/1.178 | 2640/2816 |
| uniform16_h1_dense / H1 | 1.473 [1.281,1.769] | 1.789 [1.444,1.997] | 1.715 [1.332,2.199] | 1.769 [1.313,2.336] | 1.214/1.031 | 2464/2640 |
| circle16_dense / H2 | 3.107 [2.769,4.140] | 4.104 [3.683,5.625] | 3.128 [2.724,3.743] | 3.678 [3.416,4.288] | 1.321/1.176 | 2816/2992 |
| circle16_h1_dense / H1 | 1.714 [1.380,2.273] | 1.927 [1.534,2.869] | 1.509 [1.293,2.784] | 2.006 [1.600,2.671] | 1.124/1.329 | 2464/2640 |
| clusters16_dense / H2 | 2.045 [1.736,2.285] | 2.199 [1.868,3.150] | 2.083 [1.884,2.793] | 2.232 [1.919,2.847] | 1.075/1.071 | 2640/2816 |
| clusters16_h1_dense / H1 | 1.673 [1.488,1.825] | 1.693 [1.426,2.092] | 1.735 [1.459,2.246] | 1.749 [1.554,1.979] | 1.012/1.008 | 2464/2640 |
| duplicates16_dense / H2 | 2.536 [2.054,3.168] | 3.068 [2.972,3.895] | 2.323 [2.112,2.577] | 3.049 [2.829,3.232] | 1.210/1.312 | 2640/2816 |
| duplicates16_h1_dense / H1 | 1.518 [1.343,1.989] | 1.844 [1.655,2.593] | 1.563 [1.276,3.362] | 1.689 [1.540,2.217] | 1.215/1.081 | 2464/2640 |
| sphere16_dense / H2 | 2.105 [1.956,2.340] | 2.823 [2.606,3.099] | 2.151 [2.019,2.378] | 2.840 [2.592,5.052] | 1.341/1.320 | 2640/2816 |
| sphere16_h1_dense / H1 | 1.513 [1.205,2.260] | 1.834 [1.454,1.927] | 1.542 [1.180,2.057] | 1.762 [1.572,3.018] | 1.212/1.142 | 2464/2640 |
| noisy_sphere16_dense / H2 | 2.325 [2.104,2.784] | 3.100 [2.851,4.468] | 2.067 [1.889,2.279] | 2.882 [2.577,3.269] | 1.333/1.394 | 2640/2816 |
| noisy_sphere16_h1_dense / H1 | 1.417 [1.245,1.725] | 1.829 [1.475,2.502] | 1.451 [1.242,2.084] | 1.817 [1.518,2.296] | 1.291/1.252 | 2464/2640 |
| nonmetric_ties16_dense / H2 | 2.831 [2.689,3.608] | 3.569 [3.161,4.165] | 2.938 [2.650,3.505] | 3.520 [3.302,4.124] | 1.261/1.198 | 2816/2992 |
| equal_clique16_dense / H2 | 2.777 [2.245,3.666] | 2.871 [2.674,3.430] | 2.612 [2.411,3.461] | 3.055 [2.604,3.630] | 1.034/1.170 | 2816/2992 |
| low_degree16_flag / H2 | 1.551 [1.225,1.697] | 1.924 [1.660,2.332] | 1.301 [1.192,1.421] | 1.979 [1.692,2.322] | 1.241/1.521 | 2464/2640 |
| low_degree16_h1_flag / H1 | 1.095 [0.899,1.355] | 1.525 [1.349,1.741] | 1.173 [0.969,1.781] | 1.538 [1.091,1.716] | 1.392/1.311 | 2464/2464 |
| geometric16_flag / H2 | 1.464 [1.296,1.746] | 2.341 [1.995,3.600] | 1.523 [1.301,2.080] | 2.333 [2.064,3.098] | 1.599/1.532 | 2640/2816 |
| geometric16_h1_flag / H1 | 1.107 [1.003,1.271] | 1.511 [1.420,1.797] | 1.192 [0.996,2.511] | 1.615 [1.477,2.113] | 1.365/1.355 | 2464/2640 |
| clique_heavy16_flag / H2 | 1.602 [1.281,2.555] | 2.357 [2.114,2.836] | 1.512 [1.390,1.737] | 2.305 [2.160,3.021] | 1.471/1.525 | 2640/2816 |
| clique_heavy16_h1_flag / H1 | 1.182 [0.929,1.526] | 1.672 [1.532,1.992] | 1.185 [1.045,1.796] | 1.574 [1.307,2.489] | 1.415/1.328 | 2464/2640 |
| high_fill16_flag / H2 | 1.682 [1.404,2.049] | 2.661 [2.233,3.394] | 2.019 [1.523,2.405] | 2.914 [2.659,3.542] | 1.582/1.444 | 2640/2816 |
| high_fill16_h1_flag / H1 | 1.155 [0.912,1.331] | 1.677 [1.503,1.911] | 1.196 [0.969,1.898] | 1.743 [1.258,2.619] | 1.452/1.458 | 2464/2640 |
| octahedron_dense / H2 | 1.523 [1.312,1.987] | 2.616 [2.144,3.463] | 1.457 [1.184,1.963] | 2.229 [1.923,2.933] | 1.718/1.530 | 2464/2640 |
| uniform32_dense / H2 | 6.277 [6.060,8.312] | 6.938 [6.595,8.203] | 6.579 [6.057,7.061] | 7.062 [6.513,8.040] | 1.105/1.073 | 3168/3344 |
| uniform32_h1_dense / H1 | 1.553 [1.276,1.843] | 1.877 [1.700,2.631] | 1.743 [1.573,2.074] | 1.905 [1.679,3.520] | 1.208/1.093 | 2640/2816 |
| circle32_dense / H2 | 28.458 [25.883,35.172] | 28.216 [26.826,39.470] | 28.378 [25.667,47.501] | 28.302 [26.527,30.171] | 0.992/0.997 | 4400/4576 |
| circle32_h1_dense / H1 | 1.757 [1.632,2.552] | 2.185 [1.759,3.034] | 1.775 [1.607,2.386] | 2.164 [1.836,2.877] | 1.243/1.219 | 2640/2816 |
| clusters32_dense / H2 | 7.853 [7.575,9.274] | 8.701 [8.240,9.337] | 8.337 [7.580,11.285] | 9.440 [8.223,11.417] | 1.108/1.132 | 3344/3520 |
| clusters32_h1_dense / H1 | 1.640 [1.215,2.303] | 1.933 [1.657,3.042] | 1.556 [1.326,2.318] | 1.847 [1.587,2.212] | 1.179/1.187 | 2640/2816 |
| duplicates32_dense / H2 | 14.942 [13.926,16.265] | 15.669 [14.622,18.971] | 15.085 [14.265,21.746] | 15.582 [14.094,18.633] | 1.049/1.033 | 3696/3872 |
| duplicates32_h1_dense / H1 | 1.857 [1.565,2.469] | 1.995 [1.800,2.825] | 1.863 [1.422,2.420] | 1.914 [1.522,2.439] | 1.075/1.027 | 2640/2816 |
| sphere32_dense / H2 | 16.164 [14.833,19.588] | 16.983 [15.324,21.500] | 16.980 [14.200,22.662] | 17.439 [14.750,21.006] | 1.051/1.027 | 3872/4048 |
| sphere32_h1_dense / H1 | 1.612 [1.296,2.221] | 1.818 [1.520,2.314] | 1.583 [1.332,2.408] | 1.819 [1.568,2.284] | 1.128/1.149 | 2640/2816 |
| noisy_sphere32_dense / H2 | 14.781 [14.215,17.065] | 14.980 [13.891,18.927] | 16.438 [15.219,19.542] | 16.302 [14.900,19.936] | 1.013/0.992 | 3872/4048 |
| noisy_sphere32_h1_dense / H1 | 1.893 [1.526,2.662] | 2.044 [1.714,3.031] | 1.906 [1.546,2.396] | 1.773 [1.488,2.173] | 1.080/0.930 | 2640/2816 |
| nonmetric_ties32_dense / H2 | 35.698 [33.913,40.149] | 36.722 [33.970,41.118] | 37.007 [34.155,38.276] | 36.015 [34.655,38.629] | 1.029/0.973 | 5104/5280 |
| equal_clique32_dense / H2 | 21.429 [20.362,24.787] | 21.906 [20.740,23.705] | 21.507 [20.134,25.914] | 21.360 [20.847,25.090] | 1.022/0.993 | 4400/4576 |
| low_degree32_flag / H2 | 1.405 [1.130,1.925] | 2.015 [1.862,2.330] | 1.382 [1.139,1.878] | 2.063 [1.816,2.245] | 1.434/1.493 | 2640/2816 |
| low_degree32_h1_flag / H1 | 1.218 [1.077,1.798] | 1.670 [1.509,2.446] | 1.226 [1.035,1.947] | 1.707 [1.450,1.954] | 1.371/1.393 | 2464/2640 |
| geometric32_flag / H2 | 3.527 [3.037,4.359] | 4.517 [3.941,6.437] | 3.223 [3.074,3.466] | 4.209 [3.761,4.941] | 1.281/1.306 | 2816/2992 |
| geometric32_h1_flag / H1 | 1.320 [1.059,1.835] | 1.724 [1.374,2.350] | 1.215 [1.149,1.819] | 1.827 [1.500,2.134] | 1.306/1.504 | 2640/2816 |
| clique_heavy32_flag / H2 | 1.623 [1.481,1.877] | 2.486 [2.141,3.220] | 1.748 [1.367,3.311] | 2.718 [2.282,3.508] | 1.532/1.555 | 2640/2816 |
| clique_heavy32_h1_flag / H1 | 1.215 [0.907,2.493] | 1.721 [1.395,2.983] | 1.193 [0.954,2.240] | 1.743 [1.520,3.236] | 1.416/1.461 | 2464/2640 |
| high_fill32_flag / H2 | 10.410 [9.609,11.980] | 11.429 [10.613,12.987] | 10.514 [9.722,12.571] | 12.147 [10.773,14.192] | 1.098/1.155 | 3520/3696 |
| high_fill32_h1_flag / H1 | 1.350 [1.212,1.813] | 1.714 [1.523,1.968] | 1.402 [1.169,1.991] | 1.794 [1.583,2.060] | 1.269/1.280 | 2640/2816 |
| uniform16_h3_dense / H3 | 2.170 [1.947,3.018] | 2.837 [2.507,3.193] | 2.199 [1.891,2.848] | 2.607 [2.349,4.012] | 1.308/1.186 | 2640/2816 |
| equal_clique16_h3_dense / H3 | 6.494 [6.106,7.505] | 6.860 [6.316,7.859] | 6.540 [6.198,7.190] | 6.920 [6.527,7.500] | 1.056/1.058 | 3344/3520 |
| low_degree16_h3_flag / H3 | 1.241 [1.132,2.728] | 2.034 [1.626,2.422] | 1.346 [1.173,1.607] | 2.121 [1.794,2.412] | 1.639/1.576 | 2464/2640 |
| high_fill16_h3_flag / H3 | 2.010 [1.729,3.329] | 2.847 [2.528,3.278] | 1.965 [1.653,3.470] | 2.818 [2.519,3.869] | 1.416/1.434 | 2640/2816 |
| cross_polytope_h3_dense / H3 | 1.537 [1.359,2.221] | 2.135 [1.797,2.627] | 1.710 [1.483,2.112] | 2.250 [1.936,2.927] | 1.389/1.316 | 2640/2816 |
| cross_polytope_h3_cutoff_dense / H3 | 1.633 [1.466,2.357] | 2.300 [1.840,3.629] | 1.767 [1.394,2.769] | 2.154 [1.798,2.751] | 1.409/1.219 | 2640/2816 |
| cross_polytope_h3_odd_dense / H3 | 1.464 [1.295,1.819] | 1.401 [1.329,1.787] | 1.576 [1.365,1.785] | 1.509 [1.385,2.116] | 0.957/0.957 | 2640/2640 |

### Reproduce the current continuation check

Use a clean checkout of the measured candidate, fixed clean reference sources,
Boost headers, Linux, Python 3.12 and Rust 1.91.0. Choose an allowed CPU; preserve
thread settings and finish other builds/tests before collecting samples.

```sh
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 tools/benchmark_flag_h1.py \
  --suite continuation --baseline a7ee9e911dcdeaea137b682bf383f7c2e1d090d5 \
  --gudhi-source target/native-sources/gudhi \
  --ripser-source target/native-sources/ripser --cpu 0 \
  --samples 12 --repeats 2 \
  --output target/benchmarks/commit-87821d5288c5/f2-h1-continuation/run-002
```

Supply `--boost-include` for a non-system header location. Raw evidence remains
local-only, with no public archive claimed. The report update is a separate
commit, not another measurement. Future ablations or routing changes require
separate committed source identities, independent controls and fresh acceptance;
this study supplies neither a retuned policy nor default adoption.

## Historical matrix H1 measured identities

| Identity | Value |
| --- | --- |
| PR | No associated PR at measurement; proposed branch `codex/f2-h1-adaptive-kernels`, stacked on #53 |
| Candidate and harness commit | `f98d93eb7eeed56916e52bcc877ed33370609dda` (clean branch head at measurement) |
| Baseline commit | `c3fee147ca50673d2221d961397fa6c281056da2` (#53 head when measured) |
| Ripser revision | `01add51ff64aaf40889483260cc5c3b7d0f2a1e7` |
| Candidate kernel fingerprint | `e9a40925e2b8859bd5c4b899db078d88d67938059e6ab88aae7a2f675a5274ab` |
| Baseline kernel fingerprint | `73094509600195806057f8ecb1c36822eafcf63cced02e5ad4f3763fca270757` |
| Candidate plus harness fingerprint | `88555940a2b018569be81abf7a893f135ab13c4d8b760c9077b2a04c68360215` |

The controller commits both implementation and harness before sampling. Kernel
fingerprints hash the sorted relative-path/file-SHA256 inventory of `Cargo.toml`,
`Cargo.lock` and `src`; the environment artifact includes every entry. Baseline
sources come from `git archive` of the full baseline SHA. Both Rust libraries use
the same candidate native worker adapter. The report is a later documentation
commit and is not relabeled as a measured revision.

## Contract and execution

This uses [cocycle-native-v1](../protocol.md): ordinary exact F2 H0/H1, precomputed
condensed dissimilarities, closed edge cutoffs, omitted zero-length pairs, complete
or explicitly censored diagrams, and no representatives or sparse approximation.
Distances and cutoffs are quantized once to float32 and stored losslessly as f64
for all three workers. Rust computes in f64; upstream Ripser uses its default f32.
This is not an equal-scalar-width comparison or a point-entry timing study.

Fixture I/O, validation and native input layout preparation precede each timer.
Construction, reduction, owned interval assembly and intermediate destruction are
inside it. JSON serialization and final result destruction follow it. Ripser's
unrestricted enclosing-radius computation and explicit-cutoff sparse construction
remain timed. Its six numeric output replacements are the existing pinned
adapter, with no algorithm changes.

All 11 workloads were declared in `workloads()` before measurement. Uniform uses
LCG seed 1729, circle uses equally spaced angles, grid is the 12×12 planar grid,
and sphere8 uses normalized cube samples on S7, not uniform sphere sampling.
Nonmetric has 64 discrete edge weights and makes no metric assumption. Bipartite
uses cross-part weights 1 and within-part weights 2; cutoff 1 retains K64,64.
The shuffled variant permutes labels with seed 7919. Equal has all weights 1.
Fixture hashes, retained edge counts, complete sample schedules and diagrams are
stored in the raw artifact. The analytic bipartite oracle and exact interval
multisets, including multiplicity and endpoint coverage, are checked.

Two repeats use the same compiled binaries and independently seeded sample order.
Each cell has 12 fresh serial processes, one computation each, no warmup, and
position-balanced backend rotations. The controller seed is 20261008. All 792
planned measurements completed and passed output validation: 22 case/repeat
records, 264 samples per backend, no exclusions, failures or discarded outliers.
Repeats are separate; they are not independent builds or separate hosts.

The host is Linux x86_64 on an AMD EPYC 7H12, with worker affinity CPU 2 from an
inherited set of CPUs 0–15. Rust 1.92.0 uses release/opt-level=3; GCC 11.4.0 uses
C++17, `-O3 -DNDEBUG`. Workers set OMP/OpenBLAS threads to one. Each process has a
30-second wall-time limit and 2048 MiB address-space cap. Frequency and host load
are uncontrolled. Short timings and changing rankings limit interpretation;
additional repetitions alone cannot correct those environmental effects.

## Historical matrix H1 complete measured matrix

Each timing is **median [minimum, maximum] milliseconds**, with 12 samples.
B = baseline Rust, C = candidate Rust, R = upstream Ripser C++. Lower is faster.
Every workload appears in both repeats; no aggregate ranking is reported.

| Workload | B repeat 1 | C repeat 1 | R repeat 1 | B repeat 2 | C repeat 2 | R repeat 2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| uniform_h1_32 | 0.336 [0.328, 0.402] | 0.300 [0.296, 0.345] | 0.301 [0.296, 0.312] | 0.336 [0.331, 0.343] | 0.296 [0.292, 0.376] | 0.302 [0.295, 0.340] |
| uniform_h1_64 | 1.468 [0.670, 1.499] | 1.257 [0.573, 1.386] | 1.367 [0.610, 1.451] | 1.119 [0.666, 1.513] | 1.105 [0.568, 1.322] | 1.024 [0.614, 1.458] |
| uniform_h1_256 | 13.222 [13.115, 13.672] | 10.950 [10.910, 15.298] | 11.680 [11.595, 11.886] | 13.234 [13.148, 20.926] | 10.936 [10.905, 15.473] | 11.650 [11.601, 11.841] |
| circle_h1_128 | 15.986 [15.896, 17.077] | 8.799 [8.759, 8.975] | 13.487 [13.380, 19.164] | 16.002 [15.877, 16.060] | 8.833 [8.784, 15.208] | 13.560 [13.462, 13.695] |
| sphere8_h1_128 | 7.711 [7.626, 7.919] | 3.903 [3.884, 4.090] | 7.868 [7.809, 11.821] | 7.719 [7.657, 8.083] | 3.910 [3.878, 3.942] | 7.864 [7.804, 7.934] |
| nonmetric_h1_128 | 229.117 [227.253, 234.806] | 27.155 [26.985, 31.125] | 243.514 [241.342, 250.939] | 227.635 [227.078, 236.152] | 27.108 [26.980, 27.404] | 243.372 [242.130, 273.863] |
| grid_h1_144 | 7.018 [6.979, 11.691] | 4.826 [4.805, 4.877] | 4.838 [4.798, 8.098] | 7.019 [6.989, 10.128] | 4.826 [4.789, 5.828] | 4.837 [4.802, 5.274] |
| bipartite_h1_128 | 21.947 [21.651, 30.541] | 19.841 [15.156, 23.222] | 22.976 [21.983, 29.334] | 21.863 [21.496, 29.722] | 15.233 [15.158, 22.748] | 22.226 [21.915, 30.661] |
| bipartite_h1_128_cutoff | 6.790 [6.021, 11.857] | 0.883 [0.846, 1.967] | 2.025 [1.553, 2.928] | 8.987 [5.966, 12.345] | 0.858 [0.845, 1.877] | 2.780 [1.547, 3.772] |
| bipartite_shuffled_h1_128_cutoff | 9.475 [6.099, 12.960] | 1.220 [0.908, 2.022] | 3.374 [2.037, 4.648] | 9.659 [6.102, 12.701] | 2.012 [0.918, 2.134] | 2.073 [2.037, 4.702] |
| equal_h1_128 | 1.265 [1.251, 2.928] | 1.182 [1.166, 1.212] | 0.840 [0.829, 0.894] | 1.260 [1.244, 2.713] | 1.181 [1.165, 1.959] | 0.848 [0.821, 1.891] |

Candidate medians are lower than baseline medians in every original cell.
Candidate/Ripser ranking changes between repeats for uniform_h1_64; near ties
and variable spreads limit individual cross-library comparisons. Repeated gains
on these matrix H1 inputs do not establish universal superiority or validate
unmeasured entry points, dimensions, sizes or hardware.

Process RSS before timing, high-water baseline, peak RSS and high-water growth are
retained separately in raw samples and CSVs. They include input, allocator and
runtime memory and are not live algorithm allocations. No memory win is claimed.
The process address-space cap does not establish a library memory budget.

## Reproduction and evidence

Run on a clean checkout of the measured candidate with Rust and the pinned Ripser
source already installed; choose an allowed CPU. Source fetching is a separate
explicit setup step.

```sh
python3 tools/benchmark_flag_h1.py \
  --baseline c3fee147ca50673d2221d961397fa6c281056da2 \
  --ripser-source /path/to/pinned/ripser --cpu 2 --samples 12 --repeats 2 \
  --output target/benchmarks/commit-f98d93eb7eee/f2-h1-adaptive/run-002
```

Retained run: `target/benchmarks/commit-f98d93eb7eee/f2-h1-adaptive/run-001/`,
UTC 2026-10-08 02:22:24–02:23:03. `environment.json` binds source inventories,
compiler/build commands, binary hashes, CPU/affinity and upstream pin.
`results.json` retains all validated diagrams, samples and schedules;
`summary-repeat-1.csv` and `summary-repeat-2.csv` retain the complete matrix.
Fixtures, binaries, baseline source archive contents and build logs remain outside
Git. Evidence is local-only; no public archive/upload is claimed. Reproduction
produces a new attempt and must not overwrite this one.

## Independent build and repeat before merge

The 2026-10-09 audit verified the original retained binary hashes and rebuilt
all three workers before a second complete study. Measured candidate/harness:
`2dcd5e3d8a306f3cb4c4f23e486cdfc5add78db1`; baseline integrated #53 head:
`323e617246d8ed2bd5096f7dcac2e5a6e17d6de0`. Both kernel fingerprints and the
candidate-plus-harness fingerprint match the original study. The current report
is a later documentation change, not a newly measured algorithm revision.

The same 11 fixtures, native-v1 timing boundary, CPU 2, Rust 1.92.0/GCC 11.4.0,
12 fresh processes per cell, two repeats and zero warmups produced 792 validated
samples with no failures, exclusions or discarded outliers. All builds completed
before measurement; no other local build/test ran during timing. Host frequency
and load remain uncontrolled. The original and new protocols are identical;
the repeats below remain separate from the original results.

B/C/R again mean baseline/candidate/Ripser; values are median [min, max] ms.

| Workload | B repeat 1 | C repeat 1 | R repeat 1 | B repeat 2 | C repeat 2 | R repeat 2 |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| uniform_h1_32 | 0.163 [0.160, 0.179] | 0.141 [0.137, 0.172] | 0.134 [0.129, 0.160] | 0.164 [0.157, 0.196] | 0.140 [0.135, 0.171] | 0.135 [0.133, 0.157] |
| uniform_h1_64 | 0.715 [0.703, 1.260] | 0.589 [0.569, 1.030] | 0.629 [0.602, 1.396] | 0.713 [0.704, 0.811] | 0.584 [0.570, 1.038] | 0.618 [0.607, 1.075] |
| uniform_h1_256 | 14.215 [14.132, 14.662] | 10.967 [10.911, 11.103] | 11.720 [11.565, 11.930] | 14.181 [14.115, 14.754] | 10.950 [10.883, 11.013] | 11.802 [11.549, 12.092] |
| circle_h1_128 | 17.719 [17.610, 18.181] | 8.820 [8.785, 9.034] | 13.557 [13.465, 13.682] | 17.706 [17.587, 18.200] | 8.807 [8.770, 9.564] | 13.576 [13.506, 14.387] |
| sphere8_h1_128 | 8.305 [8.271, 8.506] | 3.927 [3.901, 4.225] | 7.907 [7.868, 8.220] | 8.324 [8.279, 8.604] | 3.938 [3.914, 4.106] | 7.904 [7.825, 8.246] |
| nonmetric_h1_128 | 245.756 [243.663, 251.079] | 27.283 [27.031, 29.852] | 244.142 [241.784, 251.926] | 248.932 [244.416, 260.024] | 27.516 [27.079, 31.323] | 251.885 [243.575, 253.430] |
| grid_h1_144 | 8.140 [8.013, 8.922] | 4.931 [4.854, 5.132] | 4.944 [4.849, 5.087] | 8.411 [7.976, 9.041] | 5.164 [4.888, 5.429] | 5.118 [4.852, 5.400] |
| bipartite_h1_128 | 24.483 [24.171, 28.136] | 15.620 [15.443, 17.916] | 22.543 [22.255, 24.943] | 24.973 [23.929, 26.418] | 15.887 [15.436, 19.380] | 22.498 [22.282, 23.474] |
| bipartite_h1_128_cutoff | 7.957 [7.488, 9.500] | 1.011 [0.919, 1.284] | 1.751 [1.647, 1.904] | 7.762 [7.502, 8.537] | 0.963 [0.930, 1.132] | 1.653 [1.615, 1.796] |
| bipartite_shuffled_h1_128_cutoff | 7.447 [7.180, 8.508] | 1.071 [1.013, 1.131] | 2.200 [2.109, 2.521] | 7.441 [7.168, 11.629] | 1.035 [1.010, 2.181] | 2.181 [2.111, 4.737] |
| equal_h1_128 | 1.334 [1.285, 1.381] | 1.244 [1.207, 1.336] | 0.893 [0.859, 0.930] | 1.369 [1.309, 1.572] | 1.284 [1.200, 1.379] | 0.918 [0.861, 0.990] |

All 22 new candidate/baseline median ratios are below one, ranging from 0.111
to 0.938. The repeated nonmetric ratio is 0.111/0.111 and ordered-label truncated
bipartite is 0.127/0.124. Relative to Ripser, equal-weight H1 remains about 39–40%
slower, uniform32 about 4–6% slower, and grid is near a tie; these are not
regressions relative to #53. Other matrix cases beat Ripser in both new repeats.
This historical native-v1 study does not establish graph/H2/H3 workflow
performance or general superiority over Ripser. The separate current pipeline
study above now covers those missing controls and retains its adverse outcome.

Reproduce with the command above using candidate `2dcd5e3d8a30`, baseline
`323e617246d8ed2bd5096f7dcac2e5a6e17d6de0`, and a fresh output directory.
Local-only evidence: `target/benchmarks/commit-2dcd5e3d8a30/merge-adaptive/run-001/`;
`environment.json`, both full CSVs and `results.json` retain identities, raw
samples, schedules and validated diagrams. The merge audit also keeps a copy
under `target/pr-merges/2026-10-09/run-001/adaptive-performance/`. No public raw
archive is claimed.
