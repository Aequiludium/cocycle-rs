# Adaptive F2 H1 kernels

[Benchmark reports](README.md)

Evidence class: comparative performance study on one host. Question: does the
adaptive H1 implementation improve the existing F2 continuation branch, and how
does its H0/H1 time to result compare with upstream Ripser C++?

Both retained studies show lower candidate H1 medians than the #53 baseline
on all 11 declared matrix workloads in both repeats. Earlier prose incorrectly
called nonmetric and truncated-bipartite speedups regressions; those claims are
withdrawn. The source-bound CSVs and table preserve the measured direction.

This remains a draft pending performance checks for direct graph entry points
and H2/H3 continuation, which also use the changed H1 preparation/handoff. The
matrix H1 study does not establish their absence of regression or a universal
win over Ripser. No algorithm threshold was retuned after either study.

## Measured identities

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

## Complete measured matrix

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
This does not replace graph/H2/H3 performance checks or establish general
superiority over Ripser.

Reproduce with the command above using candidate `2dcd5e3d8a30`, baseline
`323e617246d8ed2bd5096f7dcac2e5a6e17d6de0`, and a fresh output directory.
Local-only evidence: `target/benchmarks/commit-2dcd5e3d8a30/merge-adaptive/run-001/`;
`environment.json`, both full CSVs and `results.json` retain identities, raw
samples, schedules and validated diagrams. The merge audit also keeps a copy
under `target/pr-merges/2026-10-09/run-001/adaptive-performance/`. No public raw
archive is claimed.
