# Adaptive F2 H1 kernels

[Benchmark reports](README.md)

Evidence class: comparative performance study on one host. Question: does the
adaptive H1 implementation improve the existing F2 continuation branch, and how
does its H0/H1 time to result compare with upstream Ripser C++?

The candidate is **not ready for a default-performance merge**. Exact diagrams
agree on every measured sample, but repeated-weight nonmetric matrices are
8.4–9.0 times slower than the baseline, and the ordered-label truncated bipartite
case is 3.2–3.4 times slower. Sphere8, uniform and equal-weight controls improve in
both repeats; several other rows change ranking between repeats. The matrix
access/cursor cost policy needs further ablations before claiming a general win.
No threshold was retuned after inspecting this run.

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

Candidate/Ripser ranking changes between repeats for the small uniform, circle
and nonmetric cases. Candidate/baseline ranking also changes for
circle, full bipartite and shuffled truncated bipartite. The clear repeated
regressions relative to the baseline must be resolved separately from variable
or nearly tied results. Gains on sphere8 and equal-weight inputs do not establish
universal superiority or justify the current gates on unmeasured sizes/hardware.

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
