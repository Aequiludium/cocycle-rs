# Integrated spectral workstreams: resources and correctness

[Benchmarks](../README.md) / [Reporting rules](../reporting.md)

## Question and scope

Does combining the nine source PRs in [PR #160](https://github.com/Aequiludium/cocycle-rs/pull/160)
preserve their results, and what resources do the integrated verification
workloads use compared with their individual original implementations?
All twelve independent/frozen replay suites pass. The resource snapshot covers
15 cells, 300 measured processes and 30 discarded warmups, all output-validated.
It measures bounded worker workflows, including parsing and output; it is not
an isolated-kernel, cross-library or general scalability comparison.

## Revisions and protocol

| Item | Recorded value |
| --- | --- |
| Measured candidate and harness | `9588668972873f4cf71860d9ccac7b0738c2cdd0`; committed integration-branch head, not GitHub's synthetic merge commit |
| Mathematical replay | `eacf0c752914619e9f3eaaab830473aa854b35dd`; identical `src/` and `examples/` in the measured candidate |
| Controller SHA-256 | `e69361f770993cf5a15ee24c99beac4b715338a1da268c0a23bb34db47f4cd06` |
| Protocol | `spectral-integration-resource-v1`, implemented by [the controller](../../tools/benchmark_spectral_integration.py) |
| Machine | Intel(R) Core(TM) Ultra 7 155H; Windows-11-10.0.26100-SP0 |
| Compiler | `rustc 1.98.1 (48a229cea 2026-09-01)`; same edition 2024 / `-C opt-level=3` library and worker builds on both sides |
| Python | 3.12.0rc3 (tags/v3.12.0rc3:b973ab3, Sep 19 2023, 09:39:54) [MSC v.1935 64 bit (AMD64)] |
| Dirty state | Clean when the run started; worker/library sources exported from immutable Git commits |
| Controls | Serial workers after all builds; one warmup and ten measured processes per side; seed 160, balanced pair order; thread environment set to one |
| Uncontrolled factors | Host load, affinity, frequency and allocator/OS state; no address-space cap; 120-second process timeout |
| Run UTC | 2026-11-24T11:24:09.382305+00:00 through 2026-11-24T11:26:12.710216+00:00; filesystem execution metadata |

| Group | Original PR baseline SHA |
| --- | --- |
| A: sparse action / LSMR (#156) | `0798360c45e63ceb22b55c3a67f4a6787043ff9e` |
| B: Grams / harmonic restriction (#159, with #151 reference) | `7527c57baad0e9fc8fa420a4808263d6b0a33b26` |
| C: relevant pairs (#154) | `ef2d763b71a3fc0971009342c32a9b670f042304` |
| D: non-branching (#153) | `430d1c8a9dbd3508c585f1c99fb3ee6b6a2f8aad` |
| E: subspace sequences (#152) | `c25e209931dd4003819ccce6a0ca650a2d172868` |
| F: sampling / moments (#157) | `4b28c5e506a85193501120eafe91c338f4ca143f` |
| G: HeCS (#158, including #155) | `42db213a3b0fc076d5048cc44a24428aec66e764` |
| Existing main Schur control | `a3c3e1132800624743ad0f159c9303090e0b3b3f` |

## Workloads and validation

A uses late-center stars with 257/1024 retained leaves and 1,000 analytically
checked actions (`I-11^T/n`), plus the frozen 801-pair/12-control LSMR input.
B uses the tetrahedron-shell Gram and the 35-dimensional/32-dimensional-kernel
positive-restriction solve. C uses real degree-zero prefix paths of size 16/32,
five grid points and every legal pair; the resulting zero/unit path Laplacians
were checked analytically before measurement. D uses the complete supported
1,026-pair non-branching fixture stream. E runs complete smooth/abrupt fixed-degree
sequences on their frozen 97 KB inputs. F replays validated batches 0000/0303,
preserving the exact probes/draws and full Gram/moment outputs. G uses frozen
20-vertex enriched-5 inputs, with direct/core selection. The main control uses
the ordinary full 1,114-pair component fixture stream.

All real calculations retain f64 identity-coordinate contracts. Fixed inputs,
worker arguments, coefficient/degree/cutoff metadata where applicable and complete
source/binary/input hashes are in `specification.json`, `environment.json` and
the member manifest. Full labels and integer identities must agree; numerical
outputs use the unchanged absolute/relative 1e-10 gate. Only declared elapsed
timing fields are omitted. Independent correctness is established separately:

| Scope | Fresh integrated evidence |
| --- | --- |
| Rust library | 353 passing tests/doctests and five explicitly ignored tests in each debug, release and Rust 1.91 run; full all-target Clippy/MSRV checks pass |
| Python tools | Linux: 230 tests, 92 optional skips; numerical reproduction tests on Windows: 108 tests, no skips/failures |
| Persistent components / non-branching | 1,114 generic pairs; 1,026 eligible specialized pairs, with 88 explicit branching exclusions; exact rational projection, action, multiplicity and spectral-union checks |
| A | 801 construction pairs, 12 solver controls; fixed actions and Schur matrices agree with exact rational oracles. LSMR replays 208 updates and five explicit failures |
| B | 15 independent reproduction groups and 29 frozen Rust parity cases; independent rational/analytic, dense/SciPy and original-operator residual controls |
| C | 7,120 C outputs plus 801 shared-A pairs; frozen R/V, events, selections and failure contracts preserved |
| D source | 90 frozen kernel intermediates, 180 standalone comparisons and three unchanged author cubical replays |
| E | 328 QR/Ritz steps, 48 full sequence solves, 36 mapping/projection controls and 15 failure/edge cases |
| F | 1,440 moments, 7,036 frozen sampling trials and 1,280 native-RNG cases; all original rank/sandwich failures retained |
| G | 34 frozen first-round parity cases and 34 independent native controls, including two expected rank-loss cases |

Local numerical references use GUDHI 3.11.0, NumPy 2.4.6, SciPy 1.17.0 and
SymPy 1.14.0; B's frozen parity also uses its original Python 3.10/NumPy 2.2.6
environment. GUDHI corroborates topology only; exact rational/analytic and
independent dense controls provide real-matrix/spectral expectations. The PR's
separate native GUDHI/Ripser/Topp CI regressions do not supply real-spectrum truth.

## Worker time

Milliseconds, median [minimum, maximum], ten fresh processes per side.
The wrapper starts at worker main and ends after flushed output: fixture I/O,
parsing, computation, checks and output are included; process startup is excluded.
The ratio is the median of ten paired **candidate / original** ratios; above one
means the integrated workflow took longer. No rows or outliers are discarded.

| Workload | Original ms | Integrated ms | Paired ratio |
| --- | --- | --- | --- |
| A-star-257-1000 | 2.363 [2.232, 2.754] | 2.433 [2.191, 2.717] | 1.014 |
| A-star-1024-1000 | 9.689 [9.001, 13.949] | 9.743 [8.988, 11.290] | 0.987 |
| A-lsmr-frozen | 14.865 [13.461, 17.950] | 14.892 [13.203, 18.189] | 1.000 |
| B-gram-shell | 0.158 [0.117, 0.296] | 0.135 [0.124, 0.371] | 1.001 |
| B-large-nullity-solve | 0.613 [0.562, 1.064] | 0.631 [0.547, 1.752] | 0.970 |
| C-prefix-path-16 | 0.797 [0.752, 1.231] | 0.870 [0.762, 1.149] | 0.997 |
| C-prefix-path-32 | 2.844 [2.601, 7.863] | 2.941 [2.706, 4.252] | 1.085 |
| D-nonbranching-suite | 68.480 [61.987, 74.276] | 69.454 [65.628, 76.894] | 1.021 |
| E-sequence-smooth | 7.758 [6.963, 9.963] | 8.267 [7.450, 9.686] | 1.034 |
| E-sequence-abrupt | 9.209 [8.298, 10.273] | 9.535 [8.231, 10.573] | 1.045 |
| F-batch-0000 | 8.124 [7.729, 9.532] | 8.352 [7.701, 9.671] | 1.014 |
| F-batch-0303 | 27.264 [23.879, 32.393] | 27.126 [25.142, 32.381] | 1.004 |
| G-hecs-enriched20-direct | 0.113 [0.107, 0.260] | 0.132 [0.114, 0.651] | 1.102 |
| G-hecs-enriched20-core | 0.119 [0.110, 0.192] | 0.130 [0.116, 0.166] | 1.084 |
| main-schur-components | 73.236 [67.319, 89.096] | 71.956 [67.746, 82.107] | 1.004 |

## Process time and memory

Process wall time runs from process creation through complete input/output
transport. Peak working set is the absolute per-process peak after completion;
it includes runtime, input, allocator retention and output, rather than live
algorithm allocation. MiB = 2^20 bytes. Values are median [minimum, maximum].

| Workload | Original wall ms | Integrated wall ms | Original peak MiB | Integrated peak MiB |
| --- | --- | --- | --- | --- |
| A-star-257-1000 | 9.043 [7.964, 10.530] | 9.359 [7.532, 10.492] | 4.580 [4.562, 4.590] | 4.590 [4.578, 4.605] |
| A-star-1024-1000 | 16.693 [14.815, 22.526] | 16.204 [14.645, 19.095] | 4.855 [4.848, 4.941] | 4.875 [4.871, 4.879] |
| A-lsmr-frozen | 22.184 [19.081, 25.696] | 22.708 [19.410, 26.728] | 4.340 [4.312, 4.348] | 4.336 [4.328, 4.344] |
| B-gram-shell | 6.373 [5.283, 11.527] | 6.572 [5.698, 10.284] | 4.336 [4.316, 4.344] | 4.340 [4.332, 4.340] |
| B-large-nullity-solve | 6.982 [6.045, 9.490] | 7.114 [5.765, 9.252] | 4.359 [4.352, 4.445] | 4.375 [4.352, 4.449] |
| C-prefix-path-16 | 8.056 [6.540, 11.328] | 7.478 [6.274, 10.468] | 4.600 [4.566, 4.609] | 4.613 [4.605, 4.621] |
| C-prefix-path-32 | 10.320 [8.253, 13.783] | 10.094 [8.190, 11.069] | 4.604 [4.590, 4.625] | 4.629 [4.609, 4.648] |
| D-nonbranching-suite | 75.639 [69.437, 83.480] | 76.616 [72.605, 84.444] | 4.672 [4.664, 4.727] | 4.713 [4.703, 4.766] |
| E-sequence-smooth | 14.536 [13.186, 18.970] | 15.651 [12.947, 17.053] | 6.121 [6.074, 6.379] | 6.121 [6.117, 6.445] |
| E-sequence-abrupt | 16.472 [14.878, 18.793] | 17.012 [14.796, 18.854] | 6.568 [6.277, 6.574] | 6.547 [6.230, 6.570] |
| F-batch-0000 | 14.749 [13.992, 16.372] | 15.327 [14.567, 16.898] | 4.602 [4.594, 4.648] | 4.602 [4.590, 4.645] |
| F-batch-0303 | 34.218 [30.345, 41.053] | 34.163 [32.127, 40.043] | 4.777 [4.633, 4.797] | 4.779 [4.602, 4.785] |
| G-hecs-enriched20-direct | 6.677 [5.437, 9.307] | 6.609 [5.525, 11.863] | 4.523 [4.516, 4.527] | 4.520 [4.504, 4.535] |
| G-hecs-enriched20-core | 6.349 [5.633, 8.321] | 6.495 [5.765, 9.072] | 4.520 [4.516, 4.523] | 4.520 [4.516, 4.527] |
| main-schur-components | 81.110 [74.736, 97.891] | 79.299 [75.691, 90.124] | 4.689 [4.664, 4.699] | 4.732 [4.691, 4.750] |

## Interpretation and retained limitations

This integration does not adopt a new default algorithm. Every measured result
agrees with its original implementation under the declared contract. The table
retains slower rows as well as faster ones. Small Gram/path cases are sensitive
to timer, I/O and OS variation, and process time may be dominated by startup.
These measurements cannot establish universal equivalence, an optimization
speedup, a total memory bound or complete-spectrum scalability.

Worker paired ratios range from 0.970 to 1.102. The existing main Schur control
is 1.004; prefix-path-32 is 1.085, and HeCS direct/core are 1.102/1.084. The latter
workflows take roughly 0.1 ms before process startup, so their relative changes
must be read beside their absolute ranges. Subspace smooth/abrupt are 1.034/1.045.
There is no consistent speedup across the declared matrix. Process peak working
sets overlap closely; this is not an allocation-equivalence claim.

Known mathematical exclusions remain: LSMR has 25 retained terminal continuation
variable discrepancies under the existing scoped breakdown comparison; basic
clustered LOBPCG can fail; non-branching inputs exclude branching pairs; sampling
retains 630 rank-loss and 981 epsilon-0.5 sandwich failures; full KID reconstruction
and a probability theorem remain unresolved. HeCS rejects selected rank loss.
These are expected/restricted outcomes, not newly successful cases.

The invalidated pre-measurement warmup at `b3af646308884fa8c2113a2b064480eef08cf13d`
mistakenly compared `TOTAL_SECONDS` as mathematical output. Its raw files remain
under that revision. Only elapsed-field parsing was corrected; mathematical
tolerances were not changed. It contributes no samples to the tables above.

## Evidence and replay

Raw performance evidence is **local-only**, under
`target/benchmarks/commit-958866897287/spectral-integration/run-001/`.
Its member-manifest SHA-256 is `5f6c6c2418c2654e127d7518a49e6086d4a6751b777102f493ff0eac745fc45b`.
The independent replays are under
`target/benchmarks/commit-eacf0c752914/integration-correctness/run-001/`.
No durable external performance archive is configured. These raw local samples
are not claimed to be publicly downloadable. The report commit is distinct from
both measured revisions.

The run retains the exact `specification.json`, all input bytes, per-sample
stdout/stderr, extracted committed sources, instrumented wrappers, compiler logs,
binary fingerprints, environment, all samples and summaries. With that frozen
specification/input bundle and Git objects for the listed SHAs, replay using:

```sh
python tools/benchmark_spectral_integration.py --specification target/integration-plan.json --output target/benchmarks/commit-CANDIDATE/spectral-integration/run-002
```

Use a clean checkout of the measured candidate and a new output directory.
Relevant public correctness/regeneration commands are maintained in the
[tools guide](../../tools/README.md#integrated-spectral-resource-snapshots) and
the included per-track reproduction guides. Hosted CI status belongs in PR #160;
CI smoke timings, the resulting main merge and a published release are distinct
from this local resource snapshot.
