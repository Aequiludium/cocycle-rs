# Integrated spectral routes versus unaccelerated algorithms

[Benchmarks](../README.md) / [Reporting rules](../reporting.md)

## Composed partial-spectrum increment

Measured clean source: `aaf4db38aaaabb51d8ca985307d12e8476510410`.
This increment connects fixed persistent actions to partial LOBPCG, adds bounded
boundary-image positive modes, single-sided Grams, compact non-branching
projection, reverse projection updates and exact-event upper-factor reuse.
All routes are explicit; existing dense constructors and complete spectra remain
available. The historical 30-cell results below retain their original revision.

Independent correctness checks cover 6,264 route/pair cases with 108 declared
rejections, Fraction real ranks, exact-kernel/SVD matrices, NumPy eigenvalues,
original-coordinate vector residuals and GUDHI 3.11.0 F11 topology corroboration.
Maximum normalized eigenvalue error is 8.172e-15; residual is 3.095e-15.
GUDHI is not a real-spectrum oracle. Debug/release/Rust 1.91 each pass 360 tests
and doctests with five explicit ignores. Strict rustdoc, seven guide doctests,
Clippy and source/docs checks pass. One numerical regression prevents
normalizing an annihilated upper column's roundoff into a range direction.

The new performance matrix compares **unaccelerated dense algorithms against
composed routes in the same committed build**. It has 13 cells, ten measured
fresh processes per route and one discarded warmup: 260 measurements and 26
warmups. All outputs pass independently frozen expectations under the declared
1e-8 normalized error contract. No previous-PR timings enter these comparisons.
Baseline and candidate have identical inputs and requested exported values;
native preparation, numerical verification and export are included. Some
iterative routes calculate tighter internal residuals than the external contract.

| Fixture / requested output | Dense worker median (ms) | Composed median (ms) | Paired worker ratio | Paired process ratio |
| --- | ---: | ---: | ---: | ---: |
| 64 vertices, rank-3 positive spectrum | 0.382 | 0.182 | 2.20x | 1.08x |
| 256 vertices, rank-3 positive spectrum | 15.368 | 0.301 | 49.95x | 3.31x |
| 64-target path, compact projection + positive mode | 4.315 | 0.236 | 18.61x | 1.68x |
| 256-target path, compact projection + positive mode | 473.423 | 0.564 | 830.78x | 64.61x |
| 64-leaf star, three largest modes | 0.797 | 0.256 | 3.09x | 1.05x |
| 128-leaf star, three largest modes | 13.461 | 0.390 | 34.31x | 3.10x |
| 64-leaf star, three smallest positive modes | 0.828 | 0.652 | 1.18x | 1.06x |
| 128-leaf star, three smallest positive modes | 13.333 | 2.670 | 4.99x | 2.14x |
| 3-by-64 Gram, positive modes/projector | 3.104 | 0.558 | 5.50x | 1.41x |
| 3-by-128 Gram, positive modes/projector | 49.974 | 2.803 | 17.16x | 5.72x |
| H1, 31 fixed-endpoint largest modes, cold vectors | 4.386 | 8.628 | 0.53x | 0.72x |
| Same fixed-endpoint requests, warm vectors | 4.529 | 8.150 | 0.55x | 0.76x |
| H1, 496 requested largest modes, relevant factors | 26.293 | 37.186 | 0.71x | 0.75x |

Ratios are medians of paired baseline/candidate times, not ratios of these
rounded medians; values below one mean slower. Worker timing starts at main
entry and includes I/O and checks; process timing also includes Windows startup
and transport. Small worker gains therefore shrink at process level.

The rank-3 graph has three disjoint edges and isolated vertices; its three
positive values are 2. Path controls retain only two vertices while the target
has 64/256 vertices. Their large gains replace dense target Schur elimination
with forced signed components, rather than demonstrating a generic 256-variable
positive solver. Star operators are I-11^T/n on retained leaves. High-rank
positive setup still retains n-by-(n-1) range storage. Gram controls use the
previous analytical Hadamard rows and return the same three values plus their
complete positive-space projector; projector export is included on both sides.
Batch controls reuse the original 53-simplex, 32-scale H1 fixture. Empty q-bases
are excluded from one-mode requests. The direct baseline computes each full
matrix and spectrum, then selects its largest value. Batch routes use eight
lookahead columns and preserve the fixed declared solver policy.

The small batch regressions remain visible: factor/vector reuse is functional,
but does not overcome setup/iteration costs here. No faster default or automatic
adoption follows. The 256-vertex rank-3 process peak decreases from 6.09 to
4.63 MiB; the compact path decreases from 6.13 to 4.65 MiB. Relevant partial
batching increases peak from 4.68 to 5.82 MiB. These are absolute process peaks,
not allocation counts, and small differences include runtime overhead.

Environment: Intel Core Ultra 7 155H, Windows 11 build 26100, rustc 1.98.1,
edition 2024 and `-C opt-level=3`, serial measurements after all builds, thread
environment one, balanced pair order and seed 160. Affinity, frequency and host
load are uncontrolled. Recorded interval: 2026-11-24T22:59:32.105894+08:00 to
2026-11-24T22:59:42.227793+08:00. This remains a local bounded resource snapshot.

Raw evidence is local-only under
`target/benchmarks/commit-aaf4db38aaaa/composed-unaccelerated/run-001/`.
Manifest SHA-256:
`abba60a6283e08ca10f94c5de2e80a67b89d5b99c83ec6787806f4901c20f622`.
Controller SHA-256:
`93b8d8da2e5e101564b69993202b02006283ddb5717dc63b480133c0ef95fa1c`.
The bundle retains every sample/range, process peak, input/oracle, source and
binary hash and compiler invocation. Frozen preflight inputs are in
`target/composed-performance-inputs-001/`; clean committed correctness evidence
is in `target/composed-partial-oracle-003/`. The report commit is later than the
measured source. Its standalone dead-code annotation fix does not change the
measured library or performance adapter. Hosted status belongs to the exact PR
head, separately from these local measurements.

General sparse QR, O(n*r) range fill, exact public Betti certification, online
endpoint insertion, production Chebyshev filters and persistent HeCS/sampling
remain outside this increment. Adaptive LSMR was not substituted for a fixed
self-adjoint projection.

## Corrected question and conclusions

This [PR #160](https://github.com/Aequiludium/cocycle-rs/pull/160) experiment
compares **unaccelerated algorithms with alternative routes inside the same
integrated build**. The previous original-PR-head comparison answered integration
regression, not acceleration, and is withdrawn as acceleration evidence. Its raw
records are preserved; none of its timings enter the tables below.

All 30 declared cells completed: 600 measured fresh processes and 60 discarded
warmups, each checked against an independently frozen expectation. Twenty-six
cells request equivalent mathematical output. Four sampling cells explicitly
compare an approximate Gram with an exact dense baseline, retaining quality
failures separately. This is a local-only bounded resource snapshot.

The fixed-projection action route achieves 1.35-14.54x on the declared star
workloads. Relevant-pair q1 matrix construction achieves 2.47x, reducing to 1.75x
when complete spectra are also requested. Enriched HeCS controls achieve
1.46-2.11x including preparation. These gains do not generalize to every row:
q0 batching, single-pair batching, the non-branching suite and several Gram,
restriction and reuse controls show losses or near-one ratios. No aggregate speedup or faster
default is inferred.

## Measured identity and environment

| Item | Recorded value |
| --- | --- |
| Measured library, workers and controller | `8994b62c8e2de902ce36974b0d63c8fa7ace28c4`, clean committed integration-branch source |
| Baseline revision | The same commit; explicit unaccelerated route selection below, never an original PR head |
| Controller SHA-256 | `93b8d8da2e5e101564b69993202b02006283ddb5717dc63b480133c0ef95fa1c` |
| Protocol | `spectral-unaccelerated-comparison-v2`; [controller](../../tools/benchmark_spectral_integration.py) |
| Compiler | `rustc 1.98.1 (48a229cea 2026-09-01)`; edition 2024 and `-C opt-level=3` for every worker and the shared library |
| Machine | Intel(R) Core(TM) Ultra 7 155H; Windows-11-10.0.26100-SP0 |
| Python | 3.12.0rc3 (tags/v3.12.0rc3:b973ab3, Sep 19 2023, 09:39:54) [MSC v.1935 64 bit (AMD64)] |
| Measured UTC | 2026-11-24T12:18:58.589565+00:00 through 2026-11-24T12:19:16.964868+00:00 |
| Repetition/order | One discarded warmup and ten fresh processes per route/cell; seed 160, shuffled cells, balanced pair order |
| Controls | All builds finish before serial measurement; thread environment one; no concurrent local tests/builds |
| Limits | 120-second process timeout; no address-space cap; frequency, affinity and host load uncontrolled |

Sources, generated timing wrappers, compiler invocations, binaries and inputs are
hashed in the raw member manifest. The report commit is distinct from the
measured source. At that historical revision, library `src/` was unchanged from
the independently replayed `eacf0c752914619e9f3eaaab830473aa854b35dd` source;
the composition increment above has its own later implementation and evidence.

## Baselines and requested work

| Track | Unaccelerated baseline | Current alternative and matched output |
| --- | --- | --- |
| A | Full dense Schur construction followed by dense actions | Fixed-projection operator preparation and the same repeated actions; consume/check the final vector |
| B Gram | Form the primal Gram and solve its complete dense eigenproblem | Current `gram_spectrum`, including both complete Grams/eigensolves; return positive eigenvalues and the original-coordinate projector |
| B harmonic | Complete dense eigenproblem, then select the positive space | Current `HarmonicRestriction::compute`, including full-matrix checks, kernel normalization/complement and reduced solve; same positive values/projector |
| C | Direct Schur separately for every requested pair | One exact relevant-pair schedule and all requested full matrices, optionally also all complete spectra |
| D | Generic direct Schur on every eligible pair | Explicit non-branching construction; same full matrices, optionally complete spectra |
| E | Cold ChFSI from the original starting block for every matrix | Reuse the preceding block in original order, or pay for descriptor ordering and restore original order; same eigenvalues/projectors |
| F | Exact dense `B B^T` construction | Current exact-score/moment/sampling prototype; approximate Gram with separately checked seed/draws and quality, not equivalent exact output |
| G | Zero-start, identity-preconditioned CGLS on the original weighted upper matrix | Public HeCS solve, including topology, selection, exact coverage, factors, CGLS and original residual checks; same minimum-norm solution |

A uses 64/256/512 retained leaves, one late center and 1/100 actions. Both sides
use an explicit 4,000,000-entry per-buffer dense cap; these are not claims that
the default dense cap supports the larger cases. The independent action is
`I-11^T/n`. Source construction and preparation are included.

B uses 32/64-coordinate Hadamard controls with analytic positive spectrum
`[1,4,9]`, rank three and the corresponding analytic projector. Harmonic kernels
are supplied to both routes; discovering a topological kernel is not measured.
These APIs remain deliberately unoptimized references. In particular, timing
only their final small solve would hide the full eigensolve and normalization
currently performed during preparation.

C uses the frozen alpha12 seed-0 complex, its complete 32-point grid (528 legal
pairs), degrees 0/1, and a separate q1 last-pair request to expose setup losses.
Expected full matrices are original rational Round-1 outputs; independent dense
spectra and GUDHI field-11 torsion-free topology controls are rechecked before
measurement. Matrix-only and complete-spectrum workflows are separate rows.

D uses all 1,026 supported pairs from the non-branching correctness fixture
stream, preserving degrees, scales and original coordinates. The 88 branching
exclusions from that suite remain exclusions; no implicit fallback is introduced.
Direct matrices are independently checked by rational projection where available
and SVD projection otherwise.

E uses the original smooth/abrupt sequences: eight 24-coordinate matrices,
block size six, three requested modes, fixed degree eight and tolerance 1e-11.
All matrices, initial blocks and fields are identical across routes. Cold/warm
original-order routes do not perform unused descriptor sorting. Existing trace
construction remains charged internally, while matched output exports only
values/projectors; every solve still verifies fresh original residuals.

F uses the complete eight-point Rips supplied boundaries at threshold 4.0:
q0 is 8x28 and q1 is 28x56, with one or 512 leverage draws, seed zero, degree
eight and 32 identical frozen Rademacher probes. The current prototype still
pays for a full dense eigensolve, exact numerical leverage preparation, moments
and dense sampled assembly. It is not an optimized sparse-application API.
Independent rational leverage and integer RNG/CDF controls supply sampled truth.

G uses the weighted shell and the frozen enriched-12/20-extra-5 direct/core
controls. All weights and right-hand sides are identical. Both outputs are
checked against independent SVD minimum-norm solutions, and each measured
original-equation residual must be <= 1e-9. The largest observed residual is
6.325e-10. Singular ambient systems are not ordinarily inverted.

## Correctness and timer contracts

Each sample must match its independent expected vector/matrix/projector payload,
not merely the other implementation. The fixed scaled Euclidean/Frobenius gates
are 1e-10 for A/B/C/F, 1e-9 for D, 1e-8 for E and the existing 1e-7 SVD solution
gate for G, with the additional 1e-9 original residual gate. Shapes and finite
values are mandatory. Approximate F outputs are checked against their own exact
draw reconstruction; this does not turn quality failures into exact agreement.

The new adapters pass 31 standalone tests under stable and actual Rust 1.91.
Four controller regressions reject original-PR baselines, identical routes,
changed/nonfinite outputs and wrong shapes. Linux tool discovery passes 234
tests with 92 optional skips; Windows numerical reproduction discovery passes
108 tests without skips/failures. Fresh affected full replays retain 34 native
HeCS controls, 328 QR/Ritz steps, 48 sequence solves and 36 basis transfers.
Prior full integration debug/release/MSRV and independent seven-track results
remain bound to their original recorded commits; they are not relabeled here.

The primary worker timer starts at main and ends after flushed output, including
fixture I/O, source construction, preparation, computation, internal checks,
cleanup and matched export, but excluding process startup. Process wall time and
absolute process peak working set are reported separately. Internal preparation
and solve counters are retained only as phase diagnostics; they are not used
instead of the complete-workflow comparison.

## Equivalent-output workflow results

Milliseconds, median [minimum, maximum], ten measured processes per side.
**Speedup = median of paired unaccelerated / alternative times**; above one is
faster. No rows, outliers, setup costs or failed-quality sampling cases are
silently removed. F's non-equivalent approximation rows are in their own table.

| Workload | Unaccelerated ms | Alternative ms | Paired speedup |
| --- | --- | --- | --- |
| A-star64-actions1 | 0.368 [0.306, 0.559] | 0.268 [0.229, 0.415] | 1.350x |
| A-star64-actions100 | 0.591 [0.492, 0.919] | 0.383 [0.304, 0.569] | 1.383x |
| A-star256-actions1 | 2.154 [1.929, 3.252] | 0.618 [0.542, 0.965] | 3.332x |
| A-star256-actions100 | 6.450 [5.370, 9.049] | 0.976 [0.770, 1.559] | 6.338x |
| A-star512-actions1 | 7.745 [6.732, 8.376] | 1.004 [0.831, 1.647] | 7.267x |
| A-star512-actions100 | 25.310 [21.492, 26.566] | 1.632 [1.367, 2.782] | 14.540x |
| B-gram32 | 0.710 [0.623, 1.208] | 0.874 [0.725, 1.277] | 0.849x |
| B-harmonic32 | 0.695 [0.610, 1.221] | 0.885 [0.817, 1.150] | 0.805x |
| B-gram64 | 3.794 [3.439, 5.123] | 3.799 [3.112, 4.511] | 0.995x |
| B-harmonic64 | 4.636 [3.536, 6.137] | 5.244 [4.864, 6.817] | 0.918x |
| C-alpha-q0-matrices | 11.070 [9.811, 14.414] | 12.986 [11.352, 16.043] | 0.861x |
| C-alpha-q0-spectra | 24.239 [21.083, 31.741] | 29.098 [24.967, 32.130] | 0.884x |
| C-alpha-q1-matrices | 23.616 [20.607, 27.790] | 9.976 [8.642, 14.106] | 2.465x |
| C-alpha-q1-spectra | 44.848 [40.454, 45.877] | 25.616 [22.695, 29.731] | 1.747x |
| C-alpha-q1-single | 0.368 [0.262, 0.556] | 0.409 [0.361, 0.824] | 0.685x |
| D-nonbranching-matrices | 7.351 [6.234, 8.794] | 7.793 [6.558, 9.119] | 0.966x |
| D-nonbranching-spectra | 10.713 [8.516, 11.678] | 10.394 [8.754, 11.242] | 1.012x |
| E-smooth-warm | 10.459 [8.962, 14.742] | 10.595 [9.354, 12.836] | 1.002x |
| E-smooth-ordered | 10.048 [8.760, 11.369] | 9.418 [7.762, 11.312] | 1.008x |
| E-abrupt-warm | 11.040 [9.569, 12.295] | 10.214 [9.220, 14.738] | 1.041x |
| E-abrupt-ordered | 11.186 [8.206, 13.286] | 10.445 [8.415, 12.151] | 0.981x |
| G-shell-weighted-direct | 0.155 [0.135, 0.259] | 0.180 [0.150, 0.288] | 0.882x |
| G-enriched-12-5-direct | 0.905 [0.785, 1.116] | 0.592 [0.481, 0.768] | 1.463x |
| G-enriched-12-5-core | 0.953 [0.807, 1.825] | 0.570 [0.458, 0.822] | 1.882x |
| G-enriched-20-5-direct | 2.743 [2.169, 4.309] | 1.571 [1.181, 2.252] | 1.747x |
| G-enriched-20-5-core | 2.488 [2.187, 5.700] | 1.210 [1.082, 1.617] | 2.113x |

## Approximation cost and quality

These are cost ratios against an exact Gram, not equivalent-output speedups.
The quality column is the independently evaluated epsilon-0.5 PSD sandwich
and kernel/rank gate. One-draw failures remain visible and cannot be ranked as
valid replacements for the exact operator. No probability theorem is asserted.

| Workload | Exact ms | Prototype ms | Exact/prototype cost ratio | Relative sandwich error | Rank loss | Quality gate |
| --- | --- | --- | --- | --- | --- | --- |
| F-rips-4.0-q0-draws1 | 0.150 [0.130, 0.231] | 0.401 [0.272, 0.534] | 0.411 | 6.0000 | True | fail |
| F-rips-4.0-q0-draws512 | 0.204 [0.125, 0.246] | 0.360 [0.291, 0.540] | 0.465 | 0.1208 | False | pass |
| F-rips-4.0-q1-draws1 | 0.448 [0.324, 0.657] | 1.572 [1.385, 2.203] | 0.279 | 20.0000 | True | fail |
| F-rips-4.0-q1-draws512 | 0.406 [0.351, 0.619] | 1.682 [1.567, 2.476] | 0.238 | 0.3069 | False | pass |

## Complete process costs

Wall time covers process creation through full transport. Peak working set is
the absolute Windows per-process peak, including runtime/input/allocator/output,
not live algorithm allocations. MiB = 2^20 bytes. Values are median [min, max].

| Workload | Baseline wall ms | Alternative wall ms | Baseline peak MiB | Alternative peak MiB |
| --- | --- | --- | --- | --- |
| A-star64-actions1 | 10.112 [9.524, 11.698] | 10.615 [9.891, 11.260] | 4.551 [4.547, 4.555] | 4.570 [4.566, 4.578] |
| A-star64-actions100 | 10.478 [9.078, 11.528] | 10.273 [8.919, 11.686] | 4.551 [4.547, 4.559] | 4.570 [4.547, 4.578] |
| A-star256-actions1 | 12.382 [10.605, 14.818] | 10.382 [8.958, 13.415] | 5.164 [5.152, 5.172] | 4.561 [4.523, 4.570] |
| A-star256-actions100 | 17.274 [14.170, 18.936] | 11.101 [9.230, 12.829] | 5.168 [5.160, 5.172] | 4.557 [4.527, 4.566] |
| A-star512-actions1 | 18.294 [16.183, 20.133] | 11.285 [9.014, 12.665] | 8.299 [8.293, 8.301] | 4.576 [4.566, 4.582] |
| A-star512-actions100 | 35.072 [31.452, 39.338] | 11.711 [10.213, 13.820] | 8.297 [8.293, 8.305] | 4.576 [4.566, 4.582] |
| B-gram32 | 10.386 [8.991, 12.454] | 10.504 [9.701, 12.361] | 4.480 [4.477, 4.484] | 4.531 [4.516, 4.535] |
| B-harmonic32 | 10.468 [9.222, 11.997] | 10.252 [8.954, 11.692] | 4.479 [4.473, 4.484] | 4.549 [4.543, 4.562] |
| B-gram64 | 13.502 [11.939, 15.495] | 13.119 [11.543, 15.638] | 4.449 [4.441, 4.469] | 4.508 [4.500, 4.527] |
| B-harmonic64 | 14.472 [11.918, 18.377] | 15.099 [14.325, 19.618] | 4.453 [4.445, 4.457] | 4.512 [4.508, 4.520] |
| C-alpha-q0-matrices | 22.814 [17.782, 28.561] | 23.691 [21.395, 28.757] | 5.902 [5.898, 5.953] | 7.367 [7.340, 7.402] |
| C-alpha-q0-spectra | 34.827 [30.593, 42.286] | 39.557 [33.496, 44.847] | 5.982 [5.953, 5.992] | 7.381 [7.363, 7.406] |
| C-alpha-q1-matrices | 33.970 [30.692, 38.998] | 20.243 [17.006, 26.524] | 6.152 [6.148, 6.250] | 6.980 [6.973, 7.004] |
| C-alpha-q1-spectra | 54.787 [50.563, 59.412] | 36.844 [31.430, 41.078] | 6.164 [6.160, 6.207] | 7.564 [7.535, 7.598] |
| C-alpha-q1-single | 10.262 [9.336, 12.070] | 10.356 [8.540, 12.722] | 4.551 [4.531, 4.555] | 4.545 [4.539, 4.570] |
| D-nonbranching-matrices | 16.568 [15.410, 19.985] | 17.936 [16.035, 20.489] | 4.668 [4.648, 4.695] | 4.691 [4.688, 4.703] |
| D-nonbranching-spectra | 20.559 [18.143, 25.333] | 19.927 [16.235, 21.609] | 4.707 [4.656, 4.711] | 4.717 [4.707, 4.828] |
| E-smooth-warm | 19.585 [17.873, 26.149] | 20.504 [18.596, 23.093] | 4.824 [4.809, 5.016] | 4.988 [4.801, 5.012] |
| E-smooth-ordered | 19.891 [18.380, 22.877] | 19.997 [17.520, 22.816] | 4.996 [4.820, 5.027] | 4.982 [4.812, 5.004] |
| E-abrupt-warm | 20.674 [18.701, 23.852] | 20.998 [18.726, 23.921] | 4.963 [4.770, 5.008] | 4.984 [4.812, 5.008] |
| E-abrupt-ordered | 21.226 [17.232, 23.454] | 20.330 [17.046, 22.670] | 4.969 [4.766, 5.016] | 4.900 [4.801, 5.008] |
| F-rips-4.0-q0-draws1 | 10.283 [8.892, 13.348] | 11.153 [8.989, 13.390] | 4.488 [4.484, 4.492] | 4.545 [4.531, 4.551] |
| F-rips-4.0-q0-draws512 | 10.835 [8.725, 12.796] | 9.942 [8.633, 13.556] | 4.486 [4.473, 4.496] | 4.543 [4.543, 4.551] |
| F-rips-4.0-q1-draws1 | 10.473 [8.854, 13.632] | 11.912 [10.162, 15.102] | 4.512 [4.504, 4.516] | 4.521 [4.512, 4.527] |
| F-rips-4.0-q1-draws512 | 10.431 [8.925, 12.528] | 11.466 [9.882, 13.859] | 4.508 [4.496, 4.520] | 4.531 [4.527, 4.535] |
| G-shell-weighted-direct | 9.784 [8.747, 11.514] | 9.882 [8.857, 11.518] | 4.336 [4.332, 4.340] | 4.520 [4.516, 4.527] |
| G-enriched-12-5-direct | 10.260 [8.928, 15.742] | 11.078 [9.857, 14.702] | 4.430 [4.391, 4.438] | 4.582 [4.578, 4.590] |
| G-enriched-12-5-core | 10.676 [9.817, 13.409] | 10.058 [8.490, 13.313] | 4.432 [4.406, 4.438] | 4.578 [4.551, 4.582] |
| G-enriched-20-5-direct | 12.671 [10.869, 16.824] | 11.813 [10.679, 14.412] | 4.598 [4.508, 4.602] | 4.566 [4.539, 4.574] |
| G-enriched-20-5-core | 13.594 [10.884, 17.079] | 11.953 [9.390, 14.431] | 4.592 [4.512, 4.602] | 4.543 [4.539, 4.555] |

## Interpretation

- A's 512-leaf/100-action row changes from 25.31 to 1.63 ms (14.54x);
  median process peak changes from 8.30 to 4.58 MiB. Single-use gains are smaller.
- C's q1 matrix-only ratio is 2.47x and complete-spectrum ratio is 1.75x. Preparation
  does not repay on q0 or the single-pair case (0.68x). Memory is a tradeoff:
  q1 complete-spectrum peak changes from 6.16 to 7.56 MiB.
- B's current both-Gram API is 0.85x/1.00x and harmonic restriction is
  0.80x/0.92x. Retaining the large eigensolve and validation prevents a consistent
  Gram/restriction acceleration claim at this stage.
- D's supported whole-suite ratios are 0.97x/1.01x, without/with spectra.
  Mathematical eligibility alone does not establish a performance win.
- E's smooth ordered route reduces counted actions from 3,184 to 2,774 and
  iterations from 56 to 49, with a 1.01x observed workflow ratio. Abrupt reuse
  retains 3,104 versus 3,074 cold actions and no iteration reduction; its ratios
  near one do not demonstrate a robust benefit.
- G's enriched cases reduce 43 to five or 63 to three CGLS updates, while full
  workflow ratios are 1.46-2.11x after preparation. The small shell is
  0.88x despite one fewer update; iteration ratios are not runtime speedups.
- F's one-draw cases lose rank and fail the quality gate. At 512 draws the
  epsilon-0.5 gate passes, but the current prototype is slower than exact dense
  assembly (cost ratios 0.46 and 0.24). These tiny dense controls do not
  support a sparsification acceleration or production-adoption claim.

Small sub-millisecond rows are sensitive to measurement and OS variation; read
their absolute ranges, not only ratios. There is no automatic default dispatch,
universal speedup, asymptotic or total-memory guarantee. The seven tracks are
exercised in the integrated build; their gains are not multiplied or added into
an imaginary all-track pipeline. The retained LSMR terminal-breakdown exclusions,
basic-solver failures, branching exclusions and RP2 HeCS rank rejection remain
correctness boundaries rather than fabricated successful acceleration cases.

## Evidence and reproduction

Raw evidence is **local-only** under
`target/benchmarks/commit-8994b62c8e2d/spectral-unaccelerated/run-001/`.
Member-manifest SHA-256: `0d411719ec8288e92ea3b8519517d523b8ae6b364ba5f8eacc33a44f1a83356b`.
All specifications, oracle payloads, frozen input bytes, source/build/binary
identities, raw stdout/stderr, warmups, samples and complete summaries are
retained. The corrected independent preflight bundle is under
`target/integration-20261124/unaccelerated-inputs-001/`. No durable external
performance archive is configured; these raw samples are not claimed publicly
downloadable. Hosted correctness artifacts are separate and expire after 14 days.

With the frozen specification/input/oracle bundle, a clean measured checkout,
and an idle local build/test environment, replay into a new directory:

```sh
python tools/benchmark_spectral_integration.py --specification target/unaccelerated-plan.json --output target/benchmarks/commit-CANDIDATE/spectral-unaccelerated/run-002
```

The fixed matrix was rerun after preserving finite-field-input validation on
unsorted sequence routes. The earlier unaccelerated run at
`a02927df3b9ebf119b406c488c345b6ba96a403d` remains under its own source identity;
its samples are not pooled with or relabeled as the final guarded-source run.

The [tools guide](../../tools/README.md#integrated-spectral-resource-snapshots)
owns route syntax and comparison checks. The original PR-head experiment remains
under its original `9588668972873f4cf71860d9ccac7b0738c2cdd0` raw path solely as
withdrawn integration-regression evidence. Its numbers do not answer the
unaccelerated-baseline question. Final hosted status and exact PR head belong in
PR #160; a later main merge or release is distinct from this measured source.
