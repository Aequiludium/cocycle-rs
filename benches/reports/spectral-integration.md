# Integrated spectral routes versus unaccelerated algorithms

[Benchmarks](../README.md) / [Reporting rules](../reporting.md)

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

The fixed-projection action route achieves 1.50-15.34x on the declared star
workloads. Relevant-pair q1 matrix construction achieves 2.68x, reducing to 1.59x
when complete spectra are also requested. Enriched HeCS controls achieve
1.48-2.10x including preparation. These gains do not generalize to every row:
q0 batching, single-pair batching, the non-branching suite and several Gram,
restriction and reuse controls do not improve. No aggregate speedup or faster
default is inferred.

## Measured identity and environment

| Item | Recorded value |
| --- | --- |
| Measured library, workers and controller | `a02927df3b9ebf119b406c488c345b6ba96a403d`, clean committed integration-branch source |
| Baseline revision | The same commit; explicit unaccelerated route selection below, never an original PR head |
| Controller SHA-256 | `93b8d8da2e5e101564b69993202b02006283ddb5717dc63b480133c0ef95fa1c` |
| Protocol | `spectral-unaccelerated-comparison-v2`; [controller](../../tools/benchmark_spectral_integration.py) |
| Compiler | `rustc 1.98.1 (48a229cea 2026-09-01)`; edition 2024 and `-C opt-level=3` for every worker and the shared library |
| Machine | Intel(R) Core(TM) Ultra 7 155H; Windows-11-10.0.26100-SP0 |
| Python | 3.12.0rc3 (tags/v3.12.0rc3:b973ab3, Sep 19 2023, 09:39:54) [MSC v.1935 64 bit (AMD64)] |
| Measured UTC | 2026-11-24T12:05:49.377432+00:00 through 2026-11-24T12:06:03.709328+00:00 |
| Repetition/order | One discarded warmup and ten fresh processes per route/cell; seed 160, shuffled cells, balanced pair order |
| Controls | All builds finish before serial measurement; thread environment one; no concurrent local tests/builds |
| Limits | 120-second process timeout; no address-space cap; frequency, affinity and host load uncontrolled |

Sources, generated timing wrappers, compiler invocations, binaries and inputs are
hashed in the raw member manifest. The report commit is distinct from the
measured source. The Rust library `src/` is unchanged from the independently
replayed `eacf0c752914619e9f3eaaab830473aa854b35dd` integration source; changes here
are comparison adapters, benchmark-only output modes, checks and documentation.

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
| A-star64-actions1 | 0.270 [0.229, 0.539] | 0.170 [0.159, 0.213] | 1.507x |
| A-star64-actions100 | 0.388 [0.352, 0.581] | 0.259 [0.214, 0.445] | 1.495x |
| A-star256-actions1 | 1.534 [1.351, 2.792] | 0.423 [0.393, 0.500] | 3.573x |
| A-star256-actions100 | 4.477 [4.081, 6.442] | 0.678 [0.570, 0.994] | 6.867x |
| A-star512-actions1 | 5.621 [5.087, 6.786] | 0.723 [0.671, 2.118] | 7.371x |
| A-star512-actions100 | 20.038 [16.913, 23.145] | 1.298 [1.036, 2.185] | 15.341x |
| B-gram32 | 0.541 [0.458, 0.679] | 0.599 [0.509, 1.176] | 0.853x |
| B-harmonic32 | 0.535 [0.465, 0.859] | 0.736 [0.588, 1.164] | 0.801x |
| B-gram64 | 3.085 [2.680, 4.905] | 2.883 [2.587, 3.335] | 1.069x |
| B-harmonic64 | 3.515 [2.648, 5.270] | 4.177 [3.474, 6.021] | 0.826x |
| C-alpha-q0-matrices | 8.892 [7.974, 10.765] | 11.322 [9.537, 13.914] | 0.843x |
| C-alpha-q0-spectra | 20.776 [17.573, 25.419] | 21.590 [20.572, 26.564] | 0.918x |
| C-alpha-q1-matrices | 21.868 [18.468, 27.973] | 8.197 [7.322, 9.760] | 2.685x |
| C-alpha-q1-spectra | 38.443 [32.616, 44.198] | 23.449 [20.823, 28.412] | 1.593x |
| C-alpha-q1-single | 0.208 [0.171, 0.248] | 0.281 [0.258, 0.331] | 0.689x |
| D-nonbranching-matrices | 5.905 [5.114, 7.626] | 6.172 [5.309, 7.769] | 0.946x |
| D-nonbranching-spectra | 9.081 [7.675, 12.552] | 9.172 [8.660, 13.873] | 0.946x |
| E-smooth-warm | 8.050 [7.745, 9.646] | 8.576 [7.682, 10.487] | 1.004x |
| E-smooth-ordered | 8.318 [7.573, 13.774] | 7.762 [6.890, 11.827] | 1.161x |
| E-abrupt-warm | 8.041 [7.396, 10.293] | 8.603 [7.346, 11.757] | 0.969x |
| E-abrupt-ordered | 8.658 [7.734, 12.120] | 8.825 [7.246, 9.496] | 1.022x |
| G-shell-weighted-direct | 0.104 [0.093, 0.196] | 0.127 [0.113, 0.146] | 0.854x |
| G-enriched-12-5-direct | 0.715 [0.618, 1.526] | 0.436 [0.371, 0.969] | 1.481x |
| G-enriched-12-5-core | 0.703 [0.581, 1.159] | 0.413 [0.351, 0.822] | 1.562x |
| G-enriched-20-5-direct | 1.969 [1.628, 2.648] | 1.044 [0.866, 2.146] | 1.851x |
| G-enriched-20-5-core | 1.826 [1.646, 3.779] | 0.857 [0.791, 1.333] | 2.102x |

## Approximation cost and quality

These are cost ratios against an exact Gram, not equivalent-output speedups.
The quality column is the independently evaluated epsilon-0.5 PSD sandwich
and kernel/rank gate. One-draw failures remain visible and cannot be ranked as
valid replacements for the exact operator. No probability theorem is asserted.

| Workload | Exact ms | Prototype ms | Exact/prototype cost ratio | Relative sandwich error | Rank loss | Quality gate |
| --- | --- | --- | --- | --- | --- | --- |
| F-rips-4.0-q0-draws1 | 0.099 [0.084, 0.229] | 0.224 [0.193, 0.337] | 0.441 | 6.0000 | True | fail |
| F-rips-4.0-q0-draws512 | 0.106 [0.087, 0.139] | 0.235 [0.221, 0.511] | 0.416 | 0.1208 | False | pass |
| F-rips-4.0-q1-draws1 | 0.275 [0.220, 0.522] | 1.192 [1.035, 2.300] | 0.223 | 20.0000 | True | fail |
| F-rips-4.0-q1-draws512 | 0.252 [0.197, 0.373] | 1.320 [1.101, 1.819] | 0.187 | 0.3069 | False | pass |

## Complete process costs

Wall time covers process creation through full transport. Peak working set is
the absolute Windows per-process peak, including runtime/input/allocator/output,
not live algorithm allocations. MiB = 2^20 bytes. Values are median [min, max].

| Workload | Baseline wall ms | Alternative wall ms | Baseline peak MiB | Alternative peak MiB |
| --- | --- | --- | --- | --- |
| A-star64-actions1 | 6.834 [5.223, 8.741] | 6.158 [5.148, 7.639] | 4.551 [4.527, 4.566] | 4.570 [4.551, 4.578] |
| A-star64-actions100 | 7.397 [5.841, 9.794] | 7.230 [5.401, 8.118] | 4.543 [4.512, 4.559] | 4.570 [4.555, 4.582] |
| A-star256-actions1 | 8.541 [6.486, 11.219] | 7.029 [5.824, 8.989] | 5.168 [5.156, 5.180] | 4.553 [4.520, 4.562] |
| A-star256-actions100 | 11.354 [10.489, 14.266] | 7.925 [6.270, 9.372] | 5.164 [5.152, 5.172] | 4.555 [4.516, 4.566] |
| A-star512-actions1 | 12.969 [11.474, 16.112] | 7.411 [6.669, 11.341] | 8.295 [8.281, 8.305] | 4.572 [4.562, 4.586] |
| A-star512-actions100 | 27.145 [23.082, 30.029] | 7.804 [6.583, 11.924] | 8.297 [8.281, 8.305] | 4.566 [4.543, 4.578] |
| B-gram32 | 7.012 [5.543, 8.464] | 7.696 [5.736, 9.501] | 4.484 [4.461, 4.488] | 4.531 [4.516, 4.539] |
| B-harmonic32 | 7.854 [5.779, 12.676] | 7.901 [6.889, 10.616] | 4.480 [4.473, 4.484] | 4.551 [4.543, 4.562] |
| B-gram64 | 9.774 [8.573, 12.017] | 9.474 [8.607, 11.818] | 4.449 [4.445, 4.473] | 4.504 [4.492, 4.512] |
| B-harmonic64 | 10.819 [8.538, 12.787] | 11.210 [9.633, 13.706] | 4.449 [4.441, 4.457] | 4.512 [4.500, 4.523] |
| C-alpha-q0-matrices | 16.123 [14.589, 20.113] | 18.675 [16.719, 23.385] | 5.910 [5.898, 5.949] | 7.371 [7.344, 7.387] |
| C-alpha-q0-spectra | 27.806 [24.549, 33.522] | 28.399 [27.328, 33.162] | 5.984 [5.930, 6.000] | 7.393 [7.352, 7.422] |
| C-alpha-q1-matrices | 28.563 [26.301, 38.377] | 15.882 [13.626, 17.960] | 6.146 [6.074, 6.223] | 6.986 [6.973, 6.992] |
| C-alpha-q1-spectra | 45.990 [38.629, 55.222] | 30.233 [27.574, 39.736] | 6.180 [6.078, 6.301] | 7.578 [7.535, 7.605] |
| C-alpha-q1-single | 6.802 [5.756, 10.033] | 6.457 [5.875, 8.626] | 4.551 [4.547, 4.562] | 4.549 [4.527, 4.570] |
| D-nonbranching-matrices | 12.894 [11.482, 15.695] | 13.517 [11.673, 18.634] | 4.678 [4.641, 4.699] | 4.693 [4.641, 4.707] |
| D-nonbranching-spectra | 15.749 [14.979, 24.811] | 16.829 [15.195, 21.672] | 4.689 [4.668, 4.719] | 4.711 [4.680, 4.730] |
| E-smooth-warm | 14.281 [13.253, 18.238] | 15.103 [13.103, 19.280] | 4.986 [4.941, 5.004] | 4.973 [4.770, 4.992] |
| E-smooth-ordered | 15.200 [12.906, 23.398] | 13.779 [12.754, 21.939] | 4.986 [4.820, 5.027] | 4.980 [4.809, 4.996] |
| E-abrupt-warm | 14.760 [13.266, 17.052] | 15.932 [13.213, 18.770] | 4.973 [4.797, 5.008] | 4.986 [4.836, 4.996] |
| E-abrupt-ordered | 15.091 [13.600, 19.047] | 15.713 [13.522, 17.328] | 4.910 [4.586, 4.992] | 4.975 [4.805, 4.980] |
| F-rips-4.0-q0-draws1 | 6.956 [5.199, 8.962] | 7.482 [5.821, 15.116] | 4.488 [4.480, 4.508] | 4.543 [4.535, 4.555] |
| F-rips-4.0-q0-draws512 | 7.535 [6.173, 16.071] | 6.736 [5.459, 10.808] | 4.488 [4.484, 4.500] | 4.549 [4.543, 4.559] |
| F-rips-4.0-q1-draws1 | 6.811 [5.999, 9.273] | 8.192 [6.095, 10.927] | 4.506 [4.492, 4.520] | 4.520 [4.500, 4.535] |
| F-rips-4.0-q1-draws512 | 6.668 [5.619, 9.447] | 7.984 [7.318, 9.983] | 4.508 [4.496, 4.527] | 4.527 [4.520, 4.535] |
| G-shell-weighted-direct | 6.860 [5.563, 10.904] | 6.900 [5.870, 8.871] | 4.336 [4.332, 4.336] | 4.521 [4.512, 4.527] |
| G-enriched-12-5-direct | 7.844 [5.875, 10.908] | 7.166 [5.909, 12.062] | 4.426 [4.375, 4.441] | 4.580 [4.527, 4.586] |
| G-enriched-12-5-core | 7.133 [5.711, 9.341] | 6.913 [6.198, 9.660] | 4.432 [4.375, 4.438] | 4.578 [4.531, 4.582] |
| G-enriched-20-5-direct | 8.646 [7.606, 11.282] | 7.535 [6.744, 12.271] | 4.594 [4.520, 4.598] | 4.561 [4.543, 4.574] |
| G-enriched-20-5-core | 7.985 [7.518, 11.327] | 7.653 [6.808, 9.880] | 4.598 [4.508, 4.602] | 4.547 [4.531, 4.555] |

## Interpretation

- A's 512-leaf/100-action row improves from about 20.04 to 1.30 ms (15.34x);
  median process peak falls from 8.30 to 4.57 MiB. Single-use gains are smaller.
- C's q1 matrix-only gain is 2.68x and complete-spectrum gain is 1.59x. Preparation
  does not repay on q0 or the single-pair case (0.69x). Memory is a tradeoff:
  q1 complete-spectrum peak rises from 6.18 to 7.58 MiB.
- B's current both-Gram API is about 0.85x/1.07x and harmonic restriction about
  0.80x/0.83x. Retaining the large eigensolve and validation prevents a consistent
  Gram/restriction acceleration claim at this stage.
- D's supported whole-suite ratios are about 0.95x, with and without spectra.
  Mathematical eligibility alone does not establish a performance win.
- E's smooth ordered route reduces counted actions from 3,184 to 2,774 and
  iterations from 56 to 49, with a 1.16x observed workflow ratio. Abrupt reuse
  retains 3,104 versus 3,074 cold actions and no iteration reduction; its ratios
  near one do not demonstrate a robust benefit.
- G's enriched cases reduce 43 to five or 63 to three CGLS updates, while full
  workflow gains are only 1.48-2.10x after preparation. The small shell is slower
  at 0.85x despite one fewer update; iteration ratios are not runtime speedups.
- F's one-draw cases lose rank and fail the quality gate. At 512 draws the
  epsilon-0.5 gate passes, but the current prototype is slower than exact dense
  assembly (cost ratios about 0.42 and 0.19). These tiny dense controls do not
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
`target/benchmarks/commit-a02927df3b9e/spectral-unaccelerated/run-001/`.
Member-manifest SHA-256: `e6237b6dc3d67d1de0e7569a2813b295de9250c6b9bf8ab1ad729f1888b70358`.
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

The [tools guide](../../tools/README.md#integrated-spectral-resource-snapshots)
owns route syntax and comparison checks. The original PR-head experiment remains
under its original `9588668972873f4cf71860d9ccac7b0738c2cdd0` raw path as
integration-regression evidence excluded from this acceleration comparison. Its numbers do not answer the
unaccelerated-baseline question. Final hosted status and exact PR head belong in
PR #160; a later main merge or release is distinct from this measured source.
