# Persistent Laplacian reference resources

Evidence class: **resource snapshot** under the [performance reporting rules](../reporting.md).

The dense-construction and fixed-action sections below have distinct measured
revisions and workloads; their timings must not be combined.

## Dense construction snapshot

PR: [#131](https://github.com/Aequiludium/cocycle-rs/pull/131).
Measured candidate: `90a861da2ca368cee4d714919cb7f1904a3dd011` (feature-branch commit, not a tested merge).
Baseline: not applicable; no previous-head speedup or external native ranking is claimed.
The report commit is separate from the measured implementation/controller commit.

Source-manifest SHA-256: `75b53257a73418f319db1b3d578a2661f9e19ab6b2a449646447f3ec712bf761`.
Worker-manifest SHA-256: `6dc151eceb2543dff5c352bbe2a9b66b3534b31c5f5fa8b287226254a8a67d62`.
These are SHA-256 of the recorded path-to-file-hash maps encoded as UTF-8 JSON
with sorted keys and compact separators; paths retain their recorded backslashes.
Controller SHA-256: `d493b2b3f01de9b8af12958b10ba12ea29e026256bbf4c0a3c58d7f1ab9b7b71`.
Specification SHA-256: `a68608cb460b25eec3b2a277c52244bb523b76d26452819893a4ac6df446dd7f`.
The controller and workers belong to the measured commit above.

Host: Intel(R) Core(TM) Ultra 7 155H, Windows-11-10.0.26100-SP0; Rust 1.98.1, Python 3.12.14.
Oracle/workflow packages: NumPy 2.2.6, SciPy 1.15.3, GUDHI 3.12.0;
learning additionally uses scikit-learn 1.7.2 and threadpoolctl 3.6.0.
OMP/OpenBLAS/MKL thread environment is one. CPU affinity, frequency and
external host load were uncontrolled. Workers ran serially without concurrent
local compilation/tests; reproduction-internal builds are part of that cell.
Recorded plan UTC: `2026-10-15T04:06:25.624132+00:00` to `2026-10-15T04:06:30.710860+00:00`.

Protocol: `spectral-resource-v1`; one retained discarded warmup and ten fresh
measured processes per cell, seeded order 73, no excluded outliers/retries.
Native workers use `rustc --edition 2024 -C opt-level=3` against the locked
release library. Reproduction uses its unchanged internal debug Rust builds.
The table gives median (minimum-maximum); times are milliseconds and memory
is MiB. Peak working set includes runtime/input/allocator/retained output,
rather than live algorithm allocation. Numerical/storage failures are retained.

| Cell | Outcomes / 10 | Solve ms | Whole process ms | Peak working set MiB |
| --- | --- | --- | --- | --- |
| path-64-starts24-direct | success: 10 | 40.606 (37.424-47.825) | 52.461 (49.806-61.707) | 4.576 (4.551-4.613) |
| path-64-starts24-stepwise | success: 10 | 0.632 (0.568-0.925) | 12.120 (10.880-13.207) | 4.529 (4.520-4.535) |
| path-128-starts8-direct | success: 10 | 98.311 (94.376-100.933) | 111.408 (106.586-114.574) | 4.674 (4.668-4.688) |
| path-128-starts8-stepwise | success: 10 | 2.781 (2.642-3.164) | 14.516 (12.539-18.303) | 4.660 (4.637-4.668) |
| path-256-starts1-direct | success: 10 | 59.595 (58.678-64.980) | 71.835 (70.076-78.193) | 5.158 (5.141-5.164) |
| path-256-starts1-stepwise | success: 10 | 22.116 (21.501-25.704) | 35.699 (33.030-37.203) | 4.807 (4.766-4.816) |
| isolates-64-starts12-direct | success: 10 | 2.530 (2.373-3.526) | 14.000 (12.575-18.544) | 4.488 (4.465-4.523) |
| isolates-64-starts12-stepwise | success: 10 | 0.287 (0.278-0.400) | 11.157 (10.345-14.098) | 4.473 (4.457-4.488) |
| path-257-starts1-direct | storage rejected: 10 | not returned | 12.166 (10.361-12.985) | 4.469 (4.457-4.480) |
| path-257-starts1-stepwise | storage rejected: 10 | not returned | 11.778 (9.902-13.383) | 4.473 (4.457-4.488) |

These are q=0 stage-pair matrix constructions, not complete-spectrum timings.
Vertex i is born at scale i; path edges are all born at scale n. The target is n.
For k starts, source scales are the distinct floors i*(n-1)/k, i=1,...,k; the
one-start case uses floor(n/2). Isolates have no edges. Direct and reverse routes
retain the same ordered matrices, with source setup excluded from the worker's
construction timer. Result allocation is included and destruction excluded.
Whole-process time additionally includes input parsing/preparation/startup/export.
Each measured supported fixture was independently checked against the Python
restricted-boundary SVD oracle and direct/reverse matrices before timing.

The example's benchmark mode now uses the library's actual default 65,536-entry
per-buffer/cumulative-output limit. Both 257-coordinate cells reject storage in
all ten processes with InvalidParameter(max_dense_entries); their process times
are rejection costs, never successful solve times. Matrices at 256 coordinates
use 65,536 entries just for the target Gram; coexisting buffers and owned source
metadata are additional memory. This remains a bounded dense reference, not an
accepted high-performance persistent spectral kernel.

Near-zero pivots cannot silently discard coupled corrections. The scalar helper
rejects nonzero coupling and the fixed-endpoint route uses the checked direct
Schur fallback under the same execution budget and tolerances. The numerical
Jacobi contract is deliberately separate from the Hodge/Dirac Ritz backend.
The general PSD counterexample is not claimed to be realizable by an unweighted
simplicial source. Generated corpus fallback flags and all adverse outcomes
remain in raw correctness/resource evidence. No timing assertion is a unit test.

The native workers use default `Execution` without cancellation/work caps;
the documented algorithm workspace, block and iteration policies still apply.
No process timeout or OS address-space cap was imposed. Peak memory is an
observation, never a hard library quota.

This complete cell matrix was rerun after the additional direct-Schur
correction guard at the measured commit above. The earlier implementation
run remains local and is not mixed into these current statistics. Fresh
correctness replay passes 1,114 pairs; 13 full/upper pairs explicitly use
direct fallback. Discarded positive-direction correction bounds are
checked in addition to actual solve residuals, with unchanged tolerances.

## Reproduction and retained evidence

The source checkout provides [the controller](../../tools/benchmark_spectral_resources.py).
Supply a JSON list whose entries declare `name`, `kind`, `command`, explicit
fixture data/parameters, and independent `expected_values`, `atol` and `rtol`
for native eigenpairs. Persistent entries instead supply the exact whitespace
input in `stdin` and expected row-sum checksum zero. Reproduction entries use
`{attempt}` as a fresh output-directory placeholder. The fixed cells and all
generator formulas/options are specified above and in the source workers.

```sh
python tools/benchmark_spectral_resources.py --specification target/specification.json --output target/benchmarks/commit-<measured-sha12>/spectral-resources/run-002 --samples 10
```

The actual local attempt is `target/benchmarks/commit-90a861da2ca3/spectral-resources/run-001`.
It retains the exact specification/fixture hash, immutable source/worker
fingerprints, binary hashes, environment, raw stdout/stderr, every warmup and
sample, complete sample matrix, summaries and file manifest. Companion plan metadata
records dependency versions, CPU, UTC start/end and generator hash.
Raw evidence is **local only**; no durable public archive is configured.
Committed source/settings allow replay, but the original samples are not
claimed publicly retrievable. CI correctness artifacts are a separate,
14-day evidence stream and are not performance baselines.

## Fixed-projection action: dense versus sparse incidence

Issue [#144](https://github.com/Aequiludium/cocycle-rs/issues/144). Measured
candidate [b1257aff7e477249bb1fdc2e24d5ff4da18391da](https://github.com/Aequiludium/cocycle-rs/commit/b1257aff7e477249bb1fdc2e24d5ff4da18391da),
baseline [b70e984657fdb3fc30ecb60a9662b6c86cf9a6a5](https://github.com/Aequiludium/cocycle-rs/commit/b70e984657fdb3fc30ecb60a9662b6c86cf9a6a5).
The baseline is the unoptimized Rust rewrite, after the
[paper reproduction](../../docs/research/persistent-constructions.md).
The candidate changes only incidence storage/traversal to sparse signed columns;
the fixed row-space projector, input metadata and numerical policy agree.
The report commit is separate from both measured commits.

The bounded admission gates passed: 801 project pair matrices/actions against
exact rational and SVD expectations, explicit failure/resource tests, unchanged
dense limits at a nontrivial 257-coordinate action, and lower incidence storage
and action times on the declared star controls. Keep this as a caller-selected
action route. Complete spectra, general sparse-factor scaling, ill-conditioned
rank decisions and end-to-end superiority are not established.

### Identity, workload and protocol

Windows `10.0.26100`, Intel64 Family 6 Model 170 Stepping 4; Rust/Cargo 1.98.1,
Windows MSVC release build, Python 3.10.11. Both workers were built using
`cargo build --locked --release --example persistent_laplacian` at their clean
committed sources. Baseline executable SHA-256:
`d38f161c12757889678f711ee3a890be9494d421b85d15fda6ff7c457d872271`.
Source-map SHA-256 (same sorted compact JSON convention as above):
`bc3889b08f92496bc04c7d81daf46aa3d0b6666c8e06515cb562fb45620c69cb`.
Controller SHA-256:
`9ed5afd99bcb91aa9965279acf118905b5a62742298cfa1089ba12477685e312`.
Specification SHA-256:
`602dd49ac7d59b64a40b39131cef7e65730b54c34a7ce63aed967deaba8ad457`.
The controller is from baseline b70e984; candidate retained the same controller.
The environment/specification retain all source, worker, binary and fixture hashes.

Protocol `spectral-resource-v1`, order seed 73; one retained warmup and ten fresh
serial processes per cell, no excluded samples. OMP/OpenBLAS/MKL thread settings
are one; no concurrent local builds/tests. CPU affinity, frequency and host load
are uncontrolled, with no timeout or OS address-space cap. Artifact creation/
completion timestamps bound execution to 2026-11-24 09:03:32-09:03:36 UTC; these
are filesystem timestamps, not separately instrumented wall-clock events.

K consists of n leaves born at zero. L adds a center born at one and n edges
born at two; q=0, pair (0,2), standard real coefficients, no weights. All n
upper columns are retained and the outside boundary has rank one. Input vector
x_i=(i mod 17)/17 in the source's ordered leaf coordinates. The independent
expectation is x-mean(x), i.e. U=I-11^T/n. Each measured action checks every
coordinate against it before accumulating the squared-norm checksum. The output
scope is one owned operator plus one or 1000 actions, with no matrix or spectrum.
Default dense cap is 65,536 for both versions; Q retains n entries. Incidence
storage falls from n*n f64 entries to n signed index/coefficient terms, excluding
Q, vectors, source metadata and container overhead.

Source creation, operator setup, actions and whole-process time are separate.
Setup includes incidence/projection allocation and numerical checks. Action time
includes output allocation, analytic validation and checksum accumulation; it
excludes setup. Whole-process time includes startup/source creation/export and
cleanup. Peak working set is absolute process memory, not live allocation peak.

### Complete measured matrix

All twelve cells passed all ten samples. Times are ms, median (minimum-maximum).
Paired ratios are candidate/baseline for the same repeat; below one means less
action time. All raw samples and warmups remain retained.

| Leaves | Actions | Dense action ms | Sparse action ms | Paired ratio median (range) | Whole process dense / sparse median ms |
| --- | --- | --- | --- | --- | --- |
| 64 | 1 | 0.0062 (0.0052-0.0148) | 0.0049 (0.0036-0.0081) | 0.702 (0.534-0.982) | 7.670 / 7.862 |
| 64 | 1000 | 3.107 (2.577-8.637) | 0.921 (0.846-2.385) | 0.320 (0.265-0.620) | 11.163 / 9.033 |
| 128 | 1 | 0.0171 (0.0141-0.0434) | 0.0050 (0.0042-0.0107) | 0.294 (0.209-0.473) | 8.533 / 8.129 |
| 128 | 1000 | 15.432 (10.832-29.474) | 1.636 (1.267-3.768) | 0.122 (0.090-0.250) | 23.686 / 10.691 |
| 256 | 1 | 0.0557 (0.0455-0.1972) | 0.0061 (0.0052-0.0146) | 0.107 (0.039-0.248) | 8.252 / 8.893 |
| 256 | 1000 | 56.897 (44.238-118.131) | 2.459 (2.174-7.018) | 0.045 (0.032-0.061) | 65.266 / 10.919 |

Setup median dense/sparse ms for one action: 64 leaves 0.0423/0.0441;
128 leaves 0.1141/0.0581; 256 leaves 0.2129/0.0850. Thus the smallest setup and
the one-action process medians at 64/256 regress. Median peak working set
dense/sparse bytes for one action: 4,745,216/4,751,360 (64),
4,732,928/4,739,072 (128), 4,956,160/4,741,120 (256). For 1000 actions:
4,743,168/4,747,264; 4,722,688/4,734,976; 4,950,016/4,728,832 respectively.
Small-instance peak memory therefore also regresses. All spreads and phase
samples are in `summary.json` / `samples.json`, including unfavorable outcomes.

The separate 257-leaf analytic control passes with 257 stored incidence terms
and 257 projection entries. Dense materialization and the original dense Schur
constructor reject the unchanged cap; rejection costs are not solve timings.
These measurements support the sparse incidence representation for this action
route only. They do not establish global speedups or comparable GUDHI/PETLS
timings, and the default dense complete-spectrum route is unchanged.

### Replay and evidence availability

Use the example's `--operator-benchmark N ACTIONS` at the two measured commits.
The existing resource controller accepts `kind: persistent-action`, `n`,
`actions`, the exact `command`, explicit fixture metadata and
`expected_checksum = ACTIONS*sum((x_i-mean(x))^2)`. Generate n=64/128/256,
ACTIONS=1/1000, both routes, using the formulas above; invoke the controller
with ten samples and a fresh output directory. Preserve the source/binary SHA
for each route in the specification.

Local attempt:
`target/benchmarks/commit-b1257aff7e47/persistent-action/run-001`.
It contains the full specification/environment, raw stdout/stderr, samples,
summary, paired analysis and hash manifest. A preliminary invocation with an
incorrect specification path executed no workers and produced no measurements;
it is excluded from this result. Evidence remains local only, with no durable
public archive. The CI construction/oracle artifacts are independently rerunnable
and retained for 14 days; they are not copies of these local timing samples.
