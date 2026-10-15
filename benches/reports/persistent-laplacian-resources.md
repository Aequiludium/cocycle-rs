# Persistent Laplacian reference resources

Evidence class: **resource snapshot** under the [performance reporting rules](../reporting.md).

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
