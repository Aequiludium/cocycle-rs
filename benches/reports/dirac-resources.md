# Dirac resources and conditioning

Evidence class: **resource snapshot** under the [performance reporting rules](../reporting.md).

PR: [#133](https://github.com/Aequiludium/cocycle-rs/pull/133).
Measured candidate: `d5d61cb103c4c0310cb97eaad6fda926b2d66dda` (feature-branch commit, not a tested merge).
Baseline: not applicable; no previous-head speedup or external native ranking is claimed.
The report commit is separate from the measured implementation/controller commit.

Source-manifest SHA-256: `a00d2a9964201c8a5ab44dd77e3ee379ed96ef521bd2cf2611c9ae391ed562ce`.
Worker-manifest SHA-256: `35d736509daa026680ebc58088f6abe1725e53593d70062f2f4a6db278f9ab28`.
These are SHA-256 of the recorded path-to-file-hash maps encoded as UTF-8 JSON
with sorted keys and compact separators; paths retain their recorded backslashes.
Controller SHA-256: `d493b2b3f01de9b8af12958b10ba12ea29e026256bbf4c0a3c58d7f1ab9b7b71`.
Specification SHA-256: `0817dfc4fed74fb778ec5152d52846aff04d699d1060d0506ebd591536d618ae`.
The controller and workers belong to the measured commit above.

Host: Intel(R) Core(TM) Ultra 7 155H, Windows-11-10.0.26100-SP0; Rust 1.98.1, Python 3.12.14.
Oracle/workflow packages: NumPy 2.2.6, SciPy 1.15.3, GUDHI 3.12.0;
learning additionally uses scikit-learn 1.7.2 and threadpoolctl 3.6.0.
OMP/OpenBLAS/MKL thread environment is one. CPU affinity, frequency and
external host load were uncontrolled. Workers ran serially without concurrent
local compilation/tests; reproduction-internal builds are part of that cell.
Recorded overall plan UTC: `2026-10-15T03:43:16.577385+00:00` to `2026-10-15T03:55:56.865792+00:00`.
Individual native-cell UTC endpoints were not separately recorded.

Protocol: `spectral-resource-v1`; one retained discarded warmup and ten fresh
measured processes per cell, seeded order 73, no excluded outliers/retries.
Native workers use `rustc --edition 2024 -C opt-level=3` against the locked
release library. Reproduction uses its unchanged internal debug Rust builds.
The table gives median (minimum-maximum); times are milliseconds and memory
is MiB. Peak working set includes runtime/input/allocator/retained output,
rather than live algorithm allocation. Numerical/storage failures are retained.

| Cell | Outcomes / 10 | Solve ms | Whole process ms | Peak working set MiB |
| --- | --- | --- | --- | --- |
| path-128-count3-block8 | success: 10 | 47.315 (45.292-50.462) | 58.570 (56.671-62.526) | 4.568 (4.527-4.582) |
| path-1024-count3-block8 | success: 10 | 2248.431 (2126.130-2401.826) | 2260.817 (2138.572-2413.307) | 5.090 (5.074-5.355) |
| star-2048-count3-block8 | success: 10 | 9.860 (9.112-10.064) | 21.440 (20.911-24.613) | 5.631 (5.609-5.648) |
| jordan-16-count16-block16 | success: 10 | 0.331 (0.297-0.426) | 10.552 (9.011-14.137) | 4.484 (4.461-4.496) |
| jordan-32-count32-block32 | numerical_failure: 10 | 1.732 (1.668-2.487) | 12.887 (10.135-17.380) | 4.492 (4.480-4.547) |

The operator requests the three largest singular modes of B1 for paths/stars,
not the smallest-magnitude global Dirac eigenvalues. All modes are requested for
the two Jordan chains: B has unit diagonal and -2 on its first superdiagonal.
Their complete NumPy SVD is an independent index/value check; the 32-coordinate
case exercises the normal map's squared conditioning and fails numerically in
all ten processes. It is a supplied signed integer chain, not a weighted
simplicial filtration. The successful 16-coordinate case is not evidence of
reliable arbitrarily small singular values. Partial results always mark index
coverage uncertified.

Path fixtures have n vertices and n-1 oriented edges; the star has n vertices
and n-1 edges incident to vertex zero. The corresponding Dirac coordinate counts
are 2n-1; Jordan chains have 2n. Each selected boundary mode retains both signed
Dirac vectors. Analytic path singular values are 2*cos(pi*k/(2*n)), k=1,2,3;
the star requests sqrt(n),1,1. Tolerances for these independent checks are
atol=1e-8, rtol=1e-7; Jordan SVD checks use atol=1e-15, rtol=1e-4, without changing
the library's default residual tolerances. The Jordan failure is retained rather
than loosening the check.

Preparation measures signed-chain freezing and validation. Solve measures the
normal-map iteration, reconstruction and actual full-Dirac residual checks,
including retained result allocations but excluding destruction. It never
materializes the full Dirac matrix. The result is a finite resource snapshot;
it does not establish general sparse scaling or resolve normal-map conditioning.

The native workers use default `Execution` without cancellation/work caps;
the documented algorithm workspace, block and iteration policies still apply.
No process timeout or OS address-space cap was imposed. Peak memory is an
observation, never a hard library quota.

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

The actual local attempt is `target/benchmarks/commit-d5d61cb103c4/spectral-resources/run-001`.
It retains the exact specification/fixture hash, immutable source/worker
fingerprints, binary hashes, environment, raw stdout/stderr, every warmup and
sample, complete sample matrix, summaries and file manifest. Companion plan metadata
records dependency versions, CPU, UTC start/end and generator hash.
Raw evidence is **local only**; no durable public archive is configured.
Committed source/settings allow replay, but the original samples are not
claimed publicly retrievable. CI correctness artifacts are a separate,
14-day evidence stream and are not performance baselines.
