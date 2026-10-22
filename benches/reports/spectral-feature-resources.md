# Persistent spectral workflow resources

Evidence class: **resource snapshot** under the [performance reporting rules](../reporting.md).

PR: [#129](https://github.com/Aequiludium/cocycle-rs/pull/129).
Measured candidate: `af53741fc9a16208dcc14711c1b4e8574b634186` (feature-branch commit, not a tested merge).
Baseline: not applicable; no previous-head speedup or external native ranking is claimed.
The report commit is separate from the measured implementation/controller commit.

Source-manifest SHA-256: `7d40965c2f85cd3a39a62a08a259516b2133468232bcd900dd7b442bb8fbd0c2`.
Worker-manifest SHA-256: `cd32358b38fc261f25670125f2790ce2520e861b676fa68e2cbb38fb53b392dc`.
These are SHA-256 of the recorded path-to-file-hash maps encoded as UTF-8 JSON
with sorted keys and compact separators; paths retain their recorded backslashes.
Controller SHA-256: `d493b2b3f01de9b8af12958b10ba12ea29e026256bbf4c0a3c58d7f1ab9b7b71`.
Specification SHA-256: `ea1c3c2b29aa0a8048e5c03f6eb44d73fb1c72341eeb2d88f8d603b1ec830e2d`.
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
| full-spectral-feature-reproduction | passed: 10 | not returned | 47527.906 (40511.351-72524.260) | 264.633 (264.273-265.379) |

The cell is the complete existing synthetic reproduction: 160 worked pair
matrix comparisons, the pinned author feature-loop audit, real Rust descriptor
workers, and two group-held-out random-forest experiments. Each experiment has
64 matched pairs, 96 training and 32 test rows, seed 20261012, and n_jobs=1.
The tree/disk controls and their shared motifs do not establish dataset-level
generalization. Runtime is a whole-workflow resource observation, not a spectral
kernel benchmark or a cross-library ranking.

The original reproduction enables Python tracemalloc; its overhead is included.
This is a traced workflow observation rather than an uninstrumented kernel cost.

Each fresh process creates fresh Rust child build directories and includes those
builds, input generation, independent numerical/topology checks, vectorization,
training, predictions and output retention. Whole-process wall time includes
startup and final serialization. Windows peak working set covers the parent
Python process only, not the simultaneous/aggregate Rust child peaks. It is
therefore not whole-process-tree memory, allocation peak or a hard memory quota.
The initial process is a retained discarded warmup, and every measured outcome
is retained. Input seed/counts and numerical dependency versions are fixed.

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

The actual local attempt is `target/benchmarks/commit-af53741fc9a1/spectral-resources/run-001`.
It retains the exact specification/fixture hash, immutable source/worker
fingerprints, binary hashes, environment, raw stdout/stderr, every warmup and
sample, complete sample matrix, summaries and file manifest. Companion plan metadata
records dependency versions, CPU, UTC start/end and generator hash.
Raw evidence is **local only**; no durable public archive is configured.
Committed source/settings allow replay, but the original samples are not
claimed publicly retrievable. CI correctness artifacts are a separate,
14-day evidence stream and are not performance baselines.
