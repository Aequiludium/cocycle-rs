# Hodge resources and small-gap convergence

Evidence class: **resource snapshot** under the [performance reporting rules](../reporting.md).

PR: [#130](https://github.com/Aequiludium/cocycle-rs/pull/130).
Measured candidate: `373774dd68ee0cea26f69e19687ff148fd7153d8` (feature-branch commit, not a tested merge).
Baseline: not applicable; no previous-head speedup or external native ranking is claimed.
The report commit is separate from the measured implementation/controller commit.

Source-manifest SHA-256: `4433ffe5c8abce617d5189f831181f818490c4f79fb01d2b143cd640702df7af`.
Worker-manifest SHA-256: `11bffa484ed598119ca45b30992616c8a44e855c36aad38fbc6943a666e4cd17`.
These are SHA-256 of the recorded path-to-file-hash maps encoded as UTF-8 JSON
with sorted keys and compact separators; paths retain their recorded backslashes.
Controller SHA-256: `d493b2b3f01de9b8af12958b10ba12ea29e026256bbf4c0a3c58d7f1ab9b7b71`.
Specification SHA-256: `e848af7b17defae896a40eb7a8c511247316a9cf141943dcdd7d2307cc58028c`.
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
| path-128-block8-iterations800 | success: 10 | 56.836 (51.855-65.166) | 69.218 (62.709-77.433) | 4.551 (4.543-4.582) |
| path-1024-block8-iterations800 | numerical_failure: 10 | 2718.215 (2615.578-2844.422) | 2733.453 (2629.435-2857.830) | 5.582 (5.543-5.664) |
| path-1024-block16-iterations800 | success: 10 | 5695.079 (5276.547-5784.468) | 5708.232 (5289.215-5798.135) | 6.107 (6.094-6.152) |
| path-1024-block8-iterations2400 | success: 10 | 3520.614 (3282.024-3702.207) | 3535.640 (3296.467-3717.071) | 5.609 (5.594-5.668) |
| star-2048-block8-iterations800 | success: 10 | 9.815 (8.786-13.331) | 27.578 (24.310-31.697) | 6.662 (6.637-6.719) |

All cells request the three smallest approximate q=0 full Hodge pairs of
unweighted graphs. Paths have n vertices/n-1 edges and analytic complete spectrum
2-2*cos(pi*j/n), j=0,...,n-1. The star has spectrum 0, repeated 1, n. Returned
prefixes are checked against these independent expectations with atol=1e-8,
rtol=1e-7; the library's 1e-9 absolute/relative residual defaults are unchanged.

The 1024-node path fails all ten default block-8/800-update attempts. Increasing
the update budget to 2400 or the block to 16 restores convergence in this fixed
case, at the measured additional cost. These controls establish a real default
convergence limit; they do not prove a general solver defect or guarantee that
larger budgets always recover the correct extremal indexes. The private solver
is identity-preconditioned; a larger block is not a new preconditioner.

Preparation includes deterministic simplicial input construction and sparse
operator validation. Solve includes iterative starts, sparse actions, projected
Ritz solves, fresh residual checks and retained result construction, excluding
result destruction. Dense reference timings are not mixed with these rows.
Independent complete-spectrum corpus checks and the missed-invariant-subspace
regression remain separate correctness evidence. Successful results always
mark extremal index coverage uncertified.

## Independent convergence diagnosis

The complete path spectrum gives lambda1=9.4123808476e-6, lambda2=3.7649434798e-5
and a gap of 4.7061461274e-5 after the requested third root. Native block 8 needs
1053 updates/25,283 operator calls with the extended budget; block 16 needs
454 updates/21,811 calls but its larger Ritz work costs more in this snapshot.

With the same SplitMix64 seed-73 starting column span, independent SciPy 1.15.3
controls at block 8, 800 updates and tol=1e-9 give the following measured
first-three maximum residuals:

| Numerical control | Maximum residual | Requested three pass | History length |
| --- | --- | --- | --- |
| Identity | 7.878e-9 | no | 780 |
| Inverse diagonal | 4.160e-10 | yes | 799 |
| Inverse of L+1e-3 I, tridiagonal solve | 9.898e-11 | yes | 18 |

The diagonal control still warns that its complete eight-vector lookahead block
did not converge. The shifted inverse has no warning. Their selected indexes
agree with the complete analytic spectrum; these are diagnostic convergence
observations, not time comparisons or index certificates for other inputs.
The distinct SciPy iteration/postprocessing/stopping semantics prevent treating
its history length as native Rust iterations. Together the controls support a
small-gap/identity-preconditioning limitation on this case; they neither prove
a general algorithmic defect nor justify a new default without broader tests.
Replay with `tools/check_hodge.py --path-convergence` in addition to the existing
required fixture arguments. The structured inverse is a reference control, not
a newly implemented Rust preconditioner.

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

The actual local attempt is `target/benchmarks/commit-373774dd68ee/spectral-resources/run-001`.
It retains the exact specification/fixture hash, immutable source/worker
fingerprints, binary hashes, environment, raw stdout/stderr, every warmup and
sample, complete sample matrix, summaries and file manifest. Companion plan metadata
records dependency versions, CPU, UTC start/end and generator hash.
Raw evidence is **local only**; no durable public archive is configured.
Committed source/settings allow replay, but the original samples are not
claimed publicly retrievable. CI correctness artifacts are a separate,
14-day evidence stream and are not performance baselines.
