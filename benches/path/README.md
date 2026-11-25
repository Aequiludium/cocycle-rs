# Matched directed-path H0/H1 comparison

This is an in-crate algorithm comparison, not a cross-library ranking. Both
native routes receive identical ordered weighted edges and compute only the
standard PersistenceResult for complete ordinary H0/H1, including context and
interval multiplicity. The generic route constructs compatible Omega bases
(with its current private endpoint grouping); the direct route constructs the
filtered one-boundary presentation. Both use the existing ordinary reducer.
Neither requests homology representatives or returns an internal certificate.

Commit implementation and harness first, then run from a clean source tree:

```sh
python tools/benchmark_path_topology.py --output target/benchmarks/commit-SHA/path-topology/run-001
```

The versioned protocol is `cocycle-path-diagram-v1`. A single release native
worker is compiled before sampling. Each sample uses a fresh process, and workers
run serially without concurrent compilation/tests. The internal timer excludes
fixture parsing and field validation, and includes public-terminal validation,
preparation, computation, canonical result assembly, equal interval payload
export and result/workspace destruction. JSON formatting/process startup are
excluded. Wall times are recorded separately and are not used for rankings.

The declared grid uses two seeds, cycles/sparse random graphs at 8/32/128/512
vertices, dense random graphs at 8/16/32/64 vertices, and three-layer multisquare
graphs at 9/18/36/72 vertices. Every graph has a static zero-birth and a filtered
three-event variant. Default algebra is regular, non-reduced F2 through H1,
requiring allowed paths/Omega through degree two. Optional field/model arguments
produce separate studies. Inputs, source/worker hashes, exact commits and limits
are retained in each fresh output directory.

Each cell has one discarded warmup per route and at least ten measured processes
per route. Backend order is position-balanced by round; cell order is seeded.
Every successful sample must return exactly the same complete interval multiset.
All raw samples, unsuccessful outcomes and warmups are retained. A failed cell is
not ranked. Limits are a 15-second process timeout, 250000 paths, 2000000 sparse
construction entries and 250000000 cooperative work units. Sparse/work caps are
not aggregate memory quotas; no OS memory cap or process-RSS claim is made.

Report median and min/max per route. The ratio is generic/direct, greater than
one favoring direct. The paired 95% bootstrap percentile interval uses 2000
resamples of matched round indices, selecting order statistics at floor(p*(N-1)).
`direct_advantage` requires a median ratio at least 1.1 and lower bound above one;
the converse regression rule requires ratio at most 0.9 and upper bound below
one. This describes the fixed observed cells, without a universal crossover or
asymptotic claim. CPU affinity, host load and frequency remain uncontrolled.

`--quick --samples 1` exercises only the smallest graph in each family and one
seed. Its timings are smoke evidence and must not select an algorithm. Generated
fixtures, binaries, logs and raw data stay outside source Git under `target/`.
