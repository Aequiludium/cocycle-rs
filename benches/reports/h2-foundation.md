# Exact F2 H2 foundation

[Benchmarks](README.md) / [Phase-3 tasks](https://github.com/Aequiludium/cocycle-rs/issues/33)

## Frozen experiment contract

T0 fixes the ordinary generic kernel at
`9e6715f4c2e0118ab738486bd747ec3e361397de`. Counter additions compile only under
`cfg(test)`; ordinary release measurements retain this kernel. The harness
revision and fingerprints are recorded separately by each controller. Associated
PR: the Phase-3 foundation PR, supplied at submission; no earlier PR is replaced.

`benchmark_rips_pipeline.py --phase3` uses the existing
[pipeline v2](../pipeline/README.md). Exact F2, H0/H1/H2, owned diagrams, inclusive
cutoffs, zero-length bars omitted, and one common float32 quantization apply.
Seed 20261004 fixes uniform Euclidean, circle, sphere/noisy sphere, clusters,
duplicates, nonmetric ties, equal clique, low-degree, geometric, clique-heavy
and high-fill supplied flag graphs. Quick uses n=16; ordinary uses n=32.
H1 controls share each geometric/graph input; size-growth/shrink experiments
use both sizes. Supplied sparse graphs are exact, without approximation blockers.

Preregistered sampling: one warmup, 12 balanced serial measured rounds, fixed
allowed CPU, fresh process per sample. Preserve medians, ranges and raw samples.
Baseline evidence is a resource snapshot, not a stable ranking. An optimization
needs an independent repeat, benefit in two families, no unexplained >10% H2
regression and a <=5% H1 median regression with a paired bootstrap 95% upper
bound <=10%; insufficient precision requires more preregistered samples.
Memory is process VmHWM, including input/runtime/allocator/results, not live
reducer bytes. No peak-memory admission limit is inferred from RSS alone.

Native GUDHI/Ripser source pins remain in [sources.json](../native/sources.json).
Dory is pinned there with its MIT license. `benchmark_dory.py` uses one thread,
F2, square CSV, cycles disabled and upstream SAVEPD with 17-digit numeric capture.
Its process timer includes startup, input and file output; compare correctness,
but do not rank that timer against pipeline elapsed_ms. Graphs/isolates, zero
distances and censored endpoints are outside this initial Dory adapter's audited
scope and remain pending adapter work, rather than claimed passing comparisons.

## Independent correctness and counter schema

The existing bitmask/forward-boundary oracle in
[rips_expansion.rs](../../tests/rips_expansion.rs) shares no production cofacet,
clearing or elimination code. It covers all {0,1,2,missing} edge assignments for
n<=4 at cutoffs 0/1/2 and 128 seeded assignments each at n=5/6. This is finite
coverage, not all 4^15 assignments at n=6. Complete supplied-graph versus
censored dense coverage is checked separately from endpoint tuples.

The ignored `simplicial::cohomology::profiling::profile_h2` test reads
COCYCLE_H2_FIXTURE. Per dimension it records processed/cleared/simplex counts,
pivot lookups/hits, column additions, reconstruction calls, cumulative stored
transform entries, largest transform and peak working-column entries.
Clique counters record every visited candidate, including rejects, edge queries
and emitted cofacets; sparse generic access scans the shortest adjacency list,
so intersection counts are inapplicable there. Apparent/emergent/virtual counters
are zero in the generic baseline; specialized H1 has its separate existing Stats.
Counters and latency builds are separate, including debug library builds.

Dimension-boundary ablation schema: joint baseline, lifetime-only, same-type
scratch, later typed generic reuse. Record shared input/result bytes, clear-next
keys and conversion scratch, old/new topology overlap, last H1 reducer use/drop,
retained nested capacity, reset cost, allocation/runtime cost, RSS samples and
process peak. Live allocator bytes and allocation/reallocation/free counts are
unavailable without an allocator profiler; use null, never RSS as their estimate.
Failure/cancellation recovery and cutoff tests remain required controls. T2 owns
complete pairing handoff; M1/M2 own later lifetime/reuse ablations.

## Results

Measurement results are added after the implementation and harness commit is
frozen. Generated sources, manifests, fixtures, commands and samples remain in
ignored target directories with distinct baseline/candidate identities.
