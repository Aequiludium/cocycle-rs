# Exact F2 H2 foundation

[Benchmarks](README.md) / [Phase-3 tasks](https://github.com/Aequiludium/cocycle-rs/issues/33)

## Frozen experiment contract

T0 freezes measured R0 at
`a1a7cb629187fbc17c9e2e638ee8ae4bf19c2656`, based on main
`9e6715f4c2e0118ab738486bd747ec3e361397de`. Counter additions compile only under
`cfg(test)`; ordinary release measurements retain the generic algorithm. The harness
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

The R0 quick pipeline run passed all 23 cases, 69 worker schedules and 897 fresh
processes (one warmup plus 12 measured samples each), with no exclusions or
failed cases. Source fingerprint:
`6b9251ef5d3b85d67dc8b7aec359748aa9f447432aeb33a3419e28b72cb5a52e`.
Kernel and harness are both the R0 SHA above. Linux Rust 1.91.0, native compiler
commands/header/binary hashes and machine/CPU identity are in environment.json.
This establishes scoped baseline measurements and repeated output agreement,
not a speed ranking. Local-only evidence root:
`target/benchmarks/commit-a1a7cb629187/phase3-baseline/run-001/`.

Dory attempt 001 failed compilation because the numeric-capture copy initially
omitted the pinned sort.h include path; its build.log is retained. The corrected
adapter adds that source include directory and records header hashes before a
fresh attempt. Generated sources, manifests, fixtures, commands and samples
remain in ignored target directories with distinct identities. No durable public
artifact store is configured; local evidence is not a public archive.

The n=32 R0 run also passed all 23 cases/897 processes under the same source
fingerprint. Its harness is `f04e6ba6340e64635b41180b5bc41d32e635ee52` and measured
kernel is R0; evidence is phase3-baseline/run-002. Counter run-001 exercised all
23 n=32 fixtures in separate Windows Rust 1.98.1 test processes. For high_fill32,
generic H1 processed 353/384 edges; H2 cleared 353/2287 triangles, processed
1934, made 408 additions and visited 56947 candidate vertices. Failed candidates
are included; these are work counts, not ordinary release timings.

Dory run-002 (harness f04e6ba above) validated five n=16 geometric families,
13 captures each. The equal clique exited 139 and the octahedron omitted [1,2)
in H2; both failures remain in results.json and invalidate that run's overall
summary. Source inspection finds strict edge < threshold selection and forced
matrix cone truncation. Those tied cases are outside the validated Dory domain,
not retroactively passing exclusions. Run-003 uses harness
`393d9699b0514d0f1f5494cc432cd650e2122d1d`, explicitly selects
uniform32/sphere32/noisy_sphere32, and passes all 39 captures. Both sphere cases
contain a positive H2 interval, so this is not solely empty-H2 validation.
Dory process timers have coarse 0.01-second resolution on these small inputs;
their raw values/RSS establish resource snapshots only. No ranking is reported.
