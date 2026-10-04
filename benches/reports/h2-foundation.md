# Exact F2 H2 foundation

[Benchmarks](README.md) / [Phase-3 tasks](https://github.com/Aequiludium/cocycle-rs/issues/33)

## T3 fixed-tuple baseline decision

The T3 prototype is correct on the exercised domain and independently profileable,
but remains test/benchmark-only. All 26 H2 rows improve against the frozen T2
implementation. One H1 control fails the preregistered median gate: duplicates32
is 1.053131 times T2, exceeding 1.05. This is a no-admission outcome for default
dispatch, not an optimized production-engine delivery. Retain this matched
baseline for M1; do not enable H2 shortcuts or capacity reuse in it.

The sections below this T3 account retain their T0/T2 identities and outcomes;
they are not fresh T3 measurements.

### Identity, validation and sampling

| Item | T3 evidence |
| --- | --- |
| Associated change | [T3 #37](https://github.com/Aequiludium/cocycle-rs/issues/37), stacked on foundation PR #50; PR link follows submission |
| Measured candidate/kernel and harness | `d6042f1895257d99258c61996eb0999f5be0792a`, clean before/after all formal measurements |
| Matched T2 kernel baseline | `f22d6207cf63a7d6552b0989a48444e3ce8cd05f`; stacked base `1bfd739dcf056ecc8f06254e0bbfc878df550bba` differs only in reporting |
| Candidate pipeline fingerprint | `a385958df0315d6c20c7b09965796cd6978b899803c5e2109f6e70992132c5c0` |
| Candidate pipeline executable SHA-256 | `cdce57c7eff4ede13609e277e1afce0f3a316458858f1aca3ce40b3de8cdc1b0`, identical for both sizes |
| Candidate build | Rust 1.91.0, Linux release, `RUSTFLAGS="--cfg cocycle_h2_bench"`; no test counters |
| Paired controller SHA-256 | `b47a008d390a81e9c6b342a30f1f811092985b9f890765511223f762888abcdf` at the measured commit |
| Machine | Intel Core Ultra 7 155H, Ubuntu/WSL2, CPU0; frequency and other load uncontrolled |
| Supported native validation | 588 cases, 1108 GUDHI/Ripser comparisons, 39 protocol checks passed; 68 explicit exclusions |
| Exact-reference fingerprint | `16b0d9741ad2e2fc3cec8af07efc92138d90c4d861a8ed66ee5bf571370fd14d` |
| Pipeline validation snapshots | `phase3-candidate/run-001` n16 and `run-002` n32: each 23 cases, 69 worker rows, 12 samples plus warmup; no exclusions/failures |
| Paired study | `localfs-paired/run-001` n16 and `run-002` n32: each 23 cases, 96 balanced paired rounds plus warmups, 4462 serial fresh processes |
| Paired execution window | 2026-10-04 11:24-11:26 UTC; exact start/end in each environment/summary |
| Memory | Per-process VmHWM, including input/runtime/allocator/results; live allocated bytes and allocation churn unavailable |

The comparison uses the previously frozen common float32 fixtures and unchanged
pipeline v2 public workflow timer, F2, H0/H1/H2, diagrams only, exact dense/supplied
flags, inclusive cutoffs and original multiplicity. Both frozen executables are
hash-verified copies on Linux localfs `/tmp`; fixture I/O stays outside the timer.
No compilation or other measurement runs concurrently. Native timings from the
DrvFS pipeline snapshots are retained separately and are not ranked against the
localfs paired rows. GUDHI/Ripser pins and capability exclusions are unchanged;
Dory is not rerun for T3.

Independent explicit vertex-subset boundary reduction and the full generic
engine agree for 1760 seeded dense/sparse cases, including empty/isolated inputs,
ties, zero distances, cutoff and finite/essential/censored endpoints. Additional
two-octahedron cases check repeated intervals, and adjacent/extreme binary64
values preserve endpoints. Independent all-edge enumeration checks every stored
test transformation R = C V, diagonal/triangular order and pivot progress.
Parity multiplicity and interruption at every counted boundary pass. These are
finite checks; fallible reservations preserve the allocation contract without
claiming that every possible allocator failure is recoverable.

### Matched timing and memory

Ratios are T3/T2. Intervals are paired bootstrap 95% intervals for the ratio of
medians (2000 replicates, seed 20261004). Times are milliseconds; RSS columns are
maximum process KiB, T2/T3. The octahedron is n6 in both independent size suites.
The table includes every H2 row, not just selected benefits.

| Case (dense unless flagged) | T2 ms | T3 ms | Ratio | 95% interval | Peak RSS T2/T3 |
| --- | --- | --- | --- | --- | --- |
| uniform16 | 0.337710 | 0.233879 | 0.6925 | 0.6786-0.7021 | 2816/2816 |
| circle16 | 1.464163 | 0.711613 | 0.4860 | 0.4729-0.4986 | 2992/2816 |
| clusters16 | 0.369084 | 0.237461 | 0.6434 | 0.6276-0.6579 | 2816/2816 |
| duplicates16 | 0.972725 | 0.506448 | 0.5206 | 0.5029-0.5352 | 2816/2816 |
| sphere16 | 0.738232 | 0.449744 | 0.6092 | 0.5937-0.6245 | 2816/2816 |
| noisy sphere16 | 0.709135 | 0.418070 | 0.5895 | 0.5704-0.6074 | 2816/2816 |
| nonmetric ties16 | 1.397150 | 0.672151 | 0.4811 | 0.4631-0.5013 | 2992/2816 |
| equal clique16 | 1.097590 | 0.488787 | 0.4453 | 0.4338-0.4549 | 2992/2816 |
| low degree16 (flag) | 0.052708 | 0.045363 | 0.8606 | 0.8254-0.8788 | 2640/2640 |
| geometric16 (flag) | 0.132584 | 0.104680 | 0.7895 | 0.7631-0.8062 | 2816/2816 |
| clique heavy16 (flag) | 0.184144 | 0.123750 | 0.6720 | 0.6478-0.6919 | 2816/2816 |
| high fill16 (flag) | 0.450611 | 0.310485 | 0.6890 | 0.6498-0.7181 | 2816/2816 |
| octahedron, n16 suite | 0.059978 | 0.045745 | 0.7627 | 0.7373-0.7860 | 2640/2640 |
| uniform32 | 5.145765 | 2.723486 | 0.5293 | 0.5112-0.5479 | 3344/2992 |
| circle32 | 29.604637 | 8.680419 | 0.2932 | 0.2869-0.3036 | 4576/3696 |
| clusters32 | 6.559097 | 3.000751 | 0.4575 | 0.4495-0.4685 | 3520/3168 |
| duplicates32 | 14.187084 | 4.912263 | 0.3462 | 0.3351-0.3572 | 3872/3344 |
| sphere32 | 14.202564 | 5.523109 | 0.3889 | 0.3805-0.3967 | 4048/3344 |
| noisy sphere32 | 16.193848 | 6.653002 | 0.4108 | 0.4006-0.4204 | 4048/3344 |
| nonmetric ties32 | 38.785048 | 10.572293 | 0.2726 | 0.2633-0.2812 | 5280/4224 |
| equal clique32 | 25.710567 | 6.793912 | 0.2642 | 0.2575-0.2704 | 4576/3696 |
| low degree32 (flag) | 0.085086 | 0.070702 | 0.8309 | 0.7922-0.8761 | 2816/2640 |
| geometric32 (flag) | 1.825708 | 1.015819 | 0.5564 | 0.5381-0.5908 | 2992/2816 |
| clique heavy32 (flag) | 0.321323 | 0.222158 | 0.6914 | 0.6717-0.7114 | 2816/2816 |
| high fill32 (flag) | 8.183316 | 3.241343 | 0.3961 | 0.3858-0.4097 | 3696/3344 |
| octahedron, n32 suite | 0.057434 | 0.044940 | 0.7825 | 0.7573-0.8140 | 2640/2640 |

All H2 upper bounds are below one across geometric, tied and exact sparse
families, without a >10% H2 regression. H1 n16 passes 10/10 controls (maximum
median 1.022324, maximum upper95 1.044573). H1 n32 passes 9/10: duplicates32
increases from 0.4056105 to 0.427161 ms, ratio 1.053131, 95% interval
1.029576-1.079797. Its upper95 passes 1.10, but median fails 1.05. The n32
paired summary remains `failed`; do not relabel it passing or remove that row.
No code-layout or load explanation has been established. No additional sampling
or candidate selection was performed after this outcome.

### Work and ownership baseline for M1

Separate Windows Rust 1.98.1 release-test processes validate all 26 emitted raw
diagrams against the frozen pipeline results. Counters are not latency evidence.
For high_fill32, T2 and T3 both enumerate 2287 triangles, clear 353, process 1934,
perform 2342 pivot lookups and 408 owner additions, and store 2462 transformation
entries with largest column 75. Both emit 33359 tetrahedra including reconstruction.
This supports the same algebraic work with a different representation, not pruning.

The generic path visits 56947 shortest-list candidates and makes 138127 edge
queries; the prototype performs 77228 three-way intersection rounds and 231684
list comparisons. These count different operations; neither counts just emitted
tetrahedra, and these figures do not establish fewer comparisons. Tuple/parity
storage avoids generic vertex-vector cloning and coefficient arithmetic. The
prototype's heap peaks at 768 entries versus the generic normalized column's
522; these are different storage units and cannot be read as live-byte ratios.
Observed RSS is equal or lower on the measured rows, but allocator/live-byte
tradeoffs remain for M1.

One high_fill32 instrumented process records H1 return with owned death tuples
at 249 microseconds (353 entries, Vec capacity 512), extracted clearing at 260
(353 entries, HashSet capacity 448), triangle construction after regenerated
edges drop at 606 (2287 entries, Vec capacity 4096), and H2 start at 606 (no stored
columns yet). Extraction into tuples occurs during H1 pairing. H1 edge positions,
cycle flags, pivot owners, transformations and heaps have left scope at return;
shared graph/access/index and owned intervals survive. The clearing set stays
alive and per-column scratch is fresh in this baseline. These object/event
observations neither establish total live allocation nor an RSS release bug.

### Checks, reproduction and retained evidence

Normal and private-config complete Rust tests/Clippy pass, as do debug/release
tests, Rust 1.91 tests/all-target checks, nine examples, construction tests,
strict rustdoc, README/guide doctests, focused persistence-reduction checks,
source/docs/artifact checks and 88 Python tool tests. Local collaboration source
checks pass; its copied historical ignored artifact links remain unavailable.
Hosted CI is recorded after submission, separately from these local checks.

Reproduce the committed prototype with `RUSTFLAGS="--cfg cocycle_h2_bench"` and
the [pipeline commands](../pipeline/README.md#private-h2-baseline), both sizes,
`--samples 12 --cpu 0 --phase3`, and the pinned native sources. Run the full
`tools/compare_rips.py` command under the same configuration. Then unset RUSTFLAGS
and compare the validated T2/T3 run directories using `tools/profile_rips.py`,
96 pairs and CPU0. Counter processes use the separate ignored test entry there.

Evidence is local-only under `target/benchmarks/commit-d6042f189525/` in the T3
checkout: `exact-reference/run-001`, `phase3-candidate/run-001/002`,
`localfs-paired/run-001/002`, and `counters/run-002`. Original T2 builds/results
remain under their own measured identity in the foundation checkout. Environments,
worker/controller hashes, exact fixtures, every raw sample, native exclusions,
summaries, build logs and the counter executable are retained. The first counter
controller attempt (`counters/run-001`) failed sorting a Python None against a
float after a successful kernel process; its failed artifact is retained and the
canonical-key retry uses the same kernel. No public artifact upload or retention
guarantee is claimed. Paired results SHA-256:

- n16: `882cc9fbf426532010e8472b9cd8dd354a3b4d84d6ebde1161f2fe2feca01f52`.
- n32: `6476d752f67cc70b4230202648632612c7210d1ece84e9455fd57e569f7e391f`.

## Frozen experiment contract

T0 freezes measured R0 at
`a1a7cb629187fbc17c9e2e638ee8ae4bf19c2656`, based on main
`9e6715f4c2e0118ab738486bd747ec3e361397de`. Counter additions compile only under
`cfg(test)`; ordinary release measurements retain the generic algorithm. The harness
revision and fingerprints are recorded separately by each controller. Associated
PR: [Phase-3 foundation #50](https://github.com/Aequiludium/cocycle-rs/pull/50).

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

## T2 candidate and admission

Measured kernel and pipeline harness:
`f22d6207cf63a7d6552b0989a48444e3ce8cd05f`; baseline is the R0 SHA above.
Pipeline source/harness fingerprint:
`ae6af5b6671ff198a6b140f392447f7bcfef2e258db37d75556b4a851c7209d3`.
This is the implementation revision, distinct from the subsequent report commit.
Both candidate pipeline runs pass 23 cases/897 processes, with no exclusions.
The fresh native correctness run passes 588 cases, 1108 comparisons and 39
protocol checks, retaining 68 explicit capability exclusions. Its independently
defined source fingerprint is
`de760580fd0dd06ae8d1e5ecd3ba02e09518b0bf7f79c7eebff15c88245f09ad`.

The first candidate, `429fe9d0d31e34c609e109c9b75ba5ffdda5d3a6`, is rejected by
the original WSL DrvFS H1 timing gate: 10/10 n=32 controls fail in 96 balanced
pairs, with median ratios 1.169-1.495 and a maximum upper bound 1.545. Keep
`commit-429fe9d0d31e/h1-paired/run-001/`, rather than replacing its result.
Fresh-process executable demand paging affects these tiny timings. A diagnostic
with byte-identical baseline/candidate workers copied to Linux local storage
isolates that effect. At the final kernel, six uniform32 process captures per
placement show 17-19 major I/O faults for R0 and 20 for T2 on DrvFS, versus zero
for both on local storage. The diagnostic remains separate from admission.

Before the independent local-storage repeat, fix executable placement to Linux
local storage for both kernels and verify their SHA-256 hashes against the
original frozen workers. Keep pipeline v2 elapsed_ms, public API/output scope,
fixtures, one warmup, fresh processes, CPU 0, balanced position order and the
original regression limits. This changes the storage environment, not the
timing boundary or requested work. At revision 429fe9d, 96 pairs leave three
n=16 H1 controls marginally outside the median gate; the preregistered 384-pair
precision extension passes n=16 but leaves low_degree32 at 1.0511. Those failed
summaries remain in localfs-paired/run-001 through run-004 under that revision.
The final kernel adds only an inline hint at the existing H1-only entry.

Final independent local-storage H1 repeats use 384 balanced pairs per control
and 2000 paired bootstrap replicates. All 20 controls pass: at n=16 the largest
median ratio is 1.045794 and largest 95% upper bound 1.064360; at n=32 they are
1.015417 and 1.021304. This is 7700 fresh processes per size. Final H2 repeats
use 96 balanced pairs per case (2522 processes per size); every multiset agrees
and all 13 cases per size show a benefit with the upper 95% bound below one.
Ratios below one favor T2; the octahedron remains six vertices in both repeats.

| H2 family | n=16 ratio [95% paired interval] | n=32 ratio [95% paired interval] | n=32 peak RSS KiB (R0 / T2) |
| --- | --- | --- | --- |
| uniform | 0.798 [0.782, 0.814] | 0.827 [0.815, 0.834] | 3520 / 3344 |
| circle | 0.811 [0.799, 0.825] | 0.908 [0.896, 0.922] | 4752 / 4576 |
| clusters | 0.794 [0.785, 0.811] | 0.826 [0.814, 0.837] | 3696 / 3520 |
| duplicates | 0.818 [0.793, 0.835] | 0.871 [0.853, 0.892] | 3872 / 3872 |
| sphere | 0.812 [0.800, 0.831] | 0.859 [0.846, 0.873] | 4224 / 4048 |
| noisy_sphere | 0.805 [0.793, 0.818] | 0.844 [0.836, 0.858] | 4048 / 4048 |
| nonmetric_ties | 0.761 [0.747, 0.775] | 0.907 [0.887, 0.917] | 5280 / 5280 |
| equal_clique | 0.807 [0.799, 0.828] | 0.901 [0.879, 0.916] | 4752 / 4576 |
| low_degree | 0.922 [0.865, 0.950] | 0.823 [0.805, 0.860] | 2816 / 2816 |
| geometric | 0.802 [0.782, 0.819] | 0.785 [0.772, 0.795] | 2992 / 2992 |
| clique_heavy | 0.818 [0.805, 0.835] | 0.785 [0.775, 0.789] | 2816 / 2816 |
| high_fill | 0.717 [0.702, 0.738] | 0.818 [0.806, 0.829] | 3872 / 3696 |
| octahedron | 0.932 [0.899, 0.967] | 0.896 [0.857, 0.931] | 2640 / 2640 |

RSS columns are maximum observed process VmHWM, not live allocation measurements
or a proof of a dimension-peak memory bound. No case shows an unexplained >10%
H2 regression in this local-storage repeat. This supports the bounded T2
continuation on these Linux workloads; it is neither a universal speedup nor a
passing DrvFS deployment claim. CPU frequency/load remain uncontrolled; paired
sampling reduces drift without proving its absence. No cross-library ranking
or new fast H2/H3 kernel is admitted by these measurements.

## Work changes, validation and retained evidence

At the final kernel, 23 n=32 fixtures each run in separate dispatch and full
generic Windows Rust 1.98.1 counter processes (46 total). Every H2 Stats row is
identical between routes. For high_fill32, both retain 2287 triangles, clear 353,
process 1934 and perform 408 column additions; the dispatch route has zero
generic H1 simplices/columns instead of 384/353. Generic triangle candidate
visits fall from 20104 to 8804, including rejected candidates; specialized H1
visits belong to its separate Stats and are not zero work. Direct dispatch tests
also prove zero generic H0/H1 processing. Clearing retains all triangle topology.
The current counter executable is copied into its immutable run directory;
R0 retained its binary hash/logs but did not separately preserve that executable.

Independent complete-death-index tests exercise ordinary/emergent/omitted
virtual pairs under three shortcut configurations. The equal K3 incomplete-key
counterexample, exhaustive/seeded boundary oracle, every-budget octahedron retry,
cutoff/coverage, odd-prime, representative and approximation contracts pass.
Windows checks pass formatting, all-target Clippy with warnings denied, complete
debug/release tests, all nine examples, construction example tests, strict rustdoc
and focused persistence-reduction checks. Linux Rust 1.91.0 passes complete tests,
all-target checks, construction example tests and README/nine guide doctests.
All 85 Python tool tests and tracked source/documentation/artifact checks pass.
The inline follow-up repeats format/Clippy/debug/release/rustdoc and MSRV checks;
the unchanged examples/tool/guide checks retain their 429fe9d identity. Hosted
CI and merging are separate submission evidence. The copied local collaboration
checker passes tracked-source checks but reports missing historical local links;
those unrelated local records are not part of this contribution.

Final evidence root is `target/benchmarks/commit-f22d6207cf63/`: exact-reference
run-001; phase3-candidate run-001/run-002; localfs-paired run-001/run-002 for H1
and run-003/run-004 for H2; counters run-001; filesystem-audit run-001. Each paired
directory retains the exact controller source/hash, original build environments,
byte-identical executable hashes, fixture hashes, schedule, raw samples and
bootstrap results. The pipeline source/harness revision is f22d6207cf63; the
archived supplemental paired controller has its own SHA-256 in environment.json.
All generated evidence remains local and ignored; there is no durable public
artifact archive. Reviewers can reproduce pipeline fixtures using --phase3
(--quick for n=16), with the pinned native setup and sampling contract above.
Allocator/live-storage ablations remain M1/M2 work.
