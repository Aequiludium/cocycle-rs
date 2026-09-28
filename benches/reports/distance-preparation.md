# Diagram-distance preparation and allocation effects

[Benchmark reports](README.md) / [Distance protocol](../distances/README.md)

## Question and status

This investigation accompanies [Issue 14](https://github.com/Aequiludium/cocycle-rs/issues/14).
The implementation removes the facade's intermediate finite/essential vectors
and prepares algorithm-owned storage from logical diagram views. Fewer copies
are not a performance guarantee. The change remains a draft pending independent
reproduction of a finite W1 regression and discussion of the appropriate fix.
No padding allocations or experimental selectors are added to the public crate.

The counterexample is a complete dimension-zero diagram with 512 finite
intervals against 64, evaluated with W1 and L-infinity ground distance.
Endpoints are f64 and multiplicity is retained. This is a before/after and
allocation-intervention experiment, not a cross-library speed ranking.

## Reproduce

Use Linux (including Ubuntu WSL), Git, Python 3.12 and Rust 1.91.0. No Python
packages or external reference libraries are required for this experiment.
Run from a clone with the baseline history available (`git fetch --unshallow`
if needed). Commit harness edits before collecting performance evidence.

```sh
python3 tools/reproduce_distance_preparation.py --mode all \
  --samples 14 --iterations 300 \
  --output "target/benchmarks/commit-$(git rev-parse --short=12 HEAD)/diagram-distance-preparation/run-001"
```

The default candidate is committed `HEAD`; `--candidate <full-sha>` selects
another revision (also use its short SHA in the artifact directory name).
The default baseline is `75330fcf7e5b9ca764641c68b987a71bbfad8498`.
`--mode public` runs ordinary baseline/candidate APIs. `--mode causal` runs the
seven interventions below in one executable. `--mode growth` uses separate
instrumentation. `--metric w2` provides a second metric control; `--size` changes
the 8:1 workload. Outputs must be fresh directories beneath ignored `target/`.

For a build/protocol smoke check use `--samples 2 --iterations 3`; such runs do
not support performance rankings. CI runs this smoke and uploads raw results,
the fixture, identities and build logs for 14 days. A CI runner's timing is not
substituted for local comparative measurements.

## Counterexample and measurement contract

The generator uses Python `random.Random(20260928 + 512)`. It generates the 512
left intervals first, then the 64 right intervals; each uses a birth uniformly
sampled from [-1, 1] and a lifetime from [1, 3]. The transport is `COCDST1\0`,
two little-endian u64 counts, followed by little-endian f64 birth/death pairs.

- Fixture SHA-256: `2d5b79549d03532405601cbc51d380f09f7b6e1f82cff0134f7230406643d99a`.
- Previously verified W1 distance: `442.10173773952255`.
- Each process builds typed diagrams before timing, executes three warmup calls,
  then measures complete public calls including preparation, solving and destruction.
- Samples are fresh serial processes pinned to one available CPU. A seeded
  shuffled order rotates within each block to balance positions. No outliers
  are discarded. Frequency and host load are uncontrolled and recorded as such.
- Per-process average call times are summarized with median, minimum/maximum
  and exclusive-method quartiles. Process peak RSS is reported separately in
  KiB and includes input construction and warmups.

The controller exports committed sources with `git archive` into artifact-only
snapshots. Experimental edits use checked source markers and fail if the source
shape has changed. Each snapshot records original/generated source fingerprints,
compiler commands, logs and worker hashes. The normal checkout is untouched.

## Runtime interventions

Each causal path reaches the same non-inlined solver body in the same binary.
During the first warmup, prepared coordinates, midpoint, half-lifetime and scale
are compared bit-for-bit with direct preparation. Distances must match across
all paths, processes and modes. The experiment exercises the unlimited public
API; it does not validate controlled work-count equivalence.

| Path | Intervention |
| --- | --- |
| direct | Prepare Point vectors directly from logical views |
| copy_keep | Prepare from copied coordinates and retain copies through solving |
| copy_drop | Prepare from copied coordinates and release copies before solving |
| unused_copy_keep | Copy coordinates, but prepare from the original logical views |
| reserve_keep | Reserve empty coordinate vectors, then prepare directly |
| reserve_drop | Reserve empty vectors, release after preparation and before solving |
| reserve_after_prepare | Prepare directly, then reserve empty vectors through solving |

For this fixture the two reserved coordinate capacities total 9,216 bytes.
The reserve-only paths never initialize or read coordinate elements. A speed
change on those paths does not require copying or consuming coordinate data.
It shows sensitivity to the added allocation operation; it does not identify
the particular allocator, address-layout or cache mechanism.

Counters 0-6 report inclusive preparation, solver, cleanup, call count, graph,
matching/reconstruction and grouping scopes. Nested timings must not be added.
Graph/grouping/matching scopes describe the ordinary ungrouped path; early
duplicate/direct returns may bypass them. Counters 7-9 in `growth` count graph
row capacity changes, changed nonempty-buffer addresses and live Edge bytes
at those changes. They are not allocator call counts, physical memcpy traffic
or RSS. Extra growth instrumentation perturbs timing, so it is never pooled
with ordinary or causal timing.

## Results and evidence status

The committed reproducer is being validated for a new, revision-bound run.
Earlier local exploration found a regression and allocation sensitivity, but
its uncommitted diagnostic sources are not presented as measurements of this
new harness. Measured identities and results will be recorded after the
committed run. The implementation's performance acceptance remains open.

The implementation at `9fb122ae9eb6e53c55f49ee38b356d87327d52ce` passed the local
supported distance comparison suite (108/108), using Topp
`ffa1da051ca7ac5e313c74cc9fb92a2bcb20c234`, GUDHI native bottleneck
`4ec34ac55e6d2e8cfd7c322e84c1b0a56d516d51` and GUDHI/POT weighted correctness
workers (GUDHI 3.11.0, NumPy 2.4.6, POT 0.9.6.post1). Existing numerical stress
exclusions remain exclusions. Those checks establish agreement on their
supported domain, not universal numerical or performance guarantees.

Those retained reference runs originally recorded `3314006b2881fe3125b4ccead7c9316879328354`
with a dirty worktree. Their source fingerprint
`b4c599319caa88c5775fd77db6283168590ce204d42eea3022031b4cc5563b1c`
was recomputed from the committed `9fb122ae9eb6e53c55f49ee38b356d87327d52ce`
tree and matched. The original metadata remains unchanged; this is a source
binding, not a fresh reference run. Local evidence is retained under
`target/issue14/correctness-002/`, `target/issue14/smoke-002/` and
`target/issue14/source-binding-verified.json` and is not publicly archived.

Those retained reference runs originally recorded `3314006b2881fe3125b4ccead7c9316879328354`
with a dirty worktree. Their source fingerprint
`b4c599319caa88c5775fd77db6283168590ce204d42eea3022031b4cc5563b1c`
was recomputed from the committed `9fb122ae9eb6e53c55f49ee38b356d87327d52ce`
tree and matched. The original metadata remains unchanged; this is a source
binding, not a fresh reference run. Local evidence is retained under
`target/issue14/correctness-002/`, `target/issue14/smoke-002/` and
`target/issue14/source-binding-verified.json` and is not publicly archived.

Raw artifacts remain outside Git. Each reproducer run contains `environment.json`,
`fixture.bin`, `raw.jsonl`, per-mode summaries and per-build `identity.json` /
`build.log`. Failed runs retain their raw records and status. Independent
replication should report compiler/libc/CPU, full SHAs, fixture hash, all samples
and both favorable and unfavorable paths before deciding on a production fix.
