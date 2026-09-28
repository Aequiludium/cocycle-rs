# Diagram-distance preparation and allocation effects

[Benchmark reports](README.md) / [Distance protocol](../distances/README.md)

## Question and status

This investigation accompanies [Issue 14](https://github.com/Aequiludium/cocycle-rs/issues/14).
Implementation and independent replication discussion:
[Draft PR 26](https://github.com/Aequiludium/cocycle-rs/pull/26).
The implementation removes the facade's intermediate finite/essential vectors
and prepares algorithm-owned storage from logical diagram views. Fewer copies
are not a performance guarantee. A fresh committed-harness comparison did not
reproduce the earlier stable W1 regression or the earlier timing ordering of
allocation interventions. It did reproduce differences in buffer relocation.
The change remains a draft pending independent replication and discussion.
No padding allocations or experimental selectors are added to the public crate.

The proposed counterexample is a complete dimension-zero diagram with 512 finite
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

Evidence class: local comparative experiment, with an inconclusive timing
ranking. All 224 planned process samples completed, none failed or were
excluded, all returned the expected distance, and all preparation bit checks
passed. The report commit is distinct from the measured revision.

| Identity | Recorded value |
| --- | --- |
| Measured candidate and harness | `ae7bbbe7792b7f71139784e3c44ddb92028b572c`; clean at measurement |
| Measured baseline | `75330fcf7e5b9ca764641c68b987a71bbfad8498` |
| Production source equivalence | Candidate `src/` is unchanged from `9fb122ae9eb6e53c55f49ee38b356d87327d52ce` |
| Run | `target/benchmarks/commit-ae7bbbe7792b7/diagram-distance-preparation/run-001/` |
| UTC start/end | 2026-09-28 15:01:28 / 15:04:34 |
| Build | Rust 1.91.0; Cargo release library; `rustc -O` worker; no extra RUSTFLAGS |
| Environment | Intel Core Ultra 7 155H; WSL2 Linux 6.6.87.2; glibc 2.39; Python 3.12.3; CPU 0 |
| Sampling | 14 fresh processes per path, 3 warmups and 300 timed calls per process |
| Source fingerprints | Original/generated per-file SHA-256 inventories in each build's `identity.json` |

Replay this measured implementation with `--candidate ae7bbbe7792b7f71139784e3c44ddb92028b572c`,
the same committed harness, and a fresh `run-002` directory. The worker,
controller and probe hashes are recorded in `environment.json`; experimental
libraries are generated from this commit, not ordinary production builds.

| Ordinary API | Median, us | Q1-Q3, us | Min-max, us | Median process peak RSS, KiB |
| --- | ---: | --- | --- | ---: |
| baseline | 402.37 | 374.55-517.34 | 361.47-688.91 | 3,168 |
| candidate | 378.66 | 358.14-446.33 | 342.80-644.82 | 3,344 |

The candidate/baseline median ratio is 0.941. The broad overlapping distributions
do not support either a stable regression claim or a stable speedup claim.
The small RSS difference is not evidence that materialized bytes increased:
RSS includes allocator/runtime state and input setup as well as the kernel.

| Causal path | Median, us | Q1-Q3, us | Min-max, us | Median process peak RSS, KiB |
| --- | ---: | --- | --- | ---: |
| direct | 439.65 | 368.61-539.61 | 343.58-793.19 | 3,344 |
| copy_keep | 451.98 | 354.79-646.99 | 346.55-897.37 | 3,344 |
| copy_drop | 414.42 | 362.94-493.47 | 350.03-632.24 | 3,256 |
| unused_copy_keep | 439.80 | 356.55-507.80 | 346.95-721.94 | 3,168 |
| reserve_keep | 483.22 | 366.51-621.52 | 347.34-774.46 | 3,344 |
| reserve_drop | 454.60 | 359.88-606.53 | 345.69-686.48 | 3,344 |
| reserve_after_prepare | 477.07 | 372.53-589.67 | 351.86-803.50 | 3,168 |

These distributions also overlap substantially. No intervention has established
a repeatable latency benefit in this run. Inclusive graph-generation medians
are 311.39 us for direct, 282.65 us for copy_keep and 330.03 us for reserve_keep;
the corresponding matching/reconstruction medians are 128.16, 135.54 and
135.10 us. Nested scopes and independent medians must not be summed or
subtracted to manufacture an explanation.

The separate growth diagnostic does establish a change in buffer behavior:

| Path | Capacity growths/call | Relocations/call | Relocated live Edge bytes/call |
| --- | ---: | ---: | ---: |
| direct | 2,523 | 1,024.00 | 163,840.43 |
| copy_keep | 2,523 | 1,534.29 | 229,157.12 |
| copy_drop | 2,523 | 1,536.00 | 229,376.00 |
| unused_copy_keep | 2,523 | 1,534.29 | 229,157.12 |
| reserve_keep | 2,523 | 1,534.29 | 229,157.12 |
| reserve_drop | 2,523 | 1,536.00 | 229,376.00 |
| reserve_after_prepare | 2,523 | 1,024.00 | 163,840.43 |

Entries are medians of per-process per-call averages. Equal prepared values
and capacity-growth counts do not imply equal relocation behavior. This
observation is narrower than proving that copying improves latency. It does
not locate a particular libc or cache mechanism.

Earlier local exploration motivated this case: ordinary calls measured
364.96 versus 458.36 us (24 processes), and a seven-path experimental worker
measured 466.34 us direct versus 373.99 us with copies and 369.74 us with only
reserved capacity (14 processes). Those exploratory workers were uncommitted;
these figures explain the hypothesis, not an accepted performance result of
the committed harness. Their artifacts remain local under
`target/benchmarks/commit-9fb122ae9eb6/diagram-distance-boundary/run-004-imbalance/`
and `target/issue14/causal-allocation/run-001/`. The new results above supersede
any claim that this ordering has been independently reproduced. Differences
in harness code generation, allocation history, host load and frequency have
not been isolated. Unfavorable and contradictory results are retained.

Structurally, on this finite W1 input the facade's two coordinate vectors
(9,216 bytes of logical payload) disappear. The two algorithm-owned Point
vectors still contain 18,432 bytes of payload; duplicate grouping and graph /
solver workspace remain separate and unchanged in design. These payload sizes
are not measured total allocator calls, allocator metadata, simultaneous peak
workspace or process RSS. The existing logical diagnostic worker reports
preparation and kernel capacity counters separately. No zero-copy or general
memory/latency improvement is claimed.

Raw formal-run data is local-only and has not been publicly archived. SHA-256:

- `raw.jsonl`: `dfe09accb88e76d02e6cbdee216d440bd9ffb97880fa3b46d8e7ba462e431468`.
- `environment.json`: `a3c108dd897403d91365c12787826b2bbcd946b0dff968c227b53153e188e97d`.
- `candidate/identity.json`: `6c0bccc8ae197ddbb2cc7225e573abce4bf5fb7f08c45aef004d72a1a72ce8b9`.

Independent replication should first compare the fixed fixture and prepared
values, then retain all time distributions and allocator controls. Until that
discussion is complete, the issue's integrated-performance acceptance is open.

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
