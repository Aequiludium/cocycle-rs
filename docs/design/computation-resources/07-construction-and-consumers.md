# Construction, result ownership and analysis consumers

[Design index](README.md) / [Development phase R6](development-plan.md)

These domains follow the same resource rules without requiring a common
execution object. This design preserves existing direct/prepared distinctions
and the negative repeated-matching evidence. A concrete profile/consumer selects
any later code change.

## Construction remains an independent producer

```text
borrowed geometry or supplied topology
  -> validated source access
  -> direct request OR reusable preparation OR explicit materialization
  -> consumer-specific computation
```

The alternatives are not mandatory stages. A direct persistence algorithm may
interleave construction/access and reduction. Explicit construction remains
callable without PH. A supplied complex, operator or diagram can enter its own
consumer without going through geometry or a persistence result.

| Current boundary | Preserve | Resource work to inspect |
| --- | --- | --- |
| Borrowed point/matrix views | Caller storage, metric/layout and numeric policies | Validation/traversal cost; borrowing input does not make the workflow allocation-free |
| Exact capped direct Rips request | Effective query cap and correct construction provenance | Avoid unused dense distances/retained edges while preserving exhaustive pair work when required |
| Reusable prepared graph/source | Its reusable construction range and immutable guarantees | Setup amortization and retained adjacency versus narrower direct work |
| Explicit complex expansion | Ordered topology, signed incidence and construction certificates | Requested materialization, validation reuse and potentially large output |
| Approximate/sampled construction | Original labels, sampling, metric assumptions, modified edges, blockers and approximation guarantees | Its own prepared content and callback effects; not an exact-storage optimization |

Use [the existing builder][builder] and [source adapter][source] as owners. Do
not add a generalized builder cache/snapshot service before identifying repeated
work not already handled there. Source validation/guarantees can eliminate
redundant consumer scans only when the existing immutable boundary establishes
the fact being reused.

## Effectful callbacks and execution controls

A stateful builder callback is consuming. Existing approximate preparation
caches unordered pair values according to its own contract. Reusing storage
does not permit invoking callbacks again, changing their order or restarting
them after a fallback. Pure borrowed coordinate access and arbitrary callbacks
must not be treated as equivalent optimizer inputs.

The terminal creates one budget across preparation and the requested consumer.
A construction cutoff is mathematical scope, while cancellation/limit exhaustion
is an execution error. Do not turn an interrupted construction into a successful
partial source or a narrower diagram. Keep construction coverage distinct from
the largest stored edge or the current query range.

## Results are separate from execution resources

[Common persistence data][data] owns logical diagram/context. `AsRef` borrows it
cheaply; `into_parts` moves common data/witnesses; cloning copies ownership.
The refactor must not attach producer scratch to every result or borrow a result
from storage the next computation can reset.

Preserve computed-dimension membership, interval multiplicity, ties, field,
filtration scale, approximation provenance and essential/censored distinctions.
An empty computed dimension is not an uncomputed dimension. The logical diagram
contract does not promise a physical slice layout for unchecked zero-copy reuse.
Preparation may read through supported access; it cannot infer layout or certify
an external result's provenance.

An owned output is intentional memory even if it dominates peak usage. Reuse
temporary capacity around it; do not delete requested witnesses, matrices,
histories or dense export merely to lower the memory number.

## Descriptor computations

Lifetime summaries already scan a diagram with constant auxiliary space. Keep
these as direct functions without mandatory workspace setup. Betti-curve event
allocation/sorting may be a repeated-grid preparation hypothesis, but requires a
real repeated consumer and build/query/drop evidence. Its grid, endpoint policy
and intended output remain part of the descriptor request.

If shared event preparation is justified, place it in the descriptor owner and
have ordinary/admitted repeated paths use the same numerical kernel. Do not make
the producing PH engine or its workspace a dependency of diagram analysis.

## Matching: operand content versus pair-local work

Observed [distance facades][distance] validate context/coverage and extract
finite/essential points per call. Private matching owners then perform their
specialized graph/solver choices. Preserve those owners and existing adaptive
rules unless a separate selector experiment changes an admitted region.

| Candidate reusable resource | What must match | Boundary |
| --- | --- | --- |
| Scratch capacity | Typed layout/capacity and non-overlapping lifetime | Can serve different pairs after initialization |
| Operand-local extracted data | Actual operand, dimension, endpoint/context interpretation | Can survive repeated use only while those dependencies stay valid |
| Pair normalization/costs | Both operands and the selected metric/numeric policy | Not automatically reusable with just one repeated operand |
| Matching graph, potentials and mutable search state | Their actual pair-specific algorithm contract | Keep local; any warm/incremental method is a separate justified change |

[B3 lifecycle evidence][b3] rejected the tested public Workspace/Batch/pool
promotion. [B5 scoped preparation][b5] found only conditional survivors and did
not ship a public prepared API. This refactor does not reverse either outcome.
A new proposal must name a different justified boundary/consumer and compare
against the current ordinary implementation, including setup, normalization,
reset, outputs, destruction and cumulative retention.

Do not reserve worst-case dense graphs for sparse/small pairs merely to avoid
growth. Test single pair, same pair, one-to-many, distinct pairs, all pairs and
mixed-size scopes. Preparation wins and scratch-pooling wins are separate claims.
A scalar Hausdorff or other constant-space geometry route, if later adopted,
does not need a reusable context just to fit this design.

## Acceptance and normal-path integration

Affected construction/PH changes retain independent GUDHI and supported Ripser
validation. Matching retains the declared GUDHI/Topp cases and independent tiny
matching expectations. Descriptor tests verify their own formula/grid contracts.
Include signed/unequal births, approximate blockers, effectful callbacks, empty
and uncomputed dimensions, essential counts, censored endpoints, repeated
intervals, extreme magnitudes, cancellation and live-output retention as relevant.

An admitted change belongs inside the existing source or consumer owner and its
ordinary entry point. A special repeated-use handle requires its own demonstrated
consumer and public compatibility review. R6 may end in a no-adoption decision
for a proposed cache/pool; preserving a cheaper direct path is an intentional
architecture outcome, not incomplete migration to a universal workspace.

[builder]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/filtration/rips/builder.rs
[source]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/persistence/rips/mod.rs
[data]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/diagram/data.rs
[distance]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/diagram_distances/mod.rs
[b3]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/benches/reports/distance-lifecycle.md
[b5]: https://github.com/Aequiludium/cocycle-rs/blob/7b9c2ab75d353b4b3737439e1c19149974716874/benches/reports/scoped-preparation.md
