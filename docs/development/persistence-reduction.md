# Contributing a persistence algorithm

[Documentation](../README.md) / Development / [Contribution paths](algorithm-contributions.md)

Start with the algorithm's mathematical input and output. A boundary-matrix
researcher can work directly on ordered columns; an implicit cohomology researcher
can work on coface access. Neither implementation calls the default Builder.
This guide covers in-crate contributions. The internal types and functions below
are reusable implementation tools, not a stable external reducer API.

## Find the relevant computation

| Mathematical work | Implementation | Input and output |
| --- | --- | --- |
| Forward column reduction | [algebra/reduction/boundary.rs](../../src/algebra/reduction/boundary.rs) | Owned columns over a prime field to reduced columns, pairings and optional transformations |
| Ordinary intervals from selected boundaries | [persistence/boundary](../../src/persistence/boundary/mod.rs) | Selected ordered columns with dimensions/values to a diagram |
| Zero-born H0 | [flag/h0.rs](../../src/persistence/flag/h0.rs) | Vertex count and weighted edges to merge intervals |
| Specialized implicit F2 H1 | [flag/cohomology](../../src/persistence/flag/cohomology/mod.rs) | Ordered edge/triangle access to ordinary intervals |
| Ordered-cursor F2 H1 | [flag/ordered.rs](../../src/persistence/flag/ordered.rs) | Dense matrix and prepared edges to ordinary intervals |
| Empty coboundary rejection | [filtration/flag/bitset.rs](../../src/filtration/flag/bitset.rs) | Bounded neighbor bitsets over retained weighted edges |
| Prime-field implicit cohomology | [simplicial/cohomology.rs](../../src/persistence/simplicial/cohomology.rs) | Zero-born coface access to ordinary intervals |

The [mathematical specification](../reference/mathematics.md#4-persistence-and-reference-boundary-reduction)
defines orientation and persistence pairing. Verify the applicable assumptions
before sharing an implementation: zero-born union-find pairing does not implement
arbitrary vertex-birth persistence, and a sparse Rips blocker cannot be replaced
by ordinary clique enumeration.

The exact matrix F2 H1 adapter prepares sorted edges once. For at least 64
vertices, it selects adjacency access when no more than approximately two thirds
of the possible edges remain at the internal stopping scale; otherwise it keeps
borrowed matrix access. This private heuristic includes graph preparation in the
operation budget and does not change original-input coverage. Both paths reuse
the prepared edge array. Sparse intersections reject disjoint sorted vertex
ranges before scanning their neighbors. Its reverse neighbor iterators preserve
descending cofacet order and one checkpoint per intersection step. Dense visits
select the matrix layout outside the vertex loop and use checked simplex index
prefixes for safe slice access.

The matrix adjacency adapter and supplied/threshold graph adapters share the
same neighbor-bitset selection. They can attach neighbor bitsets when bounded probes
find many empty intersections that sorted-range rejection cannot detect. The
cache has at most one u64 word per graph edge, and is considered only for 64 to
4096 vertices. A zero bitwise intersection certifies an empty retained
coboundary. Nonempty intersections use the original sorted traversal, preserving
cofacet order and weights. Probe and cache work share the caller's budget;
the heuristic selects an exact strategy and never approximates topology.

For at least 64 vertices, approximately 90% retained edges and at least three
distinct retained edge weights, the matrix adapter selects ordered coboundary
cursors. Sorted neighbor tables let each source column produce its triangles in
the existing value/colex order. The working heap holds one head per live cursor;
equal heads cancel by F2 parity. Stored owners can seek past prefixes whose XOR
is zero. See the [cursor invariants](../reference/mathematics.md#ordered-coboundary-cursors).
Table preparation shares the operation budget. This is also a private heuristic;
it does not change f64 values or coverage, and preparation must be included in
performance measurements.

The cost rules have named private constants in
[flag/selection.rs](../../src/persistence/flag/selection.rs) and
[flag/bitset.rs](../../src/filtration/flag/bitset.rs). The matrix gates retain the
run-003/run-004 local native experiments' policy; the bitset probe avoided a
measured regression from always building the cache. These are empirical choices,
not geometry classifications or correctness assumptions. They do not recognize
dataset names, and no single threshold is claimed optimal on every machine.

| Cost rule | Current policy | Purpose |
| --- | --- | --- |
| Adaptive matrix access | At least 64 vertices | Amortize representation preparation |
| Sparse matrix access | m <= P - floor(P/3), P = n(n-1)/2 | Limit adjacency overhead on dense inputs |
| Ordered matrix access | m >= P - floor(P/10), at least three retained weights | Amortize sorted cursor tables and avoid degenerate filtrations |
| Neighbor cache | 64 <= n <= 4096, n ceil(n/64) <= graph edge count | Bound auxiliary words by O(n+m) |
| Cache probe | At most 32 strided edges; at least half have costly empty intersections | Avoid caches on cheap scans or overlapping neighborhoods |
| Costly scan | Minimum endpoint degree at least twice the row word count, overlapping vertex ranges | Compare bitset work to sorted scan work |

All m counts above use the internal retained scale, not the original matrix's
nominal completeness. Cutoff and coverage remain separate. Changing this policy
requires correctness checks and comparable before/after measurements under the
[native protocol](../../benches/protocol.md) and
[reporting rules](../../benches/reporting.md), including unfavorable workloads.

Both F2 H1 engines use [flag/edges.rs](../../src/persistence/flag/edges.rs) for
forward H0 classification and typed edge/column positions. Engine-specific heaps
and traversal remain separate. Cursor positions and their single-row sentinel
are local to the ordered engine.

For F2 requests above H1, dispatch runs the selected H0/H1 engine and continues
with generic cohomology only from H2. [flag/clearing.rs](../../src/persistence/flag/clearing.rs)
retains every H1 death triangle, including zero-lifetime and omitted apparent
pairs. The handoff decodes IDs into ordered vertex tuples for the same access,
cutoff and filtration order. Starting H2 with an empty clearing set when H1 has
nonzero rank would create false H2 births. No H1 pivots certify the absence of
triangles and allow the higher part to return empty without repeating preparation.
Higher dimensions retain ordinary clearing. Odd primes keep the generic path;
high-dimensional requests also retain the tuple-based fallback if compact H1
triangle IDs cannot fit usize. Composition and handoff share one controlled
budget, and failures never return the already computed lower-dimensional part.
No public backend or tuning option is added.

Exact F2 H1 reducers specialize the operation budget with a const generic.
An operation with neither a limit nor a cancellation flag uses the no-op
specialization. Any real control keeps the original shared budget, including
preparation, probes and reduction. Ordered cursors also reuse the checked simplex
index prefixes for safe reads in all three matrix layouts, preserving canonical
f64 values without repeated fallible pair-count calculations.

The ordinary F2 H1 working heap stores canonical finite nonnegative values as their exact
f64 bits. Its equality test preserves F2 multiplicity, and its boolean comparisons
use a finite-value projection of the shared filtration order. Reconstructed
cofacets are collected in a reusable fallibly allocated buffer and extended into
the heap in batches, with cancellation checks around heap maintenance. Stored
transformations and mathematical reduction order are unchanged.

## Work directly on boundary columns

The algebra reducer consumes its columns. Its pivots and transforms belong to
that invocation. It knows no point cloud, source context, diagram type or Builder.
It accepts a checkpoint callback; maintainers connect that callback to the
operation's existing work budget instead of starting another budget.

`reduce_pairs` retains reduced columns and death partners without basis
transformations. `reduce` additionally returns V with D V = R. Nonzero reduced
columns have unique normalized pivots. The `deaths[i]` entry associates birth
column i with its death column. An empty reduced column is a birth candidate;
later columns may pair it, so do not declare it essential during the forward scan.

For ordinary interval assembly, `boundary::diagram` takes owned `BoundaryInput`,
maximum dimension, field, established coverage and the current budget. Readers
select the q+1 skeleton first. `BoundaryInput` holds ordered IDs, dimensions,
values and boundary columns; its shared append operation keeps these arrays in
step. It contains no source certificate or Builder settings. This helper declares
the contiguous output domain `0..=q` used by current default adapters.

The two existing readers establish their guarantees differently:

- [filtered.rs](../../src/persistence/filtered.rs) validates external cell IDs,
  ordering, selected-field incidence and the boundary-square condition.
- [simplicial/input.rs](../../src/persistence/simplicial/input.rs) reuses frozen
  simplicial construction guarantees and converts stored oriented incidence.

An algorithm with a different column representation, traversal, preparation or
workspace can own it. It can use the algebra reducer without using `BoundaryInput`,
or implement another reducer. Shared ordinary output types do not require shared
internal storage. For H1-only or gapped outputs, declare `ComputedDimensions`
explicitly rather than using a contiguous maximum as an availability claim.

## Derive a small case before integrating

Consider three vertices born at different signed scales, then their edges and a
delayed face. Increasing vertex order fixes orientation:

| Column | Cell | Value | Boundary |
| --- | --- | --- | --- |
| 0 | a | -4 | 0 |
| 1 | b | -3 | 0 |
| 2 | c | -2 | 0 |
| 3 | ab | 0 | b - a |
| 4 | ac | 1 | c - a |
| 5 | bc | 2 | c - b |
| 6 | abc | 3 | bc - ac + ab |

The complete diagram has H0 intervals `[-4, infinity)`, `[-3, 0)` and `[-2, 1)`,
and one H1 interval `[2, 3)`. These endpoints follow from component merges and
the delayed filling; they are not obtained from a second production solver.

[boundary/tests.rs](../../src/persistence/boundary/tests.rs) constructs these
columns directly and calls the computation without a source adapter. It also
checks D V = R using independent dense modular arithmetic, normalized unique
pivots, field-dependent cellular incidence and recovery after work exhaustion.
The existing test-only [reference reducer](../../src/persistence/reference/mod.rs)
remains independent; do not reuse production elimination to generate its answers.

Once the local mathematics works, this public composition checks source integration:

```rust
use cocycle::{Result, algebra::PrimeField, complex::{Simplex, SimplicialComplex}};
use cocycle::diagram::{Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval};
use cocycle::persistence::PersistenceExt;
let cells = [
    (vec![10], -4.), (vec![20], -3.), (vec![30], -2.),
    (vec![10, 20], 0.), (vec![10, 30], 1.), (vec![20, 30], 2.),
    (vec![10, 20, 30], 3.),
].into_iter().map(|(vertices, value)| Simplex::new(vertices, value))
 .collect::<Result<Vec<_>>>()?;
let complex = SimplicialComplex::new(cells)?;
let result = complex.persistence().field(PrimeField::new(3)?).compute()?;
let expected = PersistenceDiagram::new(1, Coverage::Complete, vec![
    PersistenceInterval::new(0, -4., IntervalEnd::Essential)?,
    PersistenceInterval::new(0, -3., IntervalEnd::Finite(0.))?,
    PersistenceInterval::new(0, -2., IntervalEnd::Finite(1.))?,
    PersistenceInterval::new(1, 2., IntervalEnd::Finite(3.))?,
])?;
assert_eq!(result.diagram(), &expected);
# Ok::<(), cocycle::Error>(())
```

## Connect a default operation when appropriate

`source.rs` handles supported sources and context. `flag/dispatch.rs` selects H0,
specialized F2 H1 or general cohomology for exact flag inputs; algorithms never
call back into it. Explicit and approximate Rips adapters call simplicial
cohomology directly. `simplicial::finish_zero_born` implements the existing choice
between an implicit diagram and materialized simplicial representatives, only
for compatible zero-born access. It is not a mandatory pipeline for new algorithms.

Keep capability selection explicit. Unsupported fields, dimensions, source
coverage or witnesses must not silently change the requested operation. Reuse
mathematical parameters when they fit. Add a specialized public entry only when
users need a meaningful choice or a different operation; an internal optimization
does not require another Builder or registry.

Maintainers assemble library-produced `PersistenceData` from the diagram and
actual source context. Compatible results can compose it with their own auxiliary
output. Existing `Representative` means persistent simplicial cycles or their
query-scale dual cocycles; arbitrary cellular witnesses need their own semantics.
Keep source certificates, requested settings and the actual computed domain
distinct. See [result composition](../design/kernel.md#results-carry-enough-information-to-be-interpreted).

Run the focused check while developing:

```sh
python3 tools/check_algorithm.py persistence-reduction
```

It runs direct persistence algorithm tests, independent Rips oracle comparisons,
filtered-input/field/representative/resource integration tests, the existing flag
example and this guide's code. It needs no C++ installation. Maintainers complete
the [full checks](../../CONTRIBUTING.md#verification) and relevant pinned native
comparisons before merge. Timing evidence follows the benchmark protocol; ordinary
tests do not assert machine-dependent performance thresholds.
