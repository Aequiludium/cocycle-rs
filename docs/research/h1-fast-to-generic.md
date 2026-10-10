# H1 fast-engine transfer to generic F2

[Documentation](../README.md) / Research / [Phase-3 T1](https://github.com/Aequiludium/cocycle-rs/issues/35)

This audit reads frozen R0 `a1a7cb629187fbc17c9e2e638ee8ae4bf19c2656`.
The [baseline report](../../benches/reports/h2-foundation.md) owns measurements,
finite coverage and source identities. Existing generic H2 support is not a new
algorithm. The decisions below authorize a minimal clearing continuation, not
production admission of a new H2 fast kernel.

## Transfer matrix

| Mechanism and source | Dependency/invariant | H2 or arbitrary-dimensional decision | Cost hypothesis and remaining obligation |
| --- | --- | --- | --- |
| [Union-find and cycle_edges](../../src/persistence/flag/cohomology/mod.rs) | Zero-born H0 edges; same forward order as H1 | H1-specific preprocessing; retain, do not generalize union-find to H2 | Linear edge scan plus inverse-Ackermann merges; eliminates H0 death columns |
| [Clearing](../../src/persistence/simplicial/cohomology.rs) | Complete preceding-dimension death simplex set | Already generic; transfer specialized H1 pairing state | O(number of deaths) handoff; preserve topology and prove key/order equality |
| [Implicit cofacets](../../src/filtration/flag/cliques.rs) | Enumerate every legal codimension-one coface within inclusive cutoff | Already generic; retain | Dense vertex scans or shortest sparse adjacency; count rejected candidates and edge tests |
| [Transform reconstruction](../../src/persistence/simplicial/cohomology.rs) | R=CV, triangular V, normalized unique pivots | Already generic, including prime fields | Stored transform versus regeneration tradeoff; fill-in can be exponential |
| [Apparent pairs](../../src/persistence/flag/cohomology/mod.rs) | Mutual earliest-cofacet/latest-facet; equal values for omission | Dimension-generic theorem; H2 implementation deferred to T4 | One facet/cofacet certificate; prove higher-dimensional traversal selects the same first row |
| [Initial emergent shortcut](../../src/persistence/flag/cohomology/mod.rs) | First equal-valued row, before additions, no stored or virtual owner | H2 transfer unknown until T4; never infer from ties alone | Avoid remaining row generation; ownership predicate and one-pass ordering must be proved |
| [Virtual apparent owners](../../src/persistence/flag/cohomology/mod.rs) | Mutual certificate plus owner strictly precedes active column in reverse order | Generic mathematical possibility; H2 storage experiment deferred to T5 | Less storage, more regeneration; complete clearing must survive omitted owner storage |
| [Parity heaps](../../src/persistence/flag/cohomology/mod.rs) | F2 multiplicities cancel in pairs, including transform entries | F2-specific, dimension-independent | Potentially smaller coefficient cost; must cancel all copies, not deduplicate |
| [One-pass initialization](../../src/persistence/flag/cohomology/mod.rs) | Enumeration in decreasing ID; first equal row is earliest | Fixed traversal property, not automatically generic | Avoid repeated prefix scans; H2 visitor order and rejected-candidate accounting need tests |
| [Dense cone cutoff](../../src/filtration/rips/mod.rs) | Complete dense zero-born flag filtration; a cone at minimum row maximum | Dimension-independent on this exact source; retain current scope | Stop above acyclic cone scale; invalid for an arbitrary sparse graph or blockers |
| [EdgePosition/ColumnPosition](../../src/persistence/flag/cohomology/mod.rs) | Positions owned by one ordered edge/transform array | H1 storage-specific; never export them | Small keys; decode before releasing owning arrays |
| [Combinatorial IDs](../../src/filtration/flag/index.rs) | Checked binomial arithmetic fits usize | H1 index is fixed-dimensional; tuple continuation avoids higher binomial dependency | Overflow-safe fallback required when selecting specialized H1 for formerly generic requests |
| [Dense traversal](../../src/filtration/flag/dense.rs) | Descending candidate vertices give decreasing triangle ID | Fixed edge-to-triangle visitor; do not presume H2 order | Reject faces/cutoff candidates; tetrahedral order must be checked independently |
| [Sparse intersection](../../src/filtration/flag/sparse.rs) | Reverse sorted intersection of two adjacency lists | H1-specific two-list traversal; generic common-neighbor analogue already exists | Degree-dependent scans; multi-list/lazy experiments deferred to T7 |
| [Execution and reservations](../../src/execution/mod.rs) | One operation budget; fallible Vec/HashMap reservations where supported | Generic contract; retain | Reusing state cannot reset budget or return partial output; BTreeMap allocations remain the existing limitation |

No transferred optimization is admitted on naming or analogy alone. H0/H1
specialization remains independent of odd-prime, representative and sparse
approximation paths. Source/literature authority for pairing is
[Bauer 2021 sections 3.2-3.5](../reference/bibliography.md#b21), rather than a new
novelty claim. Dory/Ripser source-specific experiments remain T6/T7 work.

## Ordering and clearing argument

For sorted vertex tuples, forward order is increasing maximum edge value, then
dimension, then decreasing colexicographic vertices. On a fixed dimension,
colex order agrees exactly with decreasing combinatorial ID. In particular,
triangles and tetrahedra use the same rule; original graph IDs are never relabeled.
Faces precede cofaces on equal values. Reverse column processing and earliest
forward coface pivots implement the reversed-transpose reduction specified in
[mathematics section 9](../reference/mathematics.md#9-implicit-rips-persistent-cohomology).

Adjacent-dimension clearing removes the coboundary columns of death triangles
paired with H1 edges. Their columns reduce to zero in H2 and must not become H2
births. The set is the union of ordinary pivots, stored initial shortcuts
(including emergent), and omitted zero apparent pairs. Public intervals omit
zero bars, so intervals alone also cannot reconstruct the set. Virtual owner
reconstruction during a later addition does not create another death simplex.

The mutual apparent certificate makes a row-column pair available without
earlier elimination. Equal values justify suppressing its public bar, not
suppressing its clearing key. Reverse-order legality requires its edge position
to be strictly later in forward order than the active edge. This is checked in
the current H1 reducer; no such unchecked H2 shortcut is introduced by T2.

An equal-weight K3 is a direct counterexample to empty clearing: H1 pairs its
cycle edge with triangle 012 at zero lifetime. With no tetrahedra, reducing 012
again as an uncleared H2 column would incorrectly produce [1,infinity).
Correct H2 is empty. The independent exhaustive boundary tests include this
case. Missing virtual pairs can therefore change public mathematics, not just
runtime. Related [#30](https://github.com/Aequiludium/cocycle-rs/issues/30) remains
independent; a dimension range alone neither fixes clearing nor proves no H1
recomputation.

Clearing is a reduction-column predicate. Every triangle, including cleared
ones, remains available for unique tetrahedron generation. Deleting cleared
topology instead would break higher-dimensional enumeration.

## Minimal representation and ownership

T2 exports an owned interval vector and Vec<[usize;3]> of sorted original death
triangle vertices. usize avoids vertex-width narrowing; H1's checked binomial
index must still fit. If index construction overflows for a higher-dimensional
request, retain the established tuple-based generic fallback. Allocation errors
and cancellation propagate, rather than being reclassified as overflow.

The existing generic H2 continuation keeps Simplex/Vec<usize> triangles and
tetrahedra. Fixed [usize;3]/[usize;4] keys are the selected safe T3 baseline,
without u32 narrowing, unsafe reinterpretation or a public Workspace API.
They are not an implemented optimized reducer in this PR.

Decode keys while H1 access is live. Once its function returns, edge arrays,
cycle flags, owners, transform columns and heaps have reached their last use;
no returned key borrows them. Shared graph/matrix, topology, result coverage and
the WorkBudget remain necessary independent inputs. Continuation enumerates
edges/triangles as topology, but never reruns generic union-find or H1 reduction.
Its hash-key conversion needs owned three-vertex vectors with fallible reserves.
M1 measures actual conversion/topology overlap and release effects later;
drop semantics alone is not a process-RSS improvement claim.

## Section 39 entry gates

| Original gate | Answer/evidence |
| --- | --- |
| H1 inventory and classification | Complete matrix above; existing mechanisms separated from proposed transfers |
| H1-to-H2 clearing | Complete union of ordinary/emergent/virtual deaths; K3 counterexample and boundary-oracle coverage |
| Triangle/tetrahedron ordering | Value, dimension, decreasing colex; usize tuples retain original IDs |
| H2 apparent invariant | Mutual earliest-cofacet/latest-facet plus equal values for omission; implementation/ablation remains T4 |
| Baseline representation | Existing generic tuples for continuation; fixed usize triples/quads for T3 |
| Workload matrix | 23 cases at n=16 and n=32, 10 H1 controls each, fixed seed and exact graph families |
| Generic and Ripser H2 baselines | Both R0 pipeline runs passed 12 measured rounds and every sample's output validation |
| Dory H2 baseline | Audited geometric square-matrix scope; tied failures retained, no full-domain oracle claim |
| Instrumentation | Test-only per-dimension reduction and rejected-candidate counters; resource/ownership schema in baseline report |

These gates support the minimal exact continuation on the GUDHI/Ripser and
independent-oracle domain. Dory's failed tied cases prohibit treating that engine
as an oracle for those inputs. New H2 shortcuts, performance admission, typed
reuse and arbitrary-dimensional fast-kernel claims remain their later tasks.
