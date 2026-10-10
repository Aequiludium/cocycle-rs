# Data modeling: first objects and stable reads

[Documentation](../README.md) / Research / [#65](https://github.com/Aequiludium/cocycle-rs/issues/65)

Research decision, 2026-10-10. The goal is to let independently validated data
serve filtration and direct single-scale consumers with the same cell meanings.
This selects the first modeling entry for
[#66](https://github.com/Aequiludium/cocycle-rs/issues/66); it does not implement
every topological domain or a spectral solver.

## Decision

Start with existing `PointCloudView`, `WeightedGraph` and `SimplicialComplex`.
Keep geometry, graph topology and supplied simplicial topology distinct.
Reuse `FilteredComplex::cells()` to enumerate the frozen complex's IDs, then
borrow `simplex(id)` and `boundary(id)`. This already avoids looking every ID
up again by its vertex tuple. No new source trait, container, identity registry,
public matrix type or iterator API is needed for the selected consumers.

Implement the missing runnable handoff in the existing construction example and
[interface guide](../guides/interface-contracts.md). Read signed adjacent-degree
boundaries without computing persistence; separately run the existing persistence
builder on that same owner. Export basis IDs alongside the small matrix and
retain its shape when a dimension is empty. The single-scale read can ignore
stored filtration values; these are not spectral weights.

The [#63 compatibility decision](../design/interface-compatibility.md) is already
on the baseline. The current investigation confirms its existing-API reuse
choice rather than inventing a new missing primitive. The pending filtration
stage/map PR #100 is independent of this entry and is not required here.

## Sources actually read

| Source | Fixed revision and inspected material | Conclusion and limit |
| --- | --- | --- |
| Cocycle | main `a7ee9e911dcdeaea137b682bf383f7c2e1d090d5`; [points](../../src/geometry/point_cloud.rs), [graphs](../../src/complex/graph/mod.rs), [simplicial storage](../../src/complex/simplicial/mod.rs), [orientation](../../src/complex/simplicial/simplex.rs), [filtered contract](../../src/complex/filtered.rs), [supplied complexes](../guides/filtered-complexes.md), [input decisions](input-contracts.md) | Borrow validated existing owners. The filtered trait exposes IDs already; vertex labels and filtration positions are different identities. |
| TopoNetX | `ce3d414b571fff95916c4cbc20b8d2c2c9d428e6`; [simplicial incidence and Hodge](https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/toponetx/classes/simplicial_complex.py), [cell storage and incidence](https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/toponetx/classes/cell_complex.py) | Retain row/column keys when exporting incidence. Its cell owner is a dynamic 2D polygon domain, including nonregular cells; it is not a universal replacement for frozen simplices. Source inspection, not a TopoNetX execution result. |
| GUDHI | `384680c0a34a749d1686d5d59c551c55578d815e`; [Alpha definition/construction](https://github.com/GUDHI/gudhi-devel/blob/384680c0a34a749d1686d5d59c551c55578d815e/src/Alpha_complex/doc/Intro_alpha_complex.h), [bitmap cubical definition/storage/input](https://github.com/GUDHI/gudhi-devel/blob/384680c0a34a749d1686d5d59c551c55578d815e/src/Bitmap_cubical_complex/doc/Gudhi_Cubical_Complex_doc.h) | Alpha has squared-radius scales and Delaunay/predicate requirements. Cubes have product-coordinate identities, shape and axis order; they need not be triangulated. Reading these definitions does not deliver a Rust constructor. |

The two TopoNetX files were rechecked against GitHub blob IDs
`8ace5fef7ad5e9a44cce02b284578444340bbae9` and
`f807495dc56b31b4cb63466f1c86b90052b5a733`. The GUDHI document blob IDs are
`4f2e0dfd39f0df6ab5519962f38337ad9f0c7612` and
`ceaf16e704d07e3b6e36a39cad2f958c72026923`. These pin research sources;
reference-runtime validation has its own package/compiler identity.

TopoNetX's studied simplicial `incidence_matrix(0)` returns an augmented one-row
matrix; Cocycle's ordinary H0 boundary has zero rows. Its simplicial incidence
implementation accepts `weight` without using it, and up/down Laplacians reject
non-None weights. Preserve these source-specific limits from
[the full input study](input-contracts.md), rather than adopting a parameter name
as evidence of weighted spectral semantics. TopoNetX's polygon boundary follows
cyclic edges; a polygon is not a simplex determined by an unordered vertex set.

## Three minimal use cases and requirements

| Scenario and exact small input | Validation | Required queries and consumer meaning |
| --- | --- | --- |
| Point measurements: row-major `(0,0),(1,0),(0,1)`, ambient dimension 2 | Positive ambient dimension, exact checked shape and finite coordinates. Duplicates remain distinct rows. | Borrow original rows and point indices. A Rips consumer chooses a metric/threshold explicitly; a direct geometric consumer reads coordinates without creating edges or PH. Coordinates alone are not topology. |
| Sparse network: four vertices, edges `(1,0):7` and `(2,1):11`, vertex 3 isolated | Valid distinct endpoints, no loops or repeated undirected edges, finite nonnegative values. Preserve zero-valued edges and isolates. | Borrow canonical edges and adjacency. A direct graph consumer uses only vertices/edges; orient `(a,b)` with `a<b` as `b-a`. A flag consumer explicitly adds cliques. Values retain filtration meaning; an application may explicitly choose a separate conductance or Gram. |
| Supplied topology: labels 10/20/30, vertices at 0, three edges at 1, face at 2 | Increasing nonempty vertex tuples, unique simplices, every face supplied and face value <= coface value. Finite signed values are supported; missing faces are rejected. | Borrow owner-local IDs, labels, degree and signed boundary. PH gives H1 `[1,2)`; at scale 1 a direct chain reader obtains a cycle, and at scale 2 the face fills it. No PH result is required to read the chain data. |

For the triangle, restrict source traversal to each degree. The vertex basis is
`((30),(20),(10))`, the edge basis is `((20,30),(10,30),(10,20))`, and the face
basis is `((10,20,30))`. The hand-derived integer boundaries are

```text
B1 = [[1,1,0],[-1,0,1],[0,-1,-1]], B2 = [[1],[-1],[1]].
B1 B2 = 0.
```

At scale 1, B2 has shape 3-by-0, not an absent query. At scale 2 it has shape
3-by-1. With explicitly chosen standard real inner products, the full triangle
has edge Hodge matrix `B1^T B1 + B2 B2^T = 3 I`. This is a finite hand calculation
for an input-reading example, not a production eigenproblem solver. Reducing
signed coefficients modulo two is explicit; replacing arbitrary F2 input by
real 0/1 entries does not establish an integer chain complex.

## Identity, geometry, orientation and weights

| Information | First-entry contract |
| --- | --- |
| Source and snapshot | A borrowed owner binds the local read. Caller-supplied dataset labels, units and coordinate associations remain caller responsibilities; the owner does not authenticate them. Keep them with exported matrices/results across storage or cache boundaries. |
| Point/vertex identity | Point rows and graph vertices are `0..n`; supplied simplex labels can have gaps. To associate coordinates with labels 10/20/30, retain an explicit label-to-row map. Never index a coordinate array by `SimplexId`. |
| Cell identity | `SimplexId` is a position in one immutable complex's filtration order. Do not transfer it between owners, even if its integer index is in range. Rebuilding/reordering a complex requires remapping by declared cell keys. |
| Matrix coordinates | Every exported row/column basis retains its `SimplexId`; coefficients use that declared order. Reordering a basis also reorders attached costs/attributes. Shape equality alone is insufficient. |
| Orientation | Increasing simplex vertices and omission sign `(-1)^i`; graph edge orientation must be declared explicitly. A reversal/permutation needs an explicit signed map. |
| Values and weights | Filtration values are not coordinate costs, conductances or inner products. Spectral consumers choose positive-definite Grams/adjoints; F2 class consumers choose their own positive coordinate costs. Neither is inferred from stored values. |
| Coverage | A supplied complex is exactly the supplied topology. Its ID traversal certifies neither completeness of a larger Rips/Alpha source nor omitted cells. Certified explicit Rips owners retain their construction context separately. |

## Frozen storage versus local reads

The selected supplied-complex path owns all simplices and stored incidence.
Construction requires explicit face closure and costs storage proportional to
cells plus codimension-one incidences. It is suitable when the same object serves
several consumers, with immutable IDs throughout those reads.

`cells()`, `simplex(id)`, `boundary(id)`, graph edges/adjacency and coordinate rows
borrow existing buffers. The selected example allocates only degree-basis handles
and the requested adjacent matrices. It does not clone the complete complex or
perform clique expansion for a local read. Selecting a degree from `cells()` scans
the existing cells; this is not an indexed-degree lookup or a sublinear claim.
Dense export is O(rows * columns) storage and is intended only for small examples.
Large consumers should iterate signed sparse terms and retain basis maps.

Implicit Rips persistence already reconstructs needed cofaces privately. Keep
that route for PH workloads that do not need explicit storage; a `build_complex`
request deliberately materializes topology and may grow exponentially. Do not
expose expensive lazy topology construction through this frozen read interface.
No new generic on-demand container is justified by these three use cases.

## Expansion order and exclusions

| Order / decision | Benefit, cost and acceptance condition |
| --- | --- |
| First: existing points, graph and frozen simplex owners | Immediate borrowed reads and existing validation; no new dependency or container. #66 must demonstrate both readers on one object, empty dimensions, owner-local IDs, orientation, values and copy boundaries. |
| Also available: supplied Alpha simplices | Existing storage accepts independently constructed face-closed Alpha data. Retain original labels, squared-radius convention, precision and source completeness externally. Import is not a native Alpha geometry constructor. |
| Next candidate: nonperiodic rectangular cubical input | Adds image/grid data that the three existing owners do not directly represent. Begin only with a concrete image workload and defined grid shape, flattening/axis order, top-cell versus vertex values, product orientation and cell keys. Validate integer boundary-square zero before PH. Prefer implicit grid incidence over automatic triangulation. |
| Later candidate: native unweighted low-dimensional Alpha construction | Reuses simplex output but requires robust Delaunay/predicates, degeneracy/duplicate policy, precision, source IDs and squared-radius scale. GUDHI's implementation uses CGAL and documents differing numerical guarantees. Select a concrete geometry implementation and dependency/MSRV/precision policy first; no new native dependency is adopted here. |
| Defer general polygon/cellular storage | Existing `FilteredComplex` can consume externally validated finite integer boundaries for PH, including the guide's square. A production owner must separately specify attaching maps, repeated traversals, cell multiplicity and nonregular identity. Vertex sets alone lose that information. A 2D TopoNetX owner does not prove an all-dimensional representation. |
| Exclude from this entry | Periodic cubes, weighted Alpha, sheaf/path/hypergraph/Mayer objects, automatic triangulation, native spectral/Dirac/projector solvers and learning integrations. Each has its own mathematical contract and concrete workload gate. |

Cubical input is the next *conditional* modeling investigation because it adds a
distinct data family without first solving a general geometry problem. There is
no delivery-date or universal performance claim. GUDHI's cubical format permits
infinity for missing cubes; Cocycle's current filtered-cell values must be finite.
An eventual adapter must explicitly omit/remap absent cells and establish face
closure rather than silently importing infinity. Periodicity changes cell
identifications and remains a separate gate.

## Acceptance and remaining evidence

The research delivers three concrete requirement rows, a direct non-PH example,
source-pinned identity/orientation rules, a reuse decision and conditional
extension order. The first runnable handoff belongs to #66 and its exact PR
validation; the research does not count proposed extensions as implemented.
The [construction example](../../examples/complex_construction.rs) and interface
guide own executable behavior and finite expectations. Local checks, external
comparison, hosted CI, merge and release must be recorded separately with the
submitted source identity. No performance improvement follows from this decision.

The initial local three-object probe passed on Windows with Rust 1.98.1: it checked
coordinate-buffer borrowing and invalid/empty point input, canonical graph
endpoints/values/isolate/invalid vertex lookup, and frozen simplex ID-to-reference
identity, signed triangle incidence and integer boundary-square zero. Source and
Markdown checks passed. This is finite object/contract evidence, not execution
of TopoNetX, Alpha geometry or a cubical constructor.

For #66, the existing construction example now reads one circle directly as
signed incidence and `B1^T B1` with standard inner products, then computes PH
from the same owner. Its six tests include adjacent-degree identity, gapped
vertex labels, signed integer chains, unchanged owner storage, F2/F3 persistence,
empty 0-by-n/n-by-0 shapes and out-of-range ID lookup. The interface guide also
executes the scale-1 triangle energy and its full-filtration H1 `[1,2)` read.
Both debug/release example tests and the three guide doctests passed; all library
debug/release tests, all ten examples, strict rustdoc, formatting and Clippy
passed. Rust 1.91 passed the full tests, all-target checks, six construction
example tests and the guide doctests. These local results precede hosted CI.

An additional Python 3.10.11 / GUDHI 3.12.0 check compared the circle, triangle,
empty complex and two isolates over F2/F3: eight exact diagram comparisons
passed, including multiplicity, all dimensions and essential endpoints. The
triangle's two H0 `[0,1)` merges are retained alongside essential H0 and H1
`[1,2)`. Its initial hand-expectation table omitted those two H0 merges; that
check failed before correcting the expectation. No production result or tolerance
was changed. The check uses `SimplexTree.persistence` with
`min_persistence=0` and `persistence_dim_max=True` on the exact supplied cells.
These are Python binding correctness checks, not native benchmark results or a
GUDHI spectral comparison. Ripser/Topp are outside this supplied-complex
input-reading change; no Rips or distance implementation is changed.

The unchanged Python tool suite was also attempted on Windows: 98 tests ran,
with one POSIX `resource` import error, one forward-slash diagnostic assertion
failure and five platform skips. This is not a local tool-suite pass. The two
failure mechanisms are the platform limitations already recorded in the
[interface compatibility report](../design/interface-compatibility.md); hosted
Linux tool checks remain separate evidence. The original collaboration checkout
also retains five missing links to historical #32 report paths; they are outside
this contribution and do not fail the new main-based documentation checks.
