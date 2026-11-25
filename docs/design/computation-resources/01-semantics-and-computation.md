# Mathematical semantics and executable computation

[Design index](README.md) / [Development plan](development-plan.md)

This design defines the boundary every implementation phase must preserve.
Current behavior is read from [main's exports][exports] and
[kernel design][kernel]. Proposed execution changes do not replace those public
mathematical responsibilities.

## Three contracts at a call boundary

| Contract | Content | Owning side |
| --- | --- | --- |
| Mathematical | Source/run identity, field, degree, basis/order, coverage, operation, precision and required output guarantees | Source, domain request and result types |
| Computational | Concrete implementation steps; immutable reads, mutable writes, preparation, materialization and failure | Algorithm function or bounded local composition |
| Resource | Typed capacity/layout, exclusive use, initialization, checked growth, retention and last use | Caller, terminal or explicit domain scope |

The first-class execution concept is a computation unit with the last two
contracts made explicit. A unit is normally a function or useful algorithm loop.
It is large enough that its caller can reuse resources or avoid repeated work;
it is not a heap object for each scalar addition. A domain execution plan is the
composition of these units, potentially expressed by static calls and branches.

A request for a full spectrum and a request for selected Ritz candidates are
different mathematical capabilities. Choosing between dense and sparse ways to
apply the same declared operator is an implementation choice. Omitting a
required dense output changes the request; avoiding an unrequested dense
intermediate changes execution. Reviewers must make this distinction before
accepting an optimization comparison.

## Concrete lowering examples

For real Hodge action in the stored orthonormal chain basis,
`L_q x = B_q^T(B_q x) + B_(q+1)(B_(q+1)^T x)`. Its computational plan can use a
lower temporary of length `dim C_(q-1)`, an output of length `dim C_q`, and streamed
upper-column dot/accumulation. The mathematical operator does not prescribe a
dense `dim C_q` square matrix or a fresh output allocation per application.

For persistence, the requested intervals do not prescribe one reducer. H0
connectivity, specialized F2 H1, implicit cohomology and explicit boundary
reduction have different access and retention needs. Resource-aware interfaces
must allow those differences. Retaining transformations follows representative
or compatible query demand, not the availability of a large workspace.

For diagram matching, the mathematical request fixes dimension, endpoint and
metric policies. Extracting operand-local points, choosing pair-local scale,
building a graph and solving matching are separate computations. Some prepared
operand data may be reusable; pair-dependent normalization does not become
operand-local by placing it in a cache.

## Domain responsibility map

| Existing owner | Semantic responsibility preserved | Computation/resource boundary to inspect |
| --- | --- | --- |
| `geometry` | Coordinate/distance views, metric and numeric interpretation | Borrowed reads, pair evaluation, optional dense conversion; no mandatory workspace for scalars |
| `complex` | Validated topology, source-local IDs and ordered signed incidence | Reuse frozen guarantees/access without globally exporting internal working indices |
| `filtration` | Exact/approximate construction, scale/coverage and prepared source identity | Direct versus prepared versus explicit construction; consuming callback effects |
| `persistence` | Intervals, selected fields/dimensions and optional witnesses | Preserve specialized algorithms; expose selected working-resource boundaries |
| `algebra` | Prime-field arithmetic and column/reduction invariants | Compatible column operations and optional retained transforms |
| `diagram` | Output identity, dimensions, context, multiplicity and ownership | Borrow common data or move outputs; do not retain producer scratch |
| `descriptors` | Descriptor formula, grid and endpoint meaning | Direct scans or justified local event preparation |
| `diagram_distances` | Matching metric, norm/order, endpoints and context | Operand preparation versus pair graph/solver state |
| `optimization` | Serial F2 critical-set formulas and supplied-source interpretation | Lazy compatible preparation and query-local residuals |
| `homology_operators` | Exact F2 windows, projections, classes/fillings and family/run bindings | Its own sparse/packed RREF and query work; not PH reduction |
| `hodge` | Signed real degree/selection and spectral/result policy | Action output/lower scratch, dense export and iterative caller |
| `dirac` | Truncated graded real operator, offsets and signed modes | Output placement, Gram intermediate and full-operator verification |
| `persistent_laplacian` | Same-source stage pairs, basis and Schur/numerical policy | Main dense work; candidate-only action/factor/sequence boundaries |
| `execution` | Cooperative controls and per-terminal work accounting | Stay separate from scratch ownership and process-memory accounting |
| Private eigensolver | Fixed real symmetric action and approximate solver guarantees | Caller-supplied action output, bounded blocks and projected work |

Existing algorithm functions remain independent of default builders/dispatch.
Source adapters interpret a request and call them; algorithms do not call back
into default selection. Construction does not invoke PH to produce a complex.
Descriptors/matching consume compatible data without rerunning its producer.

## Identity, views and preparation

Distinguish semantic identity from physical storage compatibility. An f64 buffer
can serve unrelated computations after its previous content is dead. A sparse
structural plan depends on its pattern/layout. Numeric factors depend on the
actual data, order and numerical policy; serial critical-set content additionally
depends on its supported decomposition form. Warm vectors can transfer to a
changed operator only through a valid coordinate map and fresh evaluation.

Use existing immutable borrows and owned snapshots as the first validity
mechanism. Where an object already binds a run through `Arc`, preserve that
binding. Do not invent a universal source ID, global content hash or zero-copy
matrix interface merely to pass resources between functions. A view must describe
the access actually needed, not pretend all mathematical objects have one layout.

## Allowed and forbidden dependency changes

An algorithm may borrow typed scratch and call compatible algebra/numerical
helpers. A terminal may assemble a local plan and own preparation. Extract a
shared helper only after two real consumers agree on arithmetic, ordering,
overflow, errors and costs. Helpers that still encode domain rules stay with
their owner.

The first refactor does not change fields, construction guarantees, tie order,
coverage, output completeness, spectral targets or error categories. Safe Rust,
one crate, existing MSRV and no new runtime dependency remain constraints.
Public API changes, if a measured consumer requires one later, get their own
compatibility statement and affected rustdoc/tests; they are not automatic
consequences of this architecture.

## Review and acceptance

A patch must identify the mathematical request, the concrete units being changed,
the resource owner and the affected callers. Tests check the owning mathematical
equations or independent expected results; storage tests check separate resource
effects. Neither substitutes for the other. If a proposed shared object requires
two domains to weaken their guarantees or adopt an unrelated representation,
keep their units separate and share only the compatible primitive.

[exports]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/lib.rs
[kernel]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/docs/design/kernel.md
