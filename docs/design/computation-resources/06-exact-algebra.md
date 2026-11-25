# Exact algebra, persistence and critical queries

[Design index](README.md) / [Preparation rules](04-preparation-and-sequences.md)

This is the R5 domain design. Numerical buffer work does not authorize replacing
exact-domain representations or forcing specialized PH through a common reducer.
The change must have a concrete repeated-query or multi-output consumer.

## Existing algebraic boundaries

| Owner | Observed contract | Consequence |
| --- | --- | --- |
| [Boundary reduction][boundary] | Consumes ordered prime-field columns; `reduce` retains transforms and `reduce_pairs` omits them | Preserve output-driven retention; columns and transformed content are not free scratch |
| [Flag dispatch][dispatch] | Chooses specialized H0/F2 H1 or compatible higher cohomology/fallback | No uniform-boundary rerouting merely to share resources |
| [Critical sets][critical] | Entire supplied complex, serial F2 ELZ/lazy support, primal R/V, lazy dual/transposes and bounded U queries | Sharing requires the same full ordered input and form, not just F2 and dimensions |
| [F2 matrices][matrix] | Domain-private sparse/packed RREF and exact generalized-inverse preparation | Keep independent of PH pivot/filtration semantics |
| [Finite operator family][family] | Homology-coordinate maps; two temporary stage bases; owned histories | Reuse temporary work without deleting requested history output or run binding |

## Define the consumer before retaining content

| Requested work | Required retained content | Release/temporary boundary |
| --- | --- | --- |
| Diagram only | Pair information needed by the chosen route | Do not retain maximum transforms for unspecified future queries |
| Diagram and representatives | Compatible transform columns and interval associations | Witnesses become owned output; unrelated reduction work can drop |
| Repeated critical queries | Serial primal R/V and completed lazy per-dimension state | Query residuals, bounded U solve and returned sets remain separate |
| Repeated exact class/filling queries | Current F2 window/run and its actual projection/RREF preparation | RHS and search state remain query-local unless explicitly reusable |
| Family barcode/features | Adjacent maps and current temporary bases; full histories only when requested | Preserve finite-stage interval semantics and deliberate output size |

Keep the current const-generic retention branch in the boundary reducer unless a
measured consumer needs another real capability. Do not add packed PH columns,
chunk/twist backends or clearing options as placeholders. Existing H1 heap reuse
is useful evidence of local resource control, not a reason to re-express H1 as
explicit boundary reduction.

## Proposed repeated-query computation contract

A bounded critical-query unit reads immutable source and valid primal/dual
content, borrows exclusive residual/working capacity, grows sparse fill-in
fallibly where supported, and produces a complete result or an error. Lazy
preparation uses the query's cumulative budget and publishes only complete
state. Query scratch never becomes a valid decomposition or silently survives
as a semantic answer for a different endpoint/target.

For an actual combined serial consumer, consider a private borrow/move boundary
for existing `BoundaryReduction` plus its source/form guarantees. Do not export
a generic `PreparedReduction` just to name that possibility. First prove that
both consumers use the same complete ordered matrix, field, transformations and
support guarantees. A selected skeleton, cutoff or implicit source expansion
can make apparently similar matrices different inputs.

| Compatibility check | Failure response |
| --- | --- |
| Same ordered full input and source interpretation | Otherwise prepare through the current consumer's route |
| Same arithmetic and retained transform meaning | No cross-field or numeric rank substitution |
| Required serial support form | Reject/avoid cleared or parallel imported content |
| Actual output capabilities retained | Recompute missing content through a valid route; do not invent witnesses |
| Existing source/run identity preserved | Do not equate same numeric cell indices across sources/runs |

Preparation sharing may reduce repeated work but increases retained memory.
Measure source binding, conversion, transforms, query work and drop together.
If the consumer only needs a diagram, the current pair-only route may remain
cheaper and should remain its default.

## Representation-specific scratch

`Column` uses ordered map storage; clearing it does not offer the same reusable
node-capacity semantics as a dense vector. Do not promise allocation elimination
through a universal scratch pool. Profile traversal, insert/remove, fill-in and
retained columns before changing the representation.

F2 operator storage already adapts between sparse and packed forms. A reusable
RHS or product buffer must state binary shape, word padding, zero initialization
and exclusive lifetime; it does not import PH ordering rules. Preserve adaptive
transition thresholds until a separate representation study supports a change.
The operator family's two-stage-basis lifetime is already deliberate; identify
remaining work before proposing another persistent family cache.

## Failure and concurrency

Critical queries currently borrow `&mut self` for lazy state. Keep exclusive
query preparation unless a concrete concurrent consumer justifies splitting
immutable factors from per-query scratch. Parallelism is not introduced by this
resource refactor. A failed lazy transpose does not erase previously complete
primal/dual state, and no incomplete result is returned as success.

Exact class/filling/affine search can remain exponential. Reusing a candidate or
residual buffer changes allocation, not the mathematical search bound. Preserve
existing cooperative limits, exact equations, source identities and complete
failure behavior.

## Independent checks and acceptance

Verify `D*V=R`, normalized pivots, pair/zero-bar associations, representative
cycle/boundary equations and serial support properties as applicable. Use
independent tiny exact reduction/enumeration or hand-derived expectations rather
than the changed reducer as its own oracle. Keep GUDHI plus supported Ripser
checks for affected persistence/VR paths under the current validation protocol.

Test odd primes and specialized PH controls even if the changed consumer is F2,
so a shared helper cannot silently narrow the existing public domain. Include
tie/order changes, sparse/packed boundaries, failed-then-valid lazy queries and
retained query/results after scratch reuse. Compare cold preparation, warm query,
multiple requested dimensions and retained high water.

R5 exits with either a small admitted owner-local change for a named consumer or
an evidence-backed no-adoption decision. It does not exit by installing a generic
carrier that no production caller consumes. The previously withdrawn
uniform-boundary route remains a regression constraint.

[boundary]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/algebra/reduction/boundary.rs
[dispatch]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/persistence/flag/dispatch.rs
[critical]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/optimization/mod.rs
[matrix]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/homology_operators/matrix.rs
[family]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/homology_operators/family.rs
