# Route integration and compatibility acceptance

[Documentation](../README.md) / Research / [#86](https://github.com/Aequiludium/cocycle-rs/issues/86)

This 2026-10-11 record separates existing pipeline acceptance from future
algorithm integration. The current issue scopes require interface preparation
for #88-#93, while keeping their later running goals open. Contract examples,
hand calculations and ordinary PH are not executions of future algorithms.

## Integrated existing pipeline: #87 and #94

The frozen pre-migration baseline for closed #87 is main
`a7ee9e911dcdeaea137b682bf383f7c2e1d090d5`. The integrated input commit is
`177790cdb123be07e2192f68700b45077d812e99`, on a new feature branch. It contains
the following original PR commits through ordinary merges; their branches and
review gates remain intact. This is a reviewable integrated candidate, not an
assertion that those PRs have merged into main.

| PR / issue | Fixed original head | Actual disposition in this version |
| --- | --- | --- |
| #104 / #72 | `ade7ff83b7e043a8eed072cc6d6e9653aaca0810` | Source-owned Rips/simplicial adapters and one ordinary interval assembly owner; unchanged sealed public paths and mathematical kernels. |
| #105 / #79 | `f7bdf0b4ada8dccb55b3db3a02a2a32fee67edf3` | Retain the established result/read owners; exercise witness-preserving moves and reads after releasing sources. |
| #107 / #82 | `466609721ffd1acde4a792b8df22e4ed54ae95f2` | Retain concrete diagram preparation/solvers and callers; preserve both result and distance disposition tables. |
| #106 / #85 | `0f1e270fdc59e0d5128b85d5c129c88b3f104ea3` | Retain descriptors; include cutoff/multiplicity regression and existing consumer examples. |
| #100 / #69 | `991e867f6897aaa79163c283a2e23551d066ed6a` | Same-owner borrowed stages/inclusions, with the two stage hooks routed to #104's relocated adapters. |
| #102 / #65 | `0cff0f0bbd2daf58d7391f31fe17bef5fdfc6963` | Source/basis construction example; no new topology container. |
| #103 / #71 | `913b4d7c438a113788bdd9694f08a6e51612225a` | Retain the migration map and current-versus-future scope correction. |
| #112 / #73-#76 stage A | `cc1cddeb8175032a4ccca9ae06788dadb54a627a` | Retain proposed family contracts; no spectral, Dirac, projection or generalized solver. |

Shared index and report conflicts preserve all relevant sections. The integrated
version does not adopt #26/#49's pending preparation API, #56's experimental
dispatch, non-diagram comparisons or a learning runtime. Final file contents,
rather than superseded intermediate commits, define the candidate's scope.

### Identity, errors and consumers

`tests/filtered_complex.rs` includes
`staged_results_keep_witnesses_and_consumers_after_source_release`: a signed
triangle with gapped vertex labels, births -2, edges at 0 and filling at 2.
Same-owner inclusions preserve integer incidence. Over F3, analyzing the original
source through stage 0 yields censored H1; explicitly importing that stage as a
supplied ring yields an essential class. Comparing the censored diagram fails.
At stage 2 the owned result has H1 [0,2), and both cycle/cocycle witnesses still
address its ordinal after every input owner and request is released. The cycle's
edge boundary is independently expanded modulo 3 using original vertex labels.
The moved result supports Betti [0,1,1,0], finite mass 2, raw self-distance
0 and comparison to computed-empty H1: bottleneck 1, W1-L-infinity 1,
W2-Euclidean sqrt(2). These are hand-derived expectations.

The complete existing suite also protects implicit point/matrix Rips and flag,
external cells, signed/unequal births, prime fields, H2+ fallback, sparse blockers,
duplicate ordinals, computed-empty/missing dimensions, source coverage and
per-terminal cancellation/budget retry. A failed terminal returns an error,
never a partial diagram/scalar. Reusing `Execution` starts a fresh budget on each
terminal. Context-aware comparison checks field/declared scale; physical units
and truthful dataset declarations remain caller responsibilities. In particular,
contextual self-comparison of this supplied signed source fails because it has
no declared edge-length convention. Explicit raw comparison is a caller choice;
equal operands do not bypass the contextual premise.

### Cost and exit boundaries

The adapter and assembly moves retain call, validation, reservation/clone,
checkpoint and workspace-lifetime order. Their source-normalized bodies match
the frozen baseline; source-specific readers still distinguish trusted frozen
owners from externally validated cells. Already conforming result, distance and
descriptor owners are retained, with no duplicate implementation. The new test
and this report add no runtime allocation, dispatch or dependency to the crate.

Stage selection remains O(log n) time and constant space; same-owner inclusion
validation is constant-time and traversal is O(selected cells), without allocation.
Requested dense exports, matching point buffers, descriptor outputs and opt-in
representative skeleton/transformation storage retain their respective costs.
This is a source-level cost audit, not a new timing/RSS baseline or speedup.
Future changes to actual allocation/routing/workspaces require matched measurements.

Local checks, exact-head hosted native comparisons, review, main integration and
release are separate facts. #94's implementation acceptance is reviewable on this
combined version; its final main exit happens only after approved integration.
Closing #94 does not close the future algorithms in #88-#93, reopen #95-#97 or
certify a release. #57's eventual exit record must retain these boundaries.

### Executed candidate checks

On the integrated candidate plus the new regression, Rust 1.98.1 debug/release
and Rust 1.91 tests each pass 224 tests/doctests, with four opt-in diagnostics
ignored. All-target MSRV, Clippy with warnings denied, crate and standalone
formatting, strict rustdoc, all ten examples, construction/critical-set example
tests in debug/release/MSRV, and all eleven maintained Markdown doctest pages
pass. Package construction and the 188-file package list pass. Source, Markdown,
local collaboration and diff checks pass. The same tree's Linux tool suite
passes all 98 tests. Exact-head hosted/native acceptance is recorded in the PR,
separately from these local checks.

Two initial harness failures are retained in local evidence: the first new test
incorrectly expected contextual distance on unspecified supplied scales, and
the first Linux invocation mishandled its successful test process's exit code.
The test now protects the actual rejection contract; direct WSL execution
returns success. Neither failure required a production-code change.

Reproduce with [CONTRIBUTING's verification commands](../../CONTRIBUTING.md#verification)
on the PR's exact head. The integration-specific command is
`cargo test --locked --test filtered_complex staged_results_keep_witnesses_and_consumers_after_source_release`.
Generated command logs remain under ignored `target/p4-c-spr/`; no raw artifacts
enter the crate or PR. Later stack entries add only interface documentation and
finite hand checks, with their own evidence and no algorithm acceptance claim.

## Single-scale spectral and Dirac handoff: #88

Stage A is an interface contract. The real Laplacian/Dirac producers remain
#73/#74, and new result readers remain #109. No numerical solver, manufactured
spectral result or mandatory persistence call is introduced here. Existing
`SimplicialComplex`, `SimplicialStage`, signed incidence and `FilteredComplex`
stay in their original owners. The runnable
[stage extraction](../guides/interface-contracts.md#read-a-signed-stage-for-a-matrix-consumer)
and construction example exercise the available source-to-basis part.

| Handoff field | Existing read / caller declaration | Future validation and consumer |
| --- | --- | --- |
| Snapshot and selection | Borrow one owner; select all stored topology or an inclusive `stage(scale)`; retain source coverage/skeleton facts | Owned output retains source/selection descriptions; counts or equal `SimplexId` indices do not authenticate another owner. |
| Ordered bases | Restrict source traversal by degree; keep position-to-original-cell keys | A columns and D rows use exactly the same current basis. Empty bases retain 0-by-n/n-by-0 shapes. |
| Domain and incidence | Read integer boundary signs; explicitly choose real/complex coefficients | Check AD=0 in that domain. PH's F2 check does not certify a naive real lift. Complex adjoints conjugate. |
| Inner products and weights | Explicit positive-definite degree Grams, including the choice G=I | Validate shape/finiteness/positivity, units and adjoints. Filtration values and F2 costs are different quantities. |
| Operation and output | Degree, full/lower/upper Hodge choice; Dirac degree range/offsets; numeric policy and requested eigenpairs/vectors | Preserve operator/eigenproblem identity, multiplicity, value/vector pairing and actual residual/orthogonality evidence. Partial eigenpairs do not certify nullity. |
| Costs and ownership | Borrow incidence/bases during extraction; materialize only a requested owned layout | Dense export costs mn+np entries, with independent shapes; sparse consumers can traverse incidence. No source clone, explicit Rips expansion or PH rerun is required by this contract. |

The [family proposals](input-contracts.md#laplacian-input-and-operation-proposal)
define `L=A_adjoint A + D D_adjoint`, weighted adjoints and descriptive calls;
the [Dirac proposal](input-contracts.md#dirac-input-grading-and-truncation-proposal)
defines the graded `K=d+d_adjoint`. Result readers use the
[spectral contract](result-contracts.md#spectral-values-and-vectors), rather than
importing a diagram or treating small residual as a complete eigensystem.

For the guide's signed triangle, the unfilled edge window has
A=[[1,1,0],[-1,0,1],[0,-1,-1]], D:3-by-0. Under standard real Grams its
edge Hodge expectation has eigenvalues (0,3,3), with kernel (1,-1,1).
Adding the face D=(1,-1,1)^T gives 3I and no H1 kernel. An ordinary edge Dirac,
in basis (v0,v1,e01), squares to the ordinary H0/H1 blocks; restricting to C1
instead gives zero Dirac, not the full L1=[2]. These are hand expectations.

The following exact integer checks reproduce the chain/kernel and Dirac-square
claims without constructing a spectral result or invoking an eigensolver:

```python
# Hand contract #88; standard real Grams, ordinary nonaugmented degrees.
a = [[1, 1, 0], [-1, 0, 1], [0, -1, -1]]
z = [1, -1, 1]
lower = [[sum(row[i] * row[j] for row in a) for j in range(3)]
         for i in range(3)]
assert [sum(x * y for x, y in zip(row, z)) for row in a] == [0, 0, 0]
assert [sum(x * y for x, y in zip(row, z)) for row in lower] == [0, 0, 0]
assert [[lower[i][j] + z[i] * z[j] for j in range(3)] for i in range(3)] == [
    [3, 0, 0], [0, 3, 0], [0, 0, 3]]
assert sum(lower[i][i] for i in range(3)) == 6
assert [[sum(lower[i][k] * lower[k][j] for k in range(3))
         for j in range(3)] for i in range(3)] == [
             [3 * x for x in row] for row in lower]
k = [[0, 0, -1], [0, 0, 1], [-1, 1, 0]]
assert [[sum(k[i][r] * k[r][j] for r in range(3)) for j in range(3)]
        for i in range(3)] == [[1, -1, 0], [-1, 1, 0], [0, 0, 2]]
assert k[2][2] == 0  # C1-only Dirac loses the adjacent-degree contribution.
print("#88 signed chain/kernel and Dirac-square hand contracts passed")
```

A computed empty degree uses explicit zero shapes and a valid zero-size Gram;
an omitted required degree is missing input. Nonfinite scales/coefficients,
unidentified basis changes, singular/indefinite Grams and unsupported complex
arithmetic must fail, not yield empty or zero spectra. Unknown source completion
cannot be promoted by a matrix export. Weighted generalized and normalized
operators retain their separate definitions; an ordinary Euclidean solver is
not valid for an arbitrary G-self-adjoint coordinate matrix.

Stage B stays open: after the actual #73/#74/#109 producers exist, run a small
direct input-to-operator-to-reader workflow, compare action and requested
representation, independently check eigenpairs/nullity/Dirac grading, and record
extraction versus solving costs. Current doctests and integer checks establish
only the source/interface contract. The existing implicit PH fast path and its
performance gates are untouched.

## Persistent spectral handoff: #89

Use #69's existing finite same-owner stages and identity inclusions as the
available input, together with #73's spectral conventions. A list of ordinary
stage spectra does not implement the pair's persistent spectrum. #109's future
reader needs both stage identities and the actual persistent operator definition.

The selected first mathematical contract is the real, ordinary nonaugmented
persistent Laplacian for a weighted simplicial inclusion, following Mémoli,
Wan and Wang, [arXiv:2012.02808v3](https://arxiv.org/html/2012.02808v3),
section 2.2 equation (4), section 2.3 Theorem 2.7 and section 5.1. We read these
definitions, including the requirement that source weights restrict from the
target; this is theoretical reading, not execution of their algorithms.
For K contained in L, use the subspace of target (q+1)-chains whose boundaries
land in Cq(K), with the induced target inner product. If E is the restricted
boundary and A is K's degree-q boundary, the operator on Cq(K) is
`E E_adjoint + A_adjoint A`. Its nullity is the real persistent Betti number.
The same-stage case reduces to ordinary Hodge. Independent stage nullities,
row-deleting a target boundary, and F2 Betti numbers do not replace this contract.

| Handoff | Requirement and current boundary |
| --- | --- |
| Pair identity | One owner/snapshot, inclusive scales s<=t, source coverage/skeleton and retained stage-to-owner keys; `inclusion_into` rejects different owners, including equal clones. |
| Degree bases/maps | Ordered source/target bases and signed inclusion in every required adjacent degree; verify both chain-map squares. An ordinary barcode does not supply these maps. |
| Domain/inner product | Explicit real coefficients; positive target degree inner products restrict to source. The restricted upper-domain basis needs its induced Gram, even when original simplex Grams are identity. Arbitrary incompatible Grams are outside the first contract. |
| Operator/result | Degree q, pair (s,t), source basis, restricted domain and adjoints, normalization, arithmetic, requested eigenpairs/vectors and numerical evidence; no single-stage identity substituted. |
| Ownership/cost | Stage/inclusion reads borrow the owner without allocation; finding the restricted subspace and its basis/Gram is new producer work. Its cost is distinct from eigen-solving or copying requested output. |

This executable source check uses only the existing stage API. K has two
vertices; L connects them by a new middle vertex and two edges. No PH or spectral
terminal is called:

```rust
use cocycle::complex::{Simplex, SimplicialComplex};
let source = SimplicialComplex::new(vec![
    Simplex::new(vec![0], 0.)?, Simplex::new(vec![2], 0.)?,
    Simplex::new(vec![1], 1.)?, Simplex::new(vec![0, 1], 1.)?,
    Simplex::new(vec![1, 2], 1.)?,
])?;
let early = source.stage(0.)?;
let later = source.stage(1.)?;
assert_eq!((early.len(), later.len()), (2, 5));
assert_eq!(early.simplices().iter().map(|s| s.vertices()).collect::<Vec<_>>(),
           [vec![2], vec![0]]); // The source's tie order is retained, not resorted.
assert_eq!(early.inclusion_into(&later)?.len(), 2);
assert!(early.inclusion_into(&later)?.all(|(a, b)| a == b));
# Ok::<(), cocycle::Error>(())
```

For the hand calculation, explicitly reorder target vertices as (0,1,2) and
edges as (01,12), with standard Grams. The middle boundary row forces c01=c12.
The upper-domain basis (1,1) has Gram 2, restricted boundary (-1,1), and
persistent L0 is one half of [[1,-1],[-1,1]], with eigenvalues (0,1), nullity 1.
Ordinary K has a zero 2-by-2 Laplacian and nullity 2; ordinary L has spectrum
(0,1,3) on three vertices. Deleting the middle row gives I2, wrongly nullity 0.
The following rational checks verify this example, including its induced Gram;
they do not construct a production persistent Laplacian:

```python
# Hand contract #89; explicitly declared coordinate order and standard Grams.
from fractions import Fraction as F
b = [[-1, 0], [1, -1], [0, 1]]
upper_basis = [1, 1]
assert sum(x * y for x, y in zip(b[1], upper_basis)) == 0
restricted = [sum(x * y for x, y in zip(b[r], upper_basis)) for r in (0, 2)]
gram = sum(x * x for x in upper_basis)
assert restricted == [-1, 1] and gram == 2
persistent = [[F(x * y, gram) for y in restricted] for x in restricted]
assert persistent == [[F(1, 2), F(-1, 2)], [F(-1, 2), F(1, 2)]]
assert [sum(row) for row in persistent] == [0, 0]
assert [row[0] - row[1] for row in persistent] == [1, -1]
naive = [[sum(b[i][k] * b[j][k] for k in range(2)) for j in (0, 2)]
         for i in (0, 2)]
assert naive == [[1, 0], [0, 1]] and naive != persistent
print("#89 restricted-domain Gram and persistent-nullity hand contract passed")
```

Empty source degree is a valid zero-size operator; absent required target degree
is missing input unless the selected finite source explicitly has no such cells.
Unknown coverage, changed weights, reverse scales, cross-owner IDs and incompatible
maps must fail rather than silently compare equal shapes. The same-source
integer inclusion can be reduced to F2, but a real spectral nullity concerns
real homology; torsion may make a prime-field PH result different.

Stage B remains unimplemented. Its future acceptance runs the restricted-domain
construction on an actual producer, compares persistent zero multiplicity to an
independent real boundary-rank/image calculation, checks the same-stage limit
and eigenpair errors, and records source extraction/subspace/solver costs.
The example above fixes a first nontrivial pair without reimplementing that
producer or blocking the existing-code #94 exit.

## F2 class and transport handoff: #90

The existing source producer is #69's same-owner stage/incidence access;
#75's window/legality/query contract is the future projection producer. The
[F2 family proposal](input-contracts.md#f2-window-legality-and-query-proposal)
and [result readouts](result-contracts.md#f2-projections-actions-and-query-records)
own those definitions. This route integrates their identities and maps, not a
second solver, reference snapshot schema or a CriticalSetWorkspace projection.

| Field or transition | Required meaning |
| --- | --- |
| Degree-k window | A:m-by-n, D:n-by-p over F2, AD=0, three ordered degree bases with source/stage keys and explicit empty dimensions. Sum signed integer incidence before reducing modulo 2. |
| Positive costs | Exactly n positive costs with units/arithmetic and a basis-aligned identity; the initial bounded proposal uses checked integers. Costs are not filtration values or real Grams. |
| Stage map | Jprev/Jcur/Jnext with `A_t Jcur=Jprev A_s` and `D_t Jnext=Jcur D_s` over F2; same-owner inclusion supplies coordinate injections, not arbitrary chain maps. |
| Query/result | A usable legal P, same source/basis/cost description, operation/solver evidence and availability; class queries require cycles. `Computed(0)` is a value; unavailable optimality is not zero. |
| Transport | For a source class representative z, `P_t Jcur z` is its chosen target representative. Validate class preservation with a target-boundary preimage; do not demand that arbitrary representative complements commute as matrices. |
| Borrow/move/cost | Borrow stages/incidence while extracting a requested window; an owned result must own enough basis/source descriptions to outlive the source. Window extraction, legality verification, solving, action and boundary-preimage checking have distinct future charges/costs. |

Legality checks all four conditions: P squared equals P, AP=0, PD=0, and
(I+P)z is in image(D) for every cycle z. P=0 passes the first three on a ring
but fails the fourth. Feasible, independently verified legal, and globally
optimal are separate claims. Preserve an independent cycle/boundary-rank check;
do not promote selected mass to minimum-class mass or ordinary PH representatives
to a validated projection family.

The first declared sample is the triangle ring at s=1 and filling at t=2 from
the interface guide. In its edge order (20,30),(10,30),(10,20), the F2 cycle is
z=(1,1,1). Ring D:3-by-0; target D=z. Choose the hand projection P_s whose
first column is z and other columns are zero; choose P_t=0. Jcur and Jprev are
identity, while Jnext maps the zero source face space into the new face space.
The two chain squares hold. Transport sends z to zero, with target preimage 1
because D_t*1=z. Thus the class dies; this is not a transport failure. Costs
(2,3,5) give source selected mass 10 and target selected mass 0. Although both
matrices are legal, `P_t Jcur` does not equal `Jcur P_s`; legal projections choose
complements, not a mandatory commuting matrix family.

This exhaustive eight-chain check validates only the declared hand matrices,
their class preservation and transport witness, without solving for a P:

```python
# Hand contract #90; declared F2 triangle ring/fill, not a solver invocation.
from itertools import product
a = [[1, 1, 0], [1, 0, 1], [0, 1, 1]]
z = (1, 1, 1)
ps = [[1, 0, 0], [1, 0, 0], [1, 0, 0]]
pt = [[0, 0, 0] for _ in range(3)]
def action(matrix, x):
    return tuple(sum(b * c for b, c in zip(row, x)) % 2 for row in matrix)
zero = (0, 0, 0)
cycles = [x for x in product((0, 1), repeat=3) if action(a, x) == zero]
assert cycles == [zero, z]
for p, boundaries in ((ps, (zero,)), (pt, (zero, z))):
    for x in product((0, 1), repeat=3):
        assert action(p, action(p, x)) == action(p, x)
        assert action(a, action(p, x)) == zero
    assert all(action(p, b) == zero for b in boundaries)
    assert all(tuple(b ^ c for b, c in zip(x, action(p, x))) in boundaries
               for x in cycles)
assert action(ps, z) == z and action(pt, z) == zero
assert tuple(b ^ c for b, c in zip(z, action(pt, z))) == z  # D_t * 1.
costs = (2, 3, 5)
assert sum(c * x for c, x in zip(costs, action(ps, z))) == 10
assert sum(c * x for c, x in zip(costs, action(pt, z))) == 0
assert z not in (zero,)  # False zero-P on the ring fails homology preservation.
print("#90 F2 legality, class death and boundary transport hand contracts passed")
```

Wrong order/shape/domain, a noncycle class query, nonpositive costs, stale cost
identity, cross-owner stages or an invalid chain square must fail. Explicit
empty current degree n=0 permits a zero-size window and no costs; missing a
required degree is a different state. Topology-equivalent or equal-sized
windows do not identify their coordinates. An explicit basis permutation must
transform A/D, all three maps, costs, actions and query keys consistently.

Stage B remains future work: run actual single-scale and two-stage components
once #75 supplies validated results, independently replay legality and boundary
transport witnesses, check query/source associations and arithmetic errors, and
separate feasible/optimal evidence and extraction/solve/action costs. The local
Rust experiment or Python reference is not adopted into this production stack.
