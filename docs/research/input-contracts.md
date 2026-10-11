# Shared input contracts for persistence, spectra and F2 operators

[Documentation](../README.md) / Research / [#60](https://github.com/Aequiludium/cocycle-rs/issues/60)

Technical report, 2026-10-09, following the [interface inventory](interface-inventory.md).
This publishes the full input decision behind the
[interface compatibility report](../design/interface-compatibility.md).
The observations below refer to the fixed revisions in the source table;
publication on a newer main does not retest those external implementations.

## Decision and scope

Retain the existing borrowed geometry/complex access and the four-method
`FilteredComplex` contract for prime-field persistence. For a concrete new
consumer, start with one operation over an existing borrowed owner and its
required ordered bases and mathematical parameters. Extract a small owned chain
window only when a matrix-based solver or external representation needs it.
Do not add a universal source trait, force a filtration onto single-scale
computations, or require an inner product for F2 persistence.

This defines input responsibilities and validates finite examples; it adds no
production spectral engine, homology operator, stage-view API, chain-map type,
identity service or public Rust signature. Public operation/result signatures
are addressed by [result contracts](result-contracts.md), the #62 reuse decision
and the #63 compatibility gate. Later mathematical implementations remain separate.

## Fixed sources and actual reading

| Source | Version and files read | Adopted lesson / limit |
| --- | --- | --- |
| Cocycle | main `9e6715f4c2e0118ab738486bd747ec3e361397de`; [filtered contract][cells], [simplex/order][simplex], [storage][storage], [incidence][incidence], [external reader][reader], [mathematics][math], [kernel design][kernel] | Reuse frozen signed incidence and stable owner-local handles. A selected-field boundary-square check does not establish integer or real validity for a custom adapter. |
| homology-operator | `a8d03d1699307ce1997b3c7e81a1946bb83f5408`; [ChainWindow][window], [input identities][identities], [family validation/inclusion][family], [interface][operator-interface], [architecture][operator-architecture], [result contract][operator-results] | A/D are explicit F2 matrices with ordered named bases and positive coordinate costs. Actual family support is nested coordinate inclusions, not arbitrary chain maps. Preserve its limits and distinguish feasibility from optimality. |
| TopoNetX | `ce3d414b571fff95916c4cbc20b8d2c2c9d428e6`; [simplicial source][toponetx], [cell source][toponetx-cell], fixed README | Read implementation rather than inferring capabilities from names. Signed boundary matrices carry row/column index maps; the studied Hodge path uses ordinary transposes. No weighted-Hodge capability is inferred from the `weight` parameter. |

The main source and homology reference were extracted into isolated ignored
directories. The user's external reference checkout was only read and remained
clean. The TopoNetX source files were fetched at the fixed SHA, with blob and
SHA256 identities recorded locally; TopoNetX was not installed or executed.

Two source details prevent blind adoption of that spectral API. Its
`incidence_matrix(0)` is a one-row augmentation, whereas Cocycle's ordinary
nonaugmented H0 uses `C[-1] = 0`. TopoNetX's Hodge rank-zero branch instead uses
`B1 B1^T`; these must not be conflated. Also, the studied simplicial
`incidence_matrix` accepts `weight` but does not read it in its implementation;
`hodge_laplacian_matrix` forwards it to that path, while up/down Laplacians
explicitly reject a non-None weight. This is source inspection at this revision,
not a runtime compatibility finding or a claim about other revisions/domains.

## What the three consumers actually read

| Consumer | Required data | Additional semantics | Data it does not require |
| --- | --- | --- | --- |
| Ordinary persistence | Ordered filtered cells/cofaces, dimensions, finite values, signed incidence reduced into the chosen validated prime field, requested dimensions/range | Face-before-coface order, tie rule, coverage, source/skeleton certificates and optional simplicial witness requests | Real inner products, geometric costs and an eigen-solver |
| Single-scale Hodge/spectral computation | Explicit stage topology, ordered bases of adjacent degrees, signed real/complex boundary maps | A declared inner product/Gram matrix in each required degree, adjoint convention, arithmetic and numerical policy | A persistence diagram or invented filtration values; F2 projection weights are not implied |
| F2 homology operator | Degree k; A:m-by-n and D:n-by-p over F2; ordered bases of lengths m/n/p; AD=0 | n positive coordinate costs with semantics, units and arithmetic; the solver later establishes a legal P and separate certification | Point coordinates, a Rips constructor, a PH diagram or a real adjoint |

Dirac will need compatible graded boundary and adjoint data, rather than the F2
operator `I+P`. Persistent spectral constructions need their own definition in
addition to stages/maps; ordinary stage spectra are insufficient. Mayer
`boundary^N=0`, sheaf restrictions and path/hypergraph boundaries are separate
theory-specific inputs, not silent weakening of `boundary^2=0`.

## Numbering, bases and orientation

1. Bind an input to its source snapshot/owner. `SimplexId::index()` is the
   filtration position in one frozen `SimplicialComplex`, not a cross-complex or
   original-vertex ID. Custom `FilteredComplex::CellId` is opaque, Copy/Eq/Hash
   and stable during the operation; serialization or integer encodings are not
   required by the existing trait.
2. A degree basis is an explicitly ordered sequence of selected cells. For the
   initial adapter, preserve the source traversal restricted to that degree.
   Keep the map from each basis position back to its source cell. Dimensions
   alone do not identify a basis. Do not silently sort into lexicographic order.
3. Keep the concrete simplex orientation: increasing original vertex labels and
   coefficient `(-1)^i` when omitting vertex i. A different orientation or basis
   order requires explicit signed permutation maps in every affected matrix.
   F2 forgets the signs, but that is a coefficient conversion, not a change to
   the stored oriented source.
4. Use A rows/current columns and D current rows/next columns in those declared
   bases. Preserve all shapes, including 0-by-n, n-by-0 and empty matrices. At
   k=0 the ordinary convention is an empty previous basis and A:0-by-n.
5. A scale selection includes every cell with value <= scale. Ties enter
   together and faces remain selected. Keep stage index distinct from scale;
   repeated scales may denote different stages in an explicitly ordered family.
   Implicit Rips/blocker source access remains private and retains its own rules.

For a direct graph consumer, the stored graph may be interpreted explicitly as
a one-dimensional cell complex: vertices and edges only. It must not silently
become its clique/flag complex. A graph triangle and its filled clique complex
have different chain spaces. Borrow the graph and derive its oriented incidence
when required; no dummy persistence request is needed.

## Coefficients, weights and source identity

Integer simplicial incidence can be reduced modulo a chosen prime, or explicitly
lifted into real/complex arithmetic with its original signs. A modular boundary
matrix is not a field embedding into the reals. The small window

```text
A = [1 1 0], D = [1 1 0]^T
```

satisfies AD=0 over F2 but its naive real lift has AD=2. An adapter containing
only modular coefficients cannot manufacture the missing integer orientation or
claim real chain validity. For arbitrary supplied integer cells, the source owns
the integer chain contract and a real consumer must validate the contract it
actually needs; passing PH's field-specific ingestion is not that certificate.

Use distinct parameters for filtration values, spectral inner products and F2
class costs. Their numerical equality does not establish equal meaning. The
single-scale graph probe stores edge values 7/11 but uses a separately declared
current-degree Gram diagonal 2/3. The triangle F2 cost vector 2/3/5 is an
abstract positive cost in its declared edge basis, not its filtration values.

For the spectral window, with Gram matrices Gprev/Gcur/Gnext, the adjoints are

```text
A_adjoint = inverse(Gcur) conjugate_transpose(A) Gprev
D_adjoint = inverse(Gnext) conjugate_transpose(D) Gcur
Lcur = A_adjoint A + D D_adjoint
```

Standard inner products are an explicit choice G=I. Real symmetric/complex
Hermitian positive-definite Grams and numerical validation are spectral
responsibilities. A weighted L need not be Euclidean-symmetric; Gcur L is
Hermitian. Do not pass it to a Euclidean-symmetric solver without the appropriate
coordinate transformation. Floating residual/tolerance/eigen-solver policies
remain in #73; the probe uses exact rational calculations and no eigen-solver.

F2 costs accept the consumer's documented arithmetic: exact integer/rational or
explicit floating positive finite values. They need one value per current basis
coordinate, declared units and weight semantics. Changing order requires moving
the costs with their coordinates; changing costs invalidates geometry-dependent
results, even if topology is unchanged. Empty current spaces have empty costs;
zero/negative/nonfinite costs and accidental bool coercion are invalid.

The integration contract must bind source namespace/snapshot, degree, field,
ordered bases/orientation, weights/arithmetic/units and actual source provenance.
Runtime owner binding may suffice for a borrowed operation; persist explicit
identities when crossing snapshots, storage or cache boundaries. A global hash
service/new public identity type is not required for this research decision.

The reference `input_id` includes k/A/D/bases but excludes `source_metadata`;
`weight_id` excludes the basis and is checked alongside `basis_id`. The probe
confirms that changing only source metadata does not alter these content hashes.
Do not copy that scheme and treat content equality as proof of the same dataset.
Namespaced stable labels and retained source facts are needed for external
adaptation. None of these identity enforcement rules has been added to Cocycle
by this report.

## Stage maps

For a map from source stage s to target stage t, explicitly identify both source
snapshots, degree bases and coefficient field. A degree-preserving chain map has
three compatible matrices Fprev/Fcur/Fnext satisfying

```text
A_t Fcur = Fprev A_s
Fcur D_s = D_t Fnext
```

Check shapes and both equations in the declared field, not just equal labels,
matching Betti numbers or a same-sized matrix. A filtration-compatible map also
respects declared stage/range/units and, for a single scale-preserving filtered
map, each nonzero image term's target value is <= its source value. Physical unit
compatibility is caller-established, not inferred from equal numeric scales.

Nested stages of one source can construct coordinate inclusions from stable cell
keys and explicitly handle reordered bases. General linear maps need their own
explicit matrices; the studied `OperatorFamily` supports coordinate inclusions
only. Algebraic chain compatibility does not imply an isometry or preserve class
mass. Inherited costs keep shared-coordinate values/arithmetic; explicit variable
costs need preserved semantics/units and separate geometry provenance. Neither
case alone defines a persistent Laplacian.

Maps compose with compatible middle identities; verify Fss=I and
Ftu Fst=Fsu where required. The finite operator probe additionally checks the
induced kernel-coordinate transport identity and three-stage composition. This
is distinct from proving a general implementation or a persistent spectrum.

## Concrete finite examples and outcomes

The Rust probe only calls main's existing public graph, frozen-complex, incidence
and persistence APIs. It borrows the owners and creates small owned exported
basis/matrix data. The Python probe imports the isolated pinned reference source
without installing it or a native extension.

| Example | Actual ordered data / expectation | Observed result |
| --- | --- | --- |
| Direct single-scale path graph | V=(0,1,2), E=((0,1),(1,2)); A columns (-1,1,0) and (0,-1,1), D:2-by-0 | A^T A=[[2,-1],[-1,2]], eigenvalues 1/3 by hand; no filtration/PH dependency in this direct calculation. Weighted Gcur L=[[3,-2],[-2,5]] is symmetric positive definite. |
| Triangle boundary stage | C0=((30),(20),(10)); C1=((20,30),(10,30),(10,20)); C2 empty. A=[[1,1,0],[-1,0,1],[0,-1,-1]] | Main cap-1 PH produces H1 birth 1 with RightCensored through 1; the stage itself has beta1=1. Standard real L1 has eigenvalues 0/3/3; F2 operator preserves cycle (1,1,1), selected mass 10. |
| Filled triangle stage | Same C0/C1, C2=((10,20,30)); D=(1,-1,1)^T | Integer AD=0. Main cap-2 PH yields H1 [1,2). Real L1=3I; F2 beta1=0 and selected cycle mass 0. The coordinate chain inclusion induces a zero class transport. |
| Abstract F2 window | A=[1,1,0], D=(1,1,0)^T; labelled dimensions 1/3/1, costs (2,3,5) | Valid without inventing simplices; beta1=1, class (0,0,1) has selected mass 5. Naive real lift fails the chain equation. |
| Empty / ordinary H0 | Explicit empty 0-by-0 maps; separately A:0-by-2 and D=(1,1)^T | Empty operator beta=0; two connected vertices have ordinary H0 beta=1. Empty chain spaces remain valid, not absent computations. |
| Reorder / change costs | Reverse edge basis, A columns and attached costs; separately change one shared cost | Explicit permutation inclusion preserves the class. Basis/input identities change on reorder; weight identity changes on cost change. Inherited variable costs reject; explicit Variable policy accepts. |
| Invalid input/maps | Duplicate basis IDs, incompatible shapes, AD!=0, zero/negative/bool/nonfinite/wrong-length costs, noncycles, broken inclusion equations, nonfinite/decreasing scales | Rejected by the corresponding fixed reference contract; duplicate scales retain separate stage identities and class transport. |

The unfilled stage's intrinsic beta/Hodge calculation does not upgrade the
original cap-1 persistence query's censored class into a known essential class.
Finite explicit sources also do not certify completeness of an omitted larger
Rips skeleton. Source coverage must remain attached to the original construction.

Validation on Windows, Rust 1.98.1 and Python 3.10.11: the public Rust probe
compiled and ran; 47 checked conditions passed, including expected rejection
conditions, hand-derived matrices, reference operator identities and transport
composition. The 47 conditions are not 47 independent fixtures or general proofs.
Source inspection of TopoNetX used Python 3.12 because its class syntax requires
that interpreter. The Python 3.10 inspection attempt's SyntaxError was a tooling
mismatch, not a library defect.

An optional GUDHI Python import in the existing external environment failed with
ModuleNotFoundError; no package was installed there. No GUDHI/Ripser/Topp, native
protocol, TopoNetX runtime, performance, MSRV or hosted CI acceptance is claimed.
The original probes changed no production code. Future production mathematical
changes still require the applicable independent/native evidence.

## Reuse versus new access: selected and rejected alternatives

| Option | Cost and decision |
| --- | --- |
| Existing borrowed graph/complex, restricted ordered basis and consumer-specific parameters | Selected for concrete first consumers. Reuses validation/incidence, source lifetime and current APIs; single-scale calculations can ignore filtration values without manufacturing a PH request. |
| Explicit owned A/D window for a matrix/external consumer | Selected only where demanded. Sparse extraction needs selected IDs/index maps and incidence; a dense m-by-n/n-by-p export has O(mn+np) storage. The probe's dense export is a tiny adapter, not the recommended high-dimensional production representation. |
| Enlarge FilteredComplex with weights, adjoints, stage maps and algorithm hooks | Rejected now: irrelevant obligations for persistence, unclear fields/ownership and no implemented consumer benefit. Keep its current static-generic role. |
| Public nonfiltered chain-view/window/map trait family | Deferred until an actual in-crate consumer establishes the narrow operation and #61 defines its result. No speculative trait or six-crate split. |
| Convert F2 rows or unsigned adjacency into real boundaries | Rejected: the explicit counterexample loses the chain equation/orientation. |
| Use SimplexId/index, equal matrix dimensions or context source-kind as a cross-stage identity | Rejected: local handles and metadata counts do not establish a common basis or source snapshot. |
| Force every single-scale input through a filtration, diagram or explicit flag expansion | Rejected: changes the input question or creates unnecessary materialization. The graph-as-1-complex choice must remain explicit. |
| Implicit weight inheritance, unit inference or silently accepting unknown coverage | Rejected: conflates algebra/geometry/filtration and upgrades unsupported guarantees. |

For #62, begin in existing owners with a concrete function, read-only owner
arguments, explicit degree/stage and only the extra borrowed parameters its
mathematics uses. Allocation, cancellation, empty shapes and source identity
checks belong to that operation. Reuse signed boundaries where valid; retain
specialized implicit persistence paths. Do not expose the probe's export format
as an approved public interchange schema.

## Accepted decision and reproducible first examples

- Three calculation families have explicit read requirements and finite examples.
- Cell identity, degree basis order, coefficients/orientation, weights/Grams,
  stage maps, provenance and borrowing/materialization responsibilities are defined.
- Legal, illegal, empty and unsupported cases are explicit; F2 persistence has no
  mandatory real inner-product requirement.
- A justified minimal reuse decision and rejected/deferred alternatives are recorded.

The #60 research decision is complete. Result interpretation is supplied by
[#61](result-contracts.md), and #62/#63 reused the existing API with an executable
guide. General stage access/maps continue under #68/#69 and consumer mathematics
under #73 (Laplacian/spectra), #74 (Dirac), #75 (F2 homology) and #76
(generalized topology). #72 adapts existing PH interfaces; #96 owns adoption
decisions. Those algorithms are not completed by this report.

The original 47-condition probe and raw logs were local research artifacts,
not a hosted acceptance suite. Their results are reported above without claiming
that this repository contains a public archive of every original run. The
following self-contained Python example reproduces the core hand-derived input
checks using only the standard library. Save it as `target/input-contract-example.py`
and run `python3 target/input-contract-example.py` from the repository root.

```python
from fractions import Fraction

def multiply(a, b, inner, columns):
    return [[sum(a[i][k] * b[k][j] for k in range(inner))
             for j in range(columns)] for i in range(len(a))]

def transpose(a, columns):
    return [[a[i][j] for i in range(len(a))] for j in range(columns)]

a = [[1, 1, 0], [-1, 0, 1], [0, -1, -1]]
d = [[1], [-1], [1]]
assert multiply(a, d, 3, 1) == [[0], [0], [0]]
lower = multiply(transpose(a, 3), a, 3, 3)
assert multiply(lower, d, 3, 1) == [[0], [0], [0]]
assert multiply(lower, lower, 3, 3) == [[3 * x for x in row] for row in lower]
assert sum(lower[i][i] for i in range(3)) == 6  # Eigenvalues 0, 3, 3.
upper = multiply(d, transpose(d, 1), 1, 3)
assert [[x + y for x, y in zip(row, other)]
        for row, other in zip(lower, upper)] == [[3, 0, 0], [0, 3, 0], [0, 0, 3]]

path = [[-1, 0], [1, -1], [0, 1]]
assert multiply(transpose(path, 2), path, 3, 2) == [[2, -1], [-1, 2]]
gprev, gcur = [1, 2, 3], [2, 3]
weighted = [[sum(Fraction(path[r][i] * gprev[r] * path[r][j], gcur[i])
                 for r in range(3)) for j in range(2)] for i in range(2)]
assert [[gcur[i] * weighted[i][j] for j in range(2)]
        for i in range(2)] == [[3, -2], [-2, 5]]
real_product = multiply([[1, 1, 0]], [[1], [1], [0]], 3, 1)
assert real_product == [[2]] and real_product[0][0] % 2 == 0
print("Signed stages, Hodge windows, Grams and the F2 lift counterexample passed.")
```

For actual Rust owner/basis extraction, run the three doctests in the
[interface guide](../guides/interface-contracts.md) with its documented maintainer
command. They exercise existing public APIs, including the signed triangle and
zero-column D shape. Their acceptance is distinct from this rational calculation.

The pinned homology reference can reproduce the ring/fill action without a
native extension. Use a fresh ignored directory, preserving any prior checkout:

```sh
git clone https://github.com/proffitteoy/homology-operator.git target/input-contract-reference
git -C target/input-contract-reference checkout --detach a8d03d1699307ce1997b3c7e81a1946bb83f5408
PYTHONPATH=target/input-contract-reference/src python3 target/operator-contract-example.py
```

Save the Python block under [result contracts](result-contracts.md#reproduce-the-reference-action)
as `target/operator-contract-example.py`. It builds these exact A/D windows,
checks selected masses 10 and 0, and demonstrates an unavailable minimum.
It reproduces representative decisions, not all 47 original conditions.

## Computation-family preparation on main a7ee9e9

This 2026-10-11 supplement supplies stage A of
[#73](https://github.com/Aequiludium/cocycle-rs/issues/73),
[#74](https://github.com/Aequiludium/cocycle-rs/issues/74),
[#75](https://github.com/Aequiludium/cocycle-rs/issues/75) and
[#76](https://github.com/Aequiludium/cocycle-rs/issues/76), under
[#70](https://github.com/Aequiludium/cocycle-rs/issues/70).
The source baseline is main `a7ee9e911dcdeaea137b682bf383f7c2e1d090d5`.
These are proposed operation contracts and independent hand expectations;
no new algorithm, public type, trait, numerical backend or runtime dependency is
implemented. Stage B of #73–#75 remains future algorithm work. The #76 decision
below selects a future research example without approving its implementation.

### Source reading and compatibility disposition

The current `src/complex/filtered.rs`, `src/complex/simplicial/{mod,simplex,incidence}.rs`,
`src/persistence/filtered.rs`, `src/execution/mod.rs`, interface guide and existing
construction example were read at that main. TopoNetX's fixed simplicial source
was reread at `ce3d414b571fff95916c4cbc20b8d2c2c9d428e6`, blob
`8ace5fef7ad5e9a44cce02b284578444340bbae9`: `incidence_matrix`,
`hodge_laplacian_matrix`, `dirac_operator_matrix` and hypergraph closure.
The homology reference's `docs/INTERFACE.md` was reread at
`a8d03d1699307ce1997b3c7e81a1946bb83f5408`, blob
`9a5423763e4e7b5ebb8b2ab10063001a3d0f573a`.
This is source inspection, not a new execution of either external package.

| Existing owner | Retain now | Future consumer's additional responsibility |
| --- | --- | --- |
| Frozen simplicial owner and integer incidence | Increasing-vertex signs, traversal order and owner-local IDs; no change to storage or `FilteredComplex` | Select degree bases, preserve position-to-cell maps and establish the actual coefficient-domain chain equation. |
| Generic filtered PH reader and implicit flag/Rips paths | Field-specific reduction and current dispatch; do not route spectra through PH or materialize implicit cliques | A concrete spectral/projector consumer extracts only the required stage/window, independently of persistence requests. |
| `Execution` and errors | Fresh terminal budget, cooperative cancellation, fallible allocation; no memory/thread quota | Define charge units and checkpoints for extraction, validation, matrix action, numerical solving and certificate replay before implementing the operation. |
| Diagram, representatives and `CriticalSetWorkspace` | Existing result/witness and optional critical-set ownership | Spectrum and homology-projector results need their own owners; critical sets do not supply a legal P. |
| Pending PRs #100/#102/#103/#104 | Preserve their stage, modeling and PH-migration scope; no branch update here | Use a merged version or an explicit future stack if required; these contracts use main's existing access and do not require those PRs to run. |

For each extraction, traverse `cells()` once, select every cell with value <= the
declared scale (or all cells for an explicitly unfiltered use), then restrict that
order by degree. Validate unique keys, faces, dimensions and summed integer
incidences; duplicate terms sum before coefficient conversion, with checked
arithmetic. Supply explicit shapes when either basis is empty. Borrow the source
during the call; an owned output retains owned basis/source descriptions rather
than dangling `SimplexId` interpretations. External identities remain caller
declarations until checked. No global identity service is required.

### Laplacian input and operation proposal

Use the three-degree window `C[k-1] <-A- C[k] <-D- C[k+1]` above, over the
declared real or complex scalar domain, with `AD=0` and an explicitly ordinary
nonaugmented H0. Signed integer simplicial data can feed either domain.
Arbitrary custom sources need their integer/domain validity established;
PH's selected-field check is insufficient. Reject a missing adjacent degree
needed by the requested operator; an explicitly absent degree in the selected
finite complex is a valid zero-size map. A supplied k-skeleton can intentionally
define a lower-only problem, but cannot claim the full Hodge operator of a
larger complex. Source cutoff and skeleton information stay attached.

Require finite Hermitian positive-definite `Gprev`, `Gcur`, `Gnext`, including an
explicit standard-inner-product choice. A zero-dimensional Gram is valid.
Singular/indefinite Grams, nonfinite coefficients, incompatible bases/shapes and
unknown weight semantics fail validation. Complex conjugation is mandatory;
ordinary transpose is valid only for real inputs. A future real-only backend
must report complex input unsupported rather than discard imaginary parts.

The proposed calls, in descriptive notation rather than compilable/public Rust,
are:

```text
laplacian_action(source_or_window, degree, selection, grams, x, execution)
    -> Result<owned coordinates L x with basis and operator identity>
laplacian_matrix(source_or_window, degree, selection, grams, layout, execution)
    -> Result<owned shaped matrix with the same identity>
laplacian_spectrum(source_or_window, degree, selection, grams,
                   eigenpair_selection, vectors_requested, numeric_policy, execution)
    -> Result<owned spectral result described in result-contracts.md>
```

Use `L=A_adjoint A + D D_adjoint`, with the adjoints defined above. Distinguish
full, lower-only and upper-only operators in identity and output. If
`Gcur=R* R`, the Hermitian coordinate problem is `H=R L inverse(R)`;
equivalently `(Gcur L)v=lambda Gcur v`. Applying a generalized solver to
`L v=lambda diag(L) v` is a different problem. Do not form inverses unnecessarily
in an eventual implementation; stable solves/factorization need their own review.
Unit declarations must make both summands commensurate; the request supplies
Gram/weight units and resulting eigenvalue units rather than inheriting edge
length units. Normalized, persistent and indefinite-inner-product variants are
separate definitions, presently unsupported by this proposal.

The [existing path and triangle hand matrices](#concrete-finite-examples-and-outcomes)
are the first fixtures. For the path with `Gprev=diag(1,2,3)`,
`Gcur=diag(2,3)` and empty next degree, `L=[[3/2,-1],[-2/3,5/3]]` and its
hand eigenvalues are `(19-sqrt(97))/12` and `(19+sqrt(97))/12`.
It is self-adjoint in G, not Euclidean symmetric. For the unfilled signed
triangle, `z=(1,-1,1)` spans the real kernel; filling adds `D=z`, making `L=3I`.
Its F2 cycle is `(1,1,1)`, which is not the same real coordinate vector.
The ordinary two-vertex edge has `L0=[[1,-1],[-1,1]]` and nullity 1.
These are expectations, not a delivered eigensolver or nullity certificate.

### Dirac input, grading and truncation proposal

Use a finite ordered direct sum `C[a] + ... + C[b]`, with explicit degree
offsets, adjacent maps `B[j]:C[j]->C[j-1]`, `B[j-1]B[j]=0` in the declared
real/complex field, and block Gram `G=diag(G[a],...,G[b])`. Empty degrees retain
their offsets. For the first proposed ordinary full complex, `a=0` and
`C[-1]=0`; neither augmentation nor cross-degree inner-product terms are implied.

Write `d` for the degree-lowering block operator and `K=d+d_adjoint` for Dirac,
to avoid confusing it with the next-boundary matrix D. The proposed operations
are `dirac_action(graded_source, selection, grams, x, execution)`,
`dirac_matrix(graded_source, selection, grams, layout, execution)` and
`dirac_spectrum(graded_source, selection, grams, eigenpair_selection,
vectors_requested, numeric_policy, execution)`, returning the corresponding
owned action/matrix/spectral information. These names are text proposals.
Each coordinate/result retains degree and local basis position; a spectral
vector generally mixes degrees and must not be assigned a single PH degree.

With compatible adjoints and `d^2=0`, `K^2` has the Hodge blocks
`B[j]_adjoint B[j] + B[j+1] B[j+1]_adjoint`. This identity applies to the same
chosen graded complex. Cutting to degrees [a,b] sets maps crossing its endpoints
to zero: the bottom block loses its lower term and the top loses its upper term.
Squaring a truncated K does not recover the full source's endpoint Laplacians.
Any future alternative boundary condition must be named and validated.
Compatible units/scalings are needed for the summed block operator and its
eigenvalues; degree-dependent units cannot silently be added.

For one ordinary edge and standard Grams, order `(v0,v1,e01)`:

```text
K   = [[0,0,-1], [0,0,1], [-1,1,0]]
K^2 = [[1,-1,0], [-1,1,0], [0,0,2]]
expected spectrum = (-sqrt(2), 0, sqrt(2)); kernel = span((1,1,0))
```

Selecting only `C1` instead gives K=0, whereas full `L1=[2]`; this detects the
truncation error. In the full ordinary finite case, grading parity anticommutes
with K, so nonzero eigenvalues occur in +/- pairs. That check requires the same
unmodified block definition; mass terms/other Dirac variants are not covered.
Laplacian and Dirac share mathematical conventions, without requiring one
solver's private implementation as a dependency of the other.

### F2 window, legality and query proposal

Extract A:m-by-n and D:n-by-p in the three named degree bases, reducing summed
integer incidence modulo 2 and preserving the original signed source. Check
`AD=0` over F2, k>=0, unique basis keys, explicit empty shapes and n strictly
positive coordinate costs. The first bounded Rust proposal uses exact positive
integer costs with checked arithmetic; rational, floating and arbitrary-precision
costs remain explicit future choices, not silently rounded inputs. Changing a
basis permutes A columns, D rows and costs together. Changing costs or source
claims changes the applicable result identity even if topology is unchanged.

The text-level entry is
`solve_homology_projection(window, costs, requested_certificate, execution)`;
direct reads are `project(x)`, `operator_action(x)`, `class_representative(z)`,
`same_class(z,y)`, `selected_mass(z)` and `support(z)` on its validated result.
Arbitrary chains are allowed for project/action; class queries require `Az=0`
for every argument. Unsupported optimality, arithmetic or representation
requests must remain unsupported; they cannot be downgraded silently to feasible.
No reference solver dispatch registry or snapshot schema is adopted.

For P:n-by-n over F2, legality requires all of

```text
P^2=P, AP=0, PD=0,
for every cycle z: (I+P)z belongs to image(D).
```

It suffices to check the last condition on an independently obtained cycle
basis, with boundary preimages/rank evidence. The zero P passes the first three
conditions even when nonzero homology exists, so those tests alone are unsound.
Legality yields `rank(P)=n-rank(A)-rank(D)` and `ker(I+P)=image(P)`, not
`ker(I+P)=ker(A)`. P chooses a representative complement; it is not a real Hodge
orthogonal projector. Explicit matrix output is optional only when an actual
factorized representation supplies equivalent actions and independently
verifiable legality; matrix-free support is not promised by a label.

For the abstract F2 window `A=[1,1,0]`, `D=(1,1,0)^T`, costs `(2,3,5)`, a
hand choice is `P=diag(0,0,1)`. Its cycle basis is `(1,1,0),(0,0,1)`, the first
is a boundary, beta=1, and `(1,1,1)` projects to `(0,0,1)` of selected mass 5.
Zero P fails homology preservation on `(0,0,1)`. `(1,0,0)` is not a cycle.
The selected representative, support and feasibility evidence refer to this P;
minimum class mass and globally minimum stretch remain separate requests.
Family transport needs the two chain-map equations in [stage maps](#stage-maps)
and compatible identities/cost policies; neither an ordinary barcode nor
CriticalSetWorkspace supplies them. #90 owns later running integration.

### Generalized-model comparison and selected first investigation

The primary readings below cover definitions and input requirements, not a
full reproduction or performance review of each paper. These are four distinct
candidate theories; no generalized model is inserted into the existing trait.

| Candidate and source read | Data missing from current owners | Decision and reason |
| --- | --- | --- |
| Cellular sheaf: Hansen/Ghrist, [arXiv:1808.01513v2](https://arxiv.org/html/1808.01513v2), sections 2.2–3.2 | Stalk spaces and ordered stalk coordinates; face-to-coface linear restrictions with composition, signed coboundaries and stalk Grams. Cell incidence alone supplies none of the restriction values. | Select a one-edge rank-one real sheaf for the next investigation: tests nonconstant restrictions with minimal extra data and a hand-derived spectral expectation. |
| Path homology: Grigoryan/Lin/Muranov/Yau, [arXiv:1207.2834v4](https://arxiv.org/pdf/1207.2834v4), section 3.2 and Example 3.9 | Directed allowed paths, regular/nonregular and reduced/ordinary convention, plus a basis of boundary-invariant chains `Omega[p]={x in allowed[p]: boundary(x) in allowed[p-1]}`. Basis vectors can be combinations of paths. | Defer: `WeightedGraph` is undirected and simplex deletion is not a path-chain constructor. For `0->1->2`, boundary(012)=12-02+01 contains forbidden 02; merely enumerating allowed paths fails. |
| Embedded hypergraph homology: Bressan/Li/Ren/Wu, [arXiv:1610.00890v5](https://arxiv.org/html/1610.00890v5), sections 2 and 3.1–3.2 | Original hyperedge membership and the chosen embedded chain subspaces; simplicial closure changes the mathematical input. | Defer: first choose the intended homology, preserve non-face-closed membership and construct its chain basis. TopoNetX's inspected closure makes a simplicial complex, which does not establish embedded-homology support. |
| Mayer N-complex: Shen/Liu/Wei, [DOI:10.3934/fods.2024032](https://doi.org/10.3934/fods.2024032), abstract's N-differential definition | N, coefficient domain (and any q/root convention), degree bases, N-nilpotent differential, and the chosen homology index/powers. | Defer: `d^N=0` for N>2 permits `d^2!=0`; ordinary PH's square-zero boundary cannot express it. Persistent versions need their own maps and definitions. |

The Mayer full-text endpoint returned a wrapper, so that row uses the primary
abstract only; no detailed Mayer algorithm reading/reproduction is claimed.

The selected sheaf investigation borrows one explicitly oriented graph edge,
sets all three stalks to R with standard Grams, and supplies restrictions
`R(v0,e)=1`, `R(v1,e)=2`. Its coboundary is `delta=[-1,2]`, so the hand
expectation is `L0=delta^T delta=[[1,-2],[-2,4]]`, kernel `span((2,1))` and
eigenvalues `(0,5)`. Unlike the constant-sheaf kernel, the vertex values need
not be equal. No composition through a 2-cell is needed in this one-edge case;
general sheaves must check every face-chain composition. Result identity includes
stalk bases and restriction maps, not merely the underlying graph identity.
This is sheaf cohomology, not a claim that an ordinary PH class has changed.

A future operation can receive `(borrowed graph, ordered stalk bases,
restriction maps, grams)` directly; no six-crate model hierarchy is justified.
Its implementation gate is actual restriction assembly/action, the constant
sheaf specialization, zero/rectangular stalk cases, independent residual and
kernel checks, resource/error behavior and a source/fixture-bound cost report.
No sheaf engine is approved or delivered here. The other candidates remain
deferred for the model-specific gaps above, not cancelled and not blockers of
old PH migration. As a separate hand counterexample, maps `C2->C1->C0` all equal
to [1] and zero outside satisfy `d^3=0` but `d^2:C2->C0=[1]`; feeding them to
`FilteredComplex` would violate its contract rather than extend it.

The [family output and resource contracts](result-contracts.md#computation-family-outputs-and-resources)
complete these input proposals. Their future algorithms and #88–#90 running
integration must establish their own actual support and evidence.

### Recheck the hand tables without a solver

Save this block as `target/family-hand-check.py` and run
`python3 target/family-hand-check.py`. It multiplies the stated tiny matrices and
checks exact polynomial/coefficient identities. It performs no eigendecomposition,
projection construction or production algorithm. The source/basis extraction
remains exercised by the existing interface-guide doctest.

```python
from fractions import Fraction as Q
from itertools import product

def mul(a, b):
    return [[sum(x * b[k][j] for k, x in enumerate(row))
             for j in range(len(b[0]))] for row in a]

def adj(a):
    return [[x.conjugate() for x in row] for row in zip(*a)]

def add(a, b):
    return [[x + y for x, y in zip(r, s)] for r, s in zip(a, b)]

path = [[-1, 0], [1, -1], [0, 1]]
assert mul(adj(path), path) == [[2, -1], [-1, 2]]
g0 = [[1, 0, 0], [0, 2, 0], [0, 0, 3]]
g1 = [[2, 0], [0, 3]]
l = mul([[Q(1, 2), 0], [0, Q(1, 3)]], mul(adj(path), mul(g0, path)))
assert l == [[Q(3, 2), -1], [Q(-2, 3), Q(5, 3)]]
assert mul(g1, l) == [[3, -2], [-2, 5]]
assert l[0][0] + l[1][1] == Q(19, 6)
assert l[0][0] * l[1][1] - l[0][1] * l[1][0] == Q(11, 6)
# Characteristic polynomial: 6 lambda^2 - 19 lambda + 11.
a = [[1, 1, 0], [-1, 0, 1], [0, -1, -1]]
d = [[1], [-1], [1]]
assert mul(a, d) == [[0], [0], [0]]
assert mul(mul(adj(a), a), d) == [[0], [0], [0]]
assert add(mul(adj(a), a), mul(d, adj(d))) == [[3, 0, 0], [0, 3, 0], [0, 0, 3]]
k = [[0, 0, -1], [0, 0, 1], [-1, 1, 0]]
assert mul(k, k) == [[1, -1, 0], [-1, 1, 0], [0, 0, 2]]
assert mul(k, [[1], [1], [0]]) == [[0], [0], [0]]
assert mul([[0]], [[0]]) != [[2]]  # Truncated C1 versus full L1.
# A complex square-zero window; using transpose instead of adjoint is wrong.
ca, cd = [[1, 1j]], [[-1j], [1]]
assert mul(ca, cd) == [[0]]
assert add(mul(adj(ca), ca), mul(cd, adj(cd))) == [[2, 0], [0, 2]]

fa, fd = [[1, 1, 0]], [[1], [1], [0]]
p = [[0, 0, 0], [0, 0, 0], [0, 0, 1]]
assert mul(p, p) == p and mul(fa, p) == [[0, 0, 0]]
assert mul(p, fd) == [[0], [0], [0]] and mul(fa, fd)[0][0] % 2 == 0
for z in product((0, 1), repeat=3):
    if (z[0] + z[1]) % 2 == 0:
        pz = [r[0] % 2 for r in mul(p, [[x] for x in z])]
        assert [(x + y) % 2 for x, y in zip(z, pz)] in ([0, 0, 0], [1, 1, 0])
assert [0, 0, 1] not in ([0, 0, 0], [1, 1, 0])  # Zero P loses this class.
assert sum(w * x for w, x in zip((2, 3, 5), (0, 0, 1))) == 5
delta = [[-1, 2]]
sheaf_l = mul(adj(delta), delta)
assert sheaf_l == [[1, -2], [-2, 4]]
assert mul(delta, [[2], [1]]) == [[0]]
assert sheaf_l[0][0] + sheaf_l[1][1] == 5
assert sheaf_l[0][0] * sheaf_l[1][1] == sheaf_l[0][1] * sheaf_l[1][0]
nd = [[0, 1, 0], [0, 0, 1], [0, 0, 0]]  # Basis order C0,C1,C2.
assert mul(nd, nd)[0][2] == 1
assert mul(mul(nd, nd), nd) == [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
print("Hand matrix, Gram, Dirac, F2 legality, sheaf and N-complex checks passed.")
```

Local preparation checks on 2026-10-11 passed with Python 3.12.13 and Rust
1.98.1: the block above, source/Markdown checks, the existing construction
example and its four tests, and the three interface-guide doctests. The primary
local collaboration checker also passed. All production source, tests, examples,
workers, Cargo files and workflows remain byte-identical to the fixed main.
No external algorithm runtime, native oracle, performance, new feature or hosted
CI result is inferred from these documentation checks.

[cells]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/complex/filtered.rs
[simplex]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/complex/simplicial/simplex.rs
[storage]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/complex/simplicial/mod.rs
[incidence]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/complex/simplicial/incidence.rs
[reader]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/persistence/filtered.rs
[math]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/docs/reference/mathematics.md
[kernel]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/docs/design/kernel.md
[window]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/src/homology_operator/chain.py
[identities]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/src/homology_operator/result.py
[family]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/src/homology_operator/family.py
[operator-interface]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/docs/INTERFACE.md
[operator-architecture]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/docs/ARCHITECTURE.md
[operator-results]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/docs/RESULT_MODEL.md
[toponetx]: https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/toponetx/classes/simplicial_complex.py
[toponetx-cell]: https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/toponetx/classes/cell_complex.py
