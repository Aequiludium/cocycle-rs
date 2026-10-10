# Finite filtration stages and source-preserving inclusions

[Documentation](../README.md) / [Stage usage](../guides/filtered-complexes.md#read-stages-and-their-inclusions)

Research and implementation decision, 2026-10-09, for
[#67](https://github.com/Aequiludium/cocycle-rs/issues/67),
[#68](https://github.com/Aequiludium/cocycle-rs/issues/68) and
[#69](https://github.com/Aequiludium/cocycle-rs/issues/69).
The source baseline is Cocycle `cedf5965b774855614b761a1ad311d9dc7992c49`.
The earlier [interface decision](../design/interface-compatibility.md) confirms
stable IDs, signed incidence and consumer-owned bases. This capability uses that
existing contract; it does not certify every later modeling or spectral framework.

## Fixed-source reading and resulting decisions

| Source inspected | Observation | Decision and remaining limit |
| --- | --- | --- |
| Cocycle baseline: `complex/simplicial`, `filtration/simplicial`, `persistence/source.rs` | Frozen simplices use value/dimension/decreasing-colex order; construction wrappers separately certify scale and dimension | Borrow a prefix with the existing IDs and incidence; preserve certificates in contextual analysis |
| [Oineus filtration](https://github.com/anigmetov/oineus/blob/e52814a1ffb5b8a81e71ff1e93b4c14194673f0f/include/oineus/filtration.h) | `sort_and_set` orders by dimension, value (possibly negated), then original ID; UID lookup and original/sorted permutations are separate; `set_values` re-sorts | Do not copy a dimension-major order into a global scale prefix. Map by original owner/identity, never by equal numeric sorted indices |
| [Oineus inclusion filtration](https://github.com/anigmetov/oineus/blob/e52814a1ffb5b8a81e71ff1e93b4c14194673f0f/include/oineus/inclusion_filtration.h) | Owns domain/codomain filtrations and locates boundary rows by UID; its size checks impose additional restrictions | Useful evidence that two orderings need explicit correspondence; not a general chain-map correctness oracle or a requirement to copy both stages |
| [Oineus mapping-cylinder draft](https://github.com/anigmetov/oineus/blob/e52814a1ffb5b8a81e71ff1e93b4c14194673f0f/include/oineus/simpl_map_filtration.h) | Takes an index-vector mapping; the inspected constructor is incomplete | Do not infer runnable arbitrary-map support from the header or README. No algorithm text is copied |
| [GUDHI Simplex_tree](https://github.com/GUDHI/gudhi-devel/blob/4ec34ac55e6d2e8cfd7c322e84c1b0a56d516d51/src/Simplex_tree/include/gudhi/Simplex_tree.h) | `filtration_simplex_range` requires a total-order convention; custom initialization for multiple parameters orders a selected scalar problem. Insertion can add subfaces/change existing values | Compare simplex/value sets independently; preserve Cocycle's reject-rather-than-repair policy. Ordering a multi-parameter object does not deliver multi-parameter persistence |

The three cached Oineus headers were matched to fixed Git blobs
`113c945d2d29b3ee854512302557df320bc7ed96`,
`e6159944bf09814d68cb3a566f5c89901f3f88a1` and
`a915c1fe80abf81d1af6b13ffde7c0f264cbf7e8`, respectively. Source reading establishes
the stated observations, not Oineus runtime correctness or feature parity.
The GUDHI runtime probe below uses Python GUDHI 3.12.0, separately from the pinned
C++ source inspection and the repository's native GUDHI/Ripser comparisons.

## Research acceptance and implementation handoff

The research exit for #68 is distinct from the implementation gate for #69.
On 2026-10-10, prerequisite #60 was accepted through merged
[PR #99](https://github.com/Aequiludium/cocycle-rs/pull/99). The following decisions
cover the research acceptance criteria and define the bounded implementation.

| Research requirement | Decision and evidence in this report |
| --- | --- |
| Scale, order, truncation and completeness | Inclusive scalar sublevels; original units and IDs; face-before-coface ties; separate scale and skeleton coverage; three hand-derived finite families |
| Valid and invalid stage maps | Same-owner identity inclusions commute with integer boundaries and compose; foreign owners, reverse scales and mixed certificates fail; an oriented-edge swap explains the signed chain-map requirement |
| Multi-parameter scope | Incomparable parameter vectors need a partial order and a separately defined invariant; a total-order projection does not inherit an ordinary barcode guarantee |

The decision is to deliver borrowed finite simplicial stages and same-source
inclusions. General chain maps, generic-cell stage owners, relabeling,
multi-parameter invariants and persistent spectral operators remain separate
research work. #69's implementation and runnable examples are in
[PR #100](https://github.com/Aequiludium/cocycle-rs/pull/100); its declared #66
modeling prerequisite has a runnable delivery in
[PR #102](https://github.com/Aequiludium/cocycle-rs/pull/102), independently awaiting
review and merge. Research acceptance does not merge either implementation or
publish a release.

## Scales, ties and what is complete

| Finite example | Value convention | Four stages and independently derived sizes |
| --- | --- | --- |
| Exact Rips on line points 0, 1, 2 | Maximum edge length; vertices enter at zero | Scales -1, 0, 1, 2: 0, 3, 5, 7 simplices |
| Lower-star four-cycle 01, 12, 23, 03 with vertex values (-2, 1, -1, 0) | Maximum supplied vertex value in each simplex | Scales -2, -1, 0, 1: 1, 2, 5, 8 simplices |
| Supplied triangle: vertex births -2, -1, 0; edges at 1; face at 2 | Caller-defined scalar, possibly unrelated to geometry | Scales -1, 0, 1, 2: 2, 3, 6, 7 simplices |

All values and query scales are finite binary64; signed zero is canonicalized.
The comparison is inclusive, so an entire tied event enters at once. Faces precede
cofaces, including at equal values. A stage before the first event is empty; an
event-free scale gap does not authorize a reverse inclusion. Rips builders still
require nonnegative distances/construction thresholds. A negative stage of an
already constructed Rips complex is a well-defined empty sublevel set.

`EdgeLength` declares a parameter convention, not metres or normalization.
Lower-star values keep the scalar field's units; supplied values use
`Unspecified`. Same-source stages share that interpretation automatically.
Cross-source units and rescaling remain caller obligations; this implementation
rejects cross-source inclusions rather than guessing those facts.

Scale coverage and dimension coverage are independent. `Coverage::Through(t)`
certifies scales only through t; stage selection beyond t fails. A requested
maximum simplex dimension is a construction cap, not proof that omitted cofaces
are absent. `source_filtration()` exposes the unchanged certificate. For example,
the complete Rips of these three points contains a triangle at 2; constructing
only edges cannot certify H1 at 2. Contextual persistence retains the existing
conservative `InsufficientSkeleton` check, even when a smaller particular stage
might admit a stronger certificate. No new scale-dependent dimension proof is
inferred. Sparse Rips stages retain the actual stored blocked topology and its
approximation metadata, rather than expanding the modified graph as an ordinary
flag complex.

A bare supplied complex is complete relative to its supplied topology. A stage
of that source is an observation cap, not automatically a new complete source.
Contextual analysis censors surviving classes while later original events remain
unobserved. Choosing the prefix via `PersistenceBuilder::from_complex` is the
explicit alternative source interpretation described in the usage guide.

## Minimum stage and map information

One `SimplicialStage` borrows the immutable owner, the optional construction
wrapper, its finite scale and a prefix length. Source-local IDs, increasing-vertex
orientation and integer signed boundary terms are reused. No new source trait,
owned stage container, weight convention or algorithm registry is needed.
Selection takes O(log n) time and O(1) storage. Degree-specific basis arrays,
matrix exports and owned copies are consumer-requested costs; the stage retains
the full source for its borrow lifetime. It does not free excluded source storage.

An inclusion from K_s to K_t requires the same owner and certificate interpretation
and s <= t. Its map is the identity on existing simplex IDs with coefficient +1;
integer boundaries commute and compositions are identities. Valid examples include
the ring at 1 into its filled triangle at 2, empty into nonempty, and equal scales.
Different owners, even clones with equal indices, reverse scales, and mixing a
certified wrapper with its bare complex are rejected before returning a map.
Individual `SimplexId` values cannot detect foreign ownership; callers must keep
them with their source or use the checked stage-pair operation.

An arbitrary degree-preserving map needs explicit bases and coefficients with
D_q(B) F_q = F_(q-1) D_q(A). For an edge oriented a to b, swapping a and b while
keeping F_1(edge)=edge is invalid over the integers: the two sides are a-b and
b-a. Choosing F_1(edge)=-edge makes the swap a valid chain map. Checking only
dimensions or only F2 loses this sign distinction. General simplicial maps also
need explicit treatment of collapsed simplices. Those operations are deferred,
not accepted by `inclusion_into` and not described as delivered.

## First capability and deferred research

Deliver finite single-parameter stages of supplied or explicitly constructed
simplicial sources and their same-source inclusions. The existing implicit Rips
paths stay available without materializing all stages. Generic custom cells are
deferred: their IDs need not encode positions, and trait implementations own
integer validity; a new checked stage owner needs a concrete consumer first.

For multiple parameters, vectors (0, 1) and (1, 0) are incomparable under
coordinatewise order. Lexicographic sorting or projecting to one scalar chooses
a different problem. Future work must define vector-valued face monotonicity,
query partial order, compatible maps/commuting diagrams and the requested
invariant before choosing storage or a solver. Ordinary barcodes, scalar stage
spectra, persistent Laplacians and general F2 class transport make different
claims; none follows from this stage access API alone.

## Runnable independent checks

The public contract tests and example cover all three finite families, signed and
tied scales, closure, integer incidence, inclusion composition, foreign owners,
coverage/skeleton errors, representatives and execution limits:

```sh
cargo test --locked --test filtered_complex --example complex_construction
cargo run --locked --example complex_construction
```

The following optional Python probe compares the actual lower-star example's four
stage cell/value sets against independent GUDHI construction. Save it under
ignored `target/` and run from the repository root with Python and GUDHI 3.12.0.
It does not use persistence to construct either filtration; ordering is compared
as simplex/value sets because the libraries' tie conventions differ.

```python
import ast
import subprocess
import gudhi

assert gudhi.__version__ == "3.12.0"
output = subprocess.check_output(
    ["cargo", "run", "--locked", "--example", "complex_construction"], text=True)
values = [-2., 1., -1., 0.]
tree = gudhi.SimplexTree()
for vertex, value in enumerate(values):
    tree.insert([vertex], filtration=value)
for edge in [(0, 1), (1, 2), (2, 3), (0, 3)]:
    tree.insert(edge, filtration=max(values[v] for v in edge))
cells = [(tuple(sorted(vertices)), value) for vertices, value in tree.get_filtration()]
count = 0
for line in output.splitlines():
    if not line.startswith("Stage "):
        continue
    label, data = line.split(": ", 1)
    scale = float(label.removeprefix("Stage "))
    actual = sorted((tuple(vertices), value) for vertices, value in ast.literal_eval(data))
    expected = sorted(cell for cell in cells if cell[1] <= scale)
    assert actual == expected, (scale, actual, expected)
    count += 1
assert count == 4
print("GUDHI lower-star stage sets: 4/4")
```

Finite fixture parity is not a proof for arbitrary user data. The PR records
the exact tested implementation, additional stage probe scope, local/hosted
checks and remaining prerequisite status. Existing native reference CI compares
the underlying supported Rips paths; it is distinct from stage-specific tests
and does not certify arbitrary maps or physical units. No performance ranking,
general no-regression timing result, merge or release is established by this report.

Local stage-specific validation on Windows used Rust 1.98.1, Python 3.10.11 and
GUDHI 3.12.0. An additional standalone public-API probe checked all three families:
12 stage cell/value sets, 24 stage/field diagrams over F2/F3 and 24 Betti readouts
agreed. The probe projects GUDHI's full-source endpoints to the same observation
cap, keeping finite, essential and censored categories distinct. Rust source
hashes and raw outputs are retained locally; they are not a published native
artifact bundle. The smaller committed probe above replays four stage-set checks.
Integer boundary/map validity is checked separately by hand-derived Rust tests.

Native Ripser at revision `01add51ff64aaf40889483260cc5c3b7d0f2a1e7` additionally
matched the line-Rips fixture's six supported stage/field diagrams (scales 0, 1,
2 over F2/F3), with explicit threshold 2. Its source SHA-256 was
`6ef9c828944316974ed62408a55b351cba896c8e548eb7193af733d005585bf0`;
GCC 13.3.0 built it under Ubuntu WSL with `USE_COEFFICIENTS`. These integer
distances are exactly representable in both Ripser f32 and Cocycle f64. The
negative empty stage was checked analytically; lower-star and arbitrary supplied
filtrations are excluded from Ripser input support. The first shell invocation
failed while passing the modulus argument; its failure is retained separately
from the successful explicit-argument runs. This small probe is not the full
native reference suite or a timing comparison.
