# F2 homology operators

[Documentation](../README.md) / Guides

Feasible F2 operators are an **Unreleased** opt-in addition in this checkout.
See the [version status](../README.md#version-and-capability-status). Use
`homology_operators` for one weighted chain window or a finite sequence of
same-source stages. The result selects class representatives; it does not solve
for the shortest representative or a globally minimum-stretch projection.

The [mathematical specification](../reference/mathematics.md#21-feasible-f2-homology-projections)
derives the projection and transport identities. The
[comparison report](https://github.com/Aequiludium/cocycle-rs/blob/codex/homology-operator-integration/benches/reports/homology-operators.md)
records exact oracle agreement, measured scopes and unfavorable results.
No Python or foreign solver is a Rust runtime dependency. This implementation
requires Rust 1.91 and follows the pinned
[homology-operator source](https://github.com/proffitteoy/homology-operator/tree/3daff9bd69033d1ee9c5cba5026bfff94df8d6ac),
with its MIT attribution in [third-party notices](../../THIRD_PARTY_NOTICES.md).

## Query an ordered chain window

```rust
use cocycle::homology_operators::{ChainWindow, F2Matrix, HomologyOperator};

let a = F2Matrix::from_rows(&[vec![1, 1, 0]], 3)?;
let d = F2Matrix::from_rows(&[vec![1], vec![1], vec![0]], 1)?;
let window = ChainWindow::new(
    1, [a, d],
    [vec!["vertex"], vec!["edge-a", "edge-b", "edge-c"], vec!["face"]],
    vec![2, 3, 5], "sample-42/snapshot-1", "abstract cost",
)?;
let operator = HomologyOperator::feasible(window)?;
assert_eq!(operator.betti(), 1);
let query = operator.class_query(&[0, 0, 1])?;
assert_eq!(query.representative(), &[0, 0, 1]);
assert_eq!(query.selected_mass(), 5);
assert_eq!(query.support().collect::<Vec<_>>(), vec![2]);
assert!(operator.is_boundary(&[1, 1, 0])?);
assert!(operator.class_query(&[1, 0, 0]).is_err()); // A noncycle.
# Ok::<(), cocycle::Error>(())
```

A maps the current degree to the previous degree; D maps the next degree to the
current degree. Shapes are m-by-n and n-by-p, and AD=0 over F2. Three explicit
ordered bases have lengths m/n/p and distinct IDs within each degree. Numeric
IDs alone do not identify an owner. Ordinary H0 has an empty previous space;
0-by-n and n-by-0 matrices retain their shapes. An empty current space is valid.

Costs are exactly n strictly positive `u64` values, in current-basis order, with
a declared unit. They are independent of filtration values and real Gram
matrices. Checked mass overflow returns `Error::NumericalFailure`; no saturated
or approximate cost is returned. A source label is a caller-declared snapshot
namespace and cannot authenticate a dataset. Rational, floating and arbitrary-
precision costs are not accepted by this port.

`F2Matrix::from_rows` rejects nonbinary values and shape mismatches.
`from_sparse_columns` accepts distinct in-range nonzero row indices per column,
including empty columns, without a dense input buffer. Internal sparse/packed
storage and elimination workspaces are private.

## Interpret and retain a query

`HomologyOperator` owns its validated window. `window()` borrows its original
bases, costs and source declaration without copying them. `project(x)` accepts
any binary chain; `class_query(z)` additionally requires Az=0. An imported
`from_projection` checks all four legality conditions against the actual A/D;
a zero P is rejected when it erases a nonzero homology class.

`selected_mass` is the weight of the chosen representative. `class_distance`
uses the chosen representative of the difference class, not a shortest-class
optimization. `stretch_with` exhaustively evaluates this P on nonzero cycles;
it does not search other legal P. L=I+P is an F2 action, not a real Hodge
Laplacian, spectral solver or eigenvalue approximation.

A `ClassQuery` owns a shared immutable window and its actual operator-run
binding. It remains readable after dropping the operator. `same_class`,
`shared_support` and `union_support` require the same run. Cloning an operator
preserves that binding; a fresh solve/import creates a new one, even if P and
costs are equal. Rebuilding topology or changing costs cannot reinterpret an
existing record. There is no Python snapshot/history serialization contract.

## Transport through borrowed stages

```rust
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::homology_operators::OperatorFamily;

let source = SimplicialComplex::new(vec![
    Simplex::new(vec![10], 0.)?, Simplex::new(vec![20], 0.)?,
    Simplex::new(vec![30], 0.)?, Simplex::new(vec![10, 20], 1.)?,
    Simplex::new(vec![10, 30], 1.)?, Simplex::new(vec![20, 30], 1.)?,
    Simplex::new(vec![10, 20, 30], 2.)?,
])?;
let family = OperatorFamily::from_stages(
    &[source.stage(1.)?, source.stage(2.)?], 1,
    vec![vec![2, 3, 5]; 2], "triangle-v1", "abstract cost",
    &Execution::default(),
)?;
let query = family.stage(0)?.class_query(&[1, 1, 1])?;
let tracked = family.track_query(&query, 0, 1)?;
assert_eq!(query.selected_mass(), 10);
assert_eq!(tracked.selected_mass(), 0);
let target = family.stage(1)?;
let filling = target.boundary_preimage(&[1, 1, 1])?.unwrap();
assert_eq!(filling, vec![1]);
assert_eq!(target.window().boundaries()[1].apply(&filling)?, vec![1, 1, 1]);
# Ok::<(), cocycle::Error>(())
```

`from_stages` checks the real `SimplicialStage::inclusion_into` relationship
before extraction. Different owners, decreasing scales or mixed certified/bare
interpretations are rejected. Equal cloned complexes and equal source labels
are not enough. Scales come from the stages; negative and repeated scales keep
their meaning. One execution budget covers extraction, validation, solving and
family checks. The owned family outlives the borrowed stages and their source.
It describes supplied windows, without certifying unseen cofaces.

For a general `FilteredComplex`, `ChainWindow::from_complex` reads signed
incidence, sums repeated terms modulo two, and preserves degree-restricted
source traversal. `through` optionally selects value <= scale. Use
`from_complex_with` for controlled extraction or pass a selected stage with
`through=None`. This validates the F2 window, not the integer chain identity or
completeness of a larger geometric source.

`OperatorFamily::new` accepts already validated operators with caller-declared
sources and nested coordinate bases. It checks both three-degree chain squares.
It cannot authenticate their physical source owner. Shared-coordinate costs
must agree; `with_weight_policy(..., WeightPolicy::Variable)` explicitly permits
changed target geometry. Such transport does not assert an isometry or mass bound.
General linear chain maps are outside the coordinate-inclusion interface.

`track_query` rejects a record from the wrong source run and returns a record
bound to the target run. `track_class` accepts raw coordinates in the declared
source basis. To check class preservation, embed the source representative in
the target basis, XOR the target representative, and solve that residual with
`boundary_preimage`. `Some` contains next-degree coefficients; `None` means no
preimage exists, while `Some(vec![])` is a valid empty witness. Legal projections
need not commute as matrices.

`transport` uses deterministic source/target homology coordinates. `barcode`
returns stage-index intervals with multiplicity. A final `death=None` means
unknown continuation beyond the observed stages. Repeated scales do not collapse
stage IDs, and the result is not implicitly converted to a complete PH diagram.

## Control work and understand storage

Controlled extraction, projection, query, tracking, filling, transport, barcode
and stretch methods use fresh `Execution` budgets. Failure or cancellation
returns no apparently successful partial result. Source access, sorting and
allocation are checked around their calls and cannot be interrupted internally.
Work units count arithmetic/coordinate visits, including query mass visits;
they are version-dependent, not bytes, deadlines or a memory quota.

Wide matrices start with ordered sparse row indices when a row needs at least
eight packed words. They switch to packed storage when total nonzeros exceeds
half the total packed-word count; this is a private representation heuristic,
not a 50% bit-density threshold. Factor actions keep vectors packed, coordinate
solves reuse one request-local row transform, and acceptance recomputes cycle/
image reductions independently of the solver inverses.

Finite barcode preparation composes adjacent homology-coordinate maps and keeps
only two stage bases during preparation. The rank table and explicit maps still
have quadratic worst-case costs. Stretch is exponential in dim ker(A), because
input mass depends on cycle representatives, not only quotient coordinates.
There is no persistent family cache, streaming barcode, historical barcode-basis
algorithm or infinite-filtration conclusion. Run
`cargo run --locked --example homology_operators` for the stage workflow; use the
[verification protocol](https://github.com/Aequiludium/cocycle-rs/blob/codex/homology-operator-integration/benches/homology-operators.md)
for external comparisons and distinct timing boundaries.
