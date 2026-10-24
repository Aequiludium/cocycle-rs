# Borrow inputs, read results and control composed calls

[Documentation](../README.md) / Usage

This guide follows the current checkout, including **Unreleased** stage calls.
See the [version status](../README.md#version-and-capability-status) for release scope
and the [architecture](../development/architecture.md#six-responsibilities-and-data-routes)
for current owners and future capability boundaries.

These examples use existing public interfaces. They connect modeling, stage
selection, persistence, result access, comparison and analysis without a common
source trait or a universal result container. The
[interface compatibility report](../design/interface-compatibility.md) records
the decisions and the remaining spectral/operator boundaries.

The [data-modeling decision](../research/data-modeling.md) selects existing
point, graph and simplex owners. The two readers below share one frozen complex:
signed single-scale boundary export and ordinary persistence. Run
`cargo run --locked --example complex_construction` for another runnable version
that constructs one lower-star circle and reads it both ways.

## Read a signed stage for a matrix consumer

Borrow a frozen complex and its `stage(scale)` with cells of value <= scale. Preserve
degree-specific basis order, owner-local IDs and increasing-vertex orientation.
The explicit conversion below returns owned dense matrices only because the
small external consumer needs them. It is tutorial code, not a new public API.

```rust
use cocycle::complex::{FilteredComplex, Simplex, SimplexId, SimplicialComplex};
use cocycle::diagram::IntervalEnd;
use cocycle::persistence::PersistenceExt;

let complex = SimplicialComplex::new(vec![
    Simplex::new(vec![10], 0.)?, Simplex::new(vec![20], 0.)?,
    Simplex::new(vec![30], 0.)?, Simplex::new(vec![10, 20], 1.)?,
    Simplex::new(vec![10, 30], 1.)?, Simplex::new(vec![20, 30], 1.)?,
    Simplex::new(vec![10, 20, 30], 2.)?,
])?;
let source = &complex;
let scale = 1.;
let stage = source.stage(scale)?;
let filled = source.stage(2.)?;
assert_eq!(stage.len(), 6);
assert_eq!(stage.inclusion_into(&filled)?.len(), 6);
let bases: [Vec<SimplexId>; 3] = std::array::from_fn(|degree| {
    source.cells()
        .filter(|&id| {
            let s = source.simplex(id).unwrap();
            s.dimension() == degree && s.value() <= scale
        }).collect()
});
let keys: Vec<_> = bases[1].iter()
    .map(|id| source.simplex(*id).unwrap().vertices()).collect();
assert_eq!(keys, [vec![20, 30], vec![10, 30], vec![10, 20]]);

// Include shapes even when a matrix has no rows or columns.
let boundary = |degree: usize| {
    let rows = bases[degree - 1].len();
    let columns = bases[degree].len();
    let mut matrix = vec![vec![0_i32; columns]; rows];
    for (column, id) in bases[degree].iter().enumerate() {
        for term in source.boundary(*id).unwrap() {
            let row = bases[degree - 1].iter()
                .position(|face| *face == term.face).unwrap();
            matrix[row][column] = i32::from(term.coefficient);
        }
    }
    ((rows, columns), matrix)
};
let (a_shape, a) = boundary(1);
let (d_shape, d) = boundary(2);
assert_eq!(a_shape, (3, 3));
assert_eq!(a, [[1, 1, 0], [-1, 0, 1], [0, -1, -1]]);
assert_eq!(d_shape, (3, 0));
assert_eq!(d, vec![Vec::<i32>::new(); 3]);

// F2 coefficients are an explicit conversion; keep the original signed data.
let a_f2: Vec<Vec<_>> = a.iter()
    .map(|row| row.iter().map(|x| x.rem_euclid(2)).collect()).collect();
assert_eq!(a_f2, [[1, 1, 0], [1, 0, 1], [0, 1, 1]]);

// A direct calculation at scale 1, with standard real inner products.
// The signed edge cycle x has B1*x=0. No persistence request was needed.
let x = [1, -1, 1];
let energy: i32 = a.iter()
    .map(|row| row.iter().zip(x).map(|(b, x)| b * x).sum::<i32>().pow(2))
    .sum();
assert_eq!(energy, 0); // x^T B1^T B1 x; C2 is empty at this scale.

// The other reader uses this exact owner and its stored filtration values.
let ids_before: Vec<_> = source.cells().collect();
let result = source.persistence().max_homology_dimension(1).compute()?;
let h1 = result.diagram().dimension(1)?.iter().next().unwrap();
assert_eq!((h1.birth(), h1.end()), (1., IntervalEnd::Finite(2.)));
assert_eq!(source.cells().collect::<Vec<_>>(), ids_before);
assert_eq!(source.simplex(bases[1][0]).unwrap().vertices(), &[20, 30]);
# Ok::<(), cocycle::Error>(())
```

The selected degree bases allocate O(number of selected cells) handles and
borrow the source's vertex labels. The dense export allocates one entry per
matrix coordinate and uses linear basis lookups; it is a deliberately small
example, not a recommended large-complex representation. Sparse consumers can
iterate signed incidence instead. The source is unchanged, and no persistence
request is needed for this conversion.

`FilteredComplex::cells()` already supplies the stored IDs without vertex-key
lookups. A direct read of the entire object can omit the scale predicate, as the
construction example does. The stage selection above includes all ties, while
the PH query reads the full filtration and obtains H1 `[1,2)`. A scale-1 truncated
PH query would instead retain a right-censored endpoint; the direct stage's
cycle does not certify essentiality of that truncated query.

For an ordinary nonaugmented H0 consumer, the previous space is empty and its
boundary has shape 0-by-n. Export that shape explicitly. Simplex IDs belong to
this owner; equal integer indices from another complex are not matching cells.
Signed simplicial incidence can be reduced into F2; arbitrary modular input
cannot be lifted into a real chain complex by replacing each bit with 0 or 1.

A single-scale graph consumer can similarly borrow `WeightedGraph::edges()`
and derive vertex-edge incidence. Its interpretation as a one-dimensional
complex is explicit: a graph triangle is not automatically a filled clique.
Filtration values, spectral Grams and F2 coordinate costs remain separate.
This conversion prepares inputs; it implements no Laplacian, Dirac or F2
projection solver. Separate candidate [Hodge](hodge-laplacian.md),
[Dirac](dirac.md) and [persistent Laplacian](persistent-laplacian.md) guides
show how their interfaces consume signed stages and preserve result coordinates.

The stage borrow checks the inclusion into the filled triangle; the matrix
reader selects the same source-local IDs with the explicit scale predicate.
Stage selection itself allocates no storage. See the
[stage guide](filtered-complexes.md#read-stages-and-their-inclusions) for checked
inclusions and for the distinction between analyzing the original source through
a stage and choosing its stored topology as a new complete supplied source.

## Move persistence data while retaining its meaning

Here the point coordinates use one caller-declared dimensionless scale. Both
comparisons concern that same source. A matching scale convention alone does
not certify physical units across unrelated datasets.

```rust
use cocycle::algebra::PrimeField;
use cocycle::descriptors::{betti_curve, finite_lifetime_summary};
use cocycle::diagram::{Coverage, PersistenceData};
use cocycle::diagram_distances::bottleneck_distance_results;
use cocycle::filtration::RipsBuilder;
use cocycle::geometry::PointCloudView;
use cocycle::persistence::{PersistenceExt, RepresentativeRequest, RepresentativeSelection};
use cocycle::Error;

let coordinates = [0., 1., 2.];
let points = PointCloudView::new(&coordinates, 3, 1)?;
let source = RipsBuilder::from_points(points);
let requests = [RepresentativeRequest::new(0, 0.5, RepresentativeSelection::Both)?];
let result = source.persistence().representatives(&requests).compute()?;
let borrowed: &PersistenceData = result.as_ref();
assert!(std::ptr::eq(borrowed.diagram(), result.diagram()));
assert!(std::ptr::eq(borrowed.context(), result.context()));
assert_eq!(borrowed.diagram().coverage(), Coverage::Complete);

let (data, representatives) = result.into_parts();
for representative in representatives.unwrap() {
    let interval = data.diagram().interval(representative.interval_index()).unwrap();
    assert_eq!(interval.dimension(), representative.dimension());
    assert!(interval.birth() <= representative.scale());
}
let h0 = data.diagram().dimension(0)?;
assert_eq!(h0.len(), 3); // Two finite bars, one essential bar.
assert_eq!(betti_curve(data.diagram(), 0, &[0., 1.])?, [3, 1]);
assert_eq!(finite_lifetime_summary(data.diagram(), 0)?.finite_count(), 2);
assert_eq!(bottleneck_distance_results(&data, &data, 0)?, 0.);

let mod3 = source.persistence().field(PrimeField::new(3)?).compute()?;
for incompatible in [
    bottleneck_distance_results(&data, &mod3, 0),
    bottleneck_distance_results(&mod3, &data, 0),
] {
    assert!(matches!(incompatible, Err(Error::IncompatibleDiagramContext { .. })));
}
let diagram = data.into_diagram(); // Explicitly discard mathematical context.
assert_eq!(diagram.dimension(0)?.len(), 3);
# Ok::<(), cocycle::Error>(())
```

`AsRef<PersistenceData>` borrows existing diagram/context without cloning or
reconstructing them. `into_parts()` moves both data and optional representatives;
their logical interval ordinals still address the moved diagram. `into_data()`
instead discards representatives. `into_diagram()` also discards context.
Copying an external array or calling `clone()` is a separate explicit operation.

For an empty computed dimension, `dimension(k)?` returns a valid empty view.
An uncomputed dimension, including a gap below `max_dimension()`, is an error.
Traverse `computed_dimensions().iter()`, and check every requested dimension
on every operand. Raw diagram distance calls do not check fields or source
conventions. Result adapters check fields and declared edge-length conventions;
physical units, normalization and dataset provenance remain caller obligations.

## Reuse controls without implying a cumulative budget

Each terminal operation creates a fresh private work budget. One persistence
builder shares it across that operation's preparation, reduction and requested
representatives; one distance call shares it across validation and matching.
Reusing an `Execution` across separate calls reuses limits and the cancellation
flag, not an already-consumed counter. Discover a sufficient limit for this
fixture rather than depending on unstable internal work counts.

```rust
use cocycle::diagram::{Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval};
use cocycle::diagram_distances::{bottleneck_distance_with, wasserstein_2_euclidean};
use cocycle::execution::Execution;
use cocycle::Error;
use std::sync::atomic::{AtomicBool, Ordering};

let first = PersistenceDiagram::new(0, Coverage::Complete, vec![
    PersistenceInterval::new(0, 0., IntervalEnd::Finite(2.))?,
])?;
let empty = PersistenceDiagram::new(0, Coverage::Complete, vec![])?;
let sufficient = (0..=4096).find(|&limit| {
    bottleneck_distance_with(&first, &empty, 0,
        &Execution::default().max_work(limit)) == Ok(1.)
}).expect("small fixture must complete within this search range");
assert!(sufficient > 0);
let controls = Execution::default().max_work(sufficient);
assert_eq!(bottleneck_distance_with(&first, &empty, 0, &controls)?, 1.);
assert_eq!(bottleneck_distance_with(&first, &empty, 0, &controls)?, 1.);
assert!(matches!(bottleneck_distance_with(&empty, &empty, 0,
    &Execution::default().max_work(0)), Err(Error::WorkLimitExceeded { .. })));

let flag = AtomicBool::new(true);
let cancellation = Execution::default().cancellation(&flag);
assert_eq!(bottleneck_distance_with(&first, &empty, 99, &cancellation),
    Err(Error::Cancelled)); // Pre-cancellation precedes dimension validation.
flag.store(false, Ordering::Relaxed);
assert_eq!(bottleneck_distance_with(&first, &empty, 0, &cancellation)?, 1.);

let essential = PersistenceDiagram::new(0, Coverage::Complete, vec![
    PersistenceInterval::new(0, 0., IntervalEnd::Essential)?,
])?;
assert_eq!(bottleneck_distance_with(&essential, &empty, 0,
    &Execution::default())?, f64::INFINITY); // A mathematical value.
let extreme = PersistenceDiagram::new(0, Coverage::Complete, vec![
    PersistenceInterval::new(0, -f64::MAX, IntervalEnd::Finite(f64::MAX))?,
])?;
assert!(matches!(wasserstein_2_euclidean(&extreme, &empty, 0),
    Err(Error::NumericalFailure { .. }))); // No scalar or partial value.
# Ok::<(), cocycle::Error>(())
```

There is no public cumulative budget for a sequence of terminals, memory quota,
deadline or pipeline-wide rollback. A failed call returns no partial result;
already completed outputs from previous calls remain valid. A future combined
operation must define and test its own shared budget rather than exposing
private reducer state or claiming that repeated calls provide that facility.
Work counts are not timing, bytes or a cross-algorithm performance score.

## Compare point sets and retain sample identities

The **Unreleased** `geometry::hausdorff_distance` compares nonempty point sets
directly in common Euclidean coordinates, without constructing topology or
computing persistence. The caller establishes matching physical units, frame
and normalization; equal dimensions alone do not establish those facts.

```rust
use cocycle::geometry::{PointCloudView, DissimilarityView,
    hausdorff_distance, hausdorff_distance_matrix};
let a = PointCloudView::new(&[0., 2.], 2, 1)?;
let b = PointCloudView::new(&[1., 3.], 2, 1)?;
let c = PointCloudView::new(&[4., 6.], 2, 1)?;
assert_eq!(hausdorff_distance(a, b)?, 1.);
let matrix = hausdorff_distance_matrix(&[("a", a), ("b", b), ("c", c)])?;
assert_eq!(matrix.sample_ids(), ["a", "b", "c"]); // Both axes, same order.
assert_eq!(matrix.values(), [1., 4., 3.]); // (1,0), (2,0), (2,1).
let nearest = (1..matrix.len()).min_by(|&i, &j|
    matrix.get(0, i).unwrap().total_cmp(&matrix.get(0, j).unwrap())).unwrap();
assert_eq!(matrix.sample_ids()[nearest], "b");
let borrowed = DissimilarityView::new(matrix.values(), matrix.len())?;
assert_eq!(borrowed.values().as_ptr(), matrix.values().as_ptr());
# Ok::<(), cocycle::Error>(())
```

The value is the worst nearest-neighbor Euclidean discrepancy, in coordinate
units. Point order and duplicate rows do not affect the underlying set. Distinct
sample IDs may therefore have zero distance. No correspondence, registration,
weights or automatic normalization is inferred. Hausdorff is sensitive to
outliers. It is a metric on finite sets, only a pseudometric on records with
duplicate rows, and does not certify positive definiteness of a derived kernel.
See the [mathematical definition](../reference/mathematics.md#25-euclidean-hausdorff-distance).

The matrix owns its exact, nonempty, unique UTF-8 IDs in input order and copies
no coordinates. It stores one finite value per unordered pair with an implicit
positive-zero diagonal. Invalid indices return `None`; zero samples and one
valid sample have distinct empty numeric buffers. Every sample is validated,
even when no pair will be evaluated. Empty point sets or unequal dimensions are
errors. Numeric failures, failed allocation, shape overflow, work exhaustion and
cancellation return no partial matrix. Missing observations must be handled
by the application before calling this interface, rather than encoded as zeros.

The `_with` variants use `Execution`; a matrix shares one fresh budget across
sample validation, ID comparisons/copies and all coordinate scans and writes.
Individual calls start separate budgets. Pre-cancellation precedes validation;
zero work rejects even a zero-sample matrix. Cancellation is checked during
coordinate and ID loops and around allocations, whose interiors are not
interruptible. Existing `PointCloudView` validation precedes the operation.
There is no memory quota or hard deadline.

Two directed scans use O(m*n*d) time and O(1) scalar scratch; no point-pair
matrix is materialized. For s samples the retained numeric storage is
s*(s-1)/2 f64 values, plus owned IDs and their vector. ID uniqueness uses
O(s*s) string comparisons. Explicit `DissimilarityView` conversion above scans
the output for validation but does not copy its values. The
[example](../../examples/hausdorff.rs) also supports `--fixtures` for independent
checks and `--measure` for exploratory resource snapshots; preparation is outside
its timers. There is no internal parallelism, cache or approximate-neighbor path.

## Compose direct comparison routes

Run `cargo run --locked --example comparison_routes` for three real compositions:

| Input and operation | Consumer and independently checked output |
| --- | --- |
| Three declared point sets in one Cartesian frame, Hausdorff sample matrix | Axes `a,b,c`, lower triangle `(1,3,2)` metres; nearest to `a` is `b`. Explicit `DissimilarityView` conversion borrows the numeric buffer. |
| Two actual Rips PH results from points `(0,2)` and `(0,3)`, context-checked H0 matching | Finite bars `[0,2)` and `[0,3)` plus one essential bar each; bottleneck-L-infinity, W1-L-infinity and W2-Euclidean each give 1 metre. A one-metre scalar threshold succeeds. |
| Two complete filtrations on the same triangle keys, application-owned L-infinity comparison | Sup difference 2 seconds, H1 bars `[3,6)` and `[3,7)`, and all three raw matching metrics 1 second. The bottleneck value satisfies the declared fixed-domain bound. |

The example's small `Sample` declaration retains ID, source snapshot, units,
alignment and normalization for both operands. Unknown or unequal declarations
fail before comparison. Different snapshots are permitted; source equality is
not a prerequisite for meaningful comparison. These strings remain caller
assertions, not authenticated provenance. No metadata registry is added to the
crate. The example preserves these fields alongside every exported scalar.

The Rips route borrows stored results for the context-checked calls, then moves
their diagrams with `into_diagram()`. Moving discards context explicitly; raw
calls afterwards agree because the application retains the declarations. PH is
not rerun for comparison or threshold analysis. Results remain usable after a
failed comparison. Unequal fields, missing degrees, incomplete coverage, cancelled
execution and work exhaustion are checked with actual computed results.

The filtration recipe compares `max_sigma |f(sigma)-g(sigma)|` only on identical
supplied simplex-key domains with known common units and label namespace. It
looks up keys through each existing owner, because filtration reordering changes
local IDs. Equal shapes or equal local ID numbers do not establish alignment.
The unique empty-domain pair has distance zero; mismatched domains and numeric
overflow fail. This is bounded application code, not a new general filtration
metric API or a comparison of arbitrary incidence arrays. General relabeling,
interleaving estimation and independent source registration remain unsupported.
The [fixed-domain stability argument](../reference/mathematics.md#8-stability-assumptions)
explains the checked bottleneck bound.

Supplied PH results have an unspecified parameter convention and are rejected
by the context-aware adapter, even when both say unspecified. This example
deliberately uses the existing raw matching API after establishing common field,
seconds and domain, rather than forging an edge-length context. No spectral or
projection comparison follows from equal array sizes or these examples.

Preparation, PH, comparison, explicit conversion and nearest-sample analysis
have separate timer fields in the output. `--measure` repeats the complete three
routes without serialization; it includes construction, PH and output handling.
The matrix owns only its IDs and values; neither geometry comparison nor the
filtered sup recipe clones the source. Diagram matching retains its existing
private allocations, while moving a diagram needs no PH recomputation. Each
library terminal starts its own execution budget; these compositions offer no
cumulative transaction, deadline or total-memory quota. The seven-cell
application filtration loop has no additional execution-control API.

Raw outputs and resource snapshots belong under `target/`. The independent
checker is a source-checkout developer tool, excluded from the Cargo package;
packaged examples and their tests remain runnable without Python dependencies.
