# Borrow inputs, read results and control composed calls

[Documentation](../README.md) / Usage

These examples use existing public interfaces. They connect modeling, stage
selection, persistence, result access, comparison and analysis without a common
source trait or a universal result container. The
[interface compatibility report](../design/interface-compatibility.md) records
the decisions and the remaining spectral/operator boundaries.

## Read a signed stage for a matrix consumer

Borrow a frozen complex and select cells with value <= scale. Preserve its
degree-specific basis order, owner-local IDs and increasing-vertex orientation.
The explicit conversion below returns owned dense matrices only because the
small external consumer needs them. It is tutorial code, not a new public API.

```rust
use cocycle::complex::{Simplex, SimplexId, SimplicialComplex};

let complex = SimplicialComplex::new(vec![
    Simplex::new(vec![10], 0.)?, Simplex::new(vec![20], 0.)?,
    Simplex::new(vec![30], 0.)?, Simplex::new(vec![10, 20], 1.)?,
    Simplex::new(vec![10, 30], 1.)?, Simplex::new(vec![20, 30], 1.)?,
    Simplex::new(vec![10, 20, 30], 2.)?,
])?;
let source = &complex;
let scale = 1.;
let bases: [Vec<SimplexId>; 3] = std::array::from_fn(|degree| {
    source.simplices().iter()
        .filter(|s| s.dimension() == degree && s.value() <= scale)
        .map(|s| source.find(s.vertices()).unwrap()).collect()
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
# Ok::<(), cocycle::Error>(())
```

The selected degree bases allocate O(number of selected cells) handles and
borrow the source's vertex labels. The dense export allocates one entry per
matrix coordinate and uses linear basis lookups; it is a deliberately small
example, not a recommended large-complex representation. Sparse consumers can
iterate signed incidence instead. The source is unchanged, and no persistence
request is needed for this conversion.

For an ordinary nonaugmented H0 consumer, the previous space is empty and its
boundary has shape 0-by-n. Export that shape explicitly. Simplex IDs belong to
this owner; equal integer indices from another complex are not matching cells.
Signed simplicial incidence can be reduced into F2; arbitrary modular input
cannot be lifted into a real chain complex by replacing each bit with 0 or 1.

A single-scale graph consumer can similarly borrow `WeightedGraph::edges()`
and derive vertex-edge incidence. Its interpretation as a one-dimensional
complex is explicit: a graph triangle is not automatically a filled clique.
Filtration values, spectral Grams and F2 coordinate costs remain separate.
This prepares inputs; it implements no Laplacian, Dirac or F2 projection solver.

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
