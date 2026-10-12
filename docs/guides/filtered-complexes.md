# Supplied filtered complexes

[Documentation](../README.md) / Guides

Use `SimplicialComplex` for a supplied, face-closed collection of simplices with
finite filtration values. Use the small `FilteredComplex` trait when a different
storage format can supply cells in filtration order and oriented boundaries.
Both paths compute ordinary persistence over a selected prime field. They support
signed scales and unequal vertex birth times. Alpha geometry and a production
cubical container are not implemented by these interfaces.

General supplied simplicial filtrations and generic filtered-cell requests use
the same boundary reduction core. Zero-born concrete simplices retain the
existing coface path, as described below. The concrete adapter reads constructor-validated stored
incidence; the generic adapter validates the external cell contract. Concrete
source metadata and construction coverage remain attached to the result. Cycle
and cocycle requests retain the separate representative computation. Implicit
flag/Rips requests continue to use their specialized dispatch.

Algorithm authors can follow the [construction walkthrough](../development/complex-construction.md)
for a tested lower-star constructor example and its focused development command.

## Terminology and types

| Name | Meaning |
| --- | --- |
| `Simplex` | One nonempty simplex, its increasing vertex IDs and filtration value |
| `SimplicialComplex` | Concrete immutable storage of simplices, values, lookup and incidence |
| `FilteredComplex` | Read-only contract: ordered cells, dimension, value and oriented boundary |
| `SimplicialFiltration` | Stored simplicial complex plus construction provenance and certified source coverage |

These names describe mathematical objects and Rust responsibilities. They do not
claim a GUDHI storage implementation. GUDHI's [Filtered Complexes documentation][g-complex]
defines a simplicial complex and its filtration; its [FilteredComplex concept][g-contract]
defines algorithm requirements. Our trait is a local contract, not a literal port
of that concept. `SimplicialComplex` uses arrays, a lookup table and stored
codimension-one incidence, not a simplex tree. `FilteredSimplicialComplex` remains
a compatibility alias for `SimplicialComplex`.

## Construct and analyze simplices

Supply every nonempty face explicitly. Vertex IDs may have gaps. Input order is
arbitrary; the constructor sorts by value, dimension and the existing colex tie
rule. It rejects duplicate simplices, missing faces, invalid vertex sequences,
and faces appearing later than cofaces. It neither inserts missing faces nor
silently repairs filtration values. `new_with` adds cooperative execution controls.

```rust
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::diagram::IntervalEnd;
use cocycle::persistence::{PersistenceExt, RepresentativeRequest, RepresentativeSelection};

// The Alpha filtration of a unit equilateral triangle, supplied analytically.
// This example supplies the filtration; it does not run an Alpha constructor.
let complex = SimplicialComplex::new(vec![
    Simplex::new(vec![10], 0.)?,
    Simplex::new(vec![20], 0.)?,
    Simplex::new(vec![30], 0.)?,
    Simplex::new(vec![10, 20], 0.25)?,
    Simplex::new(vec![10, 30], 0.25)?,
    Simplex::new(vec![20, 30], 0.25)?,
    Simplex::new(vec![10, 20, 30], 1. / 3.)?,
])?;
let requests = [RepresentativeRequest::new(1, 0.3, RepresentativeSelection::Both)?];
let result = complex.persistence().representatives(&requests).compute()?;
let h1 = result.diagram().dimension(1)?.iter().next().unwrap();
assert_eq!(h1.birth(), 0.25);
assert_eq!(h1.end(), IntervalEnd::Finite(1. / 3.));
assert_eq!(complex.max_filtration_value(), Some(1. / 3.));
# Ok::<(), cocycle::Error>(())
```

The largest edge value is not the last event of a general filtration. All simplex
values contribute to its scale coverage. Source units are caller-defined for a
supplied complex; the library does not label these numbers as edge lengths.
`ComputationContext::filtration()` exposes the shared source context and its scale
convention. Rips-specific metadata lives in the `FiltrationSource::Rips` variant.

Negative vertex births, negative analysis cutoffs, negative representative queries
and signed Betti-curve grids are supported. Rips construction still validates
nonnegative distances and edge thresholds. The legacy `PersistenceOptions::new`
and direct Rips builder analysis retain their nonnegative cutoff contract;
supplied-complex and contextual explicit analysis use signed filtration cutoffs.

## Implement another storage representation

The four-method trait requires neither vertex tuples nor contiguous integer IDs.
Here a single square is a 2-cell with four oriented boundary edges; it is not
converted to triangles.

```rust
use cocycle::complex::FilteredComplex;
use cocycle::persistence::PersistenceBuilder;
use cocycle::diagram::IntervalEnd;

struct Square;
impl FilteredComplex for Square {
    type CellId = usize;
    fn cells(&self) -> impl Iterator<Item = usize> + '_ { 0..9 }
    fn dimension(&self, cell: usize) -> usize {
        match cell { 0..=3 => 0, 4..=7 => 1, _ => 2 }
    }
    fn filtration_value(&self, cell: usize) -> f64 { self.dimension(cell) as f64 }
    fn boundary(&self, cell: usize) -> impl Iterator<Item = (usize, i32)> + '_ {
        const BOUNDARY: [&[(usize, i32)]; 9] = [
            &[], &[], &[], &[],
            &[(0, -1), (1, 1)], &[(1, -1), (2, 1)],
            &[(2, -1), (3, 1)], &[(3, -1), (0, 1)],
            &[(4, 1), (5, 1), (6, 1), (7, 1)],
        ];
        BOUNDARY[cell].iter().copied()
    }
}
let result = PersistenceBuilder::from_complex(&Square).compute()?;
let h1 = result.diagram().dimension(1)?.iter().next().unwrap();
assert_eq!((h1.birth(), h1.end()), (1., IntervalEnd::Finite(2.)));
# Ok::<(), cocycle::Error>(())
```

Cell handles remain stable during a computation. Values are finite and ordered;
boundary cells precede cofaces, including ties. Boundaries decrease dimension by
one and satisfy boundary-of-boundary zero over the integers. Coefficients are
signed `i32` incidence numbers, reduced into the selected prime field. Accessors
and iteration borrow stable storage and must honor this mathematical contract.

The generic computation checks unique IDs, finite ordered values, and the
retained skeleton's boundary references, dimensions and boundary-square zero in
the selected field. This last check does not certify integer validity or validity
in every other field. No unsafe operation relies on the trait. Metadata is read
for all cells to establish coverage; only the requested q+1 skeleton's boundary
matrix is retained. Diagram-only reduction omits basis transformations.

`from_complex` returns an ordinary configurable `PersistenceBuilder` and an owned
`PersistenceResult`, but does not support vertex-labelled representatives: the
cell contract supplies no vertex lists. A nonempty representative request returns
an error. For simplicial representatives use `complex.persistence()` instead.
Iteration/validation/reduction share the execution budget; controls cannot
interrupt the interior of a custom accessor or iterator call.

## Read stages and their inclusions

`SimplicialComplex::stage(scale)` and `SimplicialFiltration::stage(scale)` borrow
the closed sublevel set of stored simplices with value <= scale. Every tie is
included. Selection uses a binary search and retains the source's simplex IDs,
vertex orientation and signed boundaries; it allocates no simplex or map storage.
Finite negative scales are valid, including empty stages before any cell enters.

```rust
use cocycle::complex::{FilteredComplex, Simplex, SimplicialComplex};
use cocycle::diagram::{Coverage, IntervalEnd};
use cocycle::persistence::PersistenceExt;
use cocycle::Error;

let complex = SimplicialComplex::new(vec![
    Simplex::new(vec![10], -2.)?, Simplex::new(vec![20], -1.)?,
    Simplex::new(vec![30], 0.)?, Simplex::new(vec![10, 20], 1.)?,
    Simplex::new(vec![10, 30], 1.)?, Simplex::new(vec![20, 30], 1.)?,
    Simplex::new(vec![10, 20, 30], 2.)?,
])?;
let stages = [-1., 0., 1., 2.].map(|scale| complex.stage(scale).unwrap());
assert_eq!(stages.map(|stage| stage.len()), [2, 3, 6, 7]);
let ring = stages[2];
let filled = stages[3];
for (id, image) in ring.inclusion_into(&filled)? {
    assert_eq!(id, image);
    assert_eq!(ring.boundary(id), filled.boundary(image));
}
let edges: Vec<_> = ring.cells()
    .filter(|id| FilteredComplex::dimension(&ring, *id) == 1).collect();
assert_eq!(edges.len(), 3);
let result = ring.persistence().compute()?;
assert_eq!(result.diagram().coverage(), Coverage::Through(1.));
assert_eq!(result.diagram().dimension(1)?.iter().next().unwrap().end(),
    IntervalEnd::RightCensored { through: 1. });
assert!(matches!(filled.inclusion_into(&ring), Err(Error::InvalidParameter { .. })));
let another = complex.clone();
assert!(matches!(ring.inclusion_into(&another.stage(2.)?),
    Err(Error::InvalidParameter { .. })));
# Ok::<(), cocycle::Error>(())
```

The [construction example](../../examples/complex_construction.rs) prints four
lower-star stages and their inclusions before any persistence computation.
Consumers can independently traverse `cells()` and signed boundaries, or borrow
`simplices()` and `simplex(id)`. Degree-specific bases and external matrices are
explicit consumer-owned conversions. Absent cells return `None`; bare IDs still
require the caller to use the correct owner.

A certified stage exposes its original `SimplicialFiltration` through
`source_filtration()`, retaining context, approximation, scale coverage and
construction-dimension sufficiency. Selecting above `Coverage::Through(t)` is an
error. A stage may be read from an insufficient skeleton, but its contextual
`.persistence()` still rejects an unsupported homology dimension. It never
completes the skeleton or changes sparse blockers. An analysis cutoff must be
finite and no greater than the stage scale; omitting it selects that scale.
Analysis and representative computation retain the existing execution controls.

`stage.persistence()` analyzes the original filtration through the stage scale;
later unobserved deaths remain censored. Explicitly using
`PersistenceBuilder::from_complex(&stage)` or `stage.source_complex()` selects
stored topology as a new supplied source, discarding the relation to later events
and construction certificates. A ring can then have an essential H1 class even
if the original filtration later fills it. This is an explicit change of source.

`inclusion_into` checks that both stages borrow the same source and certificate
interpretation and that their scales do not decrease. It yields identity pairs
without copying or renumbering; its coefficients are +1 over the integers.
Cloned owners, reverse scales and certified/bare interpretation mixtures are
rejected, even if the integer IDs or cell sets happen to agree. General chain
maps and multi-parameter filtrations are outside this first capability. The
[source comparison and decisions](../research/filtration-stages.md) explain those
limits, units, ordering and independent validation.

## Supplied topology versus certified construction

A supplied complex is the mathematical source in its entirety. A triangle's
unfilled boundary has an essential H1 class. In contrast, an expansion of a Rips
source only through dimension one may be missing a triangle that kills that class.
`SimplicialFiltration` retains the source's scale coverage and dimension sufficiency
certificate and rejects such insufficient Rips skeletons.

Calling `filtration.complex().persistence()` explicitly chooses the stored complex
as the source and discards the larger source's certificate. Use
`filtration.persistence()` to retain the construction's meaning. Ordinary analysis
does not silently switch between these interpretations.

## Extension boundary

An Alpha constructor can compute a triangulation and filtration values in private
mutable working storage, then validate/freeze a `SimplicialComplex`. Adding a public
mutable simplex tree is not a prerequisite. A future Alpha source must also declare
scale units, vertex mapping and construction coverage; an arbitrary supplied
complex cannot manufacture those certificates.

Rips implicit coface algorithms retain a private, explicitly zero-born contract.
They are distinct from the public boundary contract. Diagram-only explicit
simplicial analysis reuses union-find and coface clearing when all vertices
are born at zero and the query cutoff is nonnegative or absent. The decision
checks actual simplex values; it also supports non-flag topology and delayed
higher-simplex values through stored cofaces. Other inputs use boundary reduction,
and representative requests retain their separate transformation work. Generic
`PersistenceBuilder::from_complex` always uses the boundary contract. Native
workers preserve their protocol; performance claims require measurements bound
to the measured source commit. The [integration report](../../benches/reports/critical-sets.md)
records the current-main gate that rejected automatic boundary rerouting.

[g-complex]: https://gudhi.inria.fr/doc/latest/group__simplex__tree.html
[g-contract]: https://gudhi.inria.fr/doc/latest/struct_filtered_complex.html
