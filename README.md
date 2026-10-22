![Cocycle — topology, computed in Rust](assets/banner.svg)

[![CI](https://github.com/Aequiludium/cocycle-rs/actions/workflows/ci.yml/badge.svg?branch=main)](https://github.com/Aequiludium/cocycle-rs/actions/workflows/ci.yml)
[![Documentation](https://github.com/Aequiludium/cocycle-rs/actions/workflows/pages.yml/badge.svg?branch=main)](https://aequiludium.github.io/cocycle-rs/)
[![Rust 1.91+](https://img.shields.io/badge/Rust-1.91%2B-102D32?style=flat-square)](Cargo.toml)
[![License: MIT](https://img.shields.io/badge/license-MIT-087D70?style=flat-square)](LICENSE)

**Find structure across scales.** Cocycle turns point clouds and pairwise
dissimilarities into persistence diagrams: which components and loops appear,
and how long they last. A general-purpose TDA library in pure Rust, with no runtime
dependencies and no unsafe code.

- **Rips persistence** — dimension-generic ordinary persistence over prime fields,
  with a specialized F₂ H₀/H₁ path. Euclidean points or borrowed matrix layouts.
- **Inspectable construction** — exact threshold graphs, supplied weighted flag
  filtrations, and frozen simplicial complexes with boundary/cofacet queries.
- **Supplied complexes** — validated simplicial storage and a small filtered-cell
  trait for boundary reduction, including signed scales and unequal vertex births.
- **Real Hodge analysis** — sparse action, explicit matrices, bounded dense
  spectra and matrix-free partial eigenpairs under the standard inner product (Unreleased).
- **Unreleased stage access** — borrowed sublevel stages and checked inclusions
  between increasing scales of the same frozen source, preserving IDs and signs.
- **Sparse approximation** — deterministic sampling, modified edges and higher-simplex
  blockers, with explicit metric hypotheses and approximation provenance.
- **Optional representatives** — owned persistent cycle bases and query-scale dual
  cocycles, associated with intervals and original vertex IDs.
- **Meaningful results** — owned diagrams preserve multiplicity and distinguish
  finite deaths, essential classes, and right-censored intervals. Explicit computed
  dimension sets distinguish empty results from uncomputed dimensions; common
  diagram/context data can be shared by specialized result types.
- **Useful summaries** — finite lifetimes, persistence entropy in nats, and Betti
  curves, computed directly from a diagram.
- **Unreleased Dirac operators** — sparse real action and bounded matrix export
  on signed chains, selected signed boundary modes, explicit degree truncation
  and source-local coordinates.
- **Diagram distances** — exact bottleneck (L-infinity), W1 (L-infinity) and
  W2 (Euclidean) matching of complete diagrams, preserving essential multiplicity.
- **Unreleased critical sets** — opt-in sparse F₂ Big Steps endpoint proposals
  with bounded Partial U; updates and differentiation remain caller-owned.
- **Unreleased persistent Laplacians** — bounded real stage-pair Schur operators,
  action, complete numerical spectra and fixed-endpoint reverse elimination.
  See the [usage and limits](docs/guides/persistent-laplacian.md).

![A square's Rips filtration: points, a loop at edge length 1, and filled triangles at sqrt(2).](assets/filtration.svg)

## Quick start

**Rust 1.91+ · MIT.** Add Cocycle to your application's `Cargo.toml`:

```toml
[dependencies]
cocycle = "0.1.1"
```

```rust
use cocycle::descriptors::betti_curve;
use cocycle::geometry::PointCloudView;
use cocycle::filtration::RipsBuilder;
use cocycle::persistence::PersistenceExt;

fn main() -> cocycle::Result<()> {
    let square = [0., 0., 1., 0., 1., 1., 0., 1.];
    let points = PointCloudView::new(&square, 4, 2)?;
    let result = RipsBuilder::from_points(points).persistence().compute()?;
    let diagram = result.diagram();
    assert_eq!(betti_curve(diagram, 1, &[0., 1., 2.])?, [0, 1, 0]);
    Ok(())
}
```

The square's H₁ interval is `[1, sqrt(2))`. Scales are **edge lengths**; a class
surviving an incomplete cutoff is censored, not dead. See the
[user guide](docs/guides/rips.md) for cutoffs, input layouts, and result semantics.
The [0.1.1 release notes](CHANGELOG.md#011) describe the latest maintenance fixes;
the [initial release notes](CHANGELOG.md#010) describe the supported scope,
limitations and migration from development Git revisions. See the
[version policy and release procedure](CONTRIBUTING.md#release-procedure) for
compatibility and publishing. The dependency above selects 0.1.1;
this checkout's Unreleased additions are identified in the
[version status](docs/README.md#version-and-capability-status).

**Unreleased stage access.** This example requires the candidate API. A stage
borrows
all stored simplices with value <= scale, including ties. An inclusion preserves
owner-local IDs and signed incidence; equal IDs from separate owners are not a map.

```rust
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::diagram::{Coverage, IntervalEnd};
use cocycle::persistence::PersistenceExt;

let complex = SimplicialComplex::new(vec![
    Simplex::new(vec![10], -1.)?, Simplex::new(vec![20], -1.)?,
    Simplex::new(vec![30], -1.)?, Simplex::new(vec![10, 20], 1.)?,
    Simplex::new(vec![10, 30], 1.)?, Simplex::new(vec![20, 30], 1.)?,
    Simplex::new(vec![10, 20, 30], 2.)?,
])?;
let ring = complex.stage(1.)?;
let filled = complex.stage(2.)?;
assert_eq!((ring.len(), filled.len()), (6, 7));
for (id, image) in ring.inclusion_into(&filled)? {
    assert_eq!(id, image);
    assert_eq!(ring.boundary(id), filled.boundary(image));
}
let result = ring.persistence().max_homology_dimension(1).compute()?;
assert_eq!(result.diagram().coverage(), Coverage::Through(1.));
assert_eq!(result.diagram().dimension(1)?.iter().next().unwrap().end(),
    IntervalEnd::RightCensored { through: 1. });
# Ok::<(), cocycle::Error>(())
```

Reading a stage does not compute persistence. Analyzing the ring as a stage of
its source retains the cutoff; its surviving loop is censored. Re-importing its
cells as a standalone supplied complex makes a different coverage claim. See the
[stage guide](docs/guides/filtered-complexes.md#read-stages-and-their-inclusions)
for certified sources, errors and the distinction from general chain maps.

## Explore

[Documentation website](https://aequiludium.github.io/cocycle-rs/) ·
[API reference](https://docs.rs/cocycle) · [Documentation map](docs/README.md) ·
[Filtered complexes](docs/guides/filtered-complexes.md) · [Rips guide](docs/guides/rips.md) · [Graph construction](docs/guides/rips-construction.md) · [Mathematics](docs/reference/mathematics.md) ·
[Architecture](docs/development/architecture.md) · [Benchmarks](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/README.md) ·
[Rips comparison](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/reports/rips-comparison.md) · [Contributing](CONTRIBUTING.md) · [Roadmap](docs/design/roadmap.md)

Contributing an algorithm? Start with the [research-area navigation](docs/development/algorithm-contributions.md)
and the [diagram-analysis](docs/development/diagram-analysis.md) or
[complex-construction](docs/development/complex-construction.md) walkthrough.
Run `cargo run --locked --example diagram_analysis` to analyze supplied intervals
without constructing a complex.
Run `cargo run --locked --example complex_construction` to construct a lower-star
filtration, inspect its simplices and compute persistence using existing APIs.
Run `cargo run --locked --example dirac` for real Dirac action on a filled triangle;
see the [Dirac guide](docs/guides/dirac.md).
The unreleased source API also provides [F2 homology operators](docs/guides/homology-operators.md):
run `cargo run --locked --example homology_operators` for class queries, finite
transport and bounded exact integer geometry queries on historical features.
Run `cargo run --locked --example hodge_laplacian` for standard real Hodge
matrices and spectra; see the [Hodge guide](docs/guides/hodge-laplacian.md).

Read the [published API reference](https://docs.rs/cocycle), or build it locally
with `cargo doc --no-deps --open`. From a repository
checkout, run the example with
`cargo run --locked --example square`. Dimension-generic prime-field persistence and
explicit complex queries are demonstrated by `cargo run --example rips_sphere`.
[Prime fields and representative bases](docs/guides/rips-representatives.md) are
demonstrated by `cargo run --example rips_representatives`.
[Sparse Rips approximation](docs/guides/sparse-rips.md), including metric hypotheses
and sampling provenance, is demonstrated by `cargo run --example sparse_rips`.

To compare two complete diagrams, run
`cargo run --locked --example diagram_distances`. The [distance example](examples/diagram_distances.rs)
demonstrates all three metrics; see the [matching contract](docs/reference/mathematics.md#16-diagram-matching-distances)
for ground costs, essential intervals and error semantics.

Phase 4 connects six responsibilities: data, filtration, computation, results,
comparison and analysis. Filtration selects scales when needed; computation
produces concrete results that comparison and analysis can borrow independently.
See the [architecture](docs/development/architecture.md#six-responsibilities-and-data-routes),
[roadmap](docs/design/roadmap.md#phase-4-priorities) and
[interface guide](docs/guides/interface-contracts.md). Standard real Hodge,
Dirac and bounded persistent Laplacian operators are available in this source
checkout. General weighted models, certified global partial spectra, scalable
persistent spectra, optimal F₂ projections, non-diagram similarity and production
learning integration remain future contracts. The
[route report](docs/research/route-integration.md) retains the earlier handoff.

Work and memory depend on the input and reduction fill-in. Rips comparisons
use native GUDHI and upstream Ripser C++; see the maintained
[Rips comparison](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/reports/rips-comparison.md) for tested scope, correctness
and performance observations. Generated measurements and logs stay outside source Git.

[Issue tracker](https://github.com/Aequiludium/cocycle-rs/issues) ·
Code and original artwork are [MIT licensed](LICENSE).
