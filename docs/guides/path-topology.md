# Directed path topology

[Documentation](../README.md) / User guides

This Unreleased API keeps directed edges and nonadjacent vertex repetitions.
It is independent of the undirected `WeightedGraph` and simplicial storage.

```rust
use cocycle::complex::{DirectedGraph, PathBoundaryMode, PathComplex, path_boundary};

let graph = DirectedGraph::new(3, vec![(1, 0), (0, 1)])?;
let paths = graph.allowed_paths(2, 7)?;
assert_eq!(paths.paths(2), [vec![0, 1, 0], vec![1, 0, 1]]);
assert_eq!(paths.paths(0).len(), 3); // Includes the isolated vertex.
assert_eq!(paths.truncated_at(), Some(2));
let explicit = PathComplex::new(2, &[vec![0], vec![1], vec![0, 1]], 3)?;
assert!(explicit.paths(2).is_empty()); // No higher paths are inferred.
let boundary = path_boundary(&[0, 1, 0], PathBoundaryMode::NonRegular)?;
assert!(boundary.contains(&(vec![0, 0], -1)));
# Ok::<(), cocycle::Error>(())
```

The positive path cap counts all retained paths, including singletons. Every
explicit path set must contain all singletons and both endpoint truncations.
Missing internal faces remain visible in the ambient integer boundary. Regular
mode discards nonregular faces; NonRegular mode retains them. Both use the
non-reduced convention, with zero singleton boundary.

The `_with` variants accept `execution::Execution`. Each terminal starts one
budget and returns an error on exhaustion or cancellation, without partial
output. Graph construction, explicit validation, enumeration and ambient boundary
construction count their own coordinate/candidate visits. Passing the same
Execution to separate terminals does not share counters. Sorting and individual
container operations are cooperative checkpoints, not interruptible operations.
Path caps and work units are not byte quotas. Vector growth uses fallible
reservations; B-tree node allocation can still abort on allocation failure.

Allowed paths need not form a chain space: their boundary may leave the allowed
set. These sources do not implement the integer `FilteredComplex` contract.
Computing H_d requires field-specific Omega spaces through d+1; a graph cutoff
must not silently certify absent higher spaces. See the
[mathematical contract](../reference/mathematics.md#25-directed-allowed-path-sources).

## Exact static homology

`persistence::compute_path_homology` constructs Omega over the selected
`algebra::PrimeField`. It returns an owned `PathHomology` with read-only field,
model, path identities, Omega bases, coordinate boundaries, ranks and homology
representatives. No modular coefficient is reinterpreted as integer incidence.

```rust
use cocycle::algebra::PrimeField;
use cocycle::complex::{DirectedGraph, PathBoundaryMode};
use cocycle::persistence::compute_path_homology;
let source = DirectedGraph::new(3, vec![(0, 1), (1, 2), (2, 0)])?
    .allowed_paths(2, 100)?;
let result = compute_path_homology(&source, 1, PrimeField::new(4_294_967_291)?,
    PathBoundaryMode::Regular, 10_000)?;
drop(source);
assert_eq!(result.betti(), [1, 1]);
assert_eq!(result.omega_bases(2).unwrap().len(), 0);
# Ok::<(), cocycle::Error>(())
```

`PathVector` is a read-only borrowed coordinate view with `terms()` and
`coefficient()`. The sparse container is private; only nonzero canonical residues
are iterated. Omega columns and homology
cycles use the result's allowed-path order; boundary columns use the preceding
Omega basis. The chain includes degrees 0 through d+1, while representatives and
Betti numbers cover 0 through d. A missing accessor degree is outside that
window. Basis choices and their order are result-local and may change between versions;
no pivot/block strategy, canonical basis, minimum basis or metric is promised.

The positive `max_entries` bounds layer/column counts, each sparse matrix and
each elimination's retained plus temporary entries, including transformations.
It does not bound all simultaneously live objects, B-tree overhead, or bytes.
`compute_path_homology_with` shares one Execution budget across environment
boundaries, kernels, coordinate solves, homology representatives and coordinate
copying. Every arithmetic/coordinate/column visit is charged; source validation
occurred before this terminal. Exhaustion returns an error without a result.

The independent dense-row certificate tests include both models on all 64
three-vertex graphs, the largest u32 prime, a higher-dimensional join and a
characteristic-sensitive Omega example. These are bounded correctness checks,
not evidence for a speedup or arbitrary-size workloads. Run
`cargo run --locked --example path_homology` for an owned-result workflow.

## Directed graph persistence

```rust
use cocycle::algebra::PrimeField;
use cocycle::complex::PathBoundaryMode;
use cocycle::persistence::compute_path_persistence;
let edges = [(0,1,1.0), (1,2,1.0), (2,0,1.0), (0,2,2.0)];
let result = compute_path_persistence(3, &edges, 1, PrimeField::new(3)?,
    PathBoundaryMode::Regular, 100, 10_000)?;
assert_eq!(cocycle::descriptors::betti_curve(result.diagram(), 1, &[0.0,1.0,2.0])?,
    [0,1,0]);
# Ok::<(), cocycle::Error>(())
```

Every vertex is born at zero. Finite nonnegative edge births control appearance
only; each path is born at its largest edge birth. Equal represented f64 values
are simultaneous, including signed zero. All graph events are included: coverage
is Complete and the computed domain is 0 through d. Zero-length bars are omitted
by the existing ordinary interval assembler. No cutoff or unequal vertex births
are offered by this operation.

`PathPersistence` composes canonical `PersistenceData`, implements
`AsRef<PersistenceData>`, and owns the compatible Omega certificate. Births index
`times()`; paths are final-graph coordinates. A basis prefix at an event spans
that event's Omega. This certificate does not claim interval-specific simplex
representatives. Source context distinguishes both path models, records the
selected prime field and declares Unspecified scale units. Context-aware diagram
distances consequently retain their existing rejection of unsupported units.
Raw-diagram comparisons remain caller-interpreted operations.

The `_with` entry shares one Execution counter across validation, graph/path
enumeration, all static stage computations, basis extensions and the existing
ordinary reducer. `max_paths` bounds each stage. Positive `max_entries` bounds
construction matrices, events/columns and Omega elimination workspaces. Ordinary
reduction fill-in follows the existing reducer's work controls rather than this
entry cap; neither is an aggregate memory quota. This recomputing baseline makes
no scalability or performance claim.

The Unreleased exhaustive `FiltrationKind` enum gains RegularDirectedPaths and
NonRegularDirectedPaths. External exhaustive matches must add these variants.
`FiltrationSource` is already non-exhaustive. Existing Rips, flag, simplicial and
filtered-cell dispatch are unchanged. Run
`cargo run --locked --example path_persistence` for result composition.
