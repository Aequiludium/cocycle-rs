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

The default terminal returns the existing `PersistenceResult`, retaining only
the ordinary diagram and context. `compute_path_persistence_certificate` and its
`_with` variant explicitly request an owned `PathPersistenceCertificate`. It
composes `PersistenceData` and implements `AsRef<PersistenceData>`. Births index
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

## Algebraic duals, projections and minimum bases

`PathHomology::projection(max_entries)` returns annotations A, projections P,
and homotopies h in the stored Omega coordinates. The construction checks
P squared = P, B P = 0, A B = 0, identity evaluation on the stored homology
basis, and I-P = B h+h B through degree d. The last h lands in Omega_(d+1).
`cohomology(max_entries)` transposes the annotations into cocycles in Omega dual
coordinates. Both have `_with` execution variants with one terminal-wide budget.
They are algebraic splittings and linear duals, without an orthogonal/harmonic
interpretation, cup product or canonical extension to all path coordinates.

`compute_minimum_path_basis` takes a DirectedGraph, one nonnegative u64 cost per
lexicographic directed edge, the model, and positive path, sparse-entry, candidate
and search-state caps. It computes a global minimum F2 H1 basis by exhaustively
enumerating elementary cycles of the underlying edge multigraph and greedily
selecting independent homology annotations. Reverse edges keep different
coordinates; their two-edge cycles remain eligible. Nonnegative costs, including
zero, justify elementary-cycle candidates. Costs and their total use checked
u128 arithmetic. The `_with` variant shares work controls across enumeration,
homology, projection, cycle search, selection and certificate checks.

This search is exponential. Candidate/support/frontier and search exhaustion are
errors; no partial basis is reported as optimal. Caps are independent counts,
not total-memory guarantees. Sparse projection and search use the same exact
field arithmetic; B-tree allocations retain the standard allocation limitation.

Omega construction groups forbidden rows by path endpoints before computing
independent kernels. Endpoint closure guarantees that only internal deletions
can be forbidden. The independent dense tests check the resulting complete spaces
in both models, including characteristic-sensitive and high-prime examples.
This grouping has no public tuning switch and no claimed speedup. The minimum
basis is independently checked against all cycle subsets of three-vertex graphs.

## Direct H0 and H1 presentation

`compute_low_dimensional_path_persistence` takes weighted edges, a PrimeField,
the model, a path cap and an entry cap. It constructs generators of partial
Omega2 directly from bigons, triangles and anchored quadrangle differences,
including the full-model returning-path constraints. Earliest completed routes
provide anchors, so the generators span the one-boundaries at every event.
They need not be independent and are not an Omega2 basis.
The terminal returns the existing `PersistenceResult`; dependent generators,
anchors and traversal counters are private working state.

The upstream adapter feeds this dependent boundary presentation to the existing
ordinary reducer, requests only H0/H1 and returns canonical PersistenceData with
the same path context as the generic route. Artificial H2 of the presentation is
not computed or exposed. The research implementation's separate forest/chord
reducer is not duplicated in production. The formulas follow
[Dey, Li and Wang, Theorem 3.1](https://arxiv.org/abs/2001.09549v1); this code does
not implement its arboricity enumeration bound.

The path cap includes vertices, edges and every visited two-path. The entry cap
bounds generators, presentation columns and initial incidence terms. The `_with`
variant shares one work counter through construction and ordinary reduction.
Ordinary fill-in follows that work control, without a byte/aggregate memory cap.
Full diagrams match the independently checked generic route for the tested
small graphs, ties and prime fields; no speedup or large-network result is claimed.

## Numerical real path operators

The separate `hodge` entries recompute Omega from integer ambient boundaries
over f64. They do not lift a finite-field basis or implement integer
FilteredComplex. Allowed paths have their standard orthonormal metric; the
returned Omega bases use two reorthogonalization passes.

```rust
use cocycle::complex::{DirectedGraph, PathBoundaryMode};
use cocycle::hodge::compute_real_path_chain;
let source = DirectedGraph::new(2, vec![(0,1)])?.allowed_paths(2, 10)?;
let chain = compute_real_path_chain(&source, 1, PathBoundaryMode::Regular, 1000, 1e-10)?;
let operator = chain.laplacian(0)?;
drop(chain);
assert_eq!(operator.apply(&[1.0,0.0])?, [1.0,-1.0]);
assert!(operator.boundary_residual() < 1e-8);
# Ok::<(), cocycle::Error>(())
```

`PathRealChain` owns allowed path identities, orthonormal Omega columns and
coordinate boundaries through d+1. A graph truncation must include that upper
degree. The finite relative tolerance lies strictly in (0,1e-3). Pivots within
the ambiguity band around the tolerance, nonfinite arithmetic and failed kernel,
coordinate-boundary or boundary-square residuals return errors. Results retain
the tolerance and largest measured absolute boundary residual; this is numerical
evidence, without a certified exact rank or arbitrary-condition guarantee.

`compute_persistent_path_laplacian` takes weighted edges, degree, an event-index
pair s<=t, model, (path cap,dense cap) and tolerance. Events are sorted zero plus
edge births. The upper domain consists of target Omega_(d+1) vectors whose full
boundary lies in source Omega_d, with its induced real orthonormal metric.
`PathRealLaplacian` owns final-graph path identities, the embedded source basis
and rectangular D/U factors, and applies D transpose D + U U transpose. The
library API builds neither square matrices nor spectra and does not use the
generic Hodge/Dirac integer-source adapter.

The positive dense cap bounds each matrix/vector, layer count, ambient word
workspace dimension squared and kernel width squared, including zero entries.
It differs from the exact sparse cap and does not aggregate live memory.
Vector reservations are fallible where supported; B-tree allocations can abort.
Construction `_with` terminals share work across requested chains, kernels
and residuals. Snapshot and action `_with` terminals each start a new counter.
Sorting/container operations are cooperative checkpoints; numerical arithmetic
and coordinate visits are charged inside the main loops.

`tests/path_real.rs` checks ambient equations, orthonormality, symmetric
positive-semidefinite action, hand-derived operators, the characteristic-sensitive
real/finite-field distinction and invalid/resource cases. The standalone worker
`tools/reference/path_real.rs` exposes factors and sampled action columns for
independent rational/SVD verification; it is not a production spectrum interface
or performance benchmark. Run `cargo run --locked --example path_laplacian`.

The Rust series does not include the research Python Grounded, Weighted, Mayer,
cup-product or large-network workstreams, nor a sparse path eigensolver.
