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

`PathVector` stores only nonzero canonical residues. Omega columns and homology
cycles use the result's allowed-path order; boundary columns use the preceding
Omega basis. The chain includes degrees 0 through d+1, while representatives and
Betti numbers cover 0 through d. A missing accessor degree is outside that
window. Bases are deterministic, without a canonical, minimum or metric claim.

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
