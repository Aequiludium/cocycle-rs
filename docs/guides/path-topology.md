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
