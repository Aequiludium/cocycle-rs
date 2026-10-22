# Diagram distances

[Documentation](../index.md) / Guides

Compare a computed homology dimension of two **complete** persistence diagrams.
The Rust implementation supports three matching metrics:

| Function in `cocycle::diagram_distances` | Ground cost | Objective |
| --- | --- | --- |
| `bottleneck_distance` | L-infinity | Minimum largest matching cost |
| `wasserstein_1_infinity` | L-infinity | Minimum sum of matching costs |
| `wasserstein_2_euclidean` | Euclidean | Square root of minimum sum of squared costs |

## Match a finite interval to the diagonal

An empty computed dimension is a valid diagram. Comparing one interval `[0, 2)`
with an empty H1 diagram matches that interval to the diagonal:

```rust
use cocycle::diagram::{Coverage, IntervalEnd, PersistenceDiagram, PersistenceInterval};
use cocycle::diagram_distances::{
    bottleneck_distance, wasserstein_1_infinity, wasserstein_2_euclidean,
};

fn main() -> cocycle::Result<()> {
    let first = PersistenceDiagram::new(
        1,
        Coverage::Complete,
        vec![PersistenceInterval::new(1, 0., IntervalEnd::Finite(2.))?],
    )?;
    let second = PersistenceDiagram::new(1, Coverage::Complete, vec![])?;
    assert_eq!(bottleneck_distance(&first, &second, 1)?, 1.);
    assert_eq!(wasserstein_1_infinity(&first, &second, 1)?, 1.);
    assert!((wasserstein_2_euclidean(&first, &second, 1)? - 2.0_f64.sqrt()).abs() < 1e-14);
    Ok(())
}
```

The L-infinity diagonal cost is half the lifetime. The Euclidean diagonal cost
is the lifetime divided by $\sqrt{2}$. The
[runnable distance example](../../examples/diagram_distances.rs) also demonstrates
bounded execution with `bottleneck_distance_with`; run
`cargo run --locked --example diagram_distances` from a source checkout.

## Establish comparability

Repeated intervals retain their multiplicity. Essential intervals match only
essential intervals; unequal essential counts produce positive infinity.
Uncomputed dimensions are errors. All `Coverage::Through` diagrams are rejected,
even if they currently contain no censored interval, because their future births
and deaths are unknown.

Raw-diagram calls leave comparable scale units to the caller. The `_results`
functions additionally check equal coefficient characteristics and declared
edge-length scales; they still cannot certify physical units or normalization.
Do not replace a censored endpoint with a finite death to make a comparison pass.

See the [matching specification](../reference/mathematics.md#16-diagram-matching-distances)
for essential matching, numerical behavior, errors and resource controls, and
the [API reference](../reference/api.md) for function contracts.
