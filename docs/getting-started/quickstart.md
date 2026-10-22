# A square's persistent loop

[Documentation](../index.md) / Getting started

Use the [installation guide](installation.md) to create an application with
Cocycle 0.1.1 and Rust 1.91 or later. Replace `src/main.rs` with:

```rust
use cocycle::descriptors::betti_curve;
use cocycle::filtration::RipsBuilder;
use cocycle::geometry::PointCloudView;
use cocycle::persistence::PersistenceExt;

fn main() -> cocycle::Result<()> {
    let square = [0., 0., 1., 0., 1., 1., 0., 1.];
    let points = PointCloudView::new(&square, 4, 2)?;
    let result = RipsBuilder::from_points(points).persistence().compute()?;
    let curve = betti_curve(result.diagram(), 1, &[0., 1., 2.])?;
    assert_eq!(curve, [0, 1, 0]);
    println!("H1 Betti numbers at scales 0, 1, 2: {curve:?}");
    Ok(())
}
```

Run `cargo run`. The output is:

```text
H1 Betti numbers at scales 0, 1, 2: [0, 1, 0]
```

The row-major coordinate buffer contains four 2D points. `PointCloudView`
borrows that buffer; the computed result owns its diagram. The default request
computes H0 and H1 over F2.

At edge length 1, the square's sides form a loop. At $\sqrt{2}$, the diagonal
edges and triangles fill it, so the H1 interval is $[1,\sqrt{2})$. Births are
included and finite deaths are excluded. This is why the curve at scales
0, 1 and 2 is `[0, 1, 0]`.

To inspect the diagram and continue with matrix inputs or cutoffs, read the
[Rips guide](../guides/rips.md). A class still present at an incomplete cutoff
is right-censored; it is not automatically essential. For a fuller runnable
workflow, use the [square example](../../examples/square.rs) with
`cargo run --locked --example square` from a source checkout.
