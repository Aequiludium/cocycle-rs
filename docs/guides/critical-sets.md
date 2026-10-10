# Big-steps critical sets

[Documentation](../README.md) / Guides

Use `optimization::CriticalSetWorkspace` to propose finite birth/death or
essential-birth moves on a frozen `SimplicialComplex` over F2. It retains sparse
primal R/V, lazily reduces the anti-transposed boundary for birth queries, and
computes only the required bounded rows of U = V inverse. This is a concrete
critical-set operation; ordinary persistence requests retain their own storage
policy. No foreign runtime or optimizer dependency is required.

The [paper framework and three algorithms](../research/big-steps.md) explain the
reduction and critical-set method. The [performance and integration report](../../benches/reports/critical-sets.md)
records measured costs and why the public API exposes one sparse workspace.

## Query a finite point

```rust
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::optimization::CriticalSetWorkspace;

let source = SimplicialComplex::new(vec![
    Simplex::new(vec![0], -2.)?,
    Simplex::new(vec![1], -1.)?,
    Simplex::new(vec![0, 1], 0.)?,
])?;
let mut workspace = CriticalSetWorkspace::new(&source)?;
let (birth, death) = workspace.finite_pairs().next().unwrap();
assert_eq!(source.simplex(birth).unwrap().vertices(), [1]);
assert_eq!(workspace.critical_set(death, 1.)?, [death]);
let targets = workspace.targets(&[(birth, -0.5), (death, -0.5)])?;
assert_eq!(targets, [(birth, -0.5), (death, -0.5)]);
# Ok::<(), cocycle::Error>(())
```

Pairs include zero-lifetime bars. Essentiality is relative to the complete
supplied complex; an omitted skeleton or truncated geometric source is not
certified. IDs belong to this source and must be looked up again after rebuilding
the filtration. The workspace borrows the source and never changes its values.

## Interpret proposals and update data

Each query returns same-dimensional simplices in source order, including ties at
the target boundary. A zero move yields an empty set. Finite death increases and
birth decreases solve bounded primal and dual U rows; the opposite directions
read V columns. Essential birth decreases use primal V and increases use dual V.
The [mathematical specification](../reference/mathematics.md#19-big-steps-critical-sets-over-f2)
states the formulas and their lazy-reduction precondition.

`targets` resolves overlaps by the largest absolute displacement from each
simplex's original value. Equal displacements keep the first input proposal.
This matches Oineus's `Max` at the
[pinned source](https://github.com/anigmetov/oineus/blob/e52814a1ffb5b8a81e71ff1e93b4c14194673f0f/include/oineus/top_optimizer.h#L1006).
It is a heuristic. It does not guarantee descent of an arbitrary diagram loss
or track a selected point through arbitrary simultaneous updates.

With targets frozen, the squared surrogate has simplex gradient `2*(f-target)`.
To optimize lower-star data, apply the chain rule to active maximum vertices and
reconstruct every simplex value from the updated vertex data. At a tie a caller
must choose a subgradient convention. For arbitrary explicit simplex updates,
include the required face/coface closure and resolve its conflicts before
constructing a valid new source. Returned sets alone do not perform this closure.

The [runnable example](../../examples/critical_sets.rs) collapses one short H0
interval toward the diagonal, accumulates a frozen-target lower-star gradient,
and rebuilds the filtration. Run `cargo run --locked --example critical_sets`.
Its small fixed step reduces the interval length from 1 to 0.25; this example
does not establish convergence or a general step-size rule.

## Storage and execution

Construction keeps sparse primal R/V for the whole supplied complex. The dual
state is built only when a birth query needs it. V transposes are cached only for
requested dimensions. U rows are solved afresh and not retained; increasing
pivots allow the solve to stop beyond the target value without forming a full
inverse. Sparse fill-in may still be quadratic and reduction work cubic.

`new_with`, `critical_set_with` and `targets_with` each start a fresh `Execution`
budget. A batch shares one budget across all queries and overlap handling.
Counted units include incidence, column visits, sparse terms and output items;
they are not bytes or stable performance counters. Controls do not impose a
memory quota. Interruption returns an error without partial output. Complete
lazy caches may survive; incomplete caches are discarded so a retry is valid.

Only the existing serial F2 lazy reducer supplies this state. Importing a
parallel or cleared decomposition would require proving or restoring the ELZ
support condition, including unit diagonals, before using critical-set formulas.
This API does not accept external decompositions, compute a full U, provide
automatic differentiation, or implement a general optimization loop.
