# Real single-scale Hodge Laplacians

[Documentation](../index.md) / Guides

This guide describes the interface under review in
[PR #130](https://github.com/Aequiludium/cocycle-rs/pull/130), checked against
[source 138db59](https://github.com/Aequiludium/cocycle-rs/tree/138db59377314272aa16a0accb5b2fc4568375a3/src/hodge).
It is absent from this documentation branch and the released 0.1.1 crate.
The example requires that candidate checkout. Rustdoc in that checkout owns the
complete signatures and errors; the [API reference](../reference/api.md#operator-interfaces-under-review)
records the interface entry points.

A Hodge operator acts on real chains of one degree in a finite
[supplied complex or stage](filtered-complexes.md). It uses signed incidence
and the standard Euclidean inner product. No persistence calculation is needed,
and filtration values do not become operator weights.

## Compute an ordinary edge spectrum

```rust,ignore
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::hodge::{HodgeLaplacian, HodgeSelection, IterativeOptions, SpectrumOptions};

let source = SimplicialComplex::new(vec![
    Simplex::new(vec![0], 0.)?,
    Simplex::new(vec![1], 0.)?,
    Simplex::new(vec![0, 1], 0.)?,
])?;
let execution = Execution::default().max_work(10_000);
let operator = HodgeLaplacian::new_with(&source, 0, HodgeSelection::Full, &execution)?;
assert_eq!(operator.matrix(&execution)?, [1., -1., -1., 1.]);
assert_eq!(operator.action(&[2., 5.], &execution)?, [-3., 3.]);

let complete = operator.spectrum(0..2, &SpectrumOptions::default(), &execution)?;
assert!(complete.values()[0].abs() < 1e-10);
assert!((complete.values()[1] - 2.).abs() < 1e-10);
assert_eq!(complete.numerical_nullity(1e-10)?, 1);

let selected = operator.smallest_eigenpairs(1, &IterativeOptions::default(), &execution)?;
assert!(selected.residuals()[0] < 1e-9);
assert!(!selected.is_index_coverage_certified());
# Ok::<(), cocycle::Error>(())
```

Run `cargo run --locked --example hodge_laplacian` in the candidate checkout for
additional full, lower and upper examples.

## Select the operator and its coordinates

Write $B_q:C_q\to C_{q-1}$ for signed incidence in the selected source. The
ordinary, nonaugmented convention has $B_0=0$:

$$
L_q=B_q^\top B_q+B_{q+1}B_{q+1}^\top.
$$

`HodgeSelection::Full`, `Lower` and `Upper` select respectively this sum,
$B_q^\top B_q$ and $B_{q+1}B_{q+1}^\top$. Only the full operator's exact nullity
equals the real Betti number. The harmonic space is
$\ker B_q\cap\ker B_{q+1}^\top$; it is not an F2 projection or a shortest-cycle
selection. See the [F2 operator guide](homology-operators.md) for that separate model.

`new` and `new_with` accept `C: FilteredComplex<CellId = Id>` with
`Id: Copy + Eq + Hash`. They validate ordering, signed incidence and the integer
chain identity in the adjacent-degree window. An identity holding only modulo
two is insufficient. Use `source.stage(scale)?` to select cells of value <= scale.

`basis()`, `previous_basis()` and `next_basis()` preserve source traversal within
their degrees. Matrix rows and action coordinates use `basis()`, with row-major
matrix storage. The operator and results own their data; IDs remain local to the
source owner. Preserve that owner description for geometric interpretation.
Absent adjacent spaces are zero spaces of the supplied topology, not evidence
that an incompletely constructed larger complex has no cofaces.

## Choose a matrix or a spectral result

| Operation | Returned data and work |
| --- | --- |
| `action(x, execution)` | Owned real coordinates from sparse incidence action; no dense matrix |
| `matrix(execution)` / `matrix_with_limit(max_dense_entries, execution)` | Owned row-major `Vec<f64>` of all $n^2$ entries |
| `spectrum(range, options, execution)` | `HodgeSpectrum<Id>`; complete dense decomposition before selecting the half-open algebraic index range |
| `smallest_eigenpairs(count, options, execution)` | `HodgeEigenpairs<Id>`; approximate ascending low Ritz pairs through sparse action |

Both spectral results expose `basis()`, `degree()`, `selection()`, `values()`,
`vectors()`, `residuals()` and `orthogonality_error()`. For $r$ returned values,
vectors are $n\times r$ row-major columns: entry `(row, column)` is at
`row*r + column`. Eigenvector signs are arbitrary; repeated eigenspaces may rotate.

The dense result additionally retains `matrix()`, `index_range()` and
`is_complete()`. `numerical_nullity(threshold)` and
`numerical_kernel_columns(threshold)` require a complete spectrum and a finite
nonnegative absolute threshold. They select $|\lambda|\le\text{threshold}$;
they do not certify exact rank or Betti numbers.

The iterative result retains no full matrix. `iterations()`,
`operator_applications()` and `peak_subspace_dimension()` report actual work.
`is_index_coverage_certified()` is always false: a small residual can coexist
with a missed zero mode, repeated root or smaller invariant subspace. A caller
requiring the true smallest prefix must establish coverage independently.

## Control work and understand numerical limits

`SpectrumOptions` defaults to dimension 128, 64 cyclic Jacobi sweeps,
absolute/relative tolerance 1e-12, and 1,048,576 entries per dense square buffer.
`matrix()` uses that same entry cap without the spectral dimension limit.
Dense storage is quadratic; a range request still pays for full decomposition.

`IterativeOptions` defaults to block size 8, 800 updates, seed 73,
absolute/relative tolerance 1e-9 and 1,048,576 entries per basis/image or projected
buffer. Identity-preconditioned block LOBPCG uses a Ritz space of at most
`min(n, 3*block_size)` columns. Larger blocks can approach quadratic storage.
Fresh sparse action checks residuals before a result is returned.

Construction and each terminal start separate `Execution` budgets. Passing the
same controls to several calls does not combine their counters. Cancellation,
work exhaustion, invalid shapes, nonfinite arithmetic and nonconvergence return
errors without a successful partial result. Buffer caps and cooperative work
limits are not total RAM quotas or elapsed-time limits.

The default 1024-vertex path exhausts its iteration budget in the retained
[resource study](https://github.com/Aequiludium/cocycle-rs/blob/138db59377314272aa16a0accb5b2fc4568375a3/benches/reports/hodge-resources.md).
Independent complete spectra and preconditioning controls diagnose that case;
they do not establish universal convergence or a new Rust preconditioner.
General Gram matrices, complex/weighted/normalized models and persistent pairs
remain outside this interface. Use the [persistent guide](persistent-laplacian.md)
for the separate two-stage operator.

## Check the candidate independently

The candidate's `tools/check_hodge.py --iterative` compares actual example output
with complete NumPy spectra, exact rational ranks, signed matrices and harmonic
projectors. Pinned TopoNetX/perslap functions and supported GUDHI field-101
checks provide additional evidence. GUDHI's finite-field Betti numbers are not
general real-rank certificates. The
[candidate mathematical specification](https://github.com/Aequiludium/cocycle-rs/blob/138db59377314272aa16a0accb5b2fc4568375a3/docs/reference/mathematics.md#22-real-single-scale-hodge-laplacians)
owns the detailed numerical policy and derivation.
