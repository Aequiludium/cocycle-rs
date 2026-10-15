# Standard real discrete Dirac operators

[Documentation](../index.md) / Guides

This guide describes the interface under review in
[PR #133](https://github.com/Aequiludium/cocycle-rs/pull/133), checked against
[source 27cd5eb](https://github.com/Aequiludium/cocycle-rs/tree/27cd5ebfb6218ddc41099e53c27daf7050f40d35/src/dirac).
It is absent from this documentation branch and the released 0.1.1 crate.
The example requires that candidate checkout. Rustdoc in that checkout owns the
complete signatures and errors; the [API reference](../reference/api.md#operator-interfaces-under-review)
records the interface entry points.

A Dirac operator couples adjacent degrees of a finite signed chain complex.
It accepts a [supplied complex or stage](filtered-complexes.md), uses the
Euclidean inner product and needs no persistence calculation. Filtration values
select stages and never become operator weights.

## Apply Dirac to one edge

```rust,ignore
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::dirac::{DiracOperator, IterativeOptions};
use cocycle::execution::Execution;

let source = SimplicialComplex::new(vec![
    Simplex::new(vec![10], -2.)?,
    Simplex::new(vec![20], -1.)?,
    Simplex::new(vec![10, 20], 3.)?,
])?;
let stage = source.stage(3.)?;
let operator = DiracOperator::new(&stage, 0..=1)?;
assert_eq!(operator.degree_offsets(), &[0, 2, 3]);
assert_eq!(operator.action(&[2., 5., 7.])?, [-7., 7., 3.]);
assert_eq!(operator.matrix(9)?, [0., 0., -1., 0., 0., 1., -1., 1., 0.]);

let execution = Execution::default().max_work(10_000);
let pairs = operator.boundary_spectrum(
    1, 1, 1e-8, &IterativeOptions::default(), &execution,
)?;
assert!((pairs.values()[0] + 2_f64.sqrt()).abs() < 1e-9);
assert!((pairs.values()[1] - 2_f64.sqrt()).abs() < 1e-9);
assert!(pairs.residuals().iter().all(|&value| value < 1e-9));
assert!(!pairs.is_index_coverage_certified());

let edges_only = DiracOperator::new(&stage, 1..=1)?;
assert!(edges_only.is_truncated());
assert_eq!(edges_only.action(&[7.])?, [0.]);
# Ok::<(), cocycle::Error>(())
```

Run `cargo run --locked --example dirac` in the candidate checkout for a filled
triangle with printed bases. The unequal vertex births above make the two
vertex coordinates explicit; in general read `basis()` rather than guessing
coordinate order from labels.

## Select degrees and interpret the operator

On $C_a\oplus\cdots\oplus C_b$, retain $B_j:C_j\to C_{j-1}$ only when
$a<j\le b$. Their degree-lowering block map $d$ gives

$$
D=d+d^\top,\qquad
D^2\big|_{C_j}=B_j^\top B_j+B_{j+1}B_{j+1}^\top.
$$

Maps crossing either selected endpoint are zero. Ordinary nonaugmented chains
have no degree-minus-one space. Thus `1..=1` removes both adjacent maps rather
than retaining the full complex's degree-one Hodge block. Comparison with
[Hodge](hodge-laplacian.md) must use the same degree truncation.

`new` and `new_with` accept `C: FilteredComplex<CellId = I>` with
`I: Copy + Eq + Hash`. They validate signed integer boundaries and $d^2=0$;
validity only modulo two is insufficient. `basis()` owns IDs grouped by degree
and source traversal order. `degree_offsets()` has one offset per degree plus
the final coordinate count; empty degrees have equal adjacent offsets.
`degrees()` returns the inclusive interval. `is_truncated()` reports excluded
cells of the supplied source, not completeness of an unseen larger complex.

The operator and spectral results outlive source storage, but their IDs remain
owner-local. There is no implicit source matching, persistence-diagram conversion
or rank certificate.

## Choose action, matrix export or signed modes

| Operation | Returned data and work |
| --- | --- |
| `new(source, degrees)` / `new_with(source, degrees, execution)` | Owned `DiracOperator<I>` with frozen signed incidence |
| `entries()` | Iterator of `(row, column, f64)` entries, including each boundary term and transpose |
| `action(x)` / `action_with(x, execution)` | Owned `Vec<f64>` in the entire graded basis; sparse action |
| `matrix(max_entries)` / `matrix_with(max_entries, execution)` | Owned row-major $N\times N$ matrix, with explicit output-entry cap |
| `boundary_spectrum(degree, count, zero_threshold, options, execution)` | `DiracBoundarySpectrum<I>` for one retained boundary map |

`boundary_spectrum` targets the `count` largest positive singular modes of
$B_j$ through matrix-free $B_j^\top B_j$. For $B_jv=\sigma u$, it reconstructs
$(u,\pm v)/\sqrt2$ with eigenvalues $\pm\sigma$ in the whole selected Dirac
basis. `zero_threshold` is finite, nonnegative and in singular-value units;
requested modes at or below it are rejected, not replaced by padding.

`singular_values()` is descending. `values()` is ascending:
$[-\sigma_{\max},\ldots,-\sigma_{\min},\sigma_{\min},\ldots,\sigma_{\max}]$.
For `count=r`, `vectors()` is $N\times2r$ row-major vector columns.
`basis()`, `degree_offsets()`, `degrees()` and `boundary_degree()` retain the
coordinate problem. `residuals()` checks the actual full truncated Dirac action;
`gram_residuals()` describes the normal-map solve separately.

This call does not select the global smallest-magnitude Dirac eigenvalues or
return the kernel. `is_index_coverage_certified()` is always false. Small
residuals do not rule out omitted invariant subspaces or multiplicities.
`iterations()`, `gram_applications()` and `peak_subspace_dimension()` expose
actual solver work, without treating Gram applications as all Dirac actions.

## Control work and understand numerical limits

`MAX_DIRAC_DEGREES` is 4096, checked before degree-offset allocation or source
access. Dense export checks $N^2\le\text{max_entries}$ before allocation;
sparse action has no such dense-output requirement. Each controlled constructor
or terminal starts a fresh `Execution` budget. Sorting, source callbacks and
allocation are cooperative checkpoints. Work units are not seconds or bytes.

The exported `IterativeOptions` shares the Hodge candidate's private solver:
block size 8, 800 updates, seed 73, absolute/relative residual tolerances 1e-9
and 1,048,576 entries per basis/image or projected buffer. It does not share the
persistent candidate's distinct Schur/Jacobi policy. Nonfinite arithmetic,
overflow, work/cancellation limits and failed convergence or paired-vector
checks return errors without an apparently successful partial result.

Normal maps square the nonzero condition number. Small singular values may be
unreliable even when action is sparse. The retained
[resource study](https://github.com/Aequiludium/cocycle-rs/blob/27cd5ebfb6218ddc41099e53c27daf7050f40d35/benches/reports/dirac-resources.md)
includes clustered graph spectra and an ill-conditioned signed integer chain
that fails explicitly. It establishes no general sparse scaling or smallest-mode
accuracy. Weighted/complex adjoints, persistent Dirac and signal filters remain
outside this interface.

## Check the candidate independently

The candidate's `tools/check_dirac.py --pairs` checks complete SVD/Dirac spectra,
requested cardinality, multiplicities, signed residuals and same-window Hodge
blocks. Supported pinned TopoNetX functions and GUDHI F3 topology checks are
separate evidence; finite-field topology does not certify real spectral indexes.
The [candidate mathematical specification](https://github.com/Aequiludium/cocycle-rs/blob/27cd5ebfb6218ddc41099e53c27daf7050f40d35/docs/reference/mathematics.md#22-standard-real-discrete-dirac)
owns the derivation and detailed contracts.
