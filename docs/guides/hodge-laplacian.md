# Real single-scale Hodge Laplacians

[Documentation](../index.md) / Guides

`HodgeEigenpairs::is_index_coverage_certified()` always returns false. Residual
convergence verifies returned Ritz pairs, not their extremal indexes, unseen zero
modes or omitted multiplicities. Independent complete NumPy spectra check the
selected prefix on finite fixtures, including a 128-node path, a 64-cycle and
repeated disconnected edges. The missed-invariant-subspace regression remains.

The [API entry points](../reference/api.md#operator-interfaces-under-review) list the callable interfaces.

This **Unreleased** capability implements #73's initial standard-real chain
window. It reads a finite [supplied complex or stage](filtered-complexes.md),
keeps its signed incidence and degree order, and computes action, matrices and
small dense spectra without a persistence request. Filtration values never
become weights. The [mathematical specification](../reference/mathematics.md#22-real-single-scale-hodge-laplacians)
defines the operator and its harmonic space.

## Compute an ordinary edge spectrum

```rust
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::hodge::{HodgeLaplacian, HodgeSelection, IterativeOptions, SpectrumOptions};

let source = SimplicialComplex::new(vec![
    Simplex::new(vec![0], 0.)?,
    Simplex::new(vec![1], 0.)?,
    Simplex::new(vec![0, 1], 0.)?,
])?;
let execution = Execution::default().max_work(10_000);
let laplacian = HodgeLaplacian::new_with(&source, 0, HodgeSelection::Full, &execution)?;
assert_eq!(laplacian.matrix(&execution)?, [1., -1., -1., 1.]);
assert_eq!(laplacian.matrix_with_limit(4, &execution)?, [1., -1., -1., 1.]);
assert_eq!(laplacian.action(&[2., 5.], &execution)?, [-3., 3.]);
let spectrum = laplacian.spectrum(0..2, &SpectrumOptions::default(), &execution)?;
assert!(spectrum.values()[0].abs() < 1e-10);
assert!((spectrum.values()[1] - 2.).abs() < 1e-10);
assert_eq!(spectrum.numerical_nullity(1e-10)?, 1);
let selected = laplacian.smallest_eigenpairs(1, &IterativeOptions::default(), &execution)?;
assert!(selected.values()[0].abs() < 1e-9);
assert!(selected.residuals()[0] < 1e-9);
for column in spectrum.numerical_kernel_columns(1e-10)? {
    let vector: Vec<_> = (0..spectrum.basis().len())
        .map(|row| spectrum.vectors()[row * spectrum.values().len() + column])
        .collect();
    assert_eq!(vector.len(), spectrum.basis().len());
    assert!(spectrum.residuals()[column] < 1e-10);
}
# Ok::<(), cocycle::Error>(())
```

## Source, coordinates and ownership

`HodgeLaplacian::new_with` scans source cell IDs, then reads only boundaries in
degrees k and k+1. It owns the ordered three-degree snapshot. Basis slices use
source traversal restricted to each degree; IDs are local to the source owner.
Preparation rejects unknown/duplicate IDs, wrong-degree/nonpreceding faces,
nonzero vertex boundaries and nonzero integer `AD`. Duplicate terms are summed;
the chain check uses checked i128 arithmetic before any f64 conversion. Summed
incidences must fit the exactly representable f64 integer range. Filtration
ordering and unrelated windows remain the source's responsibility.

`new` delegates to `new_with` with unlimited default controls. Controls supplied
only after `new` cannot limit completed preparation. The example uses `new_with`
so preparation is explicitly bounded. The same immutable `Execution` value is
reused for later terminals, each of which resets its counter; there is no shared
prepare-plus-spectrum limit. Spectrum itself shares one budget across matrix
conversion, eigensolving and error measurement. Source construction above is a
separate operation outside the Hodge preparation budget.

`action` returns an owned coordinate vector; `matrix` returns an n by n
row-major buffer. Interpret both using the prepared operator's `basis()`.
`matrix` checks n² against `DEFAULT_MAX_DENSE_ENTRIES` (1,048,576 entries,
8 MiB for f64) before dense allocation, even with unlimited `Execution`.
Use `matrix_with_limit(max_dense_entries, execution)` to set a smaller cap or
explicitly permit a larger trusted computation. Zero permits only n=0.
`HodgeSpectrum` owns its domain basis handles, actual matrix, degree, summands
and numeric policy alongside the values and vector columns. These define its
coordinate eigenproblem after source release. Keep your source descriptions to
interpret IDs geometrically; matching numeric IDs do not identify different
owners. Vectors supplied to action must already be in the declared basis:
length validation cannot detect a caller's undeclared permutation.

Missing degrees in this finite source define explicit zero-size maps. At k=0,
the lower boundary has zero rows; no augmented vertex sum is added. Top-degree
and empty-degree computations are supported. A supplied skeleton is its own
finite problem; it does not certify the full operator of a larger source with
omitted cofaces. Use an explicit stage to choose a scale. This is an ordinary
stage operator, not a persistent Laplacian between two stages.

`Full`, `Lower` and `Upper` select `A^T A + D D^T`, `A^T A` and `D D^T`.
Only the full operator's harmonic nullity has the ordinary real Betti meaning.
The API fixes real f64 coefficients and the standard inner product. General
SPD/Hermitian Grams, complex coefficients and normalized operators cannot be
expressed through this API; no input is implicitly converted to those variants.

## Matrix-free partial eigenpairs

`smallest_eigenpairs(count, options, execution)` targets a small number of low
eigenpairs directly through sparse action. It uses identity-preconditioned block
LOBPCG with the local search space spanned by current vectors X, residuals R and
previous directions P. Two-pass orthogonalization drops relative directions
below 1e-13; only the small symmetric Rayleigh-Ritz matrix is diagonalized.
The private real symmetric core is shared with the Dirac boundary-spectrum
implementation; it introduces no runtime dependency or operator normalization.

`IterativeOptions` defaults to an 8-vector block, 800 updates, seed 73 and
absolute/relative residual tolerances 1e-9. Require count <= min(n,block_size).
With b=min(n,block_size), the search dimension is <=min(n,3*b). Each basis/image
block and projected square buffer must fit `max_basis_entries` (default
1,048,576 entries), checked before starting. Several buffers coexist, so this
is not a total memory quota. Choose b much smaller than n: if b approaches n,
the subspace itself becomes quadratic. No dense Hodge matrix is constructed.

`HodgeEigenpairs` owns basis handles, degree, summands and policy, plus ascending
Ritz values and n-by-count row-major vectors. Final residuals use fresh sparse
actions and must satisfy `absolute_tolerance + relative_tolerance*|lambda|`.
The result reports orthogonality, block updates, operator applications and peak
Ritz dimension. Exceeding the iteration/workspace/work limit returns an error,
not an apparently converged subset. Count zero makes no operator applications.
There is no fallback to full dense materialization.

Residual convergence does not certify extremal index coverage. A retained
synthetic test uses A=I-u*u^T with u orthogonal to the fixed starting block:
returned eigenvalues near one have tiny residuals while a zero mode was missed.
This is an internal numerical counterexample, not an extension of the public
integer-boundary input. Change the seed/block or use an independent reference
when coverage matters. This result deliberately offers no nullity/Betti reader.
Identity preconditioning can converge slowly on small-gap paths; no universal
speedup or convergence guarantee is asserted.

The three-term block iteration follows [Knyazev (2001)](https://doi.org/10.1137/S1064827500366124).
[SciPy 1.15.3's LOBPCG](https://github.com/scipy/scipy/blob/v1.15.3/scipy/sparse/linalg/_eigen/lobpcg/lobpcg.py)
is an independent numerical reference, not copied production code. Weighted and
normalized higher-order Laplacians in [Horak and Jost, arXiv:1105.2712v3](https://arxiv.org/abs/1105.2712v3)
change the mathematical model and remain explicit future capabilities.

## Spectrum and numerical evidence

The dense reference path uses scaled cyclic Jacobi rotations. It preserves the
crate's absence of runtime dependencies. `SpectrumOptions` sets a checked
dimension limit (default 128), `max_dense_entries` (default 1,048,576),
maximum sweeps (64) and finite nonnegative
absolute/relative tolerances (both 1e-12). Convergence uses the off-diagonal
Frobenius norm in scaled coordinates, equivalent to an unscaled tolerance of
`absolute_tolerance + relative_tolerance * ||L||_F`. Failed convergence is an
error. This bounded path has not been evaluated as a large sparse solver.

The requested index range selects from the full algebraically sorted dense
decomposition. Multiplicity is preserved. `vectors()` contains an n by r
row-major matrix with value j paired to vector column j; signs and repeated
eigenspace bases are arbitrary. Small negative eigenvalues are retained.
Both the dimension and n² entry limits are checked before dense allocation,
including for an empty or partial eigenpair selection. Raising max_dimension
alone does not bypass the shared dense entry cap.

Results report actual absolute residuals from sparse action, relative residuals
using `||L||_F`, returned-column Frobenius orthogonality error, sweeps and policy.
A zero relative denominator gives `None` and retains the absolute residual.
An empty selection on a nonempty domain remains partial.
`numerical_nullity(threshold)` requires all eigenpairs and uses an explicit,
independent absolute zero threshold. It provides floating-point evidence, not
an exact real Betti certificate. Prime-field Betti numbers cannot generally be
substituted for real homology. Repeated eigenvectors should be checked by
subspace projectors rather than elementwise equality.

`numerical_kernel_columns(threshold)` borrows the selected vector-column indexes
without allocation or another solve and requires all eigenpairs. Index j pairs
`values()[j]`, `residuals()[j]` and coordinate `vectors()[row*r+j]`, with
r=`values().len()`. For `Full` this reads the numerical harmonic space. The
same read on `Lower`/`Upper` refers only to that summand's kernel.

Convergence, small residual and agreement with exact nullity are three separate
claims. A backend regression uses the exact-rank-three block matrix
`diag([[5e-13,-5e-13],[-5e-13,5e-13]], 1e12, 2e12)`. Under the default global
stopping policy it can pass without a sweep, with residual 5e-13, yet classify
zero kernel columns at threshold 1e-14 despite exact nullity one. A separately
declared stricter policy resolves this particular small block. Neither result
provides a general nullity certificate. This is a synthetic solver test, not
an expansion of the public integer-boundary API to arbitrary Gram matrices.

## What is reproduced from the reference implementations

The fixed [ICML 2023 source](https://github.com/tomogwen/persistentlaplaciandatascience/blob/b36b080f739653d84284f3f7ca9ccea5fad286ae/src/perslap.py)
is executed as a reference, including its original `spectrum_hodgelap`.
The correspondence is precise rather than an assertion that every internal
algorithm matches:

| Reference operation | Actual source behavior | Native implementation and verification |
| --- | --- | --- |
| `perslap.boundary` | Its positive-degree matrices are the negative of increasing-vertex boundaries on the closed simplicial fixtures. | Retain existing source orientation; verify that degree-wise overall sign change explicitly. It cancels in both Hodge products. |
| `perslap.hodge_laplacian` | Constructs `B_down.T @ B_down + B_up @ B_up.T`; degree zero lower summand is zero. | Same mathematical operator; sparse action first, dense conversion on request. Compare full matrices in supported degrees; hand-derived Rust tests check lower/upper selections separately. |
| `perslap.spectrum_hodgelap` | Calls NumPy `eig` for every degree below its stored maximum; returns matrices and eigenvalue arrays, not eigenvectors. | Execute the original function and compare sorted real eigenvalues with multiplicity. The Rust eigenvector/residual owner adds information the reference does not return. |
| TopoNetX signed incidence/Hodge | Standard omission signs and sparse products; degree-zero incidence alone is augmented. | Compare signed Hodge matrices under ordinary H0; do not use its augmented incidence-zero map in the Rust operator. |
| Lim §6.2 | Computes a harmonic component as least-squares residual after projecting onto `im A^T` and `im D`. | Independently compute that residual with NumPy `lstsq` and compare it to the returned numerical harmonic-column projector applied to x. |

The Rust scaled cyclic Jacobi routine is an independently selected bounded
solver. It is **not** a port of the paper's NumPy eigensolver or a numerical
method prescribed by Lim. It is retained as a small native reference path with
explicit convergence policy and the adverse regression above; mature numerical
libraries remain the independent oracle. Its tests cover near-repeated roots,
a positive eigenvalue of order 1e-14, global scales 1e-180/1e180, and a 128-node
path with analytic spectrum `2-2*cos(pi*j/128)`. They do not establish accuracy
for every ill-conditioned matrix up to dimension 128.

## Cost, cancellation and independent reproduction

Preparation retains sparse incidence and source ID metadata; exact chain checks
can cost more than reading incidence. Action needs O(n[k-1]+n[k]) scratch and
work proportional to the adjacent incidence plus coordinates. Matrix conversion
uses O(n²) output and n sparse actions. The dense solver retains O(n²) storage
and takes O(max_sweeps * n³) work, including full decomposition for a partial
selection. Checked shapes, byte capacities and fallible reservations precede
dense allocation; not every standard-library allocation is recoverable.
The entry cap bounds each dense square buffer individually. Spectra retain
several such buffers plus selected vectors and scratch; the cap is not a total
allocation or process-memory quota. Sparse action has no dense entry limit.

Each terminal uses a fresh `Execution` budget. Constructor traversal, incidence
and integer chain products; action coordinates and incidence; dense entries,
rotations and residual/orthogonality products are charged. Allocation, sorting
and source iterator internals are cooperative checkpoints, not interruptible
instructions. Work units are not elapsed time or a memory quota. No partial
result is returned on cancellation, overflow, nonfinite arithmetic or exhaustion.

From a source checkout, run:

```sh
cargo test --locked --test hodge
cargo run --locked --example hodge_laplacian > target/hodge-rust.jsonl
python tools/check_hodge.py --rust target/hodge-rust.jsonl --output target/hodge-oracle.json --iterative
```

The optional Python oracle requires NumPy; external probes additionally require
SciPy and GUDHI. Its input is the actual Rust example output. It independently
builds signed simplicial boundaries, uses NumPy `eigh`, computes exact rational
ranks and checks actions, residuals, orthogonality and eigenspace projectors.
Fixtures include the edge, unfilled/filled triangle, relabelled path, isolated
vertices, tetrahedron surface/solid, disconnected and empty sources, and twelve
fixed six-vertex complexes. Added fixtures include the repository's six-vertex
real projective plane and the 128-node path. RP2 has real beta1=0 while its F2
beta1=1; the real Hodge producer is checked independently of the F2 reducer.
Rust tests also cover generic cells, orientation,
stages, malformed integer chains, partial spectra and execution errors.

For optional external probes, download the exact files and provide
`--toponetx PATH --perslap PATH --gudhi` to the same command:

- [TopoNetX ce3d414, simplicial source](https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/toponetx/classes/simplicial_complex.py).
  SHA-256 is checked. The unchanged incidence/Hodge function bodies run with a
  skeleton adapter and SciPy; this does not validate the entire package.
  [Its license](https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/LICENSE.txt)
  is MIT. Degree-zero sources without edges and absent degrees are explicit exclusions.
- [ICML 2023 perslap b36b080](https://github.com/tomogwen/persistentlaplaciandatascience/blob/b36b080f739653d84284f3f7ca9ccea5fad286ae/src/perslap.py).
  Git blob identity is checked. Unchanged selected functions run with a minimal
  row-membership adapter; the original spectrum function runs NumPy `eig`.
  Highest/absent/empty degrees are explicit exclusions;
  no synthetic padding is added to force support. The file header declares MIT;
  a standalone root license was not established. No source is copied into Rust.
- GUDHI computes ordinary Betti numbers over prime field 101, with the top
  degree enabled. This is separate finite-field evidence; independent rational
  ranks supply the real Betti expectation for every window.

The oracle also checks the documented weighted-path Gram matrix and analytic
eigenvalues as a future-extension handoff. That check does not implement weighted
Rust Hodge, Dirac (#74), persistent Laplacians (#124), or downstream #88/#109.
