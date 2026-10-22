# Real persistent Laplacians

[Documentation](../index.md) / Guides

This guide describes the interface under review in
[PR #131](https://github.com/Aequiludium/cocycle-rs/pull/131), checked against
[source 867d9ac](https://github.com/Aequiludium/cocycle-rs/tree/867d9acfddde3cd68f9ad24a8a17a36574220acc/src/persistent_laplacian).
It is absent from this documentation branch and the released 0.1.1 crate.
The example requires that candidate checkout. Rustdoc in that checkout owns the
complete signatures and errors; the [API reference](../reference/api.md#operator-interfaces-under-review)
records the interface entry points.

A persistent Laplacian belongs to an inclusion $K_s\subseteq L_t$, rather than
one ordinary stage. It uses real coefficients, standard inner products and
ordinary nonaugmented chains. The candidate computes a bounded dense matrix,
action and complete numerical eigenvalues while preserving the pair's bases
and source metadata.

## Compute a two-stage operator

```rust,ignore
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::persistent_laplacian::{
    LaplacianComponent, PersistentLaplacian, SpectralOptions, fixed_endpoint_laplacians,
};

let source = SimplicialComplex::new(vec![
    Simplex::new(vec![0], 0.)?, Simplex::new(vec![2], 0.)?,
    Simplex::new(vec![1], 1.)?, Simplex::new(vec![0, 1], 2.)?,
    Simplex::new(vec![1, 2], 2.)?,
])?;
let k = source.stage(0.)?;
let l = source.stage(2.)?;
let execution = Execution::default();
let options = SpectralOptions::default();
let operator = PersistentLaplacian::compute(&k, &l, 0, options, &execution)?;
assert_eq!(operator.matrix(), &[0.5, -0.5, -0.5, 0.5]);
assert_eq!(operator.action(&[1., -1.], &execution)?, [1., -1.]);
let spectrum = operator.spectrum(&execution)?;
assert!(spectrum.eigenvalues[0].abs() < 1e-10);
assert!((spectrum.eigenvalues[1] - 1.).abs() < 1e-10);

let down = PersistentLaplacian::compute_down(&k, 0, options, &execution)?;
assert_eq!(down.component(), LaplacianComponent::Down);
assert_eq!(down.scales(), (0., 0.));
assert_eq!(down.matrix(), &[0.; 4]);
let starts = [source.stage(-1.)?, k, l];
let reverse = fixed_endpoint_laplacians(&starts, &l, 0, options, &execution)?;
assert_eq!(reverse[1].matrix(), operator.matrix());
# Ok::<(), cocycle::Error>(())
```

Run `cargo run --locked --example persistent_laplacian` in the candidate
checkout for additional pair and component output.

## Establish the inclusion and operator

Inputs are borrowed `SimplicialStage` values from the same original owner,
with source scale <= target scale and checked inclusion. Equal cloned complexes
or equal integer IDs do not establish common ownership. Certified constructions
must contain the required skeleton: full/upper need degree $q+1$ unless expansion
exhausted the source; lower-only needs degree $q$. Bare supplied complexes define
their entire mathematical source. Missing certified cofaces are errors.

Reorder the target upper Gram by the retained $q$-basis of $K$ and the eliminated
coordinates of $L$. With blocks $A,E,C$, the standard real operator is

$$
L_q^{s,t}=\underbrace{A-EC^+E^\top}_{U_q^{s,t}}
          +\underbrace{B_q(K)^\top B_q(K)}_{D_q^s}.
$$

The direct implementation solves $CX=E^\top$ in minimum norm before forming
$A-EX$; it does not materialize $C^+$. When $K=L$, the result is the ordinary
[full Hodge operator](hodge-laplacian.md). Exact full nullity has the persistent
Betti interpretation. Numerical threshold counts and component nullities do not
supply that certificate.

## Choose components and read results

All constructors below take a degree, `SpectralOptions` by value and an
`&Execution`, and return `Result<PersistentLaplacian>` unless a vector is shown.

| Operation | Input and returned component |
| --- | --- |
| `PersistentLaplacian::compute(k, l, degree, options, execution)` | Two stages; `LaplacianComponent::Full` |
| `compute_up(k, l, degree, options, execution)` | Two stages; `Up`, without the source lower term |
| `compute_down(k, degree, options, execution)` | Source only; `Down`, with recorded scales `(s,s)` |
| `fixed_endpoint_laplacians(starts, target, degree, options, execution)` | Ordered source stages and one target; `Vec<PersistentLaplacian>` in input order |
| `fixed_endpoint_up_laplacians(starts, target, degree, options, execution)` | Same fixed-target route, retaining upper components only |

The result owns `basis()`, `target_basis()`, `degree()`, `scales()`,
`source_context()`, `source_coverage()` and `construction_dimension()`.
`matrix()` borrows the owned row-major matrix; `action(x, execution)` returns
owned coordinates in `basis()` order. Dropping source storage does not invalidate
the operator. Metadata describes that source/pair, not cross-source equivalence.

`spectrum(execution)` returns `LaplacianSpectrum` with public fields
`component`, `eigenvalues`, `residual` and `zero_threshold`. It returns all
eigenvalues in increasing order with multiplicity and without clipping small
negative roundoff. Eigenvectors and partial index selection are not exposed.
Retain the originating operator's metadata when storing this result: its component
tag alone does not authenticate a pair or its coordinates.

For standard real chains, $UD=DU=0$. The nonzero full spectrum is the multiset
union of the component nonzero spectra; full zero multiplicity is
$z(U)+z(D)-\dim C_q(K)$. Concatenating both complete spectra doubles the domain
incorrectly. This interface does not infer exact ranks or recombine thresholded
spectra automatically.

## Control elimination and understand numerical limits

`SpectralOptions` defaults to `relative_tolerance=1e-12`,
`residual_tolerance=1e-9`, `max_sweeps=64` and `max_dense_entries=65_536`.
The last bounds each dense buffer and the cumulative fixed-endpoint output
entries. A target Gram therefore has at most 256 coordinates by default;
coexisting buffers are additional memory. Storage is quadratic, and each cyclic
Jacobi sweep has cubic work. This is a bounded dense reference implementation.

Fixed-endpoint construction reuses the target upper term, eliminates coordinates
in reverse filtration order and adds each source's own lower term. It preserves
ties, repeated scales and empty stages. A near-zero pivot may be skipped only
when retained coupling is exactly zero. Otherwise the route uses a checked direct
Schur fallback under the same cumulative `Execution` budget, tolerances and
storage policy. `used_direct_fallback()` records that choice.

Direct solving checks discarded positive-direction correction bounds as well as
solve and symmetry residuals. A small solve residual alone cannot justify dropping
a large Schur correction. Unsafe cases return `NumericalFailure`; increasing a
tolerance is not a correctness repair. Realizability of the general PSD
counterexamples by the current unweighted simplicial inputs has not been established.

`solver_residual()` describes normalized Schur/range/correction checks, whereas
`LaplacianSpectrum::residual` measures the selected matrix's eigensystem.
`reference_scale()` retains the larger assembled/input-Gram scale for PSD and
near-zero guidance. These are numerical diagnostics, not forward-error or exact
rank certificates. The persistent Jacobi solver retains its own policy separately
from the Hodge/Dirac block solver.

Construction and fixed-endpoint work use their documented cumulative budgets;
later action/spectrum calls start fresh budgets. Cancellation, work exhaustion,
invalid owners/scales, missing skeletons, storage overflow and numerical failures
return errors without successful partial output. Buffer limits are not a process
RAM quota. Weighted models, general chain maps, dynamic endpoint updates and
large sparse spectral solving remain outside this interface.

## Check the candidate independently

The candidate's comparison tool checks 1,114 finite pairs using exact Fraction
restricted-domain matrices/ranks and independent NumPy Schur/SVD routes, with
pinned author-function and supported GUDHI topology checks. Its
[mathematical specification](https://github.com/Aequiludium/cocycle-rs/blob/867d9acfddde3cd68f9ad24a8a17a36574220acc/docs/reference/mathematics.md#22-real-persistent-laplacians)
gives the derivation. The [resource study](https://github.com/Aequiludium/cocycle-rs/blob/867d9acfddde3cd68f9ad24a8a17a36574220acc/benches/reports/persistent-laplacian-resources.md)
retains repeated construction costs and default storage rejections. These checks
do not establish a scalable persistent spectral kernel or production integration
with [spectral-feature learning](../research/persistent-spectral-features.md).
