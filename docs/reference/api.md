# API reference

[Documentation](../index.md) / Reference

Rustdoc owns individual signatures, errors and API contracts. Use the
[0.1.1 reference](https://docs.rs/cocycle/0.1.1/cocycle/) when your dependency is
`cocycle = "0.1.1"`; [docs.rs](https://docs.rs/cocycle) follows the published
crate and may differ from this site's main guides.

For the current checkout, run `cargo doc --locked --no-deps --open` from the
repository root. The generated entry is `target/doc/cocycle/index.html`.
Source-only capabilities are identified in the
[version status](../README.md#version-and-capability-status).

| API area | Start with |
| --- | --- |
| Borrowed points and dissimilarities; Rips and flag construction | [Rips construction guide](../guides/rips-construction.md) |
| Supplied complexes, signed boundaries, stages and inclusions | [Filtered complexes](../guides/filtered-complexes.md) |
| Persistence requests, fields and representative bases | [Rips guide](../guides/rips.md), [representatives](../guides/rips-representatives.md) |
| Diagrams, endpoints, coverage and result ownership | [Interface contracts](../guides/interface-contracts.md) |
| Diagram summaries and Betti curves | [Diagram analysis](../development/diagram-analysis.md) |
| Bottleneck, W1 and W2 matching | [Diagram distances](../guides/diagram-distances.md) |
| F2 critical sets and partial U | [Critical sets](../guides/critical-sets.md) |
| Feasible F2 projections and finite transport | [Homology operators](../guides/homology-operators.md) |

The [mathematical specification](mathematics.md) owns shared definitions;
[architecture](../development/architecture.md) locates their implementation.
Design and research sketches do not establish callable APIs.

## Operator interfaces under review

The following entries describe separate implementation candidates, absent from
this documentation branch and released 0.1.1 rustdoc. Their guides identify
the checked source revision; build rustdoc in that candidate checkout for exact
generic bounds, lifetimes, error variants and all accessors. The namespace is
specific to each operator, with no universal spectral result container.

| Candidate namespace and input | Construction / computation | Result and interpretation |
| --- | --- | --- |
| `cocycle::hodge`; `&C` where `C: FilteredComplex<CellId = Id>` | `HodgeLaplacian::new(source, degree, selection)` or `new_with(..., execution)` | `Result<HodgeLaplacian<Id>>`; `HodgeSelection::{Full, Lower, Upper}` chooses summands |
| Prepared `HodgeLaplacian<Id>` | `spectrum(range, &SpectrumOptions, execution)` | `Result<HodgeSpectrum<Id>>`; full dense decomposition before index selection |
| Prepared `HodgeLaplacian<Id>` | `smallest_eigenpairs(count, &IterativeOptions, execution)` | `Result<HodgeEigenpairs<Id>>`; approximate ascending Ritz pairs, index coverage uncertified |
| `cocycle::dirac`; `&C` where `C: FilteredComplex<CellId = I>` | `DiracOperator::new(source, degrees: RangeInclusive<usize>)` or `new_with(..., execution)` | `Result<DiracOperator<I>>` with graded basis and offsets |
| Prepared `DiracOperator<I>` | `boundary_spectrum(degree, count, zero_threshold, &IterativeOptions, execution)` | `Result<DiracBoundarySpectrum<I>>`; descending selected singular modes and ascending signed Dirac pairs |
| `cocycle::persistent_laplacian`; two `&SimplicialStage` values | `PersistentLaplacian::compute(k, l, degree, options, execution)` / `compute_up(...)` | `Result<PersistentLaplacian>`; owned pair matrix and Full/Up tag |
| One `&SimplicialStage` | `PersistentLaplacian::compute_down(k, degree, options, execution)` | Same owned type with Down tag and source scales `(s,s)` |
| Ordered stage slice and one target | `fixed_endpoint_laplacians(starts, target, degree, options, execution)` / `fixed_endpoint_up_laplacians(...)` | `Result<Vec<PersistentLaplacian>>` in input order |
| Prepared `PersistentLaplacian` | `spectrum(execution)` | `Result<LaplacianSpectrum>`; public component/eigenvalues/residual/zero-threshold fields, no eigenvectors |

All execution arguments above are `&Execution<'_>`. Persistent constructors take
`SpectralOptions` by value; Hodge/Dirac spectral calls borrow their options.
The Hodge and Dirac candidates re-export the same `IterativeOptions` policy.
The persistent candidate's `SpectralOptions` has distinct Schur/Jacobi semantics.

Use the [Hodge](../guides/hodge-laplacian.md), [Dirac](../guides/dirac.md) and
[persistent](../guides/persistent-laplacian.md) guides for coordinate order,
matrix/vector layout, result ownership, thresholds, defaults and failure behavior.
In particular, the two iterative result types return false from
`is_index_coverage_certified()`. Residual-small pairs cannot automatically serve
as a guaranteed smallest prefix for a downstream consumer.

The [spectral-feature reproduction](../research/persistent-spectral-features.md)
has a Python function/schema and CLI, with no new `cocycle::descriptors` or
classifier API. Its array contract must not be presented as a callable Rust
interface or a completed persistent-kernel integration.
