# Changelog

## Unreleased

- Point the package documentation URL to the generated API reference on docs.rs
  and link it from the README and documentation index alongside the usage guides.

## 0.1.0

Initial release of Cocycle, a pure Rust topological data analysis library. The
crate is named `cocycle`; the source repository is `Aequiludium/cocycle-rs`.
Requires **Rust 1.91 or newer**, uses edition 2024, and is **MIT licensed**.
There are no runtime dependencies, foreign TDA backends or unsafe Rust.

### Analysis workflow

- Borrow Euclidean point clouds or lower-triangular, upper-triangular and square
  dissimilarity matrices; construct exact threshold graphs from points, matrices
  or distance callbacks. Supplied weighted graphs define flag filtrations.
- Use `RipsBuilder` to configure a filtration, then `.persistence().compute()` to
  compute ordinary persistent homology over a validated prime field. F2 is the
  default; specialized H0/H1 computation and dimension-generic implicit
  cohomology are selected internally.
- Use `build_complex` when explicit topology is needed. Frozen simplicial
  complexes support simplex lookup, oriented boundaries and cofacet queries.
  Construction dimension, analysis dimension and scale coverage remain distinct.
- Analyze supplied simplicial or filtered-cell complexes through
  `PersistenceBuilder::from_complex`. The four-method `FilteredComplex` contract
  supports signed filtration values and unequal vertex births without imposing
  Rips edge rules on other algorithms.
- Construct deterministic sparse Rips approximations with explicit metric
  hypotheses, insertion-radius provenance, modified edges and higher-simplex
  blockers. Both direct persistence and explicit construction are available.
- Request persistent cycle bases and query-scale dual cocycles on simplicial
  sources. Representatives retain original vertex IDs and interval identities.

### Diagrams and analysis

- Owned diagrams retain interval multiplicity and distinguish finite deaths,
  essential classes and right-censored intervals. Rips scales are edge lengths;
  surviving a cutoff does not imply an essential class.
- Explicit computed-dimension sets distinguish computed-empty dimensions from
  uncomputed dimensions, including gaps. Queries and consumers validate the
  selected dimension rather than inferring membership from the maximum.
- Derive finite lifetime summaries, persistence entropy in nats and Betti curves.
- Compare complete diagrams with exact bottleneck distance using L-infinity
  costs, W1 using L-infinity costs, or W2 using Euclidean costs. Matching preserves
  essential multiplicity. Context-aware operations validate fields and declared
  scale conventions; callers still establish common units and normalization.
- Use cooperative execution limits and cancellation on Builder operations and
  controlled distance entry points. Distance `_with` variants share one budget
  across preparation, matching and accumulation. Existing uncontrolled calls
  remain available; limits are not hard process-memory or wall-clock caps.

### Documentation and verification

- Runnable examples cover point clouds, threshold graphs, supplied complexes,
  higher-dimensional Rips, prime fields, representatives, sparse approximation,
  diagram descriptors and matching distances.
- Contributor walkthroughs separate complex construction, persistence reduction
  and diagram analysis so that new algorithms can reuse appropriate data and
  validation without adopting unrelated algorithms' internal representations.
- Mathematical tests include analytic fixtures, independent boundary/rank and
  matching oracles, interval multiplicity, coverage and resource-failure checks.
  Native comparison suites use GUDHI and upstream Ripser C++, not Python wrappers.
- Performance reports bind results to measured commits and protocols. They are
  evidence for those revisions and workloads, not a speed claim for every 0.1.0
  workflow. Raw measurements, reference-source downloads and logs are excluded
  from the repository and crate payload.

### Limits and next steps

This release supports a Rust workflow from point clouds or supplied complexes to
persistence diagrams, summaries and diagram distances. It does **not** claim full
GUDHI feature parity. Alpha and cubical constructors, mutable simplex-tree
editing, reduced homology, zigzag and multiparameter persistence, non-prime
coefficient rings, language bindings and dataframe integration are not included.
The lower-star construction example is a contributor tutorial, not a production
constructor API.

Exact point-cloud graph construction scans point pairs without a spatial index.
Explicit simplex expansion can grow exponentially; implicit computation still
faces combinatorial enumeration and reduction fill-in. Requested representatives
need additional computation and storage. Sparse approximation requires its stated
metric hypotheses. Diagram distances currently require complete coverage and
reject censored inputs. Numerical results use finite `f64` filtration values,
not exact real arithmetic. See the [mathematical specification](docs/reference/mathematics.md)
for precise contracts and the [roadmap](docs/design/roadmap.md) for future work.

### Moving from development Git revisions

There is no earlier crates.io release to migrate from. Users of older Git
snapshots should note the pre-publication API changes:

- Prefer `RipsBuilder` and `.persistence().compute()` for new code. Legacy Rips
  functions remain supported in 0.1.0; their documented execution scope and
  dimension limits still apply.
- `PersistenceDiagram::intervals()` yields values in canonical logical order.
  Replace slice-based access with `len`, `is_empty`, `interval(index)` and
  `dimension(k)?.iter()`, or explicitly collect an owned buffer. Storage layout
  and input-buffer reuse are not public contracts.
- `max_dimension()` is only the greatest computed dimension. Use
  `computed_dimensions()` or `dimension(k)` for membership. Default builders
  continue to compute every dimension through the requested maximum.
- Common diagram/context data lives in `PersistenceData`. Context-aware distance
  functions accept compatible wrappers through `AsRef<PersistenceData>`;
  explicit function-item types may need updating. Representative indices resolve
  through `diagram.interval(rep.interval_index())`.

Future breaking API changes will use a new minor version while the crate remains
below 1.0. Version-specific release notes will describe the migration.
