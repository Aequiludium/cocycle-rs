# Descriptor contracts and external handoff information

[Documentation](../README.md) / Research / [#84](https://github.com/Aequiludium/cocycle-rs/issues/84) and [#85](https://github.com/Aequiludium/cocycle-rs/issues/85)

**Historical baseline and future handoff record.** The `a7ee9e9` source/status
observations below remain dated evidence. Current descriptors are described by
their guide/rustdoc and the [architecture](../development/architecture.md#six-responsibilities-and-data-routes).
The [version status](../README.md#version-and-capability-status) owns release scope;
the [route record](route-integration.md) links future learning/feedback handoffs.

Current-scope decision, 2026-10-11. Baseline main is
`a7ee9e911dcdeaea137b682bf383f7c2e1d090d5`. The current issues cover existing
descriptors and interface information only. #61 and #71 are closed;
#87 supplies the pre-migration baseline below. New features, serialization,
bindings, training and external learning runtime integration are outside this
contribution. Future #92 and #93 work does not block this decision.

## Existing descriptor contract (#84)

`betti_curve(&PersistenceDiagram, dimension, &[f64])` and
`finite_lifetime_summary(&PersistenceDiagram, dimension)` borrow an existing
owned result. They neither compute PH nor inspect reducer storage. A caller can
pass result.diagram() or data.diagram(); extracting these references does not
clone intervals or restore discarded context. See the
[result interpretation](result-contracts.md) and
[descriptor tutorial](../development/diagram-analysis.md).

| Input / output | Current behavior and independent regression |
| --- | --- |
| Computed-empty versus missing degree | Computed-empty gives zero counts, finite_count=0, total=0, max/entropy=None, exclusions=0. Missing degree returns DimensionNotComputed, including for an empty grid. `descriptors_distinguish_gaps_from_computed_empty_dimensions` protects gaps. |
| Grid | Finite and strictly increasing, including signed scales; duplicate +/-0, NaN, infinity and descending values error. Empty grid is empty after degree validation. `grid_validation_rejects_invalid_order_values_and_uncomputed_dimensions` protects the contract. |
| Complete / Through coverage | Complete may be queried past endpoints. Through(t) rejects every query >t, even for an empty selected degree; t itself is included. `betti_counts_births_deaths_and_the_censoring_cutoff_correctly` and the new cutoff-birth regression cover boundaries. |
| Endpoints and multiplicity | Count birth <= query; finite death excludes the death scale. Essential classes remain alive, censored classes remain alive through their cutoff, repeated intervals count separately. Public finite zero lifetimes are rejected by PersistenceInterval; birth=through is legal for censoring. |
| Finite summaries | Finite positive lifetimes only; dimension-specific essential/censored exclusion counts. Compensated total, optional maximum and Shannon entropy in nats; no normalization by ln(count). One finite bar has entropy zero; no finite bars gives None, including censored-only input. |
| Numerical failures | Finite endpoints can have overflowing lifetime/total; NumericalFailure is an error, not infinity or empty output. Extreme probabilities may round to zero and contribute their limiting zero entropy. Existing independent large/small-lifetime tests remain. |
| Units and provenance | Total/max lifetimes use the supplied scale unit; entropy/counts are dimensionless. Raw diagrams do not authenticate samples, units, metric hypotheses or fields. One descriptor query does not automatically compare unrelated sources. Preserve context separately when needed. |

The existing supplied-diagram example has repeated H1 [-1,1) bars and an
essential class born at zero. At [-2,-1,0,1,2] its counts are [0,2,3,1,1], finite
sum 4 and entropy ln(2). The truncated counterpart excludes its censored class
from finite summaries and counts it at the cutoff. These are existing examples,
not a new feature algorithm. The example now also asserts the entropy/empty
summary contracts explicitly, without changing output or invoking PH.

## Source disposition and consumers (#85)

[#71's fixed migration map](https://github.com/Aequiludium/cocycle-rs/blob/4599b12c83756a7feb76f25d8898e5bc2cfba402/docs/research/interface-inventory.md)
selects the following concrete retention decisions:

| File / symbol | Actual consumers / tests | Disposition |
| --- | --- | --- |
| `src/descriptors/mod.rs` exports | src/lib.rs, README, diagram-analysis tutorial | Retain the established public paths. |
| `src/descriptors/betti_curve.rs::betti_curve` | examples/diagram_analysis.rs, square.rs, rips_graph.rs; tests/descriptors.rs, euclidean.rs, filtered_complex.rs; README/Rips/interface doctests | Retain the logical degree view and existing event sorting. No stale caller requires migration. |
| `src/descriptors/lifetime_statistics.rs::{LifetimeSummary,finite_lifetime_summary}` | diagram_analysis.rs, square.rs; tests/descriptors.rs | Retain borrowed scans, compensated sum, exclusions and None semantics. No second summary owner. |
| `src/diagram/{dimensions,persistence_diagram,interval}.rs`, src/error.rs | Descriptor, contract and distance consumers | Retain availability/coverage/endpoint validation and structured errors; do not duplicate them in an analysis adapter. |
| `src/optimization/`, shared src/algebra/ | examples/critical_sets.rs, tests/critical_sets.rs, private reduction tests | Retain opt-in CriticalSetWorkspace with frozen complex/algebra computation. It is not #75's homology projector and belongs outside a training layer. |

Import/call inspection finds no geometry/reducer dependency, stale consumer or
proven obsolete private duplicate. All production src and Cargo files remain
byte-identical to baseline. The disposition is a supported no-op for production
migration, with targeted example/contract regressions; the issue does not require
moving already correctly owned code. Pending #104's source adapters and #105's
result readers are preserved, without adopting their trees.

For #94, finite summary retains O(m) time / O(1) auxiliary space; Betti retains
O(m log m + g) time / O(m + g) storage. Borrowing does not imply zero allocation
for Betti curves. No production sorting, copy, allocation, dispatch or source
lifetime changes occur; no new timing/resource comparison or speedup is claimed.

## Future external input information (#84 to #92)

This is a document-level handoff, not a format/API implementation. Retain a
sample ID and source revision alongside all future data; row position alone is
insufficient for a separate label table. Labels join by ID with unique keys and
complete matches, never by accidental order. Keep declared sample order and an
explicit feature shape; missing/uncomputed/unsupported/failed is not numeric zero.

| Candidate input | Required information before an external consumer can interpret it |
| --- | --- |
| Points / raw graph | Original point/vertex IDs, coordinate shape and units, edge list/direction, weight meaning, isolated vertices and source/sample identity. No PH or diagram is required. A graph is not automatically its clique complex. |
| Frozen topology / selected filtration scale | Source-local cell basis, original labels, degree, oriented boundary, selected finite scale/coverage, coefficient interpretation and feature rows in exactly that basis order. IDs from different owners cannot be interchanged. |
| Diagram feature table | Computed dimensions, endpoint policy, grid in shared units, fixed column order/length, essential/censored exclusion counts and empty optional fields. Learned bounds/scaling must be fit on training samples only. |
| Future spectrum / F2 operator output | Degree and ordered basis, field/inner product or Gram, units, selected range, residual/multiplicity or feasibility/optimality certificate, separate source and cost meaning. No existing descriptor certifies a spectral solver or projector. |

For a document-level #92 example, a frozen triangle graph may hand off vertex
features shaped 3-by-c0 and edge features shaped 3-by-c1, with signed incidence
rows/columns named by the same vertex/edge bases. It needs no filtration sweep,
PH, diagram or distance. A filled triangle additionally needs one face row and
its boundary. The current graph/frozen owner can expose the topology, but this
crate has no tensor conversion, Laplacian solver or external model execution.
Source-local simplex order need not equal lexicographic vertex feature order;
keep an explicit row map and orientations. Do not present this sketch as a
working TDL integration.

| Fixed source actually read | Conclusion / selected boundary | Limit |
| --- | --- | --- |
| [Persim PersistenceImager.fit/transform](https://github.com/scikit-tda/persim/blob/096a25fe38644d31c6aabb880a24a05046c03559/persim/images.py) | fit derives image bounds from diagrams; transform returns fixed-resolution arrays in sample order. Future learned bounds require train-only fitting. | Empty transform and per-diagram min/max in fit need an explicit application policy. No persistence image or Persim runtime is integrated. |
| [TopoModelX SCN2.forward](https://github.com/pyt-team/TopoModelX/blob/8cfdfbcb12eb8db6e0abd1816a1555e6f5cb604a/topomodelx/nn/simplicial/scn2.py) | Takes degree-specific features and Laplacians directly, demonstrating that a diagram/distance is not an input requirement. | Bases/channels/operators must align; no Torch/model execution is delivered. |
| [TopoNetX simplicial adjacency/coadjacency](https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/toponetx/classes/simplicial_complex.py) | index=True retains matrix/simplex correspondence; the read version rejects non-None weights. | Preserve row identities and signed incidence; do not reinterpret filtration values as spectral weights. No package integration is delivered. |

Encoding remains undecided. A simple file could expose identity and units but
materializes/serializes arrays; a binding can avoid text parsing but adds ABI,
packaging and lifetime contracts; an independent adapter has explicit copies
and its own dependency/version owner. Record those alternatives until an actual
consumer determines the minimal format. No new serialization or ML dependency
is selected. Training, visualizations, new feature algorithms and #93 feedback
remain separate future work.

## Pre-migration baseline (#87) and current validation

The baseline record uses main a7ee9e9, with unchanged production kernels and
existing hand fixtures. Full tests cover old free/builder Rips, supplied flag,
external/frozen cells, signed births, non-flag complexes, optional representatives,
computed-empty/missing degrees, multiplicity, endpoints, context, cancellation,
budgets and retry behavior. Independent private explicit/bitmask/field/duality
oracles remain distinct from production kernels.

The three existing consumption choices are independently checked on the unit
square and a square scaled by two: direct H1 Betti [0,1,0,0] and [0,0,1,0] at
[0,1,2,3]; contextual bottleneck distance sqrt(2)-1; and an application reading
that scalar to check it is below 0.5. This last step is ordinary scalar consumption,
not a new analysis/learning algorithm. Fields and declared scale conventions
match; coordinate units are chosen consistently by this fixture's caller.

Reproduce existing checks from a source checkout:

```sh
cargo test --locked --all-features
cargo test --locked --release --all-features
cargo test --locked --test descriptors --test contracts --test critical_sets
cargo run --locked --example square
cargo run --locked --example rips_graph
cargo run --locked --example diagram_analysis
cargo run --locked --example diagram_distances
python tools/check_algorithm.py diagram-analysis
```

Full CONTRIBUTING verification includes MSRV, all examples, strict rustdoc,
Markdown doctests and package checks. Local checks, hosted CI and actual merge
remain separate. Native performance studies need not be repeated for unchanged
kernels; raw historical evidence is never relabeled as a new measurement.


## Local acceptance for this scope

Rust 1.98.1 / Python 3.10.11 on Windows, 2026-10-11:

- Exact archived main baseline: 216 tests/doctests pass, four diagnostics ignored.
  The updated five-doctest tutorial passes against that baseline's own library,
  including all three existing consumer choices. This completes the reusable
  pre-migration record required by #87; it does not require future algorithms.
- Candidate debug/release and Rust 1.91 suites: 217 pass, four diagnostics ignored.
  The extra regression counts repeated censored classes born exactly at the
  cutoff without manufacturing finite lifetimes, and preserves coverage errors
  even in computed-empty dimensions. All-target MSRV check, Clippy, formatting,
  focused descriptor checks and strict rustdoc pass.
- All ten existing examples and eleven Markdown doctest files pass; the changed
  tutorial additionally runs against both base and candidate. Package verification
  and the packaged diagram-analysis example pass. No runtime/CI/manifest change.
- Independent GUDHI 3.12.0 reproduces the example's repeated H1 [-1,1), [-1,1),
  [0,infinity) intervals, five Betti counts and finite sum 4, from two triangles
  filled at 1 plus one unfilled loop born at 0. Full H0/H1 output is retained
  locally. Kernels/VR/distances are unchanged; no new native performance claim.
- Linux Python tools: 98 pass. Windows retains an unchanged POSIX-resource import
  error, a slash-specific diagnostic failure and five skips; this is not a passing
  Windows tool-suite result. Source/docs/index/diff checks pass.
- Primary local collaboration checks pass. Copied ignored research records lack
  their original target artifacts; that relocation check is recorded separately.

Commands and raw local logs remain in the contribution worktree's ignored
`target/analysis-evidence/`. Earlier exploratory graph-export/learning files and
results are local historical records and are absent from the final contribution
against main. Hosted CI belongs to the final head, with review/merge/release
recorded separately; older-head checks do not certify this revision.
