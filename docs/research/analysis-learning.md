# Analysis and optional learning decisions

[Documentation](../README.md) / Research / [#84](https://github.com/Aequiludium/cocycle-rs/issues/84) and [#85](https://github.com/Aequiludium/cocycle-rs/issues/85)

Decision record, 2026-10-11. Baseline main is
`a7ee9e911dcdeaea137b682bf383f7c2e1d090d5`. #61, #63 and #71 are closed;
[#71's migration map](https://github.com/Aequiludium/cocycle-rs/blob/4599b12c83756a7feb76f25d8898e5bc2cfba402/docs/research/interface-inventory.md)
selects retention of the existing descriptor and optional optimization owners.
This delivery does not consume new result, distance, stage or spectral interfaces.
#79-L and #82 therefore introduce no additional prerequisite for these examples.

## Read sources and select the first workloads

| Fixed source actually read | Conclusion and decision | Remaining limit |
| --- | --- | --- |
| [Persim images.py](https://github.com/scikit-tda/persim/blob/096a25fe38644d31c6aabb880a24a05046c03559/persim/images.py), `PersistenceImager.fit/transform` | fit derives image bounds from supplied diagrams; transform preserves sample order and returns a fixed-resolution image. Use train-only fitting for learned bounds. Retain existing Betti/lifetime features for the first statistical example instead of adding images. | Empty transform returns zeros, while fit uses per-diagram min/max; finite endpoints and empty training input need an explicit application policy. Persim is source-inspected, not executed here. |
| [TopoModelX SCN2.forward](https://github.com/pyt-team/TopoModelX/blob/8cfdfbcb12eb8db6e0abd1816a1555e6f5cb604a/topomodelx/nn/simplicial/scn2.py) at `8cfdfbcb12eb8db6e0abd1816a1555e6f5cb604a` | Requires degree-specific feature tensors and Laplacians directly; no diagram or distance is required. Select a graph-feature classifier first to exercise that optional boundary with fewer dependencies. | No SCN2, Torch training or TDL integration is claimed. Matrix bases and channel order must align before any future TDL adapter. |
| [TopoNetX simplicial matrices](https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/toponetx/classes/simplicial_complex.py), adjacency/coadjacency | `index=True` retains matrix-to-simplex correspondence; the inspected version rejects non-None weights. Preserve basis identity and signed incidence explicitly in future matrix exchange. | Graph degree features below are unweighted; filtration values are not silently spectral weights. TopoNetX is source-inspected, not executed. |
| [scikit-learn Pipeline.fit](https://github.com/scikit-learn/scikit-learn/blob/25dee604bae18205b01548348388baf7a1cdfe0e/sklearn/pipeline.py) and [StandardScaler](https://github.com/scikit-learn/scikit-learn/blob/25dee604bae18205b01548348388baf7a1cdfe0e/sklearn/preprocessing/_data.py), source tag and executable version 1.7.2 | A tabular sample-by-feature matrix suffices. Fit a scaler and one-level tree only on the declared training rows; predict held-out rows in file order. | A tiny hand-labeled fixture establishes exchange and fitting behavior, not predictive benefit on independent application data. |

Statistics: reuse the supplied repeated H1 bars `[-1,1)`, `[-1,1)` and an
essential class born at zero. At grid `[-2,-1,0,1,2]` the expected Betti counts
are `[0,2,3,1,1]`; total finite lifetime is 4 and entropy is ln(2) nats.
The truncated counterpart counts its censored class at the cutoff and excludes
it from finite-lifetime statistics. No topology construction or PH is needed.

Learning: supply ten actual WeightedGraph objects: paths and rings on
3, 4, 6 (train) and 5, 7 (test) vertices. For n vertices a path has n-1 edges,
mean degree 2(n-1)/n and leaf fraction 2/n; a ring has n edges, mean degree 2
and leaf fraction 0. Labels come from the fixture construction, not PH or a
learned model. The classifier uses leaf fraction only. Vertex/edge counts and
mean degree remain exported audit columns, and labels/IDs/split are never features.
Both examples have hand-derived expectations and run without diagram distances.

## Existing source disposition and consumers (#85-L)

| File / symbol | Live consumers and regressions | Disposition |
| --- | --- | --- |
| `src/descriptors/mod.rs` exports | `src/lib.rs`, README and diagram-analysis walkthrough | Retain public paths; no second descriptor API. |
| `src/descriptors/betti_curve.rs::betti_curve` | `examples/{diagram_analysis,square,rips_graph}.rs`, `tests/{descriptors,euclidean,filtered_complex}.rs`, Rips/interface guides | Retain borrowed degree view and event sorting. Computed-empty is zeros; missing degree errors even for an empty grid. Finite deaths are excluded; births and censored cutoff are included. |
| `src/descriptors/lifetime_statistics.rs::{LifetimeSummary,finite_lifetime_summary}` | `diagram_analysis`, `square`, `tests/descriptors.rs` | Retain compensated finite sum, nats, dimension-specific essential/censored exclusions, None for empty max/entropy. |
| `src/diagram/{dimensions,persistence_diagram,interval}.rs` and `src/error.rs` | Descriptor, contract and distance tests | Retain availability, coverage and errors. Zero finite lifetimes are rejected by PersistenceInterval; they never reach descriptors. A class born at its censored cutoff is legal. |
| `src/optimization/` and shared `src/algebra/` | `examples/critical_sets.rs`, `tests/critical_sets.rs`, private independent reduction tests | Retain opt-in CriticalSetWorkspace with complex/algebra computation. It is neither a learning-layer object nor #75's homology projector; differentiation and application feedback remain #93. |

The descriptor import/call search finds no reducer/geometry dependency, stale
caller or proven obsolete private duplicate. All these files remain byte-identical
to baseline. The existing result.diagram() consumer is already the intended
borrowed read; moving it would add no capability. #85-L is therefore a justified
no-op, supported by source inspection and baseline/candidate consumer regression.
#87 G0/G1 is exercised by full existing route tests and square/rips_graph examples;
this does not close #87's wider future-route acceptance.

For #94, finite summary retains O(m) time / O(1) auxiliary space; Betti retains
O(m log m + g) time / O(m + g) space. No production read, copy, sorting,
allocation or source lifetime changes, so no new performance comparison is
needed or claimed. Graph feature extraction visits n adjacency lengths in O(n)
time / O(1) auxiliary space, borrows the graph and returns four scalars. Graph
construction retains its existing O(n + m log m) cost. CSV and Python conversion
materialize ten rows; scaling/classification costs belong to the external process.

## Exchange contract and selected scope (#85-N)

The existing [Rust example](../../examples/diagram_analysis.rs) retains its default
statistical mode. `--export PATH` selects a separate graph-only branch and emits
UTF-8 CSV with the exact header:

```text
sample_id,split,label,source,vertex_count,edge_count,mean_degree,leaf_fraction
```

Each row is one graph in declared fixture order, with unique `path-N`/`cycle-N`
sample ID, train/test split and `source=unweighted-graph-v1`. All numeric features
are dimensionless. Rust f64 shortest-roundtrip decimal output is parsed as Python
float64; there is no quantization or learned normalization in Rust. The four
features have fixed order, independent of vertex labels. For empty graphs the
feature function returns None; isolated nonempty graphs have zero degree/leaf
features. These inputs are tested but excluded from the two-class teaching set.

The [external example](../../examples/analysis_learning.py) validates the full
fixture, unique IDs, identity/label/count consistency, finite numbers, split and
feature formulas before importing the optional ML runtime. It rejects empty
files, missing columns/dimensions, missing/mismatched samples, duplicate IDs,
unknown sources and nonfinite features; no missing value is replaced with zero.
Rows may be reordered together: IDs/labels/splits travel with their features,
and prediction IDs retain the selected row order. There is no separate positional
label join. Float formula checks use absolute 1e-14 tolerance for values <= 2,
well above one f64 rounding step and far below the fixture's distinct values.

Pipeline.fit sees six training rows only; four test rows never set scaler
statistics, thresholds or hyperparameters. Fixed max_depth=1/random_state=0
are chosen in advance. The fitted scaler's sample count and mean are checked
against training-only expectations. Predictions are compared to the hand labels.
Do not tune this example after looking at held-out outcomes.

Run from a source checkout; Cargo creates target before the export step:

```sh
cargo run --locked --example diagram_analysis
cargo run --locked --example diagram_analysis -- --export target/graph-features.csv
python -m venv target/learning-venv
# POSIX: target/learning-venv/bin/python; Windows: target/learning-venv/Scripts/python
# Substitute that interpreter for python in the following two commands.
python -m pip install scikit-learn==1.7.2 numpy==2.2.6 scipy==1.15.3 joblib==1.5.3 threadpoolctl==3.6.0
python examples/analysis_learning.py target/graph-features.csv
```

CSV is chosen over a binding to keep Cargo dependency-free and make identities
and conversion visible. A binding could avoid text parsing but adds ABI, packaging
and array/lifetime contracts. A separate in-memory adapter may suit repeated
large exchanges after measured need. CSV is not a universal matrix/operator format.
For spectra, select degree, basis, Gram, physical units, eigenvalue range, residual
and multiplicity before exporting statistics. F2 projection features need class,
basis, cost units and feasibility/optimality certificates. These production
features await their actual computation owners; no endpoint-derived proxy is used.

Continue the two selected examples. Defer landscapes/images, spectra/operator
features, matrix TDL, visualization and differentiable feedback until concrete
inputs, algorithms and independent acceptance exist. No complete ML platform,
real-world accuracy, speedup, new kernel or release claim is made.


## Local acceptance at the prepared source

On 2026-10-11, Rust 1.98.1 / Python 3.10.11 on Windows:

- Before edits, descriptor/contract/critical-set integration tests: 38 passed.
  After edits, production src, Cargo manifest and lock are unchanged; the full
  debug/release and Rust 1.91 suites each pass 216 tests/doctests, four ignored
  diagnostics. All-target MSRV check, Clippy, formatting and strict rustdoc pass.
- All ten default examples, analysis/construction/critical-set example tests in
  debug/release, Rust 1.91 analysis tests and all eleven Markdown doctest files
  pass. The focused diagram-analysis command passes. Package verification and
  the packaged analysis example and its tests pass; optional Python files are
  included under the existing examples package pattern.
- The actual CSV exporter and scikit-learn 1.7.2 run pass with NumPy 2.2.6,
  SciPy 1.15.3, joblib 1.5.3 and threadpoolctl 3.6.0. Six train/four test rows;
  training scaler mean 0.25; test IDs path-5, cycle-5, path-7, cycle-7 match
  the four construction labels. This checks the stated teaching set only.
- Exported fixture SHA-256:
  `06985147a4ccfa1228e04c46f76278cdd122e6a4f4692dc37d6e142e6140ef3f`.
  A separate GUDHI 3.12.0 SimplexTree check reproduces all ten degree/count
  features without clique expansion. Two disjoint triangles with edges at -1
  and faces at 1, plus an unfilled triangle with edges at 0 (vertices at -2),
  independently give H1 [-1,1), [-1,1), [0,infinity), the five stated Betti
  counts and total finite lifetime 4. The full H0/H1 oracle output is retained.
- Four new standard-library exchange regressions pass; the whole Linux Python
  tool suite passes 102 tests. Windows runs 102 with one existing POSIX resource
  import error, one existing slash-specific diagnostic failure and five skips;
  it is not a passing Windows tool-suite result.
- Source/docs/diff/index checks pass. Primary local collaboration checks pass;
  copied ignored research records lack their original target artifacts, so that
  relocation check is recorded separately. No local record or artifact is staged.

Raw outputs and the independent oracle script stay under ignored
`target/analysis-evidence/` in the contribution worktree. They are local evidence,
not a durable hosted artifact archive. CI separately tests the optional external
example with pinned packages and logs the fixture hash. Hosted CI, review, merge
and release must be recorded separately at the final commit; none is inferred
from these local checks. No new native speed/performance comparison is claimed,
and unchanged VR/distance algorithms are covered by existing CI reference jobs.
