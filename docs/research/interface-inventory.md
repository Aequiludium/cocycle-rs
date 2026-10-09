# Interface inventory and ownership costs

[Documentation](../README.md) / Research / [#59](https://github.com/Aequiludium/cocycle-rs/issues/59)

Technical report, 2026-10-09. This is the complete baseline inventory behind
[input contracts](input-contracts.md), [result contracts](result-contracts.md)
and the [#62/#63 compatibility decision](../design/interface-compatibility.md).
It distinguishes public production operations, pending changes and mathematical
capabilities that still need implementation. The six frameworks classify
responsibilities; they do not require six crates, modules or public traits.

## Revision and inspection method

The original inventory used the complete source archive of main
[`9e6715f4c2e0118ab738486bd747ec3e361397de`](https://github.com/Aequiludium/cocycle-rs/commit/9e6715f4c2e0118ab738486bd747ec3e361397de).
Public entry points were checked against implementations, rustdoc, design pages,
examples and tests. All eight open PRs were inspected through their actual
base/head/state/draft, changed files, fixed base-to-head diffs and hosted checks.
The separate #32 report PR and both dependency stacks were included.

#12/#13/#15 had already reached main through #20/#25/#21; their merge commits
were verified as ancestors of this baseline. #14/#30/#31 still owned their
existing pending work. Research for #18/#19 in the then-unmerged #32 did not
establish production delivery. The local #49 experiment was not used as main
or as a substitute for the newer remote #49 head.

**Historical snapshot:** sections below describe that original baseline. During
publication, main had advanced to
[`cedf5965b774855614b761a1ad311d9dc7992c49`](https://github.com/Aequiludium/cocycle-rs/commit/cedf5965b774855614b761a1ad311d9dc7992c49).
The publication delta is recorded after the historical PR table. In particular,
the old H2 dispatch observation must not be read as the current main behavior.
Original measurements and check counts are retained with their original heads;
merging a report does not remeasure a later implementation.

## Six-framework capability map at the original baseline

| Framework | Production entry and behavior | Examples / tests | Missing capability and responsible follow-up |
| --- | --- | --- | --- |
| Data and topology modeling | [Geometry][geometry] provides borrowed `PointCloudView`, `DissimilarityView` and `DissimilarityMatrixView` over f64 buffers, with square/lower/upper layouts. `MetricPolicy` and `validate_metric` distinguish unchecked, caller-declared and checked metric evidence. [WeightedGraph][graph] owns edges/adjacency; [SimplicialComplex][simplicial] owns face-closed topology and signed incidence. [FilteredComplex][filtered-contract] has four statically generic methods: cells/dimension/value/boundary. | `square`, `rips_graph`, `complex_construction`; `matrix`, `euclidean`, `graph`, `filtered_complex` tests | No Alpha construction, production cubical container or general sheaf/path/hypergraph/Mayer model. The lower-star example composes existing APIs; it is not a new library constructor. Owner-local IDs do not establish sample/basis identity across results. #60 and #65/#66 own bases, weights and new source definitions. |
| Filtration and stage access | [RipsBuilder / ApproximateRipsBuilder][filtration] can `prepare()` an owned source or `build_complex(q)` a `SimplicialFiltration`. `FlagFiltration` describes clique filtration on a supplied graph. Supplied simplicial/cell inputs support signed values and unequal vertex births. Construction cap, analysis cap, coverage and simplex dimension have separate meanings. | `rips_graph`, `flag_persistence`, `sparse_rips`, `complex_construction`; `rips_construction`, `rips_expansion`, `sparse_rips`, `filtered_complex` tests | Expanded topology permits borrowed boundary/cofacet reads; unexpanded source simplex access remains private. Filtering an existing owner by scale is not a public identified stage/map API or a persistent operator. #60 and #68/#69 establish stage and chain-map requirements. |
| Topological computation | Sealed [PersistenceExt / PersistenceBuilder][persistence-source], including `from_complex`, compute ordinary PH over a validated prime-u32 field, default F2. Exact flag/Rips H0 uses union-find; F2 H0/H1 uses specialized implicit cohomology. At this baseline H2+ or odd-prime higher requests use generic cohomology; non-zero-born supplied inputs use boundary reduction. Simplicial sources can request persistent cycles and scale-specific dual cocycles. | `square`, `rips_sphere`, `rips_representatives`, `flag_persistence`; `prime_fields`, `rips_expansion`, `filtered_complex`, `execution` tests and independent private reducer | [Baseline dispatch][dispatch] contains `dimension > 1 || field.characteristic() != 2`. #50's continuation was then pending; see publication delta. #51 fixed-tuple H2 is experimental/private-cfg; #56 had failed performance gates. No Laplacian/eigensolver, Dirac, F2 homology operator or generalized-topology engine. #71 reuses PH; #72–#75 investigate new algorithms; #76 owns adoption. |
| Results and representations | [PersistenceDiagram / ComputedDimensions][diagram] preserve multiset, Finite/Essential/RightCensored and computed-empty versus absent dimensions. [PersistenceData][data] owns diagram/context; [PersistenceResult][result] combines optional representatives. Reads and `AsRef` borrow; `into_data` / `into_parts` move. `DiagramDimension` is a logical borrowed view yielding Copy interval values. | `diagram_analysis`, `rips_representatives`; `contracts`, `descriptors`, `diagram_distances` tests | Preserve the already-delivered #12/#13 semantics. Default builders compute contiguous `0..=q`; manually constructed diagrams can have gaps. Witness ordinals belong to their owning result. No full external result/context constructor or spectral/projector result exists. #61 and #78/#79 own those concrete decisions. |
| Distances and similarity | [Diagram distances][distances]: bottleneck/L-infinity, W1/L-infinity, W2/Euclidean; raw, `_results` and `_with(Execution)` calls. Require complete coverage and a computed selected dimension; preserve essential multiplicity. Unequal essential counts give +infinity; numeric failure gives Error. Context adapters also check field and declared edge-length convention. | `diagram_distances`; `diagram_distances`, `execution` tests and tiny exhaustive matching oracle | Every `Coverage::Through` is rejected, even without a censored bar. Raw calls do not check physical units/source; contextual calls cannot authenticate them either. Sparse-diagram comparisons do not certify original-data errors. #26/#49 own pending preparation; #81/#82 define other concrete similarities. |
| Analysis and optional learning | [betti_curve / finite_lifetime_summary][descriptors] read diagrams without rerunning PH. Summary exposes finite count, total/max persistence and natural-log entropy, with excluded essential/censored counts. | `diagram_analysis`; `descriptors` tests and diagram-analysis walkthrough | Two descriptor families, not a complete statistics/learning engine. No landscapes/images, trainer, visualization engine or ML runtime. Betti grids must be finite and strictly increasing, within coverage. Empty computed degree is valid; absent degree errors. #84/#85 and optional direct data-to-learning #92 own new consumers. |

Raw diagrams can already enter distance/analysis without modeling, filtration or
PH. A graph may enter a concrete single-scale calculation without a PH request.
That optional route does not implement a new spectral/operator engine. Likewise,
navigation between frameworks is not evidence that all routes exist.

## Borrowing, materialization and scratch boundaries

| Boundary | Observed ownership and cost | Decision / source |
| --- | --- | --- |
| Coordinates/matrix to view/builder | Views validate and borrow input slices. Builders hold views/options; `.persistence()` creates a request without computing it. | Retain [geometry][geometry], [builder][rips-builder] and [source][persistence-source]. Borrowing does not imply an allocation-free pipeline. |
| Points to exact PH | Without a finite cap, the direct point path creates an owned condensed distance buffer. With a finite cap it evaluates all pairs and retains threshold-graph edges inside the effective analysis range. | [source::exact][persistence-source]: separate O(n²) pair visits from retained-edge storage; preserve original construction range in context. |
| Input/callback to `prepare()` | Allocates owned edges, sorted adjacency and metadata, without copying the supplied matrix buffer. Callback terminal consumes the callback builder and calls each unordered pair once; later prepared analyses do not replay it. `WeightedGraph::new` takes the edge Vec; `FlagFiltration::new` moves the graph. | Explicit conversion cost in [exact][rips-exact], [builder][rips-builder] and [graph][graph]. |
| Source to explicit complex | Allocates simplices, vertex tuples, lookup and boundaries/cofacets; output may grow exponentially. Frozen queries subsequently borrow storage. | [storage][simplicial], [filtration result][simplicial-filtration]. Do not require expansion for diagram-only consumers. |
| Approximate source and context | Stores modified edges, permutation, insertion radii, retained original labels and metric evidence. Result assembly may clone owned approximation metadata. | [metadata][approximation], [builder][approx-builder], [assembly][approx-compute]. Epsilon >= 1 / unchecked metrics have no multiplicative bound; positive-radius subsampling changes the target and cannot inherit the same original-input bound silently. |
| Borrowed complex to reduction | Does not mutate the owner. Generic cells read all metadata, validate selected q+1 skeleton and create ID maps, columns, coefficients and reduction state. Frozen simplicial input reuses construction guarantees but still creates selected owned columns; zero-born diagram-only paths can use implicit cofaces. | [generic reader][filtered-reader], [frozen reader][simplicial-reader]. A field-specific boundary-square check is not an integer/all-fields certificate. |
| Source to representatives | Opt-in skeleton/transform storage, cycle/dual solves and owned vertex-labelled terms. Generic `from_complex` rejects simplex-labelled representatives. | [architecture][architecture], [representatives][representatives], [reader][filtered-reader]. Witnesses are not a free diagram field or generic cellular representation. |
| Result to borrowed data / moved parts | Diagram/context/representative access and `AsRef` borrow. `into_parts` moves both owners; `into_data` drops representatives; `into_diagram` drops context too. Logical interval iteration copies small values without allocating a second interval array or exposing rows. Explicit `clone()` copies owned buffers. | [data][data], [result][result], [diagram][diagram]. Preserve ordinal association on moves and record losses on conversion. |
| Baseline diagram to scalar distance | Each call's `points()` allocates finite `Vec<[f64;2]>` and essential `Vec<f64>`; kernels prepare additional storage. Borrowing a result adapter does not remove that cost. | [baseline facade][distances]. Existing #14/#26 owns this boundary; do not recreate it as Phase-4 work. |
| #26 / #49 preparation | #26 feeds logical degree views into private algorithm-specific preparation and streams essential births. #49 explicitly retains prepared operand arrays, borrowing an immutable source and fixing metric/dimension/context. Query budgets are fresh; solver/graph/normalization scratch remains pair-local. | [#26 source][pr26-distance], [#49 source][pr49-distance], both pending at publication. Arrays still allocate/copy; old #49 measurements do not certify the newer stacked head. |
| Diagram to descriptors | Lifetime summary scans with O(1) auxiliary storage. Betti curves allocate/sort events and an output Vec: O(m log m + g) time, O(m + g) storage. | [summary][lifetime-summary], [curve][betti-curve]. Borrowed reads do not imply zero allocations for every descriptor. |
| Execution / workspace | Reusable Execution is configuration. Each terminal gets a fresh private WorkBudget; a persistence builder covers preparation/reduction/requested representatives, a controlled distance covers validation/preparation/solve. Legacy free-function scope remains persistence-only. | [Execution][execution]. Counts are not time, bytes, RSS or threads. Accessor/callback/sort/allocation internals are not forcibly interruptible. #52 delivered no public Workspace; #32 established no hard-memory quota or production routing. |

## All eight pending PRs at the original snapshot

These exact heads were OPEN and unmerged at inventory time. Check counts are
GitHub check runs, including duplicated push/PR jobs; they are not counts of
independent mathematical validations. Mergeability labels establish neither
correctness nor performance admission.

| PR | Fixed head / actual base | Contribution and retained limit | Hosted snapshot |
| --- | --- | --- | --- |
| [#26](https://github.com/Aequiludium/cocycle-rs/pull/26) | `26866ac9223a7bd5cc7c35286e1090dc5759e509` / main `9e6715f4c2e0118ab738486bd747ec3e361397de` | #14 logical preparation preserves public metrics, coverage, multiplicity and controls; work counts can change. Conflicting measurements remain; no universal speedup. | 18 SUCCESS |
| [#49](https://github.com/Aequiludium/cocycle-rs/pull/49) | `5bb77b4b624580fe3cede4c8d538d06518de1a5d` / #26 `26866ac9223a7bd5cc7c35286e1090dc5759e509` | #31 opt-in PreparedDiagram preserves source lifetime, metric/dimension/raw-versus-context rules and failure recovery. Earlier `b60dc873` timings do not admit this head relative to #26; distinct-operand controls include regressions. | 18 SUCCESS |
| [#32](https://github.com/Aequiludium/cocycle-rs/pull/32) | `a0a034027f9fa5ccb89123bcf199cf8f7f0dfc06` / main `9e6715f4c2e0118ab738486bd747ec3e361397de` | Five Direction-B reports; changes confined to benches/reports, with no source/manifest/tool/workflow modification. Ancestry repair protects #26. Research no-go/conditional reuse is not production routing or Workspace delivery. | 16 SUCCESS |
| [#50](https://github.com/Aequiludium/cocycle-rs/pull/50) | `1bfd739dcf056ecc8f06254e0bbfc878df550bba` / main `9e6715f4c2e0118ab738486bd747ec3e361397de` | T0/T1/T2 continuation owns complete H1 death-triangle clearing, including omitted virtual pairs; generic computation starts at H2. Unsupported specialization/ID overflow retains fallback. H1-only requests do not collect clearing keys. | 16 SUCCESS |
| [#51](https://github.com/Aequiludium/cocycle-rs/pull/51) | `c5441935250145d0a37874689a0dbe3f4bb13d0e` / #50 `1bfd739dcf056ecc8f06254e0bbfc878df550bba` | Fixed-tuple H2 is tests/private `--cfg cocycle_h2_bench`; ordinary builds retain #50 generic continuation. T3 duplicates32 H1 median ratio 1.053131 exceeds 1.05, failing production admission. Favorable H2 rows do not establish a default speedup. | 16 SUCCESS |
| [#52](https://github.com/Aequiludium/cocycle-rs/pull/52) | `231ffbed0fea0dc4467bd9412a9ddf1a2a10c5f9` / #51 `c5441935250145d0a37874689a0dbe3f4bb13d0e` | M1 lifetime/scratch audit, test tracing/retry and complete source/harness reuse guards. H1 state already ends before H2 allocation. Generated ablations not adopted; no allocator-live-byte proof, public Workspace or completed M2. | 16 SUCCESS |
| [#53](https://github.com/Aequiludium/cocycle-rs/pull/53) | `c3fee147ca50673d2221d961397fa6c281056da2` / #52 `231ffbed0fea0dc4467bd9412a9ddf1a2a10c5f9` | Serial research closeout, M2 no-go, independent H0–H3 sphere regression and six native fixtures. No default reducer change relative to its parent. Finite H3 evidence is not arbitrary-dimensional proof or a parallel engine. | 16 SUCCESS; 8 CANCELLED retained separately |
| [#56](https://github.com/Aequiludium/cocycle-rs/pull/56) | `2dcd5e3d8a306f3cb4c4f23e486cdfc5add78db1` / #53 `c3fee147ca50673d2221d961397fa6c281056da2` | Draft adaptive F2/H1, bounded bitsets and ordered cursors retain H1-to-H2 handoff. Recorded ratios include 8.4–9.0× nonmetric-matrix and 3.2–3.4× truncated-bipartite regressions relative to #53. Correctness CI does not remove performance gates. | 16 SUCCESS |

The historical dependency chains were `main -> #26 -> #49` and
`main -> #50 -> #51 -> #52 -> #53 -> #56`, with #32 independently based on main.
#53's successful runs [37713713477](https://github.com/Aequiludium/cocycle-rs/actions/runs/37713713477)
and [37713709199](https://github.com/Aequiludium/cocycle-rs/actions/runs/37713709199)
are separate from canceled [37713713055](https://github.com/Aequiludium/cocycle-rs/actions/runs/37713713055).
The inventory did not rerun or recertify old performance/oracle archives.

### Publication delta on main cedf596

Main now contains merge commits #32 `7dc517d`, #98 `03e896d`, #50 `168015b`,
#51 `eb3e291`, #52 `5ce4f13` and #53 `cedf596`, verified through main's first-parent
history. #98 carries the compatibility report and executable guide. The merged
[H1 continuation report](h1-fast-to-generic.md) and
[benchmark report index](../../benches/reports/README.md) own the new source/research facts.

The current [dispatch](../../src/persistence/flag/dispatch.rs) runs specialized
F2 H1 and continues generically from H2 with clearing where supported, retaining
odd-prime and overflow fallback. The private H2 benchmark cfg still does not
make the prototype the ordinary production path. Research reports, failed
admission gates and no-go decisions retain their original mathematical meaning
after merge.

The publication-time open PR list is #26 and #49 at the same heads/base
relationship as above, plus draft #56 at
`b08144d717e2de42b02088433b0887bd3a6b1d1f` targeting main. The original #56 ratios
belong to `2dcd5e3`, not this newer head. This docs publication provides no new
performance admission for it. #62/#63 closed through #98; #59–#61 remain open
until their full reports are accepted. All these states are dated observations.

## Keep, supplement or defer

| Decision | Evidence and reason | Owner |
| --- | --- | --- |
| Keep borrowed inputs, frozen owners and four-method FilteredComplex | [Kernel design][kernel] already separates static open cells, sealed default facade and private implicit access. No consumer needs a universal registry. | #60 |
| Keep logical diagrams, dimension/endpoint/context and typed simplex witnesses | Existing #12/#13/#15 are delivered; reuse cheap common-data borrowing and explicit moves rather than repeating them. | #61 |
| Preserve pending optimization ownership | Logical preparation, scoped reuse and low/high-dimensional continuation belong to #14/#26, #31/#49 and #30/#50. Re-query heads before using their evidence. | Original issues / PRs |
| Supplement input responsibilities | Local IDs and signed incidence do not define nonfiltered degree bases, real/complex adjoints, positive cost/Gram semantics, chain-map endpoints or source identity. Edge filtration values are not inner products. | #60; #68/#69 |
| Supplement concrete result interpretation | Common-data construction prevents minting internal source context, but cannot represent every spectrum/projector/query certificate. Empty, absent, unsupported and failed need distinct interpretation. AsRef cannot create provenance. | #61; #78/#79 |
| Recognize algorithm gaps as algorithm work | Alpha/cubical, Laplacian/persistent Laplacian, Dirac, legal F2 projection/min-stretch/transport and generalized-topology mathematics do not arise from interface sketches. | #65/#66; #72–#76 |
| Define consumer domains before new similarities/features | Distances currently accept complete diagrams; descriptors cover two operation families. New spectrum/operator similarity or ML/statistics needs units, fields, identity, arithmetic/error semantics and real conversion costs. | #81/#82; #84/#85; #91/#92 |
| Defer speculative access/import/map shapes | Borrow an existing owner for the first single-scale consumer; narrow additions need demonstrated benefit. Do not create six empty crates, a universal Distance/result container or an identity service. | #60/#61 decisions; #62/#63 gate |
| Retain optimization admission limits | Failed #51/#56 gates, #52/#53 no-go results and #32 resource proxies remain. Old #49 reuse ratios do not establish current-head admission. | Existing research / #96 |

Do not equate a real Hodge Laplacian with F2 `I+P`, force Mayer
`boundary^N=0` into the ordinary `boundary^2=0` contract, or treat stage spectra
as a persistent Laplacian. Those rejected shortcuts lack the required
mathematical definitions and implemented consumers.

## Follow-up source and validation checklist

| Task | Read before implementation | Establish / preserve |
| --- | --- | --- |
| #60 input contracts | [FilteredComplex][filtered-contract], [storage][simplicial], [source context][context], [kernel][kernel]; fixed continuation/access diffs | Owner-local IDs, ordered bases, field/weights/Grams, signed orientation, single-scale/empty/invalid windows, both chain-map equations, source identity and actual conversion cost. |
| #61 result contracts | [diagram][diagram], [data][data], [result][result], [witnesses][representatives], [distance facade][distances]; fixed #26/#49 diffs | Four concrete result families; gaps and computed-empty degrees, endpoint coverage, units/sources, borrow/move/drop, external declarations versus established guarantees, incompatible examples. |
| #62 minimum change | Accepted #60/#61 decisions, [conventions][conventions], [contribution paths][contributions], [CONTRIBUTING][contributing]; fresh main/PR refs | Only demonstrated missing operations, directly affected rustdoc/examples/tests, specialized persistence and legacy execution scope. Its accepted outcome is existing-API reuse with a runnable guide. |
| #63 compatibility | [testing][testing], existing contracts/examples, [reporting][reporting] and relevant native protocols | Construction/query caps, dimensions, signed/zero-born sources, fields, multiplicity, source certificates, witness equations/ordinals and controls. Mathematical changes trigger GUDHI plus Ripser/Topp; performance adoption needs same-version/fixture/worker evidence. |
| Framework and optional routes | The capability map above, [architecture][architecture], individual #64–#94 definitions | Do not repeat delivered PH/diagram/distance/descriptor work. New spectra/operators/similarities need definitions and first examples. Learning remains an optional consumer; a compatibility pass is not framework completion. |

## Original validation and portable reproduction

The original isolated main archive was checked on Windows with Rust/Cargo
1.98.1 and Python 3.10.11, using a separate build directory. It did not switch or
modify the user's feature checkout.

| Check | Original result | Evidence boundary |
| --- | --- | --- |
| `cargo test --locked --all-features` | 196 passed, including 11 doctests; two profiling tests ignored | Existing baseline contracts/oracles run locally; not a new MSRV/release/native/performance acceptance. |
| Nine `cargo run --locked --example ...` calls | square, rips_graph, flag_persistence, rips_sphere, rips_representatives, sparse_rips, diagram_analysis, complex_construction and diagram_distances succeeded | Existing compositions run; missing frameworks remain missing. |
| `cargo test --locked --example complex_construction` | Four passed | Hand circle, ties/orientation, invalid, empty/isolates and multiple-field cases. |
| Source / docs checks | 244 source and 40 Markdown files passed | Baseline source/local-link policy; no external URL or algorithm certification. |
| Hosted checks | Read the exact historical records above | Observed CI, without rerunning it or treating canceled jobs as successes. |
| Local collaboration checker | Five historical missing #32 report links in the unrelated old checkout | Preserved failure, outside tracked baseline docs. Not a failure of the production examples. |

The first existing example gives H0 `[-2,+infinity)`, `[-1,0)` and H1
`[1,+infinity)` from signed/non-zero-born supplied topology. Diagram analysis
gives Betti curve `[0,2,3,1,1]`, finite total persistence 4 and entropy ln(2),
and distinguishes computed-empty H0 from an H1-only domain. Diagram distances
give bottleneck=1, W1=1 and W2=1.4142135623730951 on supplied diagrams.

To reproduce the inventory's production examples, use a fresh checkout of the
fixed commit and run the listed commands. The shell loop below runs all nine:

```sh
cargo test --locked --all-features
cargo test --locked --example complex_construction
for example in square rips_graph flag_persistence rips_sphere rips_representatives sparse_rips diagram_analysis complex_construction diagram_distances; do
    cargo run --locked --example "$example" || exit 1
done
python3 tools/check_source.py
python3 tools/check_docs.py
```

Later main runs can have different test/file counts; bind each result to its
actual commit. Original archives, remote snapshots and raw logs were local
research artifacts, not a public evidence bundle. The fixed source links and
commands expose the public inspection/reproduction surface. This inventory
adds no production algorithm and makes no new native GUDHI/Ripser/Topp,
benchmark, release/MSRV or hosted-CI claim.

[geometry]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/geometry/mod.rs
[graph]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/complex/graph/mod.rs
[simplicial]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/complex/simplicial/mod.rs
[filtered-contract]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/complex/filtered.rs
[filtered-reader]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/persistence/filtered.rs
[simplicial-reader]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/persistence/simplicial/input.rs
[filtration]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/filtration/mod.rs
[rips-builder]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/filtration/rips/builder.rs
[rips-exact]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/filtration/rips/exact.rs
[approx-builder]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/filtration/rips/approximation/builder.rs
[approximation]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/filtration/rips/approximation/metadata.rs
[approx-compute]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/persistence/rips/approximation.rs
[simplicial-filtration]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/filtration/simplicial/mod.rs
[context]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/filtration/context.rs
[persistence-source]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/persistence/source.rs
[dispatch]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/persistence/flag/dispatch.rs
[diagram]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/diagram/persistence_diagram.rs
[data]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/diagram/data.rs
[result]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/diagram/computation.rs
[representatives]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/diagram/representative.rs
[distances]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/diagram_distances/mod.rs
[pr26-distance]: https://github.com/Aequiludium/cocycle-rs/blob/26866ac9223a7bd5cc7c35286e1090dc5759e509/src/diagram_distances/mod.rs
[pr49-distance]: https://github.com/Aequiludium/cocycle-rs/blob/5bb77b4b624580fe3cede4c8d538d06518de1a5d/src/diagram_distances/mod.rs
[descriptors]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/descriptors/mod.rs
[betti-curve]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/descriptors/betti_curve.rs
[lifetime-summary]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/descriptors/lifetime_statistics.rs
[execution]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/execution/mod.rs
[architecture]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/docs/development/architecture.md
[kernel]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/docs/design/kernel.md
[conventions]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/docs/development/conventions.md
[contributions]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/docs/development/algorithm-contributions.md
[contributing]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/CONTRIBUTING.md
[testing]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/docs/development/testing.md
[reporting]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/benches/reporting.md
