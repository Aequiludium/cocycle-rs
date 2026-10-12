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
| Topological computation | Sealed [PersistenceExt / PersistenceBuilder][persistence-source], including `from_complex`, compute ordinary PH over a validated prime-u32 field, default F2. Exact flag/Rips H0 uses union-find; F2 H0/H1 uses specialized implicit cohomology. At this baseline H2+ or odd-prime higher requests use generic cohomology; non-zero-born supplied inputs use boundary reduction. Simplicial sources can request persistent cycles and scale-specific dual cocycles. | `square`, `rips_sphere`, `rips_representatives`, `flag_persistence`; `prime_fields`, `rips_expansion`, `filtered_complex`, `execution` tests and independent private reducer | [Baseline dispatch][dispatch] contains `dimension > 1 || field.characteristic() != 2`. #50's continuation was then pending; see publication delta. #51 fixed-tuple H2 is experimental/private-cfg; #56 remains draft pending direct-graph and H2/H3 performance checks (see the corrected report in the publication delta). No Laplacian/eigensolver, Dirac, F2 homology operator or generalized-topology engine. #72 adapts existing PH interfaces; #73 owns Laplacian/spectra, #74 Dirac, #75 F2 homology and #76 generalized topology; #96 owns adoption decisions. |
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
| [#56](https://github.com/Aequiludium/cocycle-rs/pull/56) | `2dcd5e3d8a306f3cb4c4f23e486cdfc5add78db1` / #53 `c3fee147ca50673d2221d961397fa6c281056da2` | Draft adaptive F2/H1, bounded bitsets and ordered cursors retain H1-to-H2 handoff. The [corrected report][pr56-corrected] withdraws the earlier regression interpretation: both retained studies have lower candidate medians than their #53 baselines on all 11 matrix H1 cases in both repeats. Direct-graph and H2/H3 performance remain unmeasured; correctness CI does not establish those gates. | 16 SUCCESS |

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
`b08144d717e2de42b02088433b0887bd3a6b1d1f` targeting main. Its
[corrected report][pr56-corrected] withdraws the earlier 8.4–9.0× nonmetric and
3.2–3.4× truncated-bipartite regression claims. The original study and independent
rebuild/repeat both show lower candidate medians on every tested matrix H1 case;
the independent study reports candidate/baseline ratios of 0.111–0.938. Its
measured candidate/harness is `2dcd5e3d8a306f3cb4c4f23e486cdfc5add78db1`, against
integrated #53 baseline `323e617246d8ed2bd5096f7dcac2e5a6e17d6de0`. The report
head is not a new measurement. #56 remains draft because direct-graph and H2/H3
performance is unmeasured, not because the measured matrix H1 cases regressed.
This docs publication reruns no performance study and supplies no new admission.
#62/#63 closed through #98; #59–#61 remain open until their full reports are accepted. All these states are dated observations.

## Keep, supplement or defer

| Decision | Evidence and reason | Owner |
| --- | --- | --- |
| Keep borrowed inputs, frozen owners and four-method FilteredComplex | [Kernel design][kernel] already separates static open cells, sealed default facade and private implicit access. No consumer needs a universal registry. | #60 |
| Keep logical diagrams, dimension/endpoint/context and typed simplex witnesses | Existing #12/#13/#15 are delivered; reuse cheap common-data borrowing and explicit moves rather than repeating them. | #61 |
| Preserve pending optimization ownership | Logical preparation, scoped reuse and low/high-dimensional continuation belong to #14/#26, #31/#49 and #30/#50. Re-query heads before using their evidence. | Original issues / PRs |
| Supplement input responsibilities | Local IDs and signed incidence do not define nonfiltered degree bases, real/complex adjoints, positive cost/Gram semantics, chain-map endpoints or source identity. Edge filtration values are not inner products. | #60; #68/#69 |
| Supplement concrete result interpretation | Common-data construction prevents minting internal source context, but cannot represent every spectrum/projector/query certificate. Empty, absent, unsupported and failed need distinct interpretation. AsRef cannot create provenance. | #61; #78/#79 |
| Recognize algorithm gaps as algorithm work | Alpha/cubical, Laplacian/persistent Laplacian, Dirac, legal F2 projection/min-stretch/transport and generalized-topology mathematics do not arise from interface sketches. | #65/#66; #73–#76 |
| Define consumer domains before new similarities/features | Distances currently accept complete diagrams; descriptors cover two operation families. New spectrum/operator similarity or ML/statistics needs units, fields, identity, arithmetic/error semantics and real conversion costs. | #81/#82; #84/#85; #91/#92 |
| Defer speculative access/import/map shapes | Borrow an existing owner for the first single-scale consumer; narrow additions need demonstrated benefit. Do not create six empty crates, a universal Distance/result container or an identity service. | #60/#61 decisions; #62/#63 gate |
| Retain optimization admission limits | The failed #51 admission gate, unmeasured #56 direct-graph/H2/H3 performance gates, #52/#53 no-go results and #32 resource proxies remain. Old #49 reuse ratios do not establish current-head admission. | Existing research / #96 |

Do not equate a real Hodge Laplacian with F2 `I+P`, force Mayer
`boundary^N=0` into the ordinary `boundary^2=0` contract, or treat stage spectra
as a persistent Laplacian. Those rejected shortcuts lack the required
mathematical definitions and implemented consumers.

## Computation-family interface supplement

At main `a7ee9e911dcdeaea137b682bf383f7c2e1d090d5`, the 2026-10-11 current
stage of #73–#76 is interface preparation and model investigation. The
[input supplement](input-contracts.md#computation-family-preparation-on-main-a7ee9e9)
records exact source reading, retained owners, narrow proposed operations,
coefficient/Gram/cost conditions and independent hand examples. The
[output supplement](result-contracts.md#computation-family-outputs-and-resources)
records ownership, selections, numeric evidence, feasibility/optimality,
errors and resource scope. No production API/kernel or pending PR changes.

Laplacian/Dirac/F2 algorithms remain unimplemented on that main. The generalized
comparison selects a one-edge real sheaf as the next investigation and defers
path, embedded hypergraph and Mayer implementations with model-specific reasons.
The candidate choice completes a research decision, not a generalized engine.
Stage-B source, independent checks, costs and running integration remain separate.

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
[pr56-corrected]: https://github.com/Aequiludium/cocycle-rs/blob/b08144d717e2de42b02088433b0887bd3a6b1d1f/benches/reports/f2-h1-adaptive.md

## Legacy migration decisions on main a7ee9e9

Research completion for [#71](https://github.com/Aequiludium/cocycle-rs/issues/71),
2026-10-11. This section extends the historical inventory above with file-level
dispositions; it does not repeat #59-#63 or adopt pending implementations.
The old plan base is `9e6715f4c2e0118ab738486bd747ec3e361397de`; the inspected,
freshly fetched implementation base is
`a7ee9e911dcdeaea137b682bf383f7c2e1d090d5`. Subsequent implementation PRs must
record their own base/head and validation. A completed migration decision is
distinct from its implementation, review, merge or release.

### M1: files, symbols, consumers and disposition

Actions below mean **retain in place**, **move internally**, **compatibility
forward**, **delete obsolete private duplicate**, or **defer**. Retention is a
boundary-specific decision with real consumers, not a blanket claim that
internal refactoring is unnecessary. No public deletion is selected.

| Files and symbols | Actual consumers / existing checks | Disposition, target owner and reason |
| --- | --- | --- |
| `geometry/{point_cloud,dissimilarity,matrix,metric}.rs`: `PointCloudView`, `DissimilarityView`, `DissimilarityMatrixView`, `validate_metric` | Rips builders, exact point/matrix computations; `tests/{contracts,matrix,euclidean,rips_api}.rs` | Retain in geometry: borrowed input validation and metric evidence already precede computation. No second input container (#65/#66). |
| `complex/graph/{mod,adjacency}.rs`: `WeightedGraph`; `complex/simplicial/{mod,simplex,incidence}.rs`: `SimplicialComplex`, `SimplexId`, `boundary` | Flag/Rips expansion, frozen PH and critical sets; `tests/{graph,filtered_complex,critical_sets}.rs`, construction example | Retain in complex: one frozen owner, oriented incidence and local IDs; graph-to-clique interpretation stays explicit. #100 owns stage additions. |
| `complex/filtered.rs`: `FilteredComplex::{cells,dimension,filtration_value,boundary}` | `PersistenceBuilder::from_complex`, frozen implementation; `tests/filtered_complex.rs` | Retain the four-method static contract. #100 owns any stage-associated increment; no universal solver/source trait. |
| `filtration/rips/{builder,exact,expansion}.rs`: `RipsBuilder`, `ThresholdRips`, `RipsExpansion`; `filtration/flag/` and `filtration/expansion/` | Persistence source extension, explicit/implicit access and construction examples; `tests/{rips_construction,rips_expansion,rips_api,flag}.rs` | Retain in filtration: construction and access do not call PH. Requested cap, retained edges and original source coverage differ. |
| `filtration/rips/approximation/{builder,metadata,access,blocker,expansion}.rs`: `ApproximateRipsBuilder`, `SparseRips`, `RipsApproximation` | Approximate PH and explicit sparse expansion; `tests/{sparse_rips,rips_api}.rs` | Retain blocker, labels and approximation hypotheses in their existing owner; a modified graph alone cannot reconstruct this topology. |
| `persistence/builder.rs`: `PersistenceBuilder::{from_complex,compute,compute_with}` | All source extensions and external callers; `tests/{rips_api,execution,filtered_complex}.rs` | Retain request validation, borrowed lifetimes and one terminal budget. Existing function-pointer hook needs no registry or wrapper type (#72). |
| `persistence/source.rs`: `PersistenceExt`, sealed implementations and `source!` | Public `.persistence()` calls; compile-fail rustdoc, `tests/rips_api.rs` | Retain consumer-owned sealed wiring here; move the four substantial adapters listed below to their existing mathematical owners (#72 batch 1). |
| `source.rs::exact` | `RipsBuilder` hook; point/matrix layouts, independent caps, cancellation and representatives in `tests/rips_api.rs` | Move body to `persistence/rips/mod.rs::compute_builder_budget`. Input preparation, effective cap and exact-Rips context belong with the legacy Rips entry points. Preserve the current calls/checkpoints verbatim. |
| `source.rs::approximate` | `ApproximateRipsBuilder` hook; `tests/{sparse_rips,rips_api}.rs` | Move to `persistence/rips/approximation.rs::compute_builder_budget`, alongside prepared-source computation. No new preparation or metadata clone. |
| `source.rs::{explicit,supplied}` | `SimplicialFiltration` / `SimplicialComplex` hooks; signed/delayed/non-flag cases in `tests/filtered_complex.rs` | Move to `persistence/simplicial/mod.rs::{compute_filtration_budget,compute_supplied_budget}`. Range/skeleton checks, compute and result assembly then have one existing owner (#72 batch 1). |
| `persistence/filtered.rs::{read,compute}` | `from_complex`; field-dependent cellular and malformed external inputs in `tests/filtered_complex.rs` | Retain validating external reader: ID/order/incidence and selected-field boundary-square checks cannot be replaced by trusted-owner assumptions (#72). |
| `persistence/simplicial/input.rs::read` | General frozen computation and `representatives/basis.rs`; low-degree selection tests | Retain trusted frozen reader: validates at construction, selects only q+1/cutoff columns and maps selected owner IDs. Do not impose external validation costs (#72). |
| `persistence/boundary/{mod,input}.rs`: `BoundaryInput`, `diagram`; `algebra/reduction/boundary.rs::{reduce_pairs,reduce}` | External/frozen readers, representative and critical-set consumers; `boundary/tests.rs`, independent D V = R check | Retain existing shared reduction. `reduce_pairs` omits V while `reduce` retains it intentionally. Neither is an obsolete duplicate (#72). |
| `persistence/flag/{dispatch,h0}.rs`, `flag/cohomology/`, `union_find.rs`: `compute_dense`, `compute_graph`, F2 H1 | Legacy/builder matrix and graph routes; dispatch/cohomology tests, `tests/{flag,prime_fields,rips_expansion,rips_resources}.rs` | Retain admitted specialized routing, full clearing keys, H1-only no-clearing mode and workspace lifetime. #56 exclusively owns adaptive candidate changes; no default adoption (#72). |
| `persistence/simplicial/cohomology.rs`: `compute`, higher-degree continuation | Frozen zero-born and implicit sparse/odd-prime/overflow fallback paths; cohomology tests, independent bitmask oracle in `tests/rips_expansion.rs` | Retain generic coface kernel and private H2 cfg boundary. Arbitrary births/signed cutoffs still route to boundaries; generic fallbacks are live (#72). |
| `persistence/rips/{mod,expanded,approximation}.rs`, `persistence/flag/mod.rs`: all re-exported `compute_*`, `rips_from_*` | Existing public free functions and source hooks; `tests/{rips,rips_api,rips_construction,sparse_rips,flag}.rs` | Retain public paths and legacy budget scope. Batch 1 adds internal adapter owners, not duplicate public terminals or a forced builder conversion (#72). |
| `persistence/simplicial/representatives/{request,basis,complex,dual}.rs`: `RepresentativeRequest`, basis extraction | Opt-in Rips/flag/frozen witnesses; `tests/{prime_fields,rips_api,filtered_complex}.rs` | Retain request semantics and independent closure/rank/duality tests. No default materialization; basis transformations have a real consumer (#72). |
| `persistence/mod.rs::assemble_diagram`, `RawIntervals` | Boundary, legacy Rips, simplicial finish, representative basis and test/reference/profile callers | Move internally in #72 batch 2 to a small private assembly owner if the post-batch-1 review confirms it: canonical zero-lifetime/coverage handling has multiple live callers and currently lives in the export facade. Migrate imports together; no second implementation. |
| `diagram/{persistence_diagram,dimensions,interval}.rs`: logical views, ordinals, `ComputedDimensions`, `IntervalEnd` | PH assembly, distances and descriptors; `tests/{contracts,descriptors,diagram_distances}.rs` | Retain result-layer reads: computed-empty/absent, multiplicity, canonical ordinals and censoring already have consumer tests (#79-L). Public `intervals()` compatibility remains supported. |
| `diagram/{data,computation,representative}.rs`: `PersistenceData`, `PersistenceResult`, `ComputationContext`, `into_parts` | PH adapters and contextual distance facade; `tests/{rips_api,prime_fields,diagram_distances,execution}.rs` | Retain borrowed data and owned result assembly; move no context to algebra. Explicit dropping conversions retain their documented losses (#79-L). |
| `diagram_distances/mod.rs`: raw / `_results` / controlled facades and `points`; bottleneck/wasserstein solvers | Distance example, `tests/{diagram_distances,execution}.rs`, exhaustive tiny matching | Retain public paths and solver ownership (#82-L); defer logical/prepared extraction changes to #14/#26 and #31/#49. Do not duplicate or adopt those paused branches. |
| `descriptors/{betti_curve,lifetime_statistics}.rs`: `betti_curve`, `finite_lifetime_summary` | Diagram-analysis example/walkthrough, `tests/descriptors.rs` | Retain diagram-only reads (#85-L): no PH call, expected allocation differs by descriptor. New features/ML adapters stay in #85-N. |
| `optimization/{mod,tests}.rs`: `CriticalSetWorkspace`; `algebra/{column,field}/` | Critical-set example, `tests/critical_sets.rs`, opt-in optimization tests | Retain sparse F2 critical sets as optional computation. Shared algebra does not supply a homology-preserving projector; #75 owns that increment. |
| `execution/mod.rs`, `persistence/execution.rs`, `error.rs`, `lib.rs` | Existing domain exports/terminals; `tests/{execution,rips_resources,contracts}.rs` | Retain current public compatibility forwards and structured errors. `WorkBudget` is fresh per terminal, not an aggregate quota, RSS limit or concurrency guarantee. |
| `persistence/reference/`, private profiling/cfg modules | Independent small oracles and explicit benchmark workers | Retain independent reference mathematics and experiment isolation. The caller search found no proven dead private computation to delete; similar formulas are not sufficient deletion evidence. |

Paths in this table are relative to `src/` unless prefixed with `tests/`.
The actual PH chain is extension/free terminal -> source-specific preparation
and range interpretation -> existing dispatch/coface or boundary kernel ->
diagram/witness assembly -> owned context result. The two boundary readers stay
separate, and a representative request remains an explicit different-cost path.

### M2: merged changes and parallel ownership

Main's first-parent history from the old plan base contains #32 (`7dc517d`),
#98 (`03e896d`), #50 (`168015b`), #51 (`eb3e291`), #52 (`5ce4f13`),
#53 (`cedf596`), #101 (`eb38e52`) and #99 (`a7ee9e9`). This is a 71-file delta:
resource reports, interface decisions, H1-to-H2 continuation, private H2 research,
workspace/no-go records, opt-in critical sets and final research publication.
It introduces no production Laplacian, Dirac, F2 projector or generalized engine.
The shared `BoundaryInput`/`diagram`/`reduce_pairs` already existed at the old base.

Live open PRs were checked again on 2026-10-11, including changed files:

| PR / sole implementation owner | Actual base SHA / head SHA | Overlap and migration rule |
| --- | --- | --- |
| #100 / stages, #68/#69 | `a7ee9e911dcdeaea137b682bf383f7c2e1d090d5` / `6f2b291067e4a834518c9f055dcc991b44c781ee` | `persistence/source.rs`, frozen/stage modules, architecture and filtered tests overlap batch 1. Move only main's four adapters; leave new stage wiring with #100. Merge whichever lands second normally and run stage tests then. No stage implementation is copied into this work. |
| #102 / modeling example, #65/#66 | same main base / `0cff0f0bbd2daf58d7391f31fe17bef5fdfc6963` | Docs/index/guide and construction example, no production source. Batch 1 need not wait; example/doc conflicts belong to ordinary integration. |
| #26 / logical distance preparation, #14 | same main base / `8dcc8d885d93c0133988a47865f0ea3aaba740b3` | Distance facade/solvers/workers/tests; no first-batch PH source overlap. Preserve existing review and cost decisions. |
| #49 / opt-in operands, #31 | `8dcc8d885d93c0133988a47865f0ea3aaba740b3` / `c30ee1ad353c144dc9b415caeca3d7ffc9fd4b9f` | Preserve `main <- #26 <- #49`; #82-L does not independently rewrite preparation. |
| #56 / draft adaptive H1 | `cedf5965b774855614b761a1ad311d9dc7992c49` / `b08144d717e2de42b02088433b0887bd3a6b1d1f` | Dispatch/cohomology, source access and contribution docs overlap later kernel work. Candidate stays draft and outside both migration batches; current main routing is the baseline. |

These are snapshots, not continuing mergeability or CI claims. The source cap,
stage owner and logical-distance preparation have one implementation authority
each. Refresh these refs/files before each submission.

### M3: behavior and cost gates before editing

The pre-edit `cargo test --locked --all-features` on main a7ee9e9 passed
216 tests including doctests; four diagnostic tests remained ignored.
Source and Markdown checks passed. This is local Windows correctness evidence,
not a fresh native performance/oracle result. Primary local collaboration checks
are reported separately; copied local artifact links are not portable evidence.

[#87](https://github.com/Aequiludium/cocycle-rs/issues/87) G0/G1 and
[#94](https://github.com/Aequiludium/cocycle-rs/issues/94) A/B use this matrix
inside each migration PR; neither entire issue must close before implementation.

| Gate | Existing fixtures / command scope | Cost boundary and required response |
| --- | --- | --- |
| Public paths, request lifetimes, per-terminal errors | `rips_api`, source compile-fail doctests, `execution`, `rips_resources`; old free functions plus builder sources | Builder preparation/reduction/witnesses share one budget. Legacy point preparation retains its documented narrower scope. No extra wrapper invocation/budget reset. |
| Exact matrix/points/graph and H0-H3 | `rips`, `flag`, `prime_fields`, `rips_expansion`, dispatch/cohomology independent oracle tests | Preserve dense/sparse/odd-prime/overflow routes, virtual death pairs, H1-only clearing suppression, H1 workspace end before H2. Dispatch/access/lifetime changes trigger same-fixture base/head pipeline measurements. |
| External cells versus frozen arbitrary births/non-flag sources | `filtered_complex`, `boundary` tests and construction example | Construction is separate from input read/validation, reduction and assembly. Trusted read skips external boundary-square validation; low-degree requests retain only selected ID maps/columns. Allocation/validation changes require separate read and terminal timings. |
| Sparse approximation and representatives | `sparse_rips`, `prime_fields`, `rips_api`, representative closure/duality checks | Retain blocker-aware access, original labels, metric/provenance and opt-in q+1/V storage. Changed metadata clones or witness preparation require representative and diagram-only controls. |
| Multiset/availability/coverage/context and downstream reads | `contracts`, `descriptors`, `diagram_distances`, `rips_api`; `square`, `diagram_analysis`, `diagram_distances` examples | Preserve finite/essential/censored, computed-empty/absent, field/context and witness ordinals. Borrow, move, discard and clone differ. Distance/descriptor consumption never reruns PH. |
| Engineering and independent scope | CONTRIBUTING debug/release, examples/guide doctests, format/Clippy/docs, MSRV and package checks | Pure body relocation needs these checks and a statement of unchanged operations/call topology. Any semantic change needs independent expectations plus GUDHI/Ripser where supported; distance changes additionally need Topp. Unsupported references are exclusions. |

For a cost-affecting batch use [pipeline boundaries](../../benches/pipeline/README.md)
and [reporting rules](../../benches/reporting.md) at the exact base/head with the
same seed/input/field/cutoff, alternating cold/warm repetitions and retained raw
failures. Isolate input preparation/read/validation, kernel, assembly and requested
witnesses. Record RSS separately from capacity/work counts. Reuse the admitted
H1/H2 fixtures and low-/high-degree controls; old timings do not certify new code.
No speedup or zero-allocation target is invented for mechanical migration.

### M4: executable small-PR queue

| Batch / owner | Modification and preserved path | Acceptance and rollback |
| --- | --- | --- |
| #72 batch 1: source adapter ownership | Move `exact`, `approximate`, `explicit`, `supplied` from `source.rs` into three existing Rips/simplicial files, renaming private hooks; update only extension wiring and affected architecture/contribution notes | Body/statement comparison against a7ee9e9, unchanged function-pointer signature/checkpoint/call sequence; full current contract and independent reducer tests, examples, lint/docs/MSRV/package. Public paths and benchmark workers stay callable. Revert that PR alone; no stored format changes. |
| #72 batch 2: interval assembly facade | Reassess `assemble_diagram` and `RawIntervals` after batch 1; if retained as a distinct assembly owner, move the single helper and all live imports together, with source-linked documentation | Existing multiset/endpoint/reference tests and standalone workers; no duplicate helper or new public type. Any assembly allocation or routing change requires measured controls; otherwise keep the body unchanged. Revert the batch independently. |
| #79-L: result consumers | Retain the existing diagram/data/result files and reads listed in M1; any demonstrated stale caller is fixed in this owner, with explicit ordinal/loss checks | `contracts`, witness moves/contextual distances; no pending preparation adoption. Source/tests currently support a boundary-specific no-op. |
| #82-L: distance consumers | Retain facade/solver ownership; unique preparation changes remain #26/#49 | Existing exhaustive/public contracts; measure extraction/solve if changed. No preparation clone is introduced by #72. |
| #85-L: descriptors / optional optimization | Retain diagram-only descriptor files and separate CriticalSetWorkspace | `descriptors`, analysis example, critical-set tests; measure changed sorting/copy/allocation only. Future learning is outside this migration; its separate owner is #111. |

No dead private duplicate has been proven. In particular, test oracles, generic
fallbacks, basis-retaining reduction and experimental cfg code retain distinct
consumers or mathematical/cost reasons. This decision does not authorize deleting
them to manufacture a refactor. All PH/common M decisions are now sufficient
to start #72; the other three owners retain their own implementation acceptance.

### N: future computation contracts outside the migration scope

The [input](input-contracts.md) and [result](result-contracts.md) reports already
fix signed bases, Grams, projection certificates and finite examples at their
listed upstream revisions. Their external runtimes/older experiments are not
rerun by this decision. Under the current [Phase-4 scope](https://github.com/Aequiludium/cocycle-rs/issues/57),
#71 delivers only the M1-M4 legacy migration decisions. The candidates below
describe future inputs, outputs and mathematical acceptance conditions, not
algorithm implementation or completion of #73-#76 in this refactoring phase.
#73-#75 prepare interface information independently from #60/#61; they do not
wait for each other's algorithms. Future implementation needs its own explicit
task and acceptance, outside #72/#79/#82/#85 and the #94 refactoring exit.

| Family / future candidate | Reusable code and missing mathematical data | Proposed input/output, independent validation and cost / future responsibility |
| --- | --- | --- |
| Ordinary PH: continue internal migration | Existing prime field, cofaces, two readers, boundary reducer, interval/witness owners are sufficient; no new backend needed | Signed delayed triangle: births -4/-3/-2, edges 0/1/2, face 3 -> H0 (-4,infinity), (-3,0), (-2,1), H1 (2,3). Existing hand-derived boundary test plus independent explicit/bitmask oracle; GUDHI/Ripser for supported route changes. Implicit diagram-only access stays implicit; representatives remain optional. #72 batches above. |
| Single-scale Laplacian: future minimal unweighted real action; no current algorithm or eigensolver adoption | Reuse graph or frozen **signed integer** incidence and owner/basis map, not modular PH coefficients. Choose ordinary nonaugmented chains and Euclidean Gram I. Need f64 action/output, overflow/nonfinite/residual policy; no current production spectral owner | A path graph as a 1-complex has B1=[[-1,0],[1,-1],[0,1]], L1=B1^T B1=[[2,-1],[-1,2]], spectra {1,3} by hand. Triangle ring has L1 spectrum {0,3,3}; filled triangle L1=3I. #73 currently documents the proposed action/matrix/spectrum interface. Borrowed sparse action, independent matrices and any numerical solver belong to separately accepted future implementation. Action is O(nnz B); tiny dense validation O(n^2), general dense spectrum O(n^3) time/O(n^2) storage. Default PH incurs none of it. |
| Dirac: future graded action with independently documented real conventions; defer persistent Dirac | Reuse ordered adjacent signed boundaries and Gram I. For a finite graded chain choose D=B+B^T on the direct sum with zero blocks; require B^2=0 over the reals and keep degree offsets | For path C0+C1, D=[[0,B1],[B1^T,0]], D^2=diag(L0,L1), spectrum {0,+/-1,+/-sqrt(3)}. Filled triangle includes C2 and square includes all Hodge blocks. Check block multiplication exactly, spectral +/- pairing and empty-degree shapes with independent matrices. #74 currently documents graded action/matrix/spectrum requirements alongside #73, without an algorithm dependency. Future sparse action is O(nnz B); full export is quadratic; any numerical scheme needs separate implementation acceptance. No fabricated filtration or F2 I+P substitution. |
| F2 homology operator: future bounded feasible projection; defer global optimality/general geometry | Existing frozen incidence/field can extract A:m-by-n and D:n-by-p with AD=0 mod 2, source-local ordered bases and **separate** positive costs. Current PH pairs/critical sets do not provide legal P/certificates. Local candidate code is not main or automatically authorized for publication | Ring edge space: A is signed triangle incidence reduced mod 2, D has zero columns; cycle 111 costs (2,3,5) -> selected mass 10. Filled ring D=111 -> beta1=0, selected class mass 0. Require AP=0, PD=0, P^2=P and z+Pz in im(D) for every z in ker(A), with rank beta=dim ker(A)-rank D. Exhaustive tiny binary linear algebra checks feasibility independently; homology-operator `a8d03d1` supplies a comparison contract, not a runtime dependency. #75 currently documents explicit budget/failure/certificate output requirements; projection and class-query solvers are future implementation work. Dense extraction O(mn+np), factorized/exhaustive solving bounded and potentially large; feasibility never implies optimality/shortest mass. |
| Generalized topology: later research candidate in #76; no current production adoption | Sheaf restrictions/stalk dimensions, directed allowed path chains, hypergraph embedded chain spaces and Mayer N-differentials are missing. Shared column arithmetic may be reused only after a theory's chain law is established; ordinary `FilteredComplex` is insufficient | Sheaf first case: two vertices/one edge with declared stalk maps -> coboundary and Laplacian; directed path first case: 0->1->2 requires checking whether its boundary stays in the allowed chain space; hypergraph case: one 3-vertex hyperedge lacks automatic face closure; Mayer case must specify N and generalized homology for d^N=0. Future #76 research would compare at least three candidates against original definitions and independent finite matrices before a separate implementation decision. Costs depend on stalk dimensions/allowed paths/N and cannot inherit PH benchmarks. No universal generalized backend is adopted. |

The generalized deferral has concrete obstructions: nontrivial sheaf restrictions
are not scalar +/-1 simplex incidence; the boundary of the allowed directed path
(0,1,2) includes (0,2), which need not be allowed; a hyperedge is not automatically
its downward-closed simplicial complex; and d^N=0 for N>2 does not satisfy the
ordinary d^2=0 reducer contract. These justify postponing integration rather than
renaming existing PH. Original-paper/full-runtime evaluation remains future #76 work.
The references in #76 are reading candidates, not sources claimed newly read here.

The PH migration queue in M1-M4 is this #71 delivery. The finite-chain examples
below only illustrate proposed contracts and independent hand-derived expectations;
they do not select, implement or adopt a production algorithm. #73-#75 retain
future algorithm goals, but only prepare interface information during the current
refactoring phase; actual implementation requires separately scoped future work.
#76 remains a later research candidate. These issues remain open independently
of #71's migration decision and cannot block #94 through unimplemented algorithms.
Parallel reducers, unsafe acceleration, new dependencies and persistent spectral
claims are not selected.

### Reproduce the first finite definitions

This standalone Python standard-library check was run on 2026-10-11. It verifies
the stated path/filled-triangle matrices, Dirac square and tiny feasible F2
projectors independently of production PH and any spectral/projector runtime.
The eigenvectors and block matrix are hand-derived expectations. Exhaustive
binary vectors check all projector identities in these three-coordinate spaces.
This is finite definition validation, not implementation or general optimality
evidence. Save the following as `target/migration-definitions.py` and run
`python target/migration-definitions.py`:

```python
from itertools import product


def transpose(a):
    return [list(row) for row in zip(*a)]


def multiply(a, b):
    return [[sum(x * y for x, y in zip(row, col))
             for col in transpose(b)] for row in a]


path = [[-1, 0], [1, -1], [0, 1]]
assert multiply(transpose(path), path) == [[2, -1], [-1, 2]]
# Check the two hand-derived eigenvectors, without running an eigensolver.
assert multiply([[2, -1], [-1, 2]], [[1], [1]]) == [[1], [1]]
assert multiply([[2, -1], [-1, 2]], [[1], [-1]]) == [[3], [-3]]
dirac = [[0] * 5 for _ in range(5)]
for i in range(3):
    for j in range(2):
        dirac[i][j + 3] = dirac[j + 3][i] = path[i][j]
assert multiply(dirac, dirac) == [
    [1, -1, 0, 0, 0], [-1, 2, -1, 0, 0], [0, -1, 1, 0, 0],
    [0, 0, 0, 2, -1], [0, 0, 0, -1, 2]]
triangle = [[-1, -1, 0], [1, 0, -1], [0, 1, 1]]
face = [[1], [-1], [1]]
assert multiply(triangle, face) == [[0], [0], [0]]
down = multiply(transpose(triangle), triangle)
up = multiply(face, transpose(face))
assert [[down[i][j] + up[i][j] for j in range(3)] for i in range(3)] == [
    [3, 0, 0], [0, 3, 0], [0, 0, 3]]


def action(a, z):
    return tuple(sum(x * y for x, y in zip(row, z)) % 2 for row in a)


cycles = {z for z in product((0, 1), repeat=3)
          if action(triangle, z) == (0, 0, 0)}
assert cycles == {(0, 0, 0), (1, 1, 1)}
ring_projection = [[0, 0, 1]] * 3
filled_projection = [[0] * 3 for _ in range(3)]
for projection, boundaries, expected in [
    (ring_projection, {(0, 0, 0)}, (1, 1, 1)),
    (filled_projection, cycles, (0, 0, 0)),
]:
    for z in product((0, 1), repeat=3):
        pz = action(projection, z)
        assert pz in cycles and action(projection, pz) == pz
        if z in cycles:
            assert tuple(x ^ y for x, y in zip(z, pz)) in boundaries
    assert all(action(projection, d) == (0, 0, 0) for d in boundaries)
    assert action(projection, (1, 1, 1)) == expected
print('Exact path/triangle/Dirac and exhaustive F2 projection checks passed.')
```
