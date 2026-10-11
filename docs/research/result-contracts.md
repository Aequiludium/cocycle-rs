# Reading results and establishing comparison scope

[Documentation](../README.md) / Research / [#61](https://github.com/Aequiludium/cocycle-rs/issues/61)

Technical report, 2026-10-09, following the [interface inventory](interface-inventory.md)
and [shared input contracts](input-contracts.md). This publishes the full result
decision behind the [interface compatibility report](../design/interface-compatibility.md).
The fixed-source findings remain dated observations, independent of publication,
subsequent main changes and later mathematical feature delivery.

## Decision and production boundary

Keep separate mathematical result owners. Reuse `PersistenceData` only for
ordinary persistence results with compatible diagram/context semantics. Read
logical intervals, existing borrowed data, scalar queries and requested actions
directly; materialize a complete representation only when the consumer needs it.
For future spectral/projection operations, define their own narrow output beside
the concrete consumer, including its established domain and evidence. Do not
introduce a universal output container, a shared `exact` flag, or a mandatory
full matrix/eigenbasis/diagram conversion.

Some descriptions can use common vocabulary: the source snapshot, method/source
revision, arithmetic, declared units and actual validation. They remain facts
owned by the source/request/operation, not a new metadata service or a fabricated
certificate. Per-result guarantees, query states and witnesses remain specific.
Each multi-input consumer validates every required operand, with a comparison
policy that says which facts must agree and which may differ.

This research introduces no production result type, spectral engine, homology
operator, serializer, identity checker or new Rust signature. The finite probes
call existing main APIs and fixed external reference code. #62 owns minimal
implementation and #63 owns compatibility acceptance; future spectral/Dirac/F2
capabilities retain their separate mathematics and validation tasks.

## Fixed sources, reading and choices

| Source | Revision and material actually read | Chosen lesson / remaining boundary |
| --- | --- | --- |
| Cocycle main | `9e6715f4c2e0118ab738486bd747ec3e361397de`; [kernel design][kernel], [diagram analysis][analysis], [mathematics][math], [diagram][diagram], [endpoints][endpoints], [common data][data], [result/context][result], [source scale][context], [distance facade][distances], [context/dimension tests][distance-tests] | Keep logical access, explicit dimension sets, endpoint coverage, owned common data and explicit context-dropping moves. Existing field/parameter checks do not authenticate physical units or dataset identity. |
| homology-operator | `a8d03d1699307ce1997b3c7e81a1946bb83f5408`; [result model][operator-model], [result source][operator-result], [readouts][operator], [solver source][solver], [solver contract][solver-contract], ChainWindow and validation source | Read one legal P with its basis/weights/run; separate availability, arithmetic, feasibility, objective and optimality. Take the semantics, not its entire Python snapshot schema or six-identity policy as a mandatory Rust design. |
| Persim | `096a25fe38644d31c6aabb880a24a05046c03559`; README, [bottleneck][persim-bottleneck] and [Wasserstein][persim-wasserstein] source | A bare scalar/matching does not retain units, field, dimensions or discarded points. Fixed-source Wasserstein uses W1 with Euclidean ground cost, whereas Cocycle W1 uses L-infinity. Do not infer a metric from the shared name. |
| TopoNetX | `ce3d414b571fff95916c4cbc20b8d2c2c9d428e6`; [spectrum implementation][spectrum] and the input report's incidence/Hodge inspection | Preserve the eigenproblem and requested/returned selection. Repeated eigenvectors, rounding and ordinary versus generalized problems need explicit checks. The probe ran two extracted functions; it did not install or accept the whole package. |

At the original research run, remote main was rechecked and equaled this baseline.
The main archive and fixed homology source are reused without edits. New external
files were fetched by fixed ref, with Git blob/SHA256 records. External source
checkouts, pending PRs and frozen experiment evidence are preserved.

## Facts that can be described together, and facts that cannot

| Fact | Authority and rule |
| --- | --- |
| Dataset namespace/snapshot and original construction | Source/caller evidence. Source kind, vertex count, matching matrix hashes or an asserted provenance string do not authenticate the same dataset. Different datasets can legitimately be compared if the operation's other premises hold. |
| Degree and requested/established domain | Operation plus result. PH has computed dimensions and coverage; spectra have a matrix degree and full/partial spectral selection; projections have their chain window; a distance retains its selected operands in the call. Never derive availability from output length alone. |
| Units and normalization | Source/model/operation. Distinguish filtration scale, Gram/adjoint convention, coordinate costs, eigenvalue units and distance metric/order. `EdgeLength` is a parameter convention, not metres or an equality-of-units certificate. `None` means unknown, not dimensionless. |
| Arithmetic and precision | Exact F2 coefficients, exact rational/integer costs and binary64 calculations describe different facts. Preserve actual rounding/tolerance/residual information; missing error estimates stay missing. |
| Method, version, configuration, execution | Record what was actually used; fixed source/version is reproducibility context, not a correctness proof. A new run can produce the same mathematical content with a different run identity. |
| Witnesses and certificates | Stay with the result they establish. Interval associations, oriented basis coordinates, eigenvectors/subspaces, projection feasibility and optimality proofs are not interchangeable. |

Only operation-relevant fields are required. Spectral consumers do not invent PH
coverage; PH results do not acquire Grams; an F2 projection is not a real
orthogonal projector. A result may outlive borrowed source inputs while retaining
owned source facts, but ownership alone establishes no provenance.

## Four concrete reading examples

### Ordinary persistence diagrams

Use the existing public owner and bind a logical view before retaining its iterator:

```rust,ignore
let data = result.as_ref(); // &PersistenceData; existing data is borrowed
let h1 = data.diagram().dimension(1)?;
for interval in h1.iter() {
    let _ = (interval.birth(), interval.end());
}
let _ = (data.diagram().computed_dimensions(), data.diagram().coverage());
let _ = (data.context().characteristic(), data.context().filtration().scale());
```

This is a reading fragment; the [interface guide](../guides/interface-contracts.md)
supplies a runnable owner and exercises borrowing and witness-preserving moves.
The diagram `{1,3}` with two copies of H1 `[-2,0)` has two logical intervals in
H1, a valid empty H3, and missing H0/H2/H4. Traverse the declared set, not
`0..=max_dimension()`. Keep multiplicity and canonical ordinals. Access gives
logical values and promises no stable physical row addresses or storage layout.

Finite intervals are half-open, essential intervals are justified only in the
complete supplied source, and censored intervals include their stated cutoff.
The truncated triangle's class at 1 remains censored even if a later stage or
an exported matrix makes another mathematical question computable. A manual
diagram constructor checks structural consistency, not the caller's mathematical
origin or truthful declaration of complete coverage.

### Spectral values and vectors

For a future consumer, first read degree, operator definition, field, bases,
Grams/adjoints, eigenproblem/normalization, eigenvalue units and selection. Then
read the requested values and optional vectors; a value-only consumer need not
request an entire eigenbasis. Match each vector to its value and source basis.

The direct path in the input report has ordinary real `L=[[2,-1],[-1,2]]`, hand-derived
eigenvalues `(1,3)` and normalized columns `(1,1)/sqrt(2)` and `(1,-1)/sqrt(2)`.
Its binary64 residual norm is 0 in this fixture and orthogonality error is about
`2.44e-16`. Those are finite numerical observations with a hand-derived answer,
not a general error certificate. A partial selection containing only 3 can
omit a zero eigenvalue; it cannot prove zero Betti number or a full spectral gap.

Signs are arbitrary for simple eigenvectors; an orthonormal basis may rotate
within a repeated eigenspace. Compare appropriate subspace projectors/principal
angles after a declared coordinate identification, rather than raw columns.
Use G-normalization for a G-inner product. An explicit basis transformation must
also transform the operator and relevant Grams; equal array shapes are insufficient.

The fixed TopoNetX small branch returns all eigenpairs even when `n_components=1`.
For `3I` in dimension 3, its dictionary keyed by rounded eigenvalues returns three
copies of the same vector: rank 1, zero residual, failed orthogonality/completeness.
In dimension 11, its large branch instead solves `L v = lambda Diag(L) v`, giving
values near 1 for `L=3I`; the full ordinary spectrum function gives eleven 3s.
Keep these fixed-function findings and problem definitions explicit. No fix to
that repository or broad claim about other versions is part of this task.

### F2 projections, actions and query records

The pinned reference supports direct reads without exporting a whole matrix:

```python
q = op.readout("selected_mass", (1, 1, 1))
assert q.state == "Computed" and q.value == 10
assert q.identity == op.identity
record = op.to_result()  # unread kernel/stretch remain NotComputed
```

For the input report's unfilled triangle and integer costs `(2,3,5)`, the class representative
is `(1,1,1)` and selected mass is 10. After filling, the same cycle has computed
zero representative, zero mass and empty support. These are available values,
not `NoClass` or failure. Class queries reject noncycles. `minimum_class_mass`
is `Unavailable`; selected mass must not be renamed as the true minimum.

`Ready` means there is a validated usable P; the default solve is `Feasible`,
not optimal. Its initially unread objective is `NotComputed`; reading current
stretch gives exact rational 1 in this finite ring, without changing the solver's
certificate to `ExactOptimal`. Exact F2 topology can coexist with nonexact
floating masses. Floating `(0.1,0.2,0.3)` uses binary64 `fsum` with no manufactured
tolerance bound. Rational costs retain exact Fraction values.

The reference binds queries to all six identities including solver_run_id.
Two independent feasible runs have equal operator content but different run IDs;
mixing their queries is rejected by that reference. Content hashes omit
`source_metadata`: changing only dataset claims can preserve operator content.
It is not source authentication. Rust's future concrete consumer should require
the identities its operation needs; this does not mandate copying that strict
six-field Python policy into every Rust result.

### Diagram matching distances

Read the scalar from the existing `Result<f64>` only after checking success and
the operation selected. A single H1 `[-2,0)` versus computed-empty H1 gives
bottleneck 1, W1-L-infinity 1 and W2-Euclidean `sqrt(2)`. The returned scalar
does not store the operand identities, coverage, units, field or method. Retain
those facts with the calling application's record when needed; current code adds
no distance-record or matching-witness public API.

Unequal essential counts yield mathematical positive infinity. Numerical
overflow, cancellation and budget exhaustion are errors and produce no partial
scalar. Zero means a computed value under the chosen distance, not success for
an unsupported or uncomputed input. Truncated coverage is rejected even with
no intervals. Comparing approximation diagrams describes those supplied diagrams;
it does not establish a distance/error bound for unknown original data.

Fixed Persim source returns bottleneck 1 and W1-Euclidean `sqrt(2)` for this
finite fixture. With essential `(0,+infinity)` versus empty it warns, removes
the point and returns 0. That is a finite-point projection, not agreement with
Cocycle's full endpoint contract. Empty-input padding also creates a zero-cost
matching row for a synthetic diagonal point; exported indices require explicit
mapping/filtering. Extra columns must not silently become coordinate features.

## Availability, empty values and failures

Use each concrete API's existing states/errors. This table is vocabulary for
interpretation, not a new common enum or permission to coerce errors to zero.

| Observation | Reading rule |
| --- | --- |
| Successful but empty | PH computed-empty degree; explicit zero-size spectral domain; valid empty batch or support; requested empty representative list. Preserve that domain and the result's units/coverage. |
| Not requested/not computed | PH absent dimension is `DimensionNotComputed`; representatives `None` means no request. An unread reference kernel/stretch is `NotComputed`. Length zero cannot substitute for these distinctions. |
| Unsupported | Current reference minimum-class-mass is `Unavailable`; future unsupported coefficient/eigenproblem/action must report that specific limitation. Do not make a success record. |
| Failed or interrupted | Preserve invalid input, numerical failure, cancellation and resource failure with their diagnostics. Current distance errors have no partial value. A future solver may explicitly retain a validated candidate while still recording interrupted search. |
| Empty mathematical domain | Reference empty cycle-domain stretch is `EmptyDomain` under its explicit value convention; nonempty cycles with beta=0 have a computed stretch 0. A domain convention is not a failed calculation. |
| Infinite mathematical value | Distinct from every failure above. Current diagram distance permits +infinity for unequal essential counts; choose an explicit tagged encoding if an external format forbids nonfinite numbers. That encoding is not yet a Cocycle schema. |

Reference `QueryResult` accepts computed `0`, `False` and empty tuple; its missing
states cannot carry a computed zero. Absence has no inferred exactness. Not every
thrown Python error automatically creates a failed OperatorResult. A solver
without a feasible projection cannot be restored as a usable operator.

## Exactness, approximation and certificates

State the claim that was established rather than using one universal boolean:

- Algorithmically exact PH/distance means no algorithmic approximation on the
  supplied finite input. Filtration coordinates and matching arithmetic can still
  be binary64. PH exactness does not make a custom real lift valid or complete
  an omitted source skeleton.
- Input approximation is a property of construction/hypotheses. Preserve sparse
  Rips provenance, field, scale and coverage; no generic scalar comparison upgrades
  it into an original-data approximation certificate.
- Numerical spectral readouts retain dtype, solver/problem, selected/returned
  count, convergence, normalization, residuals and tolerances actually available.
  Small residual alone proves neither completeness nor orthogonality. An
  eigenvalue-error claim needs the relevant self-adjoint/generalized hypotheses
  and rounding analysis; a thresholded zero is not an exact nullity certificate.
- F2 projection legality requires idempotence, cycle image, annihilation of
  boundaries and preservation of cycle homology. The first three alone admit
  false zero projectors. Current objective, valid global bounds and independently
  replayable optimality evidence are separate. Feasible and exact coordinate
  arithmetic do not imply minimum stretch or minimum class mass.

Keep a certificate bound to its actual input, basis, costs, P, objective and
supported proof domain. Unknown or unsupported evidence stays unknown/unsupported.
Do not promote labels, JSON booleans or mere equality of reported bounds into
independent mathematical validation.

## Multi-input comparison checklist and error examples

| Check every required operand | Agreement rule / explicit failure |
| --- | --- |
| Available selected degree/domain | Selected PH H1 may compare `{1}` against `{1,3}`; requested H3 fails in both orders. Whole-domain consumers require equal domains or an explicit selected domain available on every operand. Do not silently intersect. |
| Mathematical meaning and coefficient field | PH result adapters require matching characteristics. F2 actions and real Hodge matrices cannot be treated as one projection family. Spectral eigenproblem, ordinary/normalized/generalized convention and selected coverage must agree for the chosen spectral comparison. |
| Units, scale and normalization | Check every operand, including empty ones. Unknown units do not equal dimensionless. Metres versus seconds, different coordinate-cost semantics and generalized versus ordinary eigenvalues fail unless an explicit justified conversion is made. Record that conversion and its changed context. |
| Source provenance | Inspect every source/snapshot and validation status. Different known datasets may compare; unknown origin cannot certify a requested original-data comparison. Source kinds/counts/content hashes do not authenticate dataset identity. |
| Coordinates when the operation needs them | Eigenvalue multiset comparisons need no equality of basis labels if each operator is valid in its own basis. Coordinate vectors, actions and merged query histories need the actual source/basis/Gram/P relationship, or a verified map/isometry. Shape equality or matching Betti number is insufficient. |
| Coverage/approximation and arithmetic | Complete distance inputs; compatible spectral selection; actual approximation hypotheses and numerical precision. Refuse censored-to-complete relabeling and missing convergence/error certificates. Comparisons can explicitly target approximate results without claiming original exact results. |
| Execution/readout identity | Reference query aggregation uses the same accepted P/costs/basis/run; scalar comparison of independent experiments can be legitimate without pretending they are one joint run. Preserve both runs and the operation's actual criterion. |

Current `_results` distance adapters enforce fields and declared edge-length
conventions, not this whole provenance/physical-unit checklist. Raw diagram calls
skip context enforcement entirely. This is an existing boundary, not a claim of
new enforcement. The Rust probe verifies both behaviors. The local spectral
example rejects changed degree, units, unknown source, operator definition and
partial coverage in both operand orders; its dictionaries/function are only a
finite contract demonstration, not an approved result schema or production gate.

## Borrow, move, clone and explicit conversion

| Action | Established behavior and cost |
| --- | --- |
| `diagram.dimension(k)?.iter()` | Borrowed logical view, yields Copy interval values; no full interval buffer is needed. Storage layout/address stability is not promised. |
| `AsRef<PersistenceData>` / `diagram()` / `context()` | Infallible borrow of already owned data; no cloning, sorting, reconstruction or extra mathematical guarantee. Checked owner pointers alias in the probe. |
| `PersistenceResult::into_parts()` | Moves common data and optional representatives; ordinals still address the moved diagram. No clone/sort. This retains the full established association. |
| `into_data()` / `into_diagram()` | The first discards representatives and keeps diagram/context. The second discards context and representatives; `PersistenceData::into_diagram()` also drops context. Lost provenance cannot be recovered from endpoints alone. |
| Explicit `clone()` or external array export | Clone copies owned buffers; export allocates/transforms according to its representation. Reordering/projection must remap ordinals/basis labels and keep a loss record; it cannot change mathematical claims. |
| Reference scalar/action readout | Applies the existing P and returns requested output without a full P export. Query records freeze/copy value/details, so this is not an unrestricted zero-copy promise. Kernel basis generation is explicit and cached, and can be large. |
| Reference snapshot/JSON/recovery | Snapshot validates; `to_dict`/JSON constructs owned wire data and explicit Matrix output can be quadratic. Compact handles have separate validated versions. Recovery revalidates P, supported certificates and official query histories; this has real cost and does not run a new solve. |
| Future spectral read | Borrow an already owned value/vector slice when provided; matrix conversion, eigen-solving, sorted pairing and optional eigenvectors have their own explicit allocation cost. No spectral Rust owner/accessor exists yet. |

Selecting H1 for a comparison does not change an owner to H1-only. Dropping
essential/censored points, rescaling, rounding, changing basis/weights or replacing
an action by its dense matrix is explicit conversion, with source meaning/loss
preserved. Logical borrowing never promises all internal algorithms allocate zero
bytes. No timing or RSS claim is inferred from these ownership observations.

## External imports cannot manufacture stronger guarantees

The existing Rust manual-diagram constructor can validate dimensions, finite
coordinates, endpoint/coverage consistency, sorting and multiplicity. It cannot
verify that a declared interval really came from the claimed data, field, complete
filtration or method. `PersistenceData` context construction is crate-private;
there is no general external builder that can mint library computation context.
Retain external facts alongside the imported diagram and record unknown fields.

A two-column array is insufficient to recover computed-empty dimensions,
coefficient field, physical units, source identity or whether `+infinity` meant
essential versus an encoded censoring cutoff. Import must get those declarations
from an actual external contract. NaN, invalid finite lifetimes and unidentified
endpoint encodings are rejected. Diagonal rows require explicit documented
handling; Cocycle's finite intervals require positive lifetime. Finite-only
Persim exports remain finite-only projections, including the warning/loss record.

Spectral arrays can validate shape, finiteness and pairing; checking residuals and
orthogonality additionally requires the actual operator, bases, Grams and numerical
policy. Missing vectors cannot provide eigenpair witnesses; a partial array cannot
provide full multiplicity/nullity. Projection import additionally needs A/D,
ordered bases/costs and independent homology-preserving legality checks. Neither
content hashes nor a serialized provenance field independently authenticate data.

The pinned OperatorResult JSON roundtrip preserves all fields, exact rational
encoding and the original run. Tampered betti readouts, duplicate JSON keys and
a `Feasible -> ExactOptimal` label change are rejected. Official readouts are
replayed; arbitrary custom query content without supported parameter binding is
not thereby mathematically certified. Supported recovery is not proof that every
external provenance claim or format is trustworthy.

## Result readers and migration decision on main a7ee9e9

This 2026-10-11 follow-up completes the result-representation research in
[#78](https://github.com/Aequiludium/cocycle-rs/issues/78) and the existing-result
review in [#79](https://github.com/Aequiludium/cocycle-rs/issues/79). Its fixed
production base is `a7ee9e911dcdeaea137b682bf383f7c2e1d090d5`, after #59-#63.
The current task covers existing code and future interface information.
Laplacian/eigensolver and F2 projection implementations remain future work;
the proposed reads below are not callable Rust APIs or computed spectral output.

The source review uses this base's `diagram/`, PH assembly, contextual distance
facades, descriptors, examples and standalone workers. It also rereads the pinned
Persim bottleneck/Wasserstein and TopoNetX spectrum functions and the pinned
homology-operator result/readout implementation identified above. The earlier
112-condition probe remains historical; it is not rerun or relabeled here.

### Four consumers and the information they require

| Output / concrete consumer | Minimum read and display | Authority, precision and retained identity |
| --- | --- | --- |
| Barcode / `betti_curve`, lifetime summary and diagram matching | Declared dimension membership, logical interval multiset, endpoint kind/value and coverage; borrow `PersistenceData` when field/source checks are required | Diagram owns dimensions/coverage; context owns field, source and construction/analysis ranges. Scales are finite binary64; declared scale convention is not physical-unit authentication. Manual diagrams have no inferred producer context. |
| Representative basis / cycle or cocycle inspector in `rips_representatives` | Request ordinal, interval ordinal, degree, field, query scale, kind and borrowed oriented simplex terms | Terms use increasing original vertex labels and exact prime-field residues. The query scale uses source units. Duplicate bars have separate local ordinals. Cycles/cocycles have different equations; neither is a shortest representative or a canonical vector across runs. |
| Spectrum / future value-only spectral summary or harmonic-subspace consumer | Operator/eigenproblem, degree, ordered coordinate basis, coefficient arithmetic, Grams/adjoints, units, requested/returned selection and values; paired vectors only when needed | Preserve actual method, source snapshot, convergence, residual/tolerance and normalization evidence. Full/partial selection and numerical zero thresholds are explicit. A residual alone proves neither completeness nor orthogonality. No production spectrum exists at this base. |
| Projection / future class-action and selected-mass consumer | Window/basis/field, selected P identity, source, coordinate costs/units, query availability/value and independently established feasibility/objective/certificate | F2 action is exact in its declared coordinates; floating costs can remain numerical. A valid zero action on a boundary is computed, not missing. Feasibility, current objective, global optimality and minimum-class mass are distinct. No production homology projector exists at this base. |

Reading owns mathematical interpretation. Plotting chooses labels/axes and
layout; feature extraction chooses a specified statistic and its domain. Neither
responsibility requires exposing reducer arrays or rerunning PH. Computation
may directly feed comparison or analysis without first exporting a display form.
The two implemented result families here are diagrams and representatives;
future spectrum/projection descriptions do not count as runnable consumers.

### Selected conversions and their losses

| Conversion / consumer | Semantics and cost | Lost information / decision |
| --- | --- | --- |
| `dimension(k)` and its iterator / one-degree reader | Borrow logical view; O(log(n+1)) lookup, O(k) traversal and O(1) auxiliary space; interval values are Copy | No full buffer copy. Retain the owning diagram/context for coverage/source. Selected degree alone does not describe the whole computed domain. Use existing API. |
| `AsRef<PersistenceData>` / contextual comparison | Borrow already stored data, O(1); no clone, sorting or context reconstruction | Does not add provenance or certificates. Keep existing data/wrapper boundary. |
| `into_parts()` / diagram plus witness consumer | Move existing data and optional witness vector; no clone or reordering | Keeps canonical interval ordinals. Borrow both moved outputs together when interpreting witnesses. Use existing API. |
| `into_data()` / context-aware consumer without witnesses | Move data and drop representatives | Witnesses and their associations are intentionally discarded. Use existing explicit operation. |
| `into_diagram()` / diagram-only statistic | Move diagram and drop context/witnesses | Field, source, approximation and request facts cannot be reconstructed from endpoints. Use only when those facts are not needed. |
| Finite birth/death array / external plot or feature tool | Explicit selected-degree export, O(k) owned entries plus caller-retained ordinal map | Dropping essential/censored points changes the represented result. Preserve endpoint tags/coverage separately; never convert censoring into a death. Defer a general exporter until an actual external consumer needs it. |
| Independent interval sorting / witness-aware external tool | An owned permutation/map costs O(k); sorting costs O(k log k) | Old interval ordinals require an explicit bijection. Deduplicating bars is not a bijection and loses witness identity. Do not add automatic reordering. |
| Spectral value slice / future value-only consumer | Borrow stored values; requested optional vectors do not force a full eigenbasis | Discarding vectors loses coordinate/subspace information; a partial spectrum cannot establish full nullity or a gap. Keep selection and problem metadata. |
| Eigenvectors to subspace projector / repeated-eigenvalue comparison | Explicit materialization in the aligned coordinate/Gram convention, generally quadratic output storage | Rotation/sign information is discarded deliberately; source basis and eigenspace selection remain necessary. Defer this helper until the concrete comparison is implemented. |
| Projection action versus full P / future class query or external matrix consumer | Prefer a requested action; dense n-by-n export allocates O(n squared) entries. Actual action cost belongs to the adopted representation | A bare matrix loses basis, costs, query states and certificate binding. Do not impose dense export, the Python snapshot schema or serialization on every Rust result. |

For example, two copies of H1 `[-2,0)` are two intervals, with independently
addressable witnesses. A set of pairs loses multiplicity. Encoding both
`Finite(1)` and `RightCensored { through: 1 }` as `(-2,1)` loses the difference
between death and last observation. Pinned Persim's removal of infinite deaths
is explicitly a finite-point projection, not a complete result conversion.

For the repeated eigenvalue of `3I`, two valid orthonormal bases may differ by
any rotation. Raw column equality is not an appropriate comparison. A dictionary
keyed only by eigenvalue can overwrite vectors, as in the pinned TopoNetX source.
Use paired values/vectors and, when required, a coordinate-aligned subspace
comparison. Ordinary `L v = lambda v` and generalized `L v = lambda G v` cannot
share an unexplained value-only label. Basis changes must also transform the
operator/Gram; coordinate-cost changes can change a projection objective even
when F2 homology is unchanged.

### Existing source disposition and caller coverage

Paths below are relative to `src/` unless otherwise stated. Each row retains a
single existing implementation; no new result container or compatibility copy
is introduced. The migration map in
[PR #103](https://github.com/Aequiludium/cocycle-rs/pull/103) permits a
source-supported no-op for these result boundaries.

| File / symbols | Actual producers and consumers / regression | Disposition and reason |
| --- | --- | --- |
| `diagram/mod.rs` exports and filtration compatibility reexports | All public callers and standalone workers | Retain paths. These are thin forwards, not duplicate ownership or an algorithm facade. |
| `diagram/persistence_diagram.rs`: `PersistenceDiagram`, `DiagramDimension`; `dimensions.rs`: `ComputedDimensions`; `interval.rs`: `PersistenceInterval`, `IntervalEnd` | PH assembly; distances/descriptors; `tests/contracts.rs`, `tests/descriptors.rs` | Retain logical reads and validation. Gaps, computed-empty, multiset order and censored endpoints already have public contracts. No physical slice/row API is added. |
| `diagram/data.rs`: `PersistenceData` | `PersistenceResult`, contextual distance adapters and compatible downstream wrappers; `tests/diagram_distances.rs` | Retain shared owned diagram/context and borrowed `AsRef`. One common implementation already exists. |
| `diagram/computation.rs`: `ComputationContext`, `PersistenceResult`, borrow/move/drop methods | `persistence/{source,filtered}.rs`, flag and Rips producers; `tests/{rips_api,diagram_distances,prime_fields}.rs` | Retain interpreted output and private constructors. Source lifetime and optional-witness distinctions remain explicit. Construction facts stay out of algebra. |
| `diagram/representative.rs`: `Representative`, `RepresentativeKind`, `RepresentativeTerm` | `persistence/simplicial/representatives/{basis,dual}.rs`; example and independent closure/duality tests in `tests/prime_fields.rs` | Retain result-local ordinals and owned oriented terms. Extend the existing duplicate-interval test to consume moved outputs after source release. |
| `persistence/mod.rs`: `assemble_diagram`, `RawIntervals`; `boundary/mod.rs` and Rips/simplicial call sites | Pairing kernels, reference tests and witness assembly | Keep assembly with computations; it interprets unpaired deaths using established coverage. Source-adapter relocation and any later assembly extraction belong to #72, not a second implementation in #79. |
| `diagram_distances/mod.rs`: raw/contextual facades and `points` | `tests/{diagram_distances,execution}.rs`, matching example and standalone distance worker | Retain current conversion and validation. `points` allocates matching buffers; only borrowing common data is allocation-free. Logical/prepared changes remain PR #26/#49. |
| `descriptors/{betti_curve,lifetime_statistics}.rs` | `diagram_analysis`, guide doctests, `tests/descriptors.rs` | Retain diagram-only consumers; no producer calls. Descriptor output/work allocations retain their own existing costs. |
| `examples/rips_representatives.rs`, representative guide | Public paths and owned-output consumer | Adapt the existing example to release source buffers, borrow data, move witnesses and call a consumer with only interpreted outputs. Preserve the original F3 square's hand-derived H1 `[1,2)`. |
| `tools/{diagram_dump,benchmark_driver}.rs`, `tools/reference/{rips_cocycle,sparse_cocycle}.rs`, `benches/{native,pipeline,distances,optimization}/cocycle.rs` | Logical interval iteration, public result/context access and explicit output protocols | Retain worker imports and formats. No stale result caller or dead private duplicate was found. Native comparisons retain their existing supported-case limits. |

Pending #104 at `e0e294a209d731db29589bfb61760c334568f021` moves source adapters;
it does not change diagram/result access. This follow-up uses main directly and
touches none of its production files. #100/#102 retain stage/modeling ownership.
No prepared-distance or adaptive-kernel candidate is adopted. Reverting this
follow-up restores only an example, test and documentation; no data migration
is required.

### Future reading signatures and explicit unavailable states

These are **interface sketches only**, to refine with #73/#75 before adopting
concrete names, storage or numeric algorithms. They add no trait, placeholder
method or universal status enum to the crate. A successful producer returns
owned interpreted output or a diagnostic error; it cannot return an empty
success to hide invalid input, unsupported requests or failed convergence.

```text
// Proposed spectrum reads, not existing cocycle symbols:
values(&self) -> &[f64]
vectors(&self) -> Option<borrowed paired vectors in the recorded ordered basis>
// Read degree, field, basis/Gram, problem, source, units, requested/returned
// selection and convergence/residual evidence from the same result.

// Proposed projection reads/actions, not existing cocycle symbols:
apply(&self, coordinates: &[u8]) -> Result<owned F2 coordinates>
selected_mass(&self, cycle: &[u8]) -> Result<value in recorded cost arithmetic>
// Bind action/readout to the same validated window, basis, P and costs.
// Read actual feasibility and objective/certificate evidence separately.
```

For a computed zero-dimensional operator, the spectral domain has size zero
and a complete spectrum can be empty. A requested empty partial selection still
records its selection; omitted vectors mean not requested, not a zero-rank
eigenspace. Missing residual estimates remain missing. A numerical eigenvalue
threshold does not certify exact homology. The solver determines supported
coefficient/eigenproblem cases and reports unsupported ones explicitly.

For F2 projection, empty previous/next bases retain the A/D shapes. A boundary
can have computed zero action and zero selected mass. Unread objectives,
unavailable minimum-class-mass queries and failed construction retain different
states in the concrete adopted interface. No class query accepts a noncycle.
`P squared=P`, `AP=0`, `PD=0` and homology preservation need independent evidence;
a feasible projection is not automatically optimal. Budget/cancellation errors
carry no result unless a future solver explicitly defines a separately validated
retained candidate. The actual producer owns that policy and its cost bounds.

### Migration acceptance and limits

[#87](https://github.com/Aequiludium/cocycle-rs/issues/87) G0/G1 and
[#94](https://github.com/Aequiludium/cocycle-rs/issues/94) A/B apply to the output,
ownership and downstream-consumption rows here. The fixed base's full debug
suite passed 216 tests/doctests with four ignored diagnostics before editing.
The follow-up changes no production source, public signature, runtime dependency,
dispatch, sorting, allocation site or reducer lifetime. The read and move bodies
are retained; the example performs one PH computation and uses the existing
borrowed/moved payload. This is a source-level cost audit, not a timing/RSS
measurement or a whole-pipeline zero-allocation claim.

Regression scope includes interval multiplicity, gaps/computed-empty,
finite/essential/censored coverage, fields/source context, independent cycle
closure/duality and moved witness ordinals. The runnable representative example
checks H1 Betti `[1,0]` at scales `[1,2]`, finite lifetime 1 and contextual
self-distance 0 after releasing the source. Relevant descriptor/distance tests
keep their own error and control semantics. New spectral/projection algorithms,
plotters, serialization and bindings remain unimplemented future work. Local
checks, hosted CI, merge and release remain separate evidence.

Pre-submission Windows validation passed the full debug/release and Rust 1.91
tests/checks (216 tests/doctests, four ignored diagnostics), Clippy, formatting,
strict rustdoc, ten examples, construction/critical-set example tests, eleven
Markdown doctest pages, source/docs/local-collaboration checks and the 183-file
package build. Pinned source hash/archive checks passed. The unchanged Python
tools ran 98 tests with two Windows platform failures (missing `resource` and
slash/backslash diagnostics) and five skips; Linux hosted tool/reference checks
remain separate. No new local native oracle or performance run is claimed for
this example/test/documentation change.

## Finite validation and acceptance

The standalone Rust probe was compiled against the fixed inventory-baseline library.
It checks logical reading/multiplicity/gaps, both operand orders, hand-derived
three-metric diagonal costs, endpoints/coverage, mathematical infinity versus
numeric/resource failures, common-data aliases, witness-preserving moves, absent
versus requested-empty representatives and raw versus contextual comparison.

The Python probe imports the unchanged pinned homology reference; it checks
typed query states, exact versus floating costs, unsupported minima, identities,
empty/dead classes, feasible versus objective/optimality, restoration/tampering
and invalid candidates. It also runs hand-derived spectral examples, the local
spectral checklist, two unchanged extracted TopoNetX functions and the two pinned
Persim modules using existing NumPy/SciPy/scikit-learn dependencies. No package
was installed and the installed Persim package is not used as the fixed source.

Final totals: 38 Rust conditions and 74 Python conditions, all passed. These are
112 finite conditions, including expected rejections and repeated operand-order
checks; they are not 112 independent datasets or general mathematical proofs.
Versions: Windows, Rust 1.98.1, Python 3.10.11, NumPy 2.2.6, SciPy 1.15.3 and
scikit-learn 1.7.2. The first Python attempt incorrectly expected a single
matching row for Persim's padded empty input; its failed log is retained, and
the corrected condition reads the actual point-to-diagonal row explicitly.

This adds no production mathematical change. No new all-feature, release/MSRV,
native GUDHI/Ripser/Topp, performance, full TopoNetX integration or hosted CI
acceptance is claimed by that probe. The inventory's all-feature/examples evidence
remains tied to its fixed main, rather than being relabeled as a result-contract run.

- Four result families have reading, units/domain and guarantee rules plus finite examples.
- Multi-input checks cover selected dimension, meaning, units, sources and relevant identities on every operand; current enforcement and caller obligations are separated.
- Exact arithmetic, algorithmic approximation, residuals, feasible candidates and optimality certificates have distinct meanings.
- Borrowing, moving, materialization and provenance loss are explicit; imports cannot improve the mathematical guarantee by changing format.
- Rejected/deferred designs have reasons: no universal container, global exactness flag, forced dense representation, synthetic source certificate or whole-domain intersection.

The #61 research decision is complete. #62/#63 subsequently reused the existing
PH borrow/move/error behavior with a runnable guide and a bounded compatibility
gate. Physical-unit/source policy and concrete spectral/projector result
signatures still require their actual operation before implementation. The
framework algorithms remain separate tasks.

The original probes and raw logs were retained as local research artifacts;
the 112-condition count is historical evidence, not a public CI suite bundled
with this report. The reproducible examples below expose the important inputs,
readouts and source-sensitive discrepancies. They do not claim a replay of all
112 original conditions or acceptance of complete external packages.

### Reproduce the reference action

Prepare the fixed homology source with the commands in
[input contracts](input-contracts.md#accepted-decision-and-reproducible-first-examples).
Save this block as `target/operator-contract-example.py` and run it with that
source on `PYTHONPATH`. It requires Python 3.10+ and no native extension.

```python
from homology_operator import (
    ChainWindow, HomologyOperator, Matrix, OperatorResult,
    ProjectionProblem, solve_projection,
)

a = Matrix.from_rows(((1, 1, 0), (1, 0, 1), (0, 1, 1)))
bases = (("triangle:v30", "triangle:v20", "triangle:v10"),
         ("triangle:e20,30", "triangle:e10,30", "triangle:e10,20"))
ring = ChainWindow(1, a, Matrix.zero(3, 0), *bases, (), (2, 3, 5),
                   arithmetic="ExactInteger", unit="cost")
filled = ChainWindow(1, a, Matrix.from_rows(((1,), (1,), (1,))),
                     *bases, ("triangle:f10,20,30",), (2, 3, 5),
                     arithmetic="ExactInteger", unit="cost")
op = HomologyOperator(ring, solve_projection(ProjectionProblem(ring)))
dead = HomologyOperator(filled, solve_projection(ProjectionProblem(filled)))
assert op.betti() == 1 and dead.betti() == 0
q = op.readout("selected_mass", (1, 1, 1))
assert q.state == "Computed" and q.value == 10 and q.identity == op.identity
assert dead.readout("selected_mass", (1, 1, 1)).value == 0
assert op.minimum_class_mass((1, 1, 1)).state == "Unavailable"
assert op.to_result().solver["certificate_level"] == "Feasible"
record = op.to_result()
assert OperatorResult.from_json(record.to_json()).to_dict() == record.to_dict()
promoted = record.to_dict()
promoted["solver"]["certificate_level"] = "ExactOptimal"
try:
    OperatorResult.from_dict(promoted)
except ValueError:
    pass
else:
    raise AssertionError("serialization upgraded a feasible certificate")
print("Ring/fill readouts, unavailable minimum and certificate rejection passed.")
```

### Reproduce the spectral and finite-point discrepancies

Use fresh ignored directories for fixed external sources:

```sh
git clone https://github.com/pyt-team/TopoNetX.git target/result-contract-toponetx
git -C target/result-contract-toponetx checkout --detach ce3d414b571fff95916c4cbc20b8d2c2c9d428e6
git clone https://github.com/scikit-tda/persim.git target/result-contract-persim
git -C target/result-contract-persim checkout --detach 096a25fe38644d31c6aabb880a24a05046c03559
python3 target/result-contract-example.py
```

Save the block below as `target/result-contract-example.py`, then run it from
the repository root. It requires NumPy, SciPy and scikit-learn in the chosen
environment; the original versions are listed above. It executes the two
unchanged AST-extracted TopoNetX functions and the fixed Persim modules directly,
without installing those packages. This isolates exactly the code inspected;
it does not test TopoNetX's classes, imports, decorators or package integration.

```python
import ast
import importlib.util
from pathlib import Path
import warnings

import numpy as np
import scipy as sp
import scipy.sparse as sparse
from scipy.sparse import diags

source = Path("target/result-contract-toponetx/toponetx/algorithms/spectrum.py")
tree = ast.parse(source.read_text(encoding="utf-8"))
names = {"hodge_laplacian_eigenvectors", "laplacian_spectrum"}
functions = [node for node in tree.body
             if isinstance(node, ast.FunctionDef) and node.name in names]
assert {node.name for node in functions} == names
namespace = {"np": np, "sp": sp, "sparse": sparse, "diags": diags}
exec(compile(ast.Module(body=functions, type_ignores=[]), str(source), "exec"), namespace)
small = sparse.csr_matrix(3 * np.eye(3))
values, vectors = namespace["hodge_laplacian_eigenvectors"](small, 1)
assert len(values) == 3 and np.linalg.matrix_rank(vectors) == 1
assert np.linalg.norm(small @ vectors - vectors * values) < 1e-12
large = sparse.csr_matrix(3 * np.eye(11))
values, _ = namespace["hodge_laplacian_eigenvectors"](large, 2)
assert np.allclose(values, 1)  # L v = lambda Diag(L) v.
assert np.allclose(namespace["laplacian_spectrum"](large), 3)  # Ordinary problem.

def load_module(name):
    path = Path("target/result-contract-persim/persim") / (name + ".py")
    spec = importlib.util.spec_from_file_location("fixed_" + name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module

bn = load_module("bottleneck").bottleneck
w1 = load_module("wasserstein").wasserstein
bar, empty = np.array([[-2., 0.]]), np.empty((0, 2))
assert bn(bar, empty) == 1
assert np.isclose(w1(bar, empty), np.sqrt(2))  # Euclidean ground cost.
with warnings.catch_warnings(record=True) as recorded:
    warnings.simplefilter("always")
    projected = bn(np.array([[0., np.inf]]), empty)
assert projected == 0 and recorded  # Essential point was removed with a warning.
print("Repeated vectors, generalized spectrum and finite-point projection reproduced.")
```

AST extraction needs an interpreter able to parse the fixed source. The original
TopoNetX class inspection required Python 3.12; only its spectrum functions are
parsed here. Preserve failures caused by a different Python/dependency version
as reproduction limitations. Cocycle's corresponding full-endpoint behavior is
covered by the Rust interface-guide doctests and existing distance tests; Persim's
finite-point result above cannot replace that acceptance.

[kernel]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/docs/design/kernel.md
[analysis]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/docs/development/diagram-analysis.md
[math]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/docs/reference/mathematics.md
[diagram]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/diagram/persistence_diagram.rs
[endpoints]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/diagram/interval.rs
[data]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/diagram/data.rs
[result]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/diagram/computation.rs
[context]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/filtration/context.rs
[distances]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/src/diagram_distances/mod.rs
[distance-tests]: https://github.com/Aequiludium/cocycle-rs/blob/9e6715f4c2e0118ab738486bd747ec3e361397de/tests/diagram_distances.rs
[operator-model]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/docs/RESULT_MODEL.md
[operator-result]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/src/homology_operator/result.py
[operator]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/src/homology_operator/operator.py
[solver]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/src/homology_operator/solver.py
[solver-contract]: https://github.com/proffitteoy/homology-operator/blob/a8d03d1699307ce1997b3c7e81a1946bb83f5408/docs/SOLVER_CONTRACT.md
[persim-bottleneck]: https://github.com/scikit-tda/persim/blob/096a25fe38644d31c6aabb880a24a05046c03559/persim/bottleneck.py
[persim-wasserstein]: https://github.com/scikit-tda/persim/blob/096a25fe38644d31c6aabb880a24a05046c03559/persim/wasserstein.py
[spectrum]: https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/toponetx/algorithms/spectrum.py
