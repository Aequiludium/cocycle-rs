# Interface compatibility and framework entry report

[Documentation](../README.md) / Design

Technical report, 2026-10-09, connecting the full #59–#63 interface preparation
delivery. This page records the #62 reuse decision and #63 compatibility gate;
the linked research reports supply the inventory and complete input/result
reasoning. This closes interface preparation, not the later mathematical
frameworks or release.

| Task | Full report and evidence |
| --- | --- |
| [#59](https://github.com/Aequiludium/cocycle-rs/issues/59) | [Interface inventory](../research/interface-inventory.md): six frameworks, ownership costs, all eight historical PRs, publication delta and follow-up validation sources |
| [#60](https://github.com/Aequiludium/cocycle-rs/issues/60) | [Input contracts](../research/input-contracts.md): bases/orientation, coefficients, costs/Grams, identities, chain maps, finite examples and rejected alternatives |
| [#61](https://github.com/Aequiludium/cocycle-rs/issues/61) | [Result contracts](../research/result-contracts.md): four result families, availability, exactness/certificates, compatibility, moves/import losses and fixed-source counterexamples |
| [#62](https://github.com/Aequiludium/cocycle-rs/issues/62) | Existing-API reuse decision below and [executable guide](../guides/interface-contracts.md) |
| [#63](https://github.com/Aequiludium/cocycle-rs/issues/63) | Dated compatibility verification below; later consumers retain their own mathematics/acceptance |

## Result and scope

The first consumers can use existing borrowed geometry, frozen signed incidence,
`FilteredComplex`, ordinary persistence data and concrete descriptor/distance
functions. No demonstrated consumer requires changing a production signature.
The implemented addition is the [executable interface guide](../guides/interface-contracts.md):
stage/basis conversion, witness-preserving result moves and controlled composed
calls, continuously exercised by the existing CI quality job.

The decision for #62 is **reuse, with no production API change**. Its required
calling examples and conversion/control explanations are supplied. The decision
for #63 is **pass the existing-interface compatibility gate and proceed to
individual framework implementation**. This permits work on concrete consumers;
it does not claim a spectral solver, homology operator or all six frameworks
already exist. Pending algorithms keep their own mathematical acceptance gates.

No source trait registry, global identity service, common result container,
new dependency or shared public workspace is introduced. Current default builders,
specialized persistence dispatch and ordinary result guarantees are preserved.
Public behavior belongs to rustdoc and the [kernel design](kernel.md); this report
records the bounded decision and verification rather than redefining those APIs.

## Original compatibility baseline and historical pending work

The original production baseline is main
[`9e6715f4c2e0118ab738486bd747ec3e361397de`](https://github.com/Aequiludium/cocycle-rs/commit/9e6715f4c2e0118ab738486bd747ec3e361397de).
The original candidate
[`ec97ebba2d69619330b1fbc436a420480a01946a`](https://github.com/Aequiludium/cocycle-rs/commit/ec97ebba2d69619330b1fbc436a420480a01946a)
left `src/`, tests, examples, Cargo manifest/lock and native workers unchanged
from that commit. Its local checks and hosted runs belong to that exact revision.
For a later main, rerun affected guide checks instead of assuming this dated
decision certifies the newer implementation.

[#98](https://github.com/Aequiludium/cocycle-rs/pull/98) subsequently merged as
`03e896da106c3a90b8585e98a23bf28e013ac428`, after a synchronization commit.
The full #59–#61 reports are published from main
`cedf5965b774855614b761a1ad311d9dc7992c49`, which also contains #32 and #50–#53.
The [inventory publication delta](../research/interface-inventory.md#publication-delta-on-main-cedf596)
separates those changes from the original research. This docs supplement changes
no production implementation; its verification does not relabel the old
196-test or finite-probe counts as new-main results.

The then-open PR stacks were rechecked before the original compatibility change:

| PR | Observed head | Dependency and interpretation |
| --- | --- | --- |
| [#26](https://github.com/Aequiludium/cocycle-rs/pull/26) | `26866ac9223a7bd5cc7c35286e1090dc5759e509` | Logical-view preparation, based on main; not merged |
| [#49](https://github.com/Aequiludium/cocycle-rs/pull/49) | `5bb77b4b624580fe3cede4c8d538d06518de1a5d` | Scoped preparation, based on #26; not a current-main dependency |
| [#50](https://github.com/Aequiludium/cocycle-rs/pull/50) | `1bfd739dcf056ecc8f06254e0bbfc878df550bba` | Phase-3 foundation, based on main; not merged |
| [#51](https://github.com/Aequiludium/cocycle-rs/pull/51) | `c5441935250145d0a37874689a0dbe3f4bb13d0e` | H2 experiment, based on #50; production adoption is separate |
| [#52](https://github.com/Aequiludium/cocycle-rs/pull/52) | `231ffbed0fea0dc4467bd9412a9ddf1a2a10c5f9` | Workspace audit, based on #51; retained no-go decisions |
| [#53](https://github.com/Aequiludium/cocycle-rs/pull/53) | `c3fee147ca50673d2221d961397fa6c281056da2` | Research closeout, based on #52; no automatic merge authority |
| [#56](https://github.com/Aequiludium/cocycle-rs/pull/56) | `2dcd5e3d8a306f3cb4c4f23e486cdfc5add78db1` | Draft adaptive-kernel proposal, based on #53; failed performance gates remain |

At that snapshot, report-only [#32](https://github.com/Aequiludium/cocycle-rs/pull/32)
was separately based on main. Preserve this historical table rather than treating
it as today's open-PR list. The guide uses existing APIs without borrowing pending
preparation calls or describing experimental paths as production. It makes no
performance claim.

## Confirmed requirements and minimal action

| Consumer requirement | Existing owner / concrete action | #62 outcome |
| --- | --- | --- |
| Borrow input geometry or supplied topology | Geometry views, `WeightedGraph`, frozen `SimplicialComplex`, statically generic four-method `FilteredComplex` | Reuse; no universal source trait or runtime registry |
| Select a concrete stage and read signed boundaries | Borrow `simplices()`, filter `value <= scale`, retain degree bases and owner-local IDs, read `boundary(id)` | New runnable guide; explicit small dense export when requested, no stage/map public type |
| Compute ordinary persistence and request witnesses | Existing `.persistence()` builder, field/range options and representative requests | Reuse defaults and private specialized paths; no new engine selection |
| Read and move ordinary results | `diagram()`, `context()`, `dimension(k)`, `AsRef<PersistenceData>`, `into_parts()` | New runnable composition example; no clone/reconstruction required for common-data access |
| Compare and analyze results | Concrete raw/contextual distances and diagram descriptors | Reuse selected-dimension checks; document physical-unit/source obligations and explicit context loss |
| Control combined work | One private budget inside each terminal; same borrowed cancellation flag across calls | New runnable budget/cancellation example; separate terminals do not share consumed work |
| Read spectral or F2 projection results | Consumer-specific future owners with degree/bases/field/weights/solver evidence | Defer signatures to the actual algorithm; no placeholders or artificial PH dependency |

The unchanged areas are deliberate no-op decisions with executable evidence,
not undocumented missing implementation. General stage maps, spectral/Dirac
definitions, true minimum-class optimization and external learning remain
separate framework tasks. No consumer established a need for a cumulative public
budget or a metadata authentication API in this preparation task.

## Inputs, identity and conversion

Cell IDs are owner-local; ordered degree bases identify matrix coordinates.
Preserve the original signed orientation before any coefficient conversion.
In the guide, the triangle edge basis is `((20,30),(10,30),(10,20))` and its
vertex basis is `((30),(20),(10))`. At scale 1 the signed boundary is

```text
A = [[1,1,0],[-1,0,1],[0,-1,-1]], D has shape 3-by-0.
```

At scale 2 the face column is `(1,-1,1)^T`, giving integer `AD=0`.
An arbitrary F2 window `A=[1,1,0]`, `D=(1,1,0)^T` instead has `AD=0` only
modulo two: its naive real lift gives 2. Field-valid input cannot manufacture
missing real orientation or a real chain certificate.

Filtration values, positive F2 coordinate costs and spectral inner products
have separate semantics/units. A real/complex consumer needs its own adjoints
and Gram convention. A single-scale graph can be a one-dimensional complex
without a PH computation; it must not silently become a filled flag complex.
Ordinary stage spectra do not define a persistent Laplacian. A general chain
map needs identified endpoints/bases and both commuting boundary equations;
equal shapes or Betti numbers are insufficient.

Borrowing common data and source labels does not copy a complete representation.
Selected basis handles and a dense external matrix do allocate; the guide states
the cost and preserves zero-size shapes. A public `clone()` or external encoding
is explicit. Moving `into_parts()` retains interval/witness associations;
`into_data()` drops witnesses and `into_diagram()` additionally drops context.
Projection, sorting or rescaling must retain a loss/conversion record and remap
any associations they actually preserve. An allocation-free accessor does not
promise that every downstream algorithm allocates nothing.

## Result interpretation and compatibility

| Family | Required interpretation | Comparison/import boundary |
| --- | --- | --- |
| Persistence | Declared computed dimensions, multiplicity, finite/essential/censored endpoints and coverage | Validate each selected dimension on every operand. Complete versus censored cannot change through encoding. Manual constructors validate structure, not truthful origin. |
| Spectrum | Operator/eigenproblem, degree, basis/Gram, value units, selection, paired vectors and actual numerical evidence | Ordinary, normalized and generalized spectra differ. Partial spectra do not certify nullity; small residual does not establish a full orthonormal basis. |
| F2 projection/action | One legal P, exact coordinate field, ordered bases/costs, query availability and separate solver certification | Computed zero/dead class differs from absence. Selected mass is not shortest-class mass; feasible is not optimal. Content hashes do not authenticate dataset claims. |
| Distance/descriptor | Selected degree, definition/metric/order, units/normalization, endpoint participation and arithmetic | Scalar values retain no input identities by themselves. Mathematical infinity differs from numerical/resource failure. Approximate inputs do not imply an original-data error bound. |

Computed-empty H3 in a diagram domain `{1,3}` is valid; H0/H2/H4 are missing.
An H1 comparison can use unequal full domains, while an H3 comparison must
fail if either operand lacks it. Whole-domain consumers cannot silently
intersect domains. Current contextual distance adapters check fields and declared
edge-length conventions, not physical units or authentic dataset snapshots.
Raw diagram calls leave those context checks to the caller.

External arrays cannot recover uncomputed/empty degrees, source completeness,
coefficient field or physical units from endpoints alone. Preserve unknown
metadata. A missing result, unsupported query or failed solve must not become
computed zero. Floating residuals, exact F2 coordinates, algorithmic exactness,
valid feasible projections and independently verified optimality are different
claims. A serialized label cannot upgrade them.

## Six framework entry checks

These are interface entry capabilities, not claims that every proposed algorithm
is implemented. The runnable guide covers the shared handoff; existing guides
and regressions cover the production operations it uses.

| Framework | First usable entry at this gate | Still owned by later work |
| --- | --- | --- |
| Data/modeling | Borrowed points or supplied graph/complex; explicit signed incidence and basis keys | New topology families, sheaf/path/hypergraph contracts |
| Filtration/stages | Existing filtration builders plus a borrowed, inclusive stage selection | General stage-view/map APIs and persistent spectral definitions |
| Computation | Native ordinary persistence; explicit A/D conversion for a concrete external consumer | Laplacian, Dirac, F2 operator and generalized-topology engines |
| Result access | Borrowed common persistence data, logical degree views and witness-preserving moves | Concrete spectral/projector owners with their actual evidence |
| Distances/similarity | Existing contextual diagram distance with explicit scalar meaning | New mathematically defined comparisons for other result families |
| Analysis/optional learning | Betti/lifetime readouts without re-running PH; caller can consume those values directly | Optional application-specific learning/feedback and derivative checks |

Proceed to concrete framework work using these confirmed entries. A later
consumer that demonstrates a missing operation should add the smallest justified
function in its existing owner and rerun the affected mathematical/compatibility
checks. This gate does not authorize speculative crates, abstract backends or
the merger of pending experimental implementations.

## Verification and evidence limits

Original local acceptance used Windows, Rust 1.98.1 and Python 3.10.11. The code/manifest
identity comparison against the fixed main confirms no production/runtime
resource change, new dependency or retained buffer. There is no new performance
measurement: byte-identical production source is an identity observation, not
a benchmark or a claim about total documentation-build cost.

The acceptance checks run the guide's three independent Rust doctests on both
the baseline and candidate, including signed basis/shape extraction, result
associations, both field-mismatch orders, per-terminal budget retries,
pre-cancellation/recovery, infinity and numerical failure. They are hand-derived
finite examples, not independent large native oracle comparisons.

Candidate regression validation passed 196 all-feature tests/doctests (two
profiling tests ignored), all nine existing examples and the construction
example's four tests. Source, Markdown, artifact and diff checks cover the
submitted files. All three new guide doctests passed against the baseline,
candidate and Rust 1.91, preserving the existing minimum toolchain. Full public API migration is
unnecessary because production signatures and semantics do not change.
Commands, raw counts, exact candidate SHA and hosted check status belong to the
PR evidence. Hosted CI, native reference jobs and merge/release status must be
reported separately; configuring the guide in CI is not a hosted pass.

The unchanged Python tool suite passed all 82 tests on Ubuntu WSL. Its Windows
run retained two existing failures: `benchmark_rips_pipeline.worker` imports
the POSIX-only `resource` module, and one documentation test expects forward
slashes in a Windows-path diagnostic. That run had one error, one failed test
and three skips; it is not a Windows tool-suite pass. No unrelated tooling fix
was adopted. At original validation, the primary checkout's local collaboration
checker retained five historical missing links to then-pending #32 reports.
That older checkout's local links were outside the submitted documentation;
#32 has since merged on main as recorded in the inventory publication delta.

Earlier A1/A2 local research used fixed
[homology-operator](https://github.com/proffitteoy/homology-operator/tree/a8d03d1699307ce1997b3c7e81a1946bb83f5408),
[Persim](https://github.com/scikit-tda/persim/tree/096a25fe38644d31c6aabb880a24a05046c03559)
and [TopoNetX](https://github.com/pyt-team/TopoNetX/tree/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6).
Those reference probes are historical research evidence, not acceptance of new
Rust mathematics. Persim's studied distance functions filter infinite deaths;
its W1 uses Euclidean cost. The studied TopoNetX eigenvector function keys vectors
by rounded eigenvalues and switches to a generalized diagonal problem for larger
matrices. Do not adopt their bare outputs without the actual endpoint/problem
contract. No external package becomes a Rust dependency through this report.

No local native GUDHI/Ripser/Topp rerun is required for this documentation-only
production diff; no new native agreement is claimed. Any future change to
shared mathematical conversions or algorithms must satisfy the applicable
[validation obligations](../../CONTRIBUTING.md#verification). Existing historical
performance failures, unsupported reference cases and unmerged results remain
outside this gate, with their original evidence boundaries preserved.
