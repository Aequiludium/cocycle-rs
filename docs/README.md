# Documentation

Cocycle is a native Rust TDA library. It computes dimension-generic ordinary Rips
persistence over prime fields, constructs exact threshold graphs and frozen
simplicial complexes, and accepts supplied flag filtrations. Sparse Rips
approximation includes blocker-aware computation and explicit metric hypotheses.
Owned results preserve coverage, field and approximation context; cycle and
cocycle bases are available on request at specified scales. Descriptors derive
lifetime statistics and Betti curves from diagrams.

Supplied simplicial and filtered-cell complexes support boundary reduction with
signed scales and unequal vertex births. Contributor walkthroughs cover diagram
analysis and explicit simplicial construction.
Start with a task below; design pages distinguish implemented decisions from
proposed APIs and remaining acceptance gates.

<a id="version-and-capability-status"></a>

Usage pages describe this checkout, including **Unreleased** borrowed stages,
critical sets, feasible F2 operators and source/interval assembly integration. The README dependency
selects [0.1.1](../CHANGELOG.md#011); see the [Unreleased notes](../CHANGELOG.md#unreleased)
and [architecture](development/architecture.md) for candidate additions.
Standard real single-scale Hodge action, matrices and bounded dense eigenpairs
are available in this checkout. Dirac solvers, general SPD/complex Hodge,
persistent spectra, optimal F2 projections, non-diagram
similarity, external learning and automatic feedback remain future contracts.
The [roadmap](design/roadmap.md#phase-4-priorities) sets priorities; input/result
reports retain the proposals.

This revision is prepared from [#119 at
`a2a1841c`](https://github.com/Aequiludium/cocycle-rs/commit/a2a1841c3d566758497f377da5a72f422b8728be),
observed open on 2026-10-11, including the later accessor-contract clarification.
The original audit used `64ab5fdc`. The [route record](research/route-integration.md#integrated-existing-pipeline-87-and-94)
lists the migration/source/contract inputs. Those hashes, PR states and checks
remain historical evidence after integration; they are not a live release status.
Main integration and release follow the [release procedure](../CONTRIBUTING.md#release-procedure).

The candidate excludes #26/#49's paused preparation and #56's experimental adaptive
H1 dispatch. Their adoption gates remain separate. Corrected #56 matrix-H1 timings
do not supply its unmeasured direct-graph and H2/H3 performance gates.

## Use the library

- [F2 homology operators](guides/homology-operators.md): owned chain windows,
  selected class costs and controlled finite transport from borrowed stages.

- [Real Hodge Laplacians](guides/hodge-laplacian.md): sparse action, explicit
  matrices, small dense spectra, numeric evidence and independent reproduction.
- [Big-steps critical sets](guides/critical-sets.md): sparse F2 endpoint proposals,
  bounded partial U and a lower-star simplification step.
- [Borrow inputs and compose results](guides/interface-contracts.md): signed stage
  conversion, result ownership and per-operation execution controls.
- [Supplied filtered complexes](guides/filtered-complexes.md): validated simplicial
  storage, borrowed stages and inclusions, signed scales and boundary reduction.
- [Rips guide](guides/rips.md): choose an input, compute persistence, interpret
  coverage and endpoints, and derive measurements.
- [Rips construction guide](guides/rips-construction.md): matrix layouts, custom
  distances, threshold graphs, sparse computation and execution controls.
- [Sparse Rips approximation](guides/sparse-rips.md): metric hypotheses, blockers,
  sampling provenance and implicit/explicit computation.
- [Fields and representatives](guides/rips-representatives.md): select a prime
  field and request owned cycle/cocycle bases associated with intervals.
- [Diagram-distance example](../examples/diagram_distances.rs): compare complete
  diagrams with bottleneck, W1 or W2; run `cargo run --locked --example diagram_distances`.
  The [matching specification](reference/mathematics.md#16-diagram-matching-distances)
  defines the metrics, supported endpoints and context requirements.
- [Runnable square example](../examples/square.rs): run
  `cargo run --locked --example square` from the repository root.
- [Published API reference](https://docs.rs/cocycle): browse rustdoc for signatures
  and error contracts. For the current checkout's [API source](../src/lib.rs),
  build locally with `cargo doc --locked --no-deps --open`.

## Understand and develop the core

| Question | Read |
| --- | --- |
| How do I contribute an algorithm in my research area? | [Algorithm contribution paths](development/algorithm-contributions.md) |
| How do I implement and test a diagram descriptor? | [Diagram analysis walkthrough](development/diagram-analysis.md) |
| How do I contribute diagram matching? | [Distance contribution path](development/diagram-analysis.md#contribute-diagram-distances) |
| How do I construct and validate explicit simplicial filtrations? | [Complex construction walkthrough](development/complex-construction.md) |
| How do I implement a persistence or reduction algorithm? | [Persistence algorithm walkthrough](development/persistence-reduction.md) |
| What mathematical and numerical guarantees apply? | [Mathematical specification](reference/mathematics.md) |
| Which papers justify the definitions and algorithms? | [Bibliography](reference/bibliography.md) |
| Where does code belong and how do the parts interact? | [Current architecture](development/architecture.md) |
| How should modules, APIs, errors, and source files be written? | [Code conventions](development/conventions.md) |
| How do we check correctness independently? | [Testing and validation](development/testing.md) |
| Which commands and review rules apply? | [Contributing](../CONTRIBUTING.md) |
| How are versions chosen and published? | [Release procedure](../CONTRIBUTING.md#release-procedure) |

For algorithm work, start with the contribution path and relevant mathematical
section. Maintainers use the architecture and validation obligations for kernel
integration. The specification keeps shared notation
and numbered derivations in one place; rustdoc owns individual API contracts.

## Plan and compare

| Document | Status and purpose |
| --- | --- |
| [Roadmap](design/roadmap.md) | Selected priorities; not a delivery schedule |
| [Rips API design](design/rips-api.md) | Implemented workflow design, execution boundaries and legacy compatibility/migration |
| [Complete Rips subsystem](design/rips.md) | Rips target, original GUDHI/Ripser comparison, API sketches and current acceptance matrix; use guides/rustdoc for callable APIs |
| [Rips acceptance audit](design/rips-acceptance.md) | R1-R10 evidence, resource boundaries and local/hosted validation distinction |
| [Rips implementation scope](design/rips-implementation.md) | Stage-specific new directories, source moves, code/tooling changes and review units |
| [Kernel design](design/kernel.md) | Implemented adapter, result and algorithm-boundary revisions; ownership, conversion costs, extension policy, API migration and primary design references |
| [Interface preparation reports](design/interface-compatibility.md) | Complete #59–#63 delivery chain, existing-interface reuse, runnable framework entries and dated compatibility acceptance |
| [Route integration](research/route-integration.md) | Same-version legacy acceptance and current handoffs versus future running goals for #86 |
| [Interface inventory](research/interface-inventory.md) | Historical source audit, borrow/move/copy costs and fixed PR stacks; use the architecture for current capabilities |
| [Shared input contracts](research/input-contracts.md) | Ordered bases, signed incidence, fields, weights/Grams, source identity, chain maps and reproducible finite examples |
| [Data modeling](research/data-modeling.md) | Three first objects, stable borrowed reads, two-consumer entry and conditional Alpha/cubical/cellular expansion order |
| [Result interpretation and comparison](research/result-contracts.md) | Four result families, consumer/conversion table, existing result migration, proposed spectral/projection reads and fixed-source discrepancies |
| [Descriptor contracts and external handoff](research/analysis-learning.md) | Existing descriptor disposition, consumer regressions and future sample/basis/scale information |
| [GUDHI C++ study](research/gudhi-cpp.md) | Pinned upstream source map and reading plan; excludes Python wrappers |
| [Filtration stages and maps](research/filtration-stages.md) | Fixed-source scale/order comparison, finite single-parameter stage decision, checked inclusions and deferred general/multi-parameter maps |
| [H1 optimization transfer](research/h1-fast-to-generic.md) | Source-linked F2 transfer matrix, complete clearing invariants and safe H2 continuation |
| [Big Steps algorithms](research/big-steps.md) | Paper framework, all three main-text pseudocodes, bounded Partial U and production scope |
| [Benchmarks](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/README.md) | Suite-specific protocols and maintained comparisons with measured revisions and evidence status |
| [Rips comparison](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/reports/rips-comparison.md) | Measured native correctness and performance; updated in place on reruns |
| [Computation-family preparation](research/input-contracts.md#computation-family-preparation-on-main-a7ee9e9) | Proposed Laplacian/Dirac/F2 inputs, generalized-model decisions and [output/resource contracts](research/result-contracts.md#computation-family-outputs-and-resources); algorithms remain future work |
| [Performance reporting rules](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/reporting.md) | Comparability, sampling, memory, source retention and report template |
| [Changelog](../CHANGELOG.md) | Unreleased implemented changes and versioned release notes |

For future capabilities, use the GUDHI study as evidence, the kernel design to
evaluate boundaries, and the roadmap to select the next operation. A source map
does not establish feature parity, and historical measurements do not describe
the performance of unmeasured code.

## Organization and maintenance

| Directory | Owns |
| --- | --- |
| `guides/` | Task-oriented usage of implemented capabilities |
| `reference/` | Mathematical contracts, derivations, and bibliography |
| `development/` | Current architecture, coding conventions, and validation strategy |
| `design/` | Design decisions, implementation records, future direction and acceptance gates |
| `research/` | Upstream investigations tied to a source revision and scope |

Add a document to the directory matching its purpose and link it from this index.
Use descriptive capability names such as `guides/rips.md`; add subdirectories
when a subject has multiple documents that need their own navigation. Do not
create empty sections for unimplemented capabilities or duplicate the source tree
outside the architecture document.

Keep current behavior, proposals, and dated observations explicit. Link to the
document that owns a fact instead of maintaining parallel versions. Project
documentation is in English. Paths in prose and shell commands are relative to
the repository root unless stated otherwise; Markdown links resolve relative to
their document. See [documentation maintenance](../CONTRIBUTING.md#documentation-ownership)
for checks to run when adding or moving pages.

Use this index for version status, the architecture for current owners/routes,
rustdoc and guides for callable behavior, the roadmap for future priorities,
input/result reports for proposed contracts, and the route report for composed
handoffs and dated acceptance. The inventory retains the original audit rather
than a second current capability table. Keep historical headings, hashes and
measurements intact; prepend a scope label and link to the current owner instead
of rewriting old evidence.

Guides and reference pages are included in the Cargo package. Benchmark reports,
controllers and repository CI configuration remain repository resources; links
to those excluded files use GitHub URLs. CI checks local targets and headings in
the actual package as well as in the checkout; external URLs are not fetched.
