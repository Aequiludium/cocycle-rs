# Roadmap

[Documentation](../README.md) / Design

This document records priorities, not a delivery schedule. Supported behavior is
listed in the [architecture](../development/architecture.md) and usage guides;
the [changelog](../../CHANGELOG.md) records versioned and Unreleased work.
Use the [version status](../README.md#version-and-capability-status) to
distinguish candidate implementation from releases and dated research.

## Direction

Build a native Rust TDA kernel with GUDHI's C++ capabilities as a reference and
specialized Rips computation where appropriate. Prioritize coherent input and
result contracts, composable analysis operations, ownership, and resource
behavior. Runnable Rust examples demonstrate the kernel; language bindings and
application frameworks are not first-stage deliverables.

Phase 4 extends this direction across data, filtration, computation, results,
comparison and analysis. Keep specialized persistence and concrete consumers;
filtration, diagrams and distances are optional when a route does not need them.
The six responsibilities do not require six crates or a runtime registry.

The [GUDHI C++ study](../research/gudhi-cpp.md) maps upstream capabilities and remaining
reading work. The [kernel design](kernel.md) defines responsibility boundaries
and distinguishes implemented changes from future extensions.
The [Rips subsystem design](rips.md) defines the selected complete Rips target,
including upstream Ripser, its API draft and acceptance gates. These are design
inputs, not claims of implemented parity or a request to copy either source layout.

## Current scope

The implemented core provides borrowed point clouds and matrix layouts, exact
threshold graph construction, supplied weighted flag filtrations, ordinary prime-field
persistence through dense/sparse access, explicit simplicial expansion and
incidence queries, requested cycle/cocycle bases, blocker-aware sparse approximation,
owned diagrams/context, basic descriptors and exact bottleneck/Wasserstein distances.
Cooperative work limits and cancellation are available on the richer compute
entry points. An independent explicit boundary implementation remains a test oracle; the
production representative path owns a separate reducer. See the
[construction guide](../guides/rips-construction.md) and the
[0.1.0 release scope](../../CHANGELOG.md#010).

The [filtered-complex boundary](../guides/filtered-complexes.md) is implemented:
validated supplied simplices, a four-method cell trait, signed-scale boundary
reduction and shared source context. This supplies an integration boundary for
future Alpha geometry; it does not implement triangulation or Alpha construction.

[Algorithm contribution paths](../development/algorithm-contributions.md) now
include descriptor, diagram-distance and explicit construction paths,
mathematical tests and focused checks exercised by CI. The lower-star constructor
is a teaching example, not an additional production API. The
[persistence algorithm walkthrough](../development/persistence-reduction.md)
now covers direct column work, implicit alternatives and default integration.

Borrowed stages and same-source inclusions, opt-in Big Steps critical sets and
source/interval assembly integration are implemented Unreleased additions. Use
the [stage guide](../guides/filtered-complexes.md#read-stages-and-their-inclusions),
[critical-set guide](../guides/critical-sets.md) and
[interface guide](../guides/interface-contracts.md) for calls and limits.
The [input](../research/input-contracts.md) and [result](../research/result-contracts.md)
reports and [route record](../research/route-integration.md) prepare future handoffs;
they do not implement spectral/Dirac, F2 projection or learning algorithms.

## Next priorities

<a id="phase-4-priorities"></a>

The [Phase-4 plan #57](https://github.com/Aequiludium/cocycle-rs/issues/57) sets the
next capability priorities. Complete reviewed integration and coherent documents,
then advance concrete families through narrow implementations and independent
acceptance. Interface preparation does not complete algorithm delivery.

1. **Data.** Extend an existing owner only for a concrete Alpha, cubical/cellular
   or generalized consumer, following the [modeling decisions](../research/data-modeling.md).
2. **Filtration.** Retain finite same-owner stages; research general chain maps and
   multiparameter transport separately under the [stage decisions](../research/filtration-stages.md).
3. **Computation.** Implement single-scale Laplacian/Dirac, true persistent spectral
   pairs and F2 homology operations in separate review units, following the
   [family inputs](../research/input-contracts.md#computation-family-preparation-on-main-a7ee9e9)
   and [output/resource contracts](../research/result-contracts.md#computation-family-outputs-and-resources).
4. **Results.** Add a concrete interpretable result owner with its actual producer;
   preserve ordered bases, availability and certificates under the [result decisions](../research/result-contracts.md).
5. **Comparison.** Define one useful non-diagram similarity, coordinate alignment,
   units and incompatible cases before adding an [API](../research/result-contracts.md#future-comparison-interface-information-only).
6. **Analysis.** Preserve existing descriptors; validate a minimal exporter and
   actual external consumer under the [handoff requirements](../research/analysis-learning.md#future-external-input-information-84-to-92).

After actual producers exist, the [route record](../research/route-integration.md)
guides their composition. Separate stage spectra are not a persistent Laplacian;
legal F2 projections are not automatically optimal or matrix-commuting. Feedback
remains application-owned and deferred until invalidation and derivative/tie
behavior have concrete evidence. The implemented API and kernel revisions below
provide the existing foundation for this work.

The [API design](rips-api.md) is implemented with a compatibility stage:
`RipsBuilder::build_complex` creates inspectable topology and `.persistence()`
configures direct analysis. Approximation, callback ownership, source context and
whole-operation execution controls follow the same contracts. Guides and native
workers use the new paths. Legacy names and functions remain available; removing
them requires a separately declared breaking release. Algorithm
optimizations and new filtration families remain separate review units.

The [aligned kernel design](kernel.md) owns the dependency, ownership, result and
extension decisions. Its [implementation sequence](kernel.md#implementation-sequence)
is included in the 0.1.0 release in four steps:

1. Reuse frozen simplicial incidence, cache structural facts,
   and limit direct exact point edge retention to the effective analysis range,
   without changing public result semantics.
2. Common result data and computed dimensions in the explicitly
   [declared result-API revision](kernel.md#result-api-migration), including analysis
   consumers and equivalent source-context assembly.
3. Independent owned-boundary computation, exact flag
   dispatch and direct simplicial calls for explicit/approximate Rips. Share
   compatible types and keep specialized public operations optional.
4. Contributor guides and examples reflect the current
   boundaries; persistence reduction has a focused check exercised by CI.

Diagrams support explicit computed-dimension sets;
default builders still compute every dimension through their requested maximum.
Current representatives remain simplicial bases. Common data supports cheap
`AsRef` borrowing. The design specifies static cell adaptation and the scope
of in-crate result construction. External result import, runtime plugins and
generic cellular witnesses are separate workstreams, not prerequisites for an
algorithm contribution. Mathematical optimizations and legacy-API removal also
remain separate changes.

Rips maintenance and evidence priorities remain within this broader plan:

1. **Complete the Rips subsystem.** Follow the
   [delivery sequence and exit gates](rips.md#delivery-sequence-and-exit-gates):
   exact sparse input through computation, explicit and higher-dimensional Rips,
   prime fields and requested representatives, then sparse approximation through
   computation. Stages 1-4 are implemented, including sparse approximation
   with hypotheses, blocker rules and provenance. Stage 5 has an [acceptance audit](rips-acceptance.md), resource
   regressions and native workflow measurements. Hosted CI for the submitted
   revision remains a release gate. Review new API contracts before extending scope.
   Each stage delivers a usable path with tests, examples and documentation;
   finishing one stage does not complete the whole target.
   The [implementation scope](rips-implementation.md) identifies exact directories,
   existing-file edits and review units for the first stage and later additions.
2. **Extend native evidence with each capability.** The
   [native suites](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/README.md) compare pinned GUDHI and upstream
   Ripser C++ workers. Exact and sparse correctness suites now cover construction,
   higher dimensions, fields and approximation; the pipeline records phase and
   end-to-end resource snapshots. Extend those fixtures with each new contract,
   following the [verification matrix](rips.md#verification-and-native-comparison)
   and [reporting rules](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/reporting.md). Stronger performance claims
   need a fresh comparative study with externally retained evidence.
3. **Retain focused H1 performance work.** The scoped sequence below can improve
   the existing path while broader capabilities arrive. Include nonmetric and
   bipartite fixtures alongside geometric inputs; performance claims need fresh
   native evidence.
   Review enumeration/reduction optimizations separately from API and storage
   refactors.

The earlier Rips-first plan remains the subsystem's design and evidence record.
Phase 4 broadens the next priorities without relaxing its ownership, resource,
reproducibility or result-interpretation requirements. See
[R1-R10](rips.md#required-capability-matrix).

Each crates.io release requires hosted CI for the exact release commit and
inspection of the packaged crate, following the
[release procedure](../../CONTRIBUTING.md#release-procedure). Release notes in the
[changelog](../../CHANGELOG.md) describe delivered capabilities; this roadmap
does not commit to dates for future work.

After reviewed code integration, a separate focused documentation consolidation
must align README, changelog, index, architecture and guides before release.
Run the fixed Markdown doctests and check the actual Cargo package's local links.
Keep #26/#49's paused preparation and #56's experimental adaptive H1 adoption gates
independent. Corrected matrix-H1 measurements do not supply #56's unmeasured
direct-graph/H2/H3 performance gates; documentation acceptance does not admit it.

Performance-only changes to existing paths must preserve pivot order, F2 parity,
multiplicity, and public coverage. New fields and generalized algorithms establish
their own documented invariants. Check each optimization against the reference and external
implementations; do not remove hard cases or use point count as a universal limit.

## H1 implementation sequence

**Historical subsystem sequence.** This table preserves the earlier Rips-focused
plan, not the Phase-4 queue. Recheck source and revision-bound reports before
resuming a step; the changelog owns subsequent production changes.

Keep maintenance refactors separate from changes to enumeration, reduction, or
storage. The current baseline has a single result-normalization path, named
transformation columns, and an isolated original-column shortcut search. Its
test-only reference and seven optimization configurations remain the correctness
checks for subsequent work.

| Order | Work | Evidence required |
| --- | --- | --- |
| 1 | Reuse cofacet enumeration when the original-column shortcut fails, including empty columns | Bipartite cutoff cases retain every censored interval; tie cases still fall back when the earliest pivot is owned; measure candidate scans as well as yielded cofacets |
| 2 | Compress pending working-heap entries by F2 parity | Preserve pivot order and odd multiplicities; measure circle heap peaks and the cost on small/easy inputs |
| 3 | Reduce repeated cofacet generation during transformation-column reconstruction | Nonmetric cases retain complete diagrams; measure regeneration counts together with added cache/storage costs |

For each step, run the reference and external comparisons, then compare the same
fixtures under the same precision and timing boundaries before and after the
change. Include uniform, circle, nonmetric, and full/truncated bipartite inputs.
Keep timeouts and regressions in the results. Candidate-vertex scans are currently
absent from the cofacet counter; zero yielded cofacets does not mean zero work.

Historical benchmark records describe the source at measurement time. A replay
against their saved diagrams checks agreement with those outputs; it is neither
a fresh external-library run nor a new performance baseline. Record new source
hashes and fresh measurements before making a speed or memory claim.

Dense computations and the legacy point entry point still retain a full distance
buffer. New threshold point construction and sparse computation avoid that buffer.
H0 edge sorting remains a separate scaling cost.
Sparse inputs and higher dimensions are now explicit Rips deliverables, with
their own API and validation gates rather than being folded into these H1
optimizations. Other mathematical capabilities remain separate work.

## Feature admission

A new capability needs a mathematical specification, clear input/output contracts,
an independent correctness check, and a reasonable ownership/error model.
Avoid speculative backend registries, unused traits, placeholder modules, and
application-specific research features. Language bindings remain separate from
the kernel's Rust API and release process.
