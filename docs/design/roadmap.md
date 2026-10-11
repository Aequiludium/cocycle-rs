# Roadmap

[Documentation](../README.md) / Design

This page owns future priorities, not a delivery schedule or API specification.
The [version status](../README.md#version-and-capability-status) separates releases,
Unreleased implementation and dated evidence. The
[architecture overview](../development/architecture.md#six-responsibilities-and-data-routes)
owns current capabilities; the [changelog](../../CHANGELOG.md) owns version notes.

## Direction

Build a native Rust TDA kernel around six usable responsibilities: data,
filtration, computation, results, comparison and analysis. Keep specialized
ordinary persistence paths and concrete consumers. Filtration, diagrams and
distances are optional when a consumer does not need them. Six responsibilities
do not imply six crates, empty framework modules or a universal backend registry.

The earlier Rips subsystem plan remains a bounded subsystem design and evidence
source. Phase 4 broadens the next priorities beyond Rips; it does not abandon
existing Rips correctness, compatibility or resource obligations.

## Current foundation

Borrowed geometry, checked graphs and signed simplicial/cell sources already
support ordinary prime-field persistence, concrete result reading, complete-diagram
distances and descriptors. The integrated candidate adds borrowed stages and
same-source inclusions and retains opt-in Big Steps critical sets. Use the
[stage guide](../guides/filtered-complexes.md#read-stages-and-their-inclusions),
[interface guide](../guides/interface-contracts.md) and
[critical-set guide](../guides/critical-sets.md) for callable behavior and limits.
Critical-set proposals are distinct from F2 homology projections.

The [interface preparation gate](interface-compatibility.md) and
[route record](../research/route-integration.md) retain their exact source and
validation scope. They deliver current calls, contracts and finite examples;
they do not implement Laplacian, Dirac, persistent spectra, F2 projection,
generalized topology, non-diagram similarity or external learning algorithms.

## Phase 4 priorities

The public [Phase-4 plan #57](https://github.com/Aequiludium/cocycle-rs/issues/57)
organizes work by capability and evidence. Complete coherent documentation and
reviewed integration first. Then advance each concrete family through its own
mathematical specification, minimal implementation and independent acceptance.
Interface preparation is a prerequisite; it does not close algorithm tasks.

| Responsibility / plan | Next bounded work | Decision owner and exit evidence |
| --- | --- | --- |
| Data, [#64](https://github.com/Aequiludium/cocycle-rs/issues/64) | Extend an existing owner only when a concrete Alpha, cubical/cellular or generalized consumer needs it. | [Data decisions](../research/data-modeling.md): validate topology, orientation, source identity and conversion costs independently. |
| Filtration, [#67](https://github.com/Aequiludium/cocycle-rs/issues/67) | Retain finite same-owner stages; research general chain maps and multiparameter semantics separately. | [Stage decisions](../research/filtration-stages.md): basis association, chain-map squares, composition and coverage; no inference from matching IDs across owners. |
| Computation, [#70](https://github.com/Aequiludium/cocycle-rs/issues/70) | Implement concrete single-scale Laplacian/Dirac, true persistent spectral pairs and F2 homology operations in separate review units; generalized models need their own theory. | [Input proposals](../research/input-contracts.md#computation-family-preparation-on-main-a7ee9e9) and [output/resource contracts](../research/result-contracts.md#computation-family-outputs-and-resources): independent matrices/actions, complete or explicitly partial solver evidence and measured resources. |
| Results, [#77](https://github.com/Aequiludium/cocycle-rs/issues/77) | Add an interpretable result owner only with its actual producer. | [Result decisions](../research/result-contracts.md): ordered bases, availability, multiplicity, residuals/certificates, moves and explicit provenance loss. |
| Comparison, [#80](https://github.com/Aequiludium/cocycle-rs/issues/80) | Define one useful non-diagram similarity before adding an API. | [Comparison proposals](../research/result-contracts.md#future-comparison-interface-information-only): mathematical meaning, coordinates, units, incompatible/undefined cases and independent expected values. |
| Analysis, [#83](https://github.com/Aequiludium/cocycle-rs/issues/83) | Preserve concrete descriptors; implement a minimal exporter/external consumer only for a validated use case. | [Handoff information](../research/analysis-learning.md#future-external-input-information-84-to-92): sample/stage/basis keys, shapes, units, conversion cost, roundtrip and an actual external run. |

After actual producers exist, [route work #86](https://github.com/Aequiludium/cocycle-rs/issues/86)
can test single-scale spectral/Dirac consumption, persistent spectral pairs,
F2 class transport, direct comparisons and external learning. The
[route report](../research/route-integration.md) owns these handoffs. Separate
stage spectra are not a persistent Laplacian; legal F2 projections are not
automatically optimal or matrix-commuting. Feedback remains application-owned
and deferred until invalidation and derivative/tie behavior have concrete evidence.

## Integration and adoption gates

Prepare one documentation consolidation after the #113-#119 code candidate is
reviewed and integrated. Its README, changelog, index, architecture, guides and
fixed doctest list must describe the same version before a release. Check the
actual Cargo package's local links as well as the checkout. Follow the
[release procedure](../../CONTRIBUTING.md#release-procedure); a passed historical
gate or open PR is not a published capability.

Keep [#26](https://github.com/Aequiludium/cocycle-rs/pull/26) and
[#49](https://github.com/Aequiludium/cocycle-rs/pull/49) distance-preparation
admission separate and paused until their own review and current-head performance
evidence. [#56](https://github.com/Aequiludium/cocycle-rs/pull/56) remains an
experimental adaptive H1 candidate outside this integration. Its corrected
matrix-H1 measurements do not supply the missing direct-graph/H2/H3 performance
gates. Documentation consolidation does not change any adoption decision.

## Rips maintenance and evidence

The [Rips subsystem design](rips.md), [API design](rips-api.md),
[implementation scope](rips-implementation.md) and [acceptance audit](rips-acceptance.md)
retain the R1-R10 requirements, source history and subsystem validation.
The [kernel design](kernel.md) owns the implemented dependency/result/extension
decisions, and the [GUDHI study](../research/gudhi-cpp.md) is fixed upstream
reading evidence rather than a feature-parity claim.

Review enumeration/reduction optimization separately from public API and storage
changes. Preserve pivot order, F2 parity, multiplicity, coverage and prime-field
semantics. Extend independent and pinned native comparisons for each capability.
Use the [reporting rules](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/reporting.md) for same-fixture, same-source,
same-boundary measurements; retain regressions, timeouts and unsupported cases.
Historical timings and saved-diagram replays do not measure a newer head.

## H1 implementation sequence

**Historical subsystem sequence.** The table below preserves the earlier
Rips-focused plan. It is not the Phase-4 queue or a current delivery status.
Subsequent production changes and experiment dispositions belong to the
changelog and revision-bound reports; recheck source before resuming a step.

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

A new capability needs a mathematical specification, narrow input/output
contracts, a runnable implementation, an independent correctness check and an
ownership/error/resource model. Start in an existing owner where possible.
Avoid speculative registries, unused traits, placeholder modules and application
dependencies in the kernel. Research or a justified no-go can complete a study
without delivering the proposed algorithm. Merge and release remain separate.
