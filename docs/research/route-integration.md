# Route integration and compatibility acceptance

[Documentation](../README.md) / Research / [#86](https://github.com/Aequiludium/cocycle-rs/issues/86)

This 2026-10-11 record separates existing pipeline acceptance from future
algorithm integration. The current issue scopes require interface preparation
for #88-#93, while keeping their later running goals open. Contract examples,
hand calculations and ordinary PH are not executions of future algorithms.

## Integrated existing pipeline: #87 and #94

The frozen pre-migration baseline for closed #87 is main
`a7ee9e911dcdeaea137b682bf383f7c2e1d090d5`. The integrated input commit is
`177790cdb123be07e2192f68700b45077d812e99`, on a new feature branch. It contains
the following original PR commits through ordinary merges; their branches and
review gates remain intact. This is a reviewable integrated candidate, not an
assertion that those PRs have merged into main.

| PR / issue | Fixed original head | Actual disposition in this version |
| --- | --- | --- |
| #104 / #72 | `ade7ff83b7e043a8eed072cc6d6e9653aaca0810` | Source-owned Rips/simplicial adapters and one ordinary interval assembly owner; unchanged sealed public paths and mathematical kernels. |
| #105 / #79 | `f7bdf0b4ada8dccb55b3db3a02a2a32fee67edf3` | Retain the established result/read owners; exercise witness-preserving moves and reads after releasing sources. |
| #107 / #82 | `466609721ffd1acde4a792b8df22e4ed54ae95f2` | Retain concrete diagram preparation/solvers and callers; preserve both result and distance disposition tables. |
| #106 / #85 | `0f1e270fdc59e0d5128b85d5c129c88b3f104ea3` | Retain descriptors; include cutoff/multiplicity regression and existing consumer examples. |
| #100 / #69 | `991e867f6897aaa79163c283a2e23551d066ed6a` | Same-owner borrowed stages/inclusions, with the two stage hooks routed to #104's relocated adapters. |
| #102 / #65 | `0cff0f0bbd2daf58d7391f31fe17bef5fdfc6963` | Source/basis construction example; no new topology container. |
| #103 / #71 | `913b4d7c438a113788bdd9694f08a6e51612225a` | Retain the migration map and current-versus-future scope correction. |
| #112 / #73-#76 stage A | `cc1cddeb8175032a4ccca9ae06788dadb54a627a` | Retain proposed family contracts; no spectral, Dirac, projection or generalized solver. |

Shared index and report conflicts preserve all relevant sections. The integrated
version does not adopt #26/#49's pending preparation API, #56's experimental
dispatch, non-diagram comparisons or a learning runtime. Final file contents,
rather than superseded intermediate commits, define the candidate's scope.

### Identity, errors and consumers

`tests/filtered_complex.rs` includes
`staged_results_keep_witnesses_and_consumers_after_source_release`: a signed
triangle with gapped vertex labels, births -2, edges at 0 and filling at 2.
Same-owner inclusions preserve integer incidence. Over F3, analyzing the original
source through stage 0 yields censored H1; explicitly importing that stage as a
supplied ring yields an essential class. Comparing the censored diagram fails.
At stage 2 the owned result has H1 [0,2), and both cycle/cocycle witnesses still
address its ordinal after every input owner and request is released. The cycle's
edge boundary is independently expanded modulo 3 using original vertex labels.
The moved result supports Betti [0,1,1,0], finite mass 2, raw self-distance
0 and comparison to computed-empty H1: bottleneck 1, W1-L-infinity 1,
W2-Euclidean sqrt(2). These are hand-derived expectations.

The complete existing suite also protects implicit point/matrix Rips and flag,
external cells, signed/unequal births, prime fields, H2+ fallback, sparse blockers,
duplicate ordinals, computed-empty/missing dimensions, source coverage and
per-terminal cancellation/budget retry. A failed terminal returns an error,
never a partial diagram/scalar. Reusing `Execution` starts a fresh budget on each
terminal. Context-aware comparison checks field/declared scale; physical units
and truthful dataset declarations remain caller responsibilities. In particular,
contextual self-comparison of this supplied signed source fails because it has
no declared edge-length convention. Explicit raw comparison is a caller choice;
equal operands do not bypass the contextual premise.

### Cost and exit boundaries

The adapter and assembly moves retain call, validation, reservation/clone,
checkpoint and workspace-lifetime order. Their source-normalized bodies match
the frozen baseline; source-specific readers still distinguish trusted frozen
owners from externally validated cells. Already conforming result, distance and
descriptor owners are retained, with no duplicate implementation. The new test
and this report add no runtime allocation, dispatch or dependency to the crate.

Stage selection remains O(log n) time and constant space; same-owner inclusion
validation is constant-time and traversal is O(selected cells), without allocation.
Requested dense exports, matching point buffers, descriptor outputs and opt-in
representative skeleton/transformation storage retain their respective costs.
This is a source-level cost audit, not a new timing/RSS baseline or speedup.
Future changes to actual allocation/routing/workspaces require matched measurements.

Local checks, exact-head hosted native comparisons, review, main integration and
release are separate facts. #94's implementation acceptance is reviewable on this
combined version; its final main exit happens only after approved integration.
Closing #94 does not close the future algorithms in #88-#93, reopen #95-#97 or
certify a release. #57's eventual exit record must retain these boundaries.

### Executed candidate checks

On the integrated candidate plus the new regression, Rust 1.98.1 debug/release
and Rust 1.91 tests each pass 224 tests/doctests, with four opt-in diagnostics
ignored. All-target MSRV, Clippy with warnings denied, crate and standalone
formatting, strict rustdoc, all ten examples, construction/critical-set example
tests in debug/release/MSRV, and all eleven maintained Markdown doctest pages
pass. Package construction and the 188-file package list pass. Source, Markdown,
local collaboration and diff checks pass. The same tree's Linux tool suite
passes all 98 tests. Exact-head hosted/native acceptance is recorded in the PR,
separately from these local checks.

Two initial harness failures are retained in local evidence: the first new test
incorrectly expected contextual distance on unspecified supplied scales, and
the first Linux invocation mishandled its successful test process's exit code.
The test now protects the actual rejection contract; direct WSL execution
returns success. Neither failure required a production-code change.

Reproduce with [CONTRIBUTING's verification commands](../../CONTRIBUTING.md#verification)
on the PR's exact head. The integration-specific command is
`cargo test --locked --test filtered_complex staged_results_keep_witnesses_and_consumers_after_source_release`.
Generated command logs remain under ignored `target/p4-c-spr/`; no raw artifacts
enter the crate or PR. Later stack entries add only interface documentation and
finite hand checks, with their own evidence and no algorithm acceptance claim.
