# Verification, measurement and default admission

[Design index](README.md) / [Development plan](development-plan.md)

This design applies the repository's existing validation and
[performance reporting rules][reporting] to computation/resource changes.
Build readiness, algorithm correctness, performance admission, merge and release
are separate claims. No new oracle run or benchmark has been performed for this
design package.

## Required evidence for each implementation unit

| Evidence | Required result | Scope of claim |
| --- | --- | --- |
| Contract/source review | All changed callers, read/write/alias/retention/failure responsibilities identified | Architecture and ownership |
| Resource behavior checks | Valid-failed-valid reuse, initialized reads, bounded growth and unchanged live outputs | Resource correctness, not numerical accuracy |
| Domain correctness | Independent equations/fixtures/oracles and current public integration paths agree | Exercised mathematical contracts |
| Resource observations | Growth/allocation, copied/zeroed bytes, live/retained capacity and process memory kept distinct | Measured mechanism at actual revision |
| End-to-end comparison | Same requested work and declared complete lifetime boundary; regressions retained | Scoped performance/adoption |
| Exact-head project/CI checks | Applicable full repository checks pass for the submitted head | Integration, not release |

An allocation counter cannot establish spectral accuracy. A residual check
cannot establish lower allocation or peak memory. A callable new route is not
evidence that normal terminals use it. Each claim needs its own observation.

## Correctness matrix

| Domain | Independent evidence | Important controls |
| --- | --- | --- |
| Hodge/Dirac action | Small signed-chain dense matrices and direct hand/rational expectations | Lower/upper/full, degree offsets, orientation and empty domains |
| Numerical eigenpairs | Current-operator residuals, orthogonality, real matrix/SVD references and existing worker comparisons | Clustered/repeated/near-zero modes, count zero, rank loss, nonconvergence and retained targets |
| Persistent candidate | Restricted-domain real/rational matrices, current lower term and original pair residuals | Positive lifting, compact eligibility, tied/repeated/unordered pairs and warm remapping |
| PH/representatives | Independent reduction/witness equations plus required GUDHI/Ripser supported cases | Fields, zero bars, births/ties, cutoffs, coverage and specialized routes |
| Critical/F2 queries | Exact equations, exhaustive tiny controls and source/run/support checks | Serial form, lazy failure/retry, sparse/packed transitions and finite histories |
| Construction/consumers | Owning semantic fixtures and required GUDHI/Topp/Ripser comparisons | Callback effects, approximate provenance, dimensions/endpoints and extreme values |

GUDHI topology checks do not certify real numerical rank or spectrum. Repeated
roots need invariant-subspace/equation checks where vector coefficients are not
unique. Never make the refactored kernel its only expected-value generator.

## Baselines and attribution

Before implementation, freeze the ordinary entry-point baseline SHA, inputs,
output/guarantees, limits, toolchain and worker boundary. Before formal measurement,
commit the implementation/harness, verify measured inputs, and record candidate
and baseline identities separately. Exploratory dirty runs stay local with their
actual fingerprints; they cannot later be relabeled as a commit measurement.

For numerical resource work, compare action output reuse alone, callback change
alone where meaningful, solver-block reuse alone, and the admitted combination.
Keep arithmetic/target/precision fixed. Where two changes cannot be separated
without changing the operation, state the coupled unit. For algorithmic
acceleration studies, also identify a matched-output disabled/unaccelerated
control separately from the current production baseline.

Do not attribute full-versus-partial solver cost to scratch reuse. Do not add
independent speedup ratios or subtract phase medians to infer missing costs.
Include selection, validation, preparation, conversions, execution, verification,
owned output, reset and declared destruction. A warm-query-only diagnostic is
reported separately from complete scope cost.

## Resource observations and fair execution

Observe allocation events only with an allowed external allocator/profiler or
existing permitted instrumentation. Safe internal capacity/growth counters are
useful but are labeled as growth/reserve events, not actual malloc counts. Do not
relax `unsafe_code = forbid`, add a global allocator or new runtime dependency
just to obtain a number.

Record requested/active capacity, retained high water, copied/zeroed bytes and
maximum simultaneously live known containers. Record process RSS separately;
allocator retention can make RSS differ from logical capacity. Include cache
content, live outputs and simultaneous callers when making a memory claim.

For formal before/after studies use the owning frozen protocol, serial workers,
balanced order and retained samples. The project reporting minimum is ten
independent measured processes per compared cell; distance uses its stricter
existing study contract. Record warmups/spread, environment, failures/timeouts
and every attempt. Do not measure while compiling or testing concurrently.
Use Linux for the initial formal resource comparisons and record Windows
results separately. State the platform and controls in the frozen protocol; do
not pool observations across different environments.

Preregister workload families, repetitions, eligible-region predicate and the
owning performance acceptance rule before measuring. Existing VR/distance
thresholds are study-specific; they are not silently assigned to numerical work.
For the resource refactor, fixed-shape hot growth must stop scaling with repeated
actions/iterations after capacity is sufficient, and the complete workload must
show useful benefit without an unaccepted regression. If evidence is noisy,
negative or lacks a defined gate, keep the current default and record no adoption.

## Adoption state and review record

```text
Reference -> Experimental -> Validated opt-in -> Candidate dispatch -> Default
```

Any stage can move to rejected/deferred with the reason retained. Track state
for each strategy and important combination, not only for its containing PR.
The record includes operation/output, changed choice, preparation/retention,
eligibility/compatibility, actual measured sources, correctness/performance
evidence, selector rule and fallback. Maintain it in the existing topic report
and benchmark metadata rather than a new runtime database.

Default selection needs a cheap prospective applicability rule, full lifetime
cost, interaction/holdout/repeat evidence under its declared gate, regression
controls and a valid same-contract fallback. Failed cancellation/limits are not
retry triggers. Keep rejected distance pools, withdrawn uniform-boundary routing,
draft adaptive H1 and slower candidate small spectral batches in the applicable
regression matrix.

## Existing commands and owning workflows

Use the affected paths in `CONTRIBUTING.md`, then full checks for shared-kernel
changes. These are future implementation checks, not checks claimed run here:

```sh
cargo fmt --all -- --check
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo test --locked --all-features
cargo test --locked --release --all-features
cargo +1.91.0 check --locked --all-targets --all-features
cargo +1.91.0 test --locked --all-features
python3 tools/check_source.py
python3 tools/check_docs.py
python3 tools/check_artifacts.py
```

| Affected path | Current reference / integration entry |
| --- | --- |
| Hodge | `examples/hodge_laplacian.rs`, `tools/check_hodge.py --rust PATH --output PATH --iterative`; full pins/flags in `.github/workflows/hodge-reference.yml` |
| Dirac | `tools/reference/dirac_cocycle.rs`, `tools/check_dirac.py --rust PATH --output PATH --pairs`; `.github/workflows/dirac-reference.yml` |
| Main persistent Laplacian | `examples/persistent_laplacian.rs`, `tools/compare_persistent_laplacians.py --worker PATH --output PATH`; `.github/workflows/persistent-laplacian-reference.yml` |
| PH / construction | `tools/check_algorithm.py persistence-reduction` or `complex-construction`, plus required native comparisons |
| Diagram analysis / matching | `tools/check_algorithm.py diagram-analysis` or `diagram-distances`, plus the declared independent reference protocol |
| Exact operators | `tools/check_homology_operators.py` and existing operator/family tests; follow its actual CLI before running |

PATH entries above are required real input/output paths, not literal ready-to-run
commands. Preserve pinned references and required flags from the owning workflow.
Do not download/install or rerun external sources merely to check local Markdown.

## Documentation and handoff checks

Run `python3 tools/check_source.py`, `python3 tools/check_docs.py` and the
Git-index artifact check after staging. Audit immutable source links at their
recorded revisions and check local anchors. Run affected Rust examples/doctests;
illustrative `rust,ignore` signatures in this proposal are not implemented APIs.

Keep machine-specific worktree paths, credentials, generated evidence and raw
runs outside the documentation contribution. The PR records checks actually
run for its head, separately from future implementation requirements above.
Publishing this proposal establishes no new algorithm-correctness, MSRV,
external-oracle, performance or default-adoption result.


[reporting]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/benches/reporting.md
