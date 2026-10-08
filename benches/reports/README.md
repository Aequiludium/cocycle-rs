# Benchmark reports

[Benchmarks](../README.md)

These maintained comparisons summarize tested capabilities, correctness and
performance. Each report uses a stable topic-based filename; reruns update the
same file and Git history preserves prior versions. The report body records the
actual measured commits, protocols, environment and evidence status. Follow the
[reporting rules](../reporting.md) and [template](../report-template.md).

| Report | Scope |
| --- | --- |
| [Cocycle, GUDHI and Ripser: Rips comparison](rips-comparison.md) | Native C++ references, Rips correctness, selected workflow and H0/H1 timings, limitations and reproduction |
| [Adaptive F2 H1 kernels](f2-h1-adaptive.md) | Paired committed Rust baseline/candidate and native Ripser; complete repeated results, default-policy regressions and merge limitations |
| [Exact F2 H2 foundation](h2-foundation.md) | Phase-3 baseline, counter schema and H1-to-H2 continuation evidence |
| [H1-to-H2 lifetime and scratch evidence](h2-workspace.md) | M1 ownership contract, lifetime-only and same-type scratch ablations, retained regressions and experimental decision |
| [Serial F2 transfer and high-dimensional decision](f2-transfer.md) | Frozen T4-T9 outcomes, M2 typed/reset no-go, fresh H0-H3/fallback validation and T11 no-follow-up |

Do not create a report per date, commit, PR or rerun. Related suites share a report
with separate timing contracts and tables. Keep concise conclusions and selected
results here; raw samples, fixtures, logs, binaries and archives stay outside Git.
Local-only evidence is explicitly identified and is not publicly archived.
Routine CI results belong in PRs and job output, not additional report files.
