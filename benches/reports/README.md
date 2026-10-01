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
| [Diagram-distance preparation and allocation effects](distance-preparation.md) | Finite W1 counterexample, ordinary API comparison and isolated allocation interventions; draft pending independent reproduction |
| [Integrated Phase-2 performance baseline](phase2-baseline.md) | Frozen Direction-B R0, prepared Rips, complete Rips workflows and Bottleneck/W1/W2 resources; local evidence and replay contract |
| [Persistence diagram assembly overhead](diagram-assembly.md) | Fixed-R0 private phase and output-storage diagnostics; focused H0 follow-up and retained timing-mode discrepancy |
| [Diagram-distance state lifetimes](distance-lifecycle.md) | Fixed-R0 persistent preparation/workspace ablations, cross-platform break-even, retention and selected worker-throughput evidence |

Do not create a report per date, commit, PR or rerun. Related suites share a report
with separate timing contracts and tables. Keep concise conclusions and selected
results here; raw samples, fixtures, logs, binaries and archives stay outside Git.
Local-only evidence is explicitly identified and is not publicly archived.
Routine CI results belong in PRs and job output, not additional report files.
