# Domain plans, strategy composition and selection

[Design index](README.md) / [Resource contracts](02-resource-contracts.md)

This design adds explicit composition inside existing domain owners. It keeps
the current direction: request/adapter selects an algorithm, algorithms execute
independently, and results are assembled by their domain. A plan need not be an
allocated object or a project-wide trait.

## Selection authority

| Owner | Owns | Does not own |
| --- | --- | --- |
| Public request | Output, field/degree, coverage, spectral target and meaningful explicit policies | A hidden weakening of requested semantics |
| Source adapter/domain dispatch | Applicable algorithm family and source interpretation | Arbitrary search across registered backends |
| Domain-local selector | Concrete implementation steps, retained state, reuse and materialization | Global mathematical policy for unrelated domains |
| Kernel | Cheap local representation/size/density branch | User output or provenance |
| Resource owner | Supply/check/reserve/reset storage | Automatically choosing a different mathematical operation |
| Maintainer admission | Which eligible alternatives can become automatic | Runtime loading of benchmark reports |

Observed examples are [flag dispatch][dispatch], [Wasserstein selection][wasserstein]
and [F2 sparse/packed storage][f2]. Preserve these existing owners. A new selector
is justified only where a real implementation choice needs composition beyond a
simple call or branch.

## What a domain plan contains

For an actual operation, resolve these facts before the repeated hot work:

1. The validated request and required output capabilities.
2. The selected implementation and its eligibility predicate.
3. Preparation to borrow, build, update or omit, with its validity dependencies.
4. Execution steps, output destinations and overlapping temporary lifetimes.
5. Content retained between steps/requests, versus content recomputed or dropped.
6. A reproducible route identifier/reason for diagnostic tooling.

Use existing request options and private function arguments when sufficient.
Introduce a small closed private enum/struct only if it resolves multiple real
choices that the caller must retain. Do not create an unused combination type,
public tuning menu or dynamic `Optimization` object. The diagnostic identifier
can be recorded by a benchmark adapter; it need not add a string allocation to
ordinary calls.

## Lowering and selection sequence

```text
validate request and source guarantees
  -> determine legal algorithm family
  -> check structural/numerical eligibility
  -> check configured resource limits
  -> apply the admitted cheap selection rule
  -> arrange preparation, scratch and output lifetimes
  -> execute compatible computation units under the same terminal budget
  -> validate and assemble the requested result
```

Use dimensions, known source metadata or already available density/rank facts
when cheap. A scan, decomposition, allocation or conversion used only to decide
the route is real work: charge it to the terminal and measure it. Never build an
expensive prepared representation merely to discover that the existing route
would have been cheaper unless that decision process itself passed a gate.

Known validated metadata can authorize legal work, but it does not establish a
performance winner. Automatic use requires source-bound evidence for an eligible
workload region. Outside that region keep the current valid route. Explicit
choices remain explicit when the evidence does not support automatic selection.

## Legal strategy axes and their coupling

| Scope | Real decisions | Coupling to preserve |
| --- | --- | --- |
| PH family | H0, F2 H1, higher cohomology or boundary route | Field, requested dimensions, access and specialized index limits |
| Boundary work | Pair-only versus retained transforms | Witness consumers and serial support form; not an independent arbitrary clearing switch |
| F2 operator | Sparse/packed storage and factorized/explicit action | Its own RREF and exact equations, not PH pairing semantics |
| Real spectral work | Fixed action, allowed range restriction, full/partial solver and owned export | Target, numerical policy, coordinate maps and current residuals |
| Persistent sequence, candidate only | Fixed endpoint, exact-event factor reuse and warm vectors | Update direction, source/basis dependencies and cumulative retained limits |
| Matching | Operand-local preparation, pair graph and solver/storage route | Pair normalization, endpoint/context contract and amortization |

Representation and kernel can be independent only where their common operations
and invariants actually agree. An axis vocabulary is not a request to implement
all representations. PHAT's algorithm/representation separation informs this
design; it does not make every Cocycle column or operator interchangeable.

## Compatibility rules used in review and experiments

| Combination | Required result |
| --- | --- |
| Specialized F2 H1 with an odd-prime request | Ineligible; select the existing supported family |
| Pairs-only state with a representative request | Ineligible; retain the needed transforms through the compatible route |
| Parallel/cleared R/V with existing critical formulas | Ineligible under current serial ELZ/lazy support requirements |
| Fixed self-adjoint sparse action with LOBPCG | Eligible only on the declared fixed real space with existing checks |
| Adaptive stopping-based LSMR action as that fixed operator | Not admitted by this refactor; a numerical-method change needs separate evidence |
| Positive restriction with full ambient spectrum output | Ineligible if it silently omits zero modes or required completeness |
| Upper factor cache with warm starts | Validate factors and coordinate mapping separately; recompute current values/residuals |
| Two scratch users with overlapping live content | Separate regions; sequential maximum-capacity reuse is ineligible |
| Repeated matching with pair-dependent normalized points cached as operand-local | Ineligible without a proven factorization of that dependency |

The executable check belongs beside the owning selector/constructor. A report
compatibility table explains it; it is not a second runtime source of truth.
Types/private constructors can exclude invalid combinations when that makes the
real API simpler. Do not encode an entire speculative capability algebra.

## Fallback and failure

A fallback is a defined branch with a valid trigger, not a catch-all retry.
Preserve the request, cumulative work already spent and allocation/retention
cost. Do not swallow cancellation, work exhaustion or arbitrary numerical
failure, silently approximate, or replay effectful callbacks. Stateful preparation
must be valid or discarded before the fallback observes it. Existing algorithmic
fallbacks retain their documented numerical/structural conditions.

## Evidence and main-path integration

Record operation/output, baseline and measured source, changed axis,
eligibility, preparation, compatible choices, failure behavior, evidence,
adoption state and actual selector rule. Benchmark both a single changed choice
and important combinations. If separating two choices changes the requested
operation, label the coupled route instead of claiming independent attribution.

Normal terminals call admitted kernels from the first implementation phase.
A new explicit benchmark-only route does not complete integration. A later
default selector change is a separate reviewable concern with regression tests
for both sides of the predicate and the unchanged output/failure contract.

[dispatch]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/persistence/flag/dispatch.rs
[wasserstein]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/diagram_distances/wasserstein.rs
[f2]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/homology_operators/matrix.rs
