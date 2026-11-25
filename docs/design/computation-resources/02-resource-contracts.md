# Internal resource contracts

[Design index](README.md) / [Numerical first implementation](05-numerical-path.md)

This is the proposed contract for a useful internal composition boundary. It
starts with ordinary functions, mutable slices and private typed vectors.
The [existing execution owner][execution] remains work controls, not a new
allocator or aggregate RSS quota.

## Required information, without a universal trait

| Item | Required declaration | Example |
| --- | --- | --- |
| Inputs | Shape, layout, basis, immutable borrow and preparation assumptions | One f64 vector in the operator's stored coordinate order |
| Outputs | Length/layout, full overwrite versus accumulation, ownership after success | Fill exactly `n` output coordinates |
| Aliasing | Which inputs, outputs and temporaries must remain distinct | Input cannot alias mutable output; live basis and images coexist |
| Scratch | Exact size, upper bound, or explicitly fallible data-dependent growth | Lower action length `m`; solver block at most `n * min(n,3*b)` |
| Initialization | Every value read was initialized during this use or is valid retained content | Clear an accumulator before adding sparse terms |
| Retention | Content lifetime versus capacity lifetime; last-use point | Reuse capacity after a column/iteration without reusing stale values |
| Failure | Output/scratch state on error and validity of retained preparation | Internal partial output discarded; immutable preparation stays usable |
| Work | One cumulative terminal budget, including initialization/preparation | Nested actions receive the caller's budget |

Document this next to the affected function. A private helper or small struct
is justified by an actual caller; the table is not a demand for eight new types.

## Initial numerical interface shape

```rust,ignore
fn action_into_budget(
    &self,
    input: &[f64],
    output: &mut [f64],
    lower_work: &mut [f64],
    budget: &mut WorkBudget<'_>,
) -> Result<()>;
```

This is illustrative, not an implemented signature. Hodge needs a lower
intermediate when the lower term is selected. Dirac's direct graded action needs
only its output, so it should not receive a meaningless scratch argument. A
small owner-local scratch struct is appropriate when several buffers genuinely
travel together. Keep requirements as checked local calculations until multiple
consumers justify a common requirements helper.

The existing owned-output public wrapper creates one budget, validates/reserves
the intentional result and required scratch, calls the internal function, and
returns the output only after success. An internal iterative caller supplies an
already reserved output region and passes its own budget. The kernel never
starts another terminal budget or returns a borrow into temporary local storage.

## Sizing, capacity and initialization

Compute products, sums and byte capacities with checked arithmetic before
reserving. Apply current per-buffer caps before allocating; an overflow is not
an allocation attempt. Use fallible reserve behavior compatible with the owning
domain's existing error policy. `Vec::capacity()` measures container capacity,
not allocator metadata or process high-water memory.

Safe slices contain initialized elements. A vector with spare capacity must be
resized/initialized before creating the working slice. Avoid introducing unsafe
uninitialized storage to remove zeroing. An accumulation buffer must be reset
for each use even if its allocation is retained. A buffer that will be completely
overwritten need not preserve its previous content, but no read may precede the
new write. Preserve finite checks and canonical-zero policies at their existing
semantic boundary.

The implementation may grow data-dependent sparse state fallibly; it must not
claim a fixed precomputed size for reduction fill-in. A reserve failure can leave
larger capacity in some already-reserved scratch vectors without producing valid
prepared content. The next operation reinitializes all required values.

## Composition and simultaneous liveness

Compatible sequential temporary requirements with disjoint live intervals can
reuse one typed backing vector sized to their maximum. Simultaneous uses require
disjoint slices or separate vectors, counted together. This includes retained
factors, caller-owned outputs and any old/new blocks needed during an update.
Do not add standalone phase peaks and call the sum the workflow peak.

Start with explicitly named buffers in the domain owner. Rotating two or three
solver blocks using moves/swaps can avoid allocation while preserving overlap
requirements. A raw-byte arena, automatic graph liveness solver and global pool
are not needed for these first changes. Resetting capacity is different from
releasing it; a long-lived scope must expose or document its high-water retention
and bounded eviction/drop behavior.

## Failure, cancellation and reuse

| Event | Internal output/scratch | Valid retained content | Public terminal |
| --- | --- | --- | --- |
| Shape or option rejection | No useful output promised | Unchanged | Structured input error |
| Allocation/growth failure | May be partially reserved/initialized | Unchanged unless a mutable update was started | Allocation/size error; no successful partial output |
| Cancellation/work exhaustion | May contain partial arithmetic | Existing immutable content survives; incomplete new cache entry is absent | Execution error |
| Nonfinite arithmetic/nonconvergence | Dirty scratch | Old valid factors survive if not mutated | Owning numerical error |
| Failed incremental factor update | Cannot be treated as complete | Consume/discard the failed state or mark invalid before retry | No partially updated valid handle |
| Success | Needed output complete | Publish only complete new content | Return the requested owned result |

The first private `_into` kernels permit unspecified output-on-error because
their caller discards incomplete output. Do not expose that behavior as a public
API without a separate contract. Returning `Err` is not permission to invalidate
an earlier owned result or leave a cache marked ready.

Each public terminal starts a fresh `WorkBudget`; its nested work shares that
budget. Retained scratch never retains used work counters. A valid-failed-valid
sequence must work with fresh terminal budgets. Work counts describe the
implementation and may change when real work is eliminated; preserve documented
limit/cancellation behavior rather than claiming exact cross-algorithm counts.

## Concurrency and callback effects

Immutable preparation can be shared when its existing ownership allows it.
Mutable scratch requires exclusive borrowing, disjoint regions or separate
per-execution storage. A single locked workspace would serialize users and is
not the initial concurrency design. No static mutable pool is introduced.

Resource reuse does not establish mathematical purity. A stateful callback
builder remains consuming. Caching, fusion or retries may not replay, duplicate
or reorder arbitrary callbacks. Preserve the access/effect order unless an
existing contract explicitly permits a change.

## Acceptance evidence

Exercise wrong shape, zero-size domains, dirty-capacity reuse, growth/overflow,
cancellation and valid-failed-valid execution. Verify old owned outputs remain
unchanged. For fixed-capacity hot buffers, allocation/growth observations after
reservation should be independent of action/iteration count; necessary output
allocation and dynamic fill-in are counted separately. Record reset/zeroing and
retention costs, then apply the [full admission protocol](08-verification-and-admission.md).

[execution]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/execution/mod.rs
