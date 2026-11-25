# Numerical action, dense export and bounded iteration

[Design index](README.md) / [Development phases R1-R4](development-plan.md)

This is the first implementation slice. It changes private storage/call contracts
while preserving the existing real operators, arithmetic order, public terminals,
numerical policies and output ownership. Main and candidate code are separate.

## Observed call boundaries on main

| Owner | Current behavior | First proposed change |
| --- | --- | --- |
| [Hodge action/export][hodge] | `action_budget` allocates output and a lower temporary; `matrix_budget` invokes it per coordinate column | Fill caller output and lower scratch; dense export retains its work vectors across columns |
| [Hodge partial terminal][hspectrum] | Calls shared LOBPCG with an action returning `Vec<f64>` | Supply an output-writing action closure with captured operator-local scratch |
| [Dirac full action][dirac] | `action_budget` allocates an output and traverses stored signed incidence | Fill caller output; no unnecessary scratch owner for an action that needs none |
| [Dirac boundary modes][dspectrum] | Gram closure allocates `B*x` and `B^T*B*x`; signed verification creates temporary vectors | Reuse Gram/lift/verification temporaries under the same terminal budget |
| [Shared eigensolver][solver] | Allocates basis/images/projected and Ritz/residual blocks during iterations | Accept action output buffers, then reuse bounded live blocks and small Ritz scratch |

Main Dirac dense export already fills entries directly. Keep that path; do not
replace it with repeated action solely to make the two domains look uniform.
Main persistent Schur has a separate solver/policy; it is not automatically
consolidated into this eigensolver refactor.

## R1: output-writing leaves in existing owners

The following private shapes are illustrative; public methods remain unchanged:

```rust,ignore
// Hodge: x/y in the stored q-basis; lower in the preceding-degree basis.
fn action_into_budget(
    &self,
    x: &[f64],
    y: &mut [f64],
    lower: &mut [f64],
    budget: &mut WorkBudget<'_>,
) -> Result<()>;

// Dirac's full graded action needs only its output.
fn action_into_budget(
    &self,
    x: &[f64],
    y: &mut [f64],
    budget: &mut WorkBudget<'_>,
) -> Result<()>;
```

Use direct slices before introducing a scratch type. If multiple real local
consumers need several buffers together, a private typed owner may reduce
argument/lifetime complexity. It does not become a public Workspace or shared
mathematical operator interface.

| Action | Output | Additional live temporary | Initialization |
| --- | --- | --- | --- |
| Hodge full/lower | `n = dim C_q` scalars | `m = dim C_(q-1)` scalars | Clear output and lower accumulation every use |
| Hodge upper | `n` scalars | Per-column scalar dot, no lower vector required | Clear output; preserve upper incidence traversal |
| Full truncated Dirac | `N = total retained coordinates` scalars | None beyond scalar arithmetic | Clear output before signed two-sided accumulation |
| Dirac selected boundary Gram | `n = selected upper-degree coordinates` scalars | `m = selected lower-degree coordinates` scalars | Clear `B*x` and Gram output; preserve offset mapping |

Validate input/output shape and finite inputs once along each route. Preserve
public rejection behavior and avoid repeating immutable incidence/chain
validation. The allocation-returning wrapper and the output-writing internal
caller must not accidentally double-scan the input; a small common validation
helper is acceptable if needed to preserve public validation/allocation order.

Initialize every value read after reuse, retain finite/canonical-zero checks and
existing accumulation order, and check the same cumulative budget. On error the
internal output may be dirty; public methods return no successful partial vector.
An earlier owned result and the immutable operator remain valid.

Hodge dense export allocates its intentional `n*n` matrix once after the existing
entry cap/overflow checks, then reuses `x`, `y` and lower scratch for all columns.
Do not bypass the dense output limit or change its retained matrix semantics.
The ordinary `action`, `matrix` and spectrum entry points use the new leaves;
normal-path integration is part of R1, not a later optional facade.

## R2a: change the shared action callback

Replace the private allocation-returning callback with this output-writing shape:

```rust,ignore
mut action: impl FnMut(
    &[f64],
    &mut [f64],
    &mut WorkBudget<'_>,
) -> Result<()>
```

The solver supplies an exactly sized output and requires complete overwrite.
Hodge captures its lower temporary; Dirac boundary work captures its Gram
temporary. Update both current callers and the eigensolver's colocated test
callbacks in one buildable change. A temporary private allocating adapter, if
needed to sequence commits, must not remain the ordinary route after R2a.

Keep count-zero early return, deterministic starting seed, rank checks, two-pass
orthogonalization, stopping tolerance, tie/sort policy and current residual checks.
The callback change does not introduce new spectral targeting or a preconditioner.

## R2b: reuse live solver blocks without changing the method

Let `b=min(block_size,n)` and `c=min(3*b,n)` with checked multiplication. Active
search dimension `m` can vary up to `c`. Existing `max_basis_entries` bounds each
specified block/projected buffer; it is not an aggregate bound on the following
whole live set.

| Buffer/content | Capacity bound | Important overlap |
| --- | --- | --- |
| Current and old Ritz vectors | Each at most `n*b` | Old vectors survive until new directions are formed |
| Residual and direction blocks | Each at most `n*b` | Search-basis formation reads both and current vectors |
| Search basis and action images | Each at most `n*c` | Both are live through Rayleigh-Ritz and cannot alias |
| Projected matrix | At most `c*c` | Its contents coexist with the basis/images |
| Jacobi working matrix and rotations | Each at most `c*c` | Inventory actual copying/output use before reusing either |
| Ritz values, errors and sort indices | Bounded by `b` or `c` | Preserve deterministic active-prefix interpretation |
| Fresh verification vector | `n` | Reused across returned-vector checks, using fresh action applications |
| Owned result | `n*count` plus requested scalar metadata | Moves out only after successful verification |

Start by retaining existing layouts and operation order. Reuse inner vector
capacity as well as outer lists: clearing a `Vec<Vec<f64>>` drops its inner
allocations, so it does not by itself solve repeated allocation. Retain slots
and an active prefix or rotate owned buffers after their last use. A private
block helper is justified only if it simplifies real repeated operations and
does not force a layout/numerical-method change.

Reserve lazily for the active bounded shape or compare that choice against
maximum reservation. Preallocating every worst-case block can regress tiny,
count-zero or early-converged solves. Growth may occur when active capacity
increases; after sufficient reservation, fixed-shape scratch allocation should
not scale with iteration count. Checked entry/byte limits and fallible growth
still apply before allocation.

Fresh final action applications must remain independent of cached Ritz-image
linear combinations. Reusing the verification buffer is legal; reusing old
residuals as current acceptance is not. Do not fuse/reorder floating operations,
change rank tolerances or merge numerical policies in this storage-only phase.

## R4: persistent action and factors, candidate dependency

After PR 160's code is available on the agreed base, [candidate action][paction]
needs output `n`, upper temporary `p`, dense projection coefficients `r`, and
lower temporary `m`, with overlapping upper temporary/coefficients during fixed
projection. The compact route has its own checked apply behavior. Preserve all
projection coefficients against the original upper temporary before subtraction.

Apply the output-writing callback to smallest/largest/positive and Gram routes
without changing targets, supplied initial vectors, zero-threshold/range policy
or lifted residuals. Separate current lower state from reusable upper content
using [the factor design](04-preparation-and-sequences.md). Do not substitute
adaptive LSMR, ordinary HeCS, sampling or research subspace methods into a fixed
action as part of this refactor.

## Acceptance and rollback

Use independent small signed-chain/real matrix expectations, the existing Hodge
and Dirac reference workers, and current original-coordinate residual checks.
Exercise empty/zero-count, lower/upper/full, truncated degrees, repeated roots,
near-zero/clustered cases, failed convergence, shape/limit/cancellation and dirty
reuse. Preserve errors rather than returning incomplete successful modes.

Compare ordinary action, dense export and partial solve against the pinned
pre-change ordinary entry points. Isolate action-output reuse from solver-block
reuse and then measure their combination, including setup/drop, copied/zeroed
bytes and live/retained capacity. Include tiny/cold, early-converged and difficult
iteration controls. Full spectrum versus selected partial output is a separate
capability comparison, not the allocation-refactor ablation.

R1 can be reviewed independently of R2; R2a can be reviewed before R2b. Revert a
failed buffer layout/retention policy within its owner while keeping compatible
validated leaves. Do not leave duplicate public routes or an unused generic
workspace after rejecting a candidate.

[hodge]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/hodge/mod.rs
[hspectrum]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/hodge/spectrum.rs
[dirac]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/dirac/mod.rs
[dspectrum]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/dirac/spectrum.rs
[solver]: https://github.com/Aequiludium/cocycle-rs/blob/a3c3e1132800624743ad0f159c9303090e0b3b3f/src/eigensolver.rs
[paction]: https://github.com/Aequiludium/cocycle-rs/blob/21bd4d9c1fbda28bf0c8107cd663a0e5e17cf851/src/persistent_laplacian/operator.rs
