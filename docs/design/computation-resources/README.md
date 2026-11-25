# Computation and resource architecture proposal

[Documentation](../../README.md) / Design

Status on 2026-11-25: proposed responsibilities and staged implementation for
review. This package does not change production code, public APIs, defaults,
dependencies or MSRV. The [current architecture](../../development/architecture.md)
and [kernel design](../kernel.md) remain authoritative for implemented behavior.

The execution-side organizing concept is a **computation unit with an explicit
resource contract**. Mathematical types define meaning and guarantees;
domain-owned plans compose executable work; typed storage and explicitly owned
preparation permit valid reuse. The proposal adds this boundary to the existing
source/algorithm/result architecture, with substantive internal changes in
selected owners rather than a replacement of all algorithms.

## Read the proposal

The [overview](overview.md) contains the current module/call map, target
architecture and lifecycle diagrams, the seven direct architecture answers,
worked plans and the existing negative evidence that constrains adoption.

| Design | Implementation question it resolves |
| --- | --- |
| [01. Semantics and computation boundaries](01-semantics-and-computation.md) | What is a mathematical object, an execution unit, and an ownership boundary in every current domain? |
| [02. Resource contracts](02-resource-contracts.md) | What must an internal kernel expose about output, scratch, aliasing, growth, budgets and failure? |
| [03. Plans and selection](03-plans-and-selection.md) | How are real strategies composed, selected, identified and rejected? |
| [04. Preparation and sequence lifetimes](04-preparation-and-sequences.md) | Which retained content is valid, who owns it, and how do updates/retries work? |
| [05. Numerical action and iteration](05-numerical-path.md) | How do Hodge, Dirac, LOBPCG and eligible persistent paths acquire explicit output/scratch contracts? |
| [06. Exact algebra and critical queries](06-exact-algebra.md) | What can PH, representatives, critical sets and F2 operators share without changing their algebraic form? |
| [07. Construction and analysis consumers](07-construction-and-consumers.md) | How do source preparation, output ownership, descriptors and diagram matching fit without mandatory workspaces? |
| [08. Verification and admission](08-verification-and-admission.md) | What evidence admits an internal boundary, a strategy combination or a default change? |
| [Development plan](development-plan.md) | What are the phases, dependencies, review units, exit conditions and rollback boundaries? |

Read 01-04 before changing a cross-module contract. Start with the numerical
vertical path in 05 and the evidence requirements in 08. Designs 06-07 are
separate domain assessments, not instructions to force every operation through
one representation. The development plan owns phase dependencies/status; the
topic pages own contracts.

## Version and evidence boundaries

The implemented reference is main at
`a3c3e1132800624743ad0f159c9303090e0b3b3f`. Candidate spectral integration
[PR 160](https://github.com/Aequiludium/cocycle-rs/pull/160) was open at
`21bd4d9c1fbda28bf0c8107cd663a0e5e17cf851` when the proposal was prepared.
Candidate-only action/factor/sequence work remains a conditional dependency in
R4; it is not described as an implemented main capability.

The overview links fixed source/report revisions. Existing studies are evidence
for their recorded workloads and measured sources, not a new measurement of the
proposed interfaces. The proposal preserves rejected distance pools, conditional
operand preparation, withdrawn uniform-boundary routing and slower small
spectral batches as design constraints.

## Adoption boundary

The initial changes are private interfaces in existing owners. Ordinary public
terminals must use admitted kernels; specialized algorithms, cheap direct
functions, mathematical guarantees and owned result lifetimes remain deliberate.
No global optimizer, mandatory public workspace, runtime registry or new crate is
required by this design.

Profiling and independent validation must still establish whether each boundary
or retention policy earns adoption. A failed gate keeps the simpler owner-local
route and records the negative result. Publishing this design is distinct from
implementing the refactor, selecting a new default, merging code or releasing it.
