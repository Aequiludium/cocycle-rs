# Cocycle

![Cocycle: topology, computed in Rust](../assets/banner.svg)

Cocycle is a general-purpose topological data analysis library in pure Rust,
with no runtime dependencies and no unsafe code. Compute persistence from point
clouds, dissimilarities or supplied complexes, then inspect diagrams, compare
them and derive summaries.

## Start with a square

Install the released crate with `cargo add cocycle@0.1.1`. The
[quickstart](getting-started/quickstart.md) computes a square's loop and its
Betti curve; the [installation guide](getting-started/installation.md) explains
the Rust toolchain and how to use a source checkout.

![A square's Rips filtration](../assets/filtration.svg)

## Choose a task

| Task | Read |
| --- | --- |
| Compute persistent homology from points or a matrix | [Rips persistence](guides/rips.md) |
| Build a graph or inspect a frozen simplicial complex | [Rips construction](guides/rips-construction.md) |
| Analyze a supplied filtered complex | [Filtered complexes](guides/filtered-complexes.md) |
| Choose a field and request representatives | [Fields and representatives](guides/rips-representatives.md) |
| Compute an explicitly labeled approximation | [Sparse Rips](guides/sparse-rips.md) |
| Compare complete persistence diagrams | [Diagram distances](guides/diagram-distances.md) |
| Reuse stage inputs and owned results | [Interface contracts](guides/interface-contracts.md) |
| Propose bounded F2 endpoint changes | [Critical sets](guides/critical-sets.md) |
| Query feasible F2 homology operators | [Homology operators](guides/homology-operators.md) |
| Review a real single-stage operator | [Hodge Laplacians](guides/hodge-laplacian.md), [Dirac operators](guides/dirac.md) (available on main) |
| Review a two-stage real operator | [Persistent Laplacians](guides/persistent-laplacian.md) (available on main) |

## Versions and guarantees

The site documents **main**, including additions labeled **Unreleased**.
`cocycle = "0.1.1"` selects the released crate; it does not include every main
API. The real operator guides describe Unreleased main interfaces. See the
[version and capability status](README.md#version-and-capability-status)
before using a source-only feature, and use the [API reference](reference/api.md)
for versioned signatures and contracts.

Rips scales are edge lengths. Finite, essential and right-censored intervals
remain distinct; a cutoff does not certify an essential class. The
[mathematical specification](reference/mathematics.md) owns definitions and
numerical guarantees. Design and research pages retain their proposal, source
revision and evidence status.

```{toctree}
:hidden:
:caption: Getting started

getting-started/installation
getting-started/quickstart
Documentation map and version status <README>
```

```{toctree}
:hidden:
:caption: User guides

guides/rips
guides/rips-construction
guides/filtered-complexes
guides/rips-representatives
guides/sparse-rips
guides/diagram-distances
guides/interface-contracts
guides/critical-sets
guides/homology-operators
```

```{toctree}
:hidden:
:caption: Real operators

guides/hodge-laplacian
guides/dirac
guides/persistent-laplacian
```

```{toctree}
:hidden:
:caption: Reference

reference/api
reference/mathematics
reference/bibliography
```

```{toctree}
:hidden:
:caption: Development

development/algorithm-contributions
development/diagram-analysis
development/complex-construction
development/persistence-reduction
development/architecture
development/conventions
development/testing
Contributing <https://github.com/Aequiludium/cocycle-rs/blob/main/CONTRIBUTING.md>
```

```{toctree}
:hidden:
:caption: Design and roadmap

design/roadmap
design/rips-api
design/rips
design/rips-acceptance
design/rips-implementation
design/kernel
design/interface-compatibility
```

```{toctree}
:hidden:
:caption: Research and evidence

research/route-integration
research/interface-inventory
research/input-contracts
research/result-contracts
research/data-modeling
research/filtration-stages
research/analysis-learning
research/persistent-spectral-features
research/gudhi-cpp
research/h1-fast-to-generic
research/big-steps
Benchmarks <https://github.com/Aequiludium/cocycle-rs/tree/main/benches>
```

```{toctree}
:hidden:
:caption: Project

Changelog <https://github.com/Aequiludium/cocycle-rs/blob/main/CHANGELOG.md>
Source and issues <https://github.com/Aequiludium/cocycle-rs>
MIT license <https://github.com/Aequiludium/cocycle-rs/blob/main/LICENSE>
Third-party notices <https://github.com/Aequiludium/cocycle-rs/blob/main/THIRD_PARTY_NOTICES.md>
```
