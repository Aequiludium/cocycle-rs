# API reference

[Documentation](../index.md) / Reference

Rustdoc owns individual signatures, errors and API contracts. Use the
[0.1.1 reference](https://docs.rs/cocycle/0.1.1/cocycle/) when your dependency is
`cocycle = "0.1.1"`; [docs.rs](https://docs.rs/cocycle) follows the published
crate and may differ from this site's main guides.

For the current checkout, run `cargo doc --locked --no-deps --open` from the
repository root. The generated entry is `target/doc/cocycle/index.html`.
Source-only capabilities are identified in the
[version status](../README.md#version-and-capability-status).

| API area | Start with |
| --- | --- |
| Borrowed points and dissimilarities; Rips and flag construction | [Rips construction guide](../guides/rips-construction.md) |
| Supplied complexes, signed boundaries, stages and inclusions | [Filtered complexes](../guides/filtered-complexes.md) |
| Persistence requests, fields and representative bases | [Rips guide](../guides/rips.md), [representatives](../guides/rips-representatives.md) |
| Diagrams, endpoints, coverage and result ownership | [Interface contracts](../guides/interface-contracts.md) |
| Diagram summaries and Betti curves | [Diagram analysis](../development/diagram-analysis.md) |
| Bottleneck, W1 and W2 matching | [Diagram distances](../guides/diagram-distances.md) |
| F2 critical sets and partial U | [Critical sets](../guides/critical-sets.md) |
| Feasible F2 projections and finite transport | [Homology operators](../guides/homology-operators.md) |

The [mathematical specification](mathematics.md) owns shared definitions;
[architecture](../development/architecture.md) locates their implementation.
Design and research sketches do not establish callable APIs.
