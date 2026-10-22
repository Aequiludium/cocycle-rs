# Bibliography

[Documentation](../README.md) / Reference

These sources identify definitions, theorems, and algorithmic ideas. They are not
runtime dependencies. Cocycle's engineering conventions are stated separately in
the [mathematical specification](mathematics.md).

## NM22

Arnur Nigmetov and Dmitriy Morozov. **Topological Optimization with Big Steps**.
[arXiv:2203.16748v2](https://arxiv.org/abs/2203.16748v2).
Theorems 9-11 and Algorithms 1-3 motivate the singleton critical sets and
maximum-displacement surrogate. The
[specification](mathematics.md#19-big-steps-critical-sets-over-f2) records
Cocycle's F2 source and lazy-reduction assumptions.

[Oineus source](https://github.com/anigmetov/oineus/tree/e52814a1ffb5b8a81e71ff1e93b4c14194673f0f)
at revision `e52814a1ffb5b8a81e71ff1e93b4c14194673f0f` provides the
implementation reference for bounded partial U rows (`decomposition.h`) and
the `Max` rule (`top_optimizer.h`). Cocycle independently implements the
algorithm using its existing sparse reducer; it does not import Oineus code
or add it as a runtime dependency. Runtime comparisons and mathematical
theorems are distinct evidence.

## CSB23

Lucille Calmon, Michael T. Schaub and Ginestra Bianconi.
**Dirac Signal Processing of Higher-Order Topological Signals**.
New Journal of Physics 25, 093013 (2023).
[arXiv:2301.10137v2](https://arxiv.org/abs/2301.10137v2), section 2 for the
standard real graded Dirac definition, square and signed spectral pairing.
The [specification](mathematics.md#24-standard-real-discrete-dirac) fixes the
ordinary-chain and degree-truncation conventions used here; signal filtering
and persistent/weighted variants are outside this implementation.

[TopoNetX signed incidence and Dirac methods](https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/toponetx/classes/simplicial_complex.py)
at `ce3d414b571fff95916c4cbc20b8d2c2c9d428e6` provide an independent matrix
reference. The reference's rank-zero incidence is augmented; Cocycle's Dirac
coupling starts at B1. Its `weight` argument is not evidence of a metric adjoint.

## ZC05

Afra Zomorodian and Gunnar Carlsson. **Computing Persistent Homology**.
Discrete & Computational Geometry 33, 249–274 (2005).
[Publisher / DOI](https://doi.org/10.1007/s00454-004-1146-y).

Background for persistence over fields and interval classification. The project's
original source review checked bibliographic information and the abstract; the
specific reduction and pairing rules use the full-text source B21 below.

## B21

Ulrich Bauer. **Ripser: efficient computation of Vietoris–Rips persistence
barcodes**. Journal of Applied and Computational Topology (2021).
[arXiv v2](https://arxiv.org/abs/1908.02518v2),
[full text](https://arxiv.org/html/1908.02518v2),
[DOI](https://doi.org/10.1007/s41468-021-00071-5).

Relevant locations: §2 for Rips and filtration ordering; §3.1, Proposition 3.1 and
Algorithm 1 for reduction and pairing; §§3.2–3.5 for clearing, cohomology, implicit
matrices, and apparent pairs; §4 for H0/union-find. Cocycle independently implements
these invariants in Rust. Its internal order must satisfy each optimization's
conditions; citing Ripser does not establish correctness automatically.

## CSO13

Frédéric Chazal, Vin de Silva and Steve Oudot. **Persistence stability for
geometric complexes**.
[arXiv v3 (2013)](https://arxiv.org/abs/1207.3885v3),
[full text](https://arxiv.org/html/1207.3885v3).

Theorem 2.3 connects interleavings and bottleneck stability. Theorem 5.2 bounds
Rips diagrams using Gromov–Hausdorff distance. The fixed-vertex perturbation
argument in Cocycle first derives filtration inclusions, then applies Theorem 2.3.

## A20

Nieves Atienza, Rocío González-Díaz and Manuel Soriano-Trigueros.
**On the stability of persistent entropy and new summary functions for
Topological Data Analysis**. Pattern Recognition 107, 107509 (2020).
[arXiv v7](https://arxiv.org/abs/1803.08304v7),
[full text](https://arxiv.org/html/1803.08304v7),
[DOI](https://doi.org/10.1016/j.patcog.2020.107509).

Definition 3.1 uses base-two logarithms; §3.4 discusses infinite intervals.
Cocycle's natural logarithms, finite-only selection, and empty-result `None` are
explicit project conventions. The paper's title is not an unconditional stability
guarantee for every entropy variant.

## HG19

Jakob Hansen and Robert Ghrist. **Toward a Spectral Theory of Cellular Sheaves**.
[arXiv:1808.01513v2](https://arxiv.org/html/1808.01513v2), sections 2.2.2 and 3.2:
restriction-valued coboundaries, adjoints and harmonic cochains. The bounded
[sheaf investigation](../research/input-contracts.md#generalized-model-comparison-and-selected-first-investigation)
uses standard real stalk coordinates; its fixture calculations are derived
explicitly and do not establish a general sheaf implementation.

## CCGIKU21

Mathieu Carriere, Frederic Chazal, Marc Glisse, Yuichi Ike, Hariprasad Kannan
and Yuhei Umeda. **Optimizing persistent homology based functions**.
ICML 2021, PMLR 139:1294-1303.
[Published full text](https://proceedings.mlr.press/v139/carriere21a/carriere21a.pdf),
section 3.2 and Propositions 3.2/3.4 for fixed-order regions and differentiability;
section 4.2 for total persistence. Theorem 4.2 has additional Lipschitz, step-size,
boundedness and noise hypotheses. The [feedback experiment](../research/route-integration.md#optional-feedback-and-invalidation-decision-93)
proves its own scalar contraction and makes no claim of general stochastic
convergence or differentiability at ties.

## Architecture references

| Project | Source | Design lesson |
| --- | --- | --- |
| GUDHI C++ | [Pinned source study](../research/gudhi-cpp.md) | Multiple representations and specialized engines; not globally tied to Simplex tree |
| PHAT | [Paper, §2](https://www.geometrie.tugraz.at/kerber/kerber_papers/phat_jsc.pdf) | Separate boundary-column storage from reduction strategy |
| Ripser | [Repository](https://github.com/Ripser/ripser), B21 | Specialized Rips computations can use implicit representations |
| giotto-tda | [API modules](https://giotto-ai.github.io/gtda-docs/latest/modules/index.html) | TDA includes capabilities beyond a single persistence pipeline |

These are design references, not a statement that Cocycle implements their APIs.

## External validation tools

Current comparisons use [pinned GUDHI and upstream Ripser C++ sources](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/native/sources.json)
through [native workers](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/native/README.md). Source versions, adapter
changes and precision are recorded with each run.

[GUDHI](https://gudhi.inria.fr/python/latest/) and
[Ripser.py](https://ripser.scikit-tda.org/en/latest/reference/stubs/ripser.ripser.html)
are optional historical development references. Pin versions and record precision, input,
coefficient field, scale, truncation, and interval conventions for each comparison.

The optional wrapper tools pin GUDHI 3.13.0 and Ripser.py 0.6.14. Ripser.py is a
Python binding/fork, distinct from the upstream C++ CLI. Its installed dense
input path converts distances to float32; the comparison suite uses shared
float32-exact inputs while retaining Cocycle/GUDHI f64 arithmetic. The
[online source](https://ripser.scikit-tda.org/en/latest/_modules/ripser/ripser.html)
helps locate this behavior but may describe a different version. See the
[benchmark protocol](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/README.md) and recorded artifacts for actual runs.

## RP2

Sonia Balagopalan. **Small Triangulations of Projective Spaces**.
[Author's conference slides](https://www.maths.tcd.ie/~hmigca-18/slides/SoniaCGA.pdf).

The displayed six-vertex triangulation has facets 123, 124, 135, 146, 156, 236,
245, 256, 345, 346. Tests subtract one from these labels and take barycentric
subdivision to obtain a flag complex. Independent modular boundary ranks verify
the field-sensitive fixture; the source identifies the triangulation.

## CJS15

Nicholas J. Cavanna, Mahmoodreza Jahanseir and Donald R. Sheehy.
**A Geometric Perspective on Sparse Filtrations**. CCCG 2015.
[Author-hosted full text](https://donsheehy.net/research/cavanna15geometric.pdf),
[arXiv](https://arxiv.org/abs/1506.03797).

Relevant locations: §2 for finite metrics embedded in the max norm, §3 for the
covering lemma, §4 Theorem 4 for sparse nerve approximation, and §5 for Rips
simplex birth and disappearance constraints. The parameter conversion and
edge-length convention used by Cocycle are recorded in the
[specification](mathematics.md#15-sparse-rips-approximation).
[GUDHI's C++ documentation](https://gudhi.inria.fr/doc/latest/group__rips__complex.html)
states its `(1, 1/(1-epsilon))` convention; its pinned
[Sparse_rips_complex.h](https://github.com/GUDHI/gudhi-devel/blob/cba915e3ab8e1f5b1fe26eb44b407285f7af4e78/src/Rips_complex/include/gudhi/Sparse_rips_complex.h)
defines the compared edge/blocker algorithm. Theorem statements do not certify
floating-point implementations; native agreement and mathematical guarantees
are separate forms of evidence.

## HO

[homology-operator](https://github.com/proffitteoy/homology-operator/tree/a11bf50b8c8320651e9ac015566778d957af9c17),
MIT source revision `a11bf50b8c8320651e9ac015566778d957af9c17`. The feasible
generalized-inverse construction and deterministic packed native route are the
implementation reference. The [specification](mathematics.md#21-feasible-f2-homology-projections)
gives the algebraic proof and finite-stage convention used in Cocycle. Source
agreement is not a global-optimality proof. Bounded affine geometry queries and
historical feature correction follow this reference's contracts; Python result
schemas and globally minimum-stretch solvers are not ported. This is a source
citation, not a claimed peer-reviewed algorithm publication or a Rust runtime
dependency.

## OB23

Stable-volume Definition 6 and the separate conditional perturbation results in
Theorems 3.1/4.1: [published source](https://doi.org/10.1007/s41468-023-00119-8).
The pinned [HO](#ho) implementation specifies its F2 constraints and finite replay
domain. Cocycle's `stable_volume` implements the bounded integer optimization
contract; it does not certify the geometric hypotheses of those theorems or
claim a perturbation/resampling guarantee.
