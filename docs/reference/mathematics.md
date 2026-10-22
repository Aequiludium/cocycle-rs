# Mathematical specification

[Documentation](../README.md) / Reference

This document defines the mathematical and numerical contracts shared by the
implementation, rustdoc, and tests. It distinguishes mathematical results from
Cocycle's engineering conventions. See [references](bibliography.md) for primary
sources and [testing](../development/testing.md) for independent checks.

## 1. Scope and notation

The current algorithms compute ordinary Vietoris–Rips persistent homology of
finite inputs over a validated prime field $\mathbb F_p$ through any requested
homology dimension, with F2 the default. The legacy `RipsOptions` path remains
limited to F2 H0/H1. Optional cycle and cocycle bases are described in section 14.
Reduced homology, non-prime coefficient rings, zigzag and multiparameter
persistence are not implemented. Nonnegative scales and zero H0 births are Rips-specific; the general
interval representation is not restricted to these dimensions or scales.
Supplied `SimplicialComplex` and `FilteredComplex` analysis also supports signed
filtrations and unequal vertex births; see [the guide](../guides/filtered-complexes.md).

Matching distances between complete diagrams are specified in section 16.

| Symbol | Meaning |
| --- | --- |
| $n,d$ | Vertex count and ambient dimension |
| $\delta_{ij}$ | Finite, nonnegative, symmetric dissimilarity; Euclidean distance for point clouds |
| $t,T$ | Filtration scale and optional cutoff |
| $q$ | Maximum requested homology dimension |
| $\sigma,f(\sigma)$ | Simplex and filtration value |
| $D,R$ | Boundary matrix and reduced matrix |
| $\mathcal D_k$ | Persistence diagram in dimension k, as a multiset |

## 2. Inputs and filtration

### 2.1 Distances

For $X=(x_0,\ldots,x_{n-1})$, with $x_i\in\mathbb R^d$,

$$\delta_{ij}=\lVert x_i-x_j\rVert_2.$$

Vertices are labeled observations. Duplicate coordinates remain distinct
vertices, so zero distance between different vertices is valid (a pseudometric).

**Engineering conventions.** Require $n\ge0$, $d>0$, row-major coordinates of
length $nd$, and finite `f64` coordinates. Dissimilarities must be finite and
nonnegative. Interpret `-0.0` as zero without modifying borrowed input; normalize
stored/output zero to positive zero. Reject unrepresentable distances instead of
substituting infinity. Scaled `hypot` avoids unnecessary square-sum overflow and
underflow. Reduction is combinatorially exact for the computed floating-point
filtration values, not exact real arithmetic.

Precomputed input is the strict lower triangle with an explicit vertex count:

$$
[\delta_{10},\delta_{20},\delta_{21},\ldots],\qquad
\operatorname{offset}(i,j)=i(i-1)/2+j\quad(i>j).
$$

The length is $n(n-1)/2$, the diagonal is implicitly zero, and size arithmetic is
checked. Zero and one vertex both need an empty buffer, distinguished by $n$.
The triangle inequality is not required: symmetric dissimilarities still define
a flag filtration, but their validation does not certify a metric.

### 2.2 Rips filtration

A simplex is a nonempty finite vertex set, with dimension $|\sigma|-1$. A
simplicial complex is closed under nonempty faces. A filtration satisfies
$\tau\subseteq\sigma\Rightarrow f(\tau)\le f(\sigma)$.

Rips and supplied flag filtrations use edge-length scales:

$$
f(\{i\})=0,\qquad
f(\sigma)=\max_{i,j\in\sigma}\delta_{ij}\quad(|\sigma|\ge2),
$$

$$
K_t=\{\sigma:f(\sigma)\le t\}\quad(t\ge0),\qquad
K_t=\varnothing\quad(t<0).
$$

The threshold is closed. The parameter is neither a ball radius nor a squared
distance; compare [B21 §2](bibliography.md#b21). A face uses a subset of the vertex
pairs, so its maximum edge cannot be larger. This proves face closure and
$K_s\subseteq K_t$ for $s\le t$ without requiring the triangle inequality.

The explicit reference sorts by `(value, dimension, sorted_vertex_ids)` ascending.
The implicit path uses the descending combinatorial IDs in section 9. Both put
faces before cofaces. Tie refinements affect internal pair IDs but must preserve
the positive-lifetime diagram multiset. Sorting never merges nearby values using
an epsilon.

## 3. Chains, boundaries, and homology

First consider the F2 specialization; section 13 gives oriented prime-field
coefficients. Let $C_k(K;\mathbb F_2)$ have the k-simplices as a basis. A chain can be represented
by a set of simplex indices, with addition given by symmetric difference.

$$
\partial_k[v_0,\ldots,v_k]
=\sum_{r=0}^{k}[v_0,\ldots,\widehat v_r,\ldots,v_k]
\quad(k\ge1),\qquad \partial_0=0.
$$

Each codimension-two face occurs twice when deleting two vertices, and cancels
over F2. Hence $\partial_{k-1}\partial_k=0$, and

$$
Z_k=\ker\partial_k,\quad B_k=\operatorname{im}\partial_{k+1},\quad H_k=Z_k/B_k,
$$

$$
\beta_k=\dim C_k-\operatorname{rank}\partial_k-
\operatorname{rank}\partial_{k+1}.
$$

This rank formula supplies an oracle independent of persistence pairing.
Computing H1 needs triangles: edges determine $\ker\partial_1$, but triangle
boundaries determine which cycles are filled. A $(q+1)$-skeleton suffices for
H0 through Hq; its unpaired higher-dimensional chains are not full Rips output.
See [ZC05](bibliography.md#zc05) for the algebraic background.

## 4. Persistence and reference boundary reduction

Inclusions $K_s\hookrightarrow K_t$ induce homology maps. A finite filtration
over a field admits an interval-multiset description. Finite intervals use
$[b,d)$: present at birth, absent at death.

For ordered simplices $\sigma_0,\ldots,\sigma_{m-1}$, set $D_{ij}=1$ exactly when
$\sigma_i$ is a codimension-one face of $\sigma_j$. Nonzero entries satisfy
$i<j$ and $D^2=0$. `low` is the largest row index of a nonempty column; an empty
column has no pivot. Row zero is valid and must not be an empty sentinel.

```text
pivot_owner = empty map
R = empty columns
for j in increasing filtration order:
    column = boundary(sigma[j])
    while column is nonempty:
        i = low(column)
        if pivot_owner has no i:
            break
        column = xor(column, R[pivot_owner[i]])
    R[j] = column
    if column is nonempty:
        i = low(column)
        pivot_owner[i] = j
        record pair (i, j)

after all columns:
    unpaired_births = {i: R[i] is empty and i is not in pivot_owner}
```

Only earlier columns are added, so $R=DV$ for an invertible upper-triangular V,
preserving every prefix boundary-column space. Each cancellation strictly lowers
the pivot, ensuring termination; nonzero reduced columns have distinct pivots.
The pairing theorem then identifies `(low(R[j]), j)` with persistence pairs and
unpaired empty columns with surviving births. Equal rank alone is not a proof of
this pairing theorem: the direct source is [B21 Proposition 3.1 and Algorithm 1](bibliography.md#b21).

A pair `(i,j)` gives dimension $\dim\sigma_i$ and endpoints $f(\sigma_i),f(\sigma_j)$.
Determine unpaired births only after reduction, and output only dimensions through q.
Internal validation retains zero-length index pairs; public diagrams omit $b=d$
because $[b,b)$ is empty. There is no additional minimum-persistence filter.
The independent F2 test oracle described here does not export representatives;
section 14 specifies the production representative path.

Sparse columns use ordered symmetric difference with reusable buffers. An empty
reduced column retains its index marker while its allocation may be reused: an
empty column cannot own a pivot and will not be read for later elimination.
This preserves order, pairs, and the $R=DV$ invariant.

## 5. Complete and truncated computations

Let $\Delta=\max_{i<j}\delta_{ij}$, with $\Delta=0$ if there are no vertex pairs.
At $t\ge\Delta$, the full flag complex of a nonempty input is one full simplex.
Its ordinary persistence has one essential H0 interval and no essential H1.
This is a property of these Rips inputs, not of arbitrary complexes.

`max_edge=None` requests the complete filtration. A finite $T\ge0$ includes all
necessary simplices with $f(\sigma)\le T$; $T\ge\Delta$ also establishes complete
coverage. A computation failure is an error, never a successful truncation.

| Endpoint | Interpretation |
| --- | --- |
| `Finite(d)` | Observed death; interval $[b,d)$ |
| `Essential` | Unpaired in a complete computation; death $+\infty$ |
| `RightCensored { through: T }` | Alive at T; full death is unknown |

In incomplete coverage, all surviving births, including H0, are conservatively
right-censored. The result retains computed dimensions and coverage. A single
optional death value cannot express these distinctions.

Computed dimensions form a nonempty declared set, independently of which interval
lists are nonempty. An omitted dimension is unknown, not a claim of zero homology.
`PersistenceDiagram::new(q, ...)` and current builders declare all dimensions
from zero through q; `with_dimensions` also represents noncontiguous domains.
Dimension-specific analysis requires membership in that set. The greatest
computed dimension does not establish membership of smaller dimensions.

One could extend the truncated filtration constantly and write infinite intervals
for its survivors. Censoring instead describes uncertainty about the original
full filtration. Surviving through T is not dying at T.

## 6. H0 by union-find

All vertices are born at zero. Process edges by `(weight, vertex_ids)` and merge
connected components using a deterministic representative rule. An edge joining
two components kills one H0 class at its weight; an edge within a component does
not change H0. Higher-dimensional simplices do not merge additional components.

For complete input, selected edges form a minimum spanning tree; repeated points
allow zero-weight edges. At a cutoff, they form a minimum spanning forest with
surviving components. Test this path against boundary reduction. See
[B21 §4](bibliography.md#b21) for the H0 algorithm context.

## 7. Diagram descriptors

For one dimension, select finite positive lifetimes $\ell_i=d_i-b_i>0$. Let their
count be N and total $L=\sum_i\ell_i$.

| Descriptor | Definition | No finite positive lifetimes |
| --- | --- | --- |
| Finite count | $N$ | 0 |
| Total persistence | $L$, power one without a root | 0 |
| Maximum persistence | $\max_i\ell_i$ | `None` |
| Persistence entropy | $-\sum_i p_i\ln p_i$, $p_i=\ell_i/L$ | `None` |

[A20 Definition 3.1](bibliography.md#a20) uses base-two logarithms. Cocycle explicitly
uses natural logarithms (nats), without division by $\ln N$. One positive interval
has entropy zero; an empty collection has no normalized lifetime distribution.

The summary reports excluded essential and censored counts. It describes observed
finite lifetimes and does not estimate censored deaths or replace them with T.
Non-finite arithmetic results are errors. Total persistence uses compensated
summation. If an extreme probability underflows to zero, its $p\ln p$ contribution
is taken as zero, its limiting value.

Betti curves use all intervals:

$$\beta_k(t)=\#\{[b,d)\in\mathcal D_k:b\le t<d\}.$$

Essential intervals count for $t\ge b$. Censored intervals count throughout the
known range $b\le t\le T$; truncated diagrams reject queries beyond T. Caller-supplied
grid values must be finite and strictly increasing; negative values are valid. Omitting
zero-length pairs does not change Betti numbers at any scale.

## 8. Stability assumptions

For finite metric spaces, Rips diagrams satisfy

$$d_B(\mathcal D_k(X),\mathcal D_k(Y))\le2d_{GH}(X,Y),$$

by [CSO13 Theorem 5.2](bibliography.md#cso13). Bottleneck distance $d_B$ allows
matching to the diagonal; $d_{GH}$ is Gromov–Hausdorff distance. Do not apply this
statement directly to other filtrations, arbitrary nonmetric inputs, or censored
diagrams without establishing the required assumptions.

**Fixed-vertex derivation used in tests.** If
$\max_{ij}|\delta_{ij}-\delta'_{ij}|\le\varepsilon$, each simplex's maximum edge
changes by at most $\varepsilon$. The complete filtrations include into each
other's $\varepsilon$ shifts, giving an $\varepsilon$-interleaving. Applying
[CSO13 Theorem 2.3](bibliography.md#cso13) yields $d_B\le\varepsilon$.
This argument also holds for finite symmetric dissimilarities on the same vertices.

If corresponding Euclidean points move by at most $\eta$, the triangle inequality
bounds each edge change by $2\eta$, giving the corresponding $2\eta$ bound.
Numerical tests additionally allow explicit rounding tolerance. This does not
imply continuity of thresholded interval counts.

## 9. Implicit Rips persistent cohomology

This section specifies the production H1 path and exact F2 H2+ handoff. The explicit boundary oracle
remains independent. Cohomology, clearing, implicit columns, and shortcut pairs
have separate conditions; see [B21 §§3.2–3.5 and §4](bibliography.md#b21).

### Numbering and order

For $a<b<c$,

$$
\operatorname{id}(a,b)=\binom b2+a,\qquad
\operatorname{id}(a,b,c)=\binom c3+\binom b2+a.
$$

Forward filtration order is maximum edge ascending, dimension ascending, and
combinatorial ID descending. Faces precede cofaces. Use descending comparison,
not unsigned negation. Exhaustive small-set tests check IDs and inversion.
Changed tie order may change index pairs but not the positive-lifetime multiset.

Cancel binomial denominators before checked multiplication. H1 requires
$\binom n3$ to fit the platform's `usize`, even when a cutoff retains few triangles.
Return `SizeOverflow` otherwise; 32-bit and 64-bit indexing limits differ.

### Coboundaries and pairing direction

The F2 coboundary of an edge consists of triangles containing it within the
cutoff. A transpose alone is insufficient: the small matrix oracle uses
$C=JD^\top J$, reversing both rows and columns. With S simplices, a C pair `(i,j)`
corresponds to the D pair `(S-1-j,S-1-i)`. Restore original dimensions and
filtration values, and also check unpaired indices and censored H1 classes.

Scan edges forward to obtain union-find merges and cycle births. Then process
edge coboundary columns in reverse, with the earliest forward triangle as pivot
(the highest reversed row). A paired edge/triangle yields an H1 birth/death;
an unpaired cycle-birth edge survives. H0 comes from union-find independently.
Clearing skips H0 merge edges: the adjacent-dimension pairing guarantees their
coboundaries reduce to zero. Tests disabling clearing still reduce them, verify
that result, and never mistake them for H1 births. Collect edges even after the
graph becomes connected.

### H1-to-H2 clearing handoff

For exact diagram-only F2 H2+ requests, specialized H0/H1 supplies owned intervals
and every H1 death triangle to the generic continuation starting at dimension 2.
Ordinary and emergent pivots and omitted zero apparent pairs all belong to this
clearing set; neither stored owners nor positive-lifetime intervals alone are
complete. Triangle keys are decoded to original vertex triples before H1 reducer
storage is released. The continuation enumerates complete topology, including
cleared triangles, without rerunning generic H0 or H1. An equal-weight K3 has
empty H2; failing to clear its virtual death triangle would invent an H2 birth.
Odd-prime, explicit/blocked sources and representatives retain their established
paths. If H1 combinatorial indexing overflows on an H2+ request, use the full
tuple-based generic path. One WorkBudget includes handoff extraction, conversion,
topology generation and continuation; allocation or interruption returns no result.

### Experimental fixed-tuple H2 baseline

The private T3 prototype is available in tests and builds explicitly using
`--cfg cocycle_h2_bench`. Normal builds keep the generic H2 continuation. The
benchmark configuration selects only exact F2 diagram-only H2 requests; H3+,
odd primes, representatives and blocked approximate sources keep their routes.
It has no apparent, emergent or virtual H2 shortcuts.

Triangles and tetrahedra use `[usize; 3]` and `[usize; 4]`, without narrowed
vertices or tetrahedron binomial indexing. The shared total order compares
maximum edge value, then decreasing colex vertices. Enumerate the complete
triangle level, including H1 death triangles; clearing skips their reduction
only. Sparse tetrahedra use a three-way sorted adjacency intersection. Dense
visits test all candidate vertices against the three triangle vertices. Both
charge failed candidates/comparisons to the same WorkBudget.

Process triangles in reverse order and pivot on the earliest forward tetrahedron.
Store every pivot owner, including zero-length pairs, and a parity-normalized
list of triangle positions for V, including its diagonal. Reconstruct C V on
each owner addition and cancel both rows and transformation positions by parity.
Source positions in V precede the active column in reverse computation order;
each elimination advances the pivot in forward order. Independent explicit
boundary matrices and all-edge tetrahedron enumeration check intervals, order,
triangularity and R = C V. This is a finite validation baseline for further
experiments, not production admission or an arbitrary-dimensional proof.

Test-only ownership landmarks record the H1 return with already-owned death
tuples, clearing-set extraction, triangle-level construction after edge release,
and H2 reduction start. H1 death extraction occurs during H1 pairing, before its
return; it is not a separate late decoding pass. Vector capacities and elapsed
offsets describe those objects, not total live allocations or process RSS. The
prototype deliberately retains its triangle clearing set and allocates fresh
per-column scratch; lifetime and capacity-reuse experiments belong to M1.

### Implicit reconstruction invariant

Maintain $R=CV$. For each pivot, store
$V_j=e_j+\sum_{k\in A_j}e_k$ rather than retaining $R_j$. To cancel a pivot,
regenerate $Ce_j$ and all $Ce_k$, adding them to the working column while adding
their edge indices to the current transformation. The edges in $A_j$ precede j
in reverse computation order, so V stays unit upper triangular.

Repeated triangles and transformation indices cancel by F2 parity. Cancel all
copies of the candidate pivot before using it; ordinary set deduplication is
incorrect. Each elimination strictly lowers the pivot in reversed row order.
Stored transformations are internal machinery, not public representative cocycles.

### Clearing and shortcut pairs

H0 pairings and H1 clearing must use the same total order. An apparent pair
$(\sigma,\tau)$ requires sigma to be tau's latest facet and tau to be sigma's
earliest cofacet. Omitting its public interval also requires equal filtration
values. Do not indiscriminately remove columns that eventually reduce to zero.
The kernel omits stored zero-lifetime apparent pairs and reconstructs them when
later columns need their pivots. Other shortcut pairs retain their pivot owners
and transformation columns.

Initialize the column and inspect shortcut candidates in one traversal, before
any column addition. Cofacets are enumerated by descending ID,
and their values are at least the edge value. The first equal-valued triangle is
therefore the original column's earliest cofacet. Retain preceding cofacets in a
temporary buffer. If no shortcut applies, continue the same traversal and
build the working heap from the complete buffer; do not enumerate the prefix or
an empty column again. Reuse the buffer only after its previous contents have
been consumed or discarded.

An emergent shortcut accepts this first equal-valued pivot only when it has no
stored or virtual owner. Ordinary reduction would accept it immediately; retain
$V_j=e_j$ without generating remaining rows. If owned, fall back to ordinary
reduction, not the next equal-valued triangle. Apparent-only mode additionally
checks the latest-facet condition. Mid-reduction emergent shortcuts are not used.

For a zero-lifetime apparent pair $(\sigma,\tau)$, the same initial traversal
certifies both conditions: tau is sigma's first equal-valued cofacet, and sigma
is tau's latest facet. Skip this pair without storing its pivot owner or
transformation column. It remains a virtual column $V_\sigma=e_\sigma$;
omitting its storage does not remove its role in later elimination.

When a working pivot has no stored owner, test for that virtual pair on demand:
take the triangle's latest facet sigma, require equal filtration values, and
check that the triangle is sigma's earliest equal-valued cofacet. Only an edge
preceding the active column in reverse computation order may serve as its
virtual owner. Add $Ce_\sigma$ to cancel the pivot and include sigma in the
active transformation, with F2 parity as for an ordinary stored owner. This
preserves $R=CV$ and triangularity; every elimination strictly lowers the pivot.
A triangle that fails these conditions cannot be treated as a virtual owner.

Zero-lifetime H1 pairs are omitted before raw interval storage; the shared result
assembler still removes zero-lifetime H0 pairs. Union-find, clearing, and the
stored or virtual reduction state still account for those pairs. In particular,
zero-lifetime emergent pairs that are
not apparent retain their pivot/transform entries. Positive-lifetime, essential,
and right-censored intervals retain their multiplicities. Test ordinary,
apparent, emergent, and combined configurations against the independent oracle.
Independent two-pass tests also check storage omission independently of
single-pass caching, including the invariant $R=CV$ on dense and sparse inputs.

### Cone stopping bound

For a nonempty finite symmetric nonnegative input, define

$$r_*=\min_v\max_u\delta(v,u).$$

Choose a minimizing vertex $v_*$. At $t\ge r_*$, all edges to $v_*$ exist.
For every Rips simplex sigma, its union with $v_*$ is also a simplex. Thus the
complex is a cone: one ordinary H0 class and no positive-dimensional homology.
This derivation uses only the maximum-edge rule, not the triangle inequality.

A complete computation may stop at $r_*$; with cutoff T, use $\min(T,r_*)$.
Include all events exactly at the stop. Handle empty input separately and set
$r_*=0$ for a singleton.

Internal stopping does not replace public coverage. No user cutoff still means
`Complete`. If $T<\Delta$, retain `Through(T)` and its censoring convention even
when the cone is reached earlier; never substitute $r_*$ for T.

## 10. Analytic complete-bipartite filtration

Partition n vertices into groups of sizes $a,b>0$, with $a+b=n$. Assign cross-group
distance 1, within-group off-diagonal distance 2, and diagonal zero.

- For $0\le t<1$: vertices only; $\beta_0=n$, $\beta_1=0$.
- For $1\le t<2$: the connected triangle-free graph $K_{a,b}$;
  $\beta_1=E-V+1=ab-a-b+1=(a-1)(b-1)$.
- For $t\ge2$: the full simplex; $\beta_0=1$, $\beta_1=0$.

The complete diagram contains n−1 H0 intervals $[0,1)$ and one $[0,\infty)$;
H1 contains $(a-1)(b-1)$ copies of $[1,2)$. At T=1, H1 survivors and the remaining
H0 class are censored. Exception: a=b=1 has actual diameter 1, so T=1 is complete,
the H0 survivor is essential, and H1 is empty.

This independent derivation checks repeated endpoints, quadratic output
multiplicity, and triangle-free enumeration. The output itself can contain
$\Theta(n^2)$ intervals; do not attribute all their storage to reduction fill-in
or deduplicate intervals to improve measured cost.

## 11. Exact graph inputs and shared flag computation

The [construction guide](../guides/rips-construction.md) describes the public APIs.
A finite graph G has vertices at zero, edge weights w >= 0, and a simplex for
each clique, born at its maximum edge weight. Missing edges never enter. This
complete flag filtration can have essential H1. A cutoff below its maximum edge
conservatively gives right-censored survivors under the existing diagram contract.

An exact threshold graph constructed from a full pairwise input instead certifies
only the original Rips through its recorded cutoff, unless every original pair
was included. Its maximum retained edge cannot establish original completeness.
Requesting computation beyond that certified range is an error. With no new
cutoff, computation uses the construction's available range.

The shared H1 engine requires forward edge order and decreasing combinatorial ID
order for triangle cofacets. For a fixed edge, inserting increasing third vertices
increases the triangle ID. Reverse intersection of sorted neighbor lists therefore
satisfies the dense enumerator's order. Every emitted triangle has all three edges,
and latest-facet lookup uses the same value/ID order. These facts preserve the
original-column apparent/emergent-pair tests on supplied flag filtrations. No
complete-distance cone bound is inferred from missing sparse edges. Independent
vertex-triple enumeration and ordinary reduction validate the sparse path.

Lower/upper/square layouts borrow original f64 values. Exact square symmetry and
zero diagonal are checked without tolerance repairs. Exact construction needs no
triangle inequality. Custom callbacks supply one orientation of each pair and
must respect the documented symmetric, stable-function contract. Work limits and
cancellation indicate execution failure; they do not change a filtration cutoff
or return a successful partial diagram.

## 12. Explicit skeletons and dimension-generic cohomology

For a weighted graph, each stored clique has increasing vertices and value equal
to its largest edge (vertices enter at zero). Production filtration order is
increasing value, increasing dimension, then decreasing colexicographic order.
For equal-sized simplices this last order agrees with decreasing combinatorial
rank from section 9, without requiring that rank to fit a machine integer.
An increasing-vertex simplex has oriented boundary obtained by omitting vertex
position i with coefficient (-1)^i. Each codimension-two face appears twice with
opposite signs, so the stored integer boundary squares to zero before reduction
modulo two.

An explicit p-skeleton suffices for Hq only when p >= q+1, or when construction
certifies that it already contains the whole threshold clique complex. Otherwise
missing cofaces can create artificial top-dimensional survivors. The expansion
probes unique extensions of its last dimension to establish exhaustion. This
certificate concerns only dimensions at the constructed scale; it does not extend
original-input scale coverage. Dimension and scale checks are independent.

For dimension-generic F2 computation, classify edges by forward union-find, then
process remaining edges in reverse order. At each dimension, the working column
is a set of cofacets; symmetric difference performs F2 addition. Its earliest
forward cofacet is the reverse coboundary pivot. An owned pivot triggers addition
of the owner's transformation column, regenerating its coboundary on demand. A
new pivot pairs the current simplex with that cofacet; an empty reduced column
records an unpaired class unless the simplex was cleared. Pivot simplices become
cleared columns in the next dimension. Cleared simplices still participate in
clique generation: clearing skips algebra, not topology. This is the dimensional
extension of the reversed-transpose duality in section 9 and [B21](bibliography.md#b21).

Only the current dimension, transformations, pivot ownership and a working
coboundary are needed on the implicit path. These can still be exponentially
large. Legacy explicit Rips entry points adapt stored incidence to this reducer.
Diagram-only explicit simplicial analysis also uses it when all vertices are
born at zero and the query cutoff is nonnegative (or absent). Stored cofaces,
not graph reconstruction, preserve non-flag topology and delayed simplex values.
Other explicit inputs use the filtered-cell boundary reduction in section 17.
Both explicit paths also pay
for the already materialized skeleton. H0/H1-only implicit requests
retain the specialized compact-index engine and its pair shortcuts.

As an independent high-dimensional fixture, partition 2r vertices into r pairs.
Give opposite partners distance 2 and all other pairs distance 1. Between scales
1 and 2 the clique complex is the join of r copies of S0, hence S^(r-1); at scale
2 it becomes a full simplex. Its positive-dimensional diagram consists of one
H_(r-1) interval [1,2) for r >= 2. Omitting opposite edges permanently instead
makes that sphere essential in the supplied graph filtration. Tests exercise
r=3,4,5, alongside independent forward boundary reduction of all small vertex
subsets and native C++ comparisons.

## 13. Prime fields and oriented reduction

`PrimeField` accepts exactly prime u32 characteristics. Arithmetic outputs residues
in 0..p; nonzero sparse entries lie in 1..p. The largest accepted characteristic
is 4294967291. Exact trial division validates primality, u64 intermediates cover
every u32 product and sum, and modular exponentiation computes inverses. Zero
has no multiplicative inverse and is rejected. These arithmetic choices are
implementation contracts, independent of filtration precision.

For increasingly oriented vertices, the coefficient on the facet omitting position
i is (-1)^i modulo p. Coboundary incidence uses the same coefficient. Normalizing
a pivot c multiplies its column and transformation by c^(-1); eliminating a
coefficient a against a normalized pivot adds -a times that column. F2 reduces
to parity cancellation. Dimension progression and clearing do not depend on p.
H0 connectivity uses union-find over every supported field.

Independent dense rank checks use separate arithmetic and enumeration. A
barycentric subdivision of the [six-vertex RP2 triangulation](bibliography.md#rp2)
is flag: vertices represent nonempty faces, and edges join comparable faces.
Pairwise comparable faces form chains, so its clique complex is the subdivision.
This fixture has one H1 and one H2 class over F2 and neither over odd primes.
The test independently checks those ranks, and native C++ comparisons check the
same weighted flag filtration with matching coefficient fields. Reference modulus
limits are recorded as exclusions, never silent narrowing.

## 14. Persistent representative bases

For a full forward boundary reduction R=DV in filtration order, a finite pair
(i,j) supplies cycle R_j with leading simplex i. It is available at birth i and
becomes a boundary at death j. An unpaired birth i supplies V_i. The active
subset of these cycles gives a homology basis at every scale; see
[B21, section 3.1](bibliography.md#b21). Using V_i for every finite pair would
not establish its death association. Production normalizes nonzero pivot columns,
which only rescales the selected basis vectors.

The caller requests representatives after all simplices at a finite query scale
t have entered. Active intervals satisfy birth <= t < death, with an unpaired
interval alive through an inclusive known cutoff. Zero-length pairs are omitted.
At t, reduced boundary columns whose source simplices have entered span B_q(K_t).
For active cycle basis z_i, solve cochain constraints

$$
\varphi_i(b)=0\quad(b\in B_q(K_t)),\qquad
\varphi_i(z_j)=\delta_{ij}.
$$

Sparse row elimination normalizes pivot constraints and then back-substitutes
with free coordinates set to zero. Annihilating boundaries is exactly the
cocycle equation; identity pairing proves nontriviality and independence modulo
coboundaries. This construction defines a dual basis for the requested scale,
not a promise that the returned cochain entries remain unchanged at another
scale. No shortest-support or canonical-across-algorithms claim is made.

Output identities are ordinals in one result's canonical logical interval sequence,
resolved through `PersistenceDiagram::interval(index)` independently of storage layout.
Repeated intervals retain distinct indices; request positions distinguish repeated
queries. Terms contain original, increasing simplex vertices and canonical
nonzero coefficients. Queries outside known coverage or computed dimensions fail.
Requested representatives own their terms and do not borrow a reduction workspace.

Only nonempty request slices materialize a representative skeleton and retain
forward transformations. All relevant exact paths produce the same diagram with
or without requests, while resource costs differ. Tests independently verify
closure, boundary/coboundary quotient ranks, identity pairings, birth support,
finite-death association, coefficient ranges, repeated intervals and censoring.
Native references validate prime-field interval multisets; these adapters do not
claim an upstream representative-vector oracle.

## 15. Sparse Rips approximation

Sparse Rips is a different filtration from exact threshold Rips. The construction
follows the pinned GUDHI C++ sparse edge/blocker convention, with deterministic
sampling instead of its random initial landmark. See [CJS15](bibliography.md#cjs15)
and the [source study](../design/rips.md#three-project-baseline).

Let `p_0,...,p_(n-1)` be a greedy farthest-point permutation. Its insertion radii
are `lambda_0 = infinity` and `lambda_j = min_(i<j) d(p_i,p_j)`. Ties choose the
smallest original vertex ID. Retain the initial vertex and the prefix of positive
radii at least the optional minimum insertion radius. A pseudometric may have
zero-distance duplicates; identifying these preserves its Rips persistence.

For retained `i < j`, write `d = d(p_i,p_j)`, `li = lambda_i`, `lj = lambda_j`.
The modified edge value, in edge-length units, is:

```text
alpha = d                    if d * epsilon <= 2 * lj
edge absent                  if d * epsilon > li + lj (and the first case failed)
alpha = 2 * (d - lj/epsilon)  otherwise
```

For `epsilon < 1`, the last branch is additionally absent if
`alpha * c > lj`, where `c = epsilon * (1-epsilon)/2`. Retain only values at most
the construction threshold. Vertices have value zero. A higher simplex is a
clique whose value `f` is its largest edge value and which satisfies
`lambda_v >= f*c` at every vertex. For epsilon at least one this additional
constraint is disabled. This is not ordinary flag expansion for epsilon below one.
The constraint is hereditary: a face has a no-larger value and no-smaller minimum
insertion radius. Therefore rejecting a simplex cannot remove an admissible coface.

The conditional metric theorem uses `0 < epsilon < 1`, exact arithmetic and a
greedy permutation. GUDHI states a `(1, 1/(1-epsilon))` interleaving with exact
Rips. In [CJS15](bibliography.md#cjs15), set the paper's parameter to
`delta = epsilon/(1-epsilon)`: the ball cap is `lambda/epsilon`, disappearance
radius is `lambda/[epsilon*(1-epsilon)]`, and the factor is `1+delta`.
GUDHI doubles radius scales to edge-length scales. CJS15 §2 embeds finite metrics
isometrically in the max norm, allowing its nerve argument to cover metric Rips;
§4 states the approximation theorem and §5 explains the higher-simplex test.
Bounds apply to the retained metric subset after positive-radius subsampling.
No unchanged multiplicative guarantee to the full input is claimed then.

The library records checked, caller-assumed or absent metric hypotheses. Exact
checking treats finite binary64 inputs as real dyadic numbers. FastTwoSum
recovers the sign of a rounded sum's error to test triangle equality correctly;
no tolerance silently changes the hypothesis. A pseudometric with zero-distance
pairs is allowed. The reported factor is a nominal binary64 evaluation of the
ideal theorem, not a certificate for accumulated construction rounding. Extreme
arithmetic that overflows, or epsilon causing a zero blocker factor by underflow,
returns a numerical error. Other operations use ordinary binary64 rounding.

Coverage describes the approximate filtration. Edges excluded by its rules never
enter; those omitted solely by a finite threshold make the range incomplete.
Dimension truncation is separate: Hq requires construction through q+1 unless
expansion certified exhaustion. All public simplex/representative vertex labels
refer to original input IDs, while the exposed graph uses a documented compact
map. Representative bases describe this approximate filtration and do not claim
a chain map into original Rips at the same scale.

## 16. Diagram matching distances

The `diagram_distances` domain compares one explicitly computed homology dimension
of two complete diagrams. Points are a multiset: repetitions remain separate
matching obligations. For finite points, allow partial bijections between the
two multisets and match every unused point to the diagonal.

| Operation | Ground cost between finite points | Diagonal cost | Objective |
| --- | --- | --- | --- |
| `bottleneck_distance` | L-infinity | `(death-birth)/2` | Minimum largest cost |
| `wasserstein_1_infinity` | L-infinity | `(death-birth)/2` | Minimum sum of costs |
| `wasserstein_2_euclidean` | Euclidean | `(death-birth)/sqrt(2)` | Square root of minimum sum of squared costs |

Only `(finite birth, positive infinity)` essential points are representable.
They match essential points in sorted birth order, with absolute birth difference
as cost. Unequal counts give positive infinity. Combine finite and essential
contributions by maximum, sum, or Euclidean norm for the three operations.
Empty computed dimensions are valid; uncomputed dimensions are errors. Reject
all `Coverage::Through` inputs, even if no current interval is censored: future
births and deaths are unknown. A cutoff must not replace a death or certify
essentiality. Diagonal points, non-finite births and other infinite endpoint
categories remain excluded by the existing interval constructors.

The `_results` functions additionally require equal coefficient characteristics
and declared edge-length scales on both inputs. An unspecified scale is rejected,
even when both sources have unspecified scales; matching enum values do not
establish comparable units. Supplied-complex users can explicitly pass raw diagrams
after establishing their common scale. Even declared edge lengths do not certify
physical units or normalization across datasets; callers must establish them
for both raw-diagram and result-based calls. In sparse Rips, the modified edge
value `2 * (d - lambda / epsilon)` still has edge-length units, but changes the
filtration. Comparing its diagram with exact Rips is a distance between the two
actual diagrams, not a distance between their input metrics or an approximation
error certificate. Vertex counts, requested cutoffs
and filtration kinds need not match when coverage is complete.
Approximate constructions retain their provenance in the borrowed results;
the returned scalar measures their actual diagrams without certifying a distance
between the original datasets. Raw-diagram calls cannot establish provenance.

Exactness excludes algorithmic approximation; it does not mean exact real
arithmetic. Use binary64, deterministic tie handling and stable cost aggregation.
W2 is recomputed from the chosen original cross/diagonal costs, rather than by
subtracting two nearly equal total savings. Unrepresentable required arithmetic
returns `NumericalFailure`, distinct from the mathematical infinity caused by
unequal essential counts. Underflow must not silently erase a required nonzero
cost. Scaling used by a solver must preserve representable input distinctions
or fail explicitly.

The native solvers are adapted from the MIT-licensed
[Topp source](https://github.com/proffitteoy/Topp/tree/ffa1da051ca7ac5e313c74cc9fb92a2bcb20c234).
Implementation-local license notices retain attribution. Their private routing
and benchmark switches are not public API. Independent exhaustive partial
matching, hand-derived cases and Topp/GUDHI comparisons validate the supported
domain; source translation alone is not independent evidence.

## 17. Supplied filtered-cell boundary contract

Generic filtered cells and general supplied simplicial filtrations share the
same boundary-pair reducer. Zero-born supplied simplices retain the coface
reduction described in section 13. The simplicial reader uses constructor-validated
oriented incidence; the generic reader checks its external contract. Caller
source certificates determine coverage, even when the retained boundary columns
alone appear complete. Representative basis extraction is a separate request.

For ordered cells c_i, let D[j,i] be the integer incidence coefficient of c_j in
the boundary of c_i. Nonzero entries require j < i, dimension(c_j) =
dimension(c_i)-1 and f(c_j) <= f(c_i). The mathematical source must satisfy D^2=0
over the integers. Reduction maps coefficients to the selected prime field and
checks this chain identity on the retained q+1 skeleton in that field. This
field-specific check does not prove integer validity of a custom adapter.

Unpaired zero reduced columns represent essential classes of the supplied complex,
or censored classes when an analysis cutoff truncates its filtration. Vertex birth
values are read from the source. Maximum edge value is not used to bound general
filtrations: an Alpha triangle can enter after all its edges. Certified Rips
expansions additionally preserve original scale coverage and skeleton sufficiency;
a bare supplied complex does not assert that relationship to a larger source.

## 18. Lower-star construction example

For a finite simplicial complex K and a finite scalar function g on its vertices,
define f(sigma) = max {g(v) : v in sigma} for every nonempty simplex. If tau is
a face of sigma, its vertex set is a subset, so f(tau) <= f(sigma). Thus sublevel
sets are face-closed. This lower-star filtration permits negative and tied values;
its units are those of g, not necessarily distance or squared distance.

The [contributor example](../../examples/complex_construction.rs) supplies all
vertices and all nonvertex simplices explicitly. It validates rather than infers
missing nonvertex faces. Its vertex labels are the indices of g, and shared
simplices must occur only once. This is a filtration assignment on given topology,
not a geometric triangulation algorithm or a public production constructor.

For the four-cycle with edges 01, 12, 23, 03 and vertex values (-2, 1, -1, 0),
two components are born at -2 and -1. They merge at 0 through vertex 3. At 1,
vertex 1 and its two edges close the loop. Omitting zero-lifetime pairs gives
H0 intervals [-2, infinity), [-1, 0) and H1 interval [1, infinity) over every
prime field. At cutoff -0.5 both components are right-censored. A filled triangle
with every vertex value 2 enters with all faces at 2 and has only the essential
H0 interval [2, infinity); it has no positive-lifetime H1 interval.

## 19. Big-steps critical sets over F2

Let D be the square boundary matrix in a frozen simplex order, R = DV its
serial lazy reduction, and U = V inverse, over F2. V is unit upper triangular.
Lazy reduction adds a negative column only to cancel the current low. Thus
an off-diagonal V[i,j] can be nonzero only when column i is negative and
low(R_i) > low(R_j), with low(0) = -1. Critical-set formulas require this
support condition, not merely R = DV. See [NM22](bibliography.md#nm22),
Theorems 9, 10 and 11 and Algorithms 1-3.

For a negative endpoint j with current value f_j, moving its death to t > f_j
uses the support of row U[j,*] within f_j <= f_i <= t. Moving it to t < f_j
uses column V[*,j] within t <= f_i <= f_j. Finite birth moves use the same
negative-column rules in the dual decomposition of J D transpose J, where J
reverses order and dual values are -f. Reverse the returned indices to recover
original source IDs. Essential birth decreases instead use primal V columns;
increases use dual V columns. All sets contain only the endpoint's dimension.
Closed value windows include ties. Cocycle defines zero displacement as the
empty proposal; it leaves both data and lazy preparation unchanged.

A requested U row solves u V = e_j transpose. Start with residual e_j and
repeatedly cancel its smallest pivot p with row p of V. Its diagonal is one and
other entries have larger indices, so pivots strictly increase. Record p in u
and stop when f_p exceeds t (or the original dual value falls below t).
Unprocessed larger pivots cannot affect the retained prefix. Only the queried
dimension of V needs transposing; neither the other U rows nor the full inverse
are required. This is the bounded partial-row strategy of
[Oineus](https://github.com/anigmetov/oineus/blob/e52814a1ffb5b8a81e71ff1e93b4c14194673f0f/include/oineus/decomposition.h#L5456).
Parallel or cleared reductions require appropriate ELZ restoration before
critical-set use. The present workspace constructs serial state directly.

For overlapping endpoint proposals, select for each simplex the target t
maximizing abs(t-f), keeping the first proposal at equal distance. This is
maximum displacement, including competing directions, rather than a numerical
maximum target. With these targets frozen, the squared surrogate sum(f-t)^2
has gradient 2(f-t). Chain rules, tie subgradients, face/coface closure and data
updates are caller responsibilities; this heuristic is not a general descent
or convergence theorem. Independent full inverses, pairing-only block swaps
and hand-derived lower-star intervals test these discrete formulas.

## 20. Finite single-parameter stages and inclusions

For a frozen filtered simplicial complex (K, f) and finite scalar t, define
K_t = {sigma in K : f(sigma) <= t}. Face monotonicity makes K_t a subcomplex.
The total filtration order places K_t in a prefix; all tied values enter
together, with faces preceding cofaces. A stage borrows that prefix, keeping
source-local IDs and increasing-vertex orientation.

For s <= t, inclusion I_q sends each oriented q-simplex in K_s to the identical
simplex in K_t with coefficient +1. Thus D_q(t) I_q = I_(q-1) D_q(s) over the
integers and every prime field; successive inclusions compose. Equal numeric
indices from different owners establish no map. Reversing stages or mixing
construction certificates is outside this identity inclusion contract.

General degree-preserving chain maps F_q require explicit bases, coefficient
arithmetic and the identity D_q(B) F_q = F_(q-1) D_q(A). Merely preserving dimensions
or giving an index permutation does not establish this. Such maps and coordinatewise
multi-parameter filtrations are deferred; scalarizing incomparable parameter
vectors defines a different single-parameter problem.

A selected stage limits observation of the original filtration. Its contextual
persistence keeps later classes right-censored whenever original source coverage
is incomplete at t. Certified dimension truncation remains independent of scale
truncation. Explicitly treating K_t as the entire supplied topology changes that
source and can instead make those classes essential. See the
[stage guide](../guides/filtered-complexes.md#read-stages-and-their-inclusions) and
[fixed-source decisions](../research/filtration-stages.md).

## 21. Feasible F2 homology projections

### Window and quotient

Fix three ordered bases over F2 and column maps A:C_k -> C_(k-1), D:C_(k+1) ->
C_k with shapes m-by-n and n-by-p. Require AD=0. Cycles Z=ker(A) and boundaries
B=im(D) then satisfy B subset Z, and H=Z/B has dimension
beta=n-rank(A)-rank(D). Ordinary H0 uses C_(-1)=0, without augmentation.
Empty maps retain their declared shapes. A supplied simplicial adapter reduces
signed integer incidence modulo two; this field check does not prove integer
boundary validity or completeness of omitted cofaces.

A legal representative-selection map P:C_k -> C_k satisfies all four conditions:

```text
P^2 = P,       A P = 0,       P D = 0,
z + P z is in im(D) for every z in ker(A).
```

The first three alone are insufficient: P=0 satisfies them on a triangle ring
but destroys its nonzero H1 class. Idempotence and cycle image are checked on
coordinate generators, boundary killing on D's columns, and class preservation
on a full cycle basis. Linearity extends these finite checks to the whole
spaces. Independent cycle/image reductions use the actual A/D rather than
trusting the solver's generalized inverses. This is feasibility, not optimality.

### Factorized construction and proof

Choose generalized inverses G and U with AGA=A and DUD=D. Define R=I+GA and
B=DU, then P=(I+B)R. Addition is XOR throughout. These relations give:

```text
R^2 = I + GAGA = I + GA = R,      A R = A + AGA = 0,
R z = z for z in ker(A),         R D = D,
B^2 = DUDU = DU = B,             B D = D,
R B = B,                        (I+B)B = 0.
```

Consequently P^2=(I+B)(R+B)R=P, AP=0, and PD=(I+B)D=0. For every cycle z,
Pz=z+Bz, so z+Pz=Bz is a boundary. Thus P kills boundaries and induces the
identity on H, while im(P) lies in Z. Its restriction to Z has kernel B:
if Pz=0, class preservation gives z in B; the converse follows from PD=0.
Hence Z=B direct-sum im(P), rank(P)=beta, and every homology class has exactly
one representative in im(P). The complement depends on G/U and the bases.

The deterministic implementation follows stable elimination of [M|I], including
pivots in the full right block. A different valid generalized inverse can select
a different complement/action; matching Betti numbers alone does not establish
exact reference agreement. L=I+P also satisfies L^2=L and ker(L)=im(P). This
finite-field operator is not a real/complex Hodge Laplacian and supplies no
positive-definite energy, real spectrum or Hodge-theoretic certificate.

### Positive costs and the limit of feasibility

For positive coordinate costs w_i, define mass_w(x)=sum of w_i over the support
of x. The selected class mass is mass_w(Pz), and the selected difference-class
cost is delta([z],[y])=mass_w(P(z+y)). It is well-defined on H because P kills
boundaries. Under exact arithmetic it is a metric for this fixed complement:
positivity and class preservation make delta=0 exactly for equal classes,
symmetry follows over F2, and mass_w(x+y) <= mass_w(x)+mass_w(y) proves the
triangle inequality. This is not the minimum mass over all representatives.
The bounded implementation may return overflow rather than a metric value.

For an explicit counterexample, let A be 0-by-2, D=(1,1) transpose and w=(9,2).
The legal P with both columns (1,0) transpose selects mass 9 for the nonzero
class, even though (0,1) represents that class with mass 2. Its stretch is 9/2;
the legal projection onto the second coordinate has stretch 1. Thus legality
alone cannot certify minimum stretch or shortest-class mass.

The implemented readout is lambda(P,w)=max mass_w(Pz)/mass_w(z) over nonzero
cycles z. It searches all 2^dim(Z)-1 inputs, not only a homology basis. When
Z=0 the maximum is absent; when Z is nonzero but P=0 the computed value is zero.
For beta>0, a nonzero fixed vector of P has ratio 1, so lambda(P,w)>=1.
`stretch_with` evaluates one P, not a minimization over legal projections.
Mass sums use checked u64 addition; comparing two u64/u64 ratios uses exact
u128 cross products. Positive weights make nonzero denominators positive.
The first maximum in the deterministic binary enumeration is retained on ties.

### Chain inclusions and a transport witness

For stages s <= t, let Jprev/Jcur/Jnext identify their three ordered bases.
Require A_t Jcur=Jprev A_s and D_t Jnext=Jcur D_s. Then Jcur maps cycles to
cycles and boundaries to boundaries, giving an induced map H_s -> H_t.
`from_stages` obtains these coordinate inclusions only after checking the actual
same-owner/certificate relation of borrowed simplicial stages. Abstract-window
source labels remain caller assertions and do not authenticate physical owners.

For a source representative z, choose the target representative z_t=P_t Jcur z.
Legality gives Jcur z+z_t=D_t c for some next-degree coordinate vector c.
`boundary_preimage` returns such a c for independently replaying this identity.
A missing preimage is a failed membership statement; an empty preimage can be a
valid witness. Costs change geometry but not the underlying F2 chain-map identity.

In the triangle edge order (20,30),(10,30),(10,20), A has rows (1,1,0), (1,0,1),
(0,1,1). The ring has D:3-by-0 and Z spanned by z=(1,1,1). Its filling has
D_t=z, hence H1 changes from dimension 1 to 0. Costs (2,3,5) give source selected
mass 10 and target mass 0; the unique target face coefficient c=(1) satisfies
D_t c=z. An independent legal source choice has first column z and other columns
zero, while P_t=0. This hand choice is not the claimed deterministic solver
matrix. P_t Jcur and Jcur P_s need not agree as chain matrices; their class-level
transport is the relevant invariant. General linear maps remain outside the
implemented coordinate-inclusion interface.

### Composition and finite stage intervals

Choose deterministic bases E_i for im(P_i). The transport matrix T(s,t) solves
E_t T(s,t)=P_t Jcur E_s in those bases. It composes: the difference between
J(s,t)z and P_t J(s,t)z is a target boundary, which the next chain map sends to
a boundary killed by P_u. Therefore T(s,u)=T(t,u)T(s,t), even when the P_i do
not commute as chain matrices. Adjacent coordinate transports suffice for ranks.

For N observed stages, let r(i,j)=rank T(i,j) for i<=j<N. The finite-dimensional
one-parameter module has interval multiplicities, with half-open stage indices:

```text
mu[b,d) = r(b,d-1) - r(b-1,d-1) - r(b,d) + r(b-1,d),  0<=b<d<=N.
```

Set r(-1,j)=0 and use r(i,N)=0 only as a finite-end bookkeeping convention.
Each interval contributes one to r(i,j) exactly when it contains both stages;
the displayed finite difference isolates its birth and death. The usual
field-valued interval classification is discussed in [ZC05](bibliography.md#zc05).
For d=N the public death is `None`, meaning unknown continuation beyond the
observed stages, not essentiality of an unobserved infinite filtration.
Repeated numeric scales retain distinct stage indices and multiplicity.

The implementation prepares coordinate solves per request, retains only two
stage bases during adjacent preparation, then forms the rank table. Explicit
maps and the rank table have quadratic worst-case storage costs; this is not an
RSS bound. No persistent family cache or streaming barcode is introduced.
The [operator guide](../guides/homology-operators.md) owns the input, lifetime,
cost and execution API; the [comparison report](https://github.com/Aequiludium/cocycle-rs/blob/main/benches/reports/homology-operators.md)
separates finite implementation checks from these algebraic derivations.

### Exact affine minima and historical generators

For a cycle z, `minimum_class_mass` minimizes mass_w(x) over x=z+D c. A basis
of im(D) enumerates every representative once. The returned original-basis
witness solves D c=z+x independently of the search coefficients. Strictly
positive integer costs make the zero chain uniquely minimal when attainable.
Otherwise exhaustive completion proves optimality in this one class. This
does not minimize stretch over legal projections and leaves P unchanged.

For a fixed cycle z and an explicit upper-degree cost window, a filling solves
D c=z. A particular solution c0 and ker(D) parameterize all feasible fillings.
An inconsistent system proves infeasibility. Enumeration of c0+ker(D) with
checked positive integer costs yields a minimum; ties use the smallest original
coordinate bitset. The zero boundary has the zero filling as its unique minimum.

Historical feature extraction processes adjacent maps in birth order. When a
generator image depends on surviving older generators, XOR those combinations
through its complete prior history so it maps to zero at death. Extend the
surviving image basis by independent target coordinate units for new births.
Raw A/D checks verify cycles, full quotient bases, adjacent transport and death
membership. Endpoint multiplicities alone cannot supply these corrected chains.

For a certified one-simplex insertion source, let b be the birth event and d
the finite death event. The pair-volume chain c has c_d=1 and may use upper
simplices with b<j<d. Its boundary has coefficient one at b and zero on every
degree-k simplex with b<j<d. Minimize the explicitly supplied upper-chain mass
under these linear F2 constraints. This definition need not preserve the chosen
historical class. Stable-volume [Definition 6](bibliography.md#ob23) instead uses the exact cutoff
scale(b)+epsilon for both allowed upper simplices and zero-boundary rows, and
does not constrain the birth coefficient. A feasible death-only chain attains
the lower bound given by its positive cost. Otherwise enumerate the affine
constraint kernel. Supported nullities are bounded by 16, except the general
class query, which relies on the caller's cooperative execution limit.

The implementation certifies these finite optimization problems only. It does
not infer a perturbation theorem from an exact minimum. Integer scale bounds,
units, unsupported weight/location domains and failure behavior are specified
in the [operator guide](../guides/homology-operators.md#exact-minima-and-historical-simplex-features).

## 22. Real single-scale Hodge Laplacians

For a finite integer chain window C[k-1] <-A- C[k] <-D- C[k+1], require
AD=0 over the integers before interpreting incidence over R. In ordered bases
with standard inner products, define

$$
L_k = A^\mathsf{T} A + D D^\mathsf{T}.
$$

For k=0, A has zero rows: this is ordinary nonaugmented H0. An explicitly
absent adjacent degree in the supplied finite topology is a zero-size map.
Lower and upper selections retain only the respective summand.

For every real vector x,

$$
x^\mathsf{T}L_kx = \|Ax\|_2^2 + \|D^\mathsf{T}x\|_2^2 \geq 0.
$$

Thus L is symmetric positive semidefinite and ker L = ker A intersect ker D^T.
Because im D is contained in ker A and (im D) perpendicular = ker D^T,
orthogonal decomposition inside ker A gives

$$
\ker A = \operatorname{im}D \oplus \ker L_k,
\qquad
\dim\ker L_k = n_k - \operatorname{rank}_{\mathbb{R}} A
                         - \operatorname{rank}_{\mathbb{R}} D
             = \beta_k(\mathbb{R}).
$$

The fundamental subspace decomposition then gives
C[k] = im A^T direct-sum ker L direct-sum im D. Orthogonal projection onto ker L
assigns each homology class its unique harmonic representative. These are Lim's
Theorems 5.2 and 5.3, *Hodge Laplacians on Graphs*, SIAM Review 62(3), 2020,
685–715, [DOI 10.1137/18M1223101](https://doi.org/10.1137/18M1223101),
[author PDF](https://www.stat.uchicago.edu/~lekheng/work/hodge-graph.pdf), pp. 700–701.
The positive-norm argument is invalid over arbitrary finite fields.
Prime-field reduction must not certify real chain validity or harmonic nullity.

In edge order (01,02,12), an unfilled triangle has harmonic vector (1,-1,1)
and spectrum (0,3,3). Adding its oriented 2-cell contributes that vector's
outer product, giving L1=3I. A single edge has L0=[[1,-1],[-1,1]] and
spectrum (0,2). These expectations are independent of filtration weights.

The dense reference spectrum is an ordinary real symmetric eigenproblem using
scaled cyclic Jacobi on explicitly bounded dense matrices. Residuals, norm,
orthogonality and tolerance are recorded; thresholded zeros provide numerical
evidence only. Index-selected spectra cannot determine full nullity.
Explicit matrix construction and spectra check n² before dense allocation;
their shared default entry limit is 1,048,576 per square buffer (8 MiB of f64),
independent of work budgets. Spectra also retain a separate dimension limit.
These configurable per-buffer bounds do not cap aggregate process memory.
The matrix-free partial path instead diagonalizes S^T L S in a block space
S=span(X,R,P), retaining approximate smallest Ritz pairs and checking actual
sparse residuals. Its storage is O(n*b+b²) for fixed small block size b;
no full coordinate matrix or full n-value spectrum is required. Residual
convergence alone cannot certify extremal index coverage or missing multiplicity.
It provides no harmonic-nullity claim. Both paths retain the same standard
inner product and full/lower/upper operator, without weighted normalization.
General SPD Grams change adjoints to A*=Gk^-1 A^T Gprev and
D*=Gnext^-1 D^T Gk, and require a transformed/generalized eigenproblem.
That extension, complex conjugate adjoints and persistent Laplacians are outside
the initial API. See the [Hodge guide](../guides/hodge-laplacian.md) for callable
operations, costs, exclusions and reproducible checks.

## 23. Real persistent Laplacians

For a same-source pair K included in L, use real ordinary chains with orthonormal
increasing-vertex simplex bases. Let B denote the signed (q+1)-boundary of L,
with rows grouped into K and its complement: B = [R; F]. The allowed domain of
the persistent boundary is ker(F), with the inner product inherited from C_(q+1)(L).
For any orthonormal kernel basis N, its upper term is R N N^T R^T. For a general
full-rank kernel basis Z it is R Z (Z^T Z)^-1 Z^T R^T, not R Z Z^T R^T.

Partition B B^T = [A E; E^T C]. Then the persistent operator is
Delta_q(K,L) = B_q(K)^T B_q(K) + A - E C+ E^T, where C+ is the Moore--Penrose
pseudoinverse. Singular C is permitted. At q=0, B_0 has zero rows: the lower
term vanishes, with no augmentation. K=L gives the ordinary Hodge operator.
The operator is self-adjoint and positive semidefinite; its exact nullity equals
rank(H_q(K;R) -> H_q(L;R)). For integer boundaries an independent exact check is
n_q(K) - rank_Q(B_q(K)) - rank_Q(B) + rank_Q(F).

For a fixed endpoint L, eliminate the last active q-coordinate ell of the upper
term M by M'_ij = M_ij - M_iell M_ellj / M_ellell when the pivot is nonzero.
When the pivot is zero, positive semidefiniteness implies the eliminated row is
zero and the leading principal block remains. Iterate to each requested start,
then add its own lower term. This is reverse elimination at a fixed endpoint,
not arbitrary online filtration updates. Tied births introduce no intermediate
requested stage. Spectra transform by signed permutation under vertex relabeling.

A floating near-zero pivot is not the exact zero-pivot case: even the PSD block
`[[1,eps],[eps,eps^2]]` has a unit Schur correction. The numerical reverse route
therefore skips only exactly uncoupled coordinates and otherwise uses the checked
direct solve, or returns a numerical error under the original execution limits.
For computed discarded positive directions, the direct solve bounds the sum of
max-entry corrections `maxabs(c)^2/lambda`, normalized by the original upper Gram
scale, in addition to its actual solve/range residuals. This prevents another
large right-hand-side component from masking a lost correction. These finite
floating checks do not supply an exact rank or general forward-error certificate.

The [bounded dense implementation](../guides/persistent-laplacian.md) uses f64,
relative cutoff, checked residuals and explicit storage limits. Numerical zeros
and residual checks do not certify exact ranks or a forward-error bound.
Independent rational projection and rank calculations cover finite fixtures;
finite-field persistence only corroborates torsion-free examples.

Source: Mémoli, Wan and Wang, [Persistent Laplacians: Properties, Algorithms
and Implications](https://arxiv.org/html/2012.02808v3), arXiv:2012.02808v3
(2022-06-08), Section 2.2, Theorem 2.7, Theorem 4.6, Algorithm 2 and Section 5.2,
equation 15 / Theorem 5.1. The implementation is independently written from
these formulas, rather than copied from reference code.

For U=A-E C+ E^T and D=B_q(K)^T B_q(K), the chain identity B_q(K) R N=0
implies DU=UD=0. Thus their images are orthogonal. In this finite-dimensional
setting, the nonzero spectral multiset of U+D is the union of those of U and D;
its zero multiplicity is n_q(K)-rank(U)-rank(D), equivalently
z(U)+z(D)-n_q(K). One must not concatenate both complete spectra, or infer these
exact ranks from arbitrary floating-point near-zero cutoffs.

The implementation exposes full/upper/lower results and separately tagged
spectra. The lower term depends only on K and requires its q-skeleton; it does
not require L or a q+1 skeleton. Full and upper operators retain the stronger
q+1 coverage check. Each result retains its owned bases and source metadata.
Numerical spectrum tags identify the component, not the pair/source identity.

Wolf, Fan and Monod, [Generalized Persistent Laplacians and their Spectral
Properties](https://arxiv.org/html/2509.20220v1), arXiv:2509.20220v1 (2025-09-24),
Section 5.2 / Theorems 5.6--5.7 and Section 5.3 / Theorem 5.11, study component
spectral functions and generalize this decomposition. Their infinite-dimensional
operator domains and stability hypotheses are not implemented here. The full
operator does not inherit all component comparison properties; this does not
rule out the fixed-K, increasing-L monotonicity, where D stays constant.

## 24. Standard real discrete Dirac

For a finite supplied chain complex with ordered bases, use ordinary,
nonaugmented chains: C[-1]=0. Signed integer boundary matrices
B[j]:C[j]->C[j-1] satisfy B[j-1]B[j]=0 over the integers, before any modular
reduction. The standard real inner product has Gram I in every degree, so the
adjoint is the transpose. Filtration values select stages and are not matrix
weights. See [CSB23](bibliography.md#csb23), section 2.

On C[a] + ... + C[b], let d be the degree-lowering block matrix retaining only
B[a+1], ..., B[b], and let K=d+d^T. Thus
`(Kx)[j] = B[j]^T x[j-1] + B[j+1] x[j+1]`, with terms outside [a,b] zero.
Since d^2=0 and (d^T)^2=0, K^2=dd^T+d^Td is block diagonal, with degree-j
block `B[j]^T B[j] + B[j+1] B[j+1]^T` under the **same endpoint truncation**.
K is real symmetric. Replacing oriented incidence by absolute values or lifting
F2 coefficients to 0/1 may break d^2=0 and this identity.

Let Gamma multiply degree-j coordinates by (-1)^j. Adjacent degrees have
opposite parity, so Gamma K=-K Gamma. Therefore each nonzero eigenvalue lambda
pairs with -lambda, with equal multiplicity. The nonzero spectrum is the union
of +/- the nonzero singular values of the retained adjacent B[j]: each singular
pair gives coupled vectors in degrees j-1 and j, and the chain condition makes
these boundary summands orthogonal. Each resulting vector retains its complete
graded coordinates. It must not be labeled as one homology degree.

Also `x^T K^2 x = ||dx||^2 + ||d^T x||^2`; hence ker K is the direct sum of
the graded harmonic spaces. On a complete ordinary finite real chain complex,
its dimension is the sum of the real Betti numbers. For a truncated interval,
the kernel belongs to that interval, and generally differs from the source's
homology. A floating eigenvalue threshold is numerical evidence, not an exact
homology certificate; F2 Betti numbers need not equal the real ones.

For one edge ordered (v0,v1,e01),

```text
K   = [[0,0,-1], [0,0,1], [-1,1,0]]
K^2 = [[1,-1,0], [-1,1,0], [0,0,2]]
spectrum = (-sqrt(2), 0, sqrt(2))
```

The unfilled triangle has two +/-sqrt(3) pairs and two zero eigenvalues;
the filled triangle has three pairs and one zero eigenvalue. Restricting the
single edge to C1 instead gives K=[0], while its full edge Laplacian is [2].
These fixed matrices, independent Hodge blocks, orientation/parity, empty
degrees and exact integer-chain rejection are exercised in
[`tests/dirac.rs`](../../tests/dirac.rs).

`DiracOperator` groups coordinates by degree while preserving source traversal
within each degree. The implementation accepts at most `MAX_DIRAC_DEGREES`
(4096) consecutive degrees, rejecting longer intervals before offset allocation
or source access even when the source is empty or execution is unlimited. This
bounds degree metadata, not source size, and does not restrict absolute labels.
Duplicate signed terms are summed in checked i64 arithmetic;
zero sums are removed. Coefficients must lie in [-2^53,2^53] for exact f64
conversion. Retained two-step sums use checked i128 arithmetic. Numeric overflow
is an error, including validation overflow; no modular or tolerance-based chain
test substitutes for integer zero. Selected cells' nonzero boundary references
must be earlier existing cells of degree one less, even when a bottom-endpoint
map is discarded. Unselected boundaries are not read or certified.

Sparse action visits incidence once and accumulates both its downward and
transpose contributions with finite f64 products/partial sums. No dense matrix
is needed. Explicit row-major export checks n*n and the caller's scalar-entry
cap before allocation. Cooperative work/cancellation errors return no partial
operator or buffer. The API supports only this real standard-inner-product
operator: no Gram conversion, complex conjugate adjoint, general global spectrum,
augmentation, persistent Dirac, mass term or signal filter is inferred.

For a retained boundary B[j] and a positive singular value sigma, let
B[j]v=sigma*u and B[j]^T u=sigma*v. Then u lies in im B[j] and v lies in
im B[j]^T. The integer chain condition gives B[j-1]u=0 and B[j+1]^T v=0
inside this same selected interval. Consequently w+=(u,v)/sqrt(2) is an
eigenvector of the entire retained K with eigenvalue sigma; Gamma*w+ has
eigenvalue -sigma, where Gamma is absolute-degree parity. Different retained
boundaries contribute orthogonal nonzero chain subspaces, so their nonzero
singular values contribute separate signed multiplicities. This recovers the
nonzero spectrum without independently diagonalizing the entire Dirac matrix.

The implemented boundary-spectrum path approximates selected largest singular
modes through matrix-free B[j]^T B[j] block LOBPCG and a small Ritz problem.
It verifies the reconstructed signed residuals with the actual truncated K.
Thresholded positive modes do not certify rank, the kernel or extremal coverage.
Normal-map tolerances have sigma² units; division by |sigma| translates the
absolute error scale to signed Dirac residual units and reflects conditioning.
