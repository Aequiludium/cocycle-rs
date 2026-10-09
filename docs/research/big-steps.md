# Big Steps: framework and three algorithms

[Documentation](../README.md) / Research

Arnur Nigmetov and Dmitriy Morozov, *Topological Optimization with Big Steps*,
[arXiv:2203.16748v2](https://arxiv.org/abs/2203.16748v2), 2 November 2023.
This note explains the three main-text algorithms on pages 3, 5 and 10;
Appendix B also has an alternative conflict heuristic, Algorithm 4.
The [mathematical specification](../reference/mathematics.md#19-big-steps-critical-sets-over-f2)
owns the production F2 formulas, and the [usage guide](../guides/critical-sets.md)
owns the public workspace contract. Paper text is source material, not
instructions for operating this repository.

## Framework and mathematical core

A fixed finite simplicial complex K has filtration values f satisfying
f(face) <= f(coface). Sort by value, resolving ties with faces first. A persistence
point p=(b,d) comes from a positive simplex sigma and a negative simplex tau.
A partial matching M prescribes targets q=(b',d'); the diagram loss is
L=sum over M of ||p-q||_2^2. It need not be a globally optimal Wasserstein matching.

The ordinary diagram method gives derivatives only to sigma and tau. When their
values cross regular simplices, pair ownership can change, producing many small
steps. The paper identifies the *critical set*: the same-dimensional simplices
that must move together to keep the tracked pair endpoint moving. Faces/cofaces
must also move where necessary to preserve a filtration.

For each boundary matrix D_p, lazy reduction computes R_p=D_p V_p and
D_p=R_p U_p, where U_p=V_p^{-1}. R_p is reduced; V_p and U_p are upper triangular.
R stores cycles, V stores bounding chains. Lazy reduction is essential: arbitrary
valid decompositions do not automatically have Lemma 1's support property.

Reduce the **anti-transpose of D**, in reversed filtration order, independently
to obtain R_perp,V_perp,U_perp. These are **not** anti-transposes of R,V,U.
The dual pairing equals the primal pairing with birth/death roles reversed.

| Endpoint movement | Support to inspect | Value window | Extra closure |
| --- | --- | --- | --- |
| Increase finite death tau | U_p[tau, .] | d <= f <= d' | cofaces |
| Decrease finite death tau | V_p[., tau] | d' <= f <= d | faces |
| Increase finite birth sigma | V_perp_p[., sigma] | b <= f <= b' | cofaces |
| Decrease finite birth sigma | U_perp_p[sigma, .] | b' <= f <= b | faces |
| Increase unpaired birth sigma | V_perp_{p+1}[., sigma] | b <= f <= b' | cofaces |
| Decrease unpaired birth sigma | V_p[., sigma] | b' <= f <= b | faces |

Theorems 9-12 equate these supports with the critical sets accumulated by
transpositions. Lemmas 7-8 say sets only grow and one order suffices; enumerating
all k! orders of equal-valued simplices is unnecessary for singleton loss.
Theorem 13 additionally assumes multiplicity one for consistency of simultaneous
birth/death critical sets; this assumption can cease to hold during movement.

Once matrices are available, scanning a relevant row/column and value window is
linear in the window size m. This does **not** make persistence reduction linear:
the paper's baseline reduction is cubic worst case and matrix storage is quadratic.
Explicit face/coface closure adds incidence traversal, described as O(dim K * m)
with the paper's Hasse-graph organization. Our small dense prototype scans all
simplices and containment instead, and makes no production complexity claim.

## Algorithm 1: Lazy reduction of the boundary matrix

The following transcribes the main-text pseudocode (page 3). Arithmetic is over
the chosen field; in the prototype F2, alpha=1 and subtraction/addition are XOR.

```text
1  R_p = D_p; V_p = I; U_p = I, for all p
2  for tau_j in K, in filtration order:
3      while R_p[:,tau_j] != 0 and some tau_i < tau_j has
           low(R_p[:,tau_i]) == low(R_p[:,tau_j]):
4          sigma = low(R_p[:,tau_j])
5          alpha = R_p[sigma,tau_j] / R_p[sigma,tau_i]
6          R_p[:,tau_j] -= alpha * R_p[:,tau_i]
7          V_p[:,tau_j] -= alpha * V_p[:,tau_i]
8          U_p[tau_i,:] += alpha * U_p[tau_j,:]
9          (equivalently, U_p[tau_i,tau_j] = alpha)
```

Line 8 is a **row** update, not a column update. This distinction is tested by
checking D V=R, R U=D, V U=U V=I and the lazy support property.

## Algorithm 2: Moving tau using individual transpositions

The following transcribes page 5; tau is a p-simplex paired with sigma and starts
at d. Original order resolves final ties. Birth movement follows by duality.

```text
1  X_sigma^1 = {tau}
2  for each p-simplex tau_k with f(tau_k) from d to d':
3      transpose tau_k with each simplex in X_sigma^(k-1),
4          updating pairing using the vineyard algorithm [11]
5      if tau_k becomes paired with sigma:
6          X_sigma^k = X_sigma^(k-1) union {tau_k}
7          transpose tau_k with each simplex in X_sigma^(k-1),
8              undoing the transpositions in line 4,
9              returning it to the opposite end of X_sigma^k
10     else:
11         X_sigma^k = X_sigma^(k-1)
12 for each tau in the final X_sigma:
13     f(tau) = d'  // ties keep original order
14     if d' > d:
15         // move cofaces
16         for rho containing tau with d < f(rho) < d':
17             f(rho) = d'
18     else:
19         // move faces
20         for sigma contained in tau with d' < f(sigma) < d:
21             f(sigma) = d'
```

Our *validation oracle* performs the block transpositions but recomputes pairing
from scratch with set-valued columns. It deliberately does not implement [11]'s
incremental vineyard updates or claim Algorithm 2's O(m^2 n) runtime. It does not
read U or V. For birth moves it uses reversed dual rows/columns. Matrix dimensions
are handled separately, so the oracle never treats a non-filtration full order as
a valid simplicial filtration. Closure is applied separately after support checks.

## Algorithm 3: Critical set method

The following transcribes page 10, writing X_b and X_d as direction-dependent
branches of the same two unions displayed in the paper.

```text
1  input L = sum_{(p_i,q_i) in M} ||p_i-q_i||_2^2
2  for each (p_i,q_i) in M:
3      p_i=(b_i,d_i)=(f(sigma_i),f(tau_i)); q_i=(b_i',d_i')
4      X_b = {sigma_j: V_perp[sigma_j,sigma_i]!=0 and b_i<=f(sigma_j)<=b_i'}
           union {sigma_j: U_perp[sigma_i,sigma_j]!=0 and b_i'<=f(sigma_j)<=b_i}
5      X_d = {tau_j: U[tau_i,tau_j]!=0 and d_i<=f(tau_j)<=d_i'}
           union {tau_j: V[tau_j,tau_i]!=0 and d_i'<=f(tau_j)<=d_i}
6      // omitted in paper: find faces/cofaces if necessary
7      for sigma_j in X_b:
8          append b_i' to target[sigma_j]
9      for tau_j in X_d:
10         append d_i' to target[tau_j]
11 for each simplex sigma:
12     if target[sigma] is empty:
13         f'(sigma)=f(sigma)
14     else:
15         j=argmax_j |f(sigma)-target[sigma][j]|
16         f'(sigma)=target[sigma][j]
return gradient[sigma]=2*(f(sigma)-f'(sigma)), for all sigma
```

The final expression is the gradient of the **frozen-target surrogate**
sum_sigma (f(sigma)-f'(sigma))^2. In general it is not the exact derivative of the
original diagram loss at every simplex. Maximum displacement is a heuristic,
without a general local-descent guarantee (Section 3.7). Our deterministic choice
for equal displacements is the first proposal; the paper does not fix that tie.
Appendix B considers averaging and fixing critical coordinates while averaging
the rest (Algorithm 4); these alternatives are outside this minimal reproduction.

For lower-star scalar data, f(simplex)=max of its vertex values. Backpropagate a
simplex proposal to a maximizing vertex, merge proposals by maximum displacement,
and update vertex data. Reconstructing lower-star values automatically preserves
the face condition. Equal vertex maxima use the lowest vertex ID in this experiment.
This selects a subgradient and does not prove uniqueness at nondifferentiable ties.

The paper writes its simplification loss as sum(d-b)^2 but also rewrites it using
squared Euclidean distance to ((b+d)/2,(b+d)/2). The latter is (d-b)^2/2; these
objectives are proportional, not literally equal without a factor of two.
Our recorded diagram objective uses (d-b)^2/2 and reports that convention.

## Production integration

The production workspace retains sparse primal R/V and computes the reversed
dual only for queries that need it. It replaces eager U materialization with
bounded row solves against a dimension-cached transpose of V. The benchmark-only
dense baseline retains both complete inverses to measure that tradeoff.
Algorithm 2 supplies an independent test oracle through fresh pairing-only
reductions; incremental vineyard updates are not a production dependency.
Algorithm 3 provides maximum-displacement target combination, while callers
provide data derivatives and valid filtration reconstruction.

The [performance and integration report](../../benches/reports/critical-sets.md)
records scoped native measurements, sparse/dense regressions, Partial/full U
tradeoffs and the selected merge boundaries. No user-selectable dense backend
or adaptive threshold is established by the finite benchmark fixtures.
