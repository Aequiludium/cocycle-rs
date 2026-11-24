# Non-branching weak reduction and induced metrics

[Documentation](../index.md) / Research

This page records the bounded reproduction and Rust rewrite for
[#147](https://github.com/Aequiludium/cocycle-rs/issues/147), independently
implemented from Botnan and Dong, *The persistent Laplacian of non-branching
complexes*, [arXiv:2408.16741v2](https://arxiv.org/abs/2408.16741v2),
Algorithms 1-3, PDF pages 8, 9 and 13. It separates the literal pseudocode from
explicit corrections and a compact Rust experiment. The first-round standalone
replay remains separate from the second-round
[explicit library constructors](../guides/persistent-laplacian.md#explicit-non-branching-construction).
Neither stage enables automatic backend selection or adds cubical storage.

## Mathematical scope

Use real, ordinary chains, identity weights and fixed ordered signed bases.
For a supplied boundary split `B = [B1; D]`, the target is
`U = B1 Z (Z^T Z)^(-1) Z^T B1^T`, where the columns of Z span ker(D).
The [persistent operator specification](../reference/mathematics.md#23-real-persistent-laplacians)
owns this metric convention. Filtration values select cells; they are not weights.

Eligibility is **pair-local**: D has entries in {-1,0,1} and at most two
nonzeros per row. B1 is also signed unit incidence, but its row occupancy is
unrestricted. These algebraic checks do not certify a simplicial boundary,
source identity or a valid inclusion. Synthetic signed controls are labeled
as algebraic data; in particular `[[1,1],[1,-1]]` is not asserted to be a
realizable simplicial boundary. Arbitrary Rips complexes need not qualify.

The unsigned graph has one vertex per nonzero column of D and one edge for
each row with two entries. A singleton row forces its whole connected
component to zero. In a remaining component, orient its first column with
sign +1 and traverse signed relations `z_j = -D_ri D_rj z_i`. An inconsistent
cycle has no nonzero real kernel direction. A consistent component contributes
one signed indicator, whose squared norm is its support size. Original zero
columns contribute independent coordinate vectors exactly once.

These directions have disjoint supports, so `Z^T Z` is diagonal. Its inverse,
not the identity, weights the unnormalized boundary factor `B1 Z`.
For the endpoint path control `D=[1,-1]`, the induced metric is [2], and
`U = [[1,-1],[-1,1]]/2`; deleting the outside row would incorrectly return I.

## Literal pseudocode and disclosed corrections

The controller retains the original Algorithms 1-2 separately from the
corrected route, including their failed identities and upper matrices.
The following are observations about this fixed printed pseudocode, not a
claim that the paper's mathematical characterization is false.

| Location | Literal behavior and finite counterexample | Explicit correction |
| --- | --- | --- |
| Algorithm 1, lines 10 and 12-19 | Both flag multiplication and mutation of D apply signs. On `D=[[1,1],[1,-1]]` it accepts a spurious kernel direction, although det(D)=-2. | Keep D immutable and apply accumulated signs once. |
| Algorithm 1, lines 6-26 | One arbitrary row scan can reject an orientable component even after removing D mutation: rows `[1,1,0,0]`, `[0,0,1,1]`, `[1,0,1,0]`. | Traverse the unsigned component before checking signed relations. This uses an explicit graph traversal, not a signed DSU replacement. |
| Algorithm 2, lines 3 and 18-20 | Zeroing only selected columns does not update other reoriented columns. For `[[1,1,0],[1,0,1]]`, the selected kernel is valid but returned R differs from D E V. | Explicitly form R = D E V, including all reoriented columns. |
| Algorithm 3, line 5 | The equal-domain branch prints `B_q B_q^T`, with the wrong degree/domain for the upper term. | Use `B_(q+1) B_(q+1)^T` in identity weights. |
| Algorithm 3, lines 10 and 14 | The displayed D superscript and component subscript differ from the preceding definitions. | Reduce the outside-row D, and map reduced component indices back to original upper columns. |

The corrected weak reduction retains E, V and R on the nonzero columns of D.
V has diagonal ones and ones in the **last** column of each regulable component;
E is its orientation diagonal. Zero columns of R select the corresponding
columns of E V. Original zero columns bypass that reduced index space. The
Python output freezes both maps, supports, signs, kernel, metric, integer
boundary factor and exact rational U.

## Reproduction command and evidence

The repository tools are
[`reproduce_nonbranching.py`](https://github.com/Aequiludium/cocycle-rs/blob/main/tools/reproduce_nonbranching.py),
[`nonbranching.rs`](https://github.com/Aequiludium/cocycle-rs/blob/main/tools/reference/nonbranching.rs)
and their focused Python tests. They require NumPy, SciPy, GUDHI and rustc;
these remain external experiment dependencies. Use a new output directory:

```sh
python -m unittest discover -s tools -p test_reproduce_nonbranching.py
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 python tools/reproduce_nonbranching.py --paper PATH_TO_2408_16741V2_PDF --author-source PATH_TO_SRC_SPARSE_UPLAP_PY --output target/nonbranching-replay --benchmark
```

The controller checks paper SHA-256
`343e6c6e076ce0a2d0022f45fab35c0741dcde0c7af795f5d253b4640296c483`
and the unchanged author's `src/sparse_uplap.py` SHA-256
`03af3d06ab5fcea93e7e7ac7dd2e269ad48ef2fd07df3b9bc8162e021158edd8`
from [commit 13d8557](https://github.com/ruidongsmile/UpPersLapCubicalComplex/tree/13d855749bddc5e4d8cca656d5e0775a56a581a1).
It retains original sources, environment/version records, input hashes,
literal and corrected intermediate matrices, build logs, executable hashes,
rejections and failures under the output directory. Downloaded author code is
not redistributed as Rust source. No author-code license grant is assumed.

Example 3.9 uses the paper's exact ten edge and five triangle labels. Its
outside matrix is irregular, has rank five, and yields a zero upper term.
Example 6.1 freezes all twelve edge and four square labels and the mask
excluding edge [52]. The kernel metric is diag(1,1,2) when original zero columns
are returned first; the paper lists the merged component first and obtains
diag(2,1,1). The resulting upper matrix is identical after this column permutation.

Author replay follows the notebook core:
`find_D -> find_root -> uf.find_loops -> uf.find_loops_roots ->
find_rel_bdry_opr -> find_diag -> find_up_persistent_Laplacian`.
It covers the paper matrix, the all-rows mask and a regenerated 2x2 cubical
lattice through the author's boundary/mask constructor. The lattice has one
edge outside K and four squares in L. Its vertical-edge orientation differs
from the paper's printed orientation; its own exact signed matrices are retained.
No source compatibility patch is needed in the tested environment. This is a
small replay of the notebook's computational core, not its chest image,
incremental edge-update loop or published timing experiment.

The dated Windows replay (2026-11-24) used Python 3.12.0rc3, NumPy 2.4.6,
SciPy 1.17.0, GUDHI 3.11.0 and rustc 1.98.1. All 90 frozen cases passed exact
corrected weak-reduction checks and 180 Rust baseline/compact comparisons.
The largest numerical relative Frobenius error was 3.93e-15. Literal pseudocode
failed R=D E V in 10 cases and the exact upper matrix in 8 cases; these
failures remain in the evidence, not converted into successful reproduction.
The result is **reproduced with disclosed corrections**, while literal general
signed/row-order correctness remains unresolved by the original pseudocode.

The independent oracle reuses the existing rational Gaussian-nullspace and
general Gram inverse checker, without component classification. Exact checks
include R=D E V, full kernel dimension, D Z=0, independent nonzero R columns,
diagonal induced metric and exact U. A separate numerical pseudoinverse
projector checks Frobenius error against `1e-10*max(1,||U||F)` and
`||D P||F/max(1,||D||F) <= 1e-9`. These tolerances concern small unit-scale
controls, not a general conditioning guarantee.

GUDHI checks path, triangle, disk, Mobius and shell inclusions over F2 and F11.
It corroborates topology only. The Mobius boundary inclusion has real and F11
image rank one but F2 image rank zero, illustrating why finite-field persistence
cannot certify the real spectrum. The focused tests also exhaust all 361
two-row, three-column signed non-branching matrices, independently checking
their projectors through B1=I.

## Rust parity and compact experiment

The first Rust route retains dense R/E/V/Z and is compared entry by entry to
the corrected first-round matrices. It does not claim parity with the failing
literal pseudocode. Only after those checks does the compact candidate omit
these dense diagnostics, retaining component supports and signs with the same
orientation traversal and metric. It is a representation experiment, not a
signed-DSU optimization or production integration.

Both routes construct `F_c = B1 z_c/sqrt(size(c))`, apply `F F^T` and export
the same dense U. Each dense input, transformation, factor and output buffer
has a checked 65,536-entry cap. The cap is per buffer, not total process memory.
Unsupported rows, nonunit entries, excessive dimensions and trailing input are
errors. No dense fallback or automatic cap increase is attempted. This standalone
worker has no library `Execution` budget/cancellation contract.

The optional benchmark freezes 32- and 128-column algebraic path-constraint
matrices, eight retained rows and a dense 8x8 output. It serially runs one
discarded warmup per route and ten fresh-process paired samples per size,
alternating order. Timed runs suppress transformation export in **both** routes
while retaining the baseline's construction cost. Timers separate input
preparation, reduction, normalized factor, action, dense export, JSON formatting
and process wall time. Python full eigvalsh timing is separate from Rust timing.
Windows process peak working set includes input, runtime, allocator and output;
on other platforms it is explicitly unavailable. Incidence/factor nonzeros,
support counts, dense diagnostics and returned entries are recorded separately.

For the same dated run, the median of ten fresh samples was:

| Columns | Route | Reduction/factor/action/dense total | Process wall | Peak working set | Dense diagnostic entries |
| --- | --- | --- | --- | --- | --- |
| 32 | Dense weak-reduction baseline | 72.65 us | 7.480 ms | 4,530,176 bytes | 3,072 |
| 32 | Compact candidate | 20.30 us | 7.390 ms | 4,509,696 bytes | 0 |
| 128 | Dense weak-reduction baseline | 3,316.75 us | 11.561 ms | 4,751,360 bytes | 49,152 |
| 128 | Compact candidate | 42.60 us | 7.353 ms | 4,593,664 bytes | 0 |

The paired process-wall ratios (baseline/compact), at the 10th percentile,
median and 90th percentile, were (0.927, 1.007, 1.154) for 32 columns and
(1.428, 1.605, 1.762) for 128. The small case has no stable process-wall gain;
unfavorable paired samples remain in the recorded distribution.

The baseline explicitly multiplies D E V and is intentionally unoptimized.
These numbers compare research representations; they are not measurements
against Cocycle's production Schur implementation. Startup dominates the small
case, and each action is too short for a useful isolated speed claim. Both
return all 64 upper entries and eight action entries. Python spectrum times
are recorded separately and do not establish a spectral-method improvement.

The dense input protocol still scans m*p entries, and the final U may be dense.
Small fixed controls cannot establish near-linear total complexity, a complete
spectral speedup or general production scalability. Algorithm 7/GKL, arbitrary
weights, larger author datasets and automatic dispatch remain outside this
delivery. The first-round decision retains the compact representation as a
research candidate. The second round below integrates the unoptimized corrected
baseline; it does not adopt that compact representation.

## Second-round Rust library rewrite

The Unreleased `PersistentLaplacian::compute_non_branching` and
`compute_up_non_branching` explicitly select the corrected weak reduction.
They use the existing source-pair validation, owned bases, coverage/skeleton
metadata, full/upper result tags, action and complete-spectrum operations.
Ordinary Schur constructors remain unchanged. The private owner is
`src/persistent_laplacian/nonbranching.rs`; production code does not depend on
the standalone first-round worker or the Python oracle.

This baseline retains dense R, diagonal E, V and Z. It builds the actual
R=D E V and checks the selected zero columns, then forms the original zero
directions exactly once. Its kernel constraint check, graph traversal, dense
arithmetic and output share one `Execution` budget. Every f64 incidence,
transform, kernel, factor, Gram and returned matrix buffer has the supplied
per-buffer cap. The original target Gram is built and released for the existing
reference-scale contract. No ineligible or oversized input triggers a fallback,
and `solver_residual()==0` denotes the absence of a Schur solve, not an exact
floating assembly or rank certificate.

To reuse a first-round evidence directory, run the original command with
`--library --frozen-round1 PRIOR_REPLAY`. It verifies that input hashes have
not changed, copies the frozen R/E/V/Z and metric outputs, compiles the actual
kernel source with the actual budget/allocator, and compares all 90 cases
including failure semantics. This matrix worker is a separate source-level
build for algebraic controls; the public wrappers are checked through the
ordinary crate separately:

```sh
cargo build --locked --release --example persistent_laplacian
python tools/compare_persistent_laplacians.py --worker target/release/examples/persistent_laplacian --output target/nonbranching-library-comparison --non-branching --gudhi
```

The dated 2026-11-24 public comparison checked 1,026 qualifying pairs, all with
independent rational matrices. It compared the full/upper/lower matrices,
original-coordinate actions, complete spectra, multiplicities and exact-rank
spectral unions to independent projection and the generic reverse-elimination
comparator. Maximum matrix error was 2.89e-15 and maximum component/spectral-union
error 4.45e-15. GUDHI corroborates only the declared topology cases. Branching
outside rows are retained in an exclusion manifest and are not counted as
successful comparisons; library tests separately verify their explicit errors.
Fallback counters in this comparison describe the generic reverse comparator,
not the non-branching constructors, which never fall back.

Public tests also cover owned metadata after source release, empty/absent
degrees, cloned and mixed-certificate sources, insufficient skeletons, degree
overflow and incidence/kernel storage rejection. Private tests deterministically
interrupt different phases with cancellation and work limits. This establishes
bounded Rust parity and interface integration, not a performance adoption or
the third-round sparse/compact optimization.
