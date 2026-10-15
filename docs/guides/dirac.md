# Standard real discrete Dirac operators

[Documentation](../README.md) / Guides

This **Unreleased** capability implements #74's standard-real Dirac action and
matrix construction on a finite [supplied complex or stage](filtered-complexes.md).
It reads signed integer boundaries directly, preserves graded coordinates and
uses the Euclidean inner product. It performs no persistence request and does
not turn filtration values into weights. `boundary_spectrum` also computes
selected nonzero signed modes from each retained boundary using a matrix-free
normal-map solve, with actual truncated-Dirac residual and chirality checks.

## Apply Dirac to one edge

```rust
use cocycle::complex::{Simplex, SimplicialComplex};
use cocycle::execution::Execution;
use cocycle::spectral::{DiracOperator, IterativeOptions};

let source = SimplicialComplex::new(vec![
    Simplex::new(vec![10], -2.)?,
    Simplex::new(vec![20], -1.)?,
    Simplex::new(vec![10, 20], 3.)?,
])?;
let stage = source.stage(3.)?;
let dirac = DiracOperator::new(&stage, 0..=1)?;
assert_eq!(dirac.degree_offsets(), &[0, 2, 3]);
assert_eq!(dirac.action(&[2., 5., 7.])?, [-7., 7., 3.]);
assert_eq!(dirac.matrix(9)?, [0., 0., -1., 0., 0., 1., -1., 1., 0.]);
assert_eq!(dirac.entries().count(), 4);
let pairs = dirac.boundary_spectrum(1, 1, 1e-8, &IterativeOptions::default(), &Execution::default())?;
assert!((pairs.values()[0] + 2_f64.sqrt()).abs() < 1e-9);
assert!((pairs.values()[1] - 2_f64.sqrt()).abs() < 1e-9);
assert!(pairs.residuals().iter().all(|&r| r < 1e-9));

// Every terminal creates a fresh cooperative budget.
let execution = Execution::default().max_work(100);
assert_eq!(dirac.action_with(&[2., 5., 7.], &execution)?, [-7., 7., 3.]);
assert_eq!(dirac.matrix_with(9, &execution)?, dirac.matrix(9)?);

// Removing C0 removes B1 and its transpose, so the edge-only operator is zero.
let edges_only = DiracOperator::new(&stage, 1..=1)?;
assert!(edges_only.is_truncated());
assert_eq!(edges_only.degree_offsets(), &[0, 1]);
assert_eq!(edges_only.action(&[7.])?, [0.]);
assert_eq!(edges_only.matrix(1)?, [0.]);
# Ok::<(), cocycle::Error>(())
```

The negative, unequal vertex births make the two vertex coordinates unambiguous
in this example. In general use `basis()`, rather than guessing coordinate order
from vertex labels. Run `cargo run --locked --example dirac` for a filled
triangle; its printed basis explains every coordinate of its action vectors.

## Mathematical definition and endpoint truncation

On the ordered direct sum `C[a] + ... + C[b]`, retain signed adjacent boundaries
`B[j]:C[j]->C[j-1]` only for `a < j <= b`. Let d be their degree-lowering block
operator. The standard real Dirac operator is K=d+d^T, with

```text
(Kx)[j] = B[j]^T x[j-1] + B[j+1] x[j+1].
```

Terms crossing a selected endpoint are zero. Ordinary nonaugmented chains have
`C[-1]=0`; no all-ones vertex sum is added. With the retained integer chain
condition d^2=0, K is symmetric and K^2 is block diagonal, with degree-j block
`B[j]^T B[j] + B[j+1] B[j+1]^T` using those same endpoint conventions. See the
[mathematical specification](../reference/mathematics.md#22-standard-real-discrete-dirac).

Thus the edge-only K above is [0], while the full source's edge Hodge Laplacian
is [2]. Squaring a truncated K does not recover the source's full endpoint
Laplacians. For an interval with more than one degree, its bottom Hodge block is
upper-only, its top is lower-only, and interior blocks contain both summands.
For a single selected degree, K and its square are zero.

## Source, coordinates and owned results

`new_with` scans all source IDs, dimensions and filtration values. It checks
unique IDs and finite nondecreasing values, then groups selected cells by
increasing degree and original traversal position. Within each degree, source
order is preserved. `degree_offsets()` includes empty degrees and a final total
coordinate count; repeated adjacent offsets represent a zero-dimensional space.
`degrees()` returns the inclusive interval that was actually selected.

Both constructors reject intervals containing more than
`spectral::MAX_DIRAC_DEGREES` (4096) degrees with `InvalidParameter` for `degrees`,
before allocating offsets or reading the source. This applies to empty sources
and unlimited execution, bounding the offset table to 4097 `usize` entries
(32,776 bytes on a 64-bit platform). The limit concerns interval length, not
absolute degree labels; a narrow interval at a high degree remains valid.
Unrepresentable offset counts still return `SizeOverflow`. This metadata bound
does not cap source-dependent allocations or process memory.

For selected cells, nonzero boundary terms must name earlier existing cells of
degree one less, including terms discarded at the bottom endpoint. Zero terms
are ignored. Duplicate signed terms are summed in checked i64 arithmetic and
zero sums removed. Summed incidences must lie in [-2^53,2^53], making conversion
to f64 exact. Retained two-step chain products are summed in checked i128
arithmetic. A nonzero integer boundary square is rejected even if it is zero
modulo two. Boundaries of unselected cells are not read or certified.

The operator owns its sparse columns, source-local basis handles and offsets.
It remains usable after source release. Keep the source description when IDs
need geometric interpretation: equal numeric IDs do not identify different
owners. `is_truncated()` reports whether any supplied cell was excluded; it does
not certify a larger construction's missing skeleton or filtration coverage.

`action` returns an owned length-n vector and `matrix` returns an owned n by n
row-major buffer. Both use this operator's basis and offsets; the buffers do not
carry a duplicate metadata object. Preserve that metadata when exporting them.
A correct length cannot detect a caller's undeclared coordinate permutation.
`entries()` streams both each retained signed boundary entry and its transpose
as `(row, column, value)` triples, without allocating a matrix.

Pass an explicit stage to choose a single scale. Passing a whole stored filtered
complex includes all its cells in the selected degrees, regardless of birth
values. This is an ordinary supplied-chain operator, not persistent Dirac between
two stages. General Grams, complex coefficients, normalization, augmentation,
mass terms and filtering cannot be expressed by this API.

## How sparse action and matrix export run

For every retained term `(row, coefficient)` in boundary column `column`, action
adds `coefficient * x[column]` to output `row`, and
`coefficient * x[row]` to output `column`. The first contribution applies d;
the second applies its transpose. Incidence is stored once. Products and partial
sums must stay finite; input NaN/infinity, overflow or wrong length is an error.
Accumulation follows the frozen sparse column order without silent rescaling.

If n is the coordinate count and m the coalesced incidence count, action takes
O(n+m) work and O(n) extra output storage. Construction also indexes every source
cell, sorts selected cells and duplicate terms, and visits two-step incidence
paths for exact chain validation; it need not be linear in m. Sparse columns
may retain capacity from duplicate aggregation.

`matrix(max_entries)` checks n*n and the caller's scalar-entry cap before any
dense allocation. It initializes a zero buffer, then writes the sparse triples
directly; it does not run n separate actions. Cost is O(n^2+m) work and O(n^2)
output storage. An empty domain exports an empty buffer with a zero cap. The cap
counts f64 scalars, not total process memory; the scalar payload alone is 8*n*n
bytes. Fallible reservations do not make every standard-library allocation
failure recoverable.

Each `_with` operation starts a fresh `Execution` budget. Construction charges
source/selected-cell visits, degree slots, incidence aggregation and two-step
terms. Action charges input/output coordinates and retained incidence; matrix
export charges dense entries and emitted signed entries. Cancellation, work
exhaustion and numerical/allocation errors return no partial operator or buffer.
The frozen operator can be reused after a failed action or export. Sorting,
allocation and source iterator internals are cooperative checkpoints, not
interruptible instructions. Work units are neither elapsed time nor a memory quota.

## Spectrum properties and the current solver boundary

Degree parity anticommutes with K, so nonzero eigenvalues occur in +/- pairs.
They are the nonzero singular values of the retained adjacent boundaries, with
both signs and multiplicity. Eigenvectors generally mix adjacent degrees; they
cannot be assigned a single PH degree. On a complete ordinary real chain complex,
the kernel dimension is the sum of real Betti numbers. On a truncated interval,
the kernel belongs to that interval and may differ from source homology.

`boundary_spectrum(j,count,zero_threshold,options,execution)` requests count
largest positive singular modes of one retained B[j], with a<j<=b. The private
real symmetric numerical core is shared with Hodge's block LOBPCG path. It uses
only B^T B actions; no full Dirac or Hodge matrix is formed, and no separate
Dirac dense decomposition is introduced. With Bv=sigma*u and B^T u=sigma*v,
construct w+=(u,v)/sqrt(2) and w-=Gamma*w+, where Gamma=(-1)^absolute_degree.
The chain condition removes other adjacent couplings. Final residuals use the
entire frozen truncated Dirac action, rather than just the selected two degrees.

`DiracBoundarySpectrum` owns basis, all degree offsets and the selected interval.
Singular values are descending; signed values and paired vector columns are
ascending: [-sigma_max,...,-sigma_min,sigma_min,...,sigma_max]. Negative column
j equals Gamma times positive column 2*count-1-j. Repeated singular modes retain
multiplicity. Other boundaries and the kernel require separate requests; this
API does not claim a complete globally sorted Dirac spectrum or exact rank.

Options default to an 8-vector block, 800 updates, seed 73, tolerances 1e-9 and
1,048,576 entries per basis/image or projected buffer. The output N-by-2*count
vector buffer has the same cap. Count must fit the block; zero requests no modes.
Limits, cancellation, nonconvergence and requested values at/below the explicit
zero threshold return errors, without zero padding or a dense fallback.
Tolerances apply in sigma² units to the normal map. The signed residual bound is
absolute_tolerance/|sigma| + relative_tolerance*|sigma|, in Dirac units. Normal
maps square conditioning, so tiny singular modes may be unreliable or rejected.
Residuals do not certify extremal index coverage; a fixed starting block can
miss an invariant eigenspace. The core retains such a numerical counterexample.

This implements the signed reconstruction in [CSB23 section 2.4](../reference/bibliography.md#csb23).
The block iteration follows [Knyazev (2001)](https://doi.org/10.1137/S1064827500366124),
using identity preconditioning. Weighted/complex, persistent, general global
spectrum selection and signal processing remain separate capabilities.

## Independent reproduction and observed results

Run the maintained public Rust checks from a source checkout:

```sh
cargo test --locked --test dirac
cargo test --locked --release --test dirac
cargo run --locked --example dirac
```

Public Rust tests cover hand matrices, analytic edge eigenvectors, Hodge blocks,
parity, both endpoint cuts, empty degrees, noncontiguous IDs, duplicate signed
terms, orientation, invalid integer chains, source errors, stages, every small
work limit, degree-span limits and a 20,000-vertex sparse path whose dense
export is rejected.
Signed-mode regressions additionally verify singular multiplicity, Gamma pairs,
the same truncated interval, cross-boundary orthogonality and numerical failures.

The source-checkout checker `tools/check_dirac.py` needs NumPy and consumes
actual Rust JSONL from `tools/reference/dirac_cocycle.rs`. Its expected signed
boundaries come independently from fixture facets. Run with a prepared Python
environment; ordinary Rust checks and standard-library tool unit tests do not
require these optional numerical packages:

```powershell
cargo build --locked
rustc --edition 2024 tools/reference/dirac_cocycle.rs --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps -o target/dirac_cocycle.exe
& ./target/dirac_cocycle.exe | Set-Content target/dirac-rust.jsonl -Encoding utf8
python tools/check_dirac.py --rust target/dirac-rust.jsonl --output target/dirac-oracle.json --pairs
```

Unix callers can use a binary without the `.exe` suffix and redirect stdout to
the same JSONL path. The checker rejects missing/duplicate fixtures, invalid
graded coordinates, nonfinite values and mismatches with a nonzero exit code.
Optional `--toponetx PATH --gudhi` enables the external checks below. SciPy is
required for TopoNetX and the GUDHI Python binding for `--gudhi`; explicitly
requested unavailable probes fail rather than silently passing.

Download the fixed `simplicial_complex.py` linked below under ignored `target/`,
then provide it through `--toponetx`; Python 3.12 is needed to parse that source.
No third-party source, raw data or generated report belongs in source Git.
The original full probe used Windows 11, Rust 1.98.1, Python 3.12.0rc3,
NumPy 2.4.6, SciPy 1.17.0 and GUDHI 3.11.0; new reports record actual versions.

All eleven fixtures passed: empty/isolated sources, one edge, edge-only and
absent-degree windows, a relabelled path, unfilled/filled triangles, the filled
triangle restricted to C1+C2, and a tetrahedron with/without its solid interior.
Analytic spectra, independent boundary SVD, signed pairing and same-window Hodge
blocks agreed. Absolute spectral tolerance was 1e-11; the largest Frobenius
eigenpair residual was 5.7913e-15 and orthogonality defect was 3.0704e-15.

Nine fixtures also agree with three unchanged AST-extracted methods from
[TopoNetX ce3d414](https://github.com/pyt-team/TopoNetX/blob/ce3d414b571fff95916c4cbc20b8d2c2c9d428e6/toponetx/classes/simplicial_complex.py).
The Git source blob `8ace5fef7ad5e9a44cce02b284578444340bbae9` is checked before
execution. A supplied skeleton adapter replaces class storage; this does not
validate the full package. Empty/isolated fixed-version matrix shapes are
excluded, while truncated cases compare the same principal degree window.
Its rank-zero incidence and `weight` argument do not change this API's ordinary,
standard-inner-product definition.

GUDHI F3 Betti numbers agree with total kernel dimension on seven complete,
torsion-free fixtures. This is finite-field topology evidence, not a GUDHI
Dirac-matrix comparison or a general equivalence with real homology. No such
claim is made for truncated windows.
