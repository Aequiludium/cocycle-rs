# Critical-set comparisons

[Benchmarks](../README.md) / [Reporting rules](../reporting.md)

`critical-sets-native-v1` compares the public sparse critical-set workspace, the
original dense research implementation, and pinned serial Oineus decompositions.
It also measures existing repository persistence paths as separate cost context.
[benchmark_critical_sets.py](../../tools/benchmark_critical_sets.py) orchestrates
native workers on Linux or WSL; Python is outside the measured computation.

## Sources and build

Use Rust 1.91+, Python 3.10+, g++ with C++20, and Boost headers. The controller
builds an ordinary optimized `src/lib.rs` rlib and the standalone
[Rust worker](cocycle.rs), without test counters or production feature changes.
The [C++ worker](oineus.cpp) uses Oineus 0.9.39 at
`e52814a1ffb5b8a81e71ff1e93b4c14194673f0f`. Obtain its source from the
[fixed archive](https://github.com/anigmetov/oineus/archive/e52814a1ffb5b8a81e71ff1e93b4c14194673f0f.zip)
and keep the extracted `include/` and vendored `extern/` tree together.
This adapter does not benchmark the complete `TopologyOptimizer` or autograd.

The dense baseline comes from the local research commit
`167026e7ee8dbaafe8abaf5e325149af6e93f660`, file
`benches/optimization/big_steps.rs`, Git blob
`53bd137d16e122c2ad00a9859188ff317a362472`. It is not guaranteed to exist in a
fresh public clone. The [benchmark-only core](dense_core.rs) bundles its
reduction/critical-set formulas without the old command-line program or tests.
The controller uses this hash-verified excerpt when the Git object is missing.
Supply `--dense-source` with either pinned source to select it explicitly;
the original input is saved in `workers/dense-original.rs`. Generated adaptations change the tie
order to production decreasing colex and precompute the endpoint birth mask.
The reduction and critical-set formulas remain the original dense implementation.

Generated reference modules reuse the repository's test-only vector-XOR F2
reducer, changing visibility for this worker. It is outside the public API.
Source, generated worker, binary and consumed native-header SHA-256 hashes and
compiler commands are retained in each run. The controller records the Oineus
pin and consumed-header hashes; archive provenance must be checked separately.

```sh
python3 tools/benchmark_critical_sets.py \
  --rustc /path/to/rustc \
  --oineus-source target/oineus-pinned \
  --boost-include /path/to/boost/include \
  --output target/benchmarks/commit-<sha12>/critical-sets/run-001
```

Add `--quick --samples 1` for adapter smoke validation. A smoke run is not a
performance ranking. Every output directory must be new. The default is one
discarded warmup plus 12 measured fresh processes per cell, executed serially.
Backends are shuffled with seed 220316749 and rotated by round; position counts
differ by at most one when the round count is not divisible by backend count.
Each worker is pinned to `--cpu` (default 0), with a 60-second timeout and 2-GiB
address-space cap. Frequency and host load are uncontrolled. Do not compile or
test concurrently with timed workers.

For a matched before/after concrete-diagram comparison, supply
`--baseline-build /path/to/prior/build` with a clean committed build's
`metadata.json` and `workers/cocycle`. Its protocol and binary hash are checked;
the baseline metadata is retained with the candidate. The additional
`concrete_baseline` rows run the old worker's `diagram` mode on the same inputs,
with the same timing boundaries and independent output checks.

## Inputs and requested outputs

Fixture seed 220316748 generates two sizes of each family: flag triangulated
grids, fans with delayed triangle values, flag-valued complete 2-skeletons, and
supplied complete 3-skeletons. Full sizes are grid 8/16, fan 64/256, clique
12/20, and 3-skeleton 8/14. All methods use the same canonical simplex order,
F2, f64 values, complete supplied coverage and interval multiplicity.
The 2-skeleton families request H0/H1; the 3-skeleton requests H0/H1/H2.
Critical queries retain zero-length pairs; exported diagrams omit zero bars
and retain essential intervals. No metric, approximation or censoring is used.

| Backend | Requested work |
| --- | --- |
| `critical` | Public `CriticalSetWorkspace`: sparse primal R/V, lazy dual, dimension-cached V transpose, bounded Partial U and merged targets |
| `dense` | Eager primal and dual R/V/U in dimension blocks, then the same endpoint proposals and maximum absolute displacement merge |
| `oineus_partial` | Serial sparse R/V, no clearing, lazy dual, dimension-cached V transpose and upstream bounded U rows |
| `oineus_full` | Same Oineus adapter with eager full U in each constructed decomposition; an ablation of Partial U |
| `reference` | Repository test-only F2 vector-XOR pairing reducer and normalized intervals |
| `diagram` | Concrete simplicial `source.persistence().compute()` and normalized intervals |
| `filtered` | Generic `PersistenceBuilder::from_complex()` validation, boundary reduction and normalized intervals |
| `representatives` | Concrete persistence plus cycles in every requested dimension at the middle simplex's value |
| `flag` | Graph-derived implicit F2 H0/H1 on flag-compatible grid/clique inputs |

Diagram-only, critical-set and representative outputs differ. Context rows
must not be ranked against critical-set workflows as equivalent work.
`flag` excludes delayed fans and supplied 3-skeletons. The Rust workers use
64-bit `usize` indices and explicit coefficients; Oineus uses int32 indices and
implicit F2 sparse columns. Prepared native representations also differ.

Three fixed synthetic critical workloads select pair `(i * 17) % pair_count`,
including possible repeats and zero pairs: `primal-U1` is one death increase;
`mixed-Q8` and `mixed-Q64` rotate death increase/decrease and birth
increase/decrease. Displacements are respectively 0.037, 0.037 and 0.371 of
the filtration span. Warm batches repeat 1024, 128 and 16 times respectively.
The merge selects maximum absolute displacement, retaining the first tie.
These are endpoint-query workloads, not optimizer convergence experiments.

## Timing and validation

Input parsing, source validation, canonical IDs and native input preparation
precede the timers. `build` includes D preparation and primal reduction; dense
also eagerly computes both decompositions and inverses. `cold_query` includes
targets and any lazy dual/transpose preparation. `cleanup` destroys the cold
workspace. `cold_workflow` sums those phases from the same process; pairing
transport is excluded and owned targets remain alive.

A separate workspace is rebuilt and queried once before timing repeated warm
batches. `warm_batch` is average time per batch including each result's
destruction; workspace construction and destruction are excluded. Context
timing includes public computation, normalized intervals and result destruction.
Representative certificate validation is outside its timer.

Every sample is checked for exact interval multiplicity against an independent
Python face/set-XOR reducer. Critical outputs also match exact pair identities,
essential births and the pinned dense merged targets. Requested representative
cycles undergo independent vertex-face XOR closure checks, term dimension/scale
checks and active-interval count checks. These are cycle certificates, not an
independent proof of homology-basis independence. A separate GUDHI correctness
preflight can validate the full fixtures; it does not supply timing rankings.

RSS at prepared input, input high-water mark and absolute whole-process peak
are retained in KiB. The peak covers cold and warm work and validation; it is
not live algorithm allocation or a library memory bound.

## Artifacts

Runs save the sampling plan before executing workers, generated inputs, build
logs, source/header/binary hashes, every warmup/sample, and a full summary with
median/min/max. Mismatched outputs remain in the sample record; failures and
timeouts invalidate that comparison cell. Raw data stays under ignored
`target/`. See the [reporting rules](../reporting.md) for evidence retention and
the [report index](../reports/README.md) for maintained summaries.
