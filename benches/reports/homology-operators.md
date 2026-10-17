# Cocycle and homology-operator: feasible F2 projections

[Benchmarks](../README.md) / [Protocol](../homology-operators.md) /
[Reporting rules](../reporting.md)

## Question and integration outcome

These comparative studies ask whether native Rust can reproduce the feasible
weighted F2 projection workflow, and whether its sparse/packed and request-local
preparation changes reduce the same Rust construction/readout costs. Exact
projection entries, noncycle actions, Betti numbers and selected cycle masses
agree on the tested fixtures. The first port and optimized core are distinct
measured versions; the stage interface in this PR is not retimed.

In the 15-cell Rust before/after study, most construction medians improve, while
general n=513 has a retained 0.12% higher median. Four finite barcode cells have
before/after median ratios 2.93-42.24. These are finite synthetic workloads,
not a universal speedup, complete upstream parity or solver/default admission.
The integrated capability remains feasible projections, positive integer class
geometry and finite coordinate inclusion; global optimality is not implemented.

## Measured revisions and environment

| Group | Immutable identity | Reference and environment |
| --- | --- | --- |
| Initial cross-library, 2026-10-09 | Clean candidate/harness `9582cd6268f4c35074cb7e28a7128a609b1fc6d6`; source/harness SHA256 `3d664b7f591fab5cc29201b918e0a04b097b5d700ba645903dd6467941ea1635` | homology-operator `2d0541c2ef3bf37a580fde5219d55ebe6968215b` |
| Initial Windows group | Run `operator-performance-windows-001`; worker identity digest `38fc3de37f64fe52a308954a1691205183a093e5b0e0bee52c83d5f92baa53ba` | Windows 11 build 26100, Rust 1.98.1/LLVM 22.1.8, CPython 3.12.0rc3, actual PyO3 0.29.3 native extension |
| Initial Linux group | Run `operator-performance-linux-003`; worker identity digest `689907012a29ff7635ff5f9a4fbe95dc73c36142cd380257d6fac181d431824b` | Ubuntu WSL2/glibc 2.39, Rust 1.91.0/LLVM 21.1.2, CPython 3.12.3; no reference native extension |
| Optimized Rust before/after, 2026-10-11 | Before source/controller digest `0c56e930870484ba5fe6c6c6720def1b948dea4475490dfde6b2b22b62bb3449`; after `57e35b17dfc660344878412fda341c14070945ceb2e26602722208523f090d68` | homology-operator `3daff9bd69033d1ee9c5cba5026bfff94df8d6ac`; Windows 11, Rust 1.98.1, CPython 3.12.14 |
| Optimized solve and family suites | `operator-sync-final-heads-001`, `operator-family-final-heads-001`; before/after worker file SHA256 `7fd6eb9f1d5cd9c54032176d0cddf53ca1bd5f3ca7b73b0f39b44992320a74f6` / `9cd7224e8380abe62cc21ff24a8936e2bac6d92b3ab1d3368a6a521fd3601bd4` | Same before/after worker source and Cargo release library builds, edition 2024 and optimized workers |

All groups use an Intel Core Ultra 7 155H with 22 logical CPUs and serial
fresh processes. There is no CPU pinning or frequency/host-load control.
Platforms, compiler versions and reference contracts are not pooled.
The initial reference Python digest is
`33736f974afa55fd50f96ecefd8fc5608f3f5bad3001aa7c50633014cfc4729b`;
the optimized reference digest is
`fb9bbb14933c7b3994ec2cca84915f87f7b38918c33d150095d8116c583addf2`.
The initial native identity digest is
`706be089516e9fb4abf3a955c34ad27acaf873a361a0881aeb77305ec719251c`;
the optimized native identity digest is
`330408c05244e501a5cd011af9515ee1eda8568a3e9a34b9b84c3172950fe292`.
File hashes and basename-plus-bytes controller digests are different identities.

The optimized experiment used immutable uncommitted snapshots on an all-open-PR
hypothetical baseline: #120 at `fe83aa8b37e702be05335b399f44104fa50bb75a`,
#26/#49 withdrawal heads `527d361971deb8129af587a3d4655a088d6b5fd5` /
`b0562bd67be8d8ef1537c9af21222397e2bcf44a`, and draft #56 at
`0a13783359ac9c65f1b8d49c00d0b89a41088c7b`. Complete source snapshots and
integration proof identify that measured code. Before replaces only the three
operator core files with the port at `2ddc9deb1b603da99f337493d70b1e8466a723d3`.
This PR excludes #56 kernels and other contributions. Their draft/no-go state is
unchanged, and this report does not assign snapshot timing to the submitted head.

## Workloads and comparison contract

Each solve call constructs a fresh accepted feasible result from a validated
ordered chain window, including its legality checks and Rust's owned input clone.
Startup, fixture parsing, window validation, explicit P export, serialization,
topology construction and later geometry readouts are outside the timer.
Python/native calls also include dispatch, Python certificate validation,
identity and provenance costs. They do not isolate equivalent arithmetic kernels
or a language effect. The actual native route selects `NativeFactorizedSolver`
with `matrix_free_output=True`; missing native support is a failure, not fallback.

General fixtures have nonzero A/D ranks, singular extra rows/columns and stable
coordinate transformations. Homology-only fixtures have zero A/D with explicit
empty shapes. Sparse fixtures retain sparse AD=0 structure. All use positive
small u64 costs, seed 20261009, shared fixtures and recorded hashes.
The initial study covers general/homology-only n=8/32/65. The optimized study
covers all three families at n=8/32/65/129/513.

Barcode cells are H0 path-merge families with repeated scales and independently
known multiplicities. Only `barcode()` on an accepted family is timed, excluding
input/operator construction, startup and export. They do not time complete PH,
window extraction, arbitrary maps or the new borrowed-stage terminal.

## Validation and coverage

The initial Windows corpus has 46 cases and 92 external comparisons, against
Python and actual native. Linux runs the same 46 cases against Python, without
claiming native-extension or repeated GUDHI acceptance. The optimized corpus has
50 cases, 100 external comparisons and 50 before/after Rust comparisons. Sparse
correctness extends through 511/512/513/1025 coordinates; multiword initial cases
include 129. Every measured solve/barcode output agrees; there are no discarded
mismatches in the ranked groups.

Tiny independent cycle/boundary-set enumeration checks quotient size and class
preservation without production elimination; Rust tests include 108 small
quotient configurations. Finite triangle transport also agrees with the pinned
reference and hand expectations. GUDHI, explicitly over F2, gives Betti stages
(3,0)/(1,1)/(1,0,0) and H1 [1,2): 3.12.0 in the initial Windows run, 3.11.0 in
the optimized comparison. It cannot certify abstract weighted P, selected mass
or global optimality. Its simplex list is checked to reject implicit repair.

The new stage/query interface is tested separately: same owner/certificates,
source-run and inherited-cost rejection, empty/negative/repeated scales,
one shared extraction/solve quota and post-source-drop readout. Actual target D
replays the triangle filling preimage [1]. These tests do not measure latency.
The submitted source has its own Rust/MSRV/Markdown/package checks and exact
operator comparison; hosted status belongs to the PR's exact head, not these
historical source digests. No merge or release acceptance is inferred.

## Measurement and sampling

Each cell/engine has 12 independent measured processes and one retained but
discarded warmup, in seeded position-balanced rounds. Solve sampling uses
20261009; family sampling uses 20261011. All runs are serial, with no concurrent
benchmark/build/test work scheduled during measurements. No outliers are removed.
The initial platforms have 360 measured calls and 30 warmups in total. The
optimized solve group has 720 measured calls and 60 warmups; family readout has
96 measured calls and eight warmups. There is one formal run per platform/group,
without a separate holdout or default-selection admission study.

Times below are milliseconds, medians with observed min-max. Before/after >1
means lower after median. These ratios divide separate medians, not paired
significance estimates. The complete optimized table includes Python/native
absolute times, under their different public-result validation contracts.

## Comparable Rust feasible construction

| Family | n | Before Rust ms | After Rust ms | Before / after | Python ms | Native + Python ms |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| General | 8 | 0.02160 (0.01890-0.07000) | 0.01715 (0.01510-0.02120) | 1.259 | 1.09070 (0.87640-2.07960) | 3.04595 (2.19880-4.39920) |
| Homology only | 8 | 0.02510 (0.02380-0.02910) | 0.01810 (0.01660-0.02070) | 1.387 | 1.08645 (0.91050-1.95280) | 2.18630 (1.89410-3.57460) |
| Sparse | 8 | 0.02185 (0.01840-0.03840) | 0.01645 (0.01540-0.01720) | 1.328 | 0.97425 (0.83080-1.31450) | 2.34735 (1.97790-3.51900) |
| General | 32 | 0.10545 (0.09520-0.20690) | 0.09460 (0.08520-0.12370) | 1.115 | 4.40460 (3.87710-5.71160) | 15.17090 (13.39420-18.27810) |
| Homology only | 32 | 0.09165 (0.07780-0.18950) | 0.07275 (0.06650-0.08550) | 1.260 | 5.24075 (3.95000-7.81610) | 10.20080 (8.72790-13.42470) |
| Sparse | 32 | 0.10490 (0.09790-1.01260) | 0.09825 (0.08590-0.19050) | 1.068 | 4.38645 (3.85750-6.08130) | 14.89895 (12.58000-17.32740) |
| General | 65 | 0.34230 (0.31820-0.55650) | 0.31080 (0.29560-0.37830) | 1.101 | 15.36560 (14.25570-19.66730) | 51.27090 (49.22260-63.08520) |
| Homology only | 65 | 0.24355 (0.22870-0.63660) | 0.20945 (0.19530-0.42940) | 1.163 | 16.19250 (13.44800-18.14560) | 31.74040 (28.51350-36.11650) |
| Sparse | 65 | 0.33000 (0.30920-0.85490) | 0.28505 (0.26420-0.34770) | 1.158 | 16.17160 (12.94120-20.15520) | 54.49410 (48.63860-68.79070) |
| General | 129 | 1.10175 (1.02090-1.57960) | 0.97335 (0.93850-1.29840) | 1.132 | 56.35045 (52.47490-64.62240) | 202.67520 (189.85950-211.56120) |
| Homology only | 129 | 0.70215 (0.63470-0.82850) | 0.66910 (0.55710-0.82200) | 1.049 | 57.13035 (50.75290-74.02340) | 112.45610 (100.38950-123.72860) |
| Sparse | 129 | 1.11240 (1.06420-1.67470) | 0.97850 (0.94570-1.68110) | 1.137 | 55.79350 (51.54860-80.39970) | 192.41565 (179.41650-255.82480) |
| General | 513 | 18.57525 (17.86420-27.87670) | 18.59785 (16.26600-20.95090) | 0.999 | 874.15790 (829.52390-1076.25190) | 3050.23500 (2783.33810-3476.27460) |
| Homology only | 513 | 8.37065 (7.82380-10.70560) | 7.30715 (6.27880-8.65540) | 1.146 | 879.99970 (803.51060-1009.00570) | 1564.18840 (1415.91960-1743.88630) |
| Sparse | 513 | 19.48650 (17.38710-30.46930) | 12.54680 (10.81470-22.82520) | 1.553 | 810.54685 (788.67480-981.80420) | 2787.53700 (2675.18210-3255.90260) |

General n=513 has a 0.12% higher after median, retained alongside improvements.
Overlapping ranges and tiny sample spread do not establish consistent improvement
or regression. Windows RSS/HWM is unavailable, not zero; no memory reduction is
claimed from this group.

## Finite barcode readout

| Coordinates / stages | Before ms (min-max) | After ms (min-max) | Before / after |
| --- | ---: | ---: | ---: |
| 8 / 4 | 0.16990 (0.16270-0.27120) | 0.05790 (0.05140-0.08050) | 2.93 |
| 16 / 8 | 1.64450 (1.54510-2.49830) | 0.20895 (0.19620-0.26210) | 7.87 |
| 32 / 16 | 24.53070 (23.20460-27.64340) | 1.14695 (1.11010-1.24580) | 21.39 |
| 65 / 16 | 168.09365 (160.18960-184.48970) | 3.97945 (3.67240-6.34300) | 42.24 |

All 96 outputs match the hand-derived bars. The rank table and explicit maps
still have quadratic worst-case costs; stretch remains exponential in cycle-space
dimension. These finite measurements establish no infinite-filtration conclusion.

## Earlier cross-library construction

The initial Windows group includes 216 measured outputs with exact agreement:

| Family | n | Rust ms | Python reference ms | Native factorized + Python ms |
| --- | ---: | ---: | ---: | ---: |
| General | 8 | 0.01895 (0.0163-0.0870) | 2.5704 (2.3872-3.7807) | 3.2990 (2.9724-4.7953) |
| Homology only | 8 | 0.02420 (0.0212-0.0458) | 2.2762 (2.1604-4.5334) | 2.5820 (2.2852-3.1298) |
| General | 32 | 0.10260 (0.0944-0.1632) | 34.8401 (32.9207-65.0954) | 44.1593 (40.6643-51.6167) |
| Homology only | 32 | 0.07845 (0.0713-0.0886) | 23.4197 (21.3950-26.6681) | 13.0131 (11.9740-15.9854) |
| General | 65 | 0.32130 (0.3099-0.6072) | 213.9040 (208.4728-217.6843) | 232.5602 (223.4099-293.1999) |
| Homology only | 65 | 0.22330 (0.2088-0.3113) | 123.9824 (116.4421-139.9927) | 40.3777 (37.8540-53.4956) |

Native is slower than Python on the general fixtures and homology-only n=8;
those unfavorable results are retained and their cause was not profiled. Rust's
smaller public-result contract prevents interpreting the ratios as isolated
kernel or universal Rust speedup. Paired per-round Rust/native median ratios and
95% percentile-bootstrap intervals (10,000 resamples, seed 20261009) are:
general n=8, 0.00613 [0.00500,0.00992]; homology n=8,
0.00960 [0.00846,0.01104]; general n=32, 0.00235 [0.00228,0.00241];
homology n=32, 0.00571 [0.00535,0.00655]; general n=65,
0.00137 [0.00133,0.00153]; homology n=65, 0.00540 [0.00526,0.00571].
They describe this run's cells and public-contract difference only.

## Process memory

Initial Linux: milliseconds and KiB, medians; all 144 measured outputs agree.
RSS/HWM is captured before explicit export, with baseline HWM equal to baseline
RSS in these cells:

| Family | n | Rust ms (min-max) | Python ms (min-max) | Rust baseline / peak KiB | Python baseline / peak KiB |
| --- | ---: | ---: | ---: | ---: | ---: |
| General | 8 | 0.01032 (0.00919-0.11392) | 2.23873 (2.01181-2.73924) | 1760 / 1936 | 21648 / 22000 |
| Homology only | 8 | 0.00812 (0.00791-0.03997) | 2.06521 (1.92572-2.40670) | 1760 / 1936 | 21824 / 22000 |
| General | 32 | 0.06540 (0.06133-0.10348) | 24.85120 (22.61611-25.79432) | 1760 / 1936 | 21824 / 22176 |
| Homology only | 32 | 0.04735 (0.04581-0.08274) | 17.20411 (16.62217-20.56748) | 1760 / 1936 | 21824 / 22176 |
| General | 65 | 0.27385 (0.24964-0.35796) | 149.76930 (143.01100-158.74693) | 1936 / 2112 | 21824 / 22352 |
| Homology only | 65 | 0.16934 (0.15654-0.21367) | 94.96489 (87.39460-114.37093) | 1760 / 2112 | 21824 / 22176 |

Absolute RSS includes the Python interpreter and imported modules. Growth is
only hundreds of KiB, at page-level granularity, and includes both a pre-timing
accepted result and the measured construction. It is neither an algorithm
allocation peak nor an isolated workspace reduction. Windows memory is unavailable;
Linux had no reference native extension. Cross-platform ratios are not pooled.

## Interpretation and retained unsuccessful attempts

The measurements support exact feasible-workflow agreement and the shown scoped
Rust before/after observations. They do not establish global optimality, complete
upstream replacement, universal acceleration or default adoption. Rational/float/
arbitrary-precision costs, Python result restoration, general maps, streaming and
historical barcode-basis algorithms are not ported. Final finite-stage survival
has unknown continuation rather than upstream's constant infinite extension.

The first initial GUDHI controller expected two coordinates for a filled triangle
where the oracle returns three; actual topology agreed. The corrected run passed,
and the failed output remains. Initial Linux attempts 001/002 stopped before
measurement due to Windows worktree CRLF/config interpretation; 003 uses process-
local `core.filemode=false` and `core.autocrlf=true`, keeping strict content checks.

Earlier optimized groups are retained without pooling: pre-#56
`operator-sync-windows-002` has sparse n=8 5.98% slower and homology n=129 7.49%
slower by median (digest `5cfa3aec0d358eb0e3d2406984c72654bd031c4622bae1af07af6cae43b13c30`).
`operator-sync-all-prs-001` uses earlier #26/#49 code (digest
`df4183a7a1c7b3b5d57d291aee7551cdd5fc8ff17138e0198d511d2f2133aa36`).
`operator-sync-windows-final` was interrupted and remains unranked. Missing debug
rlib, incomplete snapshots/bench files, CRLF document conflicts and package-link
failures are not promoted to passing evidence. The final optimized controller
metadata records CPython 3.12.14, correcting the earlier prose typo 3.12.13.

Historical Windows tooling retains Unix `resource` and draft #56 slash-path
failures plus five skips; Linux tools passed separately. The initial native ABI
attempt during stage preparation used Python 3.10 against cp312 and failed; a
fresh directory with matching Python 3.12 passed. The publication's first release
build failed from disk exhaustion and was rerun after space recovery. These are
environment/controller failures, not successful solver comparisons.

## Reproduction and evidence

Use the [protocol](../homology-operators.md) and fixed sources above. Correctness
uses a fresh output directory with `--samples 0`. For the historical optimized
experiment, rebuild both complete recorded snapshots with the same worker and run:

```sh
python3 tools/check_homology_operators.py --reference /path/to/pinned-reference \
  --worker /path/to/after-worker --baseline-worker /path/to/before-worker \
  --oracle-python /path/to/gudhi-python --native-dir /path/to/pinned-extension \
  --samples 12 --sizes 8,32,65,129,513 --output target/fresh-operator-run
```

Run `family-bench n stages 1` for (8,4), (16,8), (32,16), (65,16), using 12 fresh
pairs plus warmups and seed 20261011. Reproduction need not produce identical
machine-dependent times; matching source/fixtures, outputs and timing boundaries
is the requirement. New stage extraction/solve/query timing remains unmeasured.

Complete raw results, metadata, source snapshots, binaries and verification logs
are preserved outside source Git. The local archive
`homology-operator-evidence-20261012.zip` has 4,103 checksum-listed payload members,
SHA256 `c80debd8355020662078bb3ba895728ddcd94317f447c7d9d9b192dc505b2f81`,
and compressed size 11,905,835 bytes. Upload of that broader raw package is not
authorized yet; it has no public artifact URL. The maintained report/protocol
and implementation are submitted through the PR. Current CI correctness artifacts,
when produced, have the repository's 14-day retention and cannot replace a durable
archive. Local preservation is not described as public raw-data reproducibility.

Initial Windows summary SHA256:
`99204a96910c8e934db2920aeedbc6c2e13997394000e3bca3a9f222b8dc8c0c`.
Initial Linux summary SHA256:
`fb9557ed324043dcfdc8b15b7e0bdef341049c8770766fd773eae72e31b17102`.
Individual worker/fixture/source identities and all failures remain in their own
group metadata. Reports are source-bound observations, not measurements of a later
documentation commit or release.

## 2026-10-17 bounded geometry-query synchronization

This candidate extends the feasible F2 domain with the geometry
queries added after the reference's PR #93. Earlier measurements and reference
pins above remain historical; this section records correctness, not a speedup.

| Item | Verified scope |
| --- | --- |
| Cocycle base | `30f1e0cd4ea46aae7ec0dbc4dc522ff8a4daac44`, fetched `origin/main` |
| Reference | `a11bf50b8c8320651e9ac015566778d957af9c17`, 23 commits after `3daff9b`; no open reference PRs at inspection |
| New queries | Exact integer class minima and raw-D witnesses; explicitly costed minimum fillings; family-bound historical generators and resolution; one-simplex provenance; unit H1 birth cycles; finite H1/H2 pair and stable volumes |
| Public contracts | [Operator guide](../../docs/guides/homology-operators.md#exact-minima-and-historical-simplex-features) and strict rustdoc; existing selected-P readouts retain their semantics |
| Feasible regression | 50 fixtures, exact P/actions/masses against Python `FeasibleSolver` and actual `NativeFactorizedSolver` |
| Native source | No `native/` source changes between `3daff9b` and `a11bf50`; the retained extension's binary/source hashes remain recorded |
| Independent topology | GUDHI 3.11.0, F2 triangle Betti stages `(3,0)/(1,1)/(1,0,0)` and H1 `[1,2)` |
| Geometry corpus | T1--T5, 11 features in the shared Cocycle source order; exact endpoints, occurrences, histories, locations, minima and birth-cycle agreement |
| Raw geometry control | Enumerate all upper-simplex subsets independently; compare unit and weighted fixed fillings, optimal volumes and epsilon=1 stable volumes, including bitset ties |
| Corpus semantic SHA-256 | `86611887bbc278a9e5c28c6a31b63570bb75c2f8754ed6f247856889b1dfa24c` |
| Corpus file SHA-256 | `babd5ec38b772c94bc0a0383dde55d4a61bb70a6c776057766b00bc0e0c44256` |
| Rust checks | Debug/release/MSRV 1.91 each 243 passed, four existing ignored tests; all-target MSRV check, Clippy, formatting and strict rustdoc passed |
| Usage/package | All eleven runnable examples, affected Markdown doctests, package verification, packaged documentation and packaged operator example passed |
| Tool checks | WSL Ubuntu 105/105 passed; Windows run retains one missing Unix `resource` error and five skips |
| Verification boundary | Local pre-submission results; CI configuration updated for the new pin and `--geometry`. Hosted checks belong to the submitted PR's exact head and are recorded there separately |

The independent Rust class test additionally covers all 64 binary 3-by-2
boundary matrices and all eight input classes with direct coset enumeration.
Hand fixtures check cheap upper fillings, zero versus infeasible results,
integer overflow, noncycles, work failure, corrected dying histories, repeated
endpoint occurrences and family-run mixing.

Reproduce using the [protocol](../homology-operators.md), with `--geometry`.
Raw evidence in this checkout is retained under `target/sync-correctness-final`
(and the earlier `target/sync-correctness-001`), with code/worker/reference
hashes, fixtures and process output. Verification logs are under
`target/sync-verification`; MSRV and platform tool logs remain under `target/`.
These paths are ignored artifacts, outside the package and source history.

This is the bounded native Rust contract. Costs are positive `u64`; stable
scales and epsilon use the exact integer subset within +/-2^53. Specialized
birth-cycle queries accept unit edges; source locations retain original vertex
IDs. Rational/float costs, Euclidean-length certification, periodic coordinates,
Python snapshot schemas and a convex-fan perturbation guarantee are not ported.
General class minimization remains exponential in boundary rank. Nonzero
fillings, birth cycles and affine volume searches retain the reference's
rank/nullity<=16 replay domain, with the positive-cost stable-anchor shortcut.
No partial candidate is labeled optimal after resource failure, and no global
minimum-stretch projection, new timing result, merge or release is claimed.


## 2026-10-17 review follow-up: owned family resources

This addresses the review's credential, repeated-elimination, boundary-test and
resource-evidence concerns. Both checkouts in the operator-reference CI job now
use `persist-credentials: false`. Feature construction maintains incremental
packed column spans, retains only adjacent temporary homology bases, and runs
independent raw A/D verification in a separate included pass. Row transforms
are prepared only for stages with dying histories. Complete output histories
and every eager stage operator remain owned; this is still a small-scale exact
reference interface rather than a general scalable persistence backend.

| Identity / boundary | Recorded value |
| --- | --- |
| Original PR baseline | `a9939fe5833d489a8d3c7f421dc7691fea9bd68d` |
| Measured candidate | `afe06f279f86e436b6c0c0d714cb0c9c7db1e8b0`, clean measured source |
| Protocol / harness | `cocycle-f2-family-resources-v1`; shared worker/controller committed in `051bb4404349f6769aac7d47522d2dbb65ce5d61` |
| Environment | WSL2 Ubuntu, Linux 6.6.87.2 x86_64; Core Ultra 7 155H; Rust 1.91.0; Python 3.12.3 |
| Build | Release rlibs; shared worker compiled twice with `rustc --edition 2024 -C opt-level=3` |
| Algebra / input | Supplied exact F2 H1, unit costs; face-first insertions, including repeated values; no point/metric construction |
| Sampling | 10 fixtures x 2 operations x 2 engines x 10 fresh measured processes = 400; 40 warmups retained/discarded; seeded position-balanced serial order |
| Limits / controls | 30 seconds and 512 MiB address space per worker; all 20 cells completed; affinity CPUs 0--21, host load/frequency uncontrolled |
| Output validation | Every stage Betti sequence and interval list against an independent raw-simplex PH reducer; raw cycle/transport/death equations and exact before/after complete histories |
| Timing | Construction includes all stage solves/validation and retention. Features start with an owned prepared family and include historical output plus independent validation. Fixture parsing/source validation/startup/JSON export excluded |
| Memory | Linux pre-call RSS/HWM, post-call retained RSS and absolute pre-export HWM; source, runtime, prepared families and temporaries are included |

Values below are median [minimum, maximum] over ten processes. Time ratio is
baseline/candidate: above one means the candidate was faster in that cell.
RSS is absolute process peak in MiB, not algorithm allocation or a memory cap.

| Input (vertices / simplices) | Operation | Baseline ms | Candidate ms | Time ratio | Baseline peak MiB | Candidate peak MiB |
| --- | --- | --- | --- | --- | --- | --- |
| fan 4 / 11 | construct | 0.046 [0.044, 0.076] | 0.052 [0.047, 0.085] | 0.89 | 2.062 [2.062, 2.062] | 2.062 [2.062, 2.062] |
| fan 4 / 11 | features | 0.133 [0.123, 0.194] | 0.027 [0.025, 0.031] | 4.93 | 2.234 [2.062, 2.406] | 2.062 [2.062, 2.234] |
| fan 8 / 27 | construct | 0.216 [0.210, 0.223] | 0.342 [0.320, 0.569] | 0.63 | 2.234 [2.234, 2.234] | 2.234 [2.234, 2.234] |
| fan 8 / 27 | features | 0.305 [0.289, 0.462] | 0.160 [0.158, 0.236] | 1.90 | 2.406 [2.234, 2.406] | 2.406 [2.234, 2.406] |
| fan 16 / 59 | construct | 1.320 [1.265, 1.460] | 1.427 [1.394, 1.580] | 0.92 | 2.234 [2.234, 2.234] | 2.234 [2.234, 2.234] |
| fan 16 / 59 | features | 1.724 [1.627, 3.576] | 1.219 [1.167, 2.573] | 1.41 | 2.492 [2.406, 2.578] | 2.406 [2.234, 2.406] |
| fan 32 / 123 | construct | 9.930 [9.533, 10.457] | 10.199 [9.718, 11.483] | 0.97 | 2.578 [2.578, 2.578] | 2.578 [2.578, 2.578] |
| fan 32 / 123 | features | 14.346 [13.457, 20.508] | 9.548 [9.240, 10.193] | 1.50 | 2.750 [2.750, 2.922] | 2.836 [2.750, 2.922] |
| fan 64 / 251 | construct | 76.548 [72.355, 80.399] | 76.622 [73.783, 79.145] | 1.00 | 3.266 [3.266, 3.266] | 3.266 [3.266, 3.266] |
| fan 64 / 251 | features | 148.035 [141.066, 152.017] | 79.953 [75.341, 82.829] | 1.85 | 4.469 [4.469, 4.641] | 4.469 [4.297, 4.469] |
| graph 4 / 10 | construct | 0.039 [0.037, 0.060] | 0.041 [0.040, 0.044] | 0.94 | 2.062 [2.062, 2.062] | 2.062 [2.062, 2.062] |
| graph 4 / 10 | features | 0.140 [0.134, 0.164] | 0.025 [0.025, 0.029] | 5.53 | 2.234 [2.062, 2.234] | 2.062 [2.062, 2.234] |
| graph 6 / 21 | construct | 0.139 [0.138, 0.151] | 0.258 [0.240, 0.352] | 0.54 | 2.234 [2.234, 2.234] | 2.234 [2.234, 2.234] |
| graph 6 / 21 | features | 0.328 [0.284, 0.477] | 0.155 [0.154, 0.164] | 2.12 | 2.406 [2.234, 2.406] | 2.406 [2.234, 2.406] |
| graph 8 / 36 | construct | 0.451 [0.441, 0.522] | 0.566 [0.545, 0.609] | 0.80 | 2.234 [2.234, 2.234] | 2.234 [2.234, 2.234] |
| graph 8 / 36 | features | 1.393 [1.279, 1.862] | 0.764 [0.723, 0.840] | 1.82 | 2.406 [2.234, 2.406] | 2.320 [2.234, 2.406] |
| graph 12 / 78 | construct | 3.756 [3.588, 4.242] | 3.850 [3.755, 4.531] | 0.98 | 2.406 [2.406, 2.406] | 2.406 [2.406, 2.406] |
| graph 12 / 78 | features | 23.182 [21.579, 25.611] | 9.224 [8.050, 11.035] | 2.51 | 2.578 [2.578, 2.750] | 2.578 [2.578, 2.750] |
| graph 16 / 136 | construct | 18.882 [18.297, 20.202] | 18.899 [18.384, 20.852] | 1.00 | 2.578 [2.578, 2.578] | 2.578 [2.578, 2.578] |
| graph 16 / 136 | features | 237.247 [214.590, 250.233] | 58.067 [53.359, 66.231] | 4.09 | 3.781 [3.609, 3.781] | 3.609 [3.438, 3.609] |

Prepared-family RSS is part of the feature-call baseline. For the largest inputs:

| Input | Baseline pre-feature RSS MiB | Candidate pre-feature RSS MiB |
| --- | --- | --- |
| fan 251 simplices | 3.266 [3.266, 3.266] | 3.266 [3.266, 3.266] |
| graph 136 simplices | 2.578 [2.578, 2.578] | 2.578 [2.578, 2.578] |

Feature medians improved by 1.41--5.53x on these ten inputs. At 251 fan simplices,
construction remained about 77 ms and the feature query fell from about 148 to
80 ms. The complete 16-vertex graph has 136 insertions and 105 surviving H1
generators; its feature query fell from about 237 to 58 ms. These observations
cover the named 10--251-stage families only. They do not establish a universal
acceptable size, an asymptotic speedup, H2 performance, production admission or
a general TDA workload guarantee.

Construction has no measured general speedup: small inputs regress, including
27-stage fans (0.216 to 0.342 ms) and 21-stage graphs (0.139 to 0.258 ms).
These unfavorable results are retained. Peak RSS does not uniformly improve;
for example, the 123-stage fan feature call increases from median 2.750 to
2.836 MiB. High-water snapshots can be dominated by preparation or allocator
retention; small quantized differences are not allocation savings.

An earlier clean experiment at `051bb44` extracted full incidence once and
copied degree prefixes. The same 20-cell protocol showed no full-construction
benefit and larger small-input regressions (the 11-stage fan rose from 0.049 to
0.213 ms). That constructor change was withdrawn in `afe06f2`; its metadata,
all raw samples and original binaries remain under the separate ignored
`target/benchmarks/commit-051bb4404349/operator-family/run-001`. It is not the
candidate used in the table above, and no unfavorable run was discarded.

The final correctness suite additionally covers 24 seeded small filtrations in
H0/H1/H2, comparing feature multiplicities with the separate rank-based barcode
and checking original-vertex subset incidence and all histories. Filling nullity
16 completes 65536 candidates; 17 is rejected. T1--T5 covers epsilon=1,2,3,
including repeated scales, H2 and multiple active classes: 15 case/epsilon
combinations and 11 distinct features, with reference and raw-subset agreement.
Final Rust debug/release/MSRV suites passed 246 tests each (four existing ignored);
Linux tools passed 107. GUDHI 3.11.0 triangle and the 50 feasible Python/native
comparisons passed. Hosted CI for subsequent report commits is separate.

Reproduce with the [owned-family protocol](../homology-operators.md#owned-family-and-history-resources).
The final raw run is retained at ignored
`target/benchmarks/commit-afe06f279f86/operator-family/run-001`, including every
fixture, warmup, stdout/stderr and environment/source/worker/harness identity.
This is local evidence, not a public durable raw-data archive. The committed
protocol, exact source revisions, table and source-bound generator support replay.

Final metadata SHA-256: `031df466c6641fe0af4ab722f6274dbc67c3eb86de12688d8704610f77c5cd93`.
Final summary SHA-256: `8bf0cf7f779998a212dd86ef92b61b732eec2a4ef7c07279424a48fe8cb32932`.
