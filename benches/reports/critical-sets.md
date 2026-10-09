# Cocycle, dense reproduction and Oineus: critical-set comparison

[Benchmarks](../README.md) / [Protocol](../optimization/README.md) /
[Reporting rules](../reporting.md)

## Question and integration outcome

This local comparative study evaluates sparse critical sets, the dense paper
reproduction, serial Oineus Partial/full U, and existing persistence paths.
The final PR integrates the opt-in sparse workspace and bounded Partial U,
while preserving main's persistence dispatch. Generic Filtered remains the
shared boundary-reduction optimization direction; implicit Flag/Rips remains
specialized and representatives remain requested output. Dense is benchmark-only.
No runtime backend selector or adaptive threshold is introduced.

Uniform boundary routing and an H2+-only change showed a stable large-skeleton
benefit but failed the matched low-degree gate. Those changes and the attempted
column shortcuts were withdrawn. Concrete retains its existing dispatch pending
a verified routing admission gate. Dense high-dimensional layout and repeated
support opportunities remain follow-up work.

These fixed endpoint batches measure cost, not optimizer convergence. Sparse has
useful cold and memory tradeoffs but loses high-dimensional mixed queries to
dense and Oineus. All formal ranking evidence is local-only. The PR supplies
source and replay instructions; exact-head CI supplies separate correctness and
smoke checks with finite artifact retention.

## Measured revision and environment

| Item | Recorded value |
| --- | --- |
| Associated PR | Pending submission; bound in a documentation-only follow-up |
| Final measured candidate | `bf0069ba4df7068d0b89e75722863e68faa671dc`, clean feature revision |
| Merged main | `cedf5965b774855614b761a1ad311d9dc7992c49`; ordinary merge `a3a234f4d3983e0b9ca0d9f9eed91e19916e96df` |
| Replayable harness | `0202534e3376360f6c5d16cbf2673e2ec4998aa1`; measured controller and adapter bytes unchanged in final candidate |
| Frozen same-source control | `0202534e3376360f6c5d16cbf2673e2ec4998aa1`; all `src/**/*.rs` hashes equal final candidate |
| Production difference from main | Only `src/lib.rs` exports and `src/optimization/`; all existing persistence/algebra bytes match main |
| Production fingerprint | SHA-256 `515835ae587106fed7a88f61ae360ece98ef0ce5ad69c99064cdbb79aa942d62`, sorted compact JSON map of relative `src/**/*.rs` paths to SHA-256 |
| Dense | `167026e7ee8dbaafe8abaf5e325149af6e93f660`, blob `53bd137d16e122c2ad00a9859188ff317a362472`; bundled core SHA-256 `5978aa5a225d4840b76175643af0c824b80aeead7a3686e0c3254aedca55225d` |
| Oineus | 0.9.39, `e52814a1ffb5b8a81e71ff1e93b4c14194673f0f`; all 47 consumed upstream files verified against pinned archive |
| Suite / attempt | `critical-sets-native-v1` / `commit-bf0069ba4df7/critical-sets/run-001` |
| UTC timing interval | `2026-10-09T12:50:41.020989+00:00` to `2026-10-09T12:52:07.085806+00:00`; build precedes timing |
| Platform | Intel Core Ultra 7 155H; Ubuntu 24.04 under WSL2, Linux 6.6.87.2, x86_64, 22 exposed vCPUs |
| Rust | 1.91.1 (`ed61e7d7e242494fb7057f2657300d9e77bb4fcb`), LLVM 21.1.2; edition 2024, `-C opt-level=3` |
| C++ | g++ 13.3.0; C++20, `-O3 -pthread -DOINEUS_DISABLE_ICECREAM`; header dependencies recorded with `-MMD` |
| Controls | Serial workers on CPU 0; 60-second timeout, 2-GiB address-space cap; frequency/host load uncontrolled; no concurrent builds/tests |

The frozen control and final candidate have identical production source.
Their ratios expose build/process variation, not an implementation speedup.
Measured source and later documentation/PR head are distinct; documentation-only
follow-ups preserve this production fingerprint.

Rust uses 64-bit indices and explicit field coefficients in ordered sparse
columns. Oineus uses int32 indices and implicit F2 sparse columns. Its adapter
retains V, disables clearing and reduces serially. Partial/full U share one
binary. The adapter applies the same absolute-displacement merge as Rust; it does
not time `TopologyOptimizer::combine_loss`, autograd or parallel ELZ restoration.
No Oineus code is copied into production Rust. Dense eagerly computes both
primal and dual R/V/U by dimension, even for a single primal query.

## Workloads and validation

All workers receive identical f64 supplied filtrations over F2, with full
coverage and preserved multiplicity. Critical queries include zero-length
pairs; exported diagrams exclude zero bars and retain essential intervals.
There is no approximation, censoring, distance evaluation or metric assumption.
Fixture seed is 220316748. Native preparation differs by backend and precedes
timing; this comparison starts from prepared source representations.

| Fixture | Vertices | Edges | Triangles | Tetrahedra | Requested homology |
| --- | ---: | ---: | ---: | ---: | --- |
| grid-8 | 64 | 161 | 98 | 0 | H0/H1 |
| grid-16 | 256 | 705 | 450 | 0 | H0/H1 |
| delayed-fan-64 | 64 | 126 | 63 | 0 | H0/H1 |
| delayed-fan-256 | 256 | 510 | 255 | 0 | H0/H1 |
| clique-12 | 12 | 66 | 220 | 0 | H0/H1 |
| clique-20 | 20 | 190 | 1140 | 0 | H0/H1 |
| tetra-skeleton-8 | 8 | 28 | 56 | 70 | H0/H1/H2 |
| tetra-skeleton-14 | 14 | 91 | 364 | 1001 | H0/H1/H2 |

`primal-U1` moves one finite death upward by 0.037 of the filtration span.
`mixed-Q8` and `mixed-Q64` rotate death up/down and birth up/down, with 8/64
requests and displacements of 0.037/0.371 of the span. Pair selection is
`(i * 17) % pair_count` and can repeat. Maximum absolute displacement wins
conflicts; equal displacement keeps the first proposal. These are fixed
synthetic requests, not distributions sampled from an application optimizer.

GUDHI 3.12.0 independently matched the interval multisets on all eight full
fixtures in the final postflight against the saved formal-run fixtures. Every
measured and warmup output matched an independent Python vertex-face/set-XOR pairing oracle. Critical
methods also matched exact pair identities, essential births and dense merged
targets. The native smoke covered 132/132 quick cells, and all 98 Linux Python
tool tests passed on unchanged harness bytes. Final Rust checks are listed below.

Representatives request cycles in every homology dimension at the middle
simplex's value. Independent face-XOR checks validate closure, term dimension,
term filtration and the active-interval count. The resulting cycle counts are
20, 62, 21, 81, 1, 1, 26 and 30 in the table's fixture order. These certificates
do not independently prove homology-basis independence.

Independent Rust full inverses and pairing-only block transpositions verify
D V=R, R U=D, V U=U V=I and the lazy support condition. Hand-derived integration
cases cover fields 2/3/65537, signed births, non-flag tetrahedra, cutoffs and source
context. Maximum absolute displacement was checked against pinned Oineus `Max`,
including competing directions and first-proposal ties. Distances are unchanged;
earlier reproduction/Topp checks are historical correctness context.

## Measurement and sampling

The [protocol](../optimization/README.md#timing-and-validation) specifies all
boundaries. Cold time is same-process build plus first targets plus workspace
destruction. Pairing transport is excluded; owned targets remain alive. Warm
time rebuilds a separate workspace, prewarms the same batch, then averages
1024/128/16 repeated batches for U1/Q8/Q64, including each result's destruction.
Warm time excludes workspace build and destruction. Context time includes
public computation, normalized intervals and result destruction; representative
certificate validation is excluded. Input parsing, prepared-source construction,
startup, JSON formatting and controller checks are outside all timers.

The saved plan precedes worker execution. One discarded warmup and 12 fresh
measured processes per cell ran serially, with shuffled/rotated backend order
(seed 220316749). Positions differ by at most one for the six-backend context
groups. All 140 cells completed: 1680/1680 measured workers and 140/140 warmups,
zero mismatches, crashes or timeouts. An additional 24 untallied dense calls
calibrated expected targets. Four incompatible `flag` context cells were
excluded by design, representing 48 unrequested measured calls.

Time tables are milliseconds, median [minimum, maximum], all 12 observations.
The following tables cover the complete U1/Q8/Q64 cold/warm matrix. Phase and
process memory tables follow. No pooled speed score is used.

## primal-U1: one primal U-row query

Cold workflow:

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 0.5937 [0.5127, 0.6624] | 0.5632 [0.5345, 0.7177] | 1.1385 [0.8440, 3.6119] | 1.3415 [1.2441, 1.9478] |
| grid-16 | 1.9151 [1.8030, 2.4577] | 13.2330 [13.0139, 14.1153] | 1.5298 [1.4159, 1.8819] | 2.0323 [1.8545, 4.1389] |
| delayed-fan-64 | 0.4675 [0.4030, 1.1900] | 0.3472 [0.3332, 0.4692] | 0.9757 [0.7810, 1.1765] | 1.2204 [1.1151, 1.7136] |
| delayed-fan-256 | 0.8161 [0.7481, 0.8845] | 4.7522 [4.5578, 4.9980] | 1.1405 [1.0666, 1.4611] | 1.4823 [1.3427, 2.5673] |
| clique-12 | 0.6452 [0.5749, 1.4888] | 0.7632 [0.7036, 1.0127] | 1.2896 [1.0196, 1.4317] | 1.4193 [1.1598, 2.2827] |
| clique-20 | 2.3729 [2.1860, 3.0551] | 28.7516 [25.4009, 38.7245] | 1.8458 [1.5695, 2.5931] | 2.4972 [2.1043, 3.8819] |
| tetra-skeleton-8 | 0.5550 [0.4932, 0.8861] | 0.2614 [0.2497, 0.3502] | 1.0108 [0.8768, 1.4677] | 1.3947 [1.1046, 2.1184] |
| tetra-skeleton-14 | 66.9687 [65.0863, 74.6062] | 88.3752 [86.0502, 97.4864] | 4.3496 [4.1185, 7.2224] | 7.2372 [6.6917, 9.1259] |

Warm repeated batch:

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 0.000445 [0.000421, 0.000556] | 0.000423 [0.000399, 0.000468] | 0.000224 [0.000216, 0.000489] | 0.000072 [0.000069, 0.000109] |
| grid-16 | 0.000209 [0.000203, 0.000230] | 0.001558 [0.001522, 0.001689] | 0.000199 [0.000194, 0.000311] | 0.000058 [0.000056, 0.000150] |
| delayed-fan-64 | 0.000185 [0.000177, 0.000634] | 0.000381 [0.000366, 0.000458] | 0.000184 [0.000178, 0.000247] | 0.000057 [0.000054, 0.000058] |
| delayed-fan-256 | 0.000203 [0.000193, 0.000258] | 0.001131 [0.001074, 0.001459] | 0.000205 [0.000195, 0.000391] | 0.000059 [0.000056, 0.000062] |
| clique-12 | 0.000652 [0.000586, 0.003327] | 0.000259 [0.000247, 0.000279] | 0.000258 [0.000244, 0.000283] | 0.000084 [0.000078, 0.000132] |
| clique-20 | 0.000854 [0.000785, 0.001024] | 0.000649 [0.000571, 0.001442] | 0.000416 [0.000365, 0.000477] | 0.000156 [0.000146, 0.000413] |
| tetra-skeleton-8 | 0.000422 [0.000400, 0.000496] | 0.000152 [0.000147, 0.000205] | 0.000232 [0.000216, 0.000352] | 0.000094 [0.000091, 0.000098] |
| tetra-skeleton-14 | 0.001353 [0.001158, 0.001861] | 0.000449 [0.000405, 0.000700] | 0.000443 [0.000426, 0.000889] | 0.000213 [0.000194, 0.000271] |

## mixed-Q8: eight endpoint requests

Cold workflow:

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 0.9215 [0.7214, 1.1468] | 0.6060 [0.5542, 1.1119] | 1.3810 [1.0900, 2.1971] | 1.6623 [1.4731, 5.7148] |
| grid-16 | 3.0868 [2.9339, 4.1083] | 13.4923 [12.6497, 18.3062] | 2.1231 [1.9662, 5.9321] | 2.9196 [2.7153, 3.3431] |
| delayed-fan-64 | 0.6651 [0.5692, 0.9346] | 0.3558 [0.3379, 0.4404] | 1.0266 [0.9485, 1.2448] | 1.5135 [1.2857, 1.9067] |
| delayed-fan-256 | 1.3254 [1.2321, 4.2135] | 4.5830 [4.4489, 5.7847] | 1.4382 [1.3033, 2.3151] | 1.9986 [1.7947, 2.1859] |
| clique-12 | 0.9193 [0.8301, 1.6805] | 0.7885 [0.7110, 1.2217] | 1.1836 [1.0686, 1.6505] | 1.6608 [1.4034, 2.0088] |
| clique-20 | 4.8313 [4.6621, 6.4729] | 25.8079 [24.6463, 29.2976] | 2.2030 [1.9965, 4.0109] | 2.8414 [2.6475, 3.1871] |
| tetra-skeleton-8 | 0.9243 [0.7777, 1.3886] | 0.2846 [0.2705, 0.3618] | 1.3020 [1.0295, 1.8206] | 1.7760 [1.4270, 2.3745] |
| tetra-skeleton-14 | 114.2323 [110.8119, 123.1218] | 89.1358 [85.8162, 124.4143] | 5.9741 [5.7136, 7.4210] | 9.5518 [9.1410, 11.0279] |

Warm repeated batch:

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 0.001766 [0.001600, 0.002576] | 0.003148 [0.003000, 0.005135] | 0.001068 [0.000974, 0.002034] | 0.000478 [0.000453, 0.000827] |
| grid-16 | 0.001826 [0.001687, 0.005141] | 0.009850 [0.009240, 0.017958] | 0.001269 [0.001122, 0.001571] | 0.000652 [0.000589, 0.000804] |
| delayed-fan-64 | 0.001098 [0.001032, 0.001161] | 0.002464 [0.002313, 0.002949] | 0.000854 [0.000817, 0.001705] | 0.000391 [0.000373, 0.000422] |
| delayed-fan-256 | 0.001119 [0.001035, 0.002705] | 0.008024 [0.007464, 0.009122] | 0.000937 [0.000920, 0.001167] | 0.000433 [0.000413, 0.001198] |
| clique-12 | 0.002115 [0.002093, 0.002755] | 0.002477 [0.002330, 0.003806] | 0.001086 [0.001034, 0.001916] | 0.000456 [0.000424, 0.001086] |
| clique-20 | 0.002149 [0.002090, 0.002907] | 0.008015 [0.007294, 0.012667] | 0.001192 [0.001129, 0.001621] | 0.000572 [0.000556, 0.001226] |
| tetra-skeleton-8 | 0.002732 [0.002215, 0.005284] | 0.001766 [0.001642, 0.002119] | 0.001131 [0.000998, 0.001655] | 0.000574 [0.000538, 0.001025] |
| tetra-skeleton-14 | 0.020608 [0.018811, 0.027797] | 0.009080 [0.008105, 0.014294] | 0.006090 [0.005821, 0.010317] | 0.004602 [0.004297, 0.006825] |

## mixed-Q64: 64 wide endpoint requests

Cold workflow:

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 0.8982 [0.8122, 1.9033] | 0.6237 [0.5976, 0.7978] | 1.2487 [1.0983, 1.4255] | 1.7416 [1.4673, 3.3948] |
| grid-16 | 3.8259 [3.7033, 4.0196] | 14.3022 [13.6907, 15.1925] | 2.2817 [2.0812, 2.7920] | 3.2724 [2.9172, 5.2749] |
| delayed-fan-64 | 0.6900 [0.5536, 0.8896] | 0.4193 [0.3712, 0.5984] | 1.1150 [0.9990, 1.7450] | 1.4719 [1.4147, 2.6691] |
| delayed-fan-256 | 1.6014 [1.4388, 2.0064] | 4.9268 [4.4196, 6.4526] | 1.5989 [1.2943, 2.0530] | 2.1284 [1.7284, 2.7902] |
| clique-12 | 1.1196 [1.0741, 1.7310] | 0.7647 [0.7379, 0.8914] | 1.2432 [1.0713, 1.5215] | 1.6550 [1.4981, 2.0262] |
| clique-20 | 5.0014 [4.7295, 5.4364] | 25.7799 [24.8766, 26.9826] | 2.2499 [2.0534, 2.7784] | 2.7142 [2.6007, 3.1961] |
| tetra-skeleton-8 | 0.9735 [0.8921, 1.1199] | 0.3200 [0.2980, 0.3791] | 1.1667 [0.9928, 1.6242] | 1.6890 [1.3270, 4.2532] |
| tetra-skeleton-14 | 122.2357 [117.3399, 145.4995] | 89.8736 [86.0745, 109.9958] | 6.7189 [6.4857, 7.9995] | 9.8042 [9.0512, 15.4430] |

Warm repeated batch:

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 0.019503 [0.018730, 0.050356] | 0.032752 [0.030724, 0.035802] | 0.011667 [0.011027, 0.012394] | 0.006991 [0.006800, 0.009937] |
| grid-16 | 0.054014 [0.051202, 0.062858] | 0.116214 [0.105732, 0.192317] | 0.022682 [0.021119, 0.032247] | 0.014662 [0.014137, 0.016467] |
| delayed-fan-64 | 0.013540 [0.012580, 0.017023] | 0.024153 [0.022870, 0.030704] | 0.010006 [0.009513, 0.011414] | 0.006391 [0.005973, 0.014848] |
| delayed-fan-256 | 0.029183 [0.026545, 0.041081] | 0.082303 [0.073723, 0.092348] | 0.018108 [0.015866, 0.022816] | 0.014203 [0.012622, 0.014944] |
| clique-12 | 0.090792 [0.088160, 0.111942] | 0.045337 [0.042694, 0.047692] | 0.039934 [0.037428, 0.045181] | 0.019839 [0.018377, 0.040788] |
| clique-20 | 0.165247 [0.149645, 0.193093] | 0.141637 [0.135203, 0.219250] | 0.072719 [0.069185, 0.154523] | 0.043916 [0.042473, 0.049633] |
| tetra-skeleton-8 | 0.055310 [0.050366, 0.062092] | 0.024304 [0.022244, 0.027790] | 0.019280 [0.017984, 0.029066] | 0.010902 [0.010151, 0.011390] |
| tetra-skeleton-14 | 1.615127 [1.534224, 1.765494] | 0.286462 [0.252809, 0.422550] | 0.303079 [0.287268, 0.510884] | 0.188085 [0.182558, 0.208221] |

## Repository persistence paths

These methods request different work from critical queries. Concrete, Filtered
and admitted Flag request equivalent diagrams; cycles add payload. The test-only
F2 reference lacks the complete public policy/context contract and is not an
interchangeable public speedup denominator.

| Fixture | Reference F2 | Concrete | Filtered | With cycles | Implicit Flag | Frozen control |
| --- | --- | --- | --- | --- | --- | --- |
| grid-8 | 0.0380 [0.0367, 0.0832] | 0.2681 [0.2167, 0.4458] | 0.3038 [0.2724, 0.5149] | 0.3969 [0.3601, 0.6046] | 0.2054 [0.1729, 0.2636] | 0.2823 [0.2436, 0.3506] |
| grid-16 | 0.1695 [0.1640, 0.1928] | 0.7094 [0.6191, 3.3653] | 0.8287 [0.7579, 1.7926] | 1.5343 [1.3884, 4.2482] | 0.4281 [0.4087, 0.5384] | 0.7006 [0.6469, 0.9162] |
| delayed-fan-64 | 0.0255 [0.0244, 0.0267] | 0.1014 [0.0914, 0.1156] | 0.1170 [0.1046, 0.2171] | 0.1582 [0.1508, 0.2152] | excluded | 0.1016 [0.0932, 0.1084] |
| delayed-fan-256 | 0.0768 [0.0701, 0.0924] | 0.3094 [0.2790, 0.6373] | 0.3422 [0.3250, 0.7786] | 0.5017 [0.4577, 0.6735] | excluded | 0.2984 [0.2676, 1.2966] |
| clique-12 | 0.0647 [0.0584, 0.0794] | 0.3192 [0.2839, 0.4635] | 0.3872 [0.3443, 0.5286] | 0.4392 [0.4108, 0.6117] | 0.1836 [0.1647, 0.2203] | 0.3223 [0.2752, 0.4039] |
| clique-20 | 0.3780 [0.3704, 0.4231] | 0.7841 [0.7020, 1.6615] | 1.5332 [1.4067, 1.6807] | 2.0858 [2.0152, 2.2473] | 0.2387 [0.2021, 0.3305] | 0.7427 [0.7077, 1.0662] |
| tetra-skeleton-8 | 0.0456 [0.0444, 0.0504] | 0.3387 [0.3152, 0.4043] | 0.3653 [0.3339, 1.2808] | 0.4150 [0.3896, 0.4528] | excluded | 0.3439 [0.3066, 0.4202] |
| tetra-skeleton-14 | 4.8506 [4.6362, 5.6905] | 63.3886 [59.5949, 74.0708] | 25.6483 [24.8352, 31.2305] | 70.5093 [65.9379, 78.8448] | excluded | 61.0769 [57.6801, 67.5997] |

Flag admits only grids/cliques; delayed triangle values and these supplied H2
workloads require another path. Generic validation and trusted concrete incidence
reading have different adapter costs. Results belong to this standalone build.
Current Concrete has a lower median than Filtered on seven of eight rows.
Filtered wins the large tetrahedron skeleton: 25.6483 versus 63.3886 ms. Flag
is lower than Concrete on all four admitted fixtures. Filtered remains the
generic optimization owner without forcing this performance crossover into a
default dispatch change.

The earlier `ceea0f4e47e4c296902549e459432eaa3bbb2684` study predates Phase-3 main:
its Filtered 7/8 result and 25.814 versus 61.024 ms skeleton row motivated the
trials but are not fresh main-admission evidence.

Same-source control, frozen/current Concrete ratio, median [min, max]. Above one
means the current build is faster; wins count paired rounds. Identical source
makes these sensitivity controls, not code-level improvement evidence.

| Fixture | Frozen/current ratio | Current wins |
| --- | --- | --- |
| grid-8 | 1.018 [0.555, 1.304] | 7/12 |
| grid-16 | 1.000 [0.192, 1.354] | 6/12 |
| delayed-fan-64 | 1.019 [0.857, 1.050] | 7/12 |
| delayed-fan-256 | 0.957 [0.592, 4.210] | 5/12 |
| clique-12 | 1.002 [0.768, 1.338] | 6/12 |
| clique-20 | 0.968 [0.426, 1.519] | 5/12 |
| tetra-skeleton-8 | 0.993 [0.837, 1.235] | 6/12 |
| tetra-skeleton-14 | 0.960 [0.848, 1.047] | 3/12 |

Control median ratios span 0.957-1.019. Retained outliers show why small results
need caution. The clear low-degree losses of the rejected default-routing trials
are not reproduced after restoring the exact main implementation.

## Phase costs

Each entry is median build / cold-query / cleanup in milliseconds. Phase medians
need not sum to workflow medians. Cold query includes lazy dual/transpose work;
U1 exposes eager state cost and mixed requests expose both decompositions.

### primal-U1 phases

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 0.4966 / 0.0636 / 0.0198 | 0.5418 / 0.0019 / 0.0149 | 0.9991 / 0.1314 / 0.0041 | 1.3350 / 0.0013 / 0.0047 |
| grid-16 | 1.3837 / 0.4337 / 0.0917 | 13.0426 / 0.0056 / 0.1725 | 1.3467 / 0.1637 / 0.0273 | 2.0025 / 0.0014 / 0.0288 |
| delayed-fan-64 | 0.4222 / 0.0318 / 0.0131 | 0.3344 / 0.0015 / 0.0109 | 0.8800 / 0.1007 / 0.0021 | 1.2144 / 0.0014 / 0.0046 |
| delayed-fan-256 | 0.6384 / 0.1254 / 0.0455 | 4.6410 / 0.0036 / 0.1103 | 1.0264 / 0.1137 / 0.0084 | 1.4651 / 0.0013 / 0.0158 |
| clique-12 | 0.5880 / 0.0297 / 0.0180 | 0.7386 / 0.0015 / 0.0235 | 1.1716 / 0.1196 / 0.0039 | 1.4117 / 0.0015 / 0.0059 |
| clique-20 | 2.2050 / 0.0881 / 0.0909 | 28.5352 / 0.0056 / 0.2133 | 1.6728 / 0.1609 / 0.0226 | 2.4624 / 0.0024 / 0.0325 |
| tetra-skeleton-8 | 0.5273 / 0.0127 / 0.0115 | 0.2532 / 0.0013 / 0.0059 | 0.8832 / 0.1032 / 0.0025 | 1.3900 / 0.0014 / 0.0034 |
| tetra-skeleton-14 | 66.5414 / 0.0535 / 0.3912 | 88.1599 / 0.0038 / 0.2214 | 4.1309 / 0.1650 / 0.0393 | 7.1814 / 0.0025 / 0.0541 |

### mixed-Q8 phases

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 0.5517 / 0.2987 / 0.0461 | 0.5829 / 0.0071 / 0.0145 | 1.0991 / 0.2488 / 0.0075 | 1.3905 / 0.2678 / 0.0181 |
| grid-16 | 1.4193 / 1.4691 / 0.1960 | 13.2896 / 0.0235 / 0.1647 | 1.3191 / 0.6894 / 0.0501 | 1.9767 / 0.8445 / 0.0728 |
| delayed-fan-64 | 0.4489 / 0.1875 / 0.0310 | 0.3396 / 0.0053 / 0.0112 | 0.8592 / 0.1740 / 0.0044 | 1.2721 / 0.2011 / 0.0056 |
| delayed-fan-256 | 0.6179 / 0.6138 / 0.0961 | 4.4700 / 0.0170 / 0.1071 | 1.0949 / 0.3393 / 0.0146 | 1.5415 / 0.4142 / 0.0264 |
| clique-12 | 0.5364 / 0.3126 / 0.0293 | 0.7529 / 0.0057 / 0.0229 | 0.9730 / 0.2024 / 0.0070 | 1.4158 / 0.2342 / 0.0169 |
| clique-20 | 2.0692 / 2.5366 / 0.2015 | 25.5283 / 0.0184 / 0.1893 | 1.5433 / 0.5562 / 0.0444 | 2.2088 / 0.5268 / 0.0789 |
| tetra-skeleton-8 | 0.6269 / 0.2525 / 0.0245 | 0.2732 / 0.0052 / 0.0066 | 1.0785 / 0.2224 / 0.0049 | 1.5099 / 0.2272 / 0.0105 |
| tetra-skeleton-14 | 67.8598 / 45.3765 / 0.7611 | 88.9032 / 0.0247 / 0.2130 | 4.0744 / 1.8067 / 0.0819 | 7.1499 / 2.2202 / 0.1259 |

### mixed-Q64 phases

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 0.5177 / 0.3358 / 0.0465 | 0.5606 / 0.0469 / 0.0143 | 0.9537 / 0.2656 / 0.0074 | 1.4367 / 0.2908 / 0.0196 |
| grid-16 | 1.4360 / 2.1437 / 0.2605 | 13.9120 / 0.2098 / 0.1763 | 1.3682 / 0.8413 / 0.0539 | 2.2182 / 0.9163 / 0.0722 |
| delayed-fan-64 | 0.4523 / 0.2155 / 0.0317 | 0.3637 / 0.0379 / 0.0119 | 0.9203 / 0.1956 / 0.0045 | 1.2559 / 0.2229 / 0.0065 |
| delayed-fan-256 | 0.6489 / 0.8156 / 0.1168 | 4.6692 / 0.1184 / 0.1027 | 1.1071 / 0.4371 / 0.0191 | 1.6558 / 0.4417 / 0.0283 |
| clique-12 | 0.5613 / 0.5273 / 0.0373 | 0.6789 / 0.0686 / 0.0165 | 0.9507 / 0.2933 / 0.0088 | 1.3429 / 0.2660 / 0.0160 |
| clique-20 | 2.1086 / 2.6723 / 0.1925 | 25.3663 / 0.2031 / 0.1800 | 1.4986 / 0.6800 / 0.0433 | 2.0786 / 0.5838 / 0.0506 |
| tetra-skeleton-8 | 0.5676 / 0.3702 / 0.0300 | 0.2683 / 0.0454 / 0.0061 | 0.9290 / 0.2414 / 0.0058 | 1.4510 / 0.2393 / 0.0074 |
| tetra-skeleton-14 | 68.9297 / 50.8950 / 1.2039 | 89.3617 / 0.3572 / 0.1185 | 4.1297 / 2.4566 / 0.0942 | 7.1056 / 2.3976 / 0.1167 |

## Process memory

Absolute whole-process peak RSS, MiB, median [min, max]. Peaks include cold/warm
work, validation, outputs and allocator retention. They are not live algorithm
allocations or a library memory bound. Index/storage widths and prepared inputs
prevent interpreting peaks as equal ownership.

### primal-U1 peak

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 2.750 [2.750, 2.922] | 2.922 [2.922, 2.922] | 3.953 [3.781, 3.953] | 4.125 [4.125, 4.297] |
| grid-16 | 3.781 [3.781, 3.781] | 7.496 [7.453, 7.500] | 4.469 [4.297, 4.469] | 4.641 [4.641, 4.812] |
| delayed-fan-64 | 2.578 [2.578, 2.750] | 2.750 [2.750, 2.750] | 3.953 [3.781, 3.953] | 3.953 [3.953, 4.125] |
| delayed-fan-256 | 3.266 [3.266, 3.438] | 5.172 [5.172, 5.203] | 4.297 [4.125, 4.297] | 4.555 [4.469, 4.641] |
| clique-12 | 2.750 [2.750, 2.922] | 2.922 [2.922, 2.922] | 3.867 [3.781, 3.953] | 4.125 [4.125, 4.297] |
| clique-20 | 3.609 [3.609, 3.781] | 8.938 [8.938, 8.938] | 4.469 [4.297, 4.469] | 4.641 [4.641, 4.812] |
| tetra-skeleton-8 | 2.578 [2.578, 2.578] | 2.578 [2.406, 2.578] | 3.953 [3.781, 3.953] | 3.953 [3.953, 4.125] |
| tetra-skeleton-14 | 4.812 [4.812, 4.812] | 8.629 [8.629, 8.633] | 4.641 [4.469, 4.641] | 5.328 [5.328, 5.500] |

### mixed-Q8 peak

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 2.922 [2.922, 3.094] | 2.922 [2.922, 2.922] | 4.125 [4.125, 4.125] | 4.125 [4.125, 4.297] |
| grid-16 | 4.297 [4.297, 4.469] | 7.496 [7.496, 7.500] | 4.812 [4.641, 4.812] | 4.984 [4.984, 5.156] |
| delayed-fan-64 | 2.750 [2.750, 2.922] | 2.750 [2.750, 2.750] | 4.125 [3.953, 4.125] | 4.125 [4.125, 4.297] |
| delayed-fan-256 | 3.609 [3.609, 3.781] | 5.172 [5.172, 5.176] | 4.469 [4.297, 4.469] | 4.641 [4.641, 4.812] |
| clique-12 | 2.750 [2.750, 2.922] | 2.922 [2.922, 2.922] | 4.125 [3.953, 4.125] | 4.125 [4.125, 4.297] |
| clique-20 | 4.469 [4.469, 4.641] | 8.938 [8.938, 8.992] | 4.812 [4.641, 4.812] | 5.070 [4.984, 5.156] |
| tetra-skeleton-8 | 2.750 [2.750, 2.922] | 2.578 [2.578, 2.578] | 3.953 [3.781, 3.953] | 4.125 [4.125, 4.297] |
| tetra-skeleton-14 | 6.359 [6.359, 6.359] | 8.801 [8.691, 8.801] | 5.328 [5.156, 5.328] | 5.672 [5.672, 5.844] |

### mixed-Q64 peak

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 2.922 [2.922, 3.094] | 2.922 [2.922, 2.922] | 4.125 [3.953, 4.125] | 4.211 [4.125, 4.297] |
| grid-16 | 4.641 [4.641, 4.812] | 7.496 [7.496, 7.500] | 4.812 [4.641, 4.812] | 4.984 [4.984, 5.156] |
| delayed-fan-64 | 2.750 [2.750, 2.922] | 2.750 [2.750, 2.750] | 4.125 [3.953, 4.125] | 4.297 [4.125, 4.297] |
| delayed-fan-256 | 3.781 [3.781, 3.953] | 5.344 [5.227, 5.348] | 4.469 [4.297, 4.469] | 4.641 [4.641, 4.812] |
| clique-12 | 2.922 [2.922, 2.922] | 2.922 [2.922, 2.922] | 4.125 [4.125, 4.125] | 4.125 [4.125, 4.297] |
| clique-20 | 4.641 [4.641, 4.641] | 8.938 [8.938, 8.938] | 4.812 [4.641, 4.812] | 4.984 [4.984, 5.156] |
| tetra-skeleton-8 | 2.750 [2.750, 2.922] | 2.578 [2.406, 2.578] | 4.125 [3.953, 4.125] | 4.125 [4.125, 4.297] |
| tetra-skeleton-14 | 7.906 [7.906, 8.078] | 8.766 [8.766, 8.766] | 5.672 [5.500, 5.672] | 5.672 [5.672, 5.844] |

### Q64 prepared input

Median RSS / input high-water mark, MiB, sampled before cold construction:

| Fixture | Sparse | Dense | Oineus Partial U | Oineus full U |
| --- | --- | --- | --- | --- |
| grid-8 | 2.406 / 2.406 | 2.406 / 2.406 | 3.438 / 3.438 | 3.438 / 3.609 |
| grid-16 | 2.750 / 2.750 | 2.750 / 2.750 | 3.609 / 3.609 | 3.609 / 3.609 |
| delayed-fan-64 | 2.406 / 2.406 | 2.406 / 2.406 | 3.438 / 3.523 | 3.438 / 3.609 |
| delayed-fan-256 | 2.578 / 2.578 | 2.750 / 2.750 | 3.609 / 3.609 | 3.609 / 3.609 |
| clique-12 | 2.406 / 2.406 | 2.406 / 2.406 | 3.438 / 3.523 | 3.438 / 3.438 |
| clique-20 | 2.750 / 2.750 | 2.750 / 2.750 | 3.609 / 3.695 | 3.609 / 3.609 |
| tetra-skeleton-8 | 2.406 / 2.406 | 2.406 / 2.406 | 3.438 / 3.609 | 3.438 / 3.438 |
| tetra-skeleton-14 | 2.750 / 2.750 | 2.922 / 2.922 | 3.609 / 3.609 | 3.609 / 3.609 |

## Sparse/dense and Partial/full U tradeoffs

Q64 paired sparse/dense ratios, median [min, max]; below one means sparse is
faster. Wins count sparse-faster process rounds out of 12.

| Fixture | Cold sparse/dense | Cold wins | Warm sparse/dense | Warm wins |
| --- | --- | --- | --- | --- |
| grid-8 | 1.457 [1.204, 2.791] | 0/12 | 0.621 [0.524, 1.579] | 11/12 |
| grid-16 | 0.270 [0.250, 0.289] | 12/12 | 0.465 [0.279, 0.584] | 12/12 |
| delayed-fan-64 | 1.641 [1.305, 1.916] | 0/12 | 0.553 [0.490, 0.732] | 12/12 |
| delayed-fan-256 | 0.328 [0.249, 0.402] | 12/12 | 0.363 [0.327, 0.455] | 12/12 |
| clique-12 | 1.461 [1.358, 2.277] | 0/12 | 2.027 [1.850, 2.572] | 0/12 |
| clique-20 | 0.190 [0.186, 0.214] | 12/12 | 1.153 [0.741, 1.378] | 1/12 |
| tetra-skeleton-8 | 3.021 [2.583, 3.295] | 0/12 | 2.270 [2.054, 2.493] | 0/12 |
| tetra-skeleton-14 | 1.365 [1.147, 1.595] | 0/12 | 5.682 [4.112, 6.561] | 0/12 |

On tetra-skeleton-14 Q64, sparse costs 122.2357 ms cold versus dense 89.8736
ms and Oineus Partial U 6.7189 ms; warm cost is 1.615127 versus 0.286462 and
0.303079 ms. All eight sparse warm medians exceed Oineus Partial U. These
regressions remain explicit limits of the public sparse workspace.

The larger grid/fan/clique cold costs and process peaks justify sparse as a
scalable base. Tetrahedron skeletons and clique repeated queries expose dense
opportunities. Small overlapping differences do not establish a selector
threshold. Phase results point to decomposition construction and support solving;
no profiler attribution to one data structure is made.

Oineus full U spends construction work to reduce repeated support lookup cost.
Partial U avoids constructing the full inverse but solves rows repeatedly.
Ratios of cold medians full/Partial and warm medians Partial/full follow for
all workloads. Neither strategy is a universal speedup.

| Fixture | U1 cold full/Partial | U1 warm Partial/full | Q8 cold full/Partial | Q8 warm Partial/full | Q64 cold full/Partial | Q64 warm Partial/full |
| --- | --- | --- | --- | --- | --- | --- |
| grid-8 | 1.178 | 3.122 | 1.204 | 2.232 | 1.395 | 1.669 |
| grid-16 | 1.329 | 3.451 | 1.375 | 1.946 | 1.434 | 1.547 |
| delayed-fan-64 | 1.251 | 3.221 | 1.474 | 2.185 | 1.320 | 1.566 |
| delayed-fan-256 | 1.300 | 3.463 | 1.390 | 2.165 | 1.331 | 1.275 |
| clique-12 | 1.101 | 3.086 | 1.403 | 2.378 | 1.331 | 2.013 |
| clique-20 | 1.353 | 2.675 | 1.290 | 2.084 | 1.206 | 1.656 |
| tetra-skeleton-8 | 1.380 | 2.450 | 1.364 | 1.970 | 1.448 | 1.768 |
| tetra-skeleton-14 | 1.664 | 2.077 | 1.599 | 1.323 | 1.459 | 1.611 |

## Rejected current-main integration trials

Each trial used the same fixtures, 140 cells, 12 fresh processes per cell,
fixed `0202534` control and independent output checks. All 1680 measured outputs
and eight GUDHI fixtures passed. Large-skeleton gains were insufficient to admit
default changes with repeated low-degree losses. The combined column shortcuts
were withdrawn; diagnostics do not isolate one container or compiler cause.
Full trial outputs remain separate and are not relabeled as final evidence.

| Trial | Measured SHA | Manifest SHA-256 |
| --- | --- | --- |
| Uniform boundary | `7f0b081010d866deb44db39476b3e7b2c0fc5e19` | `4df05f984b106fce0aaa3e0451b2dcc031670c789f3b008e011ac3c0bf552f02` |
| Uniform + shortcuts | `620ec7594b8003d577ccddf2f1bd3245657acdaf` | `cf72525a8adec9fc03bc7ff5d7c5a2fcd9cbae68dc66cc8b00fed32afeb411bb` |
| H2+ + shortcuts | `097d36621e4250dfb1019b54494785a2c6a44882` | `77d4c73228e94be47c0fc63850ec3ca5e9c2e311569cb24902e72f3ae1b551ac` |
| H2+ original algebra | `2b1a48e43d9106c584494c831f7af1958cf03713` | `39f83fb10b097fb9e8e91c05e368d8ca06a36b1f3f8a5ff1dcb1d9a6c6899f80` |

| Fixture | Uniform boundary ratio; wins | Uniform + shortcuts ratio; wins | H2+ + shortcuts ratio; wins | H2+ original algebra ratio; wins |
| --- | --- | --- | --- | --- |
| grid-8 | 0.721 [0.418, 1.367]; 2/12 | 0.935 [0.667, 1.339]; 4/12 | 0.622 [0.480, 0.855]; 0/12 | 0.503 [0.428, 0.754]; 0/12 |
| grid-16 | 0.888 [0.657, 1.035]; 3/12 | 1.033 [0.905, 1.431]; 8/12 | 0.822 [0.637, 1.060]; 1/12 | 0.685 [0.636, 0.889]; 0/12 |
| delayed-fan-64 | 0.475 [0.314, 0.625]; 0/12 | 1.179 [0.898, 1.293]; 11/12 | 0.237 [0.193, 0.274]; 0/12 | 0.442 [0.161, 0.689]; 0/12 |
| delayed-fan-256 | 0.730 [0.566, 0.882]; 0/12 | 1.123 [0.906, 1.541]; 10/12 | 0.467 [0.370, 2.385]; 1/12 | 0.728 [0.384, 1.905]; 1/12 |
| clique-12 | 0.715 [0.608, 0.948]; 0/12 | 0.849 [0.661, 1.102]; 2/12 | 0.757 [0.559, 1.109]; 1/12 | 0.560 [0.179, 0.715]; 0/12 |
| clique-20 | 0.528 [0.451, 0.665]; 0/12 | 0.620 [0.493, 0.728]; 0/12 | 0.977 [0.682, 1.193]; 4/12 | 0.645 [0.369, 1.121]; 1/12 |
| tetra-skeleton-8 | 0.782 [0.564, 0.914]; 0/12 | 0.935 [0.734, 1.297]; 4/12 | 0.704 [0.578, 1.327]; 1/12 | 0.914 [0.633, 1.788]; 5/12 |
| tetra-skeleton-14 | 2.379 [2.112, 2.578]; 12/12 | 2.767 [2.313, 3.081]; 12/12 | 2.722 [2.456, 2.876]; 12/12 | 2.349 [2.199, 2.642]; 12/12 |

Ratios are frozen/current Concrete as above. Single-sample dirty pilots and
three copied-source column ablations remain diagnostic only, outside rankings.

## Merge scope and follow-up

Current main is ordinary-merged; no family of user-selectable engines is added:

| Capability | Disposition | Final code outcome |
| --- | --- | --- |
| Generic Filtered | Retain; priority optimization owner | Existing validated boundary core; automatic Concrete rerouting deferred after failed gates |
| Concrete interface | Retain compatibility | All production dispatch bytes match main; established boundary/coface owners reused |
| Implicit Flag/Rips | Retain specialization | Existing F2 low-degree acceleration and high-degree continuation unchanged |
| Representatives | Retain separate output | Transform storage only for requested witnesses |
| Critical sets | Integrate sparse | Primal R/V, lazy dual, dimension V transpose and bounded Partial U; no full inverse cache |
| Dense critical sets | Benchmark only | Hash-verified excerpt excluded from package/runtime; dimension blocks/repeated-query ideas deferred |
| Internal F2 | Test/reference only | No public engine option |
| Uniform/H2+ routing and column shortcuts | Reject this integration | All attempted changes withdrawn from final production diff |
| Adaptive selection, parallel ELZ, autograd | Defer | No established threshold or completed optimizer workflow |

No existing caller migration is required. Use `CriticalSetWorkspace::new` or
`new_with` on the complete frozen F2 simplicial source. Rebuild after changing
filtration values/order; IDs belong to the source. Face/coface closure, data
derivatives and tie conventions remain caller responsibilities. The maximum
absolute-displacement heuristic has no general descent/convergence guarantee.

The guide/example explain the new API; the paper note includes all three
main-text pseudocodes. The specification records lazy/Partial U invariants.
Rustdoc, errors, controls, independent tests, examples, changelog and CI paths
are updated. Dense is excluded from the Cargo package. No runtime dependency,
unsafe code, MSRV change or public backend switch is introduced.

Review the opt-in workspace separately from the default-routing decision. A
follow-up should isolate high-dimensional decomposition/support costs, use
clean same-build controls and retain unfavorable paired cases. Eight synthetic
fixtures cannot establish an adaptive size threshold. If a future route regresses,
restore private dispatch while preserving public APIs and rerun field/coverage
and independent-oracle checks. This PR applies that rollback already.

## Verification and retained evidence

Final local debug/release/MSRV suites each passed 216 tests, with four ignored
diagnostics. All-target Clippy, ten examples, example tests, eleven strict
Markdown doctest pages, strict rustdoc, source/docs and package verification
passed. The packaged critical-set example and its hand-derived simplification
test passed. The 98 Linux tool tests and bundled native 132-cell smoke passed on
unchanged harness bytes. Exact-head hosted checks are separate; local checks do
not establish remote merge or release.

All final samples, warmups, plans, fixtures, hashes, generated sources, builds,
worker/header hashes, GUDHI results, frozen control and postflight script remain:

`target/benchmarks/commit-bf0069ba4df7/critical-sets/run-001/`

`summary.json` contains all 140 cells; `analysis.json` contains paired ratios and
memory spreads; `oracle/gudhi.json` retains independent fixture results.

`manifest.json` SHA-256 is `901bb093138abba06b29d2a95c0f11954c7d1a8c1bd967ee01f78d3e4491e383`.
Pinned Oineus archive SHA-256 is
`614cc6972ad8b9c3a8edf78bd38560aa08ee37d3c4b4e253d6a07be1bb3ce23a`;
`upstream-verification.json` proves all 47 consumed upstream files matched it.

The formal run is local-only with no public artifact URL or durable retention.
Raw results, binaries and archives are not committed. CI's native reference
artifact retains a separate quick critical smoke for 14 days, identified by
run/attempt/head SHA; it is not this formal ranking run. Historical `ceea0f4`,
failed setup attempts, dirty diagnostics and rejected clean trials are preserved
in their original ignored locations.

## Replay

Use a clean Linux/WSL checkout of the measured candidate and the recorded
compiler versions, Python, g++ and Boost. Prepare pinned Oineus `include/` and
`extern/` together per the protocol. Bundled dense source removes dependence on
the unpublished old Git blob. Different hardware/compiler/boundaries need new
measurements. With Oineus at `target/oineus-pinned`:

```sh
repo="$(pwd)"
git worktree add --detach "$repo/target/critical-baseline"   0202534e3376360f6c5d16cbf2673e2ec4998aa1
python3 "$repo/target/critical-baseline/tools/benchmark_critical_sets.py"   --build-only --rustc "$(command -v rustc)"   --oineus-source "$repo/target/oineus-pinned" --boost-include /usr/include   --dense-source "$repo/benches/optimization/dense_core.rs"   --output "$repo/target/replay-baseline-run-001"
python3 tools/benchmark_critical_sets.py   --rustc "$(command -v rustc)"   --oineus-source target/oineus-pinned --boost-include /usr/include   --dense-source benches/optimization/dense_core.rs   --baseline-build target/replay-baseline-run-001   --samples 12 --cpu 0 --timeout 60   --output target/benchmarks/commit-bf0069ba4df7/critical-sets/run-002
```

Every output directory must be new. On WSL use explicit compiler/runtime paths
if Rust is not on PATH. Run GUDHI 3.12.0 separately on saved full fixtures for
the external diagram check. The controller checks every sample against its
independent pairing oracle and dense targets. Pinned native CI GUDHI/Ripser
suites add supported VR correctness, not rankings for these supplied filtrations.
