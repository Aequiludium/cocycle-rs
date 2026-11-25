# Directed path H0/H1 comparison

[Benchmarks](../README.md) / [Protocol](../path/README.md)

## Question and observed scope

Compare the current generic compatible-Omega route with the direct one-boundary presentation on identical inputs and the same ordinary PersistenceResult output. Both include validation, the same cooperative work policy, canonical interval/context assembly, payload export and cleanup. The generic baseline includes its current endpoint grouping and omits unused static representatives.

All 64 declared cells completed and matched their complete interval multisets in every sample. Direct met the predeclared advantage rule in all 64 cells: generic/direct median ratios range from 3.076 to 14.233. These are scoped measurements of regular, non-reduced F2 H0/H1 diagrams, without representatives or certificate output.

Observed advantageous vertex counts were 8/32/128/512 for cycles and sparse random digraphs, 8/16/32/64 for dense random digraphs, and 9/18/36/72 for three-layer multisquare graphs. Both zero-birth static and three-event filtered inputs, with two seeds, passed. The smallest tested sizes already favored direct; this experiment does not locate a universal crossover. Density and two-path count matter alongside vertex count.

## Measured identities and environment

- PR: [integration #168](https://github.com/Aequiludium/cocycle-rs/pull/168), with the direct route and harness contributed through [#166](https://github.com/Aequiludium/cocycle-rs/pull/166).
- Both measured candidate (direct) and baseline (generic): `84556cf532a77435af0ae1c6504c5567554fbefb`, an actual tested integration merge commit. The report commit is a later documentation update.
- Protocol/controller: `cocycle-path-diagram-v1`, `tools/benchmark_path_topology.py`, `benches/path/cocycle.rs` and [the protocol](../path/README.md), all committed in that measured tree.
- SHA-256 of canonical sorted source-hash JSON: `e46026d407bc5f14857cff8590433e415a6e219e8607fdac520172f14b6796d0`. Worker SHA-256: `0b4d6e2807ef654fccf7d793df7afcdcb1d5ab4045ba3132c57ffc9d39a90f08`.
- Host: Windows-2025Server-10.0.26100-SP0; AMD64 Family 25 Model 1 Stepping 1, AuthenticAMD; 4 logical CPUs; Python 3.12.10; Rust 1.98.1 / x86_64-pc-windows-msvc / LLVM 22.1.8.
- Cargo release library and native worker `-C opt-level=3`. Workers ran serially after compilation. Affinity, host load and CPU frequency were uncontrolled; no RSS or OS memory-cap claim is made.
- Formal job interval: 2026-11-25 14:32:28-14:34:17 UTC, including setup/build/artifact upload. The native timer excludes fixture parsing, field validation, process startup and JSON formatting. It includes the shared work-budget checks and counters.

## Sampling and correctness

There are 64 cells, two routes and ten fresh measured processes per route/cell (1280 measured samples), plus one discarded warmup per route/cell (128 warmups). All 1408 outputs and generated input hashes are retained. Round order was position-balanced, cell order seeded, and no outlier was discarded. There were zero failures, timeouts or exclusions. Each successful output was checked for protocol/field/model/domain/coverage and exact interval multiplicity agreement. Independent exact small-graph certificates and the real-path oracle are correctness evidence; they are not additional timing references.

Both routes used max_paths=250000, max_entries=2000000, max_work=250000000 and a 15-second process timeout. Construction-entry caps and cooperative work are distinct from aggregate memory quotas. Default output is only diagram/context for both routes; internal working bases or relations are released before the terminal returns.

The advantage rule requires generic/direct median >=1.1 and the lower paired 95% bootstrap percentile bound >1. The bootstrap uses 2000 paired round-index resamples and the order-statistic quantiles declared in the protocol. Ratios describe individual cells, with no pooled speedup. The tables below show the range of medians/ratios across the two seeds; every row has two independently sampled cells. Full per-cell min/max, intervals and all samples remain in the artifact.

## Static zero-birth diagrams

| Family | Vertices | Edges, seeds range | Two-paths, seeds range | Generic median ms, seeds range | Direct median ms, seeds range | Generic/direct, seeds range |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cycle | 8 | 8 | 8 | 0.165-0.200 | 0.053-0.054 | 3.076-3.671 |
| cycle | 32 | 32 | 32 | 0.436-0.446 | 0.094-0.096 | 4.653-4.657 |
| cycle | 128 | 128 | 128 | 1.442-1.444 | 0.223-0.227 | 6.372-6.477 |
| cycle | 512 | 512 | 512 | 5.286-5.298 | 0.665-0.676 | 7.814-7.971 |
| sparse | 8 | 35-40 | 147-197 | 1.250-1.627 | 0.282-0.328 | 4.437-4.957 |
| sparse | 32 | 123-134 | 481-580 | 2.832-3.224 | 0.566-0.659 | 4.894-5.003 |
| sparse | 128 | 484-533 | 1858-2194 | 8.088-9.040 | 1.168-1.334 | 6.778-6.925 |
| sparse | 512 | 2028-2070 | 8059-8299 | 33.695-34.373 | 4.137-4.412 | 7.792-8.145 |
| dense | 8 | 13-18 | 20-34 | 0.258-0.363 | 0.079-0.109 | 3.264-3.342 |
| dense | 16 | 61-64 | 225-250 | 1.512-1.738 | 0.308-0.344 | 4.913-5.057 |
| dense | 32 | 258-266 | 2027-2258 | 14.520-16.421 | 3.518-3.846 | 4.128-4.269 |
| dense | 64 | 1018-1046 | 16044-16962 | 215.086-240.084 | 27.578-31.537 | 7.613-7.799 |
| multisquare | 9 | 18 | 27 | 0.375-0.403 | 0.102-0.103 | 3.688-3.932 |
| multisquare | 18 | 72 | 216 | 1.877-1.949 | 0.357-0.358 | 5.235-5.458 |
| multisquare | 36 | 288 | 1728 | 13.616-13.901 | 2.109-2.196 | 6.330-6.456 |
| multisquare | 72 | 1152 | 13824 | 204.929-256.197 | 15.124-18.000 | 13.550-14.233 |

## Three-event persistence diagrams

| Family | Vertices | Edges, seeds range | Two-paths, seeds range | Generic median ms, seeds range | Direct median ms, seeds range | Generic/direct, seeds range |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cycle | 8 | 8 | 8 | 0.250-0.381 | 0.056-0.077 | 4.476-4.955 |
| cycle | 32 | 32 | 32 | 0.686-0.699 | 0.103-0.112 | 6.223-6.677 |
| cycle | 128 | 128 | 128 | 2.596-2.619 | 0.266-0.272 | 9.561-9.860 |
| cycle | 512 | 512 | 512 | 10.159-10.350 | 0.807-0.814 | 12.486-12.818 |
| sparse | 8 | 35-40 | 147-197 | 1.569-2.168 | 0.301-0.380 | 5.210-5.703 |
| sparse | 32 | 123-134 | 481-580 | 3.767-4.543 | 0.620-0.786 | 5.779-6.072 |
| sparse | 128 | 484-533 | 1858-2194 | 12.691-13.041 | 1.323-1.433 | 9.098-9.596 |
| sparse | 512 | 2028-2070 | 8059-8299 | 49.965-54.113 | 4.476-4.754 | 11.164-11.383 |
| dense | 8 | 13-18 | 20-34 | 0.420-0.559 | 0.083-0.110 | 5.033-5.063 |
| dense | 16 | 61-64 | 225-250 | 2.341-2.369 | 0.372-0.406 | 5.766-6.361 |
| dense | 32 | 258-266 | 2027-2258 | 21.482-25.543 | 6.460-7.270 | 3.325-3.513 |
| dense | 64 | 1018-1046 | 16044-16962 | 335.446-389.823 | 93.878-110.169 | 3.538-3.573 |
| multisquare | 9 | 18 | 27 | 0.556-0.557 | 0.111-0.116 | 4.800-5.017 |
| multisquare | 18 | 72 | 216 | 2.633-2.795 | 0.454-0.461 | 5.709-6.162 |
| multisquare | 36 | 288 | 1728 | 18.067-20.485 | 3.255-4.116 | 4.977-5.550 |
| multisquare | 72 | 1152 | 13824 | 261.969-280.817 | 33.557-33.870 | 7.807-8.291 |

## Reproduction and retained evidence

[Actions run 45634397071](https://github.com/Aequiludium/cocycle-rs/actions/runs/45634397071) completed successfully. [Performance artifact](https://github.com/Aequiludium/cocycle-rs/actions/runs/45634397071/artifacts/15032473564) contains the full 64-cell summary, every raw sample, warmups, fixtures and the native worker. Its server archive identity is `sha256:33007aa35177c3ed45949fdbed620b8bee5c801227aaf4c8ee14aef283630dae`, expiring at 2026-12-09T14:34:11Z. This 14-day Actions artifact is not permanent external storage; local extracted evidence remains under ignored `target/path-performance-hosted/`.

```sh
git checkout 84556cf532a77435af0ae1c6504c5567554fbefb
python tools/benchmark_path_topology.py --output target/benchmarks/commit-84556cf532a7/path-topology/run-001 --samples 10
```

The same [workflow](../../.github/workflows/path-topology.yml) can be dispatched with performance=true. A new machine/run produces new measurements; it does not reproduce the exact nanosecond samples by assertion. Preserve fresh output directories and the emitted identities.

## Independent real verification

The same fixed commit passed 652 self-contained rational/SVD real-operator certificates across 149 graph/model cases and rejected 16 invalid requests. Largest recorded error was 8.600e-16 against the declared 1e-8 acceptance tolerance. This is finite numerical evidence, with exact rational kernels/persistent ranks for the exercised inputs. The [numerical artifact](https://github.com/Aequiludium/cocycle-rs/actions/runs/45634397071/artifacts/15032772599) includes all certificates and versions. The standalone code and pinned dependencies are in this repository; no local research checkout is required.

## Limits

Conclusions apply to the measured regular F2 H0/H1 diagram terminals, graph families, two seeds, scales and Windows host. Higher dimensions, other fields/models, representatives/certificates, real operators and real-network distributions require separate comparisons. The general route remains useful for higher-dimensional and explicit finite path complexes; these timings alone do not define a dispatch threshold. Scalar counts do not certify a memory bound. Timings in CI smoke runs are not substituted for this formal study.
