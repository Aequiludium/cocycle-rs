# Contributing to Cocycle

Cocycle is a small, general-purpose Rust library. Contributions should improve a
concrete mathematical capability, correctness, usability, or measured performance.
Application workflows and language bindings are outside this crate.

## Development environment

Use Rust 1.91 or later, rustfmt, and Clippy. The core has no external runtime or
test dependencies. Python's standard library runs documentation checks and tool
tests. Native comparisons need a C++17 compiler, Boost headers, and pinned GUDHI
and Ripser sources; see the [native setup](benches/native/README.md). Diagram-distance
comparisons additionally require a C++20 Topp adapter, CGAL headers and pinned
Python GUDHI/NumPy/POT packages; see the [distance setup](benches/distances/README.md#prepare-the-reference-environment).
Other Python TDA wrappers remain [optional checks](tools/legacy-benchmarks.md).
Algorithm authors can start with the [contribution paths](docs/development/algorithm-contributions.md)
and the focused commands below. For wider kernel work, start with
`cargo test --locked` and the [architecture](docs/development/architecture.md).

Use the [issue tracker](https://github.com/Aequiludium/cocycle-rs/issues) for
reproducible bugs and substantial API or algorithm proposals. Small, focused fixes
can be reviewed directly through pull requests. Do not add placeholder implementations for future work.

## Code and API rules

Follow the [code conventions](docs/development/conventions.md) for file boundaries,
naming, reuse, public contracts, errors, numeric behavior, language, and style.
Keep changes focused. New dependencies need a concrete benefit and license/MSRV
review; Python tools and C++ comparison workers are outside the Rust runtime.

Public API includes documented numeric, ordering, and error semantics, not only
Rust signatures. Preserve compatibility within a 0.x minor line and describe
incompatible changes and migration in the changelog. Follow the
[release procedure](#release-procedure) for version selection and MSRV changes.

## Focused algorithm checks

Choose the command for your contribution from a source checkout:

| Path | Tutorial |
| --- | --- |
| Diagram statistics, curves and features | [Diagram analysis](docs/development/diagram-analysis.md) |
| Explicit simplicial construction | [Complex construction](docs/development/complex-construction.md) |
| Diagram matching distances | [Distance contribution path](docs/development/diagram-analysis.md#contribute-diagram-distances) |
| Persistence and reduction | [Persistence algorithm walkthrough](docs/development/persistence-reduction.md) |

```sh
python3 tools/check_algorithm.py diagram-analysis
python3 tools/check_algorithm.py complex-construction
python3 tools/check_algorithm.py diagram-distances
python3 tools/check_algorithm.py persistence-reduction
```

Each command checks source/documentation hygiene, domain formatting, Clippy for
the library and selected targets, domain tests, example tests, the example and
tutorial doctests. Diagram analysis does not run persistent homology; construction
checks also verify persistence of hand-derived complexes. Distance checks run
public contracts, private matching oracles and the distance example; context tests
also compute small supplied and Rips filtrations. Reduction checks run direct
algorithm and independent-oracle tests, plus filtered-source, field, representative
and resource integration tests. They use the existing flag example. The lower-star example's
colocated tests are explicitly run with `cargo test --example complex_construction`;
ordinary `cargo test` alone does not execute them. No native C++ setup is needed.
Python invokes the local Rust toolchain; generated files remain in Cargo's target
directory. New test files or tutorial pages in this path must also be added to
the focused check.

The author supplies mathematical assumptions, implementation and independent
tests. Maintainers help with public exports, errors, allocation/execution policies
and integration. A focused pass is local feedback; maintainers and CI complete
the applicable full verification before merge. Contributions touching other
domains use the checks below as well.

## Verification

Run these from the repository root for Rust changes. CI distributes them across
quality, test, and MSRV jobs:

```sh
cargo fmt --all -- --check
rustfmt --edition 2024 --check tools/diagram_dump.rs tools/benchmark_driver.rs benches/native/cocycle.rs tools/reference/rips_cocycle.rs tools/reference/sparse_cocycle.rs benches/pipeline/cocycle.rs benches/distances/cocycle.rs
rustfmt --edition 2024 --config skip_children=true --check benches/optimization/cocycle.rs benches/optimization/dense_core.rs
cargo clippy --locked --all-targets --all-features -- -D warnings
cargo test --locked --all-features
cargo test --locked --release --all-features
cargo run --locked --example square
cargo run --locked --example rips_graph
cargo run --locked --example flag_persistence
cargo run --locked --example rips_sphere
cargo run --locked --example rips_representatives
cargo run --locked --example sparse_rips
cargo run --locked --example diagram_analysis
cargo run --locked --example complex_construction
cargo test --locked --example complex_construction
cargo test --locked --release --example complex_construction
cargo run --locked --example diagram_distances
cargo run --locked --example critical_sets
cargo test --locked --example critical_sets
cargo test --locked --release --example critical_sets
RUSTDOCFLAGS="-D warnings" cargo doc --locked --no-deps
cargo +1.91.0 test --locked --all-features
cargo +1.91.0 check --locked --all-targets --all-features
cargo +1.91.0 test --locked --example complex_construction
cargo +1.91.0 test --locked --example critical_sets
python3 tools/check_source.py
python3 tools/check_artifacts.py
python3 tools/check_docs.py
python3 -m unittest discover -s tools -p 'test_*.py'
cargo build --locked
rustdoc --edition 2024 --test README.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
rustdoc --edition 2024 --test docs/guides/rips.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
rustdoc --edition 2024 --test docs/guides/rips-construction.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
rustdoc --edition 2024 --test docs/guides/rips-representatives.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
rustdoc --edition 2024 --test docs/guides/sparse-rips.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
rustdoc --edition 2024 --test docs/guides/filtered-complexes.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
rustdoc --edition 2024 --test docs/guides/critical-sets.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
rustdoc --edition 2024 --test docs/guides/interface-contracts.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
rustdoc --edition 2024 --test docs/development/diagram-analysis.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
rustdoc --edition 2024 --test docs/development/complex-construction.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
rustdoc --edition 2024 --test docs/development/persistence-reduction.md --extern cocycle=target/debug/libcocycle.rlib -L dependency=target/debug/deps
```

Select additional checks by the changed contract:

| Change | Required evidence |
| --- | --- |
| Documentation only | Source and Markdown checks; run affected Rust examples/doctests |
| Rust implementation or public API | Commands above; update contract tests and relevant rustdoc |
| Mathematical algorithm | Independent expectation/property and relevant native comparisons, in addition to Rust checks |
| Diagram-distance algorithm | Independent tiny matching oracle, pinned Topp and GUDHI comparisons under the [distance protocol](benches/distances/README.md) |
| Python checks or controllers | Source and Markdown checks, all `test_*.py`; exercise the changed command on a small case |
| Native adapters, builder, or benchmark protocol | Tool tests, standalone Rust formatting when affected, and the native smoke command below |
| Performance | Comparable before/after measurements under the applicable native suite and reporting rules; keep unfavorable results |
| File layout or packaging | Relevant checks above, package file-list review, package build and packaged example |

Source and documentation checks need only Python's standard library.
The artifact check also reads the Git index; run it after staging. It rejects
generated outputs regardless of age, including forced additions. CI checks the
same index policy. For
mathematical or algorithm changes, update the [specification](docs/reference/mathematics.md), add an
independent expected result or property, and run the relevant external comparisons
in [tools/README.md](tools/README.md). Performance changes need the same fixtures,
precision, and timing boundaries before and after; keep unfavorable results.
Ordinary tests must not assert machine-dependent timing thresholds.

Benchmark changes must follow the [reporting rules](benches/reporting.md) and the
affected execution contract: [H0/H1 native](benches/protocol.md) or
[Rips pipeline](benches/pipeline/README.md). Run the affected native smoke suite
when changing workers, controllers or measurement semantics. GUDHI and Ripser comparisons use C++
executables; Python TDA wrappers remain optional interface checks, outside native rankings. Tool
changes also need `python3 -m unittest discover -s tools -p 'test_*.py'`.

After the native setup, use a new output directory for each run:

```sh
python3 tools/benchmark_native.py --quick --samples 1 --output target/native-smoke
```

Supply `--boost-include` when Boost headers are outside system include paths.
For exact graph/matrix changes also run
`python3 tools/compare_rips.py --output target/rips-reference` with a fresh output
directory and the same optional Boost setting. For sparse approximation, also run
`python3 tools/compare_sparse_rips.py --output target/sparse-rips-reference`
with a fresh output directory. For workflow timing/resource changes, also run
`python3 tools/benchmark_rips_pipeline.py --quick --samples 1 --output target/rips-pipeline-smoke`.
Treat unsupported reference inputs as documented exclusions, not successful
cross-library comparisons. See the native guide for platform requirements.

For diagram-distance implementation changes, run the full supported suite after
the distance setup, with a fresh output directory:

```sh
python3 tools/compare_distances.py --gudhi-python target/distance-oracle-venv/bin/python --output target/distance-correctness-full
```

Record the tested commit and summary in the PR. CI's `--quick` distance check
covers a smaller suite. Worker/instrumentation changes also need the small
[distance resource smoke](benches/distances/README.md#check-correctness-and-the-harness),
including `--profile-rust` when changing generated timing hooks.

Documentation-only changes need source/documentation checks and any affected examples;
they do not require rerunning large performance experiments. See
[testing](docs/development/testing.md) for coverage and CI responsibilities.

## Pull requests

Describe the problem, resulting behavior, mathematical basis if relevant, and
checks actually run. Identify incomplete checks and limitations. Update the
unreleased changelog for user-visible changes. Avoid unrelated formatting or
speculative abstractions. The pull request template is a guide, not a requirement
to add irrelevant sections.

Review focuses on correctness, public contracts, independent evidence, and whether
the change fits the library's scope. Discuss technical choices respectfully and
make feedback specific and actionable.

The repository policy is recorded in the [CI ruleset](.github/rulesets/main.json)
and [review ruleset](.github/rulesets/review.json). These files must be applied
through GitHub's repository settings or API; committing them does not activate
the rules. CI checks, an up-to-date branch, and protection against force pushes
and deletion apply to everyone. Pull requests normally require one approval
from another [code owner](.github/CODEOWNERS); new commits dismiss stale approvals,
and review conversations must be resolved. Repository administrators may
explicitly bypass the review rules on a pull request, but cannot bypass CI.

GitHub Actions are pinned to full commit SHAs. Dependabot checks for updates
weekly; maintainers review these updates through the same pull request process.
Keep the required check names in the CI ruleset aligned when renaming CI jobs.

## Documentation ownership

The [documentation index](docs/README.md#organization-and-maintenance) defines
where usage, reference, development, design, and research pages belong. It also
links each authoritative document. API behavior belongs in rustdoc next to the
item; coding rules belong in the code conventions; verification commands and
contribution procedures belong in this file. Link to the owner rather than
duplicating its source tree, commands, or detailed rules.

When adding or moving a document:

1. Update the index, relative links, heading fragments, and literal paths in
   commands or tooling. Clearly label proposals and dated source observations.
2. Run `python3 tools/check_docs.py`. It recursively checks `docs/` and selected
   repository Markdown for inline/reference local targets, headings, and CJK
   text. It does not fetch external URLs or establish English prose quality.
   Run affected Rust examples using the commands above.
3. Keep CI's guide doctest path aligned with the guide. If changing the checker,
   run `python3 -m unittest discover -s tools -p 'test_check_docs.py'`.
4. Do not stage logs, benchmark outputs, fixture collections or archives. Keep
   local runs under `target/` and publish selected evidence through CI artifacts
   or external storage; concise reports link that evidence.
5. If changing directories or package inclusion, inspect
   `cargo package --locked --allow-dirty --list` and ensure nested docs remain in
   the package. Raw benchmark artifacts remain outside the crate payload.

Performance reports use stable topic-based paths; update the existing report on
reruns and let Git retain its history. Each report binds to the measured commit
and associated PR, with a full SHA, source fingerprint and protocol. Dates are
execution metadata. The report commit and measured commit are distinct; uncommitted runs remain drafts. Do not
rewrite its raw data after a refactor or treat old timing as a new measurement.
The [reporting rules](benches/reporting.md) own evidence classification,
comparability, sampling and artifact retention; use the
[report template](benches/report-template.md) for new experiments.

## Release procedure

The crate is published as [cocycle](https://crates.io/crates/cocycle). Releases are
explicit maintainer actions after a reviewed release PR and successful CI.
GitHub CI verifies changes; it does not publish packages or select versions.

### Version and notes

`Cargo.toml` is the version source. While the crate is `0.x` with `x >= 1`, use a
patch release for compatible additions, fixes and performance improvements, and
a new minor version for incompatible changes. Compatibility includes documented
mathematical semantics, defaults and errors as well as Rust signatures; see
[Cargo's compatibility guidance](https://doc.rust-lang.org/cargo/reference/semver.html).
Treat an increased MSRV as an incompatible change for this project, update the
manifest, CI and user documentation together, and state it in the release notes.
Use a suffix such as `-rc.1` only when a pre-release is needed.

[CHANGELOG.md](CHANGELOG.md) owns release notes. Contributions add user-visible
changes to `Unreleased`; a release PR moves them into a version section and leaves
`Unreleased` ready for the next changes. Describe behavior, limitations and any
migration steps. GitHub Release reuses that section rather than a separate history.

### Prepare, publish and verify

1. Open a release PR updating `Cargo.toml`, the corresponding `Cargo.lock` entry,
   the changelog and README version references. Running `cargo check` after the
   version edit refreshes the lockfile. Review the version and notes, then merge
   through the normal PR process.
2. Use a clean checkout of that merged commit. Confirm `git status --short` is
   empty and all required CI checks succeeded for its exact SHA, including the
   package job. That job builds the archive, runs its `square`, `diagram_analysis`,
   `diagram_distances` and `complex_construction` examples, and tests the packaged
   construction example. Do not substitute checks from an earlier PR commit.
3. Inspect `cargo package --locked --list` for unintended files and run
   `cargo publish --locked --registry crates-io --dry-run` from that checkout.
   Authenticate as a crates.io owner using `cargo login` if needed. Then run
   `cargo publish --locked --registry crates-io`. Keep credentials outside Git.
4. In a temporary project outside this checkout, depend on `cocycle` with the
   exact published version (`=X.Y.Z`), copy the README Rust example and run it.
   Use the registry dependency, without a path dependency or local patch. Verify
   the version's API page at `https://docs.rs/cocycle/X.Y.Z/cocycle/` once the
   documentation build completes.
5. Create an annotated `vX.Y.Z` tag on the verified source commit and push that
   tag. Create the GitHub Release from this existing tag, copying the corresponding
   changelog section. Mark pre-releases as such in GitHub. Never move a published
   version's tag to different source.

If publication times out, check crates.io before retrying: the upload may already
have succeeded. If it did, finish verification and the GitHub Release for that
same source. A published version cannot be replaced; corrections need a new
version. See [Cargo's publishing command](https://doc.rust-lang.org/cargo/commands/cargo-publish.html).

`cargo package --allow-dirty` is only for local review of uncommitted work, not
publication. API documentation is generated from rustdoc on docs.rs; guides stay
in `docs/`. Build output and release scratch files stay in `target/` or outside
the repository.

Contributions are made under the project's [MIT License](LICENSE). Attribute
external code and check its license before incorporating it; citing an algorithm
paper does not license copying an implementation.
