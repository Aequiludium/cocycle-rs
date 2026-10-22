# Installation

[Documentation](../index.md) / Getting started

## Released crate

Cocycle requires **Rust 1.91 or later** and uses edition 2024. Install the Rust
toolchain from [rustup](https://rustup.rs/), then create an application:

```sh
cargo new topology-demo
cd topology-demo
cargo add cocycle@0.1.1
```

This adds the `cocycle` crate, distributed under the MIT license, with no runtime
dependencies. See the [quickstart](quickstart.md) for `src/main.rs` and the first
calculation. The [published API](https://docs.rs/cocycle/0.1.1/cocycle/) describes
this release.

## Current source

The site also describes main APIs marked **Unreleased**, including borrowed
stages, critical sets and feasible F2 homology operators. They require a source
checkout or a Git dependency pinned to a revision containing the feature. Read
the [version status](../README.md#version-and-capability-status) and
[changelog](../../CHANGELOG.md#unreleased) before choosing a revision.

```sh
git clone https://github.com/Aequiludium/cocycle-rs.git
cd cocycle-rs
cargo run --locked --example square
cargo doc --locked --no-deps --open
```

These commands build the local checkout. Record `git rev-parse HEAD` when
reporting results; a moving main branch is not a reproducible version pin.

## Development tools

Rustfmt and Clippy check Rust changes. Python runs repository checks; Sphinx and
its documentation dependencies build this website. Native GUDHI, Ripser and
Topp comparisons have separate optional setup requirements and do not become
Rust runtime dependencies. Follow the [development environment](../../CONTRIBUTING.md#development-environment)
and [verification commands](../../CONTRIBUTING.md#verification) for the affected
kind of change.
