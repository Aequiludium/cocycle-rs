"""Reproduce the finite-diagram preparation counterexample in isolated snapshots.

Linux, Python standard library, Git and Rust 1.91+ only. See the maintained
report in benches/reports/distance-preparation.md for interpretation limits.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import os
from pathlib import Path
import platform
import random
import statistics
import struct
import subprocess
import sys
import tarfile

ROOT = Path(__file__).resolve().parents[1]
BASELINE = "75330fcf7e5b9ca764641c68b987a71bbfad8498"
VARIANTS = ("direct", "copy_keep", "copy_drop", "unused_copy_keep",
            "reserve_keep", "reserve_drop", "reserve_after_prepare")
HARNESS = ("tools/reproduce_distance_preparation.py",
           "benches/distances/preparation.rs", "benches/distances/preparation_probe.rs")
COUNTERS = ("preparation_ns", "solver_ns", "copy_cleanup_ns", "calls",
            "graph_ns", "matching_reconstruction_ns", "grouping_ns",
            "capacity_growths", "relocations", "relocated_live_bytes")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def git(*args):
    return subprocess.check_output(["git", *args], cwd=ROOT)


def fixture(size):
    rng = random.Random(20260928 + size)
    first = [(birth := rng.uniform(-1, 1), birth + rng.uniform(1, 3))
             for _ in range(size)]
    second = [(birth := rng.uniform(-1, 1), birth + rng.uniform(1, 3))
              for _ in range(max(1, size // 8))]
    return (b"COCDST1\0" + struct.pack("<QQ", len(first), len(second))
            + b"".join(struct.pack("<dd", *point) for point in first + second))


def replace_once(text, old, new):
    if text.count(old) != 1:
        raise ValueError(f"Expected exactly one source marker: {old!r}")
    return text.replace(old, new, 1)


def section(text, begin, end):
    if text.count(begin) != 1:
        raise ValueError(f"Expected exactly one section: {begin!r}")
    start = text.index(begin)
    stop = text.index(end, start)
    return text[start:stop]


def extract_snapshot(archive, destination):
    """Extract only regular files/directories, without traversal or overwrites."""
    destination.mkdir()
    with tarfile.open(fileobj=io.BytesIO(archive)) as stream:
        for member in stream.getmembers():
            path = (destination / member.name).resolve()
            if not path.is_relative_to(destination.resolve()):
                raise ValueError(f"Archive path escapes snapshot: {member.name}")
            if member.isdir():
                path.mkdir(parents=True, exist_ok=True)
            elif member.isfile():
                path.parent.mkdir(parents=True, exist_ok=True)
                with path.open("xb") as output:
                    output.write(stream.extractfile(member).read())
            else:
                raise ValueError(f"Unsupported archive entry: {member.name}")


def instrument(source, template, growth):
    """Checked edits apply only to the generated candidate source snapshot."""
    directory = source / "src/diagram_distances"
    numeric_path = directory / "wasserstein/numeric.rs"
    numeric = numeric_path.read_text(encoding="utf-8")
    copied = section(numeric, "pub(super) fn prepare_dimensions", "\nfn collect_points")
    copied = replace_once(copied, "prepare_dimensions", "prepare_coordinates")
    if copied.count("&DiagramDimension<'_>") != 2:
        raise ValueError("Unexpected preparation arguments")
    copied = copied.replace("&DiagramDimension<'_>", "&[[f64; 2]]")
    for side in ("first", "second"):
        copied = replace_once(copied, f"{side}.iter().map(|i| (i.birth(), i.end()))",
                              f"{side}.iter().map(|p| (p[0], IntervalEnd::Finite(p[1])))")
    numeric_path.write_text(numeric + "\n" + copied, encoding="utf-8")
    path = directory / "wasserstein.rs"
    text = path.read_text(encoding="utf-8")
    original = section(text, "pub(crate) fn from_dimensions", "\n#[cfg(test)]")
    adapter, growth_adapter = template.split(
        "// Growth diagnostics: inserted only for the separately timed growth mode.")
    text = replace_once(text, original, adapter + (growth_adapter if growth else ""))
    text = replace_once(text, "fn solve_prepared<const CONTROLLED: bool>(",
                        "#[inline(never)]\nfn solve_prepared<const CONTROLLED: bool>(")
    text = replace_once(text, "    let original_pairs = product(first.len(), second.len())?;",
                        "    let grouping_started = std::time::Instant::now();\n"
                        "    let original_pairs = product(first.len(), second.len())?;")
    text = replace_once(text, "    let graph = generate(&first, &second, metric, _stats, budget)?;",
                        "    PHASES[6].fetch_add(grouping_started.elapsed().as_nanos() as u64, Ordering::Relaxed);\n"
                        "    let graph_started = std::time::Instant::now();\n"
                        "    let graph = generate(&first, &second, metric, _stats, budget)?;\n"
                        "    PHASES[4].fetch_add(graph_started.elapsed().as_nanos() as u64, Ordering::Relaxed);\n"
                        "    let matching_started = std::time::Instant::now();")
    text = replace_once(text, "    restore_scale(\n"
                        "        from_matching(&first, &second, matching, metric, budget)?,\n"
                        "        scale,\n    )\n}",
                        "    let answer = restore_scale(\n"
                        "        from_matching(&first, &second, matching, metric, budget)?, scale,\n    );\n"
                        "    PHASES[5].fetch_add(matching_started.elapsed().as_nanos() as u64, Ordering::Relaxed);\n"
                        "    answer\n}")
    path.write_text(text, encoding="utf-8")
    facade = directory / "mod.rs"
    with facade.open("a", encoding="utf-8") as stream:
        stream.write("\n/// Isolated preparation experiment selector.\n"
                     "pub fn select_preparation_experiment(variant: usize, verify: bool) {\n"
                     "    wasserstein::select_experiment(variant, verify);\n}\n"
                     "/// Inclusive diagnostic counters; see the reproducer protocol.\n"
                     "pub fn preparation_experiment_phases() -> [u64; 10] {\n"
                     "    wasserstein::experiment_phases()\n}\n")
    if growth:
        path = directory / "wasserstein/graph.rs"
        text = path.read_text(encoding="utf-8")
        original = section(text, "pub(super) fn generate<const CONTROLLED: bool>(",
                           "\npub(super) struct Component")
        if original.count("                    push(") != 2:
            raise ValueError("Unexpected candidate graph push sites")
        replacement = original.replace("                    push(",
                                       "                    super::trace_edge_push(")
        path.write_text(replace_once(text, original, replacement), encoding="utf-8")


def inventory(directory):
    return {p.relative_to(directory).as_posix(): sha(p.read_bytes())
            for p in sorted(directory.rglob("*")) if p.is_file()}


def build(output, name, revision, experimental, growth):
    directory = output / name
    directory.mkdir()
    source = directory / "source"
    extract_snapshot(git("archive", revision, "Cargo.toml", "Cargo.lock", "src",
                         "benches/rips.rs"), source)
    original = inventory(source)
    if experimental:
        instrument(source, (ROOT / HARNESS[2]).read_text(encoding="utf-8"), growth)
    worker_source = directory / "worker.rs"
    worker_source.write_bytes((ROOT / HARNESS[1]).read_bytes())
    build_dir = directory / "build"
    worker = directory / "worker"
    commands = [["cargo", "build", "--release", "--lib", "--locked", "--offline",
                 "--manifest-path", str(source / "Cargo.toml")],
                ["rustc", "--edition", "2024", "-O", str(worker_source),
                 "--extern", f"cocycle={build_dir}/release/libcocycle.rlib",
                 "-L", f"dependency={build_dir}/release/deps", "-o", str(worker)]]
    if experimental:
        # Only the worker gets this cfg. The library's other stats stay disabled.
        commands[1].extend(["--cfg", "cocycle_distance_bench"])
    environment = dict(os.environ, CARGO_TARGET_DIR=str(build_dir))
    identity = {"revision": revision, "experimental": experimental, "growth": growth,
                "commands": commands, "original_sources": original,
                "generated_sources": inventory(source), "worker_source_sha256": sha(worker_source.read_bytes())}
    write_json(directory / "identity.json", identity)
    with (directory / "build.log").open("w", encoding="utf-8") as log:
        for command in commands:
            log.write(json.dumps(command) + "\n")
            log.flush()
            subprocess.run(command, cwd=ROOT, env=environment, stdout=log,
                           stderr=subprocess.STDOUT, check=True, timeout=300)
    identity["binary_sha256"] = sha(worker.read_bytes())
    write_json(directory / "identity.json", identity)
    return worker


def distribution(values):
    return {"median": statistics.median(values), "min": min(values), "max": max(values),
            "quartiles_exclusive": statistics.quantiles(values, n=4) if len(values) > 1 else None}


def measure(output, mode, workers, args, expected):
    variants = list(workers)
    order = list(variants)
    rng = random.Random(140931)
    records = {variant: [] for variant in variants}
    for round_index in range(args.samples):
        if round_index % len(order) == 0:
            rng.shuffle(order)
        offset = round_index % len(order)
        for variant in order[offset:] + order[:offset]:
            worker, selector = workers[variant]
            command = [str(worker), str(output / "fixture.bin"), args.metric,
                       str(args.iterations), str(selector)]
            record = {"mode": mode, "variant": variant, "round": round_index, "command": command}
            try:
                result = subprocess.run(command, capture_output=True, text=True, timeout=120)
                record.update(exit_code=result.returncode, stdout=result.stdout, stderr=result.stderr)
                if result.returncode == 0:
                    record.update(json.loads(result.stdout))
            except subprocess.TimeoutExpired as error:
                record.update(error=str(error), stdout=(error.stdout or b"").decode(errors="replace"),
                              stderr=(error.stderr or b"").decode(errors="replace"))
            except ValueError as error:
                record["error"] = str(error)
            with (output / "raw.jsonl").open("a", encoding="utf-8") as stream:
                stream.write(json.dumps(record) + "\n")
            if record.get("exit_code") != 0 or "error" in record:
                raise RuntimeError(f"Worker failed; see {output / 'raw.jsonl'}")
            if expected is None:
                expected = record["value"]
            if record["value"] != expected:
                raise RuntimeError("Distance outputs differ; raw records preserved")
            if mode != "public" and record["counters"][3] != args.iterations:
                raise RuntimeError("Unexpected instrumented call count")
            records[variant].append(record)
    summary = {"mode": mode, "value": expected, "outputs_identical": True,
               "samples_per_variant": args.samples, "variants": {}}
    for variant, rows in records.items():
        summary["variants"][variant] = {
            "ns_per_call": distribution([r["elapsed_ns"] / args.iterations for r in rows]),
            "peak_rss_kib": distribution([r["peak_rss_kib"] for r in rows]),
            "counters_per_call": {name: statistics.median(r["counters"][i] / args.iterations for r in rows)
                                  for i, name in enumerate(COUNTERS)}}
    write_json(output / f"{mode}-summary.json", summary)
    print(json.dumps(summary), flush=True)
    return expected


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--mode", choices=("public", "causal", "growth", "all"), default="all")
    parser.add_argument("--baseline", default=BASELINE)
    parser.add_argument("--candidate", default="HEAD")
    parser.add_argument("--metric", choices=("w1", "w2"), default="w1")
    parser.add_argument("--size", type=int, default=512)
    parser.add_argument("--samples", type=int, default=14)
    parser.add_argument("--iterations", type=int, default=300)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    if sys.platform != "linux":
        parser.error("Run on Linux (including Ubuntu WSL); /proc and CPU affinity are required")
    if min(args.samples, args.iterations, args.size) < 1:
        parser.error("samples, iterations and size must be positive")
    # Environment overrides can silently enable unrelated instrumentation or codegen.
    if any(os.environ.get(name) for name in ("RUSTFLAGS", "CARGO_ENCODED_RUSTFLAGS")):
        parser.error("Unset RUSTFLAGS and CARGO_ENCODED_RUSTFLAGS for comparable builds")
    output = args.output.resolve()
    if not output.is_relative_to((ROOT / "target").resolve()):
        parser.error("output must be a fresh directory inside this repository's ignored target/")
    harness_revision = git("rev-parse", "HEAD").decode().strip()
    dirty = [name for name in HARNESS
             if subprocess.run(["git", "diff", "--quiet", "HEAD", "--", name], cwd=ROOT).returncode
             or subprocess.run(["git", "ls-files", "--error-unmatch", name], cwd=ROOT,
                               stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL).returncode]
    if args.samples >= 10 and dirty:
        parser.error("Commit the harness before a formal run; use samples < 10 for development smoke")
    baseline = git("rev-parse", "--verify", f"{args.baseline}^{{commit}}").decode().strip()
    candidate = git("rev-parse", "--verify", f"{args.candidate}^{{commit}}").decode().strip()
    output.mkdir(parents=True, exist_ok=False)
    payload = fixture(args.size)
    (output / "fixture.bin").write_bytes(payload)
    cpu = min(os.sched_getaffinity(0))
    metadata = {"class": "smoke only" if args.samples < 10 else "local comparative experiment",
                "baseline": baseline, "candidate": candidate, "harness_revision": harness_revision,
                "harness_dirty_paths": dirty, "harness_sha256": {p: sha((ROOT / p).read_bytes()) for p in HARNESS},
                "fixture_sha256": sha(payload), "seed": 20260928 + args.size,
                "size": [args.size, max(1, args.size // 8)], "metric": args.metric,
                "iterations": args.iterations, "samples": args.samples, "warmup_calls": 3,
                "counter_order": COUNTERS, "cpu_affinity": cpu, "platform": platform.platform(),
                "python": sys.version, "libc": platform.libc_ver(),
                "rustc": subprocess.check_output(["rustc", "-Vv"], text=True),
                "cargo": subprocess.check_output(["cargo", "-V"], text=True),
                "cpuinfo": Path("/proc/cpuinfo").read_text().split("\n\n")[0],
                "frequency_and_host_load": "uncontrolled",
                "sampling": "Fresh serial processes; rotate seeded shuffled order each block; no outlier filtering",
                "timing": "Warmed public calls including preparation, solving and destruction; fixture and diagram construction excluded",
                "started": datetime.now(timezone.utc).isoformat()}
    write_json(output / "environment.json", metadata)
    try:
        modes = ("public", "causal", "growth") if args.mode == "all" else (args.mode,)
        binaries = {}
        for mode in modes:
            if mode == "public":
                binaries[mode] = {name: (build(output, name, revision, False, False), 0)
                                  for name, revision in (("baseline", baseline), ("candidate", candidate))}
            else:
                worker = build(output, mode, candidate, True, mode == "growth")
                binaries[mode] = {name: (worker, i) for i, name in enumerate(VARIANTS)}
        os.sched_setaffinity(0, {cpu})
        expected = 442.10173773952255 if args.size == 512 and args.metric == "w1" else None
        for mode in modes:
            expected = measure(output, mode, binaries[mode], args, expected)
        metadata["status"] = "passed"
    except Exception as error:
        metadata.update(status="failed", error=str(error))
        raise
    finally:
        metadata["finished"] = datetime.now(timezone.utc).isoformat()
        write_json(output / "environment.json", metadata)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
