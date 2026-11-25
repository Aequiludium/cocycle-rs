"""Serial native matched-output path H0/H1 comparison; raw runs stay in target/."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import time

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = "cocycle-path-diagram-v1"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def command(args):
    return subprocess.check_output(args, cwd=ROOT, text=True, encoding="utf-8").strip()


def fingerprint():
    paths = [*ROOT.glob("src/**/*.rs"), ROOT / "Cargo.toml", ROOT / "Cargo.lock",
             ROOT / "benches/path/cocycle.rs", Path(__file__), ROOT / "benches/path/README.md"]
    return {p.relative_to(ROOT).as_posix(): digest(p) for p in sorted(paths)}


def fixtures(quick, seeds):
    sizes = {"cycle": [8, 32, 128, 512], "sparse": [8, 32, 128, 512],
             "dense": [8, 16, 32, 64], "multisquare": [9, 18, 36, 72]}
    for seed in seeds:
        for family, scale in sizes.items():
            for n in scale[:1] if quick else scale:
                rng = random.Random(f"{seed}:{family}:{n}")
                if family == "cycle":
                    edges = [(i, (i + 1) % n) for i in range(n)]
                elif family == "multisquare":
                    m = n // 3
                    edges = [(a, b) for a in range(m) for b in range(m, 2*m)]
                    edges += [(a, b) for a in range(m, 2*m) for b in range(2*m, n)]
                else:
                    probability = min(1.0, 4/(n-1)) if family == "sparse" else 0.25
                    edges = [(a, b) for a in range(n) for b in range(n) if a != b and rng.random() < probability]
                incoming, outgoing = [0]*n, [0]*n
                for a, b in edges:
                    outgoing[a] += 1
                    incoming[b] += 1
                for filtration in ("static", "filtered"):
                    weights = [0 if filtration == "static" else rng.randrange(1, 4) for _ in edges]
                    yield {"id": f"{family}-{n}-{seed}-{filtration}", "family": family, "vertices": n,
                           "seed": seed, "filtration": filtration, "edges": edges, "weights": weights,
                           "two_paths": sum(a*b for a, b in zip(incoming, outgoing))}


def build(output):
    subprocess.run(["cargo", "build", "--offline", "--locked", "--release"], cwd=ROOT, check=True)
    worker = output / ("path-worker.exe" if os.name == "nt" else "path-worker")
    subprocess.run(["rustc", "--edition", "2024", "-D", "warnings", "-C", "opt-level=3",
                    "benches/path/cocycle.rs", "--extern", "cocycle=target/release/libcocycle.rlib",
                    "-L", "dependency=target/release/deps", "-o", str(worker)], cwd=ROOT, check=True)
    return worker


def invoke(worker, engine, fixture, args):
    values = [engine, fixture["vertices"], args.field, args.mode, args.max_paths, args.max_entries, args.max_work]
    values += [v for edge, birth in zip(fixture["edges"], fixture["weights"]) for v in (*edge, birth)]
    started = time.perf_counter_ns()
    try:
        run = subprocess.run([str(worker)], input=" ".join(map(str, values)), capture_output=True,
                             text=True, encoding="utf-8", timeout=args.timeout)
    except subprocess.TimeoutExpired:
        return {"status": "timeout", "timeout_seconds": args.timeout}
    wall = time.perf_counter_ns() - started
    if run.returncode:
        return {"status": "error", "exit_code": run.returncode, "stderr": run.stderr, "stdout": run.stdout}
    try:
        result = json.loads(run.stdout)
    except ValueError:
        return {"status": "protocol_error", "stdout": run.stdout, "stderr": run.stderr}
    if (result.get("protocol"), result.get("engine"), result.get("field"), result.get("mode"),
        result.get("dimension"), result.get("coverage")) != (PROTOCOL, engine, args.field, args.mode, 1, "complete"):
        return {"status": "protocol_error", "result": result}
    return {"status": "ok", "wall_ns": wall, **result}


def summarize(samples, seed):
    g, d = samples["generic"], samples["direct"]
    ratio = statistics.median(g)/statistics.median(d)
    rng = random.Random(seed)
    bootstrap = []
    for _ in range(2000):
        indices = [rng.randrange(len(g)) for _ in g]
        bootstrap.append(statistics.median([g[i] for i in indices])/statistics.median([d[i] for i in indices]))
    bootstrap.sort()
    lo, hi = bootstrap[int(0.025*(len(bootstrap)-1))], bootstrap[int(0.975*(len(bootstrap)-1))]
    return {"samples_per_engine": len(g), "generic_ns": {"median": statistics.median(g), "min": min(g), "max": max(g)},
            "direct_ns": {"median": statistics.median(d), "min": min(d), "max": max(d)},
            "generic_over_direct": ratio, "paired_bootstrap_95_percentile": [lo, hi],
            "classification": "direct_advantage" if ratio >= 1.1 and lo > 1 else
                              "direct_regression" if ratio <= 0.9 and hi < 1 else "inconclusive"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--samples", type=int, default=10)
    parser.add_argument("--seeds", type=int, nargs="+", default=[20261125, 20261126])
    parser.add_argument("--field", type=int, default=2)
    parser.add_argument("--mode", choices=["regular", "full"], default="regular")
    parser.add_argument("--timeout", type=float, default=15)
    parser.add_argument("--max-paths", type=int, default=250000)
    parser.add_argument("--max-entries", type=int, default=2000000)
    parser.add_argument("--max-work", type=int, default=250000000)
    args = parser.parse_args()
    if args.samples < (1 if args.quick else 10) or min(args.timeout, args.max_paths, args.max_entries, args.max_work) <= 0:
        parser.error("positive budgets and at least ten formal measured processes are required")
    before = fingerprint()
    head = command(["git", "rev-parse", "HEAD"])
    dirty = command(["git", "status", "--porcelain", "--", *before])
    if dirty and not args.quick:
        parser.error("formal runs require committed clean implementation and harness inputs")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "fixtures").mkdir()
    worker = build(output)
    protocol = {"protocol": PROTOCOL, "evidence": "smoke" if args.quick else "comparative-performance-study",
                "candidate_commit": head, "baseline_commit": head, "candidate": "direct", "baseline": "generic",
                "dirty": bool(dirty), "sources": before, "worker_sha256": digest(worker), "samples": args.samples,
                "environment": {"os": platform.platform(), "cpu": platform.processor(), "logical_cpus": os.cpu_count(),
                                "python": platform.python_version(), "rustc": command(["rustc", "-vV"]),
                                "affinity": "inherited/uncontrolled", "host_load_and_frequency": "uncontrolled"},
                "parameters": {k: v for k, v in vars(args).items() if k != "output"}, "cells": []}
    work = list(fixtures(args.quick, args.seeds[:1] if args.quick else args.seeds))
    random.Random(20261125).shuffle(work)
    with (output / "samples.jsonl").open("w", encoding="utf-8") as raw:
        for fixture in work:
            encoded = json.dumps(fixture, sort_keys=True).encode()
            (output / "fixtures" / (fixture["id"] + ".json")).write_bytes(encoded + b"\n")
            cell = {k: v for k, v in fixture.items() if k not in ("edges", "weights")}
            cell.update({"edge_count": len(fixture["edges"]), "input_sha256": hashlib.sha256(encoded).hexdigest(), "status": "ok"})
            expected, samples = None, {"generic": [], "direct": []}
            for round_index in range(-1, args.samples):
                for engine in (("generic", "direct") if round_index % 2 == 0 else ("direct", "generic")):
                    result = invoke(worker, engine, fixture, args)
                    raw.write(json.dumps({"case": fixture["id"], "round": round_index, "engine": engine, **result}) + "\n")
                    raw.flush()
                    if result["status"] != "ok":
                        cell["status"] = "incomplete"
                        cell.setdefault("failures", []).append({"engine": engine, "round": round_index, **result})
                        continue
                    if expected is None:
                        expected = result["intervals"]
                    elif result["intervals"] != expected:
                        cell["status"] = "mismatch"
                    if round_index >= 0:
                        samples[engine].append(result["elapsed_ns"])
                if cell["status"] != "ok":
                    break
            if cell["status"] == "ok":
                cell.update(summarize(samples, fixture["id"]))
            protocol["cells"].append(cell)
            protocol["sources_stable"] = fingerprint() == before
            (output / "summary.json").write_text(json.dumps(protocol, indent=2) + "\n", encoding="utf-8")
            print(json.dumps({"case": fixture["id"], "status": cell["status"], "ratio": cell.get("generic_over_direct")}), flush=True)
    if not protocol["sources_stable"] or any(c["status"] == "mismatch" for c in protocol["cells"]):
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
