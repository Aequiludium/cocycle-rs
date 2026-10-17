"""Source-bound Linux/WSL measurements of owned F2 families and full histories."""

import argparse
from itertools import combinations
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import time

from check_homology_operators import git

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = "cocycle-f2-family-resources-v1"


def fixture(kind, vertices):
    if vertices < 3:
        raise ValueError("at least three vertices required")
    cells = [(0, (i,)) for i in range(vertices)]
    if kind == "graph":
        cells += [(1, edge) for edge in combinations(range(vertices), 2)]
    elif kind == "fan":
        cells += [(1, (0, i)) for i in range(1, vertices)]
        cells += [(1, (i, i + 1)) for i in range(1, vertices - 1)]
        cells += [(2, (0, i, i + 1)) for i in range(1, vertices - 1)]
    else:
        raise ValueError("unknown family fixture")
    cells.sort(key=lambda x: (x[0], len(x[1]), tuple(-v for v in reversed(x[1]))))
    return cells


def wire(cells):
    values = [1, len(cells)]
    for level, cell in cells:
        values += [level, len(cell), *cell]
    return " ".join(map(str, values)) + "\n"


def insert(pivots, vector):
    while vector:
        pivot = vector.bit_length() - 1
        if pivot not in pivots:
            pivots[pivot] = vector
            return
        vector ^= pivots[pivot]


def contains(pivots, vector):
    while vector:
        pivot = vector.bit_length() - 1
        if pivot not in pivots:
            return False
        vector ^= pivots[pivot]
    return True


def oracle(cells):
    """Independent raw-simplex PH and boundary spans; never a runtime result."""
    names = [cell for _, cell in cells]
    indices = {cell: i for i, cell in enumerate(names)}
    boundaries = [sum(1 << indices[face] for face in combinations(cell, len(cell)-1))
                  if len(cell) > 1 else 0 for cell in names]
    reduced, pivots, births, deaths, betti, active = [], {}, [], {}, [], 0
    spaces, span = [], {}
    for j, ((_, cell), boundary) in enumerate(zip(cells, boundaries)):
        column = boundary
        while column and column.bit_length()-1 in pivots:
            column ^= reduced[pivots[column.bit_length()-1]]
        reduced.append(column)
        if not column and len(cell) == 2:
            births.append(j)
            active += 1
        elif column:
            pivot = column.bit_length()-1
            pivots[pivot] = j
            if len(cell) == 3:
                deaths[pivot] = j
                active -= 1
        if len(cell) == 3:
            insert(span, boundary)
        spaces.append(span.copy())
        betti.append(active)
    bars = [[b, deaths.get(b)] for b in births]
    bars.sort(key=lambda x: (x[0], len(cells) if x[1] is None else x[1]))
    return betti, bars, boundaries, spaces


def validate_result(result, cells, mode):
    if result.get("state") != "Computed" or result.get("protocol") != PROTOCOL or result.get("mode") != mode:
        raise ValueError("worker did not complete the requested protocol")
    betti, bars, boundaries, spaces = oracle(cells)
    if result["stages"] != len(cells) or result["betti"] != betti:
        raise ValueError("stage Betti sequence differs from independent PH")
    if mode == "construct":
        if result["bars"] or result["histories"]:
            raise ValueError("construction unexpectedly requested feature output")
        return
    if result["bars"] != bars or len(result["histories"]) != len(bars):
        raise ValueError("historical intervals differ from independent PH")
    for (birth, death), history in zip(bars, result["histories"]):
        end = len(cells) if death is None else death
        if len(history) != end-birth:
            raise ValueError("incomplete historical output")
        previous = None
        for stage, vector in enumerate(history, birth):
            basis = [j for j, (_, cell) in enumerate(cells[:stage+1]) if len(cell) == 2]
            if len(vector) != len(basis) or any(bit not in (0, 1) for bit in vector):
                raise ValueError("invalid historical coordinates")
            chain, boundary = 0, 0
            for j, bit in zip(basis, vector):
                if bit:
                    chain ^= 1 << j
                    boundary ^= boundaries[j]
            if boundary or (previous is not None and not contains(spaces[stage], chain ^ previous)):
                raise ValueError("history violates raw cycle or transport equations")
            previous = chain
        if death is not None and not contains(spaces[death], previous):
            raise ValueError("history does not die in raw boundary space")


def fingerprint(root):
    paths = list((root / "src").rglob("*.rs")) + [root / "Cargo.toml", root / "Cargo.lock"]
    value = hashlib.sha256()
    for path in sorted(paths):
        value.update(path.relative_to(root).as_posix().encode())
        value.update(path.read_bytes())
    return value.hexdigest()


def main(args):
    if platform.system() != "Linux":
        raise ValueError("family resource measurements require Linux/WSL")
    import resource
    if args.output.exists() or args.samples < 1:
        raise ValueError("positive sample count and fresh output directory required")
    args.output.mkdir(parents=True)
    roots = {"candidate": ROOT, "baseline": args.baseline_root.resolve()}
    workers = {"candidate": args.worker.resolve(), "baseline": args.baseline_worker.resolve()}
    metadata = {"protocol": PROTOCOL, "samples": args.samples, "seed": args.seed,
                "python": platform.python_version(), "platform": platform.platform(),
                "cpu": Path("/proc/cpuinfo").read_text().split("model name", 1)[-1].splitlines()[0],
                "affinity": sorted(os.sched_getaffinity(0)), "frequency_and_host_load": "uncontrolled",
                "timeout_s": args.timeout, "address_space_mib": args.memory_mib,
                "compiler": subprocess.check_output(["rustc", "--version"], text=True).strip(),
                "build": "cargo build --locked --release; rustc --edition 2024 -C opt-level=3",
                "started_unix": time.time(), "implementations": {}}
    for engine, root in roots.items():
        git(root, "diff", "--exit-code", "HEAD", "--", "src", "Cargo.toml", "Cargo.lock")
        metadata["implementations"][engine] = {"sha": git(root, "rev-parse", "HEAD").strip(),
            "source_sha256": fingerprint(root), "worker_sha256": hashlib.sha256(workers[engine].read_bytes()).hexdigest()}
    metadata["harness_sha256"] = {str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in (Path(__file__).resolve(), ROOT / "tools/reference/homology_operator.rs")}
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def call(engine, mode, cells, path):
        command = [str(workers[engine]), "family-resources", mode]
        limit = args.memory_mib * 1024 * 1024
        try:
            proc = subprocess.run(command, input=wire(cells), text=True, capture_output=True,
                timeout=args.timeout, preexec_fn=lambda: resource.setrlimit(resource.RLIMIT_AS, (limit, limit)))
            raw = {"command": command, "returncode": proc.returncode, "stdout": proc.stdout, "stderr": proc.stderr}
            path.write_text(json.dumps(raw, indent=2) + "\n")
            if proc.returncode:
                raise RuntimeError(f"failed worker retained at {path}")
            result = json.loads(proc.stdout)
            validate_result(result, cells, mode)
            if any(result[key] is None for key in ("baseline_rss_kib", "baseline_hwm_kib", "after_rss_kib", "peak_rss_kib")):
                raise ValueError("Linux RSS/HWM measurements required")
            return result
        except subprocess.TimeoutExpired as error:
            path.write_text(json.dumps({"command": command, "state": "Timeout", "timeout": error.timeout}) + "\n")
            raise

    summary = {"protocol": PROTOCOL, "classification": "ComparativePerformance" if args.samples >= 10 else "ResourceSnapshot", "cells": []}
    for kind, sizes in (("fan", args.fan_sizes), ("graph", args.graph_sizes)):
        for vertices in map(int, sizes.split(",")):
            cells = fixture(kind, vertices)
            text = wire(cells)
            stem = f"{kind}-{vertices}"
            (args.output / f"{stem}.fixture.txt").write_text(text)
            for mode in ("construct", "features"):
                records = {engine: [] for engine in workers}
                for engine in workers:
                    call(engine, mode, cells, args.output / f"{stem}-{mode}-{engine}-warmup.json")
                for sample in range(args.samples):
                    order = list(workers)
                    random.Random(args.seed + sample//2).shuffle(order)
                    if sample%2:
                        order.reverse()
                    outputs = {}
                    for engine in order:
                        outputs[engine] = call(engine, mode, cells, args.output / f"{stem}-{mode}-{engine}-{sample}.json")
                        records[engine].append(outputs[engine])
                    for key in ("stages", "betti", "bars", "histories"):
                        if outputs["candidate"][key] != outputs["baseline"][key]:
                            raise ValueError(f"before/after mismatch in {key}")
                row = {"family": kind, "vertices": vertices, "simplices": len(cells), "mode": mode,
                       "fixture_sha256": hashlib.sha256(text.encode()).hexdigest(), "engines": {}}
                for engine, values in records.items():
                    row["engines"][engine] = {key: {"median": statistics.median(v[key] for v in values),
                        "min": min(v[key] for v in values), "max": max(v[key] for v in values)}
                        for key in ("elapsed_ns", "baseline_rss_kib", "baseline_hwm_kib", "after_rss_kib", "peak_rss_kib")}
                summary["cells"].append(row)
                (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
                print(f"{stem} {mode}: {len(cells)} stages; exact before/after and raw oracle agreement", flush=True)
    metadata["finished_unix"] = time.time()
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    parser.add_argument("--baseline-worker", type=Path, required=True)
    parser.add_argument("--baseline-root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=10)
    parser.add_argument("--seed", type=int, default=20261017)
    parser.add_argument("--timeout", type=float, default=30)
    parser.add_argument("--memory-mib", type=int, default=512)
    parser.add_argument("--fan-sizes", default="4,8,16,32,64")
    parser.add_argument("--graph-sizes", default="4,6,8,12,16")
    main(parser.parse_args())
