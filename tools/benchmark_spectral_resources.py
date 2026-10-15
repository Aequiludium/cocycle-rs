"""Repeated serial Windows spectral resource snapshots; no external timing ranking.

The JSON specification is a list of commands and explicit fixture descriptions.
Native workers report preparation and solve separately. Process wall time and
Windows peak working set include startup, parsing, runtime and retained outputs.
Generated specifications, samples, logs and manifests stay outside source Git.
"""
import argparse
import ctypes
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import time


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def peak_working_set(process):
    from ctypes import wintypes
    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("faults", wintypes.DWORD)] + [
            (name, ctypes.c_size_t) for name in ("peak_ws", "ws", "peak_paged",
             "paged", "peak_nonpaged", "nonpaged", "pagefile", "peak_pagefile")]
    counter = Counters()
    counter.cb = ctypes.sizeof(counter)
    query = ctypes.windll.psapi.GetProcessMemoryInfo
    query.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD]
    if not query(process._handle, ctypes.byref(counter), counter.cb):
        raise OSError("Windows process peak working set unavailable")
    return counter.peak_ws


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=10)
    args = parser.parse_args()
    if sys.platform != "win32" or args.samples < 10:
        parser.error("this snapshot requires Windows and at least ten measured processes per cell")
    root = Path(__file__).resolve().parents[1]
    changed = subprocess.check_output(["git", "diff", "HEAD", "--name-only"], cwd=root, text=True)
    if any(name.endswith((".rs", ".py", ".toml", ".lock", ".yml")) for name in changed.splitlines()):
        raise ValueError("commit measured code and harness before measuring")
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    specs = json.loads(args.specification.read_text(encoding="utf-8"))
    (output / "specification.json").write_bytes(args.specification.read_bytes())
    identity = {
        "evidence_class": "resource_snapshot", "protocol": "spectral-resource-v1",
        "candidate_sha": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip(),
        "source_fingerprints": {str(p.relative_to(root)): digest(p) for p in root.glob("src/**/*.rs")},
        "worker_sources": {str(p.relative_to(root)): digest(p) for p in
                           [*root.glob("tools/reference/*.rs"), *root.glob("examples/*.rs")]},
        "controller_sha256": digest(__file__), "fixture_sha256": digest(args.specification),
        "platform": platform.platform(), "processor": platform.processor(), "python": sys.version,
        "toolchains": {name: subprocess.check_output([name, "--version"], text=True).strip() for name in ("cargo", "rustc")},
        "samples_per_cell": args.samples, "discarded_warmups": 1, "order_seed": 73,
        "build_flags": {spec["name"]: spec.get("build_flags", "rustc --edition 2024 -C opt-level=3; cargo --locked --release") for spec in specs},
        "memory": "absolute process peak working set bytes; includes runtime/input/allocator/output, not allocation peak",
        "controls": "serial child processes; no concurrent build/tests; affinity/frequency/host load uncontrolled",
        "environment": {key: os.environ.get(key) for key in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "commands": specs,
    }
    (output / "environment.json").write_text(json.dumps(identity, indent=2) + "\n", encoding="utf-8")
    rows = []
    rng = random.Random(73)
    for repeat in range(-1, args.samples):
        order = list(range(len(specs)))
        rng.shuffle(order)
        for position, index in enumerate(order):
            spec = specs[index]
            command = [str(value).replace("{attempt}", str(output / f"attempt-{index}-{repeat}"))
                       for value in spec["command"]]
            begin = time.perf_counter_ns()
            process = subprocess.Popen(command, cwd=spec.get("cwd", str(root)),
                                       stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                       stderr=subprocess.PIPE, text=True)
            stdout, stderr = process.communicate(spec.get("stdin", ""))
            wall_ns = time.perf_counter_ns() - begin
            peak = peak_working_set(process)
            stem = output / f"cell-{index}-sample-{repeat}"
            stem.with_suffix(".stdout").write_text(stdout, encoding="utf-8")
            stem.with_suffix(".stderr").write_text(stderr, encoding="utf-8")
            row = {"cell": spec["name"], "repeat": repeat, "position": position,
                   "warmup": repeat == -1, "process_wall_ns": wall_ns,
                   "peak_working_set_bytes": peak, "exit_code": process.returncode,
                   "executable_sha256": digest(command[0]), "status": "process_failure"}
            if process.returncode == 0:
                if spec["kind"] == "persistent":
                    words = stdout.split()
                    row.update(status="success", solve_ns=float(words[1])*1e9,
                               checksum=float(words[3]))
                    row["direct_fallbacks"] = int(words[5])
                    if abs(row["checksum"] - spec["expected_checksum"]) > 1e-7:
                        row["status"] = "validation_failure"
                elif spec["kind"] == "reproduction":
                    value = json.loads(stdout.splitlines()[-1])
                    row.update(status=value["status"])
                    if value["core_comparisons"] != 160:
                        row["status"] = "validation_failure"
                else:
                    value = json.loads(stdout)
                    row.update(value)
                    if row["status"] == "success":
                        import numpy as np
                        expected = np.asarray(spec["expected_values"])
                        actual = np.asarray(row["values"])
                        if (actual.shape != expected.shape or not np.isfinite(actual).all()
                                or not np.allclose(actual, expected, atol=spec["atol"], rtol=spec["rtol"])
                                or row["index_coverage_certified"]):
                            row["status"] = "validation_failure"
            rows.append(row)
            (output / "samples.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
            print(json.dumps({key: row[key] for key in ("cell", "repeat", "status", "process_wall_ns", "peak_working_set_bytes")}), flush=True)
    summary = []
    for spec in specs:
        selected = [row for row in rows if row["cell"] == spec["name"] and not row["warmup"]]
        item = {"cell": spec["name"], "samples": len(selected),
                "outcomes": {status: sum(row["status"] == status for row in selected)
                             for status in sorted(set(row["status"] for row in selected))}}
        for field in ("process_wall_ns", "prepare_ns", "solve_ns", "peak_working_set_bytes"):
            values = [row[field] for row in selected if field in row]
            if values:
                item[field] = {"median": statistics.median(values), "min": min(values), "max": max(values)}
        summary.append(item)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output / "manifest.json").write_text(json.dumps({p.name: digest(p) for p in output.iterdir() if p.is_file()}, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
