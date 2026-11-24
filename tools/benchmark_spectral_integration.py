"""Compare unaccelerated and alternative spectral workflows in one committed build.

A specification supplies two explicit routes per cell, frozen inputs and an
independent oracle for each route. Original-PR-versus-integrated comparisons are
rejected: those answer an integration-regression question, not acceleration.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import io
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import statistics
import subprocess
import sys
import tarfile
import time

from benchmark_spectral_resources import peak_working_set


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def validate_values(actual, expected, tolerance):
    """Use the independently fixed Frobenius/Euclidean error contract."""
    if not isinstance(actual, list) or len(actual) != len(expected):
        return False
    if not all(type(v) in (float, int) and math.isfinite(v) for v in [*actual, *expected]):
        return False
    error = math.hypot(*(a - b for a, b in zip(actual, expected)))
    scale = max(1., math.hypot(*expected))
    return math.isfinite(error) and error / scale <= tolerance


def validate_specification(spec):
    if not spec["cells"] or len({c["name"] for c in spec["cells"]}) != len(spec["cells"]):
        raise ValueError("the specification needs nonempty, uniquely named cells")
    for cell in spec["cells"]:
        if "baseline" in cell or set(cell.get("routes", {})) != {"baseline", "candidate"}:
            raise ValueError("declare unaccelerated and alternative routes in the same build; PR-head baselines are not accepted")
        if cell["routes"]["baseline"] == cell["routes"]["candidate"]:
            raise ValueError("the two algorithm routes must differ")
        if not cell.get("baseline_description") or not cell.get("candidate_description"):
            raise ValueError("both algorithm routes need explicit descriptions")
        if not 0 < cell.get("tolerance", 1e-10) <= 1e-7:
            raise ValueError("declare the independent error contract before measuring")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=10)
    args = parser.parse_args()
    if sys.platform != "win32" or args.samples < 10:
        parser.error("Windows and at least ten paired measured processes per cell are required")
    spec = json.loads(args.specification.read_text(encoding="utf-8"))
    validate_specification(spec)
    root = Path(__file__).resolve().parents[1]
    if subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=root):
        raise ValueError("commit measured sources and controller before measuring")
    candidate = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / "specification.json").write_bytes(args.specification.read_bytes())
    directory = output / "build"
    directory.mkdir()
    archive = subprocess.check_output(["git", "archive", candidate, "src", "examples", "tools/reference"], cwd=root)
    with tarfile.open(fileobj=io.BytesIO(archive)) as bundle:
        for member in bundle:
            if member.isfile():
                path = (directory / member.name).resolve()
                if not path.is_relative_to(directory):
                    raise ValueError("source archive escapes the output directory")
                path.parent.mkdir(parents=True, exist_ok=True)
                path.write_bytes(bundle.extractfile(member).read())
    library = directory / "libcocycle.rlib"
    command = ["rustc", "--edition", "2024", "--crate-name", "cocycle", "--crate-type", "rlib",
               "-C", "opt-level=3", str(directory / "src/lib.rs"), "-o", str(library)]
    with (directory / "library-build.log").open("w") as log:
        subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT, check=True)
    builds = {}

    def build(source, cfg):
        key = (source, tuple(cfg))
        if key not in builds:
            original = directory / source
            text = original.read_text(encoding="utf-8")
            text = re.sub(r'(#\[path\s*=\s*")([^"]+)("\])',
                          lambda m: m[1] + (original.parent / m[2]).resolve().as_posix() + m[3], text)
            text, count = re.subn(r"\bfn main\(", "fn workload_main(", text)
            if count != 1:
                raise ValueError("worker must have exactly one main function")
            text += '\nfn main() -> std::result::Result<(), Box<dyn std::error::Error>> {\n'
            text += '    let start = std::time::Instant::now();\n    let result = workload_main();\n'
            text += '    std::io::Write::flush(&mut std::io::stdout())?;\n'
            text += '    eprintln!("INTEGRATION_TIMING_NS {}", start.elapsed().as_nanos());\n    result\n}\n'
            stem = hashlib.sha256(repr(key).encode()).hexdigest()[:16]
            instrumented, worker = directory / (stem + ".rs"), directory / (stem + ".exe")
            instrumented.write_text(text, encoding="utf-8", newline="\n")
            command = ["rustc", "--crate-name", "spectral_comparison_worker", "--edition", "2024", "-C", "opt-level=3",
                       str(instrumented), "--extern", "cocycle=" + str(library), "-o", str(worker)]
            for name in cfg:
                command.extend(["--cfg", name])
            with (directory / (stem + "-build.log")).open("w") as log:
                subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT, check=True)
            builds[key] = dict(worker=str(worker), source=source, original_sha256=digest(original),
                               instrumented_sha256=digest(instrumented), executable_sha256=digest(worker), command=command)
        return builds[key]

    cells, oracles = [], {}
    for item in spec["cells"]:
        cell = json.loads(json.dumps(item))
        oracle = Path(cell["oracle_file"])
        frozen = output / (cell["name"] + "-oracle.json")
        frozen.write_bytes(oracle.read_bytes())
        cell.update(oracle_file=str(frozen), oracle_sha256=digest(frozen))
        oracles[cell["name"]] = json.loads(frozen.read_text(encoding="utf-8"))
        for route, data in cell["routes"].items():
            data["build"] = build(data["source"], data.get("cfg", []))
            for field in ["stdin_file", "argument_file"]:
                if field in data:
                    frozen = output / f"{cell['name']}-{route}-{field}.txt"
                    frozen.write_bytes(Path(data[field]).read_bytes())
                    data[field], data[field + "_sha256"] = str(frozen), digest(frozen)
        cells.append(cell)
    metadata = dict(candidate_sha=candidate, controller_sha256=digest(__file__), cells=cells,
                    protocol="spectral-unaccelerated-comparison-v2", evidence_class="resource_snapshot",
                    platform=platform.platform(), processor=platform.processor(), python=sys.version,
                    rustc=subprocess.check_output(["rustc", "--version"], text=True).strip(),
                    samples=args.samples, warmups=1, order_seed=160,
                    controls="serial processes after builds; affinity/frequency/host load uncontrolled",
                    worker_timer="main entry through fixture I/O, parsing, complete workflow, checks and flushed output; excludes startup",
                    process_timer="Popen through communicate; includes startup and input/output transport",
                    memory="absolute Windows process peak working set, not allocated/live algorithm bytes",
                    measurement_started_utc=datetime.now(timezone.utc).isoformat(),
                    thread_environment={k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")})
    (output / "environment.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    rows, rng = [], random.Random(160)
    offsets = {cell["name"]: rng.randrange(2) for cell in cells}
    for repeat in range(-1, args.samples):
        order = list(cells)
        rng.shuffle(order)
        for cell in order:
            routes = ["baseline", "candidate"]
            if (repeat + offsets[cell["name"]]) % 2:
                routes.reverse()
            for route in routes:
                data = cell["routes"][route]
                command = [data["build"]["worker"]]
                if "argument_file" in data:
                    command.append(data["argument_file"])
                command.extend(data.get("args", []))
                stdin = Path(data["stdin_file"]).read_text(encoding="utf-8") if "stdin_file" in data else ""
                begin = time.perf_counter_ns()
                process = subprocess.Popen(command, cwd=root, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                timed_out = False
                try:
                    stdout, stderr = process.communicate(stdin, timeout=120)
                except subprocess.TimeoutExpired:
                    timed_out = True
                    process.kill()
                    stdout, stderr = process.communicate()
                wall, peak = time.perf_counter_ns() - begin, peak_working_set(process)
                stem = output / f"{cell['name']}-{repeat}-{route}"
                stem.with_suffix(".stdout").write_text(stdout, encoding="utf-8")
                stem.with_suffix(".stderr").write_text(stderr, encoding="utf-8")
                match = re.findall(r"^INTEGRATION_TIMING_NS (\d+)$", stderr, re.MULTILINE)
                row = dict(cell=cell["name"], repeat=repeat, route=route, warmup=repeat == -1,
                           process_wall_ns=wall, peak_working_set_bytes=peak, exit_code=process.returncode,
                           timed_out=timed_out, worker_ns=int(match[0]) if len(match) == 1 else None)
                try:
                    result = json.loads(stdout)
                    expected = oracles[cell["name"]][route]["values"]
                    row["validated"] = (not timed_out and process.returncode == 0 and len(match) == 1
                                        and validate_values(result["values"], expected, cell.get("tolerance", 1e-10)))
                    for name in ["prepare_ns", "solve_ns", "actions", "iterations"]:
                        if name in result:
                            row[name] = result[name]
                    system = oracles[cell["name"]].get("linear_system")
                    if system is not None:
                        residual = [sum(a*b for a,b in zip(matrix_row, result["values"])) - rhs
                                    for matrix_row, rhs in zip(system["matrix"], system["rhs"])]
                        row["relative_residual"] = math.sqrt(sum(v*v for v in residual)) / max(1., math.sqrt(sum(v*v for v in system["rhs"])))
                        row["validated"] &= row["relative_residual"] <= 1e-9
                except (ValueError, KeyError, TypeError):
                    row["validated"] = False
                rows.append(row)
                (output / "samples.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
                if not row["validated"]:
                    raise ValueError(f"retained worker failure or independent-oracle mismatch: {row}")
        print(f"completed round {repeat}", flush=True)
    summary = []
    for cell in cells:
        item = dict(cell=cell["name"], samples_per_route=args.samples, comparison=cell.get("comparison", "equivalent"),
                    baseline=cell["baseline_description"], candidate=cell["candidate_description"],
                    quality=oracles[cell["name"]].get("quality"))
        selected = [r for r in rows if r["cell"] == cell["name"] and not r["warmup"]]
        for field in ["worker_ns", "process_wall_ns", "peak_working_set_bytes", "prepare_ns", "solve_ns", "actions", "iterations"]:
            if not all(field in r for r in selected):
                continue
            item[field] = {route: dict(median=statistics.median(r[field] for r in selected if r["route"] == route),
                                      min=min(r[field] for r in selected if r["route"] == route),
                                      max=max(r[field] for r in selected if r["route"] == route)) for route in ["baseline", "candidate"]}
            if all(r[field] > 0 for r in selected):
                item[field]["paired_median_baseline_over_candidate"] = statistics.median(
                    next(r[field] for r in selected if r["repeat"] == i and r["route"] == "baseline") /
                    next(r[field] for r in selected if r["repeat"] == i and r["route"] == "candidate") for i in range(args.samples))
        summary.append(item)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    metadata["measurement_finished_utc"] = datetime.now(timezone.utc).isoformat()
    (output / "environment.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    (output / "manifest.json").write_text(json.dumps({str(p.relative_to(output)): digest(p) for p in output.rglob("*") if p.is_file()}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(cells=len(cells), measured_processes=len(cells)*args.samples*2, output=str(output))))


if __name__ == "__main__":
    main()
