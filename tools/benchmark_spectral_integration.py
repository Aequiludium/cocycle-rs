"""Paired Windows resource snapshots of committed spectral verification workers.

Supply previously validated inputs in a JSON specification with a ``cells`` list.
Each cell names a baseline commit, worker source, arguments and optional stdin or
argument file. Both revisions use identical compiler flags and input bytes.
"""
import argparse
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


def equivalent(expected, actual):
    """Preserve labels, integer identities and numerical outputs; omit timings."""
    def check(a, b):
        if isinstance(a, dict) and isinstance(b, dict):
            a = {k: v for k, v in a.items() if not k.endswith(("_ns", "_seconds"))}
            b = {k: v for k, v in b.items() if not k.endswith(("_ns", "_seconds"))}
            return a.keys() == b.keys() and all(check(a[k], b[k]) for k in a)
        if isinstance(a, list) and isinstance(b, list):
            return len(a) == len(b) and all(check(x, y) for x, y in zip(a, b))
        if type(a) is float or type(b) is float:
            return (type(a) in (int, float) and type(b) in (int, float)
                    and math.isfinite(a) and math.isfinite(b)
                    and math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-10))
        return type(a) is type(b) and a == b
    try:
        return check(json.loads(expected), json.loads(actual))
    except json.JSONDecodeError:
        # Line-oriented worker protocols have stable labels and numeric tokens.
        number = r"[-+]?\d+(?:\.\d*)?(?:[eE][-+]?\d+)?"
        if re.split(number, expected) != re.split(number, actual):
            return False
        a, b = re.findall(number, expected), re.findall(number, actual)
        return len(a) == len(b) and all(
            x == y if re.fullmatch(r"[-+]?\d+", x) and re.fullmatch(r"[-+]?\d+", y)
            else math.isclose(float(x), float(y), rel_tol=1e-10, abs_tol=1e-10)
            for x, y in zip(a, b))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--specification", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--samples", type=int, default=10)
    args = parser.parse_args()
    if sys.platform != "win32" or args.samples < 10:
        parser.error("Windows and at least ten paired measured processes per cell are required")
    root = Path(__file__).resolve().parents[1]
    if subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], cwd=root):
        raise ValueError("commit the candidate and controller before measuring")
    candidate = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    spec = json.loads(args.specification.read_text(encoding="utf-8"))
    if not spec["cells"] or len({c["name"] for c in spec["cells"]}) != len(spec["cells"]):
        raise ValueError("the specification needs nonempty, uniquely named cells")
    (output / "specification.json").write_bytes(args.specification.read_bytes())
    builds, sources = {}, {}

    def build(revision, source, cfg):
        revision = subprocess.check_output(["git", "rev-parse", revision + "^{commit}"], cwd=root, text=True).strip()
        if revision not in sources:
            directory = output / "build" / revision
            directory.mkdir(parents=True)
            archive = subprocess.check_output(["git", "archive", revision, "src", "examples", "tools/reference"], cwd=root)
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
            sources[revision] = directory
        key = (revision, source, tuple(cfg))
        if key not in builds:
            directory = sources[revision]
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
            command = ["rustc", "--crate-name", "spectral_integration_worker", "--edition", "2024", "-C", "opt-level=3",
                       str(instrumented), "--extern", "cocycle=" + str(directory / "libcocycle.rlib"), "-o", str(worker)]
            for name in cfg:
                command.extend(["--cfg", name])
            with (directory / (stem + "-build.log")).open("w") as log:
                subprocess.run(command, cwd=root, stdout=log, stderr=subprocess.STDOUT, check=True)
            builds[key] = dict(revision=revision, worker=str(worker), source=source,
                               original_sha256=digest(original), instrumented_sha256=digest(instrumented),
                               executable_sha256=digest(worker), command=command)
        return builds[key]

    cells = []
    for item in spec["cells"]:
        cell = dict(item)
        cell["workers"] = {route: build(revision, cell["source"], cell.get("cfg", []))
                           for route, revision in [("baseline", cell["baseline"]), ("candidate", candidate)]}
        for field in ["stdin_file", "argument_file"]:
            if field in cell:
                path = Path(cell[field]).resolve()
                frozen = output / (cell["name"] + "-" + field + ".txt")
                frozen.write_bytes(path.read_bytes())
                cell[field] = str(frozen)
                cell[field + "_sha256"] = digest(frozen)
        cells.append(cell)
    metadata = dict(candidate_sha=candidate, controller_sha256=digest(__file__), cells=cells,
                    protocol="spectral-integration-resource-v1", evidence_class="resource_snapshot",
                    platform=platform.platform(), processor=platform.processor(), python=sys.version,
                    rustc=subprocess.check_output(["rustc", "--version"], text=True).strip(),
                    samples=args.samples, warmups=1, order_seed=160,
                    controls="serial processes; no concurrent build/tests; affinity/frequency/host load uncontrolled",
                    worker_timer="main entry through fixture I/O, parsing, computation, verification and flushed output; excludes process startup",
                    process_timer="Popen through communicate; includes startup and input/output transport",
                    memory="absolute Windows process peak working set, not allocated/live algorithm bytes",
                    tolerances="identical labels/integers; float abs/rel 1e-10; JSON *_ns/*_seconds omitted",
                    thread_environment={k: os.environ.get(k) for k in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")})
    (output / "environment.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    rows, expected = [], {}
    rng = random.Random(160)
    offsets = {cell["name"]: rng.randrange(2) for cell in cells}
    for repeat in range(-1, args.samples):
        order = list(cells)
        rng.shuffle(order)
        for cell in order:
            routes = ["baseline", "candidate"]
            if (repeat + offsets[cell["name"]]) % 2:
                routes.reverse()
            for route in routes:
                command = [cell["workers"][route]["worker"], *cell.get("args", [])]
                if "argument_file" in cell:
                    command.append(cell["argument_file"])
                data = Path(cell["stdin_file"]).read_text(encoding="utf-8") if "stdin_file" in cell else ""
                begin = time.perf_counter_ns()
                process = subprocess.Popen(command, cwd=root, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                timed_out = False
                try:
                    stdout, stderr = process.communicate(data, timeout=120)
                except subprocess.TimeoutExpired:
                    timed_out = True
                    process.kill()
                    stdout, stderr = process.communicate()
                wall = time.perf_counter_ns() - begin
                peak = peak_working_set(process)
                stem = output / f"{cell['name']}-{repeat}-{route}"
                stem.with_suffix(".stdout").write_text(stdout, encoding="utf-8")
                stem.with_suffix(".stderr").write_text(stderr, encoding="utf-8")
                match = re.findall(r"^INTEGRATION_TIMING_NS (\d+)$", stderr, re.MULTILINE)
                row = dict(cell=cell["name"], repeat=repeat, route=route, warmup=repeat == -1,
                           process_wall_ns=wall, peak_working_set_bytes=peak, exit_code=process.returncode, timed_out=timed_out,
                           worker_ns=int(match[0]) if len(match) == 1 else None)
                if cell["name"] not in expected:
                    expected[cell["name"]] = stdout
                row["validated"] = not timed_out and process.returncode == 0 and len(match) == 1 and equivalent(expected[cell["name"]], stdout)
                rows.append(row)
                (output / "samples.json").write_text(json.dumps(rows, indent=2) + "\n", encoding="utf-8")
                if not row["validated"]:
                    raise ValueError(f"retained worker failure or output mismatch: {row}")
        print(f"completed round {repeat}", flush=True)
    summary = []
    for cell in cells:
        item = dict(cell=cell["name"], baseline=cell["baseline"], samples_per_route=args.samples)
        selected = [r for r in rows if r["cell"] == cell["name"] and not r["warmup"]]
        for field in ["worker_ns", "process_wall_ns", "peak_working_set_bytes"]:
            item[field] = {route: dict(median=statistics.median(r[field] for r in selected if r["route"] == route),
                                      min=min(r[field] for r in selected if r["route"] == route),
                                      max=max(r[field] for r in selected if r["route"] == route)) for route in ["baseline", "candidate"]}
            item[field]["paired_median_candidate_over_baseline"] = statistics.median(
                next(r[field] for r in selected if r["repeat"] == i and r["route"] == "candidate") /
                next(r[field] for r in selected if r["repeat"] == i and r["route"] == "baseline") for i in range(args.samples))
        summary.append(item)
    (output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (output / "manifest.json").write_text(json.dumps({str(p.relative_to(output)): digest(p) for p in output.rglob("*") if p.is_file()}, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(dict(cells=len(cells), measured_processes=len(cells)*args.samples*2, output=str(output))))


if __name__ == "__main__":
    main()
