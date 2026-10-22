"""Check bounded sheaf and application-feedback examples independently.

GUDHI checks ordinary topology and recomputed persistence, never general sheaf
restrictions. No new Python package is needed unless --gudhi is requested.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import sys
import time


def check_sheaf(records):
    inputs = {
        "nonconstant": ([1, 1], 1, [[1], [2]]),
        "constant": ([1, 1], 1, [[1], [1]]),
        "rectangular": ([2, 1], 1, [[1, 0], [2]]),
        "zero_edge": ([1, 1], 0, [[], []]),
        "zero_vertices": ([0, 0], 1, [[], []]),
        "empty": ([0, 0], 0, [[], []]),
    }
    expected = {
        "nonconstant": ([1., -2., -2., 4.], [0., 5.]),
        "constant": ([1., -1., -1., 1.], [0., 2.]),
        "rectangular": ([1., 0., -2., 0., 0., 0., -2., 0., 4.], [0., 0., 5.]),
        "zero_edge": ([0.] * 4, [0., 0.]),
        "zero_vertices": ([], []), "empty": ([], []),
    }
    if len(records) != len(expected) or {r["name"] for r in records} != set(expected):
        raise ValueError("missing or repeated sheaf fixtures")
    for record in records:
        matrix, values = expected[record["name"]]
        if (record["vertices"], record["edge"], record["restrictions"]) != inputs[record["name"]]:
            raise ValueError("sheaf fixture identity differs")
        if record["kind"] != "sheaf" or record["matrix"] != matrix:
            raise ValueError("sheaf matrix differs from hand expectation")
        if len(values) != len(record["values"]) or any(
                not abs(a-b) <= 1e-10 for a, b in zip(record["values"], values)):
            raise ValueError("sheaf eigenvalues differ")
        n = len(values)
        if len(record["x"]) != n or len(record["action"]) != n:
            raise ValueError("sheaf action has the wrong shape")
        action = [sum(matrix[i*n+j] * record["x"][j] for j in range(n))
                  for i in range(n)]
        if record["action"] != action or len(record["residuals"]) != n or any(
                not 0 <= r < 1e-10 for r in record["residuals"]):
            raise ValueError("sheaf action or residual failed")
    return {"cases": len(records), "state": "Passed", "oracle": "hand matrices and spectra"}


def check_feedback(records, gudhi=None):
    updates = [r for r in records if r["kind"] == "feedback"]
    ties = [r for r in records if r["kind"] == "tie"]
    if len(records) != 12 or len(updates) != 9 or len(ties) != 3:
        raise ValueError("feedback output is incomplete")
    if len({r["snapshot"] for r in updates}) != 9:
        raise ValueError("feedback reused a snapshot key")
    for step, record in enumerate(updates):
        theta = 1 + 2. ** -step
        if (record["iteration"], record["theta"], record["mass"], record["finite_count"],
            record["loss"], record["gradient"], record["difference"], record["h"]) != (
                step, theta, theta, 1, 0.5 * 4. ** -step, theta-1, theta-1, 1/1024):
            raise ValueError("feedback differs from the exact contraction")
    h = 1/1024
    if [r["theta"] for r in ties] != [-h, 0., h]:
        raise ValueError("missing one-sided tie controls")
    for record in ties:
        mass = max(0., record["theta"])
        if (record["mass"], record["finite_count"], record["loss"]) != (
                mass, int(mass > 0), 0.5 * (mass-1)**2):
            raise ValueError("tie persistence or loss differs")
    if gudhi is not None:
        for record in records:
            tree = gudhi.SimplexTree()
            tree.insert([0], filtration=0.)
            tree.insert([1], filtration=0.)
            tree.insert([0, 1], filtration=max(0., record["theta"]))
            tree.compute_persistence(homology_coeff_field=2, min_persistence=0,
                                     persistence_dim_max=True)
            bars = tree.persistence_intervals_in_dimension(0)
            finite = [(b, d) for b, d in bars if d != float("inf")]
            if len(finite) != record["finite_count"] or sum(d-b for b, d in finite) != record["mass"]:
                raise ValueError("GUDHI feedback persistence differs")
        tree = gudhi.SimplexTree()
        tree.insert([0, 1], filtration=0.)
        tree.compute_persistence(homology_coeff_field=2, persistence_dim_max=True)
        if tree.betti_numbers() != [1, 0]:
            raise ValueError("constant-sheaf ordinary topology control differs")
    return {"updates": 9, "tie_cases": 3, "state": "Passed",
            "gudhi": gudhi.__version__ if gudhi is not None else "NotRun"}


def run(worker, output):
    started = time.perf_counter_ns()
    result = subprocess.run([str(worker.resolve())], capture_output=True, text=True,
                            check=False, timeout=30)
    record = {"command": [str(worker.resolve())], "exit_code": result.returncode,
              "stdout": result.stdout, "stderr": result.stderr,
              "process_elapsed_ns": time.perf_counter_ns()-started,
              "worker_sha256": hashlib.sha256(worker.read_bytes()).hexdigest()}
    output.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    if result.returncode:
        raise ValueError(f"worker failed: {worker.name}")
    return [json.loads(line) for line in result.stdout.splitlines()]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sheaf", type=Path, required=True)
    parser.add_argument("--feedback", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gudhi", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    root = Path(__file__).resolve().parents[1]
    source = hashlib.sha256()
    paths = list((root/"src").rglob("*.rs")) + [Path(__file__),
        root/"examples/sheaf_probe.rs", root/"examples/filtration_feedback.rs",
        root/"tests/homology_operators.rs"]
    for path in sorted(paths):
        source.update(path.relative_to(root).as_posix().encode())
        source.update(path.read_bytes())
    gudhi = None
    if args.gudhi:
        import gudhi
    summary = {"candidate_sha": subprocess.check_output(
                   ["git", "-C", str(root), "rev-parse", "HEAD"], text=True).strip(),
               "source_sha256": source.hexdigest(), "python": sys.version,
               "platform": platform.platform(),
               "sheaf": check_sheaf(run(args.sheaf, args.output/"sheaf.json")),
               "feedback": check_feedback(run(args.feedback, args.output/"feedback.json"), gudhi),
               "limits": ["bounded examples; no general sheaf API or automatic gradient API",
                          "GUDHI does not verify nonconstant restrictions",
                          "process duration includes startup; no comparative performance claim"]}
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
