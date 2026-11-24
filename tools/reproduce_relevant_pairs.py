#!/usr/bin/env python3
"""Bounded Kirchner--Vittorietti Algorithms 1-2 reproduction over Q.

Optional NumPy/GUDHI dependencies belong to this external experiment only.
Every legal pair, including internal empty stages, is checked against the
independent restricted-kernel oracle in compare_persistent_laplacians.py.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import platform
import random
import re
import shutil
import subprocess
import statistics
import sys
import time

import compare_persistent_laplacians as reference


PAPER_SHA256 = "81731c460bbeb12a210fd5745fa862919d4f9992e30e23c07e08419c483009cc"


def reduce_boundary(matrix, columns):
    """Algorithm 1: exact signed columns, no clearing or finite-field pivots."""
    reduced = [[Fraction(int(row[j])) for row in matrix] for j in range(columns)]
    transform = [[Fraction(i == j) for i in range(columns)] for j in range(columns)]
    owners = {}
    for j, column in enumerate(reduced):
        while True:
            low = next((i for i in reversed(range(len(column))) if column[i]), None)
            if low is None or low not in owners:
                break
            i = owners[low]
            factor = column[low] / reduced[i][low]
            column[:] = [a - factor * b for a, b in zip(column, reduced[i])]
            transform[j] = [a - factor * b for a, b in zip(transform[j], transform[i])]
        if low is not None:
            owners[low] = j
    return reduced, transform


def events(rows, columns, cells, stages, reduced):
    index = {value: i + 1 for i, value in enumerate(stages)}
    result = []
    for simplex, column in zip(columns, reduced):
        low = next((i for i in reversed(range(len(column))) if column[i]), None)
        result.append((0 if low is None else index[cells[rows[low]]], index[cells[simplex]]))
    return result


def selected_pairs(labels, count):
    """Algorithm 2's AND rule: its witnesses can be different buckets."""
    return [(s, t) for t in range(count + 1) for s in range(t + 1)
            if any(a <= s and b == t for a, b in labels)
            and any(a == s and b <= t for a, b in labels)]


def predecessor(selected, s, t):
    return max((p for p in selected if p[0] <= s and p[1] <= t),
               key=lambda p: (p[1], p[0]), default=None)


def chain_upper(boundary, transform, labels, s, t, size):
    """Exact Gram--Schmidt, retaining squared norms to avoid irrational roots.

    This equals B Q Q^T B^T. All selected columns are orthogonalized together
    in the FINAL (q+1)-simplex coordinates, including zero-boundary columns.
    """
    orthogonal = []
    for vector, (a, b) in zip(transform, labels):
        if a > s or b > t:
            continue
        vector = vector.copy()
        for other, squared in orthogonal:
            factor = sum(x * y for x, y in zip(vector, other)) / squared
            vector = [x - factor * y for x, y in zip(vector, other)]
        squared = sum(x * x for x in vector)
        if not squared:
            raise ArithmeticError("dependent transformation columns")
        orthogonal.append((vector, squared))
    result = [[Fraction(0) for _ in range(size)] for _ in range(size)]
    for vector, squared in orthogonal:
        image = [sum(Fraction(int(x)) * y for x, y in zip(row, vector))
                 for row in boundary[:size]]
        for i in range(size):
            for j in range(size):
                result[i][j] += image[i] * image[j] / squared
    return result


def controls():
    cases = reference.fixtures()
    cases.append(("joint-events-delayed-vertex", {
        (0,): 0, (1,): 1, (2,): 0, (3,): 0, (4,): 4,
        (0, 1): 2, (2, 3): 3}))
    cases.append(("tetrahedral-shell", {
        s: len(s) - 1 for length in range(1, 4)
        for s in itertools.combinations(range(4), length)}))
    return cases


def paper_style():
    import gudhi
    import numpy as np
    for seed in range(3):
        rng = np.random.default_rng(seed)
        points = rng.random((12, 2))
        tree = gudhi.AlphaComplex(points=points).create_simplex_tree()
        yield f"alpha12-seed{seed}", {tuple(s): v for s, v in tree.get_filtration()}, points.tolist()
        points = rng.random((8, 3))
        tree = gudhi.RipsComplex(points=points, max_edge_length=0.3).create_simplex_tree(max_dimension=4)
        yield f"rips8-radius0.3-seed{seed}", {tuple(s): v for s, v in tree.get_filtration()}, points.tolist()


def serialize(value):
    if isinstance(value, Fraction):
        return str(value)
    raise TypeError(type(value).__name__)


def check_case(name, cells, degree, output, benchmark=False, worker=None):
    import numpy as np
    begin = time.perf_counter()
    stages = sorted(set(cells.values()))
    scales = [float("-inf")] + stages  # internal sentinel only
    rows = reference.ordered(cells, float("inf"), degree)
    columns = reference.ordered(cells, float("inf"), degree + 1)
    boundary = reference.boundary(rows, columns)
    preparation_seconds = time.perf_counter() - begin
    begin = time.perf_counter()
    reduced, transform = reduce_boundary(boundary, len(columns))
    reduction_seconds = time.perf_counter() - begin
    # Check R=B V and invertible upper-triangular V independently of pivot logic.
    for j, column in enumerate(transform):
        assert column[j] == 1 and all(v == 0 for v in column[j + 1:])
        assert [sum(Fraction(int(x)) * y for x, y in zip(row, column))
                for row in boundary] == reduced[j]
    begin = time.perf_counter()
    labels = events(rows, columns, cells, stages, reduced)
    selected = selected_pairs(labels, len(stages))
    sizes = [sum(cells[row] <= value for row in rows) for value in scales]
    scheduling_seconds = time.perf_counter() - begin
    solved = {pair: chain_upper(boundary, transform, labels, *pair, sizes[pair[0]])
              for pair in selected}
    records = []
    gudhi_checks = 0
    import gudhi
    tree = gudhi.SimplexTree()
    for simplex, value in cells.items():
        tree.insert(simplex, filtration=value)
    assert tree.num_simplices() == len(cells)
    tree.compute_persistence(homology_coeff_field=11, persistence_dim_max=True)
    intervals = tree.persistence_intervals_in_dimension(degree)
    for t in range(len(scales)):
        for s in range(t + 1):
            n = sizes[s]
            representative = predecessor(selected, s, t)
            restored = [[Fraction(0) for _ in range(n)] for _ in range(n)]
            if representative is not None:
                old = solved[representative]
                for i, row in enumerate(old):
                    restored[i][:len(row)] = row
            end_rows = reference.ordered(cells, scales[t], degree)
            end_columns = reference.ordered(cells, scales[t], degree + 1)
            pair_boundary = reference.boundary(end_rows, end_columns)
            lower = reference.boundary(reference.ordered(cells, scales[s], degree - 1), rows[:n])
            # Independent row-RREF/nullspace oracle, not Algorithm 1 or V.
            expected = reference.exact_projection(pair_boundary[:n], pair_boundary[n:],
                                                   np.zeros((0, n), dtype=int))
            assert restored == expected, (name, degree, s, t, representative, restored, expected)
            down = lower.T @ lower
            full = np.array([[float(x) for x in row] for row in restored]).reshape(n, n) + down
            # GUDHI corroborates topology only; it is not a real-spectrum oracle.
            beta = n - len(reference.rref(lower.tolist(), n)[1])
            beta -= len(reference.rref(pair_boundary.tolist(), len(end_columns))[1])
            beta += len(reference.rref(pair_boundary[n:].tolist(), len(end_columns))[1])
            nullity = int(np.count_nonzero(np.abs(np.linalg.eigvalsh(full)) < 1e-8))
            assert nullity == beta, (name, degree, s, t, nullity, beta)
            field_beta = sum(b <= scales[s] and d > scales[t] for b, d in intervals)
            assert field_beta == beta, (name, degree, s, t, field_beta, beta)
            gudhi_checks += 1
            exact_full = [[restored[i][j] + int(down[i, j]) for j in range(n)] for i in range(n)]
            records.append({"pair": [s, t], "size": n, "representative": representative,
                            "embedding": list(range(0 if representative is None else sizes[representative[0]])),
                            "upper": restored, "full": exact_full, "real_nullity": beta})
    payload = {"name": name, "degree": degree, "stages": stages, "row_basis": rows,
               "column_basis": columns, "R_columns": reduced, "V_columns": transform,
               "events": labels, "selected": selected, "pairs": records}
    (output / f"{name}-q{degree}.json").write_text(
        json.dumps(payload, default=serialize, separators=(",", ":")) + "\n", encoding="utf-8")
    rust = compare_rust(worker, name, cells, degree, stages, records, output) if worker else None
    timing = time_routes(boundary, transform, labels, sizes, selected,
                         preparation_seconds, reduction_seconds, scheduling_seconds) if benchmark else None
    field_control = None
    if name == "rp2-field-counterexample":
        tree.compute_persistence(homology_coeff_field=2, persistence_dim_max=True)
        field2 = sum(b <= scales[-1] and d > scales[-1]
                     for b, d in tree.persistence_intervals_in_dimension(degree))
        real = records[-1]["real_nullity"]
        if degree in (1, 2):
            assert field2 == 1 and real == 0
        field_control = {"final_field2_betti": int(field2), "final_real_betti": real}
    return {"name": name, "degree": degree, "possible_pairs": len(stages) * (len(stages) + 1) // 2,
            "requested_pairs": len(stages) * (len(stages) + 1) // 2,
            "checked_pairs_including_sentinel": len(records),
            "selected_pairs": sum(s > 0 for s, _ in selected),
            "actually_solved_pairs_including_sentinel": len(solved), "gudhi_field11_checks": gudhi_checks,
            "exact_matrix_equalities": len(records), "max_matrix_error": 0.0,
            "preparation_seconds": preparation_seconds, "reduction_seconds": reduction_seconds,
            "scheduling_seconds": scheduling_seconds, "timings": timing, "rust": rust,
            "field_control": field_control}


def compare_rust(worker, name, cells, degree, stages, records, output):
    import numpy as np
    data = reference.protocol(cells, degree, stages, stages[-1] if stages else 0.0)
    run = subprocess.run([str(worker), "--relevant"], input=data, text=True, capture_output=True)
    (output / f"{name}-q{degree}-rust-stdout.txt").write_text(run.stdout, encoding="utf-8")
    (output / f"{name}-q{degree}-rust-stderr.txt").write_text(run.stderr, encoding="utf-8")
    run.check_returncode()
    lines = iter(run.stdout.splitlines())
    batch = next(lines).split()
    expected_records = [r for r in records if r["pair"][0] > 0]
    assert int(batch[1]) == len(expected_records)
    max_error = 0.0
    for expected in expected_records:
        header = next(lines)
        fields = [int(x) for x in re.findall(r"\d+", header)]
        s, t, n = fields[:3]
        assert [s + 1, t + 1] == expected["pair"] and n == expected["size"]
        p = expected["representative"]
        expected_p = None if p is None or p[0] == 0 else [p[0] - 1, p[1] - 1]
        assert (fields[3:] or None) == expected_p, (header, expected_p)
        actual = np.array([float(x) for x in next(lines).split()[1:]]).reshape(n, n)
        reference_matrix = np.array([[float(Fraction(x)) for x in row] for row in expected["full"]]).reshape(n, n)
        error = np.linalg.norm(actual - reference_matrix)
        assert error <= 1e-10 * max(1.0, np.linalg.norm(reference_matrix)), (name, degree, s, t, error)
        max_error = max(max_error, float(error))
        spectrum = np.array([float(x) for x in next(lines).split()[1:]])
        np.testing.assert_allclose(spectrum, np.linalg.eigvalsh(reference_matrix), rtol=1e-10, atol=1e-10)
        residual = float(next(lines).split()[1])
        assert residual <= 1e-9
        assert int(next(lines).split()[1]) == (0 if expected_p is None else len(expected["embedding"]))
        basis = json.loads(next(lines).split(" ", 1)[1])
        assert basis == [list(simplex) for simplex in reference.ordered(cells, stages[s], degree)]
        context = next(lines).split()
        assert [int(context[1]), float(context[2]), float(context[3])] == [degree, stages[s], stages[t]]
    assert next(lines, None) is None
    return {"pairs": len(expected_records), "actually_solved_pairs": int(batch[2]),
            "construction_seconds": float(batch[3]), "max_matrix_frobenius_error": max_error}


def time_routes(boundary, transform, labels, sizes, selected,
                preparation_seconds, reduction_seconds, scheduling_seconds):
    """Exploratory serial phase timing with identical full upper-matrix output.

    These are single-run diagnostics, not formal performance acceptance.
    Filtration, exact oracle, serialization and process startup are excluded.
    """
    import numpy as np
    stages = len(sizes) - 1
    pairs = [(s, t) for t in range(1, stages + 1) for s in range(1, t + 1)]
    end_counts = [sum(b <= t for _, b in labels) for t in range(stages + 1)]
    results = []
    for route in ["exhaustive-schur", "selective-schur", "new-chains"]:
        start = time.perf_counter()
        evaluated = pairs if route == "exhaustive-schur" else selected
        matrices = {}
        qr_seconds = 0.0
        for s, t in evaluated:
            n = sizes[s]
            if route == "new-chains":
                begin = time.perf_counter()
                indices = [j for j, (a, b) in enumerate(labels) if a <= s and b <= t]
                v = np.array([[float(x) for x in transform[j]] for j in indices]).reshape(len(indices), len(labels)).T
                q, _ = np.linalg.qr(v, mode="reduced")
                qr_seconds += time.perf_counter() - begin
                factor = boundary[:n] @ q
                matrix = factor @ factor.T
            else:
                b = boundary[:sizes[t], :end_counts[t]].astype(float)
                gram = b @ b.T
                matrix = gram[:n, :n].copy()
                if sizes[t] > n:
                    e, c = gram[:n, n:], gram[n:, n:]
                    matrix -= e @ np.linalg.lstsq(c, e.T, rcond=1e-12)[0]
            matrices[s, t] = matrix
        construction_seconds = time.perf_counter() - start
        begin = time.perf_counter()
        restored = []
        for s, t in pairs:
            p = (s, t) if route == "exhaustive-schur" else predecessor(selected, s, t)
            matrix = np.zeros((sizes[s], sizes[s]))
            if p is not None:
                old = matrices[p]
                matrix[:old.shape[0], :old.shape[1]] = old
            restored.append(matrix)
        reconstruction_seconds = time.perf_counter() - begin
        begin = time.perf_counter()
        spectra = [np.linalg.eigvalsh(m) for m in restored]
        spectrum_seconds = time.perf_counter() - begin
        preparation = preparation_seconds
        if route != "exhaustive-schur":
            preparation += reduction_seconds + scheduling_seconds
        results.append({"route": route, "solves": len(evaluated), "preparation_seconds": preparation,
                        "construction_seconds": construction_seconds, "orthogonalization_seconds": qr_seconds,
                        "reconstruction_seconds": reconstruction_seconds, "spectrum_seconds": spectrum_seconds,
                        "total_seconds": preparation + construction_seconds + reconstruction_seconds + spectrum_seconds,
                        "checksum": float(sum(np.sum(v) for v in spectra)),
                        "returned_entries": sum(m.size for m in restored)})
        if route == "exhaustive-schur":
            baseline = restored
        else:
            for expected, actual in zip(baseline, restored):
                error = np.linalg.norm(expected - actual)
                assert error <= 1e-10 * max(1.0, np.linalg.norm(expected)), (route, error)
    return results


def benchmark_workers(baseline, candidate, cases, output, baseline_revision):
    """Ten serial fresh-process paired samples plus one retained warmup.

    Full serialized matrices/spectra/metadata must match, not just a checksum.
    Kernel construction excludes spectra; process wall/peak include all output.
    """
    if sys.platform != "win32" or not re.fullmatch(r"[0-9a-f]{40}", baseline_revision):
        raise ValueError("resource snapshot requires Windows and a full baseline commit SHA")
    directory = output / "benchmark"
    directory.mkdir()
    root = Path(__file__).resolve().parents[1]
    changed = subprocess.check_output(["git", "diff", "HEAD", "--name-only"], cwd=root, text=True)
    if any(name.endswith((".rs", ".py", ".toml", ".lock")) for name in changed.splitlines()):
        raise ValueError("commit implementation and harness before resource measurements")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root, text=True).strip()
    manifest = {"evidence_class": "resource_snapshot", "candidate_revision": revision,
                "baseline_revision": baseline_revision, "samples": 10, "warmups": 1,
                "order_seed": 146, "platform": platform.platform(), "python": sys.version,
                "memory": "absolute process peak working set, including input/runtime/retained output",
                "controls": "serial processes; affinity/frequency/host load uncontrolled",
                "source_fingerprints": {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
                                        for p in root.glob("src/**/*.rs")},
                "harness_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                "executables": {label: hashlib.sha256(path.read_bytes()).hexdigest()
                                for label, path in [("baseline", baseline), ("candidate", candidate)]}}
    (directory / "environment.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    rng = random.Random(146)
    samples = []
    summaries = []
    workloads = []
    for name, cells, points in cases:
        if name in {"path", "joint-events-delayed-vertex", "tetrahedral-shell", "tetrahedron",
                    "alpha12-seed0", "rips8-radius0.3-seed0"}:
            degrees = [0, 1] if name == "alpha12-seed0" else [2 if name == "tetrahedron" else 1 if name == "tetrahedral-shell" else 0]
            for degree in degrees:
                workloads.append((name, cells, degree, sorted(set(cells.values()))))
            if name == "alpha12-seed0":
                workloads.append((name + "-single-pair", cells, 1, [max(cells.values())]))
    for name, cells, degree, grid in workloads:
        data = reference.protocol(cells, degree, grid, grid[-1])
        (directory / f"{name}-q{degree}-input.txt").write_text(data, encoding="utf-8")
        expected = None
        first = rng.randrange(2)
        for repeat in range(-1, 10):
            labels = ["baseline", "candidate"] if (repeat + first) % 2 else ["candidate", "baseline"]
            for label in labels:
                executable = baseline if label == "baseline" else candidate
                begin = time.perf_counter()
                process = subprocess.Popen([str(executable), "--relevant"], stdin=subprocess.PIPE,
                                           stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                stdout, stderr = process.communicate(data)
                wall = time.perf_counter() - begin
                peak = reference.peak_memory(process)
                stem = directory / f"{name}-q{degree}-{repeat}-{label}"
                stem.with_suffix(".stdout.txt").write_text(stdout, encoding="utf-8")
                stem.with_suffix(".stderr.txt").write_text(stderr, encoding="utf-8")
                if process.returncode:
                    raise RuntimeError(f"benchmark worker failed: {stem}, exit {process.returncode}")
                if peak is None:
                    raise RuntimeError("process peak working set unavailable")
                batch, *body = stdout.splitlines()
                fields = batch.split()
                actual = (fields[:3], body)
                if expected is None:
                    expected = actual
                assert actual == expected, (name, degree, repeat, label)
                samples.append({"name": name, "degree": degree, "repeat": repeat, "route": label,
                                "warmup": repeat == -1, "construction_seconds": float(fields[3]),
                                "process_seconds": wall, "peak_working_set_bytes": peak,
                                "input_sha256": hashlib.sha256(data.encode()).hexdigest(),
                                "requested_pairs": int(fields[1]), "actually_solved_pairs": int(fields[2])})
                (directory / "samples.json").write_text(json.dumps(samples, indent=2) + "\n", encoding="utf-8")
        rows = [r for r in samples if r["name"] == name and r["degree"] == degree and not r["warmup"]]
        summary = {"name": name, "degree": degree, "requested_pairs": rows[0]["requested_pairs"],
                   "actually_solved_pairs": rows[0]["actually_solved_pairs"]}
        for field in ["construction_seconds", "process_seconds", "peak_working_set_bytes"]:
            summary[field] = {label: {"median": statistics.median(r[field] for r in rows if r["route"] == label),
                                     "min": min(r[field] for r in rows if r["route"] == label),
                                     "max": max(r[field] for r in rows if r["route"] == label)}
                              for label in ["baseline", "candidate"]}
            summary[field]["paired_median_ratio"] = statistics.median(
                next(r[field] for r in rows if r["repeat"] == i and r["route"] == "candidate") /
                next(r[field] for r in rows if r["repeat"] == i and r["route"] == "baseline") for i in range(10))
        summaries.append(summary)
        (directory / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n", encoding="utf-8")
        print(name, degree, "paired resource snapshot passed", flush=True)


def main():
    if not __debug__:
        raise RuntimeError("correctness reproduction requires Python without -O")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--paper-style", action="store_true")
    parser.add_argument("--rust-worker", type=Path)
    parser.add_argument("--benchmark-baseline", type=Path)
    parser.add_argument("--baseline-revision")
    args = parser.parse_args()
    if args.benchmark_baseline and (not args.rust_worker or not args.baseline_revision):
        parser.error("benchmark requires --rust-worker and --baseline-revision")
    if hashlib.sha256(args.paper.read_bytes()).hexdigest() != PAPER_SHA256:
        raise ValueError("paper does not match the frozen arXiv v1 bytes")
    args.output.mkdir(parents=True, exist_ok=False)
    import gudhi
    import numpy as np
    manifest = {"paper_url": "https://arxiv.org/pdf/2611.09610v1", "paper_sha256": PAPER_SHA256,
                "author_code": "not verified; independent mathematical transcription",
                "python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
                "gudhi": gudhi.__version__, "command": sys.argv, "scripts": {}}
    for source in [Path(__file__), Path(reference.__file__)]:
        shutil.copy2(source, args.output / source.name)
        manifest["scripts"][source.name] = hashlib.sha256(source.read_bytes()).hexdigest()
    worker = None
    if args.rust_worker:
        worker = (args.output / args.rust_worker.name).resolve()
        shutil.copy2(args.rust_worker, worker)
        manifest["rust_worker_sha256"] = hashlib.sha256(worker.read_bytes()).hexdigest()
    (args.output / "environment.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    cases = [(name, cells, None) for name, cells in controls()]
    if args.paper_style:
        cases.extend(paper_style())
    summaries = []
    for name, cells, points in cases:
        fixture = json.dumps({"cells": [[list(s), v] for s, v in cells.items()], "points": points}, separators=(",", ":"))
        (args.output / f"{name}-input.json").write_text(fixture + "\n", encoding="utf-8")
        for degree in range(4):
            result = check_case(name, cells, degree, args.output, benchmark=points is not None, worker=worker)
            result["input_sha256"] = hashlib.sha256((fixture + "\n").encode()).hexdigest()
            summaries.append(result)
            (args.output / "summary.json").write_text(json.dumps(summaries, indent=2) + "\n", encoding="utf-8")
            print(name, degree, result["exact_matrix_equalities"], "exact pairs passed", flush=True)
    print("PASS", sum(r["exact_matrix_equalities"] for r in summaries), "exact pair matrices", flush=True)
    if args.benchmark_baseline:
        baseline = (args.output / ("baseline-" + args.benchmark_baseline.name)).resolve()
        shutil.copy2(args.benchmark_baseline, baseline)
        benchmark_workers(baseline, worker, cases, args.output, args.baseline_revision)


if __name__ == "__main__":
    main()
