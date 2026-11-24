#!/usr/bin/env python3
"""Bounded ordinary q=1 HeCS reproduction; never a persistent solver backend.

NumPy/SciPy/GUDHI are optional experiment dependencies, imported only on use.
The Julia adapter includes unchanged, hash-checked author files. Rust parity is
opt-in and follows a successful paper/oracle run on the same frozen fixtures.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import time

AUTHOR_COMMIT = "cd4d23e5d53900ece219fc0cdafc5c5e1eb81bf3"
SEED = 150
MATRIX_TOL = 1e-10
RESIDUAL_TOL = 1e-9
ROOT = Path(__file__).resolve().parents[1]
PAPER_SHA256 = "20f7fa2ba91d308a97b1ad98f9d688a4a6ae00ce630bd457dd9e71f1e6c6d537"
AUTHOR_HASHES = {
    "README.md": "5d57de3bf9b4ecdac2abbe523506c90b0c7a36902402aae9aae9bfa06da83b7a",
    "src/collapse.jl": "4be26d62090ecce1c7fc1ff81033791931ee9ebf81b65af686cda52d7fb550ca",
    "src/cgls.jl": "e5c5d6b44065ea3d395fd6798520dd8b29ea2c84764b98ca4381d11f70aa2c61",
    "src/simplicialComplex.jl": "4e17b4872210e4f24a7d6c59973f38a41b9608fc0fe9057db455b6e2cb938e60",
    "src/utils.jl": "c8250f73c8794459f5673e6960742992fadfe98f391258ac228b4635d4c8ba16",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def boundary(edges, faces):
    """Oriented simplicial boundary, independent of the collapse routine."""
    import numpy as np
    result = np.zeros((len(edges), len(faces)), dtype=np.int64)
    index = {tuple(edge): i for i, edge in enumerate(edges)}
    for j, face in enumerate(faces):
        for removed in range(3):
            edge = tuple(face[:removed] + face[removed + 1:])
            result[index[edge], j] = (-1) ** removed
    return result


def rational_rank(matrix):
    rows = [[Fraction(int(x)) for x in row] for row in matrix]
    width = len(rows[0]) if rows else 0
    rank = 0
    for col in range(width):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        value = rows[rank][col]
        rows[rank] = [x / value for x in rows[rank]]
        for i in range(rank + 1, len(rows)):
            value = rows[i][col]
            rows[i] = [x - value * y for x, y in zip(rows[i], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def collapse(edge_faces, selected):
    """Algorithm 5.1, smallest original edge ID resolves the permitted pop tie."""
    remaining = set(selected)
    pairs = []
    while remaining:
        free = next(((e, sorted(remaining.intersection(adjacent)))
                     for e, adjacent in enumerate(edge_faces)
                     if len(remaining.intersection(adjacent)) == 1), None)
        if free is None:
            break
        e, adjacent = free
        t = adjacent[0]
        pairs.append([e, t])
        remaining.remove(t)
    return pairs, sorted(remaining)


def heavy(edge_faces, weights, peel=False):
    """Algorithm 6.1; optionally the Section 6.2 preliminary 2-core peeling."""
    initial, core = collapse(edge_faces, range(len(weights)))
    selected = [t for _, t in initial] if peel else []
    candidates = core if peel else list(range(len(weights)))
    decisions = []
    for triangle in sorted(candidates, key=lambda t: (-weights[t], t)):
        pairs, remaining = collapse(edge_faces, selected + [triangle])
        accept = not remaining
        decisions.append({"triangle": triangle, "accepted": accept,
                          "pairs": pairs, "remaining": remaining})
        if accept:
            selected.append(triangle)
    pairs, remaining = collapse(edge_faces, selected)
    assert not remaining
    return {"selected": selected, "pairs": pairs, "decisions": decisions,
            "initial_pairs": initial, "core": core}


def replay(edges, faces, selected, pairs):
    """Check a certificate directly from vertex sets, without collapse adjacency."""
    remaining = set(selected)
    used_edges = set()
    for e, t in pairs:
        assert e not in used_edges and t in remaining
        adjacent = [j for j in remaining if set(edges[e]) <= set(faces[j])]
        assert adjacent == [t], (e, t, adjacent)
        used_edges.add(e)
        remaining.remove(t)
    return not remaining


def residual(matrix, x, rhs):
    import numpy as np
    return float(np.linalg.norm(matrix @ x - rhs) / max(1., np.linalg.norm(rhs)))


def cgls(matrix, rhs, original, original_rhs, lift, tol=RESIDUAL_TOL, maxit=4000):
    """Zero-start CGLS; stop on a fresh residual in the ORIGINAL equation.

    Unlike the author demonstration, zero RHS and stationary inconsistent
    least squares have explicit outcomes. No ridge or dense fallback is used.
    """
    import numpy as np
    x = np.zeros(matrix.shape[1])
    r = rhs.copy()
    s = matrix.T @ r
    p = s.copy()
    gamma = float(s @ s)
    trace = []
    for iteration in range(maxit + 1):
        fresh = residual(original, lift @ x, original_rhs)
        trace.append({"iteration": iteration, "x": x.tolist(),
                      "residual": float(np.linalg.norm(r)),
                      "normal_residual": math.sqrt(gamma),
                      "original_residual": fresh})
        if not np.isfinite(fresh) or not np.isfinite(gamma):
            return {"status": "nonfinite", "x": x.tolist(), "trace": trace}
        if fresh <= tol:
            return {"status": "converged", "x": x.tolist(), "trace": trace}
        if iteration == maxit:
            return {"status": "iteration_limit", "x": x.tolist(), "trace": trace}
        q = matrix @ p
        denominator = float(q @ q)
        if gamma == 0 or denominator == 0:
            return {"status": "breakdown", "x": x.tolist(), "trace": trace}
        alpha = gamma / denominator
        x += alpha * p
        r -= alpha * q
        s = matrix.T @ r
        next_gamma = float(s @ s)
        beta = next_gamma / gamma
        trace[-1].update(alpha=alpha, beta=beta, gamma=gamma)
        p = s + beta * p
        gamma = next_gamma
    raise AssertionError("unreachable")


def fixtures():
    import numpy as np
    from scipy.spatial import Delaunay
    cases = []

    def add(name, n, edges, faces, weights, provenance="project control", points=None):
        b = boundary(edges, faces)
        u = (b * np.asarray(weights)) @ b.T
        x = np.array([((7 * i + 3) % 11 - 5) / 5 for i in range(len(edges))])
        cases.append({"name": name, "vertices": list(range(n)), "edges": edges,
                      "faces": faces, "weights": weights, "rhs": (u @ x).tolist(),
                      "points": points, "provenance": provenance})

    add("empty", 0, [], [], [])
    add("isolated-skeleton", 4, [[0, 1], [1, 2]], [], [])
    add("triangle", 3, [[0, 1], [0, 2], [1, 2]], [[0, 1, 2]], [1])
    add("triangle-zero-rhs", 3, [[0, 1], [0, 2], [1, 2]], [[0, 1, 2]], [1])
    cases[-1]["rhs"] = [0., 0., 0.]
    edges4 = [list(e) for e in itertools.combinations(range(4), 2)]
    faces4 = [list(t) for t in itertools.combinations(range(4), 3)]
    add("shell-unit", 4, edges4, faces4, [1, 1, 1, 1])
    add("shell-weighted", 4, edges4, faces4, [1, 2, 4, 8])
    add("disk-with-extra-edge", 5, edges4 + [[3, 4]], faces4[1:], [2, 4, 8])
    # A torsion example: real rank is not inferred from F2 topology.
    rp2 = [[0, 1, 2], [0, 1, 3], [0, 2, 4], [0, 3, 5], [0, 4, 5],
           [1, 2, 5], [1, 3, 4], [1, 4, 5], [2, 3, 4], [2, 3, 5]]
    add("rp2-rank-loss", 6, [list(e) for e in itertools.combinations(range(6), 2)],
        rp2, list(range(1, 11)))
    rng = np.random.default_rng(SEED)
    for n in (8, 12, 20):
        points = np.vstack(([[0, 0], [1, 0], [1, 1], [0, 1]], rng.random((n - 4, 2))))
        triangulation = Delaunay(points)
        edge_set = {tuple(sorted(map(int, e))) for t in triangulation.simplices
                    for e in itertools.combinations(t, 2)}
        outer = {(0, 1), (1, 2), (2, 3), (0, 3)}
        interior = sorted(edge_set - outer)
        removed = rng.choice(len(interior), size=2, replace=False)
        for i in removed:
            edge_set.remove(interior[i])
        absent = sorted(set(itertools.combinations(range(n), 2)) - edge_set)
        order = rng.permutation(len(absent))
        for extra in (0, 2, 5):
            enriched = edge_set | {absent[i] for i in order[:extra]}
            edges = sorted(enriched)
            faces = [list(t) for t in itertools.combinations(range(n), 3)
                     if all(e in enriched for e in itertools.combinations(t, 2))]
            edge_weights = {e: float(w) for e, w in zip(edges, abs(rng.standard_normal(len(edges))))}
            weights = [min(edge_weights[e] for e in itertools.combinations(t, 2)) for t in faces]
            add(f"enriched-{n}-{extra}", n, [list(e) for e in edges], faces, weights,
                "regenerated Section 7.3 style; SciPy Delaunay, not author data", points.tolist())
    return cases


def protocol(case, peel=False):
    lines = [f"{len(case['vertices'])} {len(case['edges'])} {len(case['faces'])} {int(peel)}"]
    lines += [" ".join(map(str, e)) for e in case["edges"]]
    lines += [" ".join(map(str, t + [w])) for t, w in zip(case["faces"], case["weights"])]
    lines.append(" ".join(map(str, case["rhs"])))
    return "\n".join(lines) + "\n"


def python_case(case, peel):
    import numpy as np
    import gudhi
    edges, faces, weights = case["edges"], case["faces"], case["weights"]
    b = boundary(edges, faces)
    rhs = np.asarray(case["rhs"])
    m, t = b.shape
    u = (b * np.asarray(weights)) @ b.T
    adjacent = [list(map(int, np.flatnonzero(row))) for row in b]
    start = time.perf_counter_ns()
    selection = heavy(adjacent, weights, peel)
    selection_ns = time.perf_counter_ns() - start
    selected = selection["selected"]
    pairs = selection["pairs"]
    assert replay(edges, faces, selected, pairs)
    accepted = [t for _, t in selection["initial_pairs"]] if peel else []
    for d in selection["decisions"]:
        trial = accepted + [d["triangle"]]
        assert replay(edges, faces, trial, d["pairs"]) == d["accepted"]
        assert set(trial) == set(d["remaining"]) | {t for _, t in d["pairs"]}
        if d["accepted"]:
            accepted.append(d["triangle"])
    rank = rational_rank(b)
    selected_rank = rational_rank(b[:, selected])
    assert selected_rank == len(selected)
    if case["name"] == "rp2-rank-loss":
        assert all(np.count_nonzero(row) == 2 for row in b)
        assert rank == 10 and selected_rank == 9
    if case["name"] == "shell-weighted":
        assert sorted(selected) == [1, 2, 3] and not selection["initial_pairs"]
    if case["name"] == "shell-unit":
        s = np.array([1, -1, 1, -1])
        assert np.array_equal(b.T @ b, 4 * np.eye(4) - np.outer(s, s))
        assert rank == 3
    # An exact rational weighted factor identity, before square roots/roundoff.
    exact_u = [[sum(Fraction(float(weights[j])) * int(b[i, j]) * int(b[k, j])
                    for j in selected) for k in range(m)] for i in range(m)]
    order = [e for e, _ in pairs] + [i for i in range(m) if i not in {e for e, _ in pairs}]
    columns = [t for _, t in pairs]
    exact_factor = [[sum(Fraction(float(weights[j])) * int(b[i, j]) * int(b[k, j])
                        for j in columns) for k in order] for i in order]
    assert exact_factor == [[exact_u[i][k] for k in order] for i in order]
    c = b[np.ix_(order, columns)] * np.sqrt(np.asarray(weights)[columns])
    r = len(columns)
    assert np.count_nonzero(np.triu(c[:r], 1)) == 0
    reference = np.array([[float(exact_u[i][k]) for k in order] for i in order]).reshape(m, m)
    factor_error = float(np.linalg.norm(c @ c.T - reference))
    assert factor_error <= MATRIX_TOL * max(1., float(np.linalg.norm(reference)))
    eig = np.linalg.eigvalsh(u).tolist()
    if case["name"] == "shell-unit":
        assert np.allclose(eig[-3:], [4, 4, 4], atol=MATRIX_TOL, rtol=0)
    st = gudhi.SimplexTree()
    for v in case["vertices"]:
        st.insert([v])
    for simplex in edges + faces:
        st.insert(simplex)
    st.compute_persistence(homology_coeff_field=2, persistence_dim_max=True)
    betti = st.betti_numbers()
    # GUDHI corroborates incidence topology only. RP2 deliberately differs over F2.
    assert st.num_simplices() == len(case["vertices"]) + m + t
    if case["name"].startswith("shell"):
        assert betti == [1, 0, 1]
    if case["name"] == "rp2-rank-loss":
        assert betti == [1, 1, 1]  # torsion; not the real Betti numbers
    oracle_x = np.linalg.lstsq(u, rhs, rcond=1e-12)[0] if m else np.zeros(0)
    solves = {}
    for name in ("identity", "diagonal"):
        start = time.perf_counter_ns()
        scale = np.ones(m)
        if name == "diagonal":
            diagonal = np.diag(u)
            scale[diagonal > 0] = 1 / np.sqrt(diagonal[diagonal > 0])
        lift = np.diag(scale)
        a = (u * scale[:, None]) * scale[None, :]
        h = scale * rhs
        setup_ns = time.perf_counter_ns() - start
        start = time.perf_counter_ns()
        solve = cgls(a, h, u, rhs, lift)
        solve["solve_ns"] = time.perf_counter_ns() - start
        solve["setup_ns"] = setup_ns
        x = lift @ np.asarray(solve["x"])
        solve["original_x"] = x.tolist()
        solve["solution_error_mod_kernel"] = float(np.linalg.norm(u @ (x - oracle_x)))
        assert solve["status"] == "converged", (case["name"], name, solve["status"])
        solves[name] = solve
    active = []
    if selected_rank == rank:
        start = time.perf_counter_ns()
        cp = np.linalg.pinv(c, rcond=1e-12) if r else np.zeros((0, m))
        up = u[np.ix_(order, order)]
        active = cp @ up @ cp.T
        lift = cp.T[np.argsort(order), :]
        h = cp @ rhs[order]
        setup_ns = time.perf_counter_ns() - start
        if r:
            assert np.linalg.eigvalsh(active)[0] > 0
        if case["name"] == "shell-weighted":
            # B0=B1-B2+B3, independent of either numerical pseudoinverse.
            z = np.array([{1: 1, 2: -1, 3: 1}[j] / math.sqrt(weights[j]) for j in columns])
            assert np.linalg.norm(active - (np.eye(3) + weights[0] * np.outer(z, z))) <= MATRIX_TOL
        start = time.perf_counter_ns()
        solve = cgls(active, h, u, rhs, lift)
        solve["solve_ns"] = time.perf_counter_ns() - start
        solve["setup_ns"] = selection_ns + setup_ns
        solve["original_x"] = (lift @ np.asarray(solve["x"])).tolist()
        assert solve["status"] == "converged", (case["name"], "hecs", solve["status"])
        assert np.linalg.norm(np.asarray(solve["original_x"]) - oracle_x) <= 1e-7 * max(1., np.linalg.norm(oracle_x))
        solves["hecs"] = solve
    else:
        solves["hecs"] = {"status": "rank_loss"}
    return {**selection, "rank": rank, "selected_rank": selected_rank,
            "edge_permutation": order, "triangle_permutation": columns + [i for i in range(t) if i not in columns],
            "factor": c.tolist(), "factor_shape": [m, r], "factor_nnz": int(np.count_nonzero(c)),
            "exact_weighted_gram": [[str(x) for x in row] for row in exact_factor],
            "weighted_factor_error": factor_error, "selection_ns": selection_ns,
            "spectrum": eig, "active": np.asarray(active).reshape(r, r).tolist() if selected_rank == rank else [],
            "gudhi_betti_f2": betti, "oracle_x": oracle_x.tolist(), "solves": solves}


def run_logged(command, text, directory, name):
    start = time.perf_counter_ns()
    result = subprocess.run(command, input=text, text=True, capture_output=True, check=False)
    (directory / f"{name}.stdout").write_text(result.stdout, encoding="utf-8")
    (directory / f"{name}.stderr").write_text(result.stderr, encoding="utf-8")
    write_json(directory / f"{name}.process.json", {"command": list(map(str, command)),
               "returncode": result.returncode, "wall_ns": time.perf_counter_ns() - start})
    if result.returncode:
        raise RuntimeError(f"{name} failed ({result.returncode}); see retained logs")
    return json.loads(result.stdout)


def verify_rust(actual, expected, case):
    import numpy as np
    for key in ("selected", "pairs", "edge_permutation", "triangle_permutation"):
        assert actual[key] == expected[key], (case["name"], key)
    for key in ("factor", "active"):
        a, e = np.asarray(actual[key]), np.asarray(expected[key])
        assert a.shape == e.shape and np.linalg.norm(a - e) <= MATRIX_TOL * max(1., np.linalg.norm(e)), (case["name"], key)
    for name, solve in actual["solves"].items():
        assert solve["status"] == expected["solves"][name]["status"], (case["name"], name, solve)
        if solve["status"] == "converged":
            b = boundary(case["edges"], case["faces"])
            u = (b * np.asarray(case["weights"])) @ b.T
            assert residual(u, np.asarray(solve["original_x"]), np.asarray(case["rhs"])) <= RESIDUAL_TOL
            if name == "hecs":
                assert np.linalg.norm(np.asarray(solve["original_x"]) - expected["oracle_x"]) <= 1e-7 * max(1., np.linalg.norm(expected["oracle_x"]))


def rust_parity(source, worker, output):
    """Replay Round 1 bytes, without regenerating fixtures or rerunning Julia."""
    output.mkdir(parents=True, exist_ok=False)
    summary = json.loads((source / "summary.json").read_text(encoding="utf-8"))
    if summary["status"] != "passed":
        raise ValueError("Rust parity requires a passed Round-1 summary")
    cases = json.loads((source / "fixtures.json").read_text(encoding="utf-8"))
    rows = []
    write_json(output / "environment.json", {"source": str(source.resolve()),
               "fixtures_sha256": digest(source / "fixtures.json"),
               "round1_summary_sha256": digest(source / "summary.json"),
               "worker_sha256": digest(worker), "rust_source_sha256": digest(ROOT / "tools/reference/hecs.rs")})
    try:
        for row in summary["cases"]:
            case = next(c for c in cases if c["name"] == row["case"])
            name = row["case"] + ("-core" if row["peel"] else "-direct")
            directory = output / name
            directory.mkdir()
            text = (source / name / "input.txt").read_text(encoding="utf-8")
            expected = json.loads((source / name / "python.json").read_text(encoding="utf-8"))
            actual = run_logged([str(worker.resolve())], text, directory, "rust")
            verify_rust(actual, expected, case)
            candidate = run_logged([str(worker.resolve()), "--adjacency"], text, directory, "adjacency")
            verify_rust(candidate, expected, case)
            rows.append({"case": name, "input_sha256": digest(source / name / "input.txt"),
                         "expected_sha256": digest(source / name / "python.json"),
                         "status": "passed", "solves": {k: v["status"] for k, v in actual["solves"].items()}})
    except Exception as error:
        write_json(output / "summary.json", {"status": "failed", "completed": rows, "error": str(error)})
        raise
    write_json(output / "summary.json", {"status": "passed", "cases": rows})
    print(json.dumps({"status": "passed", "rust_parity_cases": len(rows), "output": str(output)}))


def compare_collapse(source, worker, output, samples):
    """One-change, fresh-process paired experiment; retains every sample."""
    output.mkdir(parents=True, exist_ok=False)
    cases = json.loads((source / "fixtures.json").read_text(encoding="utf-8"))
    if json.loads((source / "summary.json").read_text(encoding="utf-8"))["status"] != "passed":
        raise ValueError("requires passed frozen Round-1 inputs")
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    changed = subprocess.check_output(["git", "status", "--porcelain", "--", "tools/reference/hecs.rs", "tools/reproduce_hecs.py"], cwd=ROOT, text=True).strip()
    write_json(output / "environment.json", {"revision": revision, "measured_files_dirty": bool(changed),
        "platform": platform.platform(), "samples": samples, "seed": SEED,
        "worker_sha256": digest(worker), "worker_source_sha256": digest(ROOT / "tools/reference/hecs.rs"),
        "controller_sha256": digest(ROOT / "tools/reproduce_hecs.py"),
        "fixtures_sha256": digest(source / "fixtures.json"),
        "scope": "same Rust process work, dense reference; startup/I/O/JSON export excluded from worker total_ns",
        "peak_process_memory": "not measured; no memory/adoption/scaling claim",
        "reuse": "not measured; no amortization claim"})
    rows = []
    rng = random.Random(SEED)
    try:
        for case in cases:
            if not case["name"].startswith("enriched-"):
                continue
            for peel in (False, True):
                name = case["name"] + ("-core" if peel else "-direct")
                directory = output / name
                directory.mkdir()
                text = (source / name / "input.txt").read_text(encoding="utf-8")
                expected = json.loads((source / name / "python.json").read_text(encoding="utf-8"))
                values = {mode: [] for mode in ("scan", "adjacency")}
                for sample in range(-1, samples):
                    order = ["scan", "adjacency"]
                    if rng.randrange(2):
                        order.reverse()
                    for mode in order:
                        command = [str(worker.resolve())] + (["--adjacency"] if mode == "adjacency" else [])
                        actual = run_logged(command, text, directory, f"{sample + 1:02}-{mode}")
                        verify_rust(actual, expected, case)
                        if sample >= 0:
                            values[mode].append({"selection_ns": actual["selection_ns"], "total_ns": actual["total_ns"]})
                ratios = {key: [a[key] / b[key] for a, b in zip(values["adjacency"], values["scan"])]
                          for key in ("selection_ns", "total_ns")}
                rows.append({"case": name, "samples": values, "paired_ratios": ratios,
                             "median_ratio": {k: statistics.median(v) for k, v in ratios.items()},
                             "range_ratio": {k: [min(v), max(v)] for k, v in ratios.items()}})
        write_json(output / "summary.json", {"status": "passed", "cases": rows,
                   "decision": "experimental only; dense small controls, no production adoption or persistent transfer"})
    except Exception as error:
        write_json(output / "summary.json", {"status": "failed", "completed": rows, "error": str(error)})
        raise
    print(json.dumps({"status": "passed", "paired_cells": len(rows), "output": str(output)}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--author-source", type=Path)
    parser.add_argument("--paper", type=Path)
    parser.add_argument("--julia", type=Path)
    parser.add_argument("--rust-worker", type=Path)
    parser.add_argument("--fixtures", type=Path, help="replay an existing frozen fixtures.json")
    parser.add_argument("--case", help="run only this named fixture")
    parser.add_argument("--parity-from", type=Path, help="Rust-only replay of a passed Round-1 directory")
    parser.add_argument("--compare-from", type=Path, help="one-change paired Rust experiment on frozen enriched inputs")
    parser.add_argument("--samples", type=int, default=10)
    args = parser.parse_args()
    if args.compare_from:
        if not args.rust_worker or args.samples < 1:
            parser.error("--compare-from requires --rust-worker and positive --samples")
        compare_collapse(args.compare_from, args.rust_worker, args.output, args.samples)
        return
    if args.parity_from:
        if not args.rust_worker:
            parser.error("--parity-from requires --rust-worker")
        rust_parity(args.parity_from, args.rust_worker, args.output)
        return
    if not all((args.author_source, args.paper, args.julia)):
        parser.error("Round 1 requires --author-source, --paper and --julia")
    args.output.mkdir(parents=True, exist_ok=False)
    import numpy as np
    import scipy
    import gudhi
    source_files = ["README.md", "src/collapse.jl", "src/cgls.jl", "src/simplicialComplex.jl", "src/utils.jl"]
    # Git blob identities pin the included source, even for an archive checkout.
    source_hashes = {name: digest(args.author_source / name) for name in source_files}
    if digest(args.paper) != PAPER_SHA256 or any(source_hashes[k] != v for k, v in AUTHOR_HASHES.items()):
        raise ValueError("paper or author source does not match the pinned version")
    write_json(args.output / "environment.json", {
        "python": sys.version, "platform": platform.platform(), "numpy": np.__version__,
        "scipy": scipy.__version__, "gudhi": gudhi.__version__, "seed": SEED,
        "paper_url": "https://arxiv.org/pdf/2401.15492v1", "paper_sha256": digest(args.paper),
        "author_commit": AUTHOR_COMMIT, "author_files": source_hashes,
        "harness": {p: digest(ROOT / p) for p in ["tools/reproduce_hecs.py", "tools/reference/hecs.jl"]},
        "matrix_tolerance": MATRIX_TOL, "original_residual": "||Ux-b||2 / max(1,||b||2) <= 1e-9",
        "timing_class": "exploratory resource snapshot; no ranking or adoption claim",
        "incomplete_cholesky": "not run; narrowed comparison", "persistent_adaptation": "not run"})
    cases = json.loads(args.fixtures.read_text(encoding="utf-8")) if args.fixtures else fixtures()
    if args.case:
        cases = [c for c in cases if c["name"] == args.case]
        if not cases:
            raise ValueError("unknown fixture name")
    write_json(args.output / "fixtures.json", cases)
    summary = []
    try:
        for case in cases:
            for peel in (False, True):
                directory = args.output / (case["name"] + ("-core" if peel else "-direct"))
                directory.mkdir()
                text = protocol(case, peel)
                (directory / "input.txt").write_text(text, encoding="utf-8")
                expected = python_case(case, peel)
                write_json(directory / "python.json", expected)
                julia = run_logged([str(args.julia.resolve()), "--startup-file=no", "--threads=1",
                                    str(ROOT / "tools/reference/hecs.jl"), str(args.author_source.resolve())],
                                   text, directory, "julia")
                assert sorted(julia["selected"]) == sorted(expected["selected"])
                assert replay(case["edges"], case["faces"], julia["selected"], julia["pairs"])
                assert julia["factor_error"] <= MATRIX_TOL
                assert julia["selected_rank"] == expected["selected_rank"]
                if "full_rank" in julia:
                    assert julia["full_rank"] == expected["rank"]
                    assert julia["paper_upper_nnz"] == 0
                if julia["status"] == "converged":
                    assert julia["original_residual"] <= RESIDUAL_TOL
                else:
                    assert julia["status"] in ("zero_rhs_source_nonfinite", "rank_loss", "empty_source_unsupported")
                if args.rust_worker:
                    rust = run_logged([str(args.rust_worker.resolve())], text, directory, "rust")
                    verify_rust(rust, expected, case)
                summary.append({"case": case["name"], "peel": peel, "rank": expected["rank"],
                                "selected_rank": expected["selected_rank"], "factor_error": expected["weighted_factor_error"],
                                "julia": julia["status"], "rust_parity": bool(args.rust_worker),
                                "solves": {k: {"status": v["status"], "iterations": len(v.get("trace", [])) - 1,
                                                "original_residual": v.get("trace", [{}])[-1].get("original_residual")}
                                           for k, v in expected["solves"].items()}})
        write_json(args.output / "summary.json", {"status": "passed", "cases": summary})
    except Exception as error:
        write_json(args.output / "summary.json", {"status": "failed", "completed": summary,
                                                  "error": str(error)})
        raise
    print(json.dumps({"status": "passed", "cases": len(summary), "output": str(args.output)}))


if __name__ == "__main__":
    main()
