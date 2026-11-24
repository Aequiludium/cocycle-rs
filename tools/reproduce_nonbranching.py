#!/usr/bin/env python3
"""Bounded identity-metric replay of Botnan-Dong Algorithms 1-3 (#147).

Literal pseudocode, corrected weak reduction, exact oracle, author replay and
Rust parity are separate evidence. Nothing here selects a production backend.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import importlib.util
import itertools
import json
import os
from pathlib import Path
import platform
import random
import shutil
import subprocess
import sys
import time

import numpy as np
import scipy

from compare_persistent_laplacians import exact_projection, peak_memory, rational_nullspace, rref

PAPER_SHA = "343e6c6e076ce0a2d0022f45fab35c0741dcde0c7af795f5d253b4640296c483"
AUTHOR_SHA = "03af3d06ab5fcea93e7e7ac7dd2e269ad48ef2fd07df3b9bc8162e021158edd8"
AUTHOR_COMMIT = "13d855749bddc5e4d8cca656d5e0775a56a581a1"
ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    def encode(x):
        if isinstance(x, np.ndarray):
            return x.tolist()
        if isinstance(x, Fraction):
            return str(x)
        if isinstance(x, np.generic):
            return x.item()
        raise TypeError(type(x).__name__)
    Path(path).write_text(json.dumps(value, default=encode, indent=2) + "\n", encoding="utf-8")


def validate(b1, d):
    if b1.ndim != 2 or d.ndim != 2 or b1.shape[1] != d.shape[1]:
        raise ValueError("incompatible boundary shapes")
    if not np.isin(b1, [-1, 0, 1]).all() or not np.isin(d, [-1, 0, 1]).all():
        raise ValueError("only signed unit incidence is supported")
    if np.any(np.count_nonzero(d, axis=1) > 2):
        raise ValueError("branching outside row")


def literal_reorientation(d, mutate=True):
    """Algorithm 1 lines 3-26, including its D mutation and flag multiplication.

    mutate=False isolates the additional row-order problem after removing only
    the double sign application. This is a diagnostic, not the corrected route.
    """
    d = d.copy()
    flag = np.ones(d.shape[1], dtype=int)
    vertices = set(range(d.shape[1]))
    for row in range(len(d)):
        edge = np.flatnonzero(d[row]).tolist()
        if len(edge) == 2:
            u, v = edge
            if abs(flag[u] * d[row, u] + flag[v] * d[row, v]) == 2:
                if v in vertices:
                    change = v
                elif u in vertices:
                    change = u
                else:
                    return False, flag
                if mutate:
                    d[:, change] *= -1
                flag[change] = -1
        vertices.difference_update(edge)
    return True, flag


def weak_reduction(b1, d, literal=False):
    """Paper component construction, with separately labeled corrections.

    Corrected: traverse each unsigned component before testing its signs; keep
    D immutable, accumulate flags once, and actually form R = D E V. No signed
    union-find is used. Original zero columns bypass reduction exactly once.
    """
    validate(b1, d)
    p = d.shape[1]
    original_zero = np.flatnonzero(~np.any(d, axis=0)).tolist()
    nonzero = np.flatnonzero(np.any(d, axis=0)).tolist()
    a = d[:, nonzero]
    adjacency = [[] for _ in nonzero]
    forced = set()
    for row in a:
        indices = np.flatnonzero(row)
        if len(indices) == 1:
            forced.add(int(indices[0]))
        elif len(indices) == 2:
            u, v = map(int, indices)
            sign = -int(row[u]) * int(row[v])
            adjacency[u].append((v, sign))
            adjacency[v].append((u, sign))
    count = len(nonzero)
    e = np.eye(count, dtype=int)
    v = np.eye(count, dtype=int)
    r = a.copy()
    visited = set()
    components = []
    zero_reduced = []
    for root in range(count):
        if root in visited:
            continue
        signs = {root: 1}
        pending = [root]
        consistent = True
        while pending:
            u = pending.pop()
            visited.add(u)
            for other, relation in adjacency[u]:
                wanted = signs[u] * relation
                if other in signs:
                    consistent &= signs[other] == wanted
                else:
                    signs[other] = wanted
                    pending.append(other)
        support = sorted(signs)
        row_singular = bool(forced.intersection(support))
        if literal and not row_singular:
            rows = np.flatnonzero(np.any(a[:, support], axis=1))
            consistent, flags = literal_reorientation(a[np.ix_(rows, support)])
            signs = dict(zip(support, map(int, flags)))
        eligible = consistent and not row_singular
        components.append({"support": [nonzero[i] for i in support],
                           "signs": [signs[i] for i in support],
                           "kind": "row-singular" if row_singular else
                           "regulable" if consistent else "irregular"})
        if eligible:
            last = support[-1]
            for i in support:
                e[i, i] = signs[i]
                v[i, last] = 1
            r[:, last] = 0
            zero_reduced.append(last)
    if not literal:
        r = a @ e @ v
    # Algorithm 3 selects zero columns of R, with C0 handled separately.
    z = np.zeros((p, len(original_zero) + len(zero_reduced)), dtype=int)
    for j, col in enumerate(original_zero):
        z[col, j] = 1
    for j, col in enumerate(zero_reduced, len(original_zero)):
        z[nonzero, j] = (e @ v)[:, col]
    metric = z.T @ z
    factor = b1 @ z
    upper = [[sum(Fraction(int(factor[i, c]) * int(factor[j, c]),
                           int(metric[c, c])) for c in range(z.shape[1]))
              for j in range(len(b1))] for i in range(len(b1))]
    return {"R": r, "E": e, "V": v, "nonzero": nonzero,
            "original_zero": original_zero, "transformed_zero": zero_reduced,
            "components": components, "Z": z, "metric": metric,
            "factor_integer": factor, "upper_exact": upper,
            "identity": bool(np.array_equal(r, a @ e @ v)),
            "kernel": bool(not np.any(d @ z))}


def fixtures():
    """Paper labels are frozen; synthetic controls are explicitly distinguished."""
    def case(name, b, d, origin="project algebraic control", **labels):
        return dict(name=name, B1=np.asarray(b, dtype=int), D=np.asarray(d, dtype=int),
                    origin=origin, **labels)
    cases = [
        case("path", [[-1, 0], [0, 1]], [[1, -1]]),
        case("same-sign", [[1, 0]], [[1, 1]]),
        case("double-sign", [[1, 0, 0]], [[1, 1, 0], [0, 1, 1]]),
        case("row-order", np.eye(4, dtype=int),
             [[1, 1, 0, 0], [0, 0, 1, 1], [1, 0, 1, 0]]),
        case("signed-consistent", np.eye(3, dtype=int),
             [[1, 1, 0], [0, 1, 1], [1, 0, -1]]),
        case("signed-contradiction", np.eye(2, dtype=int), [[1, 1], [1, -1]]),
        case("singleton-and-free", [[1, 0, 1], [0, 1, -1]], [[1, 0, 0], [0, 0, 0]]),
        case("all-zero", np.eye(3, dtype=int), np.zeros((2, 3), dtype=int)),
        case("equal-domain", np.eye(3, dtype=int), np.zeros((0, 3), dtype=int)),
        case("no-upper-cells", np.zeros((2, 0), dtype=int), np.zeros((1, 0), dtype=int)),
        case("empty-source", np.zeros((0, 2), dtype=int), [[1, -1]]),
    ]
    mobius = np.array([[1, 0, 0, 0, 0], [-1, -1, 0, 0, 0], [1, 0, 0, 0, 1],
                       [0, 1, 1, 0, 0], [0, 1, 0, 0, 0], [0, 0, 1, 0, 0],
                       [0, 0, -1, -1, 0], [0, 0, 0, 1, 1], [0, 0, 0, 1, 0],
                       [0, 0, 0, 0, -1]])
    mask = np.array([1, 0, 0, 0, 1, 1, 0, 0, 1, 1], dtype=bool)
    cases.append(case("paper-3.9", mobius[mask], mobius[~mask], "paper Example 3.9",
                      B=mobius, mask=mask, degree=1,
                      rows=["34", "14", "13", "24", "12", "45", "25", "35", "23", "15"],
                      columns=["134", "124", "245", "235", "135"]))
    cubical = np.array([[-1, 0, 0, 0], [0, -1, 0, 0], [-1, 0, 0, 0],
                        [1, -1, 0, 0], [0, 1, 0, 0], [1, 0, -1, 0],
                        [0, 1, 0, -1], [0, 0, -1, 0], [0, 0, 1, -1],
                        [0, 0, 0, 1], [0, 0, 1, 0], [0, 0, 0, 1]])
    mask = np.arange(12) != 3
    cases.append(case("paper-6.1", cubical[mask], cubical[~mask], "paper Example 6.1",
                      B=cubical, mask=mask, degree=1,
                      rows=["12", "23", "41", "52", "63", "45", "56", "74", "85", "96", "78", "89"],
                      columns=["4512", "5623", "7845", "8956"]))
    # Controls from #143: triangle upper and tetrahedral shell.
    cases.append(case("triangle", [[1], [-1], [1]], np.zeros((0, 1), dtype=int)))
    edges = list(itertools.combinations(range(4), 2))
    faces = list(itertools.combinations(range(4), 3))
    shell = np.zeros((6, 4), dtype=int)
    for j, face in enumerate(faces):
        for removed in range(3):
            shell[edges.index(face[:removed] + face[removed + 1:]), j] = (-1) ** removed
    cases.append(case("shell", shell, np.zeros((0, 4), dtype=int)))
    disk = np.array([[1, 0], [-1, 1], [0, -1], [1, 0], [0, 1]])
    cases.append(case("disk", disk[[0, 2, 3, 4]], disk[[1]],
                      "project triangulated disk", degree=1))
    cases.append(case("r-reorientation", np.eye(3, dtype=int), [[1, 1, 0], [1, 0, 1]]))
    rng = random.Random(147)
    for number in range(60):
        p = rng.randrange(1, 9)
        d = np.zeros((rng.randrange(0, 10), p), dtype=int)
        for row in d:
            for j in rng.sample(range(p), rng.randrange(min(2, p) + 1)):
                row[j] = rng.choice([-1, 1])
        b = np.array([[rng.choice([-1, 0, 1]) for _ in range(p)] for _ in range(3)])
        cases.append(case(f"seed147-{number}", b, d))
    # Signed basis and row permutations of the worked and synthetic examples.
    for base in cases[:7] + cases[11:17]:
        b, d = base["B1"], base["D"]
        perm = list(reversed(range(d.shape[1])))
        signs = np.array([(-1) ** j for j in range(d.shape[1])])
        cases.append(case(base["name"] + "-permuted-oriented", b[:, perm] * signs,
                          d[::-1, perm] * signs))
    return cases


def verify(case, result):
    b, d = case["B1"], case["D"]
    n, p = b.shape
    oracle = exact_projection(b, d, np.zeros((0, n), dtype=int))
    assert result["identity"] and result["kernel"], case["name"]
    assert len(rational_nullspace(d.tolist(), p)) == result["Z"].shape[1], case["name"]
    nonzero_r = np.count_nonzero(np.any(result["R"], axis=0))
    assert len(rref(result["R"].tolist(), len(result["nonzero"]))[1]) == nonzero_r
    assert result["upper_exact"] == oracle, case["name"]
    z = result["Z"]
    assert np.array_equal(z.T @ z, np.diag(np.diag(result["metric"]))), case["name"]
    u = np.array(oracle, dtype=float).reshape(n, n)
    projector = np.eye(p) - np.linalg.pinv(d.astype(float), rcond=1e-12) @ d
    numerical = b @ projector @ b.T
    error = float(np.linalg.norm(u - numerical) / max(1, np.linalg.norm(u)))
    assert error <= 1e-10, (case["name"], error)
    action_residual = float(np.linalg.norm(d @ projector) / max(1, np.linalg.norm(d)))
    assert action_residual <= 1e-9
    return {"relative_frobenius": error, "projector_residual": action_residual,
            "nullity": z.shape[1], "spectrum": np.linalg.eigvalsh(u)}


def author_replay(path, case):
    """Execute unchanged pinned source on frozen paper matrices and row masks."""
    if digest(path) != AUTHOR_SHA:
        raise ValueError("author source SHA-256 mismatch")
    spec = importlib.util.spec_from_file_location("pinned_cubical", path)
    author = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(author)
    from scipy.sparse import csr_matrix
    b, mask = case["B"], case["mask"]
    assert mask.dtype == bool and mask.shape == (len(b),)
    op = author.UpPersistentLaplacianSparse(csr_matrix(b), filt_idx_1_2=mask)
    op.find_D()
    op.find_root()
    op.uf.find_loops()
    op.uf.find_loops_roots()
    op.find_rel_bdry_opr()
    op.find_diag()
    op.find_up_persistent_Laplacian()
    return {"D": op.D_mtx.toarray(), "mask": mask, "G": op.G.toarray(),
            "factor_integer": op.rel_bdry_opr.toarray(), "inverse_metric": op.diag_array,
            "upper": op.up_lap.toarray()}


def author_grid(path):
    """Regenerated 2x2 lattice, not the author's chest image data."""
    if digest(path) != AUTHOR_SHA:
        raise ValueError("author source SHA-256 mismatch")
    spec = importlib.util.spec_from_file_location("pinned_cubical_grid", path)
    author = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(author)
    x = np.zeros((5, 5), dtype=int)
    x[1::2, 1::2] = 2
    x[1, 2] = 1  # the only edge outside K; all four squares enter at 2
    boundary = author.BoundaryMatrix(x, filt_1=0, filt_2=2)
    boundary.chains_count()
    boundary.find_boundary_opr()
    boundary.find_cell_idx()
    b, mask = boundary.bdry_opr.toarray(), boundary.filt_idx_1_2
    assert b.shape == (12, 4) and mask.shape == (12,)
    assert np.flatnonzero(~mask).tolist() == [3]
    case = dict(name="regenerated-cubical-grid", B=b, mask=mask, B1=b[mask], D=b[~mask])
    return case, {"X": x, "source_shape": b.shape, "mask": mask,
                  "origin": "project regenerated paper-style 2x2 grid"}


def topology_checks():
    """GUDHI checks declared simplicial topology, never real spectra."""
    import gudhi
    from compare_persistent_laplacians import boundary, ordered
    examples = [
        ("path", {(0,): 0, (2,): 0, (1,): 1, (0, 1): 2, (1, 2): 2}, 0, 2, 0),
    ]
    for name, faces, source_edges in [
        ("triangle", [(0, 1, 2)], [(0, 1), (0, 2), (1, 2)]),
        ("disk", [(0, 1, 2), (0, 2, 3)], [(0, 1), (0, 3), (1, 2), (2, 3)]),
        ("mobius", [(1, 3, 4), (1, 2, 4), (2, 4, 5), (2, 3, 5), (1, 3, 5)],
         [(3, 4), (1, 2), (4, 5), (2, 3), (1, 5)]),
        ("shell", list(itertools.combinations(range(4), 3)), list(itertools.combinations(range(4), 2))),
    ]:
        cells = {simplex: (0 if len(simplex) == 1 or simplex in source_edges else 2)
                 for face in faces for length in range(1, 4)
                 for simplex in itertools.combinations(face, length)}
        examples.append((name, cells, 0, 2, 1))
    rows = []
    for name, cells, s, t, q in examples:
        k, l = ordered(cells, s, q), ordered(cells, t, q)
        up = boundary(l, ordered(cells, t, q + 1))
        lower = boundary(ordered(cells, s, q - 1), k)
        retained = [l.index(simplex) for simplex in k]
        outside = [i for i in range(len(l)) if i not in retained]
        result = weak_reduction(up[retained], up[outside])
        upper_rank = len(rref(result["upper_exact"], len(k))[1])
        real_beta = len(k) - len(rref(lower.tolist(), len(k))[1]) - upper_rank
        for field in [2, 11]:
            st = gudhi.SimplexTree()
            for simplex, birth in cells.items():
                st.insert(simplex, filtration=birth)
            st.persistence(homology_coeff_field=field, persistence_dim_max=True)
            count = sum(birth <= s and death > t for birth, death in st.persistence_intervals_in_dimension(q))
            # The Mobius boundary map has degree two: its image vanishes in F2.
            expected = 0 if name == "mobius" and field == 2 else real_beta
            assert count == expected, (name, field, count, expected)
            rows.append(dict(name=name, degree=q, pair=[s, t], field=field,
                             persistent_rank=int(count), real_nullity=real_beta, version=gudhi.__version__,
                             spectral_oracle=False, cells=[dict(simplex=k, birth=v) for k, v in cells.items()]))
    return rows


def protocol(case):
    b, d = case["B1"], case["D"]
    return " ".join(map(str, [len(b), len(d), b.shape[1], *b.flat, *d.flat])) + "\n"


def worker_run(worker, case, mode, repeat=1, summary=False):
    start = time.perf_counter_ns()
    command = [str(worker), mode, str(repeat)] + (["--summary"] if summary else [])
    with subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                          stderr=subprocess.PIPE, text=True) as proc:
        stdout, stderr = proc.communicate(protocol(case))
        memory = peak_memory(proc)
    wall = time.perf_counter_ns() - start
    if proc.returncode:
        raise RuntimeError(f"worker exit {proc.returncode}: {stderr}")
    result = json.loads(stdout)
    result["process_wall_ns"] = wall
    result["peak_working_set_bytes"] = memory
    return result


def benchmark_summary(measurements):
    summaries = []
    for p in sorted({row["workload"] for row in measurements}):
        groups = {mode: [r for r in measurements if r["workload"] == p and r["mode"] == mode]
                  for mode in ["baseline", "compact"]}
        for field in ["total_ns", "process_wall_ns", "peak_working_set_bytes"]:
            values = {mode: [r[field] for r in rows] for mode, rows in groups.items()}
            if any(x is None for rows in values.values() for x in rows):
                continue
            paired = [values["baseline"][i] / values["compact"][i] for i in range(10)]
            summaries.append(dict(workload=p, field=field,
                                  baseline_p10_median_p90=np.quantile(values["baseline"], [.1, .5, .9]),
                                  compact_p10_median_p90=np.quantile(values["compact"], [.1, .5, .9]),
                                  paired_baseline_over_compact=paired,
                                  paired_p10_median_p90=np.quantile(paired, [.1, .5, .9])))
    return summaries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--paper", type=Path, required=True)
    parser.add_argument("--author-source", type=Path, required=True)
    parser.add_argument("--benchmark", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    out = args.output.resolve()
    try:
        if digest(args.paper) != PAPER_SHA or digest(args.author_source) != AUTHOR_SHA:
            raise ValueError("pinned paper/author identity mismatch")
        shutil.copy2(args.paper, out / "paper.pdf")
        shutil.copy2(args.author_source, out / "author-original.py")
        save(out / "environment.json", {"python": sys.version, "numpy": np.__version__,
             "scipy": scipy.__version__, "platform": platform.platform(),
             "paper_sha256": PAPER_SHA, "author_sha256": AUTHOR_SHA,
             "author_commit": AUTHOR_COMMIT, "command": sys.argv,
             "thread_environment": {k: os.environ.get(k) for k in
                                    ["OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"]}})
        for name in ["reproduce_nonbranching.py", "compare_persistent_laplacians.py"]:
            shutil.copy2(ROOT / "tools" / name, out / name)
        shutil.copy2(ROOT / "tools/reference/nonbranching.rs", out / "nonbranching.rs")
        cases = fixtures()
        save(out / "inputs.json", cases)
        rows = []
        frozen = []
        for case in cases:
            result = weak_reduction(case["B1"], case["D"])
            checks = verify(case, result)
            literal = weak_reduction(case["B1"], case["D"], literal=True)
            literal["upper_matches_exact"] = literal["upper_exact"] == result["upper_exact"]
            frozen.append({"name": case["name"], "corrected": result, "literal": literal})
            rows.append(dict(name=case["name"], **checks,
                             literal_identity=literal["identity"], literal_kernel=literal["kernel"],
                             literal_upper=literal["upper_matches_exact"]))
        save(out / "round1.json", frozen)
        save(out / "checks.json", rows)
        # Replay author code only on oriented cubical inputs: paper and equal-domain.
        author_cases = [next(c for c in cases if c["name"] == "paper-6.1")]
        cube = author_cases[0]
        author_cases.append(dict(cube, name="cubical-equal-domain", mask=np.ones(12, dtype=bool),
                                 B1=cube["B"], D=np.zeros((0, 4), dtype=int)))
        grid, grid_metadata = author_grid(args.author_source)
        save(out / "author-grid.json", dict(grid_metadata, **grid))
        author_cases.append(grid)
        author_results = []
        for case in author_cases:
            result = author_replay(args.author_source, case)
            expected = weak_reduction(case["B1"], case["D"])
            np.testing.assert_allclose(result["upper"], np.array(expected["upper_exact"], dtype=float),
                                       rtol=1e-10, atol=1e-10)
            author_results.append(dict(name=case["name"], **result))
        save(out / "author-replay.json", author_results)
        save(out / "gudhi-topology.json", topology_checks())
        save(out / "author-compatibility.json", {"patch": "none", "entry": "src/sparse_uplap.py",
             "notebook": "experiments/compute_time_1.ipynb",
             "scope": "unchanged notebook core call chain on frozen Example 6.1 matrices/masks; no chest image benchmark"})
        save(out / "corrections.json", [
            {"algorithm": 1, "lines": "10,12-19", "literal": "multiply flag and mutate D",
             "corrected": "immutable D; apply accumulated orientation once",
             "counterexample": "double-sign; signed-contradiction"},
            {"algorithm": 1, "lines": "6-26", "literal": "one arbitrary row scan",
             "corrected": "traverse the unsigned component, then check every signed relation",
             "counterexample": "row-order, even with D mutation disabled"},
            {"algorithm": 2, "lines": "3,18-20", "literal": "R=D with selected columns zeroed",
             "corrected": "form R=D E V, including other reoriented columns",
             "counterexample": "r-reorientation (valid kernel but R differs)"},
            {"algorithm": 3, "lines": "5,10,14", "literal": "B_q B_q^T; D_(q+1)^L; I_j",
             "corrected": "B_(q+1) B_(q+1)^T; outside-row D; mapped component C_j",
             "scope": "identity weights; C0 appended once; upper factor uses inverse metric"}])
        # Algorithm 3 cannot literally multiply the printed B_q at q=0 and
        # return an operator on C_q. Retain the concrete dimensional audit.
        save(out / "algorithm3-index-audit.json", {
            "control": "same-stage ordinary path 0-1-2, q=0",
            "B_q_shape": [0, 3], "B_q_plus_1": [[-1, 0], [1, -1], [0, 1]],
            "literal_line5_output_shape": [0, 0], "required_output_shape": [3, 3],
            "corrected_upper": [[1, -1, 0], [-1, 2, -1], [0, -1, 1]],
            "line10": "use outside-row D_(q+1)^(L,K), with nonzero-column index map",
            "line14": "component support I_j is mapped to original upper columns C_j"})
        source = ROOT / "tools/reference/nonbranching.rs"
        worker = out / ("nonbranching.exe" if os.name == "nt" else "nonbranching")
        command = ["rustc", "--edition=2024", "-O", "-Dwarnings", str(source), "-o", str(worker)]
        build = subprocess.run(command, text=True, capture_output=True)
        (out / "build.log").write_text(build.stdout + build.stderr, encoding="utf-8")
        build.check_returncode()
        save(out / "rust-build.json", {"command": command, "worker_sha256": digest(worker),
             "source_sha256": digest(source), "version": subprocess.check_output(["rustc", "--version"], text=True)})
        parity = []
        for case, expected in zip(cases, frozen):
            for mode in ["baseline", "compact"]:
                result = worker_run(worker, case, mode)
                ref = expected["corrected"]
                n, p = case["B1"].shape
                u = np.array(ref["upper_exact"], dtype=float).reshape(n, n)
                np.testing.assert_allclose(np.array(result["upper"]).reshape(n, n), u, rtol=1e-10, atol=1e-10)
                assert result["original_zero"] == ref["original_zero"]
                assert result["components"] == ref["components"]
                assert result["metric"] == np.diag(ref["metric"]).tolist()
                normalized_factor = ref["factor_integer"] / np.sqrt(np.diag(ref["metric"]))
                np.testing.assert_allclose(result["factor"], normalized_factor.ravel(), rtol=1e-10, atol=1e-10)
                np.testing.assert_allclose(result["action"], u @ np.ones(n), rtol=1e-10, atol=1e-10)
                if mode == "baseline":
                    for key in ["R", "E", "V", "Z"]:
                        assert result[key] == ref[key].ravel().tolist(), (case["name"], key)
                parity.append(dict(name=case["name"], mode=mode, **result))
        save(out / "rust-parity.json", parity)
        rejection = []
        for name, data in [("branching", "1 1 3 1 0 0 1 1 1"),
                           ("nonunit", "1 1 1 1 2"),
                           ("dense-cap", "257 0 257"),
                           ("trailing-input", "0 0 0 7")]:
            for mode in ["baseline", "compact"]:
                proc = subprocess.run([str(worker), mode, "1"], input=data, text=True, capture_output=True)
                assert proc.returncode != 0, name
                rejection.append(dict(name=name, mode=mode, stderr=proc.stderr, returncode=proc.returncode))
        save(out / "rejections.json", rejection)
        if args.benchmark:
            measurements = []
            for p in [32, 128]:
                d = np.zeros((p - 1, p), dtype=int)
                for i in range(p - 1):
                    d[i, i:i + 2] = [1, -1]
                case = dict(name=f"path-{p}", B1=np.eye(8, p, dtype=int), D=d)
                save(out / f"benchmark-input-{p}.json", case)
                for mode in ["baseline", "compact"]:
                    worker_run(worker, case, mode, summary=True)  # discarded warmup
                for sample in range(10):
                    for mode in (["baseline", "compact"] if sample % 2 == 0 else ["compact", "baseline"]):
                        result = worker_run(worker, case, mode, summary=True)
                        np.testing.assert_allclose(np.array(result["upper"]).reshape(8, 8),
                                                   np.ones((8, 8)) / p, atol=1e-10)
                        spectrum_start = time.perf_counter_ns()
                        values = np.linalg.eigvalsh(np.array(result["upper"]).reshape(8, 8))
                        spectrum_ns = time.perf_counter_ns() - spectrum_start
                        measurements.append(dict(workload=p, sample=sample, mode=mode,
                                                 spectrum_python_ns=spectrum_ns, spectrum=values, **result))
            save(out / "benchmark.json", measurements)
            save(out / "benchmark-summary.json", benchmark_summary(measurements))
        save(out / "computation-record.json", {"status": "passed", "evidence_level": "bounded-check",
             "exact_checks": len(cases), "rust_comparisons": len(parity),
             "literal_identity_failures": sum(not r["literal_identity"] for r in rows),
             "literal_upper_failures": sum(not r["literal_upper"] for r in rows),
             "input_sha256": digest(out / "inputs.json"),
             "exclusions": ["Algorithm 7 GKL", "arbitrary weights", "production dispatch",
                            "full chest X-ray notebook", "unbounded scalability"],
             "author_replay": len(author_results)})
        print(json.dumps({"output": str(out), "cases": len(cases), "rust_parity": len(parity),
                          "literal_identity_failures": sum(not r["literal_identity"] for r in rows)}))
    except BaseException as error:
        save(out / "failure.json", {"type": type(error).__name__, "message": str(error)})
        raise


if __name__ == "__main__":
    main()
