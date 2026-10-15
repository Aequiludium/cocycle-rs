#!/usr/bin/env python3
"""Independent real persistent Laplacian reproduction, not a PH runtime backend.

Requires NumPy; --gudhi additionally checks torsion-free fixtures using GUDHI.
The worker is the release persistent_laplacian example. New output directories
retain fixtures, raw worker output, numerical comparisons and adverse timings.
"""
from __future__ import annotations

import argparse
import ast
import ctypes
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time


def fixtures():
    cases = [
        ("path", {(0,): 0, (2,): 0, (1,): 1, (0, 1): 2, (1, 2): 2}),
        ("triangle", {(0,): 0, (1,): 0, (2,): 0, (0, 1): 1,
                      (0, 2): 1, (1, 2): 1, (0, 1, 2): 2}),
        ("singular", {(0,): 0, (1,): 1, (2,): 1, (1, 2): 2}),
        ("zero", {(0,): 0, (1,): 1}),
        ("empty", {}),
    ]
    rng = random.Random(124)
    for number in range(8):
        n = 4 + number % 3
        cells = {(v,): rng.randrange(2) for v in range(n)}
        for length in range(2, 5):
            for simplex in itertools.combinations(range(n), length):
                faces = [simplex[:i] + simplex[i + 1:] for i in range(length)]
                if all(f in cells for f in faces) and rng.random() < 0.85:
                    cells[simplex] = max(cells[f] for f in faces) + rng.randrange(3)
        cases.append((f"seed124-{number}", cells))
    facets = [(0, 1, 2), (0, 1, 3), (0, 2, 4), (0, 3, 5), (0, 4, 5),
              (1, 2, 5), (1, 3, 4), (1, 4, 5), (2, 3, 4), (2, 3, 5)]
    rp2 = {s: len(s) - 1 for facet in facets for length in range(1, 4)
           for s in itertools.combinations(facet, length)}
    cases.append(("rp2-field-counterexample", rp2))
    # All simplex dimensions enter at separate scales: a torsion-free sphere
    # followed by its tetrahedral filling. It also exercises q=2 and q=3.
    cases.append(("tetrahedron", {s: len(s) - 1 for length in range(1, 5)
                                  for s in itertools.combinations(range(4), length)}))
    return cases


def ordered(cells, scale, degree):
    return sorted((s for s, v in cells.items() if v <= scale and len(s) == degree + 1),
                  key=lambda s: (cells[s], tuple(-v for v in reversed(s))))


def boundary(rows, columns):
    import numpy as np
    matrix = np.zeros((len(rows), len(columns)), dtype=int)
    index = {s: i for i, s in enumerate(rows)}
    for j, s in enumerate(columns):
        if len(s) == 1:
            continue  # ordinary nonaugmented B0: 0 x n
        for removed in range(len(s)):
            matrix[index[s[:removed] + s[removed + 1:]], j] = (-1) ** removed
    return matrix


def rref(matrix, columns=None):
    a = [[Fraction(int(x)) if not isinstance(x, Fraction) else x for x in row]
         for row in matrix]
    width = len(a[0]) if a else (columns or 0)
    pivots = []
    for col in range(width):
        pivot = next((i for i in range(len(pivots), len(a)) if a[i][col]), None)
        if pivot is None:
            continue
        row = len(pivots)
        a[row], a[pivot] = a[pivot], a[row]
        factor = a[row][col]
        a[row] = [x / factor for x in a[row]]
        for i in range(len(a)):
            if i != row and a[i][col]:
                factor = a[i][col]
                a[i] = [x - factor * y for x, y in zip(a[i], a[row])]
        pivots.append(col)
        if len(pivots) == len(a):
            break
    return a, pivots


def rational_nullspace(matrix, columns):
    a, pivots = rref(matrix, columns)
    free = [j for j in range(columns) if j not in pivots]
    vectors = []
    for j in free:
        vector = [Fraction(0) for _ in range(columns)]
        vector[j] = Fraction(1)
        for i, p in enumerate(pivots):
            vector[p] = -a[i][j]
        vectors.append(vector)
    return vectors


def exact_projection(r, f, lower):
    """RZ (Z^T Z)^-1 (RZ)^T + lower^T lower, without Schur or SVD."""
    import numpy as np
    n, columns = r.shape
    vectors = rational_nullspace(f.tolist(), columns)
    rank = len(vectors)
    result = [[Fraction(int(x)) for x in row] for row in (lower.T @ lower).tolist()]
    if not rank:
        return result
    gram = [[sum(a * b for a, b in zip(vectors[i], vectors[j]))
             for j in range(rank)] for i in range(rank)]
    augmented = [row + [Fraction(i == j) for j in range(rank)]
                 for i, row in enumerate(gram)]
    inverse_rows, pivots = rref(augmented)
    assert pivots[:rank] == list(range(rank))
    inverse = [row[rank:] for row in inverse_rows]
    rz = [[sum(Fraction(int(x)) * y for x, y in zip(r[i], z)) for z in vectors]
          for i in range(n)]
    product = [[sum(rz[i][a] * inverse[a][b] for a in range(rank))
                for b in range(rank)] for i in range(n)]
    for i in range(n):
        for j in range(n):
            result[i][j] += sum(product[i][b] * rz[j][b] for b in range(rank))
    return result


def oracle(cells, s, t, q):
    import numpy as np
    k = ordered(cells, s, q)
    l = ordered(cells, t, q)
    lower = boundary(ordered(cells, s, q - 1), k)
    up = boundary(l, ordered(cells, t, q + 1))
    retained = [l.index(simplex) for simplex in k]
    removed = [i for i in range(len(l)) if i not in retained]
    r, f = up[retained, :], up[removed, :]
    a, e, c = r @ r.T, r @ f.T, f @ f.T
    direct = a.astype(float)
    if removed:
        direct -= e @ np.linalg.lstsq(c, e.T, rcond=1e-12)[0]
    direct += lower.T @ lower
    if up.shape[1]:
        _, singular, vh = np.linalg.svd(f, full_matrices=True)
        threshold = 1e-12 * (singular[0] if singular.size else 0)
        rank_f = int(np.count_nonzero(singular > threshold))
        projected = r @ vh[rank_f:].T
        projection = projected @ projected.T + lower.T @ lower
    else:
        projection = (lower.T @ lower).astype(float)
    # Integer incidence ranks over Q equal ranks over R. This is independent
    # of both floating Schur and the numerical nullspace construction.
    beta = (len(k) - len(rref(lower.tolist(), len(k))[1])
            - len(rref(up.tolist(), up.shape[1])[1])
            + len(rref(f.tolist(), f.shape[1])[1]))
    exact = exact_projection(r, f, lower) if up.shape[1] <= 18 else None
    return k, direct, projection, beta, exact, up, retained, lower


def protocol(cells, q, starts, target):
    words = [str(q), str(len(starts)), str(len(cells)), str(target)]
    words += [str(s) for s in starts]
    for simplex, value in cells.items():
        words += [str(value), str(len(simplex)), *map(str, simplex)]
    return " ".join(words) + "\n"


def parse_worker(text):
    import numpy as np
    lines = iter(text.splitlines())
    for line in lines:
        header = line.split()
        assert header[0] == "PAIR"
        q, s, t, n = int(header[1]), float(header[2]), float(header[3]), int(header[4])
        basis_line = iter(next(lines).split())
        assert next(basis_line) == "BASIS"
        basis = []
        for length in basis_line:
            basis.append(tuple(int(next(basis_line)) for _ in range(int(length))))
        record = {"q": q, "s": s, "t": t, "basis": basis}
        for name in ["MATRIX", "STEPWISE", "SPECTRUM", "ACTION", "RESIDUAL"]:
            values = next(lines).split()
            assert values[0] == name
            record[name] = np.array([float(x) for x in values[1:]])
        record["MATRIX"] = record["MATRIX"].reshape(n, n)
        record["STEPWISE"] = record["STEPWISE"].reshape(n, n)
        yield record


def peak_memory(process):
    if sys.platform != "win32":
        return None  # use an isolated resource-measurement harness on Linux
    from ctypes import wintypes
    class Counters(ctypes.Structure):
        _fields_ = [("cb", wintypes.DWORD), ("faults", wintypes.DWORD)] + [
            (name, ctypes.c_size_t) for name in ["peak_ws", "ws", "peak_paged",
             "paged", "peak_nonpaged", "nonpaged", "pagefile", "peak_pagefile"]]
    counter = Counters()
    counter.cb = ctypes.sizeof(counter)
    query = ctypes.windll.psapi.GetProcessMemoryInfo
    query.argtypes = [wintypes.HANDLE, ctypes.c_void_p, wintypes.DWORD]
    if query(process._handle, ctypes.byref(counter), counter.cb):
        return counter.peak_ws
    return None


def benchmark(worker, output):
    rows = []
    for label, count, steps in [("path-many-starts", 48, 24),
                                ("path-one-start", 48, 1),
                                ("path-same-stage", 48, 1),
                                ("isolated", 48, 12)]:
        cells = {(i,): i for i in range(count)}
        if label != "isolated":
            cells.update({(i, i + 1): count for i in range(count - 1)})
        starts = sorted(set([i * (count - 1) // steps for i in range(1, steps + 1)]))
        if steps == 1:
            starts = [count // 2]
        if label == "path-same-stage":
            starts = [count]
        data = protocol(cells, 0, starts, count)
        (output / f"benchmark-{label}.txt").write_text(data, encoding="utf-8")
        # Correctness is checked on the very same workload before timing it.
        verified = subprocess.run([str(worker), "--oracle"], input=data,
                                  capture_output=True, text=True, check=True)
        (output / f"benchmark-{label}-verified.txt").write_text(verified.stdout, encoding="utf-8")
        import numpy as np
        for record in parse_worker(verified.stdout):
            _, _, projection, beta, _, _, _, _ = oracle(cells, record["s"], count, 0)
            np.testing.assert_allclose(record["MATRIX"], projection, atol=2e-9, rtol=1e-9)
            np.testing.assert_allclose(record["STEPWISE"], projection, atol=2e-9, rtol=1e-9)
            assert int(np.count_nonzero(np.abs(record["SPECTRUM"]) < 1e-8)) == beta
        for trial in range(3):
            routes = ["direct", "stepwise"] if trial % 2 == 0 else ["stepwise", "direct"]
            for route in routes:
                begin = time.perf_counter()
                child = subprocess.Popen([str(worker), "--benchmark", route],
                                         stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                         stderr=subprocess.PIPE, text=True)
                stdout, stderr = child.communicate(data)
                wall = time.perf_counter() - begin
                if child.returncode:
                    raise RuntimeError(stderr)
                timings = stdout.split()
                rows.append({"fixture": label, "route": route, "trial": trial,
                             "kernel_seconds": float(timings[1]), "process_seconds": wall,
                             "peak_working_set_bytes": peak_memory(child),
                             "checksum": float(timings[3]), "stderr": stderr})
    return rows


def compare_matlab_example(directory, worker, output):
    """Replay author-provided data/formulas, not execute MATLAB or Dionysus.

    Fixed Git blobs establish input identity. Original lower-boundary columns
    recover triangle labels/orientation so comparison is in matched coordinates,
    not merely eigenvalue equality. The original q=2 example is 35 points,
    s=1, t=1.5. GUDHI checks the rebuilt finite Rips skeleton separately.
    """
    import numpy as np
    import scipy
    from scipy import io, linalg
    import gudhi
    blobs = {
        "SchurComp.m": "96ed4f2e67f5c6e01f846e85b495bcb3dc08beaa",
        "persistLap.m": "82558d8e00587a4eb2bec4fbdc1563a2dd85e941",
        "boundary.py": "123dd341f018aa11a5200253e4a740cccf66f336",
        "VR_persistLap.m": "3e1d28c43a893421d7d2d59ff903edf5c84d9849",
        "Point.txt": "c1ad11a0dbca8d47e541331915b3110264c9ee00",
        "q1Boundary.mat": "17962a4735c28328d7d28432afd91aca3917cea3",
        "qBoundary.mat": "0554a062b5996cbbc52436ed265e004589f6ca1d",
    }
    hashes = {}
    for name, blob in blobs.items():
        raw = (directory / name).read_bytes()
        actual = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        if actual != blob:
            raise ValueError(f"fixed MATLAB reference blob mismatch: {name}")
        hashes[name] = hashlib.sha256(raw).hexdigest()
    points = np.loadtxt(directory / "Point.txt")
    distance = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=2)
    cells = {}
    for length in range(1, 5):
        for simplex in itertools.combinations(range(len(points)), length):
            value = max((distance[a, b] for a, b in itertools.combinations(simplex, 2)), default=0.0)
            if value <= 1.5:
                cells[simplex] = value
    bq = io.loadmat(directory / "qBoundary.mat")["Bq"].astype(float)
    bup = io.loadmat(directory / "q1Boundary.mat")["Bq1"].astype(float)
    keep = bq.shape[1]
    assert (bq.shape, bup.shape) == ((203, 580), (1606, 5013))
    edges = sorted(ordered(cells, 1.0, 1), key=lambda s: distance[s])
    original_basis = []
    for column in range(keep):
        rows = np.flatnonzero(bq[:, column])
        triangle = tuple(sorted({v for row in rows for v in edges[row]}))
        assert len(triangle) == 3
        assert set(edges[row] for row in rows) == set(itertools.combinations(triangle, 2))
        original_basis.append(triangle)
    assert len(set(original_basis)) == keep
    np.testing.assert_array_equal(boundary(edges, original_basis), bq)
    assert set(original_basis) == set(ordered(cells, 1.0, 2))
    # Check integer boundary-square zero on columns entirely contained in K.
    internal = np.all(bup[keep:] == 0, axis=0)
    np.testing.assert_array_equal(bq @ bup[:keep, internal], 0)
    print("Author sample: verified fixed blobs and 580 oriented K-triangle labels", flush=True)
    data = protocol(cells, 2, [1.0], 1.5)
    (output / "author-fixture.txt").write_text(data, encoding="utf-8")
    run = subprocess.run([str(worker), "--oracle"], input=data, text=True, capture_output=True)
    (output / "author-worker-stdout.txt").write_text(run.stdout, encoding="utf-8")
    (output / "author-worker-stderr.txt").write_text(run.stderr, encoding="utf-8")
    run.check_returncode()
    record, = list(parse_worker(run.stdout))
    gram = bup @ bup.T
    c = gram[keep:, keep:]
    e = gram[:keep, keep:]
    solution, _, rank, _ = linalg.lstsq(c, e.T, cond=1e-12, lapack_driver="gelsd")
    reference = gram[:keep, :keep] - e @ solution + bq.T @ bq
    permutation = [original_basis.index(simplex) for simplex in record["basis"]]
    reference = reference[np.ix_(permutation, permutation)]
    for name in ["MATRIX", "STEPWISE"]:
        np.testing.assert_allclose(record[name], reference, atol=2e-9, rtol=1e-9)
    eigenvalues = linalg.eigvalsh(reference)
    np.testing.assert_allclose(record["SPECTRUM"], eigenvalues, atol=2e-9, rtol=1e-9)
    st = gudhi.SimplexTree()
    for simplex, value in cells.items():
        st.insert(simplex, filtration=value)
    assert st.num_simplices() == len(cells)
    st.compute_persistence(homology_coeff_field=11, persistence_dim_max=True)
    count = sum(b <= 1.0 and d > 1.5 for b, d in st.persistence_intervals_in_dimension(2))
    assert count == 0
    assert np.count_nonzero(np.abs(eigenvalues) < 1e-8) == 0
    return {"scope": "fixed author data and least-norm formula replay; MATLAB not executed",
            "commit": "ebff6e3ad6f5eb90eb8058e4f5567755b62f0c3e", "sha256": hashes,
            "point_count": len(points), "degree": 2, "scales": [1.0, 1.5],
            "boundary_shapes": [list(bq.shape), list(bup.shape)], "schur_rank": int(rank),
            "direct_matrix_max_error": float(np.max(np.abs(record["MATRIX"] - reference))),
            "stepwise_matrix_max_error": float(np.max(np.abs(record["STEPWISE"] - reference))),
            "reference_solve_residual": float(np.max(np.abs(c @ solution - e.T)) / np.max(np.abs(e))),
            "rust_residuals": record["RESIDUAL"].tolist(),
            "smallest_ten_eigenvalues": eigenvalues[:10].tolist(),
            "gudhi_persistent_beta2_field11": int(count), "scipy": scipy.__version__}


def main():
    import numpy as np
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gudhi", action="store_true")
    parser.add_argument("--reference-python", type=Path)
    parser.add_argument("--reference-matlab", type=Path,
                        help="Fixed author source/data directory; needs SciPy and GUDHI")
    parser.add_argument("--benchmark", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    source_worker = args.worker.resolve()
    # Freeze the executable used for this run. In addition to reproducibility,
    # this keeps Windows builds from trying to overwrite an executing example.
    worker = (args.output / source_worker.name).resolve()
    shutil.copyfile(source_worker, worker)
    fixture_set = fixtures()
    groups = []
    for label, cells in fixture_set:
        scales = sorted(set([-1, *cells.values()]))
        top = max((len(s) - 1 for s in cells), default=0)
        for q in range(top + 2):
            for t in scales:
                groups.append((label, cells, q, [s for s in scales if s <= t], t))
    data = "".join(protocol(cells, q, starts, t) for _, cells, q, starts, t in groups)
    (args.output / "fixtures.txt").write_text(data, encoding="utf-8")
    run = subprocess.run([str(worker), "--oracle"], input=data, capture_output=True, text=True)
    (args.output / "worker-stdout.txt").write_text(run.stdout, encoding="utf-8")
    (args.output / "worker-stderr.txt").write_text(run.stderr, encoding="utf-8")
    run.check_returncode()
    external = None
    if args.reference_python:
        source = args.reference_python.read_text(encoding="utf-8")
        tree = ast.parse(source)
        selected = [n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == "schur_complement"]
        assert len(selected) == 1
        namespace = {"np": np}
        exec(compile(ast.Module(body=selected, type_ignores=[]), str(args.reference_python), "exec"), namespace)
        external = namespace["schur_complement"]
    records = iter(parse_worker(run.stdout))
    checks = []
    for label, cells, q, starts, t in groups:
        for s in starts:
            record = next(records)
            assert (record["q"], record["s"], record["t"]) == (q, s, t)
            basis, direct, projection, beta, exact, up, retained, lower = oracle(cells, s, t, q)
            assert record["basis"] == basis
            n = len(basis)
            errors = {}
            for name, matrix in [("schur", direct), ("nullspace", projection),
                                 ("stepwise", record["STEPWISE"])]:
                errors[name] = float(np.max(np.abs(record["MATRIX"] - matrix), initial=0))
                np.testing.assert_allclose(record["MATRIX"], matrix, atol=2e-9, rtol=1e-9)
            if exact is not None:
                rational = np.array([[float(x) for x in row] for row in exact]).reshape(n, n)
                np.testing.assert_allclose(record["MATRIX"], rational, atol=2e-9, rtol=1e-9)
            if external and 0 < n < len(up):
                permutation = retained + [i for i in range(len(up)) if i not in retained]
                gram = (up @ up.T)[np.ix_(permutation, permutation)]
                reference = external(gram, n, len(up)) + lower.T @ lower
                np.testing.assert_allclose(record["MATRIX"], reference, atol=2e-9, rtol=1e-9)
                errors["fixed_python_reference"] = float(np.max(np.abs(reference-direct), initial=0))
            eigenvalues = np.linalg.eigvalsh(projection)
            np.testing.assert_allclose(record["SPECTRUM"], eigenvalues, atol=2e-9, rtol=1e-9)
            assert int(np.count_nonzero(np.abs(eigenvalues) < 1e-8)) == beta
            assert np.min(eigenvalues, initial=0) >= -2e-9
            np.testing.assert_allclose(record["ACTION"], projection @ ((np.arange(n)+1)/7), atol=2e-9, rtol=1e-9)
            checks.append({"fixture": label, "degree": q, "s": s, "t": t,
                           "basis": basis, "matrix": record["MATRIX"].tolist(),
                           "spectrum": eigenvalues.tolist(), "real_persistent_rank": beta,
                           "exact_matrix": [[str(x) for x in row] for row in exact] if exact is not None else None,
                           "errors": errors, "residuals": record["RESIDUAL"].tolist()})
    assert next(records, None) is None
    topology = []
    if args.gudhi:
        import gudhi
        for label, cells in fixture_set[:4] + [fixture_set[-1]]:
            st = gudhi.SimplexTree()
            for simplex, value in cells.items():
                st.insert(simplex, filtration=value)
            st.compute_persistence(homology_coeff_field=11, min_persistence=-1,
                                   persistence_dim_max=True)
            for check in (c for c in checks if c["fixture"] == label):
                intervals = st.persistence_intervals_in_dimension(check["degree"])
                count = sum(b <= check["s"] and d > check["t"] for b, d in intervals)
                assert count == check["real_persistent_rank"]
            topology.append({"fixture": label, "version": gudhi.__version__,
                             "field": 11, "scope": "torsion-free topology corroboration only"})
        rp2 = dict(fixture_set)["rp2-field-counterexample"]
        field_counts = {}
        for field in [2, 11]:
            st = gudhi.SimplexTree()
            for simplex, value in rp2.items():
                st.insert(simplex, filtration=value)
            st.compute_persistence(homology_coeff_field=field, persistence_dim_max=True)
            field_counts[field] = st.betti_numbers()[1]
        assert field_counts == {2: 1, 11: 0}
        real = next(c for c in checks if c["fixture"] == "rp2-field-counterexample"
                    and c["degree"] == 1 and c["s"] == 2 and c["t"] == 2)
        assert real["real_persistent_rank"] == 0
        topology.append({"fixture": "rp2-field-counterexample", "version": gudhi.__version__,
                         "betti1_by_field": field_counts, "exact_real_rank": 0,
                         "scope": "counterexample to using F2 as real truth"})
    timings = benchmark(worker, args.output) if args.benchmark else []
    author_sample = (compare_matlab_example(args.reference_matlab, worker, args.output)
                     if args.reference_matlab else None)
    result = {"python": sys.version, "numpy": np.__version__, "worker": str(worker),
              "source_worker": str(source_worker),
              "worker_sha256": hashlib.sha256(worker.read_bytes()).hexdigest(),
              "pairs": len(checks), "exact_matrices": sum(c["exact_matrix"] is not None for c in checks),
              "max_matrix_error": max((max(c["errors"].values()) for c in checks), default=0),
              "checks": checks, "gudhi": topology, "benchmarks": timings,
              "matlab_author_sample": author_sample,
              "limitations": ["Numerical evidence is not an exact Betti certificate.",
                  "Peak working set includes process baseline, not allocation peak.",
                  "Bounded dense reproduction; no sparse scalability or production fast-path adoption."]}
    (args.output / "results.json").write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: result[k] for k in ["pairs", "exact_matrices", "max_matrix_error", "gudhi", "matlab_author_sample"]}))


if __name__ == "__main__":
    main()
