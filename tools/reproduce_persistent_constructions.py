#!/usr/bin/env python3
"""SPSC-A: paper Algorithms 1/2 and undamped LSMR, outside the Rust runtime.

NumPy/SciPy are optional reference dependencies. Inputs and all intermediate
results are frozen in a new output directory; no expected values come from Rust.
"""
from __future__ import annotations

import argparse
from fractions import Fraction
import hashlib
import itertools
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys

import compare_persistent_laplacians as reference


MATRIX_TOLERANCE = 1e-10
RESIDUAL_TOLERANCE = 1e-9
MATLAB_REVISION = "ebff6e3ad6f5eb90eb8058e4f5567755b62f0c3e"
PAPERS = {
    "persistent-laplacians-v3.pdf": (
        "https://arxiv.org/pdf/2012.02808v3",
        "a4fff3f81ca41d78a43d3ac88e718540b21337e97f026a324f2942bf9751aebd"),
    "lsmr.pdf": (
        "https://web.stanford.edu/group/SOL/software/lsmr/LSMR-SISC-2011.pdf",
        "d456214a23c51a8506a5e4f4cbc935a1f05ea0f2559755eba9847927e9974b82"),
}


def column_kernel(matrix):
    """Algorithm 1, lines 9-14: exact column reduction with its transform Y.

    Lowest nonzero row is the pivot. This differs from the independent oracle's
    row RREF/free-variable kernel. Fractions make pivot decisions exact here.
    """
    import numpy as np
    rows, columns = matrix.shape
    reduced = [[Fraction(int(x)) for x in matrix[:, j]] for j in range(columns)]
    transform = [[Fraction(i == j) for i in range(columns)] for j in range(columns)]
    pivots = {}
    zero = []
    for j in range(columns):
        while True:
            low = next((i for i in range(rows - 1, -1, -1) if reduced[j][i]), None)
            if low is None:
                zero.append(j)
                break
            if low not in pivots:
                pivots[low] = j
                break
            previous = pivots[low]
            factor = reduced[j][low] / reduced[previous][low]
            reduced[j] = [x - factor * y for x, y in zip(reduced[j], reduced[previous])]
            transform[j] = [x - factor * y for x, y in zip(transform[j], transform[previous])]
    for j in zero:
        assert all(sum(Fraction(int(matrix[i, k])) * transform[j][k]
                       for k in range(columns)) == 0 for i in range(rows))
    z_exact = [transform[j] for j in zero]
    z = np.array(z_exact, dtype=float).reshape(len(zero), columns).T
    return z, len(pivots), {"Y_columns": transform, "R_columns": reduced,
                            "zero_columns": zero, "Z_columns_exact": z_exact}


def constructions(b1, b2, a):
    """Identity-weight translations of MMW Algorithms 1 and 2."""
    import numpy as np
    z, rank, reduction = column_kernel(b2)
    gram = z.T @ z
    t = b1 @ z
    down = a.T @ a
    up = t @ np.linalg.solve(gram, t.T) if z.shape[1] else np.zeros((len(b1), len(b1)))
    # Algorithm 2: retain all upper columns, form and split the original Gram.
    b = np.vstack((b1, b2))
    g = b @ b.T
    n = len(b1)
    g11, g12, g21, g22 = g[:n, :n], g[:n, n:], g[n:, :n], g[n:, n:]
    solution = np.linalg.lstsq(g22, g21, rcond=1e-12)[0]
    schur = g11 - g12 @ solution
    return up + down, schur + down, {
        **reduction, "Z": z, "induced_gram": gram, "T": t,
        "G11": g11, "G12": g12, "G21": g21, "G22": g22,
        "minimum_norm_solution": solution, "rank_B2_exact": rank,
        "U": up, "D": down, "Delta": up + down,
        "schur_delta": schur + down,
    }


def lsmr(matrix, rhs, tolerance=1e-13, max_iterations=None):
    """Independent transcription of Fong--Saunders Section 2.8, steps 1-6.

    No damping or reorthogonalization. Fresh residuals supply a project stopping
    rule (not the paper's Section 3 estimators). An exhausted limit is an error.
    """
    import numpy as np
    m = np.asarray(matrix, dtype=float)
    y = np.asarray(rhs, dtype=float)
    if m.ndim != 2 or y.shape != (m.shape[0],):
        raise ValueError("LSMR shape mismatch")
    if not np.isfinite(m).all() or not np.isfinite(y).all():
        raise ValueError("LSMR requires finite data")
    if not 0 < tolerance < 1:
        raise ValueError("invalid LSMR tolerance")
    norm = np.linalg.norm
    beta = norm(y)
    u = y / beta if beta else np.zeros_like(y)
    v = m.T @ u
    alpha = norm(v)
    v = v / alpha if alpha else np.zeros_like(v)
    x = np.zeros(m.shape[1])
    h, hbar = v.copy(), np.zeros_like(v)
    alphabar, zetabar = alpha, alpha * beta
    rho, rhobar, cbar, sbar = 1., 1., 1., 0.
    trace = []
    normal_scale = norm(m, "fro") * beta
    if alpha == 0 or beta == 0:
        return x, trace
    limit = max_iterations if max_iterations is not None else 4 * min(m.shape) + 4
    for iteration in range(limit):
        old_u, old_v, old_alpha = u.copy(), v.copy(), alpha
        raw_u = m @ v - alpha * u
        beta = norm(raw_u)
        u = raw_u / beta if beta else np.zeros_like(raw_u)
        raw_v = m.T @ u - beta * v
        alpha = norm(raw_v)
        v = raw_v / alpha if alpha else np.zeros_like(raw_v)
        old_rho, old_rhobar = rho, rhobar
        rho = float(np.hypot(alphabar, beta))
        if rho == 0:
            raise ArithmeticError("LSMR zero rotation before convergence")
        c, s = alphabar / rho, beta / rho
        theta, alphabar = s * alpha, c * alpha
        thetabar = sbar * rho
        rhotemp = cbar * rho
        rhobar = float(np.hypot(rhotemp, theta))
        if rhobar == 0:
            raise ArithmeticError("LSMR zero secondary rotation")
        cbar, sbar = rhotemp / rhobar, theta / rhobar
        zeta, zetabar = cbar * zetabar, -sbar * zetabar
        hbar = h - (thetabar * rho / (old_rho * old_rhobar)) * hbar
        x = x + (zeta / (rho * rhobar)) * hbar
        h = v - (theta / rho) * h
        residual = y - m @ x
        normal = m.T @ residual
        trace.append({"iteration": iteration + 1, "alpha": alpha, "beta": beta,
                      "u": u.copy(), "v": v.copy(), "x": x.copy(),
                      "rho": rho, "rhobar": rhobar, "zetabar": zetabar,
                      "residual": residual, "normal_residual": normal,
                      "bidiagonal_forward_error": norm(m @ old_v - old_alpha * old_u - beta * u),
                      "bidiagonal_transpose_error": norm(m.T @ u - beta * old_v - alpha * v)})
        if not np.isfinite(x).all():
            raise ArithmeticError("nonfinite LSMR iterate")
        if norm(normal) <= tolerance * normal_scale:
            return x, trace
    raise ArithmeticError(f"LSMR iteration limit {limit}")


def row_projector(b2, tolerance=1e-12):
    """Project reformulation: fixed row-space basis, two MGS passes.

    The Rust action prototype uses this policy, not adaptive inner LSMR. A
    nonzero residual below the rank cutoff is ambiguous and rejected.
    """
    import numpy as np
    basis = []
    for row in b2.astype(float):
        v = row.copy()
        for _ in range(2):
            for q in basis:
                v -= np.dot(q, v) * q
        size = np.linalg.norm(v)
        if size <= tolerance * np.linalg.norm(row):
            if size > 64 * np.finfo(float).eps * max(1, np.linalg.norm(row)):
                raise ArithmeticError("ambiguous row-space rank")
            continue
        basis.append(v / size)
    return np.array(basis).reshape(len(basis), b2.shape[1])


def controls():
    cases = reference.fixtures()
    shell = {s: len(s) - 1 for length in range(1, 4)
             for s in itertools.combinations(range(4), length)}
    return cases + [("tetrahedron-shell", shell)]


def blocks(cells, s, t, q):
    import numpy as np
    k, l = reference.ordered(cells, s, q), reference.ordered(cells, t, q)
    upper = reference.ordered(cells, t, q + 1)
    lower = reference.ordered(cells, s, q - 1)
    b = reference.boundary(l, upper)
    keep = [l.index(cell) for cell in k]
    outside = [i for i in range(len(l)) if i not in keep]
    return b[keep], b[outside], reference.boundary(lower, k), {
        "K_basis": k, "L_basis": l, "upper_basis": upper, "lower_basis": lower,
        "retained_rows": keep, "outside_rows": outside,
        "B": np.vstack((b[keep], b[outside])),
    }


def serializable(value):
    import numpy as np
    if isinstance(value, Fraction):
        return str(value)
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, dict):
        return {key: serializable(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [serializable(item) for item in value]
    if isinstance(value, np.generic):
        return value.item()
    return value


def write_json(path, value):
    path.write_text(json.dumps(serializable(value), indent=2, allow_nan=False) + "\n",
                    encoding="utf-8")


def check_close(candidate, expected, label):
    import numpy as np
    error = float(np.linalg.norm(candidate - expected))
    if error > MATRIX_TOLERANCE * max(1., float(np.linalg.norm(expected))):
        raise AssertionError(f"{label}: Frobenius error {error}")
    return error


def check_pair(cells, s, t, q):
    import numpy as np
    from scipy.sparse.linalg import lsmr as scipy_lsmr
    b1, b2, a, bases = blocks(cells, s, t, q)
    basis, schur, detail = constructions(b1, b2, a)
    exact = reference.exact_projection(b1, b2, a)
    expected = np.array(exact, dtype=float).reshape(basis.shape)
    # Independent SVD projector and rational RREF ranks are NOT column reduction.
    _, singular, vh = np.linalg.svd(b2, full_matrices=True)
    rank = len(reference.rref(b2.tolist(), b2.shape[1])[1])
    q_kernel = vh[rank:].T
    projection = b1 @ q_kernel @ q_kernel.T @ b1.T + a.T @ a
    errors = [check_close(matrix, expected, label) for matrix, label in
              [(basis, "Algorithm 1"), (schur, "Algorithm 2"), (projection, "SVD")]]
    assert detail["rank_B2_exact"] == rank
    np.testing.assert_allclose(b2 @ detail["Z"], 0., atol=1e-14, rtol=0)
    rows = row_projector(b2)
    assert len(rows) == rank
    x = (np.arange(len(b1)) + 1.) / 7.
    v = np.cos(np.arange(len(b1)) + .3)
    y = b1.T @ x
    w, trace = lsmr(b2.T, y)
    z = y - b2.T @ w
    svd_z = q_kernel @ (q_kernel.T @ y)
    scipy_w = scipy_lsmr(b2.T, y, atol=1e-14, btol=1e-14,
                         maxiter=4 * min(b2.T.shape) + 4)[0]
    check_close(z, svd_z, "LSMR projection")
    check_close(w, np.linalg.lstsq(b2.T, y, rcond=1e-12)[0], "LSMR minimum norm")
    check_close(b2.T @ scipy_w, b2.T @ w, "SciPy LSMR fit")
    def fixed_action(vector):
        rhs = b1.T @ vector
        return b1 @ (rhs - rows.T @ (rows @ rhs)) + a.T @ (a @ vector)
    def inner_action(vector, tolerance):
        rhs = b1.T @ vector
        solution, _ = lsmr(b2.T, rhs, tolerance=tolerance)
        return b1 @ (rhs - b2.T @ solution) + a.T @ (a @ vector)
    action = inner_action(x, 1e-11)
    tighter = inner_action(x, 1e-14)
    checks = {
        "fixed_action_error": check_close(fixed_action(x), expected @ x, "fixed action"),
        "inner_action_error": check_close(action, expected @ x, "inexact action"),
        "fresh_tighter_error": check_close(tighter, expected @ x, "tighter action"),
        "repeat_error": check_close(inner_action(x, 1e-11), action, "repeat"),
        "linearity_error": check_close(inner_action(2*x-v, 1e-11),
                                       2*action-inner_action(v, 1e-11), "diagnostic linearity"),
        "symmetry_error": abs(float(v @ action - x @ inner_action(v, 1e-11))),
        "quadratic_form": float(x @ action),
        "projector_idempotence_error": check_close(z-b2.T@lsmr(b2.T, z)[0], z, "idempotence"),
        "kernel_residual": float(np.linalg.norm(b2 @ z)),
    }
    if checks["quadratic_form"] < -MATRIX_TOLERANCE or checks["symmetry_error"] > MATRIX_TOLERANCE:
        raise AssertionError(checks)
    eigenvalues, vectors = np.linalg.eigh(expected)
    fresh_residuals = []
    for value, vector in zip(eigenvalues, vectors.T):
        # Original exact-oracle action supplies a fresh residual, never Ritz symmetrization.
        residual = np.linalg.norm(expected @ vector - value * vector)
        normalized = residual / max(1., np.linalg.norm(expected, 2) + abs(value))
        if normalized > RESIDUAL_TOLERANCE:
            raise AssertionError("fresh eigenpair residual")
        fresh_residuals.append(float(normalized))
    return {**bases, **detail, "B1": b1, "B2": b2, "A": a,
            "exact_delta": exact, "singular_B2": singular,
            "row_basis": rows, "rank_A_exact": len(reference.rref(a.tolist(), a.shape[1])[1]),
            "eigenvalues": eigenvalues, "eigenvectors": vectors,
            "fresh_residuals": fresh_residuals, "matrix_errors": errors,
            "frozen_x": x, "frozen_v": v, "lsmr_trace": trace,
            "lsmr_solution": w, "lsmr_projection": z, "composition_diagnostics": checks}


def solver_controls():
    import numpy as np
    from scipy.sparse.linalg import lsmr as scipy_lsmr
    matrices = [np.array([[1., 0.], [0., 2.], [0., 0.]]),
                np.array([[1., 2.], [2., 4.], [0., 0.]]),
                np.zeros((3, 2)), np.eye(2), np.zeros((0, 2)), np.zeros((2, 0))]
    records = []
    for matrix in matrices:
        for y in [matrix @ np.ones(matrix.shape[1]), np.arange(matrix.shape[0]) + 1.]:
            x, trace = lsmr(matrix, y)
            expected = np.linalg.lstsq(matrix, y, rcond=1e-12)[0]
            check_close(x, expected, "solver control")
            check_close(x, scipy_lsmr(matrix, y, atol=1e-14, btol=1e-14)[0], "SciPy control")
            for step in trace:
                check_close(step["x"], scipy_lsmr(matrix, y, atol=0, btol=0, conlim=0,
                                                  maxiter=step["iteration"])[0], "SciPy iterate")
                assert max(step["bidiagonal_forward_error"], step["bidiagonal_transpose_error"]) < 1e-12
            records.append({"M": matrix, "y": y, "x": x, "expected": expected, "trace": trace})
    return records


def orientation_controls():
    import numpy as np
    rng = np.random.default_rng(144)
    records = []
    for label, cells in controls():
        for q in range(3):
            b1, b2, a, _ = blocks(cells, 0, max(cells.values(), default=0), q)
            baseline, _, _ = constructions(b1, b2, a)
            indices = rng.permutation(len(b1))
            signs = rng.choice([-1, 1], size=len(b1))
            columns = rng.permutation(b1.shape[1])
            column_signs = rng.choice([-1, 1], size=b1.shape[1])
            new_b1 = b1[indices][:, columns] * signs[:, None] * column_signs
            new_b2 = b2[:, columns] * column_signs
            new_a = a[:, indices] * signs
            expected = baseline[np.ix_(indices, indices)] * np.outer(signs, signs)
            changed, schur, _ = constructions(new_b1, new_b2, new_a)
            check_close(changed, expected, "signed/permuted basis construction")
            check_close(schur, expected, "signed/permuted Schur construction")
            records.append({"fixture": label, "q": q, "row_permutation": indices,
                            "row_signs": signs, "column_permutation": columns,
                            "column_signs": column_signs, "Delta": changed})
    return records


def run(output, sources, worker, gudhi):
    import numpy as np
    import scipy
    # A new directory is mandatory, including for failed runs.
    output.mkdir(parents=True, exist_ok=False)
    versions = {"python": sys.version, "numpy": np.__version__, "scipy": scipy.__version__}
    identities = {}
    for name, (url, digest) in PAPERS.items():
        data = (sources / name).read_bytes()
        actual = hashlib.sha256(data).hexdigest()
        if actual != digest:
            raise ValueError(f"paper hash mismatch: {name}: {actual}")
        identities[name] = {"url": url, "sha256": actual}
        shutil.copyfile(sources / name, output / name)
    matlab_blobs = {"persistLap.m": "82558d8e00587a4eb2bec4fbdc1563a2dd85e941",
                    "SchurComp.m": "96ed4f2e67f5c6e01f846e85b495bcb3dc08beaa"}
    for name, blob in matlab_blobs.items():
        raw = (sources / name).read_bytes()
        actual = hashlib.sha1(b"blob " + str(len(raw)).encode() + b"\0" + raw).hexdigest()
        if actual != blob:
            raise ValueError(f"MATLAB blob mismatch: {name}")
        identities[name] = {"revision": MATLAB_REVISION, "git_blob": actual,
                            "sha256": hashlib.sha256(raw).hexdigest()}
        shutil.copyfile(sources / name, output / name)
    # Freeze inputs BEFORE any candidate execution.
    cases = controls()
    write_json(output / "inputs.json", [{"name": label, "cells": list(cells.items())}
                                         for label, cells in cases])
    write_json(output / "policy.json", {"matrix_frobenius_relative": MATRIX_TOLERANCE,
               "fresh_residual_normalization": "norm(Delta*v-lambda*v)/max(1,norm(Delta,2)+abs(lambda))",
               "fresh_residual_tolerance": RESIDUAL_TOLERANCE,
               "rank": "exact Fraction on integer incidence; numerical cutoff 1e-12",
               "inner_stop": "norm(M.T*r) <= tolerance*norm(M,Fro)*norm(y)",
               "scope": "unit-scale finite project controls; no ill-conditioned acceptance"})
    shutil.copyfile(__file__, output / Path(__file__).name)
    shutil.copyfile(reference.__file__, output / Path(reference.__file__).name)
    write_json(output / "environment.json", {"versions": versions, "platform": platform.platform(),
               "processor": platform.processor(), "sources": identities,
               "git_head": subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip(),
               "git_status": subprocess.check_output(["git", "status", "--short"], text=True),
               "command": sys.argv, "author_runtime": {
                   "matlab": shutil.which("matlab"), "octave": shutil.which("octave"),
                   "status": "not run: no author runtime invocation; Python formula translation only"},
               "PETLS": "not run; optional external comparator, not an oracle"})
    rows, failures = [], []
    maximum_matrix, maximum_residual = 0., 0.
    protocols = []
    for label, cells in cases:
        scales = sorted({-1, *cells.values()})
        for q in range(3):
            for t in scales:
                protocols.append(reference.protocol(cells, q, [s for s in scales if s <= t], t))
                for s in scales:
                    if s > t:
                        continue
                    identity = {"fixture": label, "q": q, "s": s, "t": t}
                    try:
                        detail = check_pair(cells, s, t, q)
                        path = f"pair-{len(rows):04}.json"
                        write_json(output / path, {**identity, **detail})
                        rows.append({**identity, "artifact": path})
                        maximum_matrix = max(maximum_matrix, *detail["matrix_errors"])
                        maximum_residual = max(maximum_residual, *detail["fresh_residuals"], 0.)
                    except Exception as error:
                        failure = {**identity, "error": repr(error)}
                        failures.append(failure)
                        write_json(output / "failures.json", failures)
    # Analytic controls in explicitly declared lexicographic bases.
    p1 = np.array([[-1, 0], [0, 1]])
    p2 = np.array([[1, -1]])
    path, _, _ = constructions(p1, p2, np.zeros((0, 2), dtype=int))
    check_close(path, .5 * np.array([[1, -1], [-1, 1]]), "analytic path")
    if np.allclose(path, p1 @ p1.T):
        raise AssertionError("row-deletion control failed")
    triangle_b = np.array([[1], [-1], [1]])
    triangle_a = reference.boundary([(0,), (1,), (2,)], [(0, 1), (0, 2), (1, 2)])
    triangle, _, _ = constructions(triangle_b, np.zeros((0, 1), dtype=int), triangle_a)
    check_close(triangle, 3 * np.eye(3), "analytic triangle multiplicity")
    shell = reference.boundary(list(itertools.combinations(range(4), 2)),
                               list(itertools.combinations(range(4), 3)))
    signs = np.array([1, -1, 1, -1])
    np.testing.assert_array_equal(shell.T @ shell, 4 * np.eye(4) - np.outer(signs, signs))
    write_json(output / "solver-controls.json", solver_controls())
    write_json(output / "orientation-controls.json", orientation_controls())
    data = "".join(protocols)
    (output / "rust-input.txt").write_text(data, encoding="utf-8")
    if worker:
        frozen = output / ("worker.exe" if sys.platform == "win32" else "worker")
        shutil.copy2(worker, frozen)
        process = subprocess.run([str(frozen.resolve()), "--oracle"], input=data,
                                 capture_output=True, text=True)
        (output / "worker-stdout.txt").write_text(process.stdout, encoding="utf-8")
        (output / "worker-stderr.txt").write_text(process.stderr, encoding="utf-8")
        process.check_returncode()
        records = list(reference.parse_worker(process.stdout))
        if len(records) != len(rows):
            raise AssertionError("Rust pair count mismatch")
        for record, row in zip(records, rows):
            detail = json.loads((output / row["artifact"]).read_text())
            expected = np.array([[float(Fraction(x)) for x in line]
                                 for line in detail["exact_delta"]]).reshape(record["MATRIX"].shape)
            check_close(record["MATRIX"], expected, "existing Rust Schur")
    topology = None
    if gudhi:
        import gudhi as gd
        versions["gudhi"] = gd.__version__
        topology = []
        # Compare persistent ranks only on these declared torsion-free controls.
        supported = {"path", "triangle", "singular", "zero", "empty", "tetrahedron", "tetrahedron-shell"}
        for label, cells in cases:
            tree = gd.SimplexTree()
            for simplex, birth in cells.items():
                tree.insert(simplex, filtration=birth)
            fields = [2, 11] if label == "rp2-field-counterexample" else [11]
            for field in fields:
                tree.compute_persistence(homology_coeff_field=field, persistence_dim_max=True)
                checks = []
                if label in supported:
                    for row in rows:
                        if row["fixture"] != label:
                            continue
                        detail = json.loads((output / row["artifact"]).read_text())
                        real = (len(detail["K_basis"]) - detail["rank_A_exact"]
                                - len(reference.rref(detail["B"], len(detail["upper_basis"]))[1])
                                + detail["rank_B2_exact"])
                        count = sum(birth <= row["s"] and death > row["t"] for birth, death
                                    in tree.persistence_intervals_in_dimension(row["q"]))
                        assert count == real
                        checks.append({"q": row["q"], "s": row["s"], "t": row["t"], "rank": real})
                if label == "rp2-field-counterexample":
                    betti = tree.betti_numbers()
                    assert betti[1] == (1 if field == 2 else 0)
                topology.append({"fixture": label, "field": field, "betti": tree.betti_numbers(),
                                 "persistent_rank_checks": checks,
                                 "scope": "torsion-free corroboration" if label in supported
                                 else "field control only; no real-spectrum claim"})
    write_json(output / "topology.json", topology)
    write_json(output / "index.json", rows)
    summary = {"pairs_passed": len(rows), "failures": failures,
               "maximum_matrix_frobenius_error": maximum_matrix,
               "maximum_fresh_eigenpair_residual": maximum_residual,
               "solver_controls": 12, "rust_schur_compared": bool(worker),
               "gudhi": topology is not None,
               "limitations": ["no MATLAB/Octave execution", "no PETLS execution",
                   "LSMR composition diagnostics do not certify a fixed self-adjoint action",
                   "no ill-conditioned rank guarantee or performance claim"]}
    write_json(output / "summary.json", summary)
    manifest = {path.name: hashlib.sha256(path.read_bytes()).hexdigest()
                for path in sorted(output.iterdir()) if path.is_file()}
    write_json(output / "manifest.json", manifest)
    print(json.dumps(summary, indent=2))
    if failures:
        raise AssertionError(f"{len(failures)} retained failures")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--sources", type=Path, required=True)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--gudhi", action="store_true")
    args = parser.parse_args()
    run(args.output, args.sources, args.worker, args.gudhi)


if __name__ == "__main__":
    main()
