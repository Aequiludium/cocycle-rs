"""Controlled SCSF/ChFSI reproduction; optional NumPy, outside the Rust runtime.

This implements the defining-polynomial correction, not a silent transcription
of the printed algorithms. See docs/research/sequence-subspaces.md for scope.
"""

import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import traceback


def numpy():
    import numpy as np
    return np


def chebyshev_filter(matrix, block, degrees, lower, upper, anchor, printed=False):
    """Apply T_m((A-cI)/e)/T_m((anchor-c)/e), one degree per column.

    printed=True retains the current-iterate subtraction printed in both papers,
    with their ambiguous column split removed solely for the polynomial control.
    Degree zero is identity. The scaling anchor must lie below the stop interval.
    """
    np = numpy()
    a, y = np.asarray(matrix, dtype=float), np.asarray(block, dtype=float)
    m = np.asarray(degrees)
    if (a.ndim != 2 or a.shape[0] != a.shape[1] or y.ndim != 2
            or y.shape[0] != a.shape[0] or m.shape != (y.shape[1],)
            or not np.issubdtype(m.dtype, np.integer) or np.any(m < 0)):
        raise ValueError("filter shapes/degrees")
    if not all(np.all(np.isfinite(x)) for x in (a, y)):
        raise ValueError("nonfinite filter input")
    if (not all(math.isfinite(x) for x in (lower, upper, anchor))
            or not anchor < lower < upper):
        raise ValueError("require finite anchor < lower < upper")
    center, radius = (lower + upper) / 2, (upper - lower) / 2
    h = (a - center * np.eye(len(a))) / radius
    sigma1 = radius / (anchor - center)
    previous, current = y.copy(), y.copy()
    active = m > 0
    current[:, active] = sigma1 * (h @ y[:, active])
    sigma = sigma1
    for degree in range(1, int(m.max(initial=0))):
        next_sigma = 1 / (2 / sigma1 - sigma)
        active = m > degree
        result = current.copy()
        subtract = current if printed else previous
        result[:, active] = (2 * next_sigma * (h @ current[:, active])
                             - next_sigma * sigma * subtract[:, active])
        previous, current, sigma = current, result, next_sigma
    if not np.all(np.isfinite(current)):
        raise ValueError("nonfinite filter output")
    return current


def fft_order(fields, truncation=2):
    """Unshifted, unnormalized FFT2 upper-left block; start 0, ties by index."""
    np = numpy()
    p = np.asarray(fields, dtype=float)
    if (p.ndim != 3 or p.shape[1] != p.shape[2] or not 0 < truncation <= p.shape[1]
            or not np.all(np.isfinite(p))):
        raise ValueError("FFT field shape/truncation/finiteness")
    descriptors = np.fft.fft2(p, axes=(-2, -1))[:, :truncation, :truncation]
    if not np.all(np.isfinite(descriptors)):
        raise ValueError("nonfinite FFT descriptor")
    order = [0] if len(p) else []
    unused = set(range(1, len(p)))
    while unused:
        selected = min(unused, key=lambda j: (
            float(np.linalg.norm(descriptors[order[-1]] - descriptors[j])), j))
        order.append(selected)
        unused.remove(selected)
    return order, descriptors


def chfsi(matrix, initial, count, lower, upper, anchor, degree=8,
          tolerance=1e-11, max_iterations=200):
    """Fixed-degree real-symmetric control with explicit active/locked QR/Ritz.

    Known control bounds replace the papers' estimated bounds/degree optimizer.
    Only a converged leading prefix is locked. Every acceptance uses A anew.
    The trace retains actual filtered/QR/Ritz/locked blocks and fresh residuals.
    """
    np = numpy()
    a, initial = np.asarray(matrix, dtype=float), np.asarray(initial, dtype=float)
    n = len(a)
    if (a.shape != (n, n) or initial.ndim != 2 or initial.shape[0] != n
            or not 0 < count <= initial.shape[1] <= n
            or not np.all(np.isfinite(a)) or not np.all(np.isfinite(initial))
            or not np.allclose(a, a.T, rtol=0, atol=1e-14)
            or degree <= 0 or tolerance <= 0 or max_iterations <= 0):
        raise ValueError("ChFSI input")
    if np.linalg.matrix_rank(initial, tol=1e-13) != initial.shape[1]:
        raise ValueError("initial rank")
    active = np.linalg.qr(initial, mode="reduced")[0]
    locked = np.empty((n, 0))
    locked_values = []
    trace = []
    actions = 0
    for iteration in range(max_iterations):
        width = active.shape[1]
        filtered = chebyshev_filter(a, active, [degree] * width, lower, upper, anchor)
        raw_filtered = filtered.copy()
        actions += degree * width
        # Normalize columns before the rank/QR check to avoid scale dominance.
        lengths = np.linalg.norm(filtered, axis=0)
        if np.any(lengths == 0):
            raise ValueError("filtered rank")
        filtered /= lengths
        joined = np.column_stack((locked, filtered))
        if np.linalg.svd(joined, compute_uv=False)[-1] < 1e-13:
            raise ValueError("filtered rank")
        q = np.linalg.qr(joined, mode="reduced")[0][:, locked.shape[1]:]
        image = a @ q
        actions += width
        quotient = q.T @ image
        values, rotation = np.linalg.eigh((quotient + quotient.T) / 2)
        ritz = q @ rotation
        fresh_image = a @ ritz
        actions += width
        errors = np.linalg.norm(fresh_image - ritz * values, axis=0)
        converged = errors <= tolerance * np.maximum(1., np.abs(values))
        take = 0
        while take < min(width, count - len(locked_values)) and converged[take]:
            take += 1
        trace.append({
            "iteration": iteration, "locked_shape": list(locked.shape),
            "active_shape": list(active.shape), "qr_shape": list(q.shape),
            "ritz_shape": list(quotient.shape), "degree": degree,
            "filter_interval": [lower, upper], "anchor": anchor,
            "filter_input": active.tolist(), "filter_output": raw_filtered.tolist(),
            "filtered": filtered.tolist(), "q": q.tolist(),
            "quotient": quotient.tolist(), "ritz": ritz.tolist(),
            "locked": locked.tolist(), "values": values.tolist(),
            "fresh_residuals": errors.tolist(), "newly_locked": take,
        })
        locked = np.column_stack((locked, ritz[:, :take]))
        locked_values.extend(values[:take].tolist())
        if len(locked_values) == count:
            # Retain lookahead as well as locked vectors for the next sample.
            retained = np.column_stack((locked, ritz[:, take:]))
            residual = np.linalg.norm(a @ locked - locked * locked_values, axis=0)
            actions += count
            if np.any(residual > tolerance * np.maximum(1., np.abs(locked_values))):
                raise ValueError("fresh final residual")
            return {"values": locked_values, "vectors": locked.tolist(),
                    "retained": retained.tolist(), "residuals": residual.tolist(),
                    "actions": actions, "trace": trace}
        active = ritz[:, take:]
    raise ValueError("ChFSI convergence")


def polynomial_controls():
    np = numpy()
    diagonal = np.array([0., 0.5, 1., 2., 3., 4.])
    a, y = np.diag(diagonal), np.arange(1., 19.).reshape(6, 3)
    lower, upper, anchor = 1., 4., 0.
    t, t0 = (diagonal - 2.5) / 1.5, -2.5 / 1.5
    results = []
    for degree in (0, 1, 2, 3, 8):
        # Independent NumPy scalar polynomial evaluation, not this recurrence.
        coefficients = [0.] * degree + [1.]
        expected = (np.polynomial.chebyshev.chebval(t, coefficients)
                    / np.polynomial.chebyshev.chebval(t0, coefficients))[:, None] * y
        corrected = chebyshev_filter(a, y, [degree] * 3, lower, upper, anchor)
        printed = chebyshev_filter(a, y, [degree] * 3, lower, upper, anchor, True)
        error = float(np.max(np.abs(expected - corrected)))
        if error > 1e-12:
            raise ValueError("defining polynomial disagreement")
        results.append({"degree": degree, "corrected_max_error": error,
                        "printed_max_error": float(np.max(np.abs(expected - printed))),
                        "expected": expected.tolist(), "corrected": corrected.tolist(),
                        "printed": printed.tolist()})
    mixed = chebyshev_filter(a, y, [0, 2, 8], lower, upper, anchor)
    for column, degree in enumerate([0, 2, 8]):
        np.testing.assert_allclose(mixed[:, column],
                                   np.array(results[[0, 1, 2, 3, 8].index(degree)]
                                            ["expected"])[:, column], atol=1e-12)
    return results


def controls(abrupt=False):
    np = numpy()
    rng = np.random.default_rng(148)
    n, count, block = 24, 3, 6
    spectrum = np.r_[0.4, 0.4, 1.1, np.linspace(3., 12., n - count)]
    base = np.linalg.qr(rng.normal(size=(n, n)))[0]
    samples = [0, 5, 1, 7, 2, 6, 3, 4]
    matrices, fields, vectors = [], [], []
    for sample in samples:
        theta = sample * 0.025
        q = base.copy()
        for i in range(block):
            c, s = math.cos(theta * (i + 1)), math.sin(theta * (i + 1))
            q[:, [i, i + block]] = q[:, [i, i + block]] @ np.array([[c, -s], [s, c]])
        if abrupt:
            # Descriptors deliberately cease predicting eigenspace similarity.
            q = np.linalg.qr(rng.normal(size=(n, n)))[0]
        matrices.append((q * spectrum) @ q.T)
        vectors.append(q[:, :count])
        x, y = np.meshgrid(np.arange(4), np.arange(4), indexing="ij")
        fields.append(1 + theta * (1 + np.cos(2 * np.pi * x / 4)
                                  + 0.5 * np.sin(2 * np.pi * y / 4)))
    initial = np.random.default_rng(73).normal(size=(n, block))
    return matrices, fields, vectors, spectrum, initial


def sequence_control(abrupt=False):
    np = numpy()
    started = time.perf_counter()
    matrices, fields, truth_vectors, spectrum, initial = controls(abrupt)
    construction = time.perf_counter() - started
    started = time.perf_counter()
    order, descriptors = fft_order(fields)
    ordering = time.perf_counter() - started
    routes = {}
    for route, indices, warm in (("cold_dense", range(len(matrices)), False),
                                 ("cold_chfsi", range(len(matrices)), False),
                                 ("warm_original", range(len(matrices)), True),
                                 ("fft_scsf", order, True)):
        previous, output, solve_seconds, verify_seconds = initial, {}, 0., 0.
        for index in indices:
            started = time.perf_counter()
            if route == "cold_dense":
                values, vectors = np.linalg.eigh(matrices[index])
                result = {"values": values[:3].tolist(), "vectors": vectors[:, :3].tolist(),
                          "actions": 0, "trace": []}
            else:
                result = chfsi(matrices[index], previous if warm else initial, 3,
                               2., 12.1, 0.3)
                previous = np.array(result["retained"])
            solve_seconds += time.perf_counter() - started
            started = time.perf_counter()
            values, vectors = np.array(result["values"]), np.array(result["vectors"])
            error = float(np.max(np.abs(values - spectrum[:3])))
            projector = float(np.linalg.norm(vectors @ vectors.T
                              - truth_vectors[index] @ truth_vectors[index].T))
            normalized = float(np.linalg.norm(matrices[index] @ vectors - vectors * values)
                               / max(1., np.linalg.norm(matrices[index], "fro")))
            orthogonality = float(np.linalg.norm(vectors.T @ vectors - np.eye(3)))
            if error > 1e-10 or normalized > 1e-9 or projector > 1e-8 or orthogonality > 1e-10:
                raise ValueError(f"{route} oracle mismatch at {index}")
            result["verification"] = {"eigenvalue_error": error, "projector_error": projector,
                                      "normalized_residual": normalized,
                                      "orthogonality_error": orthogonality}
            output[index] = result
            verify_seconds += time.perf_counter() - started
        routes[route] = {"solve_seconds": solve_seconds, "verify_seconds": verify_seconds,
                         "ordering_seconds": ordering if route == "fft_scsf" else 0.,
                         "output": [output[i] for i in range(len(matrices))]}
    return {"family": "abrupt" if abrupt else "smooth", "order": order,
            "descriptors_real": descriptors.real.tolist(),
            "descriptors_imag": descriptors.imag.tolist(),
            "construction_seconds": construction, "ordering_seconds": ordering,
            "matrices": [a.tolist() for a in matrices],
            "fields": [a.tolist() for a in fields], "routes": routes}


def model_metric_control():
    """Eq. (2), 2x2 interior Dirichlet grid, h=1, positive nonconstant field.

    A=D L, S=D^(1/2) L D^(1/2), u=D^(1/2)v. Eigenvalues are ordered
    by absolute value; this test does not substitute algebraically smallest.
    """
    np = numpy()
    stencil = np.array([[-4., 1., 1., 0.], [1., -4., 0., 1.],
                        [1., 0., -4., 1.], [0., 1., 1., -4.]])
    d = np.array([1., 2., 3., 4.])
    a = d[:, None] * stencil
    symmetric = np.sqrt(d)[:, None] * stencil * np.sqrt(d)[None, :]
    values, v = np.linalg.eigh(symmetric)
    order = np.argsort(np.abs(values), kind="stable")
    values, u = values[order], np.sqrt(d)[:, None] * v[:, order]
    residual = float(np.linalg.norm(a @ u - u * values))
    metric_error = float(np.linalg.norm(u.T @ (u / d[:, None]) - np.eye(4)))
    direct = np.linalg.eigvals(a)
    np.testing.assert_allclose(values, direct[np.argsort(np.abs(direct))], atol=1e-12)
    if residual > 1e-10 or metric_error > 1e-12:
        raise ValueError("paper-model metric control")
    return {"matrix": a.tolist(), "symmetric": symmetric.tolist(),
            "values_by_magnitude": values.tolist(), "original_residual": residual,
            "metric_orthogonality_error": metric_error,
            "euclidean_symmetry_defect": float(np.linalg.norm(a - a.T))}


def source_identity(root):
    paths = ["tools/reproduce_subspaces.py", "tools/reference/subspace_cocycle.rs",
             "tools/compare_persistent_laplacians.py", "src/eigensolver.rs",
             "src/execution/mod.rs", "src/lib.rs"]
    return {"commit": subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=root,
                                               text=True).strip(),
            "status": subprocess.check_output(["git", "status", "--porcelain"], cwd=root,
                                               text=True),
            "sha256": {p: hashlib.sha256((root / p).read_bytes()).hexdigest()
                       for p in paths if (root / p).exists()}}


def rust_call(worker, path, header, matrix, initial=None, expected_error=None):
    np = numpy()
    values = np.asarray(matrix).ravel().tolist()
    if initial is not None:
        values += np.asarray(initial).ravel().tolist()
    path.write_text(header + "\n" + " ".join(format(x, ".17g") for x in values) + "\n")
    process = subprocess.run([str(worker), str(path)], capture_output=True, text=True)
    path.with_suffix(".stdout").write_text(process.stdout)
    path.with_suffix(".stderr").write_text(process.stderr)
    if expected_error is not None:
        if not process.returncode or expected_error not in process.stderr:
            raise ValueError(f"Rust failure semantics differ: {path}")
        return {"returncode": process.returncode, "error": expected_error}
    if process.returncode:
        raise RuntimeError(f"Rust worker failed ({process.returncode}): {path}")
    return json.loads(process.stdout)


def rust_parity(worker, output, sequences):
    np = numpy()
    output.mkdir()
    maximum = 0.
    calls = 0
    # Replay the actual frozen inputs from *every* filter step in all routes.
    for family in sequences:
        for route, data in family["routes"].items():
            for index, sample in enumerate(data["output"]):
                a = np.array(family["matrices"][index])
                for step in sample["trace"]:
                    y = np.array(step["filter_input"])
                    lower, upper = step["filter_interval"]
                    header = (f"filter {len(a)} {y.shape[1]} {lower} {upper} {step['anchor']} "
                              + " ".join([str(step["degree"])] * y.shape[1]))
                    actual = rust_call(worker, output / f"filter-{calls}.txt", header, a, y)
                    actual = np.array(actual["filtered"]).reshape(y.shape)
                    expected = np.array(step["filter_output"])
                    error = float(np.linalg.norm(actual - expected) / max(1., np.linalg.norm(expected)))
                    if error > 1e-10:
                        raise ValueError("Rust frozen-filter parity")
                    maximum = max(maximum, error)
                    calls += 1
    # Include every preregistered degree, identity and mixed column degrees.
    a, y = np.diag([0., 0.5, 1., 2., 3., 4.]), np.arange(1., 19.).reshape(6, 3)
    for degrees in [[d] * 3 for d in [0, 1, 2, 3, 8]] + [[8, 0, 2]]:
        header = "filter 6 3 1 4 0 " + " ".join(map(str, degrees))
        actual = rust_call(worker, output / f"filter-{calls}.txt", header, a, y)
        expected = chebyshev_filter(a, y, degrees, 1., 4., 0.)
        np.testing.assert_allclose(np.array(actual["filtered"]).reshape(y.shape), expected,
                                   atol=1e-12, rtol=1e-12)
        calls += 1
    rows = []
    for abrupt in (False, True):
        matrices, _, _, spectrum, initial = controls(abrupt)
        previous = initial
        for i, a in enumerate(matrices):
            for warm in (0, 1):
                result = rust_call(worker, output / f"lobpcg-{abrupt}-{i}-{warm}.txt",
                                   f"lobpcg 24 6 3 {warm}", a, previous if warm else None)
                np.testing.assert_allclose(result["values"], spectrum[:3], atol=1e-10, rtol=0)
                vectors = np.array(result["vectors"]).reshape(24, 3)
                values, truth = np.linalg.eigh(a)
                projector = float(np.linalg.norm(vectors @ vectors.T - truth[:, :3] @ truth[:, :3].T))
                if projector > 1e-8:
                    raise ValueError("Rust repeated-cluster projector")
                if warm:
                    # Fill lookahead deterministically; only the preceding vectors
                    # are retained, not its values, residuals or whole history.
                    previous = np.column_stack((vectors, initial[:, 3:]))
                rows.append({"family": "abrupt" if abrupt else "smooth", "index": i,
                             "warm": bool(warm), "projector_error": projector,
                             "actions": result["actions"], "iterations": result["iterations"],
                             "residuals": result["residuals"]})
    return {"filter_replays": calls, "max_filter_relative_error": maximum, "lobpcg": rows}


def rust_sequence_parity(worker, output, sequences):
    """Rewrite admission: immutable inputs and invariant intermediate comparisons.

    QR/Ritz signs and bases within repeated clusters need not be byte-identical.
    Compare active projectors, lifted quotients, ordered values, orthogonality
    and independently recomputed residuals. Action/locking counts are diagnostics.
    """
    np = numpy()
    output.mkdir()
    step_count = 0
    max_projector, max_lifted, max_values = 0., 0., 0.
    rows = []
    for family in sequences:
        fields = np.array(family["fields"])
        samples, p, _ = fields.shape
        sorted_result = rust_call(worker, output / f"order-{family['family']}.txt",
                                  f"order {samples} {p} 2", fields)
        descriptors = np.array(family["descriptors_real"]) + 1j * np.array(family["descriptors_imag"])
        actual = (np.array(sorted_result["descriptors_real"])
                  + 1j * np.array(sorted_result["descriptors_imag"])).reshape(descriptors.shape)
        np.testing.assert_allclose(actual, descriptors, atol=1e-10, rtol=0)
        if sorted_result["order"] != family["order"]:
            raise ValueError("Rust ordering permutation differs")
        matrices = np.array(family["matrices"])
        first = family["routes"]["cold_chfsi"]["output"][0]["trace"][0]
        initial = np.array(first["filter_input"])
        n, block = initial.shape
        for route, code in (("cold_chfsi", 0), ("warm_original", 1), ("fft_scsf", 2)):
            result = rust_call(worker, output / f"sequence-{family['family']}-{code}.txt",
                               f"sequence {n} {block} 3 8 200 1e-11 2 12.1 0.3 {samples} {p} 2 {code}",
                               np.r_[fields.ravel(), matrices.ravel(), initial.ravel()])
            if result["order"] != (family["order"] if code == 2 else list(range(samples))):
                raise ValueError("Rust sequence order")
            if len(result["output"]) != samples:
                raise ValueError("Rust restored sample count")
            for index, solved in enumerate(result["output"]):
                expected = family["routes"][route]["output"][index]
                vectors = np.array(solved["vectors"]).reshape(n, 3)
                reference = np.array(expected["vectors"])
                values = np.array(solved["values"])
                np.testing.assert_allclose(values, expected["values"], atol=1e-10, rtol=0)
                projector = float(np.linalg.norm(vectors @ vectors.T - reference @ reference.T))
                residual = float(np.linalg.norm(matrices[index] @ vectors - vectors * values)
                                 / max(1., np.linalg.norm(matrices[index], "fro")))
                if (projector > 1e-8 or residual > 1e-9
                        or np.linalg.norm(vectors.T @ vectors - np.eye(3)) > 1e-10):
                    raise ValueError("Rust end-to-end subspace mismatch")
                retained = np.array(solved["retained"]).reshape(n, block)
                if np.linalg.norm(retained.T @ retained - np.eye(block)) > 1e-10:
                    raise ValueError("Rust retained lookahead orthogonality")
                rows.append({"family": family["family"], "route": route, "index": index,
                             "projector_error": projector, "normalized_residual": residual,
                             "rust_actions": solved["actions"], "python_actions": expected["actions"],
                             "rust_iterations": len(solved["trace"]),
                             "python_iterations": len(expected["trace"])})
        for data in family["routes"].values():
            for index, sample in enumerate(data["output"]):
                a = matrices[index]
                for step in sample["trace"]:
                    active, locked = np.array(step["filter_input"]), np.array(step["locked"])
                    l, width = step["locked_shape"][1], step["active_shape"][1]
                    lower, upper = step["filter_interval"]
                    result = rust_call(worker, output / f"step-{step_count}.txt",
                                       f"step {n} {l + width} {l} {step['degree']} 200 1e-11 {lower} {upper} {step['anchor']}",
                                       np.r_[a.ravel(), active.ravel(), locked.ravel()])
                    for key in ("locked_shape", "active_shape", "qr_shape", "ritz_shape"):
                        if result[key] != step[key]:
                            raise ValueError("Rust frozen active/locked dimensions")
                    q, reference_q = np.array(result["q"]).reshape(n, width), np.array(step["q"])
                    quotient = np.array(result["quotient"]).reshape(width, width)
                    projector = float(np.linalg.norm(q @ q.T - reference_q @ reference_q.T))
                    lifted = float(np.linalg.norm(q @ quotient @ q.T
                                   - reference_q @ np.array(step["quotient"]) @ reference_q.T)
                                   / max(1., np.linalg.norm(a, "fro")))
                    eigenvalue_error = float(np.max(np.abs(np.array(result["values"]) - step["values"])))
                    ritz = np.array(result["ritz"]).reshape(n, width)
                    fresh = np.linalg.norm(a @ ritz - ritz * result["values"], axis=0)
                    np.testing.assert_allclose(fresh, result["fresh_residuals"], atol=1e-12, rtol=1e-9)
                    if (projector > 1e-9 or lifted > 1e-10 or eigenvalue_error > 1e-10
                            or np.linalg.norm(ritz.T @ ritz - np.eye(width)) > 1e-10):
                        raise ValueError("Rust frozen QR/Ritz parity")
                    max_projector = max(max_projector, projector)
                    max_lifted = max(max_lifted, lifted)
                    max_values = max(max_values, eigenvalue_error)
                    step_count += 1
    return {"qr_ritz_replays": step_count, "max_active_projector_error": max_projector,
            "max_lifted_quotient_error": max_lifted, "max_ritz_value_error": max_values,
            "solves": rows}


def rust_sequence_failures(worker, output, sequences):
    """Replay first-round rank/convergence controls plus explicit input exclusions."""
    np = numpy()
    output.mkdir()
    a = np.array(sequences[0]["matrices"][0])
    initial = np.array(sequences[0]["routes"]["cold_chfsi"]["output"][0]["trace"][0]["filter_input"])
    n, b = initial.shape
    rows = []
    for tag, matrix, start, iterations, degree, expected in (
            ("initial-rank", a, np.ones_like(initial), 200, 8, "initial rank"),
            ("nonconvergence", a, initial, 1, 8, "ChFSI convergence"),
            ("nonfinite", a, np.full_like(initial, np.nan), 200, 8, "ChFSI input"),
            ("zero-filter", 7.05 * np.eye(n), initial, 200, 1, "filtered rank")):
        try:
            chfsi(matrix, start, 3, 2., 12.1, 0.3, degree=degree, max_iterations=iterations)
        except ValueError as error:
            if expected not in str(error):
                raise
        else:
            raise ValueError("Python failure control unexpectedly succeeded")
        result = rust_call(worker, output / (tag + ".txt"),
                           f"chfsi {n} {b} 3 {degree} {iterations} 1e-11 2 12.1 0.3",
                           matrix, start, expected_error=expected)
        rows.append({"case": tag, **result})
    for tag, fields in (("ties-and-large-distance", np.array([np.full((4, 4), v)
                                                             for v in [0., 10000., 10000., 30000.]])),
                        ("empty-order", np.empty((0, 4, 4)))):
        result = rust_call(worker, output / (tag + ".txt"), f"order {len(fields)} 4 2", fields)
        if result["order"] != fft_order(fields)[0]:
            raise ValueError("Rust ordering edge control")
        rows.append({"case": tag, "order": result["order"]})
    return rows


def rust_transfer_call(worker, path, block, old, new, width, max_entries, projector=None,
                       expected_error=None):
    """Keep simplex identities as integer text, never round them through f64."""
    np = numpy()
    header = (f"transfer {old['source'].encode().hex() or '-'} {new['source'].encode().hex() or '-'} "
              f"{old['degree']} {new['degree']} {old['metric']} {new['metric']} "
              f"{len(old['basis'])} {len(new['basis'])} {block.shape[1]} {width} {max_entries} "
              f"{int(projector is not None)}")
    for context in (old, new):
        for simplex, sign in context["basis"]:
            header += " " + " ".join(map(str, (sign, *simplex)))
    # Freeze the first prototype's PCG64-generated candidates as worker inputs;
    # the Rust mathematical rewrite does not claim NumPy RNG implementation parity.
    fresh = np.random.default_rng(73).normal(size=(width, len(new["basis"])))
    payload = np.r_[block.ravel(), fresh.ravel(),
                    np.asarray(projector).ravel() if projector is not None else []]
    return rust_call(worker, path, header, payload, expected_error=expected_error)


def rust_basis_parity(worker, output, frozen_path):
    """Replay original pair matrices and the saved preceding Rust solved vectors."""
    np = numpy()
    output.mkdir()
    frozen = json.loads(frozen_path.read_text())
    previous, old, group = None, None, None
    rows, raw_hashes = [], {}
    for pair in frozen["changing_basis"]["pairs"]:
        current_group = (pair["source"], pair["degree"], pair["fixed"])
        if current_group != group:
            previous, old, group = None, None, current_group
        context = {"source": pair["source"], "degree": pair["degree"], "metric": "euclidean",
                   "basis": [(tuple(s), 1) for s in pair["basis"]]}
        n = len(context["basis"])
        a = np.array(pair["matrix"]).reshape(n, n)
        width = min(3, n)
        tag = f"{pair['source']}-{pair['degree']}-{pair['fixed']}-{pair['start']}-{pair['target']}"
        if previous is not None:
            values, vectors = np.linalg.eigh(a)
            positive = vectors[:, values > 1e-10]
            for projected in (False, True):
                wanted = min(width, positive.shape[1]) if projected else width
                projector = positive @ positive.T if projected else None
                expected, rank = transfer_block(previous, old, context, wanted, 4096, projector)
                result = rust_transfer_call(worker, output / f"{tag}-{projected}.txt",
                                            previous, old, context, wanted, 4096, projector)
                mapped = np.array(result["mapped"]).reshape(n, wanted)
                error = float(np.linalg.norm(mapped @ mapped.T - expected @ expected.T))
                if result["rank"] != rank or error > 1e-10:
                    raise ValueError("Rust frozen basis-transfer parity")
                row = {"source": pair["source"], "degree": pair["degree"],
                       "fixed": pair["fixed"], "start": pair["start"], "target": pair["target"],
                       "projected": projected, "width": wanted, "rank": rank,
                       "projector_error": error, "retained_bytes": n * wanted * 8}
                if not projected:
                    solved = rust_call(worker, output / f"{tag}-solve.txt",
                                       f"lobpcg {n} {width} {width} 1", a, mapped)
                    v = np.array(solved["vectors"]).reshape(n, width)
                    np.testing.assert_allclose(solved["values"], pair["spectrum"][:width], atol=1e-10, rtol=0)
                    row["new_operator_residual"] = float(np.linalg.norm(a @ v - v * solved["values"]))
                    if row["new_operator_residual"] > 1e-9:
                        raise ValueError("Rust transferred new-operator residual")
                rows.append(row)
        raw = frozen_path.parent / "changing-basis" / (tag + "-1.stdout")
        raw_hashes[raw.name] = hashlib.sha256(raw.read_bytes()).hexdigest()
        saved = json.loads(raw.read_text())
        previous = np.array(saved["vectors"]).reshape(n, len(saved["values"]))
        old = context
    old = {"source": "triangle", "degree": 1, "metric": "euclidean",
           "basis": [((0, 1), 1), ((0, 2), 1), ((1, 2), 1)]}
    new = {**old, "basis": [((1, 2), -1), ((0, 1), 1)]}
    block = np.array([[1., 0.], [0., 1.], [2., 0.]])
    edges = []
    for tag, target, width, cap, p, expected_error in (
            ("orientation", new, 1, 20, None, None),
            ("source", {**new, "source": "another-owner"}, 1, 20, None, "incompatible"),
            ("degree", {**new, "degree": 0}, 1, 20, None, "incompatible"),
            ("metric", {**new, "metric": "weighted"}, 1, 20, None, "incompatible"),
            ("invalid-sign", {**new, "basis": [((0, 1), 0)]}, 1, 20, None, "invalid oriented"),
            ("storage", new, 1, 1, None, "storage"),
            ("projected-rank", new, 2, 20, np.diag([1., 0.]), "rank insufficient"),
            ("zero-width", new, 0, 20, None, None)):
        if expected_error:
            try:
                transfer_block(block, old, target, width, cap, p)
            except ValueError:
                pass
            else:
                raise ValueError("Python transfer failure unexpectedly succeeded")
        result = rust_transfer_call(worker, output / (tag + ".txt"), block, old,
                                    target, width, cap, p, expected_error)
        if not expected_error:
            expected, rank = transfer_block(block, old, target, width, cap, p)
            np.testing.assert_allclose(np.array(result["mapped"]).reshape(len(target["basis"]), width),
                                       expected, atol=1e-12, rtol=0)
            if result["rank"] != rank:
                raise ValueError("Rust edge transfer rank")
        edges.append({"case": tag, "error": expected_error})
    large = (1 << 53) + 1
    old = {"source": "owner with space", "degree": 0, "metric": "euclidean",
           "basis": [((large,), 1)]}
    new = {**old, "basis": [((large + 1,), 1), ((large,), -1)]}
    result = rust_transfer_call(worker, output / "large-integer-identity.txt", np.array([[1.]]), old, new, 1, 20)
    np.testing.assert_allclose(result["mapped"], [0., -1.], atol=1e-12, rtol=0)
    edges.append({"case": "large-integer-identity", "error": None})
    return {"transfers": rows, "edge_controls": edges, "saved_vector_sha256": raw_hashes,
            "fresh_candidates": "explicit frozen NumPy PCG64(73) inputs, not a Rust RNG port"}


def transfer_block(block, old, new, width, max_entries, projector=None):
    """Project experiment: exact source-local simplex identity and sign mapping.

    Basis entries are (canonical sorted vertex tuple, orientation sign). Ordinary
    Euclidean real coordinates only; caller supplies the NEW positive-range
    projector when appropriate. The dense projector is an oracle experiment.
    """
    np = numpy()
    if (old["source"] != new["source"] or old["degree"] != new["degree"]
            or old["metric"] != "euclidean" or new["metric"] != "euclidean"):
        raise ValueError("incompatible source/degree/inner product")
    for context in (old, new):
        seen = set()
        for simplex, sign in context["basis"]:
            if (len(simplex) != context["degree"] + 1 or tuple(sorted(set(simplex))) != simplex
                    or simplex in seen or sign not in (-1, 1)):
                raise ValueError("invalid oriented simplex basis")
            seen.add(simplex)
    x = np.asarray(block, dtype=float)
    n = len(new["basis"])
    if (x.ndim != 2 or x.shape[0] != len(old["basis"]) or not np.all(np.isfinite(x))
            or not 0 <= width <= n or n * max(width, x.shape[1]) > max_entries):
        raise ValueError("transfer shape/finiteness/storage")
    if width == 0:
        return np.empty((n, 0)), 0
    lookup = {s: (i, sign) for i, (s, sign) in enumerate(old["basis"])}
    mapped = np.zeros((n, x.shape[1]))
    for j, (simplex, sign) in enumerate(new["basis"]):
        if simplex in lookup:
            i, old_sign = lookup[simplex]
            mapped[j] = sign * old_sign * x[i]
    if projector is not None:
        p = np.asarray(projector)
        if (p.shape != (n, n) or n * n > max_entries or not np.all(np.isfinite(p))
                or not np.allclose(p, p.T, rtol=0, atol=1e-10)
                or not np.allclose(p @ p, p, rtol=0, atol=1e-10)):
            raise ValueError("new range projector")
        mapped = p @ mapped
    basis = []
    retained_rank = 0
    # Keep only the previous block; fresh directions never create history.
    rng = np.random.default_rng(73)
    candidates = itertools.chain(mapped.T, rng.normal(size=(width, n)),
                                 (np.eye(1, n, i).ravel() for i in range(n)))
    for number, vector in enumerate(candidates):
        v = vector.copy()
        if projector is not None:
            v = p @ v
        length = np.linalg.norm(v)
        if length:
            v /= length
        for _ in range(2):
            for q in basis:
                v -= (q @ v) * q
        length = np.linalg.norm(v)
        if length > 1e-13:
            basis.append(v / length)
            if number < x.shape[1]:
                retained_rank += 1
        if len(basis) == width:
            break
    if len(basis) != width:
        raise ValueError("projected rank insufficient; caller must choose cold/full route")
    return np.column_stack(basis) if width else np.empty((n, 0)), retained_rank


def changing_basis_control(worker, output, use_gudhi=False):
    # Reuse the existing independent exact-rational oracle, never production
    # Schur/eigensolver code to generate expectations.
    import compare_persistent_laplacians as oracle
    np = numpy()
    output.mkdir()
    rows = []
    for name, cells in oracle.fixtures():
        if name not in ("path", "triangle", "tetrahedron"):
            continue
        stages = sorted(set(cells.values()))
        st = None
        if use_gudhi:
            import gudhi
            st = gudhi.SimplexTree()
            for simplex, value in cells.items():
                st.insert(simplex, filtration=value)
            st.persistence(homology_coeff_field=101, persistence_dim_max=True)
        for degree in (0, 1, 2):
            for fixed in ("start", "target"):
                previous, old = None, None
                for start, target in ((stages[0], t) for t in stages) if fixed == "start" else (
                        (s, stages[-1]) for s in stages):
                    basis = oracle.ordered(cells, start, degree)
                    all_rows = oracle.ordered(cells, target, degree)
                    outside = [s for s in all_rows if s not in basis]
                    columns = oracle.ordered(cells, target, degree + 1)
                    full = oracle.boundary(basis + outside, columns)
                    lower = oracle.boundary(oracle.ordered(cells, start, degree - 1), basis)
                    exact = oracle.exact_projection(full[:len(basis)], full[len(basis):], lower)
                    a = np.array(exact, dtype=float).reshape(len(basis), len(basis))
                    context = {"source": name, "degree": degree, "metric": "euclidean",
                               "basis": [(s, 1) for s in basis]}
                    if not basis:
                        continue
                    values, vectors = np.linalg.eigh(a)
                    nullity = int(np.sum(np.abs(values) < 1e-10))
                    if st is not None:
                        betti = st.persistent_betti_numbers(start, target)
                        if nullity != (betti[degree] if degree < len(betti) else 0):
                            raise ValueError("torsion-free GUDHI persistent nullity")
                    width = min(3, len(basis))
                    count = width
                    initial, rank = (None, 0) if previous is None else transfer_block(
                        previous, old, context, width, 4096)
                    positive = vectors[:, values > 1e-10]
                    positive_width = min(width, positive.shape[1])
                    positive_rank = 0
                    if previous is not None and positive_width:
                        projected, positive_rank = transfer_block(
                            previous, old, context, positive_width, 4096, positive @ positive.T)
                        if np.linalg.norm((np.eye(len(basis)) - positive @ positive.T) @ projected) > 1e-10:
                            raise ValueError("new positive-range projection")
                    pair = {"source": name, "degree": degree, "fixed": fixed,
                            "start": start, "target": target, "basis": basis,
                            "matrix": a.tolist(), "spectrum": values.tolist(),
                            "mapped_rank": rank, "positive_mapped_rank": positive_rank,
                            "retained_bytes": len(basis) * width * 8}
                    if worker:
                        for warm in (0, 1):
                            supplied = initial if warm and initial is not None else None
                            tag = f"{name}-{degree}-{fixed}-{start}-{target}-{warm}"
                            actual = rust_call(worker, output / (tag + ".txt"),
                                               f"lobpcg {len(basis)} {width} {count} {int(supplied is not None)}",
                                               a, supplied)
                            np.testing.assert_allclose(actual["values"], values[:count], atol=1e-10, rtol=0)
                            v = np.array(actual["vectors"]).reshape(len(basis), count)
                            error = float(np.linalg.norm(a @ v - v * actual["values"]))
                            if error > 1e-9:
                                raise ValueError("new-operator mapped residual")
                            if warm:
                                solved = v
                            pair["warm" if warm else "cold"] = {
                                "actions": actual["actions"], "iterations": actual["iterations"],
                                "fresh_residual": error}
                    # Transfer actual solved vectors when Rust is selected. The
                    # NumPy-only control transfers its independent dense solve.
                    previous, old = solved if worker else vectors[:, :width], context
                    rows.append(pair)
    return {"pairs": rows, "gudhi_field": 101 if use_gudhi else None,
            "scope": "tiny full-spectrum transfer; dense oracle positive projection separately"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--worker", type=Path, help="fresh standalone Rust parity worker")
    parser.add_argument("--gudhi", action="store_true", help="require torsion-free field-101 topology control")
    parser.add_argument("--papers", type=Path, help="directory containing scsf-v2.pdf and chfsi-v2.pdf")
    parser.add_argument("--rust-sequences", action="store_true", help="require second-round sequence rewrite parity")
    parser.add_argument("--frozen", type=Path, help="first-round results.json; reuse exact saved inputs and outputs")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    try:
        return reproduce(args)
    except Exception:
        (args.output / "failure.txt").write_text(traceback.format_exc())
        raise


def reproduce(args):
    root = Path(__file__).resolve().parents[1]
    np = numpy()
    manifest = {"source": source_identity(root), "python": sys.version,
                "numpy": np.__version__, "platform": platform.platform(),
                "github": {key: os.environ.get(key) for key in
                           ("GITHUB_SHA", "GITHUB_RUN_ID", "GITHUB_RUN_ATTEMPT", "GITHUB_REPOSITORY")},
                "threads": {key: os.environ.get(key) for key in
                            ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
                "command": sys.argv, "seed": 148,
                "scope": "corrected fixed-degree real-symmetric controls; not author datasets",
                "evidence": "numerical-evidence; diagnostic timings include trace export preparation"}
    if args.papers:
        manifest["papers_sha256"] = {name: hashlib.sha256((args.papers / name).read_bytes()).hexdigest()
                                     for name in ("scsf-v2.pdf", "chfsi-v2.pdf")}
        if manifest["papers_sha256"] != {
                "scsf-v2.pdf": "eee5de671463c0a98ef7887696cc8ae5051fe4c4798707c88d50aedb22d44617",
                "chfsi-v2.pdf": "464a08f4f812986d3eb40ad99028a4e083de39ceef559199883be781a7e845e2"}:
            raise ValueError("pinned paper PDF hashes differ")
    if args.worker:
        args.worker = args.worker.resolve(strict=True)
        manifest["worker_sha256"] = hashlib.sha256(args.worker.read_bytes()).hexdigest()
    if args.rust_sequences and not args.worker:
        raise ValueError("--rust-sequences requires --worker")
    if args.frozen:
        manifest["frozen_sha256"] = hashlib.sha256(args.frozen.read_bytes()).hexdigest()
    if args.gudhi:
        import gudhi
        manifest["gudhi"] = gudhi.__version__
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    if args.frozen:
        frozen = json.loads(args.frozen.read_text())
        for family in frozen["sequences"]:
            payload = json.dumps({"matrices": family["matrices"], "fields": family["fields"]}, sort_keys=True).encode()
            if hashlib.sha256(payload).hexdigest() != family["input_sha256"]:
                raise ValueError("frozen sequence input hash differs")
        result = {key: frozen[key] for key in ("polynomials", "metric", "sequences")}
    else:
        result = {"polynomials": polynomial_controls(), "metric": model_metric_control(),
                  "sequences": [sequence_control(False), sequence_control(True)]}
    (args.output / "controls.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if args.worker:
        result["rust_parity"] = rust_parity(args.worker, args.output / "rust", result["sequences"])
    if args.rust_sequences:
        result["rust_sequences"] = rust_sequence_parity(args.worker, args.output / "rust-sequences", result["sequences"])
        result["rust_sequence_failures"] = rust_sequence_failures(args.worker, args.output / "rust-sequence-failures", result["sequences"])
    result["changing_basis"] = changing_basis_control(args.worker, args.output / "changing-basis", args.gudhi)
    if args.rust_sequences:
        # A caller can supply the historical first-round bundle. CI first emits
        # the same prototype bundle, then freezes that bundle for the rewrite.
        frozen_path = args.frozen or args.output / "prototype-results.json"
        if not args.frozen:
            frozen_path.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
        result["rust_basis"] = rust_basis_parity(args.worker, args.output / "rust-basis", frozen_path)
    for family in result["sequences"]:
        inputs = json.dumps({"matrices": family["matrices"], "fields": family["fields"]},
                            sort_keys=True).encode()
        family["input_sha256"] = hashlib.sha256(inputs).hexdigest()
    (args.output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    summary = {family["family"]: {
        route: {"actions": sum(r["actions"] for r in data["output"]),
                "iterations": sum(len(r["trace"]) for r in data["output"]),
                "max_residual": max(r["verification"]["normalized_residual"]
                                    for r in data["output"])}
        for route, data in family["routes"].items()} for family in result["sequences"]}
    if args.rust_sequences:
        rewrite = result["rust_sequences"]
        summary["rust_rewrite"] = {
            "qr_ritz_replays": rewrite["qr_ritz_replays"],
            "sequence_solves": len(rewrite["solves"]),
            "max_active_projector_error": rewrite["max_active_projector_error"],
            "max_lifted_quotient_error": rewrite["max_lifted_quotient_error"],
            "basis_transfers": len(result["rust_basis"]["transfers"]),
            "sequence_edge_controls": len(result["rust_sequence_failures"]),
            "basis_edge_controls": len(result["rust_basis"]["edge_controls"])}
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
