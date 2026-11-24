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


def rust_call(worker, path, header, matrix, initial=None):
    np = numpy()
    values = np.asarray(matrix).ravel().tolist()
    if initial is not None:
        values += np.asarray(initial).ravel().tolist()
    path.write_text(header + "\n" + " ".join(format(x, ".17g") for x in values) + "\n")
    process = subprocess.run([str(worker), str(path)], capture_output=True, text=True)
    path.with_suffix(".stdout").write_text(process.stdout)
    path.with_suffix(".stderr").write_text(process.stderr)
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
    if args.gudhi:
        import gudhi
        manifest["gudhi"] = gudhi.__version__
    (args.output / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    result = {"polynomials": polynomial_controls(), "metric": model_metric_control(),
              "sequences": [sequence_control(False), sequence_control(True)]}
    (args.output / "controls.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    if args.worker:
        result["rust_parity"] = rust_parity(args.worker, args.output / "rust", result["sequences"])
    result["changing_basis"] = changing_basis_control(args.worker, args.output / "changing-basis", args.gudhi)
    for family in result["sequences"]:
        inputs = json.dumps({"matrices": family["matrices"], "fields": family["fields"]},
                            sort_keys=True).encode()
        family["input_sha256"] = hashlib.sha256(inputs).hexdigest()
    (args.output / "results.json").write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(json.dumps({family["family"]: {
        route: {"actions": sum(r["actions"] for r in data["output"]),
                "iterations": sum(len(r["trace"]) for r in data["output"]),
                "max_residual": max(r["verification"]["normalized_residual"]
                                    for r in data["output"])}
        for route, data in family["routes"].items()} for family in result["sequences"]}, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
