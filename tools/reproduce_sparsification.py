"""Bounded ordinary-Hodge sampling and KID reproduction, outside the Rust core.

The printed Algorithm 4.1 and Section 4.3 recurrence are intentionally separate.
The analytic-bin integral is a project diagnostic, not the unspecified paper
histogramIntegral or a probabilistic guarantee. Dependencies are imported lazily.
"""

import argparse
from collections import Counter
from datetime import datetime, timezone
from fractions import Fraction
import hashlib
import importlib.metadata
import itertools
import json
import math
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import tracemalloc


PAPERS = {
    "2502.07558v2.pdf": "373640de5a96c8fce16784aecc072ef240cdacb96476fffde7ef3263426fa8f8",
    "1708.08436v3.pdf": "1470ef2d8785dc9872f525479253b00e6d7d84b6b5b08ef8804b981e0983305e",
}
TOLERANCE = 1e-10


def numpy():
    import numpy as np
    return np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n",
                          encoding="utf-8")


def rref(rows, columns):
    """Independent exact-rational row reduction; no numerical rank cutoff."""
    matrix = [[Fraction(value) for value in row] for row in rows]
    if any(len(row) != columns for row in matrix):
        raise ValueError("ragged exact matrix")
    pivots = []
    for column in range(columns):
        pivot = next((i for i in range(len(pivots), len(matrix))
                      if matrix[i][column]), None)
        if pivot is None:
            continue
        index = len(pivots)
        matrix[index], matrix[pivot] = matrix[pivot], matrix[index]
        scale = matrix[index][column]
        matrix[index] = [x / scale for x in matrix[index]]
        for i, row in enumerate(matrix):
            if i != index and row[column]:
                scale = row[column]
                matrix[i] = [x - scale * y for x, y in zip(row, matrix[index])]
        pivots.append(column)
    return matrix, pivots


def exact_kernel(rows, columns):
    reduced, pivots = rref(rows, columns)
    result = []
    for free in sorted(set(range(columns)) - set(pivots)):
        vector = [Fraction(0)] * columns
        vector[free] = Fraction(1)
        for index, pivot in enumerate(pivots):
            vector[pivot] = -reduced[index][free]
        result.append(vector)
    return len(pivots), result


def exact_scores(boundary):
    """diag of the rational row-space projector, independently of pinv(G)."""
    np = numpy()
    boundary = np.asarray(boundary)
    if (boundary.ndim != 2 or not np.isfinite(boundary).all()
            or np.any(boundary != np.floor(boundary))):
        raise ValueError("exact oracle requires a finite integer boundary")
    rows, columns = boundary.shape
    reduced, pivots = rref(boundary.astype(int).tolist(), columns)
    rank = len(pivots)
    if not rank:
        return [Fraction(0)] * columns, 0
    basis = reduced[:rank]
    gram = [[sum(x * y for x, y in zip(a, b)) for b in basis] for a in basis]
    augmented = [row + [Fraction(i == j) for j in range(rank)]
                 for i, row in enumerate(gram)]
    inverse, inverse_pivots = rref(augmented, 2 * rank)
    if inverse_pivots != list(range(rank)):
        raise ArithmeticError("row-space Gram is singular")
    inverse = [row[rank:] for row in inverse]
    scores = [sum(basis[i][j] * inverse[i][k] * basis[k][j]
                  for i in range(rank) for k in range(rank))
              for j in range(columns)]
    if any(x < 0 for x in scores) or sum(scores) != rank:
        raise ArithmeticError("exact leverage invariant failed")
    _, kernel = exact_kernel(boundary.astype(int).tolist(), columns)
    if any(sum(Fraction(int(boundary[i, j])) * vector[j]
               for j in range(columns)) != 0
           for vector in kernel for i in range(rows)):
        raise ArithmeticError("exact right kernel failed")
    return scores, rank


def boundary_matrix(lower, upper):
    np = numpy()
    matrix = np.zeros((len(lower), len(upper)), dtype=int)
    indices = {tuple(cell): index for index, cell in enumerate(lower)}
    for column, cell in enumerate(upper):
        for index in range(len(cell)):
            face = tuple(cell[:index] + cell[index + 1:])
            matrix[indices[face], column] = (-1)**index
    return matrix


def numerical_scores(boundary):
    np = numpy()
    gram = boundary @ boundary.T
    inverse = np.linalg.pinv(gram, rcond=1e-12, hermitian=True)
    return np.einsum("ij,ij->j", boundary, inverse @ boundary)


def sample_gram(boundary, probabilities, count, seed):
    """Keep original coordinates and aggregate Gram weights, not their roots."""
    np = numpy()
    boundary = np.asarray(boundary, dtype=float)
    probabilities = np.asarray(probabilities, dtype=float)
    if (boundary.ndim != 2 or not np.isfinite(boundary).all()
            or type(count) is not int or count < 1):
        raise ValueError("finite matrix and positive integer draw count required")
    rows, columns = boundary.shape
    if (probabilities.shape != (columns,) or not np.isfinite(probabilities).all()
            or (probabilities < 0).any()):
        raise ValueError("invalid probability vector")
    if not np.any(boundary):
        if np.any(probabilities):
            raise ValueError("zero-rank input requires zero probabilities")
        return np.zeros((rows, rows)), [], np.zeros(columns), []
    if (abs(float(probabilities.sum()) - 1.0) > 1e-12
            or np.any((np.linalg.norm(boundary, axis=0) > 0) & (probabilities == 0))):
        raise ValueError("probabilities must sum to one and cover every nonzero column")
    draws = np.random.Generator(np.random.PCG64(seed)).choice(
        columns, size=count, replace=True, p=probabilities)
    counts = np.bincount(draws, minlength=columns)
    weights = np.zeros(columns)
    selected = counts > 0
    weights[selected] = counts[selected] / (count * probabilities[selected])
    sampled = (boundary * weights) @ boundary.T
    # Independent per-draw construction checks aggregation and sqrt convention.
    factor = boundary[:, draws] / np.sqrt(count * probabilities[draws])
    if not all(np.isfinite(array).all() for array in (weights, sampled, factor)):
        raise ArithmeticError("nonfinite sampling result")
    if np.linalg.norm(sampled - factor @ factor.T) > TOLERANCE * max(1, np.linalg.norm(sampled)):
        raise ArithmeticError("sampling factor/weight disagreement")
    return sampled, draws.tolist(), weights, counts.tolist()


def literal_moments(up_laplacian, probes, degree):
    """PDF Algorithm 4.1 lines 3-11, including its indexing and indentation."""
    np = numpy()
    if up_laplacian.shape[0] != probes.shape[0]:
        raise ValueError("printed up-Laplacian/probe shape mismatch")
    moments = np.zeros((degree + 1, probes.shape[0]))
    first, second = probes, up_laplacian @ probes
    for m in range(3, degree + 1):
        current = 2 * up_laplacian @ second - first
        if m % 2:
            moments[m] = np.mean(current * probes, axis=1)
            first, second = second, current
    return moments


def kid_moments(dual_gram_scaled, probes, degree):
    """Section 4.3: update every degree, collect 2 diag(T_m(H)) at odd m."""
    np = numpy()
    h = np.asarray(dual_gram_scaled, dtype=float)
    probes = np.asarray(probes, dtype=float)
    if (h.ndim != 2 or h.shape[0] != h.shape[1] or not np.isfinite(h).all()
            or probes.ndim != 2 or probes.shape[0] != h.shape[0]
            or not probes.shape[1] or not np.isfinite(probes).all()
            or type(degree) is not int or degree < 1):
        raise ValueError("invalid moment matrix, probes or degree")
    moments = np.zeros((degree + 1, h.shape[0]))
    previous, current = probes, h @ probes
    for m in range(1, degree + 1):
        if m % 2:
            moments[m] = 2 * np.mean(probes * current, axis=1)
        previous, current = current, 2 * h @ current - previous
    return moments


def direct_moments(h, degree):
    """Independent scalar cos(m arccos(x)) via dense symmetric eigendecomposition."""
    np = numpy()
    values, vectors = np.linalg.eigh(h)
    if np.min(values, initial=0) < -1e-12 or np.max(values, initial=0) > 1 + 1e-12:
        raise ValueError("diagnostic spectrum outside [0,1]")
    # Retain unmodified values; arccos accepts one-ulp eigensolver overshoot only
    # through a separately recorded roundoff projection, never score clipping.
    angles = np.arccos(np.minimum(1, np.maximum(0, values)))
    result = np.zeros((degree + 1, h.shape[0]))
    for m in range(1, degree + 1, 2):
        result[m] = 2 * (vectors**2) @ np.cos(m * angles)
    return result


def analytic_bin_scores(moments, bins=64):
    """Project-only undamped Chebyshev bin masses on (0,1], width 1/bins.

    For odd m, the exact bin integral of 2 T_m(x)/(pi sqrt(1-x*x))
    is 2[sin(m acos(a))-sin(m acos(b))]/(pi*m). No midpoint rule,
    negative-score clipping, Jackson kernel or claimed paper parity is used.
    """
    np = numpy()
    if type(bins) is not int or bins < 1:
        raise ValueError("positive bin count required")
    edges = np.linspace(0, 1, bins + 1)
    masses = np.zeros((bins, moments.shape[1]))
    for m in range(1, moments.shape[0], 2):
        integrals = 2 * (np.sin(m * np.arccos(edges[:-1]))
                         - np.sin(m * np.arccos(edges[1:]))) / (math.pi * m)
        masses += integrals[:, None] * moments[m]
    return masses.sum(axis=0), masses


def comparison_setup(boundary, rank):
    np = numpy()
    gram = boundary @ boundary.T
    values, vectors = np.linalg.eigh(gram)
    if rank:
        active = vectors[:, -rank:]
        whitening = active / np.sqrt(values[-rank:])
    else:
        whitening = np.zeros((gram.shape[0], 0))
    exact_rank, rational_kernel = exact_kernel(boundary.T.astype(int).tolist(), gram.shape[0])
    if exact_rank != rank:
        raise ArithmeticError("left/right exact ranks disagree")
    kernel = np.array(rational_kernel, dtype=float).T.reshape(gram.shape[0], -1)
    if kernel.shape[1]:
        kernel = np.linalg.qr(kernel)[0]
    return gram, whitening, kernel


def assess(gram, sampled, whitening, kernel, rank):
    """Check the whole positive-space generalized spectrum and nullspace."""
    np = numpy()
    generalized = np.linalg.eigvalsh(whitening.T @ sampled @ whitening)
    kernel_residual = float(np.linalg.norm(sampled @ kernel) / max(1, np.linalg.norm(sampled)))
    sampled_values = np.linalg.eigvalsh(sampled)
    sampled_rank = int(np.count_nonzero(sampled_values > 1e-9 * max(1, np.linalg.norm(gram, 2))))
    error = float(np.max(np.abs(generalized - 1), initial=0))
    rank_loss = sampled_rank != rank
    return {
        "generalized_eigenvalues": generalized.tolist(),
        "minimum_generalized": float(generalized[0]) if rank else None,
        "maximum_generalized": float(generalized[-1]) if rank else None,
        "relative_sandwich_error": error,
        "kernel_residual": kernel_residual,
        "numeric_sampled_rank": sampled_rank,
        "rank_loss": rank_loss,
        "kernel_equal_numerically": not rank_loss and kernel_residual <= TOLERANCE,
        "psd_sandwich_pass": {str(eps): bool(error <= eps + TOLERANCE
                                            and not rank_loss and kernel_residual <= TOLERANCE)
                              for eps in (0.25, 0.5)},
        "frobenius_relative_error": float(np.linalg.norm(sampled - gram) / max(1, np.linalg.norm(gram))),
        "sampled_matrix_nnz": int(np.count_nonzero(sampled)),
    }


def controls():
    np = numpy()
    edges = list(itertools.combinations(range(4), 2))
    faces = list(itertools.combinations(range(4), 3))
    b = boundary_matrix(edges, faces)
    s = np.array([1, -1, 1, -1])
    if not np.array_equal(b.T @ b, 4 * np.eye(4) - np.outer(s, s)):
        raise ArithmeticError("tetrahedron analytic Gram failed")
    scores, rank = exact_scores(b)
    if rank != 3 or scores != [Fraction(3, 4)] * 4:
        raise ArithmeticError("tetrahedron exact scores failed")
    h = b.T @ b / 4
    z = np.array(list(itertools.product((-1, 1), repeat=4)), dtype=float).T
    moments = kid_moments(h, z, 8)
    for degree in (1, 3, 5, 7):
        if not np.allclose(moments[degree], 1.5, atol=1e-14, rtol=0):
            raise ArithmeticError("enumerated Rademacher moment failed")
    if not np.allclose(moments, direct_moments(h, 8), atol=1e-13, rtol=0):
        raise ArithmeticError("direct polynomial moment failed")
    shape_failure = None
    try:
        literal_moments(b @ b.T / 4, z, 8)
    except ValueError as error:
        shape_failure = str(error)
    if shape_failure is None:
        raise ArithmeticError("printed shape mismatch was not exposed")
    literal = literal_moments(h, z, 8)
    # The dimension-only diagnostic still uses the printed indices and update.
    if not np.allclose(literal[3], 0.5, atol=1e-14, rtol=0):
        raise ArithmeticError("literal degree/index diagnostic changed")
    setup = comparison_setup(b, rank)
    sampled, draws, weights, counts = sample_gram(b, np.full(4, 0.25), 1, 0)
    failure = assess(setup[0], sampled, setup[1], setup[2], rank)
    if not failure["rank_loss"] or failure["psd_sandwich_pass"]["0.5"]:
        raise ArithmeticError("one-face rank-loss control was lost")
    zeros, zero_draws, _, _ = sample_gram(np.zeros((3, 2)), np.zeros(2), 4, 0)
    if np.any(zeros) or zero_draws:
        raise ArithmeticError("zero-rank control failed")
    # Exact enumeration proves the k=1 expectation, independently of RNG.
    expected = sum(np.outer(b[:, j], b[:, j]) for j in range(4))
    if not np.array_equal(expected, b @ b.T):
        raise ArithmeticError("unbiased expectation control failed")
    return {
        "boundary": b.tolist(), "edge_basis": edges, "face_basis": faces,
        "exact_rank": rank, "exact_scores": [str(x) for x in scores],
        "probabilities": [0.25] * 4, "scaled_dual_gram": h.tolist(),
        "all_16_probes": z.tolist(), "section43_moments": moments.tolist(),
        "literal_shape_failure": shape_failure,
        "literal_with_dual_gram_only": literal.tolist(),
        "one_face_failure": {**failure, "draws": draws, "weights": weights.tolist(), "counts": counts},
        "zero_rank": "explicit zero; no normalization or draws",
        "histogram_to_probability_paper_status": "unresolved; analytic-bin branch is a project diagnostic",
    }


def fixtures():
    np = numpy()
    generator = np.random.Generator(np.random.PCG64(20261124))
    points = generator.normal(0, 0.3, (8, 2))
    points[4:] += 2
    distances = np.linalg.norm(points[:, None] - points[None, :], axis=2)
    result = []
    for scale in (0.6, 1.0, 2.8, 4.0):
        cells = [[(i,) for i in range(8)]]
        for size in (2, 3):
            cells.append([cell for cell in itertools.combinations(range(8), size)
                          if all(distances[i, j] <= scale for i, j in itertools.combinations(cell, 2))])
        for degree in (0, 1):
            b = boundary_matrix(cells[degree], cells[degree + 1])
            result.append({"name": f"rips-{scale}-q{degree}", "degree": degree,
                           "scale": scale, "points": points.tolist(), "cells": cells,
                           "boundary": b.tolist(), "boundary_shape": list(b.shape),
                           "classification": "project two-cluster Rips control; not author data"})
    return result


def topology_check(fixture):
    import gudhi
    from ripser import ripser
    np = numpy()
    points = np.array(fixture["points"])
    tree = gudhi.RipsComplex(points=points, max_edge_length=fixture["scale"]).create_simplex_tree(max_dimension=2)
    expected = {tuple(cell) for cells in fixture["cells"] for cell in cells}
    observed = {tuple(cell) for cell, _ in tree.get_filtration()}
    if observed != expected:
        raise ArithmeticError("GUDHI input simplex mismatch")
    ranks = [len(rref(boundary_matrix(fixture["cells"][d - 1], fixture["cells"][d]).tolist(),
                      len(fixture["cells"][d]))[1]) for d in (1, 2)]
    betti = [len(fixture["cells"][0]) - ranks[0],
             len(fixture["cells"][1]) - sum(ranks)]
    tree.persistence(homology_coeff_field=2, persistence_dim_max=True)
    gudhi_betti = tree.betti_numbers()[:2]
    ripser_result = ripser(points, maxdim=1, thresh=fixture["scale"], coeff=2)
    ripser_betti = [int(np.count_nonzero(np.isinf(diagram[:, 1]))) for diagram in ripser_result["dgms"]]
    if betti != gudhi_betti or betti != ripser_betti:
        raise ArithmeticError("exact/GUDHI/Ripser Betti mismatch")
    return {"real_exact_betti_0_1": betti, "gudhi_f2_betti_0_1": gudhi_betti,
            "ripser_f2_betti_0_1": ripser_betti,
            "scope": "input simplices and field-2 Betti0/1 corroboration; no real-spectrum proof"}


def worker_parity(worker, output, name, b, z, degree, p, count, seed):
    """Reuse frozen arrays/draws, comparing only the verified Rust subroutines."""
    np = numpy()
    sampled, draws, weights, _ = sample_gram(b, p, count, seed)
    if not draws:
        draws = [0] * count  # Explicitly ignored placeholders for zero-rank input.
    header = [*b.shape, degree, z.shape[1], count]
    tokens = [*header, *b.astype(int).ravel().tolist(), *z.ravel().tolist(), *p.tolist(), *draws]
    request = " ".join(map(str, tokens)) + "\n"
    (output / f"{name}.stdin").write_text(request, encoding="utf-8")
    result = subprocess.run([str(worker)], input=request, capture_output=True, text=True, timeout=60)
    (output / f"{name}.stdout").write_text(result.stdout, encoding="utf-8")
    (output / f"{name}.stderr").write_text(result.stderr, encoding="utf-8")
    if result.returncode:
        raise ArithmeticError(f"Rust worker failed on {name}: {result.stderr}")
    actual = json.loads(result.stdout)
    scores, rank = exact_scores(b)
    gram = b.T @ b
    scale = float(np.linalg.eigvalsh(gram)[-1]) if rank else 0.
    h = gram / scale if rank else gram
    expected = {"scores": np.array([float(x) for x in scores]),
                "moments": kid_moments(h, z, degree).ravel(),
                "weights": weights, "sampled": sampled.ravel()}
    errors = {}
    if actual["rank"] != rank or abs(actual["scale"] - scale) > TOLERANCE * max(1, scale):
        raise ArithmeticError(f"Rust rank/scale mismatch on {name}")
    for key, value in expected.items():
        candidate = np.array(actual[key])
        if candidate.shape != value.shape:
            raise ArithmeticError(f"Rust {key} shape mismatch on {name}")
        error = float(np.linalg.norm(candidate - value))
        errors[key] = error
        if error > TOLERANCE * max(1, float(np.linalg.norm(value))):
            raise ArithmeticError(f"Rust {key} mismatch on {name}: {error}")
    return {"name": name, "rank": rank, "errors": errors,
            "scope": "numerical exact-score route, fixed-draw Gram assembly, corrected moments; no histogram/RNG parity"}


def worker_rejections(worker, output):
    # n=p=probes=degree=count=1, B=Z=p=1 and draw=0 is the valid base.
    base = [1, 1, 1, 1, 1, 1, 1, 1, 0]
    cases = {"truncated": base[:-1], "trailing": base + [0]}
    for name, index, value in (("cap", 0, 129), ("probe", 6, 0),
                               ("nan", 7, "NaN"), ("negative", 7, -1),
                               ("sum", 7, 0.5), ("draw", 8, 1), ("count", 4, 0)):
        changed = list(base)
        changed[index] = value
        cases[name] = changed
    for name, tokens in cases.items():
        request = " ".join(map(str, tokens))
        result = subprocess.run([str(worker)], input=request, capture_output=True, text=True, timeout=60)
        write_json(output / f"reject-{name}.json", {"stdin": request, "returncode": result.returncode,
                                                  "stdout": result.stdout, "stderr": result.stderr})
        if not result.returncode or result.stdout:
            raise ArithmeticError(f"Rust worker accepted invalid {name}")
    return len(cases)


def run(output, sources, quick=False, worker=None):
    np = numpy()
    if any(not (sources / name).is_file() or digest(sources / name) != fingerprint
           for name, fingerprint in PAPERS.items()):
        raise ValueError("missing paper PDFs or SHA-256 mismatch")
    output.mkdir(parents=True, exist_ok=False)
    started = datetime.now(timezone.utc).isoformat()
    tracemalloc.start()
    start = time.perf_counter()
    environment = {name: importlib.metadata.version(name)
                   for name in ("numpy", "scipy", "gudhi", "ripser")}
    root = Path(__file__).resolve().parents[1]
    git = lambda *args: subprocess.check_output(["git", *args], cwd=root, text=True).strip()
    metadata = {
        "started_utc": started, "python": sys.version, "platform": platform.platform(),
        "processor": platform.processor(), "packages": environment, "command": sys.argv,
        "script_sha256": digest(__file__), "git_head": git("rev-parse", "HEAD"),
        "git_status": git("status", "--porcelain"), "papers": PAPERS,
        "author_code": "not verified in pinned paper links; independent mathematical implementation",
        "threads": {name: os.environ.get(name) for name in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS")},
        "contract": {"field": "real; ordinary unweighted signed chains", "dtype": "float64",
                     "rank_oracle": "exact Fraction row reduction; spectral rank diagnostic threshold 1e-9*max(1,||G||2)",
                     "matrix_tolerance": TOLERANCE, "moment_degrees": [8, 16, 32],
                     "probes": [32, 128, 512], "seeds": list(range(1 if quick else 20)),
                     "rng": "NumPy PCG64; nested probes, same draw seed across routes",
                     "integration": "64 positive bins of width 1/64; exact Chebyshev bin mass; no damping",
                     "kernel": "rectangular bin average; not the paper's unspecified convolution kernel",
                     "epsilon_targets": [0.25, 0.5], "delta_target": 0.05,
                     "bound_status": "empirical only; no sample-count or joint-probability theorem"},
        "evidence_class": "correctness and single-process phase diagnostics; no performance ranking",
        "rust_parity": "verified subroutines only" if worker else "not requested",
        "worker": {"path": str(worker), "sha256": digest(worker)} if worker else None,
        "source_sha256": {str(path.relative_to(root)): digest(path)
                          for path in (Path(__file__).resolve(), root / "tools/reference/sparsification.rs",
                                       root / "Cargo.toml", root / "Cargo.lock")},
    }
    write_json(output / "metadata.json", metadata)
    control = controls()
    write_json(output / "controls.json", control)
    parity = []
    if worker:
        b = np.array(control["boundary"])
        parity.append(worker_parity(worker, output, "rust-shell", b,
                                    np.array(control["all_16_probes"]), 8, np.full(4, 0.25), 1, 0))
        for rows, columns in ((0, 0), (3, 0), (3, 2)):
            parity.append(worker_parity(worker, output, f"rust-zero-{rows}-{columns}",
                                        np.zeros((rows, columns)), np.ones((columns, 1)),
                                        3, np.zeros(columns), 1, 0))
        rejections = worker_rejections(worker, output)
    inputs = fixtures()
    write_json(output / "inputs.json", inputs)
    input_hash = digest(output / "inputs.json")
    totals = Counter()
    phase_seconds = Counter()
    seeds = range(1 if quick else 20)
    with (output / "trials.jsonl").open("w", encoding="utf-8") as trials, \
            (output / "moments.jsonl").open("w", encoding="utf-8") as moment_file:
        for fixture in inputs:
            b = np.array(fixture["boundary"], dtype=float).reshape(fixture["boundary_shape"])
            stamp = time.perf_counter()
            rational_scores, rank = exact_scores(b)
            tau = np.array([float(x) for x in rational_scores])
            numerical = numerical_scores(b)
            if not np.allclose(tau, numerical, atol=TOLERANCE, rtol=0):
                raise ArithmeticError("pinv and independent exact row projector disagree")
            phase_seconds["exact_scores_and_independent_validation"] += time.perf_counter() - stamp
            stamp = time.perf_counter()
            gram, whitening, kernel = comparison_setup(b, rank)
            phase_seconds["comparison_projection_setup"] += time.perf_counter() - stamp
            topology = topology_check(fixture)
            name = fixture["name"]
            write_json(output / f"{name}-oracle.json", {
                "rank": rank, "rational_scores": [str(x) for x in rational_scores],
                "pinv_scores": numerical.tolist(), "score_sum": float(tau.sum()),
                "topology": topology, "gram": gram.tolist(), "kernel_basis": kernel.tolist(),
                "input_boundary_sha256": hashlib.sha256(b.astype("<f8").tobytes()).hexdigest(),
                "boundary_nnz": int(np.count_nonzero(b)), "dense_gram_entries": gram.size,
                "gram_eigenvalues": np.linalg.eigvalsh(gram).tolist(),
            })
            totals["fixtures"] += 1
            if not rank:
                continue
            exact_p = tau / rank
            counts = sorted(set((1, rank, 4 * rank, 16 * rank)))

            def record(route, p, seed, count, extra=None):
                stamp = time.perf_counter()
                sampled, draws, weights, draw_counts = sample_gram(b, p, count, seed)
                phase_seconds["sampling_and_gram_fill"] += time.perf_counter() - stamp
                stamp = time.perf_counter()
                diagnostics = assess(gram, sampled, whitening, kernel, rank)
                phase_seconds["full_generalized_spectrum_and_kernel_validation"] += time.perf_counter() - stamp
                row = {"fixture": name, "route": route, "seed": seed, "k": count,
                       "probabilities": p.tolist(), "draws": draws, "counts": draw_counts,
                       "gram_weights": weights.tolist(), "sqrt_weights": np.sqrt(weights).tolist(),
                       "distinct_columns": int(np.count_nonzero(weights)), **diagnostics, **(extra or {})}
                trials.write(json.dumps(row, allow_nan=False) + "\n")
                totals[f"{route}_trials"] += 1
                totals[f"{route}_rank_loss"] += diagnostics["rank_loss"]
                totals[f"{route}_sandwich_failure_eps0.5"] += not diagnostics["psd_sandwich_pass"]["0.5"]

            for seed in seeds:
                for count in counts:
                    record("exact", exact_p, seed, count)
                    record("uniform", np.full(b.shape[1], 1 / b.shape[1]), seed, count)
            h_raw = b.T @ b
            scale = float(np.linalg.eigvalsh(h_raw)[-1])
            h = h_raw / scale
            direct = direct_moments(h, 32)
            for seed in seeds:
                probes = np.random.Generator(np.random.PCG64(seed)).choice((-1.0, 1.0), (b.shape[1], 512))
                np.save(output / f"{name}-seed{seed}-probes.npy", probes)
                for probe_count in (32, 128, 512):
                    stamp = time.perf_counter()
                    all_moments = kid_moments(h, probes[:, :probe_count], 32)
                    phase_seconds["kid_moment_recurrence"] += time.perf_counter() - stamp
                    for degree in (8, 16, 32):
                        moments = all_moments[:degree + 1]
                        if worker and seed == 0:
                            parity.append(worker_parity(worker, output,
                                                        f"rust-{name}-m{degree}-z{probe_count}",
                                                        b, probes[:, :probe_count], degree,
                                                        exact_p, counts[-1], seed))
                        stamp = time.perf_counter()
                        scores, bin_masses = analytic_bin_scores(moments)
                        phase_seconds["project_analytic_bin_integration"] += time.perf_counter() - stamp
                        raw_sum = float(scores.sum())
                        p = scores / raw_sum if raw_sum > 0 else np.zeros_like(scores)
                        valid = bool(np.isfinite(p).all() and (scores >= 0).all()
                                     and raw_sum > 0 and not np.any((tau > 0) & (p <= 0)))
                        error = float(np.max(np.abs(p - exact_p))) if raw_sum > 0 else None
                        moment_file.write(json.dumps({
                            "fixture": name, "seed": seed, "degree": degree, "probes": probe_count,
                            "scale": scale, "moments": moments.tolist(),
                            "moment_linf_error_vs_eigh": float(np.max(np.abs(moments - direct[:degree + 1]))),
                            "raw_scores": scores.tolist(), "raw_sum": raw_sum,
                            "raw_sum_over_exact_rank": raw_sum / rank, "negative_scores": int(np.count_nonzero(scores < 0)),
                            "negative_bin_masses": int(np.count_nonzero(bin_masses < 0)),
                            "probabilities": p.tolist(), "probability_linf_error": error,
                            "normalized_sum": float(p.sum()), "valid_for_sampling": valid,
                            "status": "project integral diagnostic; not full paper KID",
                        }, allow_nan=False) + "\n")
                        totals["kid_moment_cells"] += 1
                        totals["kid_negative_score_cells"] += bool((scores < 0).any())
                        if not valid:
                            totals["kid_rejected_probability_cells"] += 1
                            continue
                        for count in counts:
                            record("kid_project_integral", p, seed, count,
                                   {"degree": degree, "probes": probe_count, "probability_linf_error": error})
    current, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    if worker:
        write_json(output / "rust-parity.json", {"cases": parity, "rejections": rejections})
    write_json(output / "summary.json", {
        "counts": dict(totals), "input_file_sha256": input_hash,
        "phase_seconds": dict(phase_seconds), "total_seconds": time.perf_counter() - start,
        "python_tracemalloc_peak_bytes": peak, "python_tracemalloc_current_bytes": current,
        "memory_scope": "Python allocations only; not process RSS or total NumPy/BLAS memory",
        "ended_utc": datetime.now(timezone.utc).isoformat(),
        "exact_sampling": "reproduced on declared controls; empirical sandwich failures retained",
        "section43_moments": "analytic controls reproduced; stochastic estimates retained",
        "full_paper_kid": "unresolved: histogramIntegral, printed shape/index/update/factor differences",
        "rust_parity_cases": len(parity),
        "admission": "defer core integration, full KID rewrite and optimization; no validated full KID or cost benefit",
    })
    return dict(totals)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sources", required=True, type=Path, help="directory with both pinned PDFs")
    parser.add_argument("--output", required=True, type=Path, help="new artifact directory; never overwritten")
    parser.add_argument("--quick", action="store_true", help="seed 0 only; full default uses seeds 0..19")
    parser.add_argument("--worker", type=Path, help="optional compiled standalone Rust parity worker")
    args = parser.parse_args()
    if args.output.exists():
        print("Output already exists; choose a new directory.", file=sys.stderr)
        return 1
    try:
        worker = args.worker.resolve() if args.worker else None
        result = run(args.output, args.sources, args.quick, worker)
    except (ValueError, ArithmeticError, OSError, ImportError, subprocess.SubprocessError) as error:
        if args.output.is_dir():
            write_json(args.output / "failure.json", {"error": str(error), "type": type(error).__name__})
        print(f"Sparsification reproduction failed: {error}", file=sys.stderr)
        return 1
    print(json.dumps(result, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
