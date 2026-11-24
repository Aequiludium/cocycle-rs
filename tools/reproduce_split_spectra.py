"""Bounded paper reproduction for #145; no production spectral API.

PETLS 2508.11560v2 Sections 4.4/4.5 and Duersch et al. 1704.07458v3
Algorithms 1/2 are independently transcribed from mathematics. NumPy, SciPy,
SymPy and optionally GUDHI are external verification dependencies only.
"""

import argparse
import hashlib
import itertools
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import time
import warnings
import zipfile


ROOT = Path(__file__).resolve().parents[1]
PETLS_REVISION = "5b6564e9e05af9d24928313b92a014696247bdd7"
PAPERS = {
    "petls-v2.pdf": "5e14d20b0230d98550be0b802d0dfe9635d232bf794a33672b85fa151706fb86",
    "lobpcg-v3.pdf": "01d447acaea3d0ec49a427e554e499791a5e23c0f38fcdb8c50e159ac21bbc07",
}
RCOND = 1e-12
MATRIX_TOLERANCE = 1e-10
RESIDUAL_TOLERANCE = 1e-9


def libraries():
    # Import lazily so standard-library tool discovery needs no scientific stack.
    import numpy as np
    import scipy.linalg as la
    return np, la


def matrix(value):
    np, _ = libraries()
    if np.iscomplexobj(value):
        raise ValueError("require a real matrix")
    result = np.asarray(value, dtype=float)
    if result.ndim != 2 or not np.isfinite(result).all():
        raise ValueError("require a finite real matrix")
    return result


def symmetric(value):
    np, _ = libraries()
    h = matrix(value)
    if h.shape[0] != h.shape[1] or not np.allclose(h, h.T, rtol=0, atol=0):
        raise ValueError("require a symmetric square matrix")
    return h


def relative_error(actual, expected):
    np, _ = libraries()
    return float(np.linalg.norm(actual - expected) / max(1., np.linalg.norm(expected)))


def check_close(actual, expected, tolerance=MATRIX_TOLERANCE):
    np, _ = libraries()
    if np.shape(actual) != np.shape(expected):
        raise AssertionError(f"shape mismatch: {np.shape(actual)} != {np.shape(expected)}")
    error = relative_error(np.asarray(actual), np.asarray(expected))
    if not np.isfinite(error) or error > tolerance:
        raise AssertionError(f"relative Frobenius error {error} > {tolerance}")
    return error


def boundary(rows, columns):
    np, _ = libraries()
    result = np.zeros((len(rows), len(columns)), dtype=int)
    positions = {s: i for i, s in enumerate(rows)}
    for j, simplex in enumerate(columns):
        if len(simplex) > 1:
            for i in range(len(simplex)):
                result[positions[simplex[:i] + simplex[i + 1:]], j] = (-1) ** i
    return result


def exact_matrix(a):
    import sympy as sp
    a = matrix(a)
    if not (a == a.astype(int)).all():
        raise ValueError("exact incidence oracle requires integral input")
    return sp.Matrix(a.shape[0], a.shape[1], list(map(int, a.flat)))


def exact_kernel(a):
    import sympy as sp
    np, _ = libraries()
    a = exact_matrix(a)
    columns = a.nullspace()
    z = sp.Matrix.hstack(*columns) if columns else sp.zeros(a.cols, 0)
    return z, np.array(z, dtype=float).reshape(z.rows, z.cols)


def gram_flip(a):
    """Section 4.4: retain both full spectra, then restore original-domain zeros.

    The threshold is numerical, not an exact rank certificate. Tiny positive
    values are retained in raw spectra and flagged when their classification
    lies within a factor of ten of the cutoff.
    """
    np, la = libraries()
    a = matrix(a)
    started = time.perf_counter()
    primal, dual = a.T @ a, a @ a.T
    construction = time.perf_counter() - started
    started = time.perf_counter()
    lp, _ = la.eigh(primal)
    ld, zd = la.eigh(dual)
    solve = time.perf_counter() - started
    threshold = RCOND * max(float(np.max(np.abs(lp), initial=0)),
                            float(np.max(np.abs(ld), initial=0)))
    if np.min(lp, initial=0) < -threshold or np.min(ld, initial=0) < -threshold:
        raise ValueError("Gram spectrum is significantly negative")
    positive = ld > threshold
    lifted = a.T @ zd[:, positive] / np.sqrt(ld[positive])
    residual = relative_error(primal @ lifted, lifted * ld[positive])
    rank = int(positive.sum())
    reconstructed = np.concatenate((np.zeros(a.shape[1] - rank), ld[positive]))
    return dict(primal=primal, dual=dual, primal_spectrum=lp, dual_spectrum=ld,
                lifted_vectors=lifted, lifted_residual=residual,
                positive=ld[positive], reconstructed=reconstructed, threshold=threshold,
                numerical_rank=rank, original_domain=a.shape[1],
                selected_side="dual" if a.shape[0] < a.shape[1] else "primal",
                ambiguous=bool(np.any((lp > 0) & (lp <= 10 * threshold))),
                seconds=dict(construction=construction, eigsolve=solve))


def harmonic_restriction(h, kernel):
    """Section 4.5 with an actual complete kernel and orthonormal complement.

    Unlike arbitrary homology representatives, supplied columns must be
    independent zero modes. Completeness and positivity are numerical checks;
    callers must retain any rank ambiguity and the unfiltered reduced spectrum.
    """
    np, la = libraries()
    started = time.perf_counter()
    h, kernel = symmetric(h), matrix(kernel)
    n = len(h)
    if kernel.shape[0] != n or kernel.shape[1] > n:
        raise ValueError("harmonic basis shape")
    scale = float(la.norm(h, 2)) if n else 0.
    if kernel.shape[1]:
        singular = la.svdvals(kernel)
        if singular[-1] <= RCOND * singular[0]:
            raise ValueError("harmonic basis is dependent")
        q, _ = la.qr(kernel, mode="economic")
    else:
        q = kernel.copy()
    leakage = float(la.norm(h @ q) / max(scale, np.finfo(float).tiny))
    if leakage > RESIDUAL_TOLERANCE:
        raise ValueError("supplied representatives are not harmonic")
    x = la.null_space(q.T, rcond=RCOND)
    setup = time.perf_counter() - started
    started = time.perf_counter()
    reduced = x.T @ h @ x
    reduced = (reduced + reduced.T) / 2
    projection = time.perf_counter() - started
    started = time.perf_counter()
    values, vectors = la.eigh(reduced)
    solve = time.perf_counter() - started
    threshold = RCOND * scale
    if len(values) and values[0] <= 0:
        raise ValueError("kernel is incomplete or restricted operator is not positive")
    original = x @ vectors
    return dict(complement=x, kernel=q, reduced=reduced, spectrum=values,
                vectors=original, threshold=threshold,
                ambiguous=bool(np.any(values <= 10 * threshold)),
                leakage=leakage, residual=relative_error(h @ original, original * values),
                orthogonality=relative_error(original.T @ original, np.eye(len(values))),
                seconds=dict(setup=setup, projection=projection, eigsolve=solve))


def rayleigh_ritz(h, search, count, compress=False):
    """Algorithm 2, including diagonal scaling and checked Cholesky.

    compress=True is a separately labeled project correction: SVD of normalized
    search columns removes dependencies before applying Algorithm 2. It is not
    the paper's robust Algorithm 8. Coefficients are mapped back to raw S.
    """
    np, la = libraries()
    h, search = symmetric(h), matrix(search)
    if search.shape[0] != len(h) or not 0 < count <= min(search.shape):
        raise ValueError("Rayleigh-Ritz shape/count")
    norms = la.norm(search, axis=0)
    mapping = np.eye(search.shape[1])
    if compress:
        keep = norms > np.finfo(float).tiny
        normalized = search[:, keep] / norms[keep]
        _, singular, vt = la.svd(normalized, full_matrices=False)
        rank = int(np.sum(singular > RCOND * singular[0])) if len(singular) else 0
        if rank < count:
            raise ValueError("insufficient search rank")
        mapping = np.zeros((search.shape[1], rank))
        mapping[keep] = vt[:rank].T / (norms[keep, None] * singular[:rank])
        s = search @ mapping
    else:
        s = search
    metric = s.T @ s
    diagonal = np.diag(metric)
    if np.any(diagonal <= 0):
        raise ValueError("zero Rayleigh-Ritz search column")
    d = 1 / np.sqrt(diagonal)
    scaled = d[:, None] * metric * d[None, :]
    condition = float(np.linalg.cond(scaled))
    if not np.isfinite(condition) or condition > 1e12:
        raise ValueError("dependent or ill-conditioned Rayleigh-Ritz search")
    r = la.cholesky(scaled, lower=False)
    transform = d[:, None] * la.solve_triangular(r, np.eye(len(r)))
    projected = transform.T @ (s.T @ h @ s) @ transform
    values, z = la.eigh((projected + projected.T) / 2)
    coefficients = mapping @ transform @ z[:, :count]
    return coefficients, values[:count], dict(metric=metric, scaled_metric=scaled,
                                              cholesky=r, projected=projected,
                                              coefficients=coefficients, condition=condition)


def basic_lobpcg(h, initial, count, max_iterations=120, tolerance=1e-9,
                 compress=False, trace=None):
    """Algorithm 1 for B=I and identity preconditioner, with fresh residuals.

    No success on rank loss or exhaustion; no smallest-index certificate. The
    initial residual can already converge (an explicit zero-update guard).
    """
    np, la = libraries()
    h, x = symmetric(h), matrix(initial).copy()
    block = x.shape[1]
    if x.shape[0] != len(h) or not 0 < count <= block <= len(h):
        raise ValueError("LOBPCG requested count/block/domain")
    if max_iterations < 0 or not np.isfinite(tolerance) or tolerance <= 0:
        raise ValueError("LOBPCG stopping policy")
    trace = [] if trace is None else trace
    c, values, rr = rayleigh_ritz(h, x, block, compress)
    x = x @ c
    p = np.empty((len(h), 0))
    for iteration in range(max_iterations + 1):
        residual = h @ x - x * values
        errors = la.norm(residual, axis=0)
        # Fixed normalization for this entire run; no post-hoc tolerance changes.
        scale = max(1., float(la.norm(h, 2)))
        trace.append(dict(iteration=iteration, X=x.copy(), P=p.copy(), R=residual,
                          eigenvalues=values, normalized_residuals=errors / scale, ritz=rr))
        if np.all(errors[:count] <= tolerance * scale):
            vectors = x[:, :count]
            return dict(values=values[:count], vectors=vectors,
                        normalized_residuals=errors[:count] / scale,
                        orthogonality=float(la.norm(vectors.T @ vectors - np.eye(count))),
                        iterations=iteration, index_coverage_certified=False)
        if iteration == max_iterations:
            raise ValueError("LOBPCG iteration limit")
        w = residual.copy()  # K=I; do not replace by SciPy's implementation.
        search = np.column_stack((x, w, p))
        trace[-1].update(W=w, S=search)
        c, values, rr = rayleigh_ritz(h, search, block, compress)
        x = search @ c
        p = search[:, block:] @ c[block:]


def jsonable(value):
    if isinstance(value, dict):
        return {str(k): jsonable(v) for k, v in value.items()}
    if isinstance(value, (tuple, list)):
        return [jsonable(v) for v in value]
    if hasattr(value, "tolist"):
        if value.ndim == 2:
            # [] alone loses the domain of a 0-by-p boundary or n-by-0 kernel.
            return dict(shape=list(value.shape), data=value.tolist())
        return value.tolist()
    return value


def write_json(path, value):
    path.write_text(json.dumps(jsonable(value), indent=2, allow_nan=False) + "\n", encoding="utf-8")


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def reference_manifest(references):
    """Verify PDFs and extracted author files without redistributing source."""
    manifest = dict(papers={}, author_runtime="not established")
    for name, expected in PAPERS.items():
        actual = sha256(references / name)
        if actual != expected:
            raise ValueError(f"unexpected paper hash: {name}")
        manifest["papers"][name] = actual
    archive = references / "petls-source.zip"
    if archive.is_file():
        prefix = f"PETLS-{PETLS_REVISION}"
        source = references / "petls-source" / prefix
        paths = ["CMakeLists.txt", "cpp/CMakeLists.txt",
                 "cpp/include/petls_headers/core/Complex.hpp",
                 "cpp/src/core/Complex.cpp",
                 "cpp/include/petls_headers/core/up_algorithms.hpp",
                 "cpp/src/eigenvalues/eigs_algorithms.cpp"]
        manifest["author_archive_sha256"] = sha256(archive)
        manifest["author_files"] = {}
        with zipfile.ZipFile(archive) as original:
            for path in paths:
                expected = hashlib.sha256(original.read(f"{prefix}/{path}")).hexdigest()
                actual = sha256(source / path)
                if actual != expected:
                    raise ValueError(f"extracted author file differs from archive: {path}")
                manifest["author_files"][path] = actual
        manifest["license"] = "no top-level license in inspected author snapshot; no source copied"
    log = references / "petls-configure.log"
    if log.is_file():
        manifest["configure_log_sha256"] = sha256(log)
    return manifest


def fixtures():
    """Project controls, not the authors' datasets or published timings."""
    np, _ = libraries()
    vertices = [(i,) for i in range(4)]
    edges = list(itertools.combinations(range(4), 2))
    faces = list(itertools.combinations(range(4), 3))
    triangle_edges = list(itertools.combinations(range(3), 2))
    return dict(tall_path=boundary(vertices, [(0, 1), (1, 2)]),
                wide_triangle=boundary([(i,) for i in range(3)], triangle_edges),
                triangle_top=boundary(triangle_edges, [(0, 1, 2)]),
                shell=boundary(edges, faces), zero=np.zeros((2, 5), dtype=int),
                empty_rows=np.zeros((0, 3), dtype=int),
                empty_columns=np.zeros((3, 0), dtype=int),
                empty=np.zeros((0, 0), dtype=int))


def pair_controls():
    """All B columns retained. Exact rational induced metric is the oracle."""
    np, la = libraries()
    import sympy as sp
    edges = list(itertools.combinations(range(3), 2))
    a = boundary([(0,), (1,), (2,)], edges)
    b = boundary(edges, [(0, 1, 2)])
    # Shared persistent control: source has isolated 0,1 and filled 234;
    # the later path 0-5-1 and filled 6,7,8 have no q=0 source columns.
    vertices = [(i,) for i in range(9)]
    target_edges = [(2, 3), (2, 4), (3, 4), (0, 5), (1, 5), (6, 7), (6, 8), (7, 8)]
    full_b = boundary(vertices, target_edges)
    for name, lower, b1, b2, expected in [
        ("filled_triangle", a, b, np.zeros((0, 1), dtype=int), [3., 3., 3.]),
        ("persistent_path_and_disconnected_triangle", np.zeros((0, 5), dtype=int),
         full_b[:5], full_b[5:], [0., 0., 1., 3., 3.]),
    ]:
        started = time.perf_counter()
        ez, z = exact_kernel(b2)
        exact_kernel_time = time.perf_counter() - started
        started = time.perf_counter()
        q, _ = la.qr(z, mode="economic")
        f = b1 @ q
        u, d = f @ f.T, lower.T @ lower
        h = u + d
        factor_time = time.perf_counter() - started
        ea, eb1, eb2 = map(exact_matrix, (lower, b1, b2))
        ef = eb1 * ez
        eu = ef * (ez.T * ez).inv() * ef.T if ez.cols else sp.zeros(b1.shape[0])
        ed = ea.T * ea
        if ea * ef != sp.zeros(ea.rows, ef.cols):
            raise AssertionError("exact chain identity A F != 0")
        if eu * ed != sp.zeros(eu.rows) or ed * eu != sp.zeros(eu.rows):
            raise AssertionError("exact component product is nonzero")
        check_close(u, np.array(eu, dtype=float))
        check_close(lower @ f, np.zeros((lower.shape[0], f.shape[1])))
        check_close(u @ d, np.zeros(h.shape))
        check_close(d @ u, np.zeros(h.shape))
        upper_flip, lower_flip = gram_flip(f.T), gram_flip(lower)
        started = time.perf_counter()
        positives = np.sort(np.concatenate((upper_flip["positive"], lower_flip["positive"])))
        rank_a, rank_f = ea.rank(), ef.rank()
        union = np.concatenate((np.zeros(len(h) - rank_a - rank_f), positives))
        recombination_time = time.perf_counter() - started
        check_close(union, np.asarray(expected))
        check_close(la.eigvalsh(h), np.asarray(expected))
        exact_h = eu + ed
        kernel_columns = exact_h.nullspace()
        kernel = sp.Matrix.hstack(*kernel_columns) if kernel_columns else sp.zeros(len(h), 0)
        restriction = harmonic_restriction(h, np.array(kernel, dtype=float).reshape(len(h), -1))
        check_close(restriction["spectrum"], positives)
        yield name, dict(A=lower, B1=b1, B2=b2, Z=z, Q=q, F=f, U=u, D=d, Delta=h,
                         exact_rank_A=rank_a, exact_rank_F=rank_f, expected=expected,
                         full_spectrum=la.eigvalsh(h), union=union, restriction=restriction,
                         upper_flip=upper_flip, lower_flip=lower_flip,
                         seconds=dict(exact_kernel=exact_kernel_time, factor=factor_time,
                                      recombination=recombination_time))


def solver_controls():
    np, _ = libraries()
    rng = np.random.default_rng(145)
    for name, diagonal, block, count in [
        ("well_separated", np.arange(1., 49.), 3, 2),
        ("clustered", np.r_[1., 1.00001, 1.00002, np.arange(2., 47.)], 3, 3),
        ("multiplicity_exceeds_block", np.r_[np.ones(5), np.arange(2., 45.)], 2, 2),
    ]:
        rotation, _ = np.linalg.qr(rng.normal(size=(len(diagonal), len(diagonal))))
        h = (rotation * diagonal) @ rotation.T
        h = (h + h.T) / 2
        yield name, h, rng.normal(size=(len(h), block)), count


def run(output, references=None, gudhi=False):
    started = time.perf_counter()
    np, la = libraries()
    import scipy
    import sympy
    from scipy.sparse.linalg import lobpcg
    output.mkdir(parents=True, exist_ok=False)
    (output / "source-reproduce_split_spectra.py").write_bytes(Path(__file__).read_bytes())
    test_source = Path(__file__).with_name("test_reproduce_split_spectra.py")
    if test_source.is_file():
        (output / "source-test_reproduce_split_spectra.py").write_bytes(test_source.read_bytes())
    environment = dict(python=sys.version, numpy=np.__version__, scipy=scipy.__version__,
                       sympy=sympy.__version__, platform=platform.platform(),
                       processor=platform.processor(),
                       thread_environment={name: os.environ.get(name) for name in
                                           ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]},
                       command=sys.argv, seed=145, classification="project controls",
                       field="R (independent incidence ranks/kernels over Q)",
                       inner_product="identity; ordinary nonaugmented signed chains",
                       normalization="Frobenius/max(1,reference Frobenius); residual/ max(1,||H||2)",
                       rcond=RCOND, matrix_tolerance=MATRIX_TOLERANCE,
                       residual_tolerance=RESIDUAL_TOLERANCE,
                       petls_revision=PETLS_REVISION, paper_sha256=PAPERS,
                       paper_identity_verified=references is not None)
    for label, args in [("head", ["rev-parse", "HEAD"]), ("dirty", ["status", "--porcelain"])]:
        result = subprocess.run(["git", *args], cwd=ROOT, text=True, capture_output=True, check=False)
        environment[label] = result.stdout.strip()
    environment["source_sha256"] = sha256(Path(__file__))
    if references:
        environment["references"] = reference_manifest(references)
        log = references / "petls-configure.log"
        if log.is_file():
            (output / "petls-configure.log").write_bytes(log.read_bytes())
    results, failures, algorithm_outcomes = [], [], []

    def record(name, operation):
        t0 = time.perf_counter()
        try:
            data = operation()
            write_json(output / f"{name}.json", data)
            if isinstance(data, dict) and "variants" in data:
                for variant, result in data["variants"].items():
                    algorithm_outcomes.append(dict(case=name, variant=variant,
                                                   status=result["status"], reason=result.get("reason"),
                                                   iterations=result.get("answer", {}).get("iterations")))
            results.append(dict(case=name, status="pass", seconds=time.perf_counter() - t0))
        except Exception as error:
            failures.append(dict(case=name, exception=type(error).__name__, error=str(error)))
            write_json(output / f"{name}-failure.json", failures[-1])

    for name, a in fixtures().items():
        write_json(output / f"input-gram-{name}.json", dict(A=a))
        def flip(a=a):
            data = gram_flip(a)
            rank = exact_matrix(a).rank()
            check_close(data["reconstructed"], data["primal_spectrum"])
            check_close(data["positive"], data["primal_spectrum"][data["primal_spectrum"] > data["threshold"]])
            if data["numerical_rank"] != rank or data["lifted_residual"] > RESIDUAL_TOLERANCE:
                raise AssertionError("Gram rank or lifted residual")
            # Compare repeated eigenspaces as whole projectors, not signed vectors.
            _, vp = la.eigh(data["primal"])
            pv = vp[:, data["primal_spectrum"] > data["threshold"]]
            data["projector_error"] = check_close(data["lifted_vectors"] @ data["lifted_vectors"].T, pv @ pv.T)
            return dict(A=a, exact_rank=rank, **data)
        record(f"gram-{name}", flip)

    def pairs():
        return dict(pair_controls())
    record("pair-controls", pairs)

    def invalid_representative():
        h = np.array([[1., -1.], [-1., 1.]])
        z = np.array([[1.], [0.]])
        try:
            harmonic_restriction(h, z)
        except ValueError as error:
            reason = str(error)
        else:
            raise AssertionError("nonharmonic homology representative accepted")
        valid = harmonic_restriction(h, np.ones((2, 1)))
        check_close(valid["spectrum"], np.array([2.]))
        return dict(Delta=h, representative=z, action=h @ z, rejection=reason, valid=valid)
    record("representative-counterexample", invalid_representative)

    for name, h, initial, count in solver_controls():
        write_json(output / f"input-solver-{name}.json", dict(H=h, initial=initial, count=count))
        def solve(h=h, initial=initial, count=count, name=name):
            dense = la.eigvalsh(h)
            variants = {}
            for compress in (False, True):
                key = "project_rank_correction" if compress else "paper_basic"
                trace = []
                solve_started = time.perf_counter()
                try:
                    answer = basic_lobpcg(h, initial, count, compress=compress, trace=trace)
                    elapsed = time.perf_counter() - solve_started
                    check_close(answer["values"], dense[:count], 1e-8)
                    variants[key] = dict(status="converged", answer=answer, solve_seconds=elapsed)
                except (ValueError, la.LinAlgError) as error:
                    variants[key] = dict(status="failed", reason=str(error),
                                         solve_seconds=time.perf_counter() - solve_started)
                write_json(output / f"solver-{name}-{key}-trace.json", trace)
            solve_started = time.perf_counter()
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                values, vectors = lobpcg(h, initial.copy(), largest=False, tol=1e-8, maxiter=300)
            scipy_solve_seconds = time.perf_counter() - solve_started
            check_close(values[:count], dense[:count], 1e-8)
            if variants["project_rank_correction"]["status"] != "converged":
                raise AssertionError("project correction did not converge")
            return dict(H=h, initial=initial, count=count, dense=dense, variants=variants,
                        scipy_values=values, scipy_warnings=[str(w.message) for w in caught],
                        scipy_residual=relative_error(h @ vectors, vectors * values),
                        scipy_solve_seconds=scipy_solve_seconds,
                        scipy_dense_fallback=len(h) < 5 * initial.shape[1])
        record(f"solver-{name}", solve)

    def boundaries():
        h = np.diag(np.r_[np.zeros(32), 1., 2., 3.])
        restriction = harmonic_restriction(h, np.eye(35)[:, :32])
        exact_start = np.eye(3)[:, :2]
        answer = basic_lobpcg(restriction["reduced"], exact_start, 2)
        check_close(answer["values"], np.array([1., 2.]))
        trace = []
        iterative = basic_lobpcg(restriction["reduced"], np.ones((3, 1)), 1, trace=trace)
        lifted = restriction["complement"] @ iterative["vectors"]
        full_residual = relative_error(h @ lifted, lifted * iterative["values"])
        if full_residual > RESIDUAL_TOLERANCE:
            raise AssertionError("restricted solve failed in original coordinates")
        count_failure = None
        try:
            basic_lobpcg(restriction["reduced"], np.eye(3), 4)
        except ValueError as error:
            count_failure = str(error)
        if not count_failure:
            raise AssertionError("excessive positive count accepted")
        tiny = harmonic_restriction(np.diag([0., 1e-14, 1.]), np.eye(3)[:, :1])
        if not tiny["ambiguous"] or tiny["spectrum"][0] <= 0:
            raise AssertionError("tiny positive value was lost or treated as unambiguous")
        incomplete_rejected = False
        try:
            harmonic_restriction(np.diag([0., 1.]), np.empty((2, 0)))
        except ValueError:
            incomplete_rejected = True
        if not incomplete_rejected:
            raise AssertionError("incomplete kernel accepted")
        return dict(large_zero_block=restriction, initial_convergence=answer,
                    restricted_iteration=iterative, trace=trace, lifted_vectors=lifted,
                    original_coordinate_residual=full_residual,
                    excessive_count_rejection=count_failure, tiny_positive=tiny,
                    incomplete_kernel_rejected=incomplete_rejected)
    record("restriction-boundaries", boundaries)

    if gudhi:
        def topology():
            import gudhi as gd
            data = {}
            for name, facets, degree, expected in [
                ("triangle", [(0, 1, 2)], 1, 0),
                ("shell", list(itertools.combinations(range(4), 3)), 2, 1),
                ("edge", [(0, 1)], 0, 1),
            ]:
                tree = gd.SimplexTree()
                for facet in facets:
                    tree.insert(facet, filtration=0.)
                tree.compute_persistence(homology_coeff_field=2, persistence_dim_max=True)
                betti = len(tree.persistence_intervals_in_dimension(degree))
                if betti != expected:
                    raise AssertionError("GUDHI torsion-free Betti control")
                data[name] = dict(degree=degree, betti=betti)
            return dict(version=gd.__version__, field=2, controls=data,
                        exclusion="GUDHI PH does not validate real spectral values or general real ranks")
        record("gudhi-topology", topology)

    summary = dict(passed=len(results), failures=failures, cases=results,
                   algorithm_outcomes=algorithm_outcomes,
                   author_runtime="not established; independent transcriptions only",
                   seconds=time.perf_counter() - started,
                   evidence="bounded exact checks and numerical evidence, not a coverage certificate")
    write_json(output / "summary.json", summary)
    environment["input_sha256"] = {p.name: sha256(p) for p in sorted(output.glob("input-*.json"))}
    write_json(output / "environment.json", environment)
    write_json(output / "artifact-hashes.json", {p.name: sha256(p) for p in sorted(output.iterdir()) if p.is_file()})
    print(json.dumps(dict(output=str(output), passed=len(results), failures=failures)))
    return 1 if failures else 0


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True, help="fresh ignored output directory")
    parser.add_argument("--references", type=Path, help="directory with the two pinned PDFs")
    parser.add_argument("--gudhi", action="store_true", help="torsion-free topology corroboration only")
    args = parser.parse_args()
    return run(args.output, args.references, args.gudhi)


if __name__ == "__main__":
    raise SystemExit(main())
