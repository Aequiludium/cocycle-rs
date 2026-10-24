"""Independent ICML 2023 persistent spectral feature reproduction (#111).

This optional Python oracle is outside the Rust runtime and package. It uses
real, unweighted, nonaugmented chains. The bounded file-exchange route consumes
the existing Rust persistent_laplacian results; learning stays outside the crate.
"""

import argparse
import ast
from fractions import Fraction
import hashlib
import json
import math
import platform
from pathlib import Path
import random
import shutil
import subprocess
import sys
import time
import tracemalloc


ROOT = Path(__file__).resolve().parents[1]
AUTHOR_REVISION = "b36b080f739653d84284f3f7ca9ccea5fad286ae"
AUTHOR_BLOB = "ef49f6caeac3c7710de84795f3ffb225d14dbc29"
RCOND = 1e-12
ZERO_TOLERANCE = 1e-10
COMPARISON_TOLERANCE = 1e-9


def numpy():
    """Keep third-party packages optional for standard-library tool discovery."""
    import numpy as np
    return np


def sampling_grid(start, end, num_points):
    """Inclusive, uniform grid; num_points counts points, not intervals."""
    np = numpy()
    if (not math.isfinite(start) or not math.isfinite(end) or start > end
            or type(num_points) is not int or num_points < 1):
        raise ValueError("invalid sampling range or point count")
    if start == end and num_points == 1:
        return (float(start),)
    if start == end or num_points < 2:
        raise ValueError("nondegenerate ranges need at least two distinct points")
    return checked_grid(np.linspace(start, end, num_points).tolist())


def checked_grid(grid):
    grid = tuple(float(value) for value in grid)
    if (not grid or any(not math.isfinite(value) for value in grid)
            or any(a >= b for a, b in zip(grid, grid[1:]))):
        raise ValueError("grid must be nonempty, finite and strictly increasing")
    return grid


def faces(cell, kind):
    if len(cell) == 1:
        return []
    if kind == "simplicial" or len(cell) == 2:
        return [(cell[:index] + cell[index + 1:], (-1)**index)
                for index in range(len(cell))]
    # A regular-grid square has row-major corners a,b,c,d. Its oriented
    # boundary is ab + bd - cd - ac, ensuring B1 B2 = 0.
    a, b, c, d = cell
    return [((a, b), 1), ((b, d), 1), ((c, d), -1), ((a, c), -1)]


def validate_complex(complex_, kind="simplicial"):
    if kind not in ("simplicial", "cubical"):
        raise ValueError("unsupported complex type")
    if kind == "cubical" and len(complex_) > 3:
        raise ValueError("oracle supports only two-dimensional regular grids")
    for degree, cells in enumerate(complex_):
        expected_size = degree + 1 if kind == "simplicial" else 2**degree
        if (len(set(cells)) != len(cells)
                or any(len(cell) != expected_size or tuple(sorted(set(cell))) != cell
                       for cell in cells)):
            raise ValueError("invalid, duplicate or noncanonical cells")
        if degree and any(face not in set(complex_[degree - 1])
                          for cell in cells for face, _ in faces(cell, kind)):
            raise ValueError("complex is not closed under boundaries")


def boundary(complex_, degree, kind="simplicial"):
    np = numpy()
    cells = complex_[degree] if degree < len(complex_) else []
    previous = complex_[degree - 1] if degree else []
    matrix = np.zeros((len(previous), len(cells)))
    indices = {cell: index for index, cell in enumerate(previous)}
    if degree:
        for column, cell in enumerate(cells):
            for face, sign in faces(cell, kind):
                matrix[indices[face], column] = sign
    return matrix


def stage(filtration, scale, kind="simplicial"):
    """Closed sublevel sets, including cells born at the maximum threshold."""
    if not math.isfinite(scale):
        raise ValueError("nonfinite stage scale")
    levels = [[], [], []]
    for cell, birth in filtration.items():
        if not math.isfinite(birth):
            raise ValueError("nonfinite filtration value")
        degree = len(cell) - 1 if kind == "simplicial" else {1: 0, 2: 1, 4: 2}[len(cell)]
        if degree > 2:
            raise ValueError("small oracle supports dimensions through two")
        for face, _ in faces(cell, kind):
            if face not in filtration or filtration[face] > birth:
                raise ValueError("filtration violates face-before-coface order")
        if birth <= scale:
            levels[degree].append(cell)
    for cells in levels:
        cells.sort()
    validate_complex(levels, kind)
    return levels


def cubical_filtration(image):
    """2D vertex lower-star filtration (pixels are vertices, as in author code)."""
    np = numpy()
    image = np.asarray(image, dtype=float)
    if image.ndim != 2 or not image.size or not np.isfinite(image).all():
        raise ValueError("image must be a nonempty finite two-dimensional array")
    height, width = image.shape
    result = {(i * width + j,): float(image[i, j])
              for i in range(height) for j in range(width)}
    for i in range(height):
        for j in range(width):
            a = i * width + j
            for b in ([a + 1] if j + 1 < width else []) + ([a + width] if i + 1 < height else []):
                result[(a, b)] = max(result[(a,)], result[(b,)])
            if i + 1 < height and j + 1 < width:
                square = (a, a + 1, a + width, a + width + 1)
                result[square] = max(result[(v,)] for v in square)
    return result


def pair_blocks(K, L, degree, kind):
    validate_complex(K, kind)
    validate_complex(L, kind)
    if type(degree) is not int or degree not in (0, 1):
        raise ValueError("small oracle computes degrees zero and one")
    if len(K) < 3 or len(L) < 3:
        raise ValueError("explicit empty adjacent degrees are required")
    if any(not set(K[q]).issubset(L[q]) for q in range(3)):
        raise ValueError("pair is not an inclusion")
    positions = {cell: index for index, cell in enumerate(L[degree])}
    selected = [positions[cell] for cell in K[degree]]
    selected_set = set(selected)
    outside = [index for index in range(len(L[degree])) if index not in selected_set]
    return boundary(K, degree, kind), boundary(L, degree + 1, kind), selected, outside


def persistent_laplacian(K, L, degree, kind="simplicial"):
    """Tiny dense Schur oracle; explicit cell mapping removes prefix assumptions."""
    np = numpy()
    down, up_boundary, selected, outside = pair_blocks(K, L, degree, kind)
    up = up_boundary @ up_boundary.T
    A = up[np.ix_(selected, selected)]
    E = up[np.ix_(selected, outside)]
    C = up[np.ix_(outside, outside)]
    matrix = down.T @ down + A
    if outside:
        matrix = matrix - E @ np.linalg.pinv(C, rcond=RCOND, hermitian=True) @ E.T
    if not np.isfinite(matrix).all():
        raise ValueError("nonfinite persistent matrix")
    if not np.allclose(matrix, matrix.T, atol=COMPARISON_TOLERANCE, rtol=0):
        raise ValueError("persistent matrix lost symmetry")
    return (matrix + matrix.T) / 2


def nullspace_laplacian(K, L, degree, kind="simplicial"):
    """Independent restricted-boundary definition, without a Schur complement."""
    from scipy.linalg import null_space
    down, up, selected, outside = pair_blocks(K, L, degree, kind)
    restricted = up[selected, :] @ null_space(up[outside, :], rcond=RCOND)
    return down.T @ down + restricted @ restricted.T


def exact_rank(matrix):
    """Fraction Gaussian elimination on integer incidence; not a spectral rank."""
    rows = [[Fraction(int(value)) for value in row] for row in matrix]
    columns = matrix.shape[1]
    rank = 0
    for column in range(columns):
        pivot = next((i for i in range(rank, len(rows)) if rows[i][column]), None)
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        factor = rows[rank][column]
        rows[rank] = [value / factor for value in rows[rank]]
        for i in range(rank + 1, len(rows)):
            factor = rows[i][column]
            rows[i] = [a - factor * b for a, b in zip(rows[i], rows[rank])]
        rank += 1
    return rank


def exact_persistent_betti(K, L, degree, kind="simplicial"):
    down, up, selected, outside = pair_blocks(K, L, degree, kind)
    # dim(im B_L intersect C_q(K)) = rank(B_L) - rank(B_L outside).
    return len(selected) - exact_rank(down) - exact_rank(up) + exact_rank(up[outside, :])


def make_sample(sample_id, source_id, filtration, grid, kind="simplicial", ordinary=False):
    np = numpy()
    grid = checked_grid(grid)
    stages = [stage(filtration, value, kind) for value in grid]
    spectra = []
    for degree in (0, 1):
        for i in range(len(grid)):
            for j in range(i, len(grid)):
                matrix = persistent_laplacian(stages[i], stages[i] if ordinary else stages[j], degree, kind)
                spectra.append({"i": i, "j": j, "degree": degree,
                                "sample_id": sample_id, "source_id": source_id,
                                "basis": stages[i][degree], "domain_size": len(stages[i][degree]),
                                "coverage": "full", "status": "computed",
                                "values": np.linalg.eigvalsh(matrix).tolist()})
    return {"sample_id": sample_id, "source_id": source_id, "grid": list(grid),
            "operator": "ordinary_hodge" if ordinary else "persistent_laplacian",
            "field": "R", "inner_product": "identity", "scale_unit": "synthetic_scale",
            "eigenvalue_unit": "dimensionless", "complex_type": kind,
            "numeric_origin": {"producer": "independent_python_oracle", "dtype": "float64",
                               "solver": "numpy.linalg.eigvalsh", "schur_rcond": RCOND,
                               "zero_tolerance": ZERO_TOLERANCE}, "spectra": spectra}


def vectorize_spectra(samples, grid, degrees=(0, 1), length=8, operator="persistent_laplacian"):
    """Owned smallest-eigenvalue embedding; pad only a known exhausted domain.

    Columns are degree, start index, end index, eigenvalue rank, in that order.
    Full spectra must match the declared domain. Partial spectra must certify
    the smallest prefix and contain min(length, domain_size) values. Failed,
    absent, complex, nonfinite and incompatible input raises ValueError.
    """
    np = numpy()
    grid = checked_grid(grid)
    degrees = tuple(degrees)
    if (not degrees or any(type(q) is not int or q < 0 for q in degrees)
            or tuple(sorted(set(degrees))) != degrees
            or type(length) is not int or length < 1):
        raise ValueError("invalid degrees or truncation length")
    keys = [(q, i, j) for q in degrees for i in range(len(grid)) for j in range(i, len(grid))]
    columns = [{"degree": q, "i": i, "j": j, "start": grid[i], "end": grid[j], "rank": rank}
               for q, i, j in keys for rank in range(length)]
    values = np.zeros((len(samples), len(columns)))
    validity = np.zeros(values.shape, dtype=bool)
    domains, coverage, ids, sources, bases = [], [], [], [], []
    common = None
    for row, sample in enumerate(samples):
        sample_id, source_id = sample["sample_id"], sample["source_id"]
        if (not isinstance(sample_id, str) or not isinstance(source_id, str)
                or not sample_id or not source_id or sample_id in ids):
            raise ValueError("sample and source identity required; sample IDs must be unique")
        if tuple(sample["grid"]) != grid or sample["operator"] != operator:
            raise ValueError("scale grid or operator mismatch")
        if sample["field"] != "R" or sample["inner_product"] != "identity":
            raise ValueError("only real unweighted Euclidean spectra are supported")
        origin = sample["numeric_origin"]
        if (origin["dtype"] != "float64" or not origin["producer"] or not origin["solver"]
                or not math.isfinite(origin["zero_tolerance"]) or origin["zero_tolerance"] < 0):
            raise ValueError("invalid numeric origin")
        contract = (sample["field"], sample["inner_product"], sample["scale_unit"],
                    sample["eigenvalue_unit"], dict(origin))
        if common is not None and contract != common:
            raise ValueError("incompatible units or numeric origin across samples")
        common = contract
        records = {}
        for record in sample["spectra"]:
            key = (record["degree"], record["i"], record["j"])
            if key in records:
                raise ValueError("duplicate spectrum")
            records[key] = record
        if set(records) != set(keys):
            raise ValueError("missing or unexpected spectrum block")
        row_domains, row_coverage, row_bases = [], [], []
        stage_bases = {}
        for block, key in enumerate(keys):
            record = records[key]
            n = record["domain_size"]
            if (record["sample_id"] != sample_id or record["source_id"] != source_id
                    or record["status"] != "computed" or type(n) is not int or n < 0
                    or len(record["basis"]) != n
                    or len(set(tuple(cell) for cell in record["basis"])) != n):
                raise ValueError("invalid status, basis or spectrum identity")
            basis = tuple(tuple(cell) for cell in record["basis"])
            stage_key = key[:2]
            if stage_key in stage_bases and set(stage_bases[stage_key]) != set(basis):
                raise ValueError("inconsistent source-stage basis across end scales")
            stage_bases[stage_key] = basis
            raw = np.asarray(record["values"])
            if raw.ndim != 1 or np.iscomplexobj(raw):
                raise ValueError("spectrum must be a real one-dimensional array")
            spectrum = np.array(raw, dtype=float, copy=True)
            if not np.isfinite(spectrum).all() or np.any(spectrum < -origin["zero_tolerance"]):
                raise ValueError("nonfinite or negative spectrum")
            selected = min(length, n)
            if record["coverage"] == "full":
                if len(spectrum) != n:
                    raise ValueError("full spectrum length differs from domain size")
            elif record["coverage"] == "smallest":
                if not selected <= len(spectrum) <= n:
                    raise ValueError("insufficient smallest-prefix coverage; cannot pad missing values")
            else:
                raise ValueError("unsupported spectrum selection")
            spectrum[np.abs(spectrum) <= origin["zero_tolerance"]] = 0
            spectrum.sort()
            offset = block * length
            values[row, offset:offset + selected] = spectrum[:selected]
            validity[row, offset:offset + selected] = True
            row_domains.append(n)
            row_coverage.append(record["coverage"])
            row_bases.append(basis)
        domains.append(row_domains)
        coverage.append(row_coverage)
        bases.append(row_bases)
        ids.append(sample_id)
        sources.append(source_id)
    return {"values": values, "validity_mask": validity, "sample_ids": ids,
            "source_ids": sources, "columns": columns, "grid": list(grid),
            "degrees": list(degrees), "length": length, "operator": operator,
            "domain_sizes": domains, "spectrum_coverage": coverage, "bases": bases,
            "numeric_contract": common}


def gudhi_diagram(filtration, kind="simplicial", image=None, max_degree=1):
    import gudhi
    if kind == "simplicial":
        complex_ = gudhi.SimplexTree()
        for cell, value in filtration.items():
            complex_.insert(cell, filtration=value)
        result = complex_.persistence(homology_coeff_field=2, persistence_dim_max=True)
    else:
        complex_ = gudhi.CubicalComplex(vertices=numpy().asarray(image, dtype=float))
        result = complex_.persistence(homology_coeff_field=2)
    bars = [(int(q), float(b), None if math.isinf(d) else float(d))
            for q, (b, d) in result if q <= max_degree]
    return sorted(bars, key=lambda bar: (bar[0], bar[1], math.inf if bar[2] is None else bar[2]))


def author_functions(path):
    """Execute only unchanged pinned algebra functions, without author imports.

    The boundary adapter below is ours. This validates author matrix algebra,
    not the original image constructor, imports or MNIST pipeline.
    """
    data = path.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(data)).encode() + b"\0" + data).hexdigest()
    if blob != AUTHOR_BLOB:
        raise ValueError("author source differs from the fixed Git blob")
    tree = ast.parse(data.decode("utf-8"))
    chosen = [node for node in tree.body if isinstance(node, ast.FunctionDef)
              and node.name in ("schur_complement", "pers_lap_pair")]
    if len(chosen) != 2:
        raise ValueError("pinned author functions are absent")
    namespace = {"np": numpy(), "boundary": lambda cells, q, complex_type="simplicial":
                 boundary(cells, q, complex_type)}
    exec(compile(ast.Module(body=chosen, type_ignores=[]), str(path), "exec"), namespace)
    return namespace


def core_validation(author_source):
    np = numpy()
    author = author_functions(author_source)
    grid = sampling_grid(0, 1, 4)
    path = {(0,): 0., (2,): 0., (1,): 1., (0, 1): 1., (1, 2): 1.}
    filled = {**path, (0, 2): 1., (0, 1, 2): 1.}
    cycle = {(0,): 0., (1,): 0., (2,): 0., (0, 1): 0., (0, 2): 0.,
             (1, 2): 0., (0, 1, 2): 1.}
    image = np.array([[0., .5, 0.], [0., 1., 0.], [0., .5, 0.]])
    cases = [("path", path, "simplicial", None), ("filled", filled, "simplicial", None),
             ("singular-complement", {**path, (3,): 1.}, "simplicial", None),
             ("cycle-filled", cycle, "simplicial", None),
             ("cubical", cubical_filtration(image), "cubical", image)]
    cases.extend((f"disk-{s['label']}", s["filtration"], "simplicial", None)
                 for s in higher_order_samples(groups=1))
    sphere = {(v,): 0. for v in range(4)}
    sphere.update({(a, b): 0. for a in range(4) for b in range(a + 1, 4)})
    sphere.update({t: 1. for t in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3))})
    cases.append(("tetrahedron-boundary-h2", sphere, "simplicial", None))
    evidence = []
    maximum_error = 0.
    for name, filtration, kind, pixels in cases:
        stages = [stage(filtration, value, kind) for value in grid]
        diagram = gudhi_diagram(filtration, kind, pixels)
        for q in (0, 1):
            for i in range(len(grid)):
                for j in range(i, len(grid)):
                    K, L = stages[i], stages[j]
                    matrix = persistent_laplacian(K, L, q, kind)
                    direct = nullspace_laplacian(K, L, q, kind)
                    # Explicitly reorder each L basis with K as its prefix, as
                    # required by the unchanged author algebra function.
                    ordered = [K[d] + [cell for cell in L[d] if cell not in K[d]] for d in range(3)]
                    reference = author["pers_lap_pair"](K, ordered, q, complex_type=kind)
                    error = max(float(np.max(np.abs(matrix - direct), initial=0)),
                                float(np.max(np.abs(matrix - reference), initial=0)))
                    maximum_error = max(maximum_error, error)
                    eigenvalues = np.linalg.eigvalsh(matrix)
                    betti = exact_persistent_betti(K, L, q, kind)
                    nullity = int(np.count_nonzero(np.abs(eigenvalues) <= ZERO_TOLERANCE))
                    gudhi_betti = sum(degree == q and birth <= grid[i]
                                      and (death is None or death > grid[j])
                                      for degree, birth, death in diagram)
                    if (error > COMPARISON_TOLERANCE or nullity != betti or betti != gudhi_betti
                            or np.any(eigenvalues < -ZERO_TOLERANCE)):
                        raise AssertionError((name, q, i, j, error, nullity, betti, gudhi_betti))
                    evidence.append({"case": name, "degree": q, "i": i, "j": j,
                                     "basis": K[q], "matrix": matrix.tolist(),
                                     "eigenvalues": eigenvalues.tolist(), "real_exact_betti": betti,
                                     "gudhi_f2_betti": gudhi_betti, "matrix_error": error})
    K, L = stage(path, 0), stage(path, 1)
    np.testing.assert_allclose(persistent_laplacian(K, L, 0), [[.5, -.5], [-.5, .5]],
                               atol=COMPARISON_TOLERANCE, rtol=0)
    np.testing.assert_allclose(persistent_laplacian(K, stage(filled, 1), 0),
                               [[1.5, -1.5], [-1.5, 1.5]], atol=COMPARISON_TOLERANCE, rtol=0)
    if gudhi_diagram(path) != gudhi_diagram(filled):
        raise AssertionError("path/filled positive-length barcodes differ")
    unsafe = author["pers_lap_pair"](K, L, 0)
    sphere_end = stage(sphere, 1)
    exact_beta2 = len(sphere_end[2]) - exact_rank(boundary(sphere_end, 2))
    sphere_diagram = gudhi_diagram(sphere, max_degree=2)
    if exact_beta2 != 1 or (2, 1., None) not in sphere_diagram:
        raise AssertionError("tetrahedron boundary H2 check failed")
    return {"comparisons": len(evidence), "maximum_matrix_error": maximum_error,
            "cases": evidence, "path_filled_same_barcode": gudhi_diagram(path),
            "author_nonprefix_actual": unsafe.tolist(),
            "author_nonprefix_expected": persistent_laplacian(K, L, 0).tolist(),
            "sphere_h2": {"real_exact_betti2": exact_beta2, "gudhi_diagram": sphere_diagram,
                          "scope": "H2 topology coverage, no q=2 spectral computation"},
            "author_scope": "unchanged Schur/pair functions; our boundary and explicit prefix adapters"}


def synthetic_samples(seed=20261012, groups=64):
    """Paired trees with identical ordinary spectra/PH, different early vertices."""
    rng = random.Random(seed)
    samples = []
    seen = set()
    for group in range(groups):
        while True:
            n = rng.randrange(5, 11)
            edges = [(rng.randrange(v), v) for v in range(1, n)]
            signature = (n, tuple(edges))
            adjacency = [set() for _ in range(n)]
            for a, b in edges:
                adjacency[a].add(b)
                adjacency[b].add(a)
            far = []
            for a in range(n):
                distances, queue = {a: 0}, [a]
                for v in queue:
                    for u in adjacency[v] - distances.keys():
                        distances[u] = distances[v] + 1
                        queue.append(u)
                far.extend((a, b) for b, d in distances.items() if a < b and d >= 3)
            if far and signature not in seen:
                seen.add(signature)
                break
        for label, early in enumerate((rng.choice(edges), rng.choice(far))):
            filtration = {(v,): (0. if v in early else 1.) for v in range(n)}
            filtration.update({tuple(sorted(edge)): 1. for edge in edges})
            samples.append({"sample_id": f"tree-{group:03d}-{label}", "group": group,
                            "source_id": f"synthetic-tree-{group:03d}-early-{'-'.join(map(str, early))}",
                            "label": label, "filtration": filtration})
    return samples


def graph_pair_geometry(filtration):
    """BFS distance and grounded graph solve, independently of Schur spectra."""
    np = numpy()
    vertices = sorted(cell[0] for cell in filtration if len(cell) == 1)
    earliest = min(filtration[(v,)] for v in vertices)
    terminals = [v for v in vertices if filtration[(v,)] == earliest]
    if len(terminals) != 2:
        raise ValueError("geometric baseline requires exactly two earliest vertices")
    adjacency = {v: set() for v in vertices}
    laplacian = np.zeros((len(vertices), len(vertices)))
    indices = {v: i for i, v in enumerate(vertices)}
    for a, b in (cell for cell in filtration if len(cell) == 2):
        adjacency[a].add(b)
        adjacency[b].add(a)
        i, j = indices[a], indices[b]
        laplacian[i, i] += 1
        laplacian[j, j] += 1
        laplacian[i, j] -= 1
        laplacian[j, i] -= 1
    source, ground = terminals
    distances, queue = {source: 0}, [source]
    for v in queue:
        for u in sorted(adjacency[v] - distances.keys()):
            distances[u] = distances[v] + 1
            queue.append(u)
    if len(distances) != len(vertices):
        raise ValueError("geometric baseline requires a connected full graph")
    kept = [i for i, v in enumerate(vertices) if v != ground]
    current = np.zeros(len(kept))
    current[kept.index(indices[source])] = 1
    potential = np.linalg.solve(laplacian[np.ix_(kept, kept)], current)
    return float(distances[ground]), float(potential[kept.index(indices[source])])


def higher_order_samples(seed=20261012, groups=64):
    """One-face versus three-face cycle fillings, with matched nuisance trees.

    L contains faces 012,013,023 of a tetrahedron boundary. K contains the
    012 cycle (label 0) or 123 cycle (label 1); all four vertices are early.
    Separate early trees remain disjoint, identically per pair. The common
    disk motif occurs in both splits; this is not out-of-family generalization.
    """
    rng = random.Random(seed)
    result, seen = [], set()
    for group in range(groups):
        if groups == 1:
            tree, n = [], 0
        else:
            while True:
                n = rng.randrange(3, 7)
                tree = [(4 + rng.randrange(v), 4 + v) for v in range(1, n)]
                signature = (n, tuple(tree))
                if signature not in seen:
                    seen.add(signature)
                    break
        late = {(v,): 0. for v in range(4 + n)}
        late.update({(a, b): 1. for a in range(4) for b in range(a + 1, 4)})
        late.update({t: 1. for t in ((0, 1, 2), (0, 1, 3), (0, 2, 3))})
        late.update({edge: 0. for edge in tree})
        for label, edges in enumerate((((0, 1), (0, 2), (1, 2)), ((1, 2), (1, 3), (2, 3)))):
            result.append({"sample_id": f"disk-{group:03d}-{label}", "group": group,
                           "source_id": f"disk-tree-{group:03d}-cycle-{label}", "label": label,
                           "filtration": {**late, **{edge: 0. for edge in edges}}})
    return result


def author_feature_audit(path):
    """Run unchanged feature-loop code with trace callbacks, without data imports."""
    np = numpy()
    author_functions(path)
    tree = ast.parse(path.read_text(encoding="utf-8"))
    chosen = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "features_perslap"]
    pairs = []

    def pair_callback(K, L, **kwargs):
        pairs.append((float(K), float(L)))
        return [], np.array([[K, L], [K, L]])

    namespace = {"np": np, "cubical_complex": lambda filtration, cutoff: cutoff,
                 "compute_pers_lap_pair": pair_callback}
    exec(compile(ast.Module(body=chosen, type_ignores=[]), str(path), "exec"), namespace)
    values, metadata = namespace["features_perslap"](np.array([0., 1.]), resolution=4)
    try:
        namespace["features_perslap"](np.array([0., 1.]), resolution=4, aggregate_function=np.sum)
    except TypeError as error:
        aggregate = {"status": "failed_as_expected", "type": type(error).__name__, "message": str(error)}
    else:
        raise AssertionError("pinned author aggregation behavior changed")
    return {"scope": "unchanged feature-loop with constructor/eigen callbacks; not original end-to-end pipeline",
            "sample_points": sorted(set(a for a, _ in pairs)), "pair_trace": pairs,
            "returned_spectrum_blocks": len(values), "returned_grid_metadata": metadata,
            "aggregate_branch": aggregate, "our_four_point_grid": list(sampling_grid(0, 1, 4)),
            "source_definition": "author function returns spectral collection; our sort/truncate/pad/flatten follows paper"}


def rust_baselines(samples, grid, output):
    """Run the existing Rust Betti/lifetime descriptors on GUDHI intervals."""
    build = subprocess.run(["cargo", "build", "--locked", "--lib", "--message-format=json",
                            "--target-dir", str(output / "rust-build")], cwd=ROOT,
                           capture_output=True, text=True)
    (output / "rust-build.log").write_text(build.stdout + build.stderr, encoding="utf-8")
    build.check_returncode()
    artifacts = [json.loads(line) for line in build.stdout.splitlines() if line.startswith("{")]
    library = next(filename for item in artifacts if item.get("reason") == "compiler-artifact"
                   and item["target"]["name"] == "cocycle" for filename in item["filenames"]
                   if filename.endswith(".rlib"))
    worker = output / ("baseline.exe" if sys.platform == "win32" else "baseline")
    command = ["rustc", "--edition", "2024", "-D", "warnings",
               str(ROOT / "tools/reference/persistent_spectral_baseline.rs"),
               "--extern", f"cocycle={library}", "-L", str(Path(library).parent), "-o", str(worker)]
    compiled = subprocess.run(command, cwd=ROOT, capture_output=True, text=True)
    (output / "rust-worker-build.log").write_text(compiled.stdout + compiled.stderr, encoding="utf-8")
    compiled.check_returncode()
    test_worker = worker.with_name(worker.stem + "-tests" + worker.suffix)
    tested = subprocess.run(command[:-2] + ["--test", "-o", str(test_worker)],
                            cwd=ROOT, capture_output=True, text=True, check=True)
    tests = subprocess.run([str(test_worker)], capture_output=True, text=True, check=True)
    (output / "rust-worker-tests.log").write_text(tested.stdout + tested.stderr + tests.stdout + tests.stderr,
                                                encoding="utf-8")
    rows = []
    for sample in samples:
        diagram = gudhi_diagram(sample["filtration"])
        request = " ".join(map(str, grid)) + "\n" + "\n".join(
            f"{q} {birth} {'inf' if death is None else death}" for q, birth, death in diagram) + "\n"
        result = subprocess.run([str(worker)], input=request, capture_output=True, text=True, check=True)
        rows.append(json.loads(result.stdout))
    return numpy().array(rows), {"worker_source": "tools/reference/persistent_spectral_baseline.rs",
                                 "binary_path": str(worker),
                                 "library_sha256": file_hash(Path(library)),
                                 "binary_sha256": file_hash(worker),
                                 "interval_source": "GUDHI F2 on torsion-free worked examples",
                                 "columns_per_degree": "four Betti values, count, total, max, entropy, essential, max_valid, entropy_valid"}


def rust_request(sample, grid):
    """The numerical worker receives no labels, split or oracle values."""
    # Input IDs are separately checked in Rust before JSON serialization.
    for key in ("sample_id", "source_id"):
        value = sample[key]
        if (not isinstance(value, str) or not 1 <= len(value) <= 128
                or any(not (c.isascii() and (c.isalnum() or c in "-_.")) for c in value)):
            raise ValueError("invalid worker identity")
    lines = [sample["sample_id"], sample["source_id"], " ".join(map(str, checked_grid(grid))),
             str(len(sample["filtration"]))]
    lines.extend(f"{birth} {len(cell)} " + " ".join(map(str, cell))
                 for cell, birth in sample["filtration"].items())
    return "\n".join(lines) + "\n"


def verify_rust_sample(record, sample, grid):
    """Check full source/pair identities and independent matrices, not just spectra."""
    np = numpy()
    if record["sample_id"] != sample["sample_id"] or record["source_id"] != sample["source_id"]:
        raise ValueError("Rust sample/source identity mismatch")
    if record["numeric_origin"]["producer"] != "cocycle.persistent_laplacian":
        raise ValueError("expected actual Rust persistent Laplacian results")
    embedding = vectorize_spectra([record], grid)
    np.testing.assert_array_equal(record["validity_mask"], embedding["validity_mask"][0])
    np.testing.assert_array_equal(record["features"], embedding["values"][0])
    maximum = 0.
    for r in record["spectra"]:
        q, i, j = r["degree"], r["i"], r["j"]
        K = stage(sample["filtration"], grid[i])
        L = stage(sample["filtration"], grid[j])
        basis = [tuple(cell) for cell in r["basis"]]
        if (set(basis) != set(K[q]) or len(r["target_basis"]) != len(L[q])
                or set(map(tuple, r["target_basis"])) != set(L[q])
                or r["basis_births"] != [sample["filtration"][cell] for cell in basis]):
            raise ValueError("Rust source or target basis differs from supplied stages")
        order = [K[q].index(cell) for cell in basis]
        expected = nullspace_laplacian(K, L, q)[np.ix_(order, order)]
        actual = np.asarray(r["matrix"]).reshape((len(basis), len(basis)))
        np.testing.assert_allclose(actual, expected, atol=COMPARISON_TOLERANCE, rtol=0)
        np.testing.assert_allclose(r["values"], np.linalg.eigvalsh(expected),
                                   atol=COMPARISON_TOLERANCE, rtol=0)
        maximum = max(maximum, float(np.max(np.abs(actual - expected), initial=0)))
    return maximum


def rust_spectral_features(samples, grid, worker, output):
    from compare_persistent_laplacians import peak_memory
    records, resources = [], []
    maximum = 0.
    raw = output / "rust-spectral"
    raw.mkdir()
    for index, sample in enumerate(samples):
        request = rust_request(sample, grid)
        started = time.perf_counter()
        child = subprocess.Popen([str(worker), "--spectral"], stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        try:
            stdout, stderr = child.communicate(request, timeout=30)
        except subprocess.TimeoutExpired:
            child.kill()
            stdout, stderr = child.communicate()
            (raw / f"{index:03d}-timeout.txt").write_text(stdout + stderr, encoding="utf-8")
            raise RuntimeError("Rust spectral worker timeout")
        process_seconds = time.perf_counter() - started
        (raw / f"{index:03d}-input.txt").write_text(request, encoding="utf-8")
        (raw / f"{index:03d}-stdout.json").write_text(stdout, encoding="utf-8")
        (raw / f"{index:03d}-stderr.txt").write_text(stderr, encoding="utf-8")
        if child.returncode:
            raise RuntimeError(f"Rust export failed: {stderr}")
        resources.append({"process_seconds": process_seconds,
                          "peak_working_set_bytes": peak_memory(child)})
        record = json.loads(stdout)
        maximum = max(maximum, verify_rust_sample(record, sample, grid))
        records.append(record)
    embedding = vectorize_spectra(records, grid)
    embedding["values"] = numpy().array([record["features"] for record in records])
    embedding["validity_mask"] = numpy().array([record["validity_mask"] for record in records], dtype=bool)
    return embedding, {
        "producer": "cocycle.persistent_laplacian", "samples": len(records),
        "pairs": len(records) * 20, "maximum_matrix_error": maximum, "processes": resources,
        "resource_boundary": "one child per sample; includes startup/export; Windows peak working set includes baseline; null elsewhere",
        "limits": "64 KiB input, 256 simplices, dimensions 0..2, 4096 entries/buffer, 20M work per construction and per solve; 30s/child"}


def fit_exported_features(input_path, output_path):
    """Separate model process; only exported arrays and sample/split IDs enter."""
    np = numpy()
    from sklearn.ensemble import RandomForestClassifier
    started = time.perf_counter()
    cpu = time.process_time()
    if input_path.stat().st_size > 4 * 1024 * 1024:
        raise ValueError("model archive exceeds 4 MiB")
    with np.load(input_path, allow_pickle=False) as data:
        ids, labels = data["sample_ids"], data["labels"]
        train, test = data["train"], data["test"]
        n = len(ids)
        if (n == 0 or n > 256 or ids.ndim != 1 or ids.dtype.kind != "U"
                or not all(ids.tolist()) or len(set(ids.tolist())) != n
                or labels.shape != (n,) or labels.dtype.kind not in "iu"
                or set(labels.tolist()) != {0, 1}
                or train.ndim != 1 or test.ndim != 1 or not len(train) or not len(test)
                or train.dtype.kind not in "iu" or test.dtype.kind not in "iu"
                or sorted(np.concatenate([train, test]).tolist()) != list(range(n))):
            raise ValueError("invalid exported sample/split manifest")
        if data["groups"].shape != (n,):
            raise ValueError("invalid group identity shape")
        if set(data["groups"][train].tolist()) & set(data["groups"][test].tolist()):
            raise ValueError("group leakage in exported split")
        branches = data["branches"].tolist()
        if (data["branches"].ndim != 1 or not 1 <= len(branches) <= 16
                or len(set(branches)) != len(branches)):
            raise ValueError("invalid feature branches")
        scores = {}
        for name in branches:
            features = data[name]
            if (features.ndim != 2 or features.shape[0] != n or not 0 < features.shape[1] <= 512
                    or features.dtype.kind not in "fbiu" or not np.isfinite(features).all()):
                raise ValueError("invalid exported feature matrix")
            params = {"n_estimators": 100, "max_depth": 4, "random_state": int(data["seed"]), "n_jobs": 1}
            begin = time.perf_counter()
            model = RandomForestClassifier(**params).fit(features[train], labels[train])
            predictions = model.predict(features[test])
            scores[name] = {"accuracy": float(np.mean(predictions == labels[test])),
                            "shape": list(features.shape), "predictions": predictions.tolist(),
                            "fit_predict_seconds": time.perf_counter() - begin,
                            "feature_bytes": features.nbytes}
        write_json(output_path, {"scores": scores, "model": params,
                   "test_sample_ids": ids[test].tolist(), "input_sha256": file_hash(input_path),
                   "wall_seconds": time.perf_counter() - started, "cpu_seconds": time.process_time() - cpu,
                   "gudhi_imported": "gudhi" in sys.modules})


def learning_run(grid, output, seed, samples=None):
    np = numpy()
    from sklearn.metrics import accuracy_score
    from sklearn.model_selection import train_test_split
    tree_case = samples is None
    samples = synthetic_samples(seed) if tree_case else samples
    started = time.perf_counter()
    persistent = [make_sample(s["sample_id"], s["source_id"], s["filtration"], grid) for s in samples]
    persistent_time = time.perf_counter() - started
    started = time.perf_counter()
    ordinary = [make_sample(s["sample_id"], s["source_id"], s["filtration"], grid, ordinary=True) for s in samples]
    ordinary_time = time.perf_counter() - started
    embedding = vectorize_spectra(persistent, grid)
    hodge = vectorize_spectra(ordinary, grid, operator="ordinary_hodge")
    # The graph baseline keeps degree zero only, with the same pair columns.
    ordinary_zero = [{**s, "spectra": [r for r in s["spectra"] if r["degree"] == 0]} for s in ordinary]
    graph = vectorize_spectra(ordinary_zero, grid, degrees=(0,), operator="ordinary_hodge")
    geometry = np.array([graph_pair_geometry(s["filtration"]) for s in samples]) if tree_case else None
    started = time.perf_counter()
    ph, rust = rust_baselines(samples, grid, output)
    baseline_time = time.perf_counter() - started
    native, native_evidence = rust_spectral_features(samples, grid, Path(rust["binary_path"]), output)
    native_evidence["maximum_feature_error"] = float(np.max(np.abs(native["values"] - embedding["values"])))
    np.testing.assert_allclose(native["values"], embedding["values"], atol=COMPARISON_TOLERANCE, rtol=0)
    np.testing.assert_array_equal(native["validity_mask"], embedding["validity_mask"])
    # Model input is the Rust export. The independent Python oracle is a control.
    q0 = native["values"][:, [c["degree"] == 0 for c in native["columns"]]]
    q1 = native["values"][:, [c["degree"] == 1 for c in native["columns"]]]
    for index in range(0, len(samples), 2):
        np.testing.assert_allclose(hodge["values"][index], hodge["values"][index + 1],
                                   atol=COMPARISON_TOLERANCE, rtol=0)
        np.testing.assert_array_equal(ph[index], ph[index + 1])
        if not tree_case:
            np.testing.assert_allclose(q0[index], q0[index + 1], atol=COMPARISON_TOLERANCE, rtol=0)
        if np.allclose(embedding["values"][index], embedding["values"][index + 1],
                       atol=COMPARISON_TOLERANCE, rtol=0):
            raise AssertionError("persistent features failed to distinguish the paired construction")
    if tree_case:
        np.testing.assert_allclose(geometry[:, 0], geometry[:, 1], atol=COMPARISON_TOLERANCE, rtol=0)
        for row, sample in enumerate(samples):
            K, L = stage(sample["filtration"], grid[0]), stage(sample["filtration"], grid[-1])
            nonzero = np.linalg.eigvalsh(persistent_laplacian(K, L, 0))[-1]
            np.testing.assert_allclose(nonzero * geometry[row, 1], 2., atol=COMPARISON_TOLERANCE, rtol=0)
    train_groups, test_groups = train_test_split(sorted({s["group"] for s in samples}), test_size=.25, random_state=seed)
    train = [i for i, s in enumerate(samples) if s["group"] in train_groups]
    test = [i for i, s in enumerate(samples) if s["group"] in test_groups]
    labels = np.array([s["label"] for s in samples])
    params = {"n_estimators": 100, "max_depth": 4, "random_state": seed, "n_jobs": 1}
    branches = [("betti_lifetime", ph), ("ordinary_graph", graph["values"]),
                ("ordinary_hodge", hodge["values"]), ("persistent_q0", q0),
                ("persistent_q1", q1), ("persistent_spectral", native["values"]),
                ("persistent_with_mask", np.concatenate([native["values"], native["validity_mask"]], axis=1)),
                ("mask_only", native["validity_mask"]), ("domain_size_only", np.array(native["domain_sizes"])),
                ("oracle_persistent", embedding["values"])]
    if tree_case:
        branches.extend((("shortest_path", geometry[:, :1]), ("effective_resistance", geometry[:, 1:])))
    model_input = output / "model-input.npz"
    np.savez(model_input, **dict(branches), labels=labels, sample_ids=np.array(native["sample_ids"]),
             groups=np.array([s["group"] for s in samples]), train=np.array(train), test=np.array(test),
             branches=np.array([name for name, _ in branches]), seed=np.array(seed))
    started = time.perf_counter()
    from compare_persistent_laplacians import peak_memory
    child = subprocess.Popen([sys.executable, str(Path(__file__).resolve()), "--fit-exported",
                              str(model_input), str(output / "model-output.json")],
                             stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        stdout, stderr = child.communicate(timeout=120)
    except subprocess.TimeoutExpired:
        child.kill()
        child.communicate()
        raise RuntimeError("external model timeout")
    model_process_seconds = time.perf_counter() - started
    (output / "model-process.log").write_text(stdout + stderr, encoding="utf-8")
    if child.returncode:
        raise RuntimeError(f"external model failed: {stderr}")
    learned = json.loads((output / "model-output.json").read_text(encoding="utf-8"))
    if (learned["input_sha256"] != file_hash(model_input)
            or learned["test_sample_ids"] != [samples[i]["sample_id"] for i in test]
            or set(learned["scores"]) != {name for name, _ in branches}):
        raise ValueError("external model output identity mismatch")
    scores = learned["scores"]
    for name, features in branches:
        score = scores[name]
        prediction = np.asarray(score["predictions"])
        if (prediction.shape != labels[test].shape or not np.isin(prediction, [0, 1]).all()
                or score["shape"] != list(features.shape)
                or score["accuracy"] != float(accuracy_score(labels[test], prediction))):
            raise ValueError("invalid model prediction readback")
    arrays = {"persistent": native["values"], "validity_mask": native["validity_mask"],
              "ordinary_hodge": hodge["values"], "ordinary_graph": graph["values"],
              "betti_lifetime": ph, "persistent_q0": q0, "persistent_q1": q1, "labels": labels}
    if tree_case:
        arrays.update(shortest_path=geometry[:, :1], effective_resistance=geometry[:, 1:])
    np.savez(output / "features.npz", **arrays)
    write_json(output / "spectra.json", persistent)
    schema = {key: value for key, value in native.items() if key not in ("values", "validity_mask")}
    write_json(output / "feature_schema.json", schema)
    write_json(output / "fixtures.json", [{**s, "filtration": [[list(cell), birth] for cell, birth in s["filtration"].items()]}
                                         for s in samples])
    return {"scores": scores, "model": params, "train_rows": train, "test_rows": test,
            "test_labels": labels[test].tolist(), "train_groups": train_groups, "test_groups": test_groups,
            "rust_baseline": rust, "persistent_seconds": persistent_time, "ordinary_seconds": ordinary_time,
            "rust_spectral": native_evidence, "model_process_seconds": model_process_seconds,
            "model_cpu_seconds": learned["cpu_seconds"],
            "model_peak_working_set_bytes": peak_memory(child),
            "baseline_seconds_including_build_and_processes": baseline_time,
            "mask_bytes": embedding["validity_mask"].nbytes,
            "copy_policy": "owned zero-filled matrices and mask; each input spectrum copied before sorting",
            "scope": "paired trees" if tree_case else "paired disk cycles with disjoint nuisance trees; graph-isomorphic filtration pairs; shared motif, not out-of-family generalization",
            "geometric_baseline": "BFS and grounded full-graph solve; label defined by distance" if tree_case else "not applicable: all vertices early; q0 spectra matched",
            "generalization_claim": False}


def file_hash(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8", newline="\n")


def source_provenance():
    """Accept a real checkout or a fully hash-verified offline source snapshot."""
    try:
        directory = subprocess.check_output(["git", "rev-parse", "--show-toplevel"], cwd=ROOT,
                                            stderr=subprocess.DEVNULL, text=True).strip()
    except (FileNotFoundError, subprocess.CalledProcessError):
        directory = ""
    if directory and Path(directory).resolve() == ROOT.resolve():
        return (subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
                subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True))
    snapshot = json.loads((ROOT / "SOURCE_SNAPSHOT.json").read_text(encoding="utf-8"))
    for relative, expected in snapshot["sources"].items():
        path = (ROOT / relative).resolve()
        path.relative_to(ROOT.resolve())
        if file_hash(path) != expected["sha256"]:
            raise ValueError(f"source snapshot hash mismatch: {relative}")
    return snapshot["base_sha"], "verified source snapshot; no Git checkout"


def main():
    if len(sys.argv) == 4 and sys.argv[1] == "--fit-exported":
        fit_exported_features(Path(sys.argv[2]), Path(sys.argv[3]))
        return 0
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--author-source", type=Path, required=True, help="unmodified pinned src/perslap.py")
    parser.add_argument("--output", type=Path, required=True, help="new output directory; never overwritten")
    parser.add_argument("--seed", type=int, default=20261012)
    args = parser.parse_args()
    import importlib.metadata
    from threadpoolctl import threadpool_info, threadpool_limits
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    shutil.copyfile(args.author_source, output / "author-perslap.py")
    revision, dirty = source_provenance()
    environment = {"python": sys.version, "platform": platform.platform(), "processor": platform.processor(),
                   "versions": {name: importlib.metadata.version(name) for name in
                                ("numpy", "scipy", "gudhi", "scikit-learn", "threadpoolctl")},
                   "base_sha": revision, "git_status": dirty,
                   "sources": {path.relative_to(ROOT).as_posix(): file_hash(path) for path in
                               (Path(__file__), ROOT / "tools/reference/persistent_spectral_baseline.rs",
                                ROOT / "tools/compare_persistent_laplacians.py", ROOT / "Cargo.toml",
                                ROOT / "Cargo.lock", *sorted((ROOT / "src").rglob("*.rs")))},
                   "author_revision": AUTHOR_REVISION, "author_blob": AUTHOR_BLOB,
                   "author_sha256": file_hash(args.author_source), "command": sys.argv,
                   "toolchains": {name: subprocess.check_output([name, "--version"], text=True).strip()
                                  for name in ("cargo", "rustc")},
                   "field": "R", "inner_product": "identity", "schur_rcond": RCOND,
                   "zero_tolerance": ZERO_TOLERANCE, "absolute_comparison_tolerance": COMPARISON_TOLERANCE,
                   "seed": args.seed, "blas_backends": threadpool_info(),
                   "evidence": "bounded-check and numerical-evidence"}
    write_json(output / "environment.json", environment)
    grid = sampling_grid(0, 1, 4)
    tracemalloc.start()
    wall = time.perf_counter()
    cpu = time.process_time()
    try:
        with threadpool_limits(limits=1):
            core = core_validation(output / "author-perslap.py")
            write_json(output / "core-validation.json", core)
            write_json(output / "author-feature-audit.json", author_feature_audit(output / "author-perslap.py"))
            learning = learning_run(grid, output, args.seed)
            higher_output = output / "higher-order"
            higher_output.mkdir()
            higher = learning_run(grid, higher_output, args.seed, higher_order_samples(args.seed))
            write_json(higher_output / "results.json", higher)
        current, peak = tracemalloc.get_traced_memory()
        summary = {"status": "passed", "core_comparisons": core["comparisons"],
                   "maximum_matrix_error": core["maximum_matrix_error"], "learning": learning, "higher_order_learning": higher,
                   "wall_seconds": time.perf_counter() - wall, "parent_cpu_seconds": time.process_time() - cpu,
                   "python_tracemalloc_current_bytes": current, "python_tracemalloc_peak_bytes": peak,
                   "resource_limits": "dense tiny oracle; tracemalloc is not OS peak RSS; CPU excludes Rust children"}
        write_json(output / "results.json", summary)
        evidence_files = [p for p in output.rglob("*") if p.is_file() and "rust-build" not in p.relative_to(output).parts]
        write_json(output / "manifest.json", {p.relative_to(output).as_posix(): {"sha256": file_hash(p), "bytes": p.stat().st_size}
                                              for p in evidence_files})
        print(json.dumps({"status": summary["status"], "core_comparisons": core["comparisons"],
                          "scores": {k: v["accuracy"] for k, v in learning["scores"].items()},
                          "higher_order_scores": {k: v["accuracy"] for k, v in higher["scores"].items()}, "output": str(output)}))
    except Exception as error:
        write_json(output / "failure.json", {"status": "failed", "type": type(error).__name__, "message": str(error)})
        raise
    finally:
        tracemalloc.stop()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
