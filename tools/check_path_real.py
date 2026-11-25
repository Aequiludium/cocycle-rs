"""Independent rational/SVD verification of native real path Laplacian action.

Only this standalone verifier uses NumPy/SciPy/SymPy. Runtime Rust has no such
dependencies. Paths, deletion boundaries, rational kernels and stage-pair ranks
are reconstructed independently of Rust path enumeration and numerical kernels.
"""
import argparse
import hashlib
import itertools
import json
import platform
from pathlib import Path
import subprocess
import sys

import numpy as np
import scipy
from scipy.linalg import svd
import sympy as sp

ROOT = Path(__file__).resolve().parents[1]
TOLERANCE = 1e-8


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def columns(vectors, rows):
    return sp.Matrix.hstack(*vectors) if vectors else sp.zeros(rows, 0)


def orthogonal(matrix):
    if not matrix.cols:
        return np.zeros((matrix.rows, 0))
    u, values, _ = svd(np.array(matrix, dtype=float), full_matrices=False)
    if min(values) <= 1e-12:
        raise ValueError("oracle rational basis has ambiguous numerical conditioning")
    return u[:, :matrix.cols]


def deletion(word, regular):
    result = {}
    if len(word) > 1:
        for j in range(len(word)):
            face = word[:j] + word[j+1:]
            if regular and any(a == b for a, b in zip(face, face[1:])):
                continue
            result[face] = result.get(face, 0) + (-1)**j
    return {k: v for k, v in result.items() if v}


def prepare(n, edges, weights, d, regular):
    if n**(d+2) > 200000:
        raise ValueError("bounded independent oracle permits at most 200000 vertex words per layer")
    births = dict(zip(map(tuple, edges), weights))
    paths = [[w for w in itertools.product(range(n), repeat=q+1)
              if all(e in births for e in zip(w, w[1:]))] for q in range(d+2)]
    words = [[deletion(w, regular) for w in layer] for layer in paths]
    boundary = []
    for q, layer in enumerate(paths):
        lower = paths[q-1] if q else []
        boundary.append(sp.Matrix(len(lower), len(layer), lambda i, j: words[q][j].get(lower[i], 0)))
    times = sorted({0, *weights})
    stages = []
    for time in times:
        active = [[i for i, w in enumerate(layer)
                   if max((births[e] for e in zip(w, w[1:])), default=0) <= time] for layer in paths]
        bases, cycles = [], []
        for q, layer in enumerate(paths):
            lower = {paths[q-1][i] for i in active[q-1]} if q else set()
            faces = sorted({f for j in active[q] for f in words[q][j]})
            forbidden = [f for f in faces if f not in lower]

            def kernel_on(rows):
                matrix = sp.Matrix(len(rows), len(active[q]), lambda i, j: words[q][active[q][j]].get(rows[i], 0))
                local = matrix.nullspace()
                embedded = sp.zeros(len(layer), len(local))
                for j, vector in enumerate(local):
                    for k, i in enumerate(active[q]):
                        embedded[i, j] = vector[k]
                return embedded

            bases.append(kernel_on(forbidden))
            cycles.append(kernel_on(faces))  # Direct full ambient kernel, not Omega reduction.
        stages.append({"S": bases, "Q": [orthogonal(s) for s in bases], "Z": cycles})
    return paths, boundary, times, stages


def reference(boundary, source, target, q):
    s, upper = source["S"][q], target["S"][q+1]
    rows = s.rows
    projector = s * (s.T*s).inv() * s.T if s.cols else sp.zeros(rows, rows)
    image = boundary[q+1] * upper
    constraint = (sp.eye(rows)-projector) * image
    domain = upper * columns(constraint.nullspace(), upper.cols)
    qs, qu = source["Q"][q], orthogonal(domain)
    down = np.array(boundary[q], dtype=float) @ qs
    up = qs.T @ np.array(boundary[q+1], dtype=float) @ qu
    laplacian = down.T @ down + up @ up.T
    # A separate floating SVD domain uses orthonormal target upper coordinates.
    y = np.array(boundary[q+1], dtype=float) @ target["Q"][q+1]
    residual = y-qs @ (qs.T @ y)
    if y.shape[1]:
        _, singular, vh = svd(residual, full_matrices=True)
        rank = int(np.count_nonzero(singular > 1e-10*max(1, np.linalg.norm(y))))
        numerical_up = qs.T @ y @ vh[rank:].T
    else:
        numerical_up = np.zeros((qs.shape[1], 0))
    numerical_laplacian = down.T @ down + numerical_up @ numerical_up.T
    z = source["Z"][q]
    rho = int(z.row_join(image).rank()-image.rank())
    return laplacian, numerical_laplacian, rho, domain.cols


def cases(quick):
    all_edges = list(itertools.permutations(range(3), 2))
    for mask in ([0, 1, 3, 63] if quick else range(64)):
        edges = [e for i, e in enumerate(all_edges) if mask & (1 << i)]
        yield f"static-{mask}", 3, edges, [0]*len(edges), 1 if quick else 2
    # Fixed, hand-maintained finite filtration cases; no research checkout imports.
    filtered = [
        ("delayed-triangle", 3, [(0,1),(1,2),(2,0),(0,2)], [1,1,1,2]),
        ("delayed-diamond", 4, [(0,1),(0,2),(1,3),(2,3),(0,3)], [1,1,1,2,3]),
        ("filled-square", 4, [(0,1),(1,2),(2,3),(0,3),(0,2)], [1,1,1,1,2]),
        ("essential-cycle", 4, [(0,1),(1,2),(2,3),(3,0)], [1,2,2,3]),
        ("reciprocal-isolate", 3, [(0,1),(1,0)], [1,2]),
        ("seeded-0", 4, [(0,2),(1,0),(1,2),(1,3),(2,3)], [3,3,2,3,3]),
        ("seeded-1", 5, [(0,1),(1,0),(1,2),(1,3),(2,0),(2,1),(2,3),(2,4),(3,1),(4,2),(4,3)], [2,1,1,2,2,1,3,2,1,3,3]),
        ("seeded-2", 4, [(0,2),(1,0),(2,3)], [3,2,1]),
    ]
    for name, n, edges, weights in filtered[:1] if quick else filtered:
        yield name, n, edges, weights, 1
    yield "empty", 0, [], [], 1
    if not quick:
        yield "static-diamond-metric", 4, [(0,1),(0,2),(1,3),(2,3)], [0]*4, 1
        edges = [(0,1),(0,2),(0,3),(1,4),(2,5),(3,6),(1,5),(2,6),(3,4),
                 (4,7),(5,8),(6,9),(4,8),(5,9),(6,7),(7,10),(8,10),(9,10),
                 (1,7),(2,8),(3,9),(1,9),(2,7),(3,8)]
        yield "characteristic-sensitive", 11, edges, [0]*len(edges), 3


def check(worker, values, regular):
    result = subprocess.run([str(worker), *([] if regular else ["--non-regular"])],
                            input=" ".join(map(str, values)), capture_output=True,
                            text=True, encoding="utf-8", timeout=30)
    if result.returncode:
        raise ValueError(result.stderr)
    return json.loads(result.stdout)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    identities = [Path(__file__), ROOT/"tools/reference/path_real.rs", *ROOT.glob("src/**/*.rs"),
                  ROOT/"Cargo.toml", ROOT/"Cargo.lock", ROOT/"benches/path/requirements-lock.txt"]
    before = {p.relative_to(ROOT).as_posix(): sha(p) for p in identities}
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    summary = {"protocol": "cocycle-path-real-reference-v1", "evidence": "numerical-evidence",
               "commit": head, "sources": before, "worker_sha256": sha(args.worker),
               "versions": {"python": platform.python_version(), "numpy": np.__version__, "scipy": scipy.__version__, "sympy": sp.__version__},
               "domain": "rational Omega and exact real persistent ranks; f64 orthonormal metric/eigenvalues",
               "tolerance": TOLERANCE, "cases": 0, "certificates": 0, "max_errors": {}, "status": "failed"}
    try:
        with (args.output/"certificates.jsonl").open("w", encoding="utf-8") as records:
            for name, n, edges, weights, d in cases(args.quick):
                for regular in (True,) if name == "characteristic-sensitive" else (True, False):
                    paths, boundary, times, stages = prepare(n, edges, weights, d, regular)
                    summary["cases"] += 1
                    for s in range(len(times)):
                        for t in range(s, len(times)):
                            for q in range(d+1):
                                values = [n, q, 10000, 1000000, 1e-10, s, t, *(v for e, w in zip(edges, weights) for v in (*e, w))]
                                result = check(args.worker, values, regular)
                                if result.get("protocol") != "cocycle-path-real-certificate-v1" or result["field"] != "real" or result["paths"] != list(map(list, paths[q])) or result["dimension"] != q or result["regular"] != regular:
                                    raise ValueError("native coordinate/model identity mismatch")
                                if (result["source_time"], result["target_time"]) != (times[s], times[t]):
                                    raise ValueError("native event identity mismatch")
                                size = len(result["basis"])
                                if size != stages[s]["S"][q].cols:
                                    raise ValueError("native real Omega dimension mismatch")
                                basis = np.array(result["basis"], dtype=float).reshape(size, len(paths[q])).T
                                actual = np.array(result["laplacian_columns"], dtype=float).reshape(size, size).T
                                expected, numerical, rho, domain = reference(boundary, stages[s], stages[t], q)
                                if name == "static-diamond-metric" and q == 1 and regular:
                                    if not np.allclose(np.linalg.eigvalsh(expected), [2,2,2,4], atol=1e-12, rtol=0):
                                        raise ValueError("hand-derived diamond metric spectrum failed")
                                qs = stages[s]["Q"][q]
                                scale = max(1, float(np.linalg.norm(expected)))
                                spectrum = np.linalg.eigvalsh(actual)
                                errors = {"orthogonality": float(np.linalg.norm(basis.T@basis-np.eye(size))),
                                          "source_space": float(np.linalg.norm(basis-qs@(qs.T@basis))),
                                          "symmetry": float(np.linalg.norm(actual-actual.T))/scale,
                                          "ambient_operator": float(np.linalg.norm(basis@actual@basis.T-qs@expected@qs.T))/scale,
                                          "spectrum": float(np.linalg.norm(spectrum-np.linalg.eigvalsh(expected)))/scale,
                                          "independent_domain": float(np.linalg.norm(expected-numerical))/scale}
                                if not np.all(np.isfinite(actual)) or any(not np.isfinite(v) or v > TOLERANCE for v in errors.values()):
                                    raise ValueError(f"numerical verification failed: {name} {regular} {(s,t,q)} {errors}")
                                if np.count_nonzero(np.abs(spectrum) <= TOLERANCE*scale) != rho or (size and min(spectrum) < -TOLERANCE*scale):
                                    raise ValueError("native zero spectrum disagrees with independent exact persistent rank")
                                request = {"case": name, "n": n, "edges": edges, "weights": weights, "q": q, "regular": regular, "pair": [s,t]}
                                records.write(json.dumps({"request": request, "input_sha256": hashlib.sha256(json.dumps(request, sort_keys=True).encode()).hexdigest(),
                                                          "certificate": result, "errors": errors, "exact_persistent_rank": rho,
                                                          "persistent_domain_dimension": domain, "reference_eigenvalues": np.linalg.eigvalsh(expected).tolist()})+"\n")
                                summary["certificates"] += 1
                                for key, value in errors.items():
                                    summary["max_errors"][key] = max(summary["max_errors"].get(key, 0), value)
        summary["rejected_inputs"] = 0
        base = [2, 1, 100, 10000, 1e-10, 0, 1, 0, 1, 1, 1, 0, 2]
        for index, value in [(2,0),(3,1),(4,0),(4,"NaN"),(5,2),(6,99),(9,-1),(9,"NaN")]:
            values = list(base)
            values[index] = value
            for regular in (True, False):
                rejected = subprocess.run([str(args.worker), *([] if regular else ["--non-regular"])],
                                          input=" ".join(map(str, values)), capture_output=True,
                                          text=True, encoding="utf-8", timeout=30)
                if not rejected.returncode or rejected.stdout.strip():
                    raise ValueError("invalid request returned a partial or successful certificate")
                summary["rejected_inputs"] += 1
        summary["sources_stable"] = before == {p.relative_to(ROOT).as_posix(): sha(p) for p in identities}
        if not summary["sources_stable"]:
            raise ValueError("measured source changed during validation")
        summary["status"] = "passed"
    except (ValueError, OSError, subprocess.SubprocessError) as error:
        summary["error"] = str(error)
    (args.output/"summary.json").write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k: v for k, v in summary.items() if k not in ("sources", "versions")}), flush=True)
    return 0 if summary["status"] == "passed" else 1


if __name__ == "__main__":
    raise SystemExit(main())
