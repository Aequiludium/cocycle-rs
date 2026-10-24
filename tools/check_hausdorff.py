"""Check actual Rust Hausdorff JSONL against exact rational/high-precision oracles.

The --references option requires NumPy, SciPy and GUDHI; all other checks use
only the standard library. GUDHI's exact cKDTree wrapper reuses SciPy, so the
independent oracle is rational arithmetic, not a second external numerical core.
"""
import argparse
from decimal import Decimal, localcontext
from fractions import Fraction
import hashlib
import json
import math
from pathlib import Path
import platform


HAND = {"single_345": 25, "asymmetric": 4, "duplicates": 0,
        "translated": 1, "outlier": 9801}
STRESS = {"stress_large", "stress_tiny", "stress_subnormal",
          "stress_irrelevant_overflow", "stress_final_overflow"}
NAMES = set(HAND) | STRESS | {f"random_{d}_{i}" for d in (1, 2, 3, 8) for i in range(20)}


def rows_from(path):
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8-sig").splitlines()
            if line.strip()]
    matrices = [r for r in rows if "matrix_ids" in r]
    pairs = [r for r in rows if "name" in r]
    if len(rows) != 91 or len(pairs) != 90 or {r["name"] for r in pairs} != NAMES:
        raise ValueError("expected exactly 90 distinct pair fixtures and one sample matrix")
    if matrices != [{"matrix_ids": ["c", "a", "b", "same"],
                     "matrix_values": [3.0, 2.0, 1.0, 3.0, 0.0, 1.0]}]:
        raise ValueError("sample identity/order or condensed matrix mismatch")
    for row in pairs:
        d = row["dimension"]
        if type(d) is not int or not 1 <= d <= 8:
            raise ValueError("invalid coordinate dimension")
        for key in ("first", "second"):
            values = row[key]
            if not isinstance(values, list) or not values or len(values) % d or len(values) > 56:
                raise ValueError("invalid point-cloud shape")
            if any(type(v) not in (int, float) or not math.isfinite(v) for v in values):
                raise ValueError("nonfinite or nonnumeric coordinate")
        value = row["distance"]
        if value != "numerical_failure" and (type(value) not in (int, float)
                                             or not math.isfinite(value) or value < 0):
            raise ValueError("invalid reported distance")
    return pairs


def exact_squared(row):
    d = row["dimension"]
    a, b = [[tuple(Fraction(v) for v in r[key][i:i+d])
             for i in range(0, len(r[key]), d)] for r, key in ((row, "first"), (row, "second"))]
    # Squaring is exact over rational representations of the original binary64
    # coordinates. No Rust norm, f64 square or nearest-neighbor engine is reused.
    distances = [[sum((x-y)**2 for x, y in zip(p, q)) for q in b] for p in a]
    return max(max(min(r) for r in distances),
               max(min(distances[i][j] for i in range(len(a))) for j in range(len(b))))


def expected_distance(squared):
    with localcontext() as context:
        context.prec = 120
        return float((Decimal(squared.numerator) / Decimal(squared.denominator)).sqrt())


def close(actual, expected):
    # No absolute tolerance: zero must never pass for a tiny nonzero distance.
    if not math.isclose(actual, expected, rel_tol=5e-14, abs_tol=0.0):
        raise ValueError(f"distance mismatch: {actual!r} != {expected!r}")


def check(path, references=False):
    pairs = rows_from(path)
    external = 0
    versions = {"python": platform.python_version()}
    if references:
        import numpy as np
        import scipy
        import gudhi
        from scipy.spatial.distance import directed_hausdorff
        from gudhi.point_cloud.knn import KNearestNeighbors
        import gudhi.point_cloud.knn as knn
        versions.update(numpy=np.__version__, scipy=scipy.__version__, gudhi=gudhi.__version__,
                        gudhi_knn_sha256=hashlib.sha256(Path(knn.__file__).read_bytes()).hexdigest())
    for row in pairs:
        squared = exact_squared(row)
        if row["name"] in HAND and squared != HAND[row["name"]]:
            raise ValueError(f"hand fixture changed: {row['name']}")
        expected = expected_distance(squared)
        if math.isinf(expected):
            if row["distance"] != "numerical_failure":
                raise ValueError("unrepresentable result must be a numerical error")
        else:
            if row["distance"] == "numerical_failure":
                raise ValueError("representable distance must not be an error")
            close(row["distance"], expected)
        if references and row["name"] not in STRESS:
            d = row["dimension"]
            a = np.asarray(row["first"], dtype=np.float64).reshape(-1, d)
            b = np.asarray(row["second"], dtype=np.float64).reshape(-1, d)
            close(max(directed_hausdorff(a, b, 0)[0], directed_hausdorff(b, a, 0)[0]), expected)
            opts = dict(k=1, return_index=False, return_distance=True,
                        implementation="ckdtree", eps=0, n_jobs=1)
            close(max(float(KNearestNeighbors(**opts).fit(a).transform(b).max()),
                      float(KNearestNeighbors(**opts).fit(b).transform(a).max())), expected)
            external += 1
    return {"pairs": len(pairs), "sample_matrices": 1, "external_pairs": external,
            "stress_reference_exclusions": sorted(STRESS) if references else [],
            "relative_tolerance": 5e-14, "absolute_tolerance": 0,
            "fixture_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "versions": versions}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSONL from the Rust hausdorff --fixtures example")
    parser.add_argument("--references", action="store_true", help="require SciPy and GUDHI on ordinary-range cases")
    args = parser.parse_args()
    print(json.dumps(check(args.input, args.references), indent=2))


if __name__ == "__main__":
    main()
