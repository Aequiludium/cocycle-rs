"""Check actual Rust Dirac fixtures against independent matrices and spectra.

NumPy is required only when running the numerical checks. Optional fixed
TopoNetX methods need SciPy and Python 3.12; --gudhi requires its Python binding.
The core crate and ordinary tool unit tests acquire no new dependencies.
"""
import argparse
import ast
import hashlib
import itertools
import json
import math
from pathlib import Path
import platform
import sys
import types

TOPONETX_COMMIT = "ce3d414b571fff95916c4cbc20b8d2c2c9d428e6"
TOPONETX_BLOB = "8ace5fef7ad5e9a44cce02b284578444340bbae9"
TOLERANCE = 1e-11
CASES = {
    "empty": ([], 0, 2),
    "isolates": ([(5,), (42,)], 0, 2),
    "edge": ([(5, 42)], 0, 1),
    "edge_only_c1": ([(5, 42)], 1, 1),
    "edge_missing_degrees": ([(5, 42)], 3, 5),
    "path": ([(2, 9), (9, 30)], 0, 1),
    "triangle_ring": ([(0, 1), (0, 2), (1, 2)], 0, 2),
    "filled_triangle": ([(0, 1, 2)], 0, 2),
    "filled_triangle_c1_c2": ([(0, 1, 2)], 1, 2),
    "tetrahedron": ([(0, 1, 2, 3)], 0, 3),
    "tetrahedron_boundary": ([(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)], 0, 2),
    "path_96": ([(i, i + 1) for i in range(95)], 0, 1),
    "cycle_64": ([(min(i, (i + 1) % 64), max(i, (i + 1) % 64)) for i in range(64)], 0, 1),
    "disconnected_edges_16": ([(2*i, 2*i+1) for i in range(16)], 0, 1),
}
ROOT2, ROOT3 = math.sqrt(2), math.sqrt(3)
ANALYTIC = {
    "empty": [], "isolates": [0, 0],
    "edge": [-ROOT2, 0, ROOT2], "edge_only_c1": [0],
    "edge_missing_degrees": [], "path": [-ROOT3, -1, 0, 1, ROOT3],
    "triangle_ring": [-ROOT3] * 2 + [0] * 2 + [ROOT3] * 2,
    "filled_triangle": [-ROOT3] * 3 + [0] + [ROOT3] * 3,
    "filled_triangle_c1_c2": [-ROOT3, 0, 0, ROOT3],
    "tetrahedron": [-2] * 7 + [0] + [2] * 7,
    "tetrahedron_boundary": [-2] * 6 + [0] * 2 + [2] * 6,
    "path_96": sorted([-2*math.sin(math.pi*k/192) for k in range(1,96)]
                      + [0] + [2*math.sin(math.pi*k/192) for k in range(1,96)]),
    "cycle_64": sorted([-2*math.sin(math.pi*k/64) for k in range(1,64)]
                       + [0, 0] + [2*math.sin(math.pi*k/64) for k in range(1,64)]),
    "disconnected_edges_16": [-ROOT2]*16 + [0]*16 + [ROOT2]*16,
}


def read_rows(path):
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8-sig").splitlines()
            if line.strip()]
    required = {"name", "basis", "offsets", "truncated", "matrix", "x",
                "action", "square_action"}
    for row in rows:
        if not isinstance(row, dict) or not required <= row.keys():
            raise ValueError("incomplete Rust Dirac record")
        if not isinstance(row["name"], str) or not isinstance(row["truncated"], bool):
            raise ValueError("invalid Rust fixture name or truncation flag")
    if len(rows) != len(CASES) or {row["name"] for row in rows} != set(CASES):
        raise ValueError("Rust output must contain each of the fourteen fixtures exactly once")
    return rows


def checked_toponetx_source(path):
    source = path.read_bytes()
    blob = hashlib.sha1(b"blob " + str(len(source)).encode() + b"\0" + source).hexdigest()
    if blob != TOPONETX_BLOB:
        raise ValueError(f"TopoNetX source blob mismatch: {blob}")
    return source


def toponetx_methods(path, np):
    source = checked_toponetx_source(path)
    from scipy.sparse import csr_matrix, dok_matrix

    tree = ast.parse(source.decode("utf-8"))
    cls = next(node for node in tree.body
               if isinstance(node, ast.ClassDef) and node.name == "SimplicialComplex")
    names = {"incidence_matrix", "dirac_operator_matrix", "hodge_laplacian_matrix"}
    methods = [node for node in cls.body
               if isinstance(node, ast.FunctionDef) and node.name in names]
    if {node.name for node in methods} != names:
        raise ValueError("fixed TopoNetX source does not contain the required methods")
    namespace = {"np": np, "csr_matrix": csr_matrix, "dok_matrix": dok_matrix}
    exec(compile(ast.Module(body=methods, type_ignores=[]), str(path), "exec"), namespace)

    class SimplexTuple(tuple):
        @property
        def elements(self):
            return tuple(self)

    class FixedFunctions:
        def __init__(self, bases):
            self.bases = bases
            self.dim = len(bases) - 1
            self._simplex_trie = types.SimpleNamespace(shape=[len(x) for x in bases])

        def skeleton(self, rank):
            return [SimplexTuple(x) for x in self.bases[rank]]

    for name in names:
        setattr(FixedFunctions, name, namespace[name])
    return FixedFunctions


def faces(facets):
    return sorted({face for cell in facets for size in range(1, len(cell) + 1)
                   for face in itertools.combinations(cell, size)})


def check(rust_path, toponetx_path=None, use_gudhi=False, check_pairs=False):
    rows = read_rows(rust_path)
    import numpy as np

    fixed_functions = toponetx_methods(toponetx_path, np) if toponetx_path else None
    gudhi = None
    if use_gudhi:
        import gudhi
    scipy_version = None
    if fixed_functions:
        import scipy
        scipy_version = scipy.__version__
    results = []
    for row in rows:
        name = row["name"]
        facets, a, b = CASES[name]
        full = faces(facets)
        expected = [face for degree in range(a, b + 1)
                    for face in full if len(face) - 1 == degree]
        basis = list(map(tuple, row["basis"]))
        assert set(basis) == set(expected) and len(basis) == len(expected), name
        assert [len(face) for face in basis] == [len(face) for face in expected], name
        offsets = [0]
        for degree in range(a, b + 1):
            offsets.append(offsets[-1] + sum(len(face) - 1 == degree for face in full))
        assert row["offsets"] == offsets, name
        assert row["truncated"] == any(not a <= len(face) - 1 <= b for face in full), name
        n = len(basis)
        positions = {face: i for i, face in enumerate(basis)}
        d = np.zeros((n, n))
        for cell, column in positions.items():
            if len(cell) - 1 == a:
                continue  # Ordinary/truncated bottom boundary is zero.
            for omitted in range(len(cell)):
                face = cell[:omitted] + cell[omitted + 1:]
                d[positions[face], column] = (-1) ** omitted
        np.testing.assert_array_equal(d @ d, np.zeros((n, n)))
        expected_k = d + d.T
        rust = np.asarray(row["matrix"], dtype=float).reshape(n, n)
        x = np.asarray(row["x"], dtype=float)
        action = np.asarray(row["action"], dtype=float)
        squared_action = np.asarray(row["square_action"], dtype=float)
        for vector in (x, action, squared_action):
            assert vector.shape == (n,) and np.isfinite(vector).all(), name
        assert np.isfinite(rust).all(), name
        np.testing.assert_array_equal(rust, expected_k)
        np.testing.assert_array_equal(rust, rust.T)
        np.testing.assert_array_equal(action, expected_k @ x)
        np.testing.assert_array_equal(squared_action, (expected_k @ expected_k) @ x)
        if check_pairs:
            spectra = row['boundary_spectra']
            expected_parts = {}
            for slot, degree in enumerate(range(a+1, b+1), start=1):
                block = d[offsets[slot-1]:offsets[slot], offsets[slot]:offsets[slot+1]]
                singular_values = np.linalg.svd(block, compute_uv=False)
                nonzero = singular_values[singular_values > 1e-8]
                if len(nonzero):
                    expected_parts[degree] = nonzero
            assert {p['degree'] for p in spectra} == set(expected_parts), name
            assert len(spectra) == len(expected_parts), name
            all_signed = []
            gamma = np.asarray([(-1.)**(len(face)-1) for face in basis])
            for part in spectra:
                complete_reference = expected_parts[part['degree']]
                count = part['count']
                assert 0 < count <= len(complete_reference)
                reference = complete_reference[:count]
                np.testing.assert_allclose(part['singular_values'], reference, atol=1e-8, rtol=1e-8)
                wanted = np.concatenate([-reference, reference[::-1]])
                np.testing.assert_allclose(part['values'], wanted, atol=1e-8, rtol=1e-8)
                q = np.asarray(part['vectors']).reshape(n,2*count)
                v = np.asarray(part['values'])
                errors = np.linalg.norm(expected_k @ q-q*v,axis=0)
                np.testing.assert_allclose(part['residuals'],errors,atol=1e-12,rtol=1e-6)
                assert np.all(errors <= 1.01e-9*(1/np.abs(v)+np.abs(v)))
                assert np.linalg.norm(q.T @ q-np.eye(2*count)) < 1e-8
                np.testing.assert_allclose(q[:,:count],gamma[:,None]*q[:,::-1][:,:count],atol=0,rtol=0)
                all_signed.extend(part['values'])
            complete = np.linalg.eigvalsh(expected_k)
            if name not in {'path_96', 'cycle_64', 'disconnected_edges_16'}:
                np.testing.assert_allclose(sorted(all_signed),complete[np.abs(complete)>1e-8],atol=1e-8,rtol=1e-8)
        hodge = np.zeros((n, n))
        singular = []
        for degree in range(a, b + 1):
            ids = [i for i, face in enumerate(basis) if len(face) - 1 == degree]
            prev = [i for i, face in enumerate(basis) if len(face) - 1 == degree - 1]
            nxt = [i for i, face in enumerate(basis) if len(face) - 1 == degree + 1]
            lower, upper = d[np.ix_(prev, ids)], d[np.ix_(ids, nxt)]
            hodge[np.ix_(ids, ids)] = lower.T @ lower + upper @ upper.T
            singular.extend(s for s in np.linalg.svd(lower, compute_uv=False) if s > 1e-10)
        np.testing.assert_array_equal(rust @ rust, hodge)
        parity = np.diag([(-1) ** (len(face) - 1) for face in basis])
        np.testing.assert_array_equal(parity @ rust, -rust @ parity)
        values, vectors = np.linalg.eigh(rust)
        svd_values = sorted([-s for s in singular] + singular + [0] * (n - 2 * len(singular)))
        np.testing.assert_allclose(values, svd_values, rtol=0, atol=TOLERANCE)
        np.testing.assert_allclose(values, ANALYTIC[name], rtol=0, atol=TOLERANCE)
        np.testing.assert_allclose(values, -values[::-1], rtol=0, atol=TOLERANCE)
        residual = float(np.linalg.norm(rust @ vectors - vectors * values))
        orthogonality = float(np.linalg.norm(vectors.T @ vectors - np.eye(n)))
        assert residual < TOLERANCE and orthogonality < TOLERANCE, name
        tnx_agreement = None
        if fixed_functions and full and max(map(len, full)) > 1:
            bases = [[face for face in full if len(face) == degree + 1]
                     for degree in range(max(map(len, full)))]
            fixed = fixed_functions(bases)
            index, tnx = fixed.dirac_operator_matrix(signed=True, index=True)
            permutation = [index[face] for face in basis]
            np.testing.assert_array_equal(rust, tnx.toarray()[np.ix_(permutation, permutation)])
            if not row["truncated"]:
                for degree, base in enumerate(bases):
                    index_l, laplacian = fixed.hodge_laplacian_matrix(degree, index=True)
                    ids = [positions[face] for face in base]
                    reorder = [index_l[face] for face in base]
                    np.testing.assert_array_equal(hodge[np.ix_(ids, ids)],
                                                  laplacian.toarray()[np.ix_(reorder, reorder)])
            tnx_agreement = True
        gudhi_betti = None
        if gudhi and full and not row["truncated"]:
            tree = gudhi.SimplexTree()
            for cell in facets:
                tree.insert(cell, filtration=0.)
            tree.compute_persistence(homology_coeff_field=3, persistence_dim_max=True)
            gudhi_betti = tree.betti_numbers()
            assert sum(gudhi_betti) == sum(value == 0 for value in ANALYTIC[name]), name
        results.append({"name": name, "basis": basis, "degree_offsets": offsets,
                        "coordinates": n, "nonzeros": int(np.count_nonzero(rust)),
                        "spectrum": values.tolist(), "residual": residual,
                        "orthogonality": orthogonality, "toponetx_functions": tnx_agreement,
                        "gudhi_full_f3_betti": gudhi_betti})
    return {
        "versions": {"python": sys.version, "numpy": np.__version__, "scipy": scipy_version,
                     "gudhi": gudhi.__version__ if gudhi else None, "platform": platform.platform()},
        "rust_output_sha256": hashlib.sha256(rust_path.read_bytes()).hexdigest(),
        "toponetx_commit": TOPONETX_COMMIT if fixed_functions else None,
        "source_blob": TOPONETX_BLOB if fixed_functions else None,
        "toponetx_scope": "Three unchanged AST-extracted methods with supplied bases; not full package integration. Empty/isolated shapes excluded." if fixed_functions else "not requested",
        "gudhi_scope": "Complete ordinary torsion-free fixtures over F3; kernel dimension only. Empty/truncated cases excluded." if gudhi else "not requested",
        "absolute_spectral_tolerance": TOLERANCE,
        "rust_signed_pair_cases": len(rows) if check_pairs else 0,
        "rust_boundary_singular_modes": sum(len(p['singular_values']) for row in rows for p in row['boundary_spectra']) if check_pairs else 0,
        "cases": results,
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rust", type=Path, required=True, help="Rust fixture JSONL")
    parser.add_argument("--output", type=Path, required=True, help="result JSON under target/")
    parser.add_argument("--toponetx", type=Path, help="fixed simplicial_complex.py at ce3d414")
    parser.add_argument("--gudhi", action="store_true", help="also check full-fixture F3 Betti numbers")
    parser.add_argument("--pairs", action="store_true", help="require Rust signed boundary eigenpairs and compare NumPy SVD/chirality")
    args = parser.parse_args(argv)
    try:
        report = check(args.rust, args.toponetx, args.gudhi, args.pairs)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    except (OSError, ValueError, AssertionError, ImportError, SyntaxError, TypeError, KeyError) as error:
        print(f"Dirac check failed: {error}", file=sys.stderr)
        return 1
    cases = report["cases"]
    print(json.dumps({"passed": len(cases),
                      "toponetx_cases": sum(row["toponetx_functions"] is True for row in cases),
                      "gudhi_cases": sum(row["gudhi_full_f3_betti"] is not None for row in cases),
                      "max_residual": max(row["residual"] for row in cases),
                      "max_orthogonality": max(row["orthogonality"] for row in cases)}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
