#!/usr/bin/env python3
"""Independent NumPy/rational Hodge check; optional pinned source and GUDHI probes.

Run the Rust hodge_laplacian example first. NumPy is an optional oracle dependency,
not part of the Rust runtime or the repository's standard-library default checks.
"""

import argparse
import ast
from fractions import Fraction
import hashlib
from itertools import combinations
import json
from pathlib import Path
import platform
import sys
import textwrap

import numpy as np


def rational_rank(matrix):
    """Exact Gaussian elimination over Q, independent of Rust/spectral ranks."""
    rows = [[Fraction(int(x)) for x in row] for row in matrix]
    m, n = matrix.shape
    pivot = 0
    for column in range(n):
        row = next((i for i in range(pivot, m) if rows[i][column]), None)
        if row is None:
            continue
        rows[pivot], rows[row] = rows[row], rows[pivot]
        value = rows[pivot][column]
        rows[pivot] = [x / value for x in rows[pivot]]
        for i in range(pivot + 1, m):
            value = rows[i][column]
            rows[i] = [x - value * y for x, y in zip(rows[i], rows[pivot])]
        pivot += 1
        if pivot == m:
            break
    return pivot


def boundary(bases, degree):
    columns = bases.get(degree, [])
    rows = bases.get(degree - 1, []) if degree else []
    matrix = np.zeros((len(rows), len(columns)), dtype=np.int64)
    row_index = {simplex: i for i, simplex in enumerate(rows)}
    if degree:
        for j, simplex in enumerate(columns):
            for i in range(len(simplex)):
                face = simplex[:i] + simplex[i + 1:]
                matrix[row_index[face], j] = (-1) ** i
    return matrix


def pinned_toponetx(path):
    from scipy.sparse import csr_matrix, dok_matrix

    expected = 'ba5b33bc58c517fa018e0499a2adaf2ae63e02735986cc2c28e2352b394d2754'
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected
    # The fixed file uses Python 3.12 class parameter syntax. Read the untouched
    # method bodies individually so the existing Python 3.10 oracle can run them.
    lines = path.read_text(encoding='utf-8').splitlines(keepends=True)
    methods = []
    for name in ['incidence_matrix', 'hodge_laplacian_matrix']:
        start = next(i for i, line in enumerate(lines) if line.startswith(f'    def {name}('))
        end = next((i for i in range(start + 1, len(lines)) if lines[i].startswith(('    def ', '    @'))), len(lines))
        methods.extend(ast.parse(textwrap.dedent(''.join(lines[start:end]))).body)
    selected = ast.Module(body=methods, type_ignores=[])
    scope = {'np': np, 'csr_matrix': csr_matrix, 'dok_matrix': dok_matrix}
    exec(compile(ast.fix_missing_locations(selected), str(path), 'exec'), scope)
    cls = type('Pinned', (), {node.name: scope[node.name] for node in methods})
    cls.skeleton = lambda self, rank: self.bases.get(rank, [])
    return cls


def pinned_perslap(path):
    raw = path.read_bytes()
    blob = hashlib.sha1(b'blob ' + str(len(raw)).encode() + b'\0' + raw).hexdigest()
    assert blob == 'ef49f6caeac3c7710de84795f3ffb225d14dbc29', blob
    tree = ast.parse(path.read_text(encoding='utf-8'))
    methods = [node for node in tree.body if isinstance(node, ast.FunctionDef)
               and node.name in {'boundary', 'hodge_laplacian', 'spectrum_hodgelap'}]
    assert len(methods) == 3
    # Minimal exact row-membership adapter. Only closed simplicial fixtures run;
    # this is not validation of the external ismember package or whole project.
    def ismember(left, right, _mode):
        lookup = {tuple(row): i for i, row in enumerate(right)}
        return [tuple(row) in lookup for row in left], [lookup[tuple(row)] for row in left]
    scope = {'np': np, 'combinations': combinations, 'ismember': ismember}
    exec(compile(ast.fix_missing_locations(ast.Module(body=methods, type_ignores=[])), str(path), 'exec'), scope)
    return scope


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--rust', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--toponetx', type=Path)
    parser.add_argument('--perslap', type=Path)
    parser.add_argument('--gudhi', action='store_true')
    args = parser.parse_args()
    top = pinned_toponetx(args.toponetx) if args.toponetx else None
    perslap = pinned_perslap(args.perslap) if args.perslap else None
    if args.gudhi:
        import gudhi
    records = [json.loads(line) for line in args.rust.read_text(encoding='utf-8-sig').splitlines() if line]
    results = []
    external = {'toponetx_pass': 0, 'perslap_pass': 0, 'perslap_spectrum_pass': 0,
                'perslap_signed_boundary_pass': 0, 'gudhi_pass': 0, 'exclusions': []}
    paper_spectra = {}
    torsion = []
    max_harmonic_difference = 0.0
    max_residual = 0.0
    max_orthogonality = 0.0
    for record in records:
        simplices = [tuple(s) for s in record['simplices']]
        k = record['degree']
        bases = {}
        for simplex in simplices:
            bases.setdefault(len(simplex) - 1, []).append(simplex)
        n = len(bases.get(k, []))
        assert record['basis'] == [list(s) for s in bases.get(k, [])]
        a, d = boundary(bases, k), boundary(bases, k + 1)
        assert np.all(a @ d == 0)
        matrix = a.T @ a + d @ d.T
        actual = np.asarray(record['matrix']).reshape(n, n)
        assert np.array_equal(actual, matrix), record['name']
        assert np.array_equal(actual, actual.T)
        values, vectors = np.linalg.eigh(matrix.astype(np.float64))
        beta = n - rational_rank(a) - rational_rank(d)
        assert n - rational_rank(matrix) == beta
        np.testing.assert_allclose(record['values'], values, atol=1e-10, rtol=1e-10)
        rust_values = np.asarray(record['values'])
        rust_vectors = np.asarray(record['vectors']).reshape(n, n)
        assert np.count_nonzero(np.abs(rust_values) <= 1e-10) == beta
        kernel_columns = record['kernel_columns']
        assert kernel_columns == np.flatnonzero(np.abs(rust_values) <= 1e-10).tolist()
        # Lim section 6.2: the residual after projecting onto im A^T and im D
        # is the harmonic component. This oracle never uses eigenvectors to
        # construct that reference component.
        x = np.asarray(record['x'], dtype=float)
        lower_component = a.T @ np.linalg.lstsq(a.T.astype(float), x, rcond=None)[0]
        upper_component = d @ np.linalg.lstsq(d.astype(float), x, rcond=None)[0]
        harmonic_lstsq = x - lower_component - upper_component
        harmonic_vectors = rust_vectors[:, kernel_columns]
        harmonic_spectral = harmonic_vectors @ harmonic_vectors.T @ x
        np.testing.assert_allclose(harmonic_spectral, harmonic_lstsq, atol=1e-9, rtol=1e-9)
        if n:
            max_harmonic_difference = max(max_harmonic_difference, float(np.linalg.norm(harmonic_spectral - harmonic_lstsq)))
        if n:
            assert values[0] >= -1e-10
            measured = np.linalg.norm(actual @ rust_vectors - rust_vectors * rust_values, axis=0)
            np.testing.assert_allclose(record['residuals'], measured, atol=1e-12, rtol=1e-8)
            assert np.max(measured) < 1e-9
            orthogonality = np.linalg.norm(rust_vectors.T @ rust_vectors - np.eye(n))
            assert orthogonality < 1e-9
            assert abs(record['orthogonality_error'] - orthogonality) < 1e-12
            max_residual = max(max_residual, float(np.max(measured)))
            max_orthogonality = max(max_orthogonality, float(orthogonality))
            # Compare projectors for each repeated eigenspace, not signs/bases.
            groups = []
            for i, value in enumerate(values):
                if not groups or abs(value - values[groups[-1][0]]) > 1e-8:
                    groups.append([])
                groups[-1].append(i)
            for group in groups:
                u, v = vectors[:, group], rust_vectors[:, group]
                np.testing.assert_allclose(u @ u.T, v @ v.T, atol=1e-9, rtol=1e-9)
        np.testing.assert_allclose(record['action'], matrix @ record['x'], atol=1e-12, rtol=1e-12)
        dimension = max(bases, default=-1)
        label = f"{record['name']}:H{k}"
        if top:
            if k <= dimension and not (k == 0 and dimension == 0):
                reference = top()
                reference.bases, reference.dim = bases, dimension
                index, lap = reference.hodge_laplacian_matrix(k, signed=True, index=True)
                permutation = [index[s] for s in bases[k]]
                np.testing.assert_array_equal(lap.toarray()[np.ix_(permutation, permutation)], matrix)
                external['toponetx_pass'] += 1
            else:
                external['exclusions'].append({'case': label, 'oracle': 'TopoNetX', 'reason': 'absent degree or degree-zero source without edges rejected at fixed version'})
        if perslap:
            if 0 <= k < dimension and bases.get(k) and bases.get(k + 1):
                sc = [np.asarray(bases.get(i, [])).reshape(-1, i + 1) for i in range(dimension + 1)]
                lap, _, _ = perslap['hodge_laplacian'](sc, k)
                np.testing.assert_array_equal(lap, matrix)
                external['perslap_pass'] += 1
                # This fixed source reverses every positive-degree boundary's
                # overall sign relative to the increasing-vertex convention.
                np.testing.assert_array_equal(perslap['boundary'](sc, k + 1), -d)
                if k:
                    np.testing.assert_array_equal(perslap['boundary'](sc, k), -a)
                external['perslap_signed_boundary_pass'] += 1
                if record['name'] not in paper_spectra:
                    paper_spectra[record['name']] = perslap['spectrum_hodgelap'](sc)
                paper_matrices, spectra = paper_spectra[record['name']]
                np.testing.assert_array_equal(paper_matrices[k], matrix)
                assert np.max(np.abs(np.asarray(spectra[k]).imag)) <= 1e-10
                np.testing.assert_allclose(np.sort(np.asarray(spectra[k]).real), values, atol=1e-10, rtol=1e-10)
                external['perslap_spectrum_pass'] += 1
            else:
                external['exclusions'].append({'case': label, 'oracle': 'perslap', 'reason': 'highest/absent/empty degree outside supported fixed-version function'})
        if args.gudhi:
            source = gudhi.SimplexTree()
            for simplex in simplices:
                source.insert(simplex, filtration=0.)
            source.persistence(homology_coeff_field=101, persistence_dim_max=True)
            betti = source.betti_numbers() if simplices else []
            assert (betti[k] if k < len(betti) else 0) == beta
            external['gudhi_pass'] += 1
            if record['name'] == 'rp2':
                source.persistence(homology_coeff_field=2, persistence_dim_max=True)
                betti_f2 = source.betti_numbers()
                f2 = betti_f2[k] if k < len(betti_f2) else 0
                expected = [1,1,1,0][k]
                assert f2 == expected
                assert beta == [1,0,0,0][k]
                torsion.append({'degree': k, 'betti_rational': beta, 'betti_f2_gudhi': f2})
        results.append({'name': record['name'], 'degree': k, 'n': n, 'betti_rational': beta, 'values': values.tolist()})
    # SPD path handoff oracle only: this is deliberately not a Rust weighted API.
    a = np.array([[-1, 0], [1, -1], [0, 1]], dtype=float)
    gram = np.diag([2., 3.])
    weighted = np.linalg.solve(gram, a.T @ np.diag([1., 2., 3.]) @ a)
    np.testing.assert_allclose(weighted, [[1.5, -1], [-2/3, 5/3]])
    np.testing.assert_allclose(gram @ weighted, weighted.T @ gram)
    root = np.diag(np.sqrt([2., 3.]))
    weighted_values = np.linalg.eigvalsh(root @ weighted @ np.linalg.inv(root))
    np.testing.assert_allclose(weighted_values, [(19-np.sqrt(97))/12, (19+np.sqrt(97))/12])
    report = {'python': sys.version, 'platform': platform.platform(), 'numpy': np.__version__,
              'gudhi': gudhi.__version__ if args.gudhi else None,
              'reference_versions': {'toponetx': 'ce3d414b571fff95916c4cbc20b8d2c2c9d428e6', 'perslap': 'b36b080f739653d84284f3f7ca9ccea5fad286ae'},
              'windows': len(results), 'cases': len({r['name'] for r in results}),
              'max_residual': max_residual, 'max_orthogonality': max_orthogonality,
              'lim_lstsq_harmonic_checks': len(results), 'max_harmonic_difference': max_harmonic_difference,
              'torsion_rp2': torsion,
              'external': external, 'weighted_path_oracle': {'matrix': weighted.tolist(), 'values': weighted_values.tolist(), 'rust_support': False},
              'evidence': 'numerical eigenpairs plus exact rational ranks on this finite corpus; GUDHI uses field 101, not real coefficients',
              'results': results,
              'reference_sha256': {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in [('toponetx', args.toponetx), ('perslap', args.perslap)] if path}}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + '\n', encoding='utf-8')
    print(json.dumps({k: report[k] for k in ['cases', 'windows', 'max_residual', 'max_orthogonality']}))
    print(json.dumps({k: v for k, v in external.items() if k != 'exclusions'}))


if __name__ == '__main__':
    main()
