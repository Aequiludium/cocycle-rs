"""Replay actual comparison-route output with exact expectations and optional references.

GUDHI recomputes the two Rips and two supplied-filtration diagrams. A native
Topp worker, when requested, must match the pinned build's recorded binary hash.
No production comparison or PH kernel is used to generate expected results.
"""
import argparse
from fractions import Fraction
import json
import math
from pathlib import Path
import tempfile

from check_hausdorff import close, exact_squared, expected_distance
from distance_common import METRICS, PINS, invoke, sha256, tiny_oracle, unpack_number, write_fixture


def require(condition, message):
    if not condition:
        raise ValueError(message)


def decode_diagram(values):
    return [(birth, math.inf if death is None else death) for birth, death in values]


def cell_map(cells):
    result = {tuple(vertices): Fraction(value) for vertices, value in cells}
    require(len(result) == len(cells), "duplicate cell keys")
    return result


def check(rows, references=False, topp=None, topp_build=None):
    require(len(rows) == 3 and {r['route'] for r in rows} == {'geometry', 'results', 'filtrations'},
            "expected each of three routes exactly once")
    records = {r['route']: r for r in rows}
    declarations = {
        'geometry': (['a', 'b', 'c'], ['measurements/a/v1', 'measurements/b/v1', 'measurements/c/v1'], 'metres', 'lab-cartesian-v1'),
        'results': (['rips-a', 'rips-b'], ['rips/a/v1', 'rips/b/v1'], 'metres', 'lab-cartesian-v1'),
        'filtrations': (['filtration-a', 'filtration-b'], ['triangle/a/v1', 'triangle/b/v1'], 'seconds', 'triangle-labels-v1'),
    }
    for route, row in records.items():
        ids, snapshots, units, alignment = declarations[route]
        require((row['ids'], row['snapshots'], row['units'], row['alignment'], row['normalization'])
                == (ids, snapshots, units, alignment, 'none'), "operand identity/units/alignment mismatch")
        for key, value in row.items():
            if key.endswith('_ns'):
                require(type(value) is int and value >= 0, "invalid resource timing")
    geometry = records['geometry']
    require(geometry['dimension'] == 2 and geometry['points'] == [[0, 0, 2, 0], [0, 1, 2, 1], [0, 3, 2, 3]],
            "geometry fixture changed")
    expected = [expected_distance(exact_squared({'dimension': 2, 'first': geometry['points'][i],
                                                'second': geometry['points'][j]}))
                for i in range(3) for j in range(i)]
    require(geometry['values'] == expected and geometry['numeric_bytes'] == 24,
            "identified matrix mismatch")
    require(geometry['nearest'] == geometry['ids'][1] and geometry['converted_alias'] is True,
            "nearest-sample or conversion mismatch")

    result = records['results']
    require(result['first_points'] == [0, 2] and result['second_points'] == [0, 3], "Rips input changed")
    require(result['first_diagram'] == [[0, 2], [0, None]]
            and result['second_diagram'] == [[0, 3], [0, None]], "hand Rips diagram mismatch")
    require(result['converted'] == result['distances'] and result['within_one_metre'] is True,
            "conversion or scalar-consumer mismatch")

    filtration = records['filtrations']
    first, second = cell_map(filtration['first_cells']), cell_map(filtration['second_cells'])
    keys = {(10,), (20,), (30,), (10, 20), (10, 30), (20, 30), (10, 20, 30)}
    require(set(first) == set(second) == keys, "filtered domains differ")
    delta = max(abs(first[key] - second[key]) for key in keys)
    require(delta == 2 and filtration['sup'] == delta, "key-aligned sup distance mismatch")
    require(filtration['first_diagram'] == [[3, 6]] and filtration['second_diagram'] == [[3, 7]],
            "hand filtered diagram mismatch")
    require(filtration['bound_holds'] is True and filtration['distances'][0] <= delta,
            "fixed-domain stability bound mismatch")

    gudhi_count, topp_count = 0, 0
    versions = {}
    if references:
        import gudhi
        versions['gudhi'] = gudhi.__version__
    if (topp is None) != (topp_build is None):
        raise ValueError("Topp executable and build metadata must be supplied together")
    if topp:
        build = json.loads(topp_build.read_text(encoding='utf-8-sig'))
        require(build['topp_revision'] == PINS['topp']['revision'], "Topp revision mismatch")
        require(build['binaries_sha256']['topp'] == sha256(topp), "Topp binary hash mismatch")
        versions.update(topp_revision=build['topp_revision'], topp_binary_sha256=sha256(topp),
                        topp_build_sha256=sha256(topp_build))
    evidence = []
    for route in ('results', 'filtrations'):
        row = records[route]
        require(row['field'] == 2 and row['degree'] == (0 if route == 'results' else 1), "field or degree mismatch")
        a, b = decode_diagram(row['first_diagram']), decode_diagram(row['second_diagram'])
        require(len(row['distances']) == 3, "missing matching metric")
        for value, metric in zip(row['distances'], METRICS):
            close(value, tiny_oracle(a, b, metric))
        if references:
            for side in ('first', 'second'):
                if route == 'results':
                    points = [[v] for v in row[side + '_points']]
                    tree = gudhi.RipsComplex(points=points).create_simplex_tree(max_dimension=1)
                else:
                    tree = gudhi.SimplexTree()
                    for vertices, value in row[side + '_cells']:
                        tree.insert(vertices, filtration=value)
                tree.compute_persistence(homology_coeff_field=2, min_persistence=0, persistence_dim_max=True)
                actual = sorted(map(tuple, tree.persistence_intervals_in_dimension(row['degree']).tolist()))
                require(actual == sorted(decode_diagram(row[side + '_diagram'])), "GUDHI recomputation mismatch")
                gudhi_count += 1
        if topp:
            with tempfile.TemporaryDirectory() as directory:
                fixture = Path(directory) / 'comparison.bin'
                write_fixture(fixture, a, b)
                for value, metric in zip(row['distances'], METRICS):
                    outcome = invoke([str(topp)], fixture, metric, 'baseline')
                    require(outcome['status'] == 'completed', f"Topp worker failed: {outcome}")
                    close(unpack_number(outcome), value)
                    evidence.append({'route': route, 'metric': metric, 'worker': outcome,
                                     'fixture_sha256': sha256(fixture)})
                    topp_count += 1
    return {'routes': 3, 'exact_matrix_pairs': 3, 'exact_matching_checks': 6,
            'gudhi_diagrams': gudhi_count, 'topp_checks': topp_count,
            'versions': versions, 'topp_evidence': evidence,
            'scope': 'Declared aligned fixtures only; no generic registration, spectral comparison or global performance claim.'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--references', action='store_true', help='require GUDHI PH recomputation')
    parser.add_argument('--topp', type=Path, help='pinned native Topp worker')
    parser.add_argument('--topp-build', type=Path, help='verified build metadata for that worker')
    args = parser.parse_args()
    rows = [json.loads(line) for line in args.input.read_text(encoding='utf-8-sig').splitlines() if line.strip()]
    report = check(rows, args.references, args.topp, args.topp_build)
    report['rust_output_sha256'] = sha256(args.input)
    print(json.dumps(report, indent=2, allow_nan=False))


if __name__ == '__main__':
    main()
