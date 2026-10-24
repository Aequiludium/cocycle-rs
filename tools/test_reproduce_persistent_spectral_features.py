"""Independent and boundary checks for the optional #111 Python reproduction."""

import copy
import importlib.util
import io
import json
import os
from pathlib import Path
import random
import subprocess
import tempfile
import unittest
from unittest.mock import patch
import warnings
import zipfile

import reproduce_persistent_spectral_features as spectral


HAS_NUMPY = importlib.util.find_spec("numpy") is not None
HAS_SCIPY = importlib.util.find_spec("scipy") is not None
HAS_SKLEARN = importlib.util.find_spec("sklearn") is not None


class SourceSnapshotTests(unittest.TestCase):
    def test_offline_snapshot_is_verified_and_corruption_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            source = root / "source.py"
            source.write_text("pass\n", encoding="utf-8")
            spectral.write_json(root / "SOURCE_SNAPSHOT.json", {"base_sha": "fixed-base",
                                "sources": {"source.py": {"sha256": spectral.file_hash(source)}}})
            with patch.object(spectral, "ROOT", root):
                self.assertEqual(spectral.source_provenance(), ("fixed-base", "verified source snapshot; no Git checkout"))
                source.write_text("changed\n", encoding="utf-8")
                with self.assertRaises(ValueError):
                    spectral.source_provenance()

    def test_snapshot_cannot_read_outside_its_root(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            spectral.write_json(root / "SOURCE_SNAPSHOT.json", {"base_sha": "fixed-base",
                                "sources": {"../outside.py": {"sha256": "unused"}}})
            with patch.object(spectral, "ROOT", root), self.assertRaises(ValueError):
                spectral.source_provenance()


def sample(grid=(0.,), values=(0., 2.), n=2):
    records = [{"degree": 0, "i": i, "j": j, "sample_id": "s", "source_id": "source",
                "status": "computed", "coverage": "full", "domain_size": n,
                "basis": [(v,) for v in range(n)], "values": list(values)}
               for i in range(len(grid)) for j in range(i, len(grid))]
    return {"sample_id": "s", "source_id": "source", "grid": list(grid),
            "operator": "persistent_laplacian", "field": "R", "inner_product": "identity",
            "scale_unit": "u", "eigenvalue_unit": "v",
            "numeric_origin": {"producer": "fixture", "solver": "exact-small-example",
                               "dtype": "float64", "zero_tolerance": 1e-10}, "spectra": records}


@unittest.skipUnless(HAS_NUMPY, "optional NumPy oracle environment required")
class VectorizationTests(unittest.TestCase):
    def vectorize(self, value, length=3, grid=(0.,)):
        return spectral.vectorize_spectra([value], grid, degrees=(0,), length=length)

    def test_four_point_uniform_sampling_and_160_columns(self):
        np = spectral.numpy()
        np.testing.assert_array_equal(spectral.sampling_grid(-3, 3, 4), [-3, -1, 1, 3])
        result = spectral.vectorize_spectra([], (0., 1., 2., 3.))
        self.assertEqual(result["values"].shape, (0, 160))
        self.assertEqual(result["validity_mask"].shape, (0, 160))

    def test_invalid_grids_and_counts(self):
        for grid in ([], [0, 0], [1, 0], [0, float("nan")], [float("inf")]):
            with self.subTest(grid=grid), self.assertRaises(ValueError):
                spectral.vectorize_spectra([], grid)
        for count in (0, -1, True, 1.5, 1):
            with self.subTest(count=count), self.assertRaises(ValueError):
                spectral.sampling_grid(0, 1, count)
        self.assertEqual(spectral.sampling_grid(3, 3, 1), (3.,))
        with self.assertRaises(ValueError):
            spectral.sampling_grid(3, 3, 4)

    def test_zero_and_padding_have_distinct_masks(self):
        result = self.vectorize(sample(values=(0.,), n=1))
        self.assertEqual(result["values"].tolist(), [[0., 0., 0.]])
        self.assertEqual(result["validity_mask"].tolist(), [[True, False, False]])

    def test_empty_domain_is_valid_but_absent_spectrum_is_not(self):
        result = self.vectorize(sample(values=(), n=0))
        self.assertEqual(result["validity_mask"].tolist(), [[False, False, False]])
        value = sample()
        value["spectra"] = []
        with self.assertRaises(ValueError):
            self.vectorize(value)

    def test_all_zero_spectrum_is_retained(self):
        result = self.vectorize(sample(values=(0., 0.), n=2))
        self.assertEqual(result["validity_mask"].tolist(), [[True, True, False]])

    def test_smallest_truncation_and_repeated_eigenvalues(self):
        result = self.vectorize(sample(values=(7., 2., 0., 2.), n=4))
        self.assertEqual(result["values"].tolist(), [[0., 2., 2.]])
        self.assertEqual(result["domain_sizes"], [[4]])

    def test_global_pair_offsets_and_degree_order(self):
        value = sample(grid=(0., 1.), values=(0.,), n=1)
        for record, number in zip(value["spectra"], [1., 2., 3.]):
            record["values"] = [number]
        extra = [{**r, "degree": 1, "values": [r["values"][0] + 10]} for r in value["spectra"]]
        value["spectra"].extend(extra)
        result = spectral.vectorize_spectra([value], (0., 1.), length=2)
        self.assertEqual(result["values"].tolist(), [[1., 0., 2., 0., 3., 0., 11., 0., 12., 0., 13., 0.]])
        self.assertEqual([(c["degree"], c["i"], c["j"], c["rank"]) for c in result["columns"]][6], (1, 0, 0, 0))

    def test_complete_and_smallest_partial_selection(self):
        value = sample(values=(0., 1., 2.), n=10)
        value["spectra"][0]["coverage"] = "smallest"
        self.assertEqual(self.vectorize(value)["values"].tolist(), [[0., 1., 2.]])
        value["spectra"][0]["values"] = [0., 1.]
        with self.assertRaises(ValueError):
            self.vectorize(value)
        value["spectra"][0]["coverage"] = "largest"
        with self.assertRaises(ValueError):
            self.vectorize(value)

    def test_full_spectrum_cannot_hide_missing_values(self):
        with self.assertRaises(ValueError):
            self.vectorize(sample(values=(0.,), n=2))

    def test_invalid_values_are_rejected(self):
        for values in ([float("nan")], [float("inf")], [-1.], [1 + 0j], [[1.]]):
            with self.subTest(values=values), self.assertRaises(ValueError):
                self.vectorize(sample(values=values, n=1))

    def test_near_zero_policy_is_explicit(self):
        result = self.vectorize(sample(values=(-1e-12, 1e-12, 2e-10), n=3))
        self.assertEqual(result["values"].tolist(), [[0., 0., 2e-10]])

    def test_wrong_status_identity_basis_and_duplicate_are_rejected(self):
        for key, bad in (("status", "failed"), ("source_id", "other"), ("sample_id", "other"),
                         ("basis", [(0,), (0,)]), ("domain_size", -1)):
            value = sample()
            value["spectra"][0][key] = bad
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.vectorize(value)
        value = sample()
        value["spectra"].append(copy.deepcopy(value["spectra"][0]))
        with self.assertRaises(ValueError):
            self.vectorize(value)

    def test_units_field_operator_and_grid_mismatch(self):
        for key, bad in (("field", "F2"), ("inner_product", "weighted"),
                         ("operator", "ordinary_hodge"), ("grid", [1.])):
            value = sample()
            value[key] = bad
            with self.subTest(key=key), self.assertRaises(ValueError):
                self.vectorize(value)
        a, b = sample(), sample()
        b["sample_id"] = "s2"
        b["spectra"][0]["sample_id"] = "s2"
        b["scale_unit"] = "incompatible"
        with self.assertRaises(ValueError):
            spectral.vectorize_spectra([a, b], (0.,), degrees=(0,))

    def test_duplicate_sample_ids_and_invalid_configuration(self):
        with self.assertRaises(ValueError):
            spectral.vectorize_spectra([sample(), sample()], (0.,), degrees=(0,))
        for degrees, length in (((), 1), ((0, 0), 1), ((1, 0), 1), ((-1,), 1), ((0,), 0)):
            with self.subTest(degrees=degrees, length=length), self.assertRaises(ValueError):
                spectral.vectorize_spectra([], (0.,), degrees=degrees, length=length)

    def test_stage_basis_identity_and_row_order_are_preserved(self):
        value = sample(grid=(0., 1.))
        value["spectra"][1]["basis"] = [(10,), (20,)]
        with self.assertRaises(ValueError):
            self.vectorize(value, grid=(0., 1.))
        a, b = sample(), sample(values=(0., 3.))
        b["sample_id"], b["source_id"] = "b", "b-source"
        b["spectra"][0]["sample_id"], b["spectra"][0]["source_id"] = "b", "b-source"
        result = spectral.vectorize_spectra([b, a], (0.,), degrees=(0,), length=2)
        self.assertEqual(result["sample_ids"], ["b", "s"])
        self.assertEqual(result["values"].tolist(), [[0., 3.], [0., 2.]])
        self.assertEqual(result["bases"][0], [((0,), (1,))])

    def test_numeric_origin_must_be_compatible(self):
        a, b = sample(), sample()
        b["sample_id"] = "b"
        b["spectra"][0]["sample_id"] = "b"
        b["numeric_origin"]["zero_tolerance"] = 1e-3
        with self.assertRaises(ValueError):
            spectral.vectorize_spectra([a, b], (0.,), degrees=(0,))

    def test_input_and_output_own_their_values(self):
        value = sample(values=(2., 0.))
        before = copy.deepcopy(value)
        result = self.vectorize(value)
        result["values"][0, 0] = 17.
        self.assertEqual(value, before)


@unittest.skipUnless(HAS_NUMPY, "optional NumPy oracle environment required")
class MathematicalTests(unittest.TestCase):
    def setUp(self):
        self.path = {(0,): 0., (2,): 0., (1,): 1., (0, 1): 1., (1, 2): 1.}
        self.K, self.L = spectral.stage(self.path, 0), spectral.stage(self.path, 1)

    def test_path_and_filled_triangle_hand_matrices(self):
        np = spectral.numpy()
        np.testing.assert_allclose(spectral.persistent_laplacian(self.K, self.L, 0),
                                   [[.5, -.5], [-.5, .5]], atol=1e-12)
        filled = spectral.stage({**self.path, (0, 2): 1., (0, 1, 2): 1.}, 1)
        matrix = spectral.persistent_laplacian(self.K, filled, 0)
        np.testing.assert_allclose(matrix, [[1.5, -1.5], [-1.5, 1.5]], atol=1e-12)
        np.testing.assert_allclose(np.linalg.eigvalsh(matrix), [0, 3], atol=1e-12)
        self.assertEqual(spectral.exact_persistent_betti(self.K, self.L, 0), 1)
        self.assertEqual(spectral.exact_persistent_betti(self.K, filled, 0), 1)

    def test_persistent_h1_cycle_is_killed_by_a_filling(self):
        np = spectral.numpy()
        cycle = [[(0,), (1,), (2,)], [(0, 1), (0, 2), (1, 2)], []]
        filled = copy.deepcopy(cycle)
        filled[2] = [(0, 1, 2)]
        np.testing.assert_allclose(spectral.persistent_laplacian(cycle, filled, 1), 3 * np.eye(3))
        self.assertEqual(spectral.exact_persistent_betti(cycle, cycle, 1), 1)
        self.assertEqual(spectral.exact_persistent_betti(cycle, filled, 1), 0)

    def test_same_stage_reduces_to_hodge(self):
        matrix = spectral.persistent_laplacian(self.L, self.L, 0)
        spectral.numpy().testing.assert_array_equal(matrix, [[1, -1, 0], [-1, 2, -1], [0, -1, 1]])

    def test_singular_complement_and_disconnected_zeros(self):
        L = spectral.stage({**self.path, (3,): 1.}, 1)
        spectral.numpy().testing.assert_allclose(spectral.persistent_laplacian(self.K, L, 0),
                                                [[.5, -.5], [-.5, .5]])
        spectral.numpy().testing.assert_array_equal(spectral.persistent_laplacian(self.K, self.K, 0),
                                                   [[0, 0], [0, 0]])

    def test_empty_adjacent_degrees_keep_rectangular_shapes(self):
        self.assertEqual(spectral.boundary(self.K, 1).shape, (2, 0))
        self.assertEqual(spectral.boundary(self.K, 0).shape, (0, 2))
        self.assertEqual(spectral.persistent_laplacian(self.K, self.L, 1).shape, (0, 0))
        self.assertEqual(spectral.persistent_laplacian([[], [], []], self.L, 0).shape, (0, 0))

    def test_arbitrary_basis_order_and_input_preservation(self):
        before = copy.deepcopy(self.L)
        reordered = [list(reversed(cells)) for cells in self.L]
        expected = spectral.persistent_laplacian(self.K, self.L, 0)
        actual = spectral.persistent_laplacian(self.K, reordered, 0)
        spectral.numpy().testing.assert_allclose(actual, expected)
        K = [list(reversed(cells)) for cells in self.K]
        spectral.numpy().testing.assert_allclose(spectral.persistent_laplacian(K, reordered, 0), expected[::-1, ::-1])
        self.assertEqual(self.L, before)

    def test_noninclusion_invalid_complex_and_bad_filtration(self):
        with self.assertRaises(ValueError):
            spectral.persistent_laplacian(self.L, self.K, 0)
        with self.assertRaises(ValueError):
            spectral.stage({(0,): 1., (1,): 0., (0, 1): 0.}, 1)
        with self.assertRaises(ValueError):
            spectral.validate_complex([[(0,), (0,)], [], []])
        with self.assertRaises(ValueError):
            spectral.stage({(0,): float("nan")}, 0)
        with self.assertRaises(ValueError):
            spectral.stage(self.path, float("nan"))

    def test_cubical_chain_condition_and_square_spectrum(self):
        np = spectral.numpy()
        filtration = spectral.cubical_filtration([[0, 0], [0, 0]])
        square = spectral.stage(filtration, 0, "cubical")
        B1, B2 = spectral.boundary(square, 1, "cubical"), spectral.boundary(square, 2, "cubical")
        np.testing.assert_array_equal(B1 @ B2, np.zeros((4, 1)))
        np.testing.assert_allclose(np.linalg.eigvalsh(spectral.persistent_laplacian(square, square, 0, "cubical")),
                                   [0, 2, 2, 4], atol=1e-12)

    def test_maximum_endpoint_is_included_for_non_28_images(self):
        filtration = spectral.cubical_filtration([[0, 1, 2], [0, 1, 2]])
        end = spectral.stage(filtration, 2, "cubical")
        self.assertEqual(list(map(len, end)), [6, 7, 2])
        with self.assertRaises(ValueError):
            spectral.cubical_filtration([[float("nan")]])

    def test_ordinary_and_ph_equal_but_persistent_spectra_differ(self):
        np = spectral.numpy()
        other = {(0,): 0., (1,): 0., (2,): 1., (0, 1): 1., (1, 2): 1.}
        grid = (0., 1.)
        a = spectral.make_sample("a", "a-source", self.path, grid)
        b = spectral.make_sample("b", "b-source", other, grid)
        persistent = spectral.vectorize_spectra([a, b], grid)["values"]
        self.assertFalse(np.allclose(persistent[0], persistent[1]))
        ordinary = [spectral.make_sample("a", "a-source", self.path, grid, ordinary=True),
                    spectral.make_sample("b", "b-source", other, grid, ordinary=True)]
        result = spectral.vectorize_spectra(ordinary, grid, operator="ordinary_hodge")["values"]
        np.testing.assert_allclose(result[0], result[1])
        self.assertEqual(spectral.exact_persistent_betti(spectral.stage(other, 0), self.L, 0), 1)

    @unittest.skipUnless(HAS_SCIPY, "optional SciPy independent oracle required")
    def test_seeded_small_filtrations_against_independent_nullspace_and_exact_rank(self):
        np = spectral.numpy()
        rng = random.Random(111)
        for _ in range(20):
            filtration = {(v,): float(rng.randrange(3)) for v in range(4)}
            for a in range(4):
                for b in range(a + 1, 4):
                    filtration[(a, b)] = max(filtration[(a,)], filtration[(b,)], rng.randrange(4))
            for cell in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
                filtration[cell] = max(filtration[face] for face, _ in spectral.faces(cell, "simplicial"))
            stages = [spectral.stage(filtration, scale) for scale in range(4)]
            for q in (0, 1):
                for i in range(4):
                    for j in range(i, 4):
                        matrix = spectral.persistent_laplacian(stages[i], stages[j], q)
                        np.testing.assert_allclose(matrix, spectral.nullspace_laplacian(stages[i], stages[j], q),
                                                   atol=spectral.COMPARISON_TOLERANCE, rtol=0)
                        self.assertEqual(int(np.count_nonzero(np.abs(np.linalg.eigvalsh(matrix)) <= spectral.ZERO_TOLERANCE)),
                                         spectral.exact_persistent_betti(stages[i], stages[j], q))

    def test_resistance_equals_tree_distance_but_differs_on_a_cycle(self):
        self.assertEqual(spectral.graph_pair_geometry(self.path), (2., 2.))
        cycle = {(0,): 0., (2,): 0., (1,): 1., (3,): 1.,
                 (0, 1): 1., (1, 2): 1., (2, 3): 1., (0, 3): 1.}
        distance, resistance = spectral.graph_pair_geometry(cycle)
        self.assertEqual(distance, 2.)
        self.assertAlmostEqual(resistance, 1.)

    def test_disconnected_geometric_baseline_fails(self):
        with self.assertRaises(ValueError):
            spectral.graph_pair_geometry({(0,): 0., (1,): 0.})

    def test_higher_order_disk_hand_matrices_and_matched_q0(self):
        np = spectral.numpy()
        a, b = spectral.higher_order_samples(groups=1)
        Ka, Kb = spectral.stage(a["filtration"], 0), spectral.stage(b["filtration"], 0)
        L = spectral.stage(a["filtration"], 1)
        self.assertEqual(L, spectral.stage(b["filtration"], 1))
        np.testing.assert_allclose(spectral.persistent_laplacian(Ka, L, 1), 3 * np.eye(3), atol=1e-12)
        expected = np.array([[7, 2, -2], [2, 7, 2], [-2, 2, 7]]) / 3
        np.testing.assert_allclose(spectral.persistent_laplacian(Kb, L, 1), expected, atol=1e-12)
        np.testing.assert_allclose(np.linalg.eigvalsh(expected), [1, 3, 3], atol=1e-12)
        np.testing.assert_array_equal(spectral.persistent_laplacian(Ka, L, 0), spectral.persistent_laplacian(Kb, L, 0))
        self.assertEqual(spectral.exact_persistent_betti(Ka, L, 1), 0)
        self.assertEqual(spectral.exact_persistent_betti(Kb, L, 1), 0)

    def test_higher_order_dataset_matched_controls_and_no_q0_signal(self):
        np = spectral.numpy()
        grid = (0., 1.)
        samples = spectral.higher_order_samples(groups=4)
        for a, b in zip(samples[::2], samples[1::2]):
            permutation = {0: 1, 1: 2, 2: 3, 3: 0}
            graph_a = {tuple(sorted(permutation.get(v, v) for v in cell)): birth
                       for cell, birth in a["filtration"].items() if len(cell) <= 2}
            graph_b = {cell: birth for cell, birth in b["filtration"].items() if len(cell) <= 2}
            self.assertEqual(graph_a, graph_b)
            both = [spectral.make_sample(s["sample_id"], s["source_id"], s["filtration"], grid) for s in (a, b)]
            result = spectral.vectorize_spectra(both, grid)
            q0 = [c["degree"] == 0 for c in result["columns"]]
            q1 = [c["degree"] == 1 for c in result["columns"]]
            np.testing.assert_allclose(result["values"][0, q0], result["values"][1, q0], atol=1e-9)
            self.assertFalse(np.allclose(result["values"][0, q1], result["values"][1, q1]))
            ordinary = [spectral.make_sample(s["sample_id"], s["source_id"], s["filtration"], grid, ordinary=True) for s in (a, b)]
            values = spectral.vectorize_spectra(ordinary, grid, operator="ordinary_hodge")["values"]
            np.testing.assert_allclose(values[0], values[1], atol=1e-9)


@unittest.skipUnless(HAS_NUMPY and HAS_SCIPY and os.environ.get("SPECTRAL_WORKER"),
                     "set SPECTRAL_WORKER to the actual Rust exporter")
class RustExportTests(unittest.TestCase):
    grid = (0., 1 / 3, 2 / 3, 1.)

    def export(self, cells):
        fixture = {"sample_id": "sample", "source_id": "supplied", "filtration": cells}
        run = subprocess.run([os.environ["SPECTRAL_WORKER"], "--spectral"],
                             input=spectral.rust_request(fixture, self.grid),
                             capture_output=True, text=True, check=True, timeout=30)
        record = json.loads(run.stdout)
        spectral.verify_rust_sample(record, fixture, self.grid)
        return fixture, record

    def test_hand_path_and_basis_permutation(self):
        np = spectral.numpy()
        cells = {(0,): 0., (2,): 0., (1,): 1., (0, 1): 1., (1, 2): 1.}
        fixture, record = self.export(cells)
        pair = next(r for r in record["spectra"] if (r["degree"], r["i"], r["j"]) == (0, 0, 3))
        np.testing.assert_allclose(pair["values"], [0, 1], atol=1e-12)
        np.testing.assert_allclose(pair["matrix"], [.5, -.5, -.5, .5], atol=1e-12)
        _, reordered = self.export(dict(reversed(list(cells.items()))))
        np.testing.assert_array_equal(record["features"], reordered["features"])
        for field, bad in (("source_id", "different"), ("features", [0.] * 160)):
            damaged = copy.deepcopy(record)
            damaged[field] = bad
            with self.assertRaises((ValueError, AssertionError)):
                spectral.verify_rust_sample(damaged, fixture, self.grid)

    def test_empty_zero_and_repeated_spectra_keep_masks(self):
        _, empty = self.export({})
        self.assertEqual(sum(empty["validity_mask"]), 0)
        _, isolates = self.export({(0,): 0., (1,): 0.})
        self.assertEqual(sum(isolates["validity_mask"]), 20)
        self.assertEqual(set(isolates["features"]), {0.})
        _, triangle = self.export({(0,): 0., (1,): 0., (2,): 0., (0, 1): 0.,
                                   (0, 2): 0., (1, 2): 0., (0, 1, 2): 1.})
        q1 = next(r for r in triangle["spectra"] if (r["degree"], r["i"], r["j"]) == (1, 0, 3))
        spectral.numpy().testing.assert_allclose(q1["values"], [3, 3, 3], atol=1e-12)

    def test_worker_limits_fail_without_successful_partial_output(self):
        for request in ("x" * 65_537, "s source 0 1 2 3 257",
                        "s source 0 1 2 3 65\n" + "\n".join(f"0 1 {v}" for v in range(65))):
            run = subprocess.run([os.environ["SPECTRAL_WORKER"], "--spectral"], input=request,
                                 capture_output=True, text=True, timeout=30)
            self.assertNotEqual(run.returncode, 0)
            self.assertEqual(run.stdout, "")

    def test_seeded_existing_fixtures_use_real_rust_spectra(self):
        for fixture in spectral.synthetic_samples(groups=3) + spectral.higher_order_samples(groups=3):
            self.export(fixture["filtration"])


@unittest.skipUnless(HAS_NUMPY and HAS_SKLEARN, "optional external model environment required")
class ExportedModelTests(unittest.TestCase):
    def archive(self, path, compressed=False, **changes):
        np = spectral.numpy()
        fields = dict(sample_ids=np.array(["a", "b", "c", "d"]), labels=np.array([0, 1, 0, 1]),
                      groups=np.array([0, 1, 2, 3]), train=np.array([0, 1]), test=np.array([2, 3]),
                      branches=np.array(["persistent_spectral"]), seed=np.array(111),
                      persistent_spectral=np.array([[0.], [2.], [0.], [2.]]))
        save = np.savez_compressed if compressed else np.savez
        save(path, **{**fields, **changes})

    def reject_before_array_load(self, source, output, message):
        with patch.object(spectral.numpy(), "load", side_effect=AssertionError("array loading must not start")) as load:
            with self.assertRaisesRegex(ValueError, message):
                spectral.fit_exported_features(source, output)
            load.assert_not_called()
        self.assertFalse(output.exists())

    def test_compressed_oversize_members_including_unused_arrays_are_rejected(self):
        np = spectral.numpy()
        large = np.zeros((1024, 1024), dtype=np.float64)  # 8 MiB, compresses to kilobytes
        for name in ("persistent_spectral", "unused"):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as directory:
                source, output = Path(directory) / "input.npz", Path(directory) / "output.json"
                self.archive(source, compressed=True, **{name: large})
                self.assertLess(source.stat().st_size, 4 * 1024 * 1024)
                self.reject_before_array_load(source, output, "uncompressed byte limits")

    def test_aggregate_uncompressed_limit_precedes_array_loading(self):
        np = spectral.numpy()
        matrix = np.zeros((256, 512), dtype=np.float64)
        names = ["persistent_spectral"] + [f"feature_{i}" for i in range(15)]
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "input.npz", Path(directory) / "output.json"
            self.archive(source, compressed=True, branches=np.array(names), **dict.fromkeys(names, matrix))
            with zipfile.ZipFile(source) as archive:
                self.assertTrue(all(m.file_size < 2 * 1024 * 1024 for m in archive.infolist()))
                self.assertGreater(sum(m.file_size for m in archive.infolist()), 16 * 1024 * 1024)
            self.reject_before_array_load(source, output, "uncompressed byte limits")

    def test_forged_npy_shape_dtype_and_truncated_payload_are_rejected(self):
        cases = [("persistent_spectral", (2 ** 40, 2), "<f8"),
                 ("persistent_spectral", (-1, 2), "<f8"),
                 ("persistent_spectral", (4, 1), "|O"),
                 ("persistent_spectral", (4, 1), [("field", "<f8")]),
                 ("sample_ids", (4,), "<U1048576"),
                 ("persistent_spectral", (4, 1), "<f8")]
        for name, shape, dtype in cases:
            with self.subTest(shape=shape, dtype=dtype), tempfile.TemporaryDirectory() as directory:
                source, output = Path(directory) / "input.npz", Path(directory) / "output.json"
                self.archive(source)
                with zipfile.ZipFile(source) as archive:
                    members = {m.filename: archive.read(m) for m in archive.infolist()}
                header = io.BytesIO()
                spectral.numpy().lib.format.write_array_header_1_0(
                    header, {"descr": dtype, "fortran_order": False, "shape": shape})
                members[name + ".npy"] = header.getvalue()  # tiny file, no claimed payload
                with zipfile.ZipFile(source, "w", zipfile.ZIP_DEFLATED) as archive:
                    for key, value in members.items():
                        archive.writestr(key, value)
                self.reject_before_array_load(source, output, "shape or dtype|payload size")

    def test_npy_header_length_is_bounded_before_reading_header(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "input.npz", Path(directory) / "output.json"
            self.archive(source)
            with zipfile.ZipFile(source) as archive:
                members = {m.filename: archive.read(m) for m in archive.infolist()}
            members["persistent_spectral.npy"] = b"\x93NUMPY\x02\x00" + (2 ** 30).to_bytes(4, "little")
            with zipfile.ZipFile(source, "w", zipfile.ZIP_DEFLATED) as archive:
                for key, value in members.items():
                    archive.writestr(key, value)
            self.reject_before_array_load(source, output, "header exceeds 4096")

    def test_duplicate_and_excess_members_are_rejected_before_array_loading(self):
        for extra in (["labels.npy"], [f"extra_{i}.npy" for i in range(16)]):
            with self.subTest(extra=extra), tempfile.TemporaryDirectory() as directory:
                source, output = Path(directory) / "input.npz", Path(directory) / "output.json"
                self.archive(source)
                with warnings.catch_warnings(), zipfile.ZipFile(source, "a") as archive:
                    warnings.simplefilter("ignore", UserWarning)  # intentional duplicate fixture
                    for name in extra:
                        archive.writestr(name, b"")
                self.reject_before_array_load(source, output, "archive members")

    def test_unsupported_zip_compression_is_rejected_before_opening_members(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "input.npz", Path(directory) / "output.json"
            self.archive(source)
            with zipfile.ZipFile(source) as archive:
                members = {m.filename: archive.read(m) for m in archive.infolist()}
            with zipfile.ZipFile(source, "w", zipfile.ZIP_BZIP2) as archive:
                for key, value in members.items():
                    archive.writestr(key, value)
            with patch.object(zipfile.ZipFile, "open", side_effect=AssertionError("must not decompress")):
                self.reject_before_array_load(source, output, "stored or DEFLATE")

    def test_valid_compressed_v2_fortran_archive_still_fits(self):
        np = spectral.numpy()
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "input.npz", Path(directory) / "output.json"
            self.archive(source, compressed=True)
            with zipfile.ZipFile(source) as archive:
                members = {m.filename: archive.read(m) for m in archive.infolist()}
            values = io.BytesIO()
            np.lib.format.write_array(values, np.asfortranarray([[0., 0.], [2., 2.], [0., 0.], [2., 2.]]), version=(2, 0))
            members["persistent_spectral.npy"] = values.getvalue()
            with zipfile.ZipFile(source, "w", zipfile.ZIP_DEFLATED) as archive:
                for key, value in members.items():
                    archive.writestr(key, value)
            spectral.fit_exported_features(source, output)
            result = json.loads(output.read_text())
            self.assertEqual(result["scores"]["persistent_spectral"]["predictions"], [0, 1])
            self.assertEqual(result["input_sha256"], spectral.file_hash(source))

    def test_model_consumes_only_export_and_returns_actual_predictions(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "input.npz", Path(directory) / "output.json"
            self.archive(source)
            with patch.object(spectral, "make_sample", side_effect=AssertionError("no oracle input")):
                spectral.fit_exported_features(source, output)
            result = json.loads(output.read_text())
            self.assertEqual(result["input_sha256"], spectral.file_hash(source))
            self.assertEqual(result["test_sample_ids"], ["c", "d"])
            self.assertEqual(result["scores"]["persistent_spectral"]["predictions"], [0, 1])

    def test_loading_and_hashing_use_the_checked_snapshot(self):
        with tempfile.TemporaryDirectory() as directory:
            source, output = Path(directory) / "input.npz", Path(directory) / "output.json"
            self.archive(source)
            payload = spectral.checked_model_archive(source)
            digest = spectral.file_hash(source)
            def replace_path(path):
                path.write_bytes(b"changed after validation")
                return payload
            with patch.object(spectral, "checked_model_archive", side_effect=replace_path):
                spectral.fit_exported_features(source, output)
            result = json.loads(output.read_text())
            self.assertEqual(result["input_sha256"], digest)
            self.assertNotEqual(result["input_sha256"], spectral.file_hash(source))
            self.assertEqual(result["scores"]["persistent_spectral"]["predictions"], [0, 1])

    def test_bad_split_identity_shape_and_nonfinite_export_are_rejected(self):
        np = spectral.numpy()
        for changes in (dict(sample_ids=np.array(["a"] * 4)), dict(test=np.array([1, 3])),
                        dict(groups=np.array([0, 1, 0, 1])),
                        dict(persistent_spectral=np.zeros((3, 1))),
                        dict(persistent_spectral=np.full((4, 1), float("nan")))):
            with tempfile.TemporaryDirectory() as directory:
                source, output = Path(directory) / "input.npz", Path(directory) / "output.json"
                self.archive(source, **changes)
                with self.assertRaises(ValueError):
                    spectral.fit_exported_features(source, output)
                self.assertFalse(output.exists())


if __name__ == "__main__":
    unittest.main()
