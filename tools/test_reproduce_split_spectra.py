"""Analytic controls and failure contracts for the optional #145 reproduction."""

import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import reproduce_split_spectra as spectral


HAS_STACK = all(importlib.util.find_spec(name) for name in ("numpy", "scipy", "sympy"))


@unittest.skipUnless(HAS_STACK, "optional NumPy/SciPy/SymPy reproduction environment")
class SplitSpectrumTests(unittest.TestCase):
    def setUp(self):
        self.np, self.la = spectral.libraries()

    def test_triangle_and_shell_analytic_grams(self):
        np = self.np
        for name, expected in [("wide_triangle", [0., 3., 3.]),
                               ("triangle_top", [3.]), ("shell", [0., 4., 4., 4.])]:
            with self.subTest(name=name):
                result = spectral.gram_flip(spectral.fixtures()[name])
                spectral.check_close(result["reconstructed"], np.array(expected))
                self.assertLess(result["lifted_residual"], 1e-12)
                self.assertEqual(result["original_domain"], len(expected))

    def test_zero_and_empty_gram_domains_are_preserved(self):
        np = self.np
        for shape in [(0, 3), (3, 0), (0, 0), (2, 5)]:
            with self.subTest(shape=shape):
                result = spectral.gram_flip(np.zeros(shape))
                self.assertEqual(result["reconstructed"].shape, (shape[1],))
                self.assertEqual(result["numerical_rank"], 0)
                self.assertEqual(result["lifted_residual"], 0.)

    def test_equal_component_roots_keep_combined_multiplicity(self):
        np = self.np
        controls = dict(spectral.pair_controls())
        triangle = controls["filled_triangle"]
        spectral.check_close(self.la.eigvalsh(triangle["U"]), np.array([0., 0., 3.]))
        spectral.check_close(self.la.eigvalsh(triangle["D"]), np.array([0., 3., 3.]))
        spectral.check_close(triangle["union"], np.array([3., 3., 3.]))
        self.assertEqual(len(triangle["union"]), 3)
        persistent = controls["persistent_path_and_disconnected_triangle"]
        spectral.check_close(persistent["union"], np.array([0., 0., 1., 3., 3.]))
        self.assertEqual(persistent["B2"].shape, (4, 8))

    def test_nonorthogonal_harmonic_basis_is_normalized(self):
        np = self.np
        h = np.diag([0., 0., 2.])
        kernel = np.array([[2., 1.], [0., 3.], [0., 0.]])
        answer = spectral.harmonic_restriction(h, kernel)
        spectral.check_close(answer["spectrum"], np.array([2.]))
        spectral.check_close(answer["kernel"].T @ answer["kernel"], np.eye(2))

    def test_nonharmonic_representative_is_rejected(self):
        np = self.np
        h = np.array([[1., -1.], [-1., 1.]])
        with self.assertRaisesRegex(ValueError, "not harmonic"):
            spectral.harmonic_restriction(h, np.array([[1.], [0.]]))
        result = spectral.harmonic_restriction(h, np.ones((2, 1)))
        spectral.check_close(result["spectrum"], np.array([2.]))
        self.assertLess(result["residual"], 1e-12)

    def test_incomplete_dependent_and_wrong_shape_kernels_fail(self):
        np = self.np
        h = np.diag([0., 1.])
        for kernel in [np.empty((2, 0)), np.array([[1., 1.], [0., 0.]]), np.ones((3, 1))]:
            with self.subTest(kernel=kernel), self.assertRaises(ValueError):
                spectral.harmonic_restriction(h, kernel)

    def test_zero_operator_restricts_to_a_genuine_empty_domain(self):
        np = self.np
        result = spectral.harmonic_restriction(np.zeros((3, 3)), np.eye(3))
        self.assertEqual(result["reduced"].shape, (0, 0))
        self.assertEqual(result["spectrum"].size, 0)
        self.assertFalse(result["ambiguous"])

    def test_tiny_positive_is_retained_and_rank_ambiguity_is_reported(self):
        np = self.np
        result = spectral.harmonic_restriction(np.diag([0., 1e-14, 1.]), np.eye(3)[:, :1])
        self.assertEqual(len(result["spectrum"]), 2)
        self.assertEqual(result["spectrum"][0], 1e-14)
        self.assertTrue(result["ambiguous"])

    def test_rayleigh_ritz_matches_the_generalized_eigenproblem(self):
        np = self.np
        h = np.diag([1., 2., 5., 7.])
        s = np.array([[1., 0.], [1., 2.], [0., 1.], [2., 1.]])
        c, values, _ = spectral.rayleigh_ritz(h, s, 2)
        spectral.check_close(c.T @ (s.T @ s) @ c, np.eye(2))
        spectral.check_close(c.T @ (s.T @ h @ s) @ c, np.diag(values))
        expected = self.la.eigvalsh(s.T @ h @ s, s.T @ s)
        spectral.check_close(values, expected)

    def test_dependent_search_fails_and_correction_is_explicit(self):
        np = self.np
        h = np.diag([1., 2., 5.])
        s = np.column_stack((np.eye(3)[:, 0], np.eye(3)[:, 0], np.eye(3)[:, 1]))
        with self.assertRaisesRegex(ValueError, "dependent|ill-conditioned"):
            spectral.rayleigh_ritz(h, s, 2)
        c, values, _ = spectral.rayleigh_ritz(h, s, 2, compress=True)
        spectral.check_close(values, np.array([1., 2.]))
        spectral.check_close((s @ c).T @ (s @ c), np.eye(2))
        with self.assertRaisesRegex(ValueError, "insufficient search rank"):
            spectral.rayleigh_ritz(h, s, 3, compress=True)

    def test_cholesky_failure_cannot_become_success(self):
        np = self.np
        with patch.object(self.la, "cholesky", side_effect=self.la.LinAlgError("failed pivot")):
            with self.assertRaises(self.la.LinAlgError):
                spectral.basic_lobpcg(np.diag([1., 2., 3.]), np.eye(3)[:, :1], 1)

    def test_basic_iteration_has_independent_known_solution_and_fresh_residuals(self):
        np = self.np
        h = np.diag([1., 1., 2., 2., 2., 2.])
        initial = np.arange(1., 7.)[:, None]
        trace = []
        result = spectral.basic_lobpcg(h, initial, 1, trace=trace)
        spectral.check_close(result["values"], np.array([1.]))
        self.assertGreater(result["iterations"], 0)
        self.assertFalse(result["index_coverage_certified"])
        for step in trace:
            spectral.check_close(step["R"], h @ step["X"] - step["X"] * step["eigenvalues"])
        self.assertIn("S", trace[0])
        self.assertIn("cholesky", trace[0]["ritz"])

    def test_exhaustion_and_excessive_count_return_errors(self):
        np = self.np
        h = np.diag([1., 2., 3., 4.])
        for count, iterations in [(2, 0), (1, 0)]:
            with self.subTest(count=count), self.assertRaises(ValueError):
                spectral.basic_lobpcg(h, np.ones((4, 1)), count, max_iterations=iterations)

    def test_invalid_matrices_are_rejected(self):
        for h in [[[1., 1.], [0., 1.]], [[float("nan")]], [[1., 2.]]]:
            with self.subTest(h=h), self.assertRaises(ValueError):
                spectral.harmonic_restriction(h, self.np.empty((len(h), 0)))
        with self.assertRaises(ValueError):
            spectral.gram_flip(self.np.array([[1 + 1j]]))

    def test_json_keeps_empty_matrix_shapes_and_output_is_never_overwritten(self):
        np = self.np
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            spectral.write_json(root / "matrix.json", np.zeros((0, 7)))
            self.assertEqual(json.loads((root / "matrix.json").read_text()),
                             {"shape": [0, 7], "data": []})
            with self.assertRaises(FileExistsError):
                spectral.run(root)


class PaperIdentityTests(unittest.TestCase):
    def test_modified_pdf_is_rejected_before_numerical_work(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            for name in spectral.PAPERS:
                (root / name).write_bytes(b"unrelated paper")
            with self.assertRaisesRegex(ValueError, "unexpected paper hash"):
                spectral.reference_manifest(root)


if __name__ == "__main__":
    unittest.main()
