"""Independent controls for the SPSC-A reproduction and failure semantics."""
import importlib.util
import unittest

import reproduce_persistent_constructions as paper


@unittest.skipUnless(importlib.util.find_spec("numpy") and importlib.util.find_spec("scipy"),
                     "optional NumPy/SciPy reference environment")
class ConstructionTests(unittest.TestCase):
    def test_nonorthogonal_kernel_requires_induced_gram(self):
        import numpy as np
        # Kernel basis vectors share a coordinate: Z Z^T is not a projector.
        b1 = np.eye(3, dtype=int)
        b2 = np.array([[1, 1, 1]])
        a = np.zeros((0, 3), dtype=int)
        result, schur, detail = paper.constructions(b1, b2, a)
        expected = np.eye(3) - np.ones((3, 3)) / 3
        np.testing.assert_allclose(result, expected, atol=1e-14)
        np.testing.assert_allclose(schur, expected, atol=1e-14)
        self.assertFalse(np.allclose(detail["Z"] @ detail["Z"].T, expected))

    def test_path_has_actual_schur_correction(self):
        import numpy as np
        result, schur, _ = paper.constructions(np.array([[-1, 0], [0, 1]]),
                                               np.array([[1, -1]]), np.zeros((0, 2), dtype=int))
        np.testing.assert_allclose(result, [[.5, -.5], [-.5, .5]])
        np.testing.assert_allclose(schur, result)
        self.assertFalse(np.allclose(result, np.eye(2)))

    def test_lsmr_each_iterate_and_breakdown(self):
        records = paper.solver_controls()
        self.assertEqual(len(records), 12)
        self.assertTrue(any(len(record["trace"]) == 1 for record in records))
        self.assertTrue(any(not record["trace"] for record in records))

    def test_lsmr_iteration_exhaustion_is_not_a_projection(self):
        import numpy as np
        with self.assertRaisesRegex(ArithmeticError, "iteration limit"):
            paper.lsmr(np.diag([1., 2., 4.]), np.ones(3), max_iterations=1)
        with self.assertRaises(ValueError):
            paper.lsmr(np.eye(2), [float("nan"), 0])

    def test_signed_permutations(self):
        self.assertEqual(len(paper.orientation_controls()), 48)

    def test_rank_ambiguity_is_explicit(self):
        import numpy as np
        with self.assertRaisesRegex(ArithmeticError, "ambiguous"):
            paper.row_projector(np.array([[1., 0.], [1., 1e-13]]))


if __name__ == "__main__":
    unittest.main()
