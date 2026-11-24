"""Independent controls for the SPSC-A reproduction and failure semantics."""
import importlib.util
import unittest
import hashlib
import json
from pathlib import Path
import tempfile

import reproduce_persistent_constructions as paper


class FrozenReplayTests(unittest.TestCase):
    def test_frozen_expected_values_cannot_be_silently_replaced(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            manifest = {}
            for name in ["index.json", "solver-controls.json"]:
                raw = b"[]\n"
                (root/name).write_bytes(raw)
                manifest[name] = hashlib.sha256(raw).hexdigest()
            (root/"manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
            self.assertEqual(paper.read_frozen_lsmr(root)[1:], ([], []))
            (root/"solver-controls.json").write_text("[1]\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "hash mismatch"):
                paper.read_frozen_lsmr(root)

    def test_worker_failure_is_retained_and_truncation_is_rejected(self):
        result = paper.parse_lsmr_replay('CASE limited\nERROR NumericalFailure { context: "LSMR iteration limit" }\nEND\n')
        self.assertIn("iteration limit", result[0]["error"])
        self.assertNotIn("solution", result[0])
        with self.assertRaisesRegex(ValueError, "unterminated"):
            paper.parse_lsmr_replay("CASE truncated\n")


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

    def test_nonfinite_and_broadcast_comparisons_are_rejected(self):
        import numpy as np
        for actual, expected in [(np.array([float("nan")]), np.array([0.])),
                                 (np.array([1.]), np.array([1.,1.]))]:
            with self.assertRaisesRegex(AssertionError, "shape mismatch or nonfinite"):
                paper.check_close(actual, expected, "worker output")


if __name__ == "__main__":
    unittest.main()
