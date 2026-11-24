"""Independent polynomial, ordering and active/locked shape controls."""

import importlib.util
import unittest

import reproduce_subspaces as subject


@unittest.skipUnless(importlib.util.find_spec("numpy"), "optional NumPy reproduction")
class SubspaceTests(unittest.TestCase):
    def test_defining_polynomial_detects_both_printed_current_iterate_recursions(self):
        rows = subject.polynomial_controls()
        self.assertLess(max(r["corrected_max_error"] for r in rows), 1e-12)
        self.assertGreater(rows[2]["printed_max_error"], 0.1)

    def test_mixed_degrees_remain_attached_to_original_columns(self):
        np = subject.numpy()
        a = np.diag([0., 1., 2., 3.])
        y = np.ones((4, 5))
        degrees = [8, 0, 3, 1, 2]
        actual = subject.chebyshev_filter(a, y, degrees, 1., 3., 0.)
        for j, degree in enumerate(degrees):
            coefficients = [0.] * degree + [1.]
            expected = np.polynomial.chebyshev.chebval(np.arange(4.) - 2., coefficients)
            expected /= np.polynomial.chebyshev.chebval(-2., coefficients)
            np.testing.assert_allclose(actual[:, j], expected, atol=1e-12)

    def test_fft_ties_first_item_and_distances_beyond_printed_sentinel(self):
        np = subject.numpy()
        fields = [np.full((4, 4), v) for v in [0., 10000., 10000., 30000.]]
        order, _ = subject.fft_order(fields)
        self.assertEqual(order, [0, 1, 2, 3])
        self.assertEqual(subject.fft_order(np.empty((0, 4, 4)))[0], [])
        with self.assertRaises(ValueError):
            subject.fft_order([np.full((4, 4), float("nan"))])

    def test_active_locked_shapes_and_degenerate_cluster_projector(self):
        np = subject.numpy()
        a = np.diag([0.4, 0.4, 1.1, 3., 4., 5., 6., 7.])
        initial = np.random.default_rng(148).normal(size=(8, 5))
        result = subject.chfsi(a, initial, 3, 2., 7.1, 0.3)
        vectors = np.array(result["vectors"])
        np.testing.assert_allclose(vectors @ vectors.T, np.diag([1., 1., 1., 0., 0., 0., 0., 0.]), atol=1e-9)
        for row in result["trace"]:
            locked, active = row["locked_shape"][1], row["active_shape"][1]
            self.assertEqual(locked + active, 5)
            self.assertEqual(row["qr_shape"], [8, active])
            self.assertEqual(row["ritz_shape"], [active, active])
        with self.assertRaisesRegex(ValueError, "rank"):
            subject.chfsi(a, np.ones((8, 5)), 3, 2., 7.1, 0.3)
        with self.assertRaisesRegex(ValueError, "convergence"):
            subject.chfsi(a, initial, 3, 2., 7.1, 0.3, max_iterations=1)

    def test_nonconstant_paper_field_requires_a_metric_and_magnitude_order(self):
        result = subject.model_metric_control()
        self.assertGreater(result["euclidean_symmetry_defect"], 1.)
        self.assertLess(result["original_residual"], 1e-10)
        self.assertEqual(sorted(result["values_by_magnitude"], key=abs), result["values_by_magnitude"])

    def test_transfer_checks_identity_orientation_projection_rank_and_caps(self):
        np = subject.numpy()
        old = {"source": "triangle", "degree": 1, "metric": "euclidean",
               "basis": [((0, 1), 1), ((0, 2), 1), ((1, 2), 1)]}
        new = {**old, "basis": [((1, 2), -1), ((0, 1), 1)]}
        block = np.array([[1., 0.], [0., 1.], [2., 0.]])
        mapped, rank = subject.transfer_block(block, old, new, 1, 20)
        expected = np.array([-2., 1.]) / np.sqrt(5.)
        np.testing.assert_allclose(mapped[:, 0], expected)
        self.assertEqual(rank, 1)
        for bad in ({**new, "source": "another-owner"}, {**new, "degree": 0},
                    {**new, "metric": "weighted"}, {**new, "basis": [((0, 1), 0)]}):
            with self.assertRaises(ValueError):
                subject.transfer_block(block, old, bad, 1, 20)
        with self.assertRaisesRegex(ValueError, "storage"):
            subject.transfer_block(block, old, new, 1, 1)
        with self.assertRaisesRegex(ValueError, "rank insufficient"):
            subject.transfer_block(block, old, new, 2, 20, np.diag([1., 0.]))
        mapped, _ = subject.transfer_block(block, old, new, 1, 20, np.diag([1., 0.]))
        self.assertLess(abs(mapped[1, 0]), 1e-12)
        self.assertEqual(subject.transfer_block(block, old, new, 0, 20)[0].shape, (2, 0))


if __name__ == "__main__":
    unittest.main()
