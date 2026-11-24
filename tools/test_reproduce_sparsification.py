"""Independent analytic and failure controls for the optional sampling tool."""

from fractions import Fraction
import importlib.util
import math
import unittest

import reproduce_sparsification as sampling


HAS_NUMPY = importlib.util.find_spec("numpy") is not None


class RationalOracleTests(unittest.TestCase):
    def test_rank_and_nonorthogonal_kernel(self):
        rank, kernel = sampling.exact_kernel([[1, 2, 3], [2, 4, 6]], 3)
        self.assertEqual(rank, 1)
        self.assertEqual(kernel, [[-2, 1, 0], [-3, 0, 1]])

    def test_empty_and_zero_rank(self):
        self.assertEqual(sampling.exact_kernel([], 0), (0, []))
        self.assertEqual(sampling.exact_kernel([], 2), (0, [[1, 0], [0, 1]]))

    def test_fractional_rank_is_exact(self):
        _, pivots = sampling.rref([[Fraction(1, 3), Fraction(2, 3)],
                                  [Fraction(1, 2), Fraction(1)]], 2)
        self.assertEqual(pivots, [0])

    def test_ragged_matrix_is_rejected(self):
        with self.assertRaises(ValueError):
            sampling.rref([[1], [1, 2]], 2)


@unittest.skipUnless(HAS_NUMPY, "optional NumPy environment required")
class SamplingControls(unittest.TestCase):
    def test_tetrahedron_and_literal_failures(self):
        result = sampling.controls()
        self.assertEqual(result["exact_scores"], ["3/4"] * 4)
        self.assertEqual(result["section43_moments"][3], [1.5] * 4)
        self.assertEqual(result["literal_with_dual_gram_only"][3], [0.5] * 4)
        self.assertTrue(result["one_face_failure"]["rank_loss"])
        self.assertAlmostEqual(result["one_face_failure"]["maximum_generalized"], 3)

    def test_nonuniform_exact_scores_and_zero_column(self):
        np = sampling.numpy()
        boundary = np.array([[1, 2, 0], [0, 0, 0]])
        tau, rank = sampling.exact_scores(boundary)
        self.assertEqual(rank, 1)
        self.assertEqual(tau, [Fraction(1, 5), Fraction(4, 5), 0])
        np.testing.assert_allclose(sampling.numerical_scores(boundary), [0.2, 0.8, 0], atol=1e-14)

    def test_no_integer_truncation_in_exact_oracle(self):
        np = sampling.numpy()
        for matrix in ([[0.5]], [[math.nan]], [[math.inf]]):
            with self.assertRaises(ValueError):
                sampling.exact_scores(np.array(matrix))

    def test_orientation_does_not_change_scores(self):
        np = sampling.numpy()
        b = np.array([[1, 2, 0], [0, 1, 1]])
        tau, rank = sampling.exact_scores(b)
        changed, changed_rank = sampling.exact_scores(b * np.array([-1, 1, -1]))
        self.assertEqual((changed, changed_rank), (tau, rank))

    def test_repeated_draws_use_weights_once(self):
        np = sampling.numpy()
        b = np.array([[1., 2.]])
        gram, draws, weights, counts = sampling.sample_gram(b, [0.2, 0.8], 30, 0)
        self.assertEqual(sum(counts), 30)
        self.assertEqual(len(draws), 30)
        self.assertAlmostEqual(gram[0, 0], 5)
        self.assertAlmostEqual(weights[0], counts[0] / 6)
        self.assertAlmostEqual(weights[1], counts[1] / 24)

    def test_zero_rank_never_normalizes(self):
        np = sampling.numpy()
        for shape in ((0, 0), (3, 0), (3, 2)):
            gram, draws, weights, _ = sampling.sample_gram(np.zeros(shape), np.zeros(shape[1]), 2, 0)
            self.assertEqual(gram.shape, (shape[0], shape[0]))
            self.assertEqual(draws, [])
            self.assertFalse(weights.any())

    def test_invalid_probabilities_and_draws_are_rejected(self):
        np = sampling.numpy()
        for probabilities in ([0., 0.], [1., 0.], [-1., 2.], [math.nan, 1.], [0.5], [1., 1.]):
            with self.subTest(probabilities=probabilities), self.assertRaises(ValueError):
                sampling.sample_gram(np.array([[1., 1.]]), probabilities, 3, 0)
        for count in (0, -1, 1.5, True):
            with self.subTest(count=count), self.assertRaises(ValueError):
                sampling.sample_gram(np.array([[1., 1.]]), [0.5, 0.5], count, 0)

    def test_polynomial_moments_against_scalar_formulas(self):
        np = sampling.numpy()
        values = np.array([0., 0.25, 0.5, 1.])
        h = np.diag(values)
        result = sampling.kid_moments(h, np.ones((4, 3)), 5)
        np.testing.assert_allclose(result[1], 2 * values, atol=1e-14)
        np.testing.assert_allclose(result[3], 2 * (4 * values**3 - 3 * values), atol=1e-14)
        np.testing.assert_allclose(result[5], 2 * (16 * values**5 - 20 * values**3 + 5 * values), atol=1e-14)
        self.assertFalse(result[2].any())

    def test_moment_invalid_shapes_fail(self):
        np = sampling.numpy()
        for h, z, degree in ((np.zeros((2, 3)), np.ones((2, 1)), 3),
                             (np.eye(2), np.ones((3, 1)), 3),
                             (np.eye(2), np.ones((2, 0)), 3),
                             (np.eye(2), np.ones((2, 1)), 0)):
            with self.assertRaises(ValueError):
                sampling.kid_moments(h, z, degree)

    def test_exact_bin_integral_is_independent_of_partition(self):
        np = sampling.numpy()
        moments = np.array([[0., 0.], [1., 2.]])
        for bins in (1, 16, 64):
            scores, _ = sampling.analytic_bin_scores(moments, bins)
            np.testing.assert_allclose(scores, [2 / math.pi, 4 / math.pi], atol=1e-14)

    def test_negative_estimates_are_retained(self):
        np = sampling.numpy()
        scores, _ = sampling.analytic_bin_scores(np.array([[0.], [-1.]]))
        self.assertLess(scores[0], 0)

    def test_equal_ordinary_spectra_do_not_prove_psd_order(self):
        np = sampling.numpy()
        b = np.diag([1, 2])
        g, whitening, kernel = sampling.comparison_setup(b, 2)
        diagnostic = sampling.assess(g, np.diag([4., 1.]), whitening, kernel, 2)
        self.assertAlmostEqual(diagnostic["minimum_generalized"], 0.25)
        self.assertAlmostEqual(diagnostic["maximum_generalized"], 4)
        self.assertFalse(diagnostic["psd_sandwich_pass"]["0.5"])


if __name__ == "__main__":
    unittest.main()
