"""Paired statistics must preserve constant ratios and reject incomplete pairs."""
import unittest
from profile_rips import paired_summary


class PairedSummaryTests(unittest.TestCase):
    def test_identical_and_scaled_rounds_preserve_ratio(self):
        for factor in (1., 2., .5):
            result = paired_summary([1., 2., 3., 4.], [factor * x for x in (1., 2., 3., 4.)], 7)
            self.assertEqual(result['median_ratio'], factor)
            self.assertEqual(result['paired_bootstrap_lower95'], factor)
            self.assertEqual(result['paired_bootstrap_upper95'], factor)

    def test_pair_resampling_is_reproducible(self):
        self.assertEqual(paired_summary([1., 2., 3.], [2., 1., 4.], 23),
                         paired_summary([1., 2., 3.], [2., 1., 4.], 23))

    def test_missing_or_nonpositive_samples_are_rejected(self):
        for a, b in (([], []), ([1.], []), ([1.], [1., 2.]), ([0.], [1.]), ([1.], [-1.]),
                     ([float('nan')], [1.]), ([1.], [float('inf')])):
            with self.assertRaises(ValueError):
                paired_summary(a, b, 7)


if __name__ == '__main__':
    unittest.main()
