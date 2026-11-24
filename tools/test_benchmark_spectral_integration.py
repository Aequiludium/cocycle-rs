"""Regression controls for baseline selection and independent output acceptance."""
import unittest

from benchmark_spectral_integration import validate_specification, validate_values


class BaselineContracts(unittest.TestCase):
    def test_original_pr_baseline_is_rejected(self):
        with self.assertRaisesRegex(ValueError, "PR-head baselines"):
            validate_specification({"cells": [{"name": "wrong", "baseline": "a" * 40}]})

    def test_identical_algorithm_routes_are_rejected(self):
        with self.assertRaisesRegex(ValueError, "must differ"):
            validate_specification({"cells": [{"name": "same", "routes": {
                "baseline": {"args": ["dense"]}, "candidate": {"args": ["dense"]}}}]})

    def test_explicit_unaccelerated_routes_are_accepted(self):
        validate_specification({"cells": [{"name": "action", "routes": {
            "baseline": {"args": ["dense"]}, "candidate": {"args": ["factorized"]}},
            "baseline_description": "dense Schur and action",
            "candidate_description": "factorized action", "oracle_file": "target/oracle.json"}]})

    def test_oracle_rejects_wrong_shape_nonfinite_and_changed_values(self):
        for actual in [None, [], [float("nan")], [float("inf")], [2.]]:
            self.assertFalse(validate_values(actual, [1.], 1e-10))
        self.assertTrue(validate_values([1. + 1e-12], [1.], 1e-10))
        self.assertFalse(validate_values([1e308], [-1e308], 1e-10))


if __name__ == "__main__":
    unittest.main()
