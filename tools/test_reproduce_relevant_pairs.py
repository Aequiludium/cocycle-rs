"""Exact regressions for Algorithm 2 and coordinate reconstruction."""
from fractions import Fraction
import importlib.util
from pathlib import Path
import tempfile
import unittest

import reproduce_relevant_pairs as reproduction


class EventTests(unittest.TestCase):
    def test_distinct_buckets_jointly_trigger_a_pair(self):
        labels = [(2, 3), (1, 4)]
        selected = reproduction.selected_pairs(labels, 5)
        self.assertIn((2, 4), selected)
        self.assertNotIn((2, 4), labels)
        self.assertEqual(reproduction.predecessor(selected, 5, 5), (2, 4))

    def test_endpoint_precedes_start_in_predecessor_order(self):
        self.assertEqual(reproduction.predecessor([(3, 3), (1, 4)], 3, 5), (1, 4))
        self.assertIsNone(reproduction.predecessor([(3, 3)], 1, 4))

    def test_zero_columns_are_retained_with_invertible_transform(self):
        r, v = reproduction.reduce_boundary([[-1, -1, 0], [1, 0, -1], [0, 1, 1]], 3)
        self.assertEqual(r[2], [0, 0, 0])
        self.assertEqual(v[2], [Fraction(1), Fraction(-1), Fraction(1)])

    def test_fractional_pivot_is_not_a_finite_field_pivot(self):
        r, v = reproduction.reduce_boundary([[1, 0, 1], [0, 2, 1]], 3)
        self.assertEqual(v[2], [Fraction(-1), Fraction(-1, 2), Fraction(1)])
        self.assertEqual(r[2], [0, 0])

    @unittest.skipUnless(importlib.util.find_spec("numpy"), "optional NumPy dependency")
    def test_route_comparison_rejects_equal_checksum_wrong_matrices_and_domains(self):
        output = ["RELEVANT 0 1 2 None", "MATRIX 1 -1 -1 1", "SPECTRUM 0 2",
                  "RESIDUAL 0", "EMBEDDING 0", "BASIS [[0], [1]]", "CONTEXT 0 0 1"]
        reproduction.compare_full_outputs(output, output)
        for index, wrong in [(1, "MATRIX 1 1 -1 -1"), (2, "SPECTRUM 2"),
                             (5, "BASIS [[1], [0]]"), (6, "CONTEXT 0 0 2")]:
            changed = output.copy()
            changed[index] = wrong
            with self.assertRaises(AssertionError):
                reproduction.compare_full_outputs(output, changed)

    @unittest.skipUnless(importlib.util.find_spec("numpy") and importlib.util.find_spec("gudhi"),
                         "optional NumPy/GUDHI dependencies")
    def test_all_pairs_include_joint_trigger_and_larger_requested_domain(self):
        import json
        name, cells = next(c for c in reproduction.controls() if c[0] == "joint-events-delayed-vertex")
        with tempfile.TemporaryDirectory() as directory:
            reproduction.check_case(name, cells, 0, Path(directory))
            payload = json.loads((Path(directory) / f"{name}-q0.json").read_text())
        final = payload["pairs"][-1]
        self.assertEqual(final["representative"], [2, 4])
        self.assertEqual(final["size"], 5)
        self.assertEqual(len(final["upper"]), 5)
        self.assertEqual(final["upper"][-1][-1], "0")


if __name__ == "__main__":
    unittest.main()
