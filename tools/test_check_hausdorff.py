"""Failure controls for the independent Hausdorff output checker."""
import copy
import json
import math
from pathlib import Path
import tempfile
import unittest

import check_hausdorff as checker


class HausdorffCheckerTests(unittest.TestCase):
    def fixtures(self):
        rows = [{"name": name, "dimension": 1, "first": [0.0],
                 "second": [0.0], "distance": 0.0} for name in sorted(checker.NAMES)]
        by_name = {r["name"]: r for r in rows}
        for name, squared in checker.HAND.items():
            by_name[name].update(second=[math.sqrt(squared)], distance=math.sqrt(squared))
        by_name["stress_large"].update(second=[1e200], distance=1e200)
        by_name["stress_tiny"].update(second=[1e-200], distance=1e-200)
        by_name["stress_subnormal"].update(second=[float.fromhex("0x0.0000000000001p-1022")],
                                          distance=float.fromhex("0x0.0000000000001p-1022"))
        by_name["stress_final_overflow"].update(first=[-float.fromhex("0x1.fffffffffffffp+1023")],
                                               second=[float.fromhex("0x1.fffffffffffffp+1023")],
                                               distance="numerical_failure")
        rows.append({"matrix_ids": ["c", "a", "b", "same"],
                     "matrix_values": [3.0, 2.0, 1.0, 3.0, 0.0, 1.0]})
        return rows

    def run_rows(self, rows):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rust.jsonl"
            path.write_text("".join(json.dumps(r) + "\n" for r in rows), encoding="utf-8")
            return checker.check(path)

    def test_complete_output_passes_and_tiny_zero_fails(self):
        rows = self.fixtures()
        self.assertEqual(self.run_rows(rows)["pairs"], 90)
        tiny = next(r for r in rows if r["name"] == "stress_tiny")
        tiny["distance"] = 0.0
        with self.assertRaisesRegex(ValueError, "distance mismatch"):
            self.run_rows(rows)

    def test_missing_duplicate_and_identity_changes_fail(self):
        original = self.fixtures()
        for change in (lambda rows: rows.pop(0),
                       lambda rows: rows.__setitem__(0, copy.deepcopy(rows[1])),
                       lambda rows: rows[-1]["matrix_ids"].reverse()):
            rows = copy.deepcopy(original)
            change(rows)
            with self.assertRaises(ValueError):
                self.run_rows(rows)

    def test_wrong_overflow_or_nonfinite_value_fails(self):
        for name, value in [("stress_final_overflow", 0.0),
                            ("stress_large", "numerical_failure"),
                            ("single_345", float("nan"))]:
            rows = self.fixtures()
            next(r for r in rows if r["name"] == name)["distance"] = value
            with self.assertRaises(ValueError):
                self.run_rows(rows)

    def test_invalid_coordinate_shape_fails(self):
        rows = self.fixtures()
        rows[0]["dimension"] = 3
        with self.assertRaisesRegex(ValueError, "shape"):
            self.run_rows(rows)


if __name__ == "__main__":
    unittest.main()
