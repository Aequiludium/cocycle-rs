"""Protect sample alignment and missing-data rejection without an ML dependency."""

import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

EXAMPLE = Path(__file__).resolve().parents[1] / "examples/analysis_learning.py"
spec = importlib.util.spec_from_file_location("analysis_learning", EXAMPLE)
example = importlib.util.module_from_spec(spec)
spec.loader.exec_module(example)


class AnalysisLearningTests(unittest.TestCase):
    def setUp(self):
        self.rows = []
        # Independently assembled path/ring teaching expectations.
        for n in (3, 4, 6, 5, 7):
            for label in ("path", "cycle"):
                self.rows.append(dict(zip(example.COLUMNS, [
                    f"{label}-{n}", "train" if n in (3, 4, 6) else "test", label,
                    "unweighted-graph-v1", n, n - (label == "path"),
                    2 - 2 / n if label == "path" else 2,
                    2 / n if label == "path" else 0])))

    def read(self, rows, columns=None):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "features.csv"
            with path.open("w", encoding="utf-8", newline="") as stream:
                writer = csv.DictWriter(stream, fieldnames=columns or example.COLUMNS)
                writer.writeheader()
                writer.writerows(rows)
            return example.read_samples(path)

    def test_reordered_rows_keep_identity_and_labels_together(self):
        rows = self.read(list(reversed(self.rows)))
        self.assertEqual([row["sample_id"] for row in rows],
                         [row["sample_id"] for row in reversed(self.rows)])
        self.assertEqual(rows[0]["label"], "cycle")
        self.assertEqual(rows[0]["split"], "test")

    def test_empty_missing_and_duplicate_samples_are_rejected(self):
        for rows in ([], self.rows[:-1], self.rows + [self.rows[0]]):
            with self.subTest(size=len(rows)), self.assertRaises(ValueError):
                self.read(rows)

    def test_mismatched_identity_split_units_and_features_are_rejected(self):
        for key, value in (("label", "cycle"), ("split", "test"),
                           ("source", "weighted"), ("sample_id", "unknown"),
                           ("vertex_count", "4"), ("edge_count", "0"),
                           ("leaf_fraction", "nan"), ("mean_degree", "inf"),
                           ("leaf_fraction", "")):
            rows = [dict(row) for row in self.rows]
            rows[0][key] = value
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                self.read(rows)

    def test_missing_feature_column_is_rejected(self):
        rows = [{key: value for key, value in row.items() if key != "leaf_fraction"}
                for row in self.rows]
        with self.assertRaises(ValueError):
            self.read(rows, example.COLUMNS[:-1])


if __name__ == "__main__":
    unittest.main()
