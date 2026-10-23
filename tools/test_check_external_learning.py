"""Identity, boundary, layout and real-process regressions for #92."""

import copy
import importlib.util
import json
import os
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import check_external_learning as learning


def hand_example():
    sample = {"sample_key": "edge-a", "sample_order": 0, "snapshot": "fixture-1",
              "vertices": [[10, -1., -3., 0.], [20, 0., 5., 1.]],
              "edges": [[10, 20, 0.]], "grid": [-2., 0.], "universe": [10, 20]}
    fixture = {"samples": [sample], "labels": [{"sample_key": "edge-a", "label": 1, "split": "train"}]}
    full = {"format": "cocycle-graph-learning-v1", "sample_key": "edge-a", "snapshot": "fixture-1",
            "sample_order": 0, "selection": {"kind": "topology", "index": None, "scale": None},
            "coverage": "complete-supplied-source", "domain": "signed-integer", "scale_unit": "fixture-time",
            "topology": "supplied-graph-no-clique-expansion", "bases": {"0": [[20], [10]], "1": [[10, 20]], "2": []},
            "boundary": {"shape": [2, 1], "terms": [[1, 0, -1], [0, 0, 1]]},
            "features": {"shape": [2, 2], "row_keys": [[20], [10]], "values": [[5., 1.], [-3., 0.]],
                         "meaning": "supplied-coordinate", "units": "dimensionless"}}
    empty = copy.deepcopy(full)
    empty.update(selection={"kind": "stage", "index": 0, "scale": -2.}, coverage="inclusive-supplied-sublevel",
                 bases={"0": [], "1": [], "2": []}, boundary={"shape": [0, 0], "terms": []})
    empty["features"].update(shape=[0, 2], row_keys=[], values=[])
    stage = copy.deepcopy(full)
    stage.update(selection={"kind": "stage", "index": 1, "scale": 0.}, coverage="inclusive-supplied-sublevel")
    return fixture, [full, empty, stage]


class HandoffTests(unittest.TestCase):
    def setUp(self):
        self.fixture, self.records = hand_example()

    def test_hand_signed_edge_empty_shape_and_application_layout(self):
        learning.validate_records(self.records, self.fixture)
        full, _ = learning.route_inputs(list(reversed(self.records)), self.fixture, "topology")
        stages, _ = learning.route_inputs(self.records, self.fixture, "stage")
        self.assertEqual(full, [[1., 1., -3., 0., 5., 1., 1.]])
        self.assertEqual(stages, [[0.] * 7 + full[0]])

    def test_complete_basis_permutation_is_accepted(self):
        for r in (self.records[0], self.records[2]):
            r["bases"]["0"].reverse()
            r["features"]["row_keys"].reverse()
            r["features"]["values"].reverse()
            r["boundary"]["terms"] = [[1 - row, col, sign] for row, col, sign in r["boundary"]["terms"]]
        self.test_hand_signed_edge_empty_shape_and_application_layout()

    def test_missing_duplicate_foreign_record_or_label_fails(self):
        for mutate in (lambda r: r.pop(), lambda r: r.append(copy.deepcopy(r[0])),
                       lambda r: r[0].update(sample_key="foreign")):
            rows = copy.deepcopy(self.records)
            mutate(rows)
            with self.subTest(rows=rows), self.assertRaises(ValueError):
                learning.validate_records(rows, self.fixture)
        for labels in ([], self.fixture["labels"] * 2, [{"sample_key": "foreign", "label": 1, "split": "train"}]):
            self.fixture["labels"] = labels
            with self.assertRaises(ValueError):
                learning.validate_records(self.records, self.fixture)

    def test_source_shape_degree_sign_units_and_nonfinite_fail(self):
        mutations = [lambda r: r.update(snapshot="wrong"), lambda r: r.update(sample_order=2),
                     lambda r: r.update(domain="F2"), lambda r: r.update(scale_unit="conductance"),
                     lambda r: r["bases"].pop("2"), lambda r: r["bases"]["1"].append([20, 10]),
                     lambda r: r["boundary"].update(shape=[1, 2]),
                     lambda r: r["boundary"]["terms"][0].__setitem__(2, 1),
                     lambda r: r["boundary"]["terms"].append([1, 0, -1]),
                     lambda r: r["features"].update(shape=[2, 1]),
                     lambda r: r["features"]["values"].reverse(),
                     lambda r: r["features"]["row_keys"].reverse(),
                     lambda r: r["features"]["values"][0].__setitem__(0, float("nan"))]
        for mutate in mutations:
            rows = copy.deepcopy(self.records)
            mutate(rows[0])
            with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                learning.validate_records(rows, self.fixture)

    def test_wrong_stage_or_empty_zero_column_shape_fails(self):
        for mutate in (lambda r: r["selection"].update(scale=1.),
                       lambda r: r["boundary"].update(shape=[0, 1]),
                       lambda r: r["features"].update(shape=[0, 0])):
            rows = copy.deepcopy(self.records)
            mutate(rows[1])
            with self.assertRaises(ValueError):
                learning.validate_records(rows, self.fixture)

    def test_ambiguous_same_size_universe_is_rejected(self):
        fixture = copy.deepcopy(self.fixture)
        s = copy.deepcopy(fixture["samples"][0])
        s.update(sample_key="edge-b", sample_order=1, universe=[20, 10])
        fixture["samples"].append(s)
        fixture["labels"].append({"sample_key": "edge-b", "label": 0, "split": "test"})
        records = copy.deepcopy(self.records)
        for r in self.records:
            record = copy.deepcopy(r)
            record.update(sample_key="edge-b", sample_order=1)
            records.append(record)
        with self.assertRaisesRegex(ValueError, "universe"):
            learning.route_inputs(records, fixture, "topology")

    def test_fixture_has_disjoint_sample_split_and_real_empty_stages(self):
        fixture = learning.fixtures()
        self.assertEqual(len(fixture["samples"]), 32)
        self.assertEqual(sum(r["split"] == "train" for r in fixture["labels"]), 24)
        self.assertNotEqual([r["sample_key"] for r in fixture["labels"]],
                            [r["sample_key"] for r in fixture["samples"]])
        self.assertTrue(all(s["grid"][0] < min(v[1] for v in s["vertices"]) for s in fixture["samples"]))

    def test_process_failure_and_timeout_leave_logs_and_exit_status(self):
        with tempfile.TemporaryDirectory() as directory:
            for name, code, timeout in [("failure", "raise SystemExit(7)", 10),
                                        ("timeout", "import time; time.sleep(10)", .05)]:
                stem = Path(directory) / name
                with self.subTest(name=name), self.assertRaises(RuntimeError):
                    learning.run_process([sys.executable, "-c", code], stem, timeout=timeout)
                measurement = learning.read_json(stem.with_suffix(".process.json"))
                self.assertNotEqual(measurement["exit_code"], 0)
                self.assertEqual(measurement["timeout"], name == "timeout")
                self.assertTrue(stem.with_suffix(".stderr").exists())

    def test_nonfinite_json_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "invalid.json"
            path.write_text('{"probability":NaN}')
            with self.assertRaises(ValueError):
                learning.read_json(path)


HAS_PACKAGES = all(importlib.util.find_spec(m) is not None for m in ("numpy", "scipy", "sklearn", "gudhi"))
WORKER = os.environ.get("COCYCLE_LEARNING_WORKER")


@unittest.skipUnless(HAS_PACKAGES and WORKER, "optional packages and COCYCLE_LEARNING_WORKER required")
class RealProcessTests(unittest.TestCase):
    def test_both_routes_and_output_integrity(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "run"
            learning.run(Path(WORKER), output, 92, Path(sys.executable))
            fixture = learning.read_json(output / "fixture.json")
            records = [json.loads(line) for line in (output / "records.jsonl").read_text().splitlines()]
            result = learning.read_json(output / "predictions.json")
            self.assertEqual(len(records), 128)
            self.assertEqual(result["routes"]["topology"]["shape"], [32, 33])
            self.assertEqual(result["routes"]["stage"]["shape"], [32, 99])
            self.assertEqual(len(result["predictions"]), 64)
            for mutate in (lambda r: r["predictions"].pop(),
                           lambda r: r["predictions"].append(copy.deepcopy(r["predictions"][0])),
                           lambda r: r["predictions"][0].update(sample_key="foreign"),
                           lambda r: r["predictions"][0].update(snapshot="wrong"),
                           lambda r: r["predictions"][0].update(probability=float("nan")),
                           lambda r: r["predictions"][0].update(probability=.5),
                           lambda r: r["routes"]["stage"].update(accuracy=.5),
                           lambda r: r["routes"]["stage"]["mean"].__setitem__(10, .123),
                           lambda r: r["routes"]["stage"].update(shape=[32, 33]),
                           lambda r: r.update(input_sha256="wrong")):
                corrupted = copy.deepcopy(result)
                mutate(corrupted)
                with self.subTest(mutate=mutate), self.assertRaises(ValueError):
                    learning.validate_predictions(corrupted, records, fixture,
                                                  output / "records.jsonl", output / "fixture.json")
            # Train-only normalization: verify independently from fixture rows.
            values, _ = learning.route_inputs(records, fixture, "topology")
            self.assertAlmostEqual(result["routes"]["topology"]["mean"][6], sum(row[6] for row in values[:24]) / 24, places=14)
            with self.assertRaises(FileExistsError):
                learning.run(Path(WORKER), output, 92, Path(sys.executable))

    def test_rust_export_rejects_invalid_inputs(self):
        with tempfile.TemporaryDirectory() as directory:
            sample = learning.fixtures()["samples"][:1]
            for name, mutate in [("duplicate", lambda s: s[0]["vertices"].append(s[0]["vertices"][0])),
                                  ("nonfinite", lambda s: s[0]["vertices"][0].__setitem__(2, float("nan"))),
                                  ("grid", lambda s: s[0].update(grid=[0., 0.])),
                                  ("closure", lambda s: s[0]["edges"].append([10, 999, 1.]))]:
                rows = copy.deepcopy(sample)
                mutate(rows)
                with self.subTest(name=name), self.assertRaises(RuntimeError):
                    learning.run_process([WORKER], Path(directory) / name, learning.input_text(rows))

    def test_failed_external_call_cannot_produce_success_report(self):
        with tempfile.TemporaryDirectory() as directory:
            output = Path(directory) / "failed"
            original = learning.run_process
            def fail_consumer(command, stem, stdin=None, timeout=120):
                if stem.name == "consumer":
                    command = [sys.executable, "-c", "raise SystemExit(7)"]
                return original(command, stem, stdin, timeout)
            with patch.object(learning, "run_process", side_effect=fail_consumer):
                with self.assertRaises(RuntimeError):
                    learning.run(Path(WORKER), output, 92, Path(sys.executable))
            report = learning.read_json(output / "report.json")
            self.assertEqual(report["status"], "failure")
            self.assertEqual(report["phase"], "external_fit_predict")
            self.assertTrue((output / "records.jsonl").exists())
            self.assertFalse((output / "predictions.json").exists())
            self.assertEqual(learning.read_json(output / "consumer.process.json")["exit_code"], 7)


if __name__ == "__main__":
    unittest.main()
