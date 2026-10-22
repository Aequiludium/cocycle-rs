"""Independent controller failures and fixture chain identities."""

import unittest
from check_homology_operators import compare, expected_span, fixture, wire
from benchmark_operator_family import PROTOCOL, fixture as family_fixture, oracle, validate_result


class OperatorProtocolTests(unittest.TestCase):
    def test_family_fixture_oracle_distinguishes_filled_disk_and_many_live_classes(self):
        disk = oracle(family_fixture("fan", 4))
        graph = oracle(family_fixture("graph", 4))
        self.assertEqual(disk[0][-1], 0)
        self.assertEqual(max(disk[0]), 2)
        self.assertEqual(sum(death is None for _, death in disk[1]), 0)
        self.assertEqual(graph[0][-1], 3)
        self.assertEqual(sum(death is None for _, death in graph[1]), 3)

    def test_family_resource_failures_and_wrong_betti_cannot_pass(self):
        cells = family_fixture("fan", 4)
        with self.assertRaises(ValueError):
            validate_result({"state": "ResourceExhausted"}, cells, "construct")
        result = {"state": "Computed", "protocol": PROTOCOL, "mode": "construct",
                  "stages": len(cells), "betti": [0] * len(cells), "bars": [], "histories": []}
        with self.assertRaises(ValueError):
            validate_result(result, cells, "construct")
    def test_generator_preserves_chain_equation_across_word_boundaries(self):
        for n in (0, 1, 5, 63, 64, 65, 129):
            data = fixture(n, 20261009)
            self.assertEqual(len(data["d"]), n)
            for row in data["a"]:
                for j in range(data["p"]):
                    self.assertEqual(sum(row[i] * data["d"][i][j] for i in range(n)) % 2, 0)
            self.assertEqual(wire(data), wire(fixture(n, 20261009)))

    def test_matching_failed_workers_cannot_pass(self):
        for state in ("ResourceExhausted", "Unavailable", "NotApplicable"):
            with self.assertRaises(RuntimeError):
                compare({"state": state}, {"state": state})

    def test_projection_and_mass_mismatches_veto_comparison(self):
        first = {"state": "Computed", "betti": 1, "projection": [[1]],
                 "actions": [[1]], "masses": [2]}
        for key, value in (("projection", [[0]]), ("masses", [3]), ("betti", 0)):
            with self.assertRaises(RuntimeError):
                compare(first, {**first, key: value})

    def test_span_oracle_rejects_homology_losing_zero_projection(self):
        data = {"n": 1, "p": 0, "a": [], "d": [[]], "vectors": [[0], [1]]}
        with self.assertRaises(RuntimeError):
            expected_span(data, {"betti": 1, "actions": [[0], [0]]})


if __name__ == "__main__":
    unittest.main()
