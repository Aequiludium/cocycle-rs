"""Failure controls for bounded research output acceptance."""
import unittest

import check_research_routes as checks


class ResearchRouteTests(unittest.TestCase):
    def records(self):
        updates = [{"kind": "feedback", "snapshot": f"edge-{i}", "iteration": i,
                    "theta": 1+2.**-i, "mass": 1+2.**-i, "finite_count": 1,
                    "loss": 0.5*4.**-i, "gradient": 2.**-i, "difference": 2.**-i,
                    "h": 1/1024} for i in range(9)]
        return updates + [{"kind": "tie", "theta": x, "mass": max(0., x),
                           "finite_count": int(x > 0), "loss": 0.5*(max(0., x)-1)**2}
                          for x in [-1/1024, 0., 1/1024]]

    def test_missing_or_reused_snapshot_cannot_pass(self):
        records = self.records()
        checks.check_feedback(records)
        with self.assertRaises(ValueError):
            checks.check_feedback(records[:-1])
        records[1]["snapshot"] = records[0]["snapshot"]
        with self.assertRaises(ValueError):
            checks.check_feedback(records)

    def test_false_tie_or_gradient_is_rejected(self):
        for field, index, value in [("finite_count", 10, 1), ("gradient", 0, 0.),
                                    ("difference", 2, 10.)]:
            records = self.records()
            records[index][field] = value
            with self.assertRaises(ValueError):
                checks.check_feedback(records)

    def test_empty_or_repeated_sheaf_output_cannot_pass(self):
        with self.assertRaises(ValueError):
            checks.check_sheaf([])
        with self.assertRaises(ValueError):
            checks.check_sheaf([{"name": "nonconstant"}]*6)


if __name__ == "__main__":
    unittest.main()
