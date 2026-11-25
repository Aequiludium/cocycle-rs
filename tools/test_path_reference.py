"""Independent known-value and deliberate-corruption checks for path tooling."""
import contextlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import benchmark_path_topology as benchmark

HAS_SCIENCE = all(importlib.util.find_spec(p) for p in ("numpy", "scipy", "sympy"))


class PathBenchmarkTests(unittest.TestCase):
    def test_generator_counts_directed_returning_two_paths_and_frozen_edges(self):
        cases = list(benchmark.fixtures(True, [20261125]))
        self.assertEqual(len(cases), 8)
        for case in cases:
            edges = set(map(tuple, case["edges"]))
            expected = sum((a,b) in edges and (b,c) in edges
                           for a in range(case["vertices"]) for b in range(case["vertices"])
                           for c in range(case["vertices"]))
            self.assertEqual(case["two_paths"], expected)
        for family in ("cycle", "sparse", "dense", "multisquare"):
            pair = [c for c in cases if c["family"] == family]
            self.assertEqual(pair[0]["edges"], pair[1]["edges"])

    def test_paired_ratio_keeps_known_proportional_samples_exact(self):
        samples = {"generic": list(range(1,11)), "direct": list(range(2,21,2))}
        summary = benchmark.summarize(samples, 1)
        self.assertEqual(summary["generic_over_direct"], 0.5)
        self.assertEqual(summary["paired_bootstrap_95_percentile"], [0.5,0.5])
        self.assertEqual(summary["classification"], "direct_regression")


@unittest.skipUnless(HAS_SCIENCE, "optional pinned path numerical oracle dependencies")
class PathNumericalTests(unittest.TestCase):
    def test_exact_diamond_domain_uses_the_induced_metric(self):
        import numpy as np
        import check_path_real as oracle
        _, boundary, _, stages = oracle.prepare(4, [(0,1),(0,2),(1,3),(2,3)], [0]*4, 1, True)
        exact, numerical, rho, dimension = oracle.reference(boundary, stages[0], stages[0], 1)
        np.testing.assert_allclose(np.linalg.eigvalsh(exact), [2,2,2,4], atol=1e-12)
        np.testing.assert_allclose(exact, numerical, atol=1e-12)
        self.assertEqual((rho, dimension), (0,1))

    def test_verifier_rejects_an_injected_wrong_action(self):
        import check_path_real as oracle
        with tempfile.TemporaryDirectory() as temporary:
            root = Path(temporary)
            worker = root/"worker"
            worker.write_bytes(b"deliberate mock worker")
            corrupted = {"protocol":"cocycle-path-real-certificate-v1", "field":"real", "paths":[[0],[1],[2]], "dimension":0, "regular":True,
                         "source_time":0, "target_time":0, "basis":[[1,0,0],[0,1,0],[0,0,1]],
                         "laplacian_columns":[[1,0,0],[0,0,0],[0,0,0]]}
            with patch.object(oracle, "check", return_value=corrupted), patch("sys.argv", ["check_path_real", "--worker",str(worker),"--output",str(root/"result"),"--quick"]), contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(oracle.main(), 1)
            summary = json.loads((root/"result/summary.json").read_text())
            self.assertEqual(summary["status"], "failed")
            self.assertIn("numerical verification failed", summary["error"])


if __name__ == "__main__":
    unittest.main()
