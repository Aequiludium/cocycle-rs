"""Independent small controls, literal-paper failures and experiment contracts."""
import importlib.util
import itertools
from pathlib import Path
import shutil
import subprocess
import tempfile
import unittest

AVAILABLE = all(importlib.util.find_spec(name) for name in ["numpy", "scipy"])
if AVAILABLE:
    import numpy as np
    import reproduce_nonbranching as replay


@unittest.skipUnless(AVAILABLE, "optional NumPy/SciPy reproduction dependencies")
class WeakReductionTests(unittest.TestCase):
    def test_literal_sign_and_row_order_failures_are_retained(self):
        cases = {case["name"]: case for case in replay.fixtures()}
        same = cases["double-sign"]
        literal = replay.weak_reduction(same["B1"], same["D"], literal=True)
        self.assertFalse(literal["identity"])
        star = cases["r-reorientation"]
        literal = replay.weak_reduction(star["B1"], star["D"], literal=True)
        self.assertTrue(literal["kernel"])
        self.assertFalse(literal["identity"])
        contradiction = cases["signed-contradiction"]
        literal = replay.weak_reduction(contradiction["B1"], contradiction["D"], literal=True)
        self.assertFalse(literal["kernel"])
        row_order = cases["row-order"]
        ok, _ = replay.literal_reorientation(row_order["D"], mutate=False)
        self.assertFalse(ok)  # removing only the D mutation is insufficient

    def test_exhaustive_two_row_three_column_exact_oracle(self):
        rows = [row for row in itertools.product([-1, 0, 1], repeat=3)
                if sum(x != 0 for x in row) <= 2]
        for first, second in itertools.product(rows, repeat=2):
            d = np.array([first, second])
            case = dict(name=str((first, second)), B1=np.eye(3, dtype=int), D=d)
            replay.verify(case, replay.weak_reduction(case["B1"], d))

    def test_metric_and_original_zeros_contribute_once(self):
        cases = {case["name"]: case for case in replay.fixtures()}
        path = cases["path"]
        result = replay.weak_reduction(path["B1"], path["D"])
        self.assertEqual(result["metric"].tolist(), [[2]])
        self.assertEqual(float(result["upper_exact"][0][0]), 0.5)
        cube = cases["paper-6.1"]
        result = replay.weak_reduction(cube["B1"], cube["D"])
        self.assertEqual(result["original_zero"], [2, 3])
        self.assertEqual(np.diag(result["metric"]).tolist(), [1, 1, 2])
        self.assertEqual(result["Z"].shape[1], 3)
        # Independent transcription of the factor displayed on PDF page 22.
        paper_factor = [[-1, 0, 0], [-1, 0, 0], [-1, 0, 0], [1, 0, 0],
                        [1, -1, 0], [1, 0, -1], [0, -1, 0], [0, 1, -1],
                        [0, 0, 1], [0, 1, 0], [0, 0, 1]]
        self.assertEqual(result["factor_integer"][:, [2, 0, 1]].tolist(), paper_factor)

    def test_mobius_is_irregular_over_reals(self):
        case = next(c for c in replay.fixtures() if c["name"] == "paper-3.9")
        result = replay.weak_reduction(case["B1"], case["D"])
        self.assertEqual(result["components"][0]["kind"], "irregular")
        self.assertEqual(result["Z"].shape, (5, 0))
        self.assertFalse(np.any(np.array(result["upper_exact"], dtype=float)))

    def test_shared_triangle_and_shell_multiplicities(self):
        cases = {case["name"]: case for case in replay.fixtures()}
        triangle = cases["triangle"]
        upper = np.array(replay.weak_reduction(triangle["B1"], triangle["D"])["upper_exact"], dtype=float)
        lower = np.array([[2, 1, -1], [1, 2, 1], [-1, 1, 2]])
        np.testing.assert_array_equal(upper + lower, 3 * np.eye(3))
        shell = cases["shell"]["B1"]
        signs = np.array([1, -1, 1, -1])
        np.testing.assert_array_equal(shell.T @ shell, 4 * np.eye(4) - np.outer(signs, signs))
        np.testing.assert_allclose(np.linalg.eigvalsh(shell @ shell.T), [0, 0, 0, 4, 4, 4], atol=1e-12)

    def test_rejects_branching_and_nonunit_data(self):
        with self.assertRaisesRegex(ValueError, "branching"):
            replay.weak_reduction(np.eye(3, dtype=int), np.ones((1, 3), dtype=int))
        with self.assertRaisesRegex(ValueError, "unit"):
            replay.weak_reduction(np.eye(1, dtype=int), np.array([[2]]))

    def test_all_frozen_cases_and_orientation_covariance(self):
        cases = replay.fixtures()
        for case in cases:
            replay.verify(case, replay.weak_reduction(case["B1"], case["D"]))
        by_name = {c["name"]: c for c in cases}
        for case in cases:
            if case["name"].endswith("-permuted-oriented"):
                original = by_name[case["name"].removesuffix("-permuted-oriented")]
                first = replay.weak_reduction(case["B1"], case["D"])
                second = replay.weak_reduction(original["B1"], original["D"])
                self.assertEqual(first["upper_exact"], second["upper_exact"])

    @unittest.skipUnless(importlib.util.find_spec("gudhi"), "optional GUDHI topology comparator")
    def test_gudhi_is_distinct_from_real_spectral_evidence(self):
        rows = replay.topology_checks()
        mobius = {row["field"]: row for row in rows if row["name"] == "mobius"}
        self.assertEqual(mobius[2]["persistent_rank"], 0)
        self.assertEqual(mobius[11]["persistent_rank"], 1)
        self.assertFalse(any(row["spectral_oracle"] for row in rows))

    @unittest.skipUnless(shutil.which("rustc"), "Rust standalone compiler")
    def test_rust_parity_and_rejection(self):
        with tempfile.TemporaryDirectory() as directory:
            worker = Path(directory) / "nonbranching-test.exe"
            subprocess.run(["rustc", "--edition=2024", "-Dwarnings",
                            str(replay.ROOT / "tools/reference/nonbranching.rs"), "-o", str(worker)],
                           check=True, capture_output=True)
            for case in replay.fixtures():
                expected = replay.weak_reduction(case["B1"], case["D"])
                for mode in ["baseline", "compact"]:
                    result = replay.worker_run(worker, case, mode)
                    np.testing.assert_allclose(result["upper"],
                                               np.array(expected["upper_exact"], dtype=float).ravel(), atol=1e-10)
                    self.assertEqual(result["components"], expected["components"])
                    if mode == "baseline":
                        for key in ["R", "E", "V", "Z"]:
                            self.assertEqual(result[key], expected[key].ravel().tolist())
            for mode in ["baseline", "compact"]:
                proc = subprocess.run([str(worker), mode, "1"], input="1 1 3 1 1 1 1 1 1",
                                      text=True, capture_output=True)
                self.assertNotEqual(proc.returncode, 0)
                self.assertIn("branching", proc.stderr)


if __name__ == "__main__":
    unittest.main()
