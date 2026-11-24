"""Small independent contracts for the optional HeCS experiment."""
import importlib.util
import itertools
import unittest

import reproduce_hecs as hecs


def exhaustive_collapsible(edges, faces, remaining):
    """Search ALL allowed first collapses directly from vertex containment."""
    if not remaining:
        return True
    for edge in edges:
        adjacent = [t for t in remaining if set(edge) <= set(faces[t])]
        if len(adjacent) == 1 and exhaustive_collapsible(edges, faces, remaining - set(adjacent)):
            return True
    return False


class CollapseTests(unittest.TestCase):
    def setUp(self):
        self.edges = list(itertools.combinations(range(4), 2))
        self.faces = list(itertools.combinations(range(4), 3))
        self.adjacency = [[j for j, t in enumerate(self.faces) if set(e) <= set(t)]
                          for e in self.edges]

    def test_all_shell_subsets_against_exhaustive_search(self):
        for bits in range(16):
            selected = [j for j in range(4) if bits & (1 << j)]
            pairs, remaining = hecs.collapse(self.adjacency, selected)
            self.assertEqual(not remaining, exhaustive_collapsible(self.edges, self.faces, set(selected)))
            self.assertTrue(hecs.replay(self.edges, self.faces, [t for _, t in pairs], pairs))

    def test_weighted_shell_and_complete_lower_skeleton(self):
        result = hecs.heavy(self.adjacency + [[]], [1, 2, 4, 8])
        self.assertEqual(result["selected"], [3, 2, 1])
        self.assertEqual(result["initial_pairs"], [])
        self.assertEqual(result["decisions"][-1]["remaining"], [0, 1, 2, 3])
        self.assertTrue(hecs.replay(self.edges, self.faces, result["selected"], result["pairs"]))

    def test_empty_and_ties(self):
        self.assertEqual(hecs.heavy([], [])["selected"], [])
        self.assertEqual(hecs.heavy(self.adjacency, [1] * 4)["selected"], [0, 1, 2])


@unittest.skipUnless(importlib.util.find_spec("numpy"), "optional NumPy experiment dependency")
class CglsTests(unittest.TestCase):
    def test_rectangular_and_zero_rhs(self):
        import numpy as np
        a = np.array([[1., 0.], [0., 2.], [0., 0.]])
        b = np.array([1., 4., 0.])
        result = hecs.cgls(a, b, a, b, np.eye(2))
        self.assertEqual(result["status"], "converged")
        np.testing.assert_allclose(result["x"], [1, 2], atol=1e-12)
        zero = hecs.cgls(a, b * 0, a, b * 0, np.eye(2))
        self.assertEqual(len(zero["trace"]), 1)

    def test_incompatible_rhs_and_iteration_limit(self):
        import numpy as np
        a = np.diag([1., 0.])
        b = np.array([0., 1.])
        self.assertEqual(hecs.cgls(a, b, a, b, np.eye(2))["status"], "breakdown")
        self.assertEqual(hecs.cgls(a, np.ones(2), a, np.ones(2), np.eye(2), maxit=0)["status"], "iteration_limit")


if __name__ == "__main__":
    unittest.main()
