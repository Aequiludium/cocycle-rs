"""Contract checks for native Dory numeric input/output capture."""
from pathlib import Path
import tempfile
import unittest

from benchmark_dory import csv_input, intervals
from benchmark_rips_pipeline import phase3_cases
from compare_rips import Fixture


class DoryTests(unittest.TestCase):
    def test_square_csv_preserves_explicit_vertex_count_and_duplicate_distance(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / 'input.csv'
            csv_input(Fixture('small', 'dense', 3, 2, None, [0., .5, 2.]), path)
            self.assertEqual([[float(v) for v in row.split(',')] for row in path.read_text().splitlines()],
                             [[0., 0., .5], [0., 0., 2.], [.5, 2., 0.]])

    def test_interval_capture_retains_multiplicity_and_unpaired_sentinel(self):
        with tempfile.TemporaryDirectory() as directory:
            prefix = Path(directory) / 'out-'
            Path(str(prefix) + 'H2_pers_data.txt').write_text('1, 2\n1, 2\n1, -1\n0, 0\n')
            self.assertEqual(intervals(prefix, 2), [[2, 1., 2.], [2, 1., 2.], [2, 1., None]])

    def test_phase3_workloads_are_reproducible_and_have_h1_controls(self):
        first, second = list(phase3_cases(True)), list(phase3_cases(True))
        self.assertEqual(first, second)
        self.assertEqual(len(first), 23)
        self.assertTrue(all(c.fixture.characteristic == 2 and not c.representatives for c in first))
        controls = [c for c in first if c.fixture.q == 1]
        self.assertEqual(len(controls), 10)
        for control in controls:
            parent = next(c for c in first if c.fixture.name == control.fixture.name.removesuffix('_h1'))
            self.assertEqual(control.fixture.data, parent.fixture.data)


if __name__ == '__main__':
    unittest.main()
