"""Failed native measurements must never become performance rankings."""

import csv
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import benchmark_flag_h1 as study


class AdmissionTests(unittest.TestCase):
    def measure(self, outputs):
        case = study.Case('test', 2, [1.])
        with tempfile.TemporaryDirectory() as directory:
            fixture = Path(directory) / 'input.bin'
            case.write(fixture)
            with patch.object(study, 'worker', side_effect=outputs):
                record = study.measure(case, fixture, {n: [n] for n in study.BACKENDS},
                                       3, 1729, {}, 1., 64)
            csv_path = Path(directory) / 'summary.csv'
            study.summarize(csv_path, [record])
            with csv_path.open() as stream:
                rows = list(csv.DictReader(stream))
        return record, rows

    @staticmethod
    def completed(bars=None):
        return {'status': 'completed', 'coverage': ['complete', None],
                'intervals': bars or [[0, 0., 'F', 1.], [0, 0., 'E', 0.]],
                'elapsed_ms': 1., 'rss_before_kib': 10, 'hwm_before_kib': 10,
                'peak_rss_kib': 20}

    def test_complete_balanced_groups_are_ranked(self):
        record, rows = self.measure([self.completed() for _ in range(9)])
        self.assertEqual(record['validation'], 'passed')
        self.assertTrue(all(row['median_ms'] == '1.0' for row in rows))
        for backend in study.BACKENDS:
            positions = [s['position'] for s in record['results'][backend]['samples']]
            self.assertEqual(sorted(positions), [0, 1, 2])

    def test_failed_backend_withholds_all_rankings(self):
        record, rows = self.measure([{'status': 'timeout'}] +
                                    [self.completed() for _ in range(6)])
        self.assertEqual(record['validation'], 'failed')
        self.assertTrue(all(row['median_ms'] == '' for row in rows))

    def test_changed_multiplicity_withholds_all_rankings(self):
        record, rows = self.measure([self.completed()] * 3 +
                                    [self.completed([[0, 0., 'E', 0.]])] +
                                    [self.completed()] * 5)
        self.assertEqual(record['validation'], 'failed')
        self.assertTrue(all(row['median_ms'] == '' for row in rows))

    def test_consistently_wrong_backend_withholds_all_rankings(self):
        schedule = [entry for entry in study.schedule(study.BACKENDS, 3, 1729)
                    if not entry['warmup']]
        outputs = [self.completed([[0, 0., 'E', 0.]])
                   if entry['backend'] == 'candidate' else self.completed()
                   for entry in schedule]
        record, rows = self.measure(outputs)
        self.assertEqual(record['validation'], 'mismatch')
        self.assertTrue(all(row['median_ms'] == '' for row in rows))


if __name__ == '__main__':
    unittest.main()
