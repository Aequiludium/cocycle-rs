"""Contracts that keep native critical-set measurements semantically comparable."""
import json
from pathlib import Path
import tempfile
import subprocess
import unittest
from unittest.mock import patch

from benchmark_critical_sets import canonical, dense_source, diagram, fixture_text, fixtures, pair_oracle, validate, workloads


class CriticalSetBenchmarkTests(unittest.TestCase):
    def test_fresh_clone_uses_the_pinned_benchmark_only_excerpt(self):
        with patch('benchmark_critical_sets.git', side_effect=subprocess.CalledProcessError(128, 'git')):
            data = dense_source()
        self.assertIn(b'fn reduce(d: &Matrix, rows: usize)', data)
        self.assertIn(b'fn critical(', data)

    def test_dense_source_rejects_a_changed_baseline(self):
        with tempfile.TemporaryDirectory() as directory:
            path=Path(directory)/'big_steps.rs'
            path.write_bytes(b'fn main() {}\n')
            with self.assertRaisesRegex(ValueError,'pinned Git blob'):
                dense_source(path)

    def test_pairing_oracle_has_hand_derived_delayed_triangle(self):
        cells=canonical([([0],0.),([1],0.),([2],0.),([0,1],1.),([0,2],2.),([1,2],3.),([0,1,2],4.)])
        pairs,essential=pair_oracle(cells)
        self.assertEqual(diagram(cells,pairs,essential,1),[[0,0.,1.],[0,0.,2.],[0,0.,None],[1,3.,4.]])

    def test_shared_fixtures_are_closed_and_only_flag_compatible_rows_admit_flag(self):
        for case in fixtures(quick=True):
            cells=case['cells']; lookup={tuple(v):f for v,f in cells}
            pairs,essential=pair_oracle(cells)
            self.assertTrue(pairs)
            self.assertTrue(essential)
            for v,f in cells:
                if len(v)>1:
                    self.assertTrue(all(lookup[tuple(v[:i]+v[i+1:])]<=f for i in range(len(v))))
                if case['flag'] and len(v)==3:
                    expected=max(lookup[(v[0],v[1])],lookup[(v[0],v[2])],lookup[(v[1],v[2])])
                    self.assertEqual(f,expected)
            for _,proposals,repeats in workloads(case):
                text=fixture_text(case,proposals,repeats)
                self.assertEqual(int(text.split()[0]),len(cells))
                finite={id for pair in pairs for id in pair}
                self.assertTrue(all(id in finite for id,_ in proposals))

    def test_validation_rejects_missing_multiplicity_wrong_targets_and_pairings(self):
        expected=[[0,0.,1.],[0,0.,1.],[0,0.,None]]
        record={'intervals':expected,'targets':[[2,1.5]],'pairs':[[0,2]],'essential':[1]}
        validate(json.loads(json.dumps(record)),expected,[[2,1.5]],[(0,2)],[1])
        for update in ({'intervals':expected[1:]},{'targets':[[2,1.6]]},{'pairs':[[1,2]]},{'essential':[]}):
            invalid=dict(record,**update)
            with self.assertRaises(ValueError):
                validate(invalid,expected,[[2,1.5]],[(0,2)],[1])


if __name__=='__main__':
    unittest.main()
