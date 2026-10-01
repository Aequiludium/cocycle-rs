"""Fail-closed lifecycle transformations, exact inputs and result contracts."""
import copy
import struct
import tempfile
import unittest
from pathlib import Path

import profile_distance_lifecycle as study


class LifecycleTests(unittest.TestCase):
    def test_changed_markers_fail(self):
        for source in ('missing', 'marker marker'):
            with self.assertRaises(ValueError):
                study.replace(source, 'marker', 'new')
        for path in (study.REPO / 'src/diagram_distances').rglob('*.rs'):
            name = path.relative_to(study.REPO / 'src/diagram_distances').as_posix()
            output = study.instrument(path.read_text(), name)
            if name == 'mod.rs':
                self.assertIn('essential != second.len() - second_finite', output)
                self.assertIn('check_context(first', output)
                self.assertIn('finite(&first, &second, (first_finite, second_finite), budget)?', output)
            if name == 'wasserstein.rs':
                self.assertIn('numeric::scale_for(first.1.max(second.1))', output)
                self.assertIn('numeric::normalize(&mut a, scale, budget)?', output)
            if name == 'wasserstein/sparse.rs':
                self.assertIn('scratch = Scratch::new(', output)
                self.assertNotIn('lifecycle::empty("heap', output)

    def test_six_workloads_preserve_orientation_and_full_matrix(self):
        a, b = [(0., 2.)], [(1., 3.)]
        for pattern in study.PATTERNS:
            diagrams, ops, fixed = study.workload(a, b, pattern, 1 if pattern == 'single' else 16)
            self.assertTrue(ops)
            self.assertEqual(diagrams[fixed], a)
            if pattern == 'one_to_many':
                self.assertTrue(all(x == fixed for x,y,_ in ops))
            if pattern == 'many_to_one':
                self.assertTrue(all(y == fixed for x,y,_ in ops))
            if pattern == 'all_pairs':
                self.assertEqual(len(ops),len(diagrams)**2)
                self.assertEqual(len(set((x,y) for x,y,_ in ops)),len(ops))
            if pattern == 'repeated_batch':
                self.assertEqual({batch for _,_,batch in ops},{0,1})
        with self.assertRaises(ValueError):
            study.workload(a,b,'single',4)

    def test_binary_preserves_float_bits_and_rejects_invalid_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            path=Path(tmp)/'input.bin'
            study.write_input(path,[[(.125,.375)]],[(0,0,0)],0)
            self.assertEqual(path.read_bytes(),b'COCLIF1\0'+struct.pack('<QQQQddQQQ',1,1,0,1,.125,.375,0,0,0))
            for diagrams,ops in [([[(0.,0.)]],[(0,0,0)]),([[(0.,float('inf'))]],[(0,0,0)]),([[(0.,1.)]],[(0,1,0)])]:
                with self.assertRaises(ValueError):
                    study.write_input(path,diagrams,ops,0)

    def test_result_identity_timing_and_phases_fail_closed(self):
        record={'type':'result','protocol':study.PROTOCOL,'metric':'w1','variant':'r2','operations':4,'total_ns':20,'prep_ns':5,'loop_ns':14,'prepared_bytes':64,'retained_bytes':0,'peak_retained_bytes':0,'output_bytes':32,'hits':0,'growths':0,'phases_ns':[0]*10,'checksum':1.}
        self.assertEqual(study.validate(record,'w1','r2',4),record)
        for key,value in [('metric','w2'),('operations',3),('total_ns',-1),('prep_ns',21),('phases_ns',[0]),('retained_bytes',True),('checksum',float('nan'))]:
            changed=copy.deepcopy(record);changed[key]=value
            with self.assertRaises(ValueError):
                study.validate(changed,'w1','r2',4)


if __name__ == '__main__':
    unittest.main()
