"""Independent trace integration and fail-closed concurrency result contracts."""
import copy
import unittest
import benchmark_distance_resources as study


class ResourceTests(unittest.TestCase):
    def test_memory_integral_has_byte_second_units(self):
        points=[{'t_ns':0,'pss_bytes':2},{'t_ns':1000000000,'pss_bytes':4},{'t_ns':3000000000,'pss_bytes':4}]
        self.assertEqual(study.integrate(points,'pss_bytes'),11)
        with self.assertRaises(ValueError):study.integrate(list(reversed(points)),'pss_bytes')

    def test_invalid_worker_results_do_not_enter_comparison(self):
        good={'type':'resource-result','protocol':study.PROTOCOL,'metric':'w1','variant':'arena','threads':2,
              'jobs_per_thread':2,'durations_ns':[[1,2],[3,4]],'native_group_ns':5,'checksum':2.}
        self.assertEqual(study.validate(good,'w1','arena',2,2),good)
        for key,value in [('protocol','wrong'),('threads',1),('durations_ns',[[1,2]]),('durations_ns',[[1,-1],[3,4]]),('checksum',float('nan'))]:
            bad=copy.deepcopy(good);bad[key]=value
            with self.assertRaises(ValueError):study.validate(bad,'w1','arena',2,2)


if __name__=='__main__':unittest.main()
