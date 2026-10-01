"""Independent trace integration and fail-closed concurrency result contracts."""
import copy
import unittest
import benchmark_distance_resources as study
import analyze_distance_resources as analysis


class ResourceTests(unittest.TestCase):
    def test_replay_preserves_residence_and_staggered_peak(self):
        # Two triangular 2-byte excursions one second apart overlap at one byte
        # each. Baseline residence is counted once and held after completion.
        trace=[(0,1),(1,3),(2,1)]
        value=analysis.replay([trace,trace],[0,1],3,10)
        self.assertEqual(value['peak'],12)
        self.assertEqual(value['integral'],34)
        with self.assertRaises(ValueError):analysis.replay([list(reversed(trace))],[0],3,10)

    def test_model_fit_rejects_holdout_and_limits_cpu_capacity(self):
        rows=[dict(type='case',split='tuning',workers=1,throughput=10,cpu=.08,service=.1),
              dict(type='case',split='tuning',workers=2,throughput=15,cpu=.1,service=2/15)]
        traces=[dict(type='case',split='tuning',workers=1,memory=.2,peak=3),
                dict(type='case',split='tuning',workers=2,memory=.3,peak=5)]
        model=analysis.fit(rows,traces,capacity=2,memory_mib=20)
        self.assertAlmostEqual(model['cpu_capacity'],1.5)
        self.assertAlmostEqual(analysis.predict(model,'case',2,1)['throughput'],18.75)
        bad=copy.deepcopy(rows);bad[1]['split']='holdout'
        with self.assertRaises(ValueError):analysis.fit(bad,traces)
        with self.assertRaises(ValueError):analysis.fit_interference(model,[dict(split='holdout')])

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
