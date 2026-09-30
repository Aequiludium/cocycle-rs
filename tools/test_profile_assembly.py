"""Assembly diagnostics must fail closed and preserve analytic multiplicities."""

import copy
import json
import unittest

import profile_assembly as assembly


class AssemblyTests(unittest.TestCase):
    def setUp(self):
        self.sample = {'intervals': [[0, 0., 1.], [0, 0., None]], 'phases_ms': [0., 0., 0., .01, 0.]}
        self.record = {'route': 'h0_union_find', 'kernel_ns': 4000, 'kernel_calls': 1,
                       'assembly_ns': 2000, 'assembly_calls': 1, 'raw_calls': 1,
                       'raw_count': 3, 'raw_capacity': 4, 'final_count': 2, 'final_capacity': 4,
                       'raw_element_bytes': 32, 'final_element_bytes': 32, 'materialization_ns': 500,
                       'diagram_validation_ns': 300, 'sorting_ns': 1000, 'capacity_growths': 1}

    def parse(self, record):
        return assembly.profile(assembly.PREFIX + json.dumps(record), self.sample)

    def test_distinct_denominators_and_overlapping_capacity(self):
        result = self.parse(self.record)
        self.assertEqual(result['assembly_fraction'], .2)
        self.assertEqual(result['kernel_assembly_fraction'], 1 / 3)
        self.assertEqual(result['public_remainder_ns'], 4000)
        self.assertEqual(result['coexisting_output_capacity_bytes'], 256)

    def test_missing_duplicate_and_malformed_records_fail(self):
        for stderr in ('', assembly.PREFIX + '[]', assembly.PREFIX + '{}',
                       (assembly.PREFIX + json.dumps(self.record) + '\n') * 2):
            with self.subTest(stderr=stderr), self.assertRaises(ValueError):
                assembly.profile(stderr, self.sample)
        for key, value in [('route', 'unknown'), ('kernel_calls', 2), ('assembly_calls', 0),
                           ('raw_calls', 0), ('raw_count', 1), ('final_count', 3),
                           ('raw_capacity', 2), ('final_capacity', 1), ('kernel_ns', 10000),
                           ('sorting_ns', 2000), ('capacity_growths', 3),
                           ('raw_element_bytes', 0), ('assembly_ns', True), ('kernel_ns', -1)]:
            changed = dict(self.record, **{key: value})
            with self.subTest(key=key, value=value), self.assertRaises(ValueError):
                self.parse(changed)
        missing = copy.copy(self.record)
        del missing['kernel_ns']
        with self.assertRaises(ValueError):
            self.parse(missing)

    def test_analytic_triangle_free_and_filled_triangle_outputs(self):
        cases = {c.fixture.name: c for c in assembly.cases(True)}
        bipartite = assembly.expected(cases['bipartite4_h1'])
        self.assertEqual(bipartite.count([1, 1., None]), 9)
        self.assertEqual(bipartite.count([0, 0., 1.]), 7)
        self.assertEqual(bipartite.count([0, 0., None]), 1)
        filled = assembly.expected(cases['filled_triangles4'])
        self.assertEqual(filled.count([0, 1., 2.]), 8)
        self.assertEqual(filled.count([0, 1., None]), 4)
        self.assertEqual(filled.count([1, 2., 3.]), 4)
        self.assertEqual(len(assembly.expected(cases['forest32_h0'])), 32)
        self.assertEqual(assembly.expected(cases['isolates32_h0']), [[0, 0., None]] * 32)

    def test_changed_markers_fail_and_decomposition_is_conditional(self):
        for source in ('missing', 'marker marker'):
            with self.assertRaises(ValueError):
                assembly.replace_once(source, 'marker', 'replacement')
        source = (assembly.REPO / 'src/persistence/mod.rs').read_text()
        coarse = assembly.instrument(source, 'persistence/mod.rs')
        detailed = assembly.instrument(source, 'persistence/mod.rs', True)
        self.assertNotIn('fn assembly_diagram_phases', coarse)
        self.assertNotIn('_previous_capacity', coarse)
        self.assertIn('_previous_capacity', detailed)
        self.assertIn('PersistenceInterval::new(dimension, birth, end)?', detailed)
        diagram = (assembly.REPO / 'src/diagram/persistence_diagram.rs').read_text()
        self.assertEqual(assembly.instrument(diagram, 'diagram/persistence_diagram.rs'), diagram)
        self.assertEqual(assembly.instrument(diagram, 'diagram/persistence_diagram.rs', True)
                         .count('intervals.sort_unstable_by(compare_intervals);'), 1)


if __name__ == '__main__':
    unittest.main()
