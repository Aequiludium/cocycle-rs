"""Failure controls for route identity, converted results and independent expectations."""
import copy
import unittest

from check_comparison_routes import check


class ComparisonRouteTests(unittest.TestCase):
    def rows(self):
        common = dict(normalization='none', units='metres', alignment='lab-cartesian-v1')
        return [dict(common, route='geometry', ids=['a', 'b', 'c'],
                     snapshots=['measurements/a/v1', 'measurements/b/v1', 'measurements/c/v1'],
                     dimension=2, points=[[0, 0, 2, 0], [0, 1, 2, 1], [0, 3, 2, 3]],
                     values=[1, 3, 2], nearest='b', converted_alias=True, numeric_bytes=24),
                dict(common, route='results', ids=['rips-a', 'rips-b'],
                     snapshots=['rips/a/v1', 'rips/b/v1'], field=2, degree=0,
                     first_points=[0, 2], second_points=[0, 3],
                     first_diagram=[[0, 2], [0, None]], second_diagram=[[0, 3], [0, None]],
                     distances=[1, 1, 1], converted=[1, 1, 1], within_one_metre=True),
                dict(route='filtrations', ids=['filtration-a', 'filtration-b'],
                     snapshots=['triangle/a/v1', 'triangle/b/v1'], units='seconds',
                     alignment='triangle-labels-v1', normalization='none', field=2, degree=1,
                     first_cells=[[[10], 0], [[20], 0], [[30], 0], [[10, 20], 1],
                                  [[10, 30], 2], [[20, 30], 3], [[10, 20, 30], 6]],
                     second_cells=[[[10], 0], [[20], 0], [[30], 0], [[10, 20], 3],
                                   [[10, 30], 1], [[20, 30], 2], [[10, 20, 30], 7]],
                     first_diagram=[[3, 6]], second_diagram=[[3, 7]], sup=2,
                     distances=[1, 1, 1], bound_holds=True)]

    def test_complete_three_routes_pass(self):
        self.assertEqual(check(self.rows())['exact_matching_checks'], 6)

    def test_missing_route_or_operand_declaration_fails(self):
        for edit in (lambda r: r.pop(), lambda r: r[0]['ids'].reverse(),
                     lambda r: r[1].update(units='centimetres'),
                     lambda r: r[2].update(alignment=None)):
            rows = self.rows()
            edit(rows)
            with self.assertRaises(ValueError):
                check(rows)

    def test_bad_result_or_conversion_or_consumer_fails(self):
        original = self.rows()
        for route, field, value in [(1, 'distances', [0, 1, 1]),
                                    (1, 'converted', [0, 0, 0]),
                                    (1, 'within_one_metre', False),
                                    (0, 'nearest', 'c'), (2, 'sup', 0)]:
            rows = copy.deepcopy(original)
            rows[route][field] = value
            with self.assertRaises(ValueError):
                check(rows)

    def test_wrong_domain_and_duplicate_cells_fail(self):
        rows = self.rows()
        rows[2]['second_cells'][0][0] = [42]
        with self.assertRaises(ValueError):
            check(rows)
        rows = self.rows()
        rows[2]['second_cells'].append(rows[2]['second_cells'][0])
        with self.assertRaises(ValueError):
            check(rows)


if __name__ == '__main__':
    unittest.main()
