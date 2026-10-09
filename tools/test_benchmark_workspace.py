"""Retained evidence and build reuse must fail closed without native execution."""
import copy
from contextlib import ExitStack
import json
import os
from pathlib import Path
import shutil
import subprocess
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import benchmark_workspace as workspace
from benchmark_rips_pipeline import Case
from compare_rips import Fixture


@unittest.skipUnless(workspace.sys.platform == 'linux', 'workspace controller requires Linux')
class WorkspaceControllerTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='workspace-controller-test-')
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.repo = self.root / 'repo'
        self.write_source(self.repo)
        for name in workspace.HARNESS_FILES:
            path = self.repo / name
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_text(name)
        native = self.root / 'native'
        (native / 'build').mkdir(parents=True)
        for name in ('gudhi', 'ripser'):
            (native / 'build' / name).write_bytes(b'not executed')
        self.native_environment = native / 'environment.json'
        self.native_environment.write_text(json.dumps({
            'pins': workspace.PINS,
            'binaries_sha256': {name: workspace.sha256(native / 'build' / name)
                                for name in ('gudhi', 'ripser')}}))
        self.identity = {'kernel_commit': 'f' * 40, 'harness_commit': 'f' * 40, 'dirty': False}
        self.case = Case(Fixture('pair', 'dense', 2, 0, None, [1.]), 'dense')
        self.sample = {'status': 'completed', 'elapsed_ms': 1., 'peak_rss_kib': 100,
                       'intervals': [[0, 0., 1.], [0, 0., None]]}
        self.trace = ('workspace_event={"vec_capacity_bytes":{"working":32}}\n'
                      'workspace_intervals=[(0, 0.0, Some(1.0)), (0, 0.0, None)]\n')

    @staticmethod
    def write_source(root):
        (root / 'src').mkdir(parents=True)
        (root / 'Cargo.toml').write_text('manifest')
        (root / 'Cargo.lock').write_text('lock')
        (root / 'src/lib.rs').write_text('source')

    def build_variant(self, name, output, local, log):
        root = output / name
        self.write_source(root)
        item = {'source_sha256': workspace.fingerprint(root), 'commands': []}
        for kind, filename in (('binary', name), ('trace', name + '-trace')):
            path = local / filename
            path.write_bytes(b'not executed')
            item[kind] = str(path)
            item[kind + '_sha256'] = workspace.sha256(path)
        return item

    def invoke(self, label, *, reuse=None, trace=None):
        output = self.root / (self.root.name + '-' + label)
        local = Path('/tmp') / ('cocycle-workspace-' + output.name + '-' + 'f' * 12)
        self.addCleanup(shutil.rmtree, local, True)
        args = SimpleNamespace(output=output, native_environment=self.native_environment,
                               samples=10, cpu=min(os.sched_getaffinity(0)), order_seed=1,
                               reuse_build=reuse)
        with ExitStack() as stack:
            stack.enter_context(patch.object(workspace, 'REPO', self.repo))
            # Isolate the hash guard from the existing clean-worktree/commit guard.
            stack.enter_context(patch.object(workspace, 'provenance', return_value=self.identity))
            stack.enter_context(patch.object(workspace, 'workspace_cases', return_value=[self.case]))
            stack.enter_context(patch.object(workspace, 'build_variant', side_effect=self.build_variant))
            stack.enter_context(patch.object(workspace.subprocess, 'check_output', return_value='rustc test'))
            if isinstance(trace, Exception):
                stack.enter_context(patch.object(workspace.subprocess, 'run', side_effect=trace))
            else:
                stack.enter_context(patch.object(workspace.subprocess, 'run', return_value=
                    SimpleNamespace(stdout=self.trace if trace is None else trace,
                                    stderr='', returncode=0)))
            measured = stack.enter_context(patch.object(workspace, 'worker', side_effect=
                lambda *args: copy.deepcopy(self.sample)))
            try:
                workspace.run(args)
            except ValueError as error:
                return output, error, measured.call_count
        return output, None, measured.call_count

    def test_failed_diagnostics_retain_samples_and_failure_summary(self):
        failures = ('', self.trace.replace('Some(1.0)', 'Some(2.0)'),
                    subprocess.TimeoutExpired('diagnostic', 30))
        for index, trace in enumerate(failures):
            with self.subTest(trace=repr(trace)):
                output, error, calls = self.invoke('trace-' + str(index), trace=trace)
                self.assertIsInstance(error, ValueError)
                result, = json.loads((output / 'results.json').read_text())
                summary = json.loads((output / 'summary.json').read_text())
                self.assertEqual(summary['status'], 'failed')
                self.assertEqual(summary['processes'], calls)
                self.assertEqual(len(result['samples']), calls)
                self.assertTrue(all(row['status'] == 'completed' for row in result['samples']))
                self.assertTrue(all('comparison_error' in row for row in result['traces'].values()))
                for row in result['summary'].values():
                    self.assertEqual(row['samples'], 10)
                    for key in ('median_ms', 'min_ms', 'max_ms', 'max_peak_rss_kib'):
                        self.assertIsNone(row[key])

    def test_reuse_accepts_matching_identity_and_rejects_changed_source_or_harness(self):
        previous, error, _ = self.invoke('initial')
        self.assertIsNone(error)
        _, error, calls = self.invoke('matching', reuse=previous)
        self.assertIsNone(error)
        self.assertGreater(calls, 0)
        paths = ['src/lib.rs', *workspace.HARNESS_FILES]
        for index, name in enumerate(paths):
            with self.subTest(path=name):
                path = self.repo / name
                original = path.read_bytes()
                path.write_bytes(original + b'changed')
                _, error, calls = self.invoke('changed-' + str(index), reuse=previous)
                path.write_bytes(original)
                self.assertRegex(str(error), 'same frozen source, harness')
                self.assertEqual(calls, 0)
        environment_path = previous / 'environment.json'
        environment = json.loads(environment_path.read_text())
        del environment['harness_files_sha256']
        environment_path.write_text(json.dumps(environment))
        _, error, calls = self.invoke('missing-identity', reuse=previous)
        self.assertRegex(str(error), 'same frozen source, harness')
        self.assertEqual(calls, 0)


if __name__ == '__main__':
    unittest.main()
