"""M1 lifetime/scratch ablations; pipeline v2 latency and separate capacity traces.

Generated variants change only explicit release or same-type heap reset sites.
They never enter the ordinary library. Run fresh serial processes on Linux.
"""
import argparse
import ast
from dataclasses import replace
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys

from benchmark_rips_pipeline import (REPO, compare_samples, phase3_cases,
                                     canonical, provenance, schedule, worker)
from build_native import PINS, sha256


def replace_once(source, old, new):
    if source.count(old) != 1:
        raise ValueError('workspace ablation no longer matches its frozen source')
    return source.replace(old, new)


def lifetime_source(source):
    marker = '        let next_level = next_dimension(access, &level, &mut || budget.step())?;'
    return replace_once(source, marker,
                        '        drop(columns);\n'
                        '        #[cfg(test)]\n'
                        '        let columns: Vec<Column<usize>> = Vec::new();\n' + marker)


def scratch_source(source):
    # All owner payloads stay owned; only two existing, identically typed heaps
    # survive column boundaries. Both logical contents are reset before use.
    declaration = ('        // Fresh scratch per column is deliberate: M1 measures reuse separately.\n'
                   '        let mut working: BinaryHeap<Reverse<TupleEntry<4>>> = BinaryHeap::new();\n'
                   '        let mut transform = BinaryHeap::new();')
    source = replace_once(source, declaration,
                          '        working.clear();\n        transform.clear();')
    marker = '    for j in (0..triangles.len()).rev() {'
    source = replace_once(source, marker,
                          '    let mut working: BinaryHeap<Reverse<TupleEntry<4>>> = BinaryHeap::new();\n'
                          '    let mut transform = BinaryHeap::new();\n' + marker)
    end = ('        "h2_reduction_end",\n        &raw,\n        &triangles,\n'
           '        &cleared,\n        &owners,\n        &columns,\n        (0, 0),')
    return replace_once(source, end, end.rsplit('(0, 0),', 1)[0]
                        + '(working.capacity() * std::mem::size_of::<Reverse<TupleEntry<4>>>(),\n'
                        + '         transform.capacity() * std::mem::size_of::<usize>()),')


def workspace_cases():
    seen = set()
    for quick in (True, False):
        for case in phase3_cases(quick):
            if case.name not in seen:
                seen.add(case.name)
                yield case
    # Existing generic continuation's old/new-level overlap; no T9 reuse claim.
    for case in phase3_cases(True):
        if case.fixture.q == 2 and any(name in case.name for name in
                                      ('uniform', 'low_degree', 'high_fill', 'equal_clique')):
            yield replace(case, fixture=replace(case.fixture, name=case.fixture.name + '_h3', q=3))


def fingerprint(root):
    digest = hashlib.sha256()
    for path in [root / 'Cargo.toml', root / 'Cargo.lock', *sorted((root / 'src').rglob('*.rs'))]:
        digest.update(path.relative_to(root).as_posix().encode() + b'\0' + path.read_bytes())
    return digest.hexdigest()


def build_variant(name, output, local, log):
    root = output / name
    root.mkdir()
    shutil.copytree(REPO / 'src', root / 'src')
    for filename in ('Cargo.toml', 'Cargo.lock'):
        shutil.copy2(REPO / filename, root / filename)
    source = root / 'src/persistence/flag/cohomology/mod.rs'
    generic = root / 'src/persistence/simplicial/cohomology.rs'
    if name == 'scratch':
        source.write_text(scratch_source(source.read_text()), encoding='utf-8')
    elif name == 'lifetime':
        generic.write_text(lifetime_source(generic.read_text()), encoding='utf-8')
    identity = fingerprint(root)
    env = dict(os.environ, RUSTFLAGS='--cfg cocycle_h2_bench')
    commands = []

    def run(command):
        commands.append(command)
        subprocess.run(command, cwd=root, env=env, stdout=log, stderr=subprocess.STDOUT, check=True)

    run(['cargo', 'build', '--locked', '--offline', '--release', '--lib'])
    binary = local / name
    run(['rustc', '--edition=2024', '-O', '-D', 'warnings',
         str(REPO / 'benches/pipeline/cocycle.rs'), '--extern',
         f'cocycle={root}/target/release/libcocycle.rlib', '-L',
         f'dependency={root}/target/release/deps', '-o', str(binary)])
    # Test-only capacity accounting is compiled separately from latency.
    run(['cargo', 'test', '--locked', '--offline', '--release', '--lib', '--no-run',
         '--message-format=json'])
    run(['cargo', 'test', '--locked', '--offline', '--release', '--lib'])
    tests = list((root / 'target/release/deps').glob('cocycle-*'))
    tests = [path for path in tests if path.is_file() and os.access(path, os.X_OK)]
    if len(tests) != 1:
        raise ValueError('expected exactly one diagnostic test executable')
    diagnostic = local / (name + '-trace')
    shutil.copy2(tests[0], diagnostic)
    if fingerprint(root) != identity:
        raise ValueError('generated kernel changed during build')
    return {'source_sha256': identity, 'binary_sha256': sha256(binary),
            'trace_sha256': sha256(diagnostic), 'commands': commands,
            'binary': str(binary), 'trace': str(diagnostic)}


def summarize(samples, names):
    result = {}
    for name in names:
        rows = [sample for sample in samples if sample['backend'] == name and not sample['warmup']]
        valid = all(sample['status'] == 'completed' and 'comparison_error' not in sample for sample in samples)
        times = [sample['elapsed_ms'] for sample in rows] if valid else []
        result[name] = {'samples': len(rows), 'median_ms': statistics.median(times) if times else None,
                        'min_ms': min(times, default=None), 'max_ms': max(times, default=None),
                        'max_peak_rss_kib': max((sample['peak_rss_kib'] for sample in rows), default=None)
                        if valid else None}
    return result


def trace_result(stdout, returncode, anchor):
    result = {'returncode': returncode, 'events': []}
    try:
        if returncode:
            raise ValueError('diagnostic process failed')
        result['events'] = [json.loads(line.split('workspace_event=', 1)[1])
                            for line in stdout.splitlines() if 'workspace_event=' in line]
        if not result['events']:
            raise ValueError('missing workspace events')
        raw, = [ast.literal_eval(line.split('workspace_intervals=', 1)[1].replace('Some(', '('))
                for line in stdout.splitlines() if 'workspace_intervals=' in line]
        intervals = [[dimension, birth, death] for dimension, birth, death in raw
                     if death is None or death != birth]
        if canonical({'intervals': intervals}) != canonical(anchor):
            raise ValueError('diagnostic interval multiset mismatch')
        for event in result['events']:
            if any(type(value) is not int or value < 0 for value in event['vec_capacity_bytes'].values()):
                raise ValueError('invalid vector capacity observation')
        result['intervals_validated'] = True
    except (ValueError, TypeError, KeyError, SyntaxError) as error:
        result['comparison_error'] = str(error)
    return result


def run(args):
    if sys.platform != 'linux' or args.samples < 10:
        raise ValueError('Linux and at least ten independent measured processes per cell required')
    if args.cpu not in os.sched_getaffinity(0):
        raise ValueError('CPU must belong to the inherited allowed affinity')
    identity = provenance('HEAD')
    initial = fingerprint(REPO)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    # Linux-local executable storage avoids DrvFS executable paging effects.
    local = Path('/tmp') / ('cocycle-workspace-' + output.name + '-' + identity['kernel_commit'][:12])
    local.mkdir(exist_ok=False)
    reference = json.loads(args.native_environment.read_text())
    if reference['pins'] != PINS:
        raise ValueError('native references differ from current source pins')
    native = args.native_environment.parent / 'build'
    for name in ('gudhi', 'ripser'):
        if sha256(native / name) != reference['binaries_sha256'][name]:
            raise ValueError('native reference binary hash mismatch')
        shutil.copy2(native / name, local / name)
    variants = {}
    if args.reuse_build:
        previous = json.loads((args.reuse_build / 'environment.json').read_text())
        if (previous['kernel_commit'] != identity['kernel_commit']
                or previous['controller_sha256'] != sha256(Path(__file__))
                or previous['native_environment_sha256'] != sha256(args.native_environment)):
            raise ValueError('repeat requires exactly the same frozen kernel, controller and references')
        for name, item in previous['variants'].items():
            if fingerprint(args.reuse_build / name) != item['source_sha256']:
                raise ValueError('generated source changed before repeat')
            for kind, key in (('binary', 'binary_sha256'), ('trace', 'trace_sha256')):
                original = Path(item[kind])
                if sha256(original) != item[key]:
                    raise ValueError('frozen executable changed before repeat')
                destination = local / (name if kind == 'binary' else name + '-trace')
                shutil.copy2(original, destination)
                item[kind] = str(destination)
            variants[name] = item
    else:
        with (output / 'build.log').open('w') as log:
            for name in ('joint', 'lifetime', 'scratch'):
                variants[name] = build_variant(name, output, local, log)
    # Finish all builds before any measured process; no compilation during runs.
    environment = {**identity, 'source_sha256': initial, 'variants': variants,
                   'native_environment_sha256': sha256(args.native_environment), 'pins': PINS,
                   'native_environment': reference,
                   'harness_files_sha256': {name: sha256(REPO / name) for name in
                                            ('benches/pipeline/cocycle.rs', 'tools/benchmark_workspace.py',
                                             'tools/benchmark_rips_pipeline.py', 'tools/benchmark_inputs.py',
                                             'tools/compare_rips.py', 'tools/build_native.py')},
                   'controller_sha256': sha256(Path(__file__)), 'samples': args.samples,
                   'order_seed': args.order_seed,
                   'reuse_build': str(args.reuse_build) if args.reuse_build else None,
                   'warmups': 1, 'cpu': args.cpu, 'allowed_cpus': sorted(os.sched_getaffinity(0)),
                   'cpuinfo': Path('/proc/cpuinfo').read_text().split('\n\n')[0],
                   'platform': platform.platform(),
                   'timeout_seconds': 30, 'address_space_mib': 2048,
                   'rustc': subprocess.check_output(['rustc', '-Vv'], text=True),
                   'started_utc': datetime.now(timezone.utc).isoformat(),
                   'protocol': 'cocycle-workspace-v1', 'worker_protocol': 'cocycle-rips-pipeline-v2',
                   'parallel_workers': False,
                   'scope': 'controlled M1 ablations, not production admission',
                   'direct_allocator_bytes': None, 'allocator_calls': None,
                   'frequency_and_host_load': 'uncontrolled'}
    (output / 'environment.json').write_text(json.dumps(environment, indent=2) + '\n')
    fixtures = output / 'fixtures'
    fixtures.mkdir()
    records = []
    names = [*variants, 'gudhi', 'ripser']
    for index, case in enumerate(workspace_cases()):
        path = fixtures / (case.name + '.txt')
        case.fixture.write(path)
        record = {'case': case.name, 'dimension': case.fixture.q,
                  'fixture_sha256': sha256(path), 'samples': [], 'traces': {}}
        records.append(record)
        anchors = {}
        for entry in schedule(names, args.samples, args.order_seed + index):
            name = entry['backend']
            sample = worker([str(local / name), str(path), case.path, case.layout,
                             str(case.epsilon), 'no'], case, 30, 2048, args.cpu)
            sample.update(entry)
            record['samples'].append(sample)
            if sample['status'] == 'completed':
                try:
                    # Native references may omit coverage or representative
                    # metadata shared by two later Rust samples.
                    for anchor in [sample, *anchors.values()]:
                        compare_samples(anchor, sample, case)
                    anchors.setdefault(name, sample)
                except ValueError as error:
                    sample['comparison_error'] = str(error)
        for name in variants:
            env = dict(os.environ, COCYCLE_H2_FIXTURE=str(path), COCYCLE_H2_ROUTE='dispatch',
                       COCYCLE_WORKSPACE_TRACE='1')
            trace_path = output / (case.name + '-' + name + '.log')
            try:
                result = subprocess.run(['taskset', '-c', str(args.cpu), variants[name]['trace'],
                                         'simplicial::cohomology::profiling::profile_h2', '--ignored',
                                         '--nocapture', '--test-threads=1'], env=env,
                                        capture_output=True, text=True, timeout=30)
                trace_path.write_text(result.stdout + result.stderr)
                record['traces'][name] = trace_result(
                    result.stdout, result.returncode, next(iter(anchors.values()), None))
            except subprocess.TimeoutExpired as error:
                trace_path.write_text(repr(error))
                record['traces'][name] = {'events': [], 'comparison_error': 'diagnostic timeout'}
            record['traces'][name]['log_sha256'] = sha256(trace_path)
        valid_traces = all(trace.get('intervals_validated') for trace in record['traces'].values())
        rows = record['samples'] if valid_traces else [*record['samples'], {'status': 'trace_error'}]
        record['summary'] = summarize(rows, names)
        (output / 'results.json').write_text(json.dumps(records, indent=2) + '\n')
        print(case.name, {name: row['median_ms'] for name, row in record['summary'].items()}, flush=True)
    passed = all(all(sample['status'] == 'completed' and 'comparison_error' not in sample
                     for sample in record['samples'])
                 and all(trace.get('intervals_validated') for trace in record['traces'].values())
                 for record in records)
    unchanged = fingerprint(REPO) == initial and provenance('HEAD') == identity
    summary = {'status': 'passed' if passed and unchanged else 'failed',
               'cases': len(records), 'processes': sum(len(record['samples']) for record in records),
               'trace_processes': len(records) * len(variants), 'source_unchanged': unchanged,
               'finished_utc': datetime.now(timezone.utc).isoformat()}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    if summary['status'] != 'passed':
        raise ValueError('workspace ablation failed; retained samples cannot establish a benefit')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--native-environment', type=Path, required=True)
    parser.add_argument('--samples', type=int, default=15)
    parser.add_argument('--cpu', type=int, default=0)
    parser.add_argument('--order-seed', type=int, default=20261004)
    parser.add_argument('--reuse-build', type=Path)
    run(parser.parse_args())
