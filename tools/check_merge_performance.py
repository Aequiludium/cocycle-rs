"""Bounded local scalar/prepared merge checks; never claim universal speedups.

Protocol distance-merge-v1: 3 discarded scalar warmups, timed preparation plus
K queries and cleanup, fresh serial processes, balanced positions, shared f64
fixtures, and separate repeat summaries. Runtime outputs must match exactly.
"""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import random
import resource
import statistics
import subprocess

import reproduce_distance_preparation as replay
from benchmark_rips_pipeline import schedule

ROOT = replay.ROOT
REVISIONS = {'main': '7dc517d2cc8127c972ad0bee0fd2c681126c4815',
             'logical': '26866ac9223a7bd5cc7c35286e1090dc5759e509',
             'scoped': '5bb77b4b624580fe3cede4c8d538d06518de1a5d'}
BACKENDS = ('main_scalar', 'logical_scalar', 'scoped_scalar', 'scoped_prepared')


def save(path, data):
    path.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--samples', type=int, default=12)
    parser.add_argument('--cpu', type=int, default=0)
    args = parser.parse_args()
    if args.samples < 12 or args.samples % 4 or args.cpu not in os.sched_getaffinity(0):
        parser.error('use >=12 samples in multiples of four and an allowed CPU')
    if replay.git('status', '--porcelain').strip():
        parser.error('commit the harness before measuring')
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=False)
    worker = ROOT / 'benches/distances/merge_check.rs'
    identities = [worker, Path(__file__), ROOT / 'tools/reproduce_distance_preparation.py',
                  ROOT / 'tools/benchmark_rips_pipeline.py']
    meta = {'protocol': 'distance-merge-v1', 'revisions': REVISIONS,
            'harness_commit': replay.git('rev-parse', 'HEAD').decode().strip(),
            'harness_hashes': {str(p.relative_to(ROOT)): replay.sha(p.read_bytes()) for p in identities},
            'samples': args.samples, 'repeats': 2, 'warmups': 3, 'cpu': args.cpu,
            'inherited_affinity': sorted(os.sched_getaffinity(0)), 'seed': 20261009,
            'frequency_and_host_load': 'uncontrolled', 'status': 'building',
            'rustc': subprocess.check_output(['rustc', '-Vv'], text=True),
            'start_utc': datetime.now(timezone.utc).isoformat(), 'commands': []}
    save(out / 'environment.json', meta)
    commands = {}

    def run(command):
        meta['commands'].append(command)
        with (out / 'build.log').open('a') as log:
            subprocess.run(command, cwd=ROOT, check=True, stdout=log, stderr=subprocess.STDOUT)

    try:
        for name, revision in REVISIONS.items():
            source = out / (name + '-source')
            replay.extract_snapshot(replay.git('archive', revision, 'src', 'Cargo.toml',
                                               'Cargo.lock', 'benches'), source)
            target = out / (name + '-build')
            run(['cargo', 'build', '--release', '--locked', '--offline', '--lib',
                 '--manifest-path', str(source / 'Cargo.toml'), '--target-dir', str(target)])
            binary = out / name
            command = ['rustc', '--edition=2024', '-O', '-D', 'warnings', str(worker),
                       '--extern', f'cocycle={target}/release/libcocycle.rlib',
                       '-L', f'dependency={target}/release/deps', '-o', str(binary)]
            if name == 'scoped': command.extend(['--cfg', 'prepared_candidate'])
            run(command)
            commands[name + '_scalar'] = [str(binary)]
            if name == 'scoped': commands['scoped_prepared'] = [str(binary)]
        meta['binary_hashes'] = {name: replay.sha(Path(cmd[0]).read_bytes()) for name, cmd in commands.items()}
        meta['status'] = 'running'
        save(out / 'environment.json', meta)
        records = []
        for size in (32, 512):
            fixture = out / f'fixture-{size}.bin'
            fixture.write_bytes(replay.fixture(size))
            for metric in ('bottleneck', 'w1', 'w2'):
                for k in (1, 8, 64):
                    for repeat in (1, 2):
                        record = {'n': size, 'm': size // 8, 'metric': metric, 'k': k,
                                  'repeat': repeat, 'fixture_sha256': replay.sha(fixture.read_bytes()),
                                  'samples': {name: [] for name in BACKENDS}, 'status': 'running'}
                        entries = [e for e in schedule(BACKENDS, args.samples,
                                   20261009 + size + k * 10 + repeat) if not e['warmup']]
                        record['schedule'] = entries
                        anchor = None
                        for entry in entries:
                            name = entry['backend']
                            command = [*commands[name], str(fixture), metric, str(k),
                                       'prepared' if name == 'scoped_prepared' else 'scalar']

                            def limits():
                                resource.setrlimit(resource.RLIMIT_AS, (2048 * 1024**2,) * 2)
                                os.sched_setaffinity(0, {args.cpu})

                            result = subprocess.run(command, capture_output=True, text=True,
                                                    timeout=30, preexec_fn=limits)
                            sample = {**entry, 'exit_code': result.returncode,
                                      'stdout': result.stdout, 'stderr': result.stderr}
                            record['samples'][name].append(sample)
                            records_snapshot = [*records, record]
                            save(out / 'results.json', records_snapshot)
                            if result.returncode:
                                raise ValueError('failed worker; raw sample retained, rankings withheld')
                            data = json.loads(result.stdout)
                            sample.update(data)
                            if not all(math.isfinite(data[x]) for x in ('value', 'elapsed_ms')) or data['elapsed_ms'] <= 0:
                                raise ValueError('invalid measured output')
                            if anchor is None: anchor = data['value']
                            if data['value'] != anchor:
                                raise ValueError('distance outputs differ; rankings withheld')
                        record['status'] = 'passed'
                        record['summary'] = {name: {
                            'median_ms': statistics.median(s['elapsed_ms'] for s in samples),
                            'min_ms': min(s['elapsed_ms'] for s in samples),
                            'max_ms': max(s['elapsed_ms'] for s in samples)}
                            for name, samples in record['samples'].items()}
                        records.append(record)
                        save(out / 'results.json', records)
                        med = {name: values['median_ms'] for name, values in record['summary'].items()}
                        print(size, metric, k, repeat, '26/main', med['logical_scalar']/med['main_scalar'],
                              '49/26', med['scoped_scalar']/med['logical_scalar'],
                              'prepared/49', med['scoped_prepared']/med['scoped_scalar'], flush=True)
        if replay.git('status', '--porcelain').strip() or any(
                replay.sha(p.read_bytes()) != meta['harness_hashes'][str(p.relative_to(ROOT))] for p in identities):
            raise ValueError('harness changed during measurement')
        meta['status'] = 'passed'
    except (OSError, ValueError, subprocess.SubprocessError) as error:
        meta.update(status='failed', error=str(error))
        raise
    finally:
        meta['end_utc'] = datetime.now(timezone.utc).isoformat()
        save(out / 'environment.json', meta)


if __name__ == '__main__':
    main()
