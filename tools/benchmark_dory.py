"""Pinned serial Dory H2 CLI resource snapshots, separate from pipeline timings.

Includes process startup, CSV parsing and PD file output. No cross-library speed
ranking is valid against the pipeline's in-process elapsed_ms boundary.
"""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import platform
import statistics
import subprocess

from benchmark_rips_pipeline import phase3_cases, provenance
from build_native import checked_source, sha256


def csv_input(fixture, path):
    if fixture.mode != 'dense':
        raise ValueError('this adapter uses square matrices with explicit vertex counts')
    def value(a, b):
        return 0. if a == b else fixture.data[max(a, b) * (max(a, b) - 1) // 2 + min(a, b)]
    path.write_text('\n'.join(','.join(format(value(a, b), '.17f') for b in range(fixture.n))
                             for a in range(fixture.n)) + '\n')


def intervals(prefix, dimension):
    path = Path(str(prefix) + f'H{dimension}_pers_data.txt')
    rows = []
    for line in path.read_text().splitlines():
        birth, death = map(float, line.split(','))
        if not math.isfinite(birth) or not math.isfinite(death):
            raise ValueError('nonfinite Dory output')
        if death != birth:
            rows.append([dimension, birth, None if death == -1 else death])
    return rows


def run(args):
    if args.samples < 1 or args.timeout <= 0:
        raise ValueError('positive samples and timeout required')
    identity = provenance('HEAD')
    source = checked_source(args.dory_source, 'dory')
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    original = source / 'Dory/Dory.c'
    generated = output / 'Dory.c'
    # Numeric capture only. SAVEPD is an upstream compile-time output option.
    text = original.read_text().replace('%.12lf', '%.17g').replace('%0.12lf', '%.17g')
    generated.write_text(text)
    binary = output / 'dory'
    command = ['gcc', str(generated), '-I', str(source / 'Dory'), '-O3', '-fopenmp', '-pthread', '-DSAVEPD', '-lm', '-o', str(binary)]
    with (output / 'build.log').open('w') as log:
        subprocess.run(command, check=True, stdout=log, stderr=subprocess.STDOUT)
    environment = {**identity, 'source_revision': 'b6d9b081f626adc60b224bb044c04e7abf5cc599',
                   'license': 'MIT', 'original_sha256': sha256(original), 'capture_sha256': sha256(generated),
                   'headers_sha256': {p.name: sha256(p) for p in sorted((source / 'Dory').glob('*.h'))},
                   'binary_sha256': sha256(binary), 'build_command': command,
                   'compiler': subprocess.check_output(['gcc', '--version'], text=True),
                   'platform': platform.platform(), 'scope': __doc__, 'threads': 1,
                   'samples': args.samples, 'warmups': 1, 'started_utc': datetime.now(timezone.utc).isoformat()}
    (output / 'environment.json').write_text(json.dumps(environment, indent=2) + '\n')
    records = []
    references = json.loads(args.pipeline_results.read_text())
    if any(r['validation'] != 'passed' for r in references):
        raise ValueError('pipeline correctness must pass before Dory measurement')
    # Square-matrix baseline covers positive distinct inputs, ties and a nontrivial
    # H2 sphere. Graph/zero-edge capabilities need a separate adapter audit.
    selected = [case.fixture for case in phase3_cases(args.quick)
                if case.fixture.q == 2 and case.fixture.mode == 'dense'
                and any(f'p3_{name}' in case.fixture.name for name in
                        ('uniform', 'circle', 'sphere', 'noisy_sphere', 'clusters', 'equal_clique'))]
    # Octahedron has a known finite H2 interval [1,2), not just an empty PD check.
    from compare_rips import Fixture
    selected.append(Fixture('octahedron', 'dense', 6, 2, None,
                            [2. if a // 2 == b // 2 else 1. for b in range(6) for a in range(b)]))
    if args.case:
        selected = [f for f in selected if f.name in args.case]
        if {f.name for f in selected} != set(args.case):
            raise ValueError('unknown case; select exact names from the frozen workload matrix')
    for fixture in selected:
        directory = output / fixture.name
        directory.mkdir()
        csv = directory / 'input.csv'
        csv_input(fixture, csv)
        record = {'name': fixture.name, 'csv_sha256': sha256(csv), 'samples': []}
        records.append(record)
        for sample in range(args.samples + 1):
            prefix = directory / f'sample-{sample}-'
            timing = directory / f'time-{sample}.txt'
            cmd = ['/usr/bin/time', '-f', '%e %M', '-o', str(timing), str(binary), str(csv),
                   str(max(fixture.data)), '0', '1', str(prefix), '2', '0', '0', '1']
            try:
                process = subprocess.run(cmd, capture_output=True, text=True, timeout=args.timeout)
                (directory / f'log-{sample}.txt').write_text(process.stdout + process.stderr)
                if process.returncode:
                    raise ValueError(f'Dory exit {process.returncode}')
                rows = intervals(prefix, 1) + intervals(prefix, 2)
                reference = next(r for r in references if r['name'] == fixture.name + '_dense_lower_diagram')
                expected = [r for r in reference['workers']['cocycle']['samples'][0]['intervals'] if r[0] > 0]
                key = lambda r: (r[0], r[1], math.inf if r[2] is None else r[2])
                if sorted(rows, key=key) != sorted(expected, key=key):
                    raise ValueError(f'H1/H2 multiset mismatch: {rows} != {expected}')
                seconds, rss = timing.read_text().split()
                record['samples'].append({'warmup': sample == 0, 'status': 'completed',
                                          'process_elapsed_seconds': float(seconds), 'peak_rss_kib': int(rss),
                                          'intervals_h1_h2': rows, 'command': cmd})
            except (ValueError, OSError, subprocess.TimeoutExpired) as error:
                record['samples'].append({'warmup': sample == 0, 'status': 'failed', 'error': str(error), 'command': cmd})
                break
        valid = [s for s in record['samples'] if s['status'] == 'completed' and not s['warmup']]
        record['status'] = 'completed' if len(valid) == args.samples else 'failed'
        record['median_process_seconds'] = statistics.median(s['process_elapsed_seconds'] for s in valid) if record['status'] == 'completed' else None
        (output / 'results.json').write_text(json.dumps(records, indent=2) + '\n')
        print(fixture.name, record['status'], flush=True)
    unchanged = provenance('HEAD') == identity
    summary = {'status': 'completed' if unchanged and all(r['status'] == 'completed' for r in records) else 'failed',
               'cases': len(records), 'sources_unchanged': unchanged,
               'validation': 'exact H1/H2 multiset comparison with validated pipeline outputs',
               'finished_utc': datetime.now(timezone.utc).isoformat()}
    (output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    if summary['status'] != 'completed':
        raise SystemExit(1)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--dory-source', type=Path, required=True)
    parser.add_argument('--pipeline-results', type=Path, required=True)
    parser.add_argument('--samples', type=int, default=12)
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('--timeout', type=float, default=30.)
    parser.add_argument('--case', action='append', help='select audited cases; omitted cases are not passes')
    run(parser.parse_args())
