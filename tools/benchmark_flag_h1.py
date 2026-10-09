"""Paired F2 H1 study: committed Rust baseline/candidate and pinned native Ripser.

Uses cocycle-native-v1 workers unchanged. Each repeat has fresh, position-balanced
processes and no warmups. A failed or incomplete case cannot receive rankings.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import io
import os
from pathlib import Path
import platform
import random
import subprocess
import tarfile

import build_native
from benchmark_inputs import Case, make_case
from benchmark_native import compare, expected, save, summarize, validate, worker
from benchmark_rips_pipeline import provenance, schedule

REPO = build_native.REPO
BACKENDS = ('baseline', 'candidate', 'ripser_cpp')


def workloads():
    cases = [make_case(family, n) for family, n in (
        ('uniform', 32), ('uniform', 64), ('uniform', 256), ('circle', 128),
        ('sphere8', 128), ('nonmetric', 128), ('grid', 144), ('bipartite', 128))]
    cases.append(make_case('bipartite', 128, 1.))
    labels = list(range(128))
    random.Random(7919).shuffle(labels)
    values = [1. if (labels[i] < 64) != (labels[j] < 64) else 2.
              for i in range(128) for j in range(i)]
    cases.append(Case('bipartite_shuffled_h1_128_cutoff', 128, values, cutoff=1.))
    cases.append(Case('equal_h1_128', 128, [1.] * (128 * 127 // 2)))
    return cases


def measure(case, fixture, commands, samples, seed, env, timeout, memory_mib):
    record = {'case': case.name, 'kind': 'benchmark', 'n': case.n,
              'homology_max': case.q, 'cutoff': case.cutoff,
              'fixture_sha256': build_native.sha256(fixture),
              'retained_edges': sum(case.cutoff is None or v <= case.cutoff
                                    for v in case.values), 'results': {}}
    # schedule() reserves round zero for a warmup; omit it entirely here.
    record['schedule'] = [entry for entry in schedule(BACKENDS, samples, seed)
                          if not entry['warmup']]
    for entry in record['schedule']:
        name = entry['backend']
        previous = record['results'].get(name)
        if previous and previous['status'] != 'completed':
            continue
        result = worker([*commands[name], str(fixture)], case, name,
                        timeout, memory_mib, env)
        if result['status'] != 'completed':
            result['samples'] = previous.get('samples', []) if previous else []
            record['results'][name] = result
            continue
        try:
            if previous:
                compare(previous, result)
        except ValueError as error:
            record['results'][name] = {'status': 'mismatch', 'error': str(error),
                                       'samples': previous['samples'],
                                       'observed': result}
            continue
        if not previous:
            previous = {key: result[key] for key in ('status', 'coverage', 'intervals')}
            previous['samples'] = []
            record['results'][name] = previous
        previous['samples'].append({**entry, **{key: value for key, value in result.items()
                                               if key not in ('intervals', 'coverage')}})
    validate(record, expected(case))
    if any(record['results'].get(name, {}).get('status') != 'completed'
           or len(record['results'][name].get('samples', [])) != samples
           for name in BACKENDS):
        record['validation'] = 'failed'
        record['error'] = 'incomplete or failed sample group; all rankings withheld'
    return record


def archive_baseline(revision, destination):
    paths = ('Cargo.toml', 'Cargo.lock', 'src', 'benches/rips.rs',
             'benches/native/cocycle.rs')
    data = subprocess.check_output(['git', 'archive', revision, *paths], cwd=REPO)
    destination.mkdir()
    with tarfile.open(fileobj=io.BytesIO(data)) as archive:
        for member in archive.getmembers():
            resolved = (destination / member.name).resolve()
            if not resolved.is_relative_to(destination.resolve()) or not (
                    member.isfile() or member.isdir()):
                raise ValueError('unexpected baseline archive entry')
        archive.extractall(destination, filter='data')


def fingerprint(root, paths):
    files = sorted(path for item in paths for path in
                   ((root / item).rglob('*') if (root / item).is_dir() else [root / item])
                   if path.is_file())
    inventory = {path.relative_to(root).as_posix(): build_native.sha256(path) for path in files}
    digest = hashlib.sha256()
    for name, value in inventory.items():
        digest.update(name.encode() + b'\0' + value.encode() + b'\n')
    return {'sha256': digest.hexdigest(), 'files': inventory}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--baseline', required=True)
    parser.add_argument('--ripser-source', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--cpu', type=int, required=True)
    parser.add_argument('--samples', type=int, default=12)
    parser.add_argument('--repeats', type=int, default=2)
    parser.add_argument('--timeout', type=float, default=30.)
    parser.add_argument('--memory-mib', type=int, default=2048)
    args = parser.parse_args()
    if args.samples < 12 or args.samples % len(BACKENDS) or args.repeats < 2:
        parser.error('use at least 12 samples in multiples of three and two repeats')
    if args.cpu not in os.sched_getaffinity(0):
        parser.error('CPU is outside inherited affinity')
    identity = provenance('HEAD')
    baseline = subprocess.check_output(['git', 'rev-parse', '--verify',
                                       args.baseline + '^{commit}'], cwd=REPO, text=True).strip()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    kernel_paths = ('Cargo.toml', 'Cargo.lock', 'src')
    harness_paths = ('tools/benchmark_flag_h1.py', 'tools/benchmark_native.py',
                     'tools/benchmark_inputs.py', 'tools/benchmark_rips_pipeline.py',
                     'tools/build_native.py', 'benches/native', 'benches/protocol.md')
    sources = fingerprint(REPO, (*kernel_paths, *harness_paths))
    metadata = {**identity, 'baseline_commit': baseline, 'sources': sources,
                'kernel': fingerprint(REPO, kernel_paths), 'protocol': 'cocycle-native-v1',
                'start_utc': datetime.now(timezone.utc).isoformat(), 'status': 'building',
                'samples_per_cell': args.samples, 'repeats': args.repeats,
                'seed': 20261008, 'warmups': 0, 'timeout_seconds': args.timeout,
                'address_space_mib': args.memory_mib, 'platform': platform.platform(),
                'cpu': args.cpu, 'inherited_affinity': sorted(os.sched_getaffinity(0)),
                'frequency_and_host_load': 'uncontrolled', 'build_commands': [],
                'rustc': subprocess.check_output(['rustc', '-Vv'], text=True),
                'cxx': subprocess.check_output(['c++', '--version'], text=True),
                'rustflags': os.environ.get('RUSTFLAGS'),
                'cpuinfo': Path('/proc/cpuinfo').read_text()}
    save(output / 'environment.json', metadata)
    env = {**os.environ, 'OMP_NUM_THREADS': '1', 'OPENBLAS_NUM_THREADS': '1'}
    metadata['thread_environment'] = {key: env[key] for key in
                                      ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS')}

    def run(command, cwd=REPO):
        metadata['build_commands'].append({'argv': command, 'cwd': str(cwd)})
        with (output / 'build.log').open('a') as log:
            subprocess.run(command, cwd=cwd, check=True, stdout=log, stderr=subprocess.STDOUT)

    try:
        baseline_root = output / 'baseline-source'
        archive_baseline(baseline, baseline_root)
        metadata['baseline_kernel'] = fingerprint(baseline_root, kernel_paths)
        binaries = {}
        for name, root in (('baseline', baseline_root), ('candidate', REPO)):
            target = output / (name + '-build')
            run(['cargo', 'build', '--release', '--locked', '--offline', '--lib',
                 '--target-dir', str(target)], root)
            binary = output / name
            # One candidate worker adapter for both Rust revisions.
            run(['rustc', '--edition=2024', '-C', 'opt-level=3', '-D', 'warnings',
                 str(REPO / 'benches/native/cocycle.rs'), '--extern',
                 f'cocycle={target}/release/libcocycle.rlib', '-L',
                 f'dependency={target}/release/deps', '-o', str(binary)])
            binaries[name] = binary
        ripser = build_native.checked_source(args.ripser_source, 'ripser')
        original = ripser / 'ripser.cpp'
        if build_native.sha256(original) != build_native.PINS['ripser']['source_sha256']:
            raise ValueError('Ripser source bytes differ from pin')
        transformed = output / 'ripser_instrumented.cpp'
        transformed.write_text(build_native.instrument_ripser(original.read_text()))
        binaries['ripser_cpp'] = output / 'ripser_cpp'
        run(['c++', '-std=c++17', '-O3', '-DNDEBUG', '-I', str(output),
             str(REPO / 'benches/native/ripser.cpp'), '-o', str(binaries['ripser_cpp'])])
        metadata['ripser_pin'] = build_native.PINS['ripser']
        metadata['instrumented_ripser_sha256'] = build_native.sha256(transformed)
        metadata['binaries_sha256'] = {name: build_native.sha256(path)
                                       for name, path in binaries.items()}
        commands = {name: ['taskset', '-c', str(args.cpu), str(path)]
                    for name, path in binaries.items()}
        metadata['status'] = 'running'
        save(output / 'environment.json', metadata)
        fixtures = output / 'fixtures'
        fixtures.mkdir()
        cases = workloads()
        records = []
        for index, case in enumerate(cases):
            fixture = fixtures / (case.name + '.bin')
            case.write(fixture)
            for repeat in range(args.repeats):
                record = measure(case, fixture, commands, args.samples,
                                 20261008 + index * 100 + repeat, env,
                                 args.timeout, args.memory_mib)
                record['repeat'] = repeat + 1
                records.append(record)
                save(output / 'results.json', records)
                print(case.name, repeat + 1, record['validation'], flush=True)
        if provenance('HEAD') != identity or fingerprint(REPO, (*kernel_paths, *harness_paths)) != sources:
            raise ValueError('sources changed during measurement; rankings invalid')
        for repeat in range(args.repeats):
            summarize(output / f'summary-repeat-{repeat + 1}.csv',
                      [r for r in records if r['repeat'] == repeat + 1])
        metadata['status'] = 'passed' if all(r['validation'] == 'passed' for r in records) else 'failed'
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        metadata['status'] = 'failed'
        metadata['error'] = str(error)
        raise
    finally:
        metadata['end_utc'] = datetime.now(timezone.utc).isoformat()
        save(output / 'environment.json', metadata)
    return 0 if metadata['status'] == 'passed' else 1


if __name__ == '__main__':
    raise SystemExit(main())
