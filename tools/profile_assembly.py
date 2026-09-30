"""Measure private persistence assembly in generated copies of a committed kernel.

cocycle-rips-assembly-v1 is a diagnostic protocol, separate from pipeline v2.
Ordinary and instrumented workers use identical inputs and public calls. Only
the generated source contains hooks; production Rust and public APIs are unchanged.
"""

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import subprocess
import sys

import benchmark_rips_pipeline as pipeline
from benchmark_inputs import distances, float32_case, performance_cases
from build_native import REPO, sha256
from compare_rips import Fixture

PROTOCOL = 'cocycle-rips-assembly-v1'
R0 = 'b2c3bd5eebf0c1193f30ea651012b8ae61b274c1'
PREFIX = 'COCYCLE_ASSEMBLY '

# This support is appended only to the generated crate's private persistence owner.
SUPPORT = r'''
std::thread_local! {
    static ASSEMBLY_RAW: std::cell::Cell<(usize, usize, usize)> = const { std::cell::Cell::new((0, 0, 0)) };
    static ASSEMBLY_KERNEL: std::cell::Cell<(&'static str, u64, usize)> = const { std::cell::Cell::new(("missing", 0, 0)) };
    static ASSEMBLY_PHASES: std::cell::Cell<(u64, u64)> = const { std::cell::Cell::new((0, 0)) };
    static ASSEMBLY_RESULT: std::cell::Cell<(u64, usize, usize, usize, u64, usize)> = const { std::cell::Cell::new((0, 0, 0, 0, 0, 0)) };
}

pub(crate) struct AssemblyKernelScope(std::time::Instant, &'static str);
impl AssemblyKernelScope {
    pub(crate) fn new(route: &'static str) -> Self { Self(std::time::Instant::now(), route) }
}
impl Drop for AssemblyKernelScope {
    fn drop(&mut self) {
        let elapsed = self.0.elapsed().as_nanos() as u64;
        ASSEMBLY_KERNEL.with(|cell| {
            let previous = cell.get();
            cell.set((self.1, previous.1 + elapsed, previous.2 + 1));
        });
    }
}

fn assembly_raw(raw: RawIntervals) -> RawIntervals {
    ASSEMBLY_RAW.with(|cell| {
        let previous = cell.get();
        cell.set((raw.len(), raw.capacity(), previous.2 + 1));
    });
    raw
}

pub(crate) fn assembly_diagram_phases(validation: u64, sorting: u64) {
    ASSEMBLY_PHASES.with(|cell| cell.set((validation, sorting)));
}

/// Emit diagnostics after all worker clocks and process-memory readings.
#[doc(hidden)]
pub fn emit_assembly_profile() {
    let (raw_count, raw_capacity, raw_calls) = ASSEMBLY_RAW.with(|cell| cell.get());
    let (route, kernel_ns, kernel_calls) = ASSEMBLY_KERNEL.with(|cell| cell.get());
    let (assembly_ns, final_count, final_capacity, assembly_calls, materialization_ns, growths) = ASSEMBLY_RESULT.with(|cell| cell.get());
    let (validation_ns, sorting_ns) = ASSEMBLY_PHASES.with(|cell| cell.get());
    eprintln!(concat!(
        "COCYCLE_ASSEMBLY {{\"route\":\"{}\",\"kernel_ns\":{},\"kernel_calls\":{},",
        "\"assembly_ns\":{},\"assembly_calls\":{},\"raw_calls\":{},\"raw_count\":{},",
        "\"raw_capacity\":{},\"final_count\":{},\"final_capacity\":{},",
        "\"raw_element_bytes\":{},\"final_element_bytes\":{},",
        "\"materialization_ns\":{},\"diagram_validation_ns\":{},\"sorting_ns\":{},\"capacity_growths\":{}}}"
    ), route, kernel_ns, kernel_calls, assembly_ns, assembly_calls, raw_calls,
       raw_count, raw_capacity, final_count, final_capacity,
       std::mem::size_of::<(usize, f64, Option<f64>)>(),
       std::mem::size_of::<PersistenceInterval>(), materialization_ns, validation_ns, sorting_ns, growths);
}
'''


def replace_once(source, before, after):
    """Reject changed or ambiguous private markers instead of guessing."""
    if source.count(before) != 1:
        raise ValueError(f'expected one assembly marker: {before[:100]}')
    return source.replace(before, after)


def scope(source, marker, route):
    if source.count(marker) != 1:
        raise ValueError(f'expected one kernel marker: {marker}')
    start = source.index('{', source.index(marker)) + 1
    return source[:start] + f'\n    let _assembly_kernel = crate::persistence::AssemblyKernelScope::new("{route}");' + source[start:]


def instrument(source, name, decompose=False):
    """Make only observational changes in the generated source copy."""
    if name == 'persistence/mod.rs':
        source = replace_once(source, '    let mut intervals = Vec::new();',
            '    let _assembly_start = std::time::Instant::now();\n'
            + ('    let mut _assembly_growths = 0usize;\n' if decompose else '    let _assembly_growths = 0usize;\n')
            + '    let mut intervals = Vec::new();')
        if decompose:
            # Record changes at the reserve operation, not allocator-call counts.
            source = replace_once(source, '        intervals\n            .try_reserve(1)',
                '        let _previous_capacity = intervals.capacity();\n        intervals\n            .try_reserve(1)')
            source = replace_once(source, '        intervals.push(interval);',
                '        if intervals.capacity() != _previous_capacity { _assembly_growths += 1; }\n        intervals.push(interval);')
        finish = '''    let _final_count = intervals.len();
    let _final_capacity = intervals.capacity();
    let _materialization_ns = MATERIALIZATION;
    let diagram = PersistenceDiagram::new(max_dimension, coverage, intervals)?;
    let _assembly_ns = _assembly_start.elapsed().as_nanos() as u64;
    ASSEMBLY_RESULT.with(|cell| {
        let previous = cell.get();
        cell.set((_assembly_ns, _final_count, _final_capacity, previous.3 + 1, _materialization_ns, _assembly_growths));
    });
    Ok(diagram)'''.replace('MATERIALIZATION', '_assembly_start.elapsed().as_nanos() as u64' if decompose else '0')
        support = SUPPORT if decompose else replace_once(SUPPORT, '''pub(crate) fn assembly_diagram_phases(validation: u64, sorting: u64) {
    ASSEMBLY_PHASES.with(|cell| cell.set((validation, sorting)));
}
''', '')
        return replace_once(source, '    PersistenceDiagram::new(max_dimension, coverage, intervals)', finish) + support
    if name == 'persistence/rips/mod.rs':
        return replace_once(source, 'assemble_diagram(options.max_dimension(), coverage, raw)',
                            'assemble_diagram(options.max_dimension(), coverage, super::assembly_raw(raw))')
    if name == 'persistence/simplicial/mod.rs':
        source = replace_once(source, '                implicit(budget)?,',
                              '                super::assembly_raw(implicit(budget)?),')
        before = '''                cohomology::compute(
                    &access,
                    options.max_homology_dimension(),
                    options.field(),
                    budget,
                )?,'''
        return replace_once(source, before, before.replace('cohomology::compute(', 'super::assembly_raw(cohomology::compute(').replace(')?,', ')?),'))
    if name == 'persistence/boundary/mod.rs':
        source = scope(source, 'fn diagram<I>(', 'boundary_reduction')
        return replace_once(source, '    super::assemble_diagram(dimension, coverage, intervals)',
            '    drop(_assembly_kernel);\n    super::assemble_diagram(dimension, coverage, super::assembly_raw(intervals))')
    routes = {'persistence/flag/h0.rs': 'h0_union_find',
              'persistence/flag/cohomology/mod.rs': 'specialized_f2_h1',
              'persistence/simplicial/cohomology.rs': 'generic_cohomology'}
    if name in routes:
        return scope(source, 'fn compute(', routes[name])
    if name == 'diagram/persistence_diagram.rs' and decompose:
        source = replace_once(source, '        let coverage = match coverage {',
            '        let _validation_start = std::time::Instant::now();\n        let coverage = match coverage {')
        return replace_once(source, '        intervals.sort_unstable_by(compare_intervals);',
            '        let _validation_ns = _validation_start.elapsed().as_nanos() as u64;\n'
            '        let _sorting_start = std::time::Instant::now();\n'
            '        intervals.sort_unstable_by(compare_intervals);\n'
            '        crate::persistence::assembly_diagram_phases(_validation_ns, _sorting_start.elapsed().as_nanos() as u64);')
    return source


def worker_source(source, profiled):
    # Supplied nonzero vertex births force the ordinary boundary path, without
    # changing dispatch or creating a library API. Each disjoint triangle fills.
    boundary = '''    let result: PersistenceResult = if args[2] == "boundary" {
        let complex = timed(&mut phases[1], || -> cocycle::Result<_> {
            use cocycle::complex::{Simplex, SimplicialComplex};
            let mut simplices = Vec::new();
            for vertex in 0..n { simplices.push(Simplex::new(vec![vertex], 1.)?); }
            for edge in edges { simplices.push(Simplex::new(edge.vertices.to_vec(), edge.value)?); }
            for a in (0..n).step_by(3) { simplices.push(Simplex::new(vec![a, a + 1, a + 2], 3.)?); }
            SimplicialComplex::new(simplices)
        })?;
        simplex_count = Some(complex.len());
        timed(&mut phases[3], || analyze(&complex, &opts, &[], &limits))?
    } else if mode == "flag" {'''
    source = replace_once(source, '    let result: PersistenceResult = if mode == "flag" {', boundary)
    if profiled:
        source = replace_once(source, '    std::hint::black_box(result);',
                              '    cocycle::persistence::emit_assembly_profile();\n    std::hint::black_box(result);')
    return source


def cases(quick=False):
    for original in performance_cases(quick):
        if original.mode == 'points':
            coordinates = original.values
            data = distances([coordinates[i:i + original.d] for i in range(0, len(coordinates), original.d)])
            yield pipeline.Case(Fixture(original.name, 'dense', original.n, original.q, original.cutoff, data,
                                        precision='f64 Euclidean points'), 'points', coordinates=coordinates)
        else:
            original = float32_case(original)
            yield pipeline.Case(Fixture(original.name, 'dense', original.n, original.q, original.cutoff,
                                        original.values), 'dense')
    for n in ([32] if quick else [10000, 100000]):
        edges = [[i - 1, i, 1. + (i % 7) / 8.] for i in range(1, n)]
        yield pipeline.Case(Fixture(f'forest{n}_h0', 'flag', n, 0, None, edges), 'flag')
        yield pipeline.Case(Fixture(f'isolates{n}_h0', 'flag', n, 0, None, []), 'flag')
    for side in ([4] if quick else [32, 64, 128]):
        edges = [[a, b, 1.] for a in range(side) for b in range(side, side * 2)]
        yield pipeline.Case(Fixture(f'bipartite{side}_h1', 'flag', side * 2, 1, None, edges), 'flag')
    for blocks in ([4] if quick else [16, 128, 1024]):
        edges = [[a + i, a + j, 2.] for a in range(0, blocks * 3, 3) for i, j in ((0, 1), (0, 2), (1, 2))]
        yield pipeline.Case(Fixture(f'filled_triangles{blocks}', 'flag', blocks * 3, 1, None, edges,
                                    characteristic=3), 'boundary')
    for case in pipeline.cases(quick):
        if case.path == 'expanded' and case.fixture.name.startswith(('nonmetric', 'sphere_h2')):
            yield case


def expected(case):
    """Independent analytic interval multisets for output-heavy controls."""
    name, n = case.fixture.name, case.fixture.n
    if name.startswith('forest'):
        return [[0, 0., edge[2]] for edge in case.fixture.data] + [[0, 0., None]]
    if name.startswith('isolates'):
        return [[0, 0., None]] * n
    if name.startswith('bipartite'):
        side = n // 2
        return [[0, 0., 1.]] * (n - 1) + [[0, 0., None]] + [[1, 1., None]] * ((side - 1) ** 2)
    if name.startswith('filled_triangles'):
        blocks = n // 3
        return [[0, 1., 2.]] * (2 * blocks) + [[0, 1., None]] * blocks + [[1, 2., 3.]] * blocks
    return None


def profile(stderr, sample):
    lines = [line[len(PREFIX):] for line in stderr.splitlines() if line.startswith(PREFIX)]
    if len(lines) != 1:
        raise ValueError('expected one assembly diagnostic record')
    result = json.loads(lines[0])
    keys = {'route', 'kernel_ns', 'kernel_calls', 'assembly_ns', 'assembly_calls', 'raw_calls',
            'raw_count', 'raw_capacity', 'final_count', 'final_capacity', 'raw_element_bytes',
            'final_element_bytes', 'materialization_ns', 'diagram_validation_ns', 'sorting_ns', 'capacity_growths'}
    if not isinstance(result, dict) or set(result) != keys:
        raise ValueError('unexpected assembly diagnostic fields')
    if result['route'] not in ('h0_union_find', 'specialized_f2_h1', 'generic_cohomology', 'boundary_reduction'):
        raise ValueError('missing or unknown producer route')
    for key, value in result.items():
        if key != 'route' and (type(value) is not int or value < 0):
            raise ValueError(f'invalid assembly counter: {key}')
    if any(result[key] != 1 for key in ('raw_calls', 'kernel_calls', 'assembly_calls')):
        raise ValueError('expected exactly one raw producer and assembly call')
    if result['final_count'] != len(sample['intervals']) or result['raw_count'] < result['final_count']:
        raise ValueError('raw/final interval counts disagree')
    if result['raw_capacity'] < result['raw_count'] or result['final_capacity'] < result['final_count']:
        raise ValueError('capacity smaller than interval count')
    analysis_ns = sample['phases_ms'][3] * 1_000_000
    private_ns = result['assembly_ns'] + result['kernel_ns']
    if analysis_ns <= 0 or private_ns <= 0 or private_ns > analysis_ns + 100:
        raise ValueError('private phases exceed the same-call public analysis duration')
    if sum(result[k] for k in ('materialization_ns', 'diagram_validation_ns', 'sorting_ns')) > result['assembly_ns']:
        raise ValueError('assembly decomposition exceeds enclosing assembly duration')
    if not result['raw_element_bytes'] or not result['final_element_bytes']:
        raise ValueError('zero interval element size')
    if result['capacity_growths'] > result['final_count']:
        raise ValueError('more capacity growths than materialized intervals')
    result['analysis_ms'] = analysis_ns / 1_000_000
    result['assembly_fraction'] = result['assembly_ns'] / analysis_ns
    result['kernel_assembly_fraction'] = result['assembly_ns'] / private_ns
    result['public_remainder_ns'] = max(0, analysis_ns - private_ns)
    result['raw_bytes'] = result['raw_capacity'] * result['raw_element_bytes']
    result['final_bytes'] = result['final_capacity'] * result['final_element_bytes']
    result['coexisting_output_capacity_bytes'] = result['raw_bytes'] + result['final_bytes']
    return result


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def fingerprint():
    return hashlib.sha256((pipeline.source_hash() + sha256(Path(__file__))).encode()).hexdigest()


def build(output, decompose):
    directory = output / 'build'
    directory.mkdir()
    generated = directory / 'profile-source'
    shutil.copytree(REPO / 'src', generated / 'src')
    for name in ('Cargo.toml', 'Cargo.lock'):
        shutil.copyfile(REPO / name, generated / name)
    # Cargo validates explicitly declared targets even for a library-only build.
    (generated / 'benches').mkdir()
    shutil.copyfile(REPO / 'benches/rips.rs', generated / 'benches/rips.rs')
    originals, modified = {}, {}
    for path in sorted((generated / 'src').rglob('*.rs')):
        name = path.relative_to(generated / 'src').as_posix()
        originals['src/' + name] = sha256(REPO / 'src' / name)
        path.write_text(instrument(path.read_text(), name, decompose), newline='\n')
        modified['src/' + name] = sha256(path)
    commands = []
    for variant in ('ordinary', 'profile'):
        worker = directory / f'{variant}.rs'
        worker.write_text(worker_source((pipeline.PIPELINE / 'cocycle.rs').read_text(), variant == 'profile'), newline='\n')
        target = directory / f'{variant}-cargo'
        manifest = (REPO if variant == 'ordinary' else generated) / 'Cargo.toml'
        commands.extend([
            ['cargo', 'build', '--release', '--locked', '--offline', '--lib', '--manifest-path', str(manifest), '--target-dir', str(target)],
            ['rustc', '--edition=2024', '-O', '-D', 'warnings', str(worker), '--extern',
             f'cocycle={target}/release/libcocycle.rlib', '-L', f'dependency={target}/release/deps', '-o', str(directory / variant)],
        ])
    save(directory / 'source.json', {'originals_sha256': originals, 'generated_sha256': modified,
                                     'workers_sha256': {v: sha256(directory / f'{v}.rs') for v in ('ordinary', 'profile')}})
    with (directory / 'build.log').open('w') as log:
        for command in commands:
            log.write(json.dumps(command) + '\n')
            log.flush()
            subprocess.run(command, cwd=REPO, check=True, stdout=log, stderr=subprocess.STDOUT, timeout=240)
    return {'commands': commands, 'binaries_sha256': {v: sha256(directory / v) for v in ('ordinary', 'profile')},
            'rustc': subprocess.check_output(['rustc', '-Vv'], text=True),
            'cargo': subprocess.check_output(['cargo', '-V'], text=True)}


def measure(case, args, directory, case_index):
    fixture = directory / 'fixtures' / (case.name + '.txt')
    case.fixture.write(fixture)
    if case.coordinates is not None:
        Path(str(fixture) + '.points').write_text(' '.join(map(str, case.coordinates)) + '\n')
    record = {'name': case.name, 'n': case.fixture.n, 'q': case.fixture.q, 'field': case.fixture.characteristic,
              'path': case.path, 'layout': case.layout, 'cutoff': case.fixture.cutoff,
              'input_edges_or_distances': len(case.fixture.data), 'fixture_sha256': sha256(fixture),
              'points_sha256': sha256(Path(str(fixture) + '.points')) if case.coordinates is not None else None,
              'precision': case.fixture.precision, 'samples': [], 'validation': 'passed'}
    anchor = None
    analytic = expected(case)
    record['oracle'] = 'analytic interval multiset' if analytic is not None else 'ordinary R0 worker; kernel acceptance recorded in baseline'
    raw_dir = directory / 'raw' / case.name
    raw_dir.mkdir(parents=True)
    for entry in pipeline.schedule(('ordinary', 'profile'), args.samples, args.order_seed + case_index):
        variant = entry['backend']
        command = [str(directory / 'build' / variant), str(fixture), case.path, case.layout, str(case.epsilon), 'no']
        sample = pipeline.worker(command, case, args.timeout, args.address_space_mib, args.cpu, retain_stderr=True)
        sample.update(entry)
        if sample['status'] == 'completed':
            try:
                if analytic is not None and pipeline.canonical(sample) != pipeline.canonical({'intervals': analytic}):
                    raise ValueError('independent analytic diagram mismatch')
                if anchor is not None:
                    pipeline.compare_samples(anchor, sample, case)
                    if anchor['intervals'] != sample['intervals']:
                        raise ValueError('canonical interval sequence changed')
                else:
                    anchor = sample
                if variant == 'profile':
                    sample['assembly'] = profile(sample['stderr'], sample)
            except (ValueError, KeyError, TypeError) as error:
                sample['comparison_error'] = str(error)
                record['validation'] = 'failed'
        else:
            record['validation'] = 'failed'
        path = raw_dir / f'{variant}-{entry["round"]:03}.json'
        path.write_text(json.dumps(sample, allow_nan=False) + '\n')
        compact = {k: v for k, v in sample.items() if k not in ('intervals', 'stderr')}
        compact.update(raw_file=path.relative_to(directory).as_posix(), raw_sha256=sha256(path))
        record['samples'].append(compact)
    valid = [s for s in record['samples'] if not s['warmup'] and s['backend'] == 'profile']
    controls = [s for s in record['samples'] if not s['warmup'] and s['backend'] == 'ordinary']
    row = {'case': case.name, 'validation': record['validation'], 'samples': len(valid)}
    if record['validation'] == 'passed':
        assert len(valid) == len(controls) == args.samples
        profiles = [s['assembly'] for s in valid]
        for key in ('route', 'raw_count', 'raw_capacity', 'final_count', 'final_capacity',
                    'raw_bytes', 'final_bytes', 'coexisting_output_capacity_bytes'):
            values = {p[key] for p in profiles}
            if len(values) != 1:
                raise ValueError(f'non-deterministic capacity/count/route: {case.name} {key}')
            row[key] = values.pop()
        metrics = {'kernel_ms': [p['kernel_ns'] / 1_000_000 for p in profiles],
                   'assembly_ms': [p['assembly_ns'] / 1_000_000 for p in profiles],
                   'analysis_ms': [p['analysis_ms'] for p in profiles],
                   'assembly_fraction': [p['assembly_fraction'] for p in profiles],
                   'kernel_assembly_fraction': [p['kernel_assembly_fraction'] for p in profiles],
                   'public_remainder_ms': [p['public_remainder_ns'] / 1_000_000 for p in profiles],
                   'ordinary_analysis_ms': [s['phases_ms'][3] for s in controls]}
        if args.decompose:
            metrics.update({k: [p[k] / 1_000_000 for p in profiles] for k in
                            ('materialization_ns', 'diagram_validation_ns', 'sorting_ns')})
            row['capacity_growths'] = sorted({p['capacity_growths'] for p in profiles})
        for key, values in metrics.items():
            row[key + '_median'] = statistics.median(values)
            row[key + '_min'] = min(values)
            row[key + '_max'] = max(values)
        row['max_peak_rss_kib'] = max(s['peak_rss_kib'] for s in valid)
        row['ordinary_max_peak_rss_kib'] = max(s['peak_rss_kib'] for s in controls)
        row['max_hwm_growth_kib'] = max(s['peak_rss_kib'] - s['hwm_before_kib'] for s in valid)
    return record, row


def run(args):
    if sys.platform != 'linux':
        raise ValueError('Linux required for per-process memory and resource controls')
    if args.samples < (1 if args.quick else 12) or args.samples % (1 if args.quick else 2):
        raise ValueError('formal study requires at least 12 balanced measured samples')
    if args.timeout <= 0 or not math.isfinite(args.timeout) or args.address_space_mib < 128:
        raise ValueError('positive timeout and at least 128 MiB address space required')
    if args.cpu is not None and args.cpu not in os.sched_getaffinity(0):
        raise ValueError('worker CPU is outside allowed affinity')
    if args.exploratory:
        identity = {'kernel_commit': args.kernel_revision, 'harness_commit': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
                    'dirty': True}
    else:
        identity = pipeline.provenance(args.kernel_revision)
    initial = fingerprint()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    (output / 'fixtures').mkdir()
    environment = {**identity, 'protocol_id': PROTOCOL, 'source_sha256': initial,
                   'started_utc': datetime.now(timezone.utc).isoformat(), 'command': [sys.executable, *sys.argv],
                   'platform': platform.platform(), 'cpuinfo': Path('/proc/cpuinfo').read_text().split('\n\n')[0],
                   'allowed_cpus': sorted(os.sched_getaffinity(0)), 'selected_cpu': args.cpu,
                   'thread_policy': 'serial fresh single-threaded Rust workers', 'samples': args.samples, 'warmups': 1,
                   'order_seed': args.order_seed, 'quick': args.quick, 'decompose': args.decompose,
                   'timeout_seconds': args.timeout, 'address_space_mib': args.address_space_mib,
                   'RUSTFLAGS': os.environ.get('RUSTFLAGS'), 'frequency_thermal_host_load_numa': 'uncontrolled',
                   'timing_contract': __doc__, 'evidence_class': 'assembly_diagnostics',
                   'memory_contract': 'owned raw/final Vec capacity overlap during conversion; excludes allocator metadata, realloc transient storage and other solver state',
                   'build': build(output, args.decompose)}
    save(output / 'environment.json', environment)
    records, rows = [], []
    selected = [case for case in cases(args.quick) if not args.cases or any(key in case.name for key in args.cases)]
    if not selected:
        raise ValueError('case filters selected no workload')
    for index, case in enumerate(selected):
        record, row = measure(case, args, output, index)
        records.append(record)
        rows.append(row)
        save(output / 'results.json', records)
        save(output / 'measurements.json', rows)
        print(case.name, record['validation'], flush=True)
    unchanged = fingerprint() == initial
    if not args.exploratory:
        unchanged = unchanged and pipeline.provenance(args.kernel_revision) == identity
    summary = {**identity, 'protocol_id': PROTOCOL, 'sources_unchanged': unchanged,
               'status': 'passed' if unchanged and all(r['validation'] == 'passed' for r in records) else 'failed',
               'cases': len(records), 'samples_per_worker': args.samples,
               'measured_processes': len(records) * args.samples * 2,
               'finished_utc': datetime.now(timezone.utc).isoformat(),
               'failed_cases': [r['name'] for r in records if r['validation'] != 'passed']}
    save(output / 'summary.json', summary)
    if summary['status'] == 'passed':
        fields = list(dict.fromkeys(key for row in rows for key in row))
        with (output / 'matrix.csv').open('w', newline='') as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
    print(json.dumps(summary))
    return int(summary['status'] != 'passed')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--kernel-revision', default=R0)
    parser.add_argument('--samples', type=int, default=12)
    parser.add_argument('--order-seed', type=int, default=1701)
    parser.add_argument('--cpu', type=int)
    parser.add_argument('--timeout', type=float, default=60.)
    parser.add_argument('--address-space-mib', type=int, default=2048)
    parser.add_argument('--quick', action='store_true')
    parser.add_argument('--exploratory', action='store_true')
    parser.add_argument('--decompose', action='store_true')
    parser.add_argument('--cases', nargs='+', help='explicit case-name substrings for conditional Stage B')
    args = parser.parse_args()
    try:
        return run(args)
    except (OSError, ValueError, subprocess.CalledProcessError) as error:
        print(str(error), file=sys.stderr)
        return 1


if __name__ == '__main__':
    raise SystemExit(main())
