"""Generate and measure benchmark-only distance lifecycle copies of fixed R0.

Ordinary R0 remains separately linked. No lifecycle type enters the library API.
"""
import argparse
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import shutil
import statistics
import struct
import subprocess
import sys
import time
import threading
from concurrent.futures import ThreadPoolExecutor

from benchmark_distances import fixtures, FAMILIES, schedule
from build_native import REPO, sha256
from benchmark_rips_pipeline import provenance
from distance_common import tiny_oracle

R0 = 'b2c3bd5eebf0c1193f30ea651012b8ae61b274c1'
PROTOCOL = 'cocycle-distance-lifecycle-v1'
VARIANTS = ('public', 'r0', 'r1', 'r2', 'r3', 'r4', 'r5', 'bounded')
KS = (1, 2, 4, 8, 16, 64, 256, 1024)
METRICS = ('bottleneck', 'w1', 'w2')
PATTERNS = ('single', 'same_pair', 'one_to_many', 'many_to_one', 'repeated_batch', 'all_pairs')


def replace(source, before, after, count=1):
    if source.count(before) != count:
        raise ValueError(f'changed lifecycle marker: {before[:100]}')
    return source.replace(before, after)


def function(source, marker):
    if source.count(marker) != 1:
        raise ValueError(f'changed lifecycle function: {marker}')
    begin = source.index(marker)
    brace = source.index('{', begin)
    depth = 1
    end = brace + 1
    while depth:
        depth += (source[end] == '{') - (source[end] == '}')
        end += 1
    return begin, brace, end


def clock(source, marker, phase):
    _, brace, _ = function(source, marker)
    return source[:brace + 1] + f'\n    let _lifecycle_clock = crate::lifecycle::Clock::new({phase});' + source[brace + 1:]


def drop_vectors(source, declaration, fields, prefix, generics=''):
    body = '\n'.join(f'crate::lifecycle::put("{prefix}-{name}", std::mem::take(&mut self.{name}));' for name in fields)
    return source + f'\nimpl{generics} Drop for {declaration} {{ fn drop(&mut self) {{ {body} }} }}\n'


def instrument(source, name):
    if name == 'bottleneck.rs':
        source = replace(source, "struct Prepared<'a> {", "pub(crate) struct Prepared<'a> {")
        source = replace(source, '    fn bytes(&self) -> usize {', '    pub(crate) fn bytes(&self) -> usize {', 2)
        begin, brace, end = function(source, 'fn solve<const CONTROLLED: bool>(')
        body = source[brace + 1:end - 1]
        marker = '    let total = first'
        split = body.index(marker)
        tail = body[split:]
        source = source[:brace + 1] + body[:split] + '\n    lifecycle_solve(&first, &second, _options, _stats, budget)\n' + source[end - 1:]
        source += '''
pub(crate) fn lifecycle_prepare(points: &[[f64; 2]]) -> Result<Prepared<'_>> {
    Prepared::new(points, &mut WorkBudget::unlimited())
}
pub(crate) fn lifecycle_prepare_with<'a, const C: bool>(points: &'a [[f64; 2]], budget: &mut WorkBudget<'_, C>) -> Result<Prepared<'a>> { Prepared::new(points, budget) }
pub(crate) fn lifecycle_solve<const CONTROLLED: bool>(first: &Prepared<'_>, second: &Prepared<'_>, _options: Options, _stats: &mut Diagnostics, budget: &mut WorkBudget<'_, CONTROLLED>) -> Result<f64> {
''' + tail + '\n}\n'
        source = clock(source, 'fn coordinates<const CONTROLLED: bool>(', 1)
        source = replace(source, '        let mut diagonals = filled(points.len(), 0.0)?;', '        let diagonal_clock = crate::lifecycle::Clock::new(2);\n        let mut diagonals = filled(points.len(), 0.0)?;')
        source = replace(source, '        // Equal-coordinate groups remain contiguous;', '        drop(diagonal_clock);\n        let order_clock = crate::lifecycle::Clock::new(3);\n        // Equal-coordinate groups remain contiguous;')
        source = replace(source, '        let mut representatives = Vec::new();', '        drop(order_clock);\n        let _group_clock = crate::lifecycle::Clock::new(4);\n        let mut representatives = Vec::new();')
        source = clock(source, 'pub(crate) fn lifecycle_solve<const CONTROLLED: bool>(', 9)
        source = clock(source, 'fn candidates<const CONTROLLED: bool>(', 7)
    elif name == 'wasserstein.rs':
        source = replace(source, 'fn solve_prepared<const CONTROLLED: bool>(', 'pub(crate) fn solve_prepared<const CONTROLLED: bool>(')
        source = replace(source, '    first: Vec<Point>,\n    second: Vec<Point>,', '    first: &[Point],\n    second: &[Point],')
        source = replace(source, '        first,\n        second,\n        scale,', '        &first,\n        &second,\n        scale,')
        source = replace(source, 'solve_prepared(first, second, scale,', 'solve_prepared(&first, &second, scale,')
        source = replace(source, 'capacity_bytes(&first).saturating_add(capacity_bytes(&second))', '(first.len() + second.len()) * std::mem::size_of::<Point>()')
        source += '''
pub(crate) use numeric::Point as LifecyclePoint;
pub(crate) fn lifecycle_prepare<const C: bool>(view: &crate::diagram::DiagramDimension<'_>, count: usize, budget: &mut WorkBudget<'_, C>) -> Result<(Vec<Point>, f64, bool)> {
    numeric::collect_points(view.iter().map(|i| (i.birth(), i.end())), count, budget)
}
pub(crate) fn lifecycle_solve<const CONTROLLED: bool>(first: &(Vec<Point>, f64, bool), second: &(Vec<Point>, f64, bool), metric: Metric, stats: &mut Stats, budget: &mut WorkBudget<'_, CONTROLLED>) -> Result<f64> {
    let scale = numeric::scale_for(first.1.max(second.1));
    if scale == 1.0 {
        if !first.2 || !second.2 { return Err(numerical()); }
        solve_prepared(&first.0, &second.0, scale, metric, Options::default(), stats, budget)
    } else {
        let mut a = first.0.clone(); let mut b = second.0.clone();
        numeric::normalize(&mut a, scale, budget)?; numeric::normalize(&mut b, scale, budget)?;
        solve_prepared(&a, &b, scale, metric, Options::default(), stats, budget)
    }
}
'''
        source = clock(source, 'pub(crate) fn solve_prepared<const CONTROLLED: bool>(', 9)
    elif name == 'wasserstein/numeric.rs':
        source = replace(source, 'pub(super) struct Point {', 'pub(crate) struct Point {')
        for marker in ('fn collect_points<const CONTROLLED: bool>(', 'fn normalize<const CONTROLLED: bool>(', 'fn scale_for('):
            source = replace(source, marker, 'pub(super) ' + marker)
        source = clock(source, 'pub(super) fn collect_points<const CONTROLLED: bool>(', 5)
        source = clock(source, 'pub(super) fn normalize<const CONTROLLED: bool>(', 6)
    elif name == 'wasserstein/graph.rs':
        source = clock(source, 'pub(super) fn generate<const CONTROLLED: bool>(', 7)
    elif name == 'mod.rs':
        source = clock(source, "fn dimension_view<'a, const CONTROLLED: bool>(", 0)
        source = replace(source, 'fn check_context(first:', 'pub(crate) fn check_context(first:')
        begin, brace, end = function(source, 'fn distance<const CONTROLLED: bool>(')
        body = source[brace + 1:end - 1]
        start = body.index('        match kind {')
        stop = body.index('\n    };', start)
        # Reuse the original coverage, dimension, essential and accumulation rules.
        body = body[:start] + '        finite(&first, &second, (first_finite, second_finite), budget)?' + body[stop:]
        source += '''
pub(crate) fn lifecycle_distance<const CONTROLLED: bool>(first: &PersistenceDiagram, second: &PersistenceDiagram, dimension: usize, kind: Kind, budget: &mut WorkBudget<'_, CONTROLLED>, finite: impl FnOnce(&DiagramDimension<'_>, &DiagramDimension<'_>, (usize, usize), &mut WorkBudget<'_, CONTROLLED>) -> Result<f64>) -> Result<f64> {
''' + body + '\n}\n'
    elif name == 'bottleneck/matching.rs':
        begin, _, end = function(source, '    pub(super) fn new(size: usize)')
        body = source[begin:end]
        body = replace(body, 'let mut stack = Vec::new();', 'let mut stack = crate::lifecycle::empty("bn-stack");')
        fields = ('left', 'right', 'order', 'degrees', 'seen')
        for field in fields:
            value = 'NONE' if field in ('left', 'right') else '0'
            body = replace(body, f'{field}: filled(size, {value})?', f'{field}: crate::lifecycle::filled("bn-{field}", size, {value})?')
        source = source[:begin] + body + source[end:]
        source = drop_vectors(source, 'Workspace', (*fields, 'stack'), 'bn')
    elif name == 'bottleneck/geometry.rs':
        for field in ('left', 'right'):
            source = replace(source, f'{field}: filled(pair.size, NONE)?', f'{field}: crate::lifecycle::filled("geo-{field}", pair.size, NONE)?')
        replacements = {'levels': ('self.pair.size', 'NONE'), 'active': ('self.pair.second.points.len()', 'false'), 'remaining': ('self.tree.nodes.len()', '0')}
        for field, (length, value) in replacements.items():
            source = replace(source, f'self.{field} = filled({length}, {value})?;', f'self.{field} = crate::lifecycle::filled("geo-{field}", {length}, {value})?;')
        for field in ('queue', 'stack', 'free_diagonal'):
            source = replace(source, f'self.{field} = Vec::new();', f'self.{field} = crate::lifecycle::empty("geo-{field}");')
        source = drop_vectors(source, "Oracle<'a, 'p, 'q>", ('left', 'right', *replacements, 'queue', 'stack', 'free_diagonal'), 'geo', "<'a, 'p, 'q>")
    elif name == 'wasserstein/dense.rs':
        fields = {'u': ('short_len', '0.0'), 'v': ('long_len', '0.0'), 'minimum': ('long_len', '0.0'), 'matched': ('long_len', '0'), 'predecessor': ('long_len', '0'), 'used': ('long_len', 'false'), 'initially_matched': ('short_len', 'false')}
        for field, (length, value) in fields.items():
            source = replace(source, f'let mut {field} = buffer({length}, {value})?;', f'let mut {field} = crate::lifecycle::vector("dense-{field}", {length}, {value})?;')
    elif name == 'wasserstein/sparse.rs':
        for field, value in {'distance': 'f64::INFINITY', 'previous_node': 'usize::MAX', 'previous_edge': 'usize::MAX'}.items():
            source = replace(source, f'{field}: buffer(nodes, {value})?', f'{field}: crate::lifecycle::filled("sparse-{field}", nodes, {value})?')
        source = drop_vectors(source, 'Scratch', ('distance', 'previous_node', 'previous_edge'), 'sparse')
        # BinaryHeap and residual graph remain pair-local in this first candidate.
    return source


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def build(output):
    output.mkdir(parents=True, exist_ok=False)
    generated = output / 'src'
    shutil.copytree(REPO / 'src', generated)
    hashes = {}
    for path in sorted((generated / 'diagram_distances').rglob('*.rs')):
        name = path.relative_to(generated / 'diagram_distances').as_posix()
        path.write_text(instrument(path.read_text(), name), newline='\n')
        hashes[name] = sha256(path)
    shutil.copyfile(REPO / 'benches/distances/lifecycle.rs', output / 'worker.rs')
    shutil.copyfile(REPO / 'benches/distances/lifecycle_pool.rs', output / 'pool.rs')
    target = output / 'cargo'
    suffix = '.exe' if sys.platform == 'win32' else ''
    binary = output / ('worker' + suffix)
    commands = [ ['cargo', 'build', '--release', '--locked', '--offline', '--lib', '--target-dir', str(target)],
                 ['rustc', '--edition=2024', '-O', '--cfg', 'cocycle_distance_bench', str(output / 'worker.rs'),
                  '--extern', 'cocycle=' + str(target / 'release/libcocycle.rlib'), '-L', 'dependency=' + str(target / 'release/deps'), '-o', str(binary)] ]
    with (output / 'build.log').open('w') as log:
        for command in commands:
            log.write(json.dumps(command) + '\n'); log.flush()
            subprocess.run(command, cwd=REPO, check=True, stdout=log, stderr=subprocess.STDOUT, timeout=240)
    metadata = {'commands': commands, 'generated_sha256': hashes, 'worker_sha256': sha256(output / 'worker.rs'),
                'pool_sha256': sha256(output / 'pool.rs'), 'binary_sha256': sha256(binary),
                'generator_sha256': sha256(Path(__file__)),
                'original_sha256': {p.relative_to(REPO).as_posix(): sha256(p) for p in sorted((REPO/'src').rglob('*.rs'))},
                'rustc': subprocess.check_output(['rustc', '-Vv'], text=True)}
    save(output / 'build.json', metadata)
    subprocess.run([str(binary), '--selftest'], check=True, stdout=(output / 'selftest.json').open('w'), timeout=60)
    return binary


def workload(first, second, pattern, k):
    """K is an operation target; all-pairs uses the next full square matrix."""
    shift = lambda points, i: [(b + i / 256, d + i / 256) for b, d in points]
    if pattern == 'single':
        if k != 1:
            raise ValueError('single is K=1 only')
        return [first, second], [(0, 1, 0)], 0
    if pattern == 'same_pair':
        return [first, second], [(0, 1, 0)] * k, 0
    if pattern in ('one_to_many', 'many_to_one'):
        diagrams = [first] + [shift(second, i) for i in range(k)]
        ops = [(0, i + 1, 0) if pattern == 'one_to_many' else (i + 1, 0, 0) for i in range(k)]
        return diagrams, ops, 0
    if pattern == 'repeated_batch':
        diagrams, ops = [], []
        # Eight operations per batch; distinct query and targets between batches.
        for base in range(0, k, 8):
            query = len(diagrams)
            diagrams.append(shift(first, base))
            for j in range(min(8, k - base)):
                diagrams.append(shift(second, base + j))
                ops.append((query, len(diagrams) - 1, base // 8))
        return diagrams, ops, 0
    if pattern == 'all_pairs':
        n = max(2, math.ceil(math.sqrt(k)))
        diagrams = [shift(first if i % 2 == 0 else second, i) for i in range(n)]
        return diagrams, [(a, b, a) for a in range(n) for b in range(n)], 0
    raise ValueError('unknown pattern')


def write_input(path, diagrams, ops, fixed):
    if not diagrams or not ops or not 0 <= fixed < len(diagrams):
        raise ValueError('empty or invalid lifecycle input')
    data = bytearray(b'COCLIF1\0' + struct.pack('<QQQ', len(diagrams), len(ops), fixed))
    for points in diagrams:
        data.extend(struct.pack('<Q', len(points)))
        for b, d in points:
            if not math.isfinite(b) or not math.isfinite(d) or b >= d:
                raise ValueError('timed fixtures require finite positive lifetimes')
            data.extend(struct.pack('<dd', b, d))
    for a, b, batch in ops:
        if not 0 <= a < len(diagrams) or not 0 <= b < len(diagrams) or batch < 0:
            raise ValueError('invalid operation index')
        data.extend(struct.pack('<QQQ', a, b, batch))
    path.write_bytes(data)


def memory(pid):
    if sys.platform == 'win32':
        import ctypes
        from ctypes import wintypes
        class Counters(ctypes.Structure):
            _fields_ = [('cb', wintypes.DWORD), ('faults', wintypes.DWORD)] + [(n, ctypes.c_size_t) for n in ('peak_rss', 'rss', 'peak_paged', 'paged', 'peak_nonpaged', 'nonpaged', 'pagefile', 'peak_pagefile')]
        kernel = ctypes.WinDLL('kernel32', use_last_error=True)
        kernel.OpenProcess.restype = wintypes.HANDLE
        handle = kernel.OpenProcess(0x410, False, pid)
        if not handle:
            raise OSError('OpenProcess for worker memory failed')
        try:
            counters = Counters(); counters.cb = ctypes.sizeof(counters)
            psapi = ctypes.WinDLL('psapi', use_last_error=True)
            psapi.GetProcessMemoryInfo.argtypes = [wintypes.HANDLE, ctypes.POINTER(Counters), wintypes.DWORD]
            if not psapi.GetProcessMemoryInfo(handle, ctypes.byref(counters), counters.cb):
                raise OSError('GetProcessMemoryInfo failed')
            psapi_cpu = ctypes.WinDLL('kernel32', use_last_error=True)
            times = [wintypes.FILETIME() for _ in range(4)]
            psapi_cpu.GetProcessTimes.argtypes = [wintypes.HANDLE] + [ctypes.POINTER(wintypes.FILETIME)]*4
            if not psapi_cpu.GetProcessTimes(handle,*[ctypes.byref(t) for t in times]):
                raise OSError('GetProcessTimes failed')
            cpu_ns = sum((t.dwHighDateTime<<32)+t.dwLowDateTime for t in times[2:])*100
            return {'rss_bytes': counters.rss, 'peak_rss_bytes': counters.peak_rss, 'cpu_ns':cpu_ns}
        finally:
            kernel.CloseHandle.argtypes = [wintypes.HANDLE]; kernel.CloseHandle(handle)
    status = Path(f'/proc/{pid}/status').read_text()
    values = {line.split(':')[0]: int(line.split()[1]) * 1024 for line in status.splitlines() if line.startswith(('VmRSS:', 'VmHWM:'))}
    stat=Path(f'/proc/{pid}/stat').read_text().rpartition(')')[2].split()
    cpu_ns=(int(stat[11])+int(stat[12]))*1_000_000_000//os.sysconf('SC_CLK_TCK')
    return {'rss_bytes': values['VmRSS'], 'peak_rss_bytes': values['VmHWM'], 'cpu_ns':cpu_ns}


def validate(record, metric, variant, operations):
    if record.get('type') != 'result' or record.get('protocol') != PROTOCOL or record.get('metric') != metric or record.get('variant') != variant or record.get('operations') != operations:
        raise ValueError('wrong lifecycle result identity')
    for key in ('total_ns', 'prep_ns', 'loop_ns', 'prepared_bytes', 'retained_bytes', 'peak_retained_bytes', 'output_bytes', 'hits', 'growths'):
        if type(record.get(key)) is not int or record[key] < 0:
            raise ValueError(f'invalid lifecycle {key}')
    if len(record.get('phases_ns', [])) != 10 or any(type(x) is not int or x < 0 for x in record['phases_ns']):
        raise ValueError('invalid phase clocks')
    if not isinstance(record.get('checksum'), (int, float)) or not math.isfinite(record['checksum']) or record['checksum'] < 0 or record['total_ns'] < record['prep_ns']:
        raise ValueError('invalid lifecycle timing/checksum')
    return record


def invoke(binary, path, metric, variant, detailed, trace, operations, cpu, log, barrier=None):
    command = [str(binary), str(path), metric, variant, str(int(detailed)), str(int(trace)), '1']
    if cpu is not None and sys.platform == 'linux':
        command = ['taskset', '-c', str(cpu), *command]
    with log.open('w') as stderr:
        proc = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=stderr, text=True)
        deadline = threading.Timer(180, proc.kill); deadline.start()
        try:
            ready = json.loads(proc.stdout.readline())
            if ready != {'type': 'ready'}:
                raise ValueError('missing worker ready handshake')
            before = memory(proc.pid)
            if barrier is not None:
                barrier.wait(timeout=180)
            started = time.perf_counter_ns()
            proc.stdin.write('ack\n'); proc.stdin.flush()
            traces = []
            while True:
                record = json.loads(proc.stdout.readline())
                observed = memory(proc.pid)
                record['observed_memory'] = observed
                proc.stdin.write('ack\n'); proc.stdin.flush()
                if record.get('type') == 'trace':
                    if not trace or record['index'] != len(traces):
                        raise ValueError('wrong trace sequence')
                    traces.append(record)
                else:
                    record = validate(record, metric, variant, operations)
                    break
            finished = time.perf_counter_ns()
            proc.stdin.close()
            if proc.wait(timeout=10) != 0 or proc.stdout.read().strip():
                raise ValueError('worker failed or emitted trailing output')
            if trace and len(traces) != operations:
                raise ValueError('missing operation trace')
            record.update(command=command, trace=traces, started_ns=started, finished_ns=finished, wall_ns=finished-started, cpu_ns=observed['cpu_ns']-before['cpu_ns'])
            return record
        finally:
            deadline.cancel()
            if proc.poll() is None:
                proc.kill(); proc.wait()


def cells(args):
    source = list(fixtures(args.families, args.sizes, False))
    source = [f for f in source if f['split'] in args.splits]
    for fixture in source:
        for pattern in args.patterns:
            for k in ([1] if pattern == 'single' else args.ks):
                for metric in args.metrics:
                    yield fixture, pattern, k, metric


def run(args):
    # Formal attribution requires a clean harness and the immutable R0 kernel.
    if subprocess.check_output(['git','diff',R0,'--','Cargo.toml','Cargo.lock','src'],text=True).strip():
        raise ValueError('kernel differs from immutable R0')
    source_identity = {'dirty': subprocess.check_output(['git','status','--porcelain','--untracked-files=normal'],text=True).strip()}
    if source_identity['dirty'] and not args.exploratory:
        raise ValueError('commit harness before formal measurement, or use --exploratory')
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=False)
    if args.worker_dir:
        build_dir = args.worker_dir.resolve()
        binary = build_dir / ('worker.exe' if sys.platform == 'win32' else 'worker')
        metadata = json.loads((build_dir / 'build.json').read_text())
        if metadata['binary_sha256'] != sha256(binary) or metadata['worker_sha256'] != sha256(REPO / 'benches/distances/lifecycle.rs') or metadata['pool_sha256'] != sha256(REPO / 'benches/distances/lifecycle_pool.rs'):
            raise ValueError('worker source or binary changed since build')
        if metadata['generator_sha256'] != sha256(Path(__file__)) or any(sha256(REPO / p) != h for p,h in metadata['original_sha256'].items()):
            raise ValueError('generator or original source changed since build')
    else:
        binary = build(output / 'build'); build_dir = output / 'build'
    subprocess.run([str(binary), '--selftest'], check=True, capture_output=True, timeout=60)
    planned = list(cells(args))
    manifest = {'protocol': PROTOCOL, 'r0': R0, 'source_identity': source_identity, 'exploratory': args.exploratory, 'revision': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(), 'platform': platform.platform(), 'python': sys.version, 'binary_sha256': sha256(binary), 'build_dir': str(build_dir), 'samples': args.samples, 'variants': args.variants, 'workers': args.workers, 'cpu': args.cpu, 'detailed': args.detailed, 'trace': args.poison, 'created_utc': datetime.now(timezone.utc).isoformat(), 'cells': []}
    for f,p,k,m in planned:
        manifest['cells'].append({'fixture': f['name'], 'pattern': p, 'k': k, 'metric': m})
    save(output / 'manifest.json', manifest)
    records = []
    for cell_index, (fixture, pattern, k, metric) in enumerate(planned):
        name = f'{fixture["name"]}-{pattern}-k{k}-{metric}'
        cell = output / name; cell.mkdir()
        if args.poison:
            huge = next(f for f in fixtures([fixture['family']], [args.huge_size], False) if f['split']==fixture['split'])
            diagrams = [fixture['first'], fixture['second'], huge['first'], huge['second']]
            ops = [(0,1,i) for i in range(3)] + [(2,3,3)] + [(0,1,i) for i in range(4,12)]
            fixed = 0
        else:
            diagrams, ops, fixed = workload(fixture['first'], fixture['second'], pattern, k)
        write_input(cell / 'input.bin', diagrams, ops, fixed)
        if args.poison:
            # Reference solves run in a different process: do not pre-poison the
            # measured worker's allocator/RSS with the huge control operation.
            reference = json.loads(subprocess.check_output([str(binary),'--reference',str(cell/'input.bin'),metric],text=True,timeout=180))
            if len(reference.get('values',[])) != len(ops) or any(not math.isfinite(v) or v<0 for v in reference['values']):
                raise ValueError('invalid separate poison reference')
            save(cell/'reference.json',reference)
            (cell/'input.bin.expected.bin').write_bytes(struct.pack('<'+'d'*len(ops),*reference['values']))
        save(cell / 'input.json', {'diagrams': diagrams, 'ops': ops, 'fixed': fixed, 'sha256': sha256(cell / 'input.bin')})
        oracle_checksum = None
        if max(map(len,diagrams)) <= 8:
            oracle_checksum = sum(tiny_oracle(diagrams[a],diagrams[b],metric) for a,b,_ in ops)
        samples = []
        for entry in schedule(args.variants, args.samples, args.order_seed + cell_index):
            label = f'r{entry["round"]}-p{entry["position"]}-{entry["worker"]}'
            if args.workers == 1:
                result = invoke(binary, cell / 'input.bin', metric, entry['worker'], args.detailed, args.poison, len(ops), args.cpu, cell / (label + '.stderr'))
                if oracle_checksum is not None and not math.isclose(result['checksum'],oracle_checksum,rel_tol=1e-10,abs_tol=1e-10):
                    save(cell / (label + '-failed.json'), result)
                    raise ValueError('independent tiny oracle mismatch')
            else:
                barrier = threading.Barrier(args.workers)
                with ThreadPoolExecutor(max_workers=args.workers) as executor:
                    futures = [executor.submit(invoke, binary, cell / 'input.bin', metric, entry['worker'], False, False, len(ops), None if args.cpu is None else args.cpu+i, cell / f'{label}-w{i}.stderr', barrier) for i in range(args.workers)]
                    children = [f.result() for f in futures]
                ipc_wall = max(c['finished_ns'] for c in children)-min(c['started_ns'] for c in children)
                wall = max(c['compute_finished_unix_ns'] for c in children)-min(c['compute_started_unix_ns'] for c in children)
                result = {'children': children, 'aggregate_rss_bytes': sum(c['observed_memory']['rss_bytes'] for c in children), 'aggregate_peak_rss_bytes': sum(c['observed_memory']['peak_rss_bytes'] for c in children), 'wall_ns': wall, 'ipc_wall_ns':ipc_wall, 'throughput_per_s': args.workers * len(ops) * 1e9 / wall, 'cpu_fraction':sum(c['cpu_ns'] for c in children)/wall/args.workers}
            result.update(entry); samples.append(result); save(cell / (label + '.json'), result)
        summaries = {}
        for variant in args.variants:
            measured = [s for s in samples if s['worker'] == variant and not s['warmup']]
            summaries[variant] = {key: statistics.median([s[key] for s in measured]) for key in (('wall_ns','throughput_per_s','aggregate_rss_bytes','aggregate_peak_rss_bytes','cpu_fraction') if args.workers > 1 else ('total_ns','prep_ns','loop_ns','prepared_bytes','retained_bytes','output_bytes','hits','growths'))}
            children = [c for s in measured for c in s.get('children',[s])]
            latencies = sorted(c['total_ns']/len(ops) for c in children)
            summaries[variant].update(median_amortized_ns=statistics.median(latencies),p95_amortized_ns=latencies[math.ceil(.95*len(latencies))-1],median_worker_rss_bytes=statistics.median(c['observed_memory']['rss_bytes'] for c in children))
        record = {'name': name, 'family': fixture['family'], 'size': fixture['size'], 'split': fixture['split'], 'pattern': pattern, 'k': k, 'metric': metric, 'operations': len(ops), 'diagrams': len(diagrams), 'independent_checksum': oracle_checksum, 'summaries': summaries, 'status': 'passed'}
        records.append(record); save(cell / 'summary.json', record); save(output / 'results.json', records)
        print(f'{cell_index+1}/{len(planned)} {name}', flush=True)
    save(output / 'summary.json', {'status': 'passed', 'cells': len(records), 'measured_processes': len(records)*len(args.variants)*args.samples*args.workers, 'warmup_processes': len(records)*len(args.variants)*args.workers})



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output', type=Path, required=True)
    parser.add_argument('--build-only', action='store_true')
    parser.add_argument('--exploratory', action='store_true')
    parser.add_argument('--worker-dir', type=Path)
    parser.add_argument('--samples', type=int, default=12)
    parser.add_argument('--variants', nargs='+', choices=VARIANTS, default=list(VARIANTS))
    parser.add_argument('--families', nargs='+', choices=FAMILIES, default=['uniform','duplicates'])
    parser.add_argument('--sizes', nargs='+', type=int, default=[32])
    parser.add_argument('--splits', nargs='+', choices=['tuning','holdout'], default=['tuning','holdout'])
    parser.add_argument('--patterns', nargs='+', choices=PATTERNS, default=list(PATTERNS))
    parser.add_argument('--ks', nargs='+', type=int, default=list(KS))
    parser.add_argument('--metrics', nargs='+', choices=METRICS, default=list(METRICS))
    parser.add_argument('--cpu', type=int)
    parser.add_argument('--order-seed', type=int, default=20261001)
    parser.add_argument('--detailed', action='store_true')
    parser.add_argument('--poison', action='store_true')
    parser.add_argument('--huge-size', type=int, default=4096)
    parser.add_argument('--workers', type=int, choices=[1,2,4], default=1)
    args = parser.parse_args()
    if args.samples < 1 or any(k < 1 for k in args.ks) or any(n < 1 for n in args.sizes):
        parser.error('positive samples, sizes and K required')
    for name in ('variants','families','sizes','splits','patterns','ks','metrics'):
        if len(getattr(args,name)) != len(set(getattr(args,name))):
            parser.error(name + ' must have unique entries')
    if args.build_only:
        build(args.output.resolve())
    else:
        run(args)


if __name__ == '__main__':
    main()
