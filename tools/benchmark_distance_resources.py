"""Separate Linux process/thread resource and throughput experiments on fixed R0."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import statistics
import subprocess
import threading
import time

from benchmark_distances import fixtures, schedule
from distance_common import ROOT, agrees, save, sha256, tiny_oracle, tolerance, write_fixture

R0 = 'b2c3bd5eebf0c1193f30ea651012b8ae61b274c1'
PROTOCOL = 'cocycle-distance-resources-v1'


def build(directory):
    directory.mkdir(parents=True, exist_ok=False)
    originals = {}
    for path in (ROOT / 'src').rglob('*.rs'):
        relative = path.relative_to(ROOT)
        target = directory / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(path.read_bytes())
        originals[relative.as_posix()] = sha256(path)
    source = (ROOT / 'benches/distances/cocycle.rs').read_text()
    marker = 'fn run() -> std::result::Result<(), Box<dyn std::error::Error>> {'
    if source.count(marker) != 1:
        raise ValueError('ordinary worker entry point changed')
    worker = directory / 'benches/distances/cocycle.rs'
    worker.parent.mkdir(parents=True)
    worker.write_text(source.split(marker)[0] + (ROOT / 'benches/distances/resources.rs').read_text())
    cargo = directory / 'cargo'
    commands = [['cargo', 'build', '--release', '--locked', '--offline', '--lib', '--target-dir', str(cargo)],
                ['rustc', '--edition=2024', '-O', '--cfg', 'cocycle_distance_bench', str(worker),
                 '--extern', 'cocycle=' + str(cargo / 'release/libcocycle.rlib'), '-L',
                 'dependency=' + str(cargo / 'release/deps'), '-o', str(directory / 'worker')]]
    with (directory / 'build.log').open('w') as log:
        for command in commands:
            subprocess.run(command, cwd=ROOT, check=True, stdout=log, stderr=subprocess.STDOUT, timeout=240)
    metadata = {'originals': originals, 'commands': commands, 'binary_sha256': sha256(directory / 'worker'),
                'generated_worker_sha256': sha256(worker), 'rustc': subprocess.check_output(['rustc', '-Vv'], text=True),
                'owners': {p: sha256(ROOT / p) for p in ('benches/distances/cocycle.rs', 'benches/distances/resources.rs',
                                                       'tools/benchmark_distance_resources.py', 'tools/benchmark_distances.py', 'tools/distance_common.py')}}
    save(directory / 'build.json', metadata)
    return directory / 'worker'


def observe(pid):
    stat = Path(f'/proc/{pid}/stat').read_text().rsplit(')', 1)[1].split()
    ticks = sum(int(x) for x in stat[11:13])
    fields = {}
    for line in Path(f'/proc/{pid}/smaps_rollup').read_text().splitlines():
        if ':' in line:
            name, value = line.split(':', 1)
            if name in ('Rss', 'Pss'):
                fields[name] = int(value.split()[0]) * 1024
    return {'cpu_seconds': ticks / os.sysconf('SC_CLK_TCK'), 'rss_bytes': fields['Rss'], 'pss_bytes': fields['Pss']}


def integrate(samples, key):
    if len(samples) < 2 or any(b['t_ns'] < a['t_ns'] for a,b in zip(samples,samples[1:])):
        raise ValueError('invalid monotonic resource trace')
    return sum((b['t_ns'] - a['t_ns']) / 1e9 * (a[key] + b[key]) / 2 for a,b in zip(samples,samples[1:]))


def validate(result, metric, variant, threads, jobs):
    if (result.get('type'), result.get('protocol'), result.get('metric'), result.get('variant')) != ('resource-result', PROTOCOL, metric, variant):
        raise ValueError('worker identity mismatch')
    values = result.get('durations_ns', [])
    if result.get('threads') != threads or result.get('jobs_per_thread') != jobs or len(values) != threads or any(len(v) != jobs for v in values):
        raise ValueError('worker job count mismatch')
    if any(type(t) != int or t < 0 for v in values for t in v) or result['native_group_ns'] <= 0 or not math.isfinite(result['checksum']):
        raise ValueError('invalid worker timing or scalar')
    return result


def group(binary, fixture, metric, variant, workers, mode, jobs, expected, error, cpus, trace, directory):
    processes = []
    logs = []
    stopped = threading.Event()
    samples = []
    sampling_errors = []
    try:
        for index in range(workers if mode == 'process' else 1):
            selected = [cpus[index]] if mode == 'process' else cpus[:workers]
            command = ['taskset', '-c', ','.join(map(str, selected)), str(binary), str(fixture), metric,
                       variant, str(1 if mode == 'process' else workers), str(jobs), repr(expected), repr(error)]
            log = (directory / f'child-{index}.stderr').open('w'); logs.append(log)
            processes.append(subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                                              stderr=log, text=True))
        readers = ThreadPoolExecutor(max_workers=len(processes))
        try:
            ready = [readers.submit(p.stdout.readline) for p in processes]
            if any(json.loads(f.result(timeout=120)).get('type') != 'ready' for f in ready):
                raise ValueError('missing ready handshake')
            baseline = [observe(p.pid) for p in processes]
            start = time.perf_counter_ns()
            def sample():
                values = [observe(p.pid) for p in processes]
                samples.append({'t_ns': time.perf_counter_ns() - start,
                                'rss_bytes': sum(v['rss_bytes'] for v in values),
                                'pss_bytes': sum(v['pss_bytes'] for v in values),
                                'cpu_seconds': sum(v['cpu_seconds'] for v in values)})
            def poll():
                try:
                    while not stopped.wait(.002): sample()
                except Exception as exc: sampling_errors.append(repr(exc))
            samples.append({'t_ns': 0, **{key: sum(v[key] for v in baseline) for key in ('rss_bytes','pss_bytes','cpu_seconds')}})
            sampler = threading.Thread(target=poll, daemon=True) if trace else None
            if sampler: sampler.start()
            for p in processes: p.stdin.write('go\n'); p.stdin.flush()
            pending = [readers.submit(p.stdout.readline) for p in processes]
            results = [validate(json.loads(f.result(timeout=120)), metric, variant, 1 if mode == 'process' else workers, jobs) for f in pending]
            elapsed = time.perf_counter_ns() - start
            stopped.set()
            if sampler: sampler.join(timeout=10)
            sample()
            if sampling_errors: raise ValueError(sampling_errors)
            for p in processes: p.stdin.write('ack\n'); p.stdin.flush()
            for p in processes:
                if p.wait(timeout=10) != 0 or p.stdout.read().strip(): raise ValueError('worker exit/trailing output')
        finally:
            # Kill our children before waiting for pipe readers on any failure.
            # A worker waiting at the ready/ack barrier must not block cleanup.
            for p in processes:
                if p.poll() is None: p.kill()
            readers.shutdown(wait=True)
        latencies = sorted(t for r in results for thread in r['durations_ns'] for t in thread)
        cpu = samples[-1]['cpu_seconds'] - samples[0]['cpu_seconds']
        return {'children': results, 'group_wall_ns': elapsed, 'throughput_per_s': len(latencies)*1e9/elapsed,
                'median_job_ns': statistics.median(latencies), 'p95_job_ns': latencies[math.ceil(.95*len(latencies))-1],
                'cpu_core_seconds': cpu, 'cpu_fraction': cpu/(elapsed/1e9)/workers,
                'endpoint_pss_bytes': samples[-1]['pss_bytes'], 'endpoint_rss_bytes': samples[-1]['rss_bytes'],
                'trace': samples if trace else [], 'pss_integral_byte_seconds': integrate(samples,'pss_bytes') if trace else None,
                'sampled_peak_pss_bytes': max(s['pss_bytes'] for s in samples) if trace else None,
                'dram_bytes': None, 'llc_misses': None}
    finally:
        stopped.set()
        for p in processes:
            if p.poll() is None: p.kill(); p.wait()
        for log in logs: log.close()


def run(args):
    if platform.system() != 'Linux': raise ValueError('resource protocol requires Linux procfs/taskset')
    if subprocess.check_output(['git','diff',R0,'--','src','Cargo.toml','Cargo.lock'],text=True).strip():
        raise ValueError('production kernel differs from R0')
    dirty = subprocess.check_output(['git','status','--porcelain'],text=True).strip()
    if dirty and not args.exploratory: raise ValueError('freeze harness before formal measurement')
    output = args.output.resolve(); output.mkdir(parents=True, exist_ok=False)
    directory = args.worker_dir.resolve() if args.worker_dir else output/'build'
    binary = directory/'worker' if args.worker_dir else build(directory)
    metadata = json.loads((directory/'build.json').read_text())
    if sha256(binary) != metadata['binary_sha256'] or any(sha256(ROOT/p) != h for p,h in {**metadata['originals'],**metadata['owners']}.items()):
        raise ValueError('built worker or source fingerprint changed')
    cpus = sorted(os.sched_getaffinity(0))
    if max(args.workers)>len(cpus): raise ValueError('worker count exceeds allowed CPUs')
    manifest = {'protocol':PROTOCOL,'kernel':R0,'harness':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
                'dirty':dirty,'exploratory':args.exploratory,'trace':args.trace,'samples':args.samples,'jobs':args.jobs,
                'workers':args.workers,'modes':args.modes,'variants':args.variants,'cpus':cpus,'platform':platform.platform(),
                'created_utc':datetime.now(timezone.utc).isoformat(),'build':str(directory),'binary_sha256':sha256(binary),
                'memory_measurement':'sum of sampled process PSS, proportional shared-page proxy; not cgroup',
                'dram_bytes':None,'llc_misses':None,'hardware_note':'no usable perf/IMC counter in this WSL environment',
                'cgroup_writable':os.access('/sys/fs/cgroup',os.W_OK),'cells':[]}
    save(output/'manifest.json',manifest)
    records=[]
    for fixture in fixtures(args.families,args.sizes,False):
        if fixture['split'] not in args.splits:continue
        for metric in args.metrics:
            for mode in args.modes:
                for workers in args.workers:
                    name=f'{fixture["name"]}-{metric}-{mode}-n{workers}'
                    cell=output/name;cell.mkdir()
                    write_fixture(cell/'input.bin',fixture['first'],fixture['second'])
                    expected=float(subprocess.check_output([str(binary),'--reference',str(cell/'input.bin'),metric],text=True,timeout=120))
                    if not math.isfinite(expected):raise ValueError('non-finite reference')
                    if max(len(fixture['first']),len(fixture['second']))<=8 and not agrees(expected,tiny_oracle(fixture['first'],fixture['second'],metric),fixture,metric):
                        raise ValueError('independent tiny oracle mismatch')
                    save(cell/'input.json',{'fixture':fixture,'sha256':sha256(cell/'input.bin'),'reference':expected})
                    manifest['cells'].append(name);save(output/'manifest.json',manifest)
                    samples=[]
                    for entry in schedule(args.variants,args.samples,20261004):
                        sample_dir=cell/f'r{entry["round"]}-p{entry["position"]}-{entry["worker"]}';sample_dir.mkdir()
                        value=group(binary,cell/'input.bin',metric,entry['worker'],workers,mode,args.jobs,expected,
                                    tolerance(fixture,metric,expected),cpus,args.trace,sample_dir)
                        value.update(entry);samples.append(value);save(sample_dir/'result.json',value)
                    summaries={v:{key:statistics.median(s[key] for s in samples if s['worker']==v and not s['warmup'])
                                  for key in ('throughput_per_s','median_job_ns','p95_job_ns','endpoint_pss_bytes','cpu_core_seconds','cpu_fraction')}
                               for v in args.variants}
                    record={'name':name,'family':fixture['family'],'size':fixture['size'],'split':fixture['split'],
                            'metric':metric,'mode':mode,'workers':workers,'jobs_per_worker':args.jobs,'summaries':summaries,'status':'passed'}
                    records.append(record);save(cell/'summary.json',record);save(output/'results.json',records)
                    print(f'{len(records)} {name}',flush=True)
    save(output/'summary.json',{'status':'passed','cells':len(records),'measured_groups':len(records)*len(args.variants)*args.samples,
                              'warmup_groups':len(records)*len(args.variants),'trace':args.trace})


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output',type=Path,required=True);p.add_argument('--worker-dir',type=Path)
    p.add_argument('--samples',type=int,default=12);p.add_argument('--jobs',type=int,default=64)
    p.add_argument('--workers',type=int,nargs='+',default=[1,2,4])
    p.add_argument('--modes',choices=['process','threads'],nargs='+',default=['process'])
    p.add_argument('--families',choices=['uniform','dense','sparse','duplicates'],nargs='+',default=['uniform','sparse'])
    p.add_argument('--sizes',type=int,nargs='+',default=[128]);p.add_argument('--splits',nargs='+',choices=['tuning','holdout'],default=['tuning','holdout'])
    p.add_argument('--metrics',choices=['w1','w2'],nargs='+',default=['w1','w2'])
    p.add_argument('--variants',choices=['baseline','vectors','arena','adaptive_arena'],nargs='+',default=['baseline','vectors','arena'])
    p.add_argument('--trace',action='store_true');p.add_argument('--exploratory',action='store_true')
    args=p.parse_args()
    if min(args.workers)<1 or args.samples<1 or args.jobs<1:p.error('positive worker/sample/job counts required')
    run(args)


if __name__=='__main__':main()
