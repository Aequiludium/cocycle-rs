"""Separate Linux process/thread resource and throughput experiments on fixed R0."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from datetime import datetime, timezone
import json
import math
import os
from pathlib import Path
import platform
import random
import sys
import statistics
import subprocess
import threading
import time

from benchmark_distances import FAMILIES, fixtures, schedule
from distance_common import ROOT, agrees, save, sha256, tiny_oracle, tolerance, write_fixture

R0 = 'b2c3bd5eebf0c1193f30ea651012b8ae61b274c1'
PROTOCOL = 'cocycle-distance-resources-v2'


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
                                                       'tools/benchmark_distance_resources.py', 'tools/analyze_distance_resources.py',
                                                       'tools/benchmark_distances.py', 'tools/distance_common.py')}}
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


def group(binary, fixture, metric, variant, workers, mode, jobs, expected, error, cpus, trace, directory,
          *, nodes=None, cold=False, stagger_ms=0, sample_ms=2):
    processes = []
    logs = []
    stopped = threading.Event()
    samples = []
    sampling_errors = []
    specs = nodes or [dict(fixture=fixture,metric=metric,variant=variant,jobs=jobs,expected=expected,error=error)
                      for _ in range(workers if mode == 'process' else 1)]
    if nodes and mode != 'process':raise ValueError('heterogeneous plans require isolated processes')
    try:
        for index,spec in enumerate(specs):
            selected = [cpus[index % len(cpus)]] if mode == 'process' else cpus[:workers]
            command = ['taskset', '-c', ','.join(map(str, selected)), str(binary), str(spec['fixture']), spec['metric'],
                       spec['variant'], str(1 if mode == 'process' else workers), str(spec['jobs']), repr(spec['expected']),
                       repr(spec['error']), 'cold' if cold else 'warm']
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
                                'cpu_seconds': sum(v['cpu_seconds'] for v in values),'nodes':values})
            def poll():
                try:
                    while not stopped.wait(sample_ms/1000): sample()
                except Exception as exc: sampling_errors.append(repr(exc))
            samples.append({'t_ns': 0, 'nodes':baseline, **{key: sum(v[key] for v in baseline) for key in ('rss_bytes','pss_bytes','cpu_seconds')}})
            sampler = threading.Thread(target=poll, daemon=True) if trace else None
            if sampler: sampler.start()
            def receive(index):
                spec=specs[index]
                value=validate(json.loads(processes[index].stdout.readline()),spec['metric'],spec['variant'],
                               1 if mode=='process' else workers,spec['jobs'])
                return value,time.perf_counter_ns()-start
            pending=[readers.submit(receive,i) for i in range(len(processes))]
            dispatched=[]
            for index,p in enumerate(processes):
                remaining=index*stagger_ms/1000-(time.perf_counter_ns()-start)/1e9
                if remaining>0:time.sleep(remaining)
                p.stdin.write('go\n'); p.stdin.flush();dispatched.append(time.perf_counter_ns()-start)
            received=[f.result(timeout=120) for f in pending]
            results=[r[0] for r in received]
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
                'dispatch_ns':dispatched,'completion_ns':[r[1] for r in received],
                'cpu_by_child':[b['cpu_seconds']-a['cpu_seconds'] for a,b in zip(baseline,samples[-1]['nodes'])],
                'cold':cold,'stagger_ms':stagger_ms,
                'dram_bytes': None, 'llc_misses': None}
    finally:
        stopped.set()
        for p in processes:
            if p.poll() is None: p.kill(); p.wait()
        for log in logs: log.close()


def calibrate(binary,fixture,metric,expected,error,cpus,directory,target_ms):
    """Choose observation duration from baseline only; never select a route."""
    path=directory/'calibration';path.mkdir()
    value=group(binary,fixture,metric,'baseline',1,'process',4,expected,error,cpus,False,path)
    jobs=max(1,min(4096,math.ceil(target_ms*1e6/max(1,value['median_job_ns']))))
    save(path/'result.json',dict(value,jobs_selected=jobs,target_ms=target_ms))
    return jobs


def run_plan(args,binary,cpus,output,manifest):
    """A finite preregistered list of heterogeneous process groups."""
    plans=json.loads(args.plan.read_text())
    if not plans or len({p['name'] for p in plans})!=len(plans):raise ValueError('empty/duplicate plan')
    cases={}
    inputs={}
    for plan in plans:
        name=plan['name']
        if not name or any(c not in 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789-_' for c in name):
            raise ValueError('unsafe plan name')
        if not plan['nodes']:raise ValueError('empty worker group')
        specs=[]
        for node in plan['nodes']:
            key=(node['family'],node['size'],node['split'],node['metric'])
            if key[0] not in FAMILIES or key[1]<1 or key[2] not in ('tuning','holdout') or key[3] not in ('bottleneck','w1','w2'):
                raise ValueError('invalid plan fixture')
            if node['variant'] not in ('baseline','adaptive_arena','vectors','arena'):raise ValueError('unknown route')
            if key not in inputs:
                fixture=next(f for f in fixtures([key[0]],[key[1]],False) if f['split']==key[2])
                path=output/'fixtures'/('-'.join(map(str,key)));path.mkdir(parents=True)
                write_fixture(path/'input.bin',fixture['first'],fixture['second'])
                expected=float(subprocess.check_output([str(binary),'--reference',str(path/'input.bin'),key[3]],text=True,timeout=120))
                if not math.isfinite(expected):raise ValueError('non-finite reference')
                if key[1]<=8 and not agrees(expected,tiny_oracle(fixture['first'],fixture['second'],key[3]),fixture,key[3]):
                    raise ValueError('independent tiny oracle mismatch')
                error=tolerance(fixture,key[3],expected)
                jobs=calibrate(binary,path/'input.bin',key[3],expected,error,cpus,path,args.target_ms) if args.target_ms else args.jobs
                inputs[key]=dict(fixture=path/'input.bin',metric=key[3],expected=expected,error=error,jobs=jobs)
                save(path/'input.json',dict(fixture=fixture,sha256=sha256(path/'input.bin'),reference=expected,jobs=jobs))
            spec=dict(inputs[key],variant=node['variant'],jobs=node.get('jobs',inputs[key]['jobs']))
            if spec['jobs']<1:raise ValueError('non-positive plan job count')
            specs.append(spec)
        cases[name]=(plan,specs)
        (output/name).mkdir()
    manifest.update(plan=plans,cells=list(cases));save(output/'manifest.json',manifest)
    rows={name:[] for name in cases}
    # Balance route positions within equivalent fixture/job/placement blocks;
    # heterogeneous configurations must not drift apart across an entire run.
    blocks={}
    for name,(plan,specs) in cases.items():
        key=tuple(sorted((str(s['fixture']),s['metric'],s['jobs']) for s in specs))+(plan.get('stagger_ms',args.stagger_ms),)
        blocks.setdefault(key,[]).append(name)
    scheduled=[schedule(names,args.samples,args.order_seed+i) for i,names in enumerate(blocks.values())]
    entries=[]
    for round_index in range(args.samples+1):
        order=list(range(len(scheduled)));random.Random(args.order_seed+round_index).shuffle(order)
        position=0
        for index in order:
            for entry in scheduled[index]:
                if entry['round']==round_index:
                    entries.append(dict(entry,block=index,global_position=position));position+=1
    manifest['plan_position_contract']='balanced equivalent fixture/job blocks; seeded shuffled block order'
    save(output/'manifest.json',manifest)
    for entry in entries:
        name=entry['worker'];plan,specs=cases[name]
        directory=output/name/f'r{entry["round"]}-p{entry["position"]}';directory.mkdir()
        value=group(binary,None,None,None,len(specs),'process',0,0,0,cpus,args.trace,directory,
                    nodes=specs,cold=args.cold,stagger_ms=plan.get('stagger_ms',args.stagger_ms),sample_ms=args.sample_ms)
        value.update(entry);rows[name].append(value);save(directory/'result.json',value)
        print(f'{entry["round"]} {name}',flush=True)
    records=[]
    for name,(plan,specs) in cases.items():
        measured=[s for s in rows[name] if not s['warmup']]
        summaries={key:statistics.median(s[key] for s in measured) for key in
                   ('throughput_per_s','median_job_ns','p95_job_ns','endpoint_pss_bytes','cpu_core_seconds','cpu_fraction')}
        record=dict(name=name,nodes=plan['nodes'],workers=len(specs),mode='process',status='passed',summaries=summaries)
        records.append(record);save(output/name/'summary.json',record)
    save(output/'results.json',records)
    save(output/'summary.json',dict(status='passed',cells=len(records),measured_groups=len(records)*args.samples,
                                    warmup_groups=len(records),trace=args.trace))


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
    if args.cpus:
        if not set(args.cpus)<=set(cpus):raise ValueError('requested CPU is outside allowed affinity')
        cpus=args.cpus
    manifest = {'protocol':PROTOCOL,'kernel':R0,'harness':subprocess.check_output(['git','rev-parse','HEAD'],text=True).strip(),
                'dirty':dirty,'exploratory':args.exploratory,'trace':args.trace,'samples':args.samples,'jobs':args.jobs,
                'workers':args.workers,'modes':args.modes,'variants':args.variants,'cpus':cpus,'platform':platform.platform(),
                'created_utc':datetime.now(timezone.utc).isoformat(),'build':str(directory),'binary_sha256':sha256(binary),
                'memory_measurement':'sum of sampled process PSS, proportional shared-page proxy; not cgroup',
                'dram_bytes':None,'llc_misses':None,'hardware_note':'no usable perf/IMC counter in this WSL environment',
                'cgroup_writable':os.access('/sys/fs/cgroup',os.W_OK),'cells':[]}
    manifest.update(command=sys.argv,cold=args.cold,sample_ms=args.sample_ms,stagger_ms=args.stagger_ms,
                    target_ms=args.target_ms,order_seed=args.order_seed,placement='round-robin individual CPU pinning; threads share CPU mask')
    save(output/'manifest.json',manifest)
    if args.plan:
        run_plan(args,binary,cpus,output,manifest)
        return
    records=[];calibrations={}
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
                    key=(fixture['name'],metric)
                    if key not in calibrations:
                        calibrations[key]=calibrate(binary,cell/'input.bin',metric,expected,tolerance(fixture,metric,expected),cpus,cell,args.target_ms) if args.target_ms else args.jobs
                    jobs=calibrations[key]
                    manifest['cells'].append(name);save(output/'manifest.json',manifest)
                    samples=[]
                    for entry in schedule(args.variants,args.samples,args.order_seed):
                        sample_dir=cell/f'r{entry["round"]}-p{entry["position"]}-{entry["worker"]}';sample_dir.mkdir()
                        value=group(binary,cell/'input.bin',metric,entry['worker'],workers,mode,jobs,expected,
                                    tolerance(fixture,metric,expected),cpus,args.trace,sample_dir,cold=args.cold,
                                    stagger_ms=args.stagger_ms,sample_ms=args.sample_ms)
                        value.update(entry);samples.append(value);save(sample_dir/'result.json',value)
                    summaries={v:{key:statistics.median(s[key] for s in samples if s['worker']==v and not s['warmup'])
                                  for key in ('throughput_per_s','median_job_ns','p95_job_ns','endpoint_pss_bytes','cpu_core_seconds','cpu_fraction')}
                               for v in args.variants}
                    record={'name':name,'family':fixture['family'],'size':fixture['size'],'split':fixture['split'],
                            'metric':metric,'mode':mode,'workers':workers,'jobs_per_worker':jobs,'summaries':summaries,'status':'passed'}
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
    p.add_argument('--families',choices=FAMILIES,nargs='+',default=['uniform','sparse'])
    p.add_argument('--sizes',type=int,nargs='+',default=[128]);p.add_argument('--splits',nargs='+',choices=['tuning','holdout'],default=['tuning','holdout'])
    p.add_argument('--metrics',choices=['bottleneck','w1','w2'],nargs='+',default=['w1','w2'])
    p.add_argument('--variants',choices=['baseline','vectors','arena','adaptive_arena'],nargs='+',default=['baseline','vectors','arena'])
    p.add_argument('--trace',action='store_true');p.add_argument('--exploratory',action='store_true')
    p.add_argument('--cold',action='store_true',help='skip in-process solver warmup for single-operation traces')
    p.add_argument('--sample-ms',type=float,default=2)
    p.add_argument('--stagger-ms',type=float,default=0)
    p.add_argument('--target-ms',type=float,default=0,help='calibrate jobs from baseline only, outside measured groups')
    p.add_argument('--order-seed',type=int,default=20261004)
    p.add_argument('--cpus',type=int,nargs='+')
    p.add_argument('--plan',type=Path,help='JSON list of named heterogeneous process groups')
    args=p.parse_args()
    if min(args.workers)<1 or args.samples<1 or args.jobs<1:p.error('positive worker/sample/job counts required')
    if args.sample_ms<=0 or args.stagger_ms<0 or args.target_ms<0:p.error('invalid observation durations')
    if 'bottleneck' in args.metrics and args.variants!=['baseline']:p.error('bottleneck uses baseline only')
    try:run(args)
    except Exception as exc:
        if args.output.exists():save(args.output/'failure.json',{'status':'failed','error':repr(exc),'command':sys.argv})
        raise


if __name__=='__main__':main()
