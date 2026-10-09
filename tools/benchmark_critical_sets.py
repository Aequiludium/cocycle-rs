"""Native sparse/dense/Oineus critical-set timing and repository persistence context."""
import argparse
from datetime import datetime, timezone
import hashlib
import itertools
import json
import os
from pathlib import Path
import random
import shutil
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
PROTOCOL = 'critical-sets-native-v1'
DENSE_SHA = '167026e7ee8dbaafe8abaf5e325149af6e93f660'
DENSE_BLOB = '53bd137d16e122c2ad00a9859188ff317a362472'
OINEUS_SHA = 'e52814a1ffb5b8a81e71ff1e93b4c14194673f0f'


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def git(*args, binary=False):
    command=['git','-C',str(ROOT)]
    entry=ROOT/'.git'
    if sys.platform=='linux' and entry.is_file():
        location=entry.read_text().strip().removeprefix('gitdir: ').replace('\\','/')
        if len(location)>2 and location[1]==':':
            location='/mnt/'+location[0].lower()+location[2:]
            command=['git',f'--git-dir={location}',f'--work-tree={ROOT}']
    output = subprocess.check_output([*command,*args], text=not binary)
    return output if binary else output.strip()


def dense_source(path=None):
    data = Path(path).read_bytes() if path else git('show',f'{DENSE_SHA}:benches/optimization/big_steps.rs',binary=True)
    blob = hashlib.sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()
    if blob != DENSE_BLOB:
        raise ValueError('dense baseline does not match the pinned Git blob')
    return data


def canonical(cells):
    return sorted(cells, key=lambda c: (c[1], len(c[0]), tuple(-v for v in reversed(c[0]))))


def pair_oracle(cells):
    """Independent pairing-only set XOR, using vertex faces rather than Rust state."""
    lookup = {tuple(v): i for i, (v, _) in enumerate(cells)}
    columns = []
    owners = {}
    for j, (vertices, _) in enumerate(cells):
        column = set() if len(vertices) == 1 else {lookup[tuple(vertices[:k] + vertices[k+1:])] for k in range(len(vertices))}
        assert all(i < j for i in column)
        while column and max(column) in owners:
            column ^= columns[owners[max(column)]]
        if column:
            owners[max(column)] = j
        columns.append(column)
    return sorted(owners.items()), [i for i, c in enumerate(columns) if not c and i not in owners]


def diagram(cells, pairs, essential, q):
    output = [[len(cells[b][0])-1, cells[b][1], cells[d][1]] for b, d in pairs
              if len(cells[b][0])-1 <= q and cells[b][1] < cells[d][1]]
    output += [[len(cells[b][0])-1, cells[b][1], None] for b in essential if len(cells[b][0])-1 <= q]
    return sorted(output, key=lambda x: (x[0], x[1], math_end(x[2])))


def math_end(value):
    return float('inf') if value is None else value


def fixtures(quick=False):
    cases = []
    rng = random.Random(220316748)
    sizes = ((6, 12), (32, 128), (10, 18), (8, 12)) if quick else ((8, 16), (64, 256), (12, 20), (8, 14))
    for side in sizes[0]:
        edges = set()
        triangles = []
        for x in range(side-1):
            for y in range(side-1):
                a=x*side+y; b=a+1; c=a+side; d=c+1
                triangles += [[a,b,d], [a,c,d]]
        for v in triangles:
            edges.update(itertools.combinations(v, 2))
        weights = {edge: 1.+rng.randrange(1, 10000)/10000 for edge in sorted(edges)}
        cells = [([i], 0.) for i in range(side*side)]
        cells += [(list(edge), f) for edge, f in weights.items()]
        cells += [(v, max(weights[edge] for edge in itertools.combinations(v,2))) for v in triangles]
        cases.append({'name':f'grid-{side}', 'cells':canonical(cells), 'q':1, 'flag':True})
    for n in sizes[1]:
        triangles = [[0,i,i+1] for i in range(1,n-1)] + [[0,1,n-1]]
        edges = set(itertools.chain.from_iterable(itertools.combinations(v,2) for v in triangles))
        cells = [([i],0.) for i in range(n)]
        cells += [(list(e),1.+rng.randrange(1,10000)/10000) for e in sorted(edges)]
        cells += [(v,3.+rng.randrange(1,10000)/10000) for v in triangles]
        cases.append({'name':f'delayed-fan-{n}', 'cells':canonical(cells), 'q':1, 'flag':False})
    for n in sizes[2]:
        weights = {e:1.+rng.randrange(1,10000)/10000 for e in itertools.combinations(range(n),2)}
        cells = [([i],0.) for i in range(n)] + [(list(e),f) for e,f in weights.items()]
        cells += [(list(v),max(weights[e] for e in itertools.combinations(v,2))) for v in itertools.combinations(range(n),3)]
        cases.append({'name':f'clique-{n}', 'cells':canonical(cells), 'q':1, 'flag':True})
    for n in sizes[3]:
        cells = [([i],0.) for i in range(n)]
        for p in (1,2,3):
            cells += [(list(v),2.*p-1.+rng.randrange(1,10000)/10000) for v in itertools.combinations(range(n),p+1)]
        cases.append({'name':f'tetra-skeleton-{n}', 'cells':canonical(cells), 'q':2, 'flag':False})
    return cases


def workloads(case):
    cells = case['cells']
    pairs, _ = pair_oracle(cells)
    pairs = [(b,d) for b,d in pairs if len(cells[b][0])-1 <= case['q']]
    assert pairs
    span = cells[-1][1] - cells[0][1]
    specs = [('primal-U1', 1, 0.037, 1024), ('mixed-Q8', 8, 0.037, 128), ('mixed-Q64', 64, 0.371, 16)]
    for name, count, fraction, repeats in specs:
        proposals = []
        for i in range(count):
            b,d = pairs[(i*17) % len(pairs)]
            # One death increase, then four finite-endpoint directions.
            id, sign = (d,1) if count==1 else ((d,1),(d,-1),(b,1),(b,-1))[i%4]
            proposals.append((id,cells[id][1]+sign*fraction*span))
        yield name, proposals, repeats


def fixture_text(case, proposals, repeats):
    cells = case['cells']
    scale = cells[len(cells)//2][1]
    lines = [f"{len(cells)} {case['q']} {len(proposals)} {repeats} {scale:.17g}"]
    lines += [f"{len(v)} {' '.join(map(str,v))} {f:.17g}" for v,f in cells]
    lines += [f'{i} {t:.17g}' for i,t in proposals]
    return '\n'.join(lines)+'\n'


def build(args, output):
    native = output/'workers'
    native.mkdir()
    original = dense_source(args.dense_source)
    (native/'dense-original.rs').write_bytes(original)
    dense = original.decode('utf-8')
    assert dense.count('.then(cells[i].vertices.cmp(&cells[j].vertices))') == 1
    dense = dense.replace('.then(cells[i].vertices.cmp(&cells[j].vertices))', '.then(cells[j].vertices.iter().rev().cmp(cells[i].vertices.iter().rev()))')
    dense = '#![allow(dead_code, unused_imports)]\n'+dense
    dense += '''
pub struct Input { cells: Vec<Cell> }
impl Input {
    pub fn new(source: &cocycle::complex::SimplicialComplex) -> Self {
        Self { cells: source.simplices().iter().map(|s|Cell { vertices:s.vertices().to_vec(), value:s.value() }).collect() }
    }
}
pub struct State<'a> { input:&'a Input, analysis:Analysis, births:Vec<bool> }
impl<'a> State<'a> {
    pub fn new(input:&'a Input)->Self {
        let analysis=analyze(&input.cells);let mut births=vec![false;input.cells.len()];
        for &(b,_) in &analysis.pairs {births[b]=true;}
        Self {input,analysis,births}
    }
    pub fn pairs(&self)->Vec<(usize,usize)> { self.analysis.pairs.clone() }
    pub fn essential(&self)->Vec<usize> { self.analysis.essential.clone() }
    pub fn targets(&self,proposals:&[(usize,f64)])->Vec<(usize,f64)> {
        let mut result:BTreeMap<usize,f64>=BTreeMap::new();
        for &(id,target) in proposals {
            if target==self.input.cells[id].value {continue;}
            let birth=self.births[id];
            for i in critical(&self.input.cells,&self.analysis,id,target,birth) {
                let value=self.input.cells[i].value;
                if let Some(previous)=result.get_mut(&i) {if (target-value).abs()>(*previous-value).abs(){*previous=target;}}
                else {result.insert(i,target);}
            }
        }
        result.into_iter().collect()
    }
}
'''
    (native/'dense.rs').write_text(dense)
    for name in ('column','reduction'):
        text = (ROOT/f'src/persistence/reference/{name}.rs').read_text()
        text = text.replace('pub(in crate::persistence)','pub(super)')
        (native/f'{name}.rs').write_text(text)
    (native/'reference.rs').write_text('''#![allow(dead_code)]
mod column;
mod reduction;
use cocycle::complex::{SimplexId,SimplicialComplex};
use crate::Result;
trait FilteredBoundary {
    fn len(&self)->usize;
    fn dimension(&self,index:usize)->usize;
    fn value(&self,index:usize)->f64;
    fn write_boundary(&self,index:usize,output:&mut Vec<usize>)->Result<()>;
}
struct Input<'a>(&'a SimplicialComplex,&'a [SimplexId]);
impl FilteredBoundary for Input<'_> {
    fn len(&self)->usize {self.0.len()}
    fn dimension(&self,index:usize)->usize {self.0.simplices()[index].dimension()}
    fn value(&self,index:usize)->f64 {self.0.simplices()[index].value()}
    fn write_boundary(&self,index:usize,output:&mut Vec<usize>)->Result<()> {
        let id=self.1[index];
        output.clear();output.extend(self.0.boundary(id).unwrap().iter().map(|t|t.face.index()));output.sort_unstable();Ok(())
    }
}
pub fn run(source:&SimplicialComplex,ids:&[SimplexId])->Result<(Vec<(usize,usize)>,Vec<usize>)> {
    let r=reduction::reduce(&Input(source,ids))?;Ok((r.pairs,r.unpaired))
}
''')
    # A generated module's default submodule lookup is beside reference.rs.
    reference_dir = native/'reference'
    reference_dir.mkdir()
    for name in ('column','reduction'):
        shutil.move(native/f'{name}.rs',reference_dir/f'{name}.rs')
    shutil.copy2(ROOT/'benches/optimization/cocycle.rs',native/'cocycle.rs')
    rust_flags = ['--edition','2024','-C','opt-level=3']
    commands = [
        [str(args.rustc),'--crate-name','cocycle','--crate-type','rlib',*rust_flags,str(ROOT/'src/lib.rs'),'-o',str(native/'libcocycle.rlib')],
        [str(args.rustc),*rust_flags,str(native/'cocycle.rs'),'--extern',f"cocycle={native/'libcocycle.rlib'}",'-L',f'dependency={native}', '-o',str(native/'cocycle')],
        ['g++','-std=c++20','-O3','-pthread','-DOINEUS_DISABLE_ICECREAM','-MMD','-MF',str(native/'oineus.d'),f'-I{args.oineus_source}/include',f'-I{args.oineus_source}/extern/taskflow',f'-I{args.oineus_source}/extern',f'-I{args.boost_include}',str(ROOT/'benches/optimization/oineus.cpp'),'-o',str(native/'oineus')],
    ]
    logs = []
    for command in commands:
        result = subprocess.run(command, cwd=ROOT, text=True, capture_output=True)
        logs.append({'command':command,'returncode':result.returncode,'stdout':result.stdout,'stderr':result.stderr})
        (output/'build.json').write_text(json.dumps(logs,indent=2)+'\n')
        if result.returncode:
            raise RuntimeError(result.stderr)
    return native


def limits(cpu):
    import resource
    os.sched_setaffinity(0,{cpu})
    resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3))


def invoke(native, mode, path, cpu, timeout):
    executable = native/('oineus' if mode.startswith('oineus') else 'cocycle')
    result = subprocess.run([str(executable),mode,str(path)],capture_output=True,text=True,
                            timeout=timeout,preexec_fn=lambda:limits(cpu))
    if result.returncode:
        raise RuntimeError(f'{mode}: exit {result.returncode}: {result.stderr}')
    data=json.loads(result.stdout)
    if data.get('protocol') != PROTOCOL:
        raise RuntimeError('wrong worker protocol')
    data['stderr']=result.stderr
    return data


def validate(output, expected, targets=None, pairs=None, essential=None):
    observed=sorted(output['intervals'],key=lambda x:(x[0],x[1],math_end(x[2])))
    if observed != expected:
        raise ValueError('interval mismatch')
    if targets is not None and output['targets'] != targets:
        raise ValueError('target mismatch')
    if pairs is not None and sorted(map(tuple,output['pairs'])) != pairs:
        raise ValueError('pair mismatch')
    if essential is not None and output['essential'] != essential:
        raise ValueError('essential mismatch')


def summarize(rows):
    groups={}
    for row in rows:
        if row['warmup'] or row['status'] != 'completed':
            continue
        key=(row['fixture'],row['workload'],row['backend'])
        groups.setdefault(key,[]).append(row['output'])
    summary=[]
    for key, samples in groups.items():
        times={name:{'median':statistics.median(s['times_ms'][name] for s in samples),
                     'min':min(s['times_ms'][name] for s in samples),
                     'max':max(s['times_ms'][name] for s in samples)} for name in samples[0]['times_ms']}
        summary.append({'fixture':key[0],'workload':key[1],'backend':key[2],'samples':len(samples),'times_ms':times,
                        'median_peak_rss_kib':statistics.median(s['peak_rss_kib'] for s in samples),
                        'median_input_rss_kib':statistics.median(s['rss_input_kib'] for s in samples),
                        'extra_payload_count':samples[0]['extra_payload_count']})
    return summary


def main(argv=None):
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,required=True)
    parser.add_argument('--rustc',type=Path,required=True)
    parser.add_argument('--oineus-source',type=Path,required=True)
    parser.add_argument('--boost-include',type=Path,required=True)
    parser.add_argument('--dense-source',type=Path,help='unchanged pinned big_steps.rs; defaults to the local research Git object')
    parser.add_argument('--samples',type=int,default=12)
    parser.add_argument('--cpu',type=int,default=0)
    parser.add_argument('--timeout',type=float,default=60.)
    parser.add_argument('--quick',action='store_true')
    parser.add_argument('--build-only',action='store_true')
    args=parser.parse_args(argv)
    if sys.platform != 'linux':
        parser.error('native timing requires Linux/WSL for all workers')
    if args.samples < 1:
        parser.error('--samples must be positive')
    output=args.output.resolve();output.mkdir(parents=True,exist_ok=False)
    native=build(args,output)
    metadata={'protocol':PROTOCOL,'candidate_commit':git('rev-parse','HEAD'),'production_commit':'a380cba76d417c62be46d50495e4d7cce5fef2b3','dense_commit':DENSE_SHA,'oineus_commit':OINEUS_SHA,
              'status':git('status','--porcelain'),'start_utc':datetime.now(timezone.utc).isoformat(),'samples':args.samples,'warmups':1,'cpu':args.cpu,'inherited_affinity':sorted(os.sched_getaffinity(0)),
              'timeout_seconds':args.timeout,'address_space_bytes':2*1024**3,'frequency_control':'none','host_load_control':'none',
              'compiler':subprocess.check_output([str(args.rustc),'--version','--verbose'],text=True),'cpp_compiler':subprocess.check_output(['g++','--version'],text=True),
              'uname':subprocess.check_output(['uname','-a'],text=True),'cpu_info':subprocess.check_output(['lscpu'],text=True),
              'sources':{str(p.relative_to(ROOT)):digest(p) for p in [ROOT/'tools/benchmark_critical_sets.py',ROOT/'benches/optimization/cocycle.rs',ROOT/'benches/optimization/oineus.cpp',ROOT/'src/optimization/mod.rs']},
              'workers':{str(p.relative_to(native)):digest(p) for p in native.rglob('*') if p.is_file()},
              'library_sources':{str(p.relative_to(ROOT)):digest(p) for p in (ROOT/'src').rglob('*.rs')}}
    dependencies=(native/'oineus.d').read_text().replace('\\\n',' ').split(':',1)[1].split()
    metadata['consumed_native_sources']={p:digest(p) for p in dependencies}
    (output/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    if args.build_only:
        return 0
    cases=fixtures(args.quick)
    (output/'fixtures.json').write_text(json.dumps(cases,indent=2)+'\n')
    fixture_dir=output/'fixtures';fixture_dir.mkdir()
    jobs=[]
    for case in cases:
        pairs,essential=pair_oracle(case['cells']);expected=diagram(case['cells'],pairs,essential,case['q'])
        for workload,proposals,repeats in workloads(case):
            path=fixture_dir/f"{case['name']}-{workload}.txt";path.write_text(fixture_text(case,proposals,repeats))
            jobs.append((case,workload,path,['dense','critical','oineus_partial','oineus_full'],expected,pairs,essential))
        path=fixture_dir/f"{case['name']}-diagram.txt";path.write_text(fixture_text(case,[],1))
        jobs.append((case,'persistence-context',path,['reference','diagram','filtered','representatives']+(['flag'] if case['flag'] else []),expected,None,None))
    plan={'jobs':[{'fixture':c['name'],'workload':w,'fixture_sha256':digest(p),'backends':b,'cells':len(c['cells']),'q':c['q']} for c,w,p,b,*_ in jobs],
          'sampling':f'one discarded warmup; {args.samples} seeded position-balanced fresh-process rounds; serial execution',
          'critical_time':'frozen input; D preparation + primal reduction, cold targets, workspace cleanup; normalized pairing transport excluded; owned targets retained',
          'warm_time':'separate fresh workspace prewarmed by same batch; repeated targets including result destruction; average batch ms',
          'context_time':'existing diagram paths including normalized interval export and cleanup; representative validation excluded; different requested outputs so no critical-set speedup ratio',
          'exclusions':'flag path excludes delayed fans and supplied 3-skeletons; finite endpoints only; no optimizer gradients or convergence ranking'}
    (output/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    rows=[]
    rng=random.Random(220316749)
    failed=False
    with (output/'samples.jsonl').open('w') as log:
        for case,workload,path,backends,expected,pairs,essential in jobs:
            baseline_targets=None
            if workload!='persistence-context':
                baseline=invoke(native,'dense',path,args.cpu,args.timeout)
                validate(baseline,expected,pairs=pairs,essential=essential)
                baseline_targets=baseline['targets']
            order=backends[:];rng.shuffle(order)
            for round_index in range(-1,args.samples):
                rotated=order[round_index%len(order):]+order[:round_index%len(order)]
                for position,backend in enumerate(rotated):
                    row={'fixture':case['name'],'workload':workload,'backend':backend,'round':round_index,'position':position,'warmup':round_index<0}
                    try:
                        observed=invoke(native,backend,path,args.cpu,args.timeout)
                        row['output']=observed
                        validate(observed,expected,baseline_targets,pairs,essential)
                        row.update(status='completed',output=observed)
                    except (RuntimeError,ValueError,subprocess.TimeoutExpired) as error:
                        row.update(status='failed',error=str(error));failed=True
                    rows.append(row);log.write(json.dumps(row)+'\n');log.flush()
            (output/'summary.json').write_text(json.dumps(summarize(rows),indent=2)+'\n')
            print(f"{case['name']} {workload}: {len(backends)*args.samples} measured calls completed",flush=True)
    metadata['end_utc']=datetime.now(timezone.utc).isoformat()
    metadata['planned_measured']=sum(len(b) for _,_,_,b,*_ in jobs)*args.samples
    metadata['completed_measured']=sum(r['status']=='completed' and not r['warmup'] for r in rows)
    metadata['failed_measured']=sum(r['status']!='completed' and not r['warmup'] for r in rows)
    (output/'metadata.json').write_text(json.dumps(metadata,indent=2)+'\n')
    return int(failed)


if __name__=='__main__':
    raise SystemExit(main())
