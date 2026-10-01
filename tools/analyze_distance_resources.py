"""Frozen CPU/PSS model contrasts and phase replay for the resource protocol."""
import argparse
import json
import math
from pathlib import Path
import statistics

from distance_common import save


def interpolate(points,t):
    """Linear interpolation, retaining endpoint residence outside the window."""
    if t<=points[0][0]:return points[0][1]
    for a,b in zip(points,points[1:]):
        if t<=b[0]:
            return b[1] if a[0]==b[0] else a[1]+(b[1]-a[1])*(t-a[0])/(b[0]-a[0])
    return points[-1][1]


def replay(traces,starts,end,baseline):
    """Replay incremental PSS above measured group-ready residence."""
    if len(traces)!=len(starts) or not traces or end<=0:raise ValueError('invalid replay')
    times={0,end}
    for points,start in zip(traces,starts):
        if not points or any(b[0]<a[0] for a,b in zip(points,points[1:])):raise ValueError('invalid trace')
        times.update(t+start for t,_ in points if 0<=t+start<=end)
    values=[]
    for t in sorted(times):
        values.append((t,baseline+sum(interpolate(p,t-s)-p[0][1] for p,s in zip(traces,starts))))
    integral=sum((b[0]-a[0])*(a[1]+b[1])/2 for a,b in zip(values,values[1:]))
    return dict(peak=max(v for _,v in values),integral=integral,points=values)


def linear_inflation(rows,key,base):
    points=[(r['workers']-1,r[key]/base-1) for r in rows if r.get(key) is not None and r['workers']>1]
    denominator=sum(x*x for x,_ in points)
    return sum(x*y for x,y in points)/denominator if denominator else 0.


def fit(rows,trace_rows,capacity=22,memory_mib=128):
    """Fit tuning only; held-out samples are rejected at this boundary."""
    if any(r['split']!='tuning' for r in rows+trace_rows):raise ValueError('holdout leakage')
    types=sorted({r['type'] for r in rows})
    cpu_capacity=min(capacity,max(r['throughput']*r['cpu'] for r in rows))
    model=dict(cpu_capacity=cpu_capacity,memory_mib=memory_mib,types={})
    for kind in types:
        cases=[r for r in rows if r['type']==kind]
        traces=[r for r in trace_rows if r['type']==kind]
        first=next(r for r in cases if r['workers']==1)
        trace=next(r for r in traces if r['workers']==1)
        baseline=dict(time=1/first['throughput'],cpu=first['cpu'],memory=trace['memory'],peak=trace['peak'])
        if baseline['cpu']<=0 or baseline['memory']<=0:raise ValueError('unresolved resource demand')
        baseline.update(time_slope=linear_inflation(cases,'service',baseline['time']),
                        cpu_slope=linear_inflation(cases,'cpu',baseline['cpu']),
                        memory_slope=linear_inflation(traces,'memory',baseline['memory']))
        model['types'][kind]=baseline
    return model


def demand(model,kind,n,level):
    b=model['types'][kind]
    if level==0:return dict(time=b['time'],cpu=b['time'],memory=b['time']*b['peak'])
    if level==1:return {k:b[k] for k in ('time','cpu','memory')}
    return {k:b[k]*max(.05,1+b[k+'_slope']*(n-1)) for k in ('time','cpu','memory')}


def predict(model,kind,n,level):
    d=demand(model,kind,n,level)
    return dict(throughput=min(n/d['time'],model['cpu_capacity']/d['cpu'],model['memory_mib']/d['memory']),**d)


def fit_interference(model,rows):
    if any(r['split']!='tuning' for r in rows):raise ValueError('holdout leakage')
    values={}
    for row in rows:
        kinds=row['types'];n=len(kinds)
        for index,kind in enumerate(kinds):
            others=set(kinds)-{kind}
            if len(others)!=1:continue
            other=others.pop();fraction=kinds.count(other)/n
            observed=row['child_service'][index]
            coefficient=(observed/demand(model,kind,n,2)['time']-1)/fraction
            values.setdefault(kind+'<-'+other,[]).append(coefficient)
    model['interference']={k:statistics.median(v) for k,v in values.items()}
    return model


def mixed_predict(model,row,level):
    kinds=row['types'];n=len(kinds);demands=[]
    for kind in kinds:
        d=demand(model,kind,n,min(level,2))
        factor=1.
        if level==3:
            factor=max(.05,1+sum(kinds.count(other)/n*model['interference'].get(kind+'<-'+other,0) for other in set(kinds)-{kind}))
        demands.append({k:v*factor for k,v in d.items()})
    jobs=row['jobs'];total=sum(jobs)
    wall=max(j*d['time'] for j,d in zip(jobs,demands))
    cpu=sum(j*d['cpu'] for j,d in zip(jobs,demands))/total
    memory=sum(j*d['memory'] for j,d in zip(jobs,demands))/total
    return min(total/wall,model['cpu_capacity']/cpu,model['memory_mib']/memory)


def load_run(directory):
    manifest=json.loads((directory/'manifest.json').read_text())
    status=json.loads((directory/'summary.json').read_text())
    if status['status']!='passed' or manifest['dirty'] or manifest['exploratory']:raise ValueError('incomplete/nonformal run')
    records=json.loads((directory/'results.json').read_text())
    rows=[]
    for cell in records:
        samples=[json.loads(p.read_text()) for p in (directory/cell['name']).glob('r*/result.json')]
        for variant in ([cell['name']] if 'nodes' in cell else manifest['variants']):
            groups=[s for s in samples if not s['warmup'] and s['worker']==variant]
            if len(groups)!=manifest['samples']:raise ValueError('incomplete groups')
            counts=lambda s:sum(len(d) for c in s['children'] for d in c['durations_ns'])
            med=lambda key:statistics.median(s[key] for s in groups)
            r=dict(run=directory.name,cell=cell['name'],workers=cell['workers'],mode=cell['mode'],trace=manifest['trace'],
                   throughput=med('throughput_per_s'),minimum=min(s['throughput_per_s'] for s in groups),maximum=max(s['throughput_per_s'] for s in groups),
                   cpu=statistics.median(s['cpu_core_seconds']/counts(s) for s in groups),
                   memory=statistics.median(s['pss_integral_byte_seconds']/1048576/counts(s) for s in groups) if manifest['trace'] else None,
                   peak=statistics.median(s['sampled_peak_pss_bytes']/1048576 for s in groups) if manifest['trace'] else med('endpoint_pss_bytes')/1048576,
                   median_us=med('median_job_ns')/1000,p95_us=med('p95_job_ns')/1000,groups=len(groups))
            if 'nodes' in cell:
                r.update(types=[n['family']+'-'+n['metric']+'-'+n['variant'] for n in cell['nodes']],split=cell['nodes'][0]['split'],
                         jobs=[c['jobs_per_thread'] for c in groups[0]['children']],
                         child_service=[statistics.median((s['completion_ns'][i]-s['dispatch_ns'][i])/1e9/s['children'][i]['jobs_per_thread'] for s in groups) for i in range(cell['workers'])])
            else:
                r.update(type=cell['family']+'-'+cell['metric']+'-'+variant,split=cell['split'],family=cell['family'],metric=cell['metric'],variant=variant,size=cell['size'],
                         service=cell['workers']/r['throughput'])
            rows.append(r)
    return rows


def evaluate(train,trace,mixed,holdout,held_mixed,output):
    model=fit(train,trace)
    fit_interference(model,mixed)
    errors=[]
    for row in holdout:
        if row['trace']:continue
        for level in range(4):
            prediction=predict(model,row['type'],row['workers'],min(level,2))
            errors.append(dict(type=row['type'],workers=row['workers'],level=level,observed=row['throughput'],predicted=prediction['throughput'],
                               throughput_error=abs(prediction['throughput']/row['throughput']-1),cpu_error=abs(prediction['cpu']/row['cpu']-1)))
    for row in held_mixed:
        if row['trace']:continue
        for level in range(4):
            predicted=mixed_predict(model,row,level)
            errors.append(dict(cell=row['cell'],types=row['types'],workers=row['workers'],level=level,observed=row['throughput'],predicted=predicted,
                               throughput_error=abs(predicted/row['throughput']-1)))
    p95=lambda values:sorted(values)[math.ceil(.95*len(values))-1]
    summary={str(level):dict(median_error=statistics.median(r['throughput_error'] for r in errors if r['level']==level),
                             p95_error=p95([r['throughput_error'] for r in errors if r['level']==level]),
                             maximum_error=max(r['throughput_error'] for r in errors if r['level']==level),
                             homogeneous_median=statistics.median(r['throughput_error'] for r in errors if r['level']==level and 'type' in r),
                             mixed_median=statistics.median(r['throughput_error'] for r in errors if r['level']==level and 'types' in r)) for level in range(4)}
    save(output/'model.json',model);save(output/'prediction-errors.json',errors);save(output/'model-summary.json',summary)
    knees=[];crossovers=[]
    performance=[r for r in holdout if not r['trace']]
    for kind in sorted({r['type'] for r in performance}):
        curve=[r for r in performance if r['type']==kind]
        observed=observed_knee(curve)
        for level in range(4):
            predicted=observed_knee([dict(workers=r['workers'],throughput=predict(model,kind,r['workers'],min(level,2))['throughput']) for r in curve])
            knees.append(dict(type=kind,level=level,observed=observed,predicted=predicted,error=abs(predicted-observed)))
    cells={}
    for r in performance:cells.setdefault((r['family'],r['metric'],r['workers']),[]).append(r)
    for key,rows in cells.items():
        winner=max(rows,key=lambda r:r['throughput'])['variant']
        for level in range(4):
            selected=max(rows,key=lambda r:predict(model,r['type'],r['workers'],min(level,2))['throughput'])['variant']
            crossovers.append(dict(cell=key,level=level,observed=winner,predicted=selected,mismatch=winner!=selected))
    save(output/'knee-errors.json',knees);save(output/'crossover-errors.json',crossovers)
    return summary


def resource_errors(model,rows):
    errors=[]
    for row in rows:
        if row['split']!='holdout' or not row['trace']:continue
        for level in range(4):
            d=demand(model,row['type'],row['workers'],min(level,2))
            peak=row['workers']*(model['types'][row['type']]['peak'] if level==0 else d['memory']/d['time'])
            errors.append(dict(type=row['type'],workers=row['workers'],level=level,
                               cpu_error=abs(d['cpu']/row['cpu']-1),memory_error=abs(d['memory']/row['memory']-1),
                               peak_error=abs(peak/row['peak']-1),predicted_peak=peak,observed_peak=row['peak'],
                               feasibility=[dict(budget_mib=b,predicted=peak<=b,observed=row['peak']<=b) for b in (2,4,8,16,32,128)],
                               predicted_cpu=d['cpu'],observed_cpu=row['cpu'],predicted_memory=d['memory'],observed_memory=row['memory']))
    return errors


def observed_knee(rows,threshold=.95):
    """Earliest measured count attaining 95% of the maximum on this finite grid."""
    peak=max(r['throughput'] for r in rows)
    return min(r['workers'] for r in rows if r['throughput']>=threshold*peak)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--train',type=Path,nargs='+',required=True)
    p.add_argument('--resource',type=Path,nargs='+',required=True)
    p.add_argument('--mixed',type=Path,nargs='+',required=True)
    p.add_argument('--holdout',type=Path,nargs='+',required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();args.output.mkdir(exist_ok=False,parents=True)
    load=lambda paths:[r for d in paths for r in load_run(d)]
    train=[r for r in load(args.train) if r['split']=='tuning' and r['mode']=='process' and not r['trace']]
    resources=load(args.resource)
    mixed=load(args.mixed)
    heldout=[r for r in load(args.holdout) if r['split']=='holdout' and r['mode']=='process']
    print(json.dumps(evaluate(train,[r for r in resources if r['split']=='tuning'],[r for r in mixed if r['split']=='tuning' and not r['trace']],heldout,[r for r in mixed if r['split']=='holdout'],args.output)))
    model=json.loads((args.output/'model.json').read_text())
    save(args.output/'resource-errors.json',resource_errors(model,resources))


if __name__=='__main__':main()
