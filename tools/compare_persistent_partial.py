"""Check composed matrix-free spectral routes against exact-rank/NumPy oracles.

Uses the existing explicit signed-complex controls, Fraction kernel oracle and
independent NumPy/SVD spectrum. GUDHI corroborates topology over F11, including
the declared RP2 field distinction; it is not a real eigensolver. No benchmark
claim is made by this correctness command.
"""
import argparse
import hashlib
import json
from pathlib import Path
import platform
import shutil
import subprocess

import compare_persistent_laplacians as reference


def main():
    import numpy as np
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--gudhi", action="store_true")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    worker = (args.output / args.worker.name).resolve()
    shutil.copy2(args.worker, worker)
    cases = reference.fixtures()
    groups, expected = [], {}
    for label, cells in cases:
        scales = sorted(set([-1, *cells.values()]))
        top = max((len(s)-1 for s in cells), default=0)
        for q in range(top+2):
            for t in scales:
                starts = [s for s in scales if s <= t]
                groups.append((label, cells, q, starts, t))
                for s in starts:
                    k, _, h, beta, exact, up, retained, _, _ = reference.oracle(cells, s, t, q)
                    n = len(k)
                    assert n <= 64, "these complete-block controls do not certify general iterative indexes"
                    if exact is not None:
                        h = np.array(exact, dtype=float).reshape(n, n)
                    values = np.linalg.eigvalsh(h)
                    outside = [i for i in range(len(up)) if i not in retained]
                    eligible = not np.any(np.count_nonzero(up[outside], axis=1) > 2)
                    expected[label, q, s, t] = dict(n=n, beta=beta, rank=n-beta,
                        matrix=h.tolist(), spectrum=values.tolist(), compact_eligible=bool(eligible))
    text = "".join(reference.protocol(cells, q, starts, t)
                   for _, cells, q, starts, t in groups)
    (args.output / "inputs.txt").write_text(text, encoding="utf-8")
    checks = []
    for mode in ["smallest", "largest", "positive", "compact", "fixed-cold", "fixed", "relevant"]:
        run = subprocess.run([str(worker), "--partial", mode], input=text,
                             text=True, capture_output=True, timeout=180)
        (args.output / f"{mode}.stdout").write_text(run.stdout, encoding="utf-8")
        (args.output / f"{mode}.stderr").write_text(run.stderr, encoding="utf-8")
        assert run.returncode == 0, (mode, run.stderr[-2000:])
        records = iter(json.loads(line) for line in run.stdout.splitlines())
        for label, _, q, starts, end in groups:
            pairs = ([(s,t) for i,t in enumerate(starts) for s in starts[:i+1]
                      if expected[label,q,s,t]["n"]] if mode == "relevant"
                     else [(s,end) for s in starts if mode not in ["fixed", "fixed-cold"]
                           or expected[label,q,s,end]["n"]])
            for s,t in pairs:
                record = next(records)
                oracle = expected[label,q,s,t]
                assert (record["degree"], record["s"], record["t"]) == (q,s,t)
                failure = ((mode == "compact" and not oracle["compact_eligible"])
                           or (mode == "positive" and oracle["n"] > 0 and oracle["rank"] == 0))
                assert bool(record.get("error")) == failure, (mode,label,q,s,t,record,oracle)
                if failure:
                    checks.append(dict(mode=mode, fixture=label, degree=q, s=s, t=t,
                                       expected_rejection=True))
                    continue
                n, rank = oracle["n"], oracle["rank"]
                values = record["values"]
                wanted = ([] if n == 0 else [oracle["spectrum"][0 if mode == "smallest"
                                                   else n-rank if mode == "positive" else -1]])
                assert len(values) == len(wanted)
                error = float(np.linalg.norm(np.array(values)-wanted) / max(1.,np.linalg.norm(wanted)))
                v = np.array(record["vectors"]).reshape(n, len(values))
                h = np.array(oracle["matrix"]).reshape(n,n)
                residual = float(np.linalg.norm(h @ v-v*np.array(values)) / max(1.,np.linalg.norm(h)))
                orthogonality = float(np.linalg.norm(v.T@v-np.eye(len(values))))
                assert error <= 1e-8 and residual <= 1e-8 and orthogonality <= 1e-8, (mode,label,q,s,t,error,residual,orthogonality)
                checks.append(dict(mode=mode, fixture=label, degree=q, s=s, t=t,
                                   eigenvalue_error=error, residual=residual, orthogonality=orthogonality))
        assert next(records, None) is None, mode
    topology = []
    if args.gudhi:
        import gudhi
        for label,cells in cases:
            st = gudhi.SimplexTree()
            for simplex,value in cells.items():
                st.insert(simplex,filtration=value)
            st.compute_persistence(homology_coeff_field=11, min_persistence=-1,
                                   persistence_dim_max=True)
            for (name,q,s,t),oracle in expected.items():
                if name != label:
                    continue
                beta = sum(b <= s and d > t for b,d in st.persistence_intervals_in_dimension(q))
                assert beta == oracle["beta"], (label,q,s,t,beta,oracle["beta"])
            topology.append(dict(fixture=label, version=gudhi.__version__, field=11,
                                 scope="topology only; real ranks use Fraction elimination"))
    source = Path(__file__).resolve().parents[1]
    sources = [Path(__file__), source/"examples/persistent_laplacian.rs",
               *sorted((source/"src/persistent_laplacian").glob("*.rs")), source/"src/eigensolver.rs"]
    summary = dict(checks=len(checks), expected_rejections=sum(c.get("expected_rejection",False) for c in checks),
        modes=7, tolerance=1e-8, max_eigenvalue_error=max((c.get("eigenvalue_error",0.) for c in checks),default=0.),
        max_residual=max((c.get("residual",0.) for c in checks),default=0.),
        topology=topology, numpy=np.__version__, platform=platform.platform(),
        source_sha=subprocess.check_output(["git","rev-parse","HEAD"],cwd=source,text=True).strip(),
        source_hashes={str(p.relative_to(source)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},
        worker_sha256=hashlib.sha256(worker.read_bytes()).hexdigest())
    (args.output/"oracles.json").write_text(json.dumps({repr(k):v for k,v in expected.items()},indent=2)+"\n",encoding="utf-8")
    (args.output/"checks.json").write_text(json.dumps(checks,indent=2)+"\n",encoding="utf-8")
    (args.output/"summary.json").write_text(json.dumps(summary,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    main()
