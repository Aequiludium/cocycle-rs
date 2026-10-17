"""Pinned F2 operator correctness and isolated public-workflow measurements.

The reference is external verification infrastructure, never a Rust dependency.
Raw fixtures, output, errors and timing samples go to a fresh ignored directory.
"""

import argparse
import hashlib
import json
from pathlib import Path
import platform
import os
import random
import statistics
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
PIN = "a11bf50b8c8320651e9ac015566778d957af9c17"


def digest(paths):
    value = hashlib.sha256()
    for path in sorted(paths):
        value.update(path.name.encode())
        value.update(path.read_bytes())
    return value.hexdigest()


def git(root, *arguments):
    """Read Windows-managed worktree metadata from WSL without rewriting it."""
    root = Path(root).resolve()
    command = ["git", "-C", str(root)]
    pointer = root / ".git"
    if sys.platform != "win32" and pointer.is_file():
        target = pointer.read_text().strip().removeprefix("gitdir: ").replace("\\", "/")
        if len(target) > 2 and target[1] == ":":
            target = "/mnt/" + target[0].lower() + target[2:]
            command = ["git", "--git-dir", target, "--work-tree", str(root)]
    environment = {**os.environ, "GIT_OPTIONAL_LOCKS": "0"}
    return subprocess.check_output(command + list(arguments), text=True, env=environment, stderr=subprocess.PIPE)


def kernel(rows, n):
    """Fixture-only byte elimination, independent of packed Rust operations."""
    rows = [list(row) for row in rows]
    pivots = []
    for col in range(n):
        pivot = next((i for i in range(len(pivots), len(rows)) if rows[i][col]), None)
        if pivot is None:
            continue
        i = len(pivots)
        rows[i], rows[pivot] = rows[pivot], rows[i]
        for j in range(len(rows)):
            if j != i and rows[j][col]:
                rows[j] = [a ^ b for a, b in zip(rows[j], rows[i])]
        pivots.append(col)
    basis = []
    for free in range(n):
        if free not in pivots:
            x = [0] * n
            x[free] = 1
            for i, pivot in enumerate(pivots):
                x[pivot] = rows[i][free]
            basis.append(x)
    return basis


def fixture(n, seed, family="general", exhaustive=False):
    rng = random.Random(seed)
    ra = n // 3 if family != "homology" else 0
    rd = n // 3 if family != "homology" else 0
    m, p = ra + (n > 0), rd + (n > 0)
    a = [[int(i == j) for j in range(n)] for i in range(ra)]
    if n:
        a.append([0] * n)
    d = [[int(j < rd and i == ra + j) for j in range(p)] for i in range(n)]
    for _ in range(0 if family == "sparse" else 5 * n):
        if n < 2:
            break
        i, j = rng.sample(range(n), 2)
        if rng.randrange(2):
            for row in a:
                row[i], row[j] = row[j], row[i]
            d[i], d[j] = d[j], d[i]
        else:
            for row in a:
                row[j] ^= row[i]
            d[i] = [x ^ y for x, y in zip(d[i], d[j])]
    cycles = kernel(a, n)
    vectors = [[0] * n]
    if exhaustive:
        vectors = [[(x >> j) & 1 for j in range(n)] for x in range(1 << n)]
    else:
        for _ in range(8):
            z = [0] * n
            for v in cycles:
                if rng.randrange(2):
                    z = [x ^ y for x, y in zip(z, v)]
            vectors.append(z)
        vectors += [[rng.randrange(2) for _ in range(n)] for _ in range(4)]
    return {"m": m, "n": n, "p": p, "a": a, "d": d,
            "weights": [rng.randrange(1, 20) for _ in range(n)],
            "vectors": vectors, "family": family, "seed": seed}


def wire(data):
    values = [data["m"], data["n"], data["p"], len(data["vectors"])]
    for key in ("a", "d"):
        values += [bit for row in data[key] for bit in row]
    values += data["weights"]
    values += [bit for row in data["vectors"] for bit in row]
    return " ".join(map(str, values)) + "\n"


def rss(field):
    try:
        for line in Path("/proc/self/status").read_text().splitlines():
            if line.startswith(field):
                return int(line.split()[1])
    except OSError:
        pass
    return None


def reference(args):
    sys.path.insert(0, str(args.reference / "src"))
    if args.native_dir:
        sys.path.insert(0, str(args.native_dir))
    from homology_operator import (ChainWindow, Matrix, HomologyOperator,
                                  ProjectionProblem, ResourceLimits, solve_projection)
    data = json.loads(args.fixture.read_text())
    a = Matrix.from_rows(data["a"], ncols=data["n"])
    d = Matrix.from_rows(data["d"], ncols=data["p"])
    window = ChainWindow(1, a, d, tuple(f"previous:{i}" for i in range(data["m"])),
                         tuple(f"current:{i}" for i in range(data["n"])),
                         tuple(f"next:{i}" for i in range(data["p"])),
                         tuple(data["weights"]), arithmetic="ExactInteger", unit="cost")
    limits = ResourceLimits(state_limit=100000000, wall_time_limit=600,
                            matrix_entry_limit=100000000)
    native = args.engine == "native"
    problem = ProjectionProblem(window, resource_limits=limits, matrix_free_output=native)

    def solve():
        solution = solve_projection(problem, "NativeFactorizedSolver" if native else "FeasibleSolver")
        if solution.status != "FeasibleOnly" or solution.certificate_level != "Feasible":
            raise RuntimeError((solution.status, solution.diagnostics))
        return HomologyOperator(window, solution, repository_revision=PIN)

    baseline_rss, baseline_hwm = rss("VmRSS:"), rss("VmHWM:")
    op = solve()
    vectors = [tuple(z) for z in data["vectors"]]
    cycles = [z for z in vectors if op.is_cycle(z)]
    started = time.perf_counter_ns()
    for _ in range(args.iterations):
        if args.mode == "solve":
            result = solve()
        else:
            result = [(op.class_representative(z), op.selected_mass(z)) for z in cycles]
    elapsed = time.perf_counter_ns() - started
    peak = rss("VmHWM:")
    del result
    columns = [op.project(tuple(int(i == j) for i in range(data["n"])))
               for j in range(data["n"])]
    projection = [[columns[j][i] for j in range(data["n"])] for i in range(data["n"])]
    print(json.dumps({"state": "Computed", "betti": op.betti(), "projection": projection,
                      "actions": [op.project(z) for z in vectors],
                      "masses": [op.selected_mass(z) if op.is_cycle(z) else None for z in vectors],
                      "elapsed_ns": elapsed, "iterations": args.iterations,
                      "baseline_rss_kib": baseline_rss, "baseline_hwm_kib": baseline_hwm,
                      "peak_rss_kib": peak}))


def gudhi_probe():
    import gudhi
    simplices = [([0], 0.), ([1], 0.), ([2], 0.),
                 ([0, 1], 1.), ([0, 2], 1.), ([1, 2], 1.), ([0, 1, 2], 2.)]
    outputs = []
    for through in (0., 1., 2.):
        tree = gudhi.SimplexTree()
        supplied = [(s, value) for s, value in simplices if value <= through]
        for simplex, value in supplied:
            tree.insert(simplex, filtration=value)
        actual = sorted((tuple(s), v) for s, v in tree.get_filtration())
        if actual != sorted((tuple(s), v) for s, v in supplied):
            raise RuntimeError("GUDHI silently changed the supplied stage")
        tree.compute_persistence(homology_coeff_field=2, min_persistence=0,
                                 persistence_dim_max=True)
        outputs.append(list(tree.betti_numbers()) + [0] * (2 - len(tree.betti_numbers())))
    tree.compute_persistence(homology_coeff_field=2, min_persistence=0, persistence_dim_max=True)
    h1 = tree.persistence_intervals_in_dimension(1).tolist()
    print(json.dumps({"version": gudhi.__version__, "betti_by_stage": outputs, "h1": h1,
                      "coefficient_field": 2, "oracle_used_for_result": False}))


def invoke(command, raw_path, input_text=None):
    try:
        proc = subprocess.run(command, input=input_text, text=True, capture_output=True,
                              timeout=600, check=False)
        record = {"command": list(map(str, command)), "exit_code": proc.returncode,
                  "stdout": proc.stdout, "stderr": proc.stderr}
        raw_path.write_text(json.dumps(record, indent=2) + "\n")
        if proc.returncode:
            raise RuntimeError(f"worker failed; retained {raw_path}")
        return json.loads(proc.stdout)
    except subprocess.TimeoutExpired as error:
        raw_path.write_text(json.dumps({"state": "Timeout", "command": list(map(str, command)),
                                       "timeout": error.timeout}) + "\n")
        raise


def family_probe(args):
    sys.path.insert(0, str(args.reference / "src"))
    from homology_operator import (ChainWindow, Matrix, HomologyOperator, OperatorFamily,
                                  ProjectionProblem, solve_projection)
    a = Matrix.from_rows(((1, 1, 0), (1, 0, 1), (0, 1, 1)))
    bases = (("v0", "v1", "v2"), ("e0", "e1", "e2"))
    ring = ChainWindow(1, a, Matrix.zero(3, 0), *bases, (), (2, 3, 5), arithmetic="ExactInteger")
    fill = ChainWindow(1, a, Matrix.from_rows(((1,), (1,), (1,))), *bases, ("f",), (2, 3, 5), arithmetic="ExactInteger")
    windows = (ring, ring, fill)
    operators = tuple(HomologyOperator(w, solve_projection(ProjectionProblem(w))) for w in windows)
    family = OperatorFamily((0., 1., 2.), windows, operators)
    ranks = [[family.transport_rank(i, j).value for j in range(i, 3)] for i in range(3)]
    bars = [[b["birth_stage"], b["death_stage"], b["multiplicity"]] for b in family.barcode().value]
    print(json.dumps({"ranks": ranks, "bars": bars}))


def expected_span(data, result):
    n = data["n"]
    if n > 5:
        return "NotApplicable: exhaustive tiny-span limit is 5"
    columns = [sum(row[j] << i for i, row in enumerate(data["d"])) for j in range(data["p"])]
    boundaries = {0}
    for v in columns:
        boundaries |= {x ^ v for x in boundaries.copy()}
    cycles = {x for x in range(1 << n)
              if all(sum(bit * ((x >> j) & 1) for j, bit in enumerate(row)) % 2 == 0
                     for row in data["a"])}
    if len(cycles) // len(boundaries) != 1 << result["betti"]:
        raise RuntimeError("independent quotient cardinality differs")
    for z, p in zip(data["vectors"], result["actions"]):
        x = sum(bit << i for i, bit in enumerate(z))
        image = sum(bit << i for i, bit in enumerate(p))
        if image not in cycles or (x in cycles and x ^ image not in boundaries):
            raise RuntimeError("independent span oracle rejected the action")
    return "Passed"


def compare(first, second):
    if first.get("state") != "Computed" or second.get("state") != "Computed":
        raise RuntimeError("process completion is not semantic completion")
    for key in ("betti", "projection", "actions", "masses"):
        if first[key] != second[key]:
            raise RuntimeError(f"exact mismatch in {key}")


def geometry_checks(args):
    """Replay pinned T1--T5 with independent raw-simplex feasible-set enumeration."""
    from itertools import combinations, product
    sys.path[:0] = [str(args.reference / "src"), str(args.reference / "tests")]
    from homology_operator import minimum_filling
    from test_stable_volume import build_case, upper_costs
    corpus_path = args.reference / "tests/fixtures/geometry_queries.json"
    corpus = json.loads(corpus_path.read_text(encoding="utf-8"))
    canonical = json.dumps(corpus["cases"], sort_keys=True, separators=(",", ":")).encode()
    if hashlib.sha256(canonical).hexdigest() != corpus["corpus_sha256"]:
        raise ValueError("geometry corpus digest differs")
    results = []
    for case in corpus["cases"]:
        # Cocycle owns the frozen filtration order, including decreasing colex
        # within tied dimensions. Give every implementation that same order.
        events = sorted(zip(case["levels"], case["simplices"]),
                        key=lambda x: (x[0], len(x[1]), tuple(-v for v in reversed(x[1]))))
        case = {**case, "levels": [x[0] for x in events], "simplices": [x[1] for x in events]}
        cells = list(map(tuple, case["simplices"]))
        wire_values = [str(case["degree"]), str(len(cells))]
        for level, cell in events:
            wire_values.extend(map(str, [level, len(cell), *cell]))
        output = invoke([str(args.worker), "geometry"], args.output / f"geometry-{case['id']}.json", " ".join(wire_values) + "\n")
        if output.get("state") != "Computed":
            raise RuntimeError("geometry worker failed")
        family = build_case(case)
        references = family.features()
        if references.state != "Computed" or len(references.value) != len(output["features"]):
            raise RuntimeError(f"feature count differs for {case['id']}")
        for record, ref in zip(output["features"], references.value):
            birth, death = record["birth"], record["death"]
            if (birth, death, record["occurrence"]) != (ref["birth_stage"], ref["death_stage"], ref["occurrence"]):
                raise RuntimeError(f"feature endpoints differ for {case['id']}")
            histories = [list(family.resolve_feature(ref, i).value["representative"])
                         for i in range(birth, len(cells) if death is None else death)]
            if record["histories"] != histories:
                raise RuntimeError(f"historical generators differ for {case['id']} birth={birth}")
            minimum = family.stage(birth).minimum_class_mass(tuple(histories[0]))
            if minimum.state != "Computed" or (record["minimum_mass"], record["minimum"]) != (minimum.value, list(minimum.details["representative"])):
                raise RuntimeError("minimum class mass differs")
            basis = [s for s in cells[:birth+1] if len(s) == case["degree"]+1]
            if record["locations"] != [list(s) for s, bit in zip(basis, histories[0]) if bit]:
                raise RuntimeError("source locations differ")
            if case["degree"] == 1:
                cycle = family.optimal_1_cycle(ref)
                if cycle.state != "Computed" or record["cycle"] != list(cycle.value["chain"]):
                    raise RuntimeError("birth-pair minimum cycle differs")
            if death is None:
                if record["volumes"]:
                    raise RuntimeError("infinite feature has a finite volume")
                continue
            upper_cells = [s for s in cells[:death+1] if len(s) == case["degree"]+2]
            lower_cells = [s for s in cells[:death+1] if len(s) == case["degree"]+1]
            for weighted, actual in enumerate(record["volumes"]):
                weights = tuple(cells.index(s) % 5 + 1 if weighted else 1 for s in upper_cells)
                upper = upper_costs(family, ref, weights if weighted else None)
                volume = family.optimal_volume(ref, upper)
                stable = family.stable_volume(ref, 1, upper, epsilon_unit="filtration time")
                z = family.inclusion(birth, death).apply(tuple(histories[0]))
                filling = minimum_filling(family.windows[death], z, upper)
                for query, key, mass_key in ((volume, "volume", "mass"), (stable, "stable", "stable_mass"), (filling, "filling", "filling_mass")):
                    if query.state != "Computed" or actual[key] != list(query.value["chain"]) or actual[mass_key] != query.value["mass"]:
                        raise RuntimeError(f"{key} differs in {case['id']} birth={birth}")
                # Raw set boundaries enumerate *every* upper chain; no reference
                # matrix, rank, projection or constraint solver is used here.
                feasible = {"volume": [], "stable": [], "filling": []}
                for bits in product((0, 1), repeat=len(upper_cells)):
                    selected = {s for s, b in zip(upper_cells, bits) if b}
                    boundary = set()
                    for s in selected:
                        boundary.symmetric_difference_update(combinations(s, len(s)-1))
                    mass = sum(w*b for w,b in zip(weights,bits))
                    key = (mass, sum(b << i for i,b in enumerate(bits)))
                    if boundary == {s for s,b in zip(lower_cells,z) if b}:
                        feasible["filling"].append(key)
                    for name, epsilon in (("volume", None), ("stable", 1)):
                        after = lambda s: (cells.index(s) > birth if epsilon is None else case["levels"][cells.index(s)] >= case["levels"][birth]+epsilon)
                        allowed = {s for s in upper_cells if cells.index(s) < death and after(s)} | {cells[death]}
                        forbidden = {s for s in lower_cells if cells.index(s) < death and after(s)}
                        if cells[death] in selected and selected <= allowed and not boundary & forbidden and (epsilon is not None or cells[birth] in boundary):
                            feasible[name].append(key)
                for name, mass_key in (("volume", "mass"), ("stable", "stable_mass"), ("filling", "filling_mass")):
                    actual_key = (actual[mass_key], sum(b << i for i,b in enumerate(actual[name])))
                    if not feasible[name] or min(feasible[name]) != actual_key:
                        raise RuntimeError(f"independent raw-simplex enumeration rejected {name}")
        results.append({"case": case["id"], "features": len(output["features"]), "state": "Passed",
                        "reference_and_raw_set_oracle": True})
    return {"corpus_sha256": corpus["corpus_sha256"], "file_sha256": hashlib.sha256(corpus_path.read_bytes()).hexdigest(), "cases": results,
            "scope": "Cocycle source order; integer costs and epsilon=1; unit and weighted upper chains; no perturbation guarantee"}


def main(args):
    if args.output.exists():
        raise ValueError("output directory must be fresh")
    args.output.mkdir(parents=True)
    engines = ["rust", "python"] + (["native"] if args.native_dir else []) + (["baseline"] if args.baseline_worker else [])
    reference_head = git(args.reference, "rev-parse", "HEAD").strip()
    if reference_head != PIN:
        raise ValueError(f"reference must be pinned to {PIN}")
    git(args.reference, "diff", "--exit-code", "HEAD")
    source_files = list((ROOT / "src").rglob("*.rs")) + [Path(__file__), ROOT / "tools/reference/homology_operator.rs"]
    metadata = {"protocol": "cocycle-f2-operator-v1", "reference_sha": reference_head,
                "candidate_sha": git(ROOT, "rev-parse", "HEAD").strip(),
                "candidate_status": git(ROOT, "status", "--short"),
                "source_sha256": digest(source_files), "reference_sha256": digest((args.reference / "src/homology_operator").glob("*.py")),
                "worker_sha256": digest([args.worker]), "python": sys.version,
                "reference_python": args.reference_python, "platform": platform.platform(),
                "processor": platform.processor(), "seed": args.seed, "samples": args.samples,
                "engines": engines, "command": sys.argv, "started_unix": time.time()}
    if args.baseline_worker:
        metadata["baseline_worker_sha256"] = digest([args.baseline_worker])
    if args.native_dir:
        native_files = [p for p in args.native_dir.rglob("*homology_native*") if p.suffix in {".pyd", ".so"}]
        if not native_files:
            raise ValueError("no homology native binary found")
        metadata["native_binary_sha256"] = digest(native_files)
        metadata["native_files"] = [str(p) for p in native_files]
        metadata["reference_native_source_sha256"] = digest((args.reference / "native/src").rglob("*.rs"))
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    def call(engine, data, path, mode="solve", iterations=1):
        fixture_path = path.with_suffix(".fixture.json")
        fixture_path.write_text(json.dumps(data, sort_keys=True) + "\n")
        if engine in {"rust", "baseline"}:
            worker = args.worker if engine == "rust" else args.baseline_worker
            return invoke([str(worker), mode, str(iterations)], path, wire(data))
        command = [args.reference_python, str(Path(__file__)), "--reference-worker", "--reference", str(args.reference),
                   "--fixture", str(fixture_path), "--engine", engine, "--mode", mode, "--iterations", str(iterations)]
        if args.native_dir:
            command += ["--native-dir", str(args.native_dir)]
        return invoke(command, path)

    correctness = []
    cases = [(n, "general", args.seed + n * 100 + s) for n in range(6) for s in range(5)]
    cases += [(n, family, args.seed + n) for n in (8, 32, 63, 64, 65, 127, 128, 129)
              for family in ("general", "homology")]
    cases += [(n, "sparse", args.seed + n) for n in (511, 512, 513, 1025)]
    for index, (n, family, seed) in enumerate(cases):
        data = fixture(n, seed, family, n <= 5)
        outputs = {engine: call(engine, data, args.output / f"correct-{index}-{engine}.json") for engine in engines}
        for engine in engines[1:]:
            compare(outputs["rust"], outputs[engine])
        correctness.append({"n": n, "family": family, "seed": seed, "state": "Passed",
                            "independent_span": expected_span(data, outputs["rust"]),
                            "fixture_sha256": hashlib.sha256(wire(data).encode()).hexdigest()})
    family = invoke([str(args.worker), "family"], args.output / "family-rust.json")
    if family != {"ranks": [[1, 1, 0], [1, 0], [0]], "bars": [[0, 2, 1]]}:
        raise RuntimeError("finite family disagrees with independent triangle expectation")
    reference_family = invoke([args.reference_python, str(Path(__file__)), "--family-worker", "--reference", str(args.reference)], args.output / "family-reference.json")
    if reference_family != family:
        raise RuntimeError("finite family differs from the pinned reference")
    if args.oracle_python:
        oracle = invoke([args.oracle_python, str(Path(__file__)), "--gudhi-worker"], args.output / "gudhi.json")
        if oracle["betti_by_stage"] != [[3, 0], [1, 1], [1, 0, 0]] or oracle["h1"] != [[1., 2.]]:
            raise RuntimeError("GUDHI triangle comparison differs")
    else:
        oracle = {"state": "NotRun", "reason": "oracle Python was not supplied"}
    summary = {"correctness": correctness, "gudhi": oracle, "family": family,
               "exclusions": ["GUDHI has no projection, weighted representative or abstract-chain oracle",
                              "globally minimum-stretch projections, rational/float costs and Python snapshot restoration are not ported",
                              "finite final-stage continuation is unknown; reference constant extension is not adopted"],
               "measurements": []}
    if args.geometry:
        summary["geometry"] = geometry_checks(args)
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    external = len(engines) - 1 - int("baseline" in engines)
    print(f"Correctness: {len(correctness)} cases, {external} external comparisons per case, "
          f"baseline checked: {'baseline' in engines}, exact agreement.", flush=True)
    for n in map(int, args.sizes.split(",")) if args.samples else []:
        for family_name in ("general", "homology", "sparse"):
            data = fixture(n, args.seed + 7000 + n, family_name)
            cell = {"n": n, "family": family_name, "mode": "solve", "samples": {},
                    "fixture_sha256": hashlib.sha256(wire(data).encode()).hexdigest()}
            # One independent warmup process per engine is retained and discarded.
            for engine in engines:
                call(engine, data, args.output / f"warmup-{n}-{family_name}-{engine}.json")
                cell["samples"][engine] = []
            for sample in range(args.samples):
                order = engines.copy()
                random.Random(args.seed + sample // len(engines)).shuffle(order)
                offset = sample % len(engines)
                order = order[offset:] + order[:offset]
                results = {}
                for engine in order:
                    result = call(engine, data, args.output / f"bench-{n}-{family_name}-{sample}-{engine}.json")
                    results[engine] = result
                    cell["samples"][engine].append(result)
                for engine in engines[1:]:
                    compare(results["rust"], results[engine])
            cell["median_ns"] = {engine: statistics.median(r["elapsed_ns"] for r in values)
                                 for engine, values in cell["samples"].items()}
            cell["spread_ns"] = {engine: [min(r["elapsed_ns"] for r in values), max(r["elapsed_ns"] for r in values)]
                                 for engine, values in cell["samples"].items()}
            cell["median_peak_rss_kib"] = {engine: statistics.median([r["peak_rss_kib"] for r in values])
                                           if all(r["peak_rss_kib"] is not None for r in values) else None
                                           for engine, values in cell["samples"].items()}
            summary["measurements"].append(cell)
            (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
            print(f"Measured {family_name} n={n}: {cell['median_ns']}", flush=True)
    metadata["finished_unix"] = time.time()
    (args.output / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def parse():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reference", type=Path)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--baseline-worker", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--reference-python", default=sys.executable)
    parser.add_argument("--oracle-python")
    parser.add_argument("--native-dir", type=Path)
    parser.add_argument("--geometry", action="store_true", help="also compare T1--T5 queries and an independent raw-simplex oracle")
    parser.add_argument("--samples", type=int, default=0)
    parser.add_argument("--sizes", default="8,32,65")
    parser.add_argument("--seed", type=int, default=20261009)
    parser.add_argument("--reference-worker", action="store_true")
    parser.add_argument("--gudhi-worker", action="store_true")
    parser.add_argument("--family-worker", action="store_true")
    parser.add_argument("--fixture", type=Path)
    parser.add_argument("--engine", choices=("python", "native"), default="python")
    parser.add_argument("--mode", choices=("solve", "query"), default="solve")
    parser.add_argument("--iterations", type=int, default=1)
    args = parser.parse_args()
    if args.samples < 0 or args.iterations < 1:
        parser.error("sample count must be nonnegative and iterations positive")
    if not args.gudhi_worker and args.reference is None:
        parser.error("--reference is required")
    if not args.gudhi_worker and not args.reference_worker and not args.family_worker and (args.worker is None or args.output is None):
        parser.error("--worker and --output are required")
    return args


if __name__ == "__main__":
    options = parse()
    if options.gudhi_worker:
        gudhi_probe()
    elif options.family_worker:
        family_probe(options)
    elif options.reference_worker:
        reference(options)
    else:
        main(options)
