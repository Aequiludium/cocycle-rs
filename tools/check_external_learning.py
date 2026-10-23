"""Run #92 graph and filtration exports through a separate scikit-learn process.

This optional, bounded application example is outside the Rust runtime. It
flattens supplied coordinates, explicit presence masks and graph adjacency;
it does not compute PH, spectra, distances or #111's persistent features.
"""

import argparse
import hashlib
import importlib.metadata
import itertools
import json
import math
from pathlib import Path
import platform
import random
import subprocess
import sys
import time


ROOT = Path(__file__).resolve().parents[1]
GRID = [-2., 0., 1.]
VERTICES = [10, 20, 35, 50, 70, 90]
LEARNING_PACKAGES = {"numpy": "2.2.6", "scipy": "1.15.3", "scikit-learn": "1.7.2",
                     "threadpoolctl": "3.6.0", "joblib": "1.6.0"}
ORACLE_PACKAGES = {"gudhi": "3.12.0"}


def package_versions(packages=LEARNING_PACKAGES):
    versions = {name: importlib.metadata.version(name) for name in packages}
    if versions != packages:
        raise ValueError("use the documented fixed optional package versions")
    return versions


def write_json(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8", newline="\n")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read_json(path):
    def invalid(value):
        raise ValueError(f"nonfinite JSON constant: {value}")
    return json.loads(path.read_text(encoding="utf-8"), parse_constant=invalid)


def fixtures(seed=92):
    """32 independent toy graphs; target is presence of edge (10,35).

    The other edges and coordinates are label-independent nuisance variables.
    Vertex 90 is an isolate born at 1; scale -2 is genuinely empty.
    """
    rng = random.Random(seed)
    samples, labels = [], []
    for order in range(32):
        name = f"graph-{order:02}"
        vertices = [[v, 1. if v == 90 else -1., rng.uniform(-1, 1),
                     rng.uniform(-1, 1)] for v in VERTICES]
        edges = [[10, 20, 0.], [20, 35, 0.]]
        if order % 2:
            edges.append([10, 35, 0.])
        if rng.random() < .5:
            edges.append([35, 50, 1.])
        if rng.random() < .5:
            edges.append([50, 70, 1.])
        sample = {"sample_key": name, "sample_order": order, "vertices": vertices,
                  "edges": edges, "grid": GRID, "universe": VERTICES}
        sample["snapshot"] = digest(json.dumps(sample, sort_keys=True).encode())
        samples.append(sample)
        labels.append({"sample_key": name, "label": order % 2,
                       "split": "test" if order >= 24 else "train"})
    rng.shuffle(labels)  # A positional label join must fail this example.
    return {"schema": "cocycle-learning-fixture-v1", "samples": samples,
            "labels": labels, "seed": seed}


def input_text(samples):
    lines = []
    for s in samples:
        lines.append(f'{s["sample_key"]} {s["snapshot"]} {s["sample_order"]} '
                     f'{len(s["vertices"])} {len(s["edges"])} {len(s["grid"])}')
        lines.append(" ".join(map(str, s["grid"])))
        lines.extend(" ".join(map(str, row)) for row in s["vertices"] + s["edges"])
    return "\n".join(lines) + "\n"


def unique(rows, key):
    result = {}
    for row in rows:
        identity = key(row)
        if identity in result:
            raise ValueError(f"duplicate key: {identity}")
        result[identity] = row
    return result


def record_key(record):
    return record["sample_key"], record["selection"]["kind"], record["selection"]["index"]


def sample_manifest(fixture):
    """Keep only ordered identities, encoding layout, labels and split metadata."""
    return {"format": "cocycle-graph-learning-manifest-v1",
            "samples": [{name: s[name] for name in ("sample_key", "snapshot", "sample_order", "universe", "grid")}
                        for s in fixture["samples"]], "labels": fixture["labels"]}


def validate_records(records, manifest):
    """Validate the self-contained wire contract without reading original data."""
    if manifest["format"] != "cocycle-graph-learning-manifest-v1":
        raise ValueError("unsupported sample manifest")
    samples = manifest["samples"]
    owners = unique(samples, lambda s: s["sample_key"])
    labels = unique(manifest["labels"], lambda s: s["sample_key"])
    if labels.keys() != owners.keys() or not owners:
        raise ValueError("labels must join every sample exactly once; empty dataset unsupported")
    if [s["sample_order"] for s in samples] != list(range(len(samples))):
        raise ValueError("sample order must be complete and explicit")
    for s in samples:
        universe, grid = s["universe"], s["grid"]
        if (any(type(v) is not int or v < 0 for v in universe)
                or universe != sorted(set(universe)) or not grid
                or any(type(x) not in (int, float) or not math.isfinite(x) for x in grid)
                or any(a >= b for a, b in zip(grid, grid[1:]))):
            raise ValueError("canonical vertex universe and finite increasing grid required")
        if universe != samples[0]["universe"] or grid != samples[0]["grid"]:
            raise ValueError("shared ordered vertex universe and grid required")
    for row in labels.values():
        if type(row["label"]) is not int or row["label"] not in (0, 1) or row["split"] not in ("train", "test"):
            raise ValueError("invalid binary label or split")
    expected = {(s["sample_key"], kind, index) for s in samples
                for kind, index in [("topology", None), *[("stage", i) for i in range(len(s["grid"]))]]}
    by_key = unique(records, record_key)
    if by_key.keys() != expected:
        raise ValueError("missing, foreign or unsupported sample/stage keys")
    for r in records:
        s = owners[r["sample_key"]]
        selection = r["selection"]
        is_stage = selection["kind"] == "stage"
        if is_stage and type(selection["index"]) is not int:
            raise ValueError("stage index must be an integer")
        scale = s["grid"][selection["index"]] if is_stage else None
        coverage = "inclusive-supplied-sublevel" if is_stage else "complete-supplied-source"
        if (r["format"] != "cocycle-graph-learning-v1" or r["snapshot"] != s["snapshot"]
                or r["sample_order"] != s["sample_order"] or selection["scale"] != scale
                or r["coverage"] != coverage or r["domain"] != "signed-integer"
                or r["scale_unit"] != "fixture-time"
                or r["topology"] != "supplied-graph-no-clique-expansion"):
            raise ValueError("source/selection/units mismatch")
        bases = r["bases"]
        vertex_keys = [tuple(v) for v in bases["0"]]
        edge_keys = [tuple(e) for e in bases["1"]]
        if (len(set(vertex_keys)) != len(vertex_keys) or len(set(edge_keys)) != len(edge_keys)
                or any(len(v) != 1 or type(v[0]) is not int or v[0] not in s["universe"] for v in vertex_keys)
                or any(len(e) != 2 or any(type(v) is not int for v in e) or e[0] >= e[1]
                       or any((v,) not in vertex_keys for v in e) for e in edge_keys)
                or bases.get("2") != [] or set(bases) != {"0", "1", "2"}):
            raise ValueError("basis/topology mismatch")
        boundary = r["boundary"]
        terms = unique(boundary["terms"], lambda t: tuple(t[:2]))
        expected_terms = {(vertex_keys.index((vertex,)), col): sign
                          for col, edge in enumerate(edge_keys)
                          for vertex, sign in zip(edge, (-1, 1))}
        if (boundary["shape"] != [len(vertex_keys), len(edge_keys)]
                or any(len(t) != 3 or any(type(x) is not int for x in t) for t in terms.values())
                or {key: t[2] for key, t in terms.items()} != expected_terms):
            raise ValueError("signed boundary/shape mismatch")
        features = r["features"]
        if (features["shape"] != [len(vertex_keys), 2] or features["row_keys"] != bases["0"]
                or features["meaning"] != "supplied-coordinate" or features["units"] != "dimensionless"
                or len(features["values"]) != len(vertex_keys) or any(len(row) != 2 for row in features["values"])
                or any(type(x) not in (int, float) or not math.isfinite(x) for row in features["values"] for x in row)):
            raise ValueError("feature row/shape/value mismatch")
    return by_key


def validate_export(records, fixture, gudhi=False):
    """Independent original-data oracle; never called by the model consumer."""
    by_key = validate_records(records, sample_manifest(fixture))
    owners = unique(fixture["samples"], lambda s: s["sample_key"])
    for r in records:
        s = owners[r["sample_key"]]
        scale = r["selection"]["scale"]
        wanted_vertices = {(v[0],) for v in s["vertices"] if scale is None or v[1] <= scale}
        wanted_edges = {tuple(e[:2]) for e in s["edges"] if scale is None or e[2] <= scale}
        if ({tuple(v) for v in r["bases"]["0"]} != wanted_vertices
                or {tuple(e) for e in r["bases"]["1"]} != wanted_edges):
            raise ValueError("export differs from original supplied topology/sublevel")
        supplied = {v[0]: v[2:] for v in s["vertices"]}
        if r["features"]["values"] != [supplied[v[0]] for v in r["bases"]["0"]]:
            raise ValueError("exported coordinates differ from original supplied values")
        if gudhi:
            import gudhi as gd
            tree = gd.SimplexTree()
            for v, birth, *_ in s["vertices"]:
                tree.insert([v], filtration=birth)
            for a, b, birth in s["edges"]:
                tree.insert([a, b], filtration=birth)
            actual = {tuple(cell) for cell, birth in tree.get_filtration() if scale is None or birth <= scale}
            if actual != wanted_vertices | wanted_edges:
                raise ValueError("GUDHI topology/sublevel disagreement")
    return by_key


def encode(record, universe):
    """Explicit application layout: presence, coordinates, upper adjacency.

    Absent vertices use zero coordinate slots WITH a separate presence mask.
    An absent edge is a computed adjacency zero, never a missing observation.
    """
    features = {v[0]: row for v, row in zip(record["features"]["row_keys"], record["features"]["values"])}
    edges = {tuple(e) for e in record["bases"]["1"]}
    return ([float(v in features) for v in universe]
            + [x for v in universe for x in features.get(v, [0., 0.])]
            + [float((a, b) in edges) for a, b in itertools.combinations(universe, 2)])


def route_inputs(records, manifest, route):
    by_key = validate_records(records, manifest)
    values, identities = [], []
    for s in manifest["samples"]:
        selections = [("topology", None)] if route == "topology" else [("stage", i) for i in range(len(s["grid"]))]
        if route not in ("topology", "stage"):
            raise ValueError("unknown route")
        rows = [by_key[(s["sample_key"], kind, index)] for kind, index in selections]
        values.append([x for r in rows for x in encode(r, s["universe"])])
        identities.append({"sample_key": s["sample_key"], "snapshot": s["snapshot"],
                           "route": route, "selections": [r["selection"] for r in rows]})
    if len({len(row) for row in values}) != 1:
        raise ValueError("ragged external model layout unsupported")
    return values, identities


def consume(input_path, manifest_path, output):
    """Actual separate-process fit/predict; third-party imports stay optional."""
    import numpy as np
    from sklearn.linear_model import LogisticRegression
    from sklearn.preprocessing import StandardScaler
    from threadpoolctl import threadpool_limits
    begin = time.perf_counter()
    versions = package_versions()
    manifest = read_json(manifest_path)
    records = [json.loads(line) for line in input_path.read_text(encoding="utf-8").splitlines()]
    labels = unique(manifest["labels"], lambda s: s["sample_key"])
    result = {"input_sha256": digest(input_path.read_bytes()), "manifest_sha256": digest(manifest_path.read_bytes()),
              "model": "StandardScaler + LogisticRegression(C=1, solver=lbfgs, max_iter=200)",
              "packages": versions, "python": sys.version,
              "controller_sha256": digest(Path(__file__).read_bytes()),
              "routes": {}, "predictions": []}
    with threadpool_limits(limits=1):
        for route in ("topology", "stage"):
            values, identities = route_inputs(records, manifest, route)
            train = [i for i, row in enumerate(identities) if labels[row["sample_key"]]["split"] == "train"]
            test = [i for i, row in enumerate(identities) if labels[row["sample_key"]]["split"] == "test"]
            y = np.array([labels[row["sample_key"]]["label"] for row in identities])
            if not train or not test or set(y[train].tolist()) != {0, 1}:
                raise ValueError("nonempty held-out split and two training classes required")
            x = np.array(values, dtype=float)
            fit_begin = time.perf_counter()
            scaler = StandardScaler().fit(x[train])
            model = LogisticRegression(C=1., solver="lbfgs", max_iter=200, random_state=92)
            model.fit(scaler.transform(x[train]), y[train])
            if model.n_iter_[0] >= 200:
                raise ValueError("external model did not converge")
            fit_seconds = time.perf_counter() - fit_begin
            predict_begin = time.perf_counter()
            probabilities = model.predict_proba(scaler.transform(x))[:, 1]
            predictions = model.predict(scaler.transform(x))
            result["routes"][route] = {"shape": list(x.shape), "train_count": len(train), "test_count": len(test),
                "accuracy": float(np.mean(predictions[test] == y[test])), "fit_seconds": fit_seconds,
                "predict_seconds": time.perf_counter() - predict_begin,
                "mean": scaler.mean_.tolist(), "scale": scaler.scale_.tolist(),
                "coefficients": model.coef_[0].tolist(), "intercept": float(model.intercept_[0])}
            result["predictions"].extend(dict(identity, probability=float(p), prediction=int(c))
                                         for identity, p, c in zip(identities, probabilities, predictions))
    result["consumer_seconds"] = time.perf_counter() - begin
    # Reverse wire order deliberately; callers must join by declared keys.
    result["predictions"].reverse()
    write_json(output, result)


def validate_predictions(result, records, manifest, input_path, manifest_path):
    if (result["packages"] != LEARNING_PACKAGES or result["controller_sha256"] != digest(Path(__file__).read_bytes())
            or result["input_sha256"] != digest(input_path.read_bytes())
            or result["manifest_sha256"] != digest(manifest_path.read_bytes())
            or set(result["routes"]) != {"topology", "stage"}):
        raise ValueError("external result input/route mismatch")
    by_key = unique(result["predictions"], lambda r: (r["route"], r["sample_key"]))
    expected_keys = {(route, s["sample_key"]) for route in ("topology", "stage") for s in manifest["samples"]}
    if by_key.keys() != expected_keys:
        raise ValueError("external output key mismatch")
    for route in ("topology", "stage"):
        values, identities = route_inputs(records, manifest, route)
        model = result["routes"][route]
        if model["shape"] != [len(values), len(values[0])]:
            raise ValueError("external shape mismatch")
        labels = unique(manifest["labels"], lambda s: s["sample_key"])
        train = [row for row, identity in zip(values, identities) if labels[identity["sample_key"]]["split"] == "train"]
        test_keys = [identity["sample_key"] for identity in identities if labels[identity["sample_key"]]["split"] == "test"]
        if model["train_count"] != len(train) or model["test_count"] != len(test_keys) or not train or not test_keys:
            raise ValueError("external split/count mismatch")
        for name in ("mean", "scale", "coefficients"):
            if len(model[name]) != len(values[0]) or any(not math.isfinite(x) for x in model[name]):
                raise ValueError("invalid external model parameters")
        if not math.isfinite(model["intercept"]) or any(x <= 0 for x in model["scale"]):
            raise ValueError("invalid external scaling")
        for column, (mean, scale) in enumerate(zip(model["mean"], model["scale"])):
            expected_mean = sum(row[column] for row in train) / len(train)
            variance = sum((row[column] - expected_mean)**2 for row in train) / len(train)
            expected_scale = math.sqrt(variance) if variance else 1.
            if abs(mean - expected_mean) > 1e-12 or abs(scale - expected_scale) > 1e-12:
                raise ValueError("normalization was not fit on training samples alone")
        for row, identity in zip(values, identities):
            actual = by_key[(route, identity["sample_key"])]
            if any(actual[name] != value for name, value in identity.items()):
                raise ValueError("external source/stage mismatch")
            score = model["intercept"] + sum((x - mean) / scale * coefficient for x, mean, scale, coefficient
                                             in zip(row, model["mean"], model["scale"], model["coefficients"]))
            probability = 1. / (1. + math.exp(-score)) if score >= 0 else math.exp(score) / (1. + math.exp(score))
            # float64 BLAS dot versus Python sum: tiny bounded 33/99-column case.
            if (not math.isfinite(actual["probability"]) or not 0 <= actual["probability"] <= 1
                    or abs(actual["probability"] - probability) > 1e-12
                    or type(actual["prediction"]) is not int or actual["prediction"] != int(score > 0)):
                raise ValueError("independent forward replay mismatch")
        expected_accuracy = sum(by_key[(route, key)]["prediction"] == labels[key]["label"] for key in test_keys) / len(test_keys)
        if model["accuracy"] != expected_accuracy:
            raise ValueError("external held-out metric mismatch")


def run_process(command, stem, stdin=None, timeout=120):
    begin = time.perf_counter()
    process = subprocess.Popen([str(x) for x in command], stdin=subprocess.PIPE,
                               stdout=subprocess.PIPE, stderr=subprocess.PIPE, cwd=ROOT, text=True)
    timed_out = False
    try:
        stdout, stderr = process.communicate(stdin, timeout=timeout)
    except subprocess.TimeoutExpired:
        timed_out = True
        process.kill()
        stdout, stderr = process.communicate()
    measurement = {"command": [str(x) for x in command], "wall_seconds": time.perf_counter() - begin,
                   "exit_code": process.returncode, "timeout": timed_out, "peak_working_set_bytes": None}
    if sys.platform == "win32":
        from benchmark_spectral_resources import peak_working_set
        measurement["peak_working_set_bytes"] = peak_working_set(process)
    stem.with_suffix(".stdout").write_text(stdout, encoding="utf-8", newline="\n")
    stem.with_suffix(".stderr").write_text(stderr, encoding="utf-8", newline="\n")
    write_json(stem.with_suffix(".process.json"), measurement)
    if timed_out or process.returncode:
        raise RuntimeError(f"external process failed; see {stem}.process.json and .stderr")
    return stdout, measurement


def run(worker, output, seed, python):
    output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    report = {"status": "running", "seed": seed, "candidate_sha": subprocess.check_output(
        ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip(),
        "base_sha": subprocess.check_output(["git", "merge-base", "HEAD", "origin/main"], cwd=ROOT, text=True).strip(),
        "working_tree": subprocess.check_output(["git", "status", "--porcelain"], cwd=ROOT, text=True),
        "platform": platform.platform(), "python": sys.version, "external_python": str(python),
        "external_python_sha256": digest(python.read_bytes()),
        "memory_scope": "Windows absolute child peak working set; null on other platforms",
        "timing_scope": "single serial smoke run, includes process startup; no performance ranking",
        "worker_sha256": digest(worker.read_bytes()), "source_sha256": {}}
    try:
        for path in [*ROOT.glob("src/**/*.rs"), ROOT / "Cargo.toml", ROOT / "Cargo.lock",
                     Path(__file__), ROOT / "tools/reference/learning_handoff.rs",
                     ROOT / "tools/benchmark_spectral_resources.py", ROOT / "tools/test_check_external_learning.py"]:
            report["source_sha256"][path.relative_to(ROOT).as_posix()] = digest(path.read_bytes())
        report["phase"] = "dependencies"
        report["oracle_packages"] = package_versions(ORACLE_PACKAGES)
        report["toolchains"] = {name: subprocess.check_output([name, "--version"], text=True).strip()
                                for name in ("cargo", "rustc")}
        fixture = fixtures(seed)
        fixture_path = output / "fixture.json"
        write_json(fixture_path, fixture)
        manifest = sample_manifest(fixture)
        manifest_path = output / "manifest.json"
        write_json(manifest_path, manifest)
        protocol = input_text(fixture["samples"])
        (output / "input.txt").write_text(protocol, encoding="utf-8", newline="\n")
        report["phase"] = "export"
        exported, report["export_process"] = run_process([worker], output / "export", protocol)
        input_path = output / "records.jsonl"
        input_path.write_text(exported, encoding="utf-8", newline="\n")
        records = [json.loads(line) for line in exported.splitlines()]
        validate_export(records, fixture, gudhi=True)
        report.update(records=len(records), samples=len(fixture["samples"]), export_bytes=len(exported.encode()))
        report["phase"] = "external_fit_predict"
        command = [python, Path(__file__).resolve(), "--consume", input_path,
                   "--manifest", manifest_path, "--result", output / "predictions.json"]
        _, report["consumer_process"] = run_process(command, output / "consumer")
        result = read_json(output / "predictions.json")
        validate_predictions(result, records, manifest, input_path, manifest_path)
        report["learning_packages"] = result["packages"]
        report.update(status="success", phase="validated", routes={route: {k: result["routes"][route][k]
            for k in ("shape", "train_count", "test_count", "accuracy", "fit_seconds", "predict_seconds")}
            for route in ("topology", "stage")}, oracle="independent signed incidence + GUDHI sublevel sets + scalar forward replay")
    except Exception as error:
        report.update(status="failure", error=f"{type(error).__name__}: {error}")
        raise
    finally:
        report["wall_seconds"] = time.perf_counter() - started
        report["artifact_sha256"] = {p.name: digest(p.read_bytes()) for p in output.iterdir() if p.is_file()}
        write_json(output / "report.json", report)
    print(json.dumps({"status": report["status"], "samples": report["samples"], "routes": report["routes"]}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--seed", type=int, default=92)
    parser.add_argument("--python", type=Path, default=Path(sys.executable))
    parser.add_argument("--consume", type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--result", type=Path)
    args = parser.parse_args()
    if args.consume:
        if not args.manifest or not args.result:
            parser.error("--consume requires --manifest and --result")
        consume(args.consume.resolve(), args.manifest.resolve(), args.result.resolve())
    else:
        if not args.worker or not args.output:
            parser.error("--worker and --output required")
        run(args.worker.resolve(), args.output.resolve(), args.seed, args.python.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
