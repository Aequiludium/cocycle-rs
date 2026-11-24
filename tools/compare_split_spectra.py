"""Round-2 Rust parity against the hash-pinned #145 Round-1 transcription.

Consumes immutable Round-1 artifacts when supplied. Otherwise the pinned
reference generates and freezes a fresh platform-specific Round-1 replay first.
The Rust worker never supplies expected values. Dependencies remain external.
"""
import argparse
import ast
import contextlib
import hashlib
import importlib.util
import json
from pathlib import Path
import platform
import shutil
import subprocess
import sys
import time


REFERENCE_COMMIT = "8ee8f49b1aed55816456a0ef46eeec2501151235"
REFERENCE_SHA256 = "9d33727de04b19db2d595c75f420e69a4a0c0457d7e0574fd86ea4d02c9f3b62"
ROOT = Path(__file__).resolve().parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_reference(path):
    if digest(path) != REFERENCE_SHA256:
        raise ValueError("Round-1 reference source hash mismatch")
    spec = importlib.util.spec_from_file_location("split_round1", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def decode(value):
    import numpy as np
    if isinstance(value, dict):
        if set(value) == {"shape", "data"}:
            return np.array(value["data"], dtype=float).reshape(value["shape"])
        return {k: decode(v) for k, v in value.items()}
    if isinstance(value, list):
        return [decode(v) for v in value]
    return value


def read(path):
    return decode(json.loads(path.read_text(encoding="utf-8")))


def verify_frozen(root):
    environment = read(root / "environment.json")
    if environment["source_sha256"] != REFERENCE_SHA256:
        raise ValueError("frozen artifacts used a different Round-1 source")
    hashes = read(root / "artifact-hashes.json")
    for name, expected in hashes.items():
        if Path(name).name != name or "/" in name or "\\" in name:
            raise ValueError("unsafe artifact path")
        if digest(root / name) != expected:
            raise ValueError(f"frozen artifact was modified: {name}")
    return environment


def protocol(operation, sizes, *arrays):
    import numpy as np
    return " ".join([operation, *(str(x) for x in sizes),
                     *(repr(float(x)) for a in arrays for x in np.asarray(a).flat)]) + "\n"


def parse(text):
    result = {}
    for line in text.splitlines():
        key, value = line.split(" ", 1)
        if value in ("true", "false"):
            result[key] = value == "true"
        else:
            result[key] = ast.literal_eval(value)
    return result


def run(worker, reference, frozen, output):
    import numpy as np
    oracle = load_reference(reference)
    if frozen and output.resolve().is_relative_to(frozen.resolve()):
        raise ValueError("output must be outside the frozen directory")
    output.mkdir(parents=True, exist_ok=False)
    started = time.perf_counter()
    if frozen:
        frozen_environment = verify_frozen(frozen)
        shutil.copytree(frozen, output / "round1")
    else:
        with (output / "round1-run.log").open("w", encoding="utf-8") as log, contextlib.redirect_stdout(log):
            code = oracle.run(output / "round1", gudhi=True)
        if code:
            raise ValueError("pinned Round-1 replay failed; see retained artifacts")
        frozen_environment = verify_frozen(output / "round1")
    data = output / "round1"
    worker = worker.resolve()
    environment = dict(reference_commit=REFERENCE_COMMIT, reference_sha256=digest(reference),
                       controller_sha256=digest(Path(__file__)), comparison_python=sys.version,
                       comparison_numpy=np.__version__, platform=platform.platform(),
                       frozen_environment=frozen_environment, worker=str(worker),
                       worker_sha256=digest(worker), command=sys.argv,
                       head=subprocess.run(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True,
                                           capture_output=True, check=True).stdout.strip(),
                       dirty=subprocess.run(["git", "status", "--porcelain"], cwd=ROOT, text=True,
                                            capture_output=True, check=True).stdout.strip(),
                       sources={str(p.relative_to(ROOT)): digest(p) for p in sorted((ROOT / "src").rglob("*.rs"))},
                       tolerance=dict(matrix=1e-10, normalized_residual=1e-9, eigenvalues=1e-8),
                       basis_policy="compare lifted matrices and subspace projectors, not complement coordinates")
    environment["sources"]["examples/split_spectra.rs"] = digest(ROOT / "examples/split_spectra.rs")
    environment["rustc"] = subprocess.run(["rustc", "--version"], text=True, capture_output=True, check=True).stdout.strip()
    oracle.write_json(output / "environment.json", environment)
    results, failures = [], []

    def invoke(name, request, arguments=(), expected_error=None):
        (output / f"{name}-input.txt").write_text(request, encoding="utf-8")
        try:
            completed = subprocess.run([str(worker), "--oracle", *arguments], input=request,
                                       text=True, capture_output=True, check=False, timeout=120)
        except subprocess.TimeoutExpired as error:
            for suffix, content in [("stdout", error.stdout), ("stderr", error.stderr)]:
                (output / f"{name}-{suffix}.txt").write_bytes(
                    content.encode("utf-8") if isinstance(content, str) else content or b"")
            raise
        (output / f"{name}-stdout.txt").write_text(completed.stdout, encoding="utf-8")
        (output / f"{name}-stderr.txt").write_text(completed.stderr, encoding="utf-8")
        if expected_error is not None:
            if completed.returncode == 0 or expected_error not in completed.stderr:
                raise AssertionError(f"expected nonzero error containing {expected_error!r}; got {completed.stderr}")
            return dict(status="expected failure", exit_code=completed.returncode)
        if completed.returncode:
            raise AssertionError(f"unexpected worker failure: {completed.stderr}")
        return parse(completed.stdout)

    def record(name, check):
        try:
            result = check()
            oracle.write_json(output / f"{name}-comparison.json", result)
            results.append(dict(case=name, status="pass"))
        except Exception as error:
            failure = dict(case=name, error=str(error), exception=type(error).__name__)
            failures.append(failure)
            oracle.write_json(output / f"{name}-failure.json", failure)

    for path in sorted(data.glob("input-gram-*.json")):
        name = path.stem.removeprefix("input-")
        a = read(path)["A"]
        expected = read(data / f"{name}.json")
        def check_gram(name=name, a=a, expected=expected):
            actual = invoke(name, protocol("GRAM", a.shape, a))
            n, rank = a.shape[1], expected["numerical_rank"]
            for field, key in [("PRIMAL", "primal"), ("DUAL", "dual"),
                               ("PRIMAL_SPECTRUM", "primal_spectrum"), ("DUAL_SPECTRUM", "dual_spectrum"),
                               ("POSITIVE", "positive"), ("RECONSTRUCTED", "reconstructed")]:
                ref = np.asarray(expected[key])
                oracle.check_close(np.array(actual[field]).reshape(ref.shape), ref)
            vectors = np.array(actual["LIFTED"]).reshape(n, rank)
            reference_vectors = expected["lifted_vectors"]
            oracle.check_close(vectors @ vectors.T, reference_vectors @ reference_vectors.T)
            fresh = oracle.relative_error(expected["primal"] @ vectors, vectors * np.array(actual["POSITIVE"]))
            if fresh > 1e-9 or np.linalg.norm(vectors.T @ vectors - np.eye(rank)) > 1e-9:
                raise AssertionError("independent lifted Gram residual/orthogonality")
            if not np.isclose(actual["THRESHOLD"], expected["threshold"], rtol=1e-10, atol=1e-30):
                raise AssertionError("Gram threshold differs")
            if actual["RANK"] != expected["exact_rank"] or actual["RESIDUAL"] > 1e-9 or actual["ORTHOGONALITY"] > 1e-9:
                raise AssertionError("Gram rank/lifted verification")
            if actual["DUAL_IS_SMALLER"] != (expected["selected_side"] == "dual"):
                raise AssertionError("Gram side metadata differs")
            return dict(rank=actual["RANK"], residual=actual["RESIDUAL"],
                        rust_ambiguous=actual["AMBIGUOUS"], reference_ambiguous=expected["ambiguous"])
        record(name, check_gram)

    def restriction_case(name, h, kernel, expected):
        n, k = kernel.shape
        actual = invoke(name, protocol("RESTRICT", (n, k), h, kernel))
        x = np.array(actual["COMPLEMENT"]).reshape(n, n - k)
        reduced = np.array(actual["MATRIX"]).reshape(n - k, n - k)
        ref_x = expected["complement"]
        oracle.check_close(x @ x.T, ref_x @ ref_x.T)
        oracle.check_close(x @ reduced @ x.T, ref_x @ expected["reduced"] @ ref_x.T)
        oracle.check_close(np.array(actual["VALUES"]), np.array(expected["spectrum"]))
        vectors = np.array(actual["VECTORS"]).reshape(n, n - k)
        oracle.check_close(vectors @ vectors.T, expected["vectors"] @ expected["vectors"].T)
        fresh = oracle.relative_error(h @ vectors, vectors * np.array(actual["VALUES"]))
        if fresh > 1e-9 or np.linalg.norm(vectors.T @ vectors - np.eye(n - k)) > 1e-9:
            raise AssertionError("independent restriction residual/orthogonality")
        if not np.isclose(actual["THRESHOLD"], expected["threshold"], rtol=1e-10, atol=1e-30):
            raise AssertionError("restriction threshold differs")
        if actual["RESIDUAL"] > 1e-9 or actual["ORTHOGONALITY"] > 1e-9 or actual["KERNEL_RESIDUAL"] > 1e-9:
            raise AssertionError("original-domain restriction verification")
        return dict(residual=actual["RESIDUAL"], ambiguity=actual["AMBIGUOUS"])

    pairs = read(data / "pair-controls.json")
    for name, pair in pairs.items():
        def split_case(name=name, pair=pair):
            a, f = pair["A"], pair["F"]
            actual = invoke(f"split-{name}", protocol("SPLIT", (a.shape[0], a.shape[1], f.shape[1]), a, f))
            oracle.check_close(np.array(actual["MATRIX"]).reshape(pair["Delta"].shape), pair["Delta"])
            oracle.check_close(np.array(actual["VALUES"]), np.array(pair["union"]))
            return dict(values=actual["VALUES"], chain_residual=actual["CHAIN_RESIDUAL"])
        record(f"split-{name}", split_case)
        record(f"restrict-{name}", lambda name=name, pair=pair: restriction_case(
            f"restrict-{name}", pair["Delta"], pair["restriction"]["kernel"], pair["restriction"]))

    edge = read(data / "representative-counterexample.json")
    record("nonharmonic-representative", lambda: invoke("nonharmonic-representative", protocol(
        "RESTRICT", (2, 1), edge["Delta"], edge["representative"]), expected_error="not harmonic"))
    record("edge-harmonic", lambda: restriction_case("edge-harmonic", edge["Delta"],
                                                     edge["valid"]["kernel"], edge["valid"]))

    for path in sorted(data.glob("input-solver-*.json")):
        name = path.stem.removeprefix("input-")
        inputs, expected = read(path), read(data / f"{name}.json")
        def check_solver(name=name, inputs=inputs, expected=expected):
            h, initial, count = inputs["H"], inputs["initial"], inputs["count"]
            n, block = initial.shape
            request = protocol("SOLVE", (n, 0, block, count, 120, 1e-9), h, np.empty((n, 0)), initial)
            basic = expected["variants"]["paper_basic"]
            if basic["status"] == "failed":
                return invoke(name, request, expected_error="Rayleigh-Ritz")
            actual = invoke(name, request)
            oracle.check_close(np.array(actual["PAIR_VALUES"]), np.array(basic["answer"]["values"]), 1e-8)
            vectors = np.array(actual["PAIR_VECTORS"]).reshape(n, count)
            fresh = np.linalg.norm(h @ vectors - vectors * np.array(actual["PAIR_VALUES"]), axis=0) / max(1., np.linalg.norm(h, 2))
            if np.max(fresh, initial=0) > 1e-9 or np.linalg.norm(vectors.T @ vectors - np.eye(count)) > 1e-9:
                raise AssertionError("independent positive residual/orthogonality")
            oracle.check_close(np.array(actual["PAIR_RESIDUALS"]), fresh, 1e-12)
            if max(actual["PAIR_RESIDUALS"], default=0.) > 1e-9 or actual["CERTIFIED"]:
                raise AssertionError("positive residual or index certification")
            if actual["PAIR_ORTHOGONALITY"] > 1e-9:
                raise AssertionError("positive orthogonality")
            return dict(values=actual["PAIR_VALUES"], residuals=actual["PAIR_RESIDUALS"],
                        rust_iterations=actual["ITERATIONS"], reference_iterations=basic["answer"]["iterations"])
        record(name, check_solver)

    limits = read(data / "restriction-boundaries.json")
    for name, key in [("large-nullity", "large_zero_block"), ("tiny-positive", "tiny_positive")]:
        case = limits[key]
        h = (case["vectors"] * np.array(case["spectrum"])) @ case["vectors"].T
        h = (h + h.T) / 2
        record(name, lambda name=name, h=h, case=case: restriction_case(name, h, case["kernel"], case))
        if key == "tiny_positive":
            record("ambiguous-positive-request", lambda h=h, case=case: invoke("ambiguous-positive-request",
                protocol("SOLVE", (3, 1, 1, 1, 120, 1e-9), h, case["kernel"], case["vectors"][:, :1]),
                expected_error="ambiguous positive restriction"))
        else:
            initial = case["complement"] @ np.ones((3, 1))
            def restricted_solve(h=h, case=case, initial=initial):
                actual = invoke("large-nullity-solve", protocol("SOLVE", (35, 32, 1, 1, 120, 1e-9),
                                                                h, case["kernel"], initial))
                oracle.check_close(np.array(actual["PAIR_VALUES"]), np.array(limits["restricted_iteration"]["values"]))
                vectors = np.array(actual["PAIR_VECTORS"]).reshape(35, 1)
                fresh = np.linalg.norm(h @ vectors - vectors * np.array(actual["PAIR_VALUES"])) / max(1., np.linalg.norm(h, 2))
                leakage = np.linalg.norm(case["kernel"].T @ vectors)
                if fresh > 1e-9 or leakage > 1e-9 or np.linalg.norm(vectors.T @ vectors - np.eye(1)) > 1e-9:
                    raise AssertionError("independent lifted restricted solve")
                return dict(values=actual["PAIR_VALUES"], residuals=actual["PAIR_RESIDUALS"])
            record("large-nullity-solve", restricted_solve)

    diagonal = np.diag([0., 1.])
    record("incomplete-kernel", lambda: invoke("incomplete-kernel", protocol("RESTRICT", (2, 0),
        diagonal, np.empty((2, 0))), expected_error="incomplete kernel"))
    record("dependent-kernel", lambda: invoke("dependent-kernel", protocol("RESTRICT", (2, 2),
        diagonal, [[1., 1.], [0., 0.]]), expected_error="dependent harmonic basis"))
    excessive = protocol("SOLVE", (2, 1, 1, 2, 120, 1e-9), edge["Delta"], edge["valid"]["kernel"], [[1.], [-1.]])
    record("excessive-count", lambda: invoke("excessive-count", excessive, expected_error="count/block/tolerance"))
    exhaustion = protocol("SOLVE", (3, 0, 1, 1, 0, 1e-9), np.diag([1., 2., 3.]), np.empty((3, 0)), np.ones((3, 1)))
    record("iteration-exhaustion", lambda: invoke("iteration-exhaustion", exhaustion, expected_error="iteration limit"))
    simple = protocol("GRAM", (2, 2), np.eye(2))
    for name, args, error in [("work-limit", ["--max-work", "1"], "WorkLimitExceeded"),
                              ("storage-limit", ["--max-entries", "1"], "max_dense_entries"),
                              ("cancellation", ["--cancel"], "Cancelled")]:
        record(name, lambda name=name, args=args, error=error: invoke(name, simple, args, error))
    record("recovery", lambda: invoke("recovery", simple))
    summary = dict(passed=len(results), failures=failures, results=results,
                   seconds=time.perf_counter() - started, reference_commit=REFERENCE_COMMIT,
                   exclusions=["no Rust SVD search correction or robust Algorithm 8", "no optimization measurement",
                               "GUDHI corroborates only frozen torsion-free topology, not real spectra"])
    oracle.write_json(output / "summary.json", summary)
    print(json.dumps(dict(output=str(output), passed=len(results), failures=failures)))
    return int(bool(failures))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", type=Path, required=True)
    parser.add_argument("--reference", type=Path, required=True, help="hash-pinned Round-1 Python source")
    parser.add_argument("--frozen", type=Path, help="existing immutable Round-1 artifact directory")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    return run(args.worker, args.reference, args.frozen, args.output)


if __name__ == "__main__":
    raise SystemExit(main())
