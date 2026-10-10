"""Profile private Rips experiments with separate diagnostic and latency builds.

Only standard-library Python is needed. Inputs are COCYCLE1 distance/H1 binary
fixtures. Two-pass initialization and virtual pairs can be measured separately
or together. Test instrumentation is enabled, so these timings are diagnostic
and must not replace the public-API cross-library benchmark. The optional H2
paired mode compares already validated pipeline release workers, with no test
instrumentation; it performs no compilation while timing.
"""

import argparse
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import shutil
import statistics
import subprocess
import tempfile
from datetime import datetime, timezone


def paired_summary(baseline, candidate, seed):
    """Balanced rounds retain pairing when resampling the ratio of medians."""
    if not baseline or len(baseline) != len(candidate) or any(
            not math.isfinite(v) or v <= 0 for v in [*baseline, *candidate]):
        raise ValueError("equal nonempty positive sample sequences required")
    rng = random.Random(seed)
    ratios = []
    for _ in range(2000):
        indices = [rng.randrange(len(baseline)) for _ in baseline]
        ratios.append(statistics.median(candidate[i] for i in indices)
                      / statistics.median(baseline[i] for i in indices))
    ratios.sort()
    return {"baseline_median_ms": statistics.median(baseline),
            "candidate_median_ms": statistics.median(candidate),
            "median_ratio": statistics.median(candidate) / statistics.median(baseline),
            "paired_bootstrap_lower95": ratios[49], "paired_bootstrap_upper95": ratios[1949]}


def paired_h2(args):
    from benchmark_rips_pipeline import PROTOCOL, phase3_cases, schedule, worker, compare_samples, provenance
    from build_native import sha256
    identity = provenance("HEAD")
    if identity["dirty"]:
        raise ValueError("commit the harness and implementation before measurement")
    directories = {"baseline": args.h2_baseline.resolve(), "candidate": args.h2_candidate.resolve()}
    environments = {name: json.loads((path / "environment.json").read_text())
                    for name, path in directories.items()}
    for name, path in directories.items():
        env = environments[name]
        if env["dirty"] or env["protocol_id"] != PROTOCOL or env["quick"] != args.quick:
            raise ValueError("clean matched pipeline protocol and size required")
        if json.loads((path / "summary.json").read_text())["status"] != "passed":
            raise ValueError("both source suites must have passed validation")
        if sha256(path / "build/cocycle") != env["binaries_sha256"]["cocycle"]:
            raise ValueError("frozen executable hash mismatch")
    args.output.mkdir(parents=True, exist_ok=False)
    metadata = {**identity, "controller_sha256": sha256(Path(__file__)),
                "started_utc": datetime.now(timezone.utc).isoformat(), "cpu": args.cpu,
                "samples": args.samples, "bootstrap_replicates": 2000, "order_seed": 20261004,
                "original_environments": environments, "quick": args.quick,
                "timing_scope": "unchanged pipeline v2 elapsed_ms; fresh serial release processes",
                "frequency_and_load_control": "not controlled",
                "h1_gate": {"median_ratio_max": 1.05, "upper95_max": 1.10},
                "h2_regression_screen": {"median_ratio_max": 1.10}}
    (args.output / "controller.py").write_bytes(Path(__file__).read_bytes())
    records = []
    def save():
        (args.output / "results.json").write_text(json.dumps(records, indent=2) + "\n")
    # Linux localfs avoids demand paging from DrvFS for these short workers.
    # Originals and hashes stay in the immutable build directories.
    with tempfile.TemporaryDirectory(prefix="cocycle-h2-", dir="/tmp") as local:
        binaries = {}
        for name, path in directories.items():
            binaries[name] = Path(local) / name
            shutil.copy2(path / "build/cocycle", binaries[name])
            assert sha256(binaries[name]) == environments[name]["binaries_sha256"]["cocycle"]
        metadata["executable_filesystems"] = subprocess.check_output(["df", "-T", local], text=True)
        (args.output / "environment.json").write_text(json.dumps(metadata, indent=2) + "\n")
        for index, case in enumerate(phase3_cases(args.quick)):
            fixture = directories["baseline"] / "fixtures" / (case.name + ".txt")
            if sha256(fixture) != sha256(directories["candidate"] / "fixtures" / fixture.name):
                raise ValueError("matched fixture hash required")
            record = {"case": case.name, "dimension": case.fixture.q, "fixture_sha256": sha256(fixture),
                      "samples": [], "validation": "pending"}
            records.append(record)
            anchor = None
            try:
                for entry in schedule(list(binaries), args.samples, 20261004 + index):
                    name = entry["backend"]
                    sample = worker([str(binaries[name]), str(fixture), case.path, case.layout,
                                     str(case.epsilon), "no"], case, 30, 2048, args.cpu)
                    record["samples"].append({**entry, **sample})
                    if sample["status"] != "completed":
                        raise ValueError("paired worker failed")
                    if anchor is None:
                        anchor = sample
                    compare_samples(anchor, sample, case)
                record["validation"] = "passed"
            except (ValueError, subprocess.SubprocessError) as error:
                record.update(validation="failed", error=str(error))
                save()
                (args.output / "summary.json").write_text(json.dumps({"status": "failed", "error": str(error)}))
                raise
            groups = {name: [s for s in record["samples"] if s["backend"] == name and not s["warmup"]]
                      for name in binaries}
            record.update(paired_summary([s["elapsed_ms"] for s in groups["baseline"]],
                                         [s["elapsed_ms"] for s in groups["candidate"]], 20261004 + index))
            record["peak_rss_kib"] = {name: max(s["peak_rss_kib"] for s in samples) for name, samples in groups.items()}
            record["passed"] = record["median_ratio"] <= (1.05 if case.fixture.q == 1 else 1.10) and (
                case.fixture.q != 1 or record["paired_bootstrap_upper95"] <= 1.10)
            save()
            print(case.name, record["median_ratio"], record["passed"], flush=True)
    summary = {"status": "passed" if all(r["passed"] for r in records) else "failed",
               "cases": len(records), "processes": sum(len(r["samples"]) for r in records),
               "finished_utc": datetime.now(timezone.utc).isoformat(),
               "regressions": [r["case"] for r in records if not r["passed"]]}
    final_identity = provenance("HEAD")
    if final_identity != identity or sha256(Path(__file__)) != metadata["controller_sha256"]:
        summary.update(status="failed", error="source identity changed during measurement")
    (args.output / "summary.json").write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("fixtures", type=Path, nargs="*")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cpu", type=int)
    parser.add_argument("--h2-baseline", type=Path)
    parser.add_argument("--h2-candidate", type=Path)
    parser.add_argument("--quick", action="store_true")
    parser.add_argument("--samples", type=int, default=96)
    args = parser.parse_args()
    if args.output.exists():
        parser.error("output must be new")
    if args.cpu is not None:
        os.sched_setaffinity(0, {args.cpu})
    if args.h2_baseline or args.h2_candidate:
        if not args.h2_baseline or not args.h2_candidate or args.fixtures or args.samples < 1:
            parser.error("paired H2 mode requires both validated runs, positive samples and no binary fixtures")
        paired_h2(args)
        return
    if not args.fixtures:
        parser.error("H1 profiling requires binary fixtures")
    repo = Path(__file__).resolve().parents[1]
    build = subprocess.check_output(
        ["cargo", "test", "--locked", "--offline", "--release", "--lib", "--no-run",
         "--message-format=json"], cwd=repo, text=True,
    )
    artifacts = [json.loads(line) for line in build.splitlines()]
    executable, = [x["executable"] for x in artifacts
                   if x.get("reason") == "compiler-artifact" and x.get("executable")]
    digest = hashlib.sha256()
    for path in sorted((repo / "src").rglob("*.rs")):
        digest.update(str(path.relative_to(repo)).encode())
        digest.update(path.read_bytes())
    record = {
        "source_sha256": digest.hexdigest(),
        "rustc": subprocess.check_output(["rustc", "-Vv"], text=True),
        "platform": platform.platform(), "cpu": args.cpu,
        "protocol": "release lib tests; instrumentation; one warmup, five samples; fresh process per stage; reference checked after measurement",
        "results": [],
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    # Failures leave a partial diagnostic record; completion is explicit.
    for fixture in args.fixtures:
        for stage in ["explicit", "clearing", "implicit", "cone", "apparent",
                      "two-pass", "emergent", "virtual-two-pass", "virtual"]:
            print(f"{fixture.stem}: {stage}", flush=True)
            env = dict(os.environ, COCYCLE_ABLATION_FIXTURE=str(fixture.resolve()),
                       COCYCLE_ABLATION_STAGE=stage)
            output = subprocess.check_output(
                [executable, "persistence::flag::cohomology::profiling::profile_stage",
                 "--exact", "--ignored", "--nocapture"],
                cwd=repo, env=env, text=True, timeout=180,
            )
            row, = [json.loads(line.split("ABLATION ", 1)[1])
                    for line in output.splitlines() if "ABLATION " in line]
            row.update(fixture=fixture.name, fixture_sha256=hashlib.sha256(fixture.read_bytes()).hexdigest())
            row["median_ms"] = statistics.median(row["samples_ms"])
            record["results"].append(row)
            args.output.write_text(json.dumps(record, indent=2) + "\n")
    record["status"] = "passed"
    args.output.write_text(json.dumps(record, indent=2) + "\n")


if __name__ == "__main__":
    main()
