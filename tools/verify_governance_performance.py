"""Run the existing integrated pipeline protocol in a fresh artifact directory."""
import json
import os
from pathlib import Path
import platform
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent.parent


def main():
    if platform.system() != "Linux":
        print("The existing pipeline protocol requires Linux /proc and resource limits.")
        return 78
    sources = {name: os.environ.get(name) for name in (
        "GOVERNANCE_GUDHI_SOURCE", "GOVERNANCE_RIPSER_SOURCE")}
    if not all(value and Path(value).is_dir() for value in sources.values()):
        print("Configure the pinned GUDHI and Ripser source directories from benches/native/README.md.")
        return 78
    parent = ROOT / "target/governance-verification"
    parent.mkdir(parents=True, exist_ok=True)
    output = Path(tempfile.mkdtemp(prefix="integrated-", dir=parent)) / "pipeline"
    command = [sys.executable, str(ROOT / "tools/benchmark_rips_pipeline.py"),
               "--quick", "--samples", "12", "--output", str(output),
               "--gudhi-source", sources["GOVERNANCE_GUDHI_SOURCE"],
               "--ripser-source", sources["GOVERNANCE_RIPSER_SOURCE"]]
    if os.environ.get("GOVERNANCE_BOOST_INCLUDE"):
        command += ["--boost-include", os.environ["GOVERNANCE_BOOST_INCLUDE"]]
    print("Executing:", command, flush=True)
    print("Evidence:", output, flush=True)
    started = time.perf_counter()
    result = subprocess.run(command, cwd=ROOT, check=False)
    elapsed = time.perf_counter() - started
    summary = {}
    if (output / "summary.json").is_file():
        summary = json.loads((output / "summary.json").read_text())
    complete = (result.returncode == 0 and summary.get("status") == "passed"
                and summary.get("sources_unchanged") is True
                and summary.get("cases", 0) > 0
                and summary.get("worker_measurements", 0) >= summary.get("cases", 0)
                and summary.get("samples_per_worker") == 12
                and summary.get("harness_commit") == summary.get("kernel_commit"))
    (output.parent / "governance-run.json").write_text(json.dumps({
        "status": "passed" if complete else "failed", "command": command,
        "pipeline_directory": str(output), "returncode": result.returncode,
        "controller_wall_seconds": elapsed, "pipeline_summary": summary,
        "controller_scope": "Full protocol invocation including build, fixtures, parsing, warmups, workers, independent comparison and artifact output; not a per-backend speed ranking.",
        "worker_scope": "The unchanged pipeline's declared input/construction/expansion/compute/export bins; excludes fixture I/O and independent controller validation.",
    }, indent=2) + "\n")
    return 0 if complete else result.returncode or 1


if __name__ == "__main__":
    raise SystemExit(main())
