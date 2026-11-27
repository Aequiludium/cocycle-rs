"""An empty or failed protocol must not become completed feature evidence."""
import json
from pathlib import Path
import tempfile
from types import SimpleNamespace
import unittest
from unittest.mock import patch

import verify_governance_performance as verify


class GovernancePerformanceTests(unittest.TestCase):
    def test_missing_and_failed_summary_cannot_pass_on_exit_zero(self):
        for summary in (None, {"status": "failed"}, {"status": "passed", "cases": 0}):
            with self.subTest(summary=summary), tempfile.TemporaryDirectory() as temporary:
                root = Path(temporary)

                def run(command, **_options):
                    output = Path(command[command.index("--output") + 1])
                    output.mkdir()
                    if summary is not None:
                        (output / "summary.json").write_text(json.dumps(summary))
                    return SimpleNamespace(returncode=0)

                with patch.object(verify, "ROOT", root), patch.object(verify.platform, "system", return_value="Linux"), patch.dict(verify.os.environ, {
                    "GOVERNANCE_GUDHI_SOURCE": temporary, "GOVERNANCE_RIPSER_SOURCE": temporary,
                }), patch.object(verify.subprocess, "run", side_effect=run):
                    self.assertEqual(verify.main(), 1)
                record = next(root.glob("target/governance-verification/*/governance-run.json"))
                self.assertEqual(json.loads(record.read_text())["status"], "failed")

    def test_unsupported_platform_is_environment_incomplete(self):
        with patch.object(verify.platform, "system", return_value="Windows"):
            self.assertEqual(verify.main(), 78)


if __name__ == "__main__":
    unittest.main()
