"""Exercise executable-worker freezing through the actual oracle controller."""

import importlib.util
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest


@unittest.skipIf(os.name == "nt", "POSIX executable permission regression")
@unittest.skipUnless(importlib.util.find_spec("numpy"), "optional NumPy oracle dependency")
class FrozenWorkerTests(unittest.TestCase):
    def test_frozen_worker_is_executable_and_retains_failure_output(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            worker = root / "worker"
            worker.write_text(
                f"#!{sys.executable}\nimport sys\n"
                "print('worker executed', sys.argv[1])\n"
                "sys.exit(37)\n", encoding="utf-8")
            worker.chmod(0o700)
            output = root / "comparison"
            controller = Path(__file__).with_name("compare_persistent_laplacians.py")
            run = subprocess.run(
                [sys.executable, str(controller), "--worker", str(worker),
                 "--output", str(output), "--components"],
                capture_output=True, text=True)
            self.assertNotEqual(run.returncode, 0)
            self.assertIn("exit status 37", run.stderr)
            self.assertEqual((output / "worker-stdout.txt").read_text().strip(),
                             "worker executed --components")
            self.assertTrue(os.access(output / "worker", os.X_OK))
            self.assertEqual((output / "worker").read_bytes(), worker.read_bytes())


if __name__ == "__main__":
    unittest.main()
