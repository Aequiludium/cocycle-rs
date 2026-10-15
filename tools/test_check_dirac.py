"""Rust-fixture coverage and fixed-source integrity guards, without NumPy."""
from contextlib import redirect_stderr
from io import StringIO
import json
from pathlib import Path
import tempfile
import unittest

import check_dirac


class DiracControllerTests(unittest.TestCase):
    def rows(self):
        return [{"name": name, "basis": [], "offsets": [], "truncated": False,
                 "matrix": [], "x": [], "action": [], "square_action": []}
                for name in check_dirac.CASES]

    def write(self, root, rows):
        path = root / "rust.jsonl"
        path.write_text("\n".join(json.dumps(row) for row in rows) + "\n", encoding="utf-8")
        return path

    def test_complete_fixture_set_is_preserved(self):
        with tempfile.TemporaryDirectory() as directory:
            rows = self.rows()
            self.assertEqual(check_dirac.read_rows(self.write(Path(directory), rows)), rows)

    def test_powershell_utf8_bom_is_accepted(self):
        with tempfile.TemporaryDirectory() as directory:
            rows = self.rows()
            path = self.write(Path(directory), rows)
            path.write_bytes(b'\xef\xbb\xbf' + path.read_bytes())
            self.assertEqual(check_dirac.read_rows(path), rows)

    def test_partial_duplicate_and_unknown_fixtures_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            rows = self.rows()
            unknown = [{**row, "name": "unknown"} if i == 0 else row
                       for i, row in enumerate(rows)]
            for invalid in (rows[:-1], rows[:-1] + [rows[0]], unknown):
                with self.assertRaisesRegex(ValueError, "fourteen fixtures"):
                    check_dirac.read_rows(self.write(Path(directory), invalid))

    def test_incomplete_or_malformed_records_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            rows = self.rows()
            del rows[0]["matrix"]
            with self.assertRaisesRegex(ValueError, "incomplete"):
                check_dirac.read_rows(self.write(root, rows))
            path = root / "rust.jsonl"
            path.write_text('{"name":', encoding="utf-8")
            with self.assertRaises(json.JSONDecodeError):
                check_dirac.read_rows(path)

    def test_wrong_reference_blob_is_rejected_without_execution(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "simplicial_complex.py"
            path.write_text("raise RuntimeError('must never execute')\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "blob mismatch"):
                check_dirac.checked_toponetx_source(path)

    def test_failed_input_does_not_overwrite_previous_evidence(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            output = root / "result.json"
            output.write_text("existing evidence\n", encoding="utf-8")
            diagnostic = StringIO()
            with redirect_stderr(diagnostic):
                status = check_dirac.main(["--rust", str(root / "missing.jsonl"),
                                           "--output", str(output)])
            self.assertEqual(status, 1)
            self.assertIn("Dirac check failed", diagnostic.getvalue())
            self.assertEqual(output.read_text(encoding="utf-8"), "existing evidence\n")


if __name__ == "__main__":
    unittest.main()
