"""Protocol and provenance checks for the optional Rust parity controller."""
import importlib.util
import json
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import compare_split_spectra as parity


class ReferenceTests(unittest.TestCase):
    def test_changed_round1_source_is_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            source = Path(directory) / "reference.py"
            source.write_text("raise RuntimeError('must not execute')\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "source hash mismatch"):
                parity.load_reference(source)

    def test_protocol_parser_retains_empty_arrays_and_false_coverage(self):
        parsed = parity.parse("VALUES []\nCERTIFIED false\nITERATIONS 0\nRESIDUAL 1e-14\n")
        self.assertEqual(parsed["VALUES"], [])
        self.assertFalse(parsed["CERTIFIED"])
        self.assertEqual(parsed["RESIDUAL"], 1e-14)


@unittest.skipUnless(importlib.util.find_spec("numpy"), "optional NumPy reference environment")
class FrozenArtifactTests(unittest.TestCase):
    def test_empty_matrix_shape_is_not_lost(self):
        result = parity.decode({"shape": [0, 7], "data": []})
        self.assertEqual(result.shape, (0, 7))

    def test_modified_frozen_input_and_path_escape_are_rejected(self):
        with tempfile.TemporaryDirectory() as directory:
            root = Path(directory)
            (root / "environment.json").write_text(json.dumps({"source_sha256": parity.REFERENCE_SHA256}), encoding="utf-8")
            (root / "input.json").write_text("{}\n", encoding="utf-8")
            hashes = {"input.json": parity.digest(root / "input.json")}
            (root / "artifact-hashes.json").write_text(json.dumps(hashes), encoding="utf-8")
            self.assertEqual(parity.verify_frozen(root)["source_sha256"], parity.REFERENCE_SHA256)
            (root / "input.json").write_text("changed", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "was modified"):
                parity.verify_frozen(root)
            (root / "artifact-hashes.json").write_text(json.dumps({"../outside": "ignored"}), encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "unsafe artifact path"):
                parity.verify_frozen(root)

    def test_existing_output_is_never_overwritten(self):
        with tempfile.TemporaryDirectory() as directory, patch.object(parity, "load_reference"):
            with self.assertRaises(FileExistsError):
                parity.run(Path("unused-worker"), Path("unused-reference"), None, Path(directory))
            with self.assertRaisesRegex(ValueError, "outside the frozen directory"):
                parity.run(Path("unused-worker"), Path("unused-reference"), Path(directory), Path(directory) / "nested")


if __name__ == "__main__":
    unittest.main()
