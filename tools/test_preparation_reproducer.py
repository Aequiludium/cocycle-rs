"""Integrity checks for the counterexample and isolated source transformations."""

import io
from pathlib import Path
import struct
import tarfile
import tempfile
import unittest

import reproduce_distance_preparation as repro


class PreparationReproducerTests(unittest.TestCase):
    def test_published_counterexample(self):
        payload = repro.fixture(512)
        self.assertEqual(repro.sha(payload),
                         "2d5b79549d03532405601cbc51d380f09f7b6e1f82cff0134f7230406643d99a")
        self.assertEqual(payload[:8], b"COCDST1\0")
        self.assertEqual(struct.unpack("<QQ", payload[8:24]), (512, 64))
        points = list(struct.iter_unpack("<dd", payload[24:]))
        self.assertEqual(len(points), 576)
        self.assertTrue(all(birth < death for birth, death in points))

    def test_patch_rejects_missing_or_ambiguous_marker(self):
        for source in ("nothing", "marker marker"):
            with self.assertRaises(ValueError):
                repro.replace_once(source, "marker", "replacement")

    def test_archive_rejects_escape_and_link(self):
        for name, kind in (("../escape", tarfile.REGTYPE), ("link", tarfile.SYMTYPE)):
            with self.subTest(name=name), tempfile.TemporaryDirectory() as temporary:
                buffer = io.BytesIO()
                with tarfile.open(fileobj=buffer, mode="w") as archive:
                    member = tarfile.TarInfo(name)
                    member.type = kind
                    member.linkname = "../escape"
                    archive.addfile(member)
                with self.assertRaises(ValueError):
                    repro.extract_snapshot(buffer.getvalue(), Path(temporary) / "source")
                self.assertFalse((Path(temporary) / "escape").exists())

    def test_snapshot_instrumentation_preserves_original(self):
        archive = repro.git("archive", "HEAD", "src")
        template = (repro.ROOT / repro.HARNESS[2]).read_text(encoding="utf-8")
        original_hashes = repro.inventory(repro.ROOT / "src")
        for growth in (False, True):
            with self.subTest(growth=growth), tempfile.TemporaryDirectory() as temporary:
                snapshot = Path(temporary) / "source"
                repro.extract_snapshot(archive, snapshot)
                repro.instrument(snapshot, template, growth)
                text = (snapshot / "src/diagram_distances/wasserstein.rs").read_text()
                self.assertIn("#[inline(never)]\nfn solve_prepared", text)
                self.assertEqual("fn trace_edge_push" in text, growth)
                with self.assertRaises(ValueError):
                    repro.instrument(snapshot, template, growth)
        self.assertEqual(repro.inventory(repro.ROOT / "src"), original_hashes)


if __name__ == "__main__":
    unittest.main()
