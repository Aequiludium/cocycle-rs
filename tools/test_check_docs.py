"""Regression checks for nested documentation discovery and relative links."""

from contextlib import redirect_stderr, redirect_stdout
from io import StringIO
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import check_docs


class DocumentationTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name).resolve()
        root_patch = patch.object(check_docs, "ROOT", self.root)
        root_patch.start()
        self.addCleanup(root_patch.stop)

    def write(self, name, text):
        path = self.root / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text, encoding="utf-8")

    def run_check(self):
        stdout, stderr = StringIO(), StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            status = check_docs.main()
        return status, stdout.getvalue(), stderr.getvalue()

    def test_nested_pages_resolve_parent_links_and_heading_fragments(self):
        self.write("docs/README.md", "[Guide](guides/rips.md)\n")
        self.write("docs/guides/rips.md",
                   "[Math](../reference/rips/mathematics.md#coverage)\n")
        self.write("docs/reference/rips/mathematics.md",
                   "# Coverage\n[Documentation](../../README.md)\n"
                   "[Again](#coverage)\n")
        status, stdout, stderr = self.run_check()
        self.assertEqual(status, 0, stderr)
        self.assertIn("3 Markdown files", stdout)

    def test_nested_errors_are_not_silently_skipped(self):
        self.write("docs/README.md", "# Documentation\n")
        self.write("docs/reference/rips/math.md",
                   "# Math\n[Missing](gone.md)\n"
                   "[Heading](../../README.md#missing)\n\u4e2d\u6587\n")
        status, _, stderr = self.run_check()
        self.assertEqual(status, 1)
        self.assertIn("missing local target gone.md", stderr)
        self.assertIn("missing heading ../../README.md#missing", stderr)
        self.assertIn("math.md:4: project documentation must be English", stderr)

    def test_code_examples_and_generated_artifacts_are_not_link_sources(self):
        self.write("docs/guides/rips.md",
                   "# Guide\n```text\n[Placeholder](missing.md)\n```\n"
                   "[External](https://example.invalid/not-fetched)\n")
        self.write("benches/results/run/notes.md", "[Old](missing.md)\n")
        self.write("target/generated.md", "[Generated](missing.md)\n")
        status, stdout, stderr = self.run_check()
        self.assertEqual(status, 0, stderr)
        self.assertIn("1 Markdown files", stdout)

    def test_benchmark_reports_and_native_guides_are_checked(self):
        self.write("benches/native/README.md", "[Missing](gone.cpp)\n")
        self.write("benches/reports/nested/report.md", "[Missing](gone.json)\n")
        self.write("benches/optimization/README.md", "[Missing](gone.rs)\n")
        status, _, stderr = self.run_check()
        self.assertEqual(status, 1)
        self.assertIn("benches/native/README.md", stderr)
        self.assertIn("benches/reports/nested/report.md", stderr)
        self.assertIn("benches/optimization/README.md", stderr)

    def test_reference_links_resolve_local_targets_and_headings(self):
        self.write("docs/reference/math.md", "# Coverage\n")
        self.write("docs/README.md",
                   '[Full][Math] [math][] [math]\n'
                   '[math]: reference/math.md#coverage "Contract"\n'
                   '[Remote][external]\n'
                   '[external]: https://example.invalid/not-fetched\n')
        status, _, stderr = self.run_check()
        self.assertEqual(status, 0, stderr)

    def test_reference_errors_include_original_line_numbers(self):
        self.write("docs/README.md",
                   '# Docs\n```text\n[Fake][missing]\n```\n'
                   '[Real][absent]\n[collapsed][]\n'
                   '[math]: missing.md\n[MATH]: #absent\n')
        status, _, stderr = self.run_check()
        self.assertEqual(status, 1)
        self.assertIn("README.md:5: undefined reference [absent]", stderr)
        self.assertIn("README.md:6: undefined reference [collapsed]", stderr)
        self.assertIn("README.md:7: missing local target missing.md", stderr)
        self.assertIn("README.md:8: duplicate reference definition [MATH]", stderr)
        self.assertIn("README.md:8: missing heading #absent", stderr)

    def test_code_spans_and_indented_longer_fences_are_not_links(self):
        self.write("docs/README.md",
                   '# Docs\n`[Inline](missing.md)`\n'
                   '``[Inline ` tick][missing]``\n'
                   '  ~~~text\n[Fake][missing]\n  ~~~~\n'
                   '[Real][absent]\n')
        status, _, stderr = self.run_check()
        self.assertEqual(status, 1)
        self.assertEqual(stderr.count("undefined reference"), 1)
        self.assertIn("README.md:7: undefined reference [absent]", stderr)
        self.assertNotIn("missing.md", stderr)

    def test_html_navigation_assets_and_encoded_fragments(self):
        self.write("site/index.html",
                   '<a href="guide/#a%20b">Guide</a><img src="art.svg">'
                   '<a href="https://example.invalid/">External</a>')
        self.write("site/guide/index.html",
                   '<a id="a b"></a><a href="../index.html">Home</a>'
                   '<a href="#a%20b">Here</a><script src="../code.js"></script>')
        self.write("site/art.svg", "<svg></svg>")
        self.write("site/code.js", "")
        self.assertEqual(check_docs.check_html(self.root / "site"), [])

    def test_html_missing_targets_anchors_and_site_escape(self):
        self.write("site/index.html",
                   '<a href="missing.html">Missing</a>\n'
                   '<a href="#absent">Anchor</a>\n'
                   '<a href="../outside.html">Escape</a>')
        self.write("outside.html", "<p>Outside</p>")
        errors = check_docs.check_html(self.root / "site")
        self.assertEqual(len(errors), 3)
        self.assertIn("index.html:1: missing HTML target missing.html", errors)
        self.assertIn("index.html:2: missing HTML anchor #absent", errors)
        self.assertIn("index.html:3: HTML target outside site ../outside.html", errors)

    def test_html_empty_or_missing_build_is_an_error(self):
        self.assertIn("No HTML pages", check_docs.check_html(self.root / "absent")[0])

    def test_angle_destinations_and_normalized_reference_labels(self):
        self.write("docs/math notes.md", "# Coverage\n")
        self.write("docs/README.md",
                   '[Inline](<math notes.md#coverage>)\n'
                   '[Reference][MaTh  Notes]\n'
                   '[math notes]: <math notes.md#coverage> "Title"\n')
        status, _, stderr = self.run_check()
        self.assertEqual(status, 0, stderr)

    def test_package_root_does_not_use_checkout_targets(self):
        self.write("README.md", "# Checkout\n")
        self.write("benches/report.md", "# Report\n")
        self.write("package/README.md", "[Report](benches/report.md)\n")
        stdout, stderr = StringIO(), StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            status = check_docs.main(self.root / "package")
        self.assertEqual(status, 1)
        self.assertIn("README.md:1: missing local target benches/report.md", stderr.getvalue())
        self.write("package/README.md", "[Outside](../benches/report.md)\n")
        with redirect_stdout(stdout), redirect_stderr(stderr):
            status = check_docs.main(self.root / "package")
        self.assertEqual(status, 1)
        self.assertIn("local target outside documentation root", stderr.getvalue())

    def test_package_root_checks_local_fragments_and_allows_repository_urls(self):
        self.write("package/README.md", "[Guide](docs/guide.md#coverage)\n"
                   "[Report](https://github.com/example/repo/blob/main/benches/report.md)\n")
        self.write("package/docs/guide.md", "# Coverage\n[Home](../README.md)\n")
        stdout, stderr = StringIO(), StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            status = check_docs.main(self.root / "package")
        self.assertEqual(status, 0, stderr.getvalue())
        self.assertIn("2 Markdown files", stdout.getvalue())

    def test_missing_or_empty_package_root_fails(self):
        stdout, stderr = StringIO(), StringIO()
        with redirect_stdout(stdout), redirect_stderr(stderr):
            self.assertEqual(check_docs.main(self.root / "missing"), 1)
            self.assertEqual(check_docs.main(self.root), 1)
        self.assertIn("not a directory", stderr.getvalue())
        self.assertIn("No Markdown documents", stderr.getvalue())


if __name__ == "__main__":
    unittest.main()
