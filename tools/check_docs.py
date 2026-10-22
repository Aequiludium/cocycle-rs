"""Check maintained Markdown or a built HTML site's local links.

Uses only the standard library. External URLs are not fetched; dated raw artifacts
and generated output are not documentation sources. Run from any directory.
Markdown checks also reject CJK project text; --html checks rendered site assets
and anchors instead.
"""

import argparse
from html.parser import HTMLParser
from pathlib import Path
import re
import sys
from urllib.parse import unquote, urlsplit

ROOT = Path(__file__).resolve().parents[1]
PATTERNS = ("*.md", "docs/**/*.md", "benches/*.md", "benches/native/**/*.md",
            "benches/reports/**/*.md", "benches/pipeline/**/*.md", "benches/distances/**/*.md",
            "benches/optimization/**/*.md",
            "tools/*.md", "assets/*.md", ".github/**/*.md")
DESTINATION = r'(<[^>\n]+>|[^\s)]+)(?:\s+"[^"\n]*")?'
LINK = re.compile(r"!?\[[^\]\n]*\]\(" + DESTINATION + r"\)")
DEFINITION = re.compile(r"^ {0,3}\[([^\]\n]+)\]:\s*" + DESTINATION + r"\s*$", re.MULTILINE)
REFERENCE = re.compile(r"!?\[([^\]\n]+)\]\[([^\]\n]*)\]")
CJK = re.compile(r"[\u3400-\u4dbf\u4e00-\u9fff]")


def prose(text):
    """Mask fenced code while preserving offsets for diagnostic line numbers."""
    result = []
    fence = None
    for line in text.splitlines(keepends=True):
        marker = re.match(r"^ {0,3}(`{3,}|~{3,})(.*)$", line)
        if fence is not None:
            result.append(blank(line))
            if (marker and marker[1][0] == fence[0] and len(marker[1]) >= len(fence)
                    and not marker[2].strip()):
                fence = None
        elif marker:
            fence = marker[1]
            result.append(blank(line))
        else:
            result.append(line)
    return "".join(result)


def blank(text):
    return re.sub(r"[^\n]", " ", text)


def reference_id(label):
    return " ".join(label.split()).casefold()


def links(text):
    """Yield (offset, destination, error) for our documented Markdown subset.

    Single-line inline links and reference definitions support quoted titles and
    angle-bracket destinations. Explicit references must resolve; bare brackets
    are shortcut links only when their label has a definition.
    """
    text = prose(text)
    text = re.sub(r"(?<!`)(`+)(?!`).*?(?<!`)\1(?!`)",
                  lambda match: blank(match[0]), text, flags=re.DOTALL)
    definitions = {}
    for match in DEFINITION.finditer(text):
        label = reference_id(match[1])
        if label in definitions:
            yield match.start(), None, f"duplicate reference definition [{match[1]}]"
        definitions[label] = match[2]
        yield match.start(), match[2], None
    text = DEFINITION.sub(lambda match: blank(match[0]), text)
    for match in LINK.finditer(text):
        yield match.start(), match[1], None
    text = LINK.sub(lambda match: blank(match[0]), text)
    for match in REFERENCE.finditer(text):
        label = match[2] or match[1]
        if reference_id(label) not in definitions:
            yield match.start(), None, f"undefined reference [{label}]"
    # Destinations were checked at their definitions, including unused ones.
    # Undefined bare brackets are ordinary prose, not necessarily broken links.


def anchors(text):
    used = {}
    result = set(re.findall(r'<a\s+(?:id|name)=["\']([^"\']+)', text))
    for heading in re.findall(r"^#{1,6}\s+(.+?)\s*#*\s*$", prose(text), re.MULTILINE):
        slug = re.sub(r"[^\w\- ]", "", heading.lower()).replace(" ", "-")
        count = used.get(slug, 0)
        used[slug] = count + 1
        result.add(f"{slug}-{count}" if count else slug)
    return result


def check(path, root=None):
    root = ROOT if root is None else root
    text = path.read_text(encoding="utf-8")
    errors = []
    for line, content in enumerate(text.splitlines(), 1):
        if CJK.search(content):
            errors.append(f"{path.relative_to(root).as_posix()}:{line}: project documentation must be English")
    for offset, destination, error in links(text):
        location = f"{path.relative_to(root).as_posix()}:{text.count(chr(10), 0, offset) + 1}"
        if error:
            errors.append(f"{location}: {error}")
            continue
        link = urlsplit(destination.strip("<>"))
        if link.scheme or link.netloc:
            continue
        target = (path.parent / unquote(link.path)).resolve() if link.path else path
        if not target.is_relative_to(root):
            errors.append(f"{location}: local target outside documentation root {destination}")
        elif not target.exists():
            errors.append(f"{location}: missing local target {destination}")
        elif link.fragment and target.suffix == ".md":
            if unquote(link.fragment) not in anchors(target.read_text(encoding="utf-8")):
                errors.append(f"{location}: missing heading {destination}")
    return errors


def main(root=None):
    root = (ROOT if root is None else root).resolve()
    if not root.is_dir():
        print(f"Documentation root is not a directory: {root}", file=sys.stderr)
        return 1
    paths = sorted({p for pattern in PATTERNS for p in root.glob(pattern) if p.is_file()})
    if not paths:
        print(f"No Markdown documents found in {root}", file=sys.stderr)
        return 1
    errors = [error for path in paths for error in check(path, root)]
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print(f"Documentation checks passed: {len(paths)} Markdown files; local links and CJK text scan.")
    return 0


class HTMLLinks(HTMLParser):
    """Collect local navigation, assets and anchors from a built page."""

    def __init__(self, text):
        super().__init__(convert_charrefs=True)
        self.ids = set()
        self.links = []
        self.feed(text)

    def handle_starttag(self, tag, attrs):
        attrs = dict(attrs)
        if attrs.get("id"):
            self.ids.add(attrs["id"])
        if tag == "a" and attrs.get("name"):
            self.ids.add(attrs["name"])
        for name in ("href", "src"):
            if attrs.get(name):
                self.links.append((self.getpos()[0], attrs[name]))


def check_html(root):
    """Check rendered local links without fetching external services."""
    root = root.resolve()
    pages = {path: HTMLLinks(path.read_text(encoding="utf-8"))
             for path in root.rglob("*.html") if path.is_file()}
    if not pages:
        return [f"No HTML pages found in {root}"]
    errors = []
    for path, page in pages.items():
        for line, destination in page.links:
            link = urlsplit(destination)
            if link.scheme or link.netloc:
                continue
            if not link.path:
                target = path
            elif link.path.startswith("/"):
                target = (root / unquote(link.path).lstrip("/")).resolve()
            else:
                target = (path.parent / unquote(link.path)).resolve()
            if target.is_dir():
                target /= "index.html"
            location = f"{path.relative_to(root).as_posix()}:{line}"
            if not target.is_relative_to(root):
                errors.append(f"{location}: HTML target outside site {destination}")
            elif not target.is_file():
                errors.append(f"{location}: missing HTML target {destination}")
            elif link.fragment and target in pages and unquote(link.fragment) not in pages[target].ids:
                errors.append(f"{location}: missing HTML anchor {destination}")
    return errors


def main_html(root):
    errors = check_html(root)
    if errors:
        print("\n".join(errors), file=sys.stderr)
        return 1
    print("HTML checks passed: local navigation, assets and anchors.")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=ROOT,
                        help="documentation root, including an unpacked Cargo package")
    parser.add_argument("--html", type=Path,
                        help="check a built HTML site instead of Markdown sources")
    args = parser.parse_args()
    raise SystemExit(main_html(args.html) if args.html is not None else main(args.root))
