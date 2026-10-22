"""Build the maintained Markdown as a Sphinx/Furo documentation site."""

from pathlib import Path
import re
import sys
import tomllib
from urllib.parse import quote, unquote, urlsplit

from docutils import nodes
from pygments.lexers import RustLexer

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
sys.path.insert(0, str(ROOT / "tools"))
from check_docs import DEFINITION, LINK, links, prose

package = tomllib.loads((ROOT / "Cargo.toml").read_text(encoding="utf-8"))["package"]

project = "Cocycle"
author = "Cocycle contributors"
copyright = "2026, Cocycle contributors"
version = package["version"]
release = f"{version} + Unreleased"
language = "en"

extensions = ["myst_parser", "sphinx.ext.mathjax", "sphinx_copybutton"]
source_suffix = {".md": "markdown"}
root_doc = "index"
exclude_patterns = ["_build", "__pycache__"]
myst_enable_extensions = ["colon_fence", "dollarmath", "fieldlist"]
myst_dmath_double_inline = True
myst_heading_anchors = 6

html_theme = "furo"
html_title = "Cocycle documentation"
html_baseurl = "https://aequiludium.github.io/cocycle-rs/"
html_theme_options = {
    "source_repository": "https://github.com/Aequiludium/cocycle-rs/",
    "source_branch": "main",
    "source_directory": "docs/",
    "light_css_variables": {
        "color-brand-primary": "#087D70",
        "color-brand-content": "#087D70",
    },
    "dark_css_variables": {
        "color-brand-primary": "#62C9B4",
        "color-brand-content": "#62C9B4",
    },
}
html_show_sourcelink = False
copybutton_prompt_text = r">>> |\.\.\. |\$ "
copybutton_prompt_is_regexp = True


def repository_links(app, docname, source):
    """Keep GitHub source links usable when only docs/ is published."""
    path = DOCS / f"{docname}.md"
    masked = prose(source[0])
    starts = {offset for offset, destination, error in links(source[0]) if destination and not error}
    replacements = []
    for pattern, group in ((LINK, 1), (DEFINITION, 2)):
        for match in pattern.finditer(masked):
            if match.start() not in starts or match[0].startswith("!"):
                continue  # Sphinx copies the local artwork into the site.
            link = urlsplit(match[group].strip("<>"))
            if link.scheme or link.netloc or not link.path:
                continue
            target = (path.parent / unquote(link.path)).resolve()
            if (target.is_relative_to(DOCS) or not target.is_relative_to(ROOT)
                    or not target.exists()):
                continue
            kind = "tree" if target.is_dir() else "blob"
            url = (f"https://github.com/Aequiludium/cocycle-rs/{kind}/main/"
                   + quote(target.relative_to(ROOT).as_posix()))
            if link.fragment:
                url += "#" + link.fragment
            replacements.append((*match.span(group), url))
    for start, end, url in sorted(replacements, reverse=True):
        source[0] = source[0][:start] + url + source[0][end:]


def record_html_anchors(app, doctree):
    """Register existing GitHub HTML anchors for MyST document fragments."""
    slugs = app.env.metadata[app.env.docname].setdefault("myst_slugs", {})
    for node in doctree.findall(nodes.raw):
        if node.get("format") == "html":
            for anchor in re.findall(r'<a\s+(?:id|name)=["\']([^"\']+)', node.astext()):
                slugs[anchor] = (node.line, anchor, anchor.replace("-", " "))


def setup(app):
    app.connect("source-read", repository_links)
    app.connect("doctree-read", record_html_anchors)
    app.add_lexer("rust,ignore", RustLexer)
    return {"parallel_read_safe": True, "parallel_write_safe": True}
