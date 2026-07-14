"""Repository documentation link contracts."""

from __future__ import annotations

import re
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INLINE_LINK_RE = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
REFERENCE_LINK_RE = re.compile(r"^\s*\[[^\]]+\]:\s*(\S+)", re.MULTILINE)
EXTERNAL_SCHEMES = {"http", "https", "mailto"}


def _markdown_files() -> list[Path]:
    roots = [
        *PROJECT_ROOT.glob("*.md"),
        *PROJECT_ROOT.glob("containers/**/*.md"),
        *PROJECT_ROOT.glob("docs/**/*.md"),
        *PROJECT_ROOT.glob("test/**/README.md"),
    ]
    return sorted({path for path in roots if path.is_file()})


def _destinations(markdown: str) -> list[str]:
    return [
        *(match.group(1) for match in INLINE_LINK_RE.finditer(markdown)),
        *(match.group(1) for match in REFERENCE_LINK_RE.finditer(markdown)),
    ]


def _path_part(destination: str) -> str | None:
    destination = destination.strip()
    if destination.startswith("<") and ">" in destination:
        destination = destination[1 : destination.index(">")]
    else:
        destination = destination.split(maxsplit=1)[0]

    if not destination or destination.startswith("#"):
        return None

    parsed = urlsplit(destination)
    if parsed.scheme.lower() in EXTERNAL_SCHEMES or parsed.netloc:
        return None

    return unquote(parsed.path) or None


@pytest.mark.parametrize(
    "document", _markdown_files(), ids=lambda path: str(path.relative_to(PROJECT_ROOT))
)
def test_relative_document_links_resolve(document: Path) -> None:
    """Every relative Markdown link must resolve to a repository path."""
    markdown = document.read_text(encoding="utf-8")
    missing: list[str] = []

    for destination in _destinations(markdown):
        path_part = _path_part(destination)
        if path_part is None:
            continue
        target = (document.parent / path_part).resolve()
        if not target.exists():
            missing.append(destination)

    assert not missing, f"broken relative link(s) in {document}: {missing}"
