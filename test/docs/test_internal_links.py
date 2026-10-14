"""Repository documentation link contracts."""

from __future__ import annotations

import re
import subprocess
from pathlib import Path
from urllib.parse import unquote, urlsplit

import pytest


PROJECT_ROOT = Path(__file__).resolve().parents[2]
INLINE_LINK_RE = re.compile(r"!?\[[^\]]*\]\(([^)]+)\)")
REFERENCE_LINK_RE = re.compile(r"^\s*\[[^\]]+\]:\s*(\S+)", re.MULTILINE)
EXTERNAL_SCHEMES = {"http", "https", "mailto"}

FENCED_CODE_RE = re.compile(r"(?ms)^(`{3,}|~{3,})[^\n]*\n.*?^\1[ \t]*$")
INLINE_CODE_RE = re.compile(r"`[^`\n]+`")


def _strip_code(markdown: str) -> str:
    """Remove fenced code blocks and inline code spans before link matching.

    Text inside code spans is not a Markdown link; matching the raw text
    mistakes spans like `[x](y)` or `ns["k"](arg)` for real links.
    """
    without_fences = FENCED_CODE_RE.sub("", markdown)
    return INLINE_CODE_RE.sub("", without_fences)


def _markdown_files() -> list[Path]:
    """Return the git-tracked Markdown files that count as repository docs.

    "Repository documentation" means tracked files: untracked local notes
    are not repository content and must not be swept into the check. When
    git is unavailable (e.g. running from an unpacked artifact) fail
    loudly instead of silently falling back to a glob.
    """
    try:
        result = subprocess.run(
            ["git", "ls-files", "-z"],
            cwd=PROJECT_ROOT,
            capture_output=True,
            check=False,
        )
    except OSError as error:
        raise RuntimeError(
            "git is required to enumerate tracked repository documentation"
        ) from error
    if result.returncode != 0:
        raise RuntimeError(
            "git ls-files failed "
            f"({result.returncode}): {result.stderr.decode('utf-8', 'replace').strip()}"
        )

    paths = []
    for relative in result.stdout.decode("utf-8").split("\0"):
        if not relative.endswith(".md"):
            continue
        parts = Path(relative).parts
        if (
            len(parts) == 1
            or parts[0] in {"containers", "docs"}
            or (parts[0] == "test" and parts[-1] == "README.md")
        ):
            paths.append(PROJECT_ROOT / relative)
    return sorted(paths)


def _destinations(markdown: str) -> list[str]:
    stripped = _strip_code(markdown)
    return [
        *(match.group(1) for match in INLINE_LINK_RE.finditer(stripped)),
        *(match.group(1) for match in REFERENCE_LINK_RE.finditer(stripped)),
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


# --- Backticked repository path references --------------------------------
#
# Bare paths in single-backtick spans (e.g. `docs/output-contract.md`) are
# not Markdown links, so the check above never sees them. The rules below
# are deliberately conservative so prose and code fragments do not turn
# into false positives: a span is only a candidate when it cannot be
# anything but a repository path reference.

ROOT_PREFIXES = (
    "src/",
    "test/",
    "docs/",
    "workflow/",
    "scripts/",
    "config/",
    "containers/",
    "frontend/",
    "profiles/",
)
TOP_LEVEL_FILES = {
    "AGENTS.md",
    "README.md",
    "CHANGELOG.md",
    "KNOWN_ISSUES.md",
    "pyproject.toml",
}
RUNTIME_PREFIXES = ("results/", "logs/", "workspace/", "benchmarks/", ".snakemake/")
SHELL_METACHARS = set("*${}|><=,'\"\\")
# Reference suffixes after a path: :42, :42-48, :73–81 (en/em dash),
# :symbol_name, and chains like :771、:824 joined by 、 or comma.
LINE_SUFFIX_RE = re.compile(r"(?:[、,]?\s*:(?:\d+(?:[-–—]\d+)?|[A-Za-z_]\w*))+$")


def _is_path_candidate(token: str) -> bool:
    if not token or any(char.isspace() for char in token):
        return False
    if any(char in SHELL_METACHARS for char in token):
        return False
    if token.startswith(("~", "/")):
        return False
    if token.startswith(RUNTIME_PREFIXES):
        return False
    if token.startswith(ROOT_PREFIXES):
        return True
    return token in TOP_LEVEL_FILES


def _backtick_path_candidates(markdown: str) -> list[str]:
    without_fences = FENCED_CODE_RE.sub("", markdown)
    # Spans inside a Markdown link's text are display labels; the link's
    # destination is already checked by test_relative_document_links_resolve.
    without_links = INLINE_LINK_RE.sub("", without_fences)
    without_links = REFERENCE_LINK_RE.sub("", without_links)
    return [
        match.group(0)[1:-1]
        for match in INLINE_CODE_RE.finditer(without_links)
        if _is_path_candidate(match.group(0)[1:-1])
    ]


def _candidate_exists(document: Path, token: str) -> bool:
    """Resolve against the repository root first, then the document's dir."""
    reference = token.split("#", 1)[0]
    reference = LINE_SUFFIX_RE.sub("", reference)
    if not reference:
        return True
    for base in (PROJECT_ROOT, document.parent):
        if (base / reference).resolve().exists():
            return True
    return False


@pytest.mark.parametrize(
    "document", _markdown_files(), ids=lambda path: str(path.relative_to(PROJECT_ROOT))
)
def test_backticked_repository_paths_resolve(document: Path) -> None:
    """Every backticked repository path reference must resolve."""
    markdown = document.read_text(encoding="utf-8")
    missing = [
        token
        for token in _backtick_path_candidates(markdown)
        if not _candidate_exists(document, token)
    ]
    assert not missing, f"broken backticked path reference(s) in {document}: {missing}"
