"""Public source-checkout commands keep the clean provenance boundary."""

from __future__ import annotations

import json
from pathlib import Path
import re

import pytest


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
CANONICAL_BOOTSTRAP = "python3 -I -S scripts/checkout_bootstrap.py --repository-root ."
SOURCE_CHECKOUT_DOCS = {
    "AGENTS.md": ("pytest", "validate"),
    "README.md": ("pytest", "validate"),
    "docs/development/harness.md": ("pytest",),
}
UNPROTECTED_COMMANDS = (
    re.compile(r"(?m)^\s*python3?\s+-m\s+pytest\b"),
    re.compile(r"(?m)^\s*python3?\s+scripts/validate_samples\.py\b"),
    re.compile(r"(?m)^\s*python3?\s+scripts/export_openapi\.py\b"),
    re.compile(r"(?m)^\s*python3?\s+scripts/run_local_platform\.py\b"),
)


def _one_line_commands(content: str) -> str:
    joined = re.sub(r"\\\s*\n\s*", " ", content)
    return re.sub(r"[ \t]+", " ", joined)


@pytest.mark.parametrize(
    ("relative_path", "commands"),
    SOURCE_CHECKOUT_DOCS.items(),
)
def test_public_source_checkout_commands_use_canonical_bootstrap(
    relative_path: str,
    commands: tuple[str, ...],
) -> None:
    content = (REPOSITORY_ROOT / relative_path).read_text(encoding="utf-8")
    for pattern in UNPROTECTED_COMMANDS:
        assert pattern.search(content) is None, relative_path

    one_line = _one_line_commands(content)
    for command in commands:
        assert f"{CANONICAL_BOOTSTRAP} {command}" in one_line


def test_documented_openapi_regeneration_delegates_to_canonical_bootstrap() -> None:
    agents = (REPOSITORY_ROOT / "AGENTS.md").read_text(encoding="utf-8")
    package = json.loads(
        (REPOSITORY_ROOT / "frontend/package.json").read_text(encoding="utf-8")
    )
    scripts = package["scripts"]

    assert "npm --prefix frontend run openapi:regenerate" in agents
    assert scripts["openapi:export"] == (
        "python3 -I -S ../scripts/checkout_bootstrap.py "
        "--repository-root .. openapi --output openapi.json"
    )
    assert scripts["openapi:regenerate"] == (
        "npm run openapi:export && npm run openapi:generate"
    )
