"""Fixed Peri runtime identity. Paths are supplied by the operator."""

from __future__ import annotations

from hashlib import sha256
from pathlib import Path


PERI_COMMIT = "c7d154c9985cbe08ba25e964cb7468d15014e6d0"
PATCH_SHA256 = "06272e1a66423a9e562509a38d044cd8506a1c10e61ef1fd3fe239107680587b"
BINARY_SHA256 = "60144bbe0919992b8913b9031696b2e30b1e9f82d95048ae3b07c887d5d85517"
SERVER_NAME = "helixweave_a1"
SCHEMA_TOOL = "get_authoring_schema"
DRAFT_TOOL = "get_draft_summary"
CANONICAL_SCHEMA_TOOL = f"mcp__{SERVER_NAME}__{SCHEMA_TOOL}"
CANONICAL_DRAFT_TOOL = f"mcp__{SERVER_NAME}__{DRAFT_TOOL}"
WIRE_TOOLS = frozenset({SCHEMA_TOOL, DRAFT_TOOL})
CANONICAL_TOOLS = frozenset({CANONICAL_SCHEMA_TOOL, CANONICAL_DRAFT_TOOL})
META_TOOLS = frozenset({"SearchExtraTools", "ExecuteExtraTool"})
MEMORY_MAX_BYTES = 512 * 1024 * 1024
PID_MAX = 128


def file_sha256(path: Path) -> str:
    digest = sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def binary_matches(path: Path) -> bool:
    try:
        resolved = path.resolve(strict=True)
    except OSError:
        return False
    return resolved.is_file() and file_sha256(resolved) == BINARY_SHA256
