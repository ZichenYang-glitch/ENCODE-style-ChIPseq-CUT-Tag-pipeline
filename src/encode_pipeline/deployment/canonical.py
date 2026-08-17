"""Canonical JSON and content identities for deployment documents."""

from __future__ import annotations

import hashlib
import json
from typing import Any


def canonical_json_bytes(value: object) -> bytes:
    """Render one JSON value with a single deterministic byte representation."""
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("utf-8")


def canonical_identity(value: object, *, scheme: str) -> str:
    """Return a domain-separated filesystem-safe SHA-256 identity."""
    digest = hashlib.sha256()
    for item in (scheme.encode("ascii"), canonical_json_bytes(value)):
        digest.update(len(item).to_bytes(8, byteorder="big", signed=False))
        digest.update(item)
    return f"sha256-{digest.hexdigest()}"


def without_key(value: dict[str, Any], key: str) -> dict[str, Any]:
    return {name: item for name, item in value.items() if name != key}


__all__ = ["canonical_identity", "canonical_json_bytes", "without_key"]
