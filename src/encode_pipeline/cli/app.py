"""Canonical public entry point for HelixWeave local commands."""

from __future__ import annotations

from collections.abc import Sequence
import sys

from encode_pipeline.cli import admin
from encode_pipeline.cli import local_platform
from encode_pipeline.deployment import cli as deployment_cli


def main(argv: Sequence[str] | None = None) -> int:
    """Dispatch supported commands without changing compatibility arguments."""
    arguments = list(sys.argv[1:] if argv is None else argv)
    if arguments[:1] == ["admin"]:
        return admin.main(arguments[1:])
    if arguments[:1] and arguments[0] in deployment_cli.COMMANDS:
        return deployment_cli.main(arguments)
    return local_platform.main(arguments)


__all__ = ["main"]
