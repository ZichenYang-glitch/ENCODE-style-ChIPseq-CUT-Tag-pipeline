"""Workflow-neutral identity for platform-managed local containers."""

from __future__ import annotations

from hashlib import sha256
import os
from pathlib import Path


MANAGED_CONTAINER_SCOPE_LABEL = "org.helixweave.workspace-scope"


def managed_container_scope(workspace: Path) -> str:
    """Return the stable, path-private scope owned by one absolute workspace."""
    if not isinstance(workspace, Path) or not workspace.is_absolute():
        raise ValueError("workspace must be an absolute pathlib.Path")
    value = str(workspace)
    if value != str(Path(value)) or any(
        character in value for character in ("\x00", "\n", "\r")
    ):
        raise ValueError("workspace must be a canonical absolute path")
    payload = b"helixweave-managed-container-scope-v1\x00" + os.fsencode(value)
    return sha256(payload).hexdigest()


def managed_container_endpoint_identity(
    executable: Path,
    unix_socket: Path,
) -> str:
    """Bind one command contract to server-owned Docker endpoint paths."""
    values: list[bytes] = []
    for name, path in (("executable", executable), ("unix_socket", unix_socket)):
        if not isinstance(path, Path) or not path.is_absolute():
            raise ValueError(f"{name} must be an absolute pathlib.Path")
        rendered = str(path)
        if rendered != str(Path(rendered)) or any(
            character in rendered for character in ("\x00", "\n", "\r")
        ):
            raise ValueError(f"{name} must be a canonical absolute path")
        values.append(os.fsencode(rendered))
    payload = b"helixweave-managed-container-endpoint-v1\x00" + b"\x00".join(values)
    return sha256(payload).hexdigest()
