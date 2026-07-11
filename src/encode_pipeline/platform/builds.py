"""Durable workflow build identity primitives."""

from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import PurePosixPath


@dataclass(frozen=True)
class WorkflowBuildIdentity:
    """Content-addressed identity for the workflow source used by one run."""

    workflow_id: str
    adapter_version: str
    scheme: str
    logical_entrypoint: str
    digest: str
    captured_at: datetime

    def __post_init__(self) -> None:
        for field_name in ("workflow_id", "adapter_version", "scheme"):
            value = getattr(self, field_name)
            if not isinstance(value, str) or not value.strip():
                raise ValueError(f"{field_name} must be a non-empty string")
            object.__setattr__(self, field_name, value.strip())

        entrypoint = self.logical_entrypoint
        if not isinstance(entrypoint, str) or not entrypoint.strip():
            raise ValueError("logical_entrypoint must be a non-empty string")
        entrypoint = entrypoint.strip()
        path = PurePosixPath(entrypoint)
        if (
            path.is_absolute()
            or entrypoint.startswith("~")
            or "\\" in entrypoint
            or any(part in {"", ".", ".."} for part in path.parts)
        ):
            raise ValueError("logical_entrypoint must be a safe relative POSIX path")
        object.__setattr__(self, "logical_entrypoint", entrypoint)

        digest = self.digest
        if (
            not isinstance(digest, str)
            or len(digest) != 64
            or digest.lower() != digest
            or any(character not in "0123456789abcdef" for character in digest)
        ):
            raise ValueError("digest must be a lowercase SHA-256 hex digest")
        if not isinstance(self.captured_at, datetime):
            raise ValueError("captured_at must be a datetime")

    def to_dict(self) -> dict[str, object]:
        """Return a serialization-ready representation."""
        return {
            "workflow_id": self.workflow_id,
            "adapter_version": self.adapter_version,
            "scheme": self.scheme,
            "logical_entrypoint": self.logical_entrypoint,
            "digest": self.digest,
            "captured_at": self.captured_at,
        }

    def matches(self, other: object) -> bool:
        """Return whether *other* identifies the same immutable build."""
        return isinstance(other, WorkflowBuildIdentity) and (
            self.workflow_id,
            self.adapter_version,
            self.scheme,
            self.logical_entrypoint,
            self.digest,
        ) == (
            other.workflow_id,
            other.adapter_version,
            other.scheme,
            other.logical_entrypoint,
            other.digest,
        )
