"""Private metadata shared by the archive publisher, storage and reader."""

from __future__ import annotations

from dataclasses import dataclass
import re


CHUNK_BYTES = 256 * 1024
SHA256 = re.compile(r"[0-9a-f]{64}")


class ArchiveStoreError(Exception):
    """A path-redacted, stable failure code suitable for API translation."""

    def __init__(self, code: str):
        self.code = code
        super().__init__(code)


@dataclass(frozen=True)
class SealedMember:
    size_bytes: int
    sha256: str
    chunks: tuple[str, ...]
    identity: tuple[int, ...]

    def __post_init__(self) -> None:
        if (
            self.size_bytes <= 0
            or SHA256.fullmatch(self.sha256) is None
            or len(self.chunks) != (self.size_bytes + CHUNK_BYTES - 1) // CHUNK_BYTES
            or any(SHA256.fullmatch(value) is None for value in self.chunks)
            or len(self.identity) != 7
        ):
            raise ArchiveStoreError("archive_integrity")


@dataclass(frozen=True)
class ArchiveVersion:
    project_id: str
    experiment_id: str
    experiment_revision_id: str
    archive_version_id: str
    track_key: str
    version_number: int
    display_name: str
    assembly: str
    track_type: str
    source: dict[str, str]
    members: dict[str, SealedMember]

    @property
    def has_index(self) -> bool:
        return self.track_type == "longrange"
