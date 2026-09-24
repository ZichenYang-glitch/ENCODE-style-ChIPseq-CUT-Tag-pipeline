"""Private, bounded cancellation diagnostics for the protected acceptance harness.

Two capture kinds live here, both in the owner-only directory resolved by
:func:`failure_diagnostics.resolve_private_diagnostics_root`:

* the raw stdout/stderr of a harness-owned worker session, kept instead of
  being discarded; and
* timestamped pre-destruction snapshots of Redis/RQ, SQLite, the worker process
  tree, and the managed container scope, taken while those records still exist.

Nothing in this module is published. Raw text stays out of the public API, the
downloadable evidence bundle, and notification mail; the public evidence tree
continues to carry only fixed codes, sizes, and digests. Every entry point is
best-effort and returns ``None`` on any problem, because a diagnostic must never
replace the execution result it was meant to explain.
"""

from __future__ import annotations

import json
import os
from pathlib import Path
import re
import tempfile
import time
from datetime import datetime, timezone

from .failure_diagnostics import resolve_private_diagnostics_root


PRIVATE_DIAGNOSTICS_ENV = "HELIXWEAVE_BULK_RNASEQ_PRIVATE_DIAGNOSTICS_DIR"
SNAPSHOT_SCHEMA_VERSION = "1.0.0"

_SNAPSHOT_LABEL = re.compile(r"[a-z][a-z0-9-]{0,63}")
_SNAPSHOT_BYTE_LIMIT = 4 * 1024 * 1024
_VALUE_LIMIT = 64 * 1024

__all__ = [
    "PRIVATE_DIAGNOSTICS_ENV",
    "SNAPSHOT_SCHEMA_VERSION",
    "WorkerStreamCapture",
    "capture_timestamp",
    "resolve_private_diagnostics_root",
    "write_cancellation_snapshot",
]


def capture_timestamp() -> str:
    """Return one bounded wall-clock token for a private capture name."""
    return datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%S.%fZ")


def _open_private_file(path: Path) -> object:
    descriptor = os.open(
        path,
        os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0),
        0o600,
    )
    # The worker child needs this descriptor across exec; dup2 clears the flag
    # only for the standard streams, so mark it inheritable explicitly.
    os.set_inheritable(descriptor, True)
    os.fchmod(descriptor, 0o600)
    return os.fdopen(descriptor, "w", encoding="utf-8", errors="replace")


def _write_private_text(path: Path, content: str) -> None:
    with _open_private_file(path) as stream:
        stream.write(content)
        stream.flush()
        os.fsync(stream.fileno())


def _write_private_json(path: Path, payload: object) -> None:
    rendered = json.dumps(
        payload,
        ensure_ascii=True,
        sort_keys=True,
        separators=(",", ":"),
        default=str,
    )
    if len(rendered.encode("utf-8")) > _SNAPSHOT_BYTE_LIMIT:
        rendered = json.dumps(
            {
                "snapshot_truncated": True,
                "rendered_bytes": len(rendered.encode("utf-8")),
                "prefix_text": rendered[:_SNAPSHOT_BYTE_LIMIT],
            },
            ensure_ascii=True,
            sort_keys=True,
            separators=(",", ":"),
        )
    _write_private_text(path, rendered + "\n")


class WorkerStreamCapture:
    """Owner-only stdout/stderr files for one harness-owned worker session."""

    def __init__(self, *, root: Path, ordinal: int) -> None:
        self.directory = Path(
            tempfile.mkdtemp(
                prefix=f"worker-{ordinal:02d}-{capture_timestamp()}-", dir=root
            )
        )
        self.stdout_path = self.directory / "worker.stdout"
        self.stderr_path = self.directory / "worker.stderr"
        self._stdout = _open_private_file(self.stdout_path)
        try:
            self._stderr = _open_private_file(self.stderr_path)
        except Exception:
            self._stdout.close()
            raise
        self._closed = False

    @property
    def stdout_handle(self):
        return self._stdout

    @property
    def stderr_handle(self):
        return self._stderr

    def close(self, *, returncode: int | None) -> None:
        """Flush and seal both streams; never raise."""
        if self._closed:
            return
        self._closed = True
        for handle in (self._stdout, self._stderr):
            try:
                handle.flush()
            except Exception:
                pass
            try:
                handle.close()
            except Exception:
                pass
        try:
            _write_private_json(
                self.directory / "worker-exit.json",
                {"returncode": returncode, "closed_at": capture_timestamp()},
            )
        except Exception:
            return


def write_cancellation_snapshot(
    *,
    root: Path,
    label: str,
    document: dict[str, object],
) -> str | None:
    """Write one timestamped private snapshot; never raise.

    Returns the capture directory name, or ``None`` when nothing was written.
    """
    if not isinstance(label, str) or _SNAPSHOT_LABEL.fullmatch(label) is None:
        return None
    try:
        destination = Path(
            tempfile.mkdtemp(prefix=f"{label}-{capture_timestamp()}-", dir=root)
        )
        payload = dict(document)
        payload["schema_version"] = SNAPSHOT_SCHEMA_VERSION
        payload["label"] = label
        payload["private_record"] = destination.name
        payload["captured_at_wall"] = capture_timestamp()
        payload["captured_at_monotonic"] = time.monotonic()
        _write_private_json(destination / "snapshot.json", payload)
        return destination.name
    except Exception:
        return None
