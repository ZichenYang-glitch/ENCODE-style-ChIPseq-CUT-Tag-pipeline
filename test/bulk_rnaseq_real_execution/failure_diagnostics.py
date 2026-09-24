"""Bounded, private failure logs for the existing protected acceptance harness."""

from __future__ import annotations

from hashlib import sha256
import os
from pathlib import Path
import re
import stat
import tempfile

from .support import _write_canonical_evidence_document


PRIVATE_DIAGNOSTICS_ENV = "HELIXWEAVE_BULK_RNASEQ_PRIVATE_DIAGNOSTICS_DIR"
_LOG_LIMIT = 64 * 1024
_TASK_LIMIT = 4096


def _read_regular(path: Path, *, tail: bool = False) -> tuple[bytes, int]:
    descriptor = os.open(path, os.O_RDONLY | os.O_NOFOLLOW | os.O_NONBLOCK)
    with os.fdopen(descriptor, "rb") as stream:
        metadata = os.fstat(stream.fileno())
        if not stat.S_ISREG(metadata.st_mode):
            raise ValueError("diagnostic input is not a regular file")
        if tail:
            stream.seek(max(0, metadata.st_size - _LOG_LIMIT))
        return stream.read(_LOG_LIMIT), metadata.st_size


def resolve_private_diagnostics_root(evidence_root: Path) -> Path:
    """Return the owner-only directory that holds raw, unpublished diagnostics.

    Every raw capture shares this rule set: the directory must be canonical,
    outside every uploadable evidence glob, outside a configured runner's
    cleanup, and owner-only. A caller that cannot satisfy the rule fails closed
    instead of writing raw text somewhere it could be published.
    """
    configured = os.environ.get(PRIVATE_DIAGNOSTICS_ENV)
    private_root = (
        Path(configured) if configured else evidence_root.parent / "private-diagnostics"
    )
    if not private_root.is_absolute() or private_root.resolve() != private_root:
        raise ValueError("private diagnostics directory must be canonical")
    if "evidence" in private_root.parts:
        raise ValueError("raw diagnostics must stay outside evidence upload globs")
    runner_temp = os.environ.get("RUNNER_TEMP")
    if configured and runner_temp and private_root.is_relative_to(Path(runner_temp)):
        raise ValueError("configured diagnostics must survive runner cleanup")
    private_root.mkdir(mode=0o700, parents=True, exist_ok=True)
    metadata = private_root.lstat()
    if (
        not stat.S_ISDIR(metadata.st_mode)
        or metadata.st_uid != os.getuid()
        or stat.S_IMODE(metadata.st_mode) != 0o700
    ):
        raise ValueError("private diagnostics directory must be owner-only")
    return private_root


def preserve_execution_failure(
    *,
    workspace: Path,
    evidence_root: Path,
    stage: str,
    reason_code: str | None = None,
    rq_exception: str | None = None,
) -> dict[str, object]:
    """Keep the first failed task's stderr before cleanup; never mask a failure.

    Only fixed codes, sizes and digests enter uploadable evidence. Raw text is
    kept separately in an owner-only directory, optionally outside RUNNER_TEMP.
    This records an observed nonzero task, not necessarily the causal first
    failure when multiple tasks fail concurrently.
    """
    files: dict[str, object] = {}
    document: dict[str, object] = {
        "stage": stage if stage in {"platform", "rapid-quant"} else "unknown",
        "reason_code": (
            reason_code
            if isinstance(reason_code, str)
            and re.fullmatch(r"[A-Z][A-Z0-9_]{0,127}", reason_code)
            else None
        ),
        "capture_status": "complete",
        "files": files,
    }
    try:
        private_root = resolve_private_diagnostics_root(evidence_root)
        destination = Path(tempfile.mkdtemp(prefix="failure-", dir=private_root))
        document["private_record"] = destination.name
        document["persistent_destination_configured"] = bool(
            os.environ.get(PRIVATE_DIAGNOSTICS_ENV)
        )

        def save(label: str, content: bytes, source_size: int) -> None:
            with (destination / label).open("xb") as stream:
                os.fchmod(stream.fileno(), 0o600)
                stream.write(content)
            files[label] = {
                "sha256": sha256(content).hexdigest(),
                "size": len(content),
                "source_size": source_size,
                "truncated": source_size > len(content),
            }

        failed_tasks: list[tuple[int, Path, int]] = []
        for count, exit_path in enumerate(
            (workspace / "engine/work").glob("*/*/.exitcode")
        ):
            if count >= _TASK_LIMIT:
                document["task_scan_truncated"] = True
                break
            if exit_path.resolve() != exit_path:
                raise ValueError("task log path is not canonical")
            content, _ = _read_regular(exit_path)
            if re.fullmatch(rb"[0-9]{1,3}\s*", content):
                code = int(content)
                if code != 0:
                    failed_tasks.append((exit_path.stat().st_mtime_ns, exit_path, code))
        if failed_tasks:
            _, exit_path, code = min(failed_tasks)
            document["task_exit_code"] = code
            for name in (".exitcode", ".command.err", ".command.out"):
                try:
                    content, size = _read_regular(exit_path.with_name(name))
                except FileNotFoundError:
                    continue
                save(name.removeprefix("."), content, size)
                if name == ".command.err" and (
                    b"error while loading shared libraries:" in content
                    and b"Permission denied" in content
                ):
                    document["task_signal"] = "SHARED_LIBRARY_PERMISSION_DENIED"
        path = workspace / "logs/nextflow.log"
        if path.resolve() != path:
            raise ValueError("workflow log path is not canonical")
        try:
            content, size = _read_regular(path, tail=True)
        except FileNotFoundError:
            pass
        else:
            save("nextflow.log", content, size)
        if rq_exception:
            content = rq_exception.encode("utf-8")
            save("rq-exception.txt", content[:_LOG_LIMIT], len(content))
    except Exception as error:
        document["capture_status"] = "incomplete"
        document["capture_error_type"] = type(error).__name__
        document["capture_errno"] = getattr(error, "errno", None)
    try:
        _write_canonical_evidence_document(
            document,
            evidence_root / "execution-failure.json",
            failure_message="execution failure evidence could not be retained",
        )
    except Exception:
        # Diagnostics must not replace the original execution result/exception.
        document["capture_status"] = "unpublished"
    return document
