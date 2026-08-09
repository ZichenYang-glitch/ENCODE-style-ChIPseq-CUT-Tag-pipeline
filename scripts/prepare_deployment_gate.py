#!/usr/bin/env python3
"""Prepare fixed-policy deployment Gate evidence locally; never dispatch it."""

from __future__ import annotations

import argparse
from collections.abc import Callable, Sequence
import json
import os
from pathlib import Path
import stat
import sys

from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.gate import (
    GateObserver,
    GatePolicy,
    GateRequest,
    UnavailableGateObserver,
    prepare_gate_request,
)
from encode_pipeline.deployment.models import BULK_RNASEQ_RUNTIME, ENCODE_RUNTIME


class _Parser(argparse.ArgumentParser):
    def error(self, message: str) -> None:
        del message
        raise fail("GATE_PREPARATION_ARGUMENT_INVALID", "Gate arguments are invalid.")


def _parser() -> argparse.ArgumentParser:
    parser = _Parser(description=__doc__)
    parser.add_argument("--head-sha", required=True)
    parser.add_argument("--release-identity", required=True)
    parser.add_argument("--encode-runtime-identity", required=True)
    parser.add_argument("--bulk-rnaseq-runtime-identity", required=True)
    parser.add_argument("--task-identity", required=True)
    return parser


def prepare(
    arguments: argparse.Namespace,
    *,
    observer: GateObserver,
    policy: GatePolicy,
    fault: Callable[[str], None] | None = None,
) -> GateRequest:
    plan, cleanup_script, request = prepare_gate_request(
        policy=policy,
        observer=observer,
        task_identity=arguments.task_identity,
        head_sha=arguments.head_sha,
        release_identity=arguments.release_identity,
        runtime_identities={
            ENCODE_RUNTIME: arguments.encode_runtime_identity,
            BULK_RNASEQ_RUNTIME: arguments.bulk_rnaseq_runtime_identity,
        },
    )
    output = policy.output_directory(arguments.task_identity)
    directory_descriptor = _open_output_directory(
        output, policy, arguments.task_identity
    )
    try:
        _write_exclusive(
            directory_descriptor,
            "cleanup-plan.json",
            canonical_json_bytes(plan.to_dict()),
            mode=0o600,
            expected_uid=policy.operator_uid,
            expected_gid=policy.operator_gid,
        )
        _write_exclusive(
            directory_descriptor,
            "cleanup.sh",
            cleanup_script,
            mode=0o500,
            expected_uid=policy.operator_uid,
            expected_gid=policy.operator_gid,
        )
        os.fsync(directory_descriptor)
        _inject(fault, "cleanup-evidence-synced")
        _write_exclusive(
            directory_descriptor,
            "gate-request.json",
            canonical_json_bytes(request.to_dict()),
            mode=0o600,
            expected_uid=policy.operator_uid,
            expected_gid=policy.operator_gid,
        )
        os.fsync(directory_descriptor)
        _inject(fault, "request-synced")
    except DeploymentError:
        raise
    except OSError:
        raise fail("GATE_OUTPUT_INVALID", "Gate output verification failed.") from None
    finally:
        os.close(directory_descriptor)
    return request


def _open_output_directory(path: Path, policy: GatePolicy, task_identity: str) -> int:
    expected = policy.output_directory(task_identity)
    if path != expected or not path.is_absolute():
        raise fail("GATE_OUTPUT_INVALID", "Gate output storage is invalid.")
    try:
        observed = path.lstat()
        flags = (
            os.O_RDONLY
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_DIRECTORY", 0)
            | getattr(os, "O_NOFOLLOW", 0)
        )
        descriptor = os.open(path, flags)
        anchored = os.fstat(descriptor)
    except OSError:
        raise fail("GATE_OUTPUT_INVALID", "Gate output storage is invalid.") from None
    if (
        not stat.S_ISDIR(observed.st_mode)
        or stat.S_ISLNK(observed.st_mode)
        or observed.st_uid != policy.operator_uid
        or observed.st_gid != policy.operator_gid
        or stat.S_IMODE(observed.st_mode) != 0o700
        or (observed.st_dev, observed.st_ino) != (anchored.st_dev, anchored.st_ino)
    ):
        os.close(descriptor)
        raise fail("GATE_OUTPUT_INVALID", "Gate output storage is invalid.")
    return descriptor


def _write_exclusive(
    directory_descriptor: int,
    name: str,
    content: bytes,
    *,
    mode: int,
    expected_uid: int,
    expected_gid: int,
) -> None:
    if name not in {"cleanup-plan.json", "cleanup.sh", "gate-request.json"}:
        raise fail("GATE_OUTPUT_INVALID", "Gate output verification failed.")
    flags = (
        os.O_WRONLY
        | os.O_CREAT
        | os.O_EXCL
        | getattr(os, "O_CLOEXEC", 0)
        | getattr(os, "O_NOFOLLOW", 0)
    )
    try:
        descriptor = os.open(name, flags, mode, dir_fd=directory_descriptor)
    except OSError:
        raise fail(
            "GATE_OUTPUT_EXISTS",
            "Gate output already exists or is unavailable.",
            recoverable=True,
        ) from None
    try:
        os.fchmod(descriptor, mode)
        offset = 0
        while offset < len(content):
            written = os.write(descriptor, content[offset:])
            if written <= 0:
                raise OSError
            offset += written
        os.fsync(descriptor)
        observed = os.fstat(descriptor)
        path_after = os.stat(name, dir_fd=directory_descriptor, follow_symlinks=False)
        if (
            observed.st_nlink != 1
            or not stat.S_ISREG(observed.st_mode)
            or observed.st_uid != expected_uid
            or observed.st_gid != expected_gid
            or stat.S_IMODE(observed.st_mode) != mode
            or observed.st_size != len(content)
            or (observed.st_dev, observed.st_ino)
            != (path_after.st_dev, path_after.st_ino)
        ):
            raise OSError
    except OSError:
        raise fail("GATE_OUTPUT_INVALID", "Gate output verification failed.") from None
    finally:
        os.close(descriptor)


def _inject(fault: Callable[[str], None] | None, point: str) -> None:
    if fault is not None:
        fault(point)


def main(
    argv: Sequence[str] | None = None,
    *,
    observer: GateObserver | None = None,
    policy: GatePolicy | None = None,
    fault: Callable[[str], None] | None = None,
) -> int:
    try:
        arguments = _parser().parse_args(argv)
        request = prepare(
            arguments,
            observer=UnavailableGateObserver() if observer is None else observer,
            policy=GatePolicy.supported() if policy is None else policy,
            fault=fault,
        )
    except DeploymentError as error:
        print(
            json.dumps(
                {
                    "schema_version": "helixweave-deployment-gate-preparation-v1",
                    "status": "error",
                    "issue": error.issue.to_dict(),
                },
                sort_keys=True,
                separators=(",", ":"),
            ),
            file=sys.stderr,
        )
        return 64 if error.issue.code.endswith("ARGUMENT_INVALID") else 70
    except Exception:
        print(
            '{"issue":{"code":"GATE_PREPARATION_FAILED",'
            '"message":"Gate preparation failed.","recoverable":false},'
            '"schema_version":"helixweave-deployment-gate-preparation-v1",'
            '"status":"error"}',
            file=sys.stderr,
        )
        return 70
    print(
        json.dumps(
            {
                "schema_version": "helixweave-deployment-gate-preparation-v1",
                "status": "prepared",
                "request_identity": request.identity,
                "task_identity": request.task_identity,
                "head_sha": request.head_sha,
                "release_identity": request.release_identity,
            },
            sort_keys=True,
            separators=(",", ":"),
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
