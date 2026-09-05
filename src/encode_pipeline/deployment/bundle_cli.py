"""Public offline producer for existing HelixWeave deployment bundles."""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
import hashlib
import json
import os
from pathlib import Path
import stat
import sys
from typing import TextIO

from encode_pipeline.deployment import cli as deployment_cli
from encode_pipeline.deployment.bundle_builder import (
    build_bulk_rnaseq_runtime_bundle,
    build_encode_runtime_bundle,
    build_platform_bundle,
)
from encode_pipeline.deployment.errors import DeploymentError, DeploymentIssue, fail
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    COMPONENTS,
    ENCODE_RUNTIME,
    PLATFORM,
)
from encode_pipeline.deployment.platform_runtime import (
    build_platform_wheel_lock,
    verify_platform_wheel_lock,
)


BUNDLE_CLI_RESULT_SCHEMA = "helixweave-deployment-bundle-cli-result-v1"
BUILD = "build"
WHEEL_LOCK_CREATE = "wheel-lock-create"
WHEEL_LOCK_VERIFY = "wheel-lock-verify"

_HELP = """usage: helixweave bundle <command> [arguments]

Build the existing HelixWeave deployment transports without sudo or downloads.

commands:
  wheel-lock create   Create a canonical platform wheel lock.
  wheel-lock verify   Verify a wheel lock against its complete wheelhouse.
  build               Build a platform, ENCODE, or Bulk RNA-seq bundle.

Run `helixweave bundle wheel-lock --help` or
`helixweave bundle build --help` for exact command forms.
"""

_WHEEL_LOCK_HELP = """usage:
  helixweave bundle wheel-lock create --wheelhouse ABS --output ABS
  helixweave bundle wheel-lock verify --wheelhouse ABS --lock ABS
"""

_BUILD_HELP = """usage:
  helixweave bundle build --component platform --wheel ABS --wheelhouse ABS --wheel-lock ABS --output ABS --scratch-root ABS
  helixweave bundle build --component encode-runtime --sdist-root ABS --micromamba ABS --archive-cache ABS --output ABS
  helixweave bundle build --component bulk-rnaseq-runtime --runtime-root ABS --output ABS
"""


@dataclass(frozen=True)
class BundleProducerRequest:
    command: str
    component: str = PLATFORM
    paths: Mapping[str, Path] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.command not in {BUILD, WHEEL_LOCK_CREATE, WHEEL_LOCK_VERIFY}:
            raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")
        if self.component not in COMPONENTS or not isinstance(self.paths, Mapping):
            raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")


def _absolute_path(value: str) -> Path:
    path = Path(value)
    if (
        not path.is_absolute()
        or ".." in path.parts
        or any(character in value for character in ("\x00", "\n", "\r"))
    ):
        raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")
    return path


def _path_pairs(values: tuple[str, ...], names: tuple[str, ...]) -> dict[str, Path]:
    if len(values) != 2 * len(names) or values[::2] != names:
        raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")
    return {
        name.removeprefix("--").replace("-", "_"): _absolute_path(value)
        for name, value in zip(names, values[1::2], strict=True)
    }


def parse_command(argv: Sequence[str]) -> BundleProducerRequest:
    """Parse the fixed producer grammar without echoing rejected values."""
    if not isinstance(argv, Sequence) or isinstance(argv, (str, bytes)):
        raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")
    values = tuple(argv)
    if any(
        not isinstance(item, str)
        or not item
        or len(item) > 4096
        or any(ord(character) < 32 or ord(character) == 127 for character in item)
        for item in values
    ):
        raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")
    if values[:2] == ("wheel-lock", "create"):
        paths = _path_pairs(
            values[2:],
            ("--wheelhouse", "--output"),
        )
        return BundleProducerRequest(WHEEL_LOCK_CREATE, paths=paths)
    if values[:2] == ("wheel-lock", "verify"):
        paths = _path_pairs(
            values[2:],
            ("--wheelhouse", "--lock"),
        )
        return BundleProducerRequest(WHEEL_LOCK_VERIFY, paths=paths)
    if len(values) >= 3 and values[:2] == ("build", "--component"):
        component = values[2]
        expected = {
            PLATFORM: (
                "--wheel",
                "--wheelhouse",
                "--wheel-lock",
                "--output",
                "--scratch-root",
            ),
            ENCODE_RUNTIME: (
                "--sdist-root",
                "--micromamba",
                "--archive-cache",
                "--output",
            ),
            BULK_RNASEQ_RUNTIME: (
                "--runtime-root",
                "--output",
            ),
        }.get(component)
        if expected is None:
            raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")
        return BundleProducerRequest(
            BUILD,
            component=component,
            paths=_path_pairs(values[3:], expected),
        )
    raise fail("DEPLOYMENT_COMMAND_INVALID", "Deployment command is invalid.")


def _write_lock(path: Path, content: bytes) -> None:
    descriptor = -1
    created = False
    try:
        descriptor = os.open(
            path,
            os.O_WRONLY
            | os.O_CREAT
            | os.O_EXCL
            | getattr(os, "O_CLOEXEC", 0)
            | getattr(os, "O_NOFOLLOW", 0),
            0o444,
        )
        created = True
        offset = 0
        while offset < len(content):
            written = os.write(descriptor, content[offset:])
            if written == 0:
                raise OSError
            offset += written
        os.fsync(descriptor)
    except OSError:
        if descriptor >= 0:
            os.close(descriptor)
            descriptor = -1
        if created:
            try:
                path.unlink()
            except OSError:
                pass
        raise fail(
            "DEPLOYMENT_BUNDLE_OUTPUT_INVALID",
            "Deployment bundle output is invalid.",
            component=PLATFORM,
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def _file_identity(path: Path) -> str:
    descriptor = -1
    flags = os.O_RDONLY | getattr(os, "O_CLOEXEC", 0) | getattr(os, "O_NOFOLLOW", 0)
    try:
        descriptor = os.open(path, flags)
        before = os.fstat(descriptor)
        if not stat.S_ISREG(before.st_mode) or before.st_nlink != 1:
            raise OSError
        digest = hashlib.sha256()
        while True:
            chunk = os.read(descriptor, 1024 * 1024)
            if not chunk:
                break
            digest.update(chunk)
        after = os.fstat(descriptor)
        if (before.st_dev, before.st_ino, before.st_size, before.st_mtime_ns) != (
            after.st_dev,
            after.st_ino,
            after.st_size,
            after.st_mtime_ns,
        ):
            raise OSError
        return f"sha256-{digest.hexdigest()}"
    except OSError:
        raise fail(
            "DEPLOYMENT_BUNDLE_OUTPUT_INVALID",
            "Deployment bundle output is invalid.",
        ) from None
    finally:
        if descriptor >= 0:
            os.close(descriptor)


def execute_command(request: BundleProducerRequest) -> dict[str, object]:
    """Execute one offline producer command through the existing builders."""
    paths = request.paths
    if request.command == WHEEL_LOCK_CREATE:
        lock = build_platform_wheel_lock(paths["wheelhouse"])
        _write_lock(paths["output"], lock.to_bytes())
        return {
            "component": PLATFORM,
            "lock_path": str(paths["output"]),
            "lock_identity": lock.identity,
            "wheel_count": len(lock.wheels),
        }
    if request.command == WHEEL_LOCK_VERIFY:
        lock = verify_platform_wheel_lock(paths["wheelhouse"], paths["lock"])
        return {
            "component": PLATFORM,
            "lock_path": str(paths["lock"]),
            "lock_identity": lock.identity,
            "wheel_count": len(lock.wheels),
        }

    if request.component == PLATFORM:
        manifest = build_platform_bundle(
            paths["wheel"],
            paths["wheelhouse"],
            paths["wheel_lock"],
            paths["output"],
            scratch_root=paths["scratch_root"],
        )
    elif request.component == ENCODE_RUNTIME:
        manifest = build_encode_runtime_bundle(
            paths["sdist_root"],
            paths["micromamba"],
            paths["archive_cache"],
            paths["output"],
        )
    else:
        manifest = build_bulk_rnaseq_runtime_bundle(
            paths["runtime_root"],
            paths["output"],
        )
    if manifest.component != request.component:
        raise fail(
            "DEPLOYMENT_BUNDLE_BUILD_FAILED",
            "Deployment bundle build failed.",
            component=request.component,
        )
    return {
        "component": request.component,
        "bundle_path": str(paths["output"]),
        "bundle_identity": _file_identity(paths["output"]),
        "manifest_identity": manifest.identity,
    }


def _write_json(stream: TextIO, value: object) -> None:
    stream.write(
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=False,
            separators=(",", ":"),
            sort_keys=True,
        )
    )
    stream.write("\n")
    stream.flush()


def _help_for(arguments: tuple[str, ...]) -> str | None:
    if arguments in {("--help",), ("-h",)}:
        return _HELP
    if arguments in {
        ("wheel-lock", "--help"),
        ("wheel-lock", "-h"),
        ("wheel-lock", "create", "--help"),
        ("wheel-lock", "verify", "--help"),
    }:
        return _WHEEL_LOCK_HELP
    if arguments in {("build", "--help"), ("build", "-h")}:
        return _BUILD_HELP
    return None


def main(
    argv: Sequence[str] | None = None,
    *,
    stdout: TextIO | None = None,
    stderr: TextIO | None = None,
) -> int:
    """Run one producer operation and emit one canonical JSON receipt."""
    arguments = tuple(sys.argv[1:] if argv is None else argv)
    output = sys.stdout if stdout is None else stdout
    errors = sys.stderr if stderr is None else stderr
    help_text = _help_for(arguments)
    if help_text is not None:
        output.write(help_text)
        output.flush()
        return 0
    request: BundleProducerRequest | None = None
    try:
        request = parse_command(arguments)
        result = execute_command(request)
        envelope = {
            "schema_version": BUNDLE_CLI_RESULT_SCHEMA,
            "command": request.command,
            "status": "ok",
            "result": result,
        }
    except DeploymentError as error:
        exit_code = deployment_cli._exit_code(error.issue)
        public_issue = deployment_cli._public_issue(error.issue, exit_code)
        _write_json(
            errors,
            {
                "schema_version": BUNDLE_CLI_RESULT_SCHEMA,
                "command": None if request is None else request.command,
                "status": "error",
                "issue": public_issue.to_dict(),
            },
        )
        return exit_code
    except Exception:
        _write_json(
            errors,
            {
                "schema_version": BUNDLE_CLI_RESULT_SCHEMA,
                "command": None if request is None else request.command,
                "status": "error",
                "issue": DeploymentIssue(
                    code="DEPLOYMENT_FAILED",
                    message="Deployment operation failed.",
                ).to_dict(),
            },
        )
        return deployment_cli.EXIT_OPERATION
    _write_json(output, envelope)
    return deployment_cli.EXIT_OK


__all__ = [
    "BUILD",
    "BUNDLE_CLI_RESULT_SCHEMA",
    "WHEEL_LOCK_CREATE",
    "WHEEL_LOCK_VERIFY",
    "BundleProducerRequest",
    "execute_command",
    "main",
    "parse_command",
]
