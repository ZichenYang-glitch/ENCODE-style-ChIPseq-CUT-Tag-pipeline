from __future__ import annotations

from io import StringIO
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

import encode_pipeline.deployment.bundle_cli as bundle_cli
from encode_pipeline.deployment.bundle_cli import (
    BUILD,
    BUNDLE_CLI_RESULT_SCHEMA,
    WHEEL_LOCK_CREATE,
    WHEEL_LOCK_VERIFY,
    main,
    parse_command,
)
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    ENCODE_RUNTIME,
    PLATFORM,
)


def _absolute(tmp_path: Path, name: str) -> str:
    return str(tmp_path / name)


def test_bundle_help_discovers_build_and_wheel_lock_commands() -> None:
    output = StringIO()

    assert main(["--help"], stdout=output) == 0

    rendered = output.getvalue()
    assert "wheel-lock create" in rendered
    assert "wheel-lock verify" in rendered
    assert "build" in rendered


def test_bundle_parser_accepts_only_the_fixed_absolute_grammars(tmp_path: Path) -> None:
    created = parse_command(
        [
            "wheel-lock",
            "create",
            "--wheelhouse",
            _absolute(tmp_path, "wheelhouse"),
            "--output",
            _absolute(tmp_path, "wheel-lock.json"),
        ]
    )
    verified = parse_command(
        [
            "wheel-lock",
            "verify",
            "--wheelhouse",
            _absolute(tmp_path, "wheelhouse"),
            "--lock",
            _absolute(tmp_path, "wheel-lock.json"),
        ]
    )
    platform = parse_command(
        [
            "build",
            "--component",
            PLATFORM,
            "--wheel",
            _absolute(tmp_path, "helixweave.whl"),
            "--wheelhouse",
            _absolute(tmp_path, "wheelhouse"),
            "--wheel-lock",
            _absolute(tmp_path, "wheel-lock.json"),
            "--output",
            _absolute(tmp_path, "platform.tar"),
            "--scratch-root",
            _absolute(tmp_path, "scratch"),
        ]
    )

    assert created.command == WHEEL_LOCK_CREATE
    assert verified.command == WHEEL_LOCK_VERIFY
    assert platform.command == BUILD
    assert platform.component == PLATFORM

    with pytest.raises(DeploymentError, match="DEPLOYMENT_COMMAND_INVALID"):
        parse_command(
            [
                "wheel-lock",
                "create",
                "--wheelhouse",
                "relative",
                "--output",
                _absolute(tmp_path, "wheel-lock.json"),
            ]
        )


def test_wheel_lock_create_and_verify_emit_canonical_receipts(
    tmp_path: Path,
    monkeypatch,
) -> None:
    output_path = tmp_path / "wheel-lock.json"
    lock = SimpleNamespace(
        identity="sha256-" + "a" * 64,
        wheels=(object(), object()),
        to_bytes=lambda: b'{"canonical":true}\n',
    )
    monkeypatch.setattr(bundle_cli, "build_platform_wheel_lock", lambda _root: lock)
    monkeypatch.setattr(
        bundle_cli,
        "verify_platform_wheel_lock",
        lambda _root, _path: lock,
    )
    created = StringIO()
    verified = StringIO()

    assert (
        main(
            [
                "wheel-lock",
                "create",
                "--wheelhouse",
                str(tmp_path / "wheelhouse"),
                "--output",
                str(output_path),
            ],
            stdout=created,
        )
        == 0
    )
    assert output_path.read_bytes() == lock.to_bytes()
    assert (
        main(
            [
                "wheel-lock",
                "verify",
                "--wheelhouse",
                str(tmp_path / "wheelhouse"),
                "--lock",
                str(output_path),
            ],
            stdout=verified,
        )
        == 0
    )

    for stream, command in (
        (created, WHEEL_LOCK_CREATE),
        (verified, WHEEL_LOCK_VERIFY),
    ):
        receipt = json.loads(stream.getvalue())
        assert receipt == {
            "schema_version": BUNDLE_CLI_RESULT_SCHEMA,
            "command": command,
            "status": "ok",
            "result": {
                "component": PLATFORM,
                "lock_path": str(output_path),
                "lock_identity": lock.identity,
                "wheel_count": 2,
            },
        }


@pytest.mark.parametrize(
    ("component", "flags", "builder_name"),
    [
        (
            PLATFORM,
            (
                "--wheel",
                "candidate.whl",
                "--wheelhouse",
                "wheelhouse",
                "--wheel-lock",
                "wheel-lock.json",
                "--output",
                "platform.tar",
                "--scratch-root",
                "scratch",
            ),
            "build_platform_bundle",
        ),
        (
            ENCODE_RUNTIME,
            (
                "--sdist-root",
                "sdist",
                "--micromamba",
                "micromamba",
                "--archive-cache",
                "archive-cache",
                "--output",
                "encode.tar",
            ),
            "build_encode_runtime_bundle",
        ),
        (
            BULK_RNASEQ_RUNTIME,
            (
                "--runtime-root",
                "bulk-runtime",
                "--output",
                "bulk.tar",
            ),
            "build_bulk_rnaseq_runtime_bundle",
        ),
    ],
)
def test_bundle_build_delegates_to_the_existing_component_builder(
    tmp_path: Path,
    monkeypatch,
    component: str,
    flags: tuple[str, ...],
    builder_name: str,
) -> None:
    arguments = ["build", "--component", component]
    named_paths: dict[str, Path] = {}
    for name, relative in zip(flags[::2], flags[1::2], strict=True):
        path = tmp_path / relative
        arguments.extend((name, str(path)))
        named_paths[name] = path
    output_path = named_paths["--output"]
    observed: list[tuple[object, ...]] = []

    def build(*values, **keywords):
        observed.append((*values, keywords))
        output_path.write_bytes(b"canonical bundle")
        return SimpleNamespace(component=component, identity="sha256-" + "b" * 64)

    monkeypatch.setattr(bundle_cli, builder_name, build)
    output = StringIO()

    assert main(arguments, stdout=output) == 0

    receipt = json.loads(output.getvalue())
    assert receipt["schema_version"] == BUNDLE_CLI_RESULT_SCHEMA
    assert receipt["command"] == BUILD
    assert receipt["status"] == "ok"
    assert receipt["result"] == {
        "component": component,
        "bundle_path": str(output_path),
        "bundle_identity": (
            "sha256-6f7567e64e88b3abe52d178d67d2e0b6382d46c276b9547e458d43b23397271a"
        ),
        "manifest_identity": "sha256-" + "b" * 64,
    }
    assert len(observed) == 1


def test_bundle_failure_is_json_and_does_not_echo_input_paths(
    tmp_path: Path,
    monkeypatch,
) -> None:
    private = tmp_path / "private" / "candidate.whl"

    def reject(*_values, **_keywords):
        raise fail(
            "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
            f"invalid private source {private}",
            component=PLATFORM,
        )

    monkeypatch.setattr(bundle_cli, "build_platform_bundle", reject)
    errors = StringIO()
    arguments = [
        "build",
        "--component",
        PLATFORM,
        "--wheel",
        str(private),
        "--wheelhouse",
        str(tmp_path / "wheelhouse"),
        "--wheel-lock",
        str(tmp_path / "wheel-lock.json"),
        "--output",
        str(tmp_path / "platform.tar"),
        "--scratch-root",
        str(tmp_path / "scratch"),
    ]

    assert main(arguments, stderr=errors) == 65

    rendered = errors.getvalue()
    assert str(tmp_path) not in rendered
    receipt = json.loads(rendered)
    assert receipt["status"] == "error"
    assert receipt["issue"] == {
        "code": "DEPLOYMENT_BUNDLE_SOURCE_INVALID",
        "message": "Deployment verification failed.",
        "recoverable": False,
        "component": PLATFORM,
    }
