"""Behavior coverage for the bounded command executor and runtime index parsing."""

from __future__ import annotations

import json
import subprocess
from types import SimpleNamespace

import pytest

import encode_pipeline.deployment.native_contracts as native_module
from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
    without_key,
)
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.native_contracts import (
    encode_runtime_index_bytes,
    parse_encode_runtime_index,
)
from encode_pipeline.deployment.operator import (
    SYSTEMCTL,
    BoundedCommandExecutor,
)

_MAX_OUTPUT = 64 * 1024


def _reject(callable_, *args, code: str, **kwargs) -> None:
    with pytest.raises(DeploymentError) as captured:
        callable_(*args, **kwargs)
    assert captured.value.issue.code == code


@pytest.mark.parametrize(
    "argv,timeout",
    [
        ((), 1.0),
        (("/bin/echo",), 1.0),
        ((str(SYSTEMCTL),), 0.0),
        ((str(SYSTEMCTL),), 16.0),
        ((str(SYSTEMCTL),), -1.0),
    ],
)
def test_bounded_executor_rejects_invalid_invocations(argv, timeout) -> None:
    executor = BoundedCommandExecutor()
    _reject(executor.run, argv, timeout=timeout, code="OPERATOR_COMMAND_INVALID")


def test_bounded_executor_captures_bounded_stdout() -> None:
    result = BoundedCommandExecutor().run((str(SYSTEMCTL), "--version"), timeout=15.0)
    assert result.returncode == 0
    assert b"systemd" in result.stdout


def test_bounded_executor_fails_closed_on_os_error(monkeypatch) -> None:
    def raising_run(*args, **kwargs):
        raise OSError("no such process")

    monkeypatch.setattr(subprocess, "run", raising_run)
    executor = BoundedCommandExecutor()
    _reject(
        executor.run,
        (str(SYSTEMCTL), "--version"),
        timeout=15.0,
        code="OPERATOR_COMMAND_FAILED",
    )


def test_bounded_executor_rejects_oversized_output(monkeypatch) -> None:
    def noisy_run(argv, *, stdin, stdout, stderr, env, close_fds, timeout, check):
        stdout.write(b"x" * (_MAX_OUTPUT + 1))
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(subprocess, "run", noisy_run)
    executor = BoundedCommandExecutor()
    _reject(
        executor.run,
        (str(SYSTEMCTL), "--version"),
        timeout=15.0,
        code="OPERATOR_COMMAND_FAILED",
    )


def _package(digest: str) -> dict[str, object]:
    filename = f"pkg-{digest[:4]}.conda"
    url = f"https://conda.anaconda.org/conda-forge/linux-64/{filename}#{'a' * 32}"
    return {
        "url": url,
        "platform": "linux-64",
        "filename": filename,
        "md5": "a" * 32,
        "archive_path": (
            f"{native_module.ENCODE_PACKAGE_ARCHIVE_ROOT}/{digest}/{filename}"
        ),
        "size_bytes": 7,
        "sha256": digest,
    }


def _index_document() -> dict[str, object]:
    package = _package("b" * 64)
    content = encode_runtime_index_bytes(
        workflow_build_identity=f"sha256-{'3' * 64}",
        micromamba={
            "path": native_module.ENCODE_MICROMAMBA_PATH,
            "size_bytes": 120,
            "sha256": "4" * 64,
        },
        packages=(package,),
        environments=(
            {
                "environment_path": "workflow/envs/runner.yml",
                "environment_sha256": "5" * 64,
                "lock_path": "workflow/envs/runner.lock",
                "lock_sha256": "6" * 64,
                "packages": [package["url"]],
            },
        ),
    )
    return json.loads(content)


def _canonical_index(document: dict[str, object]) -> bytes:
    document["identity"] = canonical_identity(
        without_key(document, "identity"),
        scheme=native_module.ENCODE_RUNTIME_INDEX_IDENTITY_SCHEME,
    )
    return canonical_json_bytes(document)


def test_encode_runtime_index_accepts_the_canonical_closure() -> None:
    index = parse_encode_runtime_index(_canonical_index(_index_document()))
    assert len(index.packages) == 1


@pytest.mark.parametrize(
    "content",
    [
        b"",
        b"x" * (native_module._MAX_ENCODE_INDEX_BYTES + 1),
        b"not json",
        b"[1]",
    ],
    ids=["empty", "over-size-limit", "not-json", "non-object"],
)
def test_encode_runtime_index_rejects_unreadable_content(content: bytes) -> None:
    with pytest.raises(native_module._NativeContractFault):
        parse_encode_runtime_index(content)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: {key: d[key] for key in d if key != "packages"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "workflow_build_identity": "bad"},
        lambda d: {**d, "workflow_build_identity": 1},
        lambda d: {**d, "micromamba": None},
        lambda d: {**d, "micromamba": {**d["micromamba"], "path": "other"}},
        lambda d: {**d, "micromamba": {**d["micromamba"], "size_bytes": 0}},
        lambda d: {**d, "micromamba": {**d["micromamba"], "sha256": "z" * 64}},
        lambda d: {**d, "packages": "pkg"},
        lambda d: {**d, "packages": []},
        lambda d: {**d, "packages": ["pkg"]},
        lambda d: {**d, "packages": [{**d["packages"][0], "extra": 1}]},
    ],
)
def test_encode_runtime_index_rejects_malformed_documents(mutation) -> None:
    document = mutation(_index_document())
    with pytest.raises(native_module._NativeContractFault):
        parse_encode_runtime_index(_canonical_index(document))


def test_encode_runtime_index_rejects_noncanonical_content() -> None:
    document = _index_document()
    content = json.dumps(document, indent=2).encode("utf-8")
    with pytest.raises(native_module._NativeContractFault):
        parse_encode_runtime_index(content)


def test_encode_runtime_index_rejects_identity_mismatch() -> None:
    document = _index_document()
    document["identity"] = f"sha256-{'f' * 64}"
    with pytest.raises(native_module._NativeContractFault):
        parse_encode_runtime_index(canonical_json_bytes(document))
