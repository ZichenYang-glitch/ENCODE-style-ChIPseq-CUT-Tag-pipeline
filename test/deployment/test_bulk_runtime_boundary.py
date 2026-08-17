from __future__ import annotations

import json
import os
from pathlib import Path
import runpy
import shutil
import subprocess
import sys

import pytest

from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.operator_action import (
    BulkRuntimePrepareReceipt,
    BulkRuntimePrepareRequest,
)


ROOT = Path(__file__).resolve().parents[2]
TEMPLATES = ROOT / "src" / "encode_pipeline" / "deployment" / "templates"
HELPER = TEMPLATES / "helixweave-bulk-runtime-prepare"
UNIT = TEMPLATES / "helixweave-bulk-runtime-prepare.service"
BOOTSTRAP = ROOT / "scripts" / "bootstrap_helixweave_operator.py"
PLATFORM_IDENTITY = f"sha256-{'a' * 64}"
BULK_IDENTITY = f"sha256-{'b' * 64}"
PRIOR_IDENTITY = f"sha256-{'c' * 64}"
CANDIDATE_IDENTITY = f"sha256-{'d' * 64}"
SERVICE_IDENTITY = f"sha256-{'e' * 64}"
CLIENT_IDENTITY = f"sha256-{'f' * 64}"
ENDPOINT_IDENTITY = f"sha256-{'1' * 64}"
RUNTIME_IDENTITY = f"sha256-{'2' * 64}"
IMAGE_SET_IDENTITY = f"sha256-{'3' * 64}"
TASK_IDENTITY = f"task-{'4' * 32}"


def _helper() -> dict[str, object]:
    return runpy.run_path(str(HELPER))


def _request(operation: str) -> BulkRuntimePrepareRequest:
    return BulkRuntimePrepareRequest.create(
        operation=operation,
        task_identity=TASK_IDENTITY,
        candidate_bulk_identity=BULK_IDENTITY,
        authority_platform_identity=PLATFORM_IDENTITY,
        prior_state_identity=PRIOR_IDENTITY,
        candidate_state_identity=CANDIDATE_IDENTITY,
        docker_service_identity=SERVICE_IDENTITY,
        docker_client_identity=CLIENT_IDENTITY,
        docker_endpoint_identity=ENDPOINT_IDENTITY,
        docker_daemon_uid=1234,
        docker_daemon_gid=1235,
    )


@pytest.mark.parametrize("operation", ("activate", "rollback", "verify"))
def test_bulk_runtime_helper_accepts_only_canonical_v2_operations(
    operation: str,
) -> None:
    helper = _helper()
    request = _request(operation)

    parsed = helper["_parse_request"](canonical_json_bytes(request.to_dict()))

    assert parsed == request.to_dict()
    with pytest.raises(ValueError):
        helper["_parse_request"](
            json.dumps(request.to_dict(), sort_keys=True).encode("utf-8")
        )


def test_bulk_runtime_helper_rejects_path_fields_duplicates_and_boolean_ids() -> None:
    helper = _helper()
    request = _request("activate").to_dict()
    request["candidate_bulk_identity"] = "../../private/runtime"
    unsigned = {name: item for name, item in request.items() if name != "identity"}
    request["identity"] = helper["_identity"](
        unsigned,
        helper["REQUEST_IDENTITY_SCHEME"],
    )
    with pytest.raises(ValueError):
        helper["_parse_request"](helper["_canonical"](request))

    request = _request("activate").to_dict()
    duplicated = canonical_json_bytes(request)
    duplicated = duplicated.replace(
        b'"component":"bulk-rnaseq-runtime",',
        b'"component":"bulk-rnaseq-runtime","component":"bulk-rnaseq-runtime",',
    )
    with pytest.raises(ValueError):
        helper["_parse_request"](duplicated)

    request = _request("activate").to_dict()
    request["docker_daemon_uid"] = True
    unsigned = {name: item for name, item in request.items() if name != "identity"}
    request["identity"] = helper["_identity"](
        unsigned,
        helper["REQUEST_IDENTITY_SCHEME"],
    )
    with pytest.raises(ValueError):
        helper["_parse_request"](helper["_canonical"](request))


def test_bulk_runtime_helper_validates_path_free_identity_bound_receipt() -> None:
    helper = _helper()
    request = _request("verify")
    receipt = BulkRuntimePrepareReceipt.create(
        request_identity=request.identity,
        candidate_bulk_identity=request.candidate_bulk_identity,
        runtime_identity=RUNTIME_IDENTITY,
        image_set_identity=IMAGE_SET_IDENTITY,
        image_count=3,
    )
    content = canonical_json_bytes(receipt.to_dict())

    helper["_validate_receipt"](content, request.to_dict())

    changed = receipt.to_dict()
    changed["runtime_identity"] = "/private/runtime"
    unsigned = {name: item for name, item in changed.items() if name != "identity"}
    changed["identity"] = helper["_identity"](
        unsigned,
        helper["RECEIPT_IDENTITY_SCHEME"],
    )
    with pytest.raises(ValueError):
        helper["_validate_receipt"](
            helper["_canonical"](changed),
            request.to_dict(),
        )


def test_bulk_runtime_helper_has_fixed_empty_process_boundary() -> None:
    content = HELPER.read_text(encoding="utf-8")

    assert (
        'REQUEST = Path("/run/helixweave/operator/bulk-runtime/request.json")'
        in content
    )
    assert (
        'RECEIPT = Path("/run/helixweave/operator/bulk-runtime/receipt.json")'
        in content
    )
    assert 'RELEASES = Path("/opt/helixweave/releases/platform")' in content
    assert "ENVIRONMENT = {}" in content
    assert '(str(launcher), "bulk-runtime-prepare")' in content
    assert "timeout=14500" in content
    assert "shell=True" not in content
    assert 'getattr(os, "O_NOFOLLOW", 0)' in content
    assert "stat.S_IMODE(before.st_mode) != 0o640" in content
    assert "before.st_uid != 0" in content
    assert "before.st_gid != os.getegid()" in content


def test_bulk_runtime_oneshot_exposes_only_fixed_runtime_socket_and_receipt() -> None:
    content = UNIT.read_text(encoding="utf-8")

    for directive in (
        "User=helixweave",
        "Group=helixweave",
        "UMask=0077",
        "ExecStart=/usr/libexec/helixweave-bulk-runtime-prepare",
        "TimeoutStartSec=14600s",
        "RestrictAddressFamilies=AF_UNIX",
        "ReadOnlyPaths=/opt/helixweave/releases/platform",
        "ReadOnlyPaths=/opt/helixweave/runtimes/bulk-rnaseq",
        "ReadOnlyPaths=/usr/bin/docker",
        "ReadOnlyPaths=/run/helixweave/operator/bulk-runtime/request.json",
        "ReadOnlyPaths=/run/helixweave/docker/docker.sock",
        "ReadWritePaths=/run/helixweave/operator/bulk-runtime/receipt.json",
        "InaccessiblePaths=/var/lib/helixweave/docker-rootless",
    ):
        assert directive in content
    assert "AF_INET" not in content
    assert "Environment=" not in content
    assert "ReadWritePaths=/opt" not in content
    assert "ReadWritePaths=/var/lib" not in content


def test_bulk_runtime_boundary_is_packaged_snapshotted_and_tmpfiles_owned() -> None:
    bootstrap = runpy.run_path(str(BOOTSTRAP))
    boundary = {item.source_name: item for item in bootstrap["BOUNDARY_SPECS"]}
    closure = {item.destination.as_posix(): item for item in bootstrap["CLOSURE_SPECS"]}

    assert boundary["helixweave-bulk-runtime-prepare"].snapshot is True
    assert boundary["helixweave-bulk-runtime-prepare"].destination == Path(
        "/usr/libexec/helixweave-bulk-runtime-prepare"
    )
    assert boundary["helixweave-bulk-runtime-prepare.service"].snapshot is False
    assert boundary["helixweave-bulk-runtime-prepare.service"].destination == Path(
        "/usr/lib/systemd/system/helixweave-bulk-runtime-prepare.service"
    )
    assert "boundary/helixweave-bulk-runtime-prepare" in closure
    assert "templates/helixweave-bulk-runtime-prepare.service" in closure
    assert "lib/encode_pipeline/deployment/bulk_docker_boundary.py" in closure
    assert "d /run/helixweave/operator/bulk-runtime 0750 root helixweave -" in (
        TEMPLATES / "helixweave.tmpfiles.conf"
    ).read_text(encoding="utf-8")


def test_bulk_runtime_helper_bytes_change_operator_closure_identity(
    tmp_path: Path,
) -> None:
    bootstrap = runpy.run_path(str(BOOTSTRAP))
    identities: list[str] = []
    for name, changed in (("original", False), ("changed", True)):
        source = tmp_path / f"{name}-source"
        shutil.copytree(TEMPLATES, source)
        if changed:
            helper = source / "helixweave-bulk-runtime-prepare"
            helper.write_bytes(helper.read_bytes() + b"\n# reviewed boundary change\n")
        root = tmp_path / name
        root.mkdir()
        backend = bootstrap["HostBootstrapBackend"](
            source_root=source,
            root_prefix=root,
            owner_uid=os.getuid(),
            owner_gid=os.getgid(),
            command_runner=lambda command: pytest.fail(
                f"unexpected command: {command}"
            ),
            sudoers_validator=lambda _path: True,
        )
        identities.append(
            backend.apply(
                operation="install",
                invoking_user="labadmin",
            ).closure_identity
        )

    assert identities[0] != identities[1]


def test_operator_closure_imports_bulk_boundary_under_isolated_python(
    tmp_path: Path,
) -> None:
    bootstrap = runpy.run_path(str(BOOTSTRAP))
    root = tmp_path / "host"
    root.mkdir()
    backend = bootstrap["HostBootstrapBackend"](
        source_root=TEMPLATES,
        root_prefix=root,
        owner_uid=os.getuid(),
        owner_gid=os.getgid(),
        command_runner=lambda command: pytest.fail(f"unexpected command: {command}"),
        sudoers_validator=lambda _path: True,
    )
    installed = backend.apply(operation="install", invoking_user="labadmin")
    closure = root / "opt/helixweave/operator" / installed.closure_identity
    script = closure / "bin/helixweave-operator"

    completed = subprocess.run(
        (sys.executable, "-I", "-S", "-E", "-B", str(script), "status"),
        check=False,
        capture_output=True,
        text=True,
    )

    assert (
        closure / "lib/encode_pipeline/deployment/bulk_docker_boundary.py"
    ).is_file()
    assert completed.returncode == 77
    assert json.loads(completed.stderr)["issue"]["code"] == (
        "OPERATOR_PRIVILEGE_REQUIRED"
    )
