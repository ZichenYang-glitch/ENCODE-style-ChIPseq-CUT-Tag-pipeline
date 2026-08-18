"""Behavior coverage for the fixed systemd action runner boundaries."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

import encode_pipeline.deployment.operator_action as operator_action_module
from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import COMPONENTS
from encode_pipeline.deployment.operator_action import (
    ACTION_UNIT,
    DATABASE_PREPARE_UNIT,
    VERIFICATION_CHECKS,
    DatabasePrepareReceipt,
    DatabasePrepareRequest,
    DeploymentActionReceipt,
    DeploymentActionRequest,
    ReadinessCheck,
    SystemdDatabasePreparer,
    SystemdDeploymentActionRunner,
)

_ID = f"sha256-{'a' * 64}"
_ID_B = f"sha256-{'b' * 64}"
_ID_C = f"sha256-{'c' * 64}"
_ID_D = f"sha256-{'d' * 64}"
_TASK = f"task-{'a' * 32}"
_SHA256 = "a" * 64
_HEAD = "20260809_13"
_UID = os.getuid()
_GID = os.getgid()


def _reject(callable_, *args, code: str, **kwargs) -> None:
    with pytest.raises(DeploymentError) as captured:
        callable_(*args, **kwargs)
    assert captured.value.issue.code == code


def _layout(tmp_path: Path, *, action_mode: int = 0o750, prepare_mode: int = 0o710):
    layout = DeploymentLayout.isolated(tmp_path / "host")
    layout.operator_action_root.mkdir(parents=True)
    layout.operator_action_root.chmod(action_mode)
    layout.database_prepare_root.mkdir(parents=True)
    layout.database_prepare_root.chmod(prepare_mode)
    return layout


def _action_runner(
    layout: DeploymentLayout, command_runner
) -> SystemdDeploymentActionRunner:
    return SystemdDeploymentActionRunner(
        layout,
        service_uid=_UID,
        service_gid=_GID,
        root_uid=_UID,
        root_gid=_GID,
        command_runner=command_runner,
    )


def _database_preparer(layout: DeploymentLayout, command_runner):
    return SystemdDatabasePreparer(
        layout,
        service_uid=_UID,
        service_gid=_GID,
        root_uid=_UID,
        command_runner=command_runner,
    )


def _action_request() -> DeploymentActionRequest:
    return DeploymentActionRequest.create(
        phase="observe",
        operation="observe",
        component="platform",
        task_identity=_TASK,
        deployment_identity=_ID,
        authority_platform_identity=_ID_B,
        prior_state_identity=_ID_C,
        candidate_state_identity=_ID_D,
        candidate_active={component: _ID for component in COMPONENTS},
    )


def _action_receipt(request: DeploymentActionRequest) -> DeploymentActionReceipt:
    return DeploymentActionReceipt.create(
        request_identity=request.identity,
        status="observed",
        compatibility="compatible",
        accepted_schema_heads=(_HEAD,),
        target_schema_heads=(_HEAD,),
        migration_inventory_identity=_ID_B,
        known_schema_revisions=("20260807_12",),
        migration_required=False,
        rollback_supported=True,
        api_contract_sha256=_SHA256,
        native_identities={component: _ID for component in COMPONENTS},
        frontend_identity=_ID_C,
        reference_compatibility_identity=_ID_D,
        readiness={
            name: ReadinessCheck(status="ready", reason_code="READY", identity=_ID)
            for name in VERIFICATION_CHECKS
        },
    )


def _database_request() -> DatabasePrepareRequest:
    return DatabasePrepareRequest.create(
        operation="activate",
        database_mode="existing-live",
        task_identity=_TASK,
        deployment_identity=_ID,
        prior_state_identity=_ID_B,
        candidate_state_identity=_ID_C,
        action_receipt_identity=_ID_D,
        backup_receipt_identity=_ID,
        target_schema_heads=(_HEAD,),
    )


def _database_receipt(request: DatabasePrepareRequest) -> DatabasePrepareReceipt:
    return DatabasePrepareReceipt.create(
        request_identity=request.identity,
        database_before_identity=_ID_B,
        database_after_identity=_ID_C,
        schema_heads=(_HEAD,),
    )


def _write_receipt(path: Path, content: bytes) -> None:
    path.write_bytes(content)
    path.chmod(0o600)


@pytest.mark.parametrize("value", [True, -1, "0"])
def test_action_runner_rejects_invalid_identity_values(tmp_path: Path, value) -> None:
    layout = _layout(tmp_path)
    _reject(
        SystemdDeploymentActionRunner,
        layout,
        service_uid=value,
        service_gid=_GID,
        root_uid=_UID,
        root_gid=_GID,
        code="OPERATOR_ACTION_BOUNDARY_INVALID",
    )


@pytest.mark.parametrize("value", [True, -1, "0"])
def test_database_preparer_rejects_invalid_identity_values(
    tmp_path: Path, value
) -> None:
    layout = _layout(tmp_path)
    _reject(
        SystemdDatabasePreparer,
        layout,
        service_uid=value,
        service_gid=_GID,
        root_uid=_UID,
        code="DATABASE_PREPARE_BOUNDARY_INVALID",
    )


def test_action_runner_requires_one_typed_request(tmp_path: Path) -> None:
    runner = _action_runner(_layout(tmp_path), lambda argv: 0)
    _reject(
        runner.run, {"operation": "observe"}, code="OPERATOR_ACTION_REQUEST_INVALID"
    )


def test_database_preparer_requires_one_typed_request(tmp_path: Path) -> None:
    preparer = _database_preparer(_layout(tmp_path), lambda argv: 0)
    _reject(
        preparer.prepare,
        {"operation": "activate"},
        code="DATABASE_PREPARE_REQUEST_INVALID",
    )


def test_action_runner_requires_the_owned_root(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    runner = _action_runner(layout, lambda argv: 0)
    _reject(runner.run, _action_request(), code="OPERATOR_ACTION_BOUNDARY_INVALID")

    layout = _layout(tmp_path / "mode", action_mode=0o700)
    runner = _action_runner(layout, lambda argv: 0)
    _reject(runner.run, _action_request(), code="OPERATOR_ACTION_BOUNDARY_INVALID")

    layout = DeploymentLayout.isolated(tmp_path / "linked" / "host")
    target = tmp_path / "linked" / "target"
    target.mkdir(parents=True)
    target.chmod(0o750)
    layout.operator_action_root.parent.mkdir(parents=True)
    layout.operator_action_root.symlink_to(target, target_is_directory=True)
    runner = _action_runner(layout, lambda argv: 0)
    _reject(runner.run, _action_request(), code="OPERATOR_ACTION_BOUNDARY_INVALID")


def test_database_preparer_requires_the_owned_root(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    preparer = _database_preparer(layout, lambda argv: 0)
    _reject(
        preparer.prepare,
        _database_request(),
        code="DATABASE_PREPARE_BOUNDARY_INVALID",
    )

    layout = _layout(tmp_path / "mode", prepare_mode=0o750)
    preparer = _database_preparer(layout, lambda argv: 0)
    _reject(
        preparer.prepare,
        _database_request(),
        code="DATABASE_PREPARE_BOUNDARY_INVALID",
    )


def test_action_runner_exchanges_canonical_documents(tmp_path: Path) -> None:
    layout = _layout(tmp_path)
    request = _action_request()
    observed_argv: list[tuple[str, ...]] = []

    def fake_runner(argv: tuple[str, ...]) -> int:
        observed_argv.append(argv)
        _write_receipt(
            layout.operator_action_receipt,
            canonical_json_bytes(_action_receipt(request).to_dict()),
        )
        return 0

    receipt = _action_runner(layout, fake_runner).run(request)

    assert receipt.request_identity == request.identity
    assert observed_argv == [
        (
            str(operator_action_module.SYSTEMCTL),
            "--no-ask-password",
            "start",
            "--",
            ACTION_UNIT,
        )
    ]
    staged = layout.operator_action_request.read_bytes()
    assert staged == canonical_json_bytes(request.to_dict())


def test_database_preparer_exchanges_canonical_documents(tmp_path: Path) -> None:
    layout = _layout(tmp_path)
    request = _database_request()
    observed_argv: list[tuple[str, ...]] = []

    def fake_runner(argv: tuple[str, ...]) -> int:
        observed_argv.append(argv)
        _write_receipt(
            layout.database_prepare_receipt,
            canonical_json_bytes(_database_receipt(request).to_dict()),
        )
        return 0

    receipt = _database_preparer(layout, fake_runner).prepare(request)

    assert receipt.request_identity == request.identity
    assert observed_argv == [
        (
            str(operator_action_module.SYSTEMCTL),
            "--no-ask-password",
            "start",
            "--",
            DATABASE_PREPARE_UNIT,
        )
    ]
    staged = layout.database_prepare_request.read_bytes()
    assert staged == canonical_json_bytes(request.to_dict())


def test_action_runner_fails_closed_on_runner_exception(tmp_path: Path) -> None:
    def raising_runner(argv: tuple[str, ...]) -> int:
        raise RuntimeError("systemctl unavailable")

    runner = _action_runner(_layout(tmp_path), raising_runner)
    _reject(runner.run, _action_request(), code="OPERATOR_ACTION_FAILED")


def test_action_runner_fails_closed_on_nonzero_return(tmp_path: Path) -> None:
    runner = _action_runner(_layout(tmp_path), lambda argv: 3)
    _reject(runner.run, _action_request(), code="OPERATOR_ACTION_FAILED")


def test_database_preparer_fails_closed_on_runner_exception(tmp_path: Path) -> None:
    def raising_runner(argv: tuple[str, ...]) -> int:
        raise RuntimeError("systemctl unavailable")

    preparer = _database_preparer(_layout(tmp_path), raising_runner)
    _reject(preparer.prepare, _database_request(), code="DATABASE_PREPARE_FAILED")


def test_database_preparer_fails_closed_on_nonzero_return(tmp_path: Path) -> None:
    preparer = _database_preparer(_layout(tmp_path), lambda argv: 1)
    _reject(preparer.prepare, _database_request(), code="DATABASE_PREPARE_FAILED")


def test_action_runner_rejects_a_missing_receipt(tmp_path: Path) -> None:
    runner = _action_runner(_layout(tmp_path), lambda argv: 0)
    _reject(runner.run, _action_request(), code="OPERATOR_ACTION_RECEIPT_INVALID")


def test_action_runner_rejects_malformed_receipt_content(tmp_path: Path) -> None:
    layout = _layout(tmp_path)

    def fake_runner(argv: tuple[str, ...]) -> int:
        _write_receipt(layout.operator_action_receipt, b"not json")
        return 0

    runner = _action_runner(layout, fake_runner)
    _reject(runner.run, _action_request(), code="OPERATOR_ACTION_RECEIPT_INVALID")


def test_action_runner_rejects_noncanonical_receipt_content(tmp_path: Path) -> None:
    import json

    layout = _layout(tmp_path)
    request = _action_request()

    def fake_runner(argv: tuple[str, ...]) -> int:
        _write_receipt(
            layout.operator_action_receipt,
            json.dumps(_action_receipt(request).to_dict(), indent=2).encode("utf-8"),
        )
        return 0

    runner = _action_runner(layout, fake_runner)
    _reject(runner.run, request, code="OPERATOR_ACTION_RECEIPT_INVALID")


def test_action_runner_rejects_a_foreign_receipt(tmp_path: Path) -> None:
    layout = _layout(tmp_path)
    other = DeploymentActionRequest.create(
        phase="observe",
        operation="observe",
        component="platform",
        task_identity=_TASK,
        deployment_identity=_ID,
        authority_platform_identity=_ID_B,
        prior_state_identity=_ID_D,
        candidate_state_identity=_ID_C,
        candidate_active={component: _ID for component in COMPONENTS},
    )

    def fake_runner(argv: tuple[str, ...]) -> int:
        _write_receipt(
            layout.operator_action_receipt,
            canonical_json_bytes(_action_receipt(other).to_dict()),
        )
        return 0

    runner = _action_runner(layout, fake_runner)
    _reject(runner.run, _action_request(), code="OPERATOR_ACTION_RECEIPT_INVALID")


def test_database_preparer_rejects_a_missing_receipt(tmp_path: Path) -> None:
    preparer = _database_preparer(_layout(tmp_path), lambda argv: 0)
    _reject(
        preparer.prepare,
        _database_request(),
        code="DATABASE_PREPARE_RECEIPT_INVALID",
    )


def test_database_preparer_rejects_an_invalid_receipt_document(tmp_path: Path) -> None:
    layout = _layout(tmp_path)

    def fake_runner(argv: tuple[str, ...]) -> int:
        _write_receipt(layout.database_prepare_receipt, b"{}")
        return 0

    preparer = _database_preparer(layout, fake_runner)
    _reject(
        preparer.prepare,
        _database_request(),
        code="DATABASE_PREPARE_RECEIPT_INVALID",
    )


def test_database_preparer_rejects_a_foreign_receipt(tmp_path: Path) -> None:
    layout = _layout(tmp_path)
    other = DatabasePrepareRequest.create(
        operation="activate",
        database_mode="existing-live",
        task_identity=_TASK,
        deployment_identity=_ID,
        prior_state_identity=_ID_C,
        candidate_state_identity=_ID_B,
        action_receipt_identity=_ID_D,
        backup_receipt_identity=_ID,
        target_schema_heads=(_HEAD,),
    )

    def fake_runner(argv: tuple[str, ...]) -> int:
        _write_receipt(
            layout.database_prepare_receipt,
            canonical_json_bytes(_database_receipt(other).to_dict()),
        )
        return 0

    preparer = _database_preparer(layout, fake_runner)
    _reject(
        preparer.prepare,
        _database_request(),
        code="DATABASE_PREPARE_RECEIPT_INVALID",
    )


def test_action_runner_systemctl_command_is_fixed(monkeypatch) -> None:
    observed: list[tuple[str, str]] = []
    monkeypatch.setattr(
        operator_action_module,
        "_run_fixed_systemctl",
        lambda argv, unit: observed.append((argv, unit)) or 0,
    )
    argv = ("/bin/systemctl", "--no-ask-password", "start", "--", ACTION_UNIT)
    assert SystemdDeploymentActionRunner._run_systemctl(argv) == 0
    assert observed == [(argv, ACTION_UNIT)]


def test_database_preparer_systemctl_command_is_fixed(monkeypatch) -> None:
    observed: list[tuple[str, str]] = []
    monkeypatch.setattr(
        operator_action_module,
        "_run_fixed_systemctl",
        lambda argv, unit: observed.append((argv, unit)) or 0,
    )
    argv = ("/bin/systemctl", "--no-ask-password", "start", "--", DATABASE_PREPARE_UNIT)
    assert SystemdDatabasePreparer._run_systemctl(argv) == 0
    assert observed == [(argv, DATABASE_PREPARE_UNIT)]
