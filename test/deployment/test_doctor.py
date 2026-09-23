from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from encode_pipeline.deployment.admission import DeferredNativeContractResolver
from encode_pipeline.deployment.doctor import (
    CHECKS,
    FAIL,
    PASS,
    WARNING,
    DeploymentDoctor,
    DeploymentSnapshot,
    DeploymentStateProbe,
    DoctorReport,
    ProbeResult,
    RuntimeProbe,
    fixed_probe,
)
from encode_pipeline.deployment.errors import DeploymentError, fail
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.manager import DeploymentManager, DeploymentOwnership
from encode_pipeline.deployment.models import ENCODE_RUNTIME, PLATFORM
from .support import manager_for, manifest_for, write_bundle


IDENTITY = f"sha256-{'a' * 64}"


def _probes(*, state: str = PASS, reason: str = "READY"):
    return {
        check_id: fixed_probe(state, reason, IDENTITY) for check_id, _category in CHECKS
    }


@pytest.mark.parametrize("check_id", ["database", "references"])
@pytest.mark.parametrize(
    "failure", [None, "unexpected", "deployment", "public-deployment"]
)
def test_doctor_probe_failure_diagnostic_preserves_fallback(
    check_id, failure, tmp_path, monkeypatch, assert_failure_diagnostics
):
    private = str(tmp_path / "PRIVATE_DOCTOR_PATH")
    monkeypatch.setenv("HELIX_TEST_PRIVATE", "PRIVATE_DOCTOR_ENV")
    calls = []

    def probe():
        calls.append(check_id)
        if failure == "unexpected":
            raise RuntimeError(private + " PRIVATE_DOCTOR_EXCEPTION")
        if failure == "deployment":
            raise fail("PRIVATE_DOCTOR_CODE", private + " PRIVATE_DOCTOR_EXCEPTION")
        if failure == "public-deployment":
            raise fail("RUNTIME_NOT_ACTIVE", private + " PRIVATE_DOCTOR_EXCEPTION")
        return ProbeResult(PASS, "READY")

    probes = _probes()
    probes[check_id] = probe
    report = DeploymentDoctor(probes).run()
    checks = {item.check_id: item for item in report.checks}
    assert calls == [check_id]
    assert checks[check_id].state == (FAIL if failure else PASS)
    assert checks[check_id].reason_code == (
        "RUNTIME_NOT_ACTIVE"
        if failure == "public-deployment"
        else "DOCTOR_CHECK_FAILED"
        if failure
        else "READY"
    )
    assert report.ready is (failure is None)
    assert len(checks) == len(CHECKS)
    assert all(item.state == PASS for key, item in checks.items() if key != check_id)
    assert_failure_diagnostics(
        [("deployment_doctor", check_id, "DOCTOR_CHECK_FAILED")] if failure else [],
        private=(
            private,
            "PRIVATE_DOCTOR_EXCEPTION",
            "PRIVATE_DOCTOR_ENV",
            "PRIVATE_DOCTOR_CODE",
        ),
    )


def test_doctor_report_has_one_stable_check_inventory_and_status() -> None:
    probes = _probes()
    probes["references"] = fixed_probe(WARNING, "REFERENCES_INCOMPLETE")

    report = DeploymentDoctor(probes).run()

    assert report.status == "degraded"
    assert report.ready is False
    assert [item["check_id"] for item in report.to_dict()["checks"]] == [
        check_id for check_id, _category in CHECKS
    ]
    assert report.to_dict()["schema_version"] == "helixweave-deployment-doctor-v1"


def test_doctor_validates_probe_results_inventory_and_healthy_status() -> None:
    with pytest.raises(DeploymentError) as invalid_result:
        ProbeResult("unknown", "READY")
    assert invalid_result.value.issue.code == "DOCTOR_RESULT_INVALID"

    with pytest.raises(DeploymentError) as invalid_report:
        DoctorReport.create(())
    assert invalid_report.value.issue.code == "DOCTOR_RESULT_INVALID"

    with pytest.raises(DeploymentError) as invalid_inventory:
        DeploymentDoctor({})
    assert invalid_inventory.value.issue.code == "DOCTOR_PROBES_INVALID"

    report = DeploymentDoctor(_probes()).run()
    assert report.status == "healthy"
    assert report.ready is True


def test_doctor_maps_private_and_wrong_type_probe_failures_to_public_reasons() -> None:
    def private_failure() -> ProbeResult:
        raise fail("PRIVATE_PROBE_FAILURE", "/private/path token=secret")

    probes = _probes()
    probes["database"] = private_failure
    probes["redis"] = lambda: object()  # type: ignore[assignment]

    report = DeploymentDoctor(probes).run()
    checks = {item.check_id: item for item in report.checks}

    assert checks["database"].state == FAIL
    assert checks["database"].reason_code == "DOCTOR_CHECK_FAILED"
    assert checks["redis"].state == FAIL
    assert checks["redis"].reason_code == "DOCTOR_RESULT_INVALID"
    assert "private" not in json.dumps(report.to_dict())
    assert "secret" not in json.dumps(report.to_dict())


def test_doctor_fails_closed_and_redacts_unexpected_probe_errors(
    tmp_path: Path,
) -> None:
    private = tmp_path / "reference" / "secret.fa"

    def broken() -> ProbeResult:
        raise RuntimeError(f"failed at {private} token=private")

    probes = _probes()
    probes["database"] = broken

    report = DeploymentDoctor(probes).run()
    rendered = json.dumps(report.to_dict())

    database = next(item for item in report.checks if item.check_id == "database")
    assert report.status == "unhealthy"
    assert database.state == FAIL
    assert database.reason_code == "DOCTOR_CHECK_FAILED"
    assert str(tmp_path) not in rendered
    assert "secret.fa" not in rendered
    assert "token=" not in rendered


def test_deployment_state_probe_is_read_only(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = DeploymentManager(
        layout,
        ownership=DeploymentOwnership(os.getuid(), os.getgid()),
    )
    manifest, payload = manifest_for(PLATFORM)
    bundle = write_bundle(tmp_path / "platform.tar", manifest, payload)
    manager.stage(
        bundle,
        expected_owner_uid=os.getuid(),
    )
    before = _tree_identity(tmp_path / "host")

    snapshot = DeploymentSnapshot(manager)
    result = DeploymentStateProbe(snapshot)()

    assert result.state == PASS
    assert result.reason_code == "DEPLOYMENT_STATE_READY"
    assert _tree_identity(tmp_path / "host") == before


def test_deployment_and_runtime_probes_share_one_lightweight_snapshot(
    tmp_path: Path,
    monkeypatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = DeploymentManager(
        layout,
        ownership=DeploymentOwnership(os.getuid(), os.getgid()),
    )
    manifest, payload = manifest_for(PLATFORM)
    manager.stage(
        write_bundle(tmp_path / "platform.tar", manifest, payload),
        expected_owner_uid=os.getuid(),
        expected_owner_gid=os.getgid(),
    )
    calls = 0
    original = manager.status

    def counted_status():
        nonlocal calls
        calls += 1
        return original()

    monkeypatch.setattr(manager, "status", counted_status)
    snapshot = DeploymentSnapshot(manager)

    DeploymentStateProbe(snapshot)()
    RuntimeProbe(snapshot, "encode-runtime")()
    RuntimeProbe(snapshot, "bulk-rnaseq-runtime")()

    assert calls == 1


def test_deployment_state_probe_reports_interrupted_state_without_live_checks() -> None:
    snapshot = SimpleNamespace(
        read=lambda: SimpleNamespace(
            interrupted=True,
            state=SimpleNamespace(identity=IDENTITY),
        )
    )

    result = DeploymentStateProbe(snapshot)  # type: ignore[arg-type]

    assert result().state == WARNING
    assert result().reason_code == "DEPLOYMENT_INTERRUPTED"
    assert result().evidence_identity == IDENTITY


def test_runtime_probe_admits_active_contract_once_and_rejects_invalid_scope(
    tmp_path: Path,
) -> None:
    with pytest.raises(DeploymentError) as invalid_component:
        RuntimeProbe(SimpleNamespace(), PLATFORM)  # type: ignore[arg-type]
    assert invalid_component.value.issue.code == "DOCTOR_PROBES_INVALID"

    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(ENCODE_RUNTIME)
    manager.stage(write_bundle(tmp_path / "runtime.tar", manifest, payload))
    manager.activate(ENCODE_RUNTIME, expected_staged_identity=manifest.identity)
    snapshot = DeploymentSnapshot(manager)

    first = snapshot.admit(ENCODE_RUNTIME)
    second = snapshot.admit(ENCODE_RUNTIME)
    result = RuntimeProbe(snapshot, ENCODE_RUNTIME)()

    assert second is first
    assert result == ProbeResult(PASS, "RUNTIME_READY", manifest.identity)

    empty_manager = manager_for(DeploymentLayout.isolated(tmp_path / "empty"))
    empty_manager.stage(
        write_bundle(tmp_path / "inactive-runtime.tar", manifest, payload)
    )
    empty = DeploymentSnapshot(empty_manager)
    with pytest.raises(DeploymentError) as inactive:
        empty.admit(ENCODE_RUNTIME)
    assert inactive.value.issue.code == "RUNTIME_NOT_ACTIVE"


def test_runtime_probe_requires_native_contract_admission(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(ENCODE_RUNTIME)
    manager.stage(write_bundle(tmp_path / "runtime.tar", manifest, payload))
    manager.activate(
        ENCODE_RUNTIME,
        expected_staged_identity=manifest.identity,
    )
    manager.contract_resolver = DeferredNativeContractResolver()
    probes = _probes()
    probes["encode-runtime"] = RuntimeProbe(DeploymentSnapshot(manager), ENCODE_RUNTIME)

    report = DeploymentDoctor(probes).run()

    runtime = next(item for item in report.checks if item.check_id == "encode-runtime")
    assert runtime.state == FAIL
    assert runtime.reason_code == "DEPLOYMENT_CONTRACT_ADMISSION_DEFERRED"


def _tree_identity(root: Path) -> tuple[tuple[object, ...], ...]:
    values: list[tuple[object, ...]] = []
    for path in sorted((root, *root.rglob("*"))):
        observed = path.lstat()
        relative = path.relative_to(root).as_posix() if path != root else "."
        if path.is_symlink():
            content_identity = f"link:{os.readlink(path)}"
        elif path.is_file():
            content_identity = hashlib.sha256(path.read_bytes()).hexdigest()
        else:
            content_identity = "directory"
        values.append(
            (
                relative,
                observed.st_mode,
                observed.st_uid,
                observed.st_gid,
                observed.st_size,
                observed.st_mtime_ns,
                content_identity,
            )
        )
    return tuple(values)
