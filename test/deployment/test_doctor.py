from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

from encode_pipeline.deployment.admission import DeferredNativeContractResolver
from encode_pipeline.deployment.doctor import (
    CHECKS,
    FAIL,
    PASS,
    WARNING,
    DeploymentDoctor,
    DeploymentSnapshot,
    DeploymentStateProbe,
    ProbeResult,
    RuntimeProbe,
    fixed_probe,
)
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.manager import DeploymentManager, DeploymentOwnership
from encode_pipeline.deployment.models import ENCODE_RUNTIME, PLATFORM
from .support import manager_for, manifest_for, write_bundle


IDENTITY = f"sha256-{'a' * 64}"


def _probes(*, state: str = PASS, reason: str = "READY"):
    return {
        check_id: fixed_probe(state, reason, IDENTITY) for check_id, _category in CHECKS
    }


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
