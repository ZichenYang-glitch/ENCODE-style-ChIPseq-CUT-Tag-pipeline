from __future__ import annotations

import hashlib
import os
from pathlib import Path
import stat
from types import SimpleNamespace

import pytest

from encode_pipeline.deployment import filesystem
from encode_pipeline.deployment.doctor import (
    FAIL,
    PASS,
    ActiveFrontendProbe,
    DatabaseProbe,
    DeploymentSnapshot,
    ProbeResult,
    frontend_probe,
)
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import (
    PLATFORM,
    BundleManifest,
    ComponentSlots,
    ContractDocument,
    ContractIdentity,
    ContractRequirement,
    DeploymentState,
    FileRecord,
)
from encode_pipeline.frontend_assets import load_packaged_frontend_assets
from .support import manager_for, manifest_for, write_bundle


CODE = "DEPLOYMENT_TEST"
IDENTITY_A = f"sha256-{'a' * 64}"
IDENTITY_B = f"sha256-{'b' * 64}"


def _raise_oserror(*args: object, **kwargs: object) -> None:
    raise OSError("injected storage failure")


# filesystem primitives


def test_require_directory_accepts_real_directory(tmp_path: Path) -> None:
    observed = filesystem.require_directory(tmp_path, code=CODE)

    assert stat.S_ISDIR(observed.st_mode)


def test_require_directory_rejects_missing_path(tmp_path: Path) -> None:
    with pytest.raises(DeploymentError) as captured:
        filesystem.require_directory(tmp_path / "missing", code=CODE)

    assert captured.value.issue.code == CODE


def test_require_directory_rejects_file_and_symlink(tmp_path: Path) -> None:
    file_path = tmp_path / "file"
    file_path.write_bytes(b"x")
    link_path = tmp_path / "link"
    link_path.symlink_to(tmp_path, target_is_directory=True)

    for rejected in (file_path, link_path):
        with pytest.raises(DeploymentError) as captured:
            filesystem.require_directory(rejected, code=CODE)
        assert captured.value.issue.code == CODE


def test_create_directory_rejects_relative_path() -> None:
    with pytest.raises(DeploymentError) as captured:
        filesystem.create_directory(Path("relative/tree"))

    assert captured.value.issue.code == "DEPLOYMENT_PATH_INVALID"


def test_create_directory_rejects_file_and_symlink_boundaries(
    tmp_path: Path,
) -> None:
    file_path = tmp_path / "file"
    file_path.write_bytes(b"x")
    link_path = tmp_path / "link"
    link_path.symlink_to(tmp_path, target_is_directory=True)

    for boundary in (file_path, link_path):
        with pytest.raises(DeploymentError) as captured:
            filesystem.create_directory(boundary / "child")
        assert captured.value.issue.code == "DEPLOYMENT_PATH_UNSAFE"


def test_create_directory_maps_mkdir_failure_to_public_code(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(Path, "mkdir", _raise_oserror)

    with pytest.raises(DeploymentError) as captured:
        filesystem.create_directory(tmp_path / "host")

    assert captured.value.issue.code == "DEPLOYMENT_STORAGE_UNAVAILABLE"


def test_create_directory_maps_lstat_failure_to_public_code(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(Path, "lstat", _raise_oserror)

    with pytest.raises(DeploymentError) as captured:
        filesystem.create_directory(tmp_path / "host")

    assert captured.value.issue.code == "DEPLOYMENT_STORAGE_UNAVAILABLE"


def test_read_regular_file_returns_content_and_stat(tmp_path: Path) -> None:
    target = tmp_path / "payload.bin"
    target.write_bytes(b"payload-bytes")

    content, observed = filesystem.read_regular_file(target, max_bytes=64, code=CODE)

    assert content == b"payload-bytes"
    assert stat.S_ISREG(observed.st_mode)


def test_read_regular_file_rejects_missing_and_symlink(tmp_path: Path) -> None:
    target = tmp_path / "payload.bin"
    target.write_bytes(b"x")
    link_path = tmp_path / "link"
    link_path.symlink_to(target)

    for rejected in (tmp_path / "missing", link_path):
        with pytest.raises(DeploymentError) as captured:
            filesystem.read_regular_file(rejected, max_bytes=64, code=CODE)
        assert captured.value.issue.code == CODE


def test_read_regular_file_rejects_directory_hardlink_and_oversize(
    tmp_path: Path,
) -> None:
    target = tmp_path / "payload.bin"
    target.write_bytes(b"x" * 16)
    hardlink = tmp_path / "hardlink"
    os.link(target, hardlink)

    for rejected, max_bytes in ((tmp_path, 64), (hardlink, 64), (target, 8)):
        with pytest.raises(DeploymentError) as captured:
            filesystem.read_regular_file(rejected, max_bytes=max_bytes, code=CODE)
        assert captured.value.issue.code == CODE


def test_read_regular_file_detects_witness_change(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "payload.bin"
    target.write_bytes(b"stable")
    witnesses = iter(((1,), (2,)))
    monkeypatch.setattr(filesystem, "_file_witness", lambda value: next(witnesses))

    with pytest.raises(DeploymentError) as captured:
        filesystem.read_regular_file(target, max_bytes=64, code=CODE)

    assert captured.value.issue.code == CODE


def test_hash_regular_file_returns_digest(tmp_path: Path) -> None:
    content = b"hash-me" * 128
    target = tmp_path / "payload.bin"
    target.write_bytes(content)

    digest, observed = filesystem.hash_regular_file(
        target, expected_size=len(content), code=CODE
    )

    assert digest == hashlib.sha256(content).hexdigest()
    assert stat.S_ISREG(observed.st_mode)


def test_hash_regular_file_rejects_missing_and_symlink(tmp_path: Path) -> None:
    target = tmp_path / "payload.bin"
    target.write_bytes(b"x")
    link_path = tmp_path / "link"
    link_path.symlink_to(target)

    for rejected in (tmp_path / "missing", link_path):
        with pytest.raises(DeploymentError) as captured:
            filesystem.hash_regular_file(rejected, expected_size=1, code=CODE)
        assert captured.value.issue.code == CODE


def test_hash_regular_file_rejects_growth_during_read() -> None:
    # /proc/self/status stats as a size-0 regular file but yields bytes when
    # read, which exercises the growth guard without mocking the syscall.
    with pytest.raises(DeploymentError) as captured:
        filesystem.hash_regular_file(
            Path("/proc/self/status"), expected_size=0, code=CODE
        )

    assert captured.value.issue.code == CODE


def test_hash_regular_file_detects_witness_change(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    target = tmp_path / "payload.bin"
    target.write_bytes(b"stable")
    witnesses = iter(((1,), (2,)))
    monkeypatch.setattr(filesystem, "_file_witness", lambda value: next(witnesses))

    with pytest.raises(DeploymentError) as captured:
        filesystem.hash_regular_file(target, expected_size=6, code=CODE)

    assert captured.value.issue.code == CODE


def test_fsync_directory_accepts_real_directory(tmp_path: Path) -> None:
    filesystem.fsync_directory(tmp_path)


def test_fsync_directory_rejects_missing_path(tmp_path: Path) -> None:
    with pytest.raises(DeploymentError) as captured:
        filesystem.fsync_directory(tmp_path / "missing")

    assert captured.value.issue.code == "DEPLOYMENT_STORAGE_UNAVAILABLE"


def test_fsync_directory_maps_fsync_failure_to_public_code(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(os, "fsync", _raise_oserror)

    with pytest.raises(DeploymentError) as captured:
        filesystem.fsync_directory(tmp_path)

    assert captured.value.issue.code == "DEPLOYMENT_STORAGE_UNAVAILABLE"


# layout contracts


def _layout_kwargs() -> dict[str, Path]:
    return {
        "immutable_root": Path("/opt/helixweave"),
        "configuration_root": Path("/etc/helixweave"),
        "data_root": Path("/var/lib/helixweave"),
        "log_root": Path("/var/log/helixweave"),
        "run_root": Path("/run/helixweave"),
    }


def test_layout_rejects_non_path_and_relative_roots() -> None:
    variants = (
        {**_layout_kwargs(), "immutable_root": Path("relative")},
        {**_layout_kwargs(), "immutable_root": "/opt/helixweave"},
    )
    for kwargs in variants:
        with pytest.raises(DeploymentError) as captured:
            DeploymentLayout(**kwargs)  # type: ignore[arg-type]
        assert captured.value.issue.code == "DEPLOYMENT_LAYOUT_INVALID"


def test_layout_rejects_duplicate_and_control_character_roots() -> None:
    duplicated = {key: Path("/same") for key in _layout_kwargs()}
    control = {**_layout_kwargs(), "run_root": Path("/run/helixweave\naudit")}

    for kwargs in (duplicated, control):
        with pytest.raises(DeploymentError) as captured:
            DeploymentLayout(**kwargs)
        assert captured.value.issue.code == "DEPLOYMENT_LAYOUT_INVALID"


def test_isolated_rejects_non_path_and_relative_root(tmp_path: Path) -> None:
    for rejected in (Path("relative"), "not-a-path"):
        with pytest.raises(DeploymentError) as captured:
            DeploymentLayout.isolated(rejected)  # type: ignore[arg-type]
        assert captured.value.issue.code == "DEPLOYMENT_LAYOUT_INVALID"

    layout = DeploymentLayout.isolated(tmp_path / "host")
    assert layout.data_root == tmp_path / "host" / "var" / "lib"


def test_component_store_rejects_unknown_component(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")

    with pytest.raises(DeploymentError) as captured:
        layout.component_store("unknown-component")

    assert captured.value.issue.code == "DEPLOYMENT_COMPONENT_INVALID"


def test_operator_action_and_database_prepare_properties(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    run_root = tmp_path / "host" / "run"

    assert layout.operator_action_root == run_root / "operator" / "action"
    assert (
        layout.operator_action_request
        == run_root / "operator" / "action" / "request.json"
    )
    assert (
        layout.operator_action_receipt
        == run_root / "operator" / "action" / "receipt.json"
    )
    assert layout.database_prepare_root == run_root / "database"
    assert layout.database_prepare_request == run_root / "database" / "prepare.json"
    assert (
        layout.database_prepare_receipt
        == run_root / "database" / "prepare-receipt.json"
    )


# model validation


def _manifest_dict() -> dict[str, object]:
    manifest, _payload = manifest_for(PLATFORM)
    return manifest.to_dict()


def test_file_record_rejects_non_object_and_non_string_path() -> None:
    manifest, _payload = manifest_for(PLATFORM)
    base = manifest.files[0].to_dict()

    for rejected in (["not-a-mapping"], {1: "x"}, {**base, "path": 123}):
        with pytest.raises(DeploymentError) as captured:
            FileRecord.from_dict(rejected)
        assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"


def test_file_record_rejects_invalid_size_and_mode() -> None:
    manifest, _payload = manifest_for(PLATFORM)
    base = manifest.files[0].to_dict()

    variants = (
        {"size_bytes": -1},
        {"size_bytes": True},
        {"size_bytes": "large"},
        {"size_bytes": 100 * 1024**3 + 1},
        {"mode": 0o777},
        {"mode": True},
    )
    for patch in variants:
        with pytest.raises(DeploymentError) as captured:
            FileRecord.from_dict({**base, **patch})
        assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"


def test_contract_identity_from_dict_round_trip_and_rejects_bad_keys() -> None:
    parsed = ContractIdentity.from_dict(
        {"contract": "helixweave.platform.distribution", "identity": IDENTITY_A}
    )
    assert parsed == ContractIdentity("helixweave.platform.distribution", IDENTITY_A)

    with pytest.raises(DeploymentError) as captured:
        ContractIdentity.from_dict({"contract": "helixweave.platform.distribution"})
    assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"


def test_contract_document_rejects_path_outside_contracts() -> None:
    with pytest.raises(DeploymentError) as captured:
        ContractDocument.from_dict(
            {
                "contract": "helixweave.platform.distribution",
                "identity": IDENTITY_A,
                "path": "payload/other/document.json",
            }
        )

    assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"


def test_contract_requirement_rejects_invalid_facts() -> None:
    variants = (
        {"contract": "Bad Contract", "accepted_identities": (IDENTITY_A,)},
        {"contract": "helixweave.platform.distribution", "accepted_identities": ()},
        {
            "contract": "helixweave.platform.distribution",
            "accepted_identities": (IDENTITY_B, IDENTITY_A),
        },
        {
            "contract": "helixweave.platform.distribution",
            "accepted_identities": (IDENTITY_A, IDENTITY_A),
        },
        {
            "contract": "helixweave.platform.distribution",
            "accepted_identities": ("sha256-not-hex",),
        },
    )
    for kwargs in variants:
        with pytest.raises(DeploymentError) as captured:
            ContractRequirement(**kwargs)
        assert captured.value.issue.code == "DEPLOYMENT_CONTRACT_FACTS_INVALID"


def test_contract_requirement_from_dict_rejects_bad_identities() -> None:
    base = {"contract": "helixweave.platform.distribution"}
    variants = (
        {**base, "accepted_identities": IDENTITY_A},
        {**base, "accepted_identities": []},
        {**base, "accepted_identities": [IDENTITY_A, IDENTITY_A]},
        {**base, "accepted_identities": [IDENTITY_B, IDENTITY_A]},
    )
    for document in variants:
        with pytest.raises(DeploymentError) as captured:
            ContractRequirement.from_dict(document)
        assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"


def test_manifest_rejects_wrong_schema_version_and_component() -> None:
    for patch in (
        {"schema_version": "other-schema"},
        {"component": "unknown-component"},
    ):
        document = {**_manifest_dict(), **patch}
        with pytest.raises(DeploymentError) as captured:
            BundleManifest.from_dict(document)
        assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"


def test_manifest_rejects_non_sequence_and_unordered_records() -> None:
    document = _manifest_dict()
    document["contracts"] = "not-a-sequence"
    with pytest.raises(DeploymentError) as captured:
        BundleManifest.from_dict(document)
    assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"

    reordered = _manifest_dict()
    reordered["files"] = list(reversed(reordered["files"]))
    with pytest.raises(DeploymentError) as captured:
        BundleManifest.from_dict(reordered)
    assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"


def test_manifest_rejects_tampered_identity() -> None:
    document = _manifest_dict()
    document["files"][0]["size_bytes"] += 1

    with pytest.raises(DeploymentError) as captured:
        BundleManifest.from_dict(document)

    assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"


def test_manifest_rejects_duplicate_provided_contract() -> None:
    manifest, _payload = manifest_for(PLATFORM)
    first = manifest.contracts[0]
    duplicate = ContractDocument(
        contract=first.contract,
        identity=first.identity,
        path="payload/contracts/zz-duplicate.json",
    )
    files = (
        *manifest.files,
        FileRecord(
            path=duplicate.path,
            size_bytes=1,
            sha256="0" * 64,
            mode=0o444,
        ),
    )

    with pytest.raises(DeploymentError) as captured:
        BundleManifest.create(
            component=PLATFORM,
            contracts=(*manifest.contracts, duplicate),
            files=files,
        )

    assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"


def test_component_slots_rejects_duplicate_slot_identities() -> None:
    with pytest.raises(DeploymentError) as captured:
        ComponentSlots.from_dict(
            {"active": IDENTITY_A, "previous": IDENTITY_A, "staged": None}
        )

    assert captured.value.issue.code == "DEPLOYMENT_STATE_INVALID"


def test_state_rejects_invalid_document_fields() -> None:
    base = DeploymentState.initial().to_dict()
    wrong_components = {
        key: value for key, value in base["components"].items() if key != PLATFORM
    }
    variants = (
        {"schema_version": "other-schema"},
        {"generation": -1},
        {"generation": True},
        {"created_at": 123},
        {"created_at": "x" * 65},
        {"created_at": "not-a-timestamp"},
        {"created_at": "2024-01-01T00:00:00"},
        {"components": wrong_components},
    )
    for patch in variants:
        with pytest.raises(DeploymentError) as captured:
            DeploymentState.from_dict({**base, **patch})
        assert captured.value.issue.code == "DEPLOYMENT_STATE_INVALID"


def test_state_rejects_tampered_identity() -> None:
    document = DeploymentState.initial().to_dict()
    document["generation"] = 7

    with pytest.raises(DeploymentError) as captured:
        DeploymentState.from_dict(document)

    assert captured.value.issue.code == "DEPLOYMENT_STATE_INVALID"


def test_state_transitions_require_staged_previous_and_known_component() -> None:
    state = DeploymentState.initial()

    with pytest.raises(DeploymentError) as staged_missing:
        state.activate(PLATFORM)
    assert staged_missing.value.issue.code == "DEPLOYMENT_STAGE_MISSING"

    with pytest.raises(DeploymentError) as previous_missing:
        state.rollback(PLATFORM)
    assert previous_missing.value.issue.code == "DEPLOYMENT_PREVIOUS_MISSING"

    for transition in (
        lambda: state.stage("unknown-component", IDENTITY_A),
        lambda: state.activate("unknown-component"),
        lambda: state.rollback("unknown-component"),
    ):
        with pytest.raises(DeploymentError) as captured:
            transition()
        assert captured.value.issue.code == "DEPLOYMENT_COMPONENT_INVALID"


# doctor probes


def _staged_platform_manager(tmp_path: Path):
    tmp_path.mkdir(parents=True, exist_ok=True)
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(PLATFORM)
    manager.stage(write_bundle(tmp_path / "platform.tar", manifest, payload))
    return manager, manifest


def _activated_platform_manager(tmp_path: Path):
    manager, manifest = _staged_platform_manager(tmp_path)
    manager.activate(PLATFORM, expected_staged_identity=manifest.identity)
    return manager, manifest


def test_active_frontend_probe_requires_active_platform(tmp_path: Path) -> None:
    manager, _manifest = _staged_platform_manager(tmp_path)

    result = ActiveFrontendProbe(DeploymentSnapshot(manager))()

    assert result == ProbeResult(FAIL, "RUNTIME_NOT_ACTIVE")


def test_active_frontend_probe_reports_frontend_identity(tmp_path: Path) -> None:
    manager, manifest = _activated_platform_manager(tmp_path)
    expected = next(
        item.identity
        for item in manifest.contracts
        if item.contract == "helixweave.platform.frontend-assets"
    )

    result = ActiveFrontendProbe(DeploymentSnapshot(manager))()

    assert result == ProbeResult(PASS, "FRONTEND_READY", expected)


def test_active_frontend_probe_rejects_facts_without_frontend_contract() -> None:
    snapshot = SimpleNamespace(
        read=lambda: SimpleNamespace(manifests={PLATFORM: {"active": object()}}),
        admit=lambda component: SimpleNamespace(contracts=()),
    )

    with pytest.raises(DeploymentError) as captured:
        ActiveFrontendProbe(snapshot)()  # type: ignore[arg-type]

    assert captured.value.issue.code == "DEPLOYMENT_CONTRACT_ADMISSION_FAILED"


def test_database_probe_reports_unavailable_and_ready(tmp_path: Path) -> None:
    staged_manager, _manifest = _staged_platform_manager(tmp_path / "staged")
    staged_result = DatabaseProbe(DeploymentSnapshot(staged_manager), staged_manager)()
    assert staged_result == ProbeResult(FAIL, "DATABASE_UNAVAILABLE")

    manager, _manifest = _activated_platform_manager(tmp_path / "active")
    result = DatabaseProbe(DeploymentSnapshot(manager), manager)()

    assert result.state == PASS
    assert result.reason_code == "DATABASE_READY"
    assert result.evidence_identity is not None


def test_database_probe_detects_schema_head_mismatch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    manager, _manifest = _activated_platform_manager(tmp_path)
    monkeypatch.setattr(
        manager,
        "observe_database_schema",
        lambda state: SimpleNamespace(heads=("schema-v9",), identity=IDENTITY_B),
    )

    result = DatabaseProbe(DeploymentSnapshot(manager), manager)()

    assert result == ProbeResult(FAIL, "DATABASE_SCHEMA_INCOMPATIBLE", IDENTITY_B)


def test_frontend_probe_reports_packaged_assets_ready() -> None:
    result = frontend_probe()

    assert result == ProbeResult(
        PASS,
        "FRONTEND_READY",
        load_packaged_frontend_assets().manifest.identity,
    )
