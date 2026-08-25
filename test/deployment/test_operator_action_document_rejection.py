"""Rejection-path coverage for operator action evidence document parsers."""

from __future__ import annotations

import pytest

from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.models import COMPONENTS
from encode_pipeline.deployment.operator_action import (
    VERIFICATION_CHECKS,
    BulkRuntimePrepareReceipt,
    BulkRuntimePrepareRequest,
    DatabasePrepareReceipt,
    DatabasePrepareRequest,
    DeploymentActionReceipt,
    DeploymentActionRequest,
    EncodeRuntimeEntry,
    EncodeRuntimeInventory,
    EncodeRuntimePrepareReceipt,
    EncodeRuntimePrepareRequest,
    ReadinessCheck,
)

_ID = f"sha256-{'a' * 64}"
_ID_B = f"sha256-{'b' * 64}"
_ID_C = f"sha256-{'c' * 64}"
_ID_D = f"sha256-{'d' * 64}"
_ID_E = f"sha256-{'e' * 64}"
_ID_F = f"sha256-{'f' * 64}"
_TASK = f"task-{'a' * 32}"
_SHA256 = "a" * 64
_HEAD = "20260809_13"
_KNOWN_HEAD = "20260807_12"
_NATIVE_IDENTITIES = {component: _ID for component in COMPONENTS}
_READINESS = {
    name: ReadinessCheck(status="ready", reason_code="READY", identity=_ID)
    for name in VERIFICATION_CHECKS
}


def _reject(parser, document: object, code: str) -> None:
    with pytest.raises(DeploymentError) as captured:
        parser(document)
    assert captured.value.issue.code == code


def _action_request() -> dict[str, object]:
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
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "phase"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "phase": "execute"},
        lambda d: {**d, "operation": "execute"},
        lambda d: {**d, "component": "other"},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "deployment_identity": "bad"},
        lambda d: {**d, "authority_platform_identity": 1},
        lambda d: {**d, "prior_state_identity": "bad"},
        lambda d: {**d, "candidate_state_identity": "bad"},
        lambda d: {**d, "candidate_active": []},
        lambda d: {**d, "identity": _ID_F},
    ],
)
def test_action_request_rejects_malformed_documents(mutation) -> None:
    _reject(
        DeploymentActionRequest.from_dict,
        mutation(_action_request()),
        "OPERATOR_ACTION_REQUEST_INVALID",
    )


@pytest.mark.parametrize(
    "document",
    [
        None,
        ["ready"],
        {"status": "ready", "reason_code": "READY"},
        {"status": "ready", "reason_code": "READY", "identity": _ID, "extra": 1},
        {"status": "unknown", "reason_code": "READY", "identity": _ID},
        {"status": "ready", "reason_code": "OPERATION_CONFLICT", "identity": _ID},
        {"status": "ready", "reason_code": "READY", "identity": None},
        {"status": "not-ready", "reason_code": "READY", "identity": _ID},
        {"status": "not-ready", "reason_code": "CONTRACT_INVALID", "identity": "bad"},
    ],
)
def test_readiness_check_rejects_malformed_documents(document) -> None:
    _reject(ReadinessCheck.from_dict, document, "OPERATOR_ACTION_RECEIPT_INVALID")


def _action_receipt() -> dict[str, object]:
    return DeploymentActionReceipt.create(
        request_identity=_ID,
        status="observed",
        compatibility="compatible",
        database_before_identity=None,
        database_after_identity=None,
        accepted_schema_heads=(_HEAD,),
        target_schema_heads=(_HEAD,),
        migration_inventory_identity=_ID_B,
        known_schema_revisions=(_KNOWN_HEAD,),
        migration_required=False,
        rollback_supported=True,
        api_contract_sha256=_SHA256,
        native_identities=dict(_NATIVE_IDENTITIES),
        frontend_identity=_ID_C,
        reference_compatibility_identity=_ID_D,
        readiness=dict(_READINESS),
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: ["status"],
        lambda d: {key: d[key] for key in d if key != "readiness"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "status": "executed"},
        lambda d: {**d, "compatibility": "unknown"},
        lambda d: {**d, "migration_required": "no"},
        lambda d: {**d, "rollback_supported": "yes"},
        lambda d: {**d, "readiness": []},
        lambda d: {**d, "readiness": {**d["readiness"], "extra": {}}},
        lambda d: {
            **d,
            "readiness": {
                name: value
                for name, value in d["readiness"].items()
                if name != "docker"
            },
        },
        lambda d: {**d, "native_identities": []},
        lambda d: {**d, "native_identities": {"platform": _ID}},
        lambda d: {
            **d,
            "native_identities": {**d["native_identities"], "platform": "bad"},
        },
        lambda d: {**d, "api_contract_sha256": "z" * 64},
        lambda d: {**d, "api_contract_sha256": 1},
        lambda d: {**d, "identity": _ID_F},
        lambda d: {**d, "request_identity": "bad"},
        lambda d: {**d, "migration_inventory_identity": "bad"},
        lambda d: {**d, "accepted_schema_heads": [_HEAD, _HEAD]},
        lambda d: {**d, "target_schema_heads": [_HEAD, _KNOWN_HEAD]},
        lambda d: {**d, "target_schema_heads": [_KNOWN_HEAD]},
        lambda d: {**d, "accepted_schema_heads": ["unknown_head"]},
        lambda d: {**d, "migration_required": True},
        lambda d: {**d, "rollback_supported": False},
        lambda d: {**d, "known_schema_revisions": [_HEAD]},
        lambda d: {**d, "database_before_identity": "bad"},
        lambda d: {**d, "database_after_identity": 1},
        lambda d: {**d, "frontend_identity": "bad"},
        lambda d: {**d, "reference_compatibility_identity": 1},
        lambda d: {
            **d,
            "readiness": {
                **d["readiness"],
                "docker": {
                    "status": "unknown",
                    "reason_code": "READY",
                    "identity": None,
                },
            },
        },
    ],
)
def test_action_receipt_rejects_malformed_documents(mutation) -> None:
    _reject(
        DeploymentActionReceipt.from_dict,
        mutation(_action_receipt()),
        "OPERATOR_ACTION_RECEIPT_INVALID",
    )


def test_action_receipt_accepts_migration_required_shape() -> None:
    receipt = DeploymentActionReceipt.create(
        request_identity=_ID,
        status="admitted",
        compatibility="incomplete",
        database_before_identity=_ID,
        database_after_identity=_ID_B,
        accepted_schema_heads=(_KNOWN_HEAD,),
        target_schema_heads=(_HEAD,),
        migration_inventory_identity=_ID_B,
        known_schema_revisions=(_KNOWN_HEAD,),
        migration_required=True,
        rollback_supported=False,
        api_contract_sha256=_SHA256,
        native_identities={component: None for component in COMPONENTS},
        frontend_identity=None,
        reference_compatibility_identity=None,
        readiness={
            name: ReadinessCheck(
                status="unavailable", reason_code="SERVICE_STOPPED", identity=None
            )
            for name in VERIFICATION_CHECKS
        },
    )
    assert receipt.migration_required is True
    assert receipt.database_after_identity == _ID_B


def _database_prepare_request() -> dict[str, object]:
    return DatabasePrepareRequest.create(
        operation="activate",
        database_mode="existing-live",
        task_identity=_TASK,
        deployment_identity=_ID,
        prior_state_identity=_ID_B,
        candidate_state_identity=_ID_C,
        action_receipt_identity=_ID_D,
        backup_receipt_identity=_ID_E,
        target_schema_heads=(_HEAD,),
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "database_mode"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "operation": "materialize"},
        lambda d: {**d, "database_mode": "other"},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "deployment_identity": "bad"},
        lambda d: {**d, "prior_state_identity": 1},
        lambda d: {**d, "candidate_state_identity": "bad"},
        lambda d: {**d, "action_receipt_identity": "bad"},
        lambda d: {**d, "backup_receipt_identity": 1},
        lambda d: {**d, "target_schema_heads": "20260809_13"},
        lambda d: {**d, "target_schema_heads": ["!!bad"]},
        lambda d: {**d, "identity": _ID_F},
    ],
)
def test_database_prepare_request_rejects_malformed_documents(mutation) -> None:
    _reject(
        DatabasePrepareRequest.from_dict,
        mutation(_database_prepare_request()),
        "DATABASE_PREPARE_REQUEST_INVALID",
    )


def _database_prepare_receipt() -> dict[str, object]:
    return DatabasePrepareReceipt.create(
        request_identity=_ID,
        database_before_identity=_ID_B,
        database_after_identity=_ID_C,
        schema_heads=(_HEAD,),
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "schema_heads"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "status": "other"},
        lambda d: {**d, "request_identity": "bad"},
        lambda d: {**d, "database_before_identity": 1},
        lambda d: {**d, "database_after_identity": "bad"},
        lambda d: {**d, "schema_heads": ["!!bad"]},
        lambda d: {**d, "identity": _ID_F},
    ],
)
def test_database_prepare_receipt_rejects_malformed_documents(mutation) -> None:
    _reject(
        DatabasePrepareReceipt.from_dict,
        mutation(_database_prepare_receipt()),
        "DATABASE_PREPARE_RECEIPT_INVALID",
    )


def _encode_file_entry() -> dict[str, object]:
    return EncodeRuntimeEntry(
        path="envs/encode.yaml",
        kind="file",
        sha256="a" * 64,
        size=12,
        mode=0o444,
        target=None,
    ).to_dict()


def _encode_symlink_entry() -> dict[str, object]:
    target = "lib/encode.so"
    import hashlib

    return EncodeRuntimeEntry(
        path="bin/encode",
        kind="symlink",
        sha256=hashlib.sha256(target.encode("utf-8")).hexdigest(),
        size=len(target.encode("utf-8")),
        mode=None,
        target=target,
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: ["path"],
        lambda d: {key: d[key] for key in d if key != "kind"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "path": 1},
        lambda d: {**d, "path": "../escape"},
        lambda d: {**d, "path": "/absolute"},
        lambda d: {**d, "path": ".helixweave-runtime-inventory.json"},
        lambda d: {**d, "kind": "directory"},
        lambda d: {**d, "sha256": 1},
        lambda d: {**d, "sha256": "z" * 64},
        lambda d: {**d, "size": True},
        lambda d: {**d, "size": -1},
        lambda d: {**d, "size": 2**63},
        lambda d: {**d, "mode": 0o644},
        lambda d: {**d, "mode": 0o555, "target": "elsewhere"},
    ],
)
def test_encode_runtime_file_entry_rejects_malformed_documents(mutation) -> None:
    _reject(
        EncodeRuntimeEntry.from_dict,
        mutation(_encode_file_entry()),
        "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
    )


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: {**d, "mode": 0o444},
        lambda d: {**d, "target": None},
        lambda d: {**d, "target": 1},
        lambda d: {**d, "target": "/absolute"},
        lambda d: {**d, "target": "../escape"},
        lambda d: {**d, "size": 1},
        lambda d: {**d, "sha256": "a" * 64},
    ],
)
def test_encode_runtime_symlink_entry_rejects_malformed_documents(mutation) -> None:
    _reject(
        EncodeRuntimeEntry.from_dict,
        mutation(_encode_symlink_entry()),
        "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
    )


def _encode_inventory() -> dict[str, object]:
    entries = (
        EncodeRuntimeEntry(
            path="conda-envs/encode.yaml",
            kind="file",
            sha256="a" * 64,
            size=12,
            mode=0o444,
            target=None,
        ),
        EncodeRuntimeEntry(
            path="runner/bin/conda",
            kind="file",
            sha256="b" * 64,
            size=24,
            mode=0o555,
            target=None,
        ),
        EncodeRuntimeEntry(
            path="runner/bin/snakemake",
            kind="file",
            sha256="c" * 64,
            size=36,
            mode=0o555,
            target=None,
        ),
    )
    return EncodeRuntimeInventory.create(entries).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "entries"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "tree_identity": _ID},
        lambda d: {**d, "entries": "envs/encode.yaml"},
        lambda d: {**d, "entries": [{**d["entries"][0], "kind": "other"}]},
        lambda d: {**d, "entries": []},
    ],
)
def test_encode_runtime_inventory_rejects_malformed_documents(mutation) -> None:
    _reject(
        EncodeRuntimeInventory.from_dict,
        mutation(_encode_inventory()),
        "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
    )


def _encode_prepare_request() -> dict[str, object]:
    return EncodeRuntimePrepareRequest.create(
        task_identity=_TASK,
        deployment_identity=_ID,
        authority_platform_identity=_ID_B,
        prior_state_identity=_ID_C,
        candidate_state_identity=_ID_D,
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "task_identity"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "operation": "observe"},
        lambda d: {**d, "component": "platform"},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "deployment_identity": "bad"},
        lambda d: {**d, "authority_platform_identity": 1},
        lambda d: {**d, "prior_state_identity": "bad"},
        lambda d: {**d, "candidate_state_identity": "bad"},
        lambda d: {**d, "identity": _ID_F},
    ],
)
def test_encode_prepare_request_rejects_malformed_documents(mutation) -> None:
    _reject(
        EncodeRuntimePrepareRequest.from_dict,
        mutation(_encode_prepare_request()),
        "ENCODE_RUNTIME_PREPARE_REQUEST_INVALID",
    )


def _encode_prepare_receipt() -> dict[str, object]:
    return EncodeRuntimePrepareReceipt.create(
        request_identity=_ID,
        deployment_identity=_ID_B,
        inventory=EncodeRuntimeInventory.from_dict(_encode_inventory()),
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "tree_identity"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "status": "other"},
        lambda d: {**d, "request_identity": "bad"},
        lambda d: {**d, "deployment_identity": 1},
        lambda d: {**d, "tree_identity": "bad"},
        lambda d: {**d, "inventory_sha256": "z" * 64},
        lambda d: {**d, "inventory_sha256": 1},
        lambda d: {**d, "inventory_size": True},
        lambda d: {**d, "inventory_size": 0},
        lambda d: {**d, "entry_count": -1},
        lambda d: {**d, "entry_count": 5},
        lambda d: {**d, "identity": _ID_F},
    ],
)
def test_encode_prepare_receipt_rejects_malformed_documents(mutation) -> None:
    _reject(
        EncodeRuntimePrepareReceipt.from_dict,
        mutation(_encode_prepare_receipt()),
        "ENCODE_RUNTIME_PREPARE_RECEIPT_INVALID",
    )


def _bulk_prepare_request() -> dict[str, object]:
    return BulkRuntimePrepareRequest.create(
        operation="verify",
        task_identity=_TASK,
        candidate_bulk_identity=_ID,
        authority_platform_identity=_ID_B,
        prior_state_identity=_ID_C,
        candidate_state_identity=_ID_D,
        docker_service_identity=_ID_E,
        docker_client_identity=_ID_F,
        docker_endpoint_identity=_ID,
        docker_daemon_uid=999,
        docker_daemon_gid=999,
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "docker_daemon_uid"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "operation": "materialize"},
        lambda d: {**d, "component": "platform"},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "candidate_bulk_identity": "bad"},
        lambda d: {**d, "authority_platform_identity": 1},
        lambda d: {**d, "prior_state_identity": "bad"},
        lambda d: {**d, "candidate_state_identity": "bad"},
        lambda d: {**d, "docker_service_identity": "bad"},
        lambda d: {**d, "docker_client_identity": "bad"},
        lambda d: {**d, "docker_endpoint_identity": "bad"},
        lambda d: {**d, "docker_daemon_uid": True},
        lambda d: {**d, "docker_daemon_uid": -1},
        lambda d: {**d, "docker_daemon_gid": "999"},
        lambda d: {**d, "identity": _ID_F},
    ],
)
def test_bulk_prepare_request_rejects_malformed_documents(mutation) -> None:
    _reject(
        BulkRuntimePrepareRequest.from_dict,
        mutation(_bulk_prepare_request()),
        "BULK_RUNTIME_PREPARE_REQUEST_INVALID",
    )


def _bulk_prepare_receipt() -> dict[str, object]:
    return BulkRuntimePrepareReceipt.create(
        request_identity=_ID,
        candidate_bulk_identity=_ID_B,
        runtime_identity=_ID_C,
        image_set_identity=_ID_D,
        image_count=3,
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "image_count"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "status": "other"},
        lambda d: {**d, "request_identity": "bad"},
        lambda d: {**d, "candidate_bulk_identity": "bad"},
        lambda d: {**d, "runtime_identity": 1},
        lambda d: {**d, "image_set_identity": "bad"},
        lambda d: {**d, "image_count": True},
        lambda d: {**d, "image_count": 0},
        lambda d: {**d, "image_count": -1},
        lambda d: {**d, "identity": _ID_F},
    ],
)
def test_bulk_prepare_receipt_rejects_malformed_documents(mutation) -> None:
    _reject(
        BulkRuntimePrepareReceipt.from_dict,
        mutation(_bulk_prepare_receipt()),
        "BULK_RUNTIME_PREPARE_RECEIPT_INVALID",
    )
