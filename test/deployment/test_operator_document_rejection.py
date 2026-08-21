"""Rejection-path coverage for operator service and observation evidence."""

from __future__ import annotations

import pytest

from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.models import COMPONENTS
from encode_pipeline.deployment.operator import (
    SERVICE_SOCKET_NAMES,
    SERVICE_UNITS,
    OperatorBoundaryObservation,
    OperatorObservation,
    ServiceIdentity,
    SocketWitness,
)

_ID = f"sha256-{'a' * 64}"
_ID_B = f"sha256-{'b' * 64}"
_TASK = f"task-{'a' * 32}"
_UNIT = "helixweave-api.service"
_SOCKET = SocketWitness(name="api-http", device=1, inode=2, kernel_inode=3)


def _reject(parser, document: object = None, code: str = "", **kwargs) -> None:
    with pytest.raises(DeploymentError) as captured:
        if kwargs:
            parser(**kwargs)
        else:
            parser(document)
    assert captured.value.issue.code == code


@pytest.mark.parametrize(
    "fields",
    [
        {"name": "other", "device": 1, "inode": 2, "kernel_inode": 3},
        {"name": "api-http", "device": True, "inode": 2, "kernel_inode": 3},
        {"name": "api-http", "device": 0, "inode": 2, "kernel_inode": 3},
        {"name": "api-http", "device": 1, "inode": 0, "kernel_inode": 3},
        {"name": "api-http", "device": 1, "inode": 2, "kernel_inode": 0},
    ],
)
def test_socket_witness_rejects_invalid_fields(fields) -> None:
    _reject(
        SocketWitness,
        code="OPERATOR_SERVICE_IDENTITY_INVALID",
        **fields,
    )


def _service_identity_document(unit: str = _UNIT) -> dict[str, object]:
    return ServiceIdentity.create(
        unit=unit,
        deployment_identity=_ID,
        task_identity=_TASK,
        main_pid=100,
        process_start_ticks=200,
        executable_device=10,
        executable_inode=20,
        cmdline_identity=_ID,
        boot_identity=_ID,
        invocation_identity=_ID,
        cgroup_identity=_ID,
        sockets=[
            SocketWitness(name=name, device=1, inode=2, kernel_inode=3)
            for name in SERVICE_SOCKET_NAMES[unit]
        ],
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: ["unit"],
        lambda d: {key: d[key] for key in d if key != "sockets"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "unit": "other.service"},
        lambda d: {**d, "sockets": "api-http"},
        lambda d: {**d, "sockets": ["api-http"]},
        lambda d: {**d, "sockets": [{"name": "api-http", "device": 1}]},
        lambda d: {**d, "sockets": []},
        lambda d: {
            **d,
            "sockets": [
                {"name": "api-http", "device": 1, "inode": 2, "kernel_inode": 3},
                {"name": "api-http", "device": 4, "inode": 5, "kernel_inode": 6},
            ],
        },
        lambda d: {**d, "deployment_identity": "bad"},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "cmdline_identity": "bad"},
        lambda d: {**d, "boot_identity": 1},
        lambda d: {**d, "invocation_identity": "bad"},
        lambda d: {**d, "cgroup_identity": "bad"},
        lambda d: {**d, "main_pid": True},
        lambda d: {**d, "main_pid": 0},
        lambda d: {**d, "process_start_ticks": 0},
        lambda d: {**d, "executable_device": 0},
        lambda d: {**d, "executable_inode": 0},
        lambda d: {**d, "identity": _ID_B},
    ],
)
def test_service_identity_rejects_malformed_documents(mutation) -> None:
    _reject(
        ServiceIdentity.from_dict,
        mutation(_service_identity_document()),
        "OPERATOR_SERVICE_IDENTITY_INVALID",
    )


def test_service_identity_round_trips_socketless_unit() -> None:
    document = _service_identity_document("helixweave-worker.service")
    identity = ServiceIdentity.from_dict(document)
    assert identity.unit == "helixweave-worker.service"
    assert identity.sockets == ()


@pytest.mark.parametrize(
    "document",
    [
        None,
        ["ready"],
        {"status": "ready", "reason_code": "READY"},
        {"status": "ready", "reason_code": "READY", "identity": _ID, "extra": 1},
        {"status": "unknown", "reason_code": "READY", "identity": _ID},
        {"status": "ready", "reason_code": "BOUNDARY_INVALID", "identity": _ID},
        {"status": "ready", "reason_code": "READY", "identity": None},
        {"status": "not-ready", "reason_code": "READY", "identity": None},
        {"status": "not-ready", "reason_code": "BOUNDARY_INVALID", "identity": _ID},
    ],
)
def test_boundary_observation_rejects_malformed_documents(document) -> None:
    _reject(
        OperatorBoundaryObservation.from_dict,
        document,
        "OPERATOR_OBSERVATION_INVALID",
    )


def _observation_document() -> dict[str, object]:
    return OperatorObservation.create(
        state_identity=_ID,
        active={component: _ID for component in COMPONENTS},
        database_schema_identity=_ID_B,
        database_schema_heads=("20260809_13",),
        services={unit: _ID for unit in SERVICE_UNITS},
        operator_boundary=OperatorBoundaryObservation.ready(_ID),
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: ["services"],
        lambda d: {key: d[key] for key in d if key != "active"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "schema_version": "other"},
        lambda d: {**d, "active": []},
        lambda d: {**d, "active": {"platform": _ID}},
        lambda d: {**d, "active": {**d["active"], "platform": "bad"}},
        lambda d: {**d, "services": []},
        lambda d: {
            **d,
            "services": {
                unit: identity
                for unit, identity in d["services"].items()
                if unit != "helixweave-redis.service"
            },
        },
        lambda d: {**d, "services": {**d["services"], "helixweave-api.service": "bad"}},
        lambda d: {**d, "database_schema_heads": "20260809_13"},
        lambda d: {**d, "database_schema_heads": ["b_head", "a_head"]},
        lambda d: {**d, "database_schema_heads": ["a_head", "a_head"]},
        lambda d: {**d, "database_schema_heads": ["!!bad"]},
        lambda d: {**d, "state_identity": "bad"},
        lambda d: {**d, "database_schema_identity": "bad"},
        lambda d: {**d, "database_schema_identity": None},
        lambda d: {**d, "identity": _ID_B},
        lambda d: {**d, "operator_pending_count": True},
        lambda d: {**d, "operator_pending_count": 2},
        lambda d: {**d, "operator_recovery_required_count": True},
        lambda d: {**d, "operator_recovery_required_count": 2},
        lambda d: {
            **d,
            "operator_pending_count": 1,
            "operator_recovery_required_count": 1,
        },
        lambda d: {**d, "operator_boundary": None},
        lambda d: {
            **d,
            "operator_boundary": {
                "status": "ready",
                "reason_code": "READY",
                "identity": None,
            },
        },
    ],
)
def test_operator_observation_rejects_malformed_documents(mutation) -> None:
    _reject(
        OperatorObservation.from_dict,
        mutation(_observation_document()),
        "OPERATOR_OBSERVATION_INVALID",
    )


def test_operator_observation_defaults_to_invalid_boundary() -> None:
    observation = OperatorObservation.create(
        state_identity=_ID,
        active={component: None for component in COMPONENTS},
        database_schema_identity=None,
        database_schema_heads=(),
        services={unit: None for unit in SERVICE_UNITS},
    )
    assert observation.operator_boundary.status == "not-ready"
