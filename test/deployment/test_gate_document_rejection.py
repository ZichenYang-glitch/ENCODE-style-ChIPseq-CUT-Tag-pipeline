"""Rejection-path coverage for Deployment Gate evidence documents."""

from __future__ import annotations

import pytest

from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.gate import (
    REQUIRED_CACHE_KINDS,
    RESOURCE_KINDS,
    RUNTIME_COMPONENTS,
    CacheEvidence,
    GateProcessIdentity,
    ResourceEvidence,
    RuntimeIdentity,
)

_ID = f"sha256-{'a' * 64}"
_ID_B = f"sha256-{'b' * 64}"
_TASK = f"task-{'a' * 32}"
_RESOURCE_KIND = RESOURCE_KINDS[0]
_CACHE_KIND = REQUIRED_CACHE_KINDS[0]
_RUNTIME_COMPONENT = RUNTIME_COMPONENTS[0]


def _reject(parser, document: object, code: str, forward: bool = False) -> None:
    with pytest.raises(DeploymentError) as captured:
        if forward:
            parser(document, code=code)
        else:
            parser(document)
    assert captured.value.issue.code == code


@pytest.mark.parametrize(
    "document",
    [
        None,
        ["component"],
        {"component": _RUNTIME_COMPONENT},
        {"component": _RUNTIME_COMPONENT, "identity": _ID, "extra": 1},
        {"component": "platform", "identity": _ID},
        {"component": _RUNTIME_COMPONENT, "identity": "bad"},
    ],
)
def test_runtime_identity_rejects_malformed_documents(document) -> None:
    _reject(
        RuntimeIdentity.from_dict,
        document,
        code="GATE_OBSERVATION_INVALID",
        forward=True,
    )


def test_runtime_identity_round_trip() -> None:
    document = {"component": _RUNTIME_COMPONENT, "identity": _ID}
    identity = RuntimeIdentity.from_dict(document, code="GATE_OBSERVATION_INVALID")
    assert identity.to_dict() == document


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "kind"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "kind": "other"},
        lambda d: {**d, "path_identity": "bad"},
        lambda d: {**d, "device": True},
        lambda d: {**d, "device": 0},
        lambda d: {**d, "inode": 0},
    ],
)
def test_resource_evidence_rejects_malformed_documents(mutation) -> None:
    document = mutation(
        ResourceEvidence(
            kind=_RESOURCE_KIND, path_identity=_ID, device=10, inode=20
        ).to_dict()
    )
    _reject(ResourceEvidence.from_dict, document, "GATE_OBSERVATION_INVALID")


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: {key: d[key] for key in d if key != "identity"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "kind": "other"},
        lambda d: {**d, "identity": "bad"},
        lambda d: {**d, "device": True},
        lambda d: {**d, "inode": -1},
    ],
)
def test_cache_evidence_rejects_malformed_documents(mutation) -> None:
    document = mutation(
        CacheEvidence(kind=_CACHE_KIND, identity=_ID, device=10, inode=20).to_dict()
    )
    _reject(CacheEvidence.from_dict, document, "GATE_OBSERVATION_INVALID")


def _gate_process(name: str = "runner") -> dict[str, object]:
    sockets: dict[str, object] = (
        {"socket_device": None, "socket_inode": None, "socket_kernel_inode": None}
        if name == "runner"
        else {"socket_device": 1, "socket_inode": 2, "socket_kernel_inode": 3}
    )
    return GateProcessIdentity.create(
        name=name,
        task_identity=_TASK,
        pid=100,
        process_start_ticks=200,
        executable_device=10,
        executable_inode=20,
        cmdline_identity=_ID,
        **sockets,
    ).to_dict()


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: None,
        lambda d: ["runner"],
        lambda d: {key: d[key] for key in d if key != "pid"},
        lambda d: {**d, "extra": 1},
        lambda d: {**d, "name": "other"},
        lambda d: {**d, "identity": _ID_B},
        lambda d: {**d, "task_identity": "bad"},
        lambda d: {**d, "pid": True},
        lambda d: {**d, "pid": 0},
        lambda d: {**d, "pid": 2**31},
        lambda d: {**d, "process_start_ticks": 0},
        lambda d: {**d, "executable_device": 0},
        lambda d: {**d, "executable_inode": True},
        lambda d: {**d, "cmdline_identity": "bad"},
        lambda d: {**d, "socket_device": 1},
        lambda d: {**d, "socket_inode": 2},
        lambda d: {**d, "socket_kernel_inode": 3},
    ],
)
def test_gate_process_identity_rejects_malformed_runner(mutation) -> None:
    document = mutation(_gate_process("runner"))
    _reject(GateProcessIdentity.from_dict, document, "GATE_PROCESS_IDENTITY_INVALID")


@pytest.mark.parametrize(
    "mutation",
    [
        lambda d: {**d, "socket_device": None},
        lambda d: {**d, "socket_inode": True},
        lambda d: {**d, "socket_kernel_inode": 0},
        lambda d: {**d, "socket_device": 2**63},
    ],
)
def test_gate_process_identity_rejects_malformed_service_sockets(mutation) -> None:
    document = mutation(_gate_process("redis"))
    _reject(GateProcessIdentity.from_dict, document, "GATE_PROCESS_IDENTITY_INVALID")


def test_gate_process_identity_round_trip() -> None:
    document = _gate_process("redis")
    identity = GateProcessIdentity.from_dict(document)
    assert identity.to_dict() == document
