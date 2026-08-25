"""Rejection-path coverage for deployment state documents and environment."""

from __future__ import annotations

from pathlib import Path

import pytest

from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import DeploymentState, PLATFORM
from encode_pipeline.deployment.state import (
    PLATFORM_ENV_KEYS,
    STATE_TRANSACTION_SCHEMA,
    StateRecoveryPlan,
    parse_platform_environment,
    render_platform_environment,
)

_ID = f"sha256-{'a' * 64}"
_ID_B = f"sha256-{'b' * 64}"
_TASK = f"task-{'a' * 32}"
_API_SHA256 = "a" * 64


def _reject(callable_, *args, code: str, **kwargs) -> None:
    with pytest.raises(DeploymentError) as captured:
        callable_(*args, **kwargs)
    assert captured.value.issue.code == code


def _transaction_document() -> bytes:
    return canonical_json_bytes(
        {
            "schema_version": STATE_TRANSACTION_SCHEMA,
            "operation": "activate",
            "prior_state": _ID,
            "target_state": _ID_B,
        }
    )


@pytest.mark.parametrize(
    "content",
    [
        b"not json",
        b"[1]",
        canonical_json_bytes({"operation": "activate"}),
        canonical_json_bytes(
            {
                "schema_version": STATE_TRANSACTION_SCHEMA,
                "operation": "activate",
                "prior_state": None,
                "target_state": _ID,
                "extra": 1,
            }
        ),
        canonical_json_bytes(
            {
                "schema_version": "other",
                "operation": "activate",
                "prior_state": None,
                "target_state": _ID,
            }
        ),
        canonical_json_bytes(
            {
                "schema_version": STATE_TRANSACTION_SCHEMA,
                "operation": 1,
                "prior_state": None,
                "target_state": _ID,
            }
        ),
        canonical_json_bytes(
            {
                "schema_version": STATE_TRANSACTION_SCHEMA,
                "operation": "",
                "prior_state": None,
                "target_state": _ID,
            }
        ),
        canonical_json_bytes(
            {
                "schema_version": STATE_TRANSACTION_SCHEMA,
                "operation": "bad operation",
                "prior_state": None,
                "target_state": _ID,
            }
        ),
        canonical_json_bytes(
            {
                "schema_version": STATE_TRANSACTION_SCHEMA,
                "operation": "activate",
                "prior_state": None,
                "target_state": "bad",
            }
        ),
        canonical_json_bytes(
            {
                "schema_version": STATE_TRANSACTION_SCHEMA,
                "operation": "activate",
                "prior_state": "bad",
                "target_state": _ID,
            }
        ),
    ],
)
def test_recovery_plan_rejects_malformed_content(content: bytes) -> None:
    _reject(
        StateRecoveryPlan.from_bytes,
        _TASK,
        content,
        code="DEPLOYMENT_STATE_INVALID",
    )


def test_recovery_plan_rejects_noncanonical_content() -> None:
    import json

    content = json.dumps(
        {
            "schema_version": STATE_TRANSACTION_SCHEMA,
            "operation": "activate",
            "prior_state": None,
            "target_state": _ID,
        },
        indent=2,
    ).encode("utf-8")
    _reject(
        StateRecoveryPlan.from_bytes,
        _TASK,
        content,
        code="DEPLOYMENT_STATE_INVALID",
    )


def test_recovery_plan_accepts_a_canonical_transaction() -> None:
    plan = StateRecoveryPlan.from_bytes(_TASK, _transaction_document())
    assert plan.operation == "activate"
    assert plan.prior_state_identity == _ID


def _active_state() -> DeploymentState:
    return DeploymentState.initial().stage(PLATFORM, _ID).activate(PLATFORM)


@pytest.mark.parametrize(
    "api_contract_sha256",
    [None, 1, "a" * 63, "a" * 65, "z" * 64, True],
)
def test_render_platform_environment_rejects_malformed_contract(
    tmp_path: Path, api_contract_sha256
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    _reject(
        render_platform_environment,
        layout,
        _active_state(),
        api_contract_sha256=api_contract_sha256,
        code="OPERATOR_CONFIGURATION_INVALID",
    )


def test_render_platform_environment_requires_an_active_platform(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    _reject(
        render_platform_environment,
        layout,
        DeploymentState.initial(),
        api_contract_sha256=_API_SHA256,
        code="OPERATOR_CONFIGURATION_INVALID",
    )


def test_platform_environment_round_trips(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    state = _active_state()
    rendered = render_platform_environment(
        layout, state, api_contract_sha256=_API_SHA256
    )
    parsed = parse_platform_environment(layout, state, rendered.content)
    assert parsed.identity == rendered.identity
    assert parsed.api_contract_sha256 == _API_SHA256
    assert len(rendered.content.splitlines()) == len(PLATFORM_ENV_KEYS)


def _mutated_environment(tmp_path: Path, mutation) -> bytes:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    state = _active_state()
    rendered = render_platform_environment(
        layout, state, api_contract_sha256=_API_SHA256
    )
    return mutation(rendered.content)


@pytest.mark.parametrize(
    "mutation",
    [
        lambda content: content.rstrip(b"\n"),
        lambda content: content + b"EXTRA=1\n",
        lambda content: content.split(b"\n", 1)[1],
        lambda content: content.replace(b"=", b":", 1),
        lambda content: content + content.split(b"\n", 1)[0] + b"\n",
        lambda content: content.replace(b"\n", b"\r\n", 1),
        lambda content: b"\n".join(reversed(content.split(b"\n"))),
        lambda content: content.replace(b"encode-runs", b"other-queue", 1),
        lambda content: b"x" * (64 * 1024 + 1) + b"\n",
    ],
)
def test_parse_platform_environment_rejects_drift(tmp_path: Path, mutation) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    state = _active_state()
    _reject(
        parse_platform_environment,
        layout,
        state,
        _mutated_environment(tmp_path, mutation),
        code="DEPLOYMENT_CONFIGURATION_INVALID",
    )
