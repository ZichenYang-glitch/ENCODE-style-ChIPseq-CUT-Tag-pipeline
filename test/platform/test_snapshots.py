from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.snapshots import (
    PAYLOAD_DIGEST_SCHEME,
    VALIDATION_EVIDENCE_OUTCOME,
    ValidatedInputSnapshot,
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)


NOW = datetime(2026, 7, 14, 8, 0, tzinfo=timezone.utc)


def _identity() -> WorkflowBuildIdentity:
    return WorkflowBuildIdentity(
        workflow_id="workflow-a",
        adapter_version="1.2.3",
        scheme="sha256-tree-v1",
        logical_entrypoint="workflow/Snakefile",
        digest="a" * 64,
        captured_at=NOW,
    )


def _snapshot(**changes: object) -> ValidatedInputSnapshot:
    inputs = WorkflowInputs(
        config={"threads": 2, "nested": {"enabled": True}},
        samples=[{"sample": "S1", "assay": "chipseq"}],
        options={"dry": False},
    )
    payload = canonical_workflow_inputs_json(inputs)
    values: dict[str, object] = {
        "snapshot_id": "vsnap_0123456789abcdef0123456789abcdef",
        "workflow_id": "workflow-a",
        "adapter_version": "1.2.3",
        "schema_version": "1.0.0",
        "schema_dialect": "https://json-schema.org/draft/2020-12/schema",
        "workflow_build_identity": _identity(),
        "canonical_payload": payload,
        "payload_digest_scheme": PAYLOAD_DIGEST_SCHEME,
        "payload_digest": build_workflow_inputs_digest(payload),
        "validation_outcome": VALIDATION_EVIDENCE_OUTCOME,
        "validation_issue_codes": ("INPUT_WARNING",),
        "validated_at": NOW,
        "expires_at": NOW + timedelta(minutes=30),
        "consumed_run_id": None,
        "consumed_at": None,
    }
    values.update(changes)
    return ValidatedInputSnapshot(**values)


def test_canonical_workflow_inputs_are_order_independent_and_framed() -> None:
    first = WorkflowInputs(
        config={"z": 1, "a": {"right": 2, "left": 1}},
        samples=[{"sample": "S1"}],
        options={},
    )
    second = WorkflowInputs(
        config={"a": {"left": 1, "right": 2}, "z": 1},
        samples=[{"sample": "S1"}],
        options={},
    )

    first_payload = canonical_workflow_inputs_json(first)
    second_payload = canonical_workflow_inputs_json(second)

    assert first_payload == second_payload
    assert build_workflow_inputs_digest(first_payload) == (
        build_workflow_inputs_digest(second_payload)
    )
    assert len(build_workflow_inputs_digest(first_payload)) == 64


@pytest.mark.parametrize("unsafe", [2**53, -(2**53), float("nan"), float("inf")])
def test_canonical_workflow_inputs_reject_unsafe_numbers(unsafe: object) -> None:
    with pytest.raises(ValueError, match="JSON-safe"):
        canonical_workflow_inputs_json(
            WorkflowInputs(config={"unsafe": unsafe}, samples=None, options={})
        )


def test_snapshot_round_trips_fresh_workflow_inputs() -> None:
    snapshot = _snapshot()

    first = snapshot.to_workflow_inputs()
    first.config["threads"] = 99
    second = snapshot.to_workflow_inputs()

    assert second.config["threads"] == 2
    assert second.samples == [{"assay": "chipseq", "sample": "S1"}]


@pytest.mark.parametrize(
    "changes",
    [
        {"snapshot_id": "not-opaque"},
        {"payload_digest": "b" * 64},
        {"canonical_payload": '{"config":{},"options":{},"samples":[] }'},
        {"workflow_id": "other"},
        {"adapter_version": "other"},
        {"expires_at": NOW},
        {"validation_outcome": "failed"},
        {"validation_issue_codes": ("private/path",)},
        {"consumed_run_id": "run-1", "consumed_at": None},
    ],
)
def test_snapshot_rejects_corrupt_or_inconsistent_evidence(changes) -> None:
    with pytest.raises(ValueError):
        _snapshot(**changes)


def test_snapshot_consumption_is_pairwise_and_immutable() -> None:
    snapshot = _snapshot()
    consumed = snapshot.with_consumption("run-1", NOW + timedelta(minutes=1))

    assert snapshot.consumed_run_id is None
    assert consumed.consumed_run_id == "run-1"
    assert consumed.consumed_at == NOW + timedelta(minutes=1)
    with pytest.raises(ValueError, match="already consumed"):
        consumed.with_consumption("run-2", NOW + timedelta(minutes=2))


def test_replacing_canonical_payload_without_digest_is_detected() -> None:
    snapshot = _snapshot()
    tampered = snapshot.canonical_payload.replace('"threads":2', '"threads":3')

    with pytest.raises(ValueError, match="digest"):
        replace(snapshot, canonical_payload=tampered)
