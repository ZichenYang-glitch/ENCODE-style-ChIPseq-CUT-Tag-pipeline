"""The IRENE-31 vendor outbound spec stays offline and matches local redaction."""

from __future__ import annotations

import json
from pathlib import Path

from encode_pipeline.services.agent_redaction import RedactionPolicy


ROOT = Path(__file__).resolve().parents[2]
CASES_PATH = ROOT / "docs" / "architecture" / "peri-p2-vendor-outbound-cases.json"
SPEC_PATH = ROOT / "docs" / "architecture" / "peri-p2-vendor-outbound.md"

OFFLINE_COMMIT = "93a7631d7387438898051949c15db2d90836191a"


def _cases() -> dict[str, object]:
    loaded = json.loads(CASES_PATH.read_text(encoding="utf-8"))
    assert isinstance(loaded, dict)
    return loaded


def test_owner_decisions_and_cost_cap_stay_unset() -> None:
    document = _cases()
    assert document["contract_version"] == "p2-transport/1"
    assert document["offline_transport_commit"] == OFFLINE_COMMIT
    decisions = document["owner_decisions"]
    assert decisions == {
        "vendor": None,
        "outbound_fields": None,
        "budget": None,
        "credential_boundary": None,
    }
    assert document["cost_cap"] == {
        "status": "owner_unset",
        "currency": None,
        "max_amount_minor": None,
        "max_requests": None,
        "max_input_tokens": None,
        "max_output_tokens": None,
    }


def test_stub_cases_cover_success_rejection_and_timeout() -> None:
    document = _cases()
    cases = document["cases"]
    assert [item["id"] for item in cases] == ["success", "rejection", "timeout"]
    by_id = {item["id"]: item for item in cases}

    success = by_id["success"]["expect"]
    assert success["deliver"] == "delivery"
    assert success["source"] == "test_transport"
    assert success["suggestion_fields"] == [
        {"token": "syn.alpha", "value": "from_response"}
    ]
    assert success["projection_fields"] == [
        {"token": "syn.alpha", "value": "from_request"}
    ]
    assert success["usage"] == {
        "attempted": 1,
        "responder_calls": 1,
        "delivered": 1,
        "adopted": 1,
    }
    assert success["provenance"]["vendor"] is None
    assert success["provenance"]["snapshot_id"] is None
    assert success["provenance"]["cost_cap_status"] == "owner_unset"
    assert success["provenance"]["grant_period"] == "grant_aaaaaaaa"
    assert success["provenance"]["delivered_at"] == 1.0
    assert success["provenance"]["adopt_deadline"] == 2.0

    rejection = by_id["rejection"]["expect"]
    assert rejection["deliver"] == "refusal"
    assert rejection["reason"] == "field_not_allowed"
    assert rejection["usage"]["responder_calls"] == 0
    refusal = rejection["reason"]
    for marker in rejection["absent_from_refusal"]:
        assert marker not in refusal

    timeout = by_id["timeout"]["expect"]
    assert timeout["reason"] == "timeout"
    assert timeout["usage"] == {
        "attempted": 1,
        "responder_calls": 1,
        "delivered": 0,
        "adopted": 0,
    }
    for marker in timeout["absent_from_refusal"]:
        assert marker not in timeout["reason"]

    projection = document["projection_example"]
    outbound = json.dumps(projection["outbound_fields"], ensure_ascii=False)
    for marker in projection["absent_from_outbound"]:
        assert marker not in outbound


def test_redaction_examples_match_the_local_policy() -> None:
    policy = RedactionPolicy()
    for example in _cases()["redaction_examples"]:
        if example["kind"] == "text":
            result = policy.redact_text(example["input"])
        else:
            result = policy.redact_value(example["input"])
        assert result.value == example["output"]
        assert result.absolute_paths_redacted == example["absolute_paths_redacted"]
        assert result.sensitive_keys_redacted == example["sensitive_keys_redacted"]


def test_spec_points_at_the_offline_contract_and_the_cases() -> None:
    text = SPEC_PATH.read_text(encoding="utf-8")
    for fragment in (
        OFFLINE_COMMIT,
        "p2-transport/1",
        "vendor_unselected",
        "outbound_fields_unset",
        "budget_unset",
        "credential_boundary_unset",
        "peri-p2-vendor-outbound-cases.json",
        "open_vendor_model_transport",
        "SyntheticSuggestionTransport",
    ):
        assert fragment in text
