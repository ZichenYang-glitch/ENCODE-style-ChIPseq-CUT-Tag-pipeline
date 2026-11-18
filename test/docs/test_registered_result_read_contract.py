from __future__ import annotations

import json
from pathlib import Path
import re
from typing import Any


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
EXAMPLES_PATH = (
    REPOSITORY_ROOT
    / "docs"
    / "architecture"
    / "examples"
    / "registered-result-read-v1.json"
)
EXAMPLES = json.loads(EXAMPLES_PATH.read_text(encoding="utf-8"))
READ_RECORD_KEYS = {
    "external_v1",
    "external_renamed",
    "external_upstream_v2",
    "downstream_after_upstream_v2",
    "missing_source",
}
RECORDS = {key: EXAMPLES[key] for key in READ_RECORD_KEYS}
DOWNLOAD_EXAMPLES = EXAMPLES["download_examples"]

PUBLIC_KEYS = {
    "project_id",
    "result_id",
    "result_revision_id",
    "revision_number",
    "revision_status",
    "display_name",
    "format",
    "schema_version",
    "origin",
    "method",
    "reference",
    "bytes",
    "created_at",
}
ID_PATTERNS = {
    "project_id": r"prj_[0-9a-f]{32}",
    "result_id": r"rreg_[0-9a-f]{32}",
    "result_revision_id": r"rrev_[0-9a-f]{32}",
}
FORBIDDEN_PUBLIC_KEYS = {
    "run_id",
    "artifact_generation",
    "artifact_revision",
    "publication_id",
    "storage_pool_id",
    "config_key",
    "relative_path",
    "absolute_path",
    "file_uri",
}


def _walk(value: Any):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from _walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk(child)


def _assert_id(value: str, field: str) -> None:
    assert re.fullmatch(ID_PATTERNS[field], value), (field, value)


def test_public_examples_have_only_bounded_revision_fields() -> None:
    assert EXAMPLES["example_schema"] == "registered-result-read-examples-v1"
    assert set(RECORDS) == READ_RECORD_KEYS
    for record in RECORDS.values():
        assert set(record) == PUBLIC_KEYS
        assert record["schema_version"] == 1
        assert record["revision_status"] in {"current", "superseded"}
        for field, pattern in ID_PATTERNS.items():
            _assert_id(record[field], field)
        assert re.fullmatch(r"[0-9a-f]{64}", record["bytes"]["content_sha256"])
        assert record["bytes"]["size_bytes"] > 0


def test_public_examples_do_not_invent_runs_or_expose_storage() -> None:
    for mapping in _walk(EXAMPLES):
        assert FORBIDDEN_PUBLIC_KEYS.isdisjoint(mapping)


def test_rename_appends_metadata_without_changing_evidence() -> None:
    original = RECORDS["external_v1"]
    renamed = RECORDS["external_renamed"]

    assert renamed["result_id"] == original["result_id"]
    assert renamed["result_revision_id"] != original["result_revision_id"]
    assert renamed["revision_number"] == original["revision_number"] + 1
    assert renamed["display_name"] != original["display_name"]
    for field in ("project_id", "format", "origin", "method", "reference", "bytes"):
        assert renamed[field] == original[field]


def test_downstream_parents_stay_pinned_after_upstream_replacement() -> None:
    original = RECORDS["external_v1"]
    upstream_v2 = RECORDS["external_upstream_v2"]
    downstream = RECORDS["downstream_after_upstream_v2"]
    parents = downstream["origin"]["parents"]

    assert downstream["origin"]["kind"] == "downstream"
    assert len(parents) == 2
    assert [parent["role"] for parent in parents] == ["primary", "control"]
    assert len({parent["result_revision_id"] for parent in parents}) == 2
    for parent in parents:
        _assert_id(parent["result_id"], "result_id")
        _assert_id(parent["result_revision_id"], "result_revision_id")
        _assert_id(parent["current_result_revision_id"], "result_revision_id")

    primary = parents[0]
    assert primary["result_revision_id"] == original["result_revision_id"]
    assert primary["display_name"] == original["display_name"]
    assert primary["revision_status"] == "superseded"
    assert primary["current_result_revision_id"] == upstream_v2["result_revision_id"]

    assert upstream_v2["result_id"] == original["result_id"]
    assert upstream_v2["bytes"] != original["bytes"]
    assert (
        upstream_v2["origin"]["source"]["accession"]
        == (original["origin"]["source"]["accession"])
    )
    assert (
        upstream_v2["origin"]["source"]["version"]
        != (original["origin"]["source"]["version"])
    )


def test_missing_source_and_other_unknowns_remain_readable_metadata() -> None:
    record = RECORDS["missing_source"]

    assert record["format"] == "h5ad"
    assert record["origin"] == {"kind": "external", "source": {"status": "unknown"}}
    assert record["method"] == {"status": "unknown"}
    assert record["reference"] == {"status": "unknown"}
    assert record["bytes"]["content_sha256"]
    assert record["bytes"]["size_bytes"] == 2048


def test_download_reauthenticates_old_links_and_audits_exact_revision() -> None:
    active = DOWNLOAD_EXAMPLES["active_member_old_revision"]
    revoked = DOWNLOAD_EXAMPLES["same_link_after_revocation"]
    anonymous = DOWNLOAD_EXAMPLES["anonymous_same_link"]
    wrong_project = DOWNLOAD_EXAMPLES["wrong_project"]
    original = RECORDS["external_v1"]

    assert active["request"] == {
        "project_id": original["project_id"],
        "result_revision_id": original["result_revision_id"],
    }
    assert active["response_status"] == 200
    assert active["response_content_sha256"] == original["bytes"]["content_sha256"]
    assert active["successful_download_audit_count"] == 1
    assert active["audit_target_identity_inputs"] == {
        "project_id": original["project_id"],
        "result_id": original["result_id"],
        "result_revision_id": original["result_revision_id"],
    }

    event = active["audit_event"]
    assert set(event) == {
        "event_id",
        "occurred_at",
        "action",
        "outcome",
        "actor_kind",
        "actor_user_id",
        "resource",
        "reason_code",
    }
    assert event["action"] == "registered_result.download"
    assert event["outcome"] == "succeeded"
    assert event["actor_user_id"] == active["principal"]["user_id"]
    assert re.fullmatch(r"usr_[0-9a-f]{32}", event["actor_user_id"])
    assert event["resource"]["kind"] == "registered_result"
    assert re.fullmatch(
        r"ares_registered_result_[0-9a-f]{64}",
        event["resource"]["resource_id"],
    )
    assert event["reason_code"] is None

    assert revoked["request"] == active["request"]
    assert revoked["principal"]["user_id"] == active["principal"]["user_id"]
    assert revoked["response_status"] == 401
    assert revoked["error_code"] == "AUTHENTICATION_REQUIRED"
    assert revoked["successful_download_audit_count"] == 0

    assert anonymous["request"] == active["request"]
    assert anonymous["response_status"] == 401
    assert anonymous["successful_download_audit_count"] == 0

    assert (
        wrong_project["request"]["result_revision_id"]
        == (active["request"]["result_revision_id"])
    )
    assert wrong_project["request"]["project_id"] != active["request"]["project_id"]
    assert wrong_project["response_status"] == 404
    assert wrong_project["error_code"] == "catalog_not_found"
    assert wrong_project["successful_download_audit_count"] == 0


def test_source_mutation_cannot_change_a_verified_response() -> None:
    scenario = EXAMPLES["source_mutation_after_verification"]
    original = RECORDS["external_v1"]

    assert scenario["request"]["result_revision_id"] == original["result_revision_id"]
    assert scenario["response_private_copy_verified"] is True
    assert scenario["source_changed_in_place_after_verification"] is True
    assert scenario["response_status"] == 200
    assert scenario["expected_content_sha256"] == original["bytes"]["content_sha256"]
    assert scenario["response_content_sha256"] == scenario["expected_content_sha256"]
