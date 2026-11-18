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
RECORDS = {key: value for key, value in EXAMPLES.items() if key != "example_schema"}

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
    assert set(RECORDS) == {
        "external_v1",
        "external_renamed",
        "external_upstream_v2",
        "downstream_after_upstream_v2",
        "missing_source",
    }
    for record in RECORDS.values():
        assert set(record) == PUBLIC_KEYS
        assert record["schema_version"] == 1
        assert record["revision_status"] in {"current", "superseded"}
        for field, pattern in ID_PATTERNS.items():
            _assert_id(record[field], field)
        assert re.fullmatch(r"[0-9a-f]{64}", record["bytes"]["content_sha256"])
        assert record["bytes"]["size_bytes"] > 0


def test_public_examples_do_not_invent_runs_or_expose_storage() -> None:
    for record in RECORDS.values():
        for mapping in _walk(record):
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
