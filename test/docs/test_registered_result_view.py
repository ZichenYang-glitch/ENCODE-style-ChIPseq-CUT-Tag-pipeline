"""Lock the IRENE-40 view candidates to the pinned registration bytes."""

from __future__ import annotations

import hashlib
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
    / "registered-result-view-v1.json"
)
SPEC_PATH = REPOSITORY_ROOT / "docs" / "architecture" / "registered-result-view-v1.md"
EXAMPLES = json.loads(EXAMPLES_PATH.read_text(encoding="utf-8"))
SPEC = SPEC_PATH.read_text(encoding="utf-8")
PINNED_COMMIT = "29d508647c18df94016f94d7da437e794b846f9e"
VIEWER_BASELINE = "f52a2c20adae6c6e61375b34e47211c84aa035f7"
CLUSTER_SHA256 = "d6cb16f73da05e9d2e6f9a70e838dc9955dafdf560a955f6f62bb1dc2907817a"
HIC_REVISION = "rrev_11111111111111111111111111111111"
FORBIDDEN_KEYS = {
    "run_id",
    "artifact_generation",
    "artifact_revision",
    "storage_pool_id",
    "absolute_path",
    "contact_value",
}


def _candidates() -> dict[str, dict[str, Any]]:
    items = EXAMPLES["candidates"]
    assert [item["action_id"] for item in items] == [
        "view-hic-region-resolution",
        "view-single-cell-qc-cluster-counts",
    ]
    return {item["action_id"]: item for item in items}


def _walk(value: Any):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from _walk(child)
    elif isinstance(value, list):
        for child in value:
            yield from _walk(child)


def test_view_examples_stay_design_only() -> None:
    assert EXAMPLES["example_schema"] == "registered-result-view-examples-v1"
    assert EXAMPLES["pinned_registration_commit"] == PINNED_COMMIT
    assert EXAMPLES["viewer_baseline_commit"] == VIEWER_BASELINE
    assert EXAMPLES["adapters_added_this_round"] == []
    assert EXAMPLES["implementation_authorized"] is False
    assert EXAMPLES["awaiting_lab_choice"] is True
    assert EXAMPLES["decisions_independent"] is True


def test_examples_do_not_invent_runs_or_contact_values() -> None:
    for mapping in _walk(EXAMPLES):
        assert FORBIDDEN_KEYS.isdisjoint(mapping)


def test_hic_candidate_reuses_the_run_track_without_claiming_a_matrix() -> None:
    hic = _candidates()["view-hic-region-resolution"]
    reader = hic["reusable_reader"]

    assert hic["required_format"] == "hic"
    assert reader["wash_u_type"] == "hic"
    assert reader["session_scope"] == "run_artifact"
    assert reader["passes_resolution"] is False
    assert reader["passes_normalization"] is False
    assert reader["can_open_registered_revision"] is False
    assert hic["request"]["result_revision_id"] == HIC_REVISION
    assert hic["request"]["region"] == "chr1:0-2000000"
    assert hic["request"]["resolution_bp"] == 1000000
    assert hic["request"]["normalization"] == "NONE"
    assert hic["registration_example"]["record"] == "external_v1"
    assert hic["registration_example"]["bytes_are_shape_illustration"] is True
    assert hic["registration_example"]["matrix_observed"] is False
    assert hic["current_outcome"] == {
        "status": "unsupported",
        "reason": "registered_revision_has_no_matrix_reader",
    }
    success = hic["success_when_implemented"]
    assert success["shows_requested_region"] is True
    assert success["shows_requested_resolution_bp"] is True
    assert success["shows_requested_normalization"] is True
    assert success["shows_one_contact_bound_to_revision_sha256"] is True
    assert success["same_contact_on_embedded_and_full_page"] is True
    assert success["unknown_reference_label"] == "未知"
    assert success["paint_only_after_sha256_match"] is True
    reasons = {item["case"]: item for item in hic["unsupported"]}
    assert reasons["cool"]["reason"] == "unsupported_contact_matrix"
    assert reasons["mcool"]["reason"] == "unsupported_contact_matrix"
    assert reasons["resolution_absent"]["page_state"] == "requested_bin_absent"
    assert reasons["normalization_absent"]["page_state"] == (
        "requested_normalization_absent"
    )
    assert reasons["bytes_unavailable"]["http_status"] == 409
    assert reasons["bytes_unavailable"]["page_state"] == "no_matrix"
    assert reasons["non_matrix_format"]["formats"] == ["tsv", "h5ad"]


def test_single_cell_candidate_is_the_verified_cluster_count_table() -> None:
    single_cell = _candidates()["view-single-cell-qc-cluster-counts"]
    payload = single_cell["bytes"]["media"].encode("utf-8")

    assert single_cell["choice"] == "qc"
    assert single_cell["choices_considered"] == ["qc", "embedding", "expression"]
    assert single_cell["required_format"] == "tsv"
    assert single_cell["required_columns"] == ["cluster_id", "cell_count"]
    assert single_cell["reusable_reader"]["accept_ranges"] == "none"
    assert single_cell["reusable_reader"]["parses_table"] is False
    assert single_cell["reusable_reader"]["can_return_these_bytes"] is True
    assert len(payload) == 50
    assert hashlib.sha256(payload).hexdigest() == CLUSTER_SHA256
    assert payload.decode("utf-8") == single_cell["bytes"]["media"]
    assert single_cell["bytes"]["content_sha256"] == CLUSTER_SHA256
    assert single_cell["bytes"]["present_on_viewer_baseline"] is False
    assert single_cell["bytes"]["live_result_revision_id_in_repository"] is None
    assert single_cell["bytes"]["method"]["label"] == "未知"
    assert single_cell["bytes"]["reference"]["label"] == "未知"
    assert single_cell["current_outcome"]["status"] == (
        "bytes_verified_table_not_shown"
    )
    assert single_cell["success_when_implemented"]["rows"] == [
        {"cluster_id": "synthetic-A", "cell_count": 3},
        {"cluster_id": "synthetic-B", "cell_count": 2},
    ]
    assert single_cell["success_when_implemented"][
        "same_rows_on_embedded_and_full_page"
    ]
    assert {item["choice"]: item["reason"] for item in single_cell["not_selected"]} == {
        "embedding": "no_reader",
        "expression": "h5ad_is_not_parsed",
    }
    unsupported = {item["case"]: item for item in single_cell["unsupported"]}
    assert unsupported["h5ad"]["page_state"] == "no_reader"
    assert unsupported["unexpected_columns"]["page_state"] == "unsupported_table"
    assert unsupported["bytes_unavailable"]["page_state"] == "no_table"


def test_spec_keeps_the_same_pins_as_the_example() -> None:
    required = (
        PINNED_COMMIT,
        VIEWER_BASELINE,
        "view-hic-region-resolution",
        "view-single-cell-qc-cluster-counts",
        "Accept-Ranges: none",
        "unsupported_contact_matrix",
        CLUSTER_SHA256,
        "synthetic-A",
        "本轮不新增 adapter",
        "registered_revision_has_no_matrix_reader",
        "bytes_verified_table_not_shown",
    )
    missing = [item for item in required if item not in SPEC]
    assert not missing
    assert re.search(r"adapters/", SPEC)
    assert "接触值 0" in SPEC
    assert "细胞数 0" in SPEC
