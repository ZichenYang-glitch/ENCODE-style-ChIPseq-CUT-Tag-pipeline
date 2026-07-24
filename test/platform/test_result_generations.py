"""Contract tests for public-safe result generation identities."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone
from decimal import Decimal

import pytest

from encode_pipeline.platform.result_generations import (
    ArtifactCursor,
    QcMetricCursor,
    RunResultState,
    build_artifact_content_revision,
    build_artifact_generation,
    build_qc_generation,
    decode_artifact_cursor,
    decode_qc_metric_cursor,
    encode_artifact_cursor,
    encode_qc_metric_cursor,
    new_result_attempt_id,
)
from encode_pipeline.platform.runs import (
    RunArtifactRef,
    RunQcMetric,
    build_qc_metric_id,
)


NOW = datetime(2026, 7, 17, tzinfo=timezone.utc)


def _artifact(content: bytes, *, revision_number: int = 1) -> RunArtifactRef:
    revision = build_artifact_content_revision(
        output_type="qc_summary",
        relative_path="results/qc/summary.tsv",
        content=content,
    )
    return RunArtifactRef(
        artifact_id="artifact-summary",
        run_id="run-1",
        artifact_type="file",
        name="summary.tsv",
        uri="run://runs/run-1/artifacts/artifact-summary",
        mime_type="text/tab-separated-values",
        produced_at=NOW,
        revision=revision,
        metadata={
            "output_type": "qc_summary",
            "relative_path": "results/qc/summary.tsv",
            "size_bytes": len(content),
        },
    )


def _metric(value: str) -> RunQcMetric:
    metric_key = "sequencing.total_reads"
    return RunQcMetric(
        metric_id=build_qc_metric_id(metric_key, "sample", "S1", None),
        run_id="run-1",
        metric_key=metric_key,
        display_name="Total reads",
        value=Decimal(value),
        unit="count",
        scope="sample",
        sample_id="S1",
        experiment_id=None,
        assay="rnaseq",
        qc_flag=None,
        source_artifact_id="artifact-summary",
        produced_at=NOW,
    )


def test_content_revision_detects_same_length_replacement_without_path_disclosure():
    first = build_artifact_content_revision(
        output_type="qc_summary",
        relative_path="results/qc/summary.tsv",
        content=b"100\n",
    )
    second = build_artifact_content_revision(
        output_type="qc_summary",
        relative_path="results/qc/summary.tsv",
        content=b"900\n",
    )

    assert first != second
    assert first.startswith("artifactrev-")
    assert len(first) == len("artifactrev-") + 64
    assert "results" not in first


def test_monotonic_revision_prevents_artifact_and_qc_generation_aba():
    artifact_a = _artifact(b"100\n")
    artifact_b = _artifact(b"900\n")
    manifest_a = (artifact_a,)
    manifest_b = (artifact_b,)

    generation_a1 = build_artifact_generation(
        run_id="run-1", revision=1, artifacts=manifest_a
    )
    generation_b = build_artifact_generation(
        run_id="run-1", revision=2, artifacts=manifest_b
    )
    generation_a2 = build_artifact_generation(
        run_id="run-1", revision=3, artifacts=manifest_a
    )

    assert len({generation_a1, generation_b, generation_a2}) == 3
    assert generation_a1.startswith("artifactgen-")

    qc_a1 = build_qc_generation(
        run_id="run-1",
        revision=1,
        artifact_generation=generation_a1,
        metrics=(_metric("100"),),
    )
    qc_a2 = build_qc_generation(
        run_id="run-1",
        revision=2,
        artifact_generation=generation_a2,
        metrics=(_metric("100"),),
    )
    assert qc_a1 != qc_a2
    assert qc_a1.startswith("qcgen-")


def test_qc_cursor_round_trip_is_canonical_and_generation_bound():
    generation = "qcgen-" + "a" * 64
    cursor = QcMetricCursor(
        run_id="run-1",
        qc_generation=generation,
        after_metric_id="qcmetric-" + "b" * 64,
    )

    encoded = encode_qc_metric_cursor(cursor)

    assert encoded.startswith("qccur_")
    assert (
        decode_qc_metric_cursor(
            encoded,
            run_id="run-1",
            qc_generation=generation,
        )
        == cursor
    )
    assert "/private" not in encoded
    with pytest.raises(ValueError):
        decode_qc_metric_cursor(
            encoded,
            run_id="run-2",
            qc_generation=generation,
        )
    with pytest.raises(ValueError):
        decode_qc_metric_cursor(
            encoded,
            run_id="run-1",
            qc_generation="qcgen-" + "c" * 64,
        )


def test_artifact_cursor_round_trip_is_canonical_and_generation_bound():
    generation = "artifactgen-" + "a" * 64
    cursor = ArtifactCursor(
        run_id="run-1",
        artifact_generation=generation,
        after_artifact_id="artifact-summary",
    )

    encoded = encode_artifact_cursor(cursor)

    assert encoded.startswith("artifactcur_")
    assert (
        decode_artifact_cursor(
            encoded,
            run_id="run-1",
            artifact_generation=generation,
        )
        == cursor
    )
    assert "/private" not in encoded
    with pytest.raises(ValueError):
        decode_artifact_cursor(
            encoded,
            run_id="run-2",
            artifact_generation=generation,
        )
    with pytest.raises(ValueError):
        decode_artifact_cursor(
            encoded,
            run_id="run-1",
            artifact_generation="artifactgen-" + "c" * 64,
        )


@pytest.mark.parametrize(
    "value",
    (
        "artifact-summary",
        "artifactcur_",
        "artifactcur_../private",
        "artifactcur_" + "a" * 1025,
    ),
)
def test_artifact_cursor_rejects_legacy_malformed_and_oversized_values(value):
    with pytest.raises(ValueError):
        decode_artifact_cursor(
            value,
            run_id="run-1",
            artifact_generation="artifactgen-" + "a" * 64,
        )


@pytest.mark.parametrize(
    "value",
    (
        "qcmetric-" + "a" * 64,
        "qccur_",
        "qccur_../private",
        "qccur_" + "a" * 1025,
    ),
)
def test_qc_cursor_rejects_legacy_malformed_and_oversized_values(value):
    with pytest.raises(ValueError):
        decode_qc_metric_cursor(
            value,
            run_id="run-1",
            qc_generation="qcgen-" + "a" * 64,
        )


def test_unbound_result_state_preserves_opaque_run_id_until_results_exist():
    state = RunResultState(run_id="..")

    assert state.run_id == ".."
    with pytest.raises(ValueError, match="result generation run_id is invalid"):
        replace(
            state,
            artifact_attempt_id=new_result_attempt_id(),
            artifact_attempt_status="pending",
        )


@pytest.mark.parametrize(
    ("changes", "message"),
    (
        (
            {"artifact_generation": "artifactgen-" + "a" * 64},
            "artifact generation and manifest digest must be paired",
        ),
        (
            {"qc_generation": "qcgen-" + "a" * 64},
            "QC generation and manifest digest must be paired",
        ),
        (
            {"qc_artifact_generation": "artifactgen-" + "a" * 64},
            "unbound QC state cannot name an artifact generation",
        ),
        (
            {
                "qc_generation": "qcgen-" + "a" * 64,
                "qc_manifest_digest": "b" * 64,
                "qc_artifact_generation": "artifactgen-" + "c" * 64,
            },
            "bound QC state requires a positive revision",
        ),
        (
            {
                "qc_attempt_artifact_generation": "artifactgen-" + "a" * 64,
            },
            "QC attempt generation requires an attempt",
        ),
        (
            {
                "qc_attempt_id": "resultattempt-" + "a" * 64,
                "qc_attempt_status": "pending",
            },
            "QC attempt requires an artifact generation",
        ),
        (
            {"artifact_attempt_id": "resultattempt-" + "a" * 64},
            "artifact attempt ID and status must be paired",
        ),
        (
            {
                "artifact_attempt_id": "resultattempt-" + "a" * 64,
                "artifact_attempt_status": "unknown",
            },
            "artifact attempt status is invalid",
        ),
        (
            {"artifact_outcome": "failed"},
            "failed artifact outcome requires a reason",
        ),
        (
            {
                "artifact_outcome": "succeeded",
                "artifact_reason_code": "STALE_REASON",
            },
            "non-failed artifact outcome cannot have a reason",
        ),
    ),
)
def test_result_state_rejects_incomplete_persisted_generation_closure(
    changes,
    message,
):
    state = RunResultState(run_id="run-1")

    with pytest.raises(ValueError, match=message):
        replace(state, **changes)
