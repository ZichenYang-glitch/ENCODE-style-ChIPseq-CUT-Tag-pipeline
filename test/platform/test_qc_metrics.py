"""Tests for workflow-neutral durable QC metric values."""

from datetime import datetime, timezone
from decimal import Decimal

import pytest

from encode_pipeline.platform.runs import (
    RunQcMetric,
    build_qc_metric_id,
    validate_qc_identifier_token,
)


def _metric(**overrides):
    values = {
        "metric_id": "qcmetric-1",
        "run_id": "run-1",
        "metric_key": "peaks.frip",
        "display_name": "Fraction of reads in peaks",
        "value": Decimal("0.125"),
        "unit": "fraction",
        "scope": "sample",
        "sample_id": "S1",
        "experiment_id": "EXP1",
        "assay": "chipseq",
        "qc_flag": "warning",
        "source_artifact_id": "artifact-1",
        "produced_at": datetime(2026, 7, 12, tzinfo=timezone.utc),
    }
    values.update(overrides)
    return RunQcMetric(**values)


def test_qc_metric_serializes_explicit_fields_without_metadata_container():
    metric = _metric()

    assert metric.to_dict() == {
        "metric_id": "qcmetric-1",
        "run_id": "run-1",
        "metric_key": "peaks.frip",
        "display_name": "Fraction of reads in peaks",
        "value": Decimal("0.125"),
        "unit": "fraction",
        "scope": "sample",
        "sample_id": "S1",
        "experiment_id": "EXP1",
        "assay": "chipseq",
        "qc_flag": "warning",
        "source_artifact_id": "artifact-1",
        "produced_at": datetime(2026, 7, 12, tzinfo=timezone.utc),
    }
    assert "metadata" not in metric.to_dict()


def test_qc_metric_rejects_non_decimal_values():
    with pytest.raises(ValueError, match="value must be a Decimal"):
        _metric(value=0.125)


def test_qc_metric_id_preserves_framed_sha256_contract():
    assert build_qc_metric_id("peaks.frip", "sample", "S1", "EXP1") == (
        "qcmetric-c44b0a8edd708eaf71ac8975390b0f6750438d00d72e021665324c01fc97a826"
    )


def test_qc_metric_id_rejects_non_string_coordinates():
    with pytest.raises(ValueError, match="coordinates"):
        build_qc_metric_id("peaks.frip", "sample", 0, "EXP1")


@pytest.mark.parametrize("value", ("-S1", ".S1", "_S1", "S1"))
def test_qc_identifier_token_accepts_canonical_safe_characters(value):
    assert validate_qc_identifier_token(value) == value


@pytest.mark.parametrize(
    "value",
    ("", ".", "..", "/private", "private\\path", "bad\nidentifier", "A" * 256),
)
def test_qc_identifier_token_rejects_unsafe_values(value):
    with pytest.raises(ValueError, match="identifier token"):
        validate_qc_identifier_token(value)
