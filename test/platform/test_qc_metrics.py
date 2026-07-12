"""Tests for workflow-neutral durable QC metric values."""

from datetime import datetime, timezone
from decimal import Decimal

import pytest

from encode_pipeline.platform.runs import RunQcMetric


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
