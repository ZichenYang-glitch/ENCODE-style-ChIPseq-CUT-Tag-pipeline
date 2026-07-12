"""Strict mappings for ENCODE machine-readable QC summary artifacts."""

from __future__ import annotations

import csv
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import io
import re

from encode_pipeline.platform.adapters import (
    ExtractedQcMetricCandidate,
    QcSourceDocument,
)


QC_SOURCE_OUTPUT_TYPES = (
    "mnase_qc_summary",
    "pooled_qc_summary",
    "qc_summary",
)

_QC_HEADER = (
    "sample",
    "assay",
    "target",
    "genome",
    "layout",
    "peak_mode",
    "use_control",
    "control_type",
    "final_bam",
    "peaks",
    "blacklist",
    "blacklist_filtered_bam",
    "blacklist_filtered_peaks",
    "total_reads",
    "reads_in_peaks",
    "frip",
    "peak_count",
    "blacklist_filtered_peak_count",
    "metrics_source",
    "unpaired_reads_examined",
    "read_pairs_examined",
    "secondary_or_supplementary_reads",
    "unmapped_reads",
    "unpaired_read_duplicates",
    "read_pair_duplicates",
    "read_pair_optical_duplicates",
    "percent_duplication",
    "estimated_library_size",
    "total_reads_examined",
    "duplicate_reads_estimate",
    "total_fragments",
    "distinct_fragments",
    "one_read_fragments",
    "two_read_fragments",
    "nrf",
    "pbc1",
    "pbc2",
)

_MNASE_HEADER = (
    "sample",
    "assay",
    "peak_mode",
    "sub_min",
    "sub_max",
    "mono_min",
    "mono_max",
    "di_min",
    "di_max",
    "dyad_min",
    "dyad_max",
    "sub_bam",
    "mono_bam",
    "di_bam",
    "dyad_bigwig",
    "mono_bigwig",
    "insert_size_metrics",
    "caller_danpos3_enabled",
    "caller_inps_enabled",
    "caller_sem_enabled",
    "sub_reads",
    "mono_reads",
    "di_reads",
)

_POOLED_HEADER = (
    "experiment",
    "assay",
    "target",
    "inferred_histone_class",
    "expected_peak_mode",
    "configured_peak_mode",
    "peak_mode_status",
    "biological_replicates",
    "biological_replicate_labels",
    "pooled_bam",
    "pooled_peaks",
    "pooled_peak_count",
    "pooled_FE_bdg",
    "pooled_ppois_bdg",
    "signal_tracks_status",
)

_HEADERS = {
    "qc_summary": _QC_HEADER,
    "mnase_qc_summary": _MNASE_HEADER,
    "pooled_qc_summary": _POOLED_HEADER,
}

_DECIMAL_PATTERN = re.compile(r"^(?:0|[1-9]\d{0,25})(?:\.\d{1,12})?$")


class EncodeQcSummaryError(ValueError):
    """A source document did not satisfy the fixed ENCODE QC contract."""


@dataclass(frozen=True)
class _MetricSpec:
    column: str
    metric_key: str
    display_name: str
    unit: str
    integral: bool = False
    maximum: Decimal | None = None
    positive: bool = False


_SAMPLE_METRICS = (
    _MetricSpec("total_reads", "sequencing.total_reads", "Total reads", "count", True),
    _MetricSpec(
        "frip",
        "peaks.frip",
        "Fraction of reads in peaks",
        "fraction",
        maximum=Decimal("1"),
    ),
    _MetricSpec("peak_count", "peaks.count", "Peak count", "count", True),
    _MetricSpec(
        "percent_duplication",
        "library.percent_duplication",
        "Percent duplication",
        "fraction",
        maximum=Decimal("1"),
    ),
    _MetricSpec(
        "estimated_library_size",
        "library.estimated_size",
        "Estimated library size",
        "count",
        True,
    ),
    _MetricSpec(
        "nrf",
        "library.nrf",
        "Non-redundant fraction",
        "fraction",
        maximum=Decimal("1"),
    ),
    _MetricSpec(
        "pbc1",
        "library.pbc1",
        "PCR bottleneck coefficient 1",
        "fraction",
        maximum=Decimal("1"),
    ),
    _MetricSpec(
        "pbc2",
        "library.pbc2",
        "PCR bottleneck coefficient 2",
        "ratio",
    ),
)

_MNASE_METRICS = (
    _MetricSpec(
        "sub_reads",
        "mnase.subnucleosomal_reads",
        "Sub-nucleosomal reads",
        "count",
        True,
    ),
    _MetricSpec(
        "mono_reads",
        "mnase.mononucleosomal_reads",
        "Mono-nucleosomal reads",
        "count",
        True,
    ),
    _MetricSpec(
        "di_reads",
        "mnase.dinucleosomal_reads",
        "Di-nucleosomal reads",
        "count",
        True,
    ),
)

_POOLED_METRICS = (
    _MetricSpec(
        "pooled_peak_count",
        "peaks.pooled_count",
        "Pooled peak count",
        "count",
        True,
    ),
    _MetricSpec(
        "biological_replicates",
        "replicates.biological_count",
        "Biological replicates",
        "count",
        True,
        positive=True,
    ),
)


def parse_encode_qc_sources(
    sources: tuple[QcSourceDocument, ...],
) -> tuple[ExtractedQcMetricCandidate, ...]:
    """Parse exact ENCODE summary contracts without accessing the filesystem."""
    if not isinstance(sources, tuple):
        raise EncodeQcSummaryError("sources must be a tuple")
    seen_ids: set[str] = set()
    seen_paths: set[str] = set()
    candidates: list[ExtractedQcMetricCandidate] = []
    for document in sources:
        if not isinstance(document, QcSourceDocument):
            raise EncodeQcSummaryError("source document type is invalid")
        source = document.source
        if source.artifact_id in seen_ids or source.relative_path in seen_paths:
            raise EncodeQcSummaryError("source artifact is duplicated")
        seen_ids.add(source.artifact_id)
        seen_paths.add(source.relative_path)
        header = _HEADERS.get(source.output_type)
        if header is None:
            raise EncodeQcSummaryError("source output type is unsupported")
        row = _read_one_row(document.content, header)
        if source.output_type == "qc_summary":
            candidates.extend(_sample_candidates(source, row, _SAMPLE_METRICS))
        elif source.output_type == "mnase_qc_summary":
            candidates.extend(_sample_candidates(source, row, _MNASE_METRICS))
        else:
            candidates.extend(_experiment_candidates(source, row))
    return tuple(
        sorted(
            candidates,
            key=lambda item: (
                item.scope,
                item.sample_id or "",
                item.experiment_id or "",
                item.metric_key,
            ),
        )
    )


def _read_one_row(content: bytes, header: tuple[str, ...]) -> dict[str, str]:
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise EncodeQcSummaryError("source is not UTF-8") from exc
    reader = csv.DictReader(io.StringIO(text, newline=""), delimiter="\t")
    if tuple(reader.fieldnames or ()) != header:
        raise EncodeQcSummaryError("source header is invalid")
    rows = list(reader)
    if len(rows) != 1 or None in rows[0]:
        raise EncodeQcSummaryError("source must contain exactly one row")
    return rows[0]


def _sample_candidates(source, row, specs):
    metadata = source.metadata
    sample_id = _required_metadata(metadata, "sample_id")
    assay = _required_metadata(metadata, "assay")
    scope = _required_metadata(metadata, "scope")
    if scope != "sample" or row.get("sample") != sample_id or row.get("assay") != assay:
        raise EncodeQcSummaryError("sample source identity is inconsistent")
    experiment_id = _optional_metadata(metadata, "experiment_id")
    return tuple(
        _candidate(
            source.artifact_id,
            spec,
            row,
            scope="sample",
            sample_id=sample_id,
            experiment_id=experiment_id,
            assay=assay,
        )
        for spec in specs
        if row.get(spec.column) != "NA"
    )


def _experiment_candidates(source, row):
    metadata = source.metadata
    experiment_id = _required_metadata(metadata, "experiment_id")
    assay = _required_metadata(metadata, "assay")
    scope = _required_metadata(metadata, "scope")
    if (
        scope != "experiment"
        or row.get("experiment") != experiment_id
        or row.get("assay") != assay
    ):
        raise EncodeQcSummaryError("experiment source identity is inconsistent")
    status = row.get("peak_mode_status")
    flags = {"ok": "pass", "mismatch": "warning", "unknown": None}
    if status not in flags:
        raise EncodeQcSummaryError("pooled QC status is invalid")
    return tuple(
        _candidate(
            source.artifact_id,
            spec,
            row,
            scope="experiment",
            experiment_id=experiment_id,
            assay=assay,
            qc_flag=flags[status] if spec.metric_key == "peaks.pooled_count" else None,
        )
        for spec in _POOLED_METRICS
        if row.get(spec.column) != "NA"
    )


def _candidate(
    source_artifact_id,
    spec,
    row,
    *,
    scope,
    sample_id=None,
    experiment_id=None,
    assay=None,
    qc_flag=None,
):
    raw = row.get(spec.column)
    if not isinstance(raw, str) or not raw or _DECIMAL_PATTERN.fullmatch(raw) is None:
        raise EncodeQcSummaryError("QC numeric value is invalid")
    try:
        value = Decimal(raw)
    except InvalidOperation as exc:
        raise EncodeQcSummaryError("QC numeric value is invalid") from exc
    if not value.is_finite() or value < 0:
        raise EncodeQcSummaryError("QC numeric value is invalid")
    if spec.integral and value != value.to_integral_value():
        raise EncodeQcSummaryError("QC count is not integral")
    if spec.maximum is not None and value > spec.maximum:
        raise EncodeQcSummaryError("QC numeric value is out of range")
    if spec.positive and value <= 0:
        raise EncodeQcSummaryError("QC numeric value must be positive")
    return ExtractedQcMetricCandidate(
        metric_key=spec.metric_key,
        display_name=spec.display_name,
        value=value,
        unit=spec.unit,
        scope=scope,
        sample_id=sample_id,
        experiment_id=experiment_id,
        assay=assay,
        qc_flag=qc_flag,
        source_artifact_id=source_artifact_id,
    )


def _required_metadata(metadata, key):
    value = metadata.get(key)
    if not isinstance(value, str) or not value:
        raise EncodeQcSummaryError("source metadata is incomplete")
    return value


def _optional_metadata(metadata, key):
    value = metadata.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value:
        raise EncodeQcSummaryError("source metadata is invalid")
    return value
