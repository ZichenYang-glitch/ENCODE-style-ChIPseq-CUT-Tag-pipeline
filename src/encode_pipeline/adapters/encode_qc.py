"""Strict mappings for ENCODE machine-readable QC artifacts.

Version 1 tables remain exactly readable.  Version 2 tables are identified by
their leading ``schema_version`` column and must match the exact v2 header;
unknown versions, duplicate columns, extra columns, wrong types and undefined
status values fail the source instead of becoming zeroes.
"""

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
from encode_pipeline.qc.contract import (
    MNASE_QC_V1_COLUMNS,
    MNASE_QC_V2_COLUMNS,
    POOLED_QC_COLUMNS,
    SAMPLE_QC_V1_COLUMNS,
    SAMPLE_QC_V2_COLUMNS,
)
from encode_pipeline.qc.policy import (
    AnalysisPolicyError,
    parse_analysis_policy_bytes,
)


QC_SOURCE_OUTPUT_TYPES = (
    "analysis_policy",
    "mnase_qc_summary",
    "pooled_qc_summary",
    "qc_summary",
)

_QC_HEADER = SAMPLE_QC_V1_COLUMNS
_QC_V2_HEADER = SAMPLE_QC_V2_COLUMNS
_MNASE_HEADER = MNASE_QC_V1_COLUMNS
_MNASE_V2_HEADER = MNASE_QC_V2_COLUMNS
_POOLED_HEADER = POOLED_QC_COLUMNS

_SAMPLE_OUTPUT_TYPES = frozenset({"qc_summary"})
_MNASE_OUTPUT_TYPES = frozenset({"mnase_qc_summary"})
_POOLED_OUTPUT_TYPES = frozenset({"pooled_qc_summary"})
_POLICY_OUTPUT_TYPES = frozenset(
    {"analysis_policy", "analysis_filtered_analysis_policy"}
)

_DECIMAL_PATTERN = re.compile(r"^(?:0|[1-9]\d{0,25})(?:\.\d{1,12})?$")
_SIGNED_DECIMAL_PATTERN = re.compile(r"^-?(?:0|[1-9]\d{0,25})(?:\.\d{1,12})?$")

_SAMPLE_STATUSES = frozenset(
    {
        "ok",
        "unresolved_contigs",
        "no_mapped_reads",
        "not_applicable",
        "missing",
        "parse_failed",
    }
)
_CROSS_CORRELATION_STATUSES = frozenset(
    {"reported", "not_enabled", "missing", "parse_failed"}
)
_MNASE_ANALYSIS_STATUSES = frozenset({"ok", "no_valid_templates", "parse_failed"})


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
    signed: bool = False


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

_SAMPLE_V2_METRICS = (
    _MetricSpec(
        "analysis_valid_fragments",
        "analysis.valid_fragments",
        "Analysis valid fragments",
        "count",
        True,
    ),
    _MetricSpec(
        "analysis_invalid_templates",
        "analysis.invalid_templates",
        "Analysis invalid templates",
        "count",
        True,
    ),
    _MetricSpec(
        "analysis_upstream_dropped_templates",
        "analysis.upstream_dropped_templates",
        "Upstream dropped templates",
        "count",
        True,
    ),
    _MetricSpec(
        "analysis_upstream_dropped_alignments",
        "analysis.upstream_dropped_alignments",
        "Upstream dropped alignments",
        "count",
        True,
    ),
    _MetricSpec("raw_peak_count", "peaks.raw_count", "Raw peak count", "count", True),
    _MetricSpec(
        "delivered_peak_count",
        "peaks.delivered_count",
        "Delivered peak count",
        "count",
        True,
    ),
    _MetricSpec("organelle_reads", "organelle.reads", "Organelle reads", "count", True),
    _MetricSpec(
        "organelle_total_mapped_reads",
        "organelle.total_mapped_reads",
        "Total mapped reads",
        "count",
        True,
    ),
    _MetricSpec(
        "organelle_fraction",
        "organelle.fraction",
        "Organelle read fraction",
        "fraction",
        maximum=Decimal("1"),
    ),
    _MetricSpec(
        "cross_correlation_nsc",
        "cross_correlation.nsc",
        "Normalized strand coefficient",
        "ratio",
        signed=True,
    ),
    _MetricSpec(
        "cross_correlation_rsc",
        "cross_correlation.rsc",
        "Relative strand correlation",
        "ratio",
        signed=True,
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

_MNASE_V2_METRICS = (
    _MetricSpec(
        "analysis_valid_templates",
        "mnase.valid_templates",
        "Valid MNase templates",
        "count",
        True,
    ),
    _MetricSpec(
        "analysis_invalid_templates",
        "mnase.invalid_templates",
        "Invalid MNase templates",
        "count",
        True,
    ),
    _MetricSpec(
        "analysis_classified_templates",
        "mnase.classified_templates",
        "Classified MNase templates",
        "count",
        True,
    ),
    _MetricSpec(
        "analysis_unclassified_templates",
        "mnase.unclassified_templates",
        "Unclassified MNase templates",
        "count",
        True,
    ),
    _MetricSpec(
        "analysis_fraction_classified",
        "mnase.fraction_classified",
        "Classified MNase template fraction",
        "fraction",
        maximum=Decimal("1"),
    ),
    _MetricSpec(
        "sub_fragments",
        "mnase.subnucleosomal_fragments",
        "Sub-nucleosomal fragments",
        "count",
        True,
    ),
    _MetricSpec(
        "mono_fragments",
        "mnase.mononucleosomal_fragments",
        "Mono-nucleosomal fragments",
        "count",
        True,
    ),
    _MetricSpec(
        "di_fragments",
        "mnase.dinucleosomal_fragments",
        "Di-nucleosomal fragments",
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

_POLICY_STAGE_METRICS = (
    (
        "analysis_input_alignments",
        "analysis.stage.input_alignments",
        "Analysis input alignments",
    ),
    (
        "analysis_kept_alignments",
        "analysis.stage.kept_alignments",
        "Analysis kept alignments",
    ),
    (
        "analysis_valid_templates",
        "analysis.stage.valid_templates",
        "Analysis valid templates",
    ),
    (
        "analysis_invalid_templates",
        "analysis.stage.invalid_templates",
        "Analysis invalid templates",
    ),
    (
        "upstream_dropped_templates",
        "analysis.stage.upstream_dropped_templates",
        "Upstream dropped templates",
    ),
    (
        "upstream_dropped_alignments",
        "analysis.stage.upstream_dropped_alignments",
        "Upstream dropped alignments",
    ),
    (
        "blacklist_dropped_templates",
        "analysis.stage.blacklist_dropped_templates",
        "Blacklist dropped templates",
    ),
    (
        "blacklist_dropped_alignments",
        "analysis.stage.blacklist_dropped_alignments",
        "Blacklist dropped alignments",
    ),
    (
        "organelle_reads",
        "analysis.stage.organelle_reads",
        "Organelle reads",
    ),
    (
        "organelle_total_mapped_reads",
        "analysis.stage.organelle_total_mapped_reads",
        "Organelle denominator reads",
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
        output_type = source.output_type
        if output_type in _SAMPLE_OUTPUT_TYPES:
            version, row = _read_versioned_row(
                document.content, _QC_HEADER, _QC_V2_HEADER
            )
            if version == "1":
                candidates.extend(_sample_candidates(source, row, _SAMPLE_METRICS))
            else:
                _validate_sample_v2(row)
                candidates.extend(
                    _sample_candidates(
                        source, row, _SAMPLE_METRICS + _SAMPLE_V2_METRICS
                    )
                )
        elif output_type in _MNASE_OUTPUT_TYPES:
            version, row = _read_versioned_row(
                document.content, _MNASE_HEADER, _MNASE_V2_HEADER
            )
            if version == "1":
                candidates.extend(_sample_candidates(source, row, _MNASE_METRICS))
            else:
                _validate_mnase_v2(row)
                candidates.extend(
                    _sample_candidates(source, row, _MNASE_METRICS + _MNASE_V2_METRICS)
                )
        elif output_type in _POOLED_OUTPUT_TYPES:
            row = _read_exact_row(document.content, _POOLED_HEADER)
            candidates.extend(_experiment_candidates(source, row))
        elif output_type in _POLICY_OUTPUT_TYPES:
            candidates.extend(_policy_candidates(source, document.content))
        else:
            raise EncodeQcSummaryError("source output type is unsupported")
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


def _read_versioned_row(
    content: bytes, v1_header: tuple[str, ...], v2_header: tuple[str, ...]
) -> tuple[str, dict[str, str]]:
    header, row = _read_header_and_row(content)
    if header == v1_header:
        version = "1"
    elif header == v2_header:
        version = "2"
    else:
        raise EncodeQcSummaryError("source header is invalid")
    return version, row


def _read_exact_row(content: bytes, header: tuple[str, ...]) -> dict[str, str]:
    actual, row = _read_header_and_row(content)
    if actual != header:
        raise EncodeQcSummaryError("source header is invalid")
    return row


def _read_header_and_row(content: bytes) -> tuple[tuple[str, ...], dict[str, str]]:
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeDecodeError as exc:
        raise EncodeQcSummaryError("source is not UTF-8") from exc
    reader = csv.reader(io.StringIO(text, newline=""), delimiter="\t")
    try:
        header = tuple(next(reader))
    except StopIteration as exc:
        raise EncodeQcSummaryError("source header is invalid") from exc
    if not header or len(set(header)) != len(header):
        raise EncodeQcSummaryError("source header is invalid")
    rows = [row for row in reader if row]
    if len(rows) != 1:
        raise EncodeQcSummaryError("source must contain exactly one row")
    if len(rows[0]) != len(header):
        raise EncodeQcSummaryError("source row width is invalid")
    return header, dict(zip(header, rows[0], strict=True))


def _validate_sample_v2(row: dict[str, str]) -> None:
    if row.get("schema_version") != "2":
        raise EncodeQcSummaryError("source schema version is unsupported")
    if row.get("analysis_input_stage") != "analysis_bam":
        raise EncodeQcSummaryError("analysis input stage is invalid")
    if row.get("duplicate_analysis") not in ("exclude", "retain"):
        raise EncodeQcSummaryError("duplicate policy is invalid")
    max_insert = row.get("resolved_max_insert", "")
    if max_insert != "not_applicable" and (
        not max_insert.isdigit() or int(max_insert) <= 0
    ):
        raise EncodeQcSummaryError("resolved max insert is invalid")
    if row.get("blacklist_mode") not in ("off", "copies", "analysis"):
        raise EncodeQcSummaryError("blacklist mode is invalid")
    if row.get("blacklist_applied") not in ("true", "false"):
        raise EncodeQcSummaryError("blacklist applied flag is invalid")
    if row.get("organelle_mode") not in ("auto", "exclude", "report", "keep"):
        raise EncodeQcSummaryError("organelle mode is invalid")
    organelle_status = row.get("organelle_status")
    if organelle_status not in _SAMPLE_STATUSES:
        raise EncodeQcSummaryError("organelle status is invalid")
    organelle_values = (
        row.get("organelle_reads", "NA"),
        row.get("organelle_total_mapped_reads", "NA"),
        row.get("organelle_fraction", "NA"),
    )
    if organelle_status == "ok":
        if any(value == "NA" for value in organelle_values):
            raise EncodeQcSummaryError("organelle values are incomplete")
    elif any(value != "NA" for value in organelle_values):
        raise EncodeQcSummaryError("organelle values are inconsistent")
    cc_status = row.get("cross_correlation_status")
    if cc_status not in _CROSS_CORRELATION_STATUSES:
        raise EncodeQcSummaryError("cross-correlation status is invalid")
    nsc = row.get("cross_correlation_nsc", "NA")
    rsc = row.get("cross_correlation_rsc", "NA")
    if cc_status == "reported":
        _parse_decimal(nsc, signed=True)
        _parse_decimal(rsc, signed=True)
    elif nsc != "NA" or rsc != "NA":
        raise EncodeQcSummaryError("cross-correlation values are inconsistent")
    revision = row.get("scientific_policy_revision", "")
    if not revision.isdigit() or int(revision) <= 0:
        raise EncodeQcSummaryError("scientific policy revision is invalid")


def _validate_mnase_v2(row: dict[str, str]) -> None:
    if row.get("schema_version") != "2":
        raise EncodeQcSummaryError("source schema version is unsupported")
    if row.get("analysis_input_stage") != "analysis_bam":
        raise EncodeQcSummaryError("analysis input stage is invalid")
    if row.get("analysis_status") not in _MNASE_ANALYSIS_STATUSES:
        raise EncodeQcSummaryError("MNase analysis status is invalid")
    if row.get("window_overlap") not in ("true", "false"):
        raise EncodeQcSummaryError("MNase window overlap flag is invalid")
    analysis_status = row.get("analysis_status")
    count_columns = (
        "analysis_total_templates",
        "analysis_valid_templates",
        "analysis_invalid_templates",
        "analysis_classified_templates",
        "analysis_unclassified_templates",
        "sub_fragments",
        "mono_fragments",
        "di_fragments",
    )
    if analysis_status == "parse_failed":
        if any(row.get(column) != "NA" for column in count_columns):
            raise EncodeQcSummaryError("MNase parse failure must use NA counts")
    elif any(row.get(column) == "NA" for column in count_columns):
        raise EncodeQcSummaryError("MNase analysis counts are incomplete")
    for column in (
        "analysis_total_templates",
        "analysis_valid_templates",
        "analysis_invalid_templates",
        "analysis_classified_templates",
        "analysis_unclassified_templates",
        "sub_fragments",
        "mono_fragments",
        "di_fragments",
    ):
        if row.get(column) != "NA":
            _parse_decimal(row.get(column, ""), integral=True)


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


def _policy_candidates(source, content: bytes):
    metadata = source.metadata
    sample_id = _required_metadata(metadata, "sample_id")
    assay = _required_metadata(metadata, "assay")
    scope = _required_metadata(metadata, "scope")
    experiment_id = _optional_metadata(metadata, "experiment_id")
    try:
        policy = parse_analysis_policy_bytes(content)
    except AnalysisPolicyError as exc:
        raise EncodeQcSummaryError("analysis policy is invalid") from exc
    if scope != "sample" or policy["sample"] != sample_id or policy["assay"] != assay:
        raise EncodeQcSummaryError("analysis policy identity is inconsistent")
    candidates = []
    for key, metric_key, display_name in _POLICY_STAGE_METRICS:
        value = policy["stage_counts"].get(key)
        if value is None:
            continue
        candidates.append(
            ExtractedQcMetricCandidate(
                metric_key=metric_key,
                display_name=display_name,
                value=Decimal(value),
                unit="count",
                scope="sample",
                sample_id=sample_id,
                experiment_id=experiment_id,
                assay=assay,
                qc_flag=None,
                source_artifact_id=source.artifact_id,
            )
        )
    return tuple(candidates)


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
    value = _parse_decimal(
        raw,
        integral=spec.integral,
        signed=spec.signed,
    )
    if not spec.signed and value < 0:
        raise EncodeQcSummaryError("QC numeric value is invalid")
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


def _parse_decimal(raw, *, integral=False, signed=False):
    if not isinstance(raw, str) or not raw:
        raise EncodeQcSummaryError("QC numeric value is invalid")
    pattern = _SIGNED_DECIMAL_PATTERN if signed else _DECIMAL_PATTERN
    if pattern.fullmatch(raw) is None:
        raise EncodeQcSummaryError("QC numeric value is invalid")
    try:
        value = Decimal(raw)
    except InvalidOperation as exc:
        raise EncodeQcSummaryError("QC numeric value is invalid") from exc
    if not value.is_finite():
        raise EncodeQcSummaryError("QC numeric value is invalid")
    if integral and value != value.to_integral_value():
        raise EncodeQcSummaryError("QC count is not integral")
    return value


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
