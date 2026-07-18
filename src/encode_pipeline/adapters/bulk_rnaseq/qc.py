"""Bounded machine-readable QC extraction for nf-core/rnaseq 3.26.0."""

from __future__ import annotations

from collections.abc import Mapping
from decimal import Context, Decimal, InvalidOperation, ROUND_HALF_EVEN, localcontext
from io import BytesIO
import json
import math
from pathlib import PurePosixPath
import re
import stat
from typing import Any
from zipfile import BadZipFile, ZipFile

from encode_pipeline.adapters.bulk_rnaseq.validation import (
    validate_bulk_rnaseq_inputs,
)
from encode_pipeline.adapters.bulk_rnaseq.results_contract import (
    effective_downstream_layout,
    effective_rseqc_modules,
    load_bulk_rnaseq_results_contract,
    trimmed_fastqc_enabled,
)
from encode_pipeline.adapters.bulk_rnaseq.status_evidence import (
    StatusEvidenceContractError,
    StatusEvidenceError,
    parse_cutadapt_processed_reads,
    parse_cutadapt_rows,
    parse_fastp_retained_reads,
    parse_fastp_summary,
    parse_fastqc_total_sequences,
    parse_star_log_final,
    parse_star_uniquely_mapped_percent,
    parse_status_table,
    reconcile_sample_status,
)
from encode_pipeline.platform.adapters import (
    ExtractedQcMetricCandidate,
    MAX_SAMPLE_ROWS,
    QcSourceDocument,
    WorkflowInputs,
)
from encode_pipeline.platform.results import Issue, Result


_ASSAY = "bulk-rnaseq"
_FASTQC_VERSION = "0.12.1"
_SALMON_VERSION = "1.10.3"
_MAX_SOURCE_FILES = MAX_SAMPLE_ROWS * 16 + 16
_MAX_SOURCE_BYTES = 16 * 1024 * 1024
_MAX_TOTAL_SOURCE_BYTES = 256 * 1024 * 1024
_MAX_TEXT_LINES = 20_000
_MAX_LINE_LENGTH = 8_192
_MAX_JSON_DEPTH = 16
_MAX_JSON_NODES = 50_000
_MAX_JSON_STRING = 8_192
_MAX_ZIP_ENTRIES = 256
_MAX_ZIP_ENTRY_BYTES = 16 * 1024 * 1024
_MAX_ZIP_TOTAL_BYTES = 64 * 1024 * 1024
_MAX_ZIP_RATIO = Decimal("200")
_MAX_METRICS = MAX_SAMPLE_ROWS * 128
_RATIO_QUANTUM = Decimal("0.000000000001")
_QC_DECIMAL_CONTEXT = Context(prec=64, rounding=ROUND_HALF_EVEN, Emin=-999, Emax=999)
_MIN_DECIMAL_ADJUSTED = -100
_MAX_DECIMAL_ADJUSTED = 25
_SAFE_ARTIFACT_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,255}$")
_SAFE_SAMPLE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_FEATURECOUNTS_STATUSES = frozenset(
    {
        "Assigned",
        "Unassigned_Unmapped",
        "Unassigned_Read_Type",
        "Unassigned_Singleton",
        "Unassigned_MappingQuality",
        "Unassigned_Chimera",
        "Unassigned_FragmentLength",
        "Unassigned_Duplicate",
        "Unassigned_MultiMapping",
        "Unassigned_Secondary",
        "Unassigned_NonSplit",
        "Unassigned_NoFeatures",
        "Unassigned_Overlapping_Length",
        "Unassigned_Ambiguity",
    }
)
_FASTQC_ROOT_FILES = frozenset(
    {"fastqc_data.txt", "fastqc_report.html", "fastqc.fo", "summary.txt"}
)
_FASTQC_ICON_FILES = frozenset(
    {"error.png", "fastqc_icon.png", "tick.png", "warning.png"}
)
_FASTQC_IMAGE_FILES = frozenset(
    {
        "adapter_content.png",
        "duplication_levels.png",
        "kmer_profiles.png",
        "overrepresented_sequences.png",
        "per_base_n_content.png",
        "per_base_quality.png",
        "per_base_sequence_content.png",
        "per_sequence_gc_content.png",
        "per_sequence_quality.png",
        "per_tile_quality.png",
        "sequence_duplication_levels.png",
        "sequence_length_distribution.png",
    }
)
_FASTQC_STATUS_MODULES = frozenset(
    {
        "Basic Statistics",
        "Per base sequence quality",
        "Per sequence quality scores",
        "Per sequence GC content",
        "Adapter Content",
    }
)


BULK_RNASEQ_QC_SOURCE_TYPES = tuple(
    sorted(
        {
            *(
                f"bulk_rnaseq.fastqc.{stage}.{role}.zip"
                for stage in ("raw", "trimmed", "filtered")
                for role in ("single", "read1", "read2")
            ),
            "bulk_rnaseq.featurecounts.summary",
            "bulk_rnaseq.multiqc.cutadapt",
            "bulk_rnaseq.multiqc.fail_mapped_samples",
            "bulk_rnaseq.multiqc.fail_trimmed_samples",
            "bulk_rnaseq.multiqc.picard_dups",
            "bulk_rnaseq.rseqc.bam_stat",
            "bulk_rnaseq.rseqc.infer_experiment",
            "bulk_rnaseq.rseqc.read_distribution",
            "bulk_rnaseq.rseqc.tin_summary",
            "bulk_rnaseq.salmon.meta_info",
            "bulk_rnaseq.star.log_final",
            "bulk_rnaseq.trim.fastp.json",
        }
    )
)


class _QcError(ValueError):
    def __init__(self, reason_code: str) -> None:
        super().__init__(reason_code)
        self.reason_code = reason_code


def extract_bulk_rnaseq_qc_metrics(
    inputs: WorkflowInputs,
    sources: tuple[QcSourceDocument, ...],
) -> Result[tuple[ExtractedQcMetricCandidate, ...]]:
    """Extract a closed numeric catalog without parsing report HTML."""
    with localcontext(_QC_DECIMAL_CONTEXT):
        return _extract_bulk_rnaseq_qc_metrics(inputs, sources)


def _extract_bulk_rnaseq_qc_metrics(
    inputs: WorkflowInputs,
    sources: tuple[QcSourceDocument, ...],
) -> Result[tuple[ExtractedQcMetricCandidate, ...]]:
    try:
        load_bulk_rnaseq_results_contract()
    except (OSError, TypeError, ValueError):
        return _failure("results_contract_invalid", path="results_contract")
    validated = validate_bulk_rnaseq_inputs(inputs)
    if validated.is_failure:
        return _failure("inputs_invalid", path="inputs")
    try:
        normalized = _normalized_mapping(validated.value)
        samples = _sample_contract(normalized)
        params = _parameter_contract(normalized)
        documents = _validate_sources(sources, samples)
        cutadapt_metrics: list[ExtractedQcMetricCandidate] = []
        trimgalore_input_evidence: Mapping[str, Decimal] = {}
        cutadapt_document = documents.get(("bulk_rnaseq.multiqc.cutadapt", ""))
        if cutadapt_document is not None:
            trimgalore_input_evidence = parse_cutadapt_processed_reads(
                cutadapt_document.content,
                expected_row_owners=_cutadapt_row_owners(samples, params),
            )
            cutadapt_metrics = _parse_multiqc_cutadapt(
                cutadapt_document,
                samples,
                params=params,
            )
        trimmed_failed, mapped_failed = _sample_status(
            documents,
            samples,
            params=params,
            trimgalore_input_evidence=trimgalore_input_evidence,
        )
        expected = _expected_sources(
            samples,
            params,
            trimmed_failed=trimmed_failed,
            mapped_failed=mapped_failed,
        )
        expected.update(
            coordinate
            for coordinate in documents
            if coordinate[0]
            in {
                "bulk_rnaseq.multiqc.fail_trimmed_samples",
                "bulk_rnaseq.multiqc.fail_mapped_samples",
            }
        )
        actual = set(documents)
        if actual != expected:
            raise _QcError("source_contract_invalid")

        metrics: list[ExtractedQcMetricCandidate] = []
        for coordinate in sorted(documents):
            output_type, sample_id = coordinate
            document = documents[coordinate]
            if output_type.startswith("bulk_rnaseq.fastqc."):
                metrics.extend(_parse_fastqc(document, sample_id, params=params))
            elif output_type == "bulk_rnaseq.multiqc.cutadapt":
                metrics.extend(cutadapt_metrics)
            elif output_type == "bulk_rnaseq.trim.fastp.json":
                metrics.extend(_parse_fastp(document, sample_id))
            elif output_type == "bulk_rnaseq.star.log_final":
                metrics.extend(_parse_star(document, sample_id))
            elif output_type == "bulk_rnaseq.salmon.meta_info":
                metrics.extend(_parse_salmon(document, sample_id))
            elif output_type == "bulk_rnaseq.featurecounts.summary":
                metrics.extend(_parse_featurecounts(document, sample_id, params=params))
            elif output_type == "bulk_rnaseq.multiqc.picard_dups":
                metrics.extend(
                    _parse_picard_multiqc(
                        document,
                        set(samples).difference(trimmed_failed, mapped_failed),
                    )
                )
            elif output_type in {
                "bulk_rnaseq.multiqc.fail_trimmed_samples",
                "bulk_rnaseq.multiqc.fail_mapped_samples",
            }:
                continue
            elif output_type == "bulk_rnaseq.rseqc.bam_stat":
                metrics.extend(
                    _parse_rseqc_bam_stat(
                        document,
                        sample_id,
                        layout=effective_downstream_layout(
                            samples[sample_id]["layout"], params
                        ),
                    )
                )
            elif output_type == "bulk_rnaseq.rseqc.infer_experiment":
                metrics.extend(
                    _parse_rseqc_infer(
                        document,
                        sample_id,
                        layout=effective_downstream_layout(
                            samples[sample_id]["layout"], params
                        ),
                    )
                )
            elif output_type == "bulk_rnaseq.rseqc.read_distribution":
                metrics.extend(_parse_rseqc_distribution(document, sample_id))
            elif output_type == "bulk_rnaseq.rseqc.tin_summary":
                metrics.extend(_parse_rseqc_tin(document, sample_id, params=params))
            else:  # pragma: no cover - protected by the closed source catalog
                raise _QcError("source_contract_invalid")

        _validate_trimmed_fastqc_counts(
            metrics,
            samples,
            params,
            trimmed_failed=trimmed_failed,
        )
        if len(metrics) > _MAX_METRICS:
            raise _QcError("metric_limit_exceeded")
        seen: set[tuple[str, str, str | None, str | None]] = set()
        for metric in metrics:
            metric_coordinate = (
                metric.metric_key,
                metric.scope,
                metric.sample_id,
                metric.experiment_id,
            )
            if metric_coordinate in seen:
                raise _QcError("metric_contract_invalid")
            seen.add(metric_coordinate)
        return Result.success(
            tuple(
                sorted(
                    metrics,
                    key=lambda metric: (
                        metric.metric_key,
                        metric.scope,
                        metric.sample_id or "",
                        metric.experiment_id or "",
                        metric.source_artifact_id,
                    ),
                )
            )
        )
    except _QcError as exc:
        return _failure(exc.reason_code)
    except (
        BadZipFile,
        InvalidOperation,
        KeyError,
        TypeError,
        UnicodeError,
        ValueError,
    ):
        return _failure("source_content_invalid")


def _normalized_mapping(value: object) -> Mapping[str, object]:
    if not isinstance(value, Mapping):
        raise _QcError("inputs_invalid")
    return value


def _sample_contract(
    normalized: Mapping[str, object],
) -> dict[str, dict[str, str]]:
    rows = normalized.get("samples")
    if not isinstance(rows, list):
        raise _QcError("inputs_invalid")
    samples: dict[str, dict[str, str]] = {}
    for row in rows:
        if not isinstance(row, Mapping):
            raise _QcError("inputs_invalid")
        sample = row.get("sample")
        layout = row.get("layout")
        strandedness = row.get("strandedness")
        if (
            not isinstance(sample, str)
            or _SAFE_SAMPLE.fullmatch(sample) is None
            or layout not in {"SE", "PE"}
            or strandedness not in {"auto", "forward", "reverse", "unstranded"}
        ):
            raise _QcError("inputs_invalid")
        previous = samples.setdefault(
            sample,
            {"layout": layout, "strandedness": strandedness},
        )
        if previous != {"layout": layout, "strandedness": strandedness}:
            raise _QcError("inputs_invalid")
    if not samples:
        raise _QcError("inputs_invalid")
    return samples


def _parameter_contract(normalized: Mapping[str, object]) -> dict[str, object]:
    params = normalized.get("nfcore_params")
    if not isinstance(params, Mapping):
        raise _QcError("inputs_invalid")
    return dict(params)


def _validate_sources(
    sources: tuple[QcSourceDocument, ...],
    samples: Mapping[str, Mapping[str, str]],
) -> dict[tuple[str, str], QcSourceDocument]:
    if not isinstance(sources, tuple) or len(sources) > _MAX_SOURCE_FILES:
        raise _QcError("source_contract_invalid")
    documents: dict[tuple[str, str], QcSourceDocument] = {}
    seen_ids: set[str] = set()
    seen_paths: set[str] = set()
    total_bytes = 0
    allowed = frozenset(BULK_RNASEQ_QC_SOURCE_TYPES)
    for document in sources:
        if not isinstance(document, QcSourceDocument):
            raise _QcError("source_contract_invalid")
        source = document.source
        if (
            source.output_type not in allowed
            or _SAFE_ARTIFACT_ID.fullmatch(source.artifact_id) is None
            or source.artifact_id in seen_ids
            or source.relative_path in seen_paths
            or not _safe_results_path(source.relative_path)
            or len(document.content) > _MAX_SOURCE_BYTES
        ):
            raise _QcError("source_contract_invalid")
        seen_ids.add(source.artifact_id)
        seen_paths.add(source.relative_path)
        total_bytes += len(document.content)
        if total_bytes > _MAX_TOTAL_SOURCE_BYTES:
            raise _QcError("source_contract_invalid")
        metadata = source.metadata
        if not isinstance(metadata, Mapping) or metadata.get("assay") != _ASSAY:
            raise _QcError("source_contract_invalid")
        if source.output_type in {
            "bulk_rnaseq.multiqc.cutadapt",
            "bulk_rnaseq.multiqc.picard_dups",
            "bulk_rnaseq.multiqc.fail_trimmed_samples",
            "bulk_rnaseq.multiqc.fail_mapped_samples",
        }:
            if metadata.get("scope") != "run" or metadata.get("sample_id") is not None:
                raise _QcError("source_contract_invalid")
            sample_id = ""
        else:
            raw_sample_id = metadata.get("sample_id")
            if (
                metadata.get("scope") != "sample"
                or not isinstance(raw_sample_id, str)
                or raw_sample_id not in samples
            ):
                raise _QcError("source_contract_invalid")
            sample_id = raw_sample_id
        coordinate = (source.output_type, sample_id)
        if coordinate in documents:
            raise _QcError("source_contract_invalid")
        documents[coordinate] = document
    return documents


def _safe_results_path(value: str) -> bool:
    if (
        not isinstance(value, str)
        or not value
        or len(value) > 2048
        or "\x00" in value
        or "\\" in value
        or value.startswith(("/", "~"))
    ):
        return False
    path = PurePosixPath(value)
    return (
        path.as_posix() == value
        and not path.is_absolute()
        and bool(path.parts)
        and path.parts[0] == "results"
        and len(path.parts) <= 32
        and all(part not in {"", ".", ".."} for part in value.split("/"))
    )


def _expected_sources(
    samples: Mapping[str, Mapping[str, str]],
    params: Mapping[str, object],
    *,
    trimmed_failed: frozenset[str],
    mapped_failed: frozenset[str],
) -> set[tuple[str, str]]:
    expected: set[tuple[str, str]] = set()
    skip_fastqc = _bool_param(params, "skip_fastqc")
    skip_trimming = _bool_param(params, "skip_trimming")
    remove_ribo = _bool_param(params, "remove_ribo_rna")
    for sample, values in samples.items():
        if sample not in trimmed_failed:
            expected.add(("bulk_rnaseq.star.log_final", sample))
            expected.add(("bulk_rnaseq.salmon.meta_info", sample))
        raw_roles = ("single",) if values["layout"] == "SE" else ("read1", "read2")
        downstream_layout = effective_downstream_layout(values["layout"], params)
        downstream_roles = (
            ("single",) if downstream_layout == "SE" else ("read1", "read2")
        )
        if not skip_fastqc:
            for role in raw_roles:
                expected.add((f"bulk_rnaseq.fastqc.raw.{role}.zip", sample))
        if trimmed_fastqc_enabled(params):
            for role in downstream_roles:
                if not (params.get("trimmer") == "fastp" and sample in trimmed_failed):
                    expected.add((f"bulk_rnaseq.fastqc.trimmed.{role}.zip", sample))
        if not skip_fastqc:
            for role in downstream_roles:
                if remove_ribo and sample not in trimmed_failed:
                    expected.add((f"bulk_rnaseq.fastqc.filtered.{role}.zip", sample))
        if not skip_trimming:
            trimmer = params.get("trimmer")
            if trimmer == "fastp":
                expected.add(("bulk_rnaseq.trim.fastp.json", sample))
            elif trimmer != "trimgalore":
                raise _QcError("source_contract_invalid")
        if (
            sample not in trimmed_failed
            and sample not in mapped_failed
            and not _bool_param(params, "skip_biotype_qc")
        ):
            expected.add(("bulk_rnaseq.featurecounts.summary", sample))
        if (
            sample not in trimmed_failed
            and sample not in mapped_failed
            and not _bool_param(params, "skip_rseqc")
        ):
            modules = frozenset(effective_rseqc_modules(params))
            for module, output_type in (
                ("bam_stat", "bulk_rnaseq.rseqc.bam_stat"),
                ("infer_experiment", "bulk_rnaseq.rseqc.infer_experiment"),
                ("read_distribution", "bulk_rnaseq.rseqc.read_distribution"),
                ("tin", "bulk_rnaseq.rseqc.tin_summary"),
            ):
                if module in modules:
                    expected.add((output_type, sample))
    if (
        not skip_trimming
        and params.get("trimmer") == "trimgalore"
        and not _bool_param(params, "skip_multiqc")
    ):
        expected.add(("bulk_rnaseq.multiqc.cutadapt", ""))
    if (
        not _bool_param(params, "skip_multiqc")
        and not _bool_param(params, "skip_markduplicates")
        and not _bool_param(params, "with_umi")
        and set(samples).difference(trimmed_failed, mapped_failed)
    ):
        expected.add(("bulk_rnaseq.multiqc.picard_dups", ""))
    return expected


def _sample_status(
    documents: Mapping[tuple[str, str], QcSourceDocument],
    samples: Mapping[str, Mapping[str, str]],
    *,
    params: Mapping[str, object],
    trimgalore_input_evidence: Mapping[str, Decimal],
) -> tuple[frozenset[str], frozenset[str]]:
    known_samples = frozenset(samples)
    trim_document = documents.get(("bulk_rnaseq.multiqc.fail_trimmed_samples", ""))
    mapped_document = documents.get(("bulk_rnaseq.multiqc.fail_mapped_samples", ""))
    if _bool_param(params, "skip_multiqc"):
        if trim_document is not None or mapped_document is not None:
            raise _QcError("source_contract_invalid")
        return frozenset(), frozenset()
    if _bool_param(params, "skip_trimming") and trim_document is not None:
        raise _QcError("source_contract_invalid")
    try:
        trim_table = parse_status_table(
            None if trim_document is None else trim_document.content,
            kind="trimmed",
            expected_header="Sample\tReads after trimming",
            known_samples=known_samples,
        )
        mapped_table = parse_status_table(
            None if mapped_document is None else mapped_document.content,
            kind="mapped",
            expected_header="Sample\tSTAR uniquely mapped reads (%)",
            known_samples=known_samples,
        )
    except StatusEvidenceContractError as exc:
        raise _QcError("source_contract_invalid") from exc
    except StatusEvidenceError as exc:
        raise _QcError("source_content_invalid") from exc

    trim_threshold = Decimal(_min_trimmed_reads(params))
    mapped_threshold = _min_mapped_percent(params)
    retained_evidence: dict[str, Decimal] = {}
    if not _bool_param(params, "skip_trimming"):
        for sample_id, sample in samples.items():
            if params.get("trimmer") == "fastp":
                document = documents.get(("bulk_rnaseq.trim.fastp.json", sample_id))
                if document is not None:
                    try:
                        retained_evidence[sample_id] = parse_fastp_retained_reads(
                            document.content
                        )
                    except StatusEvidenceError as exc:
                        raise _QcError("source_content_invalid") from exc
            elif params.get("trimmer") == "trimgalore":
                downstream_layout = effective_downstream_layout(
                    sample["layout"], params
                )
                roles = ("single",) if downstream_layout == "SE" else ("read1", "read2")
                counts: list[Decimal] = []
                for role in roles:
                    document = documents.get(
                        (f"bulk_rnaseq.fastqc.trimmed.{role}.zip", sample_id)
                    )
                    if document is None:
                        counts = []
                        break
                    root, filename, _path = _fastqc_identity(
                        sample_id,
                        "trimmed",
                        role,
                        params=params,
                    )
                    try:
                        counts.append(
                            parse_fastqc_total_sequences(
                                document.content,
                                expected_root=root,
                                expected_filename=filename,
                            )
                        )
                    except StatusEvidenceError as exc:
                        raise _QcError("source_content_invalid") from exc
                if counts:
                    if len(set(counts)) != 1:
                        raise _QcError("source_content_invalid")
                    retained_evidence[sample_id] = counts[0]
            elif params.get("trimmer") not in {"fastp", "trimgalore"}:
                raise _QcError("source_contract_invalid")

    mapped_evidence: dict[str, Decimal] = {}
    for sample_id in samples:
        document = documents.get(("bulk_rnaseq.star.log_final", sample_id))
        if document is None:
            continue
        try:
            mapped_evidence[sample_id] = parse_star_uniquely_mapped_percent(
                document.content
            )
        except StatusEvidenceError as exc:
            raise _QcError("source_content_invalid") from exc
    try:
        status = reconcile_sample_status(
            known_samples=known_samples,
            trimming_enabled=not _bool_param(params, "skip_trimming"),
            trimming_tool=str(params.get("trimmer")),
            trim_threshold=trim_threshold,
            mapped_threshold=mapped_threshold,
            trim_table=trim_table,
            mapped_table=mapped_table,
            retained_read_evidence=retained_evidence,
            trimgalore_input_read_evidence=trimgalore_input_evidence,
            uniquely_mapped_percent_evidence=mapped_evidence,
        )
    except StatusEvidenceContractError as exc:
        raise _QcError("source_contract_invalid") from exc
    except StatusEvidenceError as exc:
        raise _QcError("source_content_invalid") from exc
    return status.trimmed_failed, status.mapped_failed


def _cutadapt_row_owners(
    samples: Mapping[str, Mapping[str, str]],
    params: Mapping[str, object],
) -> Mapping[str, str]:
    owners: dict[str, str] = {}
    for sample_id, sample in samples.items():
        layout = effective_downstream_layout(sample["layout"], params)
        row_names = (
            (sample_id,) if layout == "SE" else (f"{sample_id}_1", f"{sample_id}_2")
        )
        for row_name in row_names:
            previous = owners.setdefault(row_name, sample_id)
            if previous != sample_id:
                raise _QcError("source_contract_invalid")
    return owners


def _bool_param(params: Mapping[str, object], name: str) -> bool:
    value = params.get(name)
    if not isinstance(value, bool):
        raise _QcError("source_contract_invalid")
    return value


def _min_trimmed_reads(params: Mapping[str, object]) -> int:
    value = params.get("min_trimmed_reads", 10_000)
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise _QcError("source_contract_invalid")
    return value


def _min_mapped_fraction(params: Mapping[str, object]) -> Decimal:
    value = _min_mapped_percent(params)
    return _fraction(value / Decimal(100))


def _min_mapped_percent(params: Mapping[str, object]) -> Decimal:
    value = _decimal(params.get("min_mapped_reads", 5))
    if value < 0 or value > 100:
        raise _QcError("source_contract_invalid")
    return value


def _parse_fastqc(
    document: QcSourceDocument,
    sample_id: str,
    *,
    params: Mapping[str, object],
) -> list[ExtractedQcMetricCandidate]:
    parts = document.source.output_type.split(".")
    if len(parts) != 5:
        raise _QcError("source_contract_invalid")
    stage, role = parts[2], parts[3]
    expected_root, expected_filename, expected_path = _fastqc_identity(
        sample_id,
        stage,
        role,
        params=params,
    )
    if document.source.relative_path != expected_path:
        raise _QcError("source_contract_invalid")
    try:
        with ZipFile(BytesIO(document.content)) as archive:
            if archive.comment:
                raise _QcError("source_content_invalid")
            infos = archive.infolist()
            if not 1 <= len(infos) <= _MAX_ZIP_ENTRIES:
                raise _QcError("source_content_invalid")
            total = 0
            summary_infos = []
            data_infos = []
            archive_members: set[str] = set()
            for info in infos:
                path = PurePosixPath(info.filename)
                mode = info.external_attr >> 16
                if (
                    not info.filename
                    or len(info.filename) > 512
                    or "\\" in info.filename
                    or path.is_absolute()
                    or len(path.parts) > 3
                    or any(part in {"", ".", ".."} for part in path.parts)
                    or info.flag_bits & 0x1
                    or info.comment
                    or info.extra
                    or stat.S_ISLNK(mode)
                    or info.file_size > _MAX_ZIP_ENTRY_BYTES
                    or info.filename in archive_members
                    or not _valid_fastqc_member(info, expected_root=expected_root)
                ):
                    raise _QcError("source_content_invalid")
                archive_members.add(info.filename)
                total += info.file_size
                if total > _MAX_ZIP_TOTAL_BYTES:
                    raise _QcError("source_content_invalid")
                compressed = max(info.compress_size, 1)
                if Decimal(info.file_size) / Decimal(compressed) > _MAX_ZIP_RATIO:
                    raise _QcError("source_content_invalid")
                if path.name == "summary.txt":
                    summary_infos.append(info)
                elif path.name == "fastqc_data.txt":
                    data_infos.append(info)
            if len(summary_infos) != 1 or len(data_infos) != 1:
                raise _QcError("source_content_invalid")
            if (
                summary_infos[0].filename.rsplit("/", 1)[0]
                != data_infos[0].filename.rsplit("/", 1)[0]
                or summary_infos[0].filename.rsplit("/", 1)[0] != expected_root
            ):
                raise _QcError("source_content_invalid")
            summary = archive.read(summary_infos[0])
            data = archive.read(data_infos[0])
    except (
        BadZipFile,
        KeyError,
        NotImplementedError,
        OSError,
        RuntimeError,
        ValueError,
    ) as exc:
        if isinstance(exc, _QcError):
            raise
        raise _QcError("source_content_invalid") from exc

    statuses = _fastqc_statuses(summary, expected_filename=expected_filename)
    basic, adapter_max, filename, data_statuses = _fastqc_data(data)
    if filename != expected_filename or data_statuses != statuses:
        raise _QcError("source_content_invalid")
    key = f"fastqc.{stage}.{role}"
    source_id = document.source.artifact_id
    metrics = [
        _metric(
            f"{key}.total_sequences",
            f"FastQC {stage} {role} total sequences",
            basic["Total Sequences"],
            "count",
            sample_id,
            source_id,
            qc_flag=statuses["Basic Statistics"],
        ),
        _metric(
            f"{key}.gc_fraction",
            f"FastQC {stage} {role} GC fraction",
            _percent(basic["%GC"]),
            "fraction",
            sample_id,
            source_id,
            qc_flag=statuses["Per sequence GC content"],
        ),
        _metric(
            f"{key}.per_base_sequence_quality_pass",
            f"FastQC {stage} {role} per-base quality pass indicator",
            Decimal(1)
            if statuses["Per base sequence quality"] == "pass"
            else Decimal(0),
            "fraction",
            sample_id,
            source_id,
            qc_flag=statuses["Per base sequence quality"],
        ),
        _metric(
            f"{key}.per_sequence_quality_pass",
            f"FastQC {stage} {role} per-sequence quality pass indicator",
            Decimal(1)
            if statuses["Per sequence quality scores"] == "pass"
            else Decimal(0),
            "fraction",
            sample_id,
            source_id,
            qc_flag=statuses["Per sequence quality scores"],
        ),
        _metric(
            f"{key}.adapter_content_max_fraction",
            f"FastQC {stage} {role} maximum adapter content fraction",
            _percent(adapter_max),
            "fraction",
            sample_id,
            source_id,
            qc_flag=statuses["Adapter Content"],
        ),
    ]
    if stage == "trimmed":
        metrics.append(
            _metric(
                f"trimming.{role}.retained_reads",
                f"Trimming {role.replace('read', 'read ')} retained reads",
                basic["Total Sequences"],
                "count",
                sample_id,
                source_id,
            )
        )
    return metrics


def _valid_fastqc_member(info: Any, *, expected_root: str) -> bool:
    path = PurePosixPath(info.filename)
    if not path.parts or path.parts[0] != expected_root:
        return False
    if info.is_dir():
        return path.parts in {
            (expected_root,),
            (expected_root, "Icons"),
            (expected_root, "Images"),
        }
    if len(path.parts) == 2:
        return path.parts[1] in _FASTQC_ROOT_FILES
    if len(path.parts) != 3:
        return False
    directory, name = path.parts[1:]
    if directory == "Icons":
        return name in _FASTQC_ICON_FILES
    if directory == "Images":
        return name in _FASTQC_IMAGE_FILES
    return False


def _fastqc_identity(
    sample_id: str,
    stage: str,
    role: str,
    *,
    params: Mapping[str, object],
) -> tuple[str, str, str]:
    role_number = {"single": "", "read1": "_1", "read2": "_2"}.get(role)
    if role_number is None:
        raise _QcError("source_contract_invalid")
    if stage == "raw":
        basename = f"{sample_id}_raw{role_number}"
        directory = "raw"
        extension = "fastq.gz"
    elif stage == "filtered":
        basename = f"{sample_id}_filtered{role_number}"
        directory = "filtered"
        extension = "fastq.gz"
    elif stage == "trimmed":
        directory = "trim"
        trimmer = params.get("trimmer")
        if trimmer == "fastp":
            basename = f"{sample_id}_trimmed{role_number}"
            extension = "fastq.gz"
        elif trimmer == "trimgalore":
            if role == "single":
                basename = f"{sample_id}_trimmed_trimmed"
            else:
                number = role_number.removeprefix("_")
                basename = f"{sample_id}_trimmed_{number}_val_{number}"
            extension = "fq.gz"
        else:
            raise _QcError("source_contract_invalid")
    else:
        raise _QcError("source_contract_invalid")
    root = f"{basename}_fastqc"
    return root, f"{basename}.{extension}", f"results/fastqc/{directory}/{root}.zip"


def _validate_trimmed_fastqc_counts(
    metrics: list[ExtractedQcMetricCandidate],
    samples: Mapping[str, Mapping[str, str]],
    params: Mapping[str, object],
    *,
    trimmed_failed: frozenset[str],
) -> None:
    if not trimmed_fastqc_enabled(params):
        return
    values = {
        (metric.sample_id, metric.metric_key): metric.value
        for metric in metrics
        if metric.metric_key.startswith("trimming.")
        and metric.metric_key.endswith(".retained_reads")
    }
    for sample_id, sample in samples.items():
        if params.get("trimmer") == "fastp" and sample_id in trimmed_failed:
            continue
        downstream_layout = effective_downstream_layout(sample["layout"], params)
        roles = ("single",) if downstream_layout == "SE" else ("read1", "read2")
        counts: list[Decimal] = []
        for role in roles:
            value = values.get((sample_id, f"trimming.{role}.retained_reads"))
            if value is None:
                raise _QcError("source_contract_invalid")
            counts.append(value)
        if len(counts) == 2 and counts[0] != counts[1]:
            raise _QcError("source_content_invalid")
        if params.get("trimmer") == "fastp":
            fastp_retained = values.get((sample_id, "trimming.retained_reads"))
            if fastp_retained is None or fastp_retained != sum(counts, Decimal(0)):
                raise _QcError("source_content_invalid")


def _fastqc_statuses(
    content: bytes,
    *,
    expected_filename: str,
) -> dict[str, str]:
    text = _decode_text(content)
    statuses: dict[str, str] = {}
    mapping = {"PASS": "pass", "WARN": "warning", "FAIL": "fail"}
    for line in text.splitlines():
        fields = line.split("\t")
        if len(fields) != 3:
            raise _QcError("source_content_invalid")
        raw_status, module, filename = fields
        if filename != expected_filename:
            raise _QcError("source_content_invalid")
        if module in _FASTQC_STATUS_MODULES:
            if module in statuses or raw_status not in mapping:
                raise _QcError("source_content_invalid")
            statuses[module] = mapping[raw_status]
    if set(statuses) != _FASTQC_STATUS_MODULES:
        raise _QcError("source_content_invalid")
    return statuses


def _fastqc_data(
    content: bytes,
) -> tuple[dict[str, Decimal], Decimal, str, dict[str, str]]:
    text = _decode_text(content)
    lines = text.splitlines()
    if not lines or lines[0] != f"##FastQC\t{_FASTQC_VERSION}":
        raise _QcError("source_content_invalid")
    basic: dict[str, Decimal] = {}
    filename: str | None = None
    adapter_max: Decimal | None = None
    module: str | None = None
    adapter_columns = 0
    statuses: dict[str, str] = {}
    status_mapping = {"pass": "pass", "warn": "warning", "fail": "fail"}
    for line in lines[1:]:
        if line.startswith(">>"):
            if line == ">>END_MODULE":
                module = None
                adapter_columns = 0
            else:
                fields = line[2:].split("\t")
                if len(fields) != 2 or fields[1] not in {"pass", "warn", "fail"}:
                    raise _QcError("source_content_invalid")
                module = fields[0]
                if module in _FASTQC_STATUS_MODULES:
                    if module in statuses:
                        raise _QcError("source_content_invalid")
                    statuses[module] = status_mapping[fields[1]]
            continue
        if module == "Basic Statistics":
            if line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) != 2:
                raise _QcError("source_content_invalid")
            if fields[0] == "Filename":
                if filename is not None:
                    raise _QcError("source_content_invalid")
                filename = fields[1]
            elif fields[0] in {"Total Sequences", "%GC"}:
                if fields[0] in basic:
                    raise _QcError("source_content_invalid")
                basic[fields[0]] = _decimal(fields[1])
        elif module == "Adapter Content":
            fields = line.split("\t")
            if line.startswith("#"):
                if len(fields) < 2:
                    raise _QcError("source_content_invalid")
                adapter_columns = len(fields) - 1
                continue
            if adapter_columns == 0 or len(fields) != adapter_columns + 1:
                raise _QcError("source_content_invalid")
            values = [_decimal(value) for value in fields[1:]]
            if any(value < 0 or value > 100 for value in values):
                raise _QcError("source_content_invalid")
            row_max = max(values)
            adapter_max = row_max if adapter_max is None else max(adapter_max, row_max)
    if (
        set(basic) != {"Total Sequences", "%GC"}
        or adapter_max is None
        or filename is None
        or set(statuses) != _FASTQC_STATUS_MODULES
    ):
        raise _QcError("source_content_invalid")
    return basic, adapter_max, filename, statuses


def _parse_multiqc_cutadapt(
    document: QcSourceDocument,
    samples: Mapping[str, Mapping[str, str]],
    *,
    params: Mapping[str, object],
) -> list[ExtractedQcMetricCandidate]:
    expected_rows: dict[str, tuple[str, str]] = {}
    for sample_id, sample in samples.items():
        downstream_layout = effective_downstream_layout(sample["layout"], params)
        roles = ("single",) if downstream_layout == "SE" else ("read1", "read2")
        for index, role in enumerate(roles, start=1):
            row_name = sample_id if role == "single" else f"{sample_id}_{index}"
            owner = (sample_id, role)
            previous = expected_rows.setdefault(row_name, owner)
            if previous != owner:
                raise _QcError("source_contract_invalid")
    try:
        rows = parse_cutadapt_rows(
            document.content,
            expected_row_owners={
                row_name: owner[0] for row_name, owner in expected_rows.items()
            },
        )
    except StatusEvidenceError as exc:
        raise _QcError("source_content_invalid") from exc

    metrics: list[ExtractedQcMetricCandidate] = []
    processed_by_sample: dict[str, list[int]] = {}
    for row in rows:
        sample_id, role = expected_rows[row.row_identity]
        reads_processed = row.reads_processed
        reads_with_adapters = row.reads_with_adapters
        bases_processed = row.bases_processed
        quality_trimmed_bases = row.quality_trimmed_bases
        bases_written = row.bases_written
        processed_by_sample.setdefault(sample_id, []).append(reads_processed)
        source_id = document.source.artifact_id
        prefix = f"trimming.{role}"
        display_role = role.replace("read", "read ")
        metrics.extend(
            (
                _metric(
                    f"{prefix}.input_reads",
                    f"Trim Galore {display_role} input reads",
                    Decimal(reads_processed),
                    "count",
                    sample_id,
                    source_id,
                ),
                _metric(
                    f"{prefix}.adapter_affected_reads",
                    f"Trim Galore {display_role} reads with adapters",
                    Decimal(reads_with_adapters),
                    "count",
                    sample_id,
                    source_id,
                ),
                _metric(
                    f"{prefix}.adapter_affected_fraction",
                    f"Trim Galore {display_role} adapter-affected fraction",
                    _divide(reads_with_adapters, reads_processed),
                    "fraction",
                    sample_id,
                    source_id,
                ),
                _metric(
                    f"{prefix}.input_bases",
                    f"Trim Galore {display_role} input bases",
                    Decimal(bases_processed),
                    "count",
                    sample_id,
                    source_id,
                ),
                _metric(
                    f"{prefix}.post_adapter_quality_bases",
                    f"Trim Galore {display_role} bases after adapter and quality trimming",
                    Decimal(bases_written),
                    "count",
                    sample_id,
                    source_id,
                ),
                _metric(
                    f"{prefix}.post_adapter_quality_base_fraction",
                    f"Trim Galore {display_role} base fraction after adapter and quality trimming",
                    _divide(bases_written, bases_processed),
                    "fraction",
                    sample_id,
                    source_id,
                ),
                _metric(
                    f"{prefix}.quality_trimmed_bases",
                    f"Trim Galore {display_role} quality-trimmed bases",
                    Decimal(quality_trimmed_bases),
                    "count",
                    sample_id,
                    source_id,
                ),
            )
        )
    for sample_id, counts in processed_by_sample.items():
        if effective_downstream_layout(
            samples[sample_id]["layout"], params
        ) == "PE" and (len(counts) != 2 or counts[0] != counts[1]):
            raise _QcError("source_content_invalid")
    return metrics


def _parse_fastp(
    document: QcSourceDocument,
    sample_id: str,
) -> list[ExtractedQcMetricCandidate]:
    try:
        evidence = parse_fastp_summary(document.content)
    except StatusEvidenceError as exc:
        raise _QcError("source_content_invalid") from exc
    values = {
        "input_reads": evidence.input_reads,
        "retained_reads": evidence.retained_reads,
        "input_bases": evidence.input_bases,
        "retained_bases": evidence.retained_bases,
    }
    return _trimming_metrics(
        sample_id,
        document.source.artifact_id,
        values,
        key_prefix="trimming",
        display_role="sample",
    )


def _trimming_metrics(
    sample_id: str,
    source_id: str,
    values: Mapping[str, int],
    *,
    key_prefix: str,
    display_role: str,
) -> list[ExtractedQcMetricCandidate]:
    input_reads = values["input_reads"]
    retained_reads = values["retained_reads"]
    input_bases = values["input_bases"]
    retained_bases = values["retained_bases"]
    if (
        input_reads <= 0
        or input_bases <= 0
        or retained_reads > input_reads
        or retained_bases > input_bases
    ):
        raise _QcError("metric_value_invalid")
    return [
        _metric(
            f"{key_prefix}.input_reads",
            f"Trimming {display_role} input reads",
            Decimal(input_reads),
            "count",
            sample_id,
            source_id,
        ),
        _metric(
            f"{key_prefix}.retained_reads",
            f"Trimming {display_role} retained reads",
            Decimal(retained_reads),
            "count",
            sample_id,
            source_id,
        ),
        _metric(
            f"{key_prefix}.read_retained_fraction",
            f"Trimming {display_role} retained read fraction",
            _divide(retained_reads, input_reads),
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            f"{key_prefix}.input_bases",
            f"Trimming {display_role} input bases",
            Decimal(input_bases),
            "count",
            sample_id,
            source_id,
        ),
        _metric(
            f"{key_prefix}.retained_bases",
            f"Trimming {display_role} retained bases",
            Decimal(retained_bases),
            "count",
            sample_id,
            source_id,
        ),
        _metric(
            f"{key_prefix}.base_retained_fraction",
            f"Trimming {display_role} retained base fraction",
            _divide(retained_bases, input_bases),
            "fraction",
            sample_id,
            source_id,
        ),
    ]


def _parse_star(
    document: QcSourceDocument,
    sample_id: str,
) -> list[ExtractedQcMetricCandidate]:
    try:
        evidence = parse_star_log_final(document.content)
    except StatusEvidenceError as exc:
        raise _QcError("source_content_invalid") from exc
    input_reads = evidence.input_templates
    counts = {
        "unique": evidence.uniquely_mapped_templates,
        "multi": evidence.accepted_multimapped_templates,
        "too_many": evidence.too_many_loci_templates,
        "too_many_mismatches": evidence.unmapped_too_many_mismatches_templates,
        "too_short": evidence.unmapped_too_short_templates,
        "other": evidence.unmapped_other_templates,
    }
    exact_fractions = {
        name: _divide(count, input_reads) for name, count in counts.items()
    }
    unique = exact_fractions["unique"]
    multi = exact_fractions["multi"]
    too_many = exact_fractions["too_many"]
    unmapped_components = {
        name: exact_fractions[name]
        for name in ("too_many_mismatches", "too_short", "other")
    }
    pure_unmapped = _divide(
        counts["too_many_mismatches"] + counts["too_short"] + counts["other"],
        input_reads,
    )
    unaccepted = _divide(
        counts["too_many"]
        + counts["too_many_mismatches"]
        + counts["too_short"]
        + counts["other"],
        input_reads,
    )
    source_id = document.source.artifact_id
    metrics = [
        _metric(
            "star.input_templates",
            "STAR input templates",
            Decimal(input_reads),
            "count",
            sample_id,
            source_id,
        ),
        _metric(
            "star.uniquely_mapped_template_fraction",
            "STAR uniquely mapped template fraction",
            unique,
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "star.accepted_multimapped_template_fraction",
            "STAR accepted multimapped template fraction",
            multi,
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "star.too_many_loci_template_fraction",
            "STAR too-many-loci template fraction",
            too_many,
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "star.pure_unmapped_template_fraction",
            "STAR pure-unmapped template fraction",
            pure_unmapped,
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "star.unaccepted_template_fraction",
            "STAR unaccepted (too-many-loci or unmapped) template fraction",
            unaccepted,
            "fraction",
            sample_id,
            source_id,
        ),
    ]
    for name, value in unmapped_components.items():
        metrics.append(
            _metric(
                f"star.unmapped_{name}_template_fraction",
                f"STAR unmapped {name.replace('_', ' ')} template fraction",
                value,
                "fraction",
                sample_id,
                source_id,
            )
        )
    return metrics


def _parse_salmon(
    document: QcSourceDocument,
    sample_id: str,
) -> list[ExtractedQcMetricCandidate]:
    payload = _strict_json(document.content)
    if (
        not isinstance(payload, Mapping)
        or payload.get("salmon_version") != _SALMON_VERSION
        or payload.get("mapping_type") != "alignment"
        or _nonnegative_int(payload.get("num_libraries")) != 1
    ):
        raise _QcError("source_content_invalid")
    processed = _nonnegative_int(payload.get("num_processed"))
    mapped = _nonnegative_int(payload.get("num_mapped"))
    percent = _percent(payload.get("percent_mapped"))
    if processed <= 0 or mapped > processed:
        raise _QcError("metric_value_invalid")
    exact_mapping_fraction = _divide(mapped, processed)
    if abs(exact_mapping_fraction - percent) > _RATIO_QUANTUM:
        raise _QcError("source_content_invalid")
    source_id = document.source.artifact_id
    return [
        _metric(
            "salmon.processed_fragments",
            "Salmon processed fragments",
            Decimal(processed),
            "count",
            sample_id,
            source_id,
        ),
        _metric(
            "salmon.mapped_fragments",
            "Salmon mapped fragments",
            Decimal(mapped),
            "count",
            sample_id,
            source_id,
        ),
        _metric(
            "salmon.mapping_fraction",
            "Salmon mapped fragment fraction",
            exact_mapping_fraction,
            "fraction",
            sample_id,
            source_id,
        ),
    ]


def _parse_featurecounts(
    document: QcSourceDocument,
    sample_id: str,
    *,
    params: Mapping[str, object],
) -> list[ExtractedQcMetricCandidate]:
    rows = _tsv_rows(document.content)
    if len(rows) < 2 or len(rows[0]) != 2 or rows[0][0] != "Status":
        raise _QcError("source_content_invalid")
    bam_name = rows[0][1]
    expected = _final_bam_name(sample_id, params)
    if bam_name != expected:
        raise _QcError("source_content_invalid")
    counts: dict[str, int] = {}
    for row in rows[1:]:
        if len(row) != 2 or row[0] in counts or row[0] not in _FEATURECOUNTS_STATUSES:
            raise _QcError("source_content_invalid")
        counts[row[0]] = _nonnegative_int(row[1])
    if set(counts) != _FEATURECOUNTS_STATUSES:
        raise _QcError("source_content_invalid")
    total = sum(counts.values())
    if total <= 0:
        raise _QcError("metric_value_invalid")
    source_id = document.source.artifact_id
    return [
        _metric(
            "featurecounts.assigned_reads",
            "featureCounts assigned reads",
            Decimal(counts["Assigned"]),
            "count",
            sample_id,
            source_id,
        ),
        _metric(
            "featurecounts.classified_reads",
            "featureCounts classified reads",
            Decimal(total),
            "count",
            sample_id,
            source_id,
        ),
        _metric(
            "featurecounts.assigned_read_fraction",
            "featureCounts assigned read fraction",
            _divide(counts["Assigned"], total),
            "fraction",
            sample_id,
            source_id,
        ),
    ]


def _final_bam_name(sample_id: str, params: Mapping[str, object]) -> str:
    if _bool_param(params, "with_umi"):
        return f"{sample_id}.umi_dedup.sorted.bam"
    if not _bool_param(params, "skip_markduplicates"):
        return f"{sample_id}.markdup.sorted.bam"
    return f"{sample_id}.sorted.bam"


def _parse_picard_multiqc(
    document: QcSourceDocument,
    samples: set[str],
) -> list[ExtractedQcMetricCandidate]:
    rows = _tsv_rows(document.content)
    if len(rows) < 2:
        raise _QcError("source_content_invalid")
    header = rows[0]
    expected_header = [
        "Sample",
        "LIBRARY",
        "UNPAIRED_READS_EXAMINED",
        "READ_PAIRS_EXAMINED",
        "SECONDARY_OR_SUPPLEMENTARY_RDS",
        "UNMAPPED_READS",
        "UNPAIRED_READ_DUPLICATES",
        "READ_PAIR_DUPLICATES",
        "READ_PAIR_OPTICAL_DUPLICATES",
        "PERCENT_DUPLICATION",
        "ESTIMATED_LIBRARY_SIZE",
    ]
    if header != expected_header:
        raise _QcError("source_content_invalid")
    indexes = {name: header.index(name) for name in header}
    metrics: list[ExtractedQcMetricCandidate] = []
    seen: set[str] = set()
    for row in rows[1:]:
        if len(row) != len(header):
            raise _QcError("source_content_invalid")
        sample = row[indexes["Sample"]]
        if sample not in samples or sample in seen:
            raise _QcError("source_contract_invalid")
        seen.add(sample)
        if not row[indexes["LIBRARY"]] or len(row[indexes["LIBRARY"]]) > 256:
            raise _QcError("source_content_invalid")
        counts = {
            name: _nonnegative_integral_decimal(row[indexes[name]])
            for name in expected_header[2:-2]
        }
        unpaired_examined = counts["UNPAIRED_READS_EXAMINED"]
        pairs_examined = counts["READ_PAIRS_EXAMINED"]
        unpaired_duplicates = counts["UNPAIRED_READ_DUPLICATES"]
        pair_duplicates = counts["READ_PAIR_DUPLICATES"]
        optical_duplicates = counts["READ_PAIR_OPTICAL_DUPLICATES"]
        denominator = unpaired_examined + 2 * pairs_examined
        if (
            denominator <= 0
            or unpaired_duplicates > unpaired_examined
            or pair_duplicates > pairs_examined
            or optical_duplicates > pair_duplicates
        ):
            raise _QcError("metric_value_invalid")
        duplication = _fraction(_decimal(row[indexes["PERCENT_DUPLICATION"]]))
        expected_duplication = _fraction(
            Decimal(unpaired_duplicates + 2 * pair_duplicates) / Decimal(denominator)
        )
        # Picard 3.4.0 serializes doubles to six fractional digits.
        if abs(duplication - expected_duplication) > Decimal("0.0000005"):
            raise _QcError("source_content_invalid")
        metrics.append(
            _metric(
                "picard.duplication_fraction",
                "Picard duplicate read fraction",
                duplication,
                "fraction",
                sample,
                document.source.artifact_id,
            )
        )
        if "ESTIMATED_LIBRARY_SIZE" in indexes and row[
            indexes["ESTIMATED_LIBRARY_SIZE"]
        ] not in {"", "-"}:
            metrics.append(
                _metric(
                    "picard.estimated_library_size",
                    "Picard estimated library size",
                    Decimal(
                        _nonnegative_integral_decimal(
                            row[indexes["ESTIMATED_LIBRARY_SIZE"]]
                        )
                    ),
                    "count",
                    sample,
                    document.source.artifact_id,
                )
            )
    if seen != samples:
        raise _QcError("source_contract_invalid")
    return metrics


def _parse_rseqc_bam_stat(
    document: QcSourceDocument,
    sample_id: str,
    *,
    layout: str,
) -> list[ExtractedQcMetricCandidate]:
    values = _rseqc_bam_stat_counts(document.content, layout=layout)
    total = values["Total records"]
    unique = values["mapq >= mapq_cut (unique)"]
    nonunique = values["mapq < mapq_cut (non-unique)"]
    unmapped = values["Unmapped reads"]
    classified = (
        values["QC failed"]
        + values["Optical/PCR duplicate"]
        + values["Non primary hits"]
        + unmapped
        + nonunique
        + unique
    )
    if total <= 0 or classified != total:
        raise _QcError("metric_value_invalid")
    source_id = document.source.artifact_id
    metrics = [
        _metric(
            "rseqc.bam_stat.total_records",
            "RSeQC total alignment records",
            Decimal(total),
            "count",
            sample_id,
            source_id,
        ),
        _metric(
            "rseqc.bam_stat.unmapped_fraction",
            "RSeQC unmapped record fraction",
            _divide(unmapped, total),
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "rseqc.bam_stat.unique_fraction",
            "RSeQC unique mapped record fraction",
            _divide(unique, total),
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "rseqc.bam_stat.nonunique_fraction",
            "RSeQC non-unique mapped record fraction",
            _divide(nonunique, total),
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "rseqc.bam_stat.nonprimary_fraction",
            "RSeQC non-primary alignment record fraction",
            _divide(values["Non primary hits"], total),
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "rseqc.bam_stat.qc_failed_fraction",
            "RSeQC QC-failed record fraction",
            _divide(values["QC failed"], total),
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "rseqc.bam_stat.duplicate_fraction",
            "RSeQC optical or PCR duplicate fraction",
            _divide(values["Optical/PCR duplicate"], total),
            "fraction",
            sample_id,
            source_id,
        ),
    ]
    if layout == "PE":
        proper = values.get("Reads mapped in proper pairs")
        if proper is None:
            raise _QcError("source_content_invalid")
        metrics.append(
            _metric(
                "rseqc.bam_stat.proper_pair_fraction",
                "RSeQC proper-pair record fraction",
                _divide(proper, total),
                "fraction",
                sample_id,
                source_id,
            )
        )
    return metrics


def _rseqc_bam_stat_counts(content: bytes, *, layout: str) -> dict[str, int]:
    text = _decode_text(content)
    patterns = {
        "Total records": r"^Total records:\s+([0-9]+)$",
        "QC failed": r"^QC failed:\s+([0-9]+)$",
        "Optical/PCR duplicate": r"^Optical/PCR duplicate:\s+([0-9]+)$",
        "Non primary hits": r"^Non primary hits\s+([0-9]+)$",
        "Unmapped reads": r"^Unmapped reads:\s+([0-9]+)$",
        "mapq < mapq_cut (non-unique)": (
            r"^mapq < mapq_cut \(non-unique\):\s+([0-9]+)$"
        ),
        "mapq >= mapq_cut (unique)": (r"^mapq >= mapq_cut \(unique\):\s+([0-9]+)$"),
        "Reads mapped in proper pairs": (r"^Reads mapped in proper pairs:\s+([0-9]+)$"),
        "Read-1": r"^Read-1:\s+([0-9]+)$",
        "Read-2": r"^Read-2:\s+([0-9]+)$",
        "Reads map to plus": r"^Reads map to '\+':\s+([0-9]+)$",
        "Reads map to minus": r"^Reads map to '-':\s+([0-9]+)$",
        "Non-splice reads": r"^Non-splice reads:\s+([0-9]+)$",
        "Splice reads": r"^Splice reads:\s+([0-9]+)$",
        "Proper-paired reads map to different chrom": (
            r"^Proper-paired reads map to different chrom:\s+([0-9]+)$"
        ),
    }
    result: dict[str, int] = {}
    for line in text.splitlines():
        if not line.strip():
            continue
        if line.startswith("#"):
            if len(line) > 512:
                raise _QcError("source_content_invalid")
            continue
        matches = [
            (name, match)
            for name, pattern in patterns.items()
            if (match := re.fullmatch(pattern, line)) is not None
        ]
        if len(matches) != 1:
            raise _QcError("source_content_invalid")
        key, match = matches[0]
        if key in result:
            raise _QcError("source_content_invalid")
        result[key] = _nonnegative_int(match.group(1))
    paired_only = {
        "Read-1",
        "Read-2",
        "Reads mapped in proper pairs",
        "Proper-paired reads map to different chrom",
    }
    if set(result) != set(patterns):
        raise _QcError("source_content_invalid")
    if layout == "SE" and any(result[key] != 0 for key in paired_only):
        raise _QcError("source_content_invalid")
    total = result["Total records"]
    if any(value > total for key, value in result.items() if key != "Total records"):
        raise _QcError("metric_value_invalid")
    unique = result["mapq >= mapq_cut (unique)"]
    for key in (
        "Read-1",
        "Read-2",
        "Reads map to plus",
        "Reads map to minus",
        "Non-splice reads",
        "Splice reads",
        "Reads mapped in proper pairs",
        "Proper-paired reads map to different chrom",
    ):
        if key in result and result[key] > unique:
            raise _QcError("metric_value_invalid")
    if (
        result["Reads map to plus"] + result["Reads map to minus"] != unique
        or result["Non-splice reads"] + result["Splice reads"] != unique
    ):
        raise _QcError("source_content_invalid")
    if layout == "PE" and result["Read-1"] + result["Read-2"] != unique:
        raise _QcError("source_content_invalid")
    if (
        result["Proper-paired reads map to different chrom"]
        > result["Reads mapped in proper pairs"]
    ):
        raise _QcError("source_content_invalid")
    return result


def _parse_rseqc_infer(
    document: QcSourceDocument,
    sample_id: str,
    *,
    layout: str,
) -> list[ExtractedQcMetricCandidate]:
    text = _decode_text(document.content)
    expected_banner = (
        "This is SingleEnd Data" if layout == "SE" else "This is PairEnd Data"
    )
    lines = [line for line in text.splitlines() if line.strip()]
    if len(lines) != 4 or lines[0].strip() != expected_banner:
        raise _QcError("source_content_invalid")
    values: dict[str, Decimal] = {}
    for line in lines[1:]:
        fields = line.rsplit(":", 1)
        if len(fields) != 2 or fields[0] in values:
            raise _QcError("source_content_invalid")
        values[fields[0]] = _fraction(_decimal(fields[1].strip()))
    failed_key = "Fraction of reads failed to determine"
    orientation_keys = (
        (
            'Fraction of reads explained by "++,--"',
            'Fraction of reads explained by "+-,-+"',
        )
        if layout == "SE"
        else (
            'Fraction of reads explained by "1++,1--,2+-,2-+"',
            'Fraction of reads explained by "1+-,1-+,2++,2--"',
        )
    )
    if set(values) != {failed_key, *orientation_keys}:
        raise _QcError("source_content_invalid")
    total = values[failed_key] + sum(
        (values[key] for key in orientation_keys), Decimal(0)
    )
    # RSeQC 5.0.4 prints each component independently to four decimals.
    if abs(total - Decimal(1)) > Decimal("0.0001"):
        raise _QcError("metric_value_invalid")
    source_id = document.source.artifact_id
    return [
        _metric(
            "rseqc.infer_experiment.undetermined_fraction",
            "RSeQC undetermined strandedness fraction",
            values[failed_key],
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "rseqc.infer_experiment.orientation_a_fraction",
            "RSeQC first orientation fraction",
            values[orientation_keys[0]],
            "fraction",
            sample_id,
            source_id,
        ),
        _metric(
            "rseqc.infer_experiment.orientation_b_fraction",
            "RSeQC second orientation fraction",
            values[orientation_keys[1]],
            "fraction",
            sample_id,
            source_id,
        ),
    ]


def _parse_rseqc_distribution(
    document: QcSourceDocument,
    sample_id: str,
) -> list[ExtractedQcMetricCandidate]:
    lines = [line for line in _decode_text(document.content).splitlines() if line]
    if len(lines) != 16:
        raise _QcError("source_content_invalid")
    preamble: dict[str, int] = {}
    for line, expected_key in zip(
        lines[:3],
        ("Total Reads", "Total Tags", "Total Assigned Tags"),
        strict=True,
    ):
        fields = line.split()
        if len(fields) < 2 or " ".join(fields[:-1]) != expected_key:
            raise _QcError("source_content_invalid")
        preamble[expected_key] = _nonnegative_int(fields[-1])
    if set(lines[3]) != {"="} or set(lines[-1]) != {"="}:
        raise _QcError("source_content_invalid")
    if lines[4].split() != ["Group", "Total_bases", "Tag_count", "Tags/Kb"]:
        raise _QcError("source_content_invalid")
    total_reads = preamble["Total Reads"]
    total_tags = preamble["Total Tags"]
    total_assigned = preamble["Total Assigned Tags"]
    if (
        total_reads <= 0
        or total_tags <= 0
        or total_assigned > total_tags
        or total_assigned <= 0
    ):
        raise _QcError("metric_value_invalid")
    names = {
        "CDS_Exons": "cds_exon",
        "5'UTR_Exons": "five_prime_utr",
        "3'UTR_Exons": "three_prime_utr",
        "Introns": "intron",
    }
    all_groups = (
        *names,
        "TSS_up_1kb",
        "TSS_up_5kb",
        "TSS_up_10kb",
        "TES_down_1kb",
        "TES_down_5kb",
        "TES_down_10kb",
    )
    metrics: list[ExtractedQcMetricCandidate] = []
    group_counts: dict[str, int] = {}
    source_id = document.source.artifact_id
    metrics.extend(
        (
            _metric(
                "rseqc.read_distribution.total_reads",
                "RSeQC read distribution total reads",
                Decimal(total_reads),
                "count",
                sample_id,
                source_id,
            ),
            _metric(
                "rseqc.read_distribution.total_tags",
                "RSeQC read distribution total tags",
                Decimal(total_tags),
                "count",
                sample_id,
                source_id,
            ),
            _metric(
                "rseqc.read_distribution.total_assigned_tags",
                "RSeQC read distribution assigned tags",
                Decimal(total_assigned),
                "count",
                sample_id,
                source_id,
            ),
        )
    )
    for line, expected_group in zip(lines[5:-1], all_groups, strict=True):
        row = line.split()
        if len(row) != 4 or row[0] != expected_group:
            raise _QcError("source_content_invalid")
        total_bases = _nonnegative_int(row[1])
        tag_count = _nonnegative_int(row[2])
        tags_per_kb = _nonnegative_decimal(row[3])
        if tag_count > total_assigned:
            raise _QcError("metric_value_invalid")
        group_counts[row[0]] = tag_count
        expected_rate = Decimal(tag_count * 1000) / Decimal(total_bases + 1)
        if abs(tags_per_kb - expected_rate) > Decimal("0.0051"):
            raise _QcError("source_content_invalid")
        if row[0] not in names:
            continue
        token = names[row[0]]
        metrics.extend(
            (
                _metric(
                    f"rseqc.read_distribution.{token}_tags",
                    f"RSeQC {token.replace('_', ' ')} tags",
                    Decimal(tag_count),
                    "count",
                    sample_id,
                    source_id,
                ),
                _metric(
                    f"rseqc.read_distribution.{token}_assigned_fraction",
                    f"RSeQC {token.replace('_', ' ')} assigned tag fraction",
                    _divide(tag_count, total_assigned),
                    "fraction",
                    sample_id,
                    source_id,
                ),
            )
        )
    if (
        group_counts["CDS_Exons"]
        + group_counts["5'UTR_Exons"]
        + group_counts["3'UTR_Exons"]
        + group_counts["Introns"]
        + group_counts["TSS_up_10kb"]
        + group_counts["TES_down_10kb"]
        != total_assigned
        or not (
            group_counts["TSS_up_1kb"]
            <= group_counts["TSS_up_5kb"]
            <= group_counts["TSS_up_10kb"]
        )
        or not (
            group_counts["TES_down_1kb"]
            <= group_counts["TES_down_5kb"]
            <= group_counts["TES_down_10kb"]
        )
    ):
        raise _QcError("source_content_invalid")
    return metrics


def _parse_rseqc_tin(
    document: QcSourceDocument,
    sample_id: str,
    *,
    params: Mapping[str, object],
) -> list[ExtractedQcMetricCandidate]:
    rows = _tsv_rows(document.content)
    expected_header = ["Bam_file", "TIN(mean)", "TIN(median)", "TIN(stdev)"]
    if len(rows) != 2 or rows[0] != expected_header or len(rows[1]) != 4:
        raise _QcError("source_content_invalid")
    if rows[1][0] != _final_bam_name(sample_id, params):
        raise _QcError("source_content_invalid")
    mean = _canonical_decimal(_nonnegative_decimal(rows[1][1]))
    median = _canonical_decimal(_nonnegative_decimal(rows[1][2]))
    standard_deviation = _canonical_decimal(_nonnegative_decimal(rows[1][3]))
    if mean > 100 or median > 100 or standard_deviation > 50:
        raise _QcError("metric_value_invalid")
    source_id = document.source.artifact_id
    return [
        _metric(
            "rseqc.tin.mean_score",
            "RSeQC mean transcript integrity number",
            mean,
            "score",
            sample_id,
            source_id,
        ),
        _metric(
            "rseqc.tin.median_score",
            "RSeQC median transcript integrity number",
            median,
            "score",
            sample_id,
            source_id,
        ),
        _metric(
            "rseqc.tin.standard_deviation",
            "RSeQC transcript integrity number standard deviation",
            standard_deviation,
            "score",
            sample_id,
            source_id,
        ),
    ]


def _strict_json(content: bytes) -> object:
    text = _decode_text(content)
    _preflight_json_text(text)

    def reject_constant(_value: str) -> None:
        raise _QcError("source_content_invalid")

    def pairs(values: list[tuple[str, object]]) -> dict[str, object]:
        result: dict[str, object] = {}
        for key, value in values:
            if key in result:
                raise _QcError("source_content_invalid")
            result[key] = value
        return result

    try:
        value = json.loads(
            text,
            parse_float=Decimal,
            parse_int=int,
            parse_constant=reject_constant,
            object_pairs_hook=pairs,
        )
    except (json.JSONDecodeError, RecursionError, TypeError, ValueError) as exc:
        if isinstance(exc, _QcError):
            raise
        raise _QcError("source_content_invalid") from exc
    _validate_json_tree(value)
    return value


def _preflight_json_text(text: str) -> None:
    depth = 0
    structural_tokens = 0
    in_string = False
    escaped = False
    string_length = 0
    for character in text:
        if in_string:
            if escaped:
                escaped = False
                string_length += 1
            elif character == "\\":
                escaped = True
            elif character == '"':
                in_string = False
            else:
                string_length += 1
            if string_length > _MAX_JSON_STRING:
                raise _QcError("source_content_invalid")
            continue
        if character == '"':
            in_string = True
            string_length = 0
        elif character in "[{":
            depth += 1
            structural_tokens += 1
            if depth > _MAX_JSON_DEPTH:
                raise _QcError("source_content_invalid")
        elif character in "]}":
            depth -= 1
            if depth < 0:
                raise _QcError("source_content_invalid")
        elif character in ",:":
            structural_tokens += 1
        if structural_tokens > _MAX_JSON_NODES:
            raise _QcError("source_content_invalid")
    if in_string or escaped or depth != 0:
        raise _QcError("source_content_invalid")


def _validate_json_tree(root: object) -> None:
    nodes = 0
    stack: list[tuple[object, int]] = [(root, 1)]
    while stack:
        value, depth = stack.pop()
        nodes += 1
        if nodes > _MAX_JSON_NODES or depth > _MAX_JSON_DEPTH:
            raise _QcError("source_content_invalid")
        if isinstance(value, str):
            if len(value) > _MAX_JSON_STRING or "\x00" in value:
                raise _QcError("source_content_invalid")
        elif isinstance(value, Mapping):
            for key, child in value.items():
                if not isinstance(key, str) or len(key) > 256 or "\x00" in key:
                    raise _QcError("source_content_invalid")
                stack.append((child, depth + 1))
        elif isinstance(value, list):
            stack.extend((child, depth + 1) for child in value)
        elif isinstance(value, Decimal):
            if not value.is_finite():
                raise _QcError("source_content_invalid")
        elif value is not None and not isinstance(value, (bool, int)):
            raise _QcError("source_content_invalid")


def _tsv_rows(content: bytes) -> list[list[str]]:
    text = _decode_text(content)
    rows = [line.split("\t") for line in text.splitlines() if line]
    if not rows or len(rows) > 5_000:
        raise _QcError("source_content_invalid")
    for row in rows:
        if not 1 <= len(row) <= 128 or any(len(cell) > 4096 for cell in row):
            raise _QcError("source_content_invalid")
    return rows


def _decode_text(content: bytes) -> str:
    if len(content) > _MAX_SOURCE_BYTES or b"\x00" in content:
        raise _QcError("source_content_invalid")
    try:
        text = content.decode("utf-8")
    except UnicodeDecodeError as exc:
        raise _QcError("source_content_invalid") from exc
    lines = text.splitlines()
    if len(lines) > _MAX_TEXT_LINES or any(
        len(line) > _MAX_LINE_LENGTH for line in lines
    ):
        raise _QcError("source_content_invalid")
    return text


def _percent(value: object) -> Decimal:
    return _fraction(_decimal(value) / Decimal(100))


def _decimal(value: object) -> Decimal:
    if isinstance(value, bool):
        raise _QcError("source_content_invalid")
    if isinstance(value, Decimal):
        result = value
    elif isinstance(value, int):
        result = Decimal(value)
    elif isinstance(value, str):
        if not value or len(value) > 128:
            raise _QcError("source_content_invalid")
        try:
            result = Decimal(value)
        except InvalidOperation as exc:
            raise _QcError("source_content_invalid") from exc
    elif isinstance(value, float) and math.isfinite(value):
        result = Decimal(str(value))
    else:
        raise _QcError("source_content_invalid")
    if not result.is_finite():
        raise _QcError("source_content_invalid")
    if result and not (
        _MIN_DECIMAL_ADJUSTED <= result.adjusted() <= _MAX_DECIMAL_ADJUSTED
    ):
        raise _QcError("source_content_invalid")
    return result


def _nonnegative_decimal(value: object) -> Decimal:
    result = _decimal(value)
    if result < 0:
        raise _QcError("metric_value_invalid")
    return result


def _nonnegative_int(value: object) -> int:
    if isinstance(value, bool):
        raise _QcError("source_content_invalid")
    if isinstance(value, int):
        result = value
    elif isinstance(value, str) and re.fullmatch(r"0|[1-9][0-9]*", value):
        result = int(value)
    else:
        raise _QcError("source_content_invalid")
    if result < 0 or result > 10**25:
        raise _QcError("metric_value_invalid")
    return result


def _nonnegative_integral_decimal(value: object) -> int:
    if (
        not isinstance(value, str)
        or re.fullmatch(r"(?:0|[1-9][0-9]*)(?:\.0+)?", value) is None
    ):
        raise _QcError("source_content_invalid")
    result = _decimal(value)
    if result < 0 or result != result.to_integral_value() or result > Decimal(10**25):
        raise _QcError("metric_value_invalid")
    return int(result)


def _fraction(value: Decimal) -> Decimal:
    if value < 0 or value > 1:
        raise _QcError("metric_value_invalid")
    return _canonical_decimal(value)


def _divide(numerator: int, denominator: int) -> Decimal:
    if denominator <= 0 or numerator < 0 or numerator > denominator:
        raise _QcError("metric_value_invalid")
    return _fraction(Decimal(numerator) / Decimal(denominator))


def _canonical_decimal(value: Decimal) -> Decimal:
    if not value.is_finite():
        raise _QcError("metric_value_invalid")
    if value and not (
        _MIN_DECIMAL_ADJUSTED <= value.adjusted() <= _MAX_DECIMAL_ADJUSTED
    ):
        raise _QcError("metric_value_invalid")
    exponent = value.as_tuple().exponent
    if not isinstance(exponent, int):
        raise _QcError("metric_value_invalid")
    if exponent < -12:
        value = value.quantize(_RATIO_QUANTUM, rounding=ROUND_HALF_EVEN)
    text = (
        format(value, "f").rstrip("0").rstrip(".")
        if "." in format(value, "f")
        else format(value, "f")
    )
    if text in {"", "-0"}:
        text = "0"
    if re.fullmatch(r"-?(?:0|[1-9][0-9]{0,25})(?:\.[0-9]{1,12})?", text) is None:
        raise _QcError("metric_value_invalid")
    return Decimal(text)


def _metric(
    key: str,
    name: str,
    value: Decimal,
    unit: str,
    sample_id: str,
    source_artifact_id: str,
    *,
    qc_flag: str | None = None,
) -> ExtractedQcMetricCandidate:
    value = _canonical_decimal(value)
    if unit == "count":
        if value < 0 or value != value.to_integral_value():
            raise _QcError("metric_value_invalid")
    elif unit == "fraction":
        value = _fraction(value)
    elif unit in {"ratio", "score"}:
        if value < 0:
            raise _QcError("metric_value_invalid")
    else:
        raise _QcError("metric_contract_invalid")
    return ExtractedQcMetricCandidate(
        metric_key=key,
        display_name=name,
        value=value,
        unit=unit,
        scope="sample",
        sample_id=sample_id,
        experiment_id=None,
        assay=_ASSAY,
        qc_flag=qc_flag,
        source_artifact_id=source_artifact_id,
    )


def _failure(reason_code: str, *, path: str = "qc_sources") -> Result[Any]:
    return Result.failure(
        [
            Issue(
                code="BULK_RNASEQ_QC_INVALID",
                message="Bulk RNA-seq QC summaries could not be extracted.",
                severity="error",
                path=path,
                source="bulk_rnaseq_adapter",
                context={"reason_code": reason_code},
            )
        ]
    )
