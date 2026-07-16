"""Closed, deterministic artifact discovery for nf-core/rnaseq 3.26.0.

The catalog in this module is deliberately a projection of the pinned
``star_salmon`` output contract.  Every accepted path is derived from validated
sample identities and normalized adapter-owned parameters.  A bounded,
single-level audit of the fixed STAR/Salmon namespaces rejects foreign sample
identities; no recursive output-tree discovery is performed.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from decimal import Decimal, InvalidOperation
import errno
import os
from pathlib import Path, PurePosixPath
import stat
import re
from typing import Any

from encode_pipeline.adapters.bulk_rnaseq.upstream import (
    NFCORE_RNASEQ_COMMIT,
    NFCORE_RNASEQ_RELEASE,
)
from encode_pipeline.adapters.bulk_rnaseq.results_contract import (
    effective_downstream_layout,
    effective_rseqc_modules,
    load_bulk_rnaseq_results_contract,
)
from encode_pipeline.adapters.bulk_rnaseq.validation import (
    validate_bulk_rnaseq_inputs,
)
from encode_pipeline.platform.adapters import (
    ExtractedArtifactCandidate,
    MAX_SAMPLE_ROWS,
    WorkflowInputs,
)
from encode_pipeline.platform.results import Issue, Result


_ASSAY = "bulk-rnaseq"
_ALIGNER = "star_salmon"
_MAX_ARTIFACT_BYTES = 1024**4
_MAX_ARTIFACT_CANDIDATES = MAX_SAMPLE_ROWS * 128 + 128
_MAX_PATH_COMPONENTS = 20
_MAX_STATUS_TABLE_BYTES = MAX_SAMPLE_ROWS * (128 + 1 + 64 + 1) + 128
_MAX_AUDITED_NAMESPACE_ENTRIES = MAX_SAMPLE_ROWS * 128 + 128
_DEFAULT_MIN_TRIMMED_READS = Decimal("10000")
_DEFAULT_MIN_MAPPED_PERCENT = Decimal("5")

# Directory names published directly below the fixed ``star_salmon`` route.
# Any other directory is a per-sample Salmon result directory and therefore
# must be one of the validated sample identities.  This is a bounded,
# single-level namespace check, not an output-tree walk.
_STAR_SALMON_FIXED_DIRECTORIES = frozenset(
    {
        "bigwig",
        "contaminants",
        "deseq2_qc",
        "dupradar",
        "featurecounts",
        "log",
        "picard_metrics",
        "preseq",
        "qualimap",
        "rseqc",
        "rustqc",
        "samtools_stats",
        "stringtie",
        "umicollapse",
        "umitools",
        "unmapped",
    }
)
_STAR_SAMPLE_FILE_SUFFIXES = (
    ".markdup.sorted.bam.bai",
    ".markdup.sorted.bam.csi",
    ".markdup.sorted.bam",
    ".umi_dedup.sorted.bam.bai",
    ".umi_dedup.sorted.bam.csi",
    ".umi_dedup.sorted.bam",
    ".umi_dedup.transcriptome.filtered.bam",
    ".umi_dedup.transcriptome.sorted.bam.bai",
    ".umi_dedup.transcriptome.sorted.bam.csi",
    ".umi_dedup.transcriptome.sorted.bam",
    ".umi_dedup.transcriptome.bam",
    ".sorted.bam.bai",
    ".sorted.bam.csi",
    ".sorted.bam",
    ".transcriptome.sorted.bam.bai",
    ".transcriptome.sorted.bam.csi",
    ".transcriptome.sorted.bam",
    ".Aligned.out.bam",
    ".Aligned.toTranscriptome.out.bam",
)
_STAR_LOG_SAMPLE_SUFFIXES = (".Log.final.out", ".SJ.out.tab")


@dataclass(frozen=True)
class _ArtifactSpec:
    output_type: str
    relative_path: str
    mime_type: str
    metadata: Mapping[str, object]
    required: bool = True


class _UnsafePathError(ValueError):
    pass


class _ArtifactLimitError(ValueError):
    pass


class _ArtifactRaceError(ValueError):
    pass


class _StatusTableError(ValueError):
    pass


class _SampleStatusUnavailable(ValueError):
    pass


class _UnknownSampleError(ValueError):
    pass


def discover_bulk_rnaseq_artifacts(
    inputs: WorkflowInputs,
    workspace: str | Path,
) -> Result[tuple[ExtractedArtifactCandidate, ...]]:
    """Discover the complete bounded artifact projection for one succeeded run."""
    try:
        load_bulk_rnaseq_results_contract()
    except (OSError, TypeError, ValueError):
        return _failure(
            "BULK_RNASEQ_ARTIFACT_CONTRACT_INVALID",
            "artifact_contract_invalid",
        )
    validated = validate_bulk_rnaseq_inputs(inputs)
    if validated.is_failure:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_INPUTS_INVALID",
            "validated_inputs_required",
        )
    try:
        normalized = _normalized_contract(validated.value)
        params = _mapping(normalized.get("nfcore_params"))
        samples = _normalized_samples(normalized.get("samples"))
        _require_supported_profile(params)
        results_descriptor = _open_results_root(workspace)
    except _UnsupportedProfileError:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_PROFILE_UNSUPPORTED",
            "unaudited_output_profile",
        )
    except _ArtifactLimitError:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_LIMIT_EXCEEDED",
            "artifact_catalog_limit",
        )
    except FileNotFoundError:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_REQUIRED_MISSING",
            "required_output_missing",
        )
    except _UnsafePathError:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_PATH_UNSAFE",
            "unsafe_artifact_path",
        )
    except (OSError, TypeError, ValueError):
        return _failure(
            "BULK_RNASEQ_ARTIFACT_CONTRACT_INVALID",
            "artifact_contract_invalid",
        )

    candidates: list[ExtractedArtifactCandidate] = []
    try:
        _reject_unknown_star_salmon_samples(results_descriptor, samples)
        trimmed_failed, mapped_failed, status_specs = _read_sample_status_tables(
            results_descriptor,
            params,
            samples,
        )
        specs, auto_bigwig = _expected_specs(
            params,
            samples,
            trimmed_failed=trimmed_failed,
            mapped_failed=mapped_failed,
        )
        specs.extend(status_specs)
        if len(specs) + len(auto_bigwig) * 3 > _MAX_ARTIFACT_CANDIDATES:
            raise _ArtifactLimitError
        for spec in specs:
            present = _verified_regular_file(results_descriptor, spec.relative_path)
            if not present:
                if spec.required:
                    if params.get("skip_multiqc") is True and _threshold_sensitive(
                        results_descriptor,
                        spec,
                        params=params,
                    ):
                        raise _SampleStatusUnavailable
                    return _failure(
                        "BULK_RNASEQ_ARTIFACT_REQUIRED_MISSING",
                        "required_output_missing",
                    )
                continue
            candidates.append(_candidate(spec))

        for sample, metadata in auto_bigwig:
            candidates.extend(
                _discover_auto_bigwig(results_descriptor, sample, metadata)
            )
    except FileNotFoundError:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_REQUIRED_MISSING",
            "required_output_missing",
        )
    except _ArtifactLimitError:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_LIMIT_EXCEEDED",
            "artifact_file_limit",
        )
    except _ArtifactRaceError:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_RACE_DETECTED",
            "artifact_identity_changed",
        )
    except _StatusTableError:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_STATUS_INVALID",
            "sample_status_invalid",
        )
    except _SampleStatusUnavailable:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_STATUS_UNAVAILABLE",
            "sample_status_unavailable",
        )
    except _UnknownSampleError:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_UNKNOWN_SAMPLE",
            "unknown_sample_output",
        )
    except (OSError, _UnsafePathError, TypeError, ValueError):
        return _failure(
            "BULK_RNASEQ_ARTIFACT_PATH_UNSAFE",
            "unsafe_artifact_path",
        )
    finally:
        os.close(results_descriptor)

    ordered = tuple(
        sorted(candidates, key=lambda item: (item.output_type, item.relative_path))
    )
    if len(ordered) > _MAX_ARTIFACT_CANDIDATES:
        return _failure(
            "BULK_RNASEQ_ARTIFACT_LIMIT_EXCEEDED",
            "artifact_catalog_limit",
        )
    return Result.success(ordered)


def _normalized_contract(value: object) -> Mapping[str, object]:
    normalized = _mapping(value)
    contract = _mapping(normalized.get("contract"))
    if (
        contract.get("workflow_id") != "bulk-rnaseq"
        or contract.get("nfcore_rnaseq_release") != NFCORE_RNASEQ_RELEASE
        or contract.get("nfcore_rnaseq_commit") != NFCORE_RNASEQ_COMMIT
    ):
        raise ValueError("unexpected normalized contract")
    return normalized


def _threshold_sensitive(
    results_descriptor: int,
    spec: _ArtifactSpec,
    *,
    params: Mapping[str, object],
) -> bool:
    sample = spec.metadata.get("sample_id")
    if not isinstance(sample, str):
        return False
    analysis_started = _verified_regular_file(
        results_descriptor,
        f"{_ALIGNER}/log/{sample}.Log.final.out",
    )
    if not analysis_started:
        return spec.output_type.startswith(
            (
                "bulk_rnaseq.star.",
                "bulk_rnaseq.salmon.",
                "bulk_rnaseq.featurecounts.",
                "bulk_rnaseq.rseqc.",
                "bulk_rnaseq.umi.",
            )
        )
    if spec.output_type.startswith(
        (
            "bulk_rnaseq.star.bigwig.",
            "bulk_rnaseq.featurecounts.",
            "bulk_rnaseq.rseqc.",
        )
    ):
        return True
    mark_duplicates = (
        params.get("skip_markduplicates") is False
        and params.get("with_umi") is not True
    )
    return mark_duplicates and spec.output_type.startswith("bulk_rnaseq.star.bam")


def _mapping(value: object) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise TypeError("expected mapping")
    return value


def _normalized_samples(value: object) -> tuple[Mapping[str, object], ...]:
    if not isinstance(value, list):
        raise TypeError("expected normalized samples")
    by_sample: dict[str, Mapping[str, object]] = {}
    for row_value in value:
        row = _mapping(row_value)
        sample = row.get("sample")
        layout = row.get("layout")
        strandedness = row.get("strandedness")
        if (
            not isinstance(sample, str)
            or not sample
            or layout not in {"SE", "PE"}
            or strandedness not in {"auto", "forward", "reverse", "unstranded"}
        ):
            raise ValueError("invalid normalized sample")
        previous = by_sample.setdefault(sample, row)
        if (
            previous.get("layout") != layout
            or previous.get("strandedness") != strandedness
        ):
            raise ValueError("sample semantics changed after validation")
    return tuple(by_sample[name] for name in sorted(by_sample))


def _reject_unknown_star_salmon_samples(
    results_descriptor: int,
    samples: tuple[Mapping[str, object], ...],
) -> None:
    """Reject sample-bearing files outside the validated sample identity set."""
    known_samples = frozenset(str(row["sample"]) for row in samples)
    known_route_files = frozenset(
        f"{sample}{suffix}"
        for sample in known_samples
        for suffix in _STAR_SAMPLE_FILE_SUFFIXES
    )
    known_log_files = frozenset(
        f"{sample}{suffix}"
        for sample in known_samples
        for suffix in _STAR_LOG_SAMPLE_SUFFIXES
    )
    route_descriptor = _open_optional_directory(results_descriptor, _ALIGNER)
    if route_descriptor is None:
        return
    try:
        with os.scandir(route_descriptor) as entries:
            for count, entry in enumerate(entries, start=1):
                if count > _MAX_AUDITED_NAMESPACE_ENTRIES:
                    raise _ArtifactLimitError
                name = entry.name
                if entry.is_dir(follow_symlinks=False):
                    if (
                        name not in known_samples
                        and name not in _STAR_SALMON_FIXED_DIRECTORIES
                    ):
                        raise _UnknownSampleError
                    continue
                if name in known_route_files:
                    continue
                sample = _sample_from_star_route_filename(name)
                if sample is not None and sample not in known_samples:
                    raise _UnknownSampleError
    finally:
        os.close(route_descriptor)

    log_descriptor = _open_optional_directory(
        results_descriptor,
        f"{_ALIGNER}/log",
    )
    if log_descriptor is None:
        return
    try:
        with os.scandir(log_descriptor) as entries:
            for count, entry in enumerate(entries, start=1):
                if count > _MAX_AUDITED_NAMESPACE_ENTRIES:
                    raise _ArtifactLimitError
                if entry.name in known_log_files:
                    continue
                sample = _sample_from_suffix(entry.name, _STAR_LOG_SAMPLE_SUFFIXES)
                if sample is not None and sample not in known_samples:
                    raise _UnknownSampleError
    finally:
        os.close(log_descriptor)


def _sample_from_star_route_filename(filename: str) -> str | None:
    return _sample_from_suffix(filename, _STAR_SAMPLE_FILE_SUFFIXES)


def _sample_from_suffix(filename: str, suffixes: tuple[str, ...]) -> str | None:
    for suffix in suffixes:
        if filename.endswith(suffix):
            sample = filename[: -len(suffix)]
            if not sample:
                raise _UnknownSampleError
            return sample
    return None


def _open_optional_directory(
    results_descriptor: int,
    relative_path: str,
) -> int | None:
    path = _validated_relative_path(relative_path)
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
    current = os.dup(results_descriptor)
    try:
        for component in path.parts:
            try:
                child = os.open(component, flags, dir_fd=current)
            except OSError as exc:
                if exc.errno == errno.ENOENT:
                    os.close(current)
                    return None
                raise _UnsafePathError from exc
            os.close(current)
            current = child
        return current
    except Exception:
        os.close(current)
        raise


def _read_sample_status_tables(
    results_descriptor: int,
    params: Mapping[str, object],
    samples: tuple[Mapping[str, object], ...],
) -> tuple[frozenset[str], frozenset[str], list[_ArtifactSpec]]:
    if params.get("skip_multiqc") is not False:
        return frozenset(), frozenset(), []
    known_samples = frozenset(str(row["sample"]) for row in samples)
    root = f"multiqc/{_ALIGNER}/multiqc_report_data"
    definitions = (
        (
            "trimmed",
            f"{root}/multiqc_fail_trimmed_samples_table.txt",
            "Sample\tReads after trimming",
            "bulk_rnaseq.multiqc.fail_trimmed_samples",
            _decimal_param(
                params,
                "min_trimmed_reads",
                _DEFAULT_MIN_TRIMMED_READS,
            ),
        ),
        (
            "mapped",
            f"{root}/multiqc_fail_mapped_samples_table.txt",
            "Sample\tSTAR uniquely mapped reads (%)",
            "bulk_rnaseq.multiqc.fail_mapped_samples",
            _decimal_param(
                params,
                "min_mapped_reads",
                _DEFAULT_MIN_MAPPED_PERCENT,
            ),
        ),
    )
    failures: dict[str, frozenset[str]] = {}
    specs: list[_ArtifactSpec] = []
    for kind, path, header, output_type, threshold in definitions:
        content = _read_optional_bounded_file(results_descriptor, path)
        if content is None:
            failures[kind] = frozenset()
            continue
        failures[kind] = _parse_status_table(
            content,
            kind=kind,
            expected_header=header,
            known_samples=known_samples,
            threshold=threshold,
        )
        specs.append(
            _ArtifactSpec(
                output_type=output_type,
                relative_path=path,
                mime_type="text/tab-separated-values",
                metadata=_experiment_metadata(),
            )
        )
    if failures["trimmed"] & failures["mapped"]:
        raise _StatusTableError
    return failures["trimmed"], failures["mapped"], specs


def _parse_status_table(
    content: bytes,
    *,
    kind: str,
    expected_header: str,
    known_samples: frozenset[str],
    threshold: Decimal,
) -> frozenset[str]:
    try:
        text = content.decode("utf-8", errors="strict")
    except UnicodeError as exc:
        raise _StatusTableError from exc
    if not text or "\x00" in text or "\r" in text or text.startswith("\ufeff"):
        raise _StatusTableError
    lines = text.split("\n")
    if lines[-1] == "":
        lines.pop()
    if (
        not lines
        or lines[0] != expected_header
        or len(lines) > len(known_samples) + 1
        or any(not line for line in lines)
    ):
        raise _StatusTableError
    failed: set[str] = set()
    observed: set[str] = set()
    for line in lines[1:]:
        fields = line.split("\t")
        if len(fields) != 2:
            raise _StatusTableError
        sample, raw_value = fields
        if (
            sample not in known_samples
            or sample in observed
            or not raw_value
            or len(raw_value) > 64
        ):
            raise _StatusTableError
        observed.add(sample)
        try:
            if re.fullmatch(r"(?:0|[1-9][0-9]*)(?:\.[0-9]+)?", raw_value) is None:
                raise _StatusTableError
            value = Decimal(raw_value)
        except InvalidOperation as exc:
            raise _StatusTableError from exc
        if not value.is_finite() or value < 0:
            raise _StatusTableError
        if kind == "trimmed" and value != value.to_integral_value():
            raise _StatusTableError
        if kind == "mapped" and value > 100:
            raise _StatusTableError
        if kind == "trimmed" and value > threshold:
            raise _StatusTableError
        if kind == "mapped" and value >= threshold:
            raise _StatusTableError
        if kind != "trimmed" or value < threshold:
            failed.add(sample)
    return frozenset(failed)


def _decimal_param(
    params: Mapping[str, object],
    name: str,
    default: Decimal,
) -> Decimal:
    raw = params.get(name)
    if raw is None:
        return default
    if isinstance(raw, bool) or not isinstance(raw, (int, float)):
        raise _StatusTableError
    try:
        value = Decimal(str(raw))
    except InvalidOperation as exc:
        raise _StatusTableError from exc
    if not value.is_finite() or value < 0:
        raise _StatusTableError
    return value


def _require_supported_profile(params: Mapping[str, object]) -> None:
    if params.get("aligner") != _ALIGNER:
        raise ValueError("unexpected product route")
    if params.get("save_reference") is True or params.get("skip_stringtie") is False:
        raise _UnsupportedProfileError


class _UnsupportedProfileError(ValueError):
    pass


def _expected_specs(
    params: Mapping[str, object],
    samples: tuple[Mapping[str, object], ...],
    *,
    trimmed_failed: frozenset[str],
    mapped_failed: frozenset[str],
) -> tuple[list[_ArtifactSpec], list[tuple[str, Mapping[str, object]]]]:
    # These two namespaces contain directory trees whose complete closed file
    # contracts have intentionally not been admitted in PR152.
    if params.get("save_reference") is True or params.get("skip_stringtie") is False:
        raise _UnsupportedProfileError

    specs: list[_ArtifactSpec] = []
    auto_bigwig: list[tuple[str, Mapping[str, object]]] = []

    def add(
        output_type: str,
        relative_path: str,
        mime_type: str,
        *,
        sample: str | None = None,
        required: bool = True,
    ) -> None:
        metadata = _sample_metadata(sample) if sample else _experiment_metadata()
        specs.append(
            _ArtifactSpec(
                output_type=output_type,
                relative_path=relative_path,
                mime_type=mime_type,
                metadata=metadata,
                required=required,
            )
        )

    add(
        "bulk_rnaseq.pipeline.software_versions",
        "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml",
        "text/yaml",
    )

    csi = params.get("bam_csi_index") is True
    with_umi = params.get("with_umi") is True
    mark_duplicates = params.get("skip_markduplicates") is False and not with_umi
    final_suffix = (
        "umi_dedup.sorted"
        if with_umi
        else "markdup.sorted"
        if mark_duplicates
        else "sorted"
    )
    publish_final_bam = (
        not with_umi
        or params.get("save_align_intermeds") is True
        or params.get("save_umi_intermeds") is True
    )

    for row in samples:
        sample = str(row["sample"])
        layout = str(row["layout"])
        downstream_layout = effective_downstream_layout(layout, params)
        strandedness = str(row["strandedness"])
        roles = _read_roles(layout)
        downstream_roles = _read_roles(downstream_layout)

        index_suffix = "csi" if csi else "bai"
        if publish_final_bam:
            add(
                "bulk_rnaseq.star.bam",
                f"{_ALIGNER}/{sample}.{final_suffix}.bam",
                "application/vnd.samtools.bam",
                sample=sample,
            )
            add(
                f"bulk_rnaseq.star.bam_index.{index_suffix}",
                f"{_ALIGNER}/{sample}.{final_suffix}.bam.{index_suffix}",
                "application/octet-stream",
                sample=sample,
            )
        if with_umi or (mark_duplicates and params.get("save_align_intermeds") is True):
            add(
                "bulk_rnaseq.star.sorted_intermediate_bam",
                f"{_ALIGNER}/{sample}.sorted.bam",
                "application/vnd.samtools.bam",
                sample=sample,
            )
            add(
                f"bulk_rnaseq.star.sorted_intermediate_bam_index.{index_suffix}",
                f"{_ALIGNER}/{sample}.sorted.bam.{index_suffix}",
                "application/octet-stream",
                sample=sample,
            )
        if with_umi:
            add(
                "bulk_rnaseq.star.transcriptome_sorted_intermediate_bam",
                f"{_ALIGNER}/{sample}.transcriptome.sorted.bam",
                "application/vnd.samtools.bam",
                sample=sample,
            )
            add(
                "bulk_rnaseq.star.transcriptome_sorted_intermediate_bam_index."
                f"{index_suffix}",
                f"{_ALIGNER}/{sample}.transcriptome.sorted.bam.{index_suffix}",
                "application/octet-stream",
                sample=sample,
            )
        add(
            "bulk_rnaseq.star.log_final",
            f"{_ALIGNER}/log/{sample}.Log.final.out",
            "text/plain",
            sample=sample,
        )
        add(
            "bulk_rnaseq.star.splice_junctions",
            f"{_ALIGNER}/log/{sample}.SJ.out.tab",
            "text/tab-separated-values",
            sample=sample,
        )
        add(
            "bulk_rnaseq.salmon.quant_transcript",
            f"{_ALIGNER}/{sample}/quant.sf",
            "text/tab-separated-values",
            sample=sample,
        )
        add(
            "bulk_rnaseq.salmon.quant_gene",
            f"{_ALIGNER}/{sample}/quant.genes.sf",
            "text/tab-separated-values",
            sample=sample,
        )
        add(
            "bulk_rnaseq.salmon.meta_info",
            f"{_ALIGNER}/{sample}/aux_info/meta_info.json",
            "application/json",
            sample=sample,
        )
        for stat_name in ("flagstat", "idxstats"):
            add(
                f"bulk_rnaseq.samtools.{stat_name}",
                f"{_ALIGNER}/samtools_stats/{sample}.{final_suffix}.bam.{stat_name}",
                "text/plain",
                sample=sample,
                required=False,
            )

        if params.get("skip_fastqc") is False:
            _add_fastqc_specs(add, sample, roles, stage="raw")
            if params.get("skip_trimming") is False:
                _add_fastqc_specs(
                    add,
                    sample,
                    downstream_roles,
                    stage="trimmed",
                    trimmer=str(params.get("trimmer")),
                )
            if params.get("remove_ribo_rna") is True:
                _add_fastqc_specs(add, sample, downstream_roles, stage="filtered")

        if params.get("skip_trimming") is False:
            _add_trimming_specs(add, params, sample, downstream_layout)

        if params.get("save_merged_fastq") is True:
            for role, suffix in _role_suffixes(layout, single="", paired=("_1", "_2")):
                add(
                    f"bulk_rnaseq.fastq.merged.{role}",
                    f"fastq/{sample}{suffix}.merged.fastq.gz",
                    "application/gzip",
                    sample=sample,
                )

        if params.get("save_unaligned") is True:
            # STAR consistently names the sole SE stream mate1 as ``_1``.
            unaligned_suffixes = (
                (("single", "_1"),)
                if downstream_layout == "SE"
                else (("read1", "_1"), ("read2", "_2"))
            )
            for role, suffix in unaligned_suffixes:
                add(
                    f"bulk_rnaseq.star.unaligned.{role}",
                    f"{_ALIGNER}/unmapped/{sample}.unmapped{suffix}.fastq.gz",
                    "application/gzip",
                    sample=sample,
                )

        if params.get("save_align_intermeds") is True:
            for suffix, output_type in (
                ("Aligned.out.bam", "bulk_rnaseq.star.aligned_genome_bam"),
                (
                    "Aligned.toTranscriptome.out.bam",
                    "bulk_rnaseq.star.aligned_transcriptome_bam",
                ),
            ):
                add(
                    output_type,
                    f"{_ALIGNER}/{sample}.{suffix}",
                    "application/vnd.samtools.bam",
                    sample=sample,
                )

        if with_umi:
            if params.get("save_umi_intermeds") is True:
                if params.get("skip_umi_extract") is False:
                    for role, suffix in _role_suffixes(
                        layout, single="", paired=("_1", "_2")
                    ):
                        add(
                            f"bulk_rnaseq.umi.extracted.{role}",
                            f"umitools/{sample}.umi_extract{suffix}.fastq.gz",
                            "application/gzip",
                            sample=sample,
                        )
            if (
                params.get("save_align_intermeds") is True
                or params.get("save_umi_intermeds") is True
            ):
                for suffix, output_type, required in (
                    (
                        "umi_dedup.transcriptome.bam",
                        "transcriptome_name_sorted_bam",
                        True,
                    ),
                    (
                        "umi_dedup.transcriptome.filtered.bam",
                        "transcriptome_filtered_bam",
                        downstream_layout == "PE",
                    ),
                    (
                        "umi_dedup.transcriptome.sorted.bam",
                        "transcriptome_sorted_bam",
                        True,
                    ),
                    (
                        "umi_dedup.transcriptome.sorted.bam.bai",
                        "transcriptome_sorted_bam_index_bai",
                        True,
                    ),
                ):
                    add(
                        f"bulk_rnaseq.umi.{output_type}",
                        f"{_ALIGNER}/{sample}.{suffix}",
                        (
                            "application/octet-stream"
                            if suffix.endswith((".bai", ".csi"))
                            else "application/vnd.samtools.bam"
                        ),
                        sample=sample,
                        required=required,
                    )
            if params.get("umitools_dedup_stats") is True:
                for target in ("genomic", "transcriptomic"):
                    stem = (
                        f"{sample}.umi_dedup.sorted"
                        if target == "genomic"
                        else f"{sample}.umi_dedup.transcriptome.sorted"
                    )
                    for statistic in (
                        "edit_distance",
                        "per_umi",
                        "per_umi_per_position",
                    ):
                        add(
                            f"bulk_rnaseq.umi.{target}_{statistic}",
                            f"{_ALIGNER}/umitools/{stem}_{statistic}.tsv",
                            "text/tab-separated-values",
                            sample=sample,
                        )

        if params.get("remove_ribo_rna") is True:
            _add_rrna_specs(add, params, sample, downstream_layout)

        if params.get("skip_bigwig") is False:
            if strandedness == "unstranded":
                add(
                    "bulk_rnaseq.star.bigwig.combined",
                    f"{_ALIGNER}/bigwig/{sample}.bigWig",
                    "application/octet-stream",
                    sample=sample,
                )
            elif strandedness in {"forward", "reverse"}:
                add(
                    "bulk_rnaseq.star.bigwig.combined",
                    f"{_ALIGNER}/bigwig/{sample}.bigWig",
                    "application/octet-stream",
                    sample=sample,
                )
                _add_directional_bigwigs(add, sample)
            else:
                auto_bigwig.append((sample, _sample_metadata(sample)))

        if params.get("skip_biotype_qc") is False:
            for suffix, output_type in (
                ("featureCounts.tsv.summary", "bulk_rnaseq.featurecounts.summary"),
                ("biotype_counts_mqc.tsv", "bulk_rnaseq.featurecounts.biotype"),
                (
                    "biotype_counts_rrna_mqc.tsv",
                    "bulk_rnaseq.featurecounts.rrna_biotype",
                ),
            ):
                add(
                    output_type,
                    f"{_ALIGNER}/featurecounts/{sample}.{suffix}",
                    "text/tab-separated-values",
                    sample=sample,
                )

        if params.get("skip_rseqc") is False:
            _add_rseqc_specs(add, params, sample, downstream_layout)

    for suffix in (
        "gene_counts",
        "gene_counts_length_scaled",
        "gene_counts_scaled",
        "gene_lengths",
        "gene_tpm",
        "transcript_counts",
        "transcript_lengths",
        "transcript_tpm",
        "tx2gene",
        "tx2gene_augmented",
    ):
        add(
            f"bulk_rnaseq.salmon.merged.{suffix}",
            f"{_ALIGNER}/salmon.merged.{suffix}.tsv",
            "text/tab-separated-values",
        )
    analysis_present = len(trimmed_failed) < len(samples)
    post_mapping_present = len(trimmed_failed | mapped_failed) < len(samples)
    if params.get("skip_multiqc") is False:
        _add_multiqc_specs(
            add,
            params,
            analysis_present=analysis_present,
            post_mapping_present=post_mapping_present,
        )

    specs = _filter_specs_for_sample_status(
        specs,
        params=params,
        trimmed_failed=trimmed_failed,
        mapped_failed=mapped_failed,
        analysis_present=analysis_present,
        post_mapping_present=post_mapping_present,
    )
    auto_bigwig = [
        item
        for item in auto_bigwig
        if item[0] not in trimmed_failed and item[0] not in mapped_failed
    ]

    _validate_specs(specs)
    return specs, auto_bigwig


def _filter_specs_for_sample_status(
    specs: list[_ArtifactSpec],
    *,
    params: Mapping[str, object],
    trimmed_failed: frozenset[str],
    mapped_failed: frozenset[str],
    analysis_present: bool,
    post_mapping_present: bool,
) -> list[_ArtifactSpec]:
    pre_alignment_prefixes = (
        "bulk_rnaseq.fastqc.raw.",
        "bulk_rnaseq.trim.",
        "bulk_rnaseq.fastq.merged.",
        "bulk_rnaseq.umi.extract",
        "bulk_rnaseq.umi.extracted.",
    ) + (
        ("bulk_rnaseq.fastqc.trimmed.",)
        if params.get("trimmer") == "trimgalore"
        else ()
    )
    post_mapping_prefixes = (
        "bulk_rnaseq.star.bigwig.",
        "bulk_rnaseq.featurecounts.",
        "bulk_rnaseq.rseqc.",
    )
    mark_duplicates = (
        params.get("skip_markduplicates") is False
        and params.get("with_umi") is not True
    )
    filtered: list[_ArtifactSpec] = []
    for spec in specs:
        sample = spec.metadata.get("sample_id")
        if isinstance(sample, str):
            if sample in trimmed_failed and not spec.output_type.startswith(
                pre_alignment_prefixes
            ):
                continue
            if sample in mapped_failed and (
                spec.output_type.startswith(post_mapping_prefixes)
                or (
                    mark_duplicates
                    and spec.output_type.startswith(
                        ("bulk_rnaseq.star.bam", "bulk_rnaseq.samtools.")
                    )
                )
            ):
                continue
        elif not isinstance(sample, str):
            if not analysis_present and spec.output_type.startswith(
                (
                    "bulk_rnaseq.salmon.merged.",
                    "bulk_rnaseq.multiqc.star",
                    "bulk_rnaseq.multiqc.salmon",
                    "bulk_rnaseq.multiqc.fastqc.filtered",
                )
            ):
                continue
            if not post_mapping_present and spec.output_type.startswith(
                (
                    "bulk_rnaseq.multiqc.picard_dups",
                    "bulk_rnaseq.multiqc.featurecounts",
                    "bulk_rnaseq.multiqc.rseqc.",
                )
            ):
                continue
        filtered.append(spec)
    return filtered


def _add_fastqc_specs(
    add, sample: str, roles: tuple[str, ...], *, stage: str, trimmer: str | None = None
) -> None:
    for index, role in enumerate(roles, start=1):
        if stage == "raw":
            stem = (
                f"{sample}_raw" + (f"_{index}" if len(roles) == 2 else "") + "_fastqc"
            )
            directory = "fastqc/raw"
        elif stage == "filtered":
            stem = (
                f"{sample}_filtered"
                + (f"_{index}" if len(roles) == 2 else "")
                + "_fastqc"
            )
            directory = "fastqc/filtered"
        else:
            if trimmer == "fastp":
                stem = f"{sample}_trimmed" + (f"_{index}" if len(roles) == 2 else "")
                stem += "_fastqc"
            elif len(roles) == 1:
                stem = f"{sample}_trimmed_trimmed_fastqc"
            else:
                stem = f"{sample}_trimmed_{index}_val_{index}_fastqc"
            directory = "fastqc/trim"
            if trimmer not in {"trimgalore", "fastp"}:
                raise ValueError("unsupported trimmer")
        for extension, mime in (("html", "text/html"), ("zip", "application/zip")):
            add(
                f"bulk_rnaseq.fastqc.{stage}.{role}.{extension}",
                f"{directory}/{stem}.{extension}",
                mime,
                sample=sample,
            )


def _add_trimming_specs(
    add, params: Mapping[str, object], sample: str, layout: str
) -> None:
    trimmer = params.get("trimmer")
    if trimmer == "trimgalore":
        # Trim Galore 2.1.0 raw reports include command-line and input filename
        # fields.  Only MultiQC's fixed cutadapt table is public; this branch
        # admits the saved filtered reads when explicitly requested.
        if params.get("save_trimmed") is True:
            trimmed_suffixes = (
                (("single", "_trimmed.fq.gz"),)
                if layout == "SE"
                else (("read1", "_1_val_1.fq.gz"), ("read2", "_2_val_2.fq.gz"))
            )
            for role, suffix in trimmed_suffixes:
                add(
                    f"bulk_rnaseq.trim.trimgalore.filtered.{role}",
                    f"trimgalore/{sample}_trimmed{suffix}",
                    "application/gzip",
                    sample=sample,
                )
    elif trimmer == "fastp":
        # The fixed JSON schema is sufficient for machine QC. The HTML report
        # is not admitted without an upstream content-safety guarantee.
        add(
            "bulk_rnaseq.trim.fastp.json",
            f"fastp/{sample}.fastp.json",
            "application/json",
            sample=sample,
        )
        if params.get("save_trimmed") is True:
            for role, suffix in _role_suffixes(
                layout, single="", paired=("_R1", "_R2")
            ):
                add(
                    f"bulk_rnaseq.trim.fastp.filtered.{role}",
                    f"fastp/{sample}{suffix}.fastp.fastq.gz",
                    "application/gzip",
                    sample=sample,
                )
    else:
        raise ValueError("unsupported trimmer")


def _add_rrna_specs(
    add, params: Mapping[str, object], sample: str, layout: str
) -> None:
    tool = params.get("ribo_removal_tool")
    if tool == "sortmerna":
        if params.get("save_non_ribo_reads") is True:
            suffixes = _role_suffixes(layout, single="", paired=("_1", "_2"))
            for role, suffix in suffixes:
                add(
                    f"bulk_rnaseq.rrna.sortmerna.filtered.{role}",
                    f"sortmerna/{sample}{suffix}.non_rRNA.fastq.gz",
                    "application/gzip",
                    sample=sample,
                )
    elif tool == "bowtie2":
        # PR150 validation allows publication only for SE.  PE still uses the
        # upstream both-mates-unmapped filtering route internally.
        if params.get("save_non_ribo_reads") is True:
            if layout != "SE":
                raise ValueError("Bowtie2 PE filtered FASTQ cannot be published")
            add(
                "bulk_rnaseq.rrna.bowtie2.filtered.single",
                f"bowtie2_rrna/{sample}.bowtie2_rrna.unmapped.fastq.gz",
                "application/gzip",
                sample=sample,
            )
    else:
        raise ValueError("unsupported rRNA removal route")


def _add_directional_bigwigs(add, sample: str) -> None:
    for direction in ("forward", "reverse"):
        add(
            f"bulk_rnaseq.star.bigwig.{direction}",
            f"{_ALIGNER}/bigwig/{sample}.{direction}.bigWig",
            "application/octet-stream",
            sample=sample,
        )


def _add_rseqc_specs(
    add, params: Mapping[str, object], sample: str, layout: str
) -> None:
    modules = _rseqc_modules(params)
    for module in modules:
        if module == "bam_stat":
            add(
                "bulk_rnaseq.rseqc.bam_stat",
                f"{_ALIGNER}/rseqc/bam_stat/{sample}.bam_stat.txt",
                "text/plain",
                sample=sample,
            )
        elif module == "infer_experiment":
            add(
                "bulk_rnaseq.rseqc.infer_experiment",
                f"{_ALIGNER}/rseqc/infer_experiment/{sample}.infer_experiment.txt",
                "text/plain",
                sample=sample,
            )
        elif module == "read_distribution":
            add(
                "bulk_rnaseq.rseqc.read_distribution",
                f"{_ALIGNER}/rseqc/read_distribution/{sample}.read_distribution.txt",
                "text/plain",
                sample=sample,
            )
        elif module == "tin":
            add(
                "bulk_rnaseq.rseqc.tin_summary",
                f"{_ALIGNER}/rseqc/tin/{sample}.summary.txt",
                "text/plain",
                sample=sample,
            )
            add(
                "bulk_rnaseq.rseqc.tin",
                f"{_ALIGNER}/rseqc/tin/{sample}.tin.xls",
                "text/tab-separated-values",
                sample=sample,
            )
        elif module == "inner_distance":
            add(
                "bulk_rnaseq.rseqc.inner_distance.distance",
                f"{_ALIGNER}/rseqc/inner_distance/txt/{sample}.inner_distance.txt",
                "text/plain",
                sample=sample,
            )
            if layout == "PE":
                for suffix in ("freq", "mean"):
                    add(
                        f"bulk_rnaseq.rseqc.inner_distance.{suffix}",
                        f"{_ALIGNER}/rseqc/inner_distance/txt/{sample}.inner_distance_{suffix}.txt",
                        "text/plain",
                        sample=sample,
                    )
        elif module == "junction_annotation":
            for directory, suffix, mime, required in (
                ("bed", "junction.bed", "text/plain", False),
                ("bed", "junction.Interact.bed", "text/plain", False),
                ("xls", "junction.xls", "text/tab-separated-values", True),
            ):
                add(
                    f"bulk_rnaseq.rseqc.junction_annotation.{suffix.lower()}",
                    f"{_ALIGNER}/rseqc/junction_annotation/{directory}/{sample}.{suffix}",
                    mime,
                    sample=sample,
                    required=required,
                )
        elif module == "junction_saturation":
            add(
                "bulk_rnaseq.rseqc.junction_saturation.plot",
                f"{_ALIGNER}/rseqc/junction_saturation/pdf/{sample}.junctionSaturation_plot.pdf",
                "application/pdf",
                sample=sample,
            )
        elif module == "read_duplication":
            for kind in ("pos", "seq"):
                add(
                    f"bulk_rnaseq.rseqc.read_duplication.{kind}",
                    f"{_ALIGNER}/rseqc/read_duplication/xls/{sample}.{kind}.DupRate.xls",
                    "text/tab-separated-values",
                    sample=sample,
                )
        else:
            raise ValueError("unexpected RSeQC module")


def _rseqc_modules(params: Mapping[str, object]) -> tuple[str, ...]:
    return effective_rseqc_modules(params)


def _add_multiqc_specs(
    add,
    params: Mapping[str, object],
    *,
    analysis_present: bool,
    post_mapping_present: bool,
) -> None:
    root = f"multiqc/{_ALIGNER}"
    data = f"{root}/multiqc_report_data"
    rseqc_modules = (
        frozenset(_rseqc_modules(params))
        if params.get("skip_rseqc") is False and post_mapping_present
        else frozenset()
    )
    for filename, output_type, enabled, required in (
        (
            "multiqc_general_stats.txt",
            "bulk_rnaseq.multiqc.general_stats",
            True,
            True,
        ),
        ("multiqc_star.txt", "bulk_rnaseq.multiqc.star", analysis_present, True),
        (
            "multiqc_salmon.txt",
            "bulk_rnaseq.multiqc.salmon",
            analysis_present,
            True,
        ),
        (
            "multiqc_cutadapt.txt",
            "bulk_rnaseq.multiqc.cutadapt",
            params.get("skip_trimming") is False
            and params.get("trimmer") == "trimgalore",
            True,
        ),
        (
            "multiqc_fastp.txt",
            "bulk_rnaseq.multiqc.fastp",
            params.get("skip_trimming") is False and params.get("trimmer") == "fastp",
            True,
        ),
        (
            "multiqc_fastqc_fastqc_raw.txt",
            "bulk_rnaseq.multiqc.fastqc.raw",
            params.get("skip_fastqc") is False,
            True,
        ),
        (
            "multiqc_fastqc_fastqc_trimmed.txt",
            "bulk_rnaseq.multiqc.fastqc.trimmed",
            params.get("skip_fastqc") is False
            and params.get("skip_trimming") is False
            and (params.get("trimmer") == "trimgalore" or analysis_present),
            True,
        ),
        (
            "multiqc_fastqc_fastqc_filtered.txt",
            "bulk_rnaseq.multiqc.fastqc.filtered",
            params.get("skip_fastqc") is False
            and params.get("remove_ribo_rna") is True
            and analysis_present,
            True,
        ),
        (
            "multiqc_picard_dups.txt",
            "bulk_rnaseq.multiqc.picard_dups",
            params.get("skip_markduplicates") is False and post_mapping_present,
            True,
        ),
        (
            "multiqc_featurecounts_biotype_plot.txt",
            "bulk_rnaseq.multiqc.featurecounts_biotype",
            params.get("skip_biotype_qc") is False and post_mapping_present,
            True,
        ),
        (
            "multiqc_rseqc_bam_stat.txt",
            "bulk_rnaseq.multiqc.rseqc.bam_stat",
            params.get("skip_rseqc") is False and "bam_stat" in rseqc_modules,
            False,
        ),
        (
            "multiqc_rseqc_infer_experiment.txt",
            "bulk_rnaseq.multiqc.rseqc.infer_experiment",
            params.get("skip_rseqc") is False and "infer_experiment" in rseqc_modules,
            False,
        ),
        (
            "multiqc_rseqc_read_distribution.txt",
            "bulk_rnaseq.multiqc.rseqc.read_distribution",
            params.get("skip_rseqc") is False and "read_distribution" in rseqc_modules,
            False,
        ),
        (
            "multiqc_rseqc_tin.txt",
            "bulk_rnaseq.multiqc.rseqc.tin",
            params.get("skip_rseqc") is False and "tin" in rseqc_modules,
            False,
        ),
    ):
        # Non-required here means a safe, exact optional machine document, not
        # an arbitrary MultiQC export. JSON, logs, sources and parquet are never
        # admitted by this catalog.
        if enabled:
            add(
                output_type,
                f"{data}/{filename}",
                "text/tab-separated-values",
                required=required,
            )


def _discover_auto_bigwig(
    results_descriptor: int,
    sample: str,
    metadata: Mapping[str, object],
) -> list[ExtractedArtifactCandidate]:
    paths = {
        "combined": f"{_ALIGNER}/bigwig/{sample}.bigWig",
        "forward": f"{_ALIGNER}/bigwig/{sample}.forward.bigWig",
        "reverse": f"{_ALIGNER}/bigwig/{sample}.reverse.bigWig",
    }
    present = {
        name: _verified_regular_file(results_descriptor, path)
        for name, path in paths.items()
    }
    if not present["combined"] or present["forward"] != present["reverse"]:
        raise FileNotFoundError
    return [
        _candidate(
            _ArtifactSpec(
                output_type=f"bulk_rnaseq.star.bigwig.{name}",
                relative_path=paths[name],
                mime_type="application/octet-stream",
                metadata=metadata,
            )
        )
        for name in ("combined", "forward", "reverse")
        if present[name]
    ]


def _validate_specs(specs: list[_ArtifactSpec]) -> None:
    seen: set[tuple[str, str]] = set()
    for spec in specs:
        path = PurePosixPath(spec.relative_path)
        if (
            not spec.output_type.startswith("bulk_rnaseq.")
            or path.as_posix() != spec.relative_path
            or path.is_absolute()
            or len(path.parts) > _MAX_PATH_COMPONENTS
            or any(part in {"", ".", ".."} for part in path.parts)
            or (spec.output_type, spec.relative_path) in seen
        ):
            raise ValueError("invalid artifact specification")
        seen.add((spec.output_type, spec.relative_path))


def _candidate(spec: _ArtifactSpec) -> ExtractedArtifactCandidate:
    return ExtractedArtifactCandidate(
        output_type=spec.output_type,
        relative_path=f"results/{spec.relative_path}",
        mime_type=spec.mime_type,
        metadata=spec.metadata,
    )


def _read_roles(layout: str) -> tuple[str, ...]:
    return ("single",) if layout == "SE" else ("read1", "read2")


def _role_suffixes(
    layout: str,
    *,
    single: str,
    paired: tuple[str, str],
) -> tuple[tuple[str, str], ...]:
    if layout == "SE":
        return (("single", single),)
    return (("read1", paired[0]), ("read2", paired[1]))


def _experiment_metadata() -> Mapping[str, object]:
    # The platform's persisted run is the experiment-level aggregate identity.
    return {"scope": "run", "assay": _ASSAY}


def _sample_metadata(sample: str) -> Mapping[str, object]:
    return {"scope": "sample", "sample_id": sample, "assay": _ASSAY}


def _open_results_root(workspace: str | Path) -> int:
    required_flags = ("O_DIRECTORY", "O_NOFOLLOW", "O_CLOEXEC", "O_PATH")
    if any(not hasattr(os, name) for name in required_flags):
        raise OSError("descriptor-safe file operations are unavailable")
    workspace_path = Path(workspace)
    if (
        not workspace_path.is_absolute()
        or not workspace_path.parts
        or any(part in {"", ".", ".."} for part in workspace_path.parts[1:])
    ):
        raise ValueError("workspace must be a normalized absolute path")
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
    current = os.open("/", flags)
    try:
        for component in (*workspace_path.parts[1:], "results"):
            try:
                child = os.open(component, flags, dir_fd=current)
            except OSError as exc:
                if exc.errno == errno.ENOENT:
                    raise FileNotFoundError from exc
                raise _UnsafePathError from exc
            os.close(current)
            current = child
            info = os.fstat(current)
            if not stat.S_ISDIR(info.st_mode):
                raise _UnsafePathError
        return current
    except Exception:
        os.close(current)
        raise


def _verified_regular_file(results_descriptor: int, relative_path: str) -> bool:
    try:
        first = _open_relative_identity(results_descriptor, relative_path)
    except FileNotFoundError:
        return False
    if not stat.S_ISREG(first[-1].st_mode):
        raise _UnsafePathError
    if first[-1].st_size < 0 or first[-1].st_size > _MAX_ARTIFACT_BYTES:
        raise _ArtifactLimitError
    try:
        second = _open_relative_identity(results_descriptor, relative_path)
    except FileNotFoundError as exc:
        raise _ArtifactRaceError from exc
    if tuple(_identity(info) for info in first) != tuple(
        _identity(info) for info in second
    ):
        raise _ArtifactRaceError
    return True


def _read_optional_bounded_file(
    results_descriptor: int,
    relative_path: str,
) -> bytes | None:
    try:
        descriptor, first = _open_relative_for_read(
            results_descriptor,
            relative_path,
        )
    except FileNotFoundError:
        return None
    try:
        file_info = first[-1]
        if not stat.S_ISREG(file_info.st_mode):
            raise _UnsafePathError
        if file_info.st_size < 0 or file_info.st_size > _MAX_STATUS_TABLE_BYTES:
            raise _ArtifactLimitError
        chunks: list[bytes] = []
        remaining = file_info.st_size + 1
        while remaining:
            chunk = os.read(descriptor, min(remaining, 64 * 1024))
            if not chunk:
                break
            chunks.append(chunk)
            remaining -= len(chunk)
        content = b"".join(chunks)
        if len(content) != file_info.st_size:
            raise _ArtifactRaceError
    finally:
        os.close(descriptor)
    try:
        second = _open_relative_identity(results_descriptor, relative_path)
    except FileNotFoundError as exc:
        raise _ArtifactRaceError from exc
    if tuple(_identity(info) for info in first) != tuple(
        _identity(info) for info in second
    ):
        raise _ArtifactRaceError
    return content


def _open_relative_for_read(
    results_descriptor: int,
    relative_path: str,
) -> tuple[int, tuple[os.stat_result, ...]]:
    path = _validated_relative_path(relative_path)
    directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
    final_flags = os.O_RDONLY | os.O_NONBLOCK | os.O_NOFOLLOW | os.O_CLOEXEC
    current = os.dup(results_descriptor)
    identities: list[os.stat_result] = []
    try:
        for index, component in enumerate(path.parts):
            final = index == len(path.parts) - 1
            try:
                child = os.open(
                    component,
                    final_flags if final else directory_flags,
                    dir_fd=current,
                )
            except OSError as exc:
                if exc.errno == errno.ENOENT:
                    raise FileNotFoundError from exc
                raise _UnsafePathError from exc
            os.close(current)
            current = child
            info = os.fstat(current)
            if not final and not stat.S_ISDIR(info.st_mode):
                raise _UnsafePathError
            identities.append(info)
        return current, tuple(identities)
    except Exception:
        os.close(current)
        raise


def _open_relative_identity(
    results_descriptor: int,
    relative_path: str,
) -> tuple[os.stat_result, ...]:
    path = _validated_relative_path(relative_path)
    directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
    final_flags = os.O_PATH | os.O_NOFOLLOW | os.O_CLOEXEC
    current = os.dup(results_descriptor)
    identities: list[os.stat_result] = []
    try:
        for index, component in enumerate(path.parts):
            final = index == len(path.parts) - 1
            try:
                child = os.open(
                    component,
                    final_flags if final else directory_flags,
                    dir_fd=current,
                )
            except OSError as exc:
                if exc.errno == errno.ENOENT:
                    raise FileNotFoundError from exc
                raise _UnsafePathError from exc
            os.close(current)
            current = child
            info = os.fstat(current)
            if not final and not stat.S_ISDIR(info.st_mode):
                raise _UnsafePathError
            identities.append(info)
        return tuple(identities)
    finally:
        os.close(current)


def _validated_relative_path(relative_path: str) -> PurePosixPath:
    path = PurePosixPath(relative_path)
    if (
        path.as_posix() != relative_path
        or path.is_absolute()
        or not path.parts
        or len(path.parts) > _MAX_PATH_COMPONENTS
        or any(part in {"", ".", ".."} for part in path.parts)
    ):
        raise _UnsafePathError
    return path


def _identity(info: os.stat_result) -> tuple[int, ...]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_nlink,
        info.st_uid,
        info.st_gid,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def _failure(code: str, reason: str) -> Result[tuple[ExtractedArtifactCandidate, ...]]:
    return Result.failure(
        [
            Issue(
                code=code,
                message="Bulk RNA-seq artifact discovery failed.",
                severity="error",
                path="results",
                source="adapter",
                context={"reason": reason},
            )
        ]
    )
