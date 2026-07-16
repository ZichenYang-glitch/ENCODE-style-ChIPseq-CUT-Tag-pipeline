"""Machine-readable QC contract tests for pinned nf-core/rnaseq 3.26.0."""

from __future__ import annotations

from decimal import localcontext
from io import BytesIO
import json
from pathlib import PurePosixPath
import stat
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo

import pytest

import encode_pipeline.adapters.bulk_rnaseq.qc as qc_module
from encode_pipeline.adapters.bulk_rnaseq.qc import (
    BULK_RNASEQ_QC_SOURCE_TYPES,
    extract_bulk_rnaseq_qc_metrics,
)
from encode_pipeline.platform.adapters import (
    MAX_SAMPLE_ROWS,
    QcSourceArtifact,
    QcSourceDocument,
    WorkflowInputs,
)


def _sample(
    sample: str,
    *,
    layout: str = "SE",
    strandedness: str = "forward",
) -> dict[str, str]:
    row = {
        "sample": sample,
        "library": "lib1",
        "lane": "L001",
        "layout": layout,
        "fastq_1": f"/data/{sample}.R1.fastq.gz",
        "strandedness": strandedness,
        "platform": "ILLUMINA",
    }
    if layout == "PE":
        row["fastq_2"] = f"/data/{sample}.R2.fastq.gz"
    return row


def _inputs(
    *,
    samples: list[dict[str, str]] | None = None,
    trimming: dict[str, object] | None = None,
    qc: dict[str, object] | None = None,
    umi: dict[str, object] | None = None,
    advanced: dict[str, object] | None = None,
) -> WorkflowInputs:
    qc_values = {
        "enabled": True,
        "fastqc": False,
        "multiqc": False,
        "rseqc": False,
        "qualimap": False,
        "dupradar": False,
        "biotype": False,
        "deseq2_pca": False,
        "preseq": False,
        "mark_duplicates": False,
        **(qc or {}),
    }
    standard = {
        "reference": {
            "reference_id": "tiny-ref",
            "fasta": "/refs/tiny.fa",
            "fasta_sha256": "a" * 64,
            "gtf": "/refs/tiny.gtf",
            "gtf_sha256": "b" * 64,
        },
        "trimming": {"enabled": False, "tool": "trimgalore"}
        if trimming is None
        else trimming,
        "qc": qc_values,
    }
    if umi is not None:
        standard["umi"] = umi
    config: dict[str, object] = {"standard": standard}
    if advanced:
        config["advanced"] = advanced
    return WorkflowInputs(
        config=config,
        samples=[_sample("S1")] if samples is None else samples,
        options={},
    )


def _source(
    output_type: str,
    content: bytes,
    *,
    sample: str | None = None,
    suffix: str = "txt",
    relative_path: str | None = None,
) -> QcSourceDocument:
    token = output_type.replace(".", "-")
    sample_token = sample or "run"
    if relative_path is None and output_type.startswith("bulk_rnaseq.fastqc."):
        _prefix, _tool, stage, role, _extension = output_type.split(".")
        number = {"single": "", "read1": "_1", "read2": "_2"}[role]
        if stage == "raw":
            root = f"{sample}_raw{number}_fastqc"
            directory = "raw"
        elif stage == "filtered":
            root = f"{sample}_filtered{number}_fastqc"
            directory = "filtered"
        elif role == "single":
            root = f"{sample}_trimmed_trimmed_fastqc"
            directory = "trim"
        else:
            index = number.removeprefix("_")
            root = f"{sample}_trimmed_{index}_val_{index}_fastqc"
            directory = "trim"
        relative_path = f"results/fastqc/{directory}/{root}.zip"
    return QcSourceDocument(
        source=QcSourceArtifact(
            artifact_id=f"artifact-{token}-{sample_token}",
            output_type=output_type,
            relative_path=(
                relative_path
                if relative_path is not None
                else f"results/qc/{sample_token}-{token}.{suffix}"
            ),
            metadata={
                "scope": "run" if sample is None else "sample",
                **({"sample_id": sample} if sample is not None else {}),
                "assay": "bulk-rnaseq",
            },
        ),
        content=content,
    )


def _star_log() -> bytes:
    return b"""\
Number of input reads | 1000
Uniquely mapped reads number | 800
Uniquely mapped reads % | 80.00%
Number of reads mapped to multiple loci | 100
% of reads mapped to multiple loci | 10.00%
% of reads mapped to too many loci | 2.00%
% of reads unmapped: too many mismatches | 1.00%
% of reads unmapped: too short | 6.00%
% of reads unmapped: other | 1.00%
"""


def _salmon_meta(*, percent: object = 75.25) -> bytes:
    return json.dumps(
        {
            "salmon_version": "1.10.3",
            "mapping_type": "alignment",
            "num_libraries": 1,
            "num_processed": 1000,
            "num_mapped": 752,
            "percent_mapped": percent,
        },
        separators=(",", ":"),
    ).encode()


def _core_sources(sample: str = "S1") -> list[QcSourceDocument]:
    return [
        _source("bulk_rnaseq.star.log_final", _star_log(), sample=sample),
        _source(
            "bulk_rnaseq.salmon.meta_info",
            _salmon_meta(),
            sample=sample,
            suffix="json",
        ),
    ]


def _featurecounts_summary(
    *,
    bam: str = "S1.sorted.bam",
    omit: str | None = None,
    extra: str | None = None,
) -> bytes:
    counts = {
        "Assigned": 750,
        "Unassigned_Unmapped": 100,
        "Unassigned_Read_Type": 0,
        "Unassigned_Singleton": 0,
        "Unassigned_MappingQuality": 0,
        "Unassigned_Chimera": 0,
        "Unassigned_FragmentLength": 0,
        "Unassigned_Duplicate": 0,
        "Unassigned_MultiMapping": 150,
        "Unassigned_Secondary": 0,
        "Unassigned_NonSplit": 0,
        "Unassigned_NoFeatures": 0,
        "Unassigned_Overlapping_Length": 0,
        "Unassigned_Ambiguity": 0,
    }
    if omit is not None:
        counts.pop(omit)
    if extra is not None:
        counts[extra] = 1
    lines = [f"Status\t{bam}", *(f"{key}\t{value}" for key, value in counts.items())]
    return ("\n".join(lines) + "\n").encode()


def _fastqc_zip(
    *,
    gc: str = "42",
    adapter_status: str = "WARN",
    data_adapter_status: str | None = None,
    adapter_fraction: str = "2.50",
    total_sequences: str = "1000",
    root: str = "S1_raw_fastqc",
    filename: str = "S1_raw.fastq.gz",
    extra_members: dict[str, bytes] | None = None,
) -> bytes:
    summary = "\n".join(
        (
            f"PASS\tBasic Statistics\t{filename}",
            f"PASS\tPer base sequence quality\t{filename}",
            f"WARN\tPer sequence quality scores\t{filename}",
            f"PASS\tPer sequence GC content\t{filename}",
            f"{adapter_status}\tAdapter Content\t{filename}",
            "",
        )
    ).encode()
    data = "\n".join(
        (
            "##FastQC\t0.12.1",
            ">>Basic Statistics\tpass",
            "#Measure\tValue",
            f"Filename\t{filename}",
            f"Total Sequences\t{total_sequences}",
            f"%GC\t{gc}",
            ">>END_MODULE",
            ">>Per base sequence quality\tpass",
            ">>END_MODULE",
            ">>Per sequence quality scores\twarn",
            ">>END_MODULE",
            ">>Per sequence GC content\tpass",
            ">>END_MODULE",
            f">>Adapter Content\t{(data_adapter_status or adapter_status).lower()}",
            "#Position\tIllumina Universal Adapter\tNextera Transposase Sequence",
            f"1\t0.00\t{adapter_fraction}",
            ">>END_MODULE",
            "",
        )
    ).encode()
    target = BytesIO()
    with ZipFile(target, "w", ZIP_DEFLATED) as archive:
        archive.writestr(f"{root}/summary.txt", summary)
        archive.writestr(f"{root}/fastqc_data.txt", data)
        for path, content in (extra_members or {}).items():
            archive.writestr(path, content)
    return target.getvalue()


def _metric_map(result) -> dict[str, object]:
    assert result.is_success, result.errors
    return {metric.metric_key: metric for metric in result.value}


def test_source_type_catalog_is_closed_versioned_and_deterministic():
    assert isinstance(BULK_RNASEQ_QC_SOURCE_TYPES, tuple)
    assert BULK_RNASEQ_QC_SOURCE_TYPES == tuple(sorted(BULK_RNASEQ_QC_SOURCE_TYPES))
    assert len(BULK_RNASEQ_QC_SOURCE_TYPES) == len(set(BULK_RNASEQ_QC_SOURCE_TYPES))
    assert "bulk_rnaseq.multiqc.report" not in BULK_RNASEQ_QC_SOURCE_TYPES
    assert "bulk_rnaseq.multiqc.data_json" not in BULK_RNASEQ_QC_SOURCE_TYPES
    assert "bulk_rnaseq.multiqc.cutadapt" in BULK_RNASEQ_QC_SOURCE_TYPES
    assert not any("trimgalore" in value for value in BULK_RNASEQ_QC_SOURCE_TYPES)
    assert qc_module._MAX_SOURCE_FILES >= MAX_SAMPLE_ROWS * 14 + 4
    assert qc_module._MAX_METRICS >= MAX_SAMPLE_ROWS * 100


def test_fastqc_zip_extracts_decimal_metrics_flags_and_source_binding():
    inputs = _inputs(qc={"fastqc": True})
    fastqc = _source(
        "bulk_rnaseq.fastqc.raw.single.zip",
        _fastqc_zip(),
        sample="S1",
        suffix="zip",
    )

    metrics = _metric_map(
        extract_bulk_rnaseq_qc_metrics(inputs, tuple([*_core_sources(), fastqc]))
    )

    total = metrics["fastqc.raw.single.total_sequences"]
    assert str(total.value) == "1000"
    assert total.unit == "count"
    assert total.scope == "sample"
    assert total.sample_id == "S1"
    assert total.experiment_id is None
    assert total.assay == "bulk-rnaseq"
    assert total.qc_flag == "pass"
    assert total.source_artifact_id == fastqc.source.artifact_id
    assert str(metrics["fastqc.raw.single.gc_fraction"].value) == "0.42"
    assert (
        str(metrics["fastqc.raw.single.adapter_content_max_fraction"].value) == "0.025"
    )
    assert (
        metrics["fastqc.raw.single.adapter_content_max_fraction"].qc_flag == "warning"
    )
    assert str(metrics["fastqc.raw.single.per_sequence_quality_pass"].value) == "0"
    assert metrics["fastqc.raw.single.per_sequence_quality_pass"].qc_flag == "warning"


@pytest.mark.parametrize(
    ("content", "relative_path"),
    (
        (_fastqc_zip(filename="UNKNOWN.fastq.gz"), None),
        (_fastqc_zip(root="UNKNOWN_fastqc"), None),
        (
            _fastqc_zip(
                extra_members={"OTHER/private.txt": b"/private/workspace TOKEN=secret"}
            ),
            None,
        ),
        (
            _fastqc_zip(extra_members={"S1_raw_fastqc\\..\\private.txt": b"secret"}),
            None,
        ),
        (
            _fastqc_zip(),
            "results/fastqc/raw/UNKNOWN_raw_fastqc.zip",
        ),
    ),
)
def test_fastqc_identity_and_closed_archive_tree_fail_closed(
    content: bytes,
    relative_path: str | None,
):
    source = _source(
        "bulk_rnaseq.fastqc.raw.single.zip",
        content,
        sample="S1",
        suffix="zip",
        relative_path=relative_path,
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(qc={"fastqc": True}),
        tuple([*_core_sources(), source]),
    )

    assert result.is_failure
    assert result.errors[0].context["reason_code"] in {
        "source_contract_invalid",
        "source_content_invalid",
    }


def test_fastqc_summary_and_data_status_conflict_fails_closed():
    source = _source(
        "bulk_rnaseq.fastqc.raw.single.zip",
        _fastqc_zip(adapter_status="WARN", data_adapter_status="PASS"),
        sample="S1",
        suffix="zip",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(qc={"fastqc": True}),
        tuple([*_core_sources(), source]),
    )

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


def test_trimmed_fastqc_provides_exact_retained_read_counts_for_pe():
    inputs = _inputs(
        samples=[_sample("S1", layout="PE")],
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={"fastqc": True},
    )
    sources = [*_core_sources()]
    for stage, count in (("raw", "1000"), ("trimmed", "800")):
        for role in ("read1", "read2"):
            index = "1" if role == "read1" else "2"
            basename = (
                f"S1_raw_{index}"
                if stage == "raw"
                else f"S1_trimmed_{index}_val_{index}"
            )
            sources.append(
                _source(
                    f"bulk_rnaseq.fastqc.{stage}.{role}.zip",
                    _fastqc_zip(
                        total_sequences=count,
                        root=f"{basename}_fastqc",
                        filename=f"{basename}.{'fastq.gz' if stage == 'raw' else 'fq.gz'}",
                    ),
                    sample="S1",
                    suffix="zip",
                )
            )

    metrics = _metric_map(extract_bulk_rnaseq_qc_metrics(inputs, tuple(sources)))

    assert str(metrics["trimming.read1.retained_reads"].value) == "800"
    assert str(metrics["trimming.read2.retained_reads"].value) == "800"
    assert metrics["trimming.read1.retained_reads"].unit == "count"
    assert "trimming.read1.retained_bases" not in metrics


def test_trimmed_fastqc_rejects_mismatched_pe_retained_read_counts():
    inputs = _inputs(
        samples=[_sample("S1", layout="PE")],
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={"fastqc": True},
    )
    sources = [*_core_sources()]
    for stage, role, count in (
        ("raw", "read1", "1000"),
        ("raw", "read2", "1000"),
        ("trimmed", "read1", "800"),
        ("trimmed", "read2", "799"),
    ):
        index = "1" if role == "read1" else "2"
        basename = (
            f"S1_raw_{index}" if stage == "raw" else f"S1_trimmed_{index}_val_{index}"
        )
        sources.append(
            _source(
                f"bulk_rnaseq.fastqc.{stage}.{role}.zip",
                _fastqc_zip(
                    total_sequences=count,
                    root=f"{basename}_fastqc",
                    filename=f"{basename}.{'fastq.gz' if stage == 'raw' else 'fq.gz'}",
                ),
                sample="S1",
                suffix="zip",
                relative_path=(
                    f"results/fastqc/{'raw' if stage == 'raw' else 'trim'}/"
                    f"{basename}_fastqc.zip"
                ),
            )
        )

    result = extract_bulk_rnaseq_qc_metrics(inputs, tuple(sources))

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


def test_pe_umi_discard_read_uses_single_end_downstream_qc_grammar():
    inputs = _inputs(
        samples=[_sample("S1", layout="PE")],
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={"fastqc": True, "rseqc": True},
        umi={
            "enabled": True,
            "mode": "read_sequence",
            "deduplication_tool": "umitools",
            "extraction_method": "string",
            "barcode_pattern": "NNNN",
            "discard_read": 2,
        },
        advanced={"rseqc_modules": "infer_experiment"},
    )
    sources = [*_core_sources()]
    for role, index in (("read1", "1"), ("read2", "2")):
        basename = f"S1_raw_{index}"
        sources.append(
            _source(
                f"bulk_rnaseq.fastqc.raw.{role}.zip",
                _fastqc_zip(
                    root=f"{basename}_fastqc",
                    filename=f"{basename}.fastq.gz",
                ),
                sample="S1",
                suffix="zip",
            )
        )
    sources.extend(
        (
            _source(
                "bulk_rnaseq.fastqc.trimmed.single.zip",
                _fastqc_zip(
                    root="S1_trimmed_trimmed_fastqc",
                    filename="S1_trimmed_trimmed.fq.gz",
                    total_sequences="800",
                ),
                sample="S1",
                suffix="zip",
            ),
            _source(
                "bulk_rnaseq.rseqc.infer_experiment",
                b"\n\nThis is SingleEnd Data\n"
                b"Fraction of reads failed to determine: 0.0500\n"
                b'Fraction of reads explained by "++,--": 0.8000\n'
                b'Fraction of reads explained by "+-,-+": 0.1500\n',
                sample="S1",
            ),
        )
    )

    metrics = _metric_map(extract_bulk_rnaseq_qc_metrics(inputs, tuple(sources)))

    assert "fastqc.raw.read1.total_sequences" in metrics
    assert "fastqc.raw.read2.total_sequences" in metrics
    assert str(metrics["trimming.single.retained_reads"].value) == "800"
    assert str(metrics["rseqc.infer_experiment.orientation_a_fraction"].value) == "0.8"


def test_fastp_low_trim_sample_omits_trimmed_fastqc_without_inventing_zero():
    inputs = _inputs(
        trimming={"enabled": True, "tool": "fastp"},
        qc={"fastqc": True, "multiqc": True},
    )
    raw_fastqc = _source(
        "bulk_rnaseq.fastqc.raw.single.zip",
        _fastqc_zip(),
        sample="S1",
        suffix="zip",
    )
    fastp = _source(
        "bulk_rnaseq.trim.fastp.json",
        json.dumps(
            {
                "fastp_version": "1.0.1",
                "summary": {
                    "before_filtering": {
                        "total_reads": 10_000,
                        "total_bases": 1_000_000,
                    },
                    "after_filtering": {
                        "total_reads": 9_999,
                        "total_bases": 900_000,
                    },
                },
                "filtering_result": {"passed_filter_reads": 9_999},
            }
        ).encode(),
        sample="S1",
        suffix="json",
    )
    status = _source(
        "bulk_rnaseq.multiqc.fail_trimmed_samples",
        b"Sample\tReads after trimming\nS1\t9999.0\n",
        suffix="tsv",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        inputs,
        (raw_fastqc, fastp, status),
    )

    metrics = _metric_map(result)
    assert str(metrics["trimming.retained_reads"].value) == "9999"
    assert not any(key.startswith("fastqc.trimmed") for key in metrics)


@pytest.mark.parametrize(
    ("layout", "umi", "downstream_roles", "retained"),
    (
        ("SE", None, ("single",), 150),
        ("PE", None, ("read1", "read2"), 300),
        (
            "PE",
            {
                "enabled": True,
                "mode": "read_sequence",
                "deduplication_tool": "umitools",
                "extraction_method": "string",
                "barcode_pattern": "NNNN",
                "discard_read": 2,
            },
            ("single",),
            150,
        ),
    ),
)
def test_fastp_and_trimmed_fastqc_retained_reads_reconcile_exactly(
    layout: str,
    umi: dict[str, object] | None,
    downstream_roles: tuple[str, ...],
    retained: int,
):
    inputs = _inputs(
        samples=[_sample("S1", layout=layout)],
        trimming={"enabled": True, "tool": "fastp"},
        qc={"fastqc": True},
        umi=umi,
    )
    sources = [*_core_sources()]
    raw_roles = ("single",) if layout == "SE" else ("read1", "read2")
    for role in raw_roles:
        number = {"single": "", "read1": "_1", "read2": "_2"}[role]
        basename = f"S1_raw{number}"
        sources.append(
            _source(
                f"bulk_rnaseq.fastqc.raw.{role}.zip",
                _fastqc_zip(
                    root=f"{basename}_fastqc",
                    filename=f"{basename}.fastq.gz",
                ),
                sample="S1",
                suffix="zip",
            )
        )
    per_file = retained // len(downstream_roles)
    for role in downstream_roles:
        number = {"single": "", "read1": "_1", "read2": "_2"}[role]
        basename = f"S1_trimmed{number}"
        sources.append(
            _source(
                f"bulk_rnaseq.fastqc.trimmed.{role}.zip",
                _fastqc_zip(
                    root=f"{basename}_fastqc",
                    filename=f"{basename}.fastq.gz",
                    total_sequences=str(per_file),
                ),
                sample="S1",
                suffix="zip",
                relative_path=f"results/fastqc/trim/{basename}_fastqc.zip",
            )
        )

    def fastp_source(retained_reads: int) -> QcSourceDocument:
        return _source(
            "bulk_rnaseq.trim.fastp.json",
            json.dumps(
                {
                    "fastp_version": "1.0.1",
                    "summary": {
                        "before_filtering": {
                            "total_reads": retained + 100,
                            "total_bases": 100_000,
                        },
                        "after_filtering": {
                            "total_reads": retained_reads,
                            "total_bases": 90_000,
                        },
                    },
                    "filtering_result": {"passed_filter_reads": retained_reads},
                }
            ).encode(),
            sample="S1",
            suffix="json",
        )

    success = extract_bulk_rnaseq_qc_metrics(
        inputs,
        tuple([*sources, fastp_source(retained)]),
    )
    assert success.is_success

    conflict = extract_bulk_rnaseq_qc_metrics(
        inputs,
        tuple([*sources, fastp_source(retained - 1)]),
    )
    assert conflict.is_failure
    assert conflict.errors[0].context == {"reason_code": "source_content_invalid"}


def test_multiqc_cutadapt_pe_has_explicit_pre_filter_semantics():
    inputs = _inputs(
        samples=[_sample("S1", layout="PE")],
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={"multiqc": True},
    )
    cutadapt = _source(
        "bulk_rnaseq.multiqc.cutadapt",
        b"Sample\tcutadapt_version\tr_processed\tr_with_adapters\tr_written\t"
        b"bp_processed\tquality_trimmed\tbp_written\tpercent_trimmed\n"
        b"S1_1\t4.0\t600\t60\t600\t60000\t1000\t54000\t10.0\n"
        b"S1_2\t4.0\t600\t120\t600\t60000\t2000\t48000\t20.0\n",
        suffix="tsv",
    )
    sources = tuple([*_core_sources(), cutadapt])

    first = extract_bulk_rnaseq_qc_metrics(inputs, sources)
    second = extract_bulk_rnaseq_qc_metrics(inputs, tuple(reversed(sources)))

    assert first == second
    metrics = _metric_map(first)
    assert str(metrics["trimming.read1.input_reads"].value) == "600"
    assert str(metrics["trimming.read1.adapter_affected_fraction"].value) == "0.1"
    assert str(metrics["trimming.read1.post_adapter_quality_bases"].value) == "54000"
    assert (
        str(metrics["trimming.read1.post_adapter_quality_base_fraction"].value) == "0.9"
    )
    assert str(metrics["trimming.read2.input_bases"].value) == "60000"
    assert str(metrics["trimming.read2.quality_trimmed_bases"].value) == "2000"
    assert (
        str(metrics["trimming.read2.post_adapter_quality_base_fraction"].value) == "0.8"
    )
    assert (
        metrics["trimming.read1.input_reads"].source_artifact_id
        == cutadapt.source.artifact_id
    )
    assert "trimming.read1.retained_reads" not in metrics


@pytest.mark.parametrize(
    "rows",
    (
        b"UNKNOWN\t4.0\t600\t60\t600\t60000\t1000\t54000\t10.0\n"
        b"S1_2\t4.0\t600\t120\t600\t60000\t2000\t48000\t20.0\n",
        b"S1_1\t4.0\t600\t60\t599\t60000\t1000\t54000\t10.0\n"
        b"S1_2\t4.0\t600\t120\t600\t60000\t2000\t48000\t20.0\n",
        b"S1_1\t4.0\t600\t60\t600\t60000\t1000\t54000\t11.0\n"
        b"S1_2\t4.0\t600\t120\t600\t60000\t2000\t48000\t20.0\n",
        b"S1_1\t4.0\t600\t60\t600\t60000\t1000\t54000\t10.0\n"
        b"S1_2\t4.0\t601\t120\t601\t60000\t2000\t48000\t20.0\n",
    ),
)
def test_multiqc_cutadapt_corrupt_or_unknown_rows_fail_closed(rows: bytes):
    inputs = _inputs(
        samples=[_sample("S1", layout="PE")],
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={"multiqc": True},
    )
    source = _source(
        "bulk_rnaseq.multiqc.cutadapt",
        b"Sample\tcutadapt_version\tr_processed\tr_with_adapters\tr_written\t"
        b"bp_processed\tquality_trimmed\tbp_written\tpercent_trimmed\n" + rows,
        suffix="tsv",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        inputs,
        tuple([*_core_sources(), source]),
    )

    assert result.is_failure
    assert result.errors[0].context["reason_code"] in {
        "source_contract_invalid",
        "source_content_invalid",
        "metric_value_invalid",
    }


def test_fastp_mixed_layout_metrics_are_not_double_counted():
    inputs = _inputs(
        samples=[_sample("S1", layout="SE"), _sample("S2", layout="PE")],
        trimming={"enabled": True, "tool": "fastp"},
    )
    sources = [*_core_sources("S1"), *_core_sources("S2")]
    for sample in ("S1", "S2"):
        payload = {
            "fastp_version": "1.0.1",
            "summary": {
                "before_filtering": {"total_reads": 200, "total_bases": 20000},
                "after_filtering": {"total_reads": 150, "total_bases": 14000},
            },
            "filtering_result": {"passed_filter_reads": 150},
        }
        sources.append(
            _source(
                "bulk_rnaseq.trim.fastp.json",
                json.dumps(payload).encode(),
                sample=sample,
                suffix="json",
            )
        )

    result = extract_bulk_rnaseq_qc_metrics(inputs, tuple(sources))

    assert result.is_success
    by_sample = {
        (metric.sample_id, metric.metric_key): metric for metric in result.value
    }
    assert str(by_sample[("S1", "trimming.input_reads")].value) == "200"
    assert str(by_sample[("S2", "trimming.input_reads")].value) == "200"
    assert str(by_sample[("S2", "trimming.read_retained_fraction")].value) == "0.75"


def test_star_and_salmon_semantics_are_distinct_and_exact_decimal():
    result = extract_bulk_rnaseq_qc_metrics(_inputs(), tuple(_core_sources()))

    metrics = _metric_map(result)
    assert str(metrics["star.uniquely_mapped_fraction"].value) == "0.8"
    assert str(metrics["star.multimapped_fraction"].value) == "0.1"
    assert str(metrics["star.unmapped_fraction"].value) == "0.1"
    assert str(metrics["salmon.mapping_fraction"].value) == "0.7525"
    assert (
        metrics["star.uniquely_mapped_fraction"].display_name
        != metrics["salmon.mapping_fraction"].display_name
    )


def test_numeric_leading_authoring_sample_identity_is_valid_for_qc():
    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(samples=[_sample("1sample")]),
        tuple(_core_sources("1sample")),
    )

    metrics = _metric_map(result)
    assert {metric.sample_id for metric in metrics.values()} == {"1sample"}


def test_decimal_results_ignore_the_ambient_context():
    sources = tuple(_core_sources())
    expected = extract_bulk_rnaseq_qc_metrics(_inputs(), sources)

    with localcontext() as context:
        context.prec = 6
        observed = extract_bulk_rnaseq_qc_metrics(_inputs(), sources)

    assert observed == expected


def test_decimal_exponent_bomb_fails_before_formatting():
    sources = _core_sources()
    sources[1] = _source(
        "bulk_rnaseq.salmon.meta_info",
        _salmon_meta(percent="1e1000000"),
        sample="S1",
        suffix="json",
    )

    result = extract_bulk_rnaseq_qc_metrics(_inputs(), tuple(sources))

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


def test_fixed_multiqc_status_tables_define_legitimate_filtered_samples():
    inputs = _inputs(
        samples=[_sample("S1"), _sample("S2")],
        qc={"multiqc": True},
    )
    trim_status = _source(
        "bulk_rnaseq.multiqc.fail_trimmed_samples",
        b"Sample\tReads after trimming\nS2\t9999.0\n",
        suffix="tsv",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        inputs,
        tuple([*_core_sources("S1"), trim_status]),
    )

    assert result.is_success
    assert {metric.sample_id for metric in result.value} == {"S1"}


def test_exact_trim_threshold_status_row_preserves_existing_analysis_metrics():
    status = _source(
        "bulk_rnaseq.multiqc.fail_trimmed_samples",
        b"Sample\tReads after trimming\nS1\t10000.0\n",
        suffix="tsv",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(qc={"multiqc": True}),
        tuple([*_core_sources(), status]),
    )

    metrics = _metric_map(result)
    assert "star.input_reads" in metrics
    assert "salmon.mapping_fraction" in metrics


def test_mapped_status_keeps_star_salmon_but_omits_downstream_qc():
    inputs = _inputs(
        samples=[_sample("S1"), _sample("S2")],
        qc={"multiqc": True},
    )
    mapped_status = _source(
        "bulk_rnaseq.multiqc.fail_mapped_samples",
        b"Sample\tSTAR uniquely mapped reads (%)\nS2\t4.5\n",
        suffix="tsv",
    )
    sources = tuple([*_core_sources("S1"), *_core_sources("S2"), mapped_status])

    result = extract_bulk_rnaseq_qc_metrics(inputs, sources)

    assert result.is_success
    assert {metric.sample_id for metric in result.value} == {"S1", "S2"}
    assert all(
        not metric.metric_key.startswith(("picard.", "rseqc.", "featurecounts."))
        for metric in result.value
    )


@pytest.mark.parametrize(
    ("advanced", "output_type", "content"),
    (
        (
            {"min_trimmed_reads": 100},
            "bulk_rnaseq.multiqc.fail_trimmed_samples",
            b"Sample\tReads after trimming\nS1\t101\n",
        ),
        (
            {"min_mapped_reads": 10},
            "bulk_rnaseq.multiqc.fail_mapped_samples",
            b"Sample\tSTAR uniquely mapped reads (%)\nS1\t10\n",
        ),
    ),
)
def test_status_table_must_prove_the_configured_filter_threshold(
    advanced: dict[str, object],
    output_type: str,
    content: bytes,
):
    status = _source(output_type, content, suffix="tsv")

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(
            qc={"multiqc": True},
            advanced=advanced,
            trimming=(
                {"enabled": True, "tool": "trimgalore"}
                if "min_trimmed_reads" in advanced
                else None
            ),
        ),
        (status,),
    )

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_contract_invalid"}


@pytest.mark.parametrize(
    "content",
    (
        b"Sample\tReads after trimming\nUNKNOWN\t1\n",
        b"Sample\tReads after trimming\nS1\t1\nS1\t2\n",
        b"Sample\tWrong\nS1\t1\n",
    ),
)
def test_status_table_unknown_duplicate_and_bad_header_fail_closed(content):
    status = _source(
        "bulk_rnaseq.multiqc.fail_trimmed_samples",
        content,
        suffix="tsv",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(qc={"multiqc": True}),
        tuple([status]),
    )

    assert result.is_failure
    assert result.errors[0].context["reason_code"] in {
        "source_contract_invalid",
        "source_content_invalid",
    }


def test_featurecounts_picard_and_rseqc_use_closed_machine_sources():
    inputs = _inputs(
        samples=[_sample("S1", layout="PE")],
        qc={
            "multiqc": True,
            "mark_duplicates": True,
            "biotype": True,
            "rseqc": True,
        },
        advanced={"rseqc_modules": "bam_stat,infer_experiment,read_distribution,tin"},
    )
    sources = [*_core_sources()]
    sources.extend(
        (
            _source(
                "bulk_rnaseq.featurecounts.summary",
                _featurecounts_summary(bam="S1.markdup.sorted.bam"),
                sample="S1",
            ),
            _source(
                "bulk_rnaseq.multiqc.picard_dups",
                b"Sample\tLIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\t"
                b"SECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\t"
                b"UNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t"
                b"READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\t"
                b"ESTIMATED_LIBRARY_SIZE\n"
                b"S1\tUnknown Library\t0\t500\t0\t0\t0\t50\t5\t0.1\t7000.0\n",
                suffix="tsv",
            ),
            _source(
                "bulk_rnaseq.rseqc.bam_stat",
                b"#Output (all numbers are read count)\n#================================\n"
                b"Total records: 1000\nQC failed: 10\n"
                b"Optical/PCR duplicate: 50\nNon primary hits 100\n"
                b"Unmapped reads: 100\nmapq < mapq_cut (non-unique): 140\n"
                b"mapq >= mapq_cut (unique): 600\nRead-1: 300\nRead-2: 300\n"
                b"Reads map to '+': 310\nReads map to '-': 290\n"
                b"Non-splice reads: 500\nSplice reads: 100\n"
                b"Reads mapped in proper pairs: 550\n"
                b"Proper-paired reads map to different chrom: 20\n",
                sample="S1",
            ),
            _source(
                "bulk_rnaseq.rseqc.infer_experiment",
                b"This is PairEnd Data\nFraction of reads failed to determine: 0.05\n"
                b'Fraction of reads explained by "1++,1--,2+-,2-+": 0.80\n'
                b'Fraction of reads explained by "1+-,1-+,2++,2--": 0.15\n',
                sample="S1",
            ),
            _source(
                "bulk_rnaseq.rseqc.read_distribution",
                b"Total Reads                   1000\n"
                b"Total Tags                    1200\n"
                b"Total Assigned Tags           800\n"
                b"============================================================\n"
                b"Group               Total_bases         Tag_count           Tags/Kb\n"
                b"CDS_Exons           999                 600                 600.00\n"
                b"5'UTR_Exons         99                  50                  500.00\n"
                b"3'UTR_Exons         99                  40                  400.00\n"
                b"Introns             499                 50                  100.00\n"
                b"TSS_up_1kb          999                 10                  10.00\n"
                b"TSS_up_5kb          999                 20                  20.00\n"
                b"TSS_up_10kb         999                 30                  30.00\n"
                b"TES_down_1kb        999                 10                  10.00\n"
                b"TES_down_5kb        999                 20                  20.00\n"
                b"TES_down_10kb       999                 30                  30.00\n"
                b"============================================================\n",
                sample="S1",
            ),
            _source(
                "bulk_rnaseq.rseqc.tin_summary",
                b"Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\n"
                b"S1.markdup.sorted.bam\t72.125\t74.5\t8.25\n",
                sample="S1",
                relative_path="results/star_salmon/rseqc/tin/S1.summary.txt",
            ),
        )
    )

    metrics = _metric_map(extract_bulk_rnaseq_qc_metrics(inputs, tuple(sources)))

    assert str(metrics["featurecounts.assigned_fraction"].value) == "0.75"
    assert str(metrics["picard.duplication_fraction"].value) == "0.1"
    assert str(metrics["picard.estimated_library_size"].value) == "7000"
    assert str(metrics["rseqc.bam_stat.unique_fraction"].value) == "0.6"
    assert str(metrics["rseqc.infer_experiment.orientation_a_fraction"].value) == "0.8"
    assert str(metrics["rseqc.read_distribution.cds_exon_tags"].value) == "600"
    assert (
        str(metrics["rseqc.read_distribution.cds_exon_assigned_fraction"].value)
        == "0.75"
    )
    assert str(metrics["rseqc.tin.mean_score"].value) == "72.125"
    assert str(metrics["rseqc.tin.median_score"].value) == "74.5"
    assert str(metrics["rseqc.tin.standard_deviation"].value) == "8.25"
    assert metrics["rseqc.tin.mean_score"].unit == "score"
    assert (
        metrics["rseqc.tin.mean_score"].source_artifact_id
        == "artifact-bulk_rnaseq-rseqc-tin_summary-S1"
    )


def test_picard_duplication_reconciles_counts_with_six_decimal_serialization():
    header = (
        b"Sample\tLIBRARY\tUNPAIRED_READS_EXAMINED\tREAD_PAIRS_EXAMINED\t"
        b"SECONDARY_OR_SUPPLEMENTARY_RDS\tUNMAPPED_READS\t"
        b"UNPAIRED_READ_DUPLICATES\tREAD_PAIR_DUPLICATES\t"
        b"READ_PAIR_OPTICAL_DUPLICATES\tPERCENT_DUPLICATION\t"
        b"ESTIMATED_LIBRARY_SIZE\n"
    )
    valid = _source(
        "bulk_rnaseq.multiqc.picard_dups",
        header + b"S1\tlib1\t3\t0\t0\t0\t1\t0\t0\t0.333333\t-\n",
        suffix="tsv",
    )
    inputs = _inputs(qc={"multiqc": True, "mark_duplicates": True})

    metrics = _metric_map(
        extract_bulk_rnaseq_qc_metrics(inputs, tuple([*_core_sources(), valid]))
    )
    assert str(metrics["picard.duplication_fraction"].value) == "0.333333"

    contradictory = _source(
        "bulk_rnaseq.multiqc.picard_dups",
        header + b"S1\tlib1\t3\t0\t0\t0\t1\t0\t0\t0.333334\t-\n",
        suffix="tsv",
    )
    result = extract_bulk_rnaseq_qc_metrics(
        inputs,
        tuple([*_core_sources(), contradictory]),
    )
    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


def test_csi_profile_omits_fixed_incompatible_rseqc_qc_sources():
    infer = _source(
        "bulk_rnaseq.rseqc.infer_experiment",
        b"\n\nThis is SingleEnd Data\n"
        b"Fraction of reads failed to determine: 0.0500\n"
        b'Fraction of reads explained by "++,--": 0.8000\n'
        b'Fraction of reads explained by "+-,-+": 0.1500\n',
        sample="S1",
    )
    inputs = _inputs(
        qc={"rseqc": True},
        advanced={
            "bam_csi_index": True,
            "rseqc_modules": "infer_experiment,read_distribution,tin",
        },
    )

    result = extract_bulk_rnaseq_qc_metrics(
        inputs,
        tuple([*_core_sources(), infer]),
    )

    metrics = _metric_map(result)
    assert "rseqc.infer_experiment.orientation_a_fraction" in metrics
    assert not any(key.startswith("rseqc.tin") for key in metrics)
    assert not any(key.startswith("rseqc.read_distribution") for key in metrics)


@pytest.mark.parametrize(
    "content",
    (
        b"Wrong\tTIN(mean)\tTIN(median)\tTIN(stdev)\nS1.sorted.bam\t1\t1\t1\n",
        b"Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\nUNKNOWN.bam\t1\t1\t1\n",
        b"Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\nS1.sorted.bam\t101\t1\t1\n",
        b"Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\nS1.sorted.bam\t1\t-1\t1\n",
        b"Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\nS1.sorted.bam\t1\t1\t101\n",
        b"Bam_file\tTIN(mean)\tTIN(median)\tTIN(stdev)\nS1.sorted.bam\tNaN\t1\t1\n",
    ),
)
def test_rseqc_tin_uses_exact_fixed_grammar_and_score_bounds(content: bytes):
    source = _source(
        "bulk_rnaseq.rseqc.tin_summary",
        content,
        sample="S1",
        relative_path="results/star_salmon/rseqc/tin/S1.summary.txt",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(qc={"rseqc": True}, advanced={"rseqc_modules": "tin"}),
        tuple([*_core_sources(), source]),
    )

    assert result.is_failure
    assert result.errors[0].context["reason_code"] in {
        "source_content_invalid",
        "metric_value_invalid",
    }


def test_rseqc_se_bam_stat_accepts_fixed_zero_paired_fields():
    source = _source(
        "bulk_rnaseq.rseqc.bam_stat",
        b"#Output (all numbers are read count)\n#================================\n"
        b"Total records: 1000\nQC failed: 10\n"
        b"Optical/PCR duplicate: 50\nNon primary hits 100\n"
        b"Unmapped reads: 100\nmapq < mapq_cut (non-unique): 140\n"
        b"mapq >= mapq_cut (unique): 600\nRead-1: 0\nRead-2: 0\n"
        b"Reads map to '+': 310\nReads map to '-': 290\n"
        b"Non-splice reads: 500\nSplice reads: 100\n"
        b"Reads mapped in proper pairs: 0\n"
        b"Proper-paired reads map to different chrom: 0\n",
        sample="S1",
    )

    metrics = _metric_map(
        extract_bulk_rnaseq_qc_metrics(
            _inputs(qc={"rseqc": True}, advanced={"rseqc_modules": "bam_stat"}),
            tuple([*_core_sources(), source]),
        )
    )

    assert str(metrics["rseqc.bam_stat.unique_fraction"].value) == "0.6"
    assert "rseqc.bam_stat.proper_pair_fraction" not in metrics


def test_rseqc_bam_stat_rejects_contradictory_unique_subtotals():
    source = _source(
        "bulk_rnaseq.rseqc.bam_stat",
        b"#Output (all numbers are read count)\n#================================\n"
        b"Total records: 1000\nQC failed: 10\n"
        b"Optical/PCR duplicate: 50\nNon primary hits 100\n"
        b"Unmapped reads: 100\nmapq < mapq_cut (non-unique): 140\n"
        b"mapq >= mapq_cut (unique): 600\nRead-1: 0\nRead-2: 0\n"
        b"Reads map to '+': 309\nReads map to '-': 290\n"
        b"Non-splice reads: 500\nSplice reads: 100\n"
        b"Reads mapped in proper pairs: 0\n"
        b"Proper-paired reads map to different chrom: 0\n",
        sample="S1",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(qc={"rseqc": True}, advanced={"rseqc_modules": "bam_stat"}),
        tuple([*_core_sources(), source]),
    )

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


def test_rseqc_read_distribution_reconciles_outer_windows_and_total():
    source = _source(
        "bulk_rnaseq.rseqc.read_distribution",
        b"Total Reads                   1000\n"
        b"Total Tags                    1200\n"
        b"Total Assigned Tags           800\n"
        b"============================================================\n"
        b"Group               Total_bases         Tag_count           Tags/Kb\n"
        b"CDS_Exons           999                 600                 600.00\n"
        b"5'UTR_Exons         99                  50                  500.00\n"
        b"3'UTR_Exons         99                  40                  400.00\n"
        b"Introns             499                 50                  100.00\n"
        b"TSS_up_1kb          999                 10                  10.00\n"
        b"TSS_up_5kb          999                 20                  20.00\n"
        b"TSS_up_10kb         999                 31                  31.00\n"
        b"TES_down_1kb        999                 10                  10.00\n"
        b"TES_down_5kb        999                 20                  20.00\n"
        b"TES_down_10kb       999                 30                  30.00\n"
        b"============================================================\n",
        sample="S1",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(
            qc={"rseqc": True},
            advanced={"rseqc_modules": "read_distribution"},
        ),
        tuple([*_core_sources(), source]),
    )

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


@pytest.mark.parametrize(
    "content",
    (
        _featurecounts_summary(omit="Unassigned_Ambiguity"),
        _featurecounts_summary(extra="Unassigned_Unknown"),
    ),
)
def test_featurecounts_requires_the_fixed_status_catalog(content: bytes):
    source = _source(
        "bulk_rnaseq.featurecounts.summary",
        content,
        sample="S1",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(qc={"biotype": True}),
        tuple([*_core_sources(), source]),
    )

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


def test_rseqc_infer_rejects_unknown_orientation_labels():
    source = _source(
        "bulk_rnaseq.rseqc.infer_experiment",
        b"This is SingleEnd Data\n"
        b"Fraction of reads failed to determine: 0.05\n"
        b'Fraction of reads explained by "++,--": 0.80\n'
        b'Fraction of reads explained by "invented": 0.15\n',
        sample="S1",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(qc={"rseqc": True}, advanced={"rseqc_modules": "infer_experiment"}),
        tuple([*_core_sources(), source]),
    )

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


@pytest.mark.parametrize(
    ("mutator", "expected_path"),
    (
        (lambda sources: sources[:-1], "qc_sources"),
        (
            lambda sources: [
                *sources,
                _source(
                    "bulk_rnaseq.salmon.meta_info",
                    _salmon_meta(),
                    sample="UNKNOWN",
                    suffix="json",
                ),
            ],
            "qc_sources",
        ),
        (lambda sources: [*sources, sources[0]], "qc_sources"),
    ),
)
def test_missing_unknown_and_duplicate_sources_fail_closed(mutator, expected_path):
    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(),
        tuple(mutator(_core_sources())),
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_QC_INVALID"
    assert result.errors[0].path == expected_path
    assert result.errors[0].context == {"reason_code": "source_contract_invalid"}


@pytest.mark.parametrize(
    "content",
    (
        b'{"salmon_version":"1.10.3","num_processed":1,"num_processed":2}',
        b'{"salmon_version":"1.10.3","percent_mapped":NaN}',
        b'{"salmon_version":"1.10.3","percent_mapped":Infinity}',
        b'{"a":{"b":{"c":{"d":{"e":{"f":{"g":{"h":{"i":{"j":{"k":{"l":{"m":{"n":{"o":{"p":{"q":1}}}}}}}}}}}}}}}}}',
    ),
)
def test_salmon_json_duplicate_non_finite_and_depth_fail_closed(content):
    sources = _core_sources()
    sources[1] = _source(
        "bulk_rnaseq.salmon.meta_info",
        content,
        sample="S1",
        suffix="json",
    )

    result = extract_bulk_rnaseq_qc_metrics(_inputs(), tuple(sources))

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


def test_deep_multiline_json_fails_with_a_stable_redacted_reason():
    content = ("[\n" * 9_999 + "0\n" + "]\n" * 9_999).encode()
    sources = _core_sources()
    sources[1] = _source(
        "bulk_rnaseq.salmon.meta_info",
        content,
        sample="S1",
        suffix="json",
    )

    result = extract_bulk_rnaseq_qc_metrics(_inputs(), tuple(sources))

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


def test_out_of_range_fraction_and_count_fail_closed():
    sources = _core_sources()
    sources[1] = _source(
        "bulk_rnaseq.salmon.meta_info",
        _salmon_meta(percent=125),
        sample="S1",
        suffix="json",
    )

    result = extract_bulk_rnaseq_qc_metrics(_inputs(), tuple(sources))

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "metric_value_invalid"}


@pytest.mark.parametrize("archive_case", ("traversal", "symlink", "bomb", "corrupt"))
def test_fastqc_archive_limits_and_paths_fail_closed(archive_case: str):
    if archive_case == "corrupt":
        content = b"not-a-zip"
    else:
        target = BytesIO()
        with ZipFile(target, "w", ZIP_DEFLATED) as archive:
            if archive_case == "traversal":
                archive.writestr("../summary.txt", b"PASS\tBasic Statistics\tx\n")
            elif archive_case == "symlink":
                info = ZipInfo("S1_fastqc/summary.txt")
                info.external_attr = (stat.S_IFLNK | 0o777) << 16
                archive.writestr(info, b"target")
            else:
                archive.writestr("S1_fastqc/huge.txt", b"0" * (65 * 1024 * 1024))
        content = target.getvalue()
    fastqc = _source(
        "bulk_rnaseq.fastqc.raw.single.zip",
        content,
        sample="S1",
        suffix="zip",
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(qc={"fastqc": True}),
        tuple([*_core_sources(), fastqc]),
    )

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_content_invalid"}


def test_source_documents_never_accept_absolute_path_metadata_or_dynamic_type():
    source = _core_sources()[0]
    unsafe = QcSourceDocument(
        source=QcSourceArtifact(
            artifact_id=source.source.artifact_id,
            output_type="bulk_rnaseq.dynamic.metric",
            relative_path=str(PurePosixPath("results/qc/unsafe.txt")),
            metadata={"scope": "sample", "sample_id": "S1", "assay": "/private"},
        ),
        content=source.content,
    )

    result = extract_bulk_rnaseq_qc_metrics(
        _inputs(),
        tuple([unsafe, _core_sources()[1]]),
    )

    assert result.is_failure
    assert result.errors[0].context == {"reason_code": "source_contract_invalid"}
