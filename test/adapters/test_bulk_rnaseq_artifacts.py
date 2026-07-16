"""Deterministic artifact discovery for pinned nf-core/rnaseq 3.26.0."""

from __future__ import annotations

from pathlib import Path
import os

import pytest

import encode_pipeline.adapters.bulk_rnaseq.artifacts as artifact_module
from encode_pipeline.adapters.bulk_rnaseq.artifacts import (
    discover_bulk_rnaseq_artifacts,
)
from encode_pipeline.platform.adapters import MAX_SAMPLE_ROWS, WorkflowInputs


_MERGED = (
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
)


def _inputs(
    *,
    layouts: tuple[str, ...] = ("SE",),
    repeated_lane: bool = False,
    trimming: dict[str, object] | None = None,
    qc: dict[str, object] | None = None,
    outputs: dict[str, object] | None = None,
    rrna: dict[str, object] | None = None,
    umi: dict[str, object] | None = None,
    advanced: dict[str, object] | None = None,
) -> WorkflowInputs:
    standard: dict[str, object] = {
        "reference": {
            "reference_id": "tiny",
            "fasta": "/inputs/ref.fa",
            "fasta_sha256": "1" * 64,
            "gtf": "/inputs/ref.gtf",
            "gtf_sha256": "2" * 64,
        },
        "trimming": trimming or {"enabled": False, "tool": "trimgalore"},
        "qc": qc or {"enabled": False},
        "outputs": {"bigwig": False, **(outputs or {})},
    }
    if rrna is not None:
        standard["ribosomal_rna_removal"] = rrna
    if umi is not None:
        standard["umi"] = umi
    rows: list[dict[str, str]] = []
    for index, layout in enumerate(layouts, start=1):
        sample = f"S{index}"
        row = {
            "sample": sample,
            "library": "lib1",
            "lane": "L001",
            "layout": layout,
            "fastq_1": (
                f"/inputs/{sample}_1.L001.fastq.gz"
                if layout == "PE"
                else f"/inputs/{sample}.L001.fastq.gz"
            ),
            "strandedness": "unstranded",
            "platform": "ILLUMINA",
        }
        if layout == "PE":
            row["fastq_2"] = f"/inputs/{sample}_2.L001.fastq.gz"
        rows.append(row)
        if repeated_lane and index == 1:
            repeated = dict(row)
            repeated["library"] = "lib2"
            repeated["lane"] = "L002"
            repeated["fastq_1"] = (
                f"/inputs/{sample}_1.L002.fastq.gz"
                if layout == "PE"
                else f"/inputs/{sample}.L002.fastq.gz"
            )
            if layout == "PE":
                repeated["fastq_2"] = f"/inputs/{sample}_2.L002.fastq.gz"
            rows.append(repeated)
    return WorkflowInputs(
        config={"standard": standard, "advanced": advanced or {}},
        samples=rows,
        options={},
    )


def _write(workspace: Path, relative: str, content: bytes = b"x") -> Path:
    path = workspace / "results" / relative
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(content)
    return path


def _write_core(
    workspace: Path, samples: tuple[str, ...], *, csi: bool = False
) -> None:
    _write(workspace, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    for sample in samples:
        for relative in (
            f"star_salmon/{sample}.sorted.bam",
            f"star_salmon/{sample}.sorted.bam.{'csi' if csi else 'bai'}",
            f"star_salmon/log/{sample}.Log.final.out",
            f"star_salmon/log/{sample}.SJ.out.tab",
            f"star_salmon/{sample}/quant.sf",
            f"star_salmon/{sample}/quant.genes.sf",
            f"star_salmon/{sample}/aux_info/meta_info.json",
        ):
            _write(workspace, relative)
    for suffix in _MERGED:
        _write(workspace, f"star_salmon/salmon.merged.{suffix}.tsv")


def _types(result) -> set[str]:
    assert result.is_success, [issue.to_dict() for issue in result.issues]
    return {candidate.output_type for candidate in result.value}


def _multiqc_only() -> dict[str, object]:
    return {
        "enabled": True,
        "fastqc": False,
        "multiqc": True,
        "rseqc": False,
        "qualimap": False,
        "dupradar": False,
        "biotype": False,
        "deseq2_pca": False,
        "preseq": False,
        "mark_duplicates": False,
    }


@pytest.mark.parametrize(
    ("layouts", "repeated_lane"),
    [(("SE",), False), (("PE",), True), (("SE", "PE"), False)],
)
def test_core_discovery_is_sample_scoped_sorted_and_deterministic(
    tmp_path: Path,
    layouts: tuple[str, ...],
    repeated_lane: bool,
):
    inputs = _inputs(layouts=layouts, repeated_lane=repeated_lane)
    samples = tuple(f"S{index}" for index in range(1, len(layouts) + 1))
    _write_core(tmp_path, samples)

    first = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    second = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert first == second
    assert first.is_success
    assert first.value == tuple(
        sorted(first.value, key=lambda item: (item.output_type, item.relative_path))
    )
    assert (
        len([item for item in first.value if item.metadata.get("sample_id") == "S1"])
        == len([item for item in first.value if item.metadata.get("sample_id") == "S2"])
        if len(layouts) == 2
        else True
    )
    assert all(item.relative_path.startswith("results/") for item in first.value)
    assert all(
        set(item.metadata) <= {"scope", "sample_id", "assay"} for item in first.value
    )
    assert not any("/inputs" in repr(item) for item in first.value)


def test_bounded_catalog_and_status_limits_cover_the_authoring_row_contract():
    assert artifact_module._MAX_ARTIFACT_CANDIDATES >= MAX_SAMPLE_ROWS * 128
    assert artifact_module._MAX_AUDITED_NAMESPACE_ENTRIES >= MAX_SAMPLE_ROWS * 128
    assert artifact_module._MAX_STATUS_TABLE_BYTES >= MAX_SAMPLE_ROWS * (
        128 + 1 + 64 + 1
    )


def test_required_missing_and_disabled_outputs_are_fail_closed(tmp_path: Path):
    inputs = _inputs()
    _write_core(tmp_path, ("S1",))
    (tmp_path / "results/star_salmon/S1/quant.sf").unlink()
    missing = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert missing.is_failure
    assert missing.errors[0].code == "BULK_RNASEQ_ARTIFACT_REQUIRED_MISSING"
    assert missing.errors[0].context == {"reason": "required_output_missing"}

    _write(tmp_path, "fastqc/raw/S1_raw_fastqc.zip")
    _write(tmp_path, "star_salmon/bigwig/S1.bigWig")
    _write(tmp_path, "not_a_sample/UNKNOWN.quant.sf")
    _write(tmp_path, "multiqc/star_salmon/multiqc_report_data/multiqc_data.json")
    _write(tmp_path, "multiqc/star_salmon/multiqc_report_data/multiqc.log")
    _write(tmp_path, "star_salmon/log/S1.Log.out")
    _write(tmp_path, "star_salmon/log/S1.Log.progress.out")
    _write(tmp_path, "star_salmon/S1/cmd_info.json")
    _write(tmp_path, "star_salmon/S1/logs/salmon_quant.log")
    _write_core(tmp_path, ("S1",))
    success = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    paths = {candidate.relative_path for candidate in success.value}
    assert "results/fastqc/raw/S1_raw_fastqc.zip" not in paths
    assert "results/star_salmon/bigwig/S1.bigWig" not in paths
    assert not any("UNKNOWN" in path for path in paths)
    assert not any(
        path.endswith(("multiqc_data.json", "multiqc.log")) for path in paths
    )
    assert not any(
        path.endswith(
            ("Log.out", "Log.progress.out", "cmd_info.json", "salmon_quant.log")
        )
        for path in paths
    )


def test_salmon_alignment_route_does_not_admit_pseudoalignment_lib_format(
    tmp_path: Path,
):
    inputs = _inputs()
    _write_core(tmp_path, ("S1",))
    _write(tmp_path, "star_salmon/S1/lib_format_counts.json")

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success
    assert "bulk_rnaseq.salmon.lib_format_counts" not in _types(result)


@pytest.mark.parametrize(
    "relative_path",
    (
        "star_salmon/UNKNOWN/quant.sf",
        "star_salmon/log/UNKNOWN.Log.final.out",
        "star_salmon/UNKNOWN.sorted.bam",
    ),
)
def test_unknown_sample_in_an_audited_namespace_fails_closed(
    tmp_path: Path,
    relative_path: str,
):
    inputs = _inputs()
    _write_core(tmp_path, ("S1",))
    _write(tmp_path, relative_path)

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_UNKNOWN_SAMPLE"
    assert result.errors[0].context == {"reason": "unknown_sample_output"}
    assert relative_path not in repr(result.errors[0].to_dict())


def test_dotted_sample_identity_is_matched_before_audited_suffixes(
    tmp_path: Path,
):
    inputs = _inputs()
    inputs.samples[0]["sample"] = "S1.batch"
    _write_core(tmp_path, ("S1.batch",))

    success = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert success.is_success
    assert {
        item.metadata.get("sample_id")
        for item in success.value
        if item.metadata.get("scope") == "sample"
    } == {"S1.batch"}

    _write(tmp_path, "star_salmon/UNKNOWN.sorted.bam")
    unknown = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert unknown.is_failure
    assert unknown.errors[0].code == "BULK_RNASEQ_ARTIFACT_UNKNOWN_SAMPLE"


def test_audited_namespace_enumeration_is_bounded(tmp_path: Path, monkeypatch):
    inputs = _inputs()
    _write_core(tmp_path, ("S1",))
    monkeypatch.setattr(artifact_module, "_MAX_AUDITED_NAMESPACE_ENTRIES", 0)

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_LIMIT_EXCEEDED"
    assert result.errors[0].context == {"reason": "artifact_file_limit"}


def test_default_star_salmon_profile_builds_a_closed_catalog_before_io(
    tmp_path: Path,
):
    inputs = WorkflowInputs(
        config={
            "standard": {
                "reference": {
                    "reference_id": "tiny",
                    "fasta": "/inputs/ref.fa",
                    "fasta_sha256": "1" * 64,
                    "gtf": "/inputs/ref.gtf",
                    "gtf_sha256": "2" * 64,
                }
            }
        },
        samples=[
            {
                "sample": "S1",
                "library": "lib1",
                "lane": "L001",
                "layout": "PE",
                "fastq_1": "/inputs/S1_1.fastq.gz",
                "fastq_2": "/inputs/S1_2.fastq.gz",
                "strandedness": "auto",
                "platform": "ILLUMINA",
            }
        ],
        options={},
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_REQUIRED_MISSING"
    assert result.errors[0].context == {"reason": "required_output_missing"}


def test_default_catalog_excludes_private_or_unadmitted_qc_namespaces():
    inputs = WorkflowInputs(
        config={
            "standard": {
                "reference": {
                    "reference_id": "tiny",
                    "fasta": "/inputs/ref.fa",
                    "fasta_sha256": "1" * 64,
                    "gtf": "/inputs/ref.gtf",
                    "gtf_sha256": "2" * 64,
                }
            }
        },
        samples=[
            {
                "sample": "S1",
                "library": "lib1",
                "lane": "L001",
                "layout": "SE",
                "fastq_1": "/inputs/S1.fastq.gz",
                "strandedness": "unstranded",
                "platform": "ILLUMINA",
            }
        ],
        options={},
    )
    validated = artifact_module.validate_bulk_rnaseq_inputs(inputs)
    assert validated.is_success
    normalized = validated.value
    specs, _auto = artifact_module._expected_specs(
        normalized["nfcore_params"],
        artifact_module._normalized_samples(normalized["samples"]),
        trimmed_failed=frozenset(),
        mapped_failed=frozenset(),
    )

    types = {spec.output_type for spec in specs}
    assert not any(
        output_type.startswith(
            (
                "bulk_rnaseq.qualimap.",
                "bulk_rnaseq.dupradar.",
                "bulk_rnaseq.deseq2_qc.",
                "bulk_rnaseq.preseq.",
            )
        )
        for output_type in types
    )
    assert "bulk_rnaseq.rseqc.tin" not in types
    assert "bulk_rnaseq.rseqc.tin_summary" not in types
    assert "bulk_rnaseq.star.log" not in types
    assert "bulk_rnaseq.salmon.log" not in types
    assert "bulk_rnaseq.trim.fastp.log" not in types


def test_reversing_sample_and_lane_order_does_not_change_candidates(tmp_path: Path):
    forward = _inputs(layouts=("PE", "SE"), repeated_lane=True)
    reversed_inputs = WorkflowInputs(
        config=forward.config,
        samples=list(reversed(forward.samples)),
        options=forward.options,
    )
    _write_core(tmp_path, ("S1", "S2"))

    first = discover_bulk_rnaseq_artifacts(forward, tmp_path)
    second = discover_bulk_rnaseq_artifacts(reversed_inputs, tmp_path)

    assert first.is_success
    assert first == second


def test_fastqc_trimgalore_and_machine_multiqc_types_are_exact(tmp_path: Path):
    qc = {
        "enabled": True,
        "fastqc": True,
        "multiqc": True,
        "rseqc": False,
        "qualimap": False,
        "dupradar": False,
        "biotype": False,
        "deseq2_pca": False,
        "preseq": False,
        "mark_duplicates": False,
    }
    inputs = _inputs(
        layouts=("PE",), trimming={"enabled": True, "tool": "trimgalore"}, qc=qc
    )
    _write_core(tmp_path, ("S1",))
    for read, trimmed in (("read1", "1_val_1"), ("read2", "2_val_2")):
        number = "1" if read == "read1" else "2"
        _write(tmp_path, f"fastqc/raw/S1_raw_{number}_fastqc.html")
        _write(tmp_path, f"fastqc/raw/S1_raw_{number}_fastqc.zip")
        _write(tmp_path, f"fastqc/trim/S1_trimmed_{trimmed}_fastqc.html")
        _write(tmp_path, f"fastqc/trim/S1_trimmed_{trimmed}_fastqc.zip")
        _write(
            tmp_path,
            f"trimgalore/S1_trimmed_{number}.fastq.gz_trimming_report.txt",
            b"Command line parameters: --paired /private/input.fastq.gz\n",
        )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report.html",
        b"<html>/private/nextflow/work/task</html>",
    )
    for name in (
        "multiqc_general_stats.txt",
        "multiqc_cutadapt.txt",
        "multiqc_fastqc_fastqc_raw.txt",
        "multiqc_fastqc_fastqc_trimmed.txt",
        "multiqc_star.txt",
        "multiqc_salmon.txt",
    ):
        _write(tmp_path, f"multiqc/star_salmon/multiqc_report_data/{name}")

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    types = _types(result)

    for stage in ("raw", "trimmed"):
        for read in ("read1", "read2"):
            assert f"bulk_rnaseq.fastqc.{stage}.{read}.zip" in types
    assert "bulk_rnaseq.trim.trimgalore.read1.report" not in types
    assert "bulk_rnaseq.trim.trimgalore.read2.report" not in types
    assert "bulk_rnaseq.multiqc.report" not in types
    assert "bulk_rnaseq.multiqc.star" in types
    assert "bulk_rnaseq.multiqc.salmon" in types


def test_disabled_multiqc_modules_are_not_cataloged_when_stale_files_exist(
    tmp_path: Path,
):
    inputs = _inputs(qc=_multiqc_only())
    _write_core(tmp_path, ("S1",))
    required = {"multiqc_general_stats.txt", "multiqc_star.txt", "multiqc_salmon.txt"}
    disabled = {
        "multiqc_cutadapt.txt",
        "multiqc_fastp.txt",
        "multiqc_fastqc_fastqc_raw.txt",
        "multiqc_fastqc_fastqc_trimmed.txt",
        "multiqc_fastqc_fastqc_filtered.txt",
        "multiqc_picard_dups.txt",
        "multiqc_featurecounts_biotype_plot.txt",
        "multiqc_rseqc_bam_stat.txt",
        "multiqc_rseqc_infer_experiment.txt",
        "multiqc_rseqc_read_distribution.txt",
        "multiqc_rseqc_tin.txt",
    }
    for name in sorted(required | disabled):
        _write(tmp_path, f"multiqc/star_salmon/multiqc_report_data/{name}")
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success
    multiqc_types = {
        candidate.output_type
        for candidate in result.value
        if candidate.output_type.startswith("bulk_rnaseq.multiqc.")
    }
    assert multiqc_types == {
        "bulk_rnaseq.multiqc.general_stats",
        "bulk_rnaseq.multiqc.salmon",
        "bulk_rnaseq.multiqc.star",
    }


def test_fastp_and_saved_reads_follow_normalized_toggles(tmp_path: Path):
    inputs = _inputs(
        layouts=("SE",),
        trimming={"enabled": True, "tool": "fastp"},
        outputs={"trimmed_reads": True, "merged_fastq": True},
    )
    _write_core(tmp_path, ("S1",))
    for relative in (
        "fastp/S1.fastp.json",
        "fastp/S1.fastp.html",
        "fastp/log/S1.fastp.log",
        "fastp/S1.fastp.fastq.gz",
        "fastq/S1.merged.fastq.gz",
    ):
        _write(tmp_path, relative)

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    types = _types(result)
    assert "bulk_rnaseq.trim.fastp.json" in types
    assert "bulk_rnaseq.trim.fastp.html" not in types
    assert "bulk_rnaseq.trim.fastp.filtered.single" in types
    assert "bulk_rnaseq.fastq.merged.single" in types
    assert "bulk_rnaseq.trim.fastp.log" not in types


def test_fastp_paired_fastqc_does_not_conflate_cat_and_stitched_outputs(
    tmp_path: Path,
):
    qc = {
        "enabled": True,
        "fastqc": True,
        "multiqc": False,
        "rseqc": False,
        "qualimap": False,
        "dupradar": False,
        "biotype": False,
        "deseq2_pca": False,
        "preseq": False,
        "mark_duplicates": False,
    }
    inputs = _inputs(
        layouts=("PE",),
        trimming={"enabled": True, "tool": "fastp"},
        qc=qc,
        outputs={"merged_fastq": True},
    )
    _write_core(tmp_path, ("S1",))
    for number in (1, 2):
        _write(tmp_path, f"fastqc/raw/S1_raw_{number}_fastqc.html")
        _write(tmp_path, f"fastqc/raw/S1_raw_{number}_fastqc.zip")
        _write(tmp_path, f"fastqc/trim/S1_trimmed_{number}_fastqc.html")
        _write(tmp_path, f"fastqc/trim/S1_trimmed_{number}_fastqc.zip")
        _write(tmp_path, f"fastq/S1_{number}.merged.fastq.gz")
    for relative in (
        "fastp/S1.fastp.json",
        "fastp/S1.fastp.html",
        "fastp/log/S1.fastp.log",
        "fastp/S1.merged.fastq.gz",
    ):
        _write(tmp_path, relative)

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    types = _types(result)
    assert "bulk_rnaseq.fastqc.trimmed.read1.zip" in types
    assert "bulk_rnaseq.fastqc.trimmed.read2.zip" in types
    # The pinned top-level workflow fixes fastp_merge=false. Standard
    # save_merged_fastq controls CAT_FASTQ only; it must not invent this file.
    assert "bulk_rnaseq.trim.fastp.merged" not in types
    assert "bulk_rnaseq.fastq.merged.read1" in types
    assert "bulk_rnaseq.fastq.merged.read2" in types


def test_star_salmon_featurecounts_rseqc_and_picard_sources_are_bound(tmp_path: Path):
    qc = {
        "enabled": True,
        "fastqc": False,
        "multiqc": True,
        "rseqc": True,
        "qualimap": False,
        "dupradar": False,
        "biotype": True,
        "deseq2_pca": False,
        "preseq": False,
        "mark_duplicates": True,
    }
    inputs = _inputs(
        qc=qc,
        advanced={"rseqc_modules": "bam_stat,infer_experiment,read_distribution,tin"},
    )
    _write_core(tmp_path, ("S1",))
    # MarkDuplicates changes the final published BAM identity.
    (tmp_path / "results/star_salmon/S1.sorted.bam").rename(
        tmp_path / "results/star_salmon/S1.markdup.sorted.bam"
    )
    (tmp_path / "results/star_salmon/S1.sorted.bam.bai").rename(
        tmp_path / "results/star_salmon/S1.markdup.sorted.bam.bai"
    )
    _write(tmp_path, "star_salmon/picard_metrics/S1.markdup.sorted.metrics.txt")
    _write(tmp_path, "star_salmon/featurecounts/S1.featureCounts.tsv")
    _write(tmp_path, "star_salmon/featurecounts/S1.featureCounts.tsv.summary")
    _write(tmp_path, "star_salmon/featurecounts/S1.biotype_counts_mqc.tsv")
    _write(tmp_path, "star_salmon/featurecounts/S1.biotype_counts_rrna_mqc.tsv")
    _write(tmp_path, "star_salmon/rseqc/bam_stat/S1.bam_stat.txt")
    _write(tmp_path, "star_salmon/rseqc/infer_experiment/S1.infer_experiment.txt")
    _write(tmp_path, "star_salmon/rseqc/read_distribution/S1.read_distribution.txt")
    _write(tmp_path, "star_salmon/rseqc/tin/S1.summary.txt")
    _write(tmp_path, "star_salmon/rseqc/tin/S1.tin.xls")
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    for name in (
        "multiqc_general_stats.txt",
        "multiqc_star.txt",
        "multiqc_salmon.txt",
        "multiqc_picard_dups.txt",
        "multiqc_featurecounts_biotype_plot.txt",
        "multiqc_rseqc_bam_stat.txt",
        "multiqc_rseqc_infer_experiment.txt",
        "multiqc_rseqc_read_distribution.txt",
    ):
        _write(tmp_path, f"multiqc/star_salmon/multiqc_report_data/{name}")

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    types = _types(result)
    assert "bulk_rnaseq.star.log_final" in types
    assert "bulk_rnaseq.salmon.meta_info" in types
    assert "bulk_rnaseq.featurecounts.summary" in types
    assert "bulk_rnaseq.featurecounts.counts" not in types
    assert "bulk_rnaseq.multiqc.picard_dups" in types
    picard_source = next(
        item
        for item in result.value
        if item.output_type == "bulk_rnaseq.multiqc.picard_dups"
    )
    assert picard_source.metadata == {"scope": "run", "assay": "bulk-rnaseq"}
    assert "bulk_rnaseq.picard.duplication_metrics" not in types
    for module in ("bam_stat", "infer_experiment", "read_distribution", "tin_summary"):
        assert f"bulk_rnaseq.rseqc.{module}" in types


def test_rseqc_inner_distance_uses_fixed_distance_and_mean_output_semantics(
    tmp_path: Path,
):
    qc = {
        "enabled": True,
        "fastqc": False,
        "multiqc": False,
        "rseqc": True,
        "qualimap": False,
        "dupradar": False,
        "biotype": False,
        "deseq2_pca": False,
        "preseq": False,
        "mark_duplicates": False,
    }
    inputs = _inputs(
        layouts=("PE",),
        qc=qc,
        advanced={"rseqc_modules": "inner_distance"},
    )
    _write_core(tmp_path, ("S1",))
    # Fixed nf-core/rnaseq 3.26.0 RSEQC_INNERDISTANCE emits these exact names:
    # main.nf labels them distance, freq, and mean respectively.
    for suffix in ("", "_freq", "_mean"):
        _write(
            tmp_path,
            f"star_salmon/rseqc/inner_distance/txt/S1.inner_distance{suffix}.txt",
        )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success
    by_type = {candidate.output_type: candidate for candidate in result.value}
    assert "bulk_rnaseq.rseqc.inner_distance.summary" not in by_type
    assert by_type["bulk_rnaseq.rseqc.inner_distance.distance"].relative_path == (
        "results/star_salmon/rseqc/inner_distance/txt/S1.inner_distance.txt"
    )
    assert by_type["bulk_rnaseq.rseqc.inner_distance.mean"].relative_path == (
        "results/star_salmon/rseqc/inner_distance/txt/S1.inner_distance_mean.txt"
    )


def test_csi_profile_omits_fixed_incompatible_rseqc_modules(tmp_path: Path):
    qc = {
        "enabled": True,
        "fastqc": False,
        "multiqc": False,
        "rseqc": True,
        "qualimap": False,
        "dupradar": False,
        "biotype": False,
        "deseq2_pca": False,
        "preseq": False,
        "mark_duplicates": False,
    }
    inputs = _inputs(
        qc=qc,
        advanced={
            "bam_csi_index": True,
            "rseqc_modules": "infer_experiment,read_distribution,tin",
        },
    )
    _write_core(tmp_path, ("S1",), csi=True)
    _write(tmp_path, "star_salmon/rseqc/infer_experiment/S1.infer_experiment.txt")

    types = _types(discover_bulk_rnaseq_artifacts(inputs, tmp_path))

    assert "bulk_rnaseq.rseqc.infer_experiment" in types
    assert "bulk_rnaseq.rseqc.read_distribution" not in types
    assert "bulk_rnaseq.rseqc.tin_summary" not in types


def test_command_bearing_samtools_stats_and_opaque_rds_are_not_public(
    tmp_path: Path,
):
    inputs = _inputs()
    _write_core(tmp_path, ("S1",))
    _write(
        tmp_path,
        "star_salmon/samtools_stats/S1.sorted.bam.stats",
        b"# This file was produced by samtools stats\n# The command line was: ...\n",
    )
    _write(tmp_path, "star_salmon/samtools_stats/S1.sorted.bam.flagstat")
    _write(tmp_path, "star_salmon/samtools_stats/S1.sorted.bam.idxstats")
    _write(tmp_path, "star_salmon/salmon.merged.gene.SummarizedExperiment.rds")

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    types = _types(result)

    assert "bulk_rnaseq.samtools.stats" not in types
    assert "bulk_rnaseq.samtools.flagstat" in types
    assert "bulk_rnaseq.samtools.idxstats" in types
    assert not any("summarized_experiment" in output_type for output_type in types)


@pytest.mark.parametrize("tool", ["sortmerna", "bowtie2"])
def test_rrna_saved_reads_are_cataloged_but_raw_logs_are_private(
    tmp_path: Path,
    tool: str,
):
    rrna = {
        "enabled": True,
        "tool": tool,
        "database_manifest": {"path": "/inputs/rrna.txt", "identity_sha256": "3" * 64},
        "save_filtered_reads": True,
    }
    inputs = _inputs(rrna=rrna)
    _write_core(tmp_path, ("S1",))
    if tool == "sortmerna":
        _write(tmp_path, "sortmerna/S1.sortmerna.log")
        _write(tmp_path, "sortmerna/S1.non_rRNA.fastq.gz")
    else:
        _write(tmp_path, "bowtie2_rrna/S1.bowtie2_rrna.bowtie2.log")
        _write(tmp_path, "bowtie2_rrna/S1.bowtie2_rrna.unmapped.fastq.gz")

    types = _types(discover_bulk_rnaseq_artifacts(inputs, tmp_path))
    assert f"bulk_rnaseq.rrna.{tool}.log" not in types
    assert f"bulk_rnaseq.rrna.{tool}.filtered.single" in types


def test_read_name_umi_route_does_not_require_extract_artifacts(tmp_path: Path):
    inputs = _inputs(
        umi={
            "enabled": True,
            "mode": "read_name",
            "deduplication_tool": "umitools",
            "read_name_separator": ":",
        }
    )
    _write_core(tmp_path, ("S1",))
    _write(tmp_path, "star_salmon/S1.transcriptome.sorted.bam")
    _write(tmp_path, "star_salmon/S1.transcriptome.sorted.bam.bai")
    _write(
        tmp_path,
        "star_salmon/umitools/genomic_dedup_log/S1.umi_dedup.sorted.log",
    )
    _write(
        tmp_path,
        (
            "star_salmon/umitools/transcriptomic_dedup_log/"
            "S1.umi_dedup.transcriptome.sorted.log"
        ),
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    types = _types(result)
    assert "bulk_rnaseq.umi.genomic_dedup_log" not in types
    assert "bulk_rnaseq.umi.transcriptomic_dedup_log" not in types
    assert "bulk_rnaseq.umi.extract_log" not in types
    assert "bulk_rnaseq.star.bam" not in types
    assert "bulk_rnaseq.star.sorted_intermediate_bam" in types
    assert "bulk_rnaseq.star.transcriptome_sorted_intermediate_bam" in types


def test_umi_alignment_intermediates_publish_transcriptome_dedup_bams(
    tmp_path: Path,
):
    inputs = _inputs(
        outputs={"alignment_intermediates": True},
        umi={
            "enabled": True,
            "mode": "read_name",
            "deduplication_tool": "umitools",
            "read_name_separator": ":",
        },
    )
    _write_core(tmp_path, ("S1",))
    for relative in (
        "star_salmon/S1.transcriptome.sorted.bam",
        "star_salmon/S1.transcriptome.sorted.bam.bai",
        "star_salmon/S1.umi_dedup.sorted.bam",
        "star_salmon/S1.umi_dedup.sorted.bam.bai",
        "star_salmon/S1.Aligned.out.bam",
        "star_salmon/S1.Aligned.toTranscriptome.out.bam",
        "star_salmon/S1.umi_dedup.transcriptome.bam",
        "star_salmon/S1.umi_dedup.transcriptome.sorted.bam",
        "star_salmon/S1.umi_dedup.transcriptome.sorted.bam.bai",
    ):
        _write(tmp_path, relative)

    types = _types(discover_bulk_rnaseq_artifacts(inputs, tmp_path))

    assert "bulk_rnaseq.umi.transcriptome_name_sorted_bam" in types
    assert "bulk_rnaseq.umi.transcriptome_sorted_bam" in types
    assert "bulk_rnaseq.umi.transcriptome_sorted_bam_index_bai" in types
    assert not any(
        output_type.startswith("bulk_rnaseq.umi.extracted") for output_type in types
    )


def test_umi_and_enabled_optional_outputs_have_exact_pinned_names(tmp_path: Path):
    inputs = _inputs(
        layouts=("PE",),
        trimming={"enabled": True, "tool": "trimgalore"},
        outputs={
            "alignment_intermediates": True,
            "bigwig": True,
            "trimmed_reads": True,
            "umi_intermediates": True,
            "unaligned_reads": True,
        },
        umi={
            "enabled": True,
            "mode": "read_sequence",
            "deduplication_tool": "umitools",
            "extraction_method": "string",
            "barcode_pattern": "NNNN",
            "emit_dedup_stats": True,
        },
        advanced={"bam_csi_index": True},
    )
    inputs.samples[0]["strandedness"] = "forward"
    _write_core(tmp_path, ("S1",), csi=True)
    for relative in (
        "trimgalore/S1_trimmed_1_val_1.fq.gz",
        "trimgalore/S1_trimmed_2_val_2.fq.gz",
        "umitools/S1.umi_extract_1.fastq.gz",
        "umitools/S1.umi_extract_2.fastq.gz",
        "star_salmon/S1.Aligned.out.bam",
        "star_salmon/S1.Aligned.toTranscriptome.out.bam",
        "star_salmon/S1.transcriptome.sorted.bam",
        "star_salmon/S1.transcriptome.sorted.bam.csi",
        "star_salmon/S1.umi_dedup.sorted.bam",
        "star_salmon/S1.umi_dedup.sorted.bam.csi",
        "star_salmon/S1.umi_dedup.transcriptome.bam",
        "star_salmon/S1.umi_dedup.transcriptome.filtered.bam",
        "star_salmon/S1.umi_dedup.transcriptome.sorted.bam",
        "star_salmon/S1.umi_dedup.transcriptome.sorted.bam.bai",
        "star_salmon/unmapped/S1.unmapped_1.fastq.gz",
        "star_salmon/unmapped/S1.unmapped_2.fastq.gz",
        "star_salmon/bigwig/S1.bigWig",
        "star_salmon/bigwig/S1.forward.bigWig",
        "star_salmon/bigwig/S1.reverse.bigWig",
    ):
        _write(tmp_path, relative)
    for target, stem in (
        ("genomic", "S1.umi_dedup.sorted"),
        ("transcriptomic", "S1.umi_dedup.transcriptome.sorted"),
    ):
        for statistic in ("edit_distance", "per_umi", "per_umi_per_position"):
            _write(
                tmp_path,
                f"star_salmon/umitools/{stem}_{statistic}.tsv",
            )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    types = _types(result)

    assert "bulk_rnaseq.star.bam_index.csi" in types
    assert "bulk_rnaseq.star.aligned_genome_bam" in types
    assert "bulk_rnaseq.star.aligned_transcriptome_bam" in types
    assert "bulk_rnaseq.star.unaligned.read1" in types
    assert "bulk_rnaseq.star.unaligned.read2" in types
    assert "bulk_rnaseq.trim.trimgalore.filtered.read1" in types
    assert "bulk_rnaseq.trim.trimgalore.filtered.read2" in types
    assert "bulk_rnaseq.umi.extracted.read1" in types
    assert "bulk_rnaseq.umi.extracted.read2" in types
    assert "bulk_rnaseq.umi.transcriptome_sorted_bam_index_bai" in types
    assert "bulk_rnaseq.star.sorted_intermediate_bam_index.csi" in types
    assert "bulk_rnaseq.star.transcriptome_sorted_intermediate_bam_index.csi" in types
    assert "bulk_rnaseq.umi.genomic_edit_distance" in types
    assert "bulk_rnaseq.umi.transcriptomic_per_umi_per_position" in types
    assert "bulk_rnaseq.star.bigwig.forward" in types
    assert "bulk_rnaseq.star.bigwig.reverse" in types
    assert "bulk_rnaseq.star.bigwig.combined" in types


def test_pe_umi_discard_read_uses_single_end_downstream_artifact_names(
    tmp_path: Path,
):
    qc = {
        "enabled": True,
        "fastqc": True,
        "multiqc": False,
        "rseqc": True,
        "qualimap": False,
        "dupradar": False,
        "biotype": False,
        "deseq2_pca": False,
        "preseq": False,
        "mark_duplicates": False,
    }
    inputs = _inputs(
        layouts=("PE",),
        trimming={"enabled": True, "tool": "trimgalore"},
        qc=qc,
        outputs={"trimmed_reads": True, "unaligned_reads": True},
        rrna={
            "enabled": True,
            "tool": "sortmerna",
            "database_manifest": {
                "path": "/inputs/rrna.txt",
                "identity_sha256": "3" * 64,
            },
            "save_filtered_reads": True,
        },
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
    _write_core(tmp_path, ("S1",))
    for relative in (
        "star_salmon/S1.transcriptome.sorted.bam",
        "star_salmon/S1.transcriptome.sorted.bam.bai",
        "fastqc/raw/S1_raw_1_fastqc.html",
        "fastqc/raw/S1_raw_1_fastqc.zip",
        "fastqc/raw/S1_raw_2_fastqc.html",
        "fastqc/raw/S1_raw_2_fastqc.zip",
        "fastqc/trim/S1_trimmed_trimmed_fastqc.html",
        "fastqc/trim/S1_trimmed_trimmed_fastqc.zip",
        "fastqc/filtered/S1_filtered_fastqc.html",
        "fastqc/filtered/S1_filtered_fastqc.zip",
        "trimgalore/S1_trimmed_trimmed.fq.gz",
        "sortmerna/S1.non_rRNA.fastq.gz",
        "star_salmon/unmapped/S1.unmapped_1.fastq.gz",
        "star_salmon/rseqc/infer_experiment/S1.infer_experiment.txt",
    ):
        _write(tmp_path, relative)

    types = _types(discover_bulk_rnaseq_artifacts(inputs, tmp_path))

    assert "bulk_rnaseq.fastqc.raw.read1.zip" in types
    assert "bulk_rnaseq.fastqc.raw.read2.zip" in types
    assert "bulk_rnaseq.fastqc.trimmed.single.zip" in types
    assert "bulk_rnaseq.fastqc.filtered.single.zip" in types
    assert "bulk_rnaseq.trim.trimgalore.filtered.single" in types
    assert "bulk_rnaseq.rrna.sortmerna.filtered.single" in types
    assert "bulk_rnaseq.star.unaligned.single" in types
    assert not any(
        output_type.endswith((".read1", ".read2"))
        for output_type in types
        if "trimmed" in output_type or "filtered" in output_type
    )


def test_trimmed_failure_table_proves_a_legal_missing_analysis_sample(
    tmp_path: Path,
):
    inputs = _inputs(
        layouts=("SE", "SE"),
        trimming={"enabled": True, "tool": "trimgalore"},
        qc=_multiqc_only(),
    )
    _write_core(tmp_path, ("S2",))
    for sample in ("S1", "S2"):
        _write(
            tmp_path,
            f"trimgalore/{sample}_trimmed.fastq.gz_trimming_report.txt",
        )
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    for name in (
        "multiqc_general_stats.txt",
        "multiqc_cutadapt.txt",
        "multiqc_star.txt",
        "multiqc_salmon.txt",
    ):
        _write(tmp_path, f"multiqc/star_salmon/multiqc_report_data/{name}")
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t9999\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    types = _types(result)
    assert "bulk_rnaseq.multiqc.fail_trimmed_samples" in types
    assert not any(
        item.metadata.get("sample_id") == "S1"
        and item.output_type.startswith(("bulk_rnaseq.star.", "bulk_rnaseq.salmon."))
        for item in result.value
    )


def test_all_trimmed_failures_do_not_require_star_salmon_or_merged_matrices(
    tmp_path: Path,
):
    inputs = _inputs(
        trimming={"enabled": True, "tool": "trimgalore"},
        qc=_multiqc_only(),
    )
    _write(tmp_path, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    _write(tmp_path, "trimgalore/S1_trimmed.fastq.gz_trimming_report.txt")
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    for name in ("multiqc_general_stats.txt", "multiqc_cutadapt.txt"):
        _write(tmp_path, f"multiqc/star_salmon/multiqc_report_data/{name}")
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t9999\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert result.is_success
    assert not any(
        item.output_type.startswith(("bulk_rnaseq.star.", "bulk_rnaseq.salmon."))
        for item in result.value
    )


def test_mapped_failure_table_keeps_star_salmon_but_not_post_mapping_outputs(
    tmp_path: Path,
):
    inputs = _inputs(layouts=("SE", "SE"), qc=_multiqc_only())
    _write_core(tmp_path, ("S1", "S2"))
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    for name in ("multiqc_general_stats.txt", "multiqc_star.txt", "multiqc_salmon.txt"):
        _write(tmp_path, f"multiqc/star_salmon/multiqc_report_data/{name}")
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_mapped_samples_table.txt",
        b"Sample\tSTAR uniquely mapped reads (%)\nS1\t4.99\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert result.is_success
    s1_types = {
        item.output_type
        for item in result.value
        if item.metadata.get("sample_id") == "S1"
    }
    assert "bulk_rnaseq.star.log_final" in s1_types
    assert "bulk_rnaseq.salmon.meta_info" in s1_types
    assert "bulk_rnaseq.star.bam" in s1_types
    assert "bulk_rnaseq.star.bam_index.bai" in s1_types
    assert not any(
        output_type.startswith(
            (
                "bulk_rnaseq.star.bigwig.",
                "bulk_rnaseq.featurecounts.",
                "bulk_rnaseq.rseqc.",
            )
        )
        for output_type in s1_types
    )
    assert "bulk_rnaseq.multiqc.fail_mapped_samples" in _types(result)


def test_mapped_failure_keeps_published_umi_and_pre_umi_bams(tmp_path: Path):
    inputs = _inputs(
        qc=_multiqc_only(),
        outputs={"umi_intermediates": True},
        umi={
            "enabled": True,
            "mode": "read_name",
            "deduplication_tool": "umitools",
            "read_name_separator": ":",
        },
    )
    _write_core(tmp_path, ("S1",))
    for relative in (
        "star_salmon/S1.transcriptome.sorted.bam",
        "star_salmon/S1.transcriptome.sorted.bam.bai",
        "star_salmon/S1.umi_dedup.sorted.bam",
        "star_salmon/S1.umi_dedup.sorted.bam.bai",
        "star_salmon/S1.umi_dedup.transcriptome.bam",
        "star_salmon/S1.umi_dedup.transcriptome.sorted.bam",
        "star_salmon/S1.umi_dedup.transcriptome.sorted.bam.bai",
    ):
        _write(tmp_path, relative)
    for name in ("multiqc_general_stats.txt", "multiqc_star.txt", "multiqc_salmon.txt"):
        _write(tmp_path, f"multiqc/star_salmon/multiqc_report_data/{name}")
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_mapped_samples_table.txt",
        b"Sample\tSTAR uniquely mapped reads (%)\nS1\t4.99\n",
    )

    types = _types(discover_bulk_rnaseq_artifacts(inputs, tmp_path))

    assert "bulk_rnaseq.star.bam" in types
    assert "bulk_rnaseq.star.sorted_intermediate_bam" in types
    assert "bulk_rnaseq.star.transcriptome_sorted_intermediate_bam" in types
    assert "bulk_rnaseq.umi.transcriptome_sorted_bam" in types


def test_exact_trim_threshold_row_does_not_hide_existing_analysis_outputs(
    tmp_path: Path,
):
    inputs = _inputs(qc=_multiqc_only())
    _write_core(tmp_path, ("S1",))
    for name in ("multiqc_general_stats.txt", "multiqc_star.txt", "multiqc_salmon.txt"):
        _write(tmp_path, f"multiqc/star_salmon/multiqc_report_data/{name}")
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t10000.0\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success
    assert "bulk_rnaseq.star.bam" in _types(result)


def test_trim_status_above_configured_threshold_fails_closed(tmp_path: Path):
    inputs = _inputs(
        trimming={"enabled": True, "tool": "trimgalore"},
        qc=_multiqc_only(),
        advanced={"min_trimmed_reads": 100},
    )
    _write_core(tmp_path, ("S1",))
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t101.0\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_STATUS_INVALID"


def test_fastp_low_trim_sample_does_not_require_post_filter_fastqc(
    tmp_path: Path,
):
    qc = {**_multiqc_only(), "fastqc": True}
    inputs = _inputs(
        trimming={"enabled": True, "tool": "fastp"},
        qc=qc,
    )
    _write(tmp_path, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    _write(tmp_path, "fastp/S1.fastp.json")
    _write(tmp_path, "fastqc/raw/S1_raw_fastqc.html")
    _write(tmp_path, "fastqc/raw/S1_raw_fastqc.zip")
    for name in (
        "multiqc_general_stats.txt",
        "multiqc_fastp.txt",
        "multiqc_fastqc_fastqc_raw.txt",
    ):
        _write(tmp_path, f"multiqc/star_salmon/multiqc_report_data/{name}")
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t9999\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success
    types = _types(result)
    assert "bulk_rnaseq.fastqc.raw.single.zip" in types
    assert "bulk_rnaseq.fastqc.trimmed.single.zip" not in types
    assert "bulk_rnaseq.multiqc.fastqc.trimmed" not in types


@pytest.mark.parametrize(
    ("trimmed", "mapped"),
    [
        (b"Sample\tReads after trimming\nOTHER\t1\n", None),
        (b"Sample\tReads after trimming\nS1\t1\nS1\t2\n", None),
        (b"wrong\theader\nS1\t1\n", None),
        (
            b"Sample\tReads after trimming\nS1\t1\n",
            b"Sample\tSTAR uniquely mapped reads (%)\nS1\t1\n",
        ),
        (
            b"Sample\tReads after trimming\n",
            b"Sample\tSTAR uniquely mapped reads (%)\nS1\t5.5\n",
        ),
    ],
)
def test_sample_failure_tables_reject_unknown_duplicate_corrupt_or_overlap(
    tmp_path: Path,
    trimmed: bytes,
    mapped: bytes | None,
):
    inputs = _inputs(qc=_multiqc_only())
    if trimmed is not None:
        _write(
            tmp_path,
            "multiqc/star_salmon/multiqc_report_data/multiqc_fail_trimmed_samples_table.txt",
            trimmed,
        )
    if mapped is not None:
        _write(
            tmp_path,
            "multiqc/star_salmon/multiqc_report_data/multiqc_fail_mapped_samples_table.txt",
            mapped,
        )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_STATUS_INVALID"
    assert result.errors[0].context == {"reason": "sample_status_invalid"}


def test_status_table_symlink_size_limit_and_replacement_are_fail_closed(
    tmp_path: Path,
    monkeypatch,
):
    inputs = _inputs(qc=_multiqc_only())
    path = _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/status-target.txt",
        b"Sample\tReads after trimming\nS1\t1\n",
    )
    status = path.with_name("multiqc_fail_trimmed_samples_table.txt")
    status.symlink_to(path.name)
    unsafe = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert unsafe.is_failure
    assert unsafe.errors[0].code == "BULK_RNASEQ_ARTIFACT_PATH_UNSAFE"

    status.unlink()
    status.write_bytes(b"Sample\tReads after trimming\nS1\t1\n")
    monkeypatch.setattr(artifact_module, "_MAX_STATUS_TABLE_BYTES", 8)
    oversized = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert oversized.is_failure
    assert oversized.errors[0].code == "BULK_RNASEQ_ARTIFACT_LIMIT_EXCEEDED"

    monkeypatch.setattr(artifact_module, "_MAX_STATUS_TABLE_BYTES", 64 * 1024)
    original = artifact_module._open_relative_for_read
    replaced = False

    def replace_after_open(descriptor: int, relative_path: str):
        nonlocal replaced
        opened = original(descriptor, relative_path)
        if (
            relative_path.endswith("multiqc_fail_trimmed_samples_table.txt")
            and not replaced
        ):
            replaced = True
            status.write_bytes(b"Sample\tReads after trimming\nS1\t2\nextra")
        return opened

    monkeypatch.setattr(
        artifact_module,
        "_open_relative_for_read",
        replace_after_open,
    )
    raced = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert raced.is_failure
    assert raced.errors[0].code == "BULK_RNASEQ_ARTIFACT_RACE_DETECTED"


def test_missing_core_output_without_a_failure_table_still_fails(tmp_path: Path):
    inputs = _inputs(qc=_multiqc_only())
    _write_core(tmp_path, ("S1",))
    (tmp_path / "results/star_salmon/S1/quant.sf").unlink()
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    for name in ("multiqc_general_stats.txt", "multiqc_star.txt", "multiqc_salmon.txt"):
        _write(tmp_path, f"multiqc/star_salmon/multiqc_report_data/{name}")

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_REQUIRED_MISSING"


def test_multiqc_disabled_low_yield_ambiguity_fails_with_explicit_status_reason(
    tmp_path: Path,
):
    inputs = _inputs()
    _write(tmp_path, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")

    low_trim = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert low_trim.is_failure
    assert low_trim.errors[0].code == "BULK_RNASEQ_ARTIFACT_STATUS_UNAVAILABLE"
    assert low_trim.errors[0].context == {"reason": "sample_status_unavailable"}

    markdup_qc = {**_multiqc_only(), "multiqc": False, "mark_duplicates": True}
    mapped_inputs = _inputs(qc=markdup_qc)
    _write_core(tmp_path, ("S1",))
    low_map = discover_bulk_rnaseq_artifacts(mapped_inputs, tmp_path)

    assert low_map.is_failure
    assert low_map.errors[0].code == "BULK_RNASEQ_ARTIFACT_STATUS_UNAVAILABLE"


def test_auto_bigwig_requires_combined_and_rejects_a_half_directional_pair(
    tmp_path: Path,
):
    inputs = _inputs(outputs={"bigwig": True})
    inputs.samples[0]["strandedness"] = "auto"
    _write_core(tmp_path, ("S1",))
    _write(tmp_path, "star_salmon/bigwig/S1.forward.bigWig")
    failed = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert failed.is_failure
    assert failed.errors[0].code == "BULK_RNASEQ_ARTIFACT_REQUIRED_MISSING"

    _write(tmp_path, "star_salmon/bigwig/S1.reverse.bigWig")
    still_missing_combined = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert still_missing_combined.is_failure

    _write(tmp_path, "star_salmon/bigwig/S1.bigWig")
    success = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert success.is_success
    assert {
        candidate.output_type
        for candidate in success.value
        if "bigwig" in candidate.output_type
    } == {
        "bulk_rnaseq.star.bigwig.combined",
        "bulk_rnaseq.star.bigwig.forward",
        "bulk_rnaseq.star.bigwig.reverse",
    }


@pytest.mark.parametrize("kind", ["symlink", "directory", "fifo", "oversize"])
def test_unsafe_required_artifact_fails_with_a_stable_redacted_issue(
    tmp_path: Path,
    monkeypatch,
    kind: str,
):
    inputs = _inputs()
    _write_core(tmp_path, ("S1",))
    target = tmp_path / "results/star_salmon/S1/quant.sf"
    target.unlink()
    if kind == "symlink":
        target.symlink_to("quant.genes.sf")
    elif kind == "directory":
        target.mkdir()
    elif kind == "fifo":
        os.mkfifo(target)
    else:
        target.write_bytes(b"xx")
        monkeypatch.setattr(artifact_module, "_MAX_ARTIFACT_BYTES", 1)

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert result.is_failure
    assert result.errors[0].code in {
        "BULK_RNASEQ_ARTIFACT_PATH_UNSAFE",
        "BULK_RNASEQ_ARTIFACT_LIMIT_EXCEEDED",
    }
    serialized = repr(result.errors[0].to_dict())
    assert str(tmp_path) not in serialized
    assert "quant.sf" not in serialized


def test_results_root_symlink_and_observable_replacement_race_fail_closed(
    tmp_path: Path,
    monkeypatch,
):
    inputs = _inputs()
    real_workspace = tmp_path / "real"
    _write_core(real_workspace, ("S1",))
    linked_workspace = tmp_path / "linked"
    linked_workspace.mkdir()
    (linked_workspace / "results").symlink_to(real_workspace / "results")

    linked = discover_bulk_rnaseq_artifacts(inputs, linked_workspace)
    assert linked.is_failure
    assert linked.errors[0].code == "BULK_RNASEQ_ARTIFACT_PATH_UNSAFE"

    original = artifact_module._open_relative_identity
    calls = 0

    def replace_after_first_open(descriptor: int, relative_path: str):
        nonlocal calls
        identity = original(descriptor, relative_path)
        if relative_path == "star_salmon/S1/quant.sf":
            calls += 1
            if calls == 1:
                _write(real_workspace, relative_path, b"replacement-with-new-identity")
        return identity

    monkeypatch.setattr(
        artifact_module,
        "_open_relative_identity",
        replace_after_first_open,
    )
    raced = discover_bulk_rnaseq_artifacts(inputs, real_workspace)
    assert raced.is_failure
    assert raced.errors[0].code == "BULK_RNASEQ_ARTIFACT_RACE_DETECTED"
    assert raced.errors[0].context == {"reason": "artifact_identity_changed"}


@pytest.mark.parametrize(
    "workspace",
    [Path("relative-workspace"), Path("/tmp/work/../escaped")],
)
def test_relative_or_lexically_traversing_workspace_is_rejected(workspace: Path):
    result = discover_bulk_rnaseq_artifacts(_inputs(), workspace)
    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_CONTRACT_INVALID"
    assert result.errors[0].context == {"reason": "artifact_contract_invalid"}


def test_unsupported_output_namespaces_fail_explicitly(tmp_path: Path):
    for output in ("reference_files", "stringtie"):
        result = discover_bulk_rnaseq_artifacts(
            _inputs(outputs={output: True}), tmp_path
        )
        assert result.is_failure
        assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_PROFILE_UNSUPPORTED"
        assert result.errors[0].context == {"reason": "unaudited_output_profile"}
