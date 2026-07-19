"""Deterministic artifact discovery for pinned nf-core/rnaseq 3.26.0."""

from __future__ import annotations

from fnmatch import fnmatchcase
from io import BytesIO
from pathlib import Path
import json
import os
from zipfile import ZIP_DEFLATED, ZipFile

import pytest

import encode_pipeline.adapters.bulk_rnaseq.artifacts as artifact_module
import encode_pipeline.adapters.bulk_rnaseq.status_evidence as status_evidence_module
from encode_pipeline.adapters.bulk_rnaseq.artifacts import (
    discover_bulk_rnaseq_artifacts,
)
from encode_pipeline.adapters.bulk_rnaseq.results_contract import (
    load_bulk_rnaseq_results_contract,
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

_STAR_LOG = b"""\
Started job on | Jul 17 00:00:00
Started mapping on | Jul 17 00:00:01
Finished on | Jul 17 00:00:02
Mapping speed, Million of reads per hour | 1800.00
Number of input reads | 1000
Average input read length | 100
UNIQUE READS:
Uniquely mapped reads number | 800
Uniquely mapped reads % | 80.00%
Average mapped length | 99.00
Number of splices: Total | 100
Number of splices: Annotated (sjdb) | 90
Number of splices: GT/AG | 80
Number of splices: GC/AG | 10
Number of splices: AT/AC | 5
Number of splices: Non-canonical | 5
Mismatch rate per base, % | 0.10%
Deletion rate per base | 0.01%
Deletion average length | 1.00
Insertion rate per base | 0.01%
Insertion average length | 1.00
MULTI-MAPPING READS:
Number of reads mapped to multiple loci | 100
% of reads mapped to multiple loci | 10.00%
Number of reads mapped to too many loci | 20
% of reads mapped to too many loci | 2.00%
UNMAPPED READS:
Number of reads unmapped: too many mismatches | 10
% of reads unmapped: too many mismatches | 1.00%
Number of reads unmapped: too short | 60
% of reads unmapped: too short | 6.00%
Number of reads unmapped: other | 10
% of reads unmapped: other | 1.00%
CHIMERIC READS:
Number of chimeric reads | 0
% of chimeric reads | 0.00%
"""


def _zero_second_star_log() -> bytes:
    return _STAR_LOG.replace(
        b"Finished on | Jul 17 00:00:02",
        b"Finished on | Jul 17 00:00:01",
    ).replace(
        b"Mapping speed, Million of reads per hour | 1800.00",
        b"Mapping speed, Million of reads per hour | inf",
    )


def test_star_status_accepts_pinned_zero_second_infinite_mapping_speed():
    assert (
        str(
            status_evidence_module.parse_star_uniquely_mapped_percent(
                _zero_second_star_log()
            )
        )
        == "80.00"
    )


def test_star_status_rejects_infinite_mapping_speed_for_nonzero_duration():
    with pytest.raises(status_evidence_module.StatusEvidenceError):
        status_evidence_module.parse_star_uniquely_mapped_percent(
            _STAR_LOG.replace(
                b"Mapping speed, Million of reads per hour | 1800.00",
                b"Mapping speed, Million of reads per hour | inf",
            )
        )


@pytest.mark.parametrize(
    "invalid_speed",
    (b"Inf", b"+inf", b"Infinity", b"NaN", b"-inf"),
)
def test_star_status_rejects_other_nonfinite_zero_second_mapping_speeds(
    invalid_speed: bytes,
):
    with pytest.raises(status_evidence_module.StatusEvidenceError):
        status_evidence_module.parse_star_uniquely_mapped_percent(
            _zero_second_star_log().replace(
                b"Mapping speed, Million of reads per hour | inf",
                b"Mapping speed, Million of reads per hour | " + invalid_speed,
            )
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
            _write(
                workspace,
                relative,
                _STAR_LOG if relative.endswith(".Log.final.out") else b"x",
            )
    for suffix in _MERGED:
        _write(workspace, f"star_salmon/salmon.merged.{suffix}.tsv")


def _write_fastp_evidence(workspace: Path, sample: str, retained: int) -> None:
    _write(
        workspace,
        f"fastp/{sample}.fastp.json",
        json.dumps(
            {
                "fastp_version": "1.0.1",
                "summary": {
                    "before_filtering": {
                        "total_reads": retained + 100,
                        "total_bases": 1000000,
                    },
                    "after_filtering": {
                        "total_reads": retained,
                        "total_bases": 900000,
                    },
                },
                "filtering_result": {"passed_filter_reads": retained},
            }
        ).encode(),
    )


def _fastqc_evidence_zip(
    root: str,
    filename: str,
    retained: int,
    *,
    include_directory_entries: bool = False,
) -> bytes:
    data = (
        "\n".join(
            (
                "##FastQC\t0.12.1",
                ">>Basic Statistics\tpass",
                "#Measure\tValue",
                f"Filename\t{filename}",
                f"Total Sequences\t{retained}",
                "%GC\t42",
                ">>END_MODULE",
                "",
            )
        )
    ).encode()
    target = BytesIO()
    with ZipFile(target, "w", ZIP_DEFLATED) as archive:
        if include_directory_entries:
            archive.writestr(f"{root}/", b"")
            archive.writestr(f"{root}/Icons/", b"")
            archive.writestr(f"{root}/Images/", b"")
        archive.writestr(f"{root}/fastqc_data.txt", data)
    return target.getvalue()


def _multiqc_stub(*identities: str) -> bytes:
    return (
        "Sample\tmetric\n" + "".join(f"{identity}\t1\n" for identity in identities)
    ).encode()


def _write_trimgalore_fastqc_false_route(
    workspace: Path,
    *,
    downstream_layout: str,
    retained: int,
    high_yield: bool,
    umi_discard: bool = False,
) -> None:
    """Write the fixed 3.26.0 Trim Galore route without separate FastQC."""
    _write(workspace, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    if downstream_layout == "SE":
        roots = (("single", "S1_trimmed_trimmed", "S1"),)
        multiqc_identities = ("S1",)
        general_identities = ("S1",)
    else:
        roots = (
            ("read1", "S1_trimmed_1_val_1", "S1_1"),
            ("read2", "S1_trimmed_2_val_2", "S1_2"),
        )
        multiqc_identities = ("S1_1", "S1_2")
        general_identities = ("S1", "S1 Read 1", "S1 Read 2")
    for _role, basename, _identity in roots:
        root = f"{basename}_fastqc"
        _write(workspace, f"fastqc/trim/{root}.html")
        _write(
            workspace,
            f"fastqc/trim/{root}.zip",
            _fastqc_evidence_zip(root, f"{basename}.fq.gz", retained),
        )
    _write(workspace, "multiqc/star_salmon/multiqc_report.html")
    _write(
        workspace,
        "multiqc/star_salmon/multiqc_report_data/multiqc_general_stats.txt",
        _multiqc_stub(*general_identities),
    )
    _write(
        workspace,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fastqc_fastqc_trimmed.txt",
        _multiqc_stub(*multiqc_identities),
    )
    cutadapt_rows = "".join(
        f"{identity}\t4.0\t20000\t0\t20000\t2000000\t0\t999900\t50.005\n"
        for identity in multiqc_identities
    )
    _write(
        workspace,
        "multiqc/star_salmon/multiqc_report_data/multiqc_cutadapt.txt",
        b"Sample\tcutadapt_version\tr_processed\tr_with_adapters\tr_written\t"
        b"bp_processed\tquality_trimmed\tbp_written\tpercent_trimmed\n"
        + cutadapt_rows.encode(),
    )
    if not high_yield:
        _write(
            workspace,
            "multiqc/star_salmon/multiqc_report_data/"
            "multiqc_fail_trimmed_samples_table.txt",
            f"Sample\tReads after trimming\nS1\t{retained}\n".encode(),
        )
        return
    if umi_discard:
        for relative in (
            "star_salmon/S1.sorted.bam",
            "star_salmon/S1.sorted.bam.bai",
            "star_salmon/S1.transcriptome.sorted.bam",
            "star_salmon/S1.transcriptome.sorted.bam.bai",
            "star_salmon/log/S1.Log.final.out",
            "star_salmon/log/S1.SJ.out.tab",
            "star_salmon/S1/quant.sf",
            "star_salmon/S1/quant.genes.sf",
            "star_salmon/S1/aux_info/meta_info.json",
            "star_salmon/samtools_stats/S1.umi_dedup.sorted.bam.flagstat",
            "star_salmon/samtools_stats/S1.umi_dedup.sorted.bam.idxstats",
        ):
            _write(
                workspace,
                relative,
                _STAR_LOG if relative.endswith(".Log.final.out") else b"x",
            )
        for suffix in _MERGED:
            _write(workspace, f"star_salmon/salmon.merged.{suffix}.tsv")
    else:
        _write_core(workspace, ("S1",))
        _write(workspace, "star_salmon/samtools_stats/S1.sorted.bam.flagstat")
        _write(workspace, "star_salmon/samtools_stats/S1.sorted.bam.idxstats")
    _write(
        workspace,
        "multiqc/star_salmon/multiqc_report_data/multiqc_star.txt",
        _multiqc_stub("S1"),
    )
    _write(
        workspace,
        "multiqc/star_salmon/multiqc_report_data/multiqc_salmon.txt",
        _multiqc_stub("S1"),
    )


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


def _valid_fastp_payload() -> dict[str, object]:
    return {
        "fastp_version": "1.0.1",
        "summary": {
            "before_filtering": {"total_reads": 100, "total_bases": 10000},
            "after_filtering": {"total_reads": 90, "total_bases": 9000},
        },
        "filtering_result": {"passed_filter_reads": 90},
    }


@pytest.mark.parametrize(
    "content",
    (
        b'{"fastp_version":"wrong","fastp_version":"1.0.1",'
        b'"summary":{"before_filtering":{"total_reads":100,"total_bases":10000},'
        b'"after_filtering":{"total_reads":90,"total_bases":9000}},'
        b'"filtering_result":{"passed_filter_reads":90}}',
        json.dumps({**_valid_fastp_payload(), "junk": "x" * 8193}).encode(),
        json.dumps(
            {**_valid_fastp_payload(), "junk": [[[[[[[[[[[[[[[[[0]]]]]]]]]]]]]]]]]},
        ).encode(),
        json.dumps(
            {**_valid_fastp_payload(), "junk": [0] * 50_001},
            separators=(",", ":"),
        ).encode(),
    ),
    ids=("duplicate-key", "overlong-string", "over-depth", "over-node-limit"),
)
def test_fastp_shared_parser_preserves_strict_json_limits(content: bytes):
    with pytest.raises(status_evidence_module.StatusEvidenceError):
        status_evidence_module.parse_fastp_summary(content)


def test_cutadapt_shared_parser_bounds_numeric_fields_before_conversion():
    header = (
        b"Sample\tcutadapt_version\tr_processed\tr_with_adapters\tr_written\t"
        b"bp_processed\tquality_trimmed\tbp_written\tpercent_trimmed\n"
    )
    content = header + b"S1\t4.0\t" + b"9" * 5000 + b"\t0\t1\t1\t0\t1\t0\n"

    with pytest.raises(status_evidence_module.StatusEvidenceError):
        status_evidence_module.parse_cutadapt_processed_reads(
            content,
            expected_row_owners={"S1": "S1"},
        )


def test_status_fastqc_zip_rejects_nonempty_directory_entry():
    root = "S1_trimmed_trimmed_fastqc"
    target = BytesIO(
        _fastqc_evidence_zip(
            root,
            "S1_trimmed_trimmed.fq.gz",
            20000,
        )
    )
    with ZipFile(target, "a", ZIP_DEFLATED) as archive:
        archive.writestr(f"{root}/", b"not-a-directory-record")

    with pytest.raises(status_evidence_module.StatusEvidenceError):
        status_evidence_module.parse_fastqc_total_sequences(
            target.getvalue(),
            expected_root=root,
            expected_filename="S1_trimmed_trimmed.fq.gz",
        )


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


@pytest.mark.parametrize(
    "relative_path",
    (
        "star_salmon/bigwig/UNKNOWN.bigWig",
        "star_salmon/featurecounts/UNKNOWN.featureCounts.tsv.summary",
        "star_salmon/samtools_stats/UNKNOWN.sorted.bam.flagstat",
        "star_salmon/umitools/UNKNOWN.umi_dedup.sorted_edit_distance.tsv",
        "star_salmon/unmapped/UNKNOWN.unmapped_1.fastq.gz",
        "star_salmon/rseqc/tin/UNKNOWN.summary.txt",
        "fastqc/raw/UNKNOWN_raw_fastqc.zip",
        "fastp/UNKNOWN.fastp.json",
        "trimgalore/UNKNOWN_trimmed.fastq.gz_trimming_report.txt",
        "fastq/UNKNOWN.merged.fastq.gz",
        "umitools/UNKNOWN.umi_extract.fastq.gz",
        "umitools/UNKNOWN.umi_extract_1.fastq.gz",
        "umitools/UNKNOWN.umi_extract_2.fastq.gz",
        "sortmerna/UNKNOWN.non_rRNA.fastq.gz",
        "bowtie2_rrna/UNKNOWN.bowtie2_rrna.unmapped.fastq.gz",
    ),
)
def test_unknown_sample_in_every_fixed_nested_namespace_fails_closed(
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


def test_recognized_nested_sample_entry_must_be_a_regular_file(tmp_path: Path):
    inputs = _inputs()
    _write_core(tmp_path, ("S1",))
    fifo = tmp_path / "results/star_salmon/bigwig/UNKNOWN.bigWig"
    fifo.parent.mkdir(parents=True)
    os.mkfifo(fifo)

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_PATH_UNSAFE"


@pytest.mark.parametrize("kind", ("directory", "fifo", "symlink"))
def test_private_recognized_nested_entry_must_remain_a_regular_file(
    tmp_path: Path,
    kind: str,
):
    inputs = _inputs()
    _write_core(tmp_path, ("S1",))
    path = tmp_path / "results/star_salmon/featurecounts/S1.featureCounts.tsv"
    path.parent.mkdir(parents=True, exist_ok=True)
    if kind == "directory":
        path.mkdir()
    elif kind == "fifo":
        os.mkfifo(path)
    else:
        target = path.with_name("private-target.txt")
        target.write_bytes(b"private")
        path.symlink_to(target.name)

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_PATH_UNSAFE"


def test_foreign_sample_added_after_namespace_audit_is_a_stable_race_failure(
    tmp_path: Path,
    monkeypatch,
):
    inputs = _inputs()
    _write_core(tmp_path, ("S1",))
    original = artifact_module._verified_regular_file
    injected = False

    def inject_foreign_sample(results_descriptor: int, relative_path: str) -> bool:
        nonlocal injected
        if not injected:
            injected = True
            _write(tmp_path, "star_salmon/bigwig/UNKNOWN.bigWig")
        return original(results_descriptor, relative_path)

    monkeypatch.setattr(
        artifact_module,
        "_verified_regular_file",
        inject_foreign_sample,
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_RACE_DETECTED"
    assert result.errors[0].context == {"reason": "artifact_identity_changed"}


def test_private_audited_entry_replacement_after_scan_is_a_race_failure(
    tmp_path: Path,
    monkeypatch,
):
    inputs = _inputs()
    _write_core(tmp_path, ("S1",))
    private = _write(
        tmp_path,
        "star_salmon/featurecounts/S1.featureCounts.tsv",
        b"first",
    )
    original = artifact_module._verified_regular_file
    replaced = False

    def replace_private_entry(results_descriptor: int, relative_path: str) -> bool:
        nonlocal replaced
        if not replaced:
            replaced = True
            private.unlink()
            private.write_bytes(b"replacement")
        return original(results_descriptor, relative_path)

    monkeypatch.setattr(
        artifact_module,
        "_verified_regular_file",
        replace_private_entry,
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_RACE_DETECTED"


def test_artifact_specs_reject_cross_type_relative_path_collision():
    metadata = {"scope": "sample", "sample_id": "A", "assay": "bulk-rnaseq"}
    specs = [
        artifact_module._ArtifactSpec(
            output_type="bulk_rnaseq.star.bigwig.forward",
            relative_path="star_salmon/bigwig/A.forward.bigWig",
            mime_type="application/octet-stream",
            metadata=metadata,
        ),
        artifact_module._ArtifactSpec(
            output_type="bulk_rnaseq.star.bigwig.combined",
            relative_path="star_salmon/bigwig/A.forward.bigWig",
            mime_type="application/octet-stream",
            metadata={**metadata, "sample_id": "A.forward"},
        ),
    ]

    with pytest.raises(ValueError, match="invalid artifact specification"):
        artifact_module._validate_specs(specs)


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


@pytest.mark.parametrize(
    ("layout", "umi", "expected_roles"),
    (
        ("SE", None, {"single"}),
        ("PE", None, {"read1", "read2"}),
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
            {"single"},
        ),
    ),
)
def test_trimgalore_bundled_fastqc_is_independent_of_raw_fastqc_toggle(
    layout: str,
    umi: dict[str, object] | None,
    expected_roles: set[str],
):
    inputs = _inputs(
        layouts=(layout,),
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={"enabled": True, "fastqc": False, "multiqc": True},
        umi=umi,
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

    assert {
        f"bulk_rnaseq.fastqc.trimmed.{role}.html" for role in expected_roles
    }.issubset(types)
    assert {
        f"bulk_rnaseq.fastqc.trimmed.{role}.zip" for role in expected_roles
    }.issubset(types)
    assert "bulk_rnaseq.multiqc.fastqc.trimmed" in types
    assert not any(".fastqc.raw." in value for value in types)


@pytest.mark.parametrize(
    ("layout", "umi", "downstream_layout", "expected_roles"),
    (
        ("SE", None, "SE", {"single"}),
        ("PE", None, "PE", {"read1", "read2"}),
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
            "SE",
            {"single"},
        ),
    ),
)
def test_trimgalore_fastqc_false_high_yield_full_discovery(
    tmp_path: Path,
    layout: str,
    umi: dict[str, object] | None,
    downstream_layout: str,
    expected_roles: set[str],
):
    inputs = _inputs(
        layouts=(layout,),
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={**_multiqc_only(), "fastqc": False},
        umi=umi,
    )
    _write_trimgalore_fastqc_false_route(
        tmp_path,
        downstream_layout=downstream_layout,
        retained=20000,
        high_yield=True,
        umi_discard=umi is not None,
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success, [issue.to_dict() for issue in result.issues]
    types = _types(result)
    assert {
        f"bulk_rnaseq.fastqc.trimmed.{role}.zip" for role in expected_roles
    }.issubset(types)
    assert "bulk_rnaseq.star.log_final" in types
    assert "bulk_rnaseq.salmon.meta_info" in types
    assert not any(".fastqc.raw." in value for value in types)


def test_trimgalore_status_accepts_standard_fastqc_directory_entries(
    tmp_path: Path,
):
    inputs = _inputs(
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={**_multiqc_only(), "fastqc": False},
    )
    _write_trimgalore_fastqc_false_route(
        tmp_path,
        downstream_layout="SE",
        retained=20000,
        high_yield=True,
    )
    root = "S1_trimmed_trimmed_fastqc"
    _write(
        tmp_path,
        f"fastqc/trim/{root}.zip",
        _fastqc_evidence_zip(
            root,
            "S1_trimmed_trimmed.fq.gz",
            20000,
            include_directory_entries=True,
        ),
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success, [issue.to_dict() for issue in result.issues]
    assert "bulk_rnaseq.fastqc.trimmed.single.zip" in _types(result)


@pytest.mark.parametrize("layout", ("SE", "PE"))
def test_trimgalore_fastqc_false_low_yield_full_discovery(
    tmp_path: Path,
    layout: str,
):
    inputs = _inputs(
        layouts=(layout,),
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={**_multiqc_only(), "fastqc": False},
    )
    _write_trimgalore_fastqc_false_route(
        tmp_path,
        downstream_layout=layout,
        retained=9999,
        high_yield=False,
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success, [issue.to_dict() for issue in result.issues]
    types = _types(result)
    assert "bulk_rnaseq.multiqc.fail_trimmed_samples" in types
    assert not any(
        value.startswith(("bulk_rnaseq.star.", "bulk_rnaseq.salmon."))
        for value in types
    )
    assert not any(".fastqc.raw." in value for value in types)


@pytest.mark.parametrize("mutation", ("missing_mate", "corrupt_zip", "mate_mismatch"))
def test_trimgalore_fastqc_false_pe_invalid_bundled_evidence_fails_closed(
    tmp_path: Path,
    mutation: str,
):
    inputs = _inputs(
        layouts=("PE",),
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={**_multiqc_only(), "fastqc": False},
    )
    _write_trimgalore_fastqc_false_route(
        tmp_path,
        downstream_layout="PE",
        retained=9999,
        high_yield=False,
    )
    read2 = tmp_path / "results/fastqc/trim/S1_trimmed_2_val_2_fastqc.zip"
    if mutation == "missing_mate":
        read2.unlink()
    elif mutation == "corrupt_zip":
        read2.write_bytes(b"not-a-zip")
    else:
        read2.write_bytes(
            _fastqc_evidence_zip(
                "S1_trimmed_2_val_2_fastqc",
                "S1_trimmed_2_val_2.fq.gz",
                9998,
            )
        )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_STATUS_INVALID"
    assert result.errors[0].context == {"reason": "sample_status_invalid"}


def test_declared_artifact_families_are_reachable_and_partition_expected_specs():
    qc = {
        "enabled": True,
        "fastqc": True,
        "multiqc": True,
        "rseqc": True,
        "qualimap": False,
        "dupradar": False,
        "biotype": True,
        "deseq2_pca": False,
        "preseq": False,
        "mark_duplicates": False,
    }
    comprehensive = _inputs(
        layouts=("PE",),
        trimming={"enabled": True, "tool": "trimgalore"},
        qc=qc,
        outputs={
            "alignment_intermediates": True,
            "bigwig": True,
            "merged_fastq": True,
            "trimmed_reads": True,
            "umi_intermediates": True,
            "unaligned_reads": True,
        },
        rrna={
            "enabled": True,
            "tool": "sortmerna",
            "save_filtered_reads": True,
            "database_manifest": {
                "path": "/inputs/rrna.txt",
                "identity_sha256": "3" * 64,
            },
        },
        umi={
            "enabled": True,
            "mode": "read_sequence",
            "deduplication_tool": "umitools",
            "extraction_method": "string",
            "barcode_pattern": "NNNN",
            "barcode_pattern_2": "NNNN",
            "emit_dedup_stats": True,
        },
        advanced={
            "rseqc_modules": (
                "bam_stat,inner_distance,infer_experiment,junction_annotation,"
                "junction_saturation,read_distribution,read_duplication,tin"
            )
        },
    )
    comprehensive.samples[0]["strandedness"] = "forward"
    fastp = _inputs(
        trimming={"enabled": True, "tool": "fastp"},
        outputs={"trimmed_reads": True},
        qc={
            "enabled": True,
            "fastqc": False,
            "multiqc": True,
            "rseqc": False,
            "qualimap": False,
            "dupradar": False,
            "biotype": False,
            "deseq2_pca": False,
            "preseq": False,
            "mark_duplicates": True,
        },
    )

    observed: set[str] = set()
    observed_multiqc_tables: set[str] = set()
    for profile in (comprehensive, fastp):
        validated = artifact_module.validate_bulk_rnaseq_inputs(profile)
        assert validated.is_success
        normalized = validated.value
        specs, auto_bigwig = artifact_module._expected_specs(
            normalized["nfcore_params"],
            artifact_module._normalized_samples(normalized["samples"]),
            trimmed_failed=frozenset(),
            mapped_failed=frozenset(),
        )
        assert not auto_bigwig
        observed.update(spec.output_type for spec in specs)
        observed_multiqc_tables.update(
            Path(spec.relative_path).name
            for spec in specs
            if spec.output_type.startswith("bulk_rnaseq.multiqc.")
        )

    contract = load_bulk_rnaseq_results_contract()
    membership = contract["artifact_family_output_types"]
    assert set(membership) == set(contract["artifact_families"])
    assert observed_multiqc_tables == set(contract["audited_published_multiqc_tables"])
    for output_type in observed:
        owners = {
            family
            for family, patterns in membership.items()
            if any(fnmatchcase(output_type, pattern) for pattern in patterns)
        }
        assert len(owners) == 1, (output_type, owners)
    for family, patterns in membership.items():
        assert any(
            fnmatchcase(output_type, pattern)
            for output_type in observed
            for pattern in patterns
        ), family
        for pattern in patterns:
            assert any(fnmatchcase(output_type, pattern) for output_type in observed), (
                family,
                pattern,
            )


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
        root = f"S1_trimmed_{trimmed}_fastqc"
        _write(
            tmp_path,
            f"fastqc/trim/{root}.zip",
            _fastqc_evidence_zip(
                root,
                f"S1_trimmed_{trimmed}.fq.gz",
                20000,
            ),
        )
        _write(
            tmp_path,
            f"trimgalore/S1_trimmed_{number}.fastq.gz_trimming_report.txt",
            b"20100 sequences processed in total\n"
            b"Sequences removed because they became shorter than the length "
            b"cutoff of 20 bp: 100\n"
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
        identities = (
            ("S1_1", "S1_2")
            if name
            in {
                "multiqc_cutadapt.txt",
                "multiqc_fastqc_fastqc_raw.txt",
                "multiqc_fastqc_fastqc_trimmed.txt",
            }
            else (
                ("S1", "S1 Read 1", "S1 Read 2")
                if name == "multiqc_general_stats.txt"
                else ("S1",)
            )
        )
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub(*identities),
        )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_cutadapt.txt",
        b"Sample\tcutadapt_version\tr_processed\tr_with_adapters\tr_written\t"
        b"bp_processed\tquality_trimmed\tbp_written\tpercent_trimmed\n"
        b"S1_1\t4.0\t20100\t0\t20100\t2010000\t0\t2010000\t0\n"
        b"S1_2\t4.0\t20100\t0\t20100\t2010000\t0\t2010000\t0\n",
    )

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
    enabled = {"multiqc_general_stats.txt", "multiqc_star.txt", "multiqc_salmon.txt"}
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
    for name in sorted(enabled | disabled):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1") if name in enabled else b"stale",
        )
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


def test_star_salmon_accepts_absent_optional_multiqc_salmon_table(
    tmp_path: Path,
):
    inputs = _inputs(qc=_multiqc_only())
    _write_core(tmp_path, ("S1",))
    root = "multiqc/star_salmon/multiqc_report_data"
    _write(tmp_path, f"{root}/multiqc_general_stats.txt", _multiqc_stub("S1"))
    _write(tmp_path, f"{root}/multiqc_star.txt", _multiqc_stub("S1"))

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success, [issue.to_dict() for issue in result.issues]
    types = _types(result)
    assert "bulk_rnaseq.multiqc.salmon" not in types
    assert "bulk_rnaseq.salmon.meta_info" in types
    assert "bulk_rnaseq.salmon.quant_transcript" in types
    assert "bulk_rnaseq.salmon.quant_gene" in types


def test_star_salmon_accepts_zero_second_infinite_mapping_speed(
    tmp_path: Path,
):
    inputs = _inputs(qc=_multiqc_only())
    _write_core(tmp_path, ("S1",))
    _write(tmp_path, "star_salmon/log/S1.Log.final.out", _zero_second_star_log())
    root = "multiqc/star_salmon/multiqc_report_data"
    _write(tmp_path, f"{root}/multiqc_general_stats.txt", _multiqc_stub("S1"))
    _write(tmp_path, f"{root}/multiqc_star.txt", _multiqc_stub("S1"))

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success, [issue.to_dict() for issue in result.issues]
    assert "bulk_rnaseq.star.log_final" in _types(result)


@pytest.mark.parametrize(
    ("output_type", "trimmer"),
    (
        ("bulk_rnaseq.multiqc.general_stats", "trimgalore"),
        ("bulk_rnaseq.multiqc.star", "trimgalore"),
        ("bulk_rnaseq.multiqc.salmon", "trimgalore"),
        ("bulk_rnaseq.multiqc.cutadapt", "trimgalore"),
        ("bulk_rnaseq.multiqc.fastp", "fastp"),
        ("bulk_rnaseq.multiqc.fastqc.raw", "trimgalore"),
        ("bulk_rnaseq.multiqc.fastqc.trimmed", "trimgalore"),
        ("bulk_rnaseq.multiqc.fastqc.trimmed", "fastp"),
        ("bulk_rnaseq.multiqc.fastqc.filtered", "trimgalore"),
        ("bulk_rnaseq.multiqc.picard_dups", "trimgalore"),
        ("bulk_rnaseq.multiqc.featurecounts_biotype", "trimgalore"),
        ("bulk_rnaseq.multiqc.rseqc.bam_stat", "trimgalore"),
        ("bulk_rnaseq.multiqc.rseqc.infer_experiment", "trimgalore"),
        ("bulk_rnaseq.multiqc.rseqc.read_distribution", "trimgalore"),
        ("bulk_rnaseq.multiqc.rseqc.tin", "trimgalore"),
    ),
)
def test_every_published_multiqc_table_policy_rejects_a_foreign_row(
    output_type: str,
    trimmer: str,
):
    samples = ({"sample": "S1", "layout": "PE", "strandedness": "unstranded"},)
    params = {
        "skip_fastqc": False,
        "skip_trimming": False,
        "trimmer": trimmer,
        "remove_ribo_rna": False,
    }
    policy = artifact_module._multiqc_row_policy(
        output_type,
        params=params,
        samples=samples,
        sample_ids=frozenset({"S1"}),
        analysis_samples=frozenset({"S1"}),
        post_mapping_samples=frozenset({"S1"}),
    )

    with pytest.raises(artifact_module._UnknownSampleError):
        artifact_module._validate_multiqc_sample_rows(
            _multiqc_stub("UNKNOWN"),
            policy,
        )


def test_multiqc_table_duplicate_and_malformed_rows_fail_closed():
    policy = artifact_module._MultiqcRowPolicy(
        owners={"S1": "S1"},
        exact_rows=frozenset({"S1"}),
        required_owners=frozenset(),
    )

    for content in (
        b"Sample\tmetric\nS1\t1\nS1\t2\n",
        b"wrong\tmetric\nS1\t1\n",
        b"Sample\tmetric\nS1\n",
        b"Sample\tmetric\textra\nS1\t1\n",
        b"Sample\tmetric\nS1\t1\textra\n",
    ):
        with pytest.raises(artifact_module._MultiqcTableError):
            artifact_module._validate_multiqc_sample_rows(content, policy)


@pytest.mark.parametrize(
    ("layout", "params", "expected"),
    (
        (
            "SE",
            {
                "skip_fastqc": False,
                "skip_trimming": False,
                "trimmer": "trimgalore",
                "remove_ribo_rna": False,
            },
            frozenset({"S"}),
        ),
        (
            "PE",
            {
                "skip_fastqc": True,
                "skip_trimming": True,
                "trimmer": "trimgalore",
                "remove_ribo_rna": False,
            },
            frozenset({"S"}),
        ),
        (
            "PE",
            {
                "skip_fastqc": False,
                "skip_trimming": True,
                "trimmer": "fastp",
                "remove_ribo_rna": False,
            },
            frozenset({"S", "S Read 1", "S Read 2"}),
        ),
        (
            "PE",
            {
                "skip_fastqc": True,
                "skip_trimming": False,
                "trimmer": "trimgalore",
                "remove_ribo_rna": False,
            },
            frozenset({"S", "S Read 1", "S Read 2"}),
        ),
        (
            "PE",
            {
                "skip_fastqc": True,
                "skip_trimming": False,
                "trimmer": "trimgalore",
                "remove_ribo_rna": False,
                "with_umi": True,
                "skip_umi_extract": False,
                "umi_discard_read": 2,
            },
            frozenset({"S"}),
        ),
        (
            "PE",
            {
                "skip_fastqc": False,
                "skip_trimming": False,
                "trimmer": "trimgalore",
                "remove_ribo_rna": False,
                "with_umi": True,
                "skip_umi_extract": False,
                "umi_discard_read": 1,
            },
            frozenset({"S", "S Read 1", "S Read 2"}),
        ),
    ),
)
def test_general_stats_uses_exact_fixed_multiqc_grouped_identities(
    layout: str,
    params: dict[str, object],
    expected: frozenset[str],
):
    samples = ({"sample": "S", "layout": layout, "strandedness": "unstranded"},)
    policy = artifact_module._multiqc_row_policy(
        "bulk_rnaseq.multiqc.general_stats",
        params=params,
        samples=samples,
        sample_ids=frozenset({"S"}),
        analysis_samples=frozenset({"S"}),
        post_mapping_samples=frozenset({"S"}),
    )

    assert policy.exact_rows == expected
    artifact_module._validate_multiqc_sample_rows(
        _multiqc_stub(*sorted(expected)),
        policy,
    )
    for invalid in (
        frozenset({"S", "S_1", "S_2"}),
        expected - {"S Read 2"},
    ):
        if invalid == expected:
            continue
        with pytest.raises(
            (artifact_module._MultiqcTableError, artifact_module._UnknownSampleError)
        ):
            artifact_module._validate_multiqc_sample_rows(
                _multiqc_stub(*sorted(invalid)),
                policy,
            )


def test_general_stats_accepts_1000_paired_samples_and_rejects_3001_rows():
    sample_ids = frozenset(f"S{index:04d}" for index in range(MAX_SAMPLE_ROWS))
    samples = tuple(
        {"sample": sample, "layout": "PE", "strandedness": "unstranded"}
        for sample in reversed(sorted(sample_ids))
    )
    params = {
        "skip_fastqc": False,
        "skip_trimming": False,
        "trimmer": "trimgalore",
        "remove_ribo_rna": False,
    }
    policy = artifact_module._multiqc_row_policy(
        "bulk_rnaseq.multiqc.general_stats",
        params=params,
        samples=samples,
        sample_ids=sample_ids,
        analysis_samples=sample_ids,
        post_mapping_samples=sample_ids,
    )
    expected = tuple(sorted(policy.exact_rows or ()))

    assert len(expected) == MAX_SAMPLE_ROWS * 3
    artifact_module._validate_multiqc_sample_rows(_multiqc_stub(*expected), policy)
    with pytest.raises(artifact_module._ArtifactLimitError):
        artifact_module._validate_multiqc_sample_rows(
            _multiqc_stub(*expected, "S0000 Read 3"),
            policy,
        )


def test_general_stats_accepts_maximum_length_paired_sample_identity():
    sample = "S" * 128
    params = {
        "skip_fastqc": False,
        "skip_trimming": True,
        "trimmer": "fastp",
        "remove_ribo_rna": False,
    }
    samples = ({"sample": sample, "layout": "PE", "strandedness": "unstranded"},)
    policy = artifact_module._multiqc_row_policy(
        "bulk_rnaseq.multiqc.general_stats",
        params=params,
        samples=samples,
        sample_ids=frozenset({sample}),
        analysis_samples=frozenset({sample}),
        post_mapping_samples=frozenset({sample}),
    )

    artifact_module._validate_multiqc_sample_rows(
        _multiqc_stub(*sorted(policy.exact_rows or ())),
        policy,
    )


def test_general_stats_mixed_layout_is_exact_and_status_independent():
    rows = (
        {"sample": "PAIR", "layout": "PE", "strandedness": "unstranded"},
        {"sample": "SINGLE", "layout": "SE", "strandedness": "unstranded"},
        # Repeated lanes normalize to the same biological owner and must not
        # create a second MultiQC identity.
        {"sample": "PAIR", "layout": "PE", "strandedness": "unstranded"},
    )
    params = {
        "skip_fastqc": False,
        "skip_trimming": False,
        "trimmer": "trimgalore",
        "remove_ribo_rna": False,
    }
    expected = frozenset({"PAIR", "PAIR Read 1", "PAIR Read 2", "SINGLE"})

    policies = tuple(
        artifact_module._multiqc_row_policy(
            "bulk_rnaseq.multiqc.general_stats",
            params=params,
            samples=ordered,
            sample_ids=frozenset({"PAIR", "SINGLE"}),
            # Simulate one trim failure and one later mapping failure: neither
            # threshold can remove the pre-threshold General Stats identities.
            analysis_samples=frozenset({"SINGLE"}),
            post_mapping_samples=frozenset(),
        )
        for ordered in (rows, tuple(reversed(rows)))
    )

    assert all(policy.exact_rows == expected for policy in policies)
    assert policies[0] == policies[1]
    artifact_module._validate_multiqc_sample_rows(
        _multiqc_stub(*reversed(sorted(expected))),
        policies[0],
    )


def test_multiqc_foreign_row_and_post_audit_replacement_are_fail_closed(
    tmp_path: Path,
    monkeypatch,
):
    inputs = _inputs(qc=_multiqc_only())
    _write_core(tmp_path, ("S1",))
    root = "multiqc/star_salmon/multiqc_report_data"
    for name in ("multiqc_general_stats.txt", "multiqc_star.txt", "multiqc_salmon.txt"):
        _write(tmp_path, f"{root}/{name}", _multiqc_stub("S1"))
    star = tmp_path / "results" / root / "multiqc_star.txt"

    star.write_bytes(_multiqc_stub("X1"))
    foreign = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert foreign.is_failure
    assert foreign.errors[0].code == "BULK_RNASEQ_ARTIFACT_UNKNOWN_SAMPLE"

    star.write_bytes(_multiqc_stub("S1"))
    original = artifact_module._verified_regular_file
    replaced = False

    def replace_after_audit(results_descriptor: int, relative_path: str) -> bool:
        nonlocal replaced
        if not replaced:
            replaced = True
            star.write_bytes(_multiqc_stub("X1"))
        return original(results_descriptor, relative_path)

    monkeypatch.setattr(
        artifact_module,
        "_verified_regular_file",
        replace_after_audit,
    )
    raced = discover_bulk_rnaseq_artifacts(inputs, tmp_path)
    assert raced.is_failure
    assert raced.errors[0].code == "BULK_RNASEQ_ARTIFACT_RACE_DETECTED"


@pytest.mark.parametrize("kind", ("directory", "fifo", "symlink"))
def test_published_multiqc_table_must_be_a_regular_nofollow_file(
    tmp_path: Path,
    kind: str,
):
    inputs = _inputs(qc=_multiqc_only())
    _write_core(tmp_path, ("S1",))
    root = tmp_path / "results/multiqc/star_salmon/multiqc_report_data"
    for name in ("multiqc_general_stats.txt", "multiqc_star.txt", "multiqc_salmon.txt"):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1"),
        )
    target = root / "multiqc_star.txt"
    target.unlink()
    if kind == "directory":
        target.mkdir()
    elif kind == "fifo":
        os.mkfifo(target)
    else:
        private = root / "private.tsv"
        private.write_bytes(_multiqc_stub("S1"))
        target.symlink_to(private.name)

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_PATH_UNSAFE"


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
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1"),
        )

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


def test_rseqc_inner_distance_single_end_placeholder_is_not_public(
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
    inputs = _inputs(qc=qc, advanced={"rseqc_modules": "inner_distance"})
    _write_core(tmp_path, ("S1",))
    _write(
        tmp_path,
        "star_salmon/rseqc/inner_distance/txt/S1.inner_distance.txt",
        b"inner_distance.py doesn't support single-end data\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success
    assert not any(
        candidate.output_type.startswith("bulk_rnaseq.rseqc.inner_distance.")
        for candidate in result.value
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
    for trimmed in ("1_val_1", "2_val_2"):
        root = f"S1_trimmed_{trimmed}_fastqc"
        _write(tmp_path, f"fastqc/trim/{root}.html")
        _write(
            tmp_path,
            f"fastqc/trim/{root}.zip",
            _fastqc_evidence_zip(
                root,
                f"S1_trimmed_{trimmed}.fq.gz",
                20000,
            ),
        )
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
        trimming={"enabled": True, "tool": "fastp"},
        qc=_multiqc_only(),
    )
    _write_core(tmp_path, ("S2",))
    _write_fastp_evidence(tmp_path, "S1", 9999)
    _write_fastp_evidence(tmp_path, "S2", 20000)
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    for name in (
        "multiqc_general_stats.txt",
        "multiqc_fastp.txt",
        "multiqc_star.txt",
        "multiqc_salmon.txt",
    ):
        rows = (
            ("S2",)
            if name in {"multiqc_star.txt", "multiqc_salmon.txt"}
            else (
                "S1",
                "S2",
            )
        )
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub(*rows),
        )
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
        trimming={"enabled": True, "tool": "fastp"},
        qc=_multiqc_only(),
    )
    _write(tmp_path, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    _write_fastp_evidence(tmp_path, "S1", 9999)
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    for name in ("multiqc_general_stats.txt", "multiqc_fastp.txt"):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1"),
        )
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


@pytest.mark.parametrize("cutadapt_reads_written", (20000, 9999))
def test_trimgalore_pe_failure_uses_only_valid_cutadapt_evidence(
    tmp_path: Path,
    cutadapt_reads_written: int,
):
    qc = {**_multiqc_only(), "fastqc": True}
    inputs = _inputs(
        layouts=("PE",),
        trimming={"enabled": True, "tool": "trimgalore"},
        qc=qc,
    )
    _write(tmp_path, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    for number, trimmed in (("1", "1_val_1"), ("2", "2_val_2")):
        _write(tmp_path, f"fastqc/raw/S1_raw_{number}_fastqc.html")
        _write(tmp_path, f"fastqc/raw/S1_raw_{number}_fastqc.zip")
        root = f"S1_trimmed_{trimmed}_fastqc"
        _write(tmp_path, f"fastqc/trim/{root}.html")
        _write(
            tmp_path,
            f"fastqc/trim/{root}.zip",
            _fastqc_evidence_zip(root, f"S1_trimmed_{trimmed}.fq.gz", 9999),
        )
    for name, identities in (
        ("multiqc_general_stats.txt", ("S1", "S1 Read 1", "S1 Read 2")),
        ("multiqc_cutadapt.txt", ("S1_1", "S1_2")),
        ("multiqc_fastqc_fastqc_raw.txt", ("S1_1", "S1_2")),
        ("multiqc_fastqc_fastqc_trimmed.txt", ("S1_1", "S1_2")),
    ):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub(*identities),
        )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_cutadapt.txt",
        b"Sample\tcutadapt_version\tr_processed\tr_with_adapters\tr_written\t"
        b"bp_processed\tquality_trimmed\tbp_written\tpercent_trimmed\n"
        + (
            f"S1_1\t4.0\t20000\t0\t{cutadapt_reads_written}\t"
            "2000000\t0\t999900\t50.005\n"
            f"S1_2\t4.0\t20000\t0\t{cutadapt_reads_written}\t"
            "2000000\t0\t999900\t50.005\n"
        ).encode(),
    )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/"
        "multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t9999\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    if cutadapt_reads_written != 20000:
        assert result.is_failure
        assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_STATUS_INVALID"
        return
    assert result.is_success
    assert "bulk_rnaseq.multiqc.fail_trimmed_samples" in _types(result)
    assert not any(
        item.output_type.startswith(("bulk_rnaseq.star.", "bulk_rnaseq.salmon."))
        for item in result.value
    )


def test_trimgalore_low_yield_remains_discoverable_when_raw_fastqc_is_disabled(
    tmp_path: Path,
):
    inputs = _inputs(
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={**_multiqc_only(), "fastqc": False},
    )
    _write(tmp_path, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    root = "S1_trimmed_trimmed_fastqc"
    _write(tmp_path, f"fastqc/trim/{root}.html")
    _write(
        tmp_path,
        f"fastqc/trim/{root}.zip",
        _fastqc_evidence_zip(root, "S1_trimmed_trimmed.fq.gz", 9999),
    )
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    for name, identities in (
        ("multiqc_general_stats.txt", ("S1",)),
        ("multiqc_fastqc_fastqc_trimmed.txt", ("S1",)),
    ):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub(*identities),
        )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_cutadapt.txt",
        b"Sample\tcutadapt_version\tr_processed\tr_with_adapters\tr_written\t"
        b"bp_processed\tquality_trimmed\tbp_written\tpercent_trimmed\n"
        b"S1\t4.0\t20000\t0\t20000\t2000000\t0\t999900\t50.005\n",
    )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/"
        "multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t9999\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success
    types = _types(result)
    assert "bulk_rnaseq.fastqc.trimmed.single.html" in types
    assert "bulk_rnaseq.fastqc.trimmed.single.zip" in types
    assert "bulk_rnaseq.multiqc.fastqc.trimmed" in types
    assert not any(".fastqc.raw." in value for value in types)
    assert not any(
        value.startswith(("bulk_rnaseq.star.", "bulk_rnaseq.salmon."))
        for value in types
    )


def test_cutadapt_same_length_replacement_after_reconciliation_is_a_race(
    tmp_path: Path,
    monkeypatch,
):
    inputs = _inputs(
        layouts=("PE",),
        trimming={"enabled": True, "tool": "trimgalore"},
        qc={**_multiqc_only(), "fastqc": True},
    )
    _write(tmp_path, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    for number, trimmed in (("1", "1_val_1"), ("2", "2_val_2")):
        _write(tmp_path, f"fastqc/raw/S1_raw_{number}_fastqc.html")
        _write(tmp_path, f"fastqc/raw/S1_raw_{number}_fastqc.zip")
        root = f"S1_trimmed_{trimmed}_fastqc"
        _write(tmp_path, f"fastqc/trim/{root}.html")
        _write(
            tmp_path,
            f"fastqc/trim/{root}.zip",
            _fastqc_evidence_zip(root, f"S1_trimmed_{trimmed}.fq.gz", 9999),
        )
    for name, identities in (
        ("multiqc_general_stats.txt", ("S1", "S1 Read 1", "S1 Read 2")),
        ("multiqc_fastqc_fastqc_raw.txt", ("S1_1", "S1_2")),
        ("multiqc_fastqc_fastqc_trimmed.txt", ("S1_1", "S1_2")),
    ):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub(*identities),
        )
    cutadapt = _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_cutadapt.txt",
        b"Sample\tcutadapt_version\tr_processed\tr_with_adapters\tr_written\t"
        b"bp_processed\tquality_trimmed\tbp_written\tpercent_trimmed\n"
        b"S1_1\t4.0\t20000\t0\t20000\t2000000\t0\t999900\t50.005\n"
        b"S1_2\t4.0\t20000\t0\t20000\t2000000\t0\t999900\t50.005\n",
    )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/"
        "multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t9999\n",
    )
    expected_specs = artifact_module._expected_specs
    replaced = False

    def replace_after_reconciliation(*args, **kwargs):
        nonlocal replaced
        result = expected_specs(*args, **kwargs)
        if not replaced:
            replaced = True
            original = cutadapt.read_bytes()
            cutadapt.write_bytes(original.replace(b"\t0\t20000\t", b"\t1\t20000\t"))
            assert cutadapt.stat().st_size == len(original)
        return result

    monkeypatch.setattr(
        artifact_module,
        "_expected_specs",
        replace_after_reconciliation,
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_RACE_DETECTED"


@pytest.mark.parametrize(
    ("field", "value"),
    (
        ("passed_filter_reads", 10000),
        ("before_total_reads", 9998),
        ("after_total_bases", 1000001),
    ),
)
def test_fastp_status_evidence_rejects_impossible_fixed_report(
    tmp_path: Path,
    field: str,
    value: int,
):
    inputs = _inputs(
        trimming={"enabled": True, "tool": "fastp"},
        qc=_multiqc_only(),
    )
    _write(tmp_path, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    _write_fastp_evidence(tmp_path, "S1", 9999)
    report = tmp_path / "results/fastp/S1.fastp.json"
    payload = json.loads(report.read_text())
    if field == "passed_filter_reads":
        payload["filtering_result"]["passed_filter_reads"] = value
    elif field == "before_total_reads":
        payload["summary"]["before_filtering"]["total_reads"] = value
    else:
        payload["summary"]["after_filtering"]["total_bases"] = value
    report.write_text(json.dumps(payload))
    for name in ("multiqc_general_stats.txt", "multiqc_fastp.txt"):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1"),
        )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/"
        "multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t9999\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_STATUS_INVALID"


def test_mapped_failure_table_keeps_star_salmon_but_not_post_mapping_outputs(
    tmp_path: Path,
):
    inputs = _inputs(
        layouts=("SE", "SE"),
        qc=_multiqc_only(),
        advanced={"min_mapped_reads": 90},
    )
    _write_core(tmp_path, ("S1", "S2"))
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    for name in (
        "multiqc_general_stats.txt",
        "multiqc_star.txt",
        "multiqc_salmon.txt",
        "multiqc_fastp.txt",
    ):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1", "S2") if name != "multiqc_fastp.txt" else b"stale",
        )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_mapped_samples_table.txt",
        b"Sample\tSTAR uniquely mapped reads (%)\nS1\t80.00\nS2\t80.00\n",
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


@pytest.mark.parametrize(
    "star_log",
    [
        b"Uniquely mapped reads % | 80.00%\n",
        _STAR_LOG.replace(
            b"Number of reads unmapped: other | 10",
            b"Number of reads unmapped: other | 11",
        ),
    ],
    ids=("truncated", "impossible-count-partition"),
)
def test_mapped_status_requires_complete_consistent_star_log(
    tmp_path: Path,
    star_log: bytes,
):
    inputs = _inputs(qc=_multiqc_only(), advanced={"min_mapped_reads": 90})
    _write_core(tmp_path, ("S1",))
    _write(tmp_path, "star_salmon/log/S1.Log.final.out", star_log)
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_mapped_samples_table.txt",
        b"Sample\tSTAR uniquely mapped reads (%)\nS1\t80.00\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_STATUS_INVALID"
    assert result.errors[0].context == {"reason": "sample_status_invalid"}


def test_mapped_failure_keeps_published_umi_and_pre_umi_bams(tmp_path: Path):
    inputs = _inputs(
        qc=_multiqc_only(),
        outputs={"umi_intermediates": True},
        advanced={"min_mapped_reads": 90},
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
    for name in (
        "multiqc_general_stats.txt",
        "multiqc_star.txt",
        "multiqc_salmon.txt",
        "multiqc_fastp.txt",
    ):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1") if name != "multiqc_fastp.txt" else b"stale",
        )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_mapped_samples_table.txt",
        b"Sample\tSTAR uniquely mapped reads (%)\nS1\t80.00\n",
    )

    types = _types(discover_bulk_rnaseq_artifacts(inputs, tmp_path))

    assert "bulk_rnaseq.star.bam" in types
    assert "bulk_rnaseq.star.sorted_intermediate_bam" in types
    assert "bulk_rnaseq.star.transcriptome_sorted_intermediate_bam" in types
    assert "bulk_rnaseq.umi.transcriptome_sorted_bam" in types


def test_exact_trim_threshold_row_does_not_hide_existing_analysis_outputs(
    tmp_path: Path,
):
    inputs = _inputs(
        trimming={"enabled": True, "tool": "fastp"},
        qc=_multiqc_only(),
    )
    _write_core(tmp_path, ("S1",))
    _write_fastp_evidence(tmp_path, "S1", 10000)
    for name in (
        "multiqc_general_stats.txt",
        "multiqc_star.txt",
        "multiqc_salmon.txt",
        "multiqc_fastp.txt",
    ):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1"),
        )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t10000.0\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_success
    assert "bulk_rnaseq.star.bam" in _types(result)


def test_trimming_disabled_rejects_stale_failure_table(tmp_path: Path):
    inputs = _inputs(qc=_multiqc_only())
    _write(tmp_path, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_general_stats.txt",
    )
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t9999\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_STATUS_INVALID"


def test_trimming_and_multiqc_disabled_reject_stale_failure_table(tmp_path: Path):
    inputs = _inputs()
    _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/"
        "multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t9999\n",
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_STATUS_INVALID"
    assert result.errors[0].context == {"reason": "sample_status_invalid"}


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
    _write_fastp_evidence(tmp_path, "S1", 9999)
    _write(tmp_path, "fastqc/raw/S1_raw_fastqc.html")
    _write(tmp_path, "fastqc/raw/S1_raw_fastqc.zip")
    for name in (
        "multiqc_general_stats.txt",
        "multiqc_fastp.txt",
        "multiqc_fastqc_fastqc_raw.txt",
    ):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1"),
        )
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


def test_status_table_same_length_replacement_after_reconciliation_is_a_race(
    tmp_path: Path,
    monkeypatch,
):
    inputs = _inputs(
        trimming={"enabled": True, "tool": "fastp"},
        qc=_multiqc_only(),
    )
    _write(tmp_path, "pipeline_info/nf_core_rnaseq_software_mqc_versions.yml")
    _write_fastp_evidence(tmp_path, "S1", 9999)
    for name in ("multiqc_general_stats.txt", "multiqc_fastp.txt"):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1"),
        )
    status = _write(
        tmp_path,
        "multiqc/star_salmon/multiqc_report_data/"
        "multiqc_fail_trimmed_samples_table.txt",
        b"Sample\tReads after trimming\nS1\t9999\n",
    )
    expected_specs = artifact_module._expected_specs
    replaced = False

    def replace_after_reconciliation(*args, **kwargs):
        nonlocal replaced
        result = expected_specs(*args, **kwargs)
        if not replaced:
            replaced = True
            status.write_bytes(b"Sample\tReads after trimming\nS1\t9998\n")
        return result

    monkeypatch.setattr(
        artifact_module,
        "_expected_specs",
        replace_after_reconciliation,
    )

    result = discover_bulk_rnaseq_artifacts(inputs, tmp_path)

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ARTIFACT_RACE_DETECTED"
    assert result.errors[0].context == {"reason": "artifact_identity_changed"}


def test_missing_core_output_without_a_failure_table_still_fails(tmp_path: Path):
    inputs = _inputs(qc=_multiqc_only())
    _write_core(tmp_path, ("S1",))
    (tmp_path / "results/star_salmon/S1/quant.sf").unlink()
    _write(tmp_path, "multiqc/star_salmon/multiqc_report.html")
    for name in ("multiqc_general_stats.txt", "multiqc_star.txt", "multiqc_salmon.txt"):
        _write(
            tmp_path,
            f"multiqc/star_salmon/multiqc_report_data/{name}",
            _multiqc_stub("S1"),
        )

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
