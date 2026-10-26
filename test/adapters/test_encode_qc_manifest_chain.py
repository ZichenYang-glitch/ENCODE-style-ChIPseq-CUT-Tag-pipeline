"""Config -> versioned files -> manifest -> adapter QC consumption chain."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from decimal import Decimal

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.manifest.make import build_manifest_rows, write_manifest_tsv
from encode_pipeline.platform.adapters import (
    QcSourceArtifact,
    QcSourceDocument,
    WorkflowInputs,
)
from encode_pipeline.qc.contract import (
    FRAGMENT_QC_V2_COLUMNS,
    MITO_QC_V2_COLUMNS,
    SAMPLE_QC_V2_COLUMNS,
)


def _write_tsv(path: Path, header, row) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerow([row.get(column, "NA") for column in header])


def _policy() -> dict:
    return {
        "schema_version": 1,
        "scientific_policy_revision": 5,
        "sample": "S1",
        "assay": "cuttag",
        "profile": "updated",
        "resolved_duplicate_analysis": "retain",
        "resolved_max_insert": 700,
        "organelle_mode": "report",
        "organelle_status": "ok",
        "matched_organelle_contigs": ["chrM"],
        "blacklist_mode": "off",
        "blacklist_applied": False,
        "input_stages": {
            "fragment_size": "analysis_bam",
            "mnase_classification": "not_applicable",
            "duplicate_metrics": "pre_duplicate_markdup_bam",
            "nrf_pbc": "pre_duplicate_mapq_bam",
            "organelle": "pre_duplicate_mapq_bam",
            "blacklist": "not_applicable",
            "peaks": "analysis_bam",
        },
        "stage_counts": {
            "analysis_input_alignments": 20,
            "analysis_kept_alignments": 18,
            "analysis_valid_templates": 9,
            "analysis_invalid_templates": 1,
            "upstream_dropped_templates": 1,
            "upstream_dropped_alignments": 2,
            "blacklist_dropped_templates": 0,
            "blacklist_dropped_alignments": 0,
            "organelle_reads": 2,
            "organelle_total_mapped_reads": 20,
        },
        "reference_summary": {"genome": "tiny", "organelle_chroms_explicit": False},
        "qc_schema_versions": {
            "sample": "2",
            "mnase": "2",
            "fragment": "2",
            "mito": "2",
        },
    }


def test_config_files_manifest_and_adapter_share_one_qc_contract(tmp_path):
    results = tmp_path / "results"
    samples = tmp_path / "samples.tsv"
    samples.write_text(
        "\t".join(
            [
                "sample",
                "fastq_1",
                "fastq_2",
                "layout",
                "assay",
                "target",
                "peak_mode",
                "genome",
                "bowtie2_index",
                "experiment",
                "biological_replicate",
            ]
        )
        + "\n"
        + "S1\tR1.fq\tR2.fq\tPE\tcuttag\tH3K4me3\tnarrow\ttiny\t/idx\tEXP1\t1\n",
        encoding="utf-8",
    )
    config = {
        "samples": str(samples),
        "outdir": str(results),
        "use_control": False,
        "multiqc": False,
        "replicate_analysis": False,
        "chipseq_idr": False,
        "analysis_profile": "updated",
        "blacklist_mode": "off",
        "qc": {
            "summary": True,
            "cuttag_fragment_size": True,
            "signal_tracks": False,
        },
        "genome_resources": {"tiny": {"effective_genome_size": 10000}},
    }

    sample_row = {column: "NA" for column in SAMPLE_QC_V2_COLUMNS}
    sample_row.update(
        {
            "schema_version": "2",
            "sample": "S1",
            "assay": "cuttag",
            "target": "H3K4me3",
            "genome": "tiny",
            "layout": "PE",
            "peak_mode": "narrow",
            "total_reads": "1000",
            "frip": "0.1",
            "peak_count": "10",
            "analysis_input_stage": "analysis_bam",
            "analysis_valid_fragments": "9",
            "analysis_invalid_templates": "1",
            "raw_peak_count": "10",
            "delivered_peak_count": "10",
            "duplicate_analysis": "retain",
            "resolved_max_insert": "700",
            "blacklist_mode": "off",
            "blacklist_applied": "false",
            "organelle_mode": "report",
            "organelle_status": "ok",
            "organelle_reads": "2",
            "organelle_total_mapped_reads": "20",
            "organelle_fraction": "0.1",
            "cross_correlation_status": "not_enabled",
            "scientific_policy_revision": "5",
        }
    )
    _write_tsv(
        results / "S1/01_qc/S1.qc_summary.tsv",
        SAMPLE_QC_V2_COLUMNS,
        sample_row,
    )
    (results / "S1/05_qc").mkdir(parents=True)
    (results / "S1/05_qc/S1.analysis_policy.json").write_text(
        json.dumps(_policy(), indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    _write_tsv(
        results / "S1/05_qc/S1.mito_fraction.tsv",
        MITO_QC_V2_COLUMNS,
        {
            "schema_version": "2",
            "sample": "S1",
            "organelle_reads": "2",
            "total_mapped_reads": "20",
            "fraction": "0.1",
            "mode": "report",
            "status": "ok",
            "input_stage": "pre_duplicate_mapq_bam",
            "matched_contigs": "chrM",
            "organelle_chroms_explicit": "false",
        },
    )
    _write_tsv(
        results / "S1/01_qc/S1.cuttag_fragment_size.tsv",
        FRAGMENT_QC_V2_COLUMNS,
        {
            "schema_version": "2",
            "sample": "S1",
            "layout": "PE",
            "fragment_count": "9",
            "status": "ok",
            "input_stage": "analysis_bam",
            "duplicate_analysis": "retain",
            "blacklist_applied": "false",
            "templates_total": "10",
            "templates_valid": "9",
            "templates_invalid": "1",
            "invalid_tlen_zero": "0",
            "invalid_orphan": "0",
            "invalid_cross_contig": "0",
            "invalid_inconsistent": "1",
            "alignments_total": "18",
            "upstream_dropped_templates": "1",
            "upstream_dropped_alignments": "2",
            "upstream_alignments_total": "20",
            "upstream_alignments_kept": "18",
        },
    )

    rows, _, _ = build_manifest_rows(config)
    manifest = tmp_path / "manifest.tsv"
    write_manifest_tsv(rows, manifest)
    with manifest.open(newline="", encoding="utf-8") as handle:
        manifest_rows = list(csv.DictReader(handle, delimiter="\t"))

    assert {
        row["output_type"] for row in manifest_rows if row["status"] == "present"
    } >= {
        "qc_summary",
        "analysis_policy",
        "mito_fraction",
        "cuttag_fragment_size",
    }

    adapter = EncodeStyleWorkflowAdapter()
    documents = []
    for row in manifest_rows:
        if row["output_type"] not in adapter.qc_source_output_types():
            continue
        path = Path(row["path"])
        metadata = {
            "scope": "sample",
            "sample_id": row["sample_id"],
            "assay": row["assay"],
        }
        if row["experiment_id"]:
            metadata["experiment_id"] = row["experiment_id"]
        documents.append(
            QcSourceDocument(
                source=QcSourceArtifact(
                    artifact_id=f"{row['output_type']}-{row['sample_id']}",
                    output_type=row["output_type"],
                    relative_path=str(path.relative_to(tmp_path)),
                    metadata=metadata,
                ),
                content=path.read_bytes(),
            )
        )

    result = adapter.extract_qc_metrics(WorkflowInputs(config={}), tuple(documents))
    assert result.is_success
    by_key = {item.metric_key: item for item in result.value}
    assert by_key["analysis.valid_fragments"].value == Decimal("9")
    assert by_key["organelle.fraction"].value == Decimal("0.1")
    assert by_key["analysis.stage.valid_templates"].value == Decimal("9")
    assert by_key["analysis.stage.blacklist_dropped_templates"].value == Decimal("0")
