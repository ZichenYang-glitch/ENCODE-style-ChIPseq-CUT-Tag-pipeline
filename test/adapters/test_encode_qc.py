"""Strict ENCODE QC summary mapping tests."""

from __future__ import annotations

import csv
from decimal import Decimal
import io

import pytest

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import (
    QcSourceArtifact,
    QcSourceDocument,
    WorkflowInputs,
)


QC_HEADER = (
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

MNASE_HEADER = (
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

POOLED_HEADER = (
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


def _tsv(header: tuple[str, ...], row: dict[str, str], *, extra_row=False) -> bytes:
    stream = io.StringIO(newline="")
    writer = csv.DictWriter(
        stream, fieldnames=header, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerow({key: row.get(key, "NA") for key in header})
    if extra_row:
        writer.writerow({key: row.get(key, "NA") for key in header})
    return stream.getvalue().encode()


def _source(
    artifact_id: str,
    output_type: str,
    content: bytes,
    *,
    sample_id: str | None = None,
    experiment_id: str | None = None,
    assay: str,
) -> QcSourceDocument:
    metadata = {"assay": assay}
    if sample_id is not None:
        metadata.update({"scope": "sample", "sample_id": sample_id})
    elif experiment_id is not None:
        metadata.update({"scope": "experiment", "experiment_id": experiment_id})
    if sample_id is not None and experiment_id is not None:
        metadata["experiment_id"] = experiment_id
    return QcSourceDocument(
        source=QcSourceArtifact(
            artifact_id=artifact_id,
            output_type=output_type,
            relative_path=f"results/{artifact_id}.tsv",
            metadata=metadata,
        ),
        content=content,
    )


def _sample_document(**overrides: str) -> QcSourceDocument:
    row = {
        "sample": "S1",
        "assay": "chipseq",
        "total_reads": "1000",
        "frip": "0.125",
        "peak_count": "50",
        "percent_duplication": "0.2",
        "estimated_library_size": "9007199254740993",
        "nrf": "0.8",
        "pbc1": "0.75",
        "pbc2": "3.0",
    }
    row.update(overrides)
    return _source(
        "artifact-sample",
        "qc_summary",
        _tsv(QC_HEADER, row),
        sample_id="S1",
        experiment_id="EXP1",
        assay="chipseq",
    )


def test_encode_qc_parser_maps_locked_sample_mnase_and_pooled_metrics():
    mnase = _source(
        "artifact-mnase",
        "mnase_qc_summary",
        _tsv(
            MNASE_HEADER,
            {
                "sample": "M1",
                "assay": "mnase",
                "sub_reads": "10",
                "mono_reads": "20",
                "di_reads": "30",
            },
        ),
        sample_id="M1",
        assay="mnase",
    )
    pooled = _source(
        "artifact-pooled",
        "pooled_qc_summary",
        _tsv(
            POOLED_HEADER,
            {
                "experiment": "EXP1",
                "assay": "chipseq",
                "peak_mode_status": "mismatch",
                "biological_replicates": "2",
                "pooled_peak_count": "75",
            },
        ),
        experiment_id="EXP1",
        assay="chipseq",
    )

    result = EncodeStyleWorkflowAdapter().extract_qc_metrics(
        WorkflowInputs(config={}),
        (_sample_document(), mnase, pooled),
    )

    assert result.is_success
    by_key = {
        (item.metric_key, item.sample_id, item.experiment_id): item
        for item in result.value
    }
    assert len(result.value) == 13
    assert by_key[("sequencing.total_reads", "S1", "EXP1")].value == Decimal("1000")
    assert by_key[("library.estimated_size", "S1", "EXP1")].value == Decimal(
        "9007199254740993"
    )
    assert by_key[("mnase.mononucleosomal_reads", "M1", None)].unit == "count"
    pooled_metric = by_key[("peaks.pooled_count", None, "EXP1")]
    assert pooled_metric.value == Decimal("75")
    assert pooled_metric.scope == "experiment"
    assert pooled_metric.qc_flag == "warning"
    assert by_key[("replicates.biological_count", None, "EXP1")].qc_flag is None


def test_encode_qc_parser_skips_explicit_na_without_failing_source():
    result = EncodeStyleWorkflowAdapter().extract_qc_metrics(
        WorkflowInputs(config={}),
        (_sample_document(frip="NA", estimated_library_size="NA"),),
    )

    assert result.is_success
    assert {item.metric_key for item in result.value} == {
        "sequencing.total_reads",
        "peaks.count",
        "library.percent_duplication",
        "library.nrf",
        "library.pbc1",
        "library.pbc2",
    }


@pytest.mark.parametrize(
    ("column", "value"),
    [
        ("total_reads", "-1"),
        ("total_reads", "1.5"),
        ("frip", "1.01"),
        ("frip", "NaN"),
        ("frip", "Infinity"),
        ("frip", "1e-3"),
        ("frip", "0.1234567890123"),
    ],
)
def test_encode_qc_parser_fails_closed_for_invalid_numbers(column, value):
    result = EncodeStyleWorkflowAdapter().extract_qc_metrics(
        WorkflowInputs(config={}),
        (_sample_document(**{column: value}),),
    )

    assert result.is_failure
    assert result.issues[0].code == "ENCODE_QC_SUMMARY_INVALID"
    assert value not in str(result.issues[0].to_dict())


@pytest.mark.parametrize(
    "document",
    [
        _source(
            "artifact-bad-header",
            "qc_summary",
            _tsv(QC_HEADER[:-1], {"sample": "S1", "assay": "chipseq"}),
            sample_id="S1",
            assay="chipseq",
        ),
        _source(
            "artifact-extra-row",
            "qc_summary",
            _tsv(QC_HEADER, {"sample": "S1", "assay": "chipseq"}, extra_row=True),
            sample_id="S1",
            assay="chipseq",
        ),
        _source(
            "artifact-bad-utf8",
            "qc_summary",
            b"\xff\xfe",
            sample_id="S1",
            assay="chipseq",
        ),
        _source(
            "artifact-wrong-sample",
            "qc_summary",
            _tsv(QC_HEADER, {"sample": "OTHER", "assay": "chipseq"}),
            sample_id="S1",
            assay="chipseq",
        ),
    ],
)
def test_encode_qc_parser_rejects_malformed_contract_or_identity(document):
    result = EncodeStyleWorkflowAdapter().extract_qc_metrics(
        WorkflowInputs(config={}),
        (document,),
    )

    assert result.is_failure
    assert result.issues[0].code == "ENCODE_QC_SUMMARY_INVALID"


def test_encode_qc_parser_rejects_duplicate_source_artifacts():
    source = _sample_document()

    result = EncodeStyleWorkflowAdapter().extract_qc_metrics(
        WorkflowInputs(config={}),
        (source, source),
    )

    assert result.is_failure


@pytest.mark.parametrize(
    ("status", "expected"),
    [("ok", "pass"), ("mismatch", "warning"), ("unknown", None)],
)
def test_encode_qc_parser_maps_only_existing_pooled_status_vocabulary(status, expected):
    source = _source(
        "artifact-pooled",
        "pooled_qc_summary",
        _tsv(
            POOLED_HEADER,
            {
                "experiment": "EXP1",
                "assay": "atac",
                "peak_mode_status": status,
                "biological_replicates": "2",
                "pooled_peak_count": "50",
            },
        ),
        experiment_id="EXP1",
        assay="atac",
    )

    result = EncodeStyleWorkflowAdapter().extract_qc_metrics(
        WorkflowInputs(config={}), (source,)
    )

    assert result.is_success
    assert (
        next(
            item for item in result.value if item.metric_key == "peaks.pooled_count"
        ).qc_flag
        == expected
    )
