"""Tests for the protected Bulk RNA-seq browser fixture projection."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import pytest

from encode_pipeline.adapters.bulk_rnaseq import BulkRnaSeqTranscriptomeBinding
from encode_pipeline.adapters.bulk_rnaseq.deployment import (
    TRANSCRIPTOME_BINDING_MANIFEST_ENV,
)
from encode_pipeline.platform.adapters import WorkflowInputs
from bulk_product_runtime import (
    BulkProductBrowserRuntime,
    prepare_bulk_product_browser_runtime,
)
from bulk_rnaseq_real_execution.support import AcceptanceFixture, GateSettings
import platform_runtime


def _fixture(tmp_path: Path) -> AcceptanceFixture:
    transcript = (tmp_path / "fixture/transcripts.fa").resolve()
    transcript.parent.mkdir()
    transcript.write_text(">TX1\nACGT\n", encoding="utf-8")
    rows = [
        {
            "sample": "PE1",
            "library": "libPE",
            "lane": "L001",
            "layout": "PE",
            "fastq_1": str((tmp_path / "fixture/PE1_R1.fastq.gz").resolve()),
            "fastq_2": str((tmp_path / "fixture/PE1_R2.fastq.gz").resolve()),
            "strandedness": "auto",
            "platform": "ILLUMINA",
        },
        {
            "sample": "SE1",
            "library": "libSE",
            "lane": "L001",
            "layout": "SE",
            "fastq_1": str((tmp_path / "fixture/SE1.fastq.gz").resolve()),
            "strandedness": "forward",
            "platform": "ILLUMINA",
        },
    ]
    return AcceptanceFixture(
        workflow_inputs=WorkflowInputs(
            config={"standard": {"analysis": {"alignment": "star"}}},
            samples=rows,
            options={"strict": True},
        ),
        transcriptome=BulkRnaSeqTranscriptomeBinding(
            reference_id="tiny",
            fasta_sha256="a" * 64,
            gtf_sha256="b" * 64,
            transcript_fasta=transcript,
            transcript_fasta_sha256="c" * 64,
        ),
        acceptance_manifest_sha256="d" * 64,
        source_manifest_sha256="e" * 64,
        source_identity_sha256="f" * 64,
        index_provenance_manifest_sha256="1" * 64,
        index_provenance_identity_sha256="2" * 64,
        required_artifact_output_types=("bulk_rnaseq.star.bam",),
        required_qc_metric_keys=("star.input_templates",),
        required_sample_ids=("PE1", "SE1"),
        required_artifact_sample_output_types=(
            ("PE1", "bulk_rnaseq.star.bam"),
            ("SE1", "bulk_rnaseq.star.bam"),
        ),
        required_qc_sample_metric_keys=(
            ("PE1", "star.input_templates"),
            ("SE1", "star.input_templates"),
        ),
        required_qc_sample_metric_values=(
            ("PE1", "star.input_templates", "384"),
            ("SE1", "star.input_templates", "128"),
        ),
    )


def test_product_projection_uses_verified_fixture_and_writes_private_binding(
    tmp_path,
):
    runtime_root = (tmp_path / "browser").resolve()
    runtime_root.mkdir()
    fixture = _fixture(tmp_path)
    fixture_manifest = (tmp_path / "fixture/acceptance.json").resolve()
    fixture_manifest.write_text("{}\n", encoding="utf-8")
    settings = GateSettings(
        runtime_root=(tmp_path / "runtime").resolve(),
        fixture_manifest=fixture_manifest,
        redis_url="redis://127.0.0.1:16379/0",
        docker_executable=Path("/usr/bin/docker"),
        docker_socket=(tmp_path / "docker.sock").resolve(),
    )
    observed_environment = None
    admitted_environment = None

    def load_settings(environ):
        nonlocal observed_environment
        observed_environment = environ
        return settings

    def admit(environ):
        nonlocal admitted_environment
        admitted_environment = dict(environ)

    projected = prepare_bulk_product_browser_runtime(
        runtime_root,
        {"gate": "enabled"},
        settings_loader=load_settings,
        fixture_loader=lambda path: fixture if path == fixture_manifest else None,
        admission_probe=admit,
    )

    assert observed_environment == {"gate": "enabled"}
    fields = projected.manifest_fields
    assert fields["bulkWorkflowId"] == "bulk-rnaseq"
    assert fields["bulkExpectedExecution"] == "available"
    assert fields["bulkConfig"] == fixture.workflow_inputs.config
    assert fields["bulkOptions"] == fixture.workflow_inputs.options
    assert fields["bulkRequiredArtifactOutputTypes"] == ["bulk_rnaseq.star.bam"]
    assert fields["bulkRequiredQcMetricKeys"] == ["star.input_templates"]
    assert fields["bulkRequiredSampleIds"] == ["PE1", "SE1"]
    assert fields["bulkRequiredArtifactSampleOutputTypes"] == [
        ["PE1", "bulk_rnaseq.star.bam"],
        ["SE1", "bulk_rnaseq.star.bam"],
    ]
    assert fields["bulkRequiredQcSampleMetricKeys"] == [
        ["PE1", "star.input_templates"],
        ["SE1", "star.input_templates"],
    ]
    assert fields["bulkRequiredQcSampleMetricValues"] == [
        ["PE1", "star.input_templates", "384"],
        ["SE1", "star.input_templates", "128"],
    ]

    samples_path = Path(fields["bulkSamplesPath"])
    with samples_path.open(encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert [row["layout"] for row in rows] == ["PE", "SE"]
    assert rows[1]["fastq_2"] == ""

    binding_path = Path(
        projected.deployment_environment[TRANSCRIPTOME_BINDING_MANIFEST_ENV]
    )
    assert admitted_environment == {
        "gate": "enabled",
        TRANSCRIPTOME_BINDING_MANIFEST_ENV: str(binding_path),
    }
    binding = json.loads(binding_path.read_text(encoding="utf-8"))
    assert binding == {
        "schema_version": "1.0.0",
        "reference_id": "tiny",
        "fasta_sha256": "a" * 64,
        "gtf_sha256": "b" * 64,
        "transcript_fasta": str(fixture.transcriptome.transcript_fasta),
        "transcript_fasta_sha256": "c" * 64,
    }
    assert TRANSCRIPTOME_BINDING_MANIFEST_ENV not in fields


def test_browser_fixture_selector_requires_explicit_real_gate_and_forwards_environment(
    tmp_path,
    monkeypatch,
):
    runtime_root = tmp_path.resolve()
    observed = None

    def prepare(_runtime_root, environ):
        nonlocal observed
        assert _runtime_root == runtime_root
        observed = environ
        return BulkProductBrowserRuntime(
            manifest_fields={"bulkExpectedExecution": "available"},
            deployment_environment={"PRIVATE_BINDING": "/task/binding.json"},
        )

    monkeypatch.setattr(
        platform_runtime,
        "prepare_bulk_product_browser_runtime",
        prepare,
    )
    source = {
        platform_runtime.PRODUCT_BULK_GATE_ENV: "1",
        "HELIXWEAVE_REQUIRE_BULK_RNASEQ_REAL_EXECUTION": "1",
        "runner_coordinate": "owned",
    }

    fields, deployment = platform_runtime.prepare_bulk_browser_fixture(
        runtime_root,
        source,
    )

    assert observed == source
    assert fields == {"bulkExpectedExecution": "available"}
    assert deployment == {"PRIVATE_BINDING": "/task/binding.json"}

    with pytest.raises(ValueError, match="must be exactly 1"):
        platform_runtime.prepare_bulk_browser_fixture(
            runtime_root,
            {platform_runtime.PRODUCT_BULK_GATE_ENV: "true"},
        )


def test_real_acceptance_flag_alone_does_not_enable_product_browser_execution(
    tmp_path,
):
    runtime_root = tmp_path.resolve()

    fields, deployment = platform_runtime.prepare_bulk_browser_fixture(
        runtime_root,
        {"HELIXWEAVE_REQUIRE_BULK_RNASEQ_REAL_EXECUTION": "1"},
    )

    assert fields["bulkExpectedExecution"] == "not_configured"
    assert deployment == {}
