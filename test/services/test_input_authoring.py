"""Durable inline-input reconstruction and workspace authoring tests."""

from __future__ import annotations

import csv
import io

import yaml

from encode_pipeline.persistence import open_run_persistence
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.services.defaults import (
    create_default_run_service,
    create_default_workflow_registry,
)
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


def test_inline_snapshot_survives_sqlite_reopen_and_materializes_canonical_files(
    tmp_path,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    row = {
        "sample": "S1",
        "fastq_1": str((tmp_path / "S1.R1.fastq.gz").resolve()),
        "fastq_2": str((tmp_path / "S1.R2.fastq.gz").resolve()),
        "layout": "PE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
        "genome": "hs",
        "bowtie2_index": str((tmp_path / "indices/hs").resolve()),
    }
    registry = create_default_workflow_registry()
    first = open_run_persistence(database_url)
    try:
        run_service = create_default_run_service(
            registry=registry,
            repository=first.repository,
        )
        created = run_service.create_run(
            WORKFLOW_ID,
            WorkflowInputs(config={}, samples=[row], options={}),
        )
    finally:
        first.close()

    second = open_run_persistence(database_url)
    try:
        run_service = create_default_run_service(
            registry=registry,
            repository=second.repository,
        )
        restored = run_service.get_run(created.run_id)
        assert restored.inputs == {
            "config": {},
            "samples": [row],
            "options": {},
        }
        assert "encode-platform-inline-samples" not in str(restored.inputs)

        base_plan = ExecutionPlanner(run_service).plan_run(created.run_id)
        assert base_plan.is_success
        workspace = (tmp_path / "workspace").resolve()
        planned = WorkspacePlanner(registry).plan_workspace(
            base_plan.value,
            base_dir=workspace,
        )
        assert planned.is_success
        materialized = WorkspaceMaterializer().materialize(
            planned.value.workspace_plan,
            workspace,
        )
        assert materialized.is_success
    finally:
        second.close()

    config = yaml.safe_load(
        (workspace / "config/config.yaml").read_text(encoding="utf-8")
    )
    samples_text = (workspace / "config/samples.tsv").read_text(encoding="utf-8")
    samples = list(csv.DictReader(io.StringIO(samples_text), delimiter="\t"))
    assert config["samples"] == "config/samples.tsv"
    assert config["outdir"] == "results"
    assert samples[0]["sample"] == "S1"
    assert samples[0]["fastq_1"] == row["fastq_1"]
    assert "encode-platform-inline-samples" not in config["samples"]
