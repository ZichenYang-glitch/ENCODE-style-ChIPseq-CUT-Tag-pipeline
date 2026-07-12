"""Safety and lifecycle tests for post-success QC summary indexing."""

from __future__ import annotations

from decimal import Decimal
import os
from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    ExtractedQcMetricCandidate,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunArtifactRef, RunStatus
from encode_pipeline.services.qc_summary_indexing import QcSummaryIndexingService
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


class QcAdapter:
    metadata = WorkflowMetadata(workflow_id="fake", name="Fake", version="1.0.0")
    capabilities = WorkflowCapabilities(supports=("qc_summary_extract",))

    def __init__(self, candidates=()):
        self.candidates = tuple(candidates)
        self.calls = 0
        self.documents = ()
        self.failure = False
        self.source_types = ("qc_summary",)

    def schema(self):
        return WorkflowSchema()

    def validate(self, inputs):
        return Result.success({})

    def preview_dag(self, inputs):
        return Result.success(DagPreview())

    def plan_workspace(self, inputs, workspace):
        return Result.success(WorkspacePlan())

    def build_command(self, plan):
        return Result.success(CommandSpec(argv=("fake",)))

    def extract_artifacts(self, inputs, workspace):
        return Result.success(())

    def qc_source_output_types(self):
        return self.source_types

    def extract_qc_metrics(self, inputs, sources):
        self.calls += 1
        self.documents = sources
        if self.failure:
            return Result.failure(
                [
                    Issue(
                        code="PRIVATE_QC_ERROR",
                        message="/private/workspace/secret",
                        technical_message="SECRET=value",
                    )
                ]
            )
        return Result.success(self.candidates)


def _project(root: Path) -> Path:
    files = {
        "pyproject.toml": "[project]\nname='fake'\n",
        "docs/architecture/artifact-inventory.yaml": "artifacts: []\n",
        "src/encode_pipeline/adapters/encode_qc.py": "CATALOG = 1\n",
        "workflow/Snakefile": "rule all:\n    input: []\n",
        "profiles/default/config.yaml": "cores: 1\n",
        "scripts/qc_writer.py": "HEADER = ('metric',)\n",
    }
    for relative, content in files.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content, encoding="utf-8")
    return root


def _candidate(**overrides):
    values = {
        "metric_key": "peaks.frip",
        "display_name": "Fraction of reads in peaks",
        "value": Decimal("0.125"),
        "unit": "fraction",
        "scope": "sample",
        "sample_id": "S1",
        "experiment_id": "EXP1",
        "assay": "chipseq",
        "qc_flag": None,
        "source_artifact_id": "artifact-1",
    }
    values.update(overrides)
    return ExtractedQcMetricCandidate(**values)


def _artifact(
    *,
    artifact_id="artifact-1",
    output_type="qc_summary",
    relative_path="results/summary.tsv",
):
    from datetime import datetime, timezone

    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id="run-1",
        artifact_type="file",
        name=Path(relative_path).name,
        uri=f"run://runs/run-1/artifacts/{artifact_id}",
        mime_type="text/tab-separated-values",
        produced_at=datetime(2026, 7, 12, tzinfo=timezone.utc),
        metadata={
            "output_type": output_type,
            "relative_path": relative_path,
            "size_bytes": 8,
            "scope": "sample",
            "sample_id": "S1",
            "experiment_id": "EXP1",
            "assay": "chipseq",
        },
    )


def _service(tmp_path, adapter, *, status=RunStatus.SUCCEEDED, artifacts=None):
    registry = WorkflowRegistry([adapter])
    run_service = RunService(registry, id_factory=lambda: "run-1")
    provider = WorkflowBuildIdentityProvider(
        registry,
        project_root=_project(tmp_path / "project"),
    )
    run_service.create_run("fake", WorkflowInputs(config={}))
    run_service.transition_run("run-1", RunStatus.VALIDATING)
    identity = provider.capture("fake").value
    run_service.complete_preflight("run-1", identity)
    if status in {RunStatus.SUCCEEDED, RunStatus.FAILED}:
        run_service.transition_run("run-1", RunStatus.QUEUED)
        run_service.transition_run("run-1", RunStatus.RUNNING)
    run_service.transition_run("run-1", status)
    workspace_root = (tmp_path / "workspaces").resolve()
    workspace = workspace_root / "run-1"
    (workspace / "results").mkdir(parents=True)
    artifact_set = () if artifacts is None else tuple(artifacts)
    if status is RunStatus.SUCCEEDED:
        run_service.replace_artifacts("run-1", artifact_set)
    service = QcSummaryIndexingService(
        run_service=run_service,
        registry=registry,
        build_identity_provider=provider,
        workspace_root=workspace_root,
    )
    return service, run_service, workspace, provider, artifact_set


def test_index_builds_deterministic_metrics_and_is_idempotent(tmp_path):
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

    first = service.index("run-1", artifacts)
    second = service.index("run-1", tuple(reversed(artifacts)))

    assert first.is_success and second.is_success
    assert first.value == second.value == run_service.list_qc_metrics("run-1")
    metric = first.value[0]
    assert metric.metric_id.startswith("qcmetric-")
    assert len(metric.metric_id) == 73
    assert metric.source_artifact_id == artifact.artifact_id
    assert str(workspace) not in str(metric.to_dict())
    assert adapter.documents[0].content == b"safe-qc\n"
    assert adapter.documents[0].source.relative_path == "results/summary.tsv"
    assert [event.event_type for event in run_service.list_events("run-1")].count(
        "qc_metrics_indexed"
    ) == 1


def test_no_declared_qc_sources_is_a_confirmed_empty_index(tmp_path):
    adapter = QcAdapter()
    artifact = _artifact(output_type="result_manifest")
    service, run_service, _workspace, _provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )

    result = service.index("run-1", artifacts)

    assert result.is_success and result.value == ()
    assert adapter.calls == 1
    assert adapter.documents == ()
    event = run_service.list_events("run-1")[-1]
    assert event.event_type == "qc_metrics_indexed"
    assert event.context == {"metric_count": 0}


def test_unsupported_adapter_fails_closed_without_reading_sources(tmp_path):
    adapter = QcAdapter((_candidate(),))
    adapter.capabilities = WorkflowCapabilities(supports=())
    artifact = _artifact()
    service, run_service, _workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()
    assert run_service.list_events("run-1")[-1].context == {
        "reason_code": "QC_INDEXING_UNSUPPORTED"
    }


def test_adapter_failure_preserves_succeeded_and_redacts_private_details(tmp_path):
    adapter = QcAdapter()
    adapter.failure = True
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )
    (workspace / "results/summary.tsv").write_text("safe\n", encoding="utf-8")

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert run_service.get_run("run-1").status is RunStatus.SUCCEEDED
    assert run_service.list_qc_metrics("run-1") == ()
    event = run_service.list_events("run-1")[-1]
    assert event.event_type == "qc_metrics_indexing_failed"
    assert event.context == {"reason_code": "QC_INDEXING_ADAPTER_FAILED"}
    assert str(tmp_path) not in str(event.to_dict())
    assert "SECRET" not in str(event.to_dict())


@pytest.mark.parametrize(
    "reason_code",
    ("/private/workspace", "QC FAILED", "qc_indexing_failed", "A" * 129),
)
def test_failure_event_rejects_non_public_reason_codes(tmp_path, reason_code):
    adapter = QcAdapter()
    _service_under_test, run_service, _workspace, _provider, _artifacts = _service(
        tmp_path,
        adapter,
    )
    events_before = run_service.list_events("run-1")

    with pytest.raises(ValueError, match="public-safe"):
        run_service.record_qc_metrics_failure("run-1", reason_code=reason_code)

    assert run_service.list_events("run-1") == events_before


@pytest.mark.parametrize("status", [RunStatus.FAILED, RunStatus.CANCELLED])
def test_non_succeeded_runs_do_not_access_adapter_or_sources(tmp_path, status):
    adapter = QcAdapter((_candidate(),))
    service, run_service, _workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        status=status,
    )

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()
    assert "qc_metrics_indexed" not in {
        event.event_type for event in run_service.list_events("run-1")
    }


@pytest.mark.parametrize(
    "relative_path",
    ["", "/tmp/escape", "results/../escape", "results\\escape", "other/file"],
)
def test_index_rejects_unsafe_persisted_source_paths(tmp_path, relative_path):
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact(relative_path=relative_path)
    service, run_service, _workspace, _provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.list_events("run-1")[-1].context == {
        "reason_code": "QC_INDEXING_SOURCE_VALIDATION_FAILED"
    }


def test_index_rejects_symlink_directory_fifo_and_oversized_sources(tmp_path):
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )
    source = workspace / "results/summary.tsv"

    source.mkdir()
    assert service.index("run-1", artifacts).is_failure
    source.rmdir()
    outside = tmp_path / "outside"
    outside.write_bytes(b"outside-secret")
    source.symlink_to(outside)
    assert service.index("run-1", artifacts).is_failure
    source.unlink()
    os.mkfifo(source)
    assert service.index("run-1", artifacts).is_failure
    source.unlink()
    source.write_bytes(b"x" * (1024 * 1024 + 1))
    assert service.index("run-1", artifacts).is_failure

    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()


def test_descriptor_race_never_passes_outside_bytes_to_adapter(
    tmp_path,
    monkeypatch,
):
    import encode_pipeline.services.qc_summary_indexing as module

    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, _run_service, workspace, _provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )
    source = workspace / "results/summary.tsv"
    source.write_bytes(b"inside-safe")
    outside = tmp_path / "outside-secret.tsv"
    outside.write_bytes(b"OUTSIDE-SENTINEL")
    real_open = module.os.open
    swapped = False

    def racing_open(path, flags, *args, **kwargs):
        nonlocal swapped
        if path == "summary.tsv" and not swapped:
            swapped = True
            source.unlink()
            source.symlink_to(outside)
        return real_open(path, flags, *args, **kwargs)

    monkeypatch.setattr(module.os, "open", racing_open)

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert swapped is True
    assert adapter.calls == 0
    assert all(
        b"OUTSIDE-SENTINEL" not in document.content for document in adapter.documents
    )


def test_index_rejects_symlinked_parent_component(tmp_path):
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    outside = tmp_path / "outside-results"
    outside.mkdir()
    (outside / "summary.tsv").write_bytes(b"OUTSIDE-SENTINEL")
    (workspace / "results").rmdir()
    (workspace / "results").symlink_to(outside, target_is_directory=True)

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()
    assert run_service.list_events("run-1")[-1].context == {
        "reason_code": "QC_INDEXING_SOURCE_VALIDATION_FAILED"
    }


@pytest.mark.parametrize(
    "candidate",
    [
        _candidate(value=Decimal("NaN")),
        _candidate(value=Decimal("Infinity")),
        _candidate(unit="/private/unit"),
        _candidate(metric_key="../escape"),
        _candidate(scope="sample", sample_id=None),
        _candidate(source_artifact_id="artifact-other"),
    ],
)
def test_index_rejects_invalid_adapter_candidates(tmp_path, candidate):
    adapter = QcAdapter((candidate,))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )
    (workspace / "results/summary.tsv").write_text("safe\n", encoding="utf-8")

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert run_service.list_qc_metrics("run-1") == ()
    assert run_service.list_events("run-1")[-1].context == {
        "reason_code": "QC_INDEXING_VALIDATION_FAILED"
    }


def test_index_rejects_duplicate_semantic_metric_without_partial_replace(tmp_path):
    adapter = QcAdapter((_candidate(), _candidate(display_name="Duplicate")))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )
    (workspace / "results/summary.tsv").write_text("safe\n", encoding="utf-8")

    assert service.index("run-1", artifacts).is_failure
    assert run_service.list_qc_metrics("run-1") == ()


def test_index_fails_closed_when_artifact_generation_or_build_changes(tmp_path):
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, run_service, workspace, provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )
    (workspace / "results/summary.tsv").write_text("safe\n", encoding="utf-8")

    stale = service.index("run-1", (_artifact(artifact_id="artifact-stale"),))
    assert stale.is_failure
    assert adapter.calls == 0
    (provider.project_root / "src/encode_pipeline/adapters/encode_qc.py").write_text(
        "CATALOG = 2\n",
        encoding="utf-8",
    )
    changed = service.index("run-1", artifacts)

    assert changed.is_failure
    assert adapter.calls == 0
    assert run_service.get_run("run-1").status is RunStatus.SUCCEEDED
