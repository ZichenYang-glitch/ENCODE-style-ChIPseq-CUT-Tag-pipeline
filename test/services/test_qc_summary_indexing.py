"""Safety and lifecycle tests for post-success QC summary indexing."""

from __future__ import annotations

from decimal import Decimal
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    ExtractedQcMetricCandidate,
    MAX_SAMPLE_ROWS,
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


class NoQcAdapter:
    metadata = WorkflowMetadata(workflow_id="fake", name="Fake", version="1.0.0")
    capabilities = WorkflowCapabilities(supports=())

    def __init__(self) -> None:
        self.calls = 0

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
    run_id="run-1",
    size_bytes=8,
    sample_id="S1",
    experiment_id="EXP1",
):
    from datetime import datetime, timezone

    return RunArtifactRef(
        artifact_id=artifact_id,
        run_id=run_id,
        artifact_type="file",
        name=Path(relative_path).name,
        uri=f"run://runs/{run_id}/artifacts/{artifact_id}",
        mime_type="text/tab-separated-values",
        produced_at=datetime(2026, 7, 12, tzinfo=timezone.utc),
        metadata={
            "output_type": output_type,
            "relative_path": relative_path,
            "size_bytes": size_bytes,
            "scope": "sample",
            "sample_id": sample_id,
            "experiment_id": experiment_id,
            "assay": "chipseq",
        },
    )


def _service(
    tmp_path,
    adapter,
    *,
    status=RunStatus.SUCCEEDED,
    artifacts=None,
    run_id="run-1",
):
    registry = WorkflowRegistry([adapter])
    run_service = RunService(registry, id_factory=lambda: run_id)
    provider = WorkflowBuildIdentityProvider(
        registry,
        project_root=_project(tmp_path / "project"),
    )
    run_service.create_run("fake", WorkflowInputs(config={}))
    run_service.transition_run(run_id, RunStatus.VALIDATING)
    identity = provider.capture("fake").value
    run_service.complete_preflight(run_id, identity)
    if status in {RunStatus.SUCCEEDED, RunStatus.FAILED}:
        run_service.transition_run(run_id, RunStatus.QUEUED)
        run_service.transition_run(run_id, RunStatus.RUNNING)
    run_service.transition_run(run_id, status)
    workspace_root = (tmp_path / "workspaces").resolve()
    workspace = workspace_root / run_id
    (workspace / "results").mkdir(parents=True)
    artifact_set = () if artifacts is None else tuple(artifacts)
    if status is RunStatus.SUCCEEDED:
        run_service.replace_artifacts(run_id, artifact_set)
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


def test_source_count_limit_covers_the_platform_authoring_row_bound():
    import encode_pipeline.services.qc_summary_indexing as module

    assert module._MAX_SOURCE_FILES >= MAX_SAMPLE_ROWS * 14 + 4
    assert module._MAX_QC_METRICS >= MAX_SAMPLE_ROWS * 100


def test_index_accepts_and_persists_a_workflow_neutral_score_unit(tmp_path):
    candidate = _candidate(
        metric_key="rseqc.tin.mean_score",
        display_name="RSeQC mean transcript integrity number",
        value=Decimal("72.125"),
        unit="score",
    )
    adapter = QcAdapter((candidate,))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

    result = service.index("run-1", artifacts)

    assert result.is_success
    assert result.value[0].unit == "score"
    assert result.value[0].value == Decimal("72.125")
    assert run_service.list_qc_metrics("run-1") == result.value


def test_digit_leading_canonical_run_id_indexes_qc(tmp_path):
    run_id = "1f72c3d4-0000-4000-8000-000000000000"
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact(run_id=run_id)
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
        run_id=run_id,
    )
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

    result = service.index(run_id, artifacts)

    assert result.is_success
    assert result.value == run_service.list_qc_metrics(run_id)
    assert len(result.value) == 1


@pytest.mark.parametrize("sample_id", ("-S1", ".S1", "_S1"))
def test_index_accepts_canonical_leading_punctuation_identifiers(
    tmp_path,
    sample_id,
):
    adapter = QcAdapter((_candidate(sample_id=sample_id),))
    artifact = _artifact(sample_id=sample_id)
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

    result = service.index("run-1", artifacts)

    assert result.is_success
    assert result.value == run_service.list_qc_metrics("run-1")
    assert result.value[0].sample_id == sample_id


@pytest.mark.parametrize("sample_id", ("/private", "..", "bad\nidentifier"))
def test_index_rejects_unsafe_source_identifiers_before_adapter(tmp_path, sample_id):
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact(sample_id=sample_id)
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()
    assert run_service.list_events("run-1")[-1].context == {
        "reason_code": "QC_INDEXING_SOURCE_VALIDATION_FAILED"
    }


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
    adapter = NoQcAdapter()
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
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert run_service.get_run("run-1").status is RunStatus.SUCCEEDED
    assert run_service.list_qc_metrics("run-1") == ()
    event = run_service.list_events("run-1")[-1]
    assert event.event_type == "qc_metrics_indexing_failed"
    assert event.context == {"reason_code": "QC_INDEXING_ADAPTER_FAILED"}
    assert str(tmp_path) not in str(event.to_dict())
    assert "SECRET" not in str(event.to_dict())


@pytest.mark.parametrize("size_bytes", (7, 9))
def test_source_size_must_match_persisted_artifact_metadata(tmp_path, size_bytes):
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact(size_bytes=size_bytes)
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.get_run("run-1").status is RunStatus.SUCCEEDED
    assert run_service.list_qc_metrics("run-1") == ()
    assert run_service.list_events("run-1")[-1].context == {
        "reason_code": "QC_INDEXING_SOURCE_VALIDATION_FAILED"
    }


@pytest.mark.parametrize(
    "size_bytes",
    (True, "8", -1, 16 * 1024 * 1024 + 1),
)
def test_source_size_metadata_must_be_a_bounded_nonnegative_integer(
    tmp_path,
    size_bytes,
):
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact(size_bytes=size_bytes)
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()


def test_fastqc_sized_machine_readable_source_is_supported(tmp_path):
    content = b"x" * (2 * 1024 * 1024)
    adapter = QcAdapter()
    artifact = _artifact(size_bytes=len(content))
    service, _run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    (workspace / "results/summary.tsv").write_bytes(content)

    result = service.index("run-1", artifacts)

    assert result.is_success
    assert adapter.calls == 1
    assert adapter.documents[0].content == content


def test_qc_source_file_count_is_bounded_before_reading(tmp_path, monkeypatch):
    import encode_pipeline.services.qc_summary_indexing as module

    monkeypatch.setattr(module, "_MAX_SOURCE_FILES", 1)
    first = _artifact(
        artifact_id="artifact-1",
        relative_path="results/one.tsv",
        size_bytes=3,
    )
    second = _artifact(
        artifact_id="artifact-2",
        relative_path="results/two.tsv",
        size_bytes=3,
    )
    adapter = QcAdapter()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(second, first),
    )
    (workspace / "results/one.tsv").write_bytes(b"one")
    (workspace / "results/two.tsv").write_bytes(b"two")

    result = service.index("run-1", tuple(reversed(artifacts)))

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.list_events("run-1")[-1].context == {
        "reason_code": "QC_INDEXING_SOURCE_VALIDATION_FAILED"
    }


def test_qc_source_total_bytes_are_bounded_before_reading(tmp_path, monkeypatch):
    import encode_pipeline.services.qc_summary_indexing as module

    monkeypatch.setattr(module, "_MAX_TOTAL_SOURCE_BYTES", 5)
    first = _artifact(
        artifact_id="artifact-1",
        relative_path="results/one.tsv",
        size_bytes=3,
    )
    second = _artifact(
        artifact_id="artifact-2",
        relative_path="results/two.tsv",
        size_bytes=3,
    )
    adapter = QcAdapter()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(first, second),
    )
    (workspace / "results/one.tsv").write_bytes(b"one")
    (workspace / "results/two.tsv").write_bytes(b"two")

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()


def test_multiple_qc_sources_are_delivered_in_deterministic_order(tmp_path):
    first = _artifact(
        artifact_id="artifact-a",
        relative_path="results/a.tsv",
        size_bytes=1,
    )
    second = _artifact(
        artifact_id="artifact-b",
        relative_path="results/b.tsv",
        size_bytes=1,
    )
    adapter = QcAdapter()
    service, _run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(second, first),
    )
    (workspace / "results/a.tsv").write_bytes(b"a")
    (workspace / "results/b.tsv").write_bytes(b"b")

    first_result = service.index("run-1", artifacts)
    first_documents = adapter.documents
    second_result = service.index("run-1", tuple(reversed(artifacts)))

    assert first_result.is_success and second_result.is_success
    assert first_documents == adapter.documents
    assert [document.source.artifact_id for document in adapter.documents] == [
        "artifact-a",
        "artifact-b",
    ]


def test_qc_source_path_component_count_is_bounded(tmp_path, monkeypatch):
    import encode_pipeline.services.qc_summary_indexing as module

    monkeypatch.setattr(module, "_MAX_SOURCE_PATH_COMPONENTS", 2)
    artifact = _artifact(
        relative_path="results/nested/summary.tsv",
        size_bytes=1,
    )
    adapter = QcAdapter()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    target = workspace / "results/nested/summary.tsv"
    target.parent.mkdir()
    target.write_bytes(b"x")

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()


def test_qc_metric_candidate_count_is_bounded_before_persistence(
    tmp_path,
    monkeypatch,
):
    import encode_pipeline.services.qc_summary_indexing as module

    monkeypatch.setattr(module, "_MAX_QC_METRICS", 1)
    adapter = QcAdapter(
        (
            _candidate(),
            _candidate(metric_key="peaks.nsc"),
        )
    )
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 1
    assert run_service.list_qc_metrics("run-1") == ()
    assert run_service.list_events("run-1")[-1].context == {
        "reason_code": "QC_INDEXING_VALIDATION_FAILED"
    }


def test_size_mismatch_preserves_previous_complete_qc_index(tmp_path):
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    source = workspace / "results/summary.tsv"
    source.write_bytes(b"safe-qc\n")
    first = service.index("run-1", artifacts)
    assert first.is_success
    prior_metrics = run_service.list_qc_metrics("run-1")

    source.write_bytes(b"safe-qc\nextra")
    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert adapter.calls == 1
    assert run_service.get_run("run-1").status is RunStatus.SUCCEEDED
    assert run_service.list_qc_metrics("run-1") == prior_metrics
    assert run_service.list_events("run-1")[-1].context == {
        "reason_code": "QC_INDEXING_SOURCE_VALIDATION_FAILED"
    }


def test_source_change_during_read_fails_before_adapter(tmp_path, monkeypatch):
    import encode_pipeline.services.qc_summary_indexing as module

    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    source = workspace / "results/summary.tsv"
    source.write_bytes(b"safe-qc\n")
    real_read = module.os.read
    changed = False

    def changing_read(descriptor, size):
        nonlocal changed
        content = real_read(descriptor, size)
        if content and not changed:
            changed = True
            with source.open("ab") as handle:
                handle.write(b"changed")
        return content

    monkeypatch.setattr(module.os, "read", changing_read)

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert changed is True
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()


def test_source_entry_replacement_during_read_fails_before_adapter(
    tmp_path,
    monkeypatch,
):
    import encode_pipeline.services.qc_summary_indexing as module

    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    source = workspace / "results/summary.tsv"
    source.write_bytes(b"safe-qc\n")
    displaced = workspace / "results/displaced.tsv"
    real_read = module.os.read
    replaced = False

    def replacing_read(descriptor, size):
        nonlocal replaced
        content = real_read(descriptor, size)
        if content and not replaced:
            replaced = True
            source.rename(displaced)
            source.write_bytes(b"other-qc")
        return content

    monkeypatch.setattr(module.os, "read", replacing_read)

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert replaced is True
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()
    event = run_service.list_events("run-1")[-1]
    assert event.context == {"reason_code": "QC_INDEXING_SOURCE_VALIDATION_FAILED"}
    assert str(workspace) not in str(event.to_dict())


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


def test_intermediate_component_swap_to_symlink_fails_closed(tmp_path, monkeypatch):
    import encode_pipeline.services.qc_summary_indexing as module

    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    source = workspace / "results/summary.tsv"
    source.write_bytes(b"inside-safe")
    outside = tmp_path / "outside-results"
    outside.mkdir()
    (outside / "summary.tsv").write_bytes(b"OUTSIDE-SENTINEL")
    original_results = workspace / "results-original"
    real_open = module.os.open
    swapped = False

    def racing_open(path, flags, *args, **kwargs):
        nonlocal swapped
        if path == "results" and not swapped:
            swapped = True
            (workspace / "results").rename(original_results)
            (workspace / "results").symlink_to(outside, target_is_directory=True)
        return real_open(path, flags, *args, **kwargs)

    monkeypatch.setattr(module.os, "open", racing_open)

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert swapped is True
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()


def test_device_mode_descriptor_is_rejected_before_read(tmp_path, monkeypatch):
    import encode_pipeline.services.qc_summary_indexing as module

    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, run_service, workspace, _provider, artifacts = _service(
        tmp_path,
        adapter,
        artifacts=(artifact,),
    )
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc")
    real_open = module.os.open
    real_fstat = module.os.fstat
    final_descriptor = None

    def tracking_open(path, flags, *args, **kwargs):
        nonlocal final_descriptor
        descriptor = real_open(path, flags, *args, **kwargs)
        if path == "summary.tsv":
            final_descriptor = descriptor
        return descriptor

    def device_fstat(descriptor):
        if descriptor == final_descriptor:
            return SimpleNamespace(st_mode=module.stat.S_IFCHR, st_size=0)
        return real_fstat(descriptor)

    monkeypatch.setattr(module.os, "open", tracking_open)
    monkeypatch.setattr(module.os, "fstat", device_fstat)

    result = service.index("run-1", artifacts)

    assert result.is_failure
    assert final_descriptor is not None
    assert adapter.calls == 0
    assert run_service.list_qc_metrics("run-1") == ()


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
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

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
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

    assert service.index("run-1", artifacts).is_failure
    assert run_service.list_qc_metrics("run-1") == ()


def test_index_fails_closed_when_artifact_generation_or_build_changes(tmp_path):
    adapter = QcAdapter((_candidate(),))
    artifact = _artifact()
    service, run_service, workspace, provider, artifacts = _service(
        tmp_path, adapter, artifacts=(artifact,)
    )
    (workspace / "results/summary.tsv").write_bytes(b"safe-qc\n")

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
