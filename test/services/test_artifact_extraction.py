from __future__ import annotations

import os
from pathlib import Path
from types import SimpleNamespace

import pytest

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    ExtractedArtifactCandidate,
    MAX_SAMPLE_ROWS,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.result_generations import validate_result_attempt_id
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.artifact_extraction import ArtifactExtractionService
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
from encode_pipeline.workers.timeouts import WorkerHardTimeout


class ArtifactAdapter:
    metadata = WorkflowMetadata(workflow_id="fake", name="Fake", version="1.0.0")
    capabilities = WorkflowCapabilities(supports=("artifact_extract",))

    def __init__(self, candidates=()):
        self.candidates = tuple(candidates)
        self.calls = 0
        self.failure = False

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
        self.calls += 1
        if self.failure:
            return Result.failure(
                [
                    Issue(
                        code="PRIVATE_FAILURE",
                        message=str(workspace),
                        context={"path": str(workspace)},
                    )
                ]
            )
        return Result.success(self.candidates)


class QcArtifactAdapter(ArtifactAdapter):
    capabilities = WorkflowCapabilities(
        supports=("artifact_extract", "qc_summary_extract")
    )

    def qc_source_output_types(self):
        return ("summary",)

    def extract_qc_metrics(self, inputs, sources):
        return Result.success(())


def _assert_failure_context(event, reason_code: str) -> None:
    assert event.context["reason_code"] == reason_code
    validate_result_attempt_id(event.context["attempt_id"])
    assert set(event.context) == {"attempt_id", "reason_code"}


def _project(root: Path) -> Path:
    files = {
        "pyproject.toml": "[project]\nname='fake'\n",
        "docs/architecture/artifact-inventory.yaml": "artifacts: []\n",
        "src/encode_pipeline/runtime.py": "VALUE = 1\n",
        "workflow/Snakefile": "rule all:\n    input: []\n",
        "profiles/default/config.yaml": "cores: 1\n",
        "scripts/tool.py": "VALUE = 1\n",
    }
    for relative, contents in files.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(contents, encoding="utf-8")
    return root


def _service(tmp_path, adapter, *, terminal=RunStatus.SUCCEEDED):
    registry = WorkflowRegistry([adapter])
    run_service = RunService(registry, id_factory=lambda: "run-1")
    provider = WorkflowBuildIdentityProvider(
        registry, project_root=_project(tmp_path / "project")
    )
    run_service.create_run("fake", WorkflowInputs(config={}))
    run_service.transition_run("run-1", RunStatus.VALIDATING)
    identity = provider.capture("fake").value
    run_service.complete_preflight("run-1", identity)
    run_service.transition_run("run-1", RunStatus.QUEUED)
    if terminal is not RunStatus.CANCELLED:
        run_service.transition_run("run-1", RunStatus.RUNNING)
    run_service.transition_run("run-1", terminal)
    workspace_root = tmp_path / "workspaces"
    (workspace_root / "run-1/results").mkdir(parents=True)
    service = ArtifactExtractionService(
        run_service=run_service,
        registry=registry,
        build_identity_provider=provider,
        workspace_root=workspace_root,
    )
    return service, run_service, workspace_root / "run-1", provider


def test_extract_builds_deterministic_opaque_refs_and_is_idempotent(tmp_path):
    candidate = ExtractedArtifactCandidate(
        output_type="summary",
        relative_path="results/summary.tsv",
        mime_type="text/tab-separated-values",
        metadata={"catalog_id": "summary", "scope": "project"},
    )
    adapter = ArtifactAdapter((candidate,))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    (workspace / "results/summary.tsv").write_text("x\n", encoding="utf-8")

    first = service.extract("run-1")
    second = service.extract("run-1")

    assert first.is_success and second.is_success
    assert first.value == second.value == run_service.list_artifacts("run-1")
    artifact = first.value[0]
    assert len(artifact.artifact_id) == 73
    assert artifact.uri == f"run://runs/run-1/artifacts/{artifact.artifact_id}"
    assert artifact.metadata["relative_path"] == "results/summary.tsv"
    assert artifact.metadata["size_bytes"] == 2
    assert str(workspace) not in str(artifact.to_dict())
    assert [event.event_type for event in run_service.list_events("run-1")].count(
        "artifacts_indexed"
    ) == 1


def test_extract_same_length_qc_source_replacement_advances_real_generation(tmp_path):
    candidate = ExtractedArtifactCandidate(
        output_type="summary",
        relative_path="results/summary.tsv",
        mime_type="text/tab-separated-values",
    )
    adapter = QcArtifactAdapter((candidate,))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    source = workspace / "results/summary.tsv"
    source.write_bytes(b"100\n")

    first = service.extract("run-1")
    source.write_bytes(b"900\n")
    second = service.extract("run-1")

    assert first.is_success and second.is_success
    assert first.value[0].metadata["size_bytes"] == 4
    assert second.value[0].metadata["size_bytes"] == 4
    assert first.value[0].artifact_id == second.value[0].artifact_id
    assert first.value[0].revision != second.value[0].revision
    generations = [
        event.context["artifact_generation"]
        for event in run_service.list_events("run-1")
        if event.event_type == "artifacts_indexed"
    ]
    assert len(generations) == 2
    assert generations[0] != generations[1]


def test_hard_timeout_after_artifact_commit_cannot_append_a_failure(
    tmp_path, monkeypatch
):
    candidate = ExtractedArtifactCandidate(
        output_type="summary",
        relative_path="results/summary.tsv",
    )
    adapter = ArtifactAdapter((candidate,))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    (workspace / "results/summary.tsv").write_bytes(b"complete\n")
    attempt_id = service.begin_attempt("run-1")
    replace_artifacts = run_service.replace_artifacts

    def commit_then_time_out(*args, **kwargs):
        replace_artifacts(*args, **kwargs)
        raise WorkerHardTimeout("deadline after durable artifact commit")

    monkeypatch.setattr(run_service, "replace_artifacts", commit_then_time_out)

    with pytest.raises(WorkerHardTimeout):
        service.extract("run-1", attempt_id=attempt_id)
    service.record_unexpected_failure("run-1", attempt_id=attempt_id)

    assert run_service.list_artifacts("run-1")
    outcome_types = [
        event.event_type
        for event in run_service.list_events("run-1")
        if event.event_type in {"artifacts_indexed", "artifact_extraction_failed"}
    ]
    assert outcome_types == ["artifacts_indexed"]


def test_candidate_limit_covers_the_platform_authoring_row_bound():
    import encode_pipeline.services.artifact_extraction as module

    assert module._MAX_ARTIFACT_CANDIDATES >= MAX_SAMPLE_ROWS * 128


@pytest.mark.parametrize(
    "relative_path",
    ["", "/tmp/escape", "results/../escape", "results\\escape", "other/file"],
)
def test_extract_rejects_unsafe_candidate_paths(tmp_path, relative_path):
    adapter = ArtifactAdapter(
        (ExtractedArtifactCandidate(output_type="bad", relative_path=relative_path),)
    )
    service, run_service, _workspace, _provider = _service(tmp_path, adapter)

    result = service.extract("run-1")

    assert result.is_failure
    assert run_service.list_artifacts("run-1") == ()
    event = run_service.list_events("run-1")[-1]
    assert event.event_type == "artifact_extraction_failed"
    _assert_failure_context(event, "ARTIFACT_EXTRACTION_VALIDATION_FAILED")
    assert str(tmp_path) not in str(event.to_dict())


def test_extract_rejects_directory_and_symlink_without_partial_replace(tmp_path):
    candidates = (
        ExtractedArtifactCandidate(
            output_type="good", relative_path="results/good.txt"
        ),
        ExtractedArtifactCandidate(output_type="bad", relative_path="results/link.txt"),
    )
    adapter = ArtifactAdapter(candidates)
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    (workspace / "results/good.txt").write_text("good", encoding="utf-8")
    (workspace / "results/target.txt").write_text("target", encoding="utf-8")
    (workspace / "results/link.txt").symlink_to(workspace / "results/target.txt")

    assert service.extract("run-1").is_failure
    assert run_service.list_artifacts("run-1") == ()

    (workspace / "results/link.txt").unlink()
    (workspace / "results/link.txt").mkdir()
    assert service.extract("run-1").is_failure
    assert run_service.list_artifacts("run-1") == ()


def test_failed_retry_preserves_the_previous_complete_artifact_set(tmp_path):
    original = ExtractedArtifactCandidate(
        output_type="summary",
        relative_path="results/summary.tsv",
    )
    adapter = ArtifactAdapter((original,))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    (workspace / "results/summary.tsv").write_text("complete\n", encoding="utf-8")
    first = service.extract("run-1")
    adapter.candidates = (
        original,
        ExtractedArtifactCandidate(
            output_type="duplicate",
            relative_path="results/summary.tsv",
        ),
    )

    retry = service.extract("run-1")

    assert first.is_success
    assert retry.is_failure
    assert run_service.list_artifacts("run-1") == first.value
    event_types = [event.event_type for event in run_service.list_events("run-1")]
    assert event_types.count("artifacts_indexed") == 1
    assert event_types[-1] == "artifact_extraction_failed"

    adapter.candidates = (original,)
    recovered = service.extract("run-1")

    assert recovered.is_success
    assert run_service.list_artifacts("run-1") == first.value
    event_types = [event.event_type for event in run_service.list_events("run-1")]
    assert event_types.count("artifacts_indexed") == 2
    assert event_types[-1] == "artifacts_indexed"


def test_adapter_failure_is_sanitized_and_keeps_succeeded(tmp_path):
    adapter = ArtifactAdapter()
    adapter.failure = True
    service, run_service, _workspace, _provider = _service(tmp_path, adapter)

    result = service.extract("run-1")

    assert result.is_failure
    assert run_service.get_run("run-1").status is RunStatus.SUCCEEDED
    event = run_service.list_events("run-1")[-1]
    _assert_failure_context(event, "ARTIFACT_EXTRACTION_ADAPTER_FAILED")
    assert "PRIVATE_FAILURE" not in str(event.to_dict())
    assert str(tmp_path) not in str(event.to_dict())


def test_unsupported_capability_fails_without_adapter_call(tmp_path):
    adapter = ArtifactAdapter()
    adapter.capabilities = WorkflowCapabilities(supports=())
    service, run_service, _workspace, _provider = _service(tmp_path, adapter)

    result = service.extract("run-1")

    assert result.is_failure
    assert adapter.calls == 0
    _assert_failure_context(
        run_service.list_events("run-1")[-1],
        "ARTIFACT_EXTRACTION_UNSUPPORTED",
    )


def test_build_identity_drift_fails_before_adapter_call(tmp_path):
    adapter = ArtifactAdapter()
    service, run_service, _workspace, provider = _service(tmp_path, adapter)
    (provider.project_root / "docs/architecture/artifact-inventory.yaml").write_text(
        "artifacts:\n- changed\n",
        encoding="utf-8",
    )

    result = service.extract("run-1")

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.get_run("run-1").status is RunStatus.SUCCEEDED
    _assert_failure_context(
        run_service.list_events("run-1")[-1],
        "ARTIFACT_EXTRACTION_BUILD_MISMATCH",
    )


def test_duplicate_paths_and_absolute_metadata_fail_closed(tmp_path):
    duplicate = ExtractedArtifactCandidate(
        output_type="one", relative_path="results/shared.txt"
    )
    adapter = ArtifactAdapter((duplicate, duplicate))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    (workspace / "results/shared.txt").write_text("x", encoding="utf-8")
    assert service.extract("run-1").is_failure
    assert run_service.list_artifacts("run-1") == ()

    adapter.candidates = (
        ExtractedArtifactCandidate(
            output_type="one",
            relative_path="results/shared.txt",
            metadata={"leak": str(workspace)},
        ),
    )
    assert service.extract("run-1").is_failure
    assert run_service.list_artifacts("run-1") == ()


@pytest.mark.parametrize(
    ("output_type", "mime_type"),
    [
        ("/private/workspace/output", "application/octet-stream"),
        (
            "C:" + chr(92) + "private" + chr(92) + "output",
            "application/octet-stream",
        ),
        ("summary\nprivate", "application/octet-stream"),
        ("summary", "file:///private/workspace/output"),
        ("summary", "text/plain\nprivate"),
    ],
)
def test_hostile_logical_fields_never_reach_persistence(
    tmp_path,
    output_type,
    mime_type,
):
    adapter = ArtifactAdapter(
        (
            ExtractedArtifactCandidate(
                output_type=output_type,
                relative_path="results/output.dat",
                mime_type=mime_type,
            ),
        )
    )
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    (workspace / "results/output.dat").write_bytes(b"output")

    result = service.extract("run-1")

    assert result.is_failure
    assert run_service.list_artifacts("run-1") == ()
    event = run_service.list_events("run-1")[-1]
    assert event.event_type == "artifact_extraction_failed"
    _assert_failure_context(event, "ARTIFACT_EXTRACTION_VALIDATION_FAILED")
    assert "/private/workspace" not in str(event.to_dict())


def test_fifo_is_rejected_without_reading_contents(tmp_path):
    adapter = ArtifactAdapter(
        (ExtractedArtifactCandidate(output_type="fifo", relative_path="results/fifo"),)
    )
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    fifo = workspace / "results/fifo"
    fifo.parent.mkdir(exist_ok=True)
    os.mkfifo(fifo)

    assert service.extract("run-1").is_failure
    assert run_service.list_artifacts("run-1") == ()


def test_candidate_collection_is_bounded_before_any_path_is_indexed(
    tmp_path,
    monkeypatch,
):
    import encode_pipeline.services.artifact_extraction as module

    monkeypatch.setattr(module, "_MAX_ARTIFACT_CANDIDATES", 1)
    adapter = ArtifactAdapter(
        (
            ExtractedArtifactCandidate(
                output_type="one",
                relative_path="results/one.txt",
            ),
            ExtractedArtifactCandidate(
                output_type="two",
                relative_path="results/two.txt",
            ),
        )
    )
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    (workspace / "results/one.txt").write_text("one", encoding="utf-8")

    result = service.extract("run-1")

    assert result.is_failure
    assert run_service.list_artifacts("run-1") == ()
    _assert_failure_context(
        run_service.list_events("run-1")[-1],
        "ARTIFACT_EXTRACTION_VALIDATION_FAILED",
    )


def test_artifact_path_component_count_is_bounded(tmp_path, monkeypatch):
    import encode_pipeline.services.artifact_extraction as module

    monkeypatch.setattr(module, "_MAX_ARTIFACT_PATH_COMPONENTS", 2)
    candidate = ExtractedArtifactCandidate(
        output_type="nested",
        relative_path="results/nested/output.txt",
    )
    adapter = ArtifactAdapter((candidate,))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    target = workspace / "results/nested/output.txt"
    target.parent.mkdir()
    target.write_text("output", encoding="utf-8")

    result = service.extract("run-1")

    assert result.is_failure
    assert run_service.list_artifacts("run-1") == ()


def test_single_artifact_size_is_bounded_without_reading_contents(
    tmp_path,
    monkeypatch,
):
    import encode_pipeline.services.artifact_extraction as module

    monkeypatch.setattr(module, "_MAX_ARTIFACT_BYTES", 4)
    candidate = ExtractedArtifactCandidate(
        output_type="oversized",
        relative_path="results/output.bin",
    )
    adapter = ArtifactAdapter((candidate,))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    (workspace / "results/output.bin").write_bytes(b"12345")

    result = service.extract("run-1")

    assert result.is_failure
    assert run_service.list_artifacts("run-1") == ()


def test_descriptor_reopen_detects_artifact_entry_replacement(
    tmp_path,
    monkeypatch,
):
    import encode_pipeline.services.artifact_extraction as module

    candidate = ExtractedArtifactCandidate(
        output_type="summary",
        relative_path="results/summary.tsv",
    )
    adapter = ArtifactAdapter((candidate,))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    source = workspace / "results/summary.tsv"
    source.write_bytes(b"first")
    displaced = workspace / "results/displaced.tsv"
    real_open = module.os.open
    final_open_count = 0

    def replacing_open(path, flags, *args, **kwargs):
        nonlocal final_open_count
        descriptor = real_open(path, flags, *args, **kwargs)
        if path == "summary.tsv":
            final_open_count += 1
            if final_open_count == 1:
                source.rename(displaced)
                source.write_bytes(b"other")
        return descriptor

    monkeypatch.setattr(module.os, "open", replacing_open)

    result = service.extract("run-1")

    assert result.is_failure
    assert final_open_count >= 2
    assert run_service.list_artifacts("run-1") == ()
    event = run_service.list_events("run-1")[-1]
    _assert_failure_context(event, "ARTIFACT_EXTRACTION_VALIDATION_FAILED")
    assert str(workspace) not in str(event.to_dict())


def test_symlinked_results_component_is_rejected(tmp_path):
    candidate = ExtractedArtifactCandidate(
        output_type="summary",
        relative_path="results/summary.tsv",
    )
    adapter = ArtifactAdapter((candidate,))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    outside = tmp_path / "outside-results"
    outside.mkdir()
    (outside / "summary.tsv").write_bytes(b"outside")
    (workspace / "results").rmdir()
    (workspace / "results").symlink_to(outside, target_is_directory=True)

    result = service.extract("run-1")

    assert result.is_failure
    assert run_service.list_artifacts("run-1") == ()


def test_device_mode_artifact_descriptor_is_rejected(tmp_path, monkeypatch):
    import encode_pipeline.services.artifact_extraction as module

    candidate = ExtractedArtifactCandidate(
        output_type="device",
        relative_path="results/output.dat",
    )
    adapter = ArtifactAdapter((candidate,))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    (workspace / "results/output.dat").write_bytes(b"safe")
    real_open = module.os.open
    real_fstat = module.os.fstat
    final_descriptor = None

    def tracking_open(path, flags, *args, **kwargs):
        nonlocal final_descriptor
        descriptor = real_open(path, flags, *args, **kwargs)
        if path == "output.dat" and final_descriptor is None:
            final_descriptor = descriptor
        return descriptor

    def device_fstat(descriptor):
        if descriptor == final_descriptor:
            return SimpleNamespace(st_mode=module.stat.S_IFCHR, st_size=0)
        return real_fstat(descriptor)

    monkeypatch.setattr(module.os, "open", tracking_open)
    monkeypatch.setattr(module.os, "fstat", device_fstat)

    result = service.extract("run-1")

    assert result.is_failure
    assert final_descriptor is not None
    assert run_service.list_artifacts("run-1") == ()


def test_symlinked_workspace_root_is_rejected(tmp_path):
    candidate = ExtractedArtifactCandidate(
        output_type="summary",
        relative_path="results/summary.tsv",
    )
    adapter = ArtifactAdapter((candidate,))
    service, run_service, workspace, _provider = _service(tmp_path, adapter)
    workspace_root = workspace.parent
    real_workspace_root = tmp_path / "real-workspaces"
    workspace_root.rename(real_workspace_root)
    workspace_root.symlink_to(real_workspace_root, target_is_directory=True)
    (workspace / "results/summary.tsv").write_text("summary\n", encoding="utf-8")

    result = service.extract("run-1")

    assert result.is_failure
    assert run_service.list_artifacts("run-1") == ()
    event = run_service.list_events("run-1")[-1]
    _assert_failure_context(event, "ARTIFACT_EXTRACTION_VALIDATION_FAILED")


@pytest.mark.parametrize("terminal", [RunStatus.FAILED, RunStatus.CANCELLED])
def test_non_succeeded_run_never_calls_adapter(tmp_path, terminal):
    adapter = ArtifactAdapter()
    service, run_service, _workspace, _provider = _service(
        tmp_path, adapter, terminal=terminal
    )

    result = service.extract("run-1")

    assert result.is_failure
    assert adapter.calls == 0
    assert run_service.list_artifacts("run-1") == ()
