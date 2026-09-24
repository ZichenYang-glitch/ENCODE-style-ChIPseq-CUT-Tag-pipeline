"""Atomic result notification policy with original services and capture-only SMTP."""

from __future__ import annotations

from pathlib import Path
import runpy
from types import SimpleNamespace

import pytest

from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.platform.adapters import ExtractedArtifactCandidate, WorkflowInputs
from encode_pipeline.platform.artifact_publications import ArtifactPublicationFilters
from encode_pipeline.platform.notifications import SmtpTerminalEmailSettings
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.artifact_extraction import ArtifactExtractionService
from encode_pipeline.services.qc_summary_indexing import QcSummaryIndexingService
from encode_pipeline.services.run_repositories import (
    InMemoryRunRepository,
    result_bundle_is_complete,
)
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.terminal_notifications import TerminalNotificationService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
from encode_pipeline.workers.jobs import _execute_claimed_run

# Reuse the adapter only: the unrelated helper creates an ENCODE profile and is
# deliberately not called. This fixture fingerprints just a synthetic inventory.
_ADAPTERS = runpy.run_path(
    str(Path(__file__).parents[1] / "services/test_artifact_extraction.py")
)
AtomicAdapter = _ADAPTERS["AtomicQcArtifactAdapter"]


class LegacyAdapter(AtomicAdapter):
    def requires_atomic_result_publication(self):
        return False


@pytest.fixture(params=("memory", "sqlite"))
def repository(request, tmp_path):
    if request.param == "memory":
        yield InMemoryRunRepository()
        return
    url = f"sqlite:///{tmp_path / 'notifications.sqlite'}"
    upgrade_database(url)
    engine = create_database_engine(url)
    try:
        yield SqlAlchemyRunRepository(create_session_factory(engine))
    finally:
        assert engine.pool.checkedout() == 0
        engine.dispose()


@pytest.fixture
def harness(tmp_path, repository):
    def create(*, atomic=True):
        adapter = (AtomicAdapter if atomic else LegacyAdapter)(
            (
                ExtractedArtifactCandidate(
                    output_type="summary",
                    relative_path="results/summary.tsv",
                    mime_type="text/tab-separated-values",
                ),
            )
        )
        registry = WorkflowRegistry([adapter])
        runs = RunService(registry, id_factory=lambda: "run-1", repository=repository)
        project = tmp_path / "project"
        inventory = project / "docs/architecture/artifact-inventory.yaml"
        inventory.parent.mkdir(parents=True)
        inventory.write_bytes(b"artifacts: []\n")
        adapter.identity_root = project
        provider = WorkflowBuildIdentityProvider(registry, project_root=project)
        runs.create_run("fake", WorkflowInputs(config={}))
        runs.transition_run("run-1", RunStatus.VALIDATING)
        runs.complete_preflight("run-1", provider.capture("fake").value)
        runs.transition_run("run-1", RunStatus.QUEUED)
        runs.transition_run("run-1", RunStatus.RUNNING)
        runs.transition_run("run-1", RunStatus.SUCCEEDED)
        workspace = tmp_path / "workspaces/run-1"
        (workspace / "results").mkdir(parents=True)
        source = workspace / "results/summary.tsv"
        source.write_bytes(b"original summary bytes\n")
        extraction = ArtifactExtractionService(
            run_service=runs,
            registry=registry,
            build_identity_provider=provider,
            workspace_root=workspace.parent,
        )
        qc = QcSummaryIndexingService(
            run_service=runs,
            registry=registry,
            build_identity_provider=provider,
            workspace_root=workspace.parent,
        )
        deliveries = []
        notifier = TerminalNotificationService(
            settings=SmtpTerminalEmailSettings(
                admin_recipients=("review@example.test",),
                sender="helixweave@example.test",
                application_base_url="http://localhost:8000",
                smtp_host="127.0.0.1",
                smtp_port=1025,
                tls_mode="local_plaintext",
            ),
            run_repository=repository,
            authentication_repository=SimpleNamespace(),
            transport=SimpleNamespace(
                send=lambda message, recipients: deliveries.append(
                    (message, recipients)
                )
            ),
        )
        runtime = SimpleNamespace(
            registry=registry,
            run_service=runs,
            local_execution_service=SimpleNamespace(
                execute=lambda *_: Result.success(object())
            ),
            artifact_extraction_service=extraction,
            qc_summary_indexing_service=qc,
            terminal_notifier=notifier,
        )
        return SimpleNamespace(
            adapter=adapter,
            runs=runs,
            runtime=runtime,
            extraction=extraction,
            source=source,
            deliveries=deliveries,
        )

    return create


def _public_snapshot(h):
    state = h.runs.get_result_state("run-1")
    return (
        h.runs.list_artifacts("run-1"),
        h.runs.list_qc_metrics("run-1"),
        state.artifact_generation,
        state.qc_generation,
        h.runs._repository.list_artifact_publications(
            filters=ArtifactPublicationFilters(), after=None, limit=100
        ),
    )


def _notification_events(h):
    return [
        event
        for event in h.runs.list_events("run-1")
        if event.event_type.startswith("terminal_email_")
    ]


@pytest.mark.parametrize("old_bundle", [False, True])
@pytest.mark.parametrize("failure", ["qc_prepare", "bundle_commit"])
def test_atomic_failure_keeps_old_results_and_sends_no_success(
    harness, monkeypatch, old_bundle, failure, capsys
):
    h = harness()
    if old_bundle:
        assert h.extraction.extract("run-1").is_success
        h.source.write_bytes(b"new rejected source bytes\n")
    before = _public_snapshot(h)
    scientific_record = h.runs.get_run("run-1")
    if failure == "qc_prepare":
        h.adapter.reject_qc = True
    else:

        def fail_commit(*_args, **_kwargs):
            raise RuntimeError("/private/bundle-commit PRIVATE_SENTINEL")

        monkeypatch.setattr(h.runs, "publish_result_bundle", fail_commit)

    _execute_claimed_run(h.runtime, "run-1", object())

    state = h.runs.get_result_state("run-1")
    assert state.artifact_attempt_status == "failed"
    assert _public_snapshot(h) == before
    assert h.runs.get_run("run-1") == scientific_record
    events = h.runs.list_events("run-1")
    assert any(event.event_type == "artifact_extraction_failed" for event in events)
    assert "PRIVATE_SENTINEL" not in repr(events)
    assert "/private/bundle-commit" not in repr(events)
    observations = [
        {
            "subject": str(message["Subject"]),
            "body": message.get_body(preferencelist=("plain",)).get_content(),
        }
        for message, _ in h.deliveries
    ]
    with capsys.disabled():
        print(
            {
                "old_bundle": old_bundle,
                "failure": failure,
                "notifications": observations,
            }
        )
    assert h.deliveries == []
    assert _notification_events(h) == []


def test_atomic_complete_bundle_sends_current_qc_once(harness):
    h = harness()
    _execute_claimed_run(h.runtime, "run-1", object())
    state = h.runs.get_result_state("run-1")
    assert result_bundle_is_complete(state, state.artifact_attempt_id)
    assert len(h.deliveries) == 1
    message, recipients = h.deliveries[0]
    assert recipients == ("review@example.test",)
    assert str(message["Subject"]) == "HelixWeave run SUCCEEDED: run-1"
    body = message.get_body(preferencelist=("plain",)).get_content()
    assert "Run status: SUCCEEDED" in body
    assert "PETs" in body and "8" in body and "sample-a" in body
    assert "original summary bytes" not in body
    events = _notification_events(h)
    assert len(events) == 1
    assert events[0].event_type == "terminal_email_sent"
    assert events[0].context["metric_count"] == 1


def test_atomic_incomplete_success_cannot_fall_back_to_separate_qc(
    harness, monkeypatch
):
    h = harness()

    def incomplete_commit(run_id, artifacts, _metrics, *, attempt_id, **_kwargs):
        # Fault at the exact commit seam: artifact-only success is not a bundle.
        h.runs.replace_artifacts(run_id, artifacts, attempt_id=attempt_id)

    monkeypatch.setattr(h.runs, "publish_result_bundle", incomplete_commit)
    _execute_claimed_run(h.runtime, "run-1", object())
    state = h.runs.get_result_state("run-1")
    assert h.deliveries == []
    assert state.qc_attempt_id is None
    assert not result_bundle_is_complete(state, state.artifact_attempt_id)


@pytest.mark.parametrize("failure", ["artifact", "qc", None])
def test_non_opt_in_preserves_existing_notification_policy(harness, failure):
    h = harness(atomic=False)
    if failure == "artifact":
        h.adapter.failure = True
    elif failure == "qc":
        h.adapter.reject_qc = True
    _execute_claimed_run(h.runtime, "run-1", object())
    assert h.runs.get_run("run-1").status is RunStatus.SUCCEEDED
    assert len(h.deliveries) == 1
    message, recipients = h.deliveries[0]
    assert recipients == ("review@example.test",)
    assert str(message["Subject"]) == "HelixWeave run SUCCEEDED: run-1"
    body = message.get_body(preferencelist=("plain",)).get_content()
    if failure is not None:
        assert "No persisted QC metrics are available." in body
    else:
        assert "PETs" in body and "sample-a" in body
    assert len(_notification_events(h)) == 1
