"""H5 F2 boundary: atomic bundle failure sends nothing over a real SMTP client.

Directed integration beside the full platform chain: original RunService,
SQLite repository, original ArtifactExtraction/QC services and the original
TerminalNotificationService, but with the production ``SmtpEmailTransport``
delivering to a loopback capture sink. This verifies the real SMTP client path
without any external email. It is NOT executed through the real RQ worker: the
F2 worker composition is covered by test/workers/test_atomic_bundle_notifications
.py, and the real-chain negative evidence (no SUCCEEDED mail on failure,
cancel or timeout) lives in test_platform_real.py. No deterministic seam exists
to fail the bundle commit inside a real worker after science success without
modifying the product, so that exact ordering stays at this layer.
"""

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
from encode_pipeline.services.run_repositories import result_bundle_is_complete
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.smtp_transport import SmtpEmailTransport
from encode_pipeline.services.terminal_notifications import TerminalNotificationService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
from encode_pipeline.workers.jobs import _execute_claimed_run

from h5_stack import CaptureSmtp

# Reuse the adapter only: the unrelated helper creates an ENCODE profile and is
# deliberately not called. This fixture fingerprints just a synthetic inventory.
_ADAPTERS = runpy.run_path(
    str(Path(__file__).parents[1] / "services/test_artifact_extraction.py")
)
AtomicAdapter = _ADAPTERS["AtomicQcArtifactAdapter"]

pytestmark = pytest.mark.real_execution


@pytest.fixture(scope="module")
def capture_smtp(tmp_path_factory):
    sink = CaptureSmtp(tmp_path_factory.mktemp("h5-f2-smtp"))
    yield sink
    sink.stop()


@pytest.fixture
def harness(tmp_path, capture_smtp):
    def create():
        adapter = AtomicAdapter(
            (
                ExtractedArtifactCandidate(
                    output_type="summary",
                    relative_path="results/summary.tsv",
                    mime_type="text/tab-separated-values",
                ),
            )
        )
        registry = WorkflowRegistry([adapter])
        url = f"sqlite:///{tmp_path / 'publication.sqlite'}"
        upgrade_database(url)
        engine = create_database_engine(url)
        repository = SqlAlchemyRunRepository(create_session_factory(engine))
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
        settings = SmtpTerminalEmailSettings(
            admin_recipients=("h5-review@example.test",),
            sender="helixweave-h5@example.test",
            application_base_url="http://127.0.0.1",
            smtp_host="127.0.0.1",
            smtp_port=capture_smtp.port,
            tls_mode="local_plaintext",
        )
        notifier = TerminalNotificationService(
            settings=settings,
            run_repository=repository,
            authentication_repository=SimpleNamespace(),
            transport=SmtpEmailTransport(settings),
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
            engine=engine,
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
def test_h5_atomic_failure_no_partial_bundle_no_delivery(
    harness, capture_smtp, monkeypatch, old_bundle, failure
):
    """First-time and post-old-bundle failures publish nothing and send nothing."""
    h = harness()
    try:
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
        assert any(e.event_type == "artifact_extraction_failed" for e in events)
        assert "PRIVATE_SENTINEL" not in repr(events)
        assert _notification_events(h) == []
        # The real SMTP client was never even asked; the sink stays empty.
        assert capture_smtp.captured() == []
    finally:
        h.engine.dispose()


def test_h5_atomic_complete_bundle_delivers_current_qc_once(harness, capture_smtp):
    """The one success notification really transits the SMTP client boundary."""
    h = harness()
    try:
        _execute_claimed_run(h.runtime, "run-1", object())
        state = h.runs.get_result_state("run-1")
        assert result_bundle_is_complete(state, state.artifact_attempt_id)
        captured = capture_smtp.captured()
        assert len(captured) == 1
        message = captured[0]
        assert message["recipients"] == ["rcpt TO:<h5-review@example.test>"]
        assert "Subject: HelixWeave run SUCCEEDED: run-1" in message["data"]
        assert "Run status: SUCCEEDED" in message["data"]
        assert "PETs" in message["data"] and "sample-a" in message["data"]
        assert "original summary bytes" not in message["data"]
        events = _notification_events(h)
        assert len(events) == 1
        assert events[0].event_type == "terminal_email_sent"
        assert events[0].context["metric_count"] == 1
    finally:
        h.engine.dispose()
