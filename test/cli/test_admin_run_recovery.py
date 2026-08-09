"""Public administrator run-recovery CLI contracts."""

from __future__ import annotations

from contextlib import contextmanager
from dataclasses import dataclass, replace
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sqlite3
from types import SimpleNamespace

import pytest

from encode_pipeline.cli import admin
from encode_pipeline.cli.local_platform import (
    RecoveryDoctorStatus,
    run_recovery_doctor,
)
from encode_pipeline.persistence import (
    SqlAlchemyRunRepository,
    create_database_engine,
    create_session_factory,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.execution import RunExecutionAssignment
from encode_pipeline.platform.run_recovery import (
    ExecutionQueueEvidenceState,
    RunExecutionQueueEvidence,
)
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.persistence.migrations import upgrade_database
from encode_pipeline.services.run_repositories import RunEventDraft
from encode_pipeline.workers import rq_queue as rq_queue_module
from encode_pipeline.workers.settings import (
    MANAGED_DOCKER_EXECUTABLE_ENV,
    MANAGED_DOCKER_SOCKET_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
    WorkerSettings,
)


DATABASE_URL = "sqlite:////tmp/helixweave-run-recovery-test.db"
RUN_ID = "run_11111111111111111111111111111111"
JOB_ID = "run-execution-" + "2" * 64


def _assignment(*, claimed: bool = True) -> RunExecutionAssignment:
    now = datetime(2026, 8, 9, 8, 0, tzinfo=timezone.utc)
    return RunExecutionAssignment(
        run_id=RUN_ID,
        job_id=JOB_ID,
        backend="rq",
        queue_name="encode-pipeline-demo",
        created_at=now,
        dispatched_at=now,
        claimed_at=now if claimed else None,
    )


def _diagnostic(**updates):
    values = {
        "run_id": RUN_ID,
        "workflow_id": "bulk-rnaseq",
        "status": RunStatus.RUNNING,
        "diagnosis_code": "RUN_RECOVERY_CALLBACK_GAP",
        "assignment": _assignment(),
        "queue_evidence": SimpleNamespace(state="finished"),
        "result_indexing": "not_started",
        "cleanup": "not_required",
        "gaps": ("callback",),
        "allowed_actions": ("fail",),
    }
    values.update(updates)
    return SimpleNamespace(**values)


@dataclass(frozen=True)
class _ActionResult:
    action: str
    run_id: str
    previous_status: RunStatus
    status: RunStatus
    reason_code: str
    changed: bool
    assignment: RunExecutionAssignment | None = None


class RecordingRecoveryService:
    def __init__(self, diagnostic=None) -> None:
        self.diagnostic = diagnostic or _diagnostic()
        self.calls: list[tuple[str, object]] = []

    def diagnose(self, run_id: str):
        self.calls.append(("diagnose", run_id))
        return self.diagnostic

    def fail_run(self, run_id, *, expected_status, expected_assignment):
        self.calls.append(("fail_run", (run_id, expected_status, expected_assignment)))
        return _ActionResult(
            action="fail",
            run_id=run_id,
            previous_status=expected_status,
            status=RunStatus.FAILED,
            reason_code="RUN_FAILED_BY_ADMIN_RECOVERY",
            changed=True,
        )

    def requeue_run(self, run_id, *, expected_status, expected_assignment):
        self.calls.append(
            ("requeue_run", (run_id, expected_status, expected_assignment))
        )
        return _ActionResult(
            action="requeue",
            run_id=run_id,
            previous_status=expected_status,
            status=expected_status,
            reason_code="RUN_REQUEUED_BY_ADMIN_RECOVERY",
            changed=True,
            assignment=SimpleNamespace(
                **{
                    **expected_assignment.__dict__,
                    "requeue_requested_at": datetime(
                        2026, 8, 9, 8, 1, tzinfo=timezone.utc
                    ),
                    "requeue_confirmed_at": datetime(
                        2026, 8, 9, 8, 2, tzinfo=timezone.utc
                    ),
                }
            ),
        )


def _factory(service, observed_urls):
    @contextmanager
    def factory(database_url: str, *, read_only: bool, queue_name: str | None = None):
        observed_urls.append((database_url, read_only, queue_name))
        yield service

    return factory


def _invoke(arguments, *, service, observed_urls):
    return admin.main(
        ["--database-url", DATABASE_URL, "run", *arguments],
        run_recovery_factory=_factory(service, observed_urls),
    )


def test_run_diagnose_emits_only_the_fixed_safe_projection(capsys) -> None:
    private_path = "/private/operator/workspaces/run-1"
    private_input = {"samples": private_path}
    diagnostic = _diagnostic(
        record=SimpleNamespace(inputs=private_input),
        private_detail=private_path,
    )
    service = RecordingRecoveryService(diagnostic)
    observed_urls = []

    assert (
        _invoke(
            ["diagnose", RUN_ID, "--queue-name", "encode-pipeline-demo"],
            service=service,
            observed_urls=observed_urls,
        )
        == 0
    )

    captured = capsys.readouterr()
    assert captured.err == ""
    assert captured.out == (
        '{"allowed_actions":["fail"],"assignment":{"backend":"rq",'
        '"cancellation_acknowledged":false,"cancellation_requested":false,'
        '"claimed":true,"dispatched":true,"job_id":"'
        + JOB_ID
        + '","queue_name":"encode-pipeline-demo","requeue_confirmed":false,'
        '"requeue_requested":false},"cleanup":"not_required",'
        '"diagnosis_code":"RUN_RECOVERY_CALLBACK_GAP","gaps":["callback"],'
        '"queue":{"state":"finished"},"result_indexing":"not_started",'
        '"run_id":"' + RUN_ID + '","status":"running","workflow_id":"bulk-rnaseq"}\n'
    )
    assert private_path not in captured.out
    assert observed_urls == [(DATABASE_URL, True, "encode-pipeline-demo")]
    assert service.calls == [("diagnose", RUN_ID)]


@pytest.mark.parametrize("command", ["fail", "requeue"])
def test_run_mutations_require_complete_explicit_cas_identity(command, capsys) -> None:
    service = RecordingRecoveryService()
    observed_urls = []

    with pytest.raises(SystemExit, match="2"):
        _invoke([command, RUN_ID], service=service, observed_urls=observed_urls)

    error = capsys.readouterr().err
    for option in ("--expected-status", "--job-id", "--backend", "--queue-name"):
        assert option in error
    assert observed_urls == []
    assert service.calls == []


def test_run_fail_uses_diagnosed_full_assignment_and_compact_result(capsys) -> None:
    service = RecordingRecoveryService()
    observed_urls = []

    assert (
        _invoke(
            [
                "fail",
                RUN_ID,
                "--expected-status",
                "running",
                "--job-id",
                JOB_ID,
                "--backend",
                "rq",
                "--queue-name",
                "encode-pipeline-demo",
            ],
            service=service,
            observed_urls=observed_urls,
        )
        == 0
    )

    assert capsys.readouterr().out == (
        '{"action":"fail","changed":true,"previous_status":"running",'
        '"reason_code":"RUN_FAILED_BY_ADMIN_RECOVERY","run_id":"'
        + RUN_ID
        + '","status":"failed"}\n'
    )
    assert service.calls == [
        ("diagnose", RUN_ID),
        ("fail_run", (RUN_ID, RunStatus.RUNNING, _assignment())),
    ]
    assert observed_urls == [(DATABASE_URL, False, "encode-pipeline-demo")]


def test_run_requeue_emits_only_safe_confirmation_markers(capsys) -> None:
    service = RecordingRecoveryService(
        _diagnostic(
            status=RunStatus.QUEUED,
            diagnosis_code="RUN_RECOVERY_QUEUE_GAP",
            assignment=_assignment(claimed=False),
            allowed_actions=("fail", "requeue"),
        )
    )
    observed_urls = []

    assert (
        _invoke(
            [
                "requeue",
                RUN_ID,
                "--expected-status",
                "queued",
                "--job-id",
                JOB_ID,
                "--backend",
                "rq",
                "--queue-name",
                "encode-pipeline-demo",
            ],
            service=service,
            observed_urls=observed_urls,
        )
        == 0
    )

    assert capsys.readouterr().out == (
        '{"action":"requeue","changed":true,"job_id":"'
        + JOB_ID
        + '","reason_code":"RUN_REQUEUED_BY_ADMIN_RECOVERY",'
        '"requeue_confirmed":true,"requeue_requested":true,"run_id":"'
        + RUN_ID
        + '","status":"queued"}\n'
    )
    assert observed_urls == [(DATABASE_URL, False, "encode-pipeline-demo")]


def test_run_mutation_identity_mismatch_is_a_stable_safe_conflict(capsys) -> None:
    service = RecordingRecoveryService()
    observed_urls = []

    exit_code = _invoke(
        [
            "fail",
            RUN_ID,
            "--expected-status",
            "running",
            "--job-id",
            "wrong-job",
            "--backend",
            "rq",
            "--queue-name",
            "encode-pipeline-demo",
        ],
        service=service,
        observed_urls=observed_urls,
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert captured.out == ""
    assert json.loads(captured.err) == {
        "error": {
            "code": "RUN_RECOVERY_CONFLICT",
            "message": "Run recovery preconditions did not match.",
        }
    }
    assert service.calls == [("diagnose", RUN_ID)]


def test_run_recovery_unexpected_failure_is_redacted(capsys) -> None:
    private_value = "/private/operator/platform.db"

    @contextmanager
    def failing_factory(
        _database_url: str,
        *,
        read_only: bool,
        queue_name: str | None = None,
    ):
        assert read_only is True
        assert queue_name is None
        raise RuntimeError(private_value)
        yield

    exit_code = admin.main(
        ["--database-url", DATABASE_URL, "run", "diagnose", RUN_ID],
        run_recovery_factory=failing_factory,
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert captured.out == ""
    assert json.loads(captured.err) == {
        "error": {
            "code": "RUN_RECOVERY_INTERNAL_ERROR",
            "message": "Run recovery command failed.",
        }
    }
    assert private_value not in captured.err


def test_real_run_diagnose_is_read_only_for_current_sqlite(
    tmp_path: Path,
    capsys,
) -> None:
    database_path = tmp_path / "platform.db"
    database_url = f"sqlite:///{database_path}"
    upgrade_database(database_url)
    before_digest = hashlib.sha256(database_path.read_bytes()).hexdigest()
    before_mtime_ns = database_path.stat().st_mtime_ns
    before_entries = tuple(sorted(path.name for path in tmp_path.iterdir()))

    exit_code = admin.main(
        ["--database-url", database_url, "run", "diagnose", "missing-run"]
    )

    captured = capsys.readouterr()
    assert exit_code == 1
    assert captured.out == ""
    assert json.loads(captured.err)["error"]["code"] == "RUN_RECOVERY_NOT_FOUND"
    assert hashlib.sha256(database_path.read_bytes()).hexdigest() == before_digest
    assert database_path.stat().st_mtime_ns == before_mtime_ns
    assert tuple(sorted(path.name for path in tmp_path.iterdir())) == before_entries


def test_real_run_diagnose_refuses_prior_schema_without_upgrade(
    tmp_path: Path,
    capsys,
) -> None:
    database_path = tmp_path / "platform.db"
    database_url = f"sqlite:///{database_path}"
    upgrade_database(database_url, revision="20260807_12")
    before_digest = hashlib.sha256(database_path.read_bytes()).hexdigest()

    exit_code = admin.main(
        ["--database-url", database_url, "run", "diagnose", "missing-run"]
    )

    assert exit_code == 1
    assert json.loads(capsys.readouterr().err)["error"]["code"] == (
        "RUN_RECOVERY_DATA_INVALID"
    )
    assert hashlib.sha256(database_path.read_bytes()).hexdigest() == before_digest
    with sqlite3.connect(f"{database_path.as_uri()}?mode=ro", uri=True) as database:
        assert database.execute(
            "SELECT version_num FROM alembic_version"
        ).fetchone() == ("20260807_12",)


def test_real_run_diagnose_missing_database_does_not_create_parent(
    tmp_path: Path,
    capsys,
) -> None:
    database_path = tmp_path / "missing" / "platform.db"

    exit_code = admin.main(
        [
            "--database-url",
            f"sqlite:///{database_path}",
            "run",
            "diagnose",
            "missing-run",
        ]
    )

    assert exit_code == 1
    assert json.loads(capsys.readouterr().err)["error"]["code"] == (
        "RUN_RECOVERY_DATA_INVALID"
    )
    assert not database_path.parent.exists()


def test_real_run_cleanup_gap_refuses_fail_without_managed_cleaner(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    database_path = tmp_path / "private" / "platform.db"
    database_path.parent.mkdir()
    database_url = f"sqlite:///{database_path}"
    private_input = str(tmp_path / "private" / "sample.csv")
    private_workspace = str(tmp_path / "private" / "workspaces")
    private_redis_url = "redis://operator:secret@private-redis.invalid:6379/0"
    now = datetime(2026, 8, 9, 8, 0, tzinfo=timezone.utc)
    validating = RunRecord(
        run_id=RUN_ID,
        workflow_id="bulk-rnaseq",
        inputs={"samples": private_input},
        status=RunStatus.VALIDATING,
        created_at=now,
        updated_at=now,
        started_at=None,
        ended_at=None,
        current_stage="validation",
        cancellation_reason=None,
        error=None,
    )
    upgrade_database(database_url)
    setup_engine = create_database_engine(database_url)
    try:
        repository = SqlAlchemyRunRepository(create_session_factory(setup_engine))
        repository.create_run(
            validating,
            RunEventDraft(
                event_type="run_created",
                message="Run created.",
                status=RunStatus.VALIDATING,
            ),
        )
        planned = replace(
            validating,
            status=RunStatus.PLANNED,
            current_stage="planning",
        )
        repository.complete_preflight(
            planned,
            WorkflowBuildIdentity(
                workflow_id="bulk-rnaseq",
                adapter_version="1.0.0",
                scheme="sha256-tree-v1",
                logical_entrypoint="main.nf",
                digest="a" * 64,
                captured_at=now,
            ),
            expected_status=RunStatus.VALIDATING,
            event=RunEventDraft(
                event_type="preflight_completed",
                message="Preflight completed.",
                status=RunStatus.PLANNED,
            ),
        )
        running = replace(
            planned,
            status=RunStatus.RUNNING,
            started_at=now,
            current_stage="execution",
        )
        repository.update_run(
            running,
            expected_status=RunStatus.PLANNED,
            event=RunEventDraft(
                event_type="status_changed",
                message="Run started.",
                status=RunStatus.RUNNING,
            ),
        )
        repository.ensure_execution_assignment(
            _assignment(),
            expected_status=RunStatus.RUNNING,
        )
    finally:
        setup_engine.dispose()

    class ExactTerminalRunQueue:
        def __init__(self, settings: WorkerSettings) -> None:
            assert settings.queue_name == "encode-pipeline-demo"

        def inspect_execution(
            self,
            _assignment: RunExecutionAssignment,
        ) -> RunExecutionQueueEvidence:
            return RunExecutionQueueEvidence(state=ExecutionQueueEvidenceState.FINISHED)

        def requeue_execution(self, _assignment: RunExecutionAssignment) -> str:
            pytest.fail("cleanup-gap recovery must not enqueue")

        def close(self) -> None:
            return None

    monkeypatch.setattr(rq_queue_module, "RqRunQueue", ExactTerminalRunQueue)
    monkeypatch.delenv(MANAGED_DOCKER_EXECUTABLE_ENV, raising=False)
    monkeypatch.delenv(MANAGED_DOCKER_SOCKET_ENV, raising=False)
    monkeypatch.setenv(WORKSPACE_ROOT_ENV, private_workspace)
    monkeypatch.setenv(REDIS_URL_ENV, private_redis_url)

    def persisted_state():
        engine = create_database_engine(database_url)
        try:
            repository = SqlAlchemyRunRepository(create_session_factory(engine))
            return (
                repository.get_run(RUN_ID),
                repository.get_execution_assignment(RUN_ID),
                repository.list_events(RUN_ID),
                repository.get_result_state(RUN_ID),
            )
        finally:
            engine.dispose()

    before = persisted_state()
    assert (
        admin.main(
            [
                "--database-url",
                database_url,
                "run",
                "diagnose",
                RUN_ID,
                "--queue-name",
                "encode-pipeline-demo",
            ]
        )
        == 0
    )
    diagnosis_output = capsys.readouterr()
    diagnosis = json.loads(diagnosis_output.out)
    assert diagnosis_output.err == ""
    assert diagnosis["queue"] == {"state": "finished"}
    assert diagnosis["cleanup"] == "blocked"
    assert diagnosis["gaps"] == ["callback", "cleanup"]
    assert diagnosis["allowed_actions"] == []

    doctor = run_recovery_doctor(
        database_url=database_url,
        redis_url=private_redis_url,
        queue_name="encode-pipeline-demo",
    )
    assert doctor.status is RecoveryDoctorStatus.ATTENTION_REQUIRED
    assert doctor.reason_code == "RUN_RECOVERY_ATTENTION_REQUIRED"
    assert doctor.counts.runs_examined == 1
    assert doctor.counts.database_gaps == 0
    assert doctor.counts.queue_gaps == 0
    assert doctor.counts.callback_gaps == 1
    assert doctor.counts.result_indexing_gaps == 0
    assert doctor.counts.cleanup_gaps == 1

    exit_code = admin.main(
        [
            "--database-url",
            database_url,
            "run",
            "fail",
            RUN_ID,
            "--expected-status",
            "running",
            "--job-id",
            JOB_ID,
            "--backend",
            "rq",
            "--queue-name",
            "encode-pipeline-demo",
        ]
    )

    failure_output = capsys.readouterr()
    assert exit_code == 1
    assert failure_output.out == ""
    assert json.loads(failure_output.err) == {
        "error": {
            "code": "RUN_RECOVERY_NOT_SAFE",
            "message": "Requested run recovery action is not safe.",
        }
    }
    combined_output = diagnosis_output.out + diagnosis_output.err + failure_output.err
    for private_value in (
        private_input,
        private_workspace,
        private_redis_url,
        "operator:secret",
        "private-redis.invalid",
        str(database_path),
    ):
        assert private_value not in combined_output
    assert persisted_state() == before
    assert all(
        event.event_type != "run_failed_by_admin_recovery"
        for event in persisted_state()[2]
    )
