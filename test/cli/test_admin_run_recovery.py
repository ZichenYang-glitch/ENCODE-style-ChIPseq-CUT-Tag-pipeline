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
from encode_pipeline.services.authentication_service import (
    AuthenticationActor,
)
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
from encode_pipeline.platform.managed_containers import (
    managed_container_endpoint_identity,
    managed_container_scope,
)
from encode_pipeline.platform.run_recovery import (
    ExecutionQueueEvidenceState,
    RunExecutionQueueEvidence,
)
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.persistence.migrations import upgrade_database
from encode_pipeline.services.run_repositories import RunEventDraft
from encode_pipeline.services import managed_containers as managed_containers_module
from encode_pipeline.services import terminal_notifications
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
FAKE_NOTIFICATIONS_ENV = Path("/tmp/helixweave-test-notifications.env")


def test_mutating_admin_recovery_composes_the_terminal_notifier(
    tmp_path,
    monkeypatch,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'composition.db'}"

    class Notifier:
        def notify_terminal_run(self, _run_id, _status, *, include_qc=False):
            del include_qc

    notifier = Notifier()
    observed = []
    notifications_environment = tmp_path / "notifications.env"
    notifications_environment.write_text(
        "\n".join(
            (
                "HELIXWEAVE_TERMINAL_EMAIL_ENABLED=true",
                "HELIXWEAVE_TERMINAL_EMAIL_ADMIN_RECIPIENTS=admin@example.test",
                "HELIXWEAVE_TERMINAL_EMAIL_FROM=helixweave@example.test",
                "HELIXWEAVE_TERMINAL_EMAIL_APPLICATION_BASE_URL=https://example.test",
                "HELIXWEAVE_SMTP_HOST=127.0.0.1",
                "HELIXWEAVE_SMTP_PORT=2525",
                "HELIXWEAVE_SMTP_TLS_MODE=local_plaintext",
                "HELIXWEAVE_SMTP_TIMEOUT_SECONDS=1",
            )
        )
        + "\n",
        encoding="utf-8",
    )
    notifications_environment.chmod(0o600)

    def compose(**kwargs):
        observed.append(kwargs)
        return notifier

    monkeypatch.setattr(
        terminal_notifications,
        "compose_terminal_run_notifier",
        compose,
    )

    with admin._open_run_recovery(
        database_url,
        read_only=False,
        notifications_environment_path=notifications_environment,
    ) as recovery:
        assert recovery._terminal_notifier is notifier

    assert len(observed) == 1
    assert set(observed[0]) == {
        "environ",
        "run_repository",
        "authentication_repository",
    }
    assert observed[0]["environ"]["HELIXWEAVE_TERMINAL_EMAIL_ENABLED"] == "true"
    assert observed[0]["environ"]["HELIXWEAVE_SMTP_HOST"] == "127.0.0.1"


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

    def fail_run(
        self,
        run_id,
        *,
        expected_status,
        expected_assignment,
        security_audit_actor=None,
    ):
        self.calls.append(
            (
                "fail_run",
                (
                    run_id,
                    expected_status,
                    expected_assignment,
                    security_audit_actor,
                ),
            )
        )
        return _ActionResult(
            action="fail",
            run_id=run_id,
            previous_status=expected_status,
            status=RunStatus.FAILED,
            reason_code="RUN_FAILED_BY_ADMIN_RECOVERY",
            changed=True,
        )

    def requeue_run(
        self,
        run_id,
        *,
        expected_status,
        expected_assignment,
        security_audit_actor=None,
    ):
        self.calls.append(
            (
                "requeue_run",
                (
                    run_id,
                    expected_status,
                    expected_assignment,
                    security_audit_actor,
                ),
            )
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
    def factory(
        database_url: str,
        *,
        read_only: bool,
        queue_name: str | None = None,
        notifications_environment_path: Path | None = None,
    ):
        del notifications_environment_path
        observed_urls.append((database_url, read_only, queue_name))
        yield service

    return factory


def _invoke(arguments, *, service, observed_urls):
    return admin.main(
        ["--database-url", DATABASE_URL, "run", *arguments],
        run_recovery_factory=_factory(service, observed_urls),
    )


def _disabled_notifications_environment(tmp_path: Path) -> Path:
    path = (tmp_path / "notifications.env").resolve()
    path.write_text(
        "HELIXWEAVE_TERMINAL_EMAIL_ENABLED=false\n",
        encoding="utf-8",
    )
    path.chmod(0o600)
    return path


def _persist_real_running_run(
    database_url: str,
    *,
    assignment: RunExecutionAssignment,
    private_input: str,
) -> None:
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
    engine = create_database_engine(database_url)
    try:
        repository = SqlAlchemyRunRepository(create_session_factory(engine))
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
            assignment,
            expected_status=RunStatus.RUNNING,
        )
    finally:
        engine.dispose()


class _ExactTerminalRunQueue:
    def __init__(self, settings: WorkerSettings) -> None:
        assert settings.queue_name == "encode-pipeline-demo"

    def inspect_execution(
        self,
        _assignment: RunExecutionAssignment,
    ) -> RunExecutionQueueEvidence:
        return RunExecutionQueueEvidence(state=ExecutionQueueEvidenceState.FINISHED)

    def requeue_execution(self, _assignment: RunExecutionAssignment) -> str:
        pytest.fail("claimed recovery must not enqueue")

    def close(self) -> None:
        return None


class _AlwaysSuccessfulCleaner:
    cleanup_scopes: list[str] = []

    def __init__(self, *, executable: Path, unix_socket: Path) -> None:
        self.endpoint_identity = managed_container_endpoint_identity(
            executable,
            unix_socket,
        )

    def cleanup(self, scope: str):
        type(self).cleanup_scopes.append(scope)
        return SimpleNamespace(is_success=True)


class _UnavailableCleaner:
    def __init__(self, *, executable: Path, unix_socket: Path) -> None:
        raise OSError(f"unavailable endpoint: {executable} {unix_socket}")


def _real_recovery_state(database_url: str):
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
    if command == "fail":
        assert "--notifications-env-file" in error
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
                "--notifications-env-file",
                str(FAKE_NOTIFICATIONS_ENV),
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
        (
            "fail_run",
            (
                RUN_ID,
                RunStatus.RUNNING,
                _assignment(),
                AuthenticationActor.local_operator(),
            ),
        ),
    ]
    assert observed_urls == [(DATABASE_URL, False, "encode-pipeline-demo")]


def test_run_fail_forwards_private_notification_environment_path(
    tmp_path,
    capsys,
) -> None:
    service = RecordingRecoveryService()
    environment_path = (tmp_path / "notifications.env").resolve()
    observed = []

    @contextmanager
    def factory(
        database_url: str,
        *,
        read_only: bool,
        queue_name: str | None = None,
        notifications_environment_path: Path | None = None,
    ):
        observed.append(
            (
                database_url,
                read_only,
                queue_name,
                notifications_environment_path,
            )
        )
        yield service

    exit_code = admin.main(
        [
            "--database-url",
            DATABASE_URL,
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
            "--notifications-env-file",
            str(environment_path),
        ],
        run_recovery_factory=factory,
    )

    assert exit_code == 0
    assert capsys.readouterr().err == ""
    assert observed == [
        (DATABASE_URL, False, "encode-pipeline-demo", environment_path),
    ]


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
            "--notifications-env-file",
            str(FAKE_NOTIFICATIONS_ENV),
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
    notifications_environment = _disabled_notifications_environment(tmp_path)
    private_input = str(tmp_path / "private" / "sample.csv")
    private_workspace = str(tmp_path / "private" / "workspaces")
    private_redis_url = "redis://operator:secret@private-redis.invalid:6379/0"
    _persist_real_running_run(
        database_url,
        assignment=_assignment(),
        private_input=private_input,
    )

    monkeypatch.setattr(rq_queue_module, "RqRunQueue", _ExactTerminalRunQueue)
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
            "--notifications-env-file",
            str(notifications_environment),
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

    requeue_exit_code = admin.main(
        [
            "--database-url",
            database_url,
            "run",
            "requeue",
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
    requeue_output = capsys.readouterr()
    assert requeue_exit_code == 1
    assert requeue_output.out == ""
    assert json.loads(requeue_output.err)["error"]["code"] == ("RUN_RECOVERY_NOT_SAFE")
    assert persisted_state() == before


@pytest.mark.parametrize("endpoint_matches", [True, False])
def test_real_cli_binds_cleanup_to_durable_scope_and_endpoint(
    endpoint_matches: bool,
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    database_path = tmp_path / "private" / "platform.db"
    database_path.parent.mkdir()
    database_url = f"sqlite:///{database_path}"
    notifications_environment = _disabled_notifications_environment(tmp_path)
    original_workspace = (tmp_path / "private" / "original-workspaces").resolve()
    original_executable = (tmp_path / "private" / "original" / "docker").resolve()
    original_socket = (tmp_path / "private" / "original.sock").resolve()
    admin_workspace = (tmp_path / "private" / "admin-workspaces").resolve()
    if endpoint_matches:
        admin_executable = original_executable
        admin_socket = original_socket
    else:
        admin_executable = (tmp_path / "private" / "admin" / "docker").resolve()
        admin_socket = (tmp_path / "private" / "admin.sock").resolve()
    private_input = str(tmp_path / "private" / "sample.csv")
    private_redis_url = "redis://operator:secret@private-redis.invalid:6379/0"
    durable_scope = managed_container_scope(original_workspace / RUN_ID)
    assignment = replace(
        _assignment(),
        managed_container_scope=durable_scope,
        managed_container_endpoint_identity=managed_container_endpoint_identity(
            original_executable,
            original_socket,
        ),
    )
    _persist_real_running_run(
        database_url,
        assignment=assignment,
        private_input=private_input,
    )
    before = _real_recovery_state(database_url)
    _AlwaysSuccessfulCleaner.cleanup_scopes = []
    monkeypatch.setattr(rq_queue_module, "RqRunQueue", _ExactTerminalRunQueue)
    monkeypatch.setattr(
        managed_containers_module,
        "ManagedContainerCleaner",
        _AlwaysSuccessfulCleaner,
    )
    monkeypatch.setenv(MANAGED_DOCKER_EXECUTABLE_ENV, str(admin_executable))
    monkeypatch.setenv(MANAGED_DOCKER_SOCKET_ENV, str(admin_socket))
    monkeypatch.setenv(WORKSPACE_ROOT_ENV, str(admin_workspace))
    monkeypatch.setenv(REDIS_URL_ENV, private_redis_url)

    diagnose_exit = admin.main(
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
    diagnose_output = capsys.readouterr()
    diagnosis = json.loads(diagnose_output.out)
    assert diagnose_exit == 0
    assert diagnose_output.err == ""
    assert diagnosis["cleanup"] == ("pending" if endpoint_matches else "blocked")
    assert diagnosis["allowed_actions"] == (["fail"] if endpoint_matches else [])
    if not endpoint_matches:
        doctor = run_recovery_doctor(
            database_url=database_url,
            redis_url=private_redis_url,
            queue_name="encode-pipeline-demo",
        )
        assert doctor.status is RecoveryDoctorStatus.ATTENTION_REQUIRED
        assert doctor.counts.cleanup_gaps == 1

    fail_exit = admin.main(
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
            "--notifications-env-file",
            str(notifications_environment),
        ]
    )
    fail_output = capsys.readouterr()
    if endpoint_matches:
        assert fail_exit == 0
        assert fail_output.err == ""
        assert json.loads(fail_output.out)["status"] == "failed"
        assert _AlwaysSuccessfulCleaner.cleanup_scopes == [durable_scope]
        after = _real_recovery_state(database_url)
        assert after[0].status is RunStatus.FAILED
        assert after[2][-1].event_type == "run_failed_by_admin_recovery"
    else:
        assert fail_exit == 1
        assert fail_output.out == ""
        assert json.loads(fail_output.err)["error"]["code"] == ("RUN_RECOVERY_NOT_SAFE")
        assert _AlwaysSuccessfulCleaner.cleanup_scopes == []
        assert _real_recovery_state(database_url) == before
        assert all(
            event.event_type != "run_failed_by_admin_recovery" for event in before[2]
        )

    assert managed_container_scope(admin_workspace / RUN_ID) != durable_scope

    combined_output = diagnose_output.out + diagnose_output.err + fail_output.err
    for private_value in (
        private_input,
        str(original_workspace),
        str(original_executable),
        str(original_socket),
        str(admin_executable),
        str(admin_socket),
        str(admin_workspace),
        private_redis_url,
        "operator:secret",
        "private-redis.invalid",
        str(database_path),
    ):
        assert private_value not in combined_output


def test_real_cli_reports_configured_but_unavailable_cleaner_as_blocked(
    tmp_path: Path,
    monkeypatch,
    capsys,
) -> None:
    database_path = tmp_path / "private" / "platform.db"
    database_path.parent.mkdir()
    database_url = f"sqlite:///{database_path}"
    notifications_environment = _disabled_notifications_environment(tmp_path)
    private_workspace = (tmp_path / "private" / "original-workspaces").resolve()
    private_executable = (tmp_path / "private" / "missing" / "docker").resolve()
    private_socket = (tmp_path / "private" / "missing.sock").resolve()
    private_input = str(tmp_path / "private" / "sample.csv")
    private_redis_url = "redis://operator:secret@private-redis.invalid:6379/0"
    assignment = replace(
        _assignment(),
        managed_container_scope=managed_container_scope(private_workspace / RUN_ID),
        managed_container_endpoint_identity=managed_container_endpoint_identity(
            private_executable,
            private_socket,
        ),
    )
    _persist_real_running_run(
        database_url,
        assignment=assignment,
        private_input=private_input,
    )
    before = _real_recovery_state(database_url)
    monkeypatch.setattr(rq_queue_module, "RqRunQueue", _ExactTerminalRunQueue)
    monkeypatch.setattr(
        managed_containers_module,
        "ManagedContainerCleaner",
        _UnavailableCleaner,
    )
    monkeypatch.setenv(MANAGED_DOCKER_EXECUTABLE_ENV, str(private_executable))
    monkeypatch.setenv(MANAGED_DOCKER_SOCKET_ENV, str(private_socket))
    monkeypatch.setenv(WORKSPACE_ROOT_ENV, str(tmp_path / "private" / "admin"))
    monkeypatch.setenv(REDIS_URL_ENV, private_redis_url)

    diagnose_exit = admin.main(
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
    diagnose_output = capsys.readouterr()
    diagnosis = json.loads(diagnose_output.out)
    assert diagnose_exit == 0
    assert diagnose_output.err == ""
    assert diagnosis["cleanup"] == "blocked"
    assert diagnosis["gaps"] == ["callback", "cleanup"]
    assert diagnosis["allowed_actions"] == []

    doctor = run_recovery_doctor(
        database_url=database_url,
        redis_url=private_redis_url,
        queue_name="encode-pipeline-demo",
    )
    assert doctor.status is RecoveryDoctorStatus.ATTENTION_REQUIRED
    assert doctor.counts.cleanup_gaps == 1

    fail_exit = admin.main(
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
            "--notifications-env-file",
            str(notifications_environment),
        ]
    )
    fail_output = capsys.readouterr()
    assert fail_exit == 1
    assert fail_output.out == ""
    assert json.loads(fail_output.err)["error"]["code"] == "RUN_RECOVERY_NOT_SAFE"
    assert _real_recovery_state(database_url) == before
    assert all(
        event.event_type != "run_failed_by_admin_recovery" for event in before[2]
    )

    combined_output = diagnose_output.out + diagnose_output.err + fail_output.err
    for private_value in (
        private_input,
        str(private_workspace),
        str(private_executable),
        str(private_socket),
        private_redis_url,
        "operator:secret",
        "private-redis.invalid",
        str(database_path),
    ):
        assert private_value not in combined_output
