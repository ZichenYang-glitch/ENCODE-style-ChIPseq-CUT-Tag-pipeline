"""Real Redis/RQ/Snakemake process-group cancellation integration."""

from __future__ import annotations

import asyncio
import json
import os
from pathlib import Path
import shlex
import shutil
import signal
import subprocess
import sys
from threading import Event, Thread
import time
from uuid import uuid4

import fastapi.routing
import httpx
import pytest
from redis import Redis
from rq import Queue
from rq.job import JobStatus
from rq.registry import FailedJobRegistry
from rq.serializers import JSONSerializer
import yaml

from encode_pipeline.api.main import create_app
from encode_pipeline.persistence.runtime import open_run_persistence
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.defaults import (
    create_default_run_service,
    create_default_workflow_registry,
)
from encode_pipeline.workers.settings import (
    QUEUE_NAME_ENV,
    REDIS_URL_ENV,
    WORKSPACE_ROOT_ENV,
)

from .process_helpers import terminate_rq_worker


REPOSITORY_ROOT = Path(__file__).resolve().parents[2]
PROFILE_ROOT = REPOSITORY_ROOT / "test" / "profiles" / "platform_worker_tiny"
TEST_REDIS_URL_ENV = "ENCODE_PIPELINE_TEST_REDIS_URL"
WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


def _request(app, method: str, url: str, **kwargs) -> httpx.Response:
    async def send() -> httpx.Response:
        transport = httpx.ASGITransport(app=app)
        async with httpx.AsyncClient(
            transport=transport,
            base_url="http://testserver",
            follow_redirects=True,
        ) as client:
            response = await client.request(method, url, **kwargs)
            await response.aread()
            return response

    return asyncio.run(send())


async def _run_in_joined_test_thread(function, *args, **kwargs):
    completed = Event()
    results: list[object] = []
    exceptions: list[BaseException] = []

    def invoke() -> None:
        try:
            results.append(function(*args, **kwargs))
        except BaseException as exc:
            exceptions.append(exc)
        finally:
            completed.set()

    thread = Thread(target=invoke)
    thread.start()
    try:
        while not completed.is_set():
            await asyncio.sleep(0.001)
    finally:
        thread.join(timeout=3)
        if thread.is_alive():  # pragma: no cover - test deadlock guard
            raise RuntimeError("test threadpool call did not terminate")
    if exceptions:
        raise exceptions[0]
    return results[0]


def test_real_worker_cancels_long_snakemake_process_group_truthfully(
    tmp_path,
    monkeypatch,
):
    redis_url = os.getenv(TEST_REDIS_URL_ENV)
    if redis_url is None:
        pytest.skip(f"{TEST_REDIS_URL_ENV} is not configured")
    real_snakemake = shutil.which("snakemake")
    assert real_snakemake is not None, (
        "snakemake must be available when the real Redis cancellation gate is enabled"
    )
    assert os.name == "posix", (
        "the enabled real Redis cancellation gate requires POSIX process groups"
    )

    project_root = tmp_path / "controlled-project"
    _create_controlled_project(project_root)
    marker_root = tmp_path / "process-markers"
    marker_root.mkdir()
    wrapper_dir = tmp_path / "test-bin"
    wrapper_dir.mkdir()
    _create_snakemake_wrapper(
        wrapper_dir / "snakemake",
        real_snakemake=real_snakemake,
    )

    queue_name = f"encode-pipeline-cancel-{uuid4().hex}"
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    workspace_root = tmp_path / "workspaces"
    monkeypatch.setenv("ENCODE_PIPELINE_DATABASE_URL", database_url)
    monkeypatch.setenv(REDIS_URL_ENV, redis_url)
    monkeypatch.setenv(QUEUE_NAME_ENV, queue_name)
    monkeypatch.setenv(WORKSPACE_ROOT_ENV, str(workspace_root))
    monkeypatch.setenv("TMPDIR", str(marker_root))
    monkeypatch.setenv(
        "PATH",
        os.pathsep.join([str(wrapper_dir), os.environ.get("PATH", os.defpath)]),
    )
    monkeypatch.setattr(
        fastapi.routing,
        "run_in_threadpool",
        _run_in_joined_test_thread,
    )

    connection = Redis.from_url(redis_url)
    queue = Queue(queue_name, connection=connection, serializer=JSONSerializer)
    worker_process: subprocess.Popen[str] | None = None
    app = create_app(
        database_url=database_url,
        workspace_root=workspace_root,
        project_root=project_root,
    )
    run_id: str | None = None
    process_group: int | None = None
    recorded_pids: tuple[int, ...] = ()
    try:
        assert connection.ping() is True
        config = yaml.safe_load(
            (PROFILE_ROOT / "config.yaml").read_text(encoding="utf-8")
        )
        samples_path = (PROFILE_ROOT / "samples.tsv").resolve()
        config["samples"] = str(samples_path)
        validation = _request(
            app,
            "POST",
            f"/api/v1/workflows/{WORKFLOW_ID}/validate",
            json={
                "config": config,
                "samples": str(samples_path),
                "options": {},
            },
        )
        assert validation.status_code == 200
        assert validation.json()["ok"] is True
        created = _request(
            app,
            "POST",
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json={"snapshot_id": validation.json()["snapshot"]["snapshot_id"]},
        )
        assert created.status_code == 201
        run_id = created.json()["run"]["run_id"]
        preflight = _request(app, "POST", f"/api/v1/runs/{run_id}/preflight")
        assert preflight.status_code == 202, preflight.text
        assert (
            _request(app, "GET", f"/api/v1/runs/{run_id}").json()["run"]["status"]
            == "planned"
        )
        started = _request(app, "POST", f"/api/v1/runs/{run_id}/start")
        assert started.status_code == 202, started.text
        assert started.json()["run"]["status"] == "queued"

        environment = dict(os.environ)
        environment.update(
            {
                "ENCODE_PIPELINE_DATABASE_URL": database_url,
                REDIS_URL_ENV: redis_url,
                QUEUE_NAME_ENV: queue_name,
                WORKSPACE_ROOT_ENV: str(workspace_root),
                "TMPDIR": str(marker_root),
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONPATH": str(project_root / "src"),
                "PATH": os.pathsep.join(
                    [str(wrapper_dir), environment.get("PATH", os.defpath)]
                ),
            }
        )
        worker_process = subprocess.Popen(
            [sys.executable, "-m", "encode_pipeline.workers.cli", "--burst"],
            cwd=project_root,
            env=environment,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            start_new_session=True,
        )

        assignment = _wait_for_running_assignment(database_url, run_id)
        _wait_for_path(marker_root / "entered")
        _wait_until(
            lambda: _logs_contain(database_url, run_id, "cancel-e2e-entered"),
            description="entered log line to be persisted before cancellation",
        )
        process_evidence = json.loads(
            (marker_root / "processes.json").read_text(encoding="utf-8")
        )
        snakemake_pid = int(
            (marker_root / "snakemake.pid").read_text(encoding="utf-8").strip()
        )
        process_group = int(process_evidence["process_group"])
        recorded_pids = (
            snakemake_pid,
            int(process_evidence["shell_pid"]),
            int(process_evidence["helper_pid"]),
            int(process_evidence["child_pid"]),
        )
        assert all(_process_exists(pid) for pid in recorded_pids)
        assert all(os.getpgid(pid) == process_group for pid in recorded_pids)

        job = queue.fetch_job(assignment.job_id)
        assert job is not None
        _wait_until(
            lambda: (
                job.get_status(refresh=True) is JobStatus.STARTED
                and bool(job.worker_name)
            ),
            description="RQ job to enter STARTED with a worker identity",
        )
        cancelled = _request(app, "POST", f"/api/v1/runs/{run_id}/cancel")
        assert cancelled.status_code in {200, 202}, cancelled.text
        cancellation_snapshot = cancelled.json()["run"]
        if cancelled.status_code == 202:
            assert cancellation_snapshot["status"] == "running"
            assert cancellation_snapshot["ended_at"] is None
            assert cancellation_snapshot["cancellation_reason"] is None
        else:
            assert cancellation_snapshot["status"] == "cancelled"
            assert cancellation_snapshot["ended_at"] is not None
            assert cancellation_snapshot["cancellation_reason"] == (
                "User requested cancellation."
            )

        stdout, stderr = worker_process.communicate(timeout=20)
        assert worker_process.returncode == 0, (stdout, stderr)
        _wait_until(
            lambda: _read_status(database_url, run_id) is RunStatus.CANCELLED,
            description="SQLite cancellation acknowledgement",
        )
        _wait_until(
            lambda: all(not _process_exists(pid) for pid in recorded_pids),
            description="Snakemake descendant processes to be reaped",
        )
        _wait_until(
            lambda: not _process_group_exists(process_group),
            description="RQ horse process group to disappear",
        )

        persistence = open_run_persistence(database_url)
        try:
            service = create_default_run_service(
                create_default_workflow_registry(),
                repository=persistence.repository,
            )
            record = service.get_run(run_id)
            persisted_assignment = service.get_execution_assignment(run_id)
            events = service.list_events(run_id, limit=100)
            logs = (
                *service.list_logs(run_id, "stdout", limit=100),
                *service.list_logs(run_id, "stderr", limit=100),
            )
        finally:
            persistence.close()

        assert record.status is RunStatus.CANCELLED
        assert record.ended_at is not None
        assert record.cancellation_reason == "User requested cancellation."
        assert persisted_assignment is not None
        assert persisted_assignment.cancellation_requested_at is not None
        assert persisted_assignment.cancellation_acknowledged_at is not None
        assert [event.event_type for event in events].count(
            "cancellation_requested"
        ) == 1
        assert [event.event_type for event in events].count(
            "cancellation_acknowledged"
        ) == 1
        terminal_events = [
            event
            for event in events
            if event.status is not None and event.status.is_terminal
        ]
        assert len(terminal_events) == 1
        assert terminal_events[0].status is RunStatus.CANCELLED
        persisted_logs = "\n".join(line for chunk in logs for line in chunk.lines)
        assert "cancel-e2e-entered" in persisted_logs

        assert not (marker_root / "helper-completed").exists()
        assert not (workspace_root / run_id / "result" / "complete.txt").exists()
        job.refresh()
        assert job.get_status() is JobStatus.STOPPED
        assert job.is_stopped
        assert not job.is_failed
        assert (
            assignment.job_id
            in FailedJobRegistry(
                queue_name,
                connection=connection,
                serializer=JSONSerializer,
            ).get_job_ids()
        )
        failure = job.latest_result()
        assert failure is not None
        assert failure.exc_string == "Job stopped by user, work-horse terminated."
    finally:
        if worker_process is not None:
            terminate_rq_worker(worker_process)
        if process_group is not None and _process_group_exists(process_group):
            try:
                os.killpg(process_group, signal.SIGKILL)
            except ProcessLookupError:
                pass
        if run_id is not None:
            persistence = open_run_persistence(database_url)
            try:
                service = create_default_run_service(
                    create_default_workflow_registry(),
                    repository=persistence.repository,
                )
                assignment = service.get_execution_assignment(run_id)
                if assignment is not None:
                    job = queue.fetch_job(assignment.job_id)
                    if job is not None:
                        job.delete()
            finally:
                persistence.close()
        queue.delete()
        connection.close()
        app.state.run_queue.close()
        app.state.persistence.close()


def _create_controlled_project(project_root: Path) -> None:
    project_root.mkdir()
    shutil.copy2(REPOSITORY_ROOT / "pyproject.toml", project_root / "pyproject.toml")
    shutil.copytree(
        REPOSITORY_ROOT / "src" / "encode_pipeline",
        project_root / "src" / "encode_pipeline",
    )
    inventory = project_root / "docs/architecture/artifact-inventory.yaml"
    inventory.parent.mkdir(parents=True)
    shutil.copy2(
        REPOSITORY_ROOT / "docs/architecture/artifact-inventory.yaml",
        inventory,
    )
    workflow = project_root / "workflow"
    profile = project_root / "profiles" / "default"
    scripts = project_root / "scripts"
    workflow.mkdir()
    profile.mkdir(parents=True)
    scripts.mkdir()
    (profile / "config.yaml").write_text("printshellcmds: true\n", encoding="utf-8")
    (workflow / "Snakefile").write_text(
        """
from pathlib import Path

HELPER = Path(workflow.basedir).parent / "scripts" / "cancellation_process.py"

rule all:
    input:
        "result/complete.txt"

rule long_running:
    output:
        "result/complete.txt"
    params:
        helper=str(HELPER)
    shell:
        "python3 {params.helper:q} {output:q}"
""".lstrip(),
        encoding="utf-8",
    )
    (scripts / "cancellation_process.py").write_text(
        """
from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys
import time

marker_root = Path(os.environ["TMPDIR"])
output = Path(sys.argv[1])
child = subprocess.Popen(
    [sys.executable, "-c", "import time; time.sleep(120)"],
)
evidence = {
    "shell_pid": os.getppid(),
    "helper_pid": os.getpid(),
    "child_pid": child.pid,
    "process_group": os.getpgrp(),
}
temporary = marker_root / "processes.json.tmp"
temporary.write_text(json.dumps(evidence), encoding="utf-8")
temporary.replace(marker_root / "processes.json")
print("cancel-e2e-entered", flush=True)
(marker_root / "entered").write_text("entered\\n", encoding="utf-8")
time.sleep(120)
child.wait(timeout=5)
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text("completed\\n", encoding="utf-8")
(marker_root / "helper-completed").write_text("completed\\n", encoding="utf-8")
""".lstrip(),
        encoding="utf-8",
    )


def _create_snakemake_wrapper(path: Path, *, real_snakemake: str) -> None:
    path.write_text(
        "#!/bin/sh\n"
        'printf "%s\\n" "$$" > "$TMPDIR/snakemake.pid"\n'
        f'exec {shlex.quote(real_snakemake)} "$@"\n',
        encoding="utf-8",
    )
    path.chmod(0o755)


def _wait_for_running_assignment(database_url: str, run_id: str):
    assignment_holder = []

    def running() -> bool:
        persistence = open_run_persistence(database_url)
        try:
            service = create_default_run_service(
                create_default_workflow_registry(),
                repository=persistence.repository,
            )
            record = service.get_run(run_id)
            assignment = service.get_execution_assignment(run_id)
            if record.status is RunStatus.RUNNING and assignment is not None:
                assignment_holder[:] = [assignment]
                return True
            if record.status.is_terminal:
                raise AssertionError(
                    f"run became terminal before cancellation: {record}"
                )
            return False
        finally:
            persistence.close()

    _wait_until(running, description="worker to claim and start the run")
    return assignment_holder[0]


def _read_status(database_url: str, run_id: str) -> RunStatus:
    persistence = open_run_persistence(database_url)
    try:
        return persistence.repository.get_run(run_id).status
    finally:
        persistence.close()


def _logs_contain(database_url: str, run_id: str, needle: str) -> bool:
    persistence = open_run_persistence(database_url)
    try:
        service = create_default_run_service(
            create_default_workflow_registry(),
            repository=persistence.repository,
        )
        chunks = (
            *service.list_logs(run_id, "stdout", limit=100),
            *service.list_logs(run_id, "stderr", limit=100),
        )
        return any(needle in line for chunk in chunks for line in chunk.lines)
    finally:
        persistence.close()


def _wait_for_path(path: Path) -> None:
    _wait_until(path.is_file, description=f"{path.name} marker")


def _wait_until(predicate, *, description: str, timeout_seconds: float = 15) -> None:
    deadline = time.monotonic() + timeout_seconds
    while time.monotonic() < deadline:
        if predicate():
            return
        time.sleep(0.05)
    raise AssertionError(f"Timed out waiting for {description}")


def _process_exists(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True


def _process_group_exists(process_group: int) -> bool:
    try:
        os.killpg(process_group, 0)
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    return True
