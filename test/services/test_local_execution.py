"""Tests for durable worker-side local execution orchestration."""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs, WorkspacePlan
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.local_execution import LocalExecutionService
from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
from encode_pipeline.services.process_runner import ProcessResult, ProcessRunner
from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.workers.timeouts import WorkerHardTimeout


WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"
EMPTY_SAMPLES = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome"
    "\tbowtie2_index\n"
)


class ControlledRunner(ProcessRunner):
    def __init__(
        self,
        *,
        exit_code: int = 0,
        fail: bool = False,
        failure_code: str = "PROCESS_RUNNER_EXECUTION_ERROR",
        callback=None,
        output_truncated: bool = False,
    ) -> None:
        super().__init__(allowed_executables=("snakemake",))
        self.exit_code = exit_code
        self.fail = fail
        self.failure_code = failure_code
        self.callback = callback
        self.output_truncated = output_truncated
        self.specs = []

    def run(self, spec, *, output_callback=None):
        self.specs.append(spec)
        if self.callback is not None:
            self.callback()
        if self.fail:
            return Result.failure(
                [
                    Issue(
                        code=self.failure_code,
                        message="Subprocess execution failed.",
                        source="process_runner",
                    )
                ]
            )
        if output_callback is not None:
            output_callback("stdout", ("first", "second"))
            output_callback("stderr", ("warning",))
        issues = []
        if self.output_truncated:
            issues.append(
                Issue(
                    code="PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED",
                    message="Subprocess output was truncated.",
                    severity="warning",
                    source="process_runner",
                )
            )
        if self.exit_code:
            issues.append(
                Issue(
                    code="PROCESS_RUNNER_NONZERO_EXIT",
                    message="Subprocess exited with a non-zero code.",
                    severity="warning",
                    source="process_runner",
                    context={"exit_code": self.exit_code},
                )
            )
        return Result.success(
            ProcessResult(exit_code=self.exit_code, stdout="", stderr=""),
            issues=issues,
        )


def _prepared_service(tmp_path: Path, runner: ProcessRunner):
    registry = create_default_workflow_registry()
    run_service = RunService(registry, id_factory=lambda: "run-1")
    samples = tmp_path / "empty-samples.tsv"
    samples.write_text(EMPTY_SAMPLES, encoding="utf-8")
    run_service.create_run(
        WORKFLOW_ID,
        WorkflowInputs(config={"samples": str(samples), "multiqc": False}),
    )
    run_service.transition_run("run-1", RunStatus.VALIDATING, stage="preflight")

    execution_planner = ExecutionPlanner(run_service)
    workspace_planner = WorkspacePlanner(registry)
    command_builder = CommandBuilder(registry)
    materializer = WorkspaceMaterializer()
    workspace_root = tmp_path / "workspaces"
    driver = LocalRunDriver(
        run_service=run_service,
        materializer=materializer,
        command_builder=command_builder,
        workspace_root=workspace_root,
        process_runner=runner,
    )
    base_plan = execution_planner.plan_run("run-1").value
    workspace_dir = driver.derive_workspace_dir("run-1")
    workspace_plan = workspace_planner.plan_workspace(base_plan, workspace_dir).value
    assert materializer.materialize(
        workspace_plan.workspace_plan,
        workspace_dir,
    ).is_success
    run_service.transition_run("run-1", RunStatus.PLANNED, stage="preflight")
    run_service.transition_run("run-1", RunStatus.QUEUED, stage="execution")

    execution_service = LocalExecutionService(
        run_service=run_service,
        execution_planner=execution_planner,
        workspace_planner=workspace_planner,
        command_builder=command_builder,
        local_run_driver=driver,
    )
    return execution_service, run_service, workspace_dir


def _service_failure(code: str = "CONTROLLED_SERVICE_FAILURE") -> Result:
    return Result.failure(
        [
            Issue(
                code=code,
                message="A controlled service boundary failed.",
                source="test",
            )
        ]
    )


def _verification_plan(
    *,
    directories: tuple[str, ...] = (),
    files: tuple[tuple[str, bytes], ...] = (),
    workspace_plan: WorkspacePlan | None = None,
) -> ExecutionPlan:
    return ExecutionPlan(
        plan_id="plan-verification",
        run_id="run-verification",
        workflow_id=WORKFLOW_ID,
        status=PlanStatus.PLANNED,
        inputs_snapshot={"config": {}},
        workspace_plan=(
            workspace_plan
            if workspace_plan is not None
            else WorkspacePlan(directories=directories, files=files)
        ),
    )


def test_local_execution_rebuilds_existing_workspace_and_succeeds(tmp_path):
    runner = ControlledRunner()
    service, run_service, _workspace = _prepared_service(tmp_path, runner)

    result = service.execute("run-1")

    assert result.is_success
    assert result.value.status is RunStatus.SUCCEEDED
    assert result.value.started_at is not None
    assert result.value.ended_at is not None
    assert len(runner.specs) == 1
    assert "-n" not in runner.specs[0].argv
    assert "--dry-run" not in runner.specs[0].argv
    assert run_service.list_logs("run-1", "stdout")[0].lines == (
        "first",
        "second",
    )
    assert run_service.list_logs("run-1", "stderr")[0].lines == ("warning",)


def test_local_execution_maps_nonzero_exit_to_failed_run(tmp_path):
    runner = ControlledRunner(exit_code=17)
    service, run_service, _workspace = _prepared_service(tmp_path, runner)

    result = service.execute("run-1")

    assert result.is_failure
    record = run_service.get_run("run-1")
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_EXECUTION_FAILED"
    assert record.error.context == {"reason_code": "LOCAL_RUN_EXECUTION_FAILED"}
    assert run_service.list_logs("run-1", "stderr")[0].lines == ("warning",)


def test_local_execution_persists_and_propagates_output_truncation_warning(tmp_path):
    runner = ControlledRunner(output_truncated=True)
    service, run_service, _workspace = _prepared_service(tmp_path, runner)

    result = service.execute("run-1")

    assert result.is_success
    assert [issue.code for issue in result.issues] == [
        "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED"
    ]
    events = run_service.list_events("run-1", limit=100)
    truncation = [
        event for event in events if event.event_type == "execution_output_truncated"
    ]
    assert len(truncation) == 1
    assert truncation[0].issue is not None
    assert truncation[0].issue.code == "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED"


def test_local_execution_maps_runner_infrastructure_failure(tmp_path):
    service, run_service, _workspace = _prepared_service(
        tmp_path,
        ControlledRunner(fail=True),
    )

    result = service.execute("run-1")

    assert result.is_failure
    record = run_service.get_run("run-1")
    assert record.status is RunStatus.FAILED
    assert record.error.context == {"reason_code": "LOCAL_RUN_PROCESS_FAILED"}


def test_local_execution_persists_the_process_runner_timeout_reason(tmp_path):
    service, run_service, _workspace = _prepared_service(
        tmp_path,
        ControlledRunner(fail=True, failure_code="PROCESS_RUNNER_TIMEOUT"),
    )

    result = service.execute("run-1")

    assert result.is_failure
    record = run_service.get_run("run-1")
    assert record.status is RunStatus.FAILED
    assert record.error is not None
    assert record.error.code == "RUN_EXECUTION_FAILED"
    assert record.error.context == {"reason_code": "PROCESS_RUNNER_TIMEOUT"}


def test_local_execution_rejects_changed_preflight_workspace(tmp_path):
    runner = ControlledRunner()
    service, run_service, workspace = _prepared_service(tmp_path, runner)
    (workspace / "config" / "config.yaml").write_text(
        "changed: true\n",
        encoding="utf-8",
    )

    result = service.execute("run-1")

    assert result.is_failure
    assert runner.specs == []
    record = run_service.get_run("run-1")
    assert record.status is RunStatus.FAILED
    assert record.error.context == {"reason_code": "LOCAL_EXECUTION_WORKSPACE_INVALID"}


def test_local_execution_does_not_start_process_when_queued_cancel_wins(
    tmp_path,
    monkeypatch,
):
    runner = ControlledRunner()
    service, run_service, _workspace = _prepared_service(tmp_path, runner)
    original_transition = run_service.transition_run
    cancelled = False

    def cancel_before_running(run_id, to_status, **kwargs):
        nonlocal cancelled
        if to_status is RunStatus.RUNNING and not cancelled:
            cancelled = True
            run_service.cancel_run(run_id, reason="race")
        return original_transition(run_id, to_status, **kwargs)

    monkeypatch.setattr(run_service, "transition_run", cancel_before_running)

    result = service.execute("run-1")

    assert result.is_failure
    assert result.issues[0].code == "LOCAL_EXECUTION_CANCELLED"
    assert runner.specs == []
    assert run_service.get_run("run-1").status is RunStatus.CANCELLED
    assert all(
        event.status is not RunStatus.FAILED
        for event in run_service.list_events("run-1", limit=100)
    )


def test_local_execution_propagates_rq_timeout_to_worker_boundary(tmp_path):
    def timeout():
        raise WorkerHardTimeout("RQ deadline reached")

    runner = ControlledRunner(callback=timeout)
    service, run_service, _workspace = _prepared_service(tmp_path, runner)

    with pytest.raises(WorkerHardTimeout, match="RQ deadline reached"):
        service.execute("run-1")

    assert run_service.get_run("run-1").status is RunStatus.RUNNING


def test_local_execution_rejects_missing_cancelled_and_nonqueued_work(
    tmp_path,
) -> None:
    service, run_service, _workspace = _prepared_service(
        tmp_path,
        ControlledRunner(),
    )

    missing = service.execute("missing")
    assert missing.is_failure
    assert missing.errors[0].code == "LOCAL_EXECUTION_RUN_NOT_FOUND"

    run_service.cancel_run("run-1", reason="cancelled before worker claim")
    cancelled = service.execute("run-1")
    assert cancelled.is_failure
    assert cancelled.errors[0].code == "LOCAL_EXECUTION_CANCELLED"

    second_root = tmp_path / "second"
    second_root.mkdir()
    second_service, second_runs, _workspace = _prepared_service(
        second_root,
        ControlledRunner(),
    )
    second_runs.transition_run("run-1", RunStatus.RUNNING)
    nonqueued = second_service.execute("run-1")
    assert nonqueued.is_failure
    assert nonqueued.errors[0].code == "LOCAL_EXECUTION_NOT_QUEUED"
    assert nonqueued.errors[0].context == {"current_status": "running"}


@pytest.mark.parametrize(
    ("race_status", "expected_code"),
    (
        pytest.param(
            RunStatus.CANCELLED,
            "LOCAL_EXECUTION_CANCELLED",
            id="cancelled-after-rebuild",
        ),
        pytest.param(
            RunStatus.FAILED,
            "LOCAL_EXECUTION_STATE_CHANGED",
            id="terminal-state-after-rebuild",
        ),
    ),
)
def test_local_execution_rechecks_lifecycle_after_rebuilding_dependencies(
    tmp_path,
    monkeypatch,
    race_status,
    expected_code,
) -> None:
    service, run_service, _workspace = _prepared_service(
        tmp_path,
        ControlledRunner(),
    )
    original_rebuild = service._rebuild_plan

    def rebuild_then_race(run_id):
        result = original_rebuild(run_id)
        assert result.is_success
        if race_status is RunStatus.CANCELLED:
            run_service.cancel_run(run_id, reason="cancelled during rebuild")
        else:
            run_service.transition_run(run_id, RunStatus.FAILED)
        return result

    monkeypatch.setattr(service, "_rebuild_plan", rebuild_then_race)

    result = service.execute("run-1")

    assert result.is_failure
    assert result.errors[0].code == expected_code
    assert run_service.get_run("run-1").status is race_status


def test_local_execution_rejects_nonrunning_transition_race(
    tmp_path,
    monkeypatch,
) -> None:
    service, run_service, _workspace = _prepared_service(
        tmp_path,
        ControlledRunner(),
    )
    original_transition = run_service.transition_run

    def fail_before_running(run_id, to_status, **kwargs):
        if to_status is RunStatus.RUNNING:
            original_transition(run_id, RunStatus.FAILED)
        return original_transition(run_id, to_status, **kwargs)

    monkeypatch.setattr(run_service, "transition_run", fail_before_running)

    result = service.execute("run-1")

    assert result.is_failure
    assert result.errors[0].code == "LOCAL_EXECUTION_STATE_CHANGED"
    assert result.errors[0].context == {"current_status": "failed"}


@pytest.mark.parametrize(
    ("race_status", "expected_code"),
    (
        pytest.param(
            RunStatus.CANCELLED,
            "LOCAL_EXECUTION_CANCELLED",
            id="cancelled-after-driver",
        ),
        pytest.param(
            RunStatus.FAILED,
            "LOCAL_EXECUTION_STATE_CHANGED",
            id="terminal-state-after-driver",
        ),
    ),
)
def test_local_execution_rechecks_lifecycle_after_driver_returns(
    tmp_path,
    monkeypatch,
    race_status,
    expected_code,
) -> None:
    service, run_service, _workspace = _prepared_service(
        tmp_path,
        ControlledRunner(),
    )
    raced = False
    original_get = run_service.get_run

    def run_then_race(_run_id, _plan):
        nonlocal raced
        if race_status is RunStatus.FAILED:
            run_service.transition_run("run-1", RunStatus.FAILED)
        else:
            raced = True
        return Result.success(object())

    def observe_race(run_id):
        record = original_get(run_id)
        if raced and record.status is RunStatus.RUNNING:
            return replace(record, status=RunStatus.CANCELLED)
        return record

    monkeypatch.setattr(service._local_run_driver, "run", run_then_race)
    monkeypatch.setattr(run_service, "get_run", observe_race)

    result = service.execute("run-1")

    assert result.is_failure
    assert result.errors[0].code == expected_code


@pytest.mark.parametrize(
    ("race_status", "expected_success", "expected_code"),
    (
        pytest.param(
            RunStatus.CANCELLED,
            False,
            "LOCAL_EXECUTION_CANCELLED",
            id="cancelled-at-success-commit",
        ),
        pytest.param(
            RunStatus.SUCCEEDED,
            True,
            None,
            id="success-commit-won",
        ),
        pytest.param(
            RunStatus.FAILED,
            False,
            "LOCAL_EXECUTION_STATE_CHANGED",
            id="other-terminal-won",
        ),
    ),
)
def test_local_execution_reconciles_terminal_success_commit_races(
    tmp_path,
    monkeypatch,
    race_status,
    expected_success,
    expected_code,
) -> None:
    service, run_service, _workspace = _prepared_service(
        tmp_path,
        ControlledRunner(),
    )
    original_transition = run_service.transition_run
    original_get = run_service.get_run
    cancelled_snapshot = False

    monkeypatch.setattr(
        service._local_run_driver,
        "run",
        lambda _run_id, _plan: Result.success(object()),
    )

    def race_success_commit(run_id, to_status, **kwargs):
        nonlocal cancelled_snapshot
        if to_status is not RunStatus.SUCCEEDED:
            return original_transition(run_id, to_status, **kwargs)
        if race_status is RunStatus.CANCELLED:
            cancelled_snapshot = True
            raise ConcurrentRunUpdateError("cancellation won")
        original_transition(run_id, race_status)
        raise ConcurrentRunUpdateError("terminal transition won")

    def observe_terminal_race(run_id):
        record = original_get(run_id)
        if cancelled_snapshot and record.status is RunStatus.RUNNING:
            return replace(record, status=RunStatus.CANCELLED)
        return record

    monkeypatch.setattr(run_service, "transition_run", race_success_commit)
    monkeypatch.setattr(run_service, "get_run", observe_terminal_race)

    result = service.execute("run-1")

    assert result.is_success is expected_success
    if expected_code is not None:
        assert result.errors[0].code == expected_code


@pytest.mark.parametrize(
    "failing_boundary",
    (
        pytest.param("execution-planner", id="planning"),
        pytest.param("workspace-planner", id="workspace"),
        pytest.param("command-builder", id="command"),
    ),
)
def test_local_execution_persists_rebuild_boundary_failures(
    tmp_path,
    monkeypatch,
    failing_boundary,
) -> None:
    service, run_service, _workspace = _prepared_service(
        tmp_path,
        ControlledRunner(),
    )
    if failing_boundary == "execution-planner":
        monkeypatch.setattr(
            service._execution_planner,
            "plan_run",
            lambda _run_id: _service_failure("EXECUTION_PLAN_REBUILD_FAILED"),
        )
    elif failing_boundary == "workspace-planner":
        monkeypatch.setattr(
            service._workspace_planner,
            "plan_workspace",
            lambda _plan, base_dir: _service_failure("WORKSPACE_REBUILD_FAILED"),
        )
    else:
        monkeypatch.setattr(
            service._command_builder,
            "build_command",
            lambda _plan, _base_dir: _service_failure("COMMAND_REBUILD_FAILED"),
        )

    result = service.execute("run-1")

    assert result.is_failure
    assert run_service.get_run("run-1").status is RunStatus.FAILED
    assert run_service.get_run("run-1").error.context == {
        "reason_code": result.errors[0].code
    }


def test_workspace_verification_rejects_each_materialized_shape_failure(
    tmp_path,
    monkeypatch,
) -> None:
    workspace = tmp_path / "workspace"
    workspace.mkdir()

    no_plan = _verification_plan(workspace_plan=None)
    no_plan = replace(no_plan, workspace_plan=None)
    assert LocalExecutionService._verify_materialized_workspace(
        no_plan,
        workspace,
    ).is_failure

    missing_root = _verification_plan()
    assert LocalExecutionService._verify_materialized_workspace(
        missing_root,
        tmp_path / "missing",
    ).is_failure

    (workspace / "not-a-directory").write_text("file", encoding="utf-8")
    bad_directory = _verification_plan(directories=("not-a-directory",))
    assert LocalExecutionService._verify_materialized_workspace(
        bad_directory,
        workspace,
    ).is_failure

    (workspace / "not-a-file").mkdir()
    bad_file = _verification_plan(files=(("not-a-file", b"expected"),))
    assert LocalExecutionService._verify_materialized_workspace(
        bad_file,
        workspace,
    ).is_failure

    readable = workspace / "readable"
    readable.write_bytes(b"expected")
    read_error = _verification_plan(files=(("readable", b"expected"),))
    original_read_bytes = Path.read_bytes

    def fail_target_read(path):
        if path == readable:
            raise OSError("controlled read failure")
        return original_read_bytes(path)

    monkeypatch.setattr(Path, "read_bytes", fail_target_read)
    assert LocalExecutionService._verify_materialized_workspace(
        read_error,
        workspace,
    ).is_failure


def test_local_execution_synthesizes_empty_driver_failure_and_retries_state_write(
    tmp_path,
    monkeypatch,
) -> None:
    service, run_service, _workspace = _prepared_service(
        tmp_path,
        ControlledRunner(),
    )
    original_transition = run_service.transition_run
    failed_writes = 0

    def reject_failure_writes(run_id, to_status, **kwargs):
        nonlocal failed_writes
        if to_status is RunStatus.FAILED:
            failed_writes += 1
            raise ConcurrentRunUpdateError("controlled persistence race")
        return original_transition(run_id, to_status, **kwargs)

    monkeypatch.setattr(run_service, "transition_run", reject_failure_writes)

    result = service._fail("run-1", ())

    assert result.is_failure
    assert result.errors[0].code == "LOCAL_EXECUTION_FAILED"
    assert failed_writes == 2
    assert original_transition is not None
    assert run_service.get_run("run-1").status is RunStatus.QUEUED


@pytest.mark.parametrize(
    ("failure_state", "expected_code"),
    (
        pytest.param("missing", "CONTROLLED_FAILURE", id="run-disappeared"),
        pytest.param("cancelled", "LOCAL_EXECUTION_CANCELLED", id="cancelled"),
        pytest.param("terminal", "CONTROLLED_FAILURE", id="already-terminal"),
        pytest.param("not-executable", "CONTROLLED_FAILURE", id="not-executable"),
    ),
)
def test_failure_persistence_never_overwrites_newer_lifecycle_state(
    tmp_path,
    monkeypatch,
    failure_state,
    expected_code,
) -> None:
    service, run_service, _workspace = _prepared_service(
        tmp_path,
        ControlledRunner(),
    )
    issue = Issue(
        code="CONTROLLED_FAILURE",
        message="Controlled failure.",
        source="test",
    )

    if failure_state == "missing":
        monkeypatch.setattr(
            run_service,
            "get_run",
            lambda _run_id: (_ for _ in ()).throw(KeyError("gone")),
        )
    elif failure_state == "cancelled":
        run_service.cancel_run("run-1", reason="cancelled before failure commit")
    elif failure_state == "terminal":
        run_service.transition_run("run-1", RunStatus.FAILED)
    else:
        queued = run_service.get_run("run-1")
        monkeypatch.setattr(
            run_service,
            "get_run",
            lambda _run_id: replace(queued, status=RunStatus.PLANNED),
        )

    result = service._fail("run-1", (issue,))

    assert result.is_failure
    assert result.errors[0].code == expected_code
