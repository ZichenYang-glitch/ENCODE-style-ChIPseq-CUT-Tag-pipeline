"""Tests for durable worker-side local execution orchestration."""

from __future__ import annotations

from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.command_builder import CommandBuilder
from encode_pipeline.services.local_execution import LocalExecutionService
from encode_pipeline.services.local_run_driver import LocalRunDriver
from encode_pipeline.services.materialization import WorkspaceMaterializer
from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
from encode_pipeline.services.process_runner import ProcessResult, ProcessRunner
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
        callback=None,
        output_truncated: bool = False,
    ) -> None:
        super().__init__(allowed_executables=("snakemake",))
        self.exit_code = exit_code
        self.fail = fail
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
                        code="PROCESS_RUNNER_EXECUTION_ERROR",
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


def test_local_execution_does_not_overwrite_concurrent_cancellation(tmp_path):
    holder = {}

    def cancel():
        holder["run_service"].cancel_run("run-1", reason="race")

    runner = ControlledRunner(callback=cancel)
    service, run_service, _workspace = _prepared_service(tmp_path, runner)
    holder["run_service"] = run_service

    result = service.execute("run-1")

    assert result.is_failure
    assert result.issues[0].code == "LOCAL_EXECUTION_CANCELLED"
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
