"""Durable local execution orchestration for worker-owned runs."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from encode_pipeline.platform.planning import WorkspacePathError, WorkspacePathPolicy
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError

if TYPE_CHECKING:
    from encode_pipeline.platform.planning import ExecutionPlan
    from encode_pipeline.services.command_builder import CommandBuilder
    from encode_pipeline.services.local_run_driver import LocalRunDriver
    from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
    from encode_pipeline.services.runs import RunService


class LocalExecutionService:
    """Rebuild and execute a preflighted workspace under durable lifecycle state."""

    def __init__(
        self,
        *,
        run_service: "RunService",
        execution_planner: "ExecutionPlanner",
        workspace_planner: "WorkspacePlanner",
        command_builder: "CommandBuilder",
        local_run_driver: "LocalRunDriver",
    ) -> None:
        from encode_pipeline.services.command_builder import CommandBuilder
        from encode_pipeline.services.local_run_driver import LocalRunDriver
        from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
        from encode_pipeline.services.runs import RunService

        if not isinstance(run_service, RunService):
            raise ValueError("LocalExecutionService requires a RunService instance")
        if not isinstance(execution_planner, ExecutionPlanner):
            raise ValueError(
                "LocalExecutionService requires an ExecutionPlanner instance"
            )
        if not isinstance(workspace_planner, WorkspacePlanner):
            raise ValueError(
                "LocalExecutionService requires a WorkspacePlanner instance"
            )
        if not isinstance(command_builder, CommandBuilder):
            raise ValueError("LocalExecutionService requires a CommandBuilder instance")
        if not isinstance(local_run_driver, LocalRunDriver):
            raise ValueError("LocalExecutionService requires a LocalRunDriver instance")

        self._run_service = run_service
        self._execution_planner = execution_planner
        self._workspace_planner = workspace_planner
        self._command_builder = command_builder
        self._local_run_driver = local_run_driver

    def execute(self, run_id: str) -> Result[RunRecord]:
        """Execute one claimed QUEUED run and persist its terminal outcome."""
        try:
            current = self._run_service.get_run(run_id)
        except KeyError:
            return Result.failure(
                [
                    Issue(
                        code="LOCAL_EXECUTION_RUN_NOT_FOUND",
                        message="Run was not found.",
                        severity="error",
                        path="run_id",
                        source="local_execution_service",
                    )
                ]
            )

        if current.status is RunStatus.CANCELLED:
            return self._cancelled_result()
        if current.status is not RunStatus.QUEUED:
            return Result.failure(
                [
                    Issue(
                        code="LOCAL_EXECUTION_NOT_QUEUED",
                        message="Run is not queued for worker execution.",
                        severity="error",
                        path="run_id",
                        source="local_execution_service",
                        context={"current_status": current.status.value},
                    )
                ]
            )

        rebuilt = self._rebuild_plan(run_id)
        if rebuilt.is_failure:
            return self._fail(run_id, rebuilt.issues)

        current = self._run_service.get_run(run_id)
        if current.status is RunStatus.CANCELLED:
            return self._cancelled_result()
        if current.status is not RunStatus.QUEUED:
            return self._state_changed_result(current)

        try:
            self._run_service.transition_run(
                run_id,
                RunStatus.RUNNING,
                stage="execution",
                message="Worker started local workflow execution.",
                context={"reason_code": "LOCAL_EXECUTION_STARTED"},
            )
        except (ConcurrentRunUpdateError, ValueError):
            current = self._run_service.get_run(run_id)
            if current.status is RunStatus.CANCELLED:
                return self._cancelled_result()
            if current.status is not RunStatus.RUNNING:
                return self._state_changed_result(current)

        execution_result = self._local_run_driver.run(run_id, rebuilt.value)
        if execution_result.is_failure:
            return self._fail(run_id, execution_result.issues)

        current = self._run_service.get_run(run_id)
        if current.status is RunStatus.CANCELLED:
            return self._cancelled_result()
        if current.status is not RunStatus.RUNNING:
            return self._state_changed_result(current)

        try:
            completed = self._run_service.transition_run(
                run_id,
                RunStatus.SUCCEEDED,
                stage="execution",
                message="Local workflow execution completed successfully.",
                context={"reason_code": "LOCAL_EXECUTION_SUCCEEDED"},
            )
        except (ConcurrentRunUpdateError, ValueError):
            current = self._run_service.get_run(run_id)
            if current.status is RunStatus.CANCELLED:
                return self._cancelled_result()
            if current.status is RunStatus.SUCCEEDED:
                return Result.success(current, issues=execution_result.issues)
            return self._state_changed_result(current)
        return Result.success(completed, issues=execution_result.issues)

    def _rebuild_plan(self, run_id: str) -> Result["ExecutionPlan"]:
        """Rebuild the command from SQLite inputs and verify existing files."""
        plan_result = self._execution_planner.plan_run(run_id)
        if plan_result.is_failure:
            return plan_result

        workspace_dir = self._local_run_driver.derive_workspace_dir(run_id)
        workspace_result = self._workspace_planner.plan_workspace(
            plan_result.value,
            base_dir=workspace_dir,
        )
        if workspace_result.is_failure:
            return workspace_result

        verification = self._verify_materialized_workspace(
            workspace_result.value,
            workspace_dir,
        )
        if verification.is_failure:
            return verification

        command_result = self._command_builder.build_command(
            workspace_result.value,
            workspace_dir,
        )
        if command_result.is_failure:
            return command_result
        return command_result

    @staticmethod
    def _verify_materialized_workspace(
        plan: "ExecutionPlan",
        workspace_dir: Path,
    ) -> Result[None]:
        """Verify that preflight files still match the reconstructed plan."""
        if plan.workspace_plan is None:
            return LocalExecutionService._workspace_failure()
        try:
            policy = WorkspacePathPolicy(base_dir=workspace_dir)
            if not workspace_dir.is_dir() or LocalExecutionService._has_symlink(
                workspace_dir
            ):
                return LocalExecutionService._workspace_failure()

            for directory in plan.workspace_plan.directories:
                resolved = policy.resolve(directory)
                if not resolved.is_dir() or LocalExecutionService._has_symlink(
                    resolved
                ):
                    return LocalExecutionService._workspace_failure()

            for relative_path, expected_contents in plan.workspace_plan.files:
                resolved = policy.resolve(relative_path)
                if not resolved.is_file() or LocalExecutionService._has_symlink(
                    resolved
                ):
                    return LocalExecutionService._workspace_failure()
                if resolved.read_bytes() != expected_contents:
                    return LocalExecutionService._workspace_failure()
        except (OSError, WorkspacePathError):
            return LocalExecutionService._workspace_failure()
        return Result.success(None)

    @staticmethod
    def _has_symlink(path: Path) -> bool:
        return any(component.is_symlink() for component in (*path.parents[::-1], path))

    @staticmethod
    def _workspace_failure() -> Result[None]:
        return Result.failure(
            [
                Issue(
                    code="LOCAL_EXECUTION_WORKSPACE_INVALID",
                    message="Preflighted workspace is missing or has changed.",
                    severity="error",
                    path="workspace",
                    source="local_execution_service",
                    hint="Create a new run and complete preflight again.",
                )
            ]
        )

    def _fail(
        self,
        run_id: str,
        issues: tuple[Issue, ...] | list[Issue],
    ) -> Result[RunRecord]:
        issue_list = list(issues)
        if not issue_list:
            issue_list.append(
                Issue(
                    code="LOCAL_EXECUTION_FAILED",
                    message="Local workflow execution failed.",
                    severity="error",
                    path="execution",
                    source="local_execution_service",
                )
            )
        reason_code = issue_list[0].code
        if reason_code == "LOCAL_RUN_PROCESS_FAILED":
            nested_reason = issue_list[0].context.get("reason_code")
            if nested_reason == "PROCESS_RUNNER_TIMEOUT":
                reason_code = nested_reason
        record_issue = Issue(
            code="RUN_EXECUTION_FAILED",
            message="Local workflow execution failed.",
            severity="error",
            path="execution",
            source="local_execution_service",
            hint="Review the persisted stdout and stderr logs.",
            context={"reason_code": reason_code},
        )

        for _attempt in range(2):
            try:
                current = self._run_service.get_run(run_id)
            except KeyError:
                return Result.failure(issue_list)
            if current.status is RunStatus.CANCELLED:
                return self._cancelled_result()
            if current.status.is_terminal:
                return Result.failure(issue_list)
            if current.status not in {RunStatus.QUEUED, RunStatus.RUNNING}:
                return Result.failure(issue_list)
            try:
                self._run_service.transition_run(
                    run_id,
                    RunStatus.FAILED,
                    stage="execution",
                    message="Local workflow execution failed.",
                    issue=record_issue,
                    context={"reason_code": reason_code},
                )
                break
            except (ConcurrentRunUpdateError, ValueError):
                continue
        return Result.failure(issue_list)

    @staticmethod
    def _cancelled_result() -> Result[RunRecord]:
        return Result.failure(
            [
                Issue(
                    code="LOCAL_EXECUTION_CANCELLED",
                    message="Run was cancelled during local execution.",
                    severity="error",
                    path="run_id",
                    source="local_execution_service",
                )
            ]
        )

    @staticmethod
    def _state_changed_result(current: RunRecord) -> Result[RunRecord]:
        return Result.failure(
            [
                Issue(
                    code="LOCAL_EXECUTION_STATE_CHANGED",
                    message="Run lifecycle state changed during local execution.",
                    severity="error",
                    path="run_id",
                    source="local_execution_service",
                    context={"current_status": current.status.value},
                )
            ]
        )
