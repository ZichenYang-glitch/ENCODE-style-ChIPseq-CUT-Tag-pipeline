"""Fail-closed local execution driver with pre-execution orchestration."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING

from encode_pipeline.platform.adapters import CommandSpec
from encode_pipeline.platform.planning import (
    PlanStatus,
    WorkspacePathError,
    WorkspacePathPolicy,
)
from encode_pipeline.platform.results import Issue, Result

if TYPE_CHECKING:
    from encode_pipeline.platform.planning import ExecutionPlan
    from encode_pipeline.services.command_builder import CommandBuilder
    from encode_pipeline.services.materialization import WorkspaceMaterializer
    from encode_pipeline.services.process_runner import ProcessResult, ProcessRunner
    from encode_pipeline.services.runs import RunService


class LocalRunDriver:
    """Fail-closed local execution driver with pre-execution orchestration.

    Validates the plan/run association, then for PENDING plans with a
    workspace_plan: derives a per-run workspace directory, materializes
    workspace files, builds the CommandSpec, executes a Snakemake dry-run,
    and returns the PLANNED ExecutionPlan on success. A PLANNED plan is run
    without dry-run flags while stdout and stderr are persisted incrementally.

    Does not transition run status.
    """

    def __init__(
        self,
        run_service: "RunService",
        materializer: "WorkspaceMaterializer",
        command_builder: "CommandBuilder",
        workspace_root: Path,
        *,
        process_runner: "ProcessRunner | None" = None,
    ) -> None:
        from encode_pipeline.services.command_builder import CommandBuilder
        from encode_pipeline.services.materialization import WorkspaceMaterializer
        from encode_pipeline.services.process_runner import ProcessRunner
        from encode_pipeline.services.runs import RunService

        if not isinstance(run_service, RunService):
            raise ValueError("LocalRunDriver requires a RunService instance")
        if not isinstance(materializer, WorkspaceMaterializer):
            raise ValueError("LocalRunDriver requires a WorkspaceMaterializer instance")
        if not isinstance(command_builder, CommandBuilder):
            raise ValueError("LocalRunDriver requires a CommandBuilder instance")
        if not isinstance(workspace_root, Path) or not workspace_root.is_absolute():
            raise ValueError("LocalRunDriver requires an absolute workspace_root Path")
        if process_runner is not None and not isinstance(process_runner, ProcessRunner):
            raise ValueError("LocalRunDriver requires a ProcessRunner instance or None")
        self._run_service = run_service
        self._materializer = materializer
        self._command_builder = command_builder
        self._workspace_root = workspace_root
        self._process_runner = (
            process_runner if process_runner is not None else ProcessRunner()
        )

    def run(self, run_id: str, plan: "ExecutionPlan") -> "Result[ExecutionPlan]":
        """Validate plan, materialize workspace, build command, and dry-run.

        For PENDING plans with a workspace_plan, materializes the workspace,
        builds the CommandSpec, executes a Snakemake dry-run, and returns the
        PLANNED ExecutionPlan on success. PLANNED plans execute through the
        controlled ProcessRunner boundary.
        """
        record = self._run_service.get_run(run_id)

        if record.run_id != plan.run_id or record.workflow_id != plan.workflow_id:
            issue = Issue(
                code="LOCAL_RUN_PLAN_MISMATCH",
                message="Plan run_id or workflow_id does not match the run record.",
                severity="error",
                path="plan",
                source="local_run_driver",
            )
            return self._refuse(issue, plan, record.run_id)

        if plan.status is PlanStatus.PENDING and plan.workspace_plan is not None:
            prepared = self._prepare(run_id, plan)
            if isinstance(prepared, Result):
                return prepared
            return Result.success(prepared)

        if plan.status is PlanStatus.PENDING:
            issue = Issue(
                code="LOCAL_RUN_MISSING_WORKSPACE_PLAN",
                message="Execution plan has no workspace plan.",
                severity="error",
                path="plan",
                source="local_run_driver",
            )
            return self._refuse(issue, plan, record.run_id)

        if not plan.can_execute:
            if plan.status is not PlanStatus.PLANNED:
                issue = Issue(
                    code="LOCAL_RUN_NOT_PLANNED",
                    message="Execution plan is not in the planned state.",
                    severity="error",
                    path="plan",
                    source="local_run_driver",
                )
            else:
                issue = Issue(
                    code="LOCAL_RUN_MISSING_COMMAND_SPEC",
                    message="Execution plan has no command spec.",
                    severity="error",
                    path="plan",
                    source="local_run_driver",
                )
            return self._refuse(issue, plan, record.run_id)

        # Dry-run flag conflict check
        if any(arg in ("-n", "--dry-run") for arg in plan.command_spec.argv):
            issue = Issue(
                code="LOCAL_RUN_DRY_RUN_FLAG_CONFLICT",
                message="CommandSpec argv already contains a dry-run flag.",
                severity="error",
                path="command_spec",
                source="local_run_driver",
            )
            return self._refuse(issue, plan, record.run_id)

        return self._execute(run_id, plan)

    def _refuse(
        self,
        issue: Issue,
        plan: "ExecutionPlan",
        run_id: str,
        *,
        additional_issues: Iterable[Issue] = (),
    ) -> "Result[ExecutionPlan]":
        """Record a runner_refused event and return a failure result."""
        self._run_service.add_event(
            run_id=run_id,
            event_type="runner_refused",
            message="Local execution refused.",
            status=None,
            context={
                "reason_code": issue.code,
                "plan_status": plan.status.value,
                "can_execute": plan.can_execute,
                "has_command_spec": plan.command_spec is not None,
            },
        )
        return Result.failure([issue] + list(additional_issues))

    def derive_workspace_dir(self, run_id: str) -> Path:
        """Return the absolute per-run workspace directory under workspace_root."""
        policy = WorkspacePathPolicy(base_dir=self._workspace_root)
        return policy.resolve(run_id)

    def _prepare(
        self,
        run_id: str,
        plan: "ExecutionPlan",
    ) -> "ExecutionPlan | Result[ExecutionPlan]":
        """Materialize workspace and build command spec for a PENDING plan.

        Returns the PLANNED ExecutionPlan on success, or a failure Result
        (already recording runner_refused) on any step failure.
        """
        # Step 1: Derive per-run workspace directory
        try:
            policy = WorkspacePathPolicy(base_dir=self._workspace_root)
            workspace_dir = policy.resolve(run_id)
        except WorkspacePathError:
            issue = Issue(
                code="LOCAL_RUN_WORKSPACE_DIR_INVALID",
                message="Could not derive a safe workspace directory for this run.",
                severity="error",
                path="run_id",
                source="local_run_driver",
            )
            return self._refuse(issue, plan, run_id)

        # Step 2: Materialize workspace files
        result = self._materializer.materialize(plan.workspace_plan, workspace_dir)
        if result.is_failure:
            return self._refuse(
                Issue(
                    code="LOCAL_RUN_MATERIALIZATION_FAILED",
                    message="Workspace materialization failed.",
                    severity="error",
                    path="plan",
                    source="local_run_driver",
                ),
                plan,
                run_id,
                additional_issues=result.issues,
            )

        self._run_service.add_event(
            run_id=run_id,
            event_type="workspace_materialized",
            message="Workspace materialized successfully.",
            status=None,
            context={
                "directory_count": len(plan.workspace_plan.directories),
                "file_count": len(plan.workspace_plan.files),
            },
        )

        # Step 3: Build command spec
        result = self._command_builder.build_command(plan, workspace_dir)
        if result.is_failure:
            return self._refuse(
                Issue(
                    code="LOCAL_RUN_COMMAND_BUILD_FAILED",
                    message="Command spec build failed.",
                    severity="error",
                    path="plan",
                    source="local_run_driver",
                ),
                plan,
                run_id,
                additional_issues=result.issues,
            )

        planned_plan = result.value
        self._run_service.add_event(
            run_id=run_id,
            event_type="command_built",
            message="Command spec built successfully.",
            status=None,
            context={
                "has_command_spec": True,
                "plan_status": planned_plan.status.value,
            },
        )

        # Dry-run flag conflict check
        original_argv = planned_plan.command_spec.argv
        if any(arg in ("-n", "--dry-run") for arg in original_argv):
            issue = Issue(
                code="LOCAL_RUN_DRY_RUN_FLAG_CONFLICT",
                message="CommandSpec argv already contains a dry-run flag.",
                severity="error",
                path="command_spec",
                source="local_run_driver",
            )
            return self._refuse(issue, planned_plan, run_id)

        # Build dry-run CommandSpec
        dry_run_spec = CommandSpec(
            argv=original_argv + ("-n",),
            cwd=planned_plan.command_spec.cwd,
            env=planned_plan.command_spec.env,
        )

        # Execute dry-run via ProcessRunner
        dry_run_result = self._process_runner.run(dry_run_spec)

        if dry_run_result.is_success:
            self._append_dry_run_logs(run_id, dry_run_result.value)

        if dry_run_result.is_success and dry_run_result.value.exit_code == 0:
            self._run_service.add_event(
                run_id=run_id,
                event_type="dry_run_completed",
                message="Snakemake dry-run completed successfully.",
                status=None,
                context={"exit_code": 0},
            )
            return planned_plan

        if dry_run_result.is_success:  # exit_code != 0
            self._run_service.add_event(
                run_id=run_id,
                event_type="dry_run_failed",
                message="Snakemake dry-run failed with a non-zero exit code.",
                status=None,
                context={
                    "reason_code": "LOCAL_RUN_DRY_RUN_FAILED",
                    "exit_code": dry_run_result.value.exit_code,
                    "issue_count": len(dry_run_result.issues),
                },
            )
            return self._refuse(
                Issue(
                    code="LOCAL_RUN_DRY_RUN_FAILED",
                    message="Snakemake dry-run failed.",
                    severity="error",
                    path="command_spec",
                    source="local_run_driver",
                ),
                planned_plan,
                run_id,
                additional_issues=dry_run_result.issues,
            )

        # ProcessRunner failure (timeout / not found / OSError)
        first_issue_code = (
            dry_run_result.issues[0].code if dry_run_result.issues else "UNKNOWN"
        )
        self._run_service.add_event(
            run_id=run_id,
            event_type="dry_run_failed",
            message="Snakemake dry-run could not be executed.",
            status=None,
            context={
                "reason_code": "LOCAL_RUN_PROCESS_FAILED",
                "process_issue_code": first_issue_code,
                "issue_count": len(dry_run_result.issues),
            },
        )
        return self._refuse(
            Issue(
                code="LOCAL_RUN_PROCESS_FAILED",
                message="ProcessRunner failed to execute the dry-run command.",
                severity="error",
                path="command_spec",
                source="local_run_driver",
            ),
            planned_plan,
            run_id,
            additional_issues=dry_run_result.issues,
        )

    def _append_dry_run_logs(
        self, run_id: str, process_result: "ProcessResult"
    ) -> None:
        """Append non-empty stdout/stderr from a dry-run ProcessResult to RunService logs."""
        if process_result.stdout:
            self._run_service.append_log(
                run_id, "stdout", process_result.stdout.splitlines()
            )
        if process_result.stderr:
            self._run_service.append_log(
                run_id, "stderr", process_result.stderr.splitlines()
            )

    def _execute(
        self,
        run_id: str,
        plan: "ExecutionPlan",
    ) -> "Result[ExecutionPlan]":
        """Execute one already-planned command and persist both output streams."""
        assert plan.command_spec is not None
        streamed: set[str] = set()

        def persist_chunk(stream_name: str, lines: tuple[str, ...]) -> None:
            if not lines:
                return
            self._run_service.append_log(run_id, stream_name, lines)
            streamed.add(stream_name)

        process_result = self._process_runner.run(
            plan.command_spec,
            output_callback=persist_chunk,
        )
        self._persist_execution_warnings(run_id, process_result.issues)
        if process_result.is_failure:
            reason_code = (
                process_result.issues[0].code
                if process_result.issues
                else "PROCESS_RUNNER_FAILED"
            )
            issue = Issue(
                code="LOCAL_RUN_PROCESS_FAILED",
                message="Workflow process could not be executed.",
                severity="error",
                path="command_spec",
                source="local_run_driver",
                context={"reason_code": reason_code},
            )
            return Result.failure([issue, *process_result.issues])

        completed = process_result.value
        # Test doubles and non-streaming ProcessRunner implementations can still
        # return captured output. Persist it once when no callback chunk for that
        # stream was observed.
        if completed.stdout and "stdout" not in streamed:
            self._run_service.append_log(
                run_id,
                "stdout",
                completed.stdout.splitlines(),
            )
        if completed.stderr and "stderr" not in streamed:
            self._run_service.append_log(
                run_id,
                "stderr",
                completed.stderr.splitlines(),
            )

        if completed.exit_code != 0:
            issue = Issue(
                code="LOCAL_RUN_EXECUTION_FAILED",
                message="Workflow process exited with a non-zero status.",
                severity="error",
                path="command_spec",
                source="local_run_driver",
                context={"exit_code": completed.exit_code},
            )
            return Result.failure([issue, *process_result.issues])

        return Result.success(plan, issues=process_result.issues)

    def _persist_execution_warnings(
        self,
        run_id: str,
        issues: Iterable[Issue],
    ) -> None:
        """Persist user-visible execution warnings before any failure return."""
        for process_issue in issues:
            if process_issue.code != "PROCESS_RUNNER_OUTPUT_LIMIT_EXCEEDED":
                continue
            self._run_service.add_event(
                run_id=run_id,
                event_type="execution_output_truncated",
                message="Workflow output exceeded the persisted log size limit.",
                status=None,
                stage="execution",
                context={"reason_code": process_issue.code},
                issue=process_issue,
            )
