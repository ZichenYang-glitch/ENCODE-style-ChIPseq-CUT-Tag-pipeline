"""Fail-closed local execution driver with pre-execution orchestration."""

from __future__ import annotations

from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING

from encode_pipeline.platform.planning import PlanStatus, WorkspacePathError, WorkspacePathPolicy
from encode_pipeline.platform.results import Issue, Result

if TYPE_CHECKING:
    from encode_pipeline.platform.adapters import CommandSpec
    from encode_pipeline.platform.planning import ExecutionPlan
    from encode_pipeline.platform.runs import RunRecord
    from encode_pipeline.services.command_builder import CommandBuilder
    from encode_pipeline.services.materialization import WorkspaceMaterializer
    from encode_pipeline.services.process_runner import ProcessRunner
    from encode_pipeline.services.runs import RunService


class LocalRunDriver:
    """Fail-closed local execution driver with pre-execution orchestration.

    Validates the plan/run association, then for PENDING plans with a
    workspace_plan: derives a per-run workspace directory, materializes
    workspace files, builds the CommandSpec, records progress events,
    and refuses execution with LOCAL_RUN_NOT_IMPLEMENTED.

    Does not invoke subprocess, Snakemake, or transition run status.
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

    def run(self, run_id: str, plan: "ExecutionPlan") -> "Result[RunRecord]":
        """Validate plan and refuse execution.

        For PENDING plans with a workspace_plan, materializes the workspace
        and builds the command spec before refusing.
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
            issue = Issue(
                code="LOCAL_RUN_NOT_IMPLEMENTED",
                message="Local execution is not implemented yet.",
                severity="error",
                path="plan",
                source="local_run_driver",
            )
            return self._refuse(issue, prepared, record.run_id)

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

        issue = Issue(
            code="LOCAL_RUN_NOT_IMPLEMENTED",
            message="Local execution is not implemented yet.",
            severity="error",
            path="plan",
            source="local_run_driver",
        )
        return self._refuse(issue, plan, record.run_id)

    def _refuse(
        self,
        issue: Issue,
        plan: "ExecutionPlan",
        run_id: str,
        *,
        additional_issues: Iterable[Issue] = (),
    ) -> "Result[RunRecord]":
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

    def _prepare(
        self,
        run_id: str,
        plan: "ExecutionPlan",
    ) -> "ExecutionPlan | Result[RunRecord]":
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

        return planned_plan
