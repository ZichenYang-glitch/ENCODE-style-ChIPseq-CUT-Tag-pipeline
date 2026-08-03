"""Local preflight orchestrator for workflow runs."""

from __future__ import annotations

from collections.abc import Mapping
from typing import TYPE_CHECKING

from encode_pipeline.platform.adapters import (
    ReferenceProfileBindingAdapter,
    WorkflowAdapter,
    WorkflowInputs,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.reference_profiles import BoundWorkflowReference
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunRecord, RunStatus
from encode_pipeline.services.run_repositories import ConcurrentRunUpdateError

if TYPE_CHECKING:
    from encode_pipeline.services.local_run_driver import LocalRunDriver
    from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
    from encode_pipeline.services.runs import RunService
    from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


class LocalPreflightService:
    """Compose planning, workspace authoring, and dry-run into a preflight path.

    The service is synchronous so it can be unit-tested without HTTP. FastAPI
    routes should call ``preflight()`` directly for the full
    ``CREATED -> VALIDATING -> PLANNED/FAILED`` path, or schedule
    ``run_preflight()`` via ``BackgroundTasks`` after they have already
    transitioned the run to ``VALIDATING``.
    """

    def __init__(
        self,
        run_service: "RunService",
        execution_planner: "ExecutionPlanner",
        workspace_planner: "WorkspacePlanner",
        local_run_driver: "LocalRunDriver",
        build_identity_provider: "WorkflowBuildIdentityProvider",
        reference_profile_resolver=None,
    ) -> None:
        if (
            getattr(build_identity_provider, "registry", None)
            is not run_service.registry
        ):
            raise ValueError(
                "build_identity_provider registry must match run_service registry"
            )
        self._run_service = run_service
        self._execution_planner = execution_planner
        self._workspace_planner = workspace_planner
        self._local_run_driver = local_run_driver
        self._build_identity_provider = build_identity_provider
        if reference_profile_resolver is not None and not callable(
            getattr(reference_profile_resolver, "resolve_run", None)
        ):
            raise ValueError("reference_profile_resolver is invalid")
        self._reference_profile_resolver = reference_profile_resolver

    def preflight(self, run_id: str) -> Result[RunRecord]:
        """Run the full local preflight path for *run_id*.

        Accepts runs only in ``CREATED`` and transitions them to
        ``VALIDATING`` before executing. Any other status is treated as a
        duplicate trigger.
        """
        try:
            record = self._run_service.get_run(run_id)
        except KeyError:
            return Result.failure(
                [
                    Issue(
                        code="PREFLIGHT_RUN_NOT_FOUND",
                        message="Run was not found.",
                        severity="error",
                        path="run_id",
                        source="preflight_service",
                    )
                ]
            )

        if record.status is not RunStatus.CREATED:
            return Result.failure(
                [
                    Issue(
                        code="PREFLIGHT_ALREADY_TRIGGERED",
                        message="Preflight has already been triggered for this run.",
                        severity="error",
                        path="run_id",
                        source="preflight_service",
                        context={"current_status": record.status.value},
                    )
                ]
            )

        try:
            self._run_service.transition_run(
                run_id,
                RunStatus.VALIDATING,
                stage="preflight",
                message="Local preflight started.",
            )
        except (ConcurrentRunUpdateError, ValueError):
            current = self._run_service.get_run(run_id)
            return Result.failure(
                [
                    Issue(
                        code="PREFLIGHT_ALREADY_TRIGGERED",
                        message="Preflight has already been triggered for this run.",
                        severity="error",
                        path="run_id",
                        source="preflight_service",
                        context={"current_status": current.status.value},
                    )
                ]
            )
        return self._run_preflight(run_id)

    def run_preflight(self, run_id: str) -> Result[RunRecord]:
        """Background-worker entry point.

        The caller must already have transitioned the run to ``VALIDATING``.
        """
        try:
            record = self._run_service.get_run(run_id)
        except KeyError:
            return Result.failure(
                [
                    Issue(
                        code="PREFLIGHT_RUN_NOT_FOUND",
                        message="Run was not found.",
                        severity="error",
                        path="run_id",
                        source="preflight_service",
                    )
                ]
            )

        if record.status is RunStatus.CANCELLED:
            return self._cancelled_result()

        if record.status is not RunStatus.VALIDATING:
            return Result.failure(
                [
                    Issue(
                        code="PREFLIGHT_ALREADY_TRIGGERED",
                        message="Preflight has already been triggered for this run.",
                        severity="error",
                        path="run_id",
                        source="preflight_service",
                        context={"current_status": record.status.value},
                    )
                ]
            )

        return self._run_preflight(run_id)

    def _run_preflight(self, run_id: str) -> Result[RunRecord]:
        try:
            current = self._run_service.get_run(run_id)
            build_before = self._capture_current_build(current)
            if build_before.is_failure:
                return self._fail(run_id, build_before.issues)

            plan_result = self._execution_planner.plan_run(run_id)
            if plan_result.is_failure:
                return self._fail(run_id, plan_result.issues)

            base_plan = plan_result.value
            workspace_dir = self._local_run_driver.derive_workspace_dir(run_id)
            workspace_result = self._workspace_planner.plan_workspace(
                base_plan,
                base_dir=workspace_dir,
                require_reference_enabled=True,
            )
            if workspace_result.is_failure:
                return self._fail(run_id, workspace_result.issues)

            prepared_plan = workspace_result.value
            run_result = self._local_run_driver.run(run_id, prepared_plan)
            if run_result.is_failure:
                return self._fail(run_id, run_result.issues)

            final_plan = run_result.value
            build_after = self._capture_current_build(self._run_service.get_run(run_id))
            if build_after.is_failure:
                return self._fail(run_id, build_after.issues)
            if not build_before.value.matches(build_after.value):
                return self._fail(
                    run_id,
                    [
                        Issue(
                            code="PREFLIGHT_WORKFLOW_BUILD_CHANGED",
                            message="Workflow source changed during preflight.",
                            severity="error",
                            source="preflight_service",
                            path="workflow",
                        )
                    ],
                )

            current = self._run_service.get_run(run_id)
            if current.status is RunStatus.CANCELLED:
                return self._cancelled_result()

            preflight_kind = final_plan.command_spec.preflight_kind
            if preflight_kind == "dry_run":
                completion_message = "Local preflight completed; dry-run succeeded."
                completion_context = {
                    "plan_status": final_plan.status.value,
                    "has_command_spec": final_plan.command_spec is not None,
                    "reason_code": "PREFLIGHT_COMPLETED",
                }
            else:
                completion_message = (
                    "Local workflow configuration preflight completed successfully."
                )
                completion_context = {
                    "plan_status": final_plan.status.value,
                    "has_command_spec": final_plan.command_spec is not None,
                    "preflight_kind": "configuration",
                    "reason_code": "PREFLIGHT_CONFIGURATION_COMPLETED",
                }
            updated = self._run_service.complete_preflight(
                run_id,
                build_after.value,
                stage="preflight",
                message=completion_message,
                context=completion_context,
            )
            return Result.success(updated)
        except Exception:
            issue = Issue(
                code="PREFLIGHT_UNEXPECTED_ERROR",
                message="An unexpected error occurred during preflight.",
                severity="error",
                source="preflight_service",
                path="preflight",
            )
            return self._fail(run_id, [issue])

    def _capture_current_build(
        self,
        record: RunRecord,
    ) -> Result[WorkflowBuildIdentity]:
        try:
            adapter = self._run_service.registry.get(record.workflow_id)
        except (KeyError, ValueError):
            return self._reference_build_failure()
        if not isinstance(adapter, ReferenceProfileBindingAdapter):
            return self._build_identity_provider.capture_executable(record.workflow_id)
        resolver = self._reference_profile_resolver
        if resolver is None:
            return self._reference_build_failure()
        config = record.inputs.get("config", {})
        options = record.inputs.get("options", {})
        if not isinstance(config, Mapping) or not isinstance(options, Mapping):
            return self._reference_build_failure()
        try:
            resolved = resolver.resolve_run(
                record.run_id,
                record.workflow_id,
                WorkflowInputs(
                    config=config,
                    samples=record.inputs.get("samples"),
                    options=options,
                ),
                require_enabled=True,
            )
        except Exception:
            return self._reference_build_failure()
        if (
            not isinstance(resolved, Result)
            or resolved.is_failure
            or not isinstance(resolved.value, BoundWorkflowReference)
            or not isinstance(resolved.value.adapter, WorkflowAdapter)
        ):
            return self._reference_build_failure()
        return self._build_identity_provider.capture_resolved_executable(
            resolved.value.adapter
        )

    @staticmethod
    def _reference_build_failure() -> Result[WorkflowBuildIdentity]:
        return Result.failure(
            [
                Issue(
                    code="REFERENCE_PROFILE_BINDING_INVALID",
                    message=(
                        "The selected Reference Profile binding could not be verified."
                    ),
                    severity="error",
                    source="reference_profile",
                    path="reference_profile_revision_id",
                )
            ]
        )

    def _fail(
        self,
        run_id: str,
        issues: tuple[Issue, ...] | list[Issue],
    ) -> Result[RunRecord]:
        issue_list = list(issues)
        reason_code = issue_list[0].code if issue_list else "PREFLIGHT_FAILED"
        current = self._run_service.get_run(run_id)

        if current.status is RunStatus.CANCELLED:
            return self._cancelled_result()

        record_issue = Issue(
            code="PREFLIGHT_FAILED",
            message="Local preflight failed.",
            severity="error",
            path="preflight",
            source="preflight_service",
            context={"reason_code": reason_code},
        )
        self._run_service.transition_run(
            run_id,
            RunStatus.FAILED,
            stage="preflight",
            message="Local preflight failed.",
            issue=record_issue,
            context={"reason_code": reason_code},
        )
        return Result.failure(issue_list)

    def _cancelled_result(self) -> Result[RunRecord]:
        return Result.failure(
            [
                Issue(
                    code="PREFLIGHT_CANCELLED",
                    message="Preflight was cancelled before it completed.",
                    severity="error",
                    path="run_id",
                    source="preflight_service",
                )
            ]
        )
