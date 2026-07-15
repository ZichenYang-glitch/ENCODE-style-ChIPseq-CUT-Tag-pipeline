"""Fail-closed local execution driver with pre-execution orchestration."""

from __future__ import annotations

from collections.abc import Iterable
import os
from pathlib import Path
import stat
from typing import TYPE_CHECKING

from encode_pipeline.platform.adapters import CommandSpec
from encode_pipeline.platform.planning import (
    PlanStatus,
    WorkspacePathError,
    WorkspacePathPolicy,
)
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.process_runner import redact_bounded_literals

if TYPE_CHECKING:
    from encode_pipeline.platform.planning import ExecutionPlan
    from encode_pipeline.services.command_builder import CommandBuilder
    from encode_pipeline.services.materialization import WorkspaceMaterializer
    from encode_pipeline.services.process_runner import ProcessResult, ProcessRunner
    from encode_pipeline.services.runs import RunService


MANAGED_LOG_MAX_BYTES = 10_000_000


class LocalRunDriver:
    """Fail-closed local execution driver with pre-execution orchestration.

    Validates the plan/run association, then for PENDING plans with a
    workspace_plan: derives a per-run workspace directory, materializes
    workspace files, builds the CommandSpec, executes its optional preflight
    command, and returns the PLANNED ExecutionPlan on success. A PLANNED plan
    is run without preflight flags while stdout and stderr are persisted
    incrementally.

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
        builds the CommandSpec, executes its optional preflight command, and
        returns the PLANNED ExecutionPlan on success. PLANNED plans execute
        through the controlled ProcessRunner boundary.
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

        # Main commands must never smuggle legacy dry-run flags into execution.
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

        preflight_argv = planned_plan.command_spec.preflight_argv
        if preflight_argv is None:
            return planned_plan

        # Build the exact adapter/legacy-declared preflight CommandSpec.
        dry_run_spec = CommandSpec(
            argv=preflight_argv,
            cwd=planned_plan.command_spec.cwd,
            env=planned_plan.command_spec.env,
            redaction_values=planned_plan.command_spec.redaction_values,
        )

        # Execute dry-run via ProcessRunner
        dry_run_result = self._process_runner.run(dry_run_spec)

        if dry_run_result.is_success:
            self._append_dry_run_logs(run_id, dry_run_result.value)

        managed_result = self._ingest_managed_logs(
            run_id=run_id,
            workspace_dir=workspace_dir,
            managed_logs=planned_plan.command_spec.preflight_managed_logs,
            redaction_values=planned_plan.command_spec.redaction_values,
            phase="preflight",
        )
        if managed_result.is_failure:
            self._record_preflight_failed(
                run_id,
                kind=planned_plan.command_spec.preflight_kind,
                failure="managed_log",
                context={
                    "reason_code": "LOCAL_RUN_MANAGED_LOG_FAILED",
                    "issue_count": len(managed_result.issues),
                },
            )
            return self._refuse(
                managed_result.issues[0],
                planned_plan,
                run_id,
                additional_issues=managed_result.issues[1:],
            )

        if dry_run_result.is_success and dry_run_result.value.exit_code == 0:
            self._record_preflight_completed(
                run_id,
                kind=planned_plan.command_spec.preflight_kind,
            )
            return planned_plan

        if dry_run_result.is_success:  # exit_code != 0
            kind = planned_plan.command_spec.preflight_kind
            self._record_preflight_failed(
                run_id,
                kind=kind,
                failure="nonzero",
                context={
                    "reason_code": (
                        "LOCAL_RUN_DRY_RUN_FAILED"
                        if kind == "dry_run"
                        else "LOCAL_RUN_PREFLIGHT_FAILED"
                    ),
                    "exit_code": dry_run_result.value.exit_code,
                    "issue_count": len(dry_run_result.issues),
                },
            )
            return self._refuse(
                Issue(
                    code=(
                        "LOCAL_RUN_DRY_RUN_FAILED"
                        if kind == "dry_run"
                        else "LOCAL_RUN_PREFLIGHT_FAILED"
                    ),
                    message=(
                        "Snakemake dry-run failed."
                        if kind == "dry_run"
                        else "Workflow configuration preflight failed."
                    ),
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
        kind = planned_plan.command_spec.preflight_kind
        self._record_preflight_failed(
            run_id,
            kind=kind,
            failure="process",
            context={
                "reason_code": "LOCAL_RUN_PROCESS_FAILED",
                "process_issue_code": first_issue_code,
                "issue_count": len(dry_run_result.issues),
            },
        )
        return self._refuse(
            Issue(
                code="LOCAL_RUN_PROCESS_FAILED",
                message=(
                    "ProcessRunner failed to execute the dry-run command."
                    if kind == "dry_run"
                    else "ProcessRunner failed to execute the preflight command."
                ),
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

    def _record_preflight_completed(self, run_id: str, *, kind: str) -> None:
        if kind == "dry_run":
            event_type = "dry_run_completed"
            message = "Snakemake dry-run completed successfully."
        else:
            event_type = "preflight_completed"
            message = "Workflow configuration preflight completed successfully."
        self._run_service.add_event(
            run_id=run_id,
            event_type=event_type,
            message=message,
            status=None,
            context={"exit_code": 0},
        )

    def _record_preflight_failed(
        self,
        run_id: str,
        *,
        kind: str,
        failure: str,
        context: dict[str, object],
    ) -> None:
        if kind == "dry_run":
            event_type = "dry_run_failed"
            messages = {
                "nonzero": "Snakemake dry-run failed with a non-zero exit code.",
                "process": "Snakemake dry-run could not be executed.",
                "managed_log": "Snakemake dry-run log could not be persisted.",
            }
        else:
            event_type = "preflight_failed"
            messages = {
                "nonzero": (
                    "Workflow configuration preflight failed with a non-zero exit code."
                ),
                "process": "Workflow configuration preflight could not be executed.",
                "managed_log": (
                    "Workflow configuration preflight log could not be persisted."
                ),
            }
        self._run_service.add_event(
            run_id=run_id,
            event_type=event_type,
            message=messages[failure],
            status=None,
            context=context,
        )

    def _ingest_managed_logs(
        self,
        *,
        run_id: str,
        workspace_dir: Path,
        managed_logs: tuple[tuple[str, str], ...],
        redaction_values: tuple[str, ...],
        phase: str,
    ) -> Result[None]:
        """Sanitize managed regular files in place and persist their lines."""
        for stream_name, path_value in managed_logs:
            try:
                sanitized = _sanitize_managed_log(
                    workspace_dir=workspace_dir,
                    path_value=path_value,
                    redaction_values=redaction_values,
                )
                lines = tuple(sanitized.splitlines())
                if lines:
                    self._run_service.append_log(run_id, stream_name, lines)
            except Exception:
                return Result.failure(
                    [
                        Issue(
                            code="LOCAL_RUN_MANAGED_LOG_FAILED",
                            message="A managed workflow log could not be persisted.",
                            severity="error",
                            path="command_spec.managed_logs",
                            source="local_run_driver",
                            context={"phase": phase},
                        )
                    ]
                )
        return Result.success(None)

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
        managed_result = self._ingest_managed_logs(
            run_id=run_id,
            workspace_dir=self.derive_workspace_dir(run_id),
            managed_logs=plan.command_spec.execution_managed_logs,
            redaction_values=plan.command_spec.redaction_values,
            phase="execution",
        )
        if managed_result.is_failure:
            return Result.failure([*managed_result.issues, *process_result.issues])
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


def _sanitize_managed_log(
    *,
    workspace_dir: Path,
    path_value: str,
    redaction_values: tuple[str, ...],
) -> str:
    """Read and rewrite one bounded workspace log through no-follow descriptors."""
    required_flags = ("O_DIRECTORY", "O_NOFOLLOW", "O_CLOEXEC", "O_NONBLOCK")
    if not all(hasattr(os, name) for name in required_flags):
        raise OSError("required no-follow flags are unavailable")
    if not isinstance(workspace_dir, Path) or not workspace_dir.is_absolute():
        raise ValueError("workspace is invalid")
    path = Path(path_value)
    try:
        relative = path.relative_to(workspace_dir)
    except ValueError as exc:
        raise ValueError("managed log escaped workspace") from exc
    if not relative.parts or any(part in {"", ".", ".."} for part in relative.parts):
        raise ValueError("managed log path is invalid")

    directory_flags = os.O_RDONLY | os.O_DIRECTORY | os.O_NOFOLLOW | os.O_CLOEXEC
    file_flags = os.O_RDWR | os.O_NOFOLLOW | os.O_NONBLOCK | os.O_CLOEXEC
    descriptors: list[int] = []
    try:
        root_fd = os.open(workspace_dir, directory_flags)
        descriptors.append(root_fd)
        root_info = os.fstat(root_fd)
        if not stat.S_ISDIR(root_info.st_mode):
            raise OSError("workspace is not a directory")

        parent_fd = root_fd
        for component in relative.parts[:-1]:
            listed = os.stat(component, dir_fd=parent_fd, follow_symlinks=False)
            child_fd = os.open(component, directory_flags, dir_fd=parent_fd)
            descriptors.append(child_fd)
            opened = os.fstat(child_fd)
            if (
                not stat.S_ISDIR(listed.st_mode)
                or not stat.S_ISDIR(opened.st_mode)
                or _inode_identity(listed) != _inode_identity(opened)
            ):
                raise OSError("managed log directory changed")
            parent_fd = child_fd

        name = relative.parts[-1]
        listed = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if not _bounded_regular_log(listed):
            raise OSError("managed log is not a bounded regular file")
        log_fd = os.open(name, file_flags, dir_fd=parent_fd)
        descriptors.append(log_fd)
        before = os.fstat(log_fd)
        if not _bounded_regular_log(before) or _stable_file_observation(
            listed
        ) != _stable_file_observation(before):
            raise OSError("managed log changed while opening")

        content = _read_bounded(log_fd)
        after_read = os.fstat(log_fd)
        listed_after_read = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if _stable_file_observation(before) != _stable_file_observation(
            after_read
        ) or _stable_file_observation(after_read) != _stable_file_observation(
            listed_after_read
        ):
            raise OSError("managed log changed while reading")

        sanitized = redact_bounded_literals(
            content.decode("utf-8", errors="replace"),
            redaction_values,
        )
        sanitized_bytes = sanitized.encode("utf-8")
        if len(sanitized_bytes) > MANAGED_LOG_MAX_BYTES:
            raise OSError("redacted managed log exceeds the bound")

        os.ftruncate(log_fd, 0)
        os.lseek(log_fd, 0, os.SEEK_SET)
        _write_all(log_fd, sanitized_bytes)
        os.ftruncate(log_fd, len(sanitized_bytes))
        os.fsync(log_fd)

        final = os.fstat(log_fd)
        listed_final = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if (
            not stat.S_ISREG(final.st_mode)
            or final.st_nlink != 1
            or final.st_size != len(sanitized_bytes)
            or _inode_identity(final) != _inode_identity(listed_final)
        ):
            raise OSError("managed log changed while rewriting")
        return sanitized
    finally:
        for descriptor in reversed(descriptors):
            try:
                os.close(descriptor)
            except OSError:
                pass


def _bounded_regular_log(info: os.stat_result) -> bool:
    return (
        stat.S_ISREG(info.st_mode)
        and info.st_nlink == 1
        and 0 <= info.st_size <= MANAGED_LOG_MAX_BYTES
    )


def _inode_identity(info: os.stat_result) -> tuple[int, int]:
    return info.st_dev, info.st_ino


def _stable_file_observation(info: os.stat_result) -> tuple[int, ...]:
    return (
        info.st_dev,
        info.st_ino,
        info.st_mode,
        info.st_nlink,
        info.st_size,
        info.st_mtime_ns,
        info.st_ctime_ns,
    )


def _read_bounded(descriptor: int) -> bytes:
    chunks: list[bytes] = []
    total = 0
    while True:
        chunk = os.read(
            descriptor,
            min(64 * 1024, MANAGED_LOG_MAX_BYTES - total + 1),
        )
        if not chunk:
            return b"".join(chunks)
        total += len(chunk)
        if total > MANAGED_LOG_MAX_BYTES:
            raise OSError("managed log exceeds the bound")
        chunks.append(chunk)


def _write_all(descriptor: int, content: bytes) -> None:
    offset = 0
    while offset < len(content):
        written = os.write(descriptor, content[offset:])
        if written <= 0:
            raise OSError("managed log rewrite made no progress")
        offset += written
