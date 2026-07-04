"""Execution planning primitives for the workflow platform."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
from pathlib import Path
from typing import Any, Mapping

from encode_pipeline.platform.adapters import CommandSpec, DagPreview, WorkspacePlan
from encode_pipeline.platform.results import Issue


class PlanStatus(str, Enum):
    """Lifecycle status for an execution plan."""

    UNSUPPORTED = "unsupported"
    PENDING = "pending"
    PLANNED = "planned"


@dataclass(frozen=True)
class ExecutionPlan:
    """Immutable envelope describing the conceptual plan for a workflow run.

    PR106 produces only ``PlanStatus.UNSUPPORTED`` plans with every conceptual
    planning slot set to ``None``. Future PRs will fill ``dag_preview``,
    ``workspace_plan``, and ``command_spec`` as adapter/workspace planning
    boundaries are implemented.
    """

    plan_id: str
    run_id: str
    workflow_id: str
    status: PlanStatus
    inputs_snapshot: Mapping[str, Any]
    dag_preview: DagPreview | None = None
    workspace_plan: WorkspacePlan | None = None
    command_spec: CommandSpec | None = None
    created_at: datetime = field(default_factory=lambda: datetime.now(timezone.utc))
    issues: tuple[Issue, ...] = ()

    def __post_init__(self) -> None:
        object.__setattr__(self, "status", _normalize_status(self.status))
        object.__setattr__(self, "inputs_snapshot", _copy_inputs(self.inputs_snapshot))
        object.__setattr__(self, "issues", _normalize_issues(self.issues))

    @property
    def can_execute(self) -> bool:
        """Return True only when the plan is planned and has a command spec."""
        return self.status is PlanStatus.PLANNED and self.command_spec is not None


def _normalize_status(value: PlanStatus | str) -> PlanStatus:
    if isinstance(value, PlanStatus):
        return value
    if isinstance(value, str):
        try:
            return PlanStatus(value)
        except ValueError as exc:
            raise ValueError(f"Invalid plan status: {value!r}") from exc
    raise ValueError(f"Invalid plan status: {value!r}")


def _copy_inputs(inputs: Mapping[str, Any]) -> dict[str, Any]:
    if not isinstance(inputs, Mapping):
        raise ValueError("ExecutionPlan inputs_snapshot must be a mapping")
    return deepcopy(dict(inputs))


def _normalize_issues(issues: tuple[Issue, ...]) -> tuple[Issue, ...]:
    issue_tuple = tuple(issues)
    for issue in issue_tuple:
        if not isinstance(issue, Issue):
            raise ValueError("ExecutionPlan issues must contain only Issue entries")
    return issue_tuple


class WorkspacePathError(ValueError):
    """Base class for workspace path policy violations."""

    code: str = "WORKSPACE_PATH_ERROR"


class WorkspaceBaseDirRelativeError(WorkspacePathError):
    """base_dir is not an absolute path."""

    code = "WORKSPACE_BASE_DIR_RELATIVE"


class WorkspacePathAbsoluteError(WorkspacePathError):
    """Planned path is absolute."""

    code = "WORKSPACE_PATH_ABSOLUTE"


class WorkspacePathTraversalError(WorkspacePathError):
    """Planned path contains '.' or '..' components."""

    code = "WORKSPACE_PATH_TRAVERSAL"


class WorkspacePathInvalidError(WorkspacePathError):
    """Planned path is otherwise invalid (empty, spaces, NUL, etc.)."""

    code = "WORKSPACE_PATH_INVALID"


class WorkspacePathEscapeError(WorkspacePathError):
    """Resolved path would leave base_dir (defense in depth)."""

    code = "WORKSPACE_PATH_ESCAPE"


@dataclass(frozen=True)
class WorkspacePathPolicy:
    """Lexical-only path safety policy for planned workspace paths."""

    base_dir: Path

    def __post_init__(self) -> None:
        base_dir = self.base_dir
        if not isinstance(base_dir, Path):
            raise WorkspaceBaseDirRelativeError("base_dir must be a Path")
        if not base_dir.is_absolute():
            raise WorkspaceBaseDirRelativeError(
                f"base_dir must be absolute: {base_dir}"
            )
        object.__setattr__(self, "base_dir", base_dir)

    def resolve(self, relative_path: str) -> Path:
        """Lexically join relative_path under base_dir and return the result.

        Does not create directories, does not call os.stat/lstat, and does not
        follow symlinks. Raises WorkspacePathError subclasses for policy violations.
        """
        if not isinstance(relative_path, str):
            raise WorkspacePathInvalidError("relative_path must be a string")

        if not relative_path:
            raise WorkspacePathInvalidError("relative_path must be non-empty")

        if "\x00" in relative_path or " " in relative_path:
            raise WorkspacePathInvalidError(
                "relative_path must not contain NUL or spaces"
            )

        if relative_path.startswith("/"):
            raise WorkspacePathAbsoluteError(
                f"relative_path must not be absolute: {relative_path}"
            )

        if relative_path.startswith("~"):
            raise WorkspacePathInvalidError(
                f"relative_path must not use home expansion: {relative_path}"
            )

        if len(relative_path) >= 2 and relative_path[1] == ":":
            raise WorkspacePathInvalidError(
                f"relative_path must not be a Windows drive path: {relative_path}"
            )

        for part in relative_path.split("/"):
            if part == "." or part == "..":
                raise WorkspacePathTraversalError(
                    f"relative_path must not contain '.' or '..' components: {relative_path}"
                )
            if not part:
                raise WorkspacePathInvalidError(
                    f"relative_path must not contain empty components: {relative_path}"
                )

        resolved = self.base_dir / relative_path
        try:
            resolved.relative_to(self.base_dir)
        except ValueError as exc:
            raise WorkspacePathEscapeError(
                f"relative_path escapes base_dir: {relative_path}"
            ) from exc

        return resolved
