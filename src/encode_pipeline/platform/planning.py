"""Execution planning primitives for the workflow platform."""

from __future__ import annotations

from copy import deepcopy
from dataclasses import dataclass, field
from datetime import datetime, timezone
from enum import Enum
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
