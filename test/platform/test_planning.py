"""Tests for execution planning platform primitives."""

from copy import deepcopy
from dataclasses import FrozenInstanceError
from datetime import datetime, timezone

import pytest

from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus


def test_plan_status_values():
    assert PlanStatus.UNSUPPORTED.value == "unsupported"
    assert PlanStatus.PENDING.value == "pending"
    assert PlanStatus.PLANNED.value == "planned"


def test_plan_status_string_roundtrip():
    assert PlanStatus("unsupported") is PlanStatus.UNSUPPORTED
    assert PlanStatus("pending") is PlanStatus.PENDING
    assert PlanStatus("planned") is PlanStatus.PLANNED


def test_execution_plan_is_frozen():
    plan = ExecutionPlan(
        plan_id="plan-1",
        run_id="run-1",
        workflow_id="wf-1",
        status=PlanStatus.UNSUPPORTED,
        inputs_snapshot={"config": {"a": 1}},
        created_at=datetime.now(timezone.utc),
    )
    with pytest.raises(FrozenInstanceError):
        plan.status = PlanStatus.PLANNED


def test_execution_plan_can_execute_false_for_unsupported():
    plan = ExecutionPlan(
        plan_id="plan-1",
        run_id="run-1",
        workflow_id="wf-1",
        status=PlanStatus.UNSUPPORTED,
        inputs_snapshot={},
        created_at=datetime.now(timezone.utc),
    )
    assert plan.can_execute is False


def test_execution_plan_can_execute_true_only_when_planned_and_command_spec():
    from encode_pipeline.platform.adapters import CommandSpec

    unsupported_with_command = ExecutionPlan(
        plan_id="plan-1",
        run_id="run-1",
        workflow_id="wf-1",
        status=PlanStatus.UNSUPPORTED,
        inputs_snapshot={},
        command_spec=CommandSpec(argv=("echo", "hi")),
        created_at=datetime.now(timezone.utc),
    )
    assert unsupported_with_command.can_execute is False

    planned_without_command = ExecutionPlan(
        plan_id="plan-2",
        run_id="run-2",
        workflow_id="wf-1",
        status=PlanStatus.PLANNED,
        inputs_snapshot={},
        created_at=datetime.now(timezone.utc),
    )
    assert planned_without_command.can_execute is False

    planned_with_command = ExecutionPlan(
        plan_id="plan-3",
        run_id="run-3",
        workflow_id="wf-1",
        status=PlanStatus.PLANNED,
        inputs_snapshot={},
        command_spec=CommandSpec(argv=("echo", "hi")),
        created_at=datetime.now(timezone.utc),
    )
    assert planned_with_command.can_execute is True


def test_execution_plan_inputs_snapshot_is_defensive_copy():
    original = {"config": {"samples": ["s1"]}}
    plan = ExecutionPlan(
        plan_id="plan-1",
        run_id="run-1",
        workflow_id="wf-1",
        status=PlanStatus.UNSUPPORTED,
        inputs_snapshot=original,
        created_at=datetime.now(timezone.utc),
    )
    assert plan.inputs_snapshot == original
    assert plan.inputs_snapshot is not original
    original["config"]["samples"].append("s2")
    assert plan.inputs_snapshot == {"config": {"samples": ["s1"]}}


def test_execution_plan_default_slots_are_none():
    plan = ExecutionPlan(
        plan_id="plan-1",
        run_id="run-1",
        workflow_id="wf-1",
        status=PlanStatus.UNSUPPORTED,
        inputs_snapshot={},
        created_at=datetime.now(timezone.utc),
    )
    assert plan.dag_preview is None
    assert plan.workspace_plan is None
    assert plan.command_spec is None


def test_execution_plan_issues_defaults_to_empty_tuple():
    plan = ExecutionPlan(
        plan_id="plan-1",
        run_id="run-1",
        workflow_id="wf-1",
        status=PlanStatus.UNSUPPORTED,
        inputs_snapshot={},
        created_at=datetime.now(timezone.utc),
    )
    assert plan.issues == ()


def test_execution_plan_rejects_non_issue_issues():
    with pytest.raises(ValueError, match="ExecutionPlan issues must contain only Issue entries"):
        ExecutionPlan(
            plan_id="plan-1",
            run_id="run-1",
            workflow_id="wf-1",
            status=PlanStatus.UNSUPPORTED,
            inputs_snapshot={},
            created_at=datetime.now(timezone.utc),
            issues=("not-an-issue",),
        )


from pathlib import Path

from encode_pipeline.platform.planning import (
    WorkspaceBaseDirRelativeError,
    WorkspacePathAbsoluteError,
    WorkspacePathPolicy,
    WorkspacePathTraversalError,
    WorkspacePathInvalidError,
)


def test_workspace_path_policy_accepts_valid_relative_paths(tmp_path):
    base_dir = tmp_path.resolve()
    policy = WorkspacePathPolicy(base_dir=base_dir)

    assert policy.resolve("logs") == base_dir / "logs"
    assert policy.resolve("results/peaks") == base_dir / "results" / "peaks"


def test_workspace_path_policy_rejects_relative_base_dir(tmp_path):
    with pytest.raises(WorkspaceBaseDirRelativeError) as exc_info:
        WorkspacePathPolicy(base_dir=Path("relative/path"))
    assert exc_info.value.code == "WORKSPACE_BASE_DIR_RELATIVE"


def test_workspace_path_policy_rejects_absolute_planned_path(tmp_path):
    base_dir = tmp_path.resolve()
    policy = WorkspacePathPolicy(base_dir=base_dir)

    with pytest.raises(WorkspacePathAbsoluteError) as exc_info:
        policy.resolve("/tmp/logs")
    assert exc_info.value.code == "WORKSPACE_PATH_ABSOLUTE"


def test_workspace_path_policy_rejects_parent_traversal(tmp_path):
    base_dir = tmp_path.resolve()
    policy = WorkspacePathPolicy(base_dir=base_dir)

    for path in ("../logs", "results/../../escape", "a/../b"):
        with pytest.raises(WorkspacePathTraversalError) as exc_info:
            policy.resolve(path)
        assert exc_info.value.code == "WORKSPACE_PATH_TRAVERSAL"


def test_workspace_path_policy_rejects_dot_component(tmp_path):
    base_dir = tmp_path.resolve()
    policy = WorkspacePathPolicy(base_dir=base_dir)

    for path in ("./logs", "results/./peaks", "."):
        with pytest.raises(WorkspacePathTraversalError) as exc_info:
            policy.resolve(path)
        assert exc_info.value.code == "WORKSPACE_PATH_TRAVERSAL"


@pytest.mark.parametrize(
    "path",
    [
        "",
        "~/logs",
        "C:" + "\\logs",
        "logs with spaces",
        "logs\x00file",
    ],
)
def test_workspace_path_policy_rejects_invalid_paths(tmp_path, path):
    base_dir = tmp_path.resolve()
    policy = WorkspacePathPolicy(base_dir=base_dir)

    with pytest.raises(WorkspacePathInvalidError) as exc_info:
        policy.resolve(path)
    assert exc_info.value.code == "WORKSPACE_PATH_INVALID"
