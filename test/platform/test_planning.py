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
