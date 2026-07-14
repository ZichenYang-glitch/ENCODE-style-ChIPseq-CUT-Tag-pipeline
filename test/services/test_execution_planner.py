"""Tests for the ExecutionPlanner service boundary."""

from pathlib import Path

import pytest

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import (
    DagPreview,
    WorkspacePlan,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
)
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.planning import ExecutionPlanner, WorkspacePlanner
from encode_pipeline.services.runs import RunService


@pytest.fixture
def fake_adapter():
    class _FakeAdapter:
        metadata = WorkflowMetadata(
            workflow_id="stub",
            name="Stub Workflow",
            version="0.0.1",
        )
        capabilities = WorkflowCapabilities(supports=("workspace_plan",))

        def schema(self) -> WorkflowSchema:
            return WorkflowSchema()

        def validate(self, inputs: WorkflowInputs) -> Result[object]:
            return Result.success(None)

        def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
            raise AssertionError("ExecutionPlanner must not call preview_dag")

        def plan_workspace(
            self,
            inputs: WorkflowInputs,
            workspace: str,
        ) -> Result[WorkspacePlan]:
            return Result.success(
                WorkspacePlan(directories=("logs", "results"), files=())
            )

        def build_command(self, plan: WorkspacePlan) -> Result[object]:
            raise AssertionError("ExecutionPlanner must not call build_command")

        def extract_artifacts(self, inputs, workspace):
            raise AssertionError("ExecutionPlanner must not call extract_artifacts")

    return _FakeAdapter()


@pytest.fixture
def run_service(fake_adapter):
    return RunService(registry=WorkflowRegistry(adapters=[fake_adapter]))


@pytest.fixture
def planner(run_service):
    from encode_pipeline.services.planning import ExecutionPlanner

    return ExecutionPlanner(run_service=run_service)


def test_plan_run_unknown_run_returns_failure():
    from encode_pipeline.services.planning import ExecutionPlanner

    run_service = RunService(registry=WorkflowRegistry(adapters=[]))
    planner = ExecutionPlanner(run_service=run_service)
    result = planner.plan_run("does-not-exist")

    assert result.is_failure is True
    assert result.value is None
    assert len(result.issues) == 1
    issue = result.issues[0]
    assert issue.code == "EXECUTION_RUN_NOT_FOUND"
    assert issue.message == "Run not found."
    assert issue.severity.value == "error"
    assert issue.path == "does-not-exist"
    assert issue.source == "execution_planner"


def test_plan_run_valid_run_returns_unsupported_plan(planner, run_service):
    inputs = WorkflowInputs(config={"genome": "hg38"}, samples=None, options={})
    record = run_service.create_run("stub", inputs)

    result = planner.plan_run(record.run_id)

    assert result.is_success is True
    plan = result.value
    assert plan is not None
    assert plan.run_id == record.run_id
    assert plan.workflow_id == record.workflow_id
    assert plan.status.value == "unsupported"
    assert plan.dag_preview is None
    assert plan.workspace_plan is None
    assert plan.command_spec is None
    assert plan.can_execute is False
    assert len(result.issues) == 1
    issue = result.issues[0]
    assert issue.code == "EXECUTION_PLANNING_UNSUPPORTED"
    assert issue.message == "Execution planning is not supported yet."
    assert issue.severity.value == "info"
    assert issue.path == "execution_plan"
    assert issue.source == "execution_planner"


def test_plan_run_does_not_change_run_status(planner, run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)
    original_status = record.status
    original_updated_at = record.updated_at

    planner.plan_run(record.run_id)

    after = run_service.get_run(record.run_id)
    assert after.status == original_status
    assert after.updated_at == original_updated_at


def test_plan_run_does_not_add_events(planner, run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)
    before_count = len(run_service.list_events(record.run_id))

    planner.plan_run(record.run_id)

    after_count = len(run_service.list_events(record.run_id))
    assert after_count == before_count


def test_plan_run_does_not_append_logs(planner, run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)

    planner.plan_run(record.run_id)

    chunks = run_service.list_logs(record.run_id, "stdout")
    assert chunks == ()


def test_plan_run_does_not_record_artifacts(planner, run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)

    planner.plan_run(record.run_id)

    artifacts = run_service.list_artifacts(record.run_id)
    assert artifacts == ()


def test_plan_run_does_not_call_adapter_planning_methods(planner, run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)

    result = planner.plan_run(record.run_id)

    assert result.is_success is True
    plan = result.value
    assert plan.status.value == "unsupported"
    assert plan.dag_preview is None
    assert plan.workspace_plan is None
    assert plan.command_spec is None


def test_plan_run_defensively_copies_inputs(planner, run_service):
    inputs = WorkflowInputs(
        config={"samples": ["s1"]},
        samples=None,
        options={},
    )
    record = run_service.create_run("stub", inputs)

    result = planner.plan_run(record.run_id)
    plan = result.value

    original_inputs = run_service.get_run(record.run_id).inputs
    original_inputs["config"]["samples"].append("s2")

    assert plan.inputs_snapshot["config"] == {"samples": ["s1"]}


def test_default_factory_returns_planner(run_service):
    from encode_pipeline.services import (
        ExecutionPlanner,
        create_default_execution_planner,
    )

    planner = create_default_execution_planner(run_service=run_service)
    assert isinstance(planner, ExecutionPlanner)


def _make_execution_plan(run_service):
    inputs = WorkflowInputs(config={}, samples=None, options={})
    record = run_service.create_run("stub", inputs)
    planner = ExecutionPlanner(run_service=run_service)
    return planner.plan_run(record.run_id).value


def test_plan_workspace_returns_pending_plan_with_workspace_plan(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(registry=run_service._registry)
    base_dir = tmp_path.resolve()

    result = workspace_planner.plan_workspace(input_plan, base_dir=base_dir)

    assert result.is_success is True
    plan = result.value
    assert plan is not None
    assert plan is not input_plan
    assert plan.status is PlanStatus.PENDING
    assert plan.workspace_plan is not None
    assert plan.workspace_plan.directories == ("logs", "results")
    assert plan.workspace_plan.files == ()
    assert plan.command_spec is None
    assert plan.can_execute is False
    assert plan.run_id == input_plan.run_id
    assert plan.workflow_id == input_plan.workflow_id
    assert plan.dag_preview == input_plan.dag_preview
    assert plan.inputs_snapshot == input_plan.inputs_snapshot


def test_workspace_planner_preserves_one_adapter_deprecation_warning(tmp_path):
    adapter = EncodeStyleWorkflowAdapter()
    row = {
        "sample": "S1",
        "fastq_1": str((tmp_path / "S1.R1.fastq.gz").resolve()),
        "fastq_2": str((tmp_path / "S1.R2.fastq.gz").resolve()),
        "layout": "PE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
        "genome": "hs",
        "bowtie2_index": str((tmp_path / "indices/hs").resolve()),
    }
    plan = ExecutionPlan(
        plan_id="plan-semantic-warning",
        run_id="run-semantic-warning",
        workflow_id=adapter.metadata.workflow_id,
        status=PlanStatus.UNSUPPORTED,
        inputs_snapshot=WorkflowInputs(
            config={
                "replicate_analysis": {"enabled": False},
                "stage4b": False,
                "chipseq_idr": {"enabled": False},
                "stage5": "FALSE",
            },
            samples=[row],
        ).to_dict(),
    )

    result = WorkspacePlanner(WorkflowRegistry([adapter])).plan_workspace(
        plan,
        base_dir=(tmp_path / "workspace").resolve(),
    )

    assert result.is_success
    assert [issue.code for issue in result.issues] == [
        "ENCODE_CONFIG_LEGACY_ALIAS_DEPRECATED",
        "ENCODE_WORKSPACE_PLANNING_COMPLETE",
    ]
    assert [issue.code for issue in result.value.issues] == [
        "ENCODE_CONFIG_LEGACY_ALIAS_DEPRECATED",
        "ENCODE_WORKSPACE_PLANNING_COMPLETE",
    ]


def test_plan_workspace_does_not_create_directories_or_files(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(registry=run_service._registry)
    base_dir = tmp_path.resolve()

    workspace_planner.plan_workspace(input_plan, base_dir=base_dir)

    assert not (base_dir / "logs").exists()
    assert not (base_dir / "results").exists()


def test_plan_workspace_does_not_mutate_input_plan(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(registry=run_service._registry)
    base_dir = tmp_path.resolve()

    workspace_planner.plan_workspace(input_plan, base_dir=base_dir)

    assert input_plan.status is PlanStatus.UNSUPPORTED
    assert input_plan.workspace_plan is None
    assert input_plan.command_spec is None


def test_plan_workspace_rejects_relative_base_dir(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(registry=run_service._registry)

    result = workspace_planner.plan_workspace(
        input_plan, base_dir=Path("relative/path")
    )

    assert result.is_failure is True
    assert result.value is None
    assert len(result.issues) == 1
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_BASE_DIR_RELATIVE"
    assert issue.path == "base_dir"


def test_plan_workspace_does_not_mutate_run_service(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    before_record = run_service.get_run(input_plan.run_id)
    before_events = run_service.list_events(input_plan.run_id)
    before_logs = run_service.list_logs(input_plan.run_id, "stdout")
    before_artifacts = run_service.list_artifacts(input_plan.run_id)

    workspace_planner = WorkspacePlanner(registry=run_service._registry)
    workspace_planner.plan_workspace(input_plan, base_dir=tmp_path.resolve())

    after_record = run_service.get_run(input_plan.run_id)
    assert after_record.status == before_record.status
    assert after_record.updated_at == before_record.updated_at
    assert run_service.list_events(input_plan.run_id) == before_events
    assert run_service.list_logs(input_plan.run_id, "stdout") == before_logs
    assert run_service.list_artifacts(input_plan.run_id) == before_artifacts


def test_plan_workspace_propagates_path_policy_failure(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(registry=run_service._registry)

    # Use a base_dir that is absolute but whose value we will monkeypatch the
    # planner to produce an invalid path. Since the default layout is fixed and
    # safe, the simplest failure path is a non-absolute base_dir.
    result = workspace_planner.plan_workspace(input_plan, base_dir=Path("not-absolute"))

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_BASE_DIR_RELATIVE"


def test_default_factory_returns_workspace_planner():
    from encode_pipeline.services import (
        WorkspacePlanner,
        create_default_workspace_planner,
    )

    planner = create_default_workspace_planner()
    assert isinstance(planner, WorkspacePlanner)


def test_workspace_planner_reconstructs_inputs_from_snapshot(run_service, tmp_path):
    from encode_pipeline.services.planning import WorkspacePlanner

    inputs = WorkflowInputs(
        config={"samples": "/abs/samples.tsv"},
        samples="/abs/samples.tsv",
        options={"strict_inputs": True},
    )
    record = run_service.create_run("stub", inputs)
    planner = ExecutionPlanner(run_service=run_service)
    plan = planner.plan_run(record.run_id).value

    workspace_planner = WorkspacePlanner(registry=run_service._registry)
    result = workspace_planner._reconstruct_inputs(plan.inputs_snapshot)

    assert result.is_success is True
    assert result.value.config == {"samples": "/abs/samples.tsv"}
    assert result.value.samples == "/abs/samples.tsv"
    assert result.value.options == {"strict_inputs": True}


def test_workspace_planner_fails_on_malformed_inputs_snapshot_missing_config(
    run_service, tmp_path
):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(registry=run_service._registry)

    result = workspace_planner.plan_workspace(
        ExecutionPlan(
            plan_id=input_plan.plan_id,
            run_id=input_plan.run_id,
            workflow_id=input_plan.workflow_id,
            status=input_plan.status,
            inputs_snapshot={"samples": None, "options": {}},
            created_at=input_plan.created_at,
        ),
        base_dir=tmp_path.resolve(),
    )

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_PLAN_MALFORMED_INPUTS"
    assert issue.path == "inputs_snapshot.config"


def test_workspace_planner_fails_on_malformed_inputs_snapshot_invalid_options_type(
    run_service, tmp_path
):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(registry=run_service._registry)

    result = workspace_planner.plan_workspace(
        ExecutionPlan(
            plan_id=input_plan.plan_id,
            run_id=input_plan.run_id,
            workflow_id=input_plan.workflow_id,
            status=input_plan.status,
            inputs_snapshot={"config": {}, "samples": None, "options": "not-a-mapping"},
            created_at=input_plan.created_at,
        ),
        base_dir=tmp_path.resolve(),
    )

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_PLAN_MALFORMED_INPUTS"
    assert issue.path == "inputs_snapshot.options"


def test_workspace_planner_fails_on_malformed_inputs_snapshot_invalid_samples_type(
    run_service, tmp_path
):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(registry=run_service._registry)

    result = workspace_planner.plan_workspace(
        ExecutionPlan(
            plan_id=input_plan.plan_id,
            run_id=input_plan.run_id,
            workflow_id=input_plan.workflow_id,
            status=input_plan.status,
            inputs_snapshot={"config": {}, "samples": 123, "options": {}},
            created_at=input_plan.created_at,
        ),
        base_dir=tmp_path.resolve(),
    )

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_PLAN_MALFORMED_INPUTS"
    assert issue.path == "inputs_snapshot.samples"


def test_workspace_planner_delegates_to_adapter_and_preserves_info_issue(
    run_service, tmp_path
):
    class _CustomAdapter:
        metadata = WorkflowMetadata(
            workflow_id="stub",
            name="Stub",
            version="0.0.1",
        )
        capabilities = WorkflowCapabilities(supports=("workspace_plan",))

        def schema(self):
            return WorkflowSchema()

        def validate(self, inputs):
            from encode_pipeline.platform.results import Result

            return Result.success(None)

        def preview_dag(self, inputs):
            from encode_pipeline.platform.results import Result

            return Result.success(DagPreview())

        def plan_workspace(self, inputs, workspace):
            from encode_pipeline.platform.results import Result, Issue

            return Result.success(
                WorkspacePlan(
                    directories=("logs",), files=(("config/config.yaml", b"x"),)
                ),
                issues=[
                    Issue(
                        code="ADAPTER_WORKSPACE_PLANNING_COMPLETE",
                        message="Adapter planned workspace.",
                        severity="info",
                        source="adapter",
                    )
                ],
            )

        def build_command(self, plan):
            from encode_pipeline.platform.results import Result

            return Result.failure([])

        def extract_artifacts(self, inputs, workspace):
            return Result.success(())

    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(
        registry=WorkflowRegistry(adapters=[_CustomAdapter()])
    )
    result = workspace_planner.plan_workspace(input_plan, base_dir=tmp_path.resolve())

    assert result.is_success is True
    updated = result.value
    assert updated.status is PlanStatus.PENDING
    assert updated.workspace_plan is not None
    assert updated.workspace_plan.directories == ("logs",)
    assert updated.workspace_plan.files == (("config/config.yaml", b"x"),)
    assert any(
        issue.code == "ADAPTER_WORKSPACE_PLANNING_COMPLETE" for issue in updated.issues
    )


def test_workspace_planner_rejects_undeclared_capability_before_adapter_call(
    run_service, tmp_path
):
    class _NoWorkspaceCapabilityAdapter:
        metadata = WorkflowMetadata(
            workflow_id="stub",
            name="Stub",
            version="0.0.1",
        )
        capabilities = WorkflowCapabilities(supports=("validation",))

        def schema(self):
            return WorkflowSchema()

        def validate(self, inputs):
            return Result.success(None)

        def preview_dag(self, inputs):
            return Result.success(DagPreview())

        def plan_workspace(self, inputs, workspace):
            raise AssertionError("undeclared workspace method was called")

        def build_command(self, plan):
            return Result.failure(
                [
                    Issue(
                        code="UNSUPPORTED",
                        message="Unsupported.",
                        source="adapter",
                    )
                ]
            )

        def extract_artifacts(self, inputs, workspace):
            return Result.failure(
                [
                    Issue(
                        code="UNSUPPORTED",
                        message="Unsupported.",
                        source="adapter",
                    )
                ]
            )

    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(
        registry=WorkflowRegistry(adapters=[_NoWorkspaceCapabilityAdapter()])
    )

    result = workspace_planner.plan_workspace(input_plan, base_dir=tmp_path.resolve())

    assert result.is_failure
    assert result.issues[0].code == "WORKSPACE_PLAN_CAPABILITY_UNSUPPORTED"
    assert result.issues[0].context == {
        "workflow_id": "stub",
        "capability": "workspace_plan",
    }


def test_workspace_planner_fails_when_workflow_not_in_registry(run_service, tmp_path):
    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(registry=WorkflowRegistry(adapters=[]))

    result = workspace_planner.plan_workspace(input_plan, base_dir=tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_PLAN_WORKFLOW_NOT_FOUND"
    assert issue.path == "workflow_id"


def test_workspace_planner_rejects_adapter_returned_absolute_path(
    run_service, tmp_path
):
    class _BadAdapter:
        metadata = WorkflowMetadata(
            workflow_id="stub",
            name="Stub",
            version="0.0.1",
        )
        capabilities = WorkflowCapabilities(supports=("workspace_plan",))

        def schema(self):
            return WorkflowSchema()

        def validate(self, inputs):
            from encode_pipeline.platform.results import Result

            return Result.success(None)

        def preview_dag(self, inputs):
            from encode_pipeline.platform.results import Result

            return Result.success(DagPreview())

        def plan_workspace(self, inputs, workspace):
            from encode_pipeline.platform.results import Result

            return Result.success(WorkspacePlan(files=(("/etc/passwd", b"x"),)))

        def build_command(self, plan):
            from encode_pipeline.platform.results import Result

            return Result.failure([])

        def extract_artifacts(self, inputs, workspace):
            return Result.success(())

    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(
        registry=WorkflowRegistry(adapters=[_BadAdapter()])
    )
    result = workspace_planner.plan_workspace(input_plan, base_dir=tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_PATH_ABSOLUTE"
    assert issue.path == "workspace_plan.files[0]"
    payload = issue.to_dict()
    for value in payload.values():
        if isinstance(value, str):
            assert "/etc/passwd" not in value
    for ctx_value in payload["context"].values():
        if isinstance(ctx_value, str):
            assert "/etc/passwd" not in ctx_value


def test_workspace_planner_rejects_adapter_returned_traversal_path(
    run_service, tmp_path
):
    class _BadAdapter:
        metadata = WorkflowMetadata(
            workflow_id="stub",
            name="Stub",
            version="0.0.1",
        )
        capabilities = WorkflowCapabilities(supports=("workspace_plan",))

        def schema(self):
            return WorkflowSchema()

        def validate(self, inputs):
            from encode_pipeline.platform.results import Result

            return Result.success(None)

        def preview_dag(self, inputs):
            from encode_pipeline.platform.results import Result

            return Result.success(DagPreview())

        def plan_workspace(self, inputs, workspace):
            from encode_pipeline.platform.results import Result

            return Result.success(WorkspacePlan(files=(("../escape", b"x"),)))

        def build_command(self, plan):
            from encode_pipeline.platform.results import Result

            return Result.failure([])

        def extract_artifacts(self, inputs, workspace):
            return Result.success(())

    input_plan = _make_execution_plan(run_service)
    workspace_planner = WorkspacePlanner(
        registry=WorkflowRegistry(adapters=[_BadAdapter()])
    )
    result = workspace_planner.plan_workspace(input_plan, base_dir=tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "WORKSPACE_PATH_TRAVERSAL"
    assert issue.path == "workspace_plan.files[0]"
