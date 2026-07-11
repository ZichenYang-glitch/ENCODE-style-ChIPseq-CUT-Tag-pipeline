"""Tests for the CommandBuilder command-spec construction boundary."""

from pathlib import Path
from uuid import uuid4

import pytest


_WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


def _make_pending_plan(
    *,
    workflow_id: str = _WORKFLOW_ID,
    inputs_snapshot: dict | None = None,
    workspace_plan_files: tuple[tuple[str, bytes], ...] = (
        ("config/config.yaml", b""),
    ),
):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus

    return ExecutionPlan(
        plan_id=str(uuid4()),
        run_id="run-1",
        workflow_id=workflow_id,
        status=PlanStatus.PENDING,
        inputs_snapshot=inputs_snapshot or {},
        workspace_plan=WorkspacePlan(files=workspace_plan_files),
    )


def test_command_builder_requires_execution_plan(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command("not-a-plan", tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "COMMAND_BUILD_INVALID_PLAN"
    assert issue.path == "plan"
    assert issue.source == "command_builder"
    assert issue.severity.value == "error"


def test_command_builder_requires_absolute_base_dir():
    from encode_pipeline.services.command_builder import CommandBuilder

    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(_make_pending_plan(), Path("relative/path"))

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "COMMAND_BUILD_BASE_DIR_RELATIVE"
    assert issue.path == "base_dir"


def test_command_builder_requires_pending_status(tmp_path):
    from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
    from encode_pipeline.services.command_builder import CommandBuilder

    for status in (PlanStatus.PLANNED, PlanStatus.UNSUPPORTED):
        plan = ExecutionPlan(
            plan_id="p1",
            run_id="run-1",
            workflow_id=_WORKFLOW_ID,
            status=status,
            inputs_snapshot={},
        )
        builder = CommandBuilder(registry=_make_registry())
        result = builder.build_command(plan, tmp_path.resolve())

        assert result.is_failure is True
        issue = result.issues[0]
        assert issue.code == "COMMAND_BUILD_INVALID_PLAN_STATUS"
        assert issue.path == "plan"


def test_command_builder_requires_workspace_plan(tmp_path):
    from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
    from encode_pipeline.services.command_builder import CommandBuilder

    plan = ExecutionPlan(
        plan_id="p1",
        run_id="run-1",
        workflow_id=_WORKFLOW_ID,
        status=PlanStatus.PENDING,
        inputs_snapshot={},
        workspace_plan=None,
    )
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "COMMAND_BUILD_MISSING_WORKSPACE_PLAN"
    assert issue.path == "workspace_plan"


def test_command_builder_requires_config_file_in_workspace(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    plan = _make_pending_plan(workspace_plan_files=(("samples.tsv", b""),))
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "COMMAND_BUILD_MISSING_CONFIG"
    assert issue.path == "workspace_plan.files"


def test_command_builder_requires_registered_workflow(tmp_path):
    from encode_pipeline.platform.adapters import WorkspacePlan
    from encode_pipeline.platform.planning import ExecutionPlan
    from encode_pipeline.platform.registry import WorkflowRegistry
    from encode_pipeline.services.command_builder import CommandBuilder

    registry = WorkflowRegistry(adapters=[])
    plan = ExecutionPlan(
        plan_id="p1",
        run_id="run-1",
        workflow_id="unknown-workflow",
        status="pending",
        inputs_snapshot={},
        workspace_plan=WorkspacePlan(files=(("config/config.yaml", b""),)),
    )
    builder = CommandBuilder(registry=registry)
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "COMMAND_BUILD_UNSUPPORTED_WORKFLOW"
    assert issue.path == "workflow"


def test_command_builder_requires_supported_engine(tmp_path):
    from encode_pipeline.platform.adapters import (
        WorkspacePlan,
        WorkflowCapabilities,
        WorkflowMetadata,
    )
    from encode_pipeline.platform.planning import ExecutionPlan
    from encode_pipeline.platform.registry import WorkflowRegistry
    from encode_pipeline.services.command_builder import CommandBuilder

    class OtherAdapter:
        metadata = WorkflowMetadata(
            workflow_id="other-engine",
            name="Other",
            version="0.0.0",
            engines=("other",),
        )
        capabilities = WorkflowCapabilities(supports=())

        def schema(self): ...
        def validate(self, inputs): ...
        def preview_dag(self, inputs): ...
        def plan_workspace(self, inputs, workspace): ...
        def build_command(self, plan): ...

    registry = WorkflowRegistry(adapters=[OtherAdapter()])
    plan = ExecutionPlan(
        plan_id="p1",
        run_id="run-1",
        workflow_id="other-engine",
        status="pending",
        inputs_snapshot={},
        workspace_plan=WorkspacePlan(files=(("config/config.yaml", b""),)),
    )
    builder = CommandBuilder(registry=registry)
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "COMMAND_BUILD_UNSUPPORTED_ENGINE"
    assert issue.path == "workflow"


def test_command_builder_rejects_invalid_cores(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    for cores in ("four", 0, -1, 3.14, True, []):
        plan = _make_pending_plan(inputs_snapshot={"options": {"cores": cores}})
        builder = CommandBuilder(registry=_make_registry())
        result = builder.build_command(plan, tmp_path.resolve())

        assert result.is_failure is True, f"cores={cores!r} should fail"
        issue = result.issues[0]
        assert issue.code == "COMMAND_BUILD_INVALID_CORES"
        assert issue.path == "plan.inputs_snapshot.options.cores"


def test_command_builder_defaults_cores_to_one(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    plan = _make_pending_plan()
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_success is True
    assert result.value.command_spec.argv[-1] == "1"


def test_command_builder_defaults_cores_when_options_missing(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    plan = _make_pending_plan(inputs_snapshot={})
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_success is True
    assert result.value.command_spec.argv[-1] == "1"


def test_command_builder_defaults_cores_when_cores_missing(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    plan = _make_pending_plan(inputs_snapshot={"options": {}})
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_success is True
    assert result.value.command_spec.argv[-1] == "1"


def test_command_builder_returns_planned_execution_plan(tmp_path):
    from encode_pipeline.platform.planning import PlanStatus
    from encode_pipeline.services.command_builder import CommandBuilder

    plan = _make_pending_plan()
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_success is True
    output = result.value
    assert output.status is PlanStatus.PLANNED
    assert output.command_spec is not None
    assert output.plan_id != plan.plan_id
    assert output.run_id == plan.run_id
    assert output.workflow_id == plan.workflow_id


def test_command_builder_command_spec_has_expected_argv(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    base_dir = tmp_path.resolve()
    plan = _make_pending_plan(inputs_snapshot={"options": {"cores": 4}})
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, base_dir)

    assert result.is_success is True
    argv = result.value.command_spec.argv
    assert argv[0] == "snakemake"
    assert argv[1] == "--snakefile"
    assert argv[3] == "--directory"
    assert argv[4] == str(base_dir)
    assert argv[5] == "--configfile"
    assert argv[6] == str(base_dir / "config" / "config.yaml")
    assert argv[7] == "--cores"
    assert argv[8] == "4"

    # Snakefile path is inside the repo
    snakefile_path = Path(argv[2])
    assert snakefile_path.is_file()
    assert snakefile_path.name == "Snakefile"
    assert "workflow" in snakefile_path.parts


def test_command_builder_preserves_inputs_and_workspace(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    plan = _make_pending_plan(inputs_snapshot={"options": {"cores": 2}})
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_success is True
    assert result.value.inputs_snapshot == plan.inputs_snapshot
    assert result.value.workspace_plan == plan.workspace_plan


def test_command_builder_command_spec_has_empty_env_and_no_cwd(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    plan = _make_pending_plan()
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_success is True
    assert result.value.command_spec.cwd is None
    assert result.value.command_spec.env == {}


def test_command_builder_refuses_missing_snakefile(tmp_path, monkeypatch):
    from encode_pipeline.services import command_builder
    from encode_pipeline.services.command_builder import CommandBuilder

    builder = CommandBuilder(registry=_make_registry())
    fake_repo = tmp_path / "repo"
    fake_repo.mkdir()
    monkeypatch.setattr(
        command_builder,
        "_bundled_snakefile_path",
        lambda: fake_repo / "workflow" / "Snakefile",
    )

    result = builder.build_command(_make_pending_plan(), tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "COMMAND_BUILD_SNAKEFILE_NOT_FOUND"
    assert issue.path == "workflow"


def test_command_builder_does_not_swallow_worker_hard_timeout(
    tmp_path,
    monkeypatch,
):
    from encode_pipeline.platform.planning import WorkspacePathPolicy
    from encode_pipeline.services.command_builder import CommandBuilder
    from encode_pipeline.workers.timeouts import WorkerHardTimeout

    def timeout(_self, _path):
        raise WorkerHardTimeout("RQ deadline reached")

    monkeypatch.setattr(WorkspacePathPolicy, "resolve", timeout)

    with pytest.raises(WorkerHardTimeout, match="RQ deadline reached"):
        CommandBuilder(registry=_make_registry()).build_command(
            _make_pending_plan(),
            tmp_path.resolve(),
        )


def test_command_builder_issues_use_logical_paths(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    plan = _make_pending_plan(inputs_snapshot={"options": {"cores": -1}})
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert str(tmp_path) not in issue.message
    assert str(tmp_path) not in (issue.technical_message or "")
    assert str(tmp_path) not in (issue.path or "")


def test_command_builder_info_issue_does_not_leak_command_or_paths(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    plan = _make_pending_plan()
    builder = CommandBuilder(registry=_make_registry())
    result = builder.build_command(plan, tmp_path.resolve())

    assert result.is_success is True
    info_issue = result.value.issues[-1]
    assert info_issue.code == "COMMAND_BUILDING_COMPLETE"
    assert info_issue.severity.value == "info"
    assert info_issue.source == "command_builder"
    assert info_issue.path == "command_spec"
    assert "snakemake" not in info_issue.message
    assert "config" not in info_issue.message
    assert str(tmp_path) not in info_issue.message


def test_command_builder_import_boundary():
    import ast
    from pathlib import Path

    source_path = (
        Path(__file__).resolve().parents[2]
        / "src/encode_pipeline/services/command_builder.py"
    )
    source = source_path.read_text(encoding="utf-8")
    tree = ast.parse(source)

    forbidden_modules = {
        "encode_pipeline.api",
        "encode_pipeline.frontend",
        "fastapi",
        "pydantic",
        "snakemake",
        "subprocess",
    }
    imported_modules: set[str] = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported_modules.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module is not None:
            imported_modules.add(node.module)

    assert not any(
        module == forbidden or module.startswith(f"{forbidden}.")
        for module in imported_modules
        for forbidden in forbidden_modules
    )


def _make_registry():
    from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
    from encode_pipeline.platform.registry import WorkflowRegistry

    return WorkflowRegistry(adapters=[EncodeStyleWorkflowAdapter()])
