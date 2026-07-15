"""Tests for the CommandBuilder command-spec construction boundary."""

from pathlib import Path
from uuid import uuid4

import pytest


_WORKFLOW_ID = "encode-style-chipseq-cuttag-atac-mnase"


class _CommandAdapter:
    def __init__(self, callback):
        from encode_pipeline.platform.adapters import (
            WorkflowCapabilities,
            WorkflowMetadata,
        )

        self.metadata = WorkflowMetadata(
            workflow_id="delegated-command",
            name="Delegated command",
            version="1.0.0",
            engines=("opaque-engine",),
        )
        self.capabilities = WorkflowCapabilities(supports=("command",))
        self._callback = callback

    def schema(self): ...
    def validate(self, inputs): ...
    def preview_dag(self, inputs): ...
    def plan_workspace(self, inputs, workspace): ...

    def build_command(self, plan, workspace):
        return self._callback(plan, workspace)

    def extract_artifacts(self, inputs, workspace): ...


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
        def extract_artifacts(self, inputs, workspace): ...

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


def test_command_builder_legacy_command_declares_exact_dry_run_preflight(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    result = CommandBuilder(registry=_make_registry()).build_command(
        _make_pending_plan(),
        tmp_path.resolve(),
    )

    assert result.is_success
    command = result.value.command_spec
    assert command is not None
    assert command.preflight_argv == command.argv + ("-n",)


def test_command_builder_delegates_command_capability_with_absolute_workspace(
    tmp_path,
):
    from encode_pipeline.platform.adapters import CommandSpec
    from encode_pipeline.platform.registry import WorkflowRegistry
    from encode_pipeline.platform.results import Result
    from encode_pipeline.services.command_builder import CommandBuilder

    observed = {}

    def build(workspace_plan, workspace):
        observed["plan"] = workspace_plan
        observed["workspace"] = workspace
        launch_dir = workspace / "engine" / "launch"
        return Result.success(
            CommandSpec(
                argv=("pinned-engine", "run"),
                cwd=str(launch_dir),
                preflight_argv=("pinned-engine", "config"),
            )
        )

    adapter = _CommandAdapter(build)
    plan = _make_pending_plan(
        workflow_id=adapter.metadata.workflow_id,
        workspace_plan_files=(),
    )
    workspace = tmp_path.resolve()

    result = CommandBuilder(registry=WorkflowRegistry([adapter])).build_command(
        plan,
        workspace,
    )

    assert result.is_success
    assert observed == {"plan": plan.workspace_plan, "workspace": workspace}
    assert result.value.command_spec == CommandSpec(
        argv=("pinned-engine", "run"),
        cwd=str(workspace / "engine" / "launch"),
        preflight_argv=("pinned-engine", "config"),
    )


def test_command_builder_accepts_managed_logs_inside_workspace(tmp_path):
    from encode_pipeline.platform.adapters import CommandSpec
    from encode_pipeline.platform.registry import WorkflowRegistry
    from encode_pipeline.platform.results import Result
    from encode_pipeline.services.command_builder import CommandBuilder

    def build(_plan, workspace):
        return Result.success(
            CommandSpec(
                argv=("pinned-engine", "run"),
                cwd=str(workspace / "engine" / "launch"),
                preflight_argv=("pinned-engine", "config"),
                preflight_kind="configuration",
                preflight_managed_logs=(
                    ("engine_preflight", str(workspace / "logs/preflight.log")),
                ),
                execution_managed_logs=(
                    ("engine", str(workspace / "logs/engine.log")),
                ),
            )
        )

    adapter = _CommandAdapter(build)
    result = CommandBuilder(registry=WorkflowRegistry([adapter])).build_command(
        _make_pending_plan(
            workflow_id=adapter.metadata.workflow_id,
            workspace_plan_files=(),
        ),
        tmp_path.resolve(),
    )

    assert result.is_success
    assert result.value.command_spec.preflight_kind == "configuration"


@pytest.mark.parametrize("case", ("outside", "planned-file", "workspace-root"))
def test_command_builder_rejects_unsafe_managed_log_path(tmp_path, case):
    from encode_pipeline.platform.adapters import CommandSpec
    from encode_pipeline.platform.registry import WorkflowRegistry
    from encode_pipeline.platform.results import Result
    from encode_pipeline.services.command_builder import CommandBuilder

    workspace = tmp_path.resolve()
    paths = {
        "outside": str(workspace.parent / "outside.log"),
        "planned-file": str(workspace / "config/planned.json"),
        "workspace-root": str(workspace),
    }
    adapter = _CommandAdapter(
        lambda _plan, _workspace: Result.success(
            CommandSpec(
                argv=("pinned-engine", "run"),
                cwd=str(workspace / "engine/launch"),
                execution_managed_logs=(("engine", paths[case]),),
            )
        )
    )
    result = CommandBuilder(registry=WorkflowRegistry([adapter])).build_command(
        _make_pending_plan(
            workflow_id=adapter.metadata.workflow_id,
            workspace_plan_files=(("config/planned.json", b"{}"),),
        ),
        workspace,
    )

    assert result.is_failure
    assert result.issues[0].code == "COMMAND_BUILD_ADAPTER_FAILED"
    assert paths[case] not in str(result.issues[0].to_dict())


def test_command_builder_sanitizes_adapter_failure(tmp_path):
    from encode_pipeline.platform.registry import WorkflowRegistry
    from encode_pipeline.platform.results import Issue, Result
    from encode_pipeline.services.command_builder import CommandBuilder

    secret = str(tmp_path / "private" / "source")
    adapter = _CommandAdapter(
        lambda _plan, _workspace: Result.failure(
            [
                Issue(
                    code="PRIVATE_FAILURE",
                    message=secret,
                    technical_message=secret,
                )
            ]
        )
    )
    result = CommandBuilder(registry=WorkflowRegistry([adapter])).build_command(
        _make_pending_plan(
            workflow_id=adapter.metadata.workflow_id,
            workspace_plan_files=(),
        ),
        tmp_path.resolve(),
    )

    assert result.is_failure
    assert [issue.code for issue in result.issues] == ["COMMAND_BUILD_ADAPTER_FAILED"]
    assert secret not in str(result.issues[0].to_dict())


def test_command_builder_sanitizes_adapter_exception(tmp_path):
    from encode_pipeline.platform.registry import WorkflowRegistry
    from encode_pipeline.services.command_builder import CommandBuilder

    secret = str(tmp_path / "private" / "source")

    def raise_private(_plan, _workspace):
        raise RuntimeError(secret)

    adapter = _CommandAdapter(raise_private)
    result = CommandBuilder(registry=WorkflowRegistry([adapter])).build_command(
        _make_pending_plan(
            workflow_id=adapter.metadata.workflow_id,
            workspace_plan_files=(),
        ),
        tmp_path.resolve(),
    )

    assert result.is_failure
    assert [issue.code for issue in result.issues] == ["COMMAND_BUILD_ADAPTER_FAILED"]
    assert secret not in str(result.issues[0].to_dict())


def test_command_builder_rejects_adapter_cwd_outside_workspace(tmp_path):
    from encode_pipeline.platform.adapters import CommandSpec
    from encode_pipeline.platform.registry import WorkflowRegistry
    from encode_pipeline.platform.results import Result
    from encode_pipeline.services.command_builder import CommandBuilder

    adapter = _CommandAdapter(
        lambda _plan, _workspace: Result.success(
            CommandSpec(argv=("pinned-engine", "run"), cwd="/outside")
        )
    )
    result = CommandBuilder(registry=WorkflowRegistry([adapter])).build_command(
        _make_pending_plan(
            workflow_id=adapter.metadata.workflow_id,
            workspace_plan_files=(),
        ),
        tmp_path.resolve(),
    )

    assert result.is_failure
    assert result.issues[0].code == "COMMAND_BUILD_ADAPTER_FAILED"


@pytest.mark.parametrize("matches", (True, False))
def test_command_builder_requires_platform_derived_container_scope(tmp_path, matches):
    from encode_pipeline.platform.adapters import CommandSpec
    from encode_pipeline.platform.managed_containers import managed_container_scope
    from encode_pipeline.platform.registry import WorkflowRegistry
    from encode_pipeline.platform.results import Result
    from encode_pipeline.services.command_builder import CommandBuilder

    workspace = tmp_path.resolve()
    scope = managed_container_scope(workspace) if matches else "f" * 64
    adapter = _CommandAdapter(
        lambda _plan, selected_workspace: Result.success(
            CommandSpec(
                argv=("pinned-engine", "run"),
                cwd=str(selected_workspace),
                managed_container_scope=scope,
                managed_container_endpoint_identity="0" * 64,
            )
        )
    )

    result = CommandBuilder(registry=WorkflowRegistry([adapter])).build_command(
        _make_pending_plan(
            workflow_id=adapter.metadata.workflow_id,
            workspace_plan_files=(),
        ),
        workspace,
    )

    assert result.is_success is matches
    if matches:
        assert result.value.command_spec.managed_container_scope == scope
    else:
        assert result.issues[0].code == "COMMAND_BUILD_ADAPTER_FAILED"
        assert scope not in str(result.issues[0].to_dict())


@pytest.mark.parametrize(
    "case",
    ("missing-cwd", "different-preflight-executable"),
)
def test_command_builder_rejects_incomplete_adapter_launch_boundary(
    tmp_path,
    case,
):
    from encode_pipeline.platform.adapters import CommandSpec
    from encode_pipeline.platform.registry import WorkflowRegistry
    from encode_pipeline.platform.results import Result
    from encode_pipeline.services.command_builder import CommandBuilder

    def build(_plan, workspace):
        if case == "missing-cwd":
            command_spec = CommandSpec(argv=("pinned-engine", "run"))
        else:
            command_spec = CommandSpec(
                argv=("pinned-engine", "run"),
                cwd=str(workspace),
                preflight_argv=("other-engine", "validate"),
            )
        return Result.success(command_spec)

    adapter = _CommandAdapter(build)
    result = CommandBuilder(registry=WorkflowRegistry([adapter])).build_command(
        _make_pending_plan(
            workflow_id=adapter.metadata.workflow_id,
            workspace_plan_files=(),
        ),
        tmp_path.resolve(),
    )

    assert result.is_failure
    assert result.issues[0].code == "COMMAND_BUILD_ADAPTER_FAILED"


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


def test_command_builder_refuses_missing_snakefile(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    fake_repo = tmp_path / "repo"
    fake_repo.mkdir()
    builder = CommandBuilder(
        registry=_make_registry(),
        project_root=fake_repo.resolve(),
    )

    result = builder.build_command(_make_pending_plan(), tmp_path.resolve())

    assert result.is_failure is True
    issue = result.issues[0]
    assert issue.code == "COMMAND_BUILD_SNAKEFILE_NOT_FOUND"
    assert issue.path == "workflow"


def test_command_builder_uses_snakefile_from_explicit_project_root(tmp_path):
    from encode_pipeline.services.command_builder import CommandBuilder

    project_root = tmp_path / "source"
    workflow_dir = project_root / "workflow"
    workflow_dir.mkdir(parents=True)
    snakefile = workflow_dir / "Snakefile"
    snakefile.write_text("rule all:\n    input: []\n", encoding="utf-8")
    builder = CommandBuilder(
        registry=_make_registry(),
        project_root=project_root.resolve(),
    )

    result = builder.build_command(_make_pending_plan(), tmp_path.resolve())

    assert result.is_success is True
    assert result.value.command_spec.argv[2] == str(snakefile.resolve())


def test_command_builder_requires_absolute_project_root():
    from encode_pipeline.services.command_builder import CommandBuilder

    with pytest.raises(
        ValueError,
        match="project_root must be an absolute pathlib.Path",
    ):
        CommandBuilder(
            registry=_make_registry(),
            project_root=Path("relative-source"),
        )


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
