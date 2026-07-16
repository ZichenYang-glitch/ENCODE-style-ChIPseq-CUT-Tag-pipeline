"""Reusable adapter-conformance coverage for bundled and test-only adapters."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone
from decimal import Decimal
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

import pytest
import yaml

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagNode,
    DagPreview,
    ExtractedArtifactCandidate,
    ExtractedQcMetricCandidate,
    QcSourceArtifact,
    QcSourceDocument,
    WorkflowAuthoringModes,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowInputModes,
    WorkflowMetadata,
    WorkflowSchema,
    WorkflowSchemaCoverage,
    WorkspacePlan,
)
from encode_pipeline.platform.planning import ExecutionPlan, PlanStatus
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.runs import RunStatus
from encode_pipeline.services.artifact_extraction import ArtifactExtractionService
from encode_pipeline.services.planning import WorkspacePlanner
from encode_pipeline.services.qc_summary_indexing import QcSummaryIndexingService
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider
from encode_pipeline.services.workflow_info import WorkflowInfoService
from encode_pipeline.testing.adapter_conformance import (
    AdapterConformanceCase,
    AdapterConformanceError,
    verify_adapter_conformance,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]


def _unsupported(method: str) -> Result:
    return Result.failure(
        [
            Issue(
                code="MINIMAL_CAPABILITY_UNSUPPORTED",
                message="The capability is not supported.",
                path=method,
                source="adapter",
                context={"method": method},
            )
        ]
    )


class MinimalConformantAdapter:
    """A deterministic workflow-neutral adapter that exists only in tests."""

    metadata = WorkflowMetadata(
        workflow_id="minimal-test-workflow",
        name="Minimal test workflow",
        version="1.0.0",
        engines=("minimal",),
        tags=("test",),
    )
    capabilities = WorkflowCapabilities(
        supports=(
            "validation",
            "dag_preview",
            "workspace_plan",
            "command",
            "input_authoring",
            "artifact_extract",
            "qc_summary_extract",
        )
    )

    def __init__(self, *, command_fails: bool = False) -> None:
        self._command_fails = command_fails

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema(
            schema_version="1.0.0",
            coverage=WorkflowSchemaCoverage(
                config="complete", samples="complete", options="complete"
            ),
            authoring_modes=WorkflowAuthoringModes(
                config=("schema_form",),
                samples=("inline_table",),
                options=("schema_form",),
            ),
            input_modes=WorkflowInputModes(
                config=("object",), samples=("inline_rows",), options=("object",)
            ),
            config_schema={
                "type": "object",
                "properties": {"label": {"type": "string"}},
                "required": ["label"],
                "additionalProperties": False,
            },
            sample_schema={
                "type": "array",
                "items": {
                    "type": "object",
                    "properties": {"sample": {"type": "string"}},
                    "required": ["sample"],
                    "additionalProperties": False,
                },
                "minItems": 1,
            },
            option_schema={
                "type": "object",
                "properties": {},
                "additionalProperties": False,
            },
        )

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        if (
            not isinstance(inputs.config.get("label"), str)
            or not isinstance(inputs.samples, list)
            or not inputs.samples
        ):
            return Result.failure(
                [
                    Issue(
                        code="MINIMAL_INPUTS_INVALID",
                        message="Workflow inputs are invalid.",
                        path="inputs",
                        source="adapter",
                    )
                ]
            )
        return Result.success(
            {
                "label": inputs.config["label"],
                "sample_count": len(inputs.samples),
            }
        )

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        return Result.success(
            DagPreview(nodes=(DagNode(id="run", label="Run workflow"),))
        )

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[WorkspacePlan]:
        validated = self.validate(inputs)
        if validated.is_failure:
            return Result.failure(validated.issues)
        payload = json.dumps(
            inputs.to_dict(), sort_keys=True, separators=(",", ":")
        ).encode("utf-8")
        return Result.success(
            WorkspacePlan(
                directories=("config", "results"),
                files=(("config/inputs.json", payload),),
            )
        )

    def build_command(
        self,
        plan: WorkspacePlan,
        workspace: str | Path,
    ) -> Result[CommandSpec]:
        if self._command_fails:
            return _unsupported("build_command")
        return Result.success(
            CommandSpec(
                argv=("minimal-workflow", "run", "config/inputs.json"),
                cwd=str(workspace),
            )
        )

    def extract_artifacts(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[tuple[ExtractedArtifactCandidate, ...]]:
        return Result.success(
            (
                ExtractedArtifactCandidate(
                    output_type="qc_summary",
                    relative_path="results/qc.tsv",
                    mime_type="text/tab-separated-values",
                    metadata={"scope": "run"},
                ),
            )
        )

    def qc_source_output_types(self) -> tuple[str, ...]:
        return ("qc_summary",)

    def extract_qc_metrics(
        self,
        inputs: WorkflowInputs,
        sources: tuple[QcSourceDocument, ...],
    ) -> Result[tuple[ExtractedQcMetricCandidate, ...]]:
        return Result.success(
            (
                ExtractedQcMetricCandidate(
                    metric_key="reads.total",
                    display_name="Total reads",
                    value=Decimal("10"),
                    unit="count",
                    scope="run",
                    source_artifact_id=sources[0].source.artifact_id,
                ),
            )
        )


def _minimal_case(tmp_path: Path, *, command_fails: bool = False):
    planning_workspace = tmp_path / "planned-workspace"
    artifact_workspace = tmp_path / "artifact-workspace"
    (artifact_workspace / "results").mkdir(parents=True)
    (artifact_workspace / "results/qc.tsv").write_text(
        "metric\tvalue\nreads.total\t10\n", encoding="utf-8"
    )
    valid_inputs = WorkflowInputs(
        config={"label": "demo"},
        samples=[{"sample": "S1"}],
        options={},
    )
    invalid_inputs = WorkflowInputs(config={}, samples=[{"sample": "S1"}])
    qc_content = b"metric\tvalue\nreads.total\t10\n"
    qc_sources = (
        QcSourceDocument(
            source=QcSourceArtifact(
                artifact_id="artifact-qc",
                output_type="qc_summary",
                relative_path="results/qc.tsv",
                metadata={"scope": "run"},
            ),
            content=qc_content,
        ),
    )
    return AdapterConformanceCase(
        adapter=MinimalConformantAdapter(command_fails=command_fails),
        valid_inputs=valid_inputs,
        invalid_inputs=invalid_inputs,
        planning_workspace=planning_workspace,
        artifact_workspace=artifact_workspace,
        qc_sources=qc_sources,
    )


def _encode_case(tmp_path: Path, monkeypatch) -> AdapterConformanceCase:
    config = yaml.safe_load(
        (PROJECT_ROOT / "test/profiles/platform_worker_tiny/config.yaml").read_text(
            encoding="utf-8"
        )
    )
    config.pop("samples")
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
    valid_inputs = WorkflowInputs(
        config=config,
        samples=[row],
        options={"strict_inputs": False},
    )
    invalid_inputs = WorkflowInputs(
        config={**config, "threads": 0},
        samples=[row],
        options={"strict_inputs": False},
    )
    adapter = EncodeStyleWorkflowAdapter()
    artifact_workspace = tmp_path / "encode-artifact-workspace"
    planned = adapter.plan_workspace(valid_inputs, artifact_workspace)
    assert planned.is_success
    for directory in planned.value.directories:
        (artifact_workspace / directory).mkdir(parents=True, exist_ok=True)
    for relative_path, content in planned.value.files:
        destination = artifact_workspace / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_bytes(content)
    monkeypatch.setattr(
        "encode_pipeline.manifest.make.build_manifest_rows",
        lambda *_args, **_kwargs: ([], 0, 0),
    )
    return AdapterConformanceCase(
        adapter=adapter,
        valid_inputs=valid_inputs,
        invalid_inputs=invalid_inputs,
        planning_workspace=tmp_path / "encode-planned-workspace",
        artifact_workspace=artifact_workspace,
        qc_sources=(),
    )


def test_minimal_adapter_passes_reusable_conformance_suite(tmp_path):
    verify_adapter_conformance(_minimal_case(tmp_path))


def test_encode_adapter_passes_same_reusable_conformance_suite(tmp_path, monkeypatch):
    verify_adapter_conformance(_encode_case(tmp_path, monkeypatch))


def test_conformance_rejects_declared_command_that_returns_failure(tmp_path):
    with pytest.raises(AdapterConformanceError, match="capability.command"):
        verify_adapter_conformance(_minimal_case(tmp_path, command_fails=True))


def test_conformance_rejects_invalid_json_schema_without_leaking_values(tmp_path):
    class InvalidSchemaAdapter(MinimalConformantAdapter):
        def schema(self) -> WorkflowSchema:
            return WorkflowSchema(
                config_schema={
                    "type": "object",
                    "properties": {"secret": {"type": "not-a-json-schema-type"}},
                }
            )

    case = replace(_minimal_case(tmp_path), adapter=InvalidSchemaAdapter())

    with pytest.raises(AdapterConformanceError) as raised:
        verify_adapter_conformance(case)

    assert str(raised.value) == ("adapter.schema: must satisfy JSON Schema 2020-12")
    assert "secret" not in str(raised.value)
    assert raised.value.__context__ is None
    assert "secret" not in "".join(traceback.format_exception(raised.value))


def test_conformance_suppresses_adapter_exception_context(tmp_path):
    class RaisingSchemaAdapter(MinimalConformantAdapter):
        def schema(self) -> WorkflowSchema:
            raise RuntimeError("/private/secret TOKEN=abc")

    case = replace(_minimal_case(tmp_path), adapter=RaisingSchemaAdapter())

    with pytest.raises(AdapterConformanceError) as raised:
        verify_adapter_conformance(case)

    assert str(raised.value) == "adapter.schema: raised an exception"
    assert raised.value.__context__ is None
    assert raised.value.__suppress_context__ is True
    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/secret" not in formatted
    assert "TOKEN=abc" not in formatted


def test_conformance_suppresses_schema_serialization_exception_context(tmp_path):
    class RaisingWorkflowSchema(WorkflowSchema):
        def to_dict(self) -> dict[str, object]:
            raise RuntimeError("/private/schema TOKEN=def")

    class RaisingSchemaDocumentAdapter(MinimalConformantAdapter):
        def schema(self) -> WorkflowSchema:
            return RaisingWorkflowSchema()

    case = replace(_minimal_case(tmp_path), adapter=RaisingSchemaDocumentAdapter())

    with pytest.raises(AdapterConformanceError) as raised:
        verify_adapter_conformance(case)

    assert str(raised.value) == "adapter.schema: raised an exception"
    assert raised.value.__context__ is None
    assert raised.value.__suppress_context__ is True
    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/schema" not in formatted
    assert "TOKEN=def" not in formatted


@pytest.mark.parametrize(
    ("attribute", "coordinate"),
    (("metadata", "adapter.metadata"), ("capabilities", "adapter.capabilities")),
)
def test_conformance_suppresses_declaration_descriptor_context(
    tmp_path,
    attribute,
    coordinate,
):
    class RaisingDeclarationAdapter(MinimalConformantAdapter):
        def __getattribute__(self, name: str):
            if name == attribute:
                raise RuntimeError("/private/declaration TOKEN=ghi")
            return super().__getattribute__(name)

    case = replace(_minimal_case(tmp_path), adapter=RaisingDeclarationAdapter())

    with pytest.raises(AdapterConformanceError) as raised:
        verify_adapter_conformance(case)

    assert str(raised.value) == f"{coordinate}: raised an exception"
    assert raised.value.__context__ is None
    assert raised.value.__suppress_context__ is True
    formatted = "".join(traceback.format_exception(raised.value))
    assert "/private/declaration" not in formatted
    assert "TOKEN=ghi" not in formatted


def test_conformance_support_import_does_not_import_pytest():
    code = """
import sys
import encode_pipeline.testing.adapter_conformance
print("pytest" in sys.modules)
"""
    environment = dict(os.environ)
    environment["PYTHONPATH"] = str(PROJECT_ROOT / "src")
    environment["PYTHONDONTWRITEBYTECODE"] = "1"

    completed = subprocess.run(
        [sys.executable, "-c", code],
        check=False,
        capture_output=True,
        text=True,
        env=environment,
    )

    assert completed.returncode == 0, completed.stderr
    assert completed.stdout.strip() == "False"


def test_minimal_adapter_uses_existing_platform_services_without_id_branches(tmp_path):
    case = _minimal_case(tmp_path)
    registry = WorkflowRegistry(adapters=[case.adapter])
    workflow_info = WorkflowInfoService(registry)
    validation = ValidationService(registry)
    workspace_planner = WorkspacePlanner(registry)
    execution_plan = ExecutionPlan(
        plan_id="plan-minimal",
        run_id="run-minimal",
        workflow_id=case.adapter.metadata.workflow_id,
        status=PlanStatus.PENDING,
        inputs_snapshot=case.valid_inputs.to_dict(),
        created_at=datetime.now(timezone.utc),
    )

    metadata = workflow_info.list_workflows()
    schema = workflow_info.get_schema(case.adapter.metadata.workflow_id)
    validated = validation.validate(
        case.adapter.metadata.workflow_id, case.valid_inputs
    )
    planned = workspace_planner.plan_workspace(
        execution_plan, tmp_path.resolve() / "service-workspace"
    )

    assert metadata == [case.adapter.metadata]
    assert schema.is_success and isinstance(schema.value, WorkflowSchema)
    assert validated.is_success
    assert planned.is_success
    assert planned.value.workspace_plan.files[0][0] == "config/inputs.json"


def _minimal_project(root: Path) -> Path:
    files = {
        "pyproject.toml": "[project]\nname='minimal'\n",
        "docs/architecture/artifact-inventory.yaml": "artifacts: []\n",
        "src/encode_pipeline/minimal.py": "VERSION = 1\n",
        "workflow/Snakefile": "rule all:\n    input: []\n",
        "profiles/default/config.yaml": "cores: 1\n",
        "scripts/minimal.py": "VALUE = 1\n",
    }
    for relative_path, content in files.items():
        destination = root / relative_path
        destination.parent.mkdir(parents=True, exist_ok=True)
        destination.write_text(content, encoding="utf-8")
    return root


def test_minimal_adapter_reuses_artifact_and_qc_platform_services(tmp_path):
    case = _minimal_case(tmp_path)
    registry = WorkflowRegistry((case.adapter,))
    run_service = RunService(registry, id_factory=lambda: "run-minimal")
    identity_provider = WorkflowBuildIdentityProvider(
        registry,
        project_root=_minimal_project(tmp_path / "minimal-project").resolve(),
    )
    run_service.create_run(case.adapter.metadata.workflow_id, case.valid_inputs)
    run_service.transition_run("run-minimal", RunStatus.VALIDATING)
    identity = identity_provider.capture(case.adapter.metadata.workflow_id)
    assert identity.is_success
    run_service.complete_preflight("run-minimal", identity.value)
    run_service.transition_run("run-minimal", RunStatus.QUEUED)
    run_service.transition_run("run-minimal", RunStatus.RUNNING)
    run_service.transition_run("run-minimal", RunStatus.SUCCEEDED)
    workspace_root = (tmp_path / "service-workspaces").resolve()
    workspace = workspace_root / "run-minimal"
    (workspace / "results").mkdir(parents=True)
    (workspace / "results/qc.tsv").write_bytes(b"metric\tvalue\nreads.total\t10\n")
    artifact_service = ArtifactExtractionService(
        run_service=run_service,
        registry=registry,
        build_identity_provider=identity_provider,
        workspace_root=workspace_root,
    )
    qc_service = QcSummaryIndexingService(
        run_service=run_service,
        registry=registry,
        build_identity_provider=identity_provider,
        workspace_root=workspace_root,
    )

    artifacts = artifact_service.extract("run-minimal")
    assert artifacts.is_success and len(artifacts.value) == 1
    metrics = qc_service.index("run-minimal", artifacts.value)

    assert metrics.is_success and len(metrics.value) == 1
    assert run_service.list_artifacts("run-minimal") == artifacts.value
    assert run_service.list_qc_metrics("run-minimal") == metrics.value
    assert metrics.value[0].metric_key == "reads.total"
