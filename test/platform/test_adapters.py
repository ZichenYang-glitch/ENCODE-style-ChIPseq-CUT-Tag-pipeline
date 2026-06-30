"""Tests for workflow-platform adapter contract primitives."""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagEdge,
    DagNode,
    DagPreview,
    WorkflowAdapter,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.results import Result


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"


@pytest.mark.parametrize(
    "kwargs",
    [
        {"workflow_id": "", "name": "Workflow", "version": "1.0.0"},
        {"workflow_id": "workflow", "name": "", "version": "1.0.0"},
        {"workflow_id": "workflow", "name": "Workflow", "version": ""},
        {"workflow_id": "   ", "name": "Workflow", "version": "1.0.0"},
        {"workflow_id": "workflow", "name": "   ", "version": "1.0.0"},
        {"workflow_id": "workflow", "name": "Workflow", "version": "   "},
    ],
)
def test_metadata_rejects_empty_required_fields(kwargs):
    with pytest.raises(ValueError):
        WorkflowMetadata(**kwargs)


def test_metadata_defaults_engine_and_normalizes_tuple_fields():
    defaulted = WorkflowMetadata(
        workflow_id="encode-chipseq",
        name="ENCODE ChIP-seq",
        version="1.0.0",
    )
    explicit = WorkflowMetadata(
        workflow_id="rnaseq",
        name="RNA-seq",
        version="2.0.0",
        engines=["nextflow"],
        tags=["rna", "expression"],
    )

    assert defaulted.engines == ("snakemake",)
    assert explicit.engines == ("nextflow",)
    assert explicit.tags == ("rna", "expression")


def test_capabilities_stores_tuple_supports_and_rejects_empty_strings():
    capabilities = WorkflowCapabilities(
        supports=["validation", "dag_preview"],
    )

    assert capabilities.supports == ("validation", "dag_preview")
    with pytest.raises(ValueError):
        WorkflowCapabilities(supports=["validation", "  "])


def test_schemas_deep_copy_nested_mappings_and_to_dict_returns_deep_copy():
    config_schema = {"properties": {"threads": {"type": "integer"}}}
    sample_schema = {"columns": ["sample", "fastq_1"]}
    option_schema = {"properties": {"trim": {"type": "boolean"}}}
    schema = WorkflowSchema(
        config_schema=config_schema,
        sample_schema=sample_schema,
        option_schema=option_schema,
    )

    config_schema["properties"]["threads"]["type"] = "string"
    sample_schema["columns"].append("fastq_2")
    option_schema["properties"]["trim"]["type"] = "string"
    serialized = schema.to_dict()
    serialized["config_schema"]["properties"]["threads"]["type"] = "number"
    serialized["sample_schema"]["columns"].append("layout")
    serialized["option_schema"]["properties"]["trim"]["type"] = "number"

    assert schema.config_schema == {
        "properties": {"threads": {"type": "integer"}},
    }
    assert schema.sample_schema == {"columns": ["sample", "fastq_1"]}
    assert schema.option_schema == {
        "properties": {"trim": {"type": "boolean"}},
    }
    assert schema.to_dict()["config_schema"] == {
        "properties": {"threads": {"type": "integer"}},
    }
    assert schema.to_dict()["sample_schema"] == {"columns": ["sample", "fastq_1"]}
    assert schema.to_dict()["option_schema"] == {
        "properties": {"trim": {"type": "boolean"}},
    }


def test_inputs_deep_copy_nested_config_options_and_to_dict_returns_deep_copy():
    config = {"samples": "samples.tsv", "resources": {"threads": 8}}
    options = {"execution": {"profile": "local"}}
    samples = [{"sample": "S1", "fastq_1": "R1.fq"}]
    inputs = WorkflowInputs(
        config=config,
        samples=samples,
        options=options,
    )

    config["resources"]["threads"] = 16
    options["execution"]["profile"] = "cluster"
    samples[0]["sample"] = "S2"
    serialized = inputs.to_dict()
    serialized["config"]["resources"]["threads"] = 32
    serialized["options"]["execution"]["profile"] = "cloud"
    serialized["samples"][0]["sample"] = "S3"

    assert inputs.config == {"samples": "samples.tsv", "resources": {"threads": 8}}
    assert inputs.options == {"execution": {"profile": "local"}}
    assert inputs.samples == [{"sample": "S1", "fastq_1": "R1.fq"}]
    assert inputs.to_dict()["config"] == {
        "samples": "samples.tsv",
        "resources": {"threads": 8},
    }
    assert inputs.to_dict()["options"] == {"execution": {"profile": "local"}}
    assert inputs.to_dict()["samples"] == [
        {"sample": "S1", "fastq_1": "R1.fq"},
    ]


def test_inputs_preserve_path_internally_and_serialize_path_as_string():
    sample_path = Path("samples.tsv")
    inputs = WorkflowInputs(config={"samples": str(sample_path)}, samples=sample_path)

    assert inputs.samples == sample_path
    assert inputs.to_dict()["samples"] == "samples.tsv"


@pytest.mark.parametrize("samples", [object(), {"sample": "S1"}, [object()]])
def test_inputs_reject_invalid_samples_object_type(samples):
    with pytest.raises(ValueError):
        WorkflowInputs(config={}, samples=samples)


def test_dag_preview_serializes_nodes_and_edges_without_tool_or_command_fields():
    preview = DagPreview(
        nodes=[
            DagNode(
                id="align",
                label="Align reads",
                inputs=["reads.fastq.gz"],
                outputs=["aligned.bam"],
            ),
            DagNode(id="qc", label="Summarize QC"),
        ],
        edges=[DagEdge(source="align", target="qc")],
    )

    assert preview.nodes[0].inputs == ("reads.fastq.gz",)
    assert preview.nodes[0].outputs == ("aligned.bam",)
    assert preview.to_dict() == {
        "nodes": [
            {
                "id": "align",
                "label": "Align reads",
                "inputs": ["reads.fastq.gz"],
                "outputs": ["aligned.bam"],
            },
            {
                "id": "qc",
                "label": "Summarize QC",
                "inputs": [],
                "outputs": [],
            },
        ],
        "edges": [{"source": "align", "target": "qc"}],
    }
    assert "tool" not in preview.to_dict()["nodes"][0]
    assert "command" not in preview.to_dict()["nodes"][0]


def test_command_spec_stores_argv_as_tuple_rejects_empty_argv_and_copies_env():
    env = {"SNAKEMAKE_PROFILE": "local"}
    command = CommandSpec(
        argv=["snakemake", "-n"],
        cwd="/workspace/run-1",
        env=env,
    )

    env["SNAKEMAKE_PROFILE"] = "cluster"

    assert command.argv == ("snakemake", "-n")
    assert command.env == {"SNAKEMAKE_PROFILE": "local"}
    assert command.to_dict() == {
        "argv": ["snakemake", "-n"],
        "cwd": "/workspace/run-1",
        "env": {"SNAKEMAKE_PROFILE": "local"},
    }
    with pytest.raises(ValueError):
        CommandSpec(argv=[])


def test_workspace_plan_stores_tuples_and_does_not_touch_filesystem(tmp_path):
    planned_dir = tmp_path / "planned"
    plan = WorkspacePlan(
        directories=[str(planned_dir)],
        files=[("config.yaml", b"samples: samples.tsv\n")],
    )

    assert plan.directories == (str(planned_dir),)
    assert plan.files == (("config.yaml", b"samples: samples.tsv\n"),)
    assert not planned_dir.exists()


def test_fake_adapter_satisfies_protocol_without_inheritance():
    class FakeAdapter:
        metadata = WorkflowMetadata(
            workflow_id="fake",
            name="Fake Workflow",
            version="1.0.0",
        )
        capabilities = WorkflowCapabilities(
            supports=("validation", "dag_preview", "workspace_plan", "command"),
        )

        def schema(self) -> WorkflowSchema:
            return WorkflowSchema(config_schema={"type": "object"})

        def validate(self, inputs: WorkflowInputs) -> Result[object]:
            return Result.success({"validated": dict(inputs.config)})

        def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
            return Result.success(
                DagPreview(nodes=[DagNode(id="step", label="Step")]),
            )

        def plan_workspace(
            self,
            inputs: WorkflowInputs,
            workspace: str | Path,
        ) -> Result[WorkspacePlan]:
            return Result.success(WorkspacePlan(directories=[str(workspace)]))

        def build_command(self, plan: WorkspacePlan) -> Result[CommandSpec]:
            return Result.success(CommandSpec(argv=["run-workflow"]))

    adapter = FakeAdapter()
    inputs = WorkflowInputs(config={"samples": "samples.tsv"})

    assert isinstance(adapter, WorkflowAdapter)
    assert isinstance(adapter.validate(inputs), Result)
    assert adapter.validate(inputs).value == {"validated": {"samples": "samples.tsv"}}
    assert isinstance(adapter.preview_dag(inputs).value, DagPreview)
    assert isinstance(adapter.plan_workspace(inputs, "/workspace").value, WorkspacePlan)
    assert isinstance(adapter.build_command(WorkspacePlan()).value, CommandSpec)


def test_platform_exports_adapter_contract_primitives():
    code = """
        from encode_pipeline.platform import (
            CommandSpec,
            DagEdge,
            DagNode,
            DagPreview,
            WorkflowAdapter,
            WorkflowCapabilities,
            WorkflowInputs,
            WorkflowMetadata,
            WorkflowSchema,
            WorkspacePlan,
        )
        print(all([
            CommandSpec,
            DagEdge,
            DagNode,
            DagPreview,
            WorkflowAdapter,
            WorkflowCapabilities,
            WorkflowInputs,
            WorkflowMetadata,
            WorkflowSchema,
            WorkspacePlan,
        ]))
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert proc.returncode == 0, proc.stderr
    assert proc.stdout.strip() == "True"


def test_importing_platform_adapters_does_not_import_workflow_specific_modules():
    code = """
        import sys
        import encode_pipeline.platform.adapters

        forbidden = [
            "encode_pipeline.config.validator",
            "encode_pipeline.samples",
            "encode_pipeline.cli.dag",
            "fastapi",
            "pydantic",
            "snakemake",
        ]
        for name in forbidden:
            print(f"{name}={name in sys.modules}")
    """
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC_ROOT)
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    proc = subprocess.run(
        [sys.executable, "-c", textwrap.dedent(code)],
        capture_output=True,
        text=True,
        env=env,
    )

    assert proc.returncode == 0, proc.stderr
    assert set(proc.stdout.splitlines()) == {
        "encode_pipeline.config.validator=False",
        "encode_pipeline.samples=False",
        "encode_pipeline.cli.dag=False",
        "fastapi=False",
        "pydantic=False",
        "snakemake=False",
    }
