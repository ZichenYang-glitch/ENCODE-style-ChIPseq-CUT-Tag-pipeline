"""Tests for workflow-platform adapter contract primitives."""

import os
import subprocess
import sys
import textwrap
from dataclasses import FrozenInstanceError
from decimal import Decimal
from math import nan
from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagEdge,
    DagNode,
    DagPreview,
    ExtractedArtifactCandidate,
    ExtractedQcMetricCandidate,
    LocalRunDriver,
    QcSourceArtifact,
    QcSourceDocument,
    QcSummaryExtractingAdapter,
    WorkflowAdapter,
    WorkflowAuthoringModes,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowInputLimits,
    WorkflowInputModes,
    WorkflowMetadata,
    WorkflowSchema,
    WorkflowSchemaCoverage,
    WorkspacePlan,
    JSON_SCHEMA_DIALECT,
    MAX_AUTHORING_REQUEST_BYTES,
    MAX_SAMPLE_CELL_LENGTH,
    MAX_SAMPLE_COLUMN_NAME_LENGTH,
    MAX_SAMPLE_COLUMNS,
    MAX_SAMPLE_ROWS,
)
from encode_pipeline.platform.results import Result


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"

INPUT_LIMIT_FIELDS = (
    ("max_request_bytes", MAX_AUTHORING_REQUEST_BYTES),
    ("max_sample_rows", MAX_SAMPLE_ROWS),
    ("max_sample_columns", MAX_SAMPLE_COLUMNS),
    ("max_sample_column_name_length", MAX_SAMPLE_COLUMN_NAME_LENGTH),
    ("max_sample_cell_length", MAX_SAMPLE_CELL_LENGTH),
)


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


def test_schema_is_top_level_frozen_and_copies_input_and_serialized_mappings():
    config_schema = {"properties": {"threads": {"type": "integer"}}}
    sample_schema = {
        "type": "array",
        "items": {"type": "object", "properties": {"sample": {"type": "string"}}},
    }
    option_schema = {"properties": {"trim": {"type": "boolean"}}}
    schema = WorkflowSchema(
        config_schema=config_schema,
        sample_schema=sample_schema,
        option_schema=option_schema,
    )

    config_schema["properties"]["threads"]["type"] = "string"
    sample_schema["items"]["properties"]["sample"]["type"] = "integer"
    option_schema["properties"]["trim"]["type"] = "string"
    serialized = schema.to_dict()
    serialized["config_schema"]["properties"]["threads"]["type"] = "number"
    serialized["sample_schema"]["items"]["properties"]["sample"]["type"] = "number"
    serialized["option_schema"]["properties"]["trim"]["type"] = "number"

    with pytest.raises(FrozenInstanceError):
        schema.schema_version = "2.0.0"

    assert schema.config_schema == {
        "$schema": JSON_SCHEMA_DIALECT,
        "type": "object",
        "properties": {"threads": {"type": "integer"}},
    }
    assert schema.sample_schema == {
        "$schema": JSON_SCHEMA_DIALECT,
        "type": "array",
        "items": {
            "type": "object",
            "properties": {"sample": {"type": "string"}},
        },
    }
    assert schema.option_schema == {
        "$schema": JSON_SCHEMA_DIALECT,
        "type": "object",
        "properties": {"trim": {"type": "boolean"}},
    }
    assert schema.to_dict()["config_schema"] == {
        "$schema": JSON_SCHEMA_DIALECT,
        "type": "object",
        "properties": {"threads": {"type": "integer"}},
    }
    assert schema.to_dict()["sample_schema"] == {
        "$schema": JSON_SCHEMA_DIALECT,
        "type": "array",
        "items": {
            "type": "object",
            "properties": {"sample": {"type": "string"}},
        },
    }
    assert schema.to_dict()["option_schema"] == {
        "$schema": JSON_SCHEMA_DIALECT,
        "type": "object",
        "properties": {"trim": {"type": "boolean"}},
    }


def test_input_limits_default_serialization_matches_platform_ceilings():
    assert WorkflowInputLimits().to_dict() == dict(INPUT_LIMIT_FIELDS)


@pytest.mark.parametrize(("field_name", "ceiling"), INPUT_LIMIT_FIELDS)
@pytest.mark.parametrize("offset", (-1, 1), ids=("lower", "upper"))
def test_input_limits_reject_values_below_or_above_platform_ceiling(
    field_name: str,
    ceiling: int,
    offset: int,
):
    with pytest.raises(ValueError, match="must equal the platform ceiling"):
        WorkflowInputLimits(**{field_name: ceiling + offset})


@pytest.mark.parametrize(("field_name", "ceiling"), INPUT_LIMIT_FIELDS)
@pytest.mark.parametrize(
    "invalid_value",
    (True, "100"),
    ids=("boolean", "string"),
)
def test_input_limits_reject_wrong_types(
    field_name: str,
    ceiling: int,
    invalid_value: object,
):
    del ceiling
    with pytest.raises(ValueError, match="must be an integer"):
        WorkflowInputLimits(**{field_name: invalid_value})


def test_schema_serializes_versioned_authoring_contract_values():
    coverage = WorkflowSchemaCoverage(
        config="partial",
        samples="complete",
        options="complete",
    )
    authoring_modes = WorkflowAuthoringModes(
        config=("schema_form", "yaml"),
        samples=("tsv_upload", "inline_table"),
        options=("schema_form",),
    )
    input_modes = WorkflowInputModes(
        config=("object",),
        samples=("inline_rows", "server_path"),
        options=("object",),
    )
    limits = WorkflowInputLimits()
    schema = WorkflowSchema(
        schema_version="1.0.0",
        schema_dialect=JSON_SCHEMA_DIALECT,
        coverage=coverage,
        authoring_modes=authoring_modes,
        input_modes=input_modes,
        limits=limits,
    )

    serialized = schema.to_dict()
    serialized["coverage"]["config"] = "complete"
    serialized["input_modes"]["samples"].append("unsafe")

    assert schema.schema_version == "1.0.0"
    assert schema.schema_dialect == JSON_SCHEMA_DIALECT
    assert schema.coverage.config == "partial"
    assert schema.input_modes.samples == ("inline_rows", "server_path")
    assert schema.to_dict()["authoring_modes"] == {
        "config": ["schema_form", "yaml"],
        "samples": ["tsv_upload", "inline_table"],
        "options": ["schema_form"],
    }
    assert schema.to_dict()["limits"] == {
        "max_request_bytes": MAX_AUTHORING_REQUEST_BYTES,
        "max_sample_rows": MAX_SAMPLE_ROWS,
        "max_sample_columns": MAX_SAMPLE_COLUMNS,
        "max_sample_column_name_length": MAX_SAMPLE_COLUMN_NAME_LENGTH,
        "max_sample_cell_length": MAX_SAMPLE_CELL_LENGTH,
    }


@pytest.mark.parametrize(
    "factory",
    [
        lambda: WorkflowSchemaCoverage(config="unknown"),
        lambda: WorkflowAuthoringModes(config=("",)),
        lambda: WorkflowInputModes(samples=("inline rows",)),
        lambda: WorkflowSchema(schema_version="not-a-version"),
        lambda: WorkflowSchema(
            schema_dialect="http://json-schema.org/draft-07/schema#"
        ),
        lambda: WorkflowSchema(config_schema={"type": "array"}),
        lambda: WorkflowSchema(sample_schema={"type": "object"}),
        lambda: WorkflowSchema(option_schema={"type": "array"}),
        lambda: WorkflowSchema(config_schema={"type": "object", "value": nan}),
        lambda: WorkflowSchema(
            config_schema={"type": "object", "value": Path("local")}
        ),
    ],
)
def test_schema_contract_rejects_invalid_or_non_json_values(factory):
    with pytest.raises(ValueError):
        factory()


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


@pytest.mark.parametrize(
    "samples",
    [
        [],
        [{f"column-{index}": "value" for index in range(MAX_SAMPLE_COLUMNS + 1)}],
        [{"x" * (MAX_SAMPLE_COLUMN_NAME_LENGTH + 1): "value"}],
        [{"sample\tname": "value"}],
        [{"sample": "x" * (MAX_SAMPLE_CELL_LENGTH + 1)}],
        [{"sample": "contains\x00nul"}],
        [{"sample": "contains\ttab"}],
        [{"sample": "contains\nnewline"}],
        [{"sample": "contains\rcarriage-return"}],
    ],
)
def test_inputs_reject_empty_or_structurally_unsafe_inline_rows(samples):
    with pytest.raises(ValueError):
        WorkflowInputs(config={}, samples=samples)


def test_inputs_reject_too_many_inline_rows():
    with pytest.raises(ValueError):
        WorkflowInputs(
            config={},
            samples=[{"sample": f"S{index}"} for index in range(MAX_SAMPLE_ROWS + 1)],
        )


@pytest.mark.parametrize(
    "samples",
    [
        [{"sample": f"S{index}"} for index in range(MAX_SAMPLE_ROWS)],
        [{f"column-{index}": "value" for index in range(MAX_SAMPLE_COLUMNS)}],
        [{"x" * MAX_SAMPLE_COLUMN_NAME_LENGTH: "value"}],
        [{"sample": "x" * MAX_SAMPLE_CELL_LENGTH}],
    ],
    ids=("rows", "columns", "column-name", "cell"),
)
def test_inputs_accept_each_exact_platform_ceiling(samples):
    inputs = WorkflowInputs(config={}, samples=samples)

    assert len(inputs.samples) == len(samples)


@pytest.mark.parametrize("samples", ["bad\x00path", "bad\npath", Path("bad\tpath")])
def test_inputs_reject_control_characters_in_sample_paths(samples):
    with pytest.raises(ValueError):
        WorkflowInputs(config={}, samples=samples)


def test_extracted_artifact_candidate_defensively_copies_metadata():
    metadata = {"catalog_id": "cpm_bigwig", "nested": {"scope": "sample"}}
    candidate = ExtractedArtifactCandidate(
        output_type="cpm_bigwig",
        relative_path="results/S1/03_bigwig/S1.CPM.bw",
        mime_type="application/octet-stream",
        metadata=metadata,
    )

    metadata["catalog_id"] = "changed"
    metadata["nested"]["scope"] = "project"

    assert candidate.metadata == {
        "catalog_id": "cpm_bigwig",
        "nested": {"scope": "sample"},
    }


def test_qc_source_values_are_immutable_and_defensively_copy_metadata():
    metadata = {"scope": "sample", "sample_id": "S1"}
    source = QcSourceArtifact(
        artifact_id="artifact-1",
        output_type="qc_summary",
        relative_path="results/S1/01_qc/S1.qc_summary.tsv",
        metadata=metadata,
    )
    document = QcSourceDocument(source=source, content=b"sample\tfrip\nS1\t0.5\n")
    candidate = ExtractedQcMetricCandidate(
        metric_key="peaks.frip",
        display_name="Fraction of reads in peaks",
        value=Decimal("0.5"),
        unit="fraction",
        scope="sample",
        sample_id="S1",
        experiment_id="EXP1",
        assay="chipseq",
        qc_flag=None,
        source_artifact_id="artifact-1",
    )

    metadata["sample_id"] = "changed"

    assert source.metadata == {"scope": "sample", "sample_id": "S1"}
    assert document.content == b"sample\tfrip\nS1\t0.5\n"
    assert candidate.value == Decimal("0.5")
    with pytest.raises(ValueError):
        QcSourceDocument(source=source, content=bytearray(b"unsafe"))
    with pytest.raises(ValueError):
        ExtractedQcMetricCandidate(
            metric_key="peaks.frip",
            display_name="FRiP",
            value=0.5,
            unit="fraction",
            scope="sample",
            sample_id="S1",
            source_artifact_id="artifact-1",
        )


def test_optional_qc_adapter_protocol_is_runtime_checkable():
    class QcAdapter:
        def qc_source_output_types(self):
            return ("qc_summary",)

        def extract_qc_metrics(self, inputs, sources):
            return Result.success(())

    class NoQcAdapter:
        pass

    adapter = QcAdapter()
    assert isinstance(adapter, QcSummaryExtractingAdapter)
    assert not isinstance(NoQcAdapter(), QcSummaryExtractingAdapter)
    assert adapter.extract_qc_metrics(WorkflowInputs(config={}), ()).value == ()


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

        def extract_artifacts(self, inputs, workspace):
            return Result.success(
                (
                    ExtractedArtifactCandidate(
                        output_type="summary",
                        relative_path="results/summary.tsv",
                    ),
                )
            )

    adapter = FakeAdapter()
    inputs = WorkflowInputs(config={"samples": "samples.tsv"})

    assert isinstance(adapter, WorkflowAdapter)
    assert isinstance(adapter.validate(inputs), Result)
    assert adapter.validate(inputs).value == {"validated": {"samples": "samples.tsv"}}
    assert isinstance(adapter.preview_dag(inputs).value, DagPreview)
    assert isinstance(adapter.plan_workspace(inputs, "/workspace").value, WorkspacePlan)
    assert isinstance(adapter.build_command(WorkspacePlan()).value, CommandSpec)
    assert (
        adapter.extract_artifacts(inputs, "/workspace").value[0].output_type
        == "summary"
    )


def test_platform_exports_adapter_contract_primitives():
    code = """
        from encode_pipeline.platform import (
            CommandSpec,
            DagEdge,
            DagNode,
            DagPreview,
            ExtractedArtifactCandidate,
            ExtractedQcMetricCandidate,
            QcSourceArtifact,
            QcSourceDocument,
            QcSummaryExtractingAdapter,
            JSON_SCHEMA_DIALECT,
            MAX_AUTHORING_REQUEST_BYTES,
            MAX_SAMPLE_CELL_LENGTH,
            MAX_SAMPLE_COLUMN_NAME_LENGTH,
            MAX_SAMPLE_COLUMNS,
            MAX_SAMPLE_ROWS,
            WorkflowAdapter,
            WorkflowAuthoringModes,
            WorkflowCapabilities,
            WorkflowInputLimits,
            WorkflowInputModes,
            WorkflowInputs,
            WorkflowMetadata,
            WorkflowSchema,
            WorkflowSchemaCoverage,
            WorkspacePlan,
        )
        print(all([
            CommandSpec,
            DagEdge,
            DagNode,
            DagPreview,
            ExtractedArtifactCandidate,
            ExtractedQcMetricCandidate,
            QcSourceArtifact,
            QcSourceDocument,
            QcSummaryExtractingAdapter,
            JSON_SCHEMA_DIALECT,
            MAX_AUTHORING_REQUEST_BYTES,
            MAX_SAMPLE_CELL_LENGTH,
            MAX_SAMPLE_COLUMN_NAME_LENGTH,
            MAX_SAMPLE_COLUMNS,
            MAX_SAMPLE_ROWS,
            WorkflowAdapter,
            WorkflowAuthoringModes,
            WorkflowCapabilities,
            WorkflowInputLimits,
            WorkflowInputModes,
            WorkflowInputs,
            WorkflowMetadata,
            WorkflowSchema,
            WorkflowSchemaCoverage,
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


def test_local_run_driver_protocol_is_runtime_checkable():
    class _Dummy:
        def run(self, run_id: str, plan: object) -> object:
            return None

    assert isinstance(_Dummy(), LocalRunDriver)
