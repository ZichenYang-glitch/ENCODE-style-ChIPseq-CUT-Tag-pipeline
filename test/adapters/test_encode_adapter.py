"""Tests for the ENCODE-style workflow adapter wrapper."""

import os
import subprocess
import sys
import textwrap
from pathlib import Path

import pytest

from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    WorkflowAdapter,
    WorkflowInputs,
    WorkspacePlan,
)
from encode_pipeline.platform.results import Issue, Result

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"


VALID_SAMPLES = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
    "S1\tR1.fq\tR2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\tidx\n"
)


def _write_samples(path: Path, content: str = VALID_SAMPLES) -> Path:
    path.write_text(content, encoding="utf-8")
    return path


def _adapter() -> EncodeStyleWorkflowAdapter:
    return EncodeStyleWorkflowAdapter()


def _first_issue(result: Result[object]) -> Issue:
    assert result.errors
    return result.errors[0]


def test_adapter_structurally_satisfies_workflow_adapter():
    adapter = _adapter()

    assert isinstance(adapter, WorkflowAdapter)


def test_metadata_and_capabilities_match_minimal_contract():
    adapter = _adapter()

    assert adapter.metadata.workflow_id == "encode-style-chipseq-cuttag-atac-mnase"
    assert adapter.metadata.name == "ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase"
    assert adapter.metadata.version
    assert adapter.metadata.engines == ("snakemake",)
    assert adapter.metadata.tags == (
        "chipseq",
        "cuttag",
        "atac",
        "mnase",
        "encode-style",
    )
    assert adapter.capabilities.supports == ("validation",)


def test_schema_returns_json_ready_hints_and_strict_inputs_option():
    schema = _adapter().schema()
    as_dict = schema.to_dict()

    assert as_dict["config_schema"]["schema_kind"] == "hints"
    assert as_dict["sample_schema"]["schema_kind"] == "hints"
    assert "samples" in as_dict["config_schema"]["properties"]
    assert as_dict["sample_schema"]["required_columns"] == [
        "sample",
        "fastq_1",
        "layout",
        "assay",
        "target",
        "peak_mode",
        "genome",
        "bowtie2_index",
    ]
    assert as_dict["option_schema"]["properties"]["strict_inputs"] == {
        "type": "boolean",
        "default": False,
        "description": "Validate FASTQ and Bowtie2 index file existence.",
    }


def test_importing_encode_adapter_does_not_import_heavy_workflow_modules():
    code = """
        import sys
        import encode_pipeline.adapters.encode

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


def test_valid_config_plus_tsv_path_returns_successful_result(tmp_path):
    samples_path = _write_samples(tmp_path / "samples.tsv")
    inputs = WorkflowInputs(
        config={"samples": str(samples_path), "use_control": False},
    )

    result = _adapter().validate(inputs)

    assert result.is_success
    assert isinstance(result.value, dict)
    assert set(result.value) == {"config", "samples"}
    assert result.value["config"]["samples"] == str(samples_path)
    assert result.value["samples"][0]["id"] == "S1"
    assert result.issues == ()


def test_workflow_inputs_samples_path_can_override_config_samples(tmp_path):
    samples_path = _write_samples(tmp_path / "samples.tsv")
    inputs = WorkflowInputs(
        config={"samples": str(tmp_path / "missing.tsv"), "use_control": False},
        samples=samples_path,
    )

    result = _adapter().validate(inputs)

    assert result.is_success
    assert result.value["config"]["samples"] == str(samples_path)


def test_adapter_validate_output_consumable_by_platform(tmp_path):
    samples_path = _write_samples(tmp_path / "samples.tsv")
    result = _adapter().validate(
        WorkflowInputs(config={"samples": str(samples_path)}),
    )

    serialized = result.to_dict(
        value_serializer=lambda value: {
            "sample_count": len(value["samples"]),
            "has_config": "config" in value,
        },
    )

    assert result.is_success
    assert all(isinstance(issue, Issue) for issue in result.issues)
    assert serialized == {
        "ok": True,
        "value": {"sample_count": 1, "has_config": True},
        "issues": [],
    }


def test_invalid_config_maps_to_config_issue(tmp_path):
    samples_path = _write_samples(tmp_path / "samples.tsv")
    result = _adapter().validate(
        WorkflowInputs(config={"samples": str(samples_path), "threads": 0}),
    )

    issue = _first_issue(result)
    assert result.is_failure
    assert issue.code == "ENCODE_CONFIG_INVALID"
    assert issue.source == "config"
    assert issue.path == "config"


def test_invalid_sample_sheet_maps_to_samples_issue(tmp_path):
    samples_path = _write_samples(
        tmp_path / "samples.tsv",
        "sample\tfastq_2\tlayout\nS1\tR2.fq\tSE\n",
    )
    result = _adapter().validate(
        WorkflowInputs(config={"samples": str(samples_path)}),
    )

    issue = _first_issue(result)
    assert result.is_failure
    assert issue.code == "ENCODE_SAMPLES_INVALID"
    assert issue.source == "samples"
    assert issue.path == "samples"


def test_invalid_strict_inputs_type_maps_to_options_issue(tmp_path):
    samples_path = _write_samples(tmp_path / "samples.tsv")
    result = _adapter().validate(
        WorkflowInputs(
            config={"samples": str(samples_path)},
            options={"strict_inputs": "yes"},
        ),
    )

    issue = _first_issue(result)
    assert result.is_failure
    assert issue.code == "ENCODE_OPTIONS_INVALID"
    assert issue.source == "adapter"
    assert issue.path == "options.strict_inputs"
    assert issue.message == "strict_inputs must be a boolean"


def test_unsupported_options_return_options_issue(tmp_path):
    samples_path = _write_samples(tmp_path / "samples.tsv")
    result = _adapter().validate(
        WorkflowInputs(
            config={"samples": str(samples_path)},
            options={"profile": "local"},
        ),
    )

    issue = _first_issue(result)
    assert result.is_failure
    assert issue.code == "ENCODE_OPTIONS_INVALID"
    assert issue.source == "adapter"
    assert issue.path == "options"
    assert issue.context == {"unsupported_options": ["profile"]}


def test_unsupported_inline_sample_rows_return_adapter_issue():
    result = _adapter().validate(
        WorkflowInputs(
            config={},
            samples=[{"sample": "S1", "fastq_1": "R1.fq"}],
        ),
    )

    issue = _first_issue(result)
    assert result.is_failure
    assert issue.code == "ENCODE_ADAPTER_UNSUPPORTED"
    assert issue.source == "adapter"
    assert issue.path == "samples"
    assert issue.context == {"feature": "inline_samples"}


def test_unsupported_methods_return_failure_without_side_effects(tmp_path):
    adapter = _adapter()
    inputs = WorkflowInputs(config={})
    workspace = tmp_path / "workspace"

    dag_result = adapter.preview_dag(inputs)
    workspace_result = adapter.plan_workspace(inputs, workspace)
    command_result = adapter.build_command(
        WorkspacePlan(
            directories=[str(workspace)],
            files=[("config.yaml", b"samples: samples.tsv\n")],
        ),
    )

    assert not workspace.exists()
    for result, expected_method, expected_type in [
        (dag_result, "preview_dag", DagPreview),
        (workspace_result, "plan_workspace", WorkspacePlan),
        (command_result, "build_command", CommandSpec),
    ]:
        assert result.is_failure
        assert result.value is None
        assert result.errors[0].code == "ENCODE_ADAPTER_UNSUPPORTED"
        assert result.errors[0].source == "adapter"
        assert result.errors[0].path == expected_method
        assert result.errors[0].context == {"method": expected_method}
        assert not isinstance(result.value, expected_type)
