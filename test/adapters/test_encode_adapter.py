"""Tests for the ENCODE-style workflow adapter wrapper."""

import csv
from concurrent.futures import ThreadPoolExecutor
import os
import subprocess
import sys
import tempfile
import textwrap
from io import StringIO
from pathlib import Path

import jsonschema
import pytest
import yaml

import encode_pipeline.adapters.encode as encode_adapter_module
from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    JSON_SCHEMA_DIALECT,
    WorkflowAdapter,
    WorkflowInputs,
    WorkspacePlan,
    QcSummaryExtractingAdapter,
)
from encode_pipeline.platform.results import Issue, Result

from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter


SRC_ROOT = Path(__file__).resolve().parents[2] / "src"


def _inline_row(root: Path, sample: str = "S1") -> dict[str, str]:
    return {
        "sample": sample,
        "fastq_1": str((root / f"{sample}.R1.fastq.gz").resolve()),
        "fastq_2": str((root / f"{sample}.R2.fastq.gz").resolve()),
        "layout": "PE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
        "genome": "hs",
        "bowtie2_index": str((root / "indices" / "hs").resolve()),
    }


def _samples_text(root: Path, sample: str = "S1") -> str:
    row = _inline_row(root, sample)
    output = StringIO()
    writer = csv.DictWriter(
        output,
        fieldnames=tuple(row),
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerow(row)
    return output.getvalue()


def _write_samples(path: Path, content: str | None = None) -> Path:
    if content is None:
        content = _samples_text(path.parent)
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
    assert isinstance(adapter, QcSummaryExtractingAdapter)


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
    assert adapter.capabilities.supports == (
        "validation",
        "workspace_plan",
        "input_authoring",
        "input_bundle_import",
        "artifact_extract",
        "qc_summary_extract",
    )


def test_qc_source_output_types_are_exact_and_do_not_include_aggregate_or_html():
    assert _adapter().qc_source_output_types() == (
        "mnase_qc_summary",
        "pooled_qc_summary",
        "qc_summary",
    )


def test_schema_returns_versioned_renderable_contract_and_strict_options():
    schema = _adapter().schema()
    as_dict = schema.to_dict()

    assert as_dict["schema_version"] == "1.1.0"
    assert as_dict["schema_dialect"] == JSON_SCHEMA_DIALECT
    assert as_dict["coverage"] == {
        "config": "partial",
        "samples": "complete",
        "options": "complete",
    }
    assert as_dict["authoring_modes"] == {
        "config": ["schema_form", "yaml"],
        "samples": ["tsv_upload", "inline_table"],
        "options": ["schema_form"],
    }
    assert as_dict["input_modes"] == {
        "config": ["object"],
        "samples": ["inline_rows", "server_path"],
        "options": ["object"],
    }
    assert "samples" not in as_dict["config_schema"].get("required", [])
    assert "samples" not in as_dict["config_schema"]["properties"]
    config_properties = as_dict["config_schema"]["properties"]
    assert "stage4b" not in config_properties
    assert "stage5" not in config_properties
    assert config_properties["replicate_analysis"] == {
        "type": "object",
        "title": "Replicate analysis",
        "description": (
            "Enable replicate-aware pooled outputs for experiments with "
            "biological replicates."
        ),
        "default": {"enabled": True},
        "properties": {
            "enabled": {
                "type": "boolean",
                "title": "Replicate analysis enabled",
                "default": True,
            }
        },
        "required": ["enabled"],
        "additionalProperties": False,
    }
    assert config_properties["chipseq_idr"] == {
        "type": "object",
        "title": "ChIP-seq IDR",
        "description": (
            "Enable narrow-peak ChIP-seq IDR for eligible experiments with "
            "exactly two biological replicates."
        ),
        "default": {"enabled": False},
        "properties": {
            "enabled": {
                "type": "boolean",
                "title": "ChIP-seq IDR enabled",
                "default": False,
            }
        },
        "required": ["enabled"],
        "additionalProperties": False,
    }
    assert as_dict["config_schema"]["properties"]["outdir"] == {
        "type": "string",
        "const": "results",
        "default": "results",
        "description": "Platform-owned run output directory.",
    }
    sample_items = as_dict["sample_schema"]["items"]
    assert sample_items["required"] == [
        "sample",
        "fastq_1",
        "layout",
        "assay",
        "target",
        "peak_mode",
        "genome",
        "bowtie2_index",
    ]
    assert sample_items["additionalProperties"] is False
    assert len(sample_items["properties"]) == 17
    assert as_dict["option_schema"]["properties"]["strict_inputs"] == {
        "type": "boolean",
        "default": False,
        "description": "Validate FASTQ and Bowtie2 index file existence.",
    }
    assert as_dict["option_schema"]["additionalProperties"] is False

    for name in ("config_schema", "sample_schema", "option_schema"):
        document = as_dict[name]
        assert document["$schema"] == JSON_SCHEMA_DIALECT
        assert document["$id"].endswith("/1.1.0")
        jsonschema.Draft202012Validator.check_schema(document)


@pytest.mark.parametrize(
    ("replicate_enabled", "chipseq_idr_enabled"),
    [(False, False), (False, True), (True, False), (True, True)],
)
def test_semantic_switches_translate_to_complete_engine_boolean_matrix(
    replicate_enabled,
    chipseq_idr_enabled,
):
    source = {
        "threads": 2,
        "replicate_analysis": {"enabled": replicate_enabled},
        "chipseq_idr": {"enabled": chipseq_idr_enabled},
    }

    result = encode_adapter_module._translate_authoring_config(source)

    assert result.is_success
    assert result.value == {
        "threads": 2,
        "stage4b": replicate_enabled,
        "stage5": chipseq_idr_enabled,
    }
    assert result.issues == ()
    assert source == {
        "threads": 2,
        "replicate_analysis": {"enabled": replicate_enabled},
        "chipseq_idr": {"enabled": chipseq_idr_enabled},
    }


@pytest.mark.parametrize("legacy", [True, "true", "TRUE", "True"])
def test_equal_legacy_true_spellings_use_existing_boolean_compatibility(legacy):
    result = encode_adapter_module._translate_authoring_config(
        {"replicate_analysis": {"enabled": True}, "stage4b": legacy}
    )

    assert result.is_success
    assert result.value == {"stage4b": True, "stage5": False}
    assert [issue.code for issue in result.issues] == [
        "ENCODE_CONFIG_LEGACY_ALIAS_DEPRECATED"
    ]


@pytest.mark.parametrize("legacy", [False, "false", "FALSE", "False"])
def test_equal_legacy_false_spellings_use_existing_boolean_compatibility(legacy):
    result = encode_adapter_module._translate_authoring_config(
        {"chipseq_idr": {"enabled": False}, "stage5": legacy}
    )

    assert result.is_success
    assert result.value == {"stage4b": True, "stage5": False}
    assert [issue.code for issue in result.issues] == [
        "ENCODE_CONFIG_LEGACY_ALIAS_DEPRECATED"
    ]


@pytest.mark.parametrize(
    ("semantic_key", "value"),
    [
        ("replicate_analysis", None),
        ("replicate_analysis", True),
        ("replicate_analysis", {}),
        ("replicate_analysis", {"enabled": "true"}),
        ("replicate_analysis", {"enabled": None}),
        ("replicate_analysis", {"enabled": True, "unknown": False}),
        ("chipseq_idr", None),
        ("chipseq_idr", {"enabled": "false"}),
    ],
)
def test_malformed_semantic_switches_fail_closed(semantic_key, value):
    result = encode_adapter_module._translate_authoring_config({semantic_key: value})

    assert result.is_failure
    assert result.errors[0].code == "ENCODE_CONFIG_SEMANTIC_INVALID"
    assert result.errors[0].path == f"config.{semantic_key}"
    assert "stage4b" not in str(result.errors[0].to_dict())
    assert "stage5" not in str(result.errors[0].to_dict())


@pytest.mark.parametrize(
    "config",
    [
        {"replicate_analysis": {"enabled": True}, "stage4b": False},
        {"replicate_analysis": {"enabled": False}, "stage4b": "TRUE"},
        {"chipseq_idr": {"enabled": True}, "stage5": "false"},
        {"chipseq_idr": {"enabled": False}, "stage5": True},
    ],
)
def test_conflicting_semantic_and_legacy_switches_fail_closed(config):
    result = encode_adapter_module._translate_authoring_config(config)

    assert result.is_failure
    issue = result.errors[0]
    assert issue.code == "ENCODE_CONFIG_SEMANTIC_CONFLICT"
    assert issue.path in {
        "config.replicate_analysis.enabled",
        "config.chipseq_idr.enabled",
    }
    assert "stage4b" not in str(issue.to_dict())
    assert "stage5" not in str(issue.to_dict())


@pytest.mark.parametrize("legacy", [None, 0, 1, "yes", " true "])
def test_invalid_legacy_value_is_not_misreported_as_semantic_conflict(
    tmp_path,
    legacy,
):
    result = _adapter().validate(
        WorkflowInputs(
            config={
                "replicate_analysis": {"enabled": True},
                "stage4b": legacy,
            },
            samples=[_inline_row(tmp_path)],
        )
    )

    assert result.is_failure
    assert result.errors[0].code == "ENCODE_CONFIG_INVALID"


def test_equal_aliases_emit_one_redacted_warning_and_remain_valid(tmp_path):
    result = _adapter().validate(
        WorkflowInputs(
            config={
                "replicate_analysis": {"enabled": False},
                "stage4b": "FALSE",
                "chipseq_idr": {"enabled": False},
                "stage5": False,
            },
            samples=[_inline_row(tmp_path)],
        )
    )

    assert result.is_success
    assert result.value["config"]["stage4b"] is False
    assert result.value["config"]["stage5"] is False
    assert [issue.code for issue in result.issues] == [
        "ENCODE_CONFIG_LEGACY_ALIAS_DEPRECATED"
    ]
    warning = result.issues[0]
    assert warning.severity.value == "warning"
    assert warning.path == "config"
    assert warning.technical_message is None
    assert warning.context == {}
    assert "stage4b" not in str(warning.to_dict())
    assert "stage5" not in str(warning.to_dict())


def test_legacy_only_switches_remain_supported_without_semantic_warning():
    result = encode_adapter_module._translate_authoring_config(
        {"stage4b": "false", "stage5": "true"}
    )

    assert result.is_success
    assert result.value == {"stage4b": "false", "stage5": "true"}
    assert result.issues == ()


def test_chipseq_idr_remains_orthogonal_to_advanced_reproducibility_policy(
    tmp_path,
):
    rows = []
    for index in (1, 2):
        row = _inline_row(tmp_path, f"S{index}")
        row.update(
            {
                "experiment": "EXP1",
                "biological_replicate": str(index),
            }
        )
        rows.append(row)

    with pytest.warns(UserWarning, match="Config contradiction"):
        result = _adapter().validate(
            WorkflowInputs(
                config={
                    "replicate_analysis": {"enabled": True},
                    "chipseq_idr": {"enabled": True},
                    "reproducibility": {
                        "enabled": True,
                        "idr": {"chipseq_narrow": False},
                    },
                },
                samples=rows,
            )
        )

    assert result.is_success
    assert result.value["config"]["stage4b"] is True
    assert result.value["config"]["stage5"] is True
    assert result.value["config"]["reproducibility"]["idr"]["chipseq_narrow"] is False


def test_schema_returns_fresh_instances_after_internal_mapping_mutation():
    adapter = _adapter()
    first = adapter.schema()
    first_properties = first.config_schema["properties"]
    assert isinstance(first_properties, dict)
    first_outdir = first_properties["outdir"]
    assert isinstance(first_outdir, dict)
    first_outdir["description"] = "mutated caller value"

    second = adapter.schema()
    second_properties = second.config_schema["properties"]
    assert isinstance(second_properties, dict)
    second_outdir = second_properties["outdir"]
    assert isinstance(second_outdir, dict)

    assert first is not second
    assert second_outdir["description"] == "Platform-owned run output directory."


def test_schema_accepts_tiny_profile_without_samples_and_inline_rows(tmp_path):
    schema = _adapter().schema()
    tiny_config = yaml.safe_load(
        (
            SRC_ROOT.parents[0] / "test/profiles/platform_worker_tiny/config.yaml"
        ).read_text(encoding="utf-8")
    )
    tiny_config.pop("samples")
    row = _inline_row(tmp_path)

    jsonschema.Draft202012Validator(schema.config_schema).validate(tiny_config)
    jsonschema.Draft202012Validator(schema.sample_schema).validate([row])
    jsonschema.Draft202012Validator(schema.option_schema).validate(
        {"strict_inputs": False}
    )
    result = _adapter().validate(
        WorkflowInputs(
            config=tiny_config, samples=[row], options={"strict_inputs": False}
        )
    )

    assert result.is_success
    assert result.value["samples"][0]["id"] == "S1"
    assert "samples" not in result.value["config"]


def test_schema_does_not_claim_custom_outdir_is_authorable():
    schema = _adapter().schema()

    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(schema.config_schema).validate(
            {"outdir": "custom-results"}
        )


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


def test_inline_sample_rows_validate_without_config_samples(tmp_path):
    result = _adapter().validate(
        WorkflowInputs(
            config={"use_control": False},
            samples=[_inline_row(tmp_path)],
        ),
    )

    assert result.is_success
    assert result.value["samples"][0]["id"] == "S1"
    assert "samples" not in result.value["config"]


def test_inline_rows_override_config_samples_without_fallback(tmp_path):
    missing_server_path = (tmp_path / "must-not-be-read.tsv").resolve()
    result = _adapter().validate(
        WorkflowInputs(
            config={"samples": str(missing_server_path), "use_control": False},
            samples=[_inline_row(tmp_path)],
        )
    )

    assert result.is_success
    assert not missing_server_path.exists()
    assert "samples" not in result.value["config"]


def test_inline_unknown_column_fails_closed_without_temp_path(tmp_path):
    row = _inline_row(tmp_path)
    row["adapter_private_column"] = "not-allowed"

    result = _adapter().validate(WorkflowInputs(config={}, samples=[row]))

    issue = _first_issue(result)
    assert result.is_failure
    assert issue.code == "ENCODE_SAMPLES_INVALID"
    assert issue.path == "samples"
    assert tempfile.gettempdir() not in str(result.to_dict())


def test_inline_scientific_failure_is_stable_and_does_not_leak_temp_path(tmp_path):
    row = _inline_row(tmp_path)
    row["layout"] = "INVALID"

    result = _adapter().validate(WorkflowInputs(config={}, samples=[row]))

    issue = _first_issue(result)
    assert result.is_failure
    assert issue.code == "ENCODE_SAMPLES_INVALID"
    assert issue.message == "Sample sheet is invalid."
    assert issue.technical_message == "Sample sheet is invalid."
    assert tempfile.gettempdir() not in str(result.to_dict())


def test_validate_reuses_external_path_policy_for_inline_rows(tmp_path):
    row = _inline_row(tmp_path)
    row["fastq_1"] = "relative.fastq.gz"

    result = _adapter().validate(WorkflowInputs(config={}, samples=[row]))

    issue = _first_issue(result)
    assert result.is_failure
    assert issue.code == "ENCODE_WORKSPACE_RELATIVE_EXTERNAL_PATH"
    assert issue.path == "samples[0].fastq_1"
    assert "relative.fastq.gz" not in str(issue.to_dict())


def test_inline_temp_directories_are_unique_and_cleaned_for_concurrent_validation(
    tmp_path, monkeypatch
):
    from encode_pipeline.adapters import encode as encode_module

    real_temporary_directory = tempfile.TemporaryDirectory
    roots: list[Path] = []

    def tracking_temporary_directory(*args, **kwargs):
        directory = real_temporary_directory(*args, **kwargs)
        roots.append(Path(directory.name))
        return directory

    monkeypatch.setattr(
        encode_module.tempfile,
        "TemporaryDirectory",
        tracking_temporary_directory,
    )
    adapter = _adapter()

    def validate(sample: str):
        return adapter.validate(
            WorkflowInputs(config={}, samples=[_inline_row(tmp_path, sample)])
        )

    with ThreadPoolExecutor(max_workers=2) as executor:
        results = list(executor.map(validate, ("S1", "S2")))

    assert all(result.is_success for result in results)
    assert {result.value["samples"][0]["id"] for result in results} == {"S1", "S2"}
    assert len(roots) == 2
    assert roots[0] != roots[1]
    assert all(not root.exists() for root in roots)


def test_unsupported_methods_return_failure_without_side_effects(tmp_path):
    adapter = _adapter()
    inputs = WorkflowInputs(config={})
    workspace = tmp_path / "workspace"

    dag_result = adapter.preview_dag(inputs)
    command_result = adapter.build_command(
        WorkspacePlan(
            directories=[str(workspace)],
            files=[("config.yaml", b"samples: samples.tsv\n")],
        ),
    )

    # plan_workspace is implemented; with invalid inputs it fails validation
    # without touching the filesystem.
    workspace_result = adapter.plan_workspace(inputs, workspace)

    assert not workspace.exists()
    for result, expected_method, expected_type in [
        (dag_result, "preview_dag", DagPreview),
        (command_result, "build_command", CommandSpec),
    ]:
        assert result.is_failure
        assert result.value is None
        assert result.errors[0].code == "ENCODE_ADAPTER_UNSUPPORTED"
        assert result.errors[0].source == "adapter"
        assert result.errors[0].path == expected_method
        assert result.errors[0].context == {"method": expected_method}
        assert not isinstance(result.value, expected_type)

    assert workspace_result.is_failure
    assert workspace_result.errors[0].code != "ENCODE_ADAPTER_UNSUPPORTED"


def _artifact_workspace(tmp_path: Path) -> Path:
    workspace = tmp_path / "workspace"
    (workspace / "config").mkdir(parents=True)
    (workspace / "results").mkdir()
    (workspace / "config/config.yaml").write_text("samples: config/samples.tsv\n")
    (workspace / "config/samples.tsv").write_text("unused\n")
    return workspace


def test_extract_artifacts_maps_present_and_dynamic_manifest_rows(
    tmp_path, monkeypatch
):
    workspace = _artifact_workspace(tmp_path)
    first = workspace / "results/S1/03_bigwig/S1.CPM.bw"
    second = workspace / "results/experiments/E1/02_align/biorep1.final.bam"
    first.parent.mkdir(parents=True)
    second.parent.mkdir(parents=True)
    first.write_bytes(b"bw")
    second.write_bytes(b"bam")
    rows = [
        {
            "output_type": "cpm_bigwig",
            "path": str(first),
            "status": "present",
            "sample_id": "S1",
        },
        {
            "output_type": "biorep1_final_bam",
            "path": str(second),
            "status": "present",
            "experiment_id": "E1",
        },
        {
            "output_type": "final_bam",
            "path": str(workspace / "results/missing.bam"),
            "status": "missing",
        },
        {"output_type": "multiqc_report", "path": "", "status": "not_applicable"},
    ]
    monkeypatch.setattr(
        "encode_pipeline.manifest.make.build_manifest_rows",
        lambda *_args, **_kwargs: (rows, 1, 1),
    )

    result = _adapter().extract_artifacts(WorkflowInputs(config={}), workspace)

    assert result.is_success
    assert [(item.output_type, item.relative_path) for item in result.value] == [
        ("biorep1_final_bam", "results/experiments/E1/02_align/biorep1.final.bam"),
        ("cpm_bigwig", "results/S1/03_bigwig/S1.CPM.bw"),
    ]
    assert result.value[0].metadata["catalog_id"] == "biorep_final_bam"


def test_extract_artifacts_skips_known_directory_aggregates(tmp_path, monkeypatch):
    workspace = _artifact_workspace(tmp_path)
    peak_dir = workspace / "results/S1/04_peaks/S1"
    peak_dir.mkdir(parents=True)
    monkeypatch.setattr(
        "encode_pipeline.manifest.make.build_manifest_rows",
        lambda *_args, **_kwargs: (
            [{"output_type": "macs3_peak", "path": str(peak_dir), "status": "present"}],
            0,
            0,
        ),
    )

    result = _adapter().extract_artifacts(WorkflowInputs(config={}), workspace)

    assert result.is_success
    assert result.value == ()


def test_extract_artifacts_fails_closed_for_unknown_vocabulary(tmp_path, monkeypatch):
    workspace = _artifact_workspace(tmp_path)
    monkeypatch.setattr(
        "encode_pipeline.manifest.make.build_manifest_rows",
        lambda *_args, **_kwargs: (
            [{"output_type": "unknown_output", "path": "", "status": "missing"}],
            1,
            0,
        ),
    )

    result = _adapter().extract_artifacts(WorkflowInputs(config={}), workspace)

    assert result.is_failure
    assert result.issues[0].code == "ENCODE_ARTIFACT_CATALOG_MISMATCH"


def test_extract_artifacts_adds_existing_result_manifest(tmp_path, monkeypatch):
    workspace = _artifact_workspace(tmp_path)
    manifest = workspace / "results/multiqc/result_manifest.tsv"
    manifest.parent.mkdir(parents=True)
    manifest.write_text("output_type\tpath\n", encoding="utf-8")
    monkeypatch.setattr(
        "encode_pipeline.manifest.make.build_manifest_rows",
        lambda *_args, **_kwargs: ([], 0, 0),
    )

    result = _adapter().extract_artifacts(WorkflowInputs(config={}), workspace)

    assert result.is_success
    assert result.value[0].output_type == "result_manifest"
    assert result.value[0].relative_path == "results/multiqc/result_manifest.tsv"
