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
    DagPreview,
    JSON_SCHEMA_DIALECT,
    WorkflowAdapter,
    WorkflowInputs,
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
        "command",
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

    assert as_dict["schema_version"] == "3.0.0"
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
    assert "genome_resources" not in config_properties
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
    ]
    assert sample_items["additionalProperties"] is False
    assert len(sample_items["properties"]) == 15
    assert as_dict["option_schema"]["properties"]["strict_inputs"] == {
        "type": "boolean",
        "default": False,
        "description": "Validate FASTQ and Bowtie2 index file existence.",
    }
    assert set(as_dict["option_schema"]["properties"]) == {"strict_inputs", "cores"}
    assert as_dict["option_schema"]["properties"]["cores"] == {
        "type": "integer",
        "minimum": 1,
        "maximum": 1024,
        "default": 1,
        "description": "Maximum CPU cores available to this workflow run.",
    }
    assert as_dict["option_schema"]["additionalProperties"] is False

    for name in ("config_schema", "sample_schema", "option_schema"):
        document = as_dict[name]
        assert document["$schema"] == JSON_SCHEMA_DIALECT
        assert document["$id"].endswith("/3.0.0")
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
        "replicate_analysis": replicate_enabled,
        "chipseq_idr": chipseq_idr_enabled,
    }
    assert result.issues == ()
    assert source == {
        "threads": 2,
        "replicate_analysis": {"enabled": replicate_enabled},
        "chipseq_idr": {"enabled": chipseq_idr_enabled},
    }


@pytest.mark.parametrize("field", ["replicate_analysis", "chipseq_idr"])
@pytest.mark.parametrize("value", [True, False, "true", "false", "TRUE", "FALSE"])
def test_scientific_scalar_switches_preserved_for_validator(field, value):
    source = {field: value}
    result = encode_adapter_module._translate_authoring_config(source)

    assert result.is_success
    expected = {"replicate_analysis": True, "chipseq_idr": False, field: value}
    assert result.value == expected
    assert result.issues == ()
    assert source == {field: value}


@pytest.mark.parametrize(
    ("semantic_key", "value"),
    [
        ("replicate_analysis", None),
        ("replicate_analysis", []),
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


@pytest.mark.parametrize("value", ["yes", " true ", "invalid"])
def test_invalid_scientific_scalar_switch_uses_config_failure(tmp_path, value):
    result = _adapter().validate(
        WorkflowInputs(
            config={"replicate_analysis": value}, samples=[_inline_row(tmp_path)]
        )
    )
    assert result.is_failure
    assert result.errors[0].code == "ENCODE_CONFIG_INVALID"


def test_semantic_switches_validate_without_compatibility_warnings(tmp_path):
    result = _adapter().validate(
        WorkflowInputs(
            config={
                "replicate_analysis": {"enabled": False},
                "chipseq_idr": {"enabled": False},
                "blacklist_mode": "off",
            },
            samples=[_inline_row(tmp_path)],
        )
    )

    assert result.is_success
    assert result.value["config"]["replicate_analysis"] is False
    assert result.value["config"]["chipseq_idr"] is False
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
                    "blacklist_mode": "off",
                },
                samples=rows,
            )
        )

    assert result.is_success
    assert result.value["config"]["replicate_analysis"] is True
    assert result.value["config"]["chipseq_idr"] is True
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


@pytest.mark.parametrize(
    "options",
    [
        {"strict_inputs": False},
        {"strict_inputs": False, "cores": 1},
        {"strict_inputs": False, "cores": 8},
        {"strict_inputs": False, "cores": 1024},
    ],
)
def test_schema_accepts_tiny_profile_authoring_config_and_inline_rows(
    tmp_path, options
):
    schema = _adapter().schema()
    tiny_config = yaml.safe_load(
        (
            SRC_ROOT.parents[0] / "test/profiles/platform_worker_tiny/config.yaml"
        ).read_text(encoding="utf-8")
    )
    tiny_config.pop("samples")
    for field in ("replicate_analysis", "chipseq_idr"):
        tiny_config[field] = {"enabled": tiny_config[field]}
    row = _inline_row(tmp_path)
    public_row = {
        key: value
        for key, value in row.items()
        if key not in {"genome", "bowtie2_index"}
    }

    jsonschema.Draft202012Validator(schema.config_schema).validate(tiny_config)
    jsonschema.Draft202012Validator(schema.sample_schema).validate([public_row])
    jsonschema.Draft202012Validator(schema.option_schema).validate(options)
    result = _adapter().validate(
        WorkflowInputs(config=tiny_config, samples=[row], options=options)
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
        config={
            "samples": str(samples_path),
            "use_control": False,
            "blacklist_mode": "off",
        },
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
        config={
            "samples": str(tmp_path / "missing.tsv"),
            "use_control": False,
            "blacklist_mode": "off",
        },
        samples=samples_path,
    )

    result = _adapter().validate(inputs)

    assert result.is_success
    assert result.value["config"]["samples"] == str(samples_path)


def test_adapter_validate_output_consumable_by_platform(tmp_path):
    samples_path = _write_samples(tmp_path / "samples.tsv")
    result = _adapter().validate(
        WorkflowInputs(
            config={
                "samples": str(samples_path),
                "blacklist_mode": "off",
            }
        ),
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


@pytest.mark.parametrize("cores", [None, True, False, "8", 3.14, 0, -1, 1025, [], {}])
def test_invalid_cores_rejected_by_schema_and_adapter(tmp_path, cores):
    adapter = _adapter()
    options = {"cores": cores}
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(adapter.schema().option_schema).validate(
            options
        )

    result = adapter.validate(
        WorkflowInputs(config={}, samples=[_inline_row(tmp_path)], options=options)
    )

    issue = _first_issue(result)
    assert result.is_failure
    assert issue.code == "ENCODE_OPTIONS_INVALID"
    assert issue.source == "adapter"
    assert issue.path == "options.cores"
    assert issue.message == "cores must be an integer between 1 and 1024"


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
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.Draft202012Validator(_adapter().schema().option_schema).validate(
            {"profile": "local"}
        )


def test_inline_sample_rows_validate_without_config_samples(tmp_path):
    result = _adapter().validate(
        WorkflowInputs(
            config={"use_control": False, "blacklist_mode": "off"},
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
            config={
                "samples": str(missing_server_path),
                "use_control": False,
                "blacklist_mode": "off",
            },
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

    result = _adapter().validate(
        WorkflowInputs(config={"blacklist_mode": "off"}, samples=[row])
    )

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
            WorkflowInputs(
                config={"blacklist_mode": "off"},
                samples=[_inline_row(tmp_path, sample)],
            )
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

    # plan_workspace is implemented; with invalid inputs it fails validation
    # without touching the filesystem.
    workspace_result = adapter.plan_workspace(inputs, workspace)

    assert not workspace.exists()
    for result, expected_method, expected_type in [
        (dag_result, "preview_dag", DagPreview),
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


def test_adapter_tss_fragment_uses_workflow_eligibility(tmp_path):
    gtf = str(Path(__file__).resolve())
    config = {
        "use_control": False,
        "multiqc": False,
        "replicate_analysis": False,
        "analysis_profile": "updated",
        "fragment_cpm": "auto",
        "tss_fragment_profile": True,
        "blacklist_mode": "off",
        "genome_resources": {"hs": {"effective_genome_size": "hs", "gtf": gtf}},
    }
    row = _inline_row(tmp_path)
    row.update({"assay": "cuttag", "target": "UNKNOWN_TARGET"})

    result = _adapter().validate(WorkflowInputs(config=config, samples=[row]))

    assert result.is_failure
    assert any(
        issue.code == "ENCODE_REFERENCE_ADMISSION_INVALID" for issue in result.issues
    )

    row["target"] = "H3K4me3"
    result = _adapter().validate(WorkflowInputs(config=config, samples=[row]))
    assert result.is_success, result.issues


def test_plan_workspace_analysis_mode_survives_snakemake_dryrun(tmp_path):
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    from _tool_resolver import resolve_tool

    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
    config = {
        "use_control": False,
        "multiqc": False,
        "replicate_analysis": False,
        "analysis_profile": "updated",
        "blacklist_mode": "analysis",
        "genome_resources": {
            "hs": {
                "effective_genome_size": "hs",
                "blacklist": str(blacklist),
            }
        },
        "qc": {
            "blacklist_filter": False,
            "frip": False,
            "library_complexity": False,
            "nrf_pbc": False,
            "signal_tracks": False,
            "summary": False,
            "cuttag_fragment_size": False,
        },
    }
    workspace = tmp_path / "workspace"
    plan = _adapter().plan_workspace(
        WorkflowInputs(config=config, samples=[_inline_row(tmp_path)]),
        workspace,
    )
    assert plan.is_success, plan.issues
    for directory in plan.value.directories:
        (workspace / directory).mkdir(parents=True)
    for relative, content in plan.value.files:
        target = workspace / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(content)

    snakemake = resolve_tool("snakemake", "SNAKEMAKE")
    assert snakemake
    result = subprocess.run(
        [
            snakemake,
            "-s",
            str(SRC_ROOT.parent / "workflow" / "Snakefile"),
            "--configfile",
            str(workspace / "config/config.yaml"),
            "--workflow-profile",
            "none",
            "-n",
            "results/.helixweave/analysis_identity.json",
        ],
        cwd=workspace,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr + result.stdout


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


def test_extract_artifacts_extracts_pooled_filtered_peak(tmp_path, monkeypatch):
    workspace = _artifact_workspace(tmp_path)
    (workspace / "config/config.yaml").write_text(
        "samples: config/samples.tsv\nblacklist_mode: analysis\n", encoding="utf-8"
    )
    peak_file = (
        workspace
        / "results/analysis_filtered/experiments/E1/04_peaks/pooled/E1_pooled_peaks_blacklist_filtered/E1_pooled_peaks.blacklist_filtered.narrowPeak"
    )
    peak_file.parent.mkdir(parents=True)
    peak_file.write_text("chr1\t10\t20\n", encoding="utf-8")
    rows = [
        {
            "output_type": "pooled_filtered_peak",
            "path": str(peak_file),
            "status": "present",
            "experiment_id": "E1",
            "assay": "chipseq",
            "genome": "hs",
        },
        {
            "output_type": "pooled_macs3_peak",
            "path": str(peak_file.parents[1] / "E1_pooled_peaks"),
            "status": "present",
            "experiment_id": "E1",
        },
    ]
    monkeypatch.setattr(
        "encode_pipeline.manifest.make.build_manifest_rows",
        lambda *_args, **_kwargs: (rows, 0, 0),
    )

    result = _adapter().extract_artifacts(
        WorkflowInputs(
            config={
                "outdir": "results",
                "analysis_outdir": "results/analysis_filtered",
            }
        ),
        workspace,
    )

    assert result.is_success
    assert len(result.value) == 1
    candidate = result.value[0]
    assert candidate.output_type == "pooled_filtered_peak"
    assert candidate.relative_path == (
        "results/analysis_filtered/experiments/E1/04_peaks/pooled/E1_pooled_peaks_blacklist_filtered/E1_pooled_peaks.blacklist_filtered.narrowPeak"
    )
    assert candidate.mime_type == "text/plain"
    assert candidate.metadata["catalog_id"] == "pooled_filtered_peak"
    assert candidate.metadata["analysis_branch"] == "analysis_filtered"


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


def test_extract_artifacts_uses_analysis_filtered_root(tmp_path):
    workspace = _artifact_workspace(tmp_path)
    blacklist = workspace / "config/blacklist.bed"
    blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
    (workspace / "config/config.yaml").write_text(
        "\n".join(
            [
                "samples: config/samples.tsv",
                "outdir: results",
                "use_control: false",
                "multiqc: false",
                "replicate_analysis: false",
                "chipseq_idr: false",
                "analysis_profile: updated",
                "blacklist_mode: analysis",
                "genome_resources:",
                "  hs:",
                '    effective_genome_size: "hs"',
                f'    blacklist: "{blacklist}"',
                "qc:",
                "  signal_tracks: false",
                "  summary: true",
                "  frip: false",
                "",
            ]
        ),
        encoding="utf-8",
    )
    _write_samples(
        workspace / "config/samples.tsv", _samples_text(workspace / "config")
    )

    raw_bam = workspace / "results/S1/02_align/S1.final.bam"
    raw_bam.parent.mkdir(parents=True)
    raw_bam.write_bytes(b"raw")
    filtered_root = workspace / "results/analysis_filtered"
    filtered_bw = filtered_root / "S1/03_bigwig/S1.CPM.bw"
    filtered_qc = filtered_root / "S1/01_qc/S1.qc_summary.tsv"
    filtered_manifest = filtered_root / "multiqc/result_manifest.tsv"
    for path in (filtered_bw, filtered_qc, filtered_manifest):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x\n", encoding="utf-8")

    result = _adapter().extract_artifacts(WorkflowInputs(config={}), workspace)

    assert result.is_success, result.issues
    paths = {item.relative_path for item in result.value}
    assert "results/S1/02_align/S1.final.bam" in paths
    assert "results/analysis_filtered/S1/03_bigwig/S1.CPM.bw" in paths
    assert "results/analysis_filtered/S1/01_qc/S1.qc_summary.tsv" in paths
    assert "results/analysis_filtered/multiqc/result_manifest.tsv" in paths

    by_path = {item.relative_path: item.metadata for item in result.value}
    filtered_meta = by_path["results/analysis_filtered/S1/03_bigwig/S1.CPM.bw"]
    assert filtered_meta["analysis_profile"] == "updated"
    assert filtered_meta["analysis_branch"] == "analysis_filtered"
    assert filtered_meta["blacklist_mode"] == "analysis"
    assert filtered_meta["coverage_unit"] == "read"
    assert filtered_meta["normalization"] == "CPM"
    assert filtered_meta["reference_assembly"] == "hs"
    raw_meta = by_path["results/S1/02_align/S1.final.bam"]
    assert raw_meta["analysis_branch"] == "raw"


def test_extract_artifacts_adds_existing_result_manifest(tmp_path, monkeypatch):
    workspace = _artifact_workspace(tmp_path)
    (workspace / "config/config.yaml").write_text(
        "samples: config/samples.tsv\nblacklist_mode: off\n", encoding="utf-8"
    )
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


def _write_workspace_config(workspace, lines):
    (workspace / "config/config.yaml").write_text(
        "\n".join(lines) + "\n", encoding="utf-8"
    )


def _metadata_by_path(workspace):
    result = _adapter().extract_artifacts(WorkflowInputs(config={}), workspace)
    assert result.is_success, result.issues
    return {item.relative_path: item.metadata for item in result.value}


@pytest.mark.parametrize("normalization", ["CPM", "RPKM", "BPM", "None"])
def test_artifact_metadata_uses_effective_bamcoverage_normalization(
    tmp_path, normalization
):
    workspace = _artifact_workspace(tmp_path)
    _write_samples(
        workspace / "config/samples.tsv", _samples_text(workspace / "config")
    )
    _write_workspace_config(
        workspace,
        [
            "samples: config/samples.tsv",
            "outdir: results",
            "blacklist_mode: off",
            "use_control: false",
            "multiqc: false",
            "replicate_analysis: false",
            "chipseq_idr: false",
            "analysis_profile: updated",
            "tool_parameters:",
            "  bamcoverage:",
            f'    normalize_using: "{normalization}"',
            "qc:",
            "  signal_tracks: false",
            "  summary: false",
            "  frip: false",
        ],
    )
    cpm = workspace / "results/S1/03_bigwig/S1.CPM.bw"
    cpm.parent.mkdir(parents=True)
    cpm.write_text("x\n", encoding="utf-8")

    metadata = _metadata_by_path(workspace)["results/S1/03_bigwig/S1.CPM.bw"]

    assert metadata["normalization"] == normalization


def test_artifact_metadata_does_not_label_fe_ppois_as_read_coverage(tmp_path):
    workspace = _artifact_workspace(tmp_path)
    _write_samples(
        workspace / "config/samples.tsv", _samples_text(workspace / "config")
    )
    chrom_sizes = workspace / "config/chrom.sizes"
    chrom_sizes.write_text("chr1\t1000\n", encoding="utf-8")
    _write_workspace_config(
        workspace,
        [
            "samples: config/samples.tsv",
            "outdir: results",
            "blacklist_mode: off",
            "use_control: false",
            "multiqc: false",
            "replicate_analysis: false",
            "chipseq_idr: false",
            "analysis_profile: updated",
            "genome_resources:",
            "  hs:",
            '    effective_genome_size: "hs"',
            f'    chrom_sizes: "{chrom_sizes}"',
            "qc:",
            "  signal_tracks: true",
            "  summary: false",
            "  frip: false",
        ],
    )
    for name in ("FE", "ppois"):
        path = workspace / f"results/S1/03_signal/S1.{name}.bw"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x\n", encoding="utf-8")

    metadata = _metadata_by_path(workspace)
    for name in ("FE", "ppois"):
        entry = metadata[f"results/S1/03_signal/S1.{name}.bw"]
        assert entry.get("coverage_unit") != "read"


def test_adapter_surfaces_use_control_false_warning(tmp_path):
    row = _inline_row(tmp_path)
    row["control_sample"] = "C1"
    result = _adapter().validate(
        WorkflowInputs(
            config={
                "use_control": False,
                "multiqc": False,
                "blacklist_mode": "off",
            },
            samples=[row],
        )
    )

    assert result.is_success, result.issues
    assert any(
        issue.code == "ENCODE_SAMPLE_WARNING" and issue.severity == "warning"
        for issue in result.issues
    )


def test_adapter_control_warning_is_request_local_under_concurrency(tmp_path):
    control_row = _inline_row(tmp_path, "S1")
    control_row["control_sample"] = "C1"
    plain_row = _inline_row(tmp_path, "S2")
    config = {
        "use_control": False,
        "multiqc": False,
        "blacklist_mode": "off",
    }

    with ThreadPoolExecutor(max_workers=2) as executor:
        control_result = executor.submit(
            _adapter().validate,
            WorkflowInputs(config=config, samples=[control_row]),
        )
        plain_result = executor.submit(
            _adapter().validate,
            WorkflowInputs(config=config, samples=[plain_row]),
        )
        control = control_result.result()
        plain = plain_result.result()

    assert control.is_success and plain.is_success
    assert any(issue.code == "ENCODE_SAMPLE_WARNING" for issue in control.issues)
    assert not any(issue.code == "ENCODE_SAMPLE_WARNING" for issue in plain.issues)


def test_artifact_branch_distinguishes_sample_named_analysis_filtered(tmp_path):
    workspace = _artifact_workspace(tmp_path)
    samples = _samples_text(workspace / "config").replace("S1", "analysis_filtered")
    _write_samples(workspace / "config/samples.tsv", samples)
    _write_workspace_config(
        workspace,
        [
            "samples: config/samples.tsv",
            "outdir: results",
            "use_control: false",
            "multiqc: false",
            "replicate_analysis: false",
            "chipseq_idr: false",
            "analysis_profile: legacy",
            "blacklist_mode: off",
            "qc:",
            "  signal_tracks: false",
            "  summary: false",
            "  frip: false",
        ],
    )
    raw_bam = (
        workspace / "results/analysis_filtered/02_align/analysis_filtered.final.bam"
    )
    raw_bam.parent.mkdir(parents=True)
    raw_bam.write_text("raw\n", encoding="utf-8")

    metadata = _metadata_by_path(workspace)
    assert (
        metadata["results/analysis_filtered/02_align/analysis_filtered.final.bam"][
            "analysis_branch"
        ]
        == "raw"
    )


def test_artifact_branch_does_not_mislabel_raw_sample_named_analysis_filtered(
    tmp_path,
):
    workspace = _artifact_workspace(tmp_path)
    samples = _samples_text(workspace / "config").replace("S1", "analysis_filtered")
    _write_samples(workspace / "config/samples.tsv", samples)
    blacklist = workspace / "config/blacklist.bed"
    blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
    _write_workspace_config(
        workspace,
        [
            "samples: config/samples.tsv",
            "outdir: results",
            "use_control: false",
            "multiqc: false",
            "replicate_analysis: false",
            "chipseq_idr: false",
            "analysis_profile: updated",
            "blacklist_mode: analysis",
            "genome_resources:",
            "  hs:",
            '    effective_genome_size: "hs"',
            f'    blacklist: "{blacklist}"',
            "qc:",
            "  signal_tracks: false",
            "  summary: false",
            "  frip: false",
        ],
    )
    raw_bam = (
        workspace / "results/analysis_filtered/02_align/analysis_filtered.final.bam"
    )
    raw_bam.parent.mkdir(parents=True)
    raw_bam.write_text("raw\n", encoding="utf-8")
    filtered_cpm = (
        workspace / "results/analysis_filtered/analysis_filtered/03_bigwig/"
        "analysis_filtered.CPM.bw"
    )
    filtered_cpm.parent.mkdir(parents=True)
    filtered_cpm.write_text("x\n", encoding="utf-8")

    metadata = _metadata_by_path(workspace)
    assert (
        metadata["results/analysis_filtered/02_align/analysis_filtered.final.bam"][
            "analysis_branch"
        ]
        == "raw"
    )
    assert (
        metadata[
            "results/analysis_filtered/analysis_filtered/03_bigwig/"
            "analysis_filtered.CPM.bw"
        ]["analysis_branch"]
        == "analysis_filtered"
    )


def _write_rows(path, rows):
    output = StringIO()
    writer = csv.DictWriter(
        output, fieldnames=tuple(rows[0]), delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        writer.writerow(row)
    path.write_text(output.getvalue(), encoding="utf-8")


def test_artifact_metadata_mnase_dyad_is_not_read_coverage(tmp_path):
    workspace = _artifact_workspace(tmp_path)
    row = _inline_row(tmp_path, "M1")
    row.update({"assay": "mnase", "target": "H3", "peak_mode": "nucleosome"})
    _write_rows(workspace / "config/samples.tsv", [row])
    chrom_sizes = workspace / "config/chrom.sizes"
    chrom_sizes.write_text("chr1\t1000\n", encoding="utf-8")
    _write_workspace_config(
        workspace,
        [
            "samples: config/samples.tsv",
            "outdir: results",
            "blacklist_mode: off",
            "use_control: false",
            "multiqc: false",
            "replicate_analysis: false",
            "chipseq_idr: false",
            "analysis_profile: updated",
            "tool_parameters:",
            "  bamcoverage:",
            '    normalize_using: "RPKM"',
            "genome_resources:",
            "  hs:",
            '    effective_genome_size: "hs"',
            f'    chrom_sizes: "{chrom_sizes}"',
            "qc:",
            "  signal_tracks: true",
            "  summary: false",
            "  frip: false",
        ],
    )
    for name in ("dyad", "mono"):
        path = workspace / f"results/M1/04_signal/M1.{name}.CPM.bw"
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x\n", encoding="utf-8")

    metadata = _metadata_by_path(workspace)
    dyad = metadata["results/M1/04_signal/M1.dyad.CPM.bw"]
    mono = metadata["results/M1/04_signal/M1.mono.CPM.bw"]
    assert "coverage_unit" not in dyad
    assert dyad["normalization"] == "RPKM"
    assert mono["coverage_unit"] == "read"
    assert mono["normalization"] == "RPKM"


def test_artifact_metadata_pooled_mnase_dyad_is_not_read_coverage(tmp_path):
    workspace = _artifact_workspace(tmp_path)
    rows = []
    for sample_id, biorep in (("M1", "1"), ("M2", "2")):
        row = _inline_row(tmp_path, sample_id)
        row.update(
            {
                "assay": "mnase",
                "target": "H3",
                "peak_mode": "nucleosome",
                "experiment": "EXP1",
                "biological_replicate": biorep,
            }
        )
        rows.append(row)
    _write_rows(workspace / "config/samples.tsv", rows)
    chrom_sizes = workspace / "config/chrom.sizes"
    chrom_sizes.write_text("chr1\t1000\n", encoding="utf-8")
    _write_workspace_config(
        workspace,
        [
            "samples: config/samples.tsv",
            "outdir: results",
            "blacklist_mode: off",
            "use_control: false",
            "multiqc: false",
            "replicate_analysis: true",
            "chipseq_idr: false",
            "analysis_profile: updated",
            "genome_resources:",
            "  hs:",
            '    effective_genome_size: "hs"',
            f'    chrom_sizes: "{chrom_sizes}"',
            "qc:",
            "  signal_tracks: true",
            "  summary: false",
            "  frip: false",
        ],
    )
    for name in ("dyad", "mono"):
        path = (
            workspace / f"results/experiments/EXP1/04_signal/EXP1.pooled.{name}.CPM.bw"
        )
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("x\n", encoding="utf-8")

    metadata = _metadata_by_path(workspace)
    dyad = metadata["results/experiments/EXP1/04_signal/EXP1.pooled.dyad.CPM.bw"]
    mono = metadata["results/experiments/EXP1/04_signal/EXP1.pooled.mono.CPM.bw"]
    assert "coverage_unit" not in dyad
    assert mono["coverage_unit"] == "read"


def test_adapter_validates_and_plans_bowtie2_max_insert_and_extra_args_migration(
    tmp_path,
):
    row = _inline_row(tmp_path, "M1")
    row.update({"assay": "mnase", "target": "H3", "peak_mode": "nucleosome"})
    config = {
        "use_control": False,
        "multiqc": False,
        "replicate_analysis": False,
        "chipseq_idr": False,
        "blacklist_mode": "off",
        "tool_parameters": {
            "bowtie2": {
                "extra_args": "-X 900 --very-sensitive",
            }
        },
    }
    first = _adapter().plan_workspace(
        WorkflowInputs(config=config, samples=[row]),
        tmp_path / "workspace-1",
    )
    assert first.is_success, first.issues
    rendered_first = yaml.safe_load(dict(first.value.files)["config/config.yaml"])
    assert rendered_first["tool_parameters"]["bowtie2"]["max_insert"] == 900
    assert (
        rendered_first["tool_parameters"]["bowtie2"]["extra_args"] == "--very-sensitive"
    )

    # Re-planning from the normalized config is byte-identical.
    second_cfg = dict(rendered_first)
    second_cfg.pop("samples", None)
    second_cfg.pop("outdir", None)
    second = _adapter().plan_workspace(
        WorkflowInputs(config=second_cfg, samples=[row]),
        tmp_path / "workspace-2",
    )
    assert second.is_success, second.issues
    assert first.value.files == second.value.files


def test_adapter_rejects_bowtie2_minins_exceeding_resolved_max_insert(tmp_path):
    row = _inline_row(tmp_path, "M1")
    row.update({"assay": "mnase", "target": "H3", "peak_mode": "nucleosome"})
    config = {
        "use_control": False,
        "multiqc": False,
        "replicate_analysis": False,
        "chipseq_idr": False,
        "tool_parameters": {
            "bowtie2": {
                "max_insert": "auto",
                "extra_args": "-I 600",
            }
        },
    }
    validated = _adapter().validate(WorkflowInputs(config=config, samples=[row]))
    assert validated.is_failure
    assert validated.errors[0].code == "ENCODE_CONFIG_INVALID"

    planned = _adapter().plan_workspace(
        WorkflowInputs(config=config, samples=[row]),
        tmp_path / "workspace",
    )
    assert planned.is_failure
    assert planned.errors[0].code == "ENCODE_CONFIG_INVALID"


def test_extract_artifacts_default_analysis_without_canonical_key(tmp_path):
    workspace = _artifact_workspace(tmp_path)
    _write_samples(
        workspace / "config/samples.tsv", _samples_text(workspace / "config")
    )
    _write_workspace_config(
        workspace,
        [
            "samples: config/samples.tsv",
            "outdir: results",
            "use_control: false",
            "multiqc: false",
            "replicate_analysis: false",
            "chipseq_idr: false",
            "analysis_profile: updated",
            "qc:",
            "  signal_tracks: false",
            "  summary: false",
            "  frip: false",
        ],
    )
    cpm = workspace / "results/analysis_filtered/S1/03_bigwig/S1.CPM.bw"
    cpm.parent.mkdir(parents=True)
    cpm.write_text("x\n", encoding="utf-8")

    result = _adapter().extract_artifacts(WorkflowInputs(config={}), workspace)

    assert result.is_success, result.issues
    by_type = {item.output_type: item for item in result.value}
    assert (
        by_type["cpm_bigwig"].relative_path
        == "results/analysis_filtered/S1/03_bigwig/S1.CPM.bw"
    )
