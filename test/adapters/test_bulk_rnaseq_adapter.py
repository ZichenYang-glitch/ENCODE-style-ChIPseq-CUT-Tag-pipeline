"""Contract and validation tests for the bulk RNA-seq adapter skeleton."""

from __future__ import annotations

from copy import deepcopy
import hashlib
import json
from pathlib import Path

import jsonschema
import pytest

from encode_pipeline import __version__
from encode_pipeline.adapters.bulk_rnaseq import BulkRnaSeqWorkflowAdapter
from encode_pipeline.adapters.bulk_rnaseq.adapter import WORKFLOW_ID
from encode_pipeline.adapters.bulk_rnaseq.authoring import SCHEMA_VERSION
from encode_pipeline.adapters.bulk_rnaseq.upstream import (
    ADVANCED_NATIVE_PARAMETERS,
    FIXED_PRODUCT_NATIVE_PARAMETERS,
    KNOWN_UPSTREAM_PARAMETERS,
    NFCORE_RNASEQ_COMMIT,
    NFCORE_RNASEQ_LICENSE,
    NFCORE_RNASEQ_RELEASE,
    PARAMETER_IMPACT_BY_NATIVE_NAME,
    PLATFORM_OWNED_NATIVE_PARAMETERS,
    RAW_ARGUMENT_NATIVE_PARAMETERS,
    STANDARD_NATIVE_PARAMETERS,
    UNSUPPORTED_NATIVE_PARAMETERS,
    UPSTREAM_LICENSE_FILE,
    UPSTREAM_LICENSE_SHA256,
    UPSTREAM_LICENSE_SIZE,
    UPSTREAM_PARAMETER_SCHEMA_FILE,
    UPSTREAM_PARAMETER_SCHEMA_SHA256,
    UPSTREAM_PARAMETER_SCHEMA_SIZE,
    UPSTREAM_SAMPLESHEET_SCHEMA_FILE,
    UPSTREAM_SAMPLESHEET_SCHEMA_SHA256,
    UPSTREAM_SAMPLESHEET_SCHEMA_SIZE,
    load_upstream_parameter_schema,
    load_upstream_samplesheet_schema,
    upstream_parameter_properties,
)
from encode_pipeline.adapters.encode import EncodeStyleWorkflowAdapter
from encode_pipeline.platform.adapters import (
    InputBundleImportingAdapter,
    QcSummaryExtractingAdapter,
    WorkflowAdapter,
    WorkflowInputs,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.testing.adapter_conformance import (
    AdapterConformanceCase,
    verify_adapter_conformance,
)


PROJECT_ROOT = Path(__file__).resolve().parents[2]
CONTRACT_ROOT = PROJECT_ROOT / "src/encode_pipeline/contracts/nfcore_rnaseq"


def _reference(*, indexes: bool = False) -> dict[str, object]:
    value: dict[str, object] = {
        "reference_id": "GRCh38-test",
        "fasta": "/refs/GRCh38.fa",
        "fasta_sha256": "a" * 64,
        "gtf": "/refs/gencode.gtf",
        "gtf_sha256": "b" * 64,
        "annotation_style": "gencode",
    }
    if indexes:
        value.update(
            {
                "star_index": {
                    "path": "/refs/star-index",
                    "identity_sha256": "c" * 64,
                },
                "salmon_index": {
                    "path": "/refs/salmon-index",
                    "identity_sha256": "d" * 64,
                },
            }
        )
    return value


def _ribosomal_rna_removal(
    *,
    tool: str = "sortmerna",
    save_filtered_reads: bool = False,
    sortmerna_index: bool = False,
) -> dict[str, object]:
    value: dict[str, object] = {
        "enabled": True,
        "tool": tool,
        "save_filtered_reads": save_filtered_reads,
        "database_manifest": {
            "path": "/refs/rrna/database-manifest.txt",
            "identity_sha256": "e" * 64,
        },
    }
    if sortmerna_index:
        value["sortmerna_index"] = {
            "path": "/refs/rrna/sortmerna-index",
            "identity_sha256": "f" * 64,
        }
    return value


def _config(**standard_updates: object) -> dict[str, object]:
    standard = {"reference": _reference(), **standard_updates}
    return {"standard": standard}


def _sample(
    *,
    sample: str = "S1",
    library: str = "lib1",
    lane: str = "L001",
    layout: str = "PE",
    strandedness: str = "auto",
) -> dict[str, str]:
    row = {
        "sample": sample,
        "library": library,
        "lane": lane,
        "layout": layout,
        "fastq_1": f"/data/{sample}.{library}.{lane}.R1.fastq.gz",
        "strandedness": strandedness,
        "platform": "ILLUMINA",
    }
    if layout == "PE":
        row["fastq_2"] = f"/data/{sample}.{library}.{lane}.R2.fastq.gz"
    return row


def _inputs(
    *,
    config: dict[str, object] | None = None,
    samples: list[dict[str, str]] | None = None,
    options: dict[str, object] | None = None,
) -> WorkflowInputs:
    return WorkflowInputs(
        config=_config() if config is None else config,
        samples=[_sample()] if samples is None else samples,
        options={} if options is None else options,
    )


def _error_code(inputs: WorkflowInputs) -> str:
    result = BulkRnaSeqWorkflowAdapter().validate(inputs)
    assert result.is_failure
    return result.errors[0].code


def test_adapter_identity_and_capabilities_are_truthful():
    adapter = BulkRnaSeqWorkflowAdapter()

    assert isinstance(adapter, WorkflowAdapter)
    assert not isinstance(adapter, QcSummaryExtractingAdapter)
    assert not isinstance(adapter, InputBundleImportingAdapter)
    assert adapter.metadata.workflow_id == "bulk-rnaseq"
    assert adapter.metadata.engines == ("nextflow",)
    assert adapter.metadata.version == __version__
    assert adapter.capabilities.supports == ("validation", "input_authoring")


def test_adapter_passes_reusable_conformance_suite(tmp_path: Path):
    valid = _inputs()
    invalid_config = _config()
    invalid_config["standard"]["reference"]["fasta_sha256"] = "not-a-digest"

    verify_adapter_conformance(
        AdapterConformanceCase(
            adapter=BulkRnaSeqWorkflowAdapter(),
            valid_inputs=valid,
            invalid_inputs=_inputs(config=invalid_config),
            planning_workspace=(tmp_path / "planning").resolve(),
            artifact_workspace=(tmp_path / "artifacts").resolve(),
        )
    )


def test_schema_is_complete_versioned_stable_and_fresh():
    adapter = BulkRnaSeqWorkflowAdapter()
    first = adapter.schema()
    second = adapter.schema()
    first_dict = first.to_dict()
    second_dict = second.to_dict()

    assert first is not second
    assert first.schema_version == SCHEMA_VERSION == "1.0.0"
    assert first.coverage.to_dict() == {
        "config": "complete",
        "samples": "complete",
        "options": "complete",
    }
    assert first.authoring_modes.to_dict() == {
        "config": ["schema_form", "yaml"],
        "samples": ["tsv_upload", "inline_table"],
        "options": ["schema_form"],
    }
    assert first.input_modes.to_dict() == {
        "config": ["object"],
        "samples": ["inline_rows"],
        "options": ["object"],
    }
    assert json.dumps(first_dict, sort_keys=True) == json.dumps(
        second_dict, sort_keys=True
    )
    for name in ("config_schema", "sample_schema", "option_schema"):
        jsonschema.Draft202012Validator.check_schema(first_dict[name])
        assert first_dict[name]["$id"].endswith("/1.0.0")

    first.config_schema["properties"].clear()
    assert "standard" in adapter.schema().config_schema["properties"]


def test_public_schema_closes_all_surfaces_and_does_not_seed_advanced_defaults():
    document = BulkRnaSeqWorkflowAdapter().schema().to_dict()
    config = document["config_schema"]
    samples = document["sample_schema"]
    options = document["option_schema"]
    advanced = config["properties"]["advanced"]
    ribosomal_rna_removal = config["properties"]["standard"]["properties"][
        "ribosomal_rna_removal"
    ]

    assert config["additionalProperties"] is False
    assert config["properties"]["standard"]["additionalProperties"] is False
    assert advanced["additionalProperties"] is False
    assert set(advanced["properties"]) == set(ADVANCED_NATIVE_PARAMETERS)
    assert all("default" not in value for value in advanced["properties"].values())
    assert ribosomal_rna_removal["additionalProperties"] is False
    assert ribosomal_rna_removal["default"] == {"enabled": False}
    assert ribosomal_rna_removal["properties"]["tool"]["enum"] == [
        "sortmerna",
        "bowtie2",
    ]
    assert "default" not in ribosomal_rna_removal["properties"]["tool"]
    assert (
        ribosomal_rna_removal["properties"]["database_manifest"]["additionalProperties"]
        is False
    )
    assert (
        ribosomal_rna_removal["properties"]["sortmerna_index"]["additionalProperties"]
        is False
    )
    assert samples["items"]["additionalProperties"] is False
    assert all(
        property_schema.get("type") == "string"
        for property_schema in samples["items"]["properties"].values()
    )
    assert options == {
        "$schema": "https://json-schema.org/draft/2020-12/schema",
        "$id": "https://helixweave.org/schemas/bulk-rnaseq/options/1.0.0",
        "title": "HelixWeave bulk RNA-seq adapter options",
        "description": "No caller-owned platform options are defined in schema 1.0.0.",
        "type": "object",
        "properties": {},
        "additionalProperties": False,
    }


def test_pinned_upstream_contract_identity_matches_original_bytes():
    expected = {
        UPSTREAM_PARAMETER_SCHEMA_FILE: (
            UPSTREAM_PARAMETER_SCHEMA_SIZE,
            UPSTREAM_PARAMETER_SCHEMA_SHA256,
        ),
        UPSTREAM_SAMPLESHEET_SCHEMA_FILE: (
            UPSTREAM_SAMPLESHEET_SCHEMA_SIZE,
            UPSTREAM_SAMPLESHEET_SCHEMA_SHA256,
        ),
        UPSTREAM_LICENSE_FILE: (UPSTREAM_LICENSE_SIZE, UPSTREAM_LICENSE_SHA256),
    }
    for filename, (size, digest) in expected.items():
        content = (CONTRACT_ROOT / filename).read_bytes()
        assert len(content) == size
        assert hashlib.sha256(content).hexdigest() == digest

    provenance = json.loads((CONTRACT_ROOT / "provenance.json").read_text())
    assert provenance["release"] == NFCORE_RNASEQ_RELEASE == "3.26.0"
    assert provenance["release_name"].endswith("Chromium Cuttlefish")
    assert provenance["draft"] is False
    assert provenance["prerelease"] is False
    assert provenance["commit"] == NFCORE_RNASEQ_COMMIT
    assert provenance["license"] == NFCORE_RNASEQ_LICENSE == "MIT"
    assert {entry["sha256"] for entry in provenance["files"]} == {
        digest for _, digest in expected.values()
    }


def test_upstream_schema_is_exact_and_parameter_policy_is_closed():
    document = load_upstream_parameter_schema()
    properties = upstream_parameter_properties()
    policy_classes = (
        STANDARD_NATIVE_PARAMETERS,
        ADVANCED_NATIVE_PARAMETERS,
        PLATFORM_OWNED_NATIVE_PARAMETERS,
        RAW_ARGUMENT_NATIVE_PARAMETERS,
        UNSUPPORTED_NATIVE_PARAMETERS,
    )

    assert document["$schema"] == "https://json-schema.org/draft/2020-12/schema"
    assert load_upstream_samplesheet_schema()["$schema"] == (
        "https://json-schema.org/draft/2020-12/schema"
    )
    assert len(properties) == len(KNOWN_UPSTREAM_PARAMETERS) == 133
    assert set().union(*policy_classes) == set(KNOWN_UPSTREAM_PARAMETERS)
    for index, group in enumerate(policy_classes):
        assert not any(
            group.intersection(other) for other in policy_classes[index + 1 :]
        )
    assert set(PARAMETER_IMPACT_BY_NATIVE_NAME) == set(
        STANDARD_NATIVE_PARAMETERS
        | ADVANCED_NATIVE_PARAMETERS
        | FIXED_PRODUCT_NATIVE_PARAMETERS
    )


@pytest.mark.parametrize("layout", ["SE", "PE"])
def test_valid_se_and_pe_inputs_normalize_to_star_salmon(layout: str):
    result = BulkRnaSeqWorkflowAdapter().validate(
        _inputs(samples=[_sample(layout=layout)])
    )

    assert result.is_success
    assert result.value["contract"] == {
        "workflow_id": WORKFLOW_ID,
        "schema_version": "1.0.0",
        "nfcore_rnaseq_release": "3.26.0",
        "nfcore_rnaseq_commit": NFCORE_RNASEQ_COMMIT,
        "upstream_parameter_schema_sha256": UPSTREAM_PARAMETER_SCHEMA_SHA256,
        "upstream_samplesheet_schema_sha256": UPSTREAM_SAMPLESHEET_SCHEMA_SHA256,
    }
    assert result.value["nfcore_params"]["aligner"] == "star_salmon"
    assert result.value["nfcore_params"]["igenomes_ignore"] is True
    assert result.value["nfcore_params"]["remove_ribo_rna"] is False
    assert result.value["nfcore_params"]["save_non_ribo_reads"] is False
    assert "ribo_removal_tool" not in result.value["nfcore_params"]
    assert "ribo_database_manifest" not in result.value["nfcore_params"]
    assert "sortmerna_index" not in result.value["nfcore_params"]
    assert result.value["nfcore_params"]["skip_alignment"] is False
    assert result.value["nfcore_params"]["skip_quantification_merge"] is False
    assert set(result.value["nfcore_params"]).issubset(PARAMETER_IMPACT_BY_NATIVE_NAME)
    expected_fastq_2 = "" if layout == "SE" else "/data/S1.lib1.L001.R2.fastq.gz"
    assert result.value["samples"][0]["fastq_2"] == expected_fastq_2


def test_repeated_sample_rows_preserve_library_and_lane_semantics():
    rows = [
        _sample(library="lib1", lane="L001"),
        _sample(library="lib1", lane="L002"),
        _sample(library="lib2", lane="L001"),
    ]

    result = BulkRnaSeqWorkflowAdapter().validate(_inputs(samples=rows))

    assert result.is_success
    assert [(row["library"], row["lane"]) for row in result.value["samples"]] == [
        ("lib1", "L001"),
        ("lib1", "L002"),
        ("lib2", "L001"),
    ]


def test_supported_fna_reference_suffix_is_accepted():
    config = _config()
    config["standard"]["reference"]["fasta"] = "/refs/GRCh38.fna.gz"

    assert BulkRnaSeqWorkflowAdapter().validate(_inputs(config=config)).is_success


def test_salmon_index_is_rejected_when_no_sample_needs_auto_strand_inference():
    config = _config(reference=_reference(indexes=True))

    assert (
        _error_code(_inputs(config=config, samples=[_sample(strandedness="reverse")]))
        == "BULK_RNASEQ_REFERENCE_CONTEXT_CONFLICT"
    )


def test_validation_and_serialization_are_deterministic_and_do_not_mutate_inputs():
    source_config = _config(
        trimming={"enabled": True, "tool": "fastp"},
        ribosomal_rna_removal=_ribosomal_rna_removal(
            save_filtered_reads=True,
            sortmerna_index=True,
        ),
        outputs={"bigwig": True, "stringtie": False},
    )
    source_samples = [_sample()]
    original_config = deepcopy(source_config)
    original_samples = deepcopy(source_samples)
    inputs = _inputs(config=source_config, samples=source_samples)
    adapter = BulkRnaSeqWorkflowAdapter()

    first = adapter.validate(inputs)
    second = adapter.validate(inputs)

    assert first.is_success and second.is_success
    assert json.dumps(first.value, sort_keys=True, separators=(",", ":")) == json.dumps(
        second.value, sort_keys=True, separators=(",", ":")
    )
    assert list(first.value["nfcore_params"]) == sorted(first.value["nfcore_params"])
    assert source_config == original_config
    assert source_samples == original_samples


def test_validation_does_not_probe_user_paths(monkeypatch: pytest.MonkeyPatch):
    forbidden = {
        "/refs/GRCh38.fa",
        "/refs/gencode.gtf",
        "/refs/rrna/database-manifest.txt",
        "/refs/rrna/sortmerna-index",
        "/data/S1.lib1.L001.R1.fastq.gz",
        "/data/S1.lib1.L001.R2.fastq.gz",
    }
    originals = {
        name: getattr(Path, name) for name in ("exists", "is_file", "stat", "open")
    }

    for method_name, original in originals.items():

        def guarded(self, *args, _name=method_name, _original=original, **kwargs):
            if str(self) in forbidden:
                raise AssertionError(f"user path probed through Path.{_name}")
            return _original(self, *args, **kwargs)

        monkeypatch.setattr(Path, method_name, guarded)

    config = _config(ribosomal_rna_removal=_ribosomal_rna_removal(sortmerna_index=True))
    assert BulkRnaSeqWorkflowAdapter().validate(_inputs(config=config)).is_success


def test_sortmerna_removal_normalizes_typed_manifest_index_and_output_selection():
    config = _config(
        ribosomal_rna_removal=_ribosomal_rna_removal(
            save_filtered_reads=True,
            sortmerna_index=True,
        )
    )

    result = BulkRnaSeqWorkflowAdapter().validate(_inputs(config=config))

    assert result.is_success
    params = result.value["nfcore_params"]
    assert params["remove_ribo_rna"] is True
    assert params["ribo_removal_tool"] == "sortmerna"
    assert params["save_non_ribo_reads"] is True
    assert params["ribo_database_manifest"] == "/refs/rrna/database-manifest.txt"
    assert params["sortmerna_index"] == "/refs/rrna/sortmerna-index"
    assert (
        result.value["reference_identity"]["ribo_database_manifest_sha256"] == "e" * 64
    )
    assert result.value["reference_identity"]["sortmerna_index_sha256"] == "f" * 64
    impacts = {
        item["native_parameter"]: item["impact"]
        for item in result.value["parameter_impacts"]
    }
    assert impacts["remove_ribo_rna"] == "route_namespace"
    assert impacts["ribo_removal_tool"] == "route_namespace"
    assert impacts["ribo_database_manifest"] == "content_only"
    assert impacts["sortmerna_index"] == "route_namespace"
    assert impacts["save_non_ribo_reads"] == "additive_artifacts"


def test_bowtie2_removal_requires_manifest_and_does_not_emit_sortmerna_index():
    config = _config(ribosomal_rna_removal=_ribosomal_rna_removal(tool="bowtie2"))

    result = BulkRnaSeqWorkflowAdapter().validate(_inputs(config=config))

    assert result.is_success
    params = result.value["nfcore_params"]
    assert params["remove_ribo_rna"] is True
    assert params["ribo_removal_tool"] == "bowtie2"
    assert params["ribo_database_manifest"] == "/refs/rrna/database-manifest.txt"
    assert params["save_non_ribo_reads"] is False
    assert "sortmerna_index" not in params
    assert (
        result.value["reference_identity"]["ribo_database_manifest_sha256"] == "e" * 64
    )
    assert "sortmerna_index_sha256" not in result.value["reference_identity"]


def test_bowtie2_can_save_filtered_reads_for_single_end_samples():
    config = _config(
        ribosomal_rna_removal=_ribosomal_rna_removal(
            tool="bowtie2",
            save_filtered_reads=True,
        )
    )

    result = BulkRnaSeqWorkflowAdapter().validate(
        _inputs(config=config, samples=[_sample(layout="SE")])
    )

    assert result.is_success
    assert result.value["nfcore_params"]["save_non_ribo_reads"] is True


@pytest.mark.parametrize(
    "samples",
    [
        [_sample(layout="PE")],
        [
            _sample(sample="S1", layout="SE"),
            _sample(sample="S2", layout="PE"),
        ],
    ],
)
def test_bowtie2_filtered_read_output_is_rejected_for_pe_or_mixed_layouts(samples):
    config = _config(
        ribosomal_rna_removal=_ribosomal_rna_removal(
            tool="bowtie2",
            save_filtered_reads=True,
        )
    )

    assert _error_code(_inputs(config=config, samples=samples)) == (
        "BULK_RNASEQ_RRNA_CONFLICT"
    )


@pytest.mark.parametrize("save_filtered_reads", [False, True])
def test_filtered_read_output_selection_maps_without_changing_tool_semantics(
    save_filtered_reads: bool,
):
    config = _config(
        ribosomal_rna_removal=_ribosomal_rna_removal(
            save_filtered_reads=save_filtered_reads
        )
    )

    result = BulkRnaSeqWorkflowAdapter().validate(_inputs(config=config))

    assert result.is_success
    assert result.value["nfcore_params"]["save_non_ribo_reads"] is (save_filtered_reads)


@pytest.mark.parametrize(
    "extra",
    [
        {"tool": "sortmerna"},
        {"save_filtered_reads": False},
        {"save_filtered_reads": True},
        {
            "database_manifest": {
                "path": "/refs/rrna/database-manifest.txt",
                "identity_sha256": "e" * 64,
            }
        },
        {
            "sortmerna_index": {
                "path": "/refs/rrna/sortmerna-index",
                "identity_sha256": "f" * 64,
            }
        },
    ],
)
def test_disabled_ribosomal_rna_removal_rejects_all_extra_configuration(extra):
    config = _config(ribosomal_rna_removal={"enabled": False, **extra})

    assert _error_code(_inputs(config=config)) == "BULK_RNASEQ_RRNA_CONFLICT"


@pytest.mark.parametrize("tool", ["sortmerna", "bowtie2"])
def test_enabled_ribosomal_rna_removal_requires_database_manifest(tool: str):
    config = _config(
        ribosomal_rna_removal={
            "enabled": True,
            "tool": tool,
            "save_filtered_reads": False,
        }
    )

    assert _error_code(_inputs(config=config)) == "BULK_RNASEQ_STANDARD_INVALID"


def test_sortmerna_index_is_rejected_for_bowtie2():
    config = _config(
        ribosomal_rna_removal=_ribosomal_rna_removal(
            tool="bowtie2",
            sortmerna_index=True,
        )
    )

    assert _error_code(_inputs(config=config)) == "BULK_RNASEQ_RRNA_CONFLICT"


@pytest.mark.parametrize(
    "ribosomal_rna_removal",
    [
        {"enabled": "true"},
        {"enabled": 1},
        {
            **_ribosomal_rna_removal(),
            "save_filtered_reads": "true",
        },
        {
            **_ribosomal_rna_removal(),
            "tool": True,
        },
        {
            "enabled": True,
            "tool": "sortmerna",
            "database_manifest": {
                "path": "/refs/rrna/database-manifest.txt",
            },
        },
        {
            **_ribosomal_rna_removal(),
            "database_manifest": {
                "path": "https://example.invalid/rrna-manifest.txt",
                "identity_sha256": "e" * 64,
            },
        },
        {
            **_ribosomal_rna_removal(),
            "sortmerna_index": {
                "path": "/refs/rrna/sortmerna-index",
                "identity_sha256": True,
            },
        },
        {
            **_ribosomal_rna_removal(),
            "unexpected": True,
        },
        {
            **_ribosomal_rna_removal(),
            "tool": "ribodetector",
        },
    ],
)
def test_ribosomal_rna_schema_rejects_string_bools_unknown_fields_and_tools(
    ribosomal_rna_removal,
):
    config = _config(ribosomal_rna_removal=ribosomal_rna_removal)

    assert _error_code(_inputs(config=config)) == "BULK_RNASEQ_STANDARD_INVALID"


def test_valid_advanced_allowlist_is_exactly_validated_and_classified():
    reference = _reference(indexes=True)
    reference["annotation_style"] = "ensembl"
    config = _config(
        reference=reference,
        outputs={"stringtie": True},
    )
    config["advanced"] = {
        "bam_csi_index": True,
        "deseq2_vst": False,
        "featurecounts_feature_type": "exon",
        "featurecounts_group_type": "gene_biotype",
        "gffread_transcript_fasta": True,
        "gtf_extra_attributes": "gene_name,transcript_id",
        "gtf_group_features": "gene_id",
        "min_mapped_reads": 10.5,
        "min_trimmed_reads": 1,
        "rseqc_modules": "bam_stat,tin",
        "skip_gtf_filter": False,
        "skip_gtf_transcript_filter": False,
        "star_ignore_sjdbgtf": True,
        "stranded_threshold": 0.7,
        "stringtie_ignore_gtf": True,
        "unstranded_threshold": 0.2,
    }

    result = BulkRnaSeqWorkflowAdapter().validate(_inputs(config=config))

    assert result.is_success
    for name, value in config["advanced"].items():
        assert result.value["nfcore_params"][name] == value
    impacts = {
        item["native_parameter"]: item["impact"]
        for item in result.value["parameter_impacts"]
    }
    assert set(config["advanced"]).issubset(impacts)
    assert impacts["bam_csi_index"] == "filename_or_format"
    assert impacts["min_trimmed_reads"] == "sample_set"
    assert impacts["rseqc_modules"] == "artifact_set"
    assert impacts["stranded_threshold"] == "route_namespace"
    assert impacts["unstranded_threshold"] == "route_namespace"


@pytest.mark.parametrize(
    ("name", "value", "code"),
    [
        ("aligner", "star_salmon", "BULK_RNASEQ_PARAMETER_SURFACE_CONFLICT"),
        ("input", "/tmp/samples.csv", "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN"),
        ("workDir", "/tmp/work", "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN"),
        ("profile", "docker", "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN"),
        ("resume", True, "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN"),
        ("plugins", ["nf-schema"], "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN"),
        ("tower", True, "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN"),
        ("wave", True, "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN"),
        ("raw_config", "process {}", "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN"),
        ("extra_star_align_args", "--runMode x", "BULK_RNASEQ_RAW_ARGUMENT_FORBIDDEN"),
        ("remove_ribo_rna", True, "BULK_RNASEQ_PARAMETER_SURFACE_CONFLICT"),
        (
            "ribo_removal_tool",
            "sortmerna",
            "BULK_RNASEQ_PARAMETER_SURFACE_CONFLICT",
        ),
        (
            "save_non_ribo_reads",
            True,
            "BULK_RNASEQ_PARAMETER_SURFACE_CONFLICT",
        ),
        (
            "ribo_database_manifest",
            "/refs/rrna/database-manifest.txt",
            "BULK_RNASEQ_PARAMETER_SURFACE_CONFLICT",
        ),
        (
            "sortmerna_index",
            "/refs/rrna/sortmerna-index",
            "BULK_RNASEQ_PARAMETER_SURFACE_CONFLICT",
        ),
        (
            "use_gpu_ribodetector",
            True,
            "BULK_RNASEQ_PLATFORM_PARAMETER_FORBIDDEN",
        ),
        ("not_an_upstream_parameter", True, "BULK_RNASEQ_ADVANCED_UNKNOWN"),
    ],
)
def test_advanced_conflicts_and_dangerous_parameters_fail_closed(
    name: str, value: object, code: str
):
    config = _config()
    config["advanced"] = {name: value}

    assert _error_code(_inputs(config=config)) == code


def test_unknown_advanced_key_is_not_reflected_in_public_issue():
    secret_key = "TOKEN_secret/private" + "x" * 4_096
    config = _config()
    config["advanced"] = {secret_key: True}

    result = BulkRnaSeqWorkflowAdapter().validate(_inputs(config=config))

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_ADVANCED_UNKNOWN"
    assert result.errors[0].path == "config.advanced"
    assert secret_key not in json.dumps(result.to_dict())


@pytest.mark.parametrize(
    ("name", "value"),
    [
        ("min_trimmed_reads", True),
        ("min_trimmed_reads", 0),
        ("min_trimmed_reads", "10"),
        ("min_mapped_reads", True),
        ("min_mapped_reads", 101),
        ("deseq2_vst", 1),
        ("gtf_group_features", "gene id"),
        ("rseqc_modules", "bam_stat,bam_stat"),
    ],
)
def test_advanced_wrong_types_and_stricter_semantics_fail_closed(
    name: str, value: object
):
    config = _config()
    config["advanced"] = {name: value}

    assert _error_code(_inputs(config=config)) == "BULK_RNASEQ_ADVANCED_INVALID"


@pytest.mark.parametrize(
    "mutate",
    [
        lambda config: config.update({"unexpected": True}),
        lambda config: config["standard"].update({"unexpected": True}),
        lambda config: config["standard"]["reference"].update(
            {"fasta_sha256": "invalid"}
        ),
        lambda config: config["standard"].update(
            {
                "umi": {
                    "enabled": True,
                    "mode": "read_sequence",
                    "deduplication_tool": "umitools",
                    "extraction_method": "regex",
                    "barcode_pattern": "(?P<umi>.{6})",
                    "discard_read": True,
                }
            }
        ),
    ],
)
def test_standard_unknown_fields_bad_digests_and_bool_integer_confusion_fail(mutate):
    config = _config()
    mutate(config)

    assert _error_code(_inputs(config=config)) in {
        "BULK_RNASEQ_CONFIG_UNKNOWN",
        "BULK_RNASEQ_STANDARD_INVALID",
    }


@pytest.mark.parametrize(
    ("rows", "code"),
    [
        (
            [
                _sample(lane="L001", strandedness="auto"),
                _sample(lane="L002", strandedness="reverse"),
            ],
            "BULK_RNASEQ_LANE_SEMANTICS_CONFLICT",
        ),
        (
            [_sample(), _sample()],
            "BULK_RNASEQ_LANE_DUPLICATE",
        ),
        (
            [
                _sample(sample="S1"),
                {
                    **_sample(sample="S2"),
                    "fastq_1": "/data/S1.lib1.L001.R1.fastq.gz",
                },
            ],
            "BULK_RNASEQ_FASTQ_DUPLICATE",
        ),
    ],
)
def test_lane_conflicts_and_duplicate_fastqs_fail_closed(rows, code):
    assert _error_code(_inputs(samples=rows)) == code


def test_layout_is_explicit_and_not_inferred_from_fastq_names():
    pe_missing_r2 = _sample(layout="PE")
    pe_missing_r2.pop("fastq_2")
    se_with_r2 = _sample(layout="SE")
    se_with_r2["fastq_2"] = "/data/S1.lib1.L001.R2.fastq.gz"
    no_layout = _sample()
    no_layout.pop("layout")

    for row in (pe_missing_r2, se_with_r2, no_layout):
        assert _error_code(_inputs(samples=[row])) == "BULK_RNASEQ_SAMPLES_INVALID"


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("fasta", "/refs/reference.txt"),
        ("fasta", "/refs/reference genome.fa"),
        ("fasta", "/refs/reference;touch.fa"),
        ("fasta", "/refs/../reference.fa"),
        ("gtf", "/refs/annotation.bed"),
        ("gtf", "/refs/gene annotation.gtf"),
    ],
)
def test_reference_paths_require_explicit_safe_fasta_and_gtf_suffixes(field, value):
    config = _config()
    config["standard"]["reference"][field] = value

    assert _error_code(_inputs(config=config)) == "BULK_RNASEQ_STANDARD_INVALID"


@pytest.mark.parametrize(
    "fastq",
    [
        "/data/read one.fastq.gz",
        "/data/read.fastq",
        "/data/read.txt",
    ],
)
def test_fastq_paths_use_controlled_gzip_contract(fastq):
    row = _sample()
    row["fastq_1"] = fastq

    assert _error_code(_inputs(samples=[row])) == "BULK_RNASEQ_SAMPLES_INVALID"


@pytest.mark.parametrize(
    "umi",
    [
        {
            "enabled": True,
            "mode": "read_name",
            "deduplication_tool": "umitools",
            "read_name_separator": ":",
            "grouping_method": "directional",
        },
        {
            "enabled": True,
            "mode": "read_sequence",
            "deduplication_tool": "umitools",
            "extraction_method": "regex",
            "barcode_pattern": "(?P<umi_1>.{6})",
            "barcode_pattern_2": "(?P<umi_2>.{6})",
            "discard_read": 2,
            "emit_dedup_stats": True,
            "primary_alignments_only": True,
        },
    ],
)
def test_typed_umi_modes_normalize_without_raw_cli_tokens(umi):
    result = BulkRnaSeqWorkflowAdapter().validate(_inputs(config=_config(umi=umi)))

    assert result.is_success
    assert result.value["nfcore_params"]["with_umi"] is True
    assert result.value["nfcore_params"]["umi_dedup_tool"] == "umitools"
    assert result.value["nfcore_params"]["skip_markduplicates"] is True
    assert "extra_star_align_args" not in result.value["nfcore_params"]


def test_string_mode_umi_pattern_uses_reviewed_symbolic_alphabet():
    umi = {
        "enabled": True,
        "mode": "read_sequence",
        "deduplication_tool": "umitools",
        "extraction_method": "string",
        "barcode_pattern": "NNNNXX",
    }

    result = BulkRnaSeqWorkflowAdapter().validate(_inputs(config=_config(umi=umi)))

    assert result.is_success
    assert result.value["nfcore_params"]["umitools_bc_pattern"] == "NNNNXX"


@pytest.mark.parametrize(
    "umi",
    [
        {
            "enabled": True,
            "mode": "read_sequence",
            "deduplication_tool": "umitools",
            "extraction_method": "regex",
            "barcode_pattern": "';touch/tmp/pwn#",
        },
        {
            "enabled": True,
            "mode": "read_sequence",
            "deduplication_tool": "umitools",
            "extraction_method": "regex",
            "barcode_pattern": "(?P<discard_1>.{6})",
        },
        {
            "enabled": True,
            "mode": "read_sequence",
            "deduplication_tool": "umitools",
            "extraction_method": "string",
            "barcode_pattern": "NN.?",
        },
        {
            "enabled": True,
            "mode": "read_name",
            "deduplication_tool": "umitools",
            "read_name_separator": "'",
        },
    ],
)
def test_umi_patterns_and_separators_reject_shell_like_or_unreviewed_strings(umi):
    assert _error_code(_inputs(config=_config(umi=umi))) in {
        "BULK_RNASEQ_STANDARD_INVALID",
        "BULK_RNASEQ_UMI_PATTERN_INVALID",
    }


@pytest.mark.parametrize(
    "standard_update",
    [
        {"umi": {"enabled": False, "grouping_method": "directional"}},
        {
            "umi": {
                "enabled": True,
                "mode": "read_name",
                "deduplication_tool": "umitools",
                "read_name_separator": ":",
                "barcode_pattern": "NNNN",
            }
        },
        {
            "umi": {
                "enabled": True,
                "mode": "read_sequence",
                "deduplication_tool": "umitools",
                "extraction_method": "string",
                "barcode_pattern": "NNNN",
                "read_name_separator": ":",
            }
        },
        {"outputs": {"umi_intermediates": True}},
        {"qc": {"enabled": False, "fastqc": True}},
        {
            "trimming": {"enabled": False, "tool": "fastp"},
            "outputs": {"trimmed_reads": True},
        },
    ],
)
def test_standard_semantic_conflicts_are_rejected(standard_update):
    assert _error_code(_inputs(config=_config(**standard_update))) in {
        "BULK_RNASEQ_UMI_CONFLICT",
        "BULK_RNASEQ_OUTPUT_CONFLICT",
        "BULK_RNASEQ_QC_CONFLICT",
    }


def test_qc_master_switch_disables_every_qc_owned_step():
    result = BulkRnaSeqWorkflowAdapter().validate(
        _inputs(config=_config(qc={"enabled": False}))
    )

    assert result.is_success
    params = result.value["nfcore_params"]
    assert params["skip_qc"] is True
    for name in (
        "skip_fastqc",
        "skip_multiqc",
        "skip_rseqc",
        "skip_qualimap",
        "skip_dupradar",
        "skip_biotype_qc",
        "skip_deseq2_qc",
        "skip_preseq",
        "skip_markduplicates",
    ):
        assert params[name] is True


@pytest.mark.parametrize(
    "samples",
    [
        [_sample(layout="SE")],
        [_sample(sample="S1", layout="PE"), _sample(sample="S2", layout="SE")],
    ],
)
@pytest.mark.parametrize("r2_field", ["barcode_pattern_2", "discard_read"])
def test_r2_specific_umi_options_require_all_samples_to_be_paired(samples, r2_field):
    umi = {
        "enabled": True,
        "mode": "read_sequence",
        "deduplication_tool": "umitools",
        "extraction_method": "regex",
        "barcode_pattern": "(?P<umi_1>.{6})",
        r2_field: "(?P<umi_2>.{6})" if r2_field == "barcode_pattern_2" else 2,
    }

    assert _error_code(_inputs(config=_config(umi=umi), samples=samples)) == (
        "BULK_RNASEQ_UMI_LAYOUT_CONFLICT"
    )


@pytest.mark.parametrize(
    ("advanced", "standard_update", "strandedness"),
    [
        (
            {"min_trimmed_reads": 10},
            {"trimming": {"enabled": False, "tool": "fastp"}},
            "auto",
        ),
        ({"star_ignore_sjdbgtf": True}, {}, "auto"),
        ({"stringtie_ignore_gtf": True}, {}, "auto"),
        (
            {"rseqc_modules": "bam_stat"},
            {"qc": {"enabled": True, "rseqc": False}},
            "auto",
        ),
        ({"deseq2_vst": True}, {"qc": {"enabled": True, "deseq2_pca": False}}, "auto"),
        ({"stranded_threshold": 0.8}, {}, "reverse"),
        ({"stranded_threshold": 0.6, "unstranded_threshold": 0.7}, {}, "auto"),
        (
            {"featurecounts_feature_type": "exon"},
            {"qc": {"enabled": True, "biotype": False}},
            "auto",
        ),
        ({"featurecounts_group_type": "gene_biotype"}, {}, "auto"),
    ],
)
def test_advanced_context_conflicts_fail_closed(
    advanced, standard_update, strandedness
):
    config = _config(**standard_update)
    config["advanced"] = advanced

    assert _error_code(
        _inputs(config=config, samples=[_sample(strandedness=strandedness)])
    ) in {
        "BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT",
        "BULK_RNASEQ_ADVANCED_INVALID",
    }


def test_options_and_server_path_samples_are_not_accepted():
    assert _error_code(_inputs(options={"strict": True})) == (
        "BULK_RNASEQ_OPTIONS_INVALID"
    )
    assert (
        _error_code(WorkflowInputs(config=_config(), samples="/tmp/samples.csv"))
        == "BULK_RNASEQ_SAMPLES_INVALID"
    )


@pytest.mark.parametrize(
    "config",
    [
        {1: "not-a-json-key"},
        {"standard": {"reference": _reference()}, "advanced": {1: True}},
        {
            "standard": {"reference": _reference()},
            "advanced": {"min_mapped_reads": float("nan")},
        },
    ],
)
def test_non_json_config_values_fail_closed_without_exceptions(config):
    assert _error_code(_inputs(config=config)) == "BULK_RNASEQ_CONFIG_INVALID"


def test_bulk_and_encode_adapters_coexist_only_in_explicit_test_registry():
    registry = WorkflowRegistry(
        (EncodeStyleWorkflowAdapter(), BulkRnaSeqWorkflowAdapter())
    )

    assert [item.workflow_id for item in registry.list_metadata()] == [
        "encode-style-chipseq-cuttag-atac-mnase",
        "bulk-rnaseq",
    ]
    assert registry.get("bulk-rnaseq").metadata.engines == ("nextflow",)


def test_default_registry_remains_encode_only():
    registry = create_default_workflow_registry()

    assert [item.workflow_id for item in registry.list_metadata()] == [
        "encode-style-chipseq-cuttag-atac-mnase"
    ]
    assert registry.has("bulk-rnaseq") is False
