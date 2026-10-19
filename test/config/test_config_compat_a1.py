"""A1 config compatibility tests.

Pin the additive compatibility surface without activating new analysis
behavior. Legacy configs must keep their parsed semantics and paths.
"""

import yaml

import pytest

from encode_pipeline.adapters.encode import (
    _prepare_workspace_config,
    _render_config_yaml,
)
from encode_pipeline.adapters.encode_authoring import (
    SCHEMA_VERSION,
    build_encode_authoring_schema,
)
from encode_pipeline.config.validate import (
    ValidationError,
    validate_config,
    validate_execution_support,
)


SAMPLES_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
    "bowtie2_index\n"
)
SAMPLES_ROW = "S1\tR1.fq\tR2.fq\tPE\tchipseq\tT\tnarrow\ths\tidx\n"


def _write_samples(tmp_path):
    path = tmp_path / "samples.tsv"
    path.write_text(SAMPLES_HEADER + SAMPLES_ROW, encoding="utf-8")
    return path


def _validate(tmp_path, **overrides):
    config = {"samples": str(_write_samples(tmp_path)), "use_control": False}
    config.update(overrides)
    return validate_config(config)


def test_omitted_new_keys_resolve_to_legacy(tmp_path):
    validated = _validate(tmp_path)

    assert validated["config_schema_version"] == "1.0.0"
    assert validated["analysis_profile"] == "legacy"
    assert validated["blacklist_mode"] == "copies"
    assert validated["qc"]["blacklist_filter"] is True
    assert validated["duplicate_analysis"] == "exclude"
    assert validated["organelle_mode"] == "keep"
    assert validated["control_mode"] == "pooled"
    assert validated["duplicate_backend"] == ""
    assert validated["reference_mode"] == "legacy"
    assert validated["external_control_policy"] == "warn"
    assert validated["fragment_cpm"] == "off"
    assert validated["tss_fragment_profile"] is False


def test_explicit_v2_resolves_to_updated_raw(tmp_path):
    validated = _validate(tmp_path, config_schema_version="2.0.0")

    assert validated["config_schema_version"] == "2.0.0"
    assert validated["analysis_profile"] == "updated"
    assert validated["blacklist_mode"] == "off"
    assert validated["qc"]["blacklist_filter"] is False
    assert validated["duplicate_analysis"] == "exclude"
    assert validated["organelle_mode"] == "auto"
    assert validated["control_mode"] == "matched"
    assert validated["duplicate_backend"] == "picard"


@pytest.mark.parametrize(
    ("qc_value", "expected_mode", "expected_alias"),
    [(True, "copies", True), (False, "off", False)],
)
def test_blacklist_alias_true_false(tmp_path, qc_value, expected_mode, expected_alias):
    validated = _validate(tmp_path, qc={"blacklist_filter": qc_value})

    assert validated["blacklist_mode"] == expected_mode
    assert validated["qc"]["blacklist_filter"] is expected_alias


def test_blacklist_alias_and_canonical_consistent(tmp_path):
    validated = _validate(
        tmp_path,
        qc={"blacklist_filter": True},
        blacklist_mode="copies",
    )
    assert validated["blacklist_mode"] == "copies"
    assert validated["qc"]["blacklist_filter"] is True


def test_blacklist_alias_and_canonical_conflict_rejected(tmp_path):
    with pytest.raises(ValidationError, match="blacklist_mode conflicts"):
        _validate(
            tmp_path,
            qc={"blacklist_filter": True},
            blacklist_mode="off",
        )


def test_updated_alias_true_without_canonical_key_maps_to_copies(tmp_path):
    validated = _validate(
        tmp_path,
        analysis_profile="updated",
        qc={"blacklist_filter": True},
    )
    assert validated["blacklist_mode"] == "copies"
    assert validated["qc"]["blacklist_filter"] is True


def test_unsupported_analysis_branch_rejected(tmp_path):
    with pytest.raises(ValidationError, match="not implemented in A1"):
        _validate(tmp_path, blacklist_mode="analysis")


def test_remove_dup_conflicts_with_duplicate_analysis(tmp_path):
    with pytest.raises(ValidationError, match="remove_dup=no conflicts"):
        _validate(
            tmp_path,
            analysis_profile="updated",
            remove_dup="no",
            duplicate_analysis="exclude",
        )
    with pytest.raises(ValidationError, match="remove_dup=yes conflicts"):
        _validate(
            tmp_path,
            analysis_profile="updated",
            remove_dup="yes",
            duplicate_analysis="retain",
        )


def test_remove_dup_auto_does_not_override_explicit_duplicate_analysis(tmp_path):
    validated = _validate(
        tmp_path,
        analysis_profile="updated",
        remove_dup="auto",
        duplicate_analysis="retain",
    )
    assert validated["duplicate_analysis"] == "retain"


def test_updated_remove_dup_no_requires_explicit_policy(tmp_path):
    with pytest.raises(
        ValidationError, match="requires an explicit duplicate_analysis"
    ):
        _validate(tmp_path, analysis_profile="updated", remove_dup="no")


def test_updated_remove_dup_no_and_retain_accepted(tmp_path):
    validated = _validate(
        tmp_path,
        analysis_profile="updated",
        remove_dup="no",
        duplicate_analysis="retain",
    )
    assert validated["duplicate_analysis"] == "retain"


def test_legacy_duplicate_policy_is_derived(tmp_path):
    validated = _validate(tmp_path, remove_dup="no")
    assert validated["duplicate_analysis"] == "retain"


@pytest.mark.parametrize(
    "overrides",
    [
        {"fragment_cpm": "on"},
        {"fragment_cpm": "auto"},
        {"tss_fragment_profile": True},
        {"analysis_profile": "legacy", "control_mode": "matched"},
        {"analysis_profile": "legacy", "duplicate_backend": "picard"},
    ],
)
def test_unimplemented_or_legacy_incompatible_modes_rejected(tmp_path, overrides):
    with pytest.raises(ValidationError):
        _validate(tmp_path, **overrides)


def test_organelle_chroms_resource_is_accepted(tmp_path):
    validated = _validate(
        tmp_path,
        genome_resources={
            "hs": {
                "effective_genome_size": "hs",
                "organelle_chroms": ["chrM", "MT"],
            }
        },
    )

    assert validated["genome_resources"]["hs"]["organelle_chroms"] == [
        "chrM",
        "MT",
    ]


@pytest.mark.parametrize(
    "organelle_chroms",
    [
        "chrM",
        ["chrM", "chrM"],
        ["chrM", ""],
        ["chrM", 1],
    ],
)
def test_organelle_chroms_resource_rejects_invalid_values(tmp_path, organelle_chroms):
    with pytest.raises(ValidationError, match="organelle_chroms"):
        _validate(
            tmp_path,
            genome_resources={
                "hs": {
                    "effective_genome_size": "hs",
                    "organelle_chroms": organelle_chroms,
                }
            },
        )


def test_workspace_render_roundtrip_is_idempotent(tmp_path):
    validated = _validate(tmp_path, config_schema_version="2.0.0")
    workspace_config = _prepare_workspace_config(validated)
    rendered = yaml.safe_load(_render_config_yaml(workspace_config))

    # The adapter rewrites samples/outdir for the workspace; restore the
    # concrete samples path before re-validating.
    rendered["samples"] = validated["samples"]
    revalidated = validate_config(rendered)

    keys = (
        "config_schema_version",
        "analysis_profile",
        "blacklist_mode",
        "remove_dup",
        "duplicate_analysis",
        "organelle_mode",
        "control_mode",
        "duplicate_backend",
        "reference_mode",
        "external_control_policy",
        "fragment_cpm",
        "tss_fragment_profile",
    )
    for key in keys:
        assert revalidated[key] == validated[key], key
    assert revalidated["qc"]["blacklist_filter"] is validated["qc"]["blacklist_filter"]


def test_validation_has_no_filesystem_side_effects(tmp_path):
    samples = _write_samples(tmp_path)
    before = sorted(p.name for p in tmp_path.iterdir())
    validate_config({"samples": str(samples), "use_control": False})
    after = sorted(p.name for p in tmp_path.iterdir())
    assert after == before


def test_authoring_schema_exposes_a1_keys_and_executable_legacy_defaults():
    schema = build_encode_authoring_schema()
    assert schema.schema_version == SCHEMA_VERSION
    properties = schema.config_schema["properties"]

    for key in (
        "config_schema_version",
        "analysis_profile",
        "blacklist_mode",
        "organelle_mode",
        "duplicate_backend",
        "duplicate_analysis",
        "control_mode",
        "reference_mode",
        "external_control_policy",
        "fragment_cpm",
        "tss_fragment_profile",
    ):
        assert key in properties

    assert properties["config_schema_version"]["default"] == "1.0.0"
    assert properties["analysis_profile"]["default"] == "legacy"
    assert properties["blacklist_mode"]["default"] == "copies"
    assert properties["duplicate_backend"]["default"] == ""
    assert properties["control_mode"]["default"] == "pooled"
    assert properties["qc"]["properties"]["blacklist_filter"]["default"] is True
    assert "" in properties["duplicate_backend"]["enum"]


def test_updated_control_mode_auto_rejected_at_execution(tmp_path):
    validated = _validate(
        tmp_path,
        analysis_profile="updated",
        control_mode="auto",
    )
    with pytest.raises(ValidationError, match="requires control_mode=matched"):
        validate_execution_support(validated)
