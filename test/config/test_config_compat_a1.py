"""Single-path config contract tests.

Omitted selectors resolve to one analysis path. analysis_profile and
remove_dup are rejected instead of being remapped.
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


def _assert_single_path_defaults(validated):
    assert "analysis_profile" not in validated
    assert "remove_dup" not in validated
    assert validated["blacklist_mode"] == "analysis"
    assert validated["qc"]["blacklist_filter"] is False
    assert validated["duplicate_analysis"] == "auto"
    assert validated["organelle_mode"] == "auto"
    assert validated["control_mode"] == "matched"
    assert validated["duplicate_backend"] == "picard"
    assert validated["reference_mode"] == "legacy"
    assert validated["external_control_policy"] == "warn"
    assert validated["fragment_cpm"] == "auto"
    assert validated["tss_fragment_profile"] is False


def test_omitted_version_uses_the_single_path(tmp_path):
    validated = _validate(tmp_path)

    assert validated["config_schema_version"] == "1.0.0"
    _assert_single_path_defaults(validated)
    validate_execution_support(validated)


def test_explicit_v2_uses_the_same_path(tmp_path):
    validated = _validate(tmp_path, config_schema_version="2.0.0")

    assert validated["config_schema_version"] == "2.0.0"
    _assert_single_path_defaults(validated)
    validate_execution_support(validated)


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


def test_blacklist_alias_without_profile_maps_to_copies(tmp_path):
    validated = _validate(tmp_path, qc={"blacklist_filter": True})
    assert validated["blacklist_mode"] == "copies"
    assert validated["qc"]["blacklist_filter"] is True


def test_analysis_blacklist_branch_is_executable(tmp_path):
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
    validated = _validate(
        tmp_path,
        blacklist_mode="analysis",
        genome_resources={
            "hs": {
                "effective_genome_size": "hs",
                "blacklist": str(blacklist),
            }
        },
    )
    validate_execution_support(validated)
    assert validated["blacklist_mode"] == "analysis"
    assert validated["qc"]["blacklist_filter"] is False


@pytest.mark.parametrize("key", ["analysis_profile", "remove_dup"])
@pytest.mark.parametrize("value", ["legacy", "updated", "auto", "yes", "no"])
def test_removed_selectors_are_rejected(tmp_path, key, value):
    with pytest.raises(ValidationError, match=key) as exc:
        _validate(tmp_path, **{key: value})
    message = str(exc.value)
    assert "duplicate_analysis" in message
    assert "auto" in message


def test_explicit_duplicate_analysis_is_idempotent(tmp_path):
    first = _validate(tmp_path, duplicate_analysis="retain")
    assert first["duplicate_analysis"] == "retain"
    assert "remove_dup" not in first
    second = validate_config(first)
    assert first == second
    validate_execution_support(second)


def test_duplicate_analysis_auto_stays_unresolved(tmp_path):
    validated = _validate(tmp_path, duplicate_analysis="auto")
    assert validated["duplicate_analysis"] == "auto"
    validate_execution_support(validated)


@pytest.mark.parametrize(
    ("key", "value"),
    [
        ("control_mode", "pooled"),
        ("control_mode", "auto"),
        ("duplicate_backend", "samtools"),
        ("duplicate_backend", ""),
    ],
)
def test_non_executable_control_and_backend_are_rejected(tmp_path, key, value):
    with pytest.raises(ValidationError, match=key):
        _validate(tmp_path, **{key: value})


def test_tss_fragment_profile_requires_fragment_cpm(tmp_path):
    validated = _validate(tmp_path, tss_fragment_profile=True)
    assert validated["tss_fragment_profile"] is True
    validate_execution_support(validated)
    with pytest.raises(ValidationError, match="requires fragment_cpm != off"):
        _validate(
            tmp_path,
            fragment_cpm="off",
            tss_fragment_profile=True,
        )


@pytest.mark.parametrize("mode", ["off", "auto", "on"])
def test_fragment_cpm_modes_accepted(tmp_path, mode):
    validated = _validate(tmp_path, fragment_cpm=mode)
    assert validated["fragment_cpm"] == mode
    if mode != "off":
        validate_execution_support(validated)


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


def test_chipseq_narrow_alias_keeps_idr_parameters_and_is_idempotent(tmp_path):
    validated = _validate(
        tmp_path,
        chipseq_idr=False,
        idr={"threshold": 0.2, "rank": "signal.value", "seed": 123},
        reproducibility={"enabled": True, "idr": {"chipseq_narrow": True}},
    )
    assert validated["chipseq_idr"] is True
    assert validated["idr"] == {
        "threshold": 0.2,
        "rank": "signal.value",
        "seed": 123,
    }
    renormalized = validate_config(validated)
    assert renormalized["chipseq_idr"] is True
    assert renormalized["idr"] == validated["idr"]

    workspace_config = _prepare_workspace_config(validated)
    rendered = yaml.safe_load(_render_config_yaml(workspace_config))
    rendered["samples"] = validated["samples"]
    revalidated = validate_config(rendered)
    assert revalidated["chipseq_idr"] is True
    assert revalidated["idr"] == validated["idr"]


def test_chipseq_narrow_alias_requires_replicate_analysis(tmp_path):
    with pytest.raises(ValidationError, match="requires replicate_analysis=true"):
        _validate(
            tmp_path,
            chipseq_idr=False,
            replicate_analysis=False,
            reproducibility={"enabled": True, "idr": {"chipseq_narrow": True}},
        )


def test_chipseq_narrow_alias_rejects_invalid_idr_parameters(tmp_path):
    with pytest.raises(ValidationError, match="idr.threshold"):
        _validate(
            tmp_path,
            chipseq_idr=False,
            idr={"threshold": 2.0},
            reproducibility={"enabled": True, "idr": {"chipseq_narrow": True}},
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
        "blacklist_mode",
        "duplicate_analysis",
        "organelle_mode",
        "control_mode",
        "duplicate_backend",
        "reference_mode",
        "external_control_policy",
        "fragment_cpm",
        "tss_fragment_profile",
    )
    assert "analysis_profile" not in revalidated
    assert "remove_dup" not in revalidated
    for key in keys:
        assert revalidated[key] == validated[key], key
    assert revalidated["qc"]["blacklist_filter"] is validated["qc"]["blacklist_filter"]


def test_analysis_mode_normalizes_and_renders_idempotently(tmp_path):
    blacklist = tmp_path / "blacklist.bed"
    blacklist.write_text("chr1\t0\t10\n", encoding="utf-8")
    resources = {
        "hs": {
            "effective_genome_size": "hs",
            "blacklist": str(blacklist),
        }
    }
    validated = _validate(
        tmp_path,
        blacklist_mode="analysis",
        genome_resources=resources,
    )
    assert validated["blacklist_mode"] == "analysis"
    assert validated["qc"]["blacklist_filter"] is False

    renormalized = validate_config(validated)
    assert renormalized["blacklist_mode"] == "analysis"
    assert renormalized["qc"]["blacklist_filter"] is False

    workspace_config = _prepare_workspace_config(validated)
    rendered = yaml.safe_load(_render_config_yaml(workspace_config))
    rendered["samples"] = validated["samples"]
    revalidated = validate_config(rendered)
    assert revalidated["blacklist_mode"] == "analysis"
    assert revalidated["qc"]["blacklist_filter"] is False


def test_validation_has_no_filesystem_side_effects(tmp_path):
    samples = _write_samples(tmp_path)
    before = sorted(p.name for p in tmp_path.iterdir())
    validate_config({"samples": str(samples), "use_control": False})
    after = sorted(p.name for p in tmp_path.iterdir())
    assert after == before


def test_authoring_schema_exposes_the_single_path():
    schema = build_encode_authoring_schema()
    assert schema.schema_version == SCHEMA_VERSION
    properties = schema.config_schema["properties"]

    for key in (
        "config_schema_version",
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
    assert "analysis_profile" not in properties
    assert "remove_dup" not in properties

    assert properties["config_schema_version"]["default"] == "1.0.0"
    assert properties["blacklist_mode"]["default"] == "analysis"
    assert properties["organelle_mode"]["default"] == "auto"
    assert properties["duplicate_backend"]["default"] == "picard"
    assert properties["duplicate_backend"]["enum"] == ["picard"]
    assert properties["duplicate_analysis"]["default"] == "auto"
    assert properties["control_mode"]["default"] == "matched"
    assert properties["control_mode"]["enum"] == ["matched"]
    assert properties["fragment_cpm"]["default"] == "auto"
    assert "blacklist_filter" not in properties["qc"]["properties"]


def test_default_template_canonical_copies_override(tmp_path):
    from pathlib import Path

    import yaml

    root = Path(__file__).resolve().parents[2]
    config = yaml.safe_load((root / "config" / "config.yaml").read_text())
    config["samples"] = str(tmp_path / "samples.tsv")
    Path(config["samples"]).write_text(
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\n"
        "S1\tR1.fq\tR2.fq\tPE\tchipseq\tCTCF\tnarrow\ths\t/idx\n",
        encoding="utf-8",
    )
    config["outdir"] = str(tmp_path / "results")
    config["blacklist_mode"] = "copies"
    validated = validate_config(config)
    assert validated["blacklist_mode"] == "copies"
    assert validated["qc"]["blacklist_filter"] is True
