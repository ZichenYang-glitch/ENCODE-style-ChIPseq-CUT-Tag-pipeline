"""Hi-TrAC authoring contracts require no scientific tools or file contents."""

from __future__ import annotations

from copy import deepcopy
import json
from pathlib import Path

from jsonschema import Draft202012Validator
import pytest

from encode_pipeline.adapters.hitrac_preprocess.authoring import (
    build_hitrac_authoring_schema,
)
from encode_pipeline.adapters.hitrac_preprocess.validation import validate_hitrac_inputs
from encode_pipeline.platform.adapters import WorkflowInputs


def _sample(sample_id="sample 1", **changes):
    row = {
        "sample_id": sample_id,
        "fastq_1": "/operator inputs/reads_1.fastq.gz",
        "fastq_2": "/operator inputs/reads_2.fastq.gz",
    }
    return {**row, **changes}


def _inputs(*, samples=None, options=None, config=None):
    return WorkflowInputs(
        config={} if config is None else config,
        samples=[_sample()] if samples is None else samples,
        options={} if options is None else options,
    )


def test_schema_declares_only_author_input_and_fresh_copies():
    schema = build_hitrac_authoring_schema()
    for field in ("config_schema", "sample_schema", "option_schema"):
        Draft202012Validator.check_schema(getattr(schema, field))
    assert schema.coverage.to_dict() == {
        "config": "complete",
        "samples": "complete",
        "options": "complete",
    }
    assert schema.input_modes.to_dict() == {
        "config": ["object"],
        "samples": ["inline_rows"],
        "options": ["object"],
    }
    assert schema.config_schema["properties"] == {}
    assert set(schema.option_schema["properties"]) == {"threads", "mapq"}
    assert set(schema.sample_schema["items"]["properties"]) == {
        "sample_id",
        "fastq_1",
        "fastq_2",
    }
    schema.option_schema["properties"]["threads"]["maximum"] = 999
    assert (
        build_hitrac_authoring_schema().option_schema["properties"]["threads"][
            "maximum"
        ]
        == 8
    )


def test_valid_multi_sample_order_defaults_and_source_preservation(monkeypatch):
    submitted = _inputs(samples=[_sample("zeta"), _sample("alpha")])
    before = deepcopy(submitted.to_dict())

    def forbidden_io(*args, **kwargs):
        raise AssertionError("validation must not inspect submitted files")

    monkeypatch.setattr(Path, "open", forbidden_io)
    monkeypatch.setattr(Path, "stat", forbidden_io)
    result = validate_hitrac_inputs(submitted)
    assert result.is_success
    assert result.value.options == {"threads": 2, "mapq": 10}
    assert [row["sample_id"] for row in result.value.samples] == ["zeta", "alpha"]
    assert submitted.to_dict() == before
    result.value.samples[0]["sample_id"] = "changed copy"
    assert submitted.to_dict() == before


@pytest.mark.parametrize(
    "options",
    [
        {},
        {"threads": 2, "mapq": 0},
        {"threads": 3, "mapq": 17},
        {"threads": 8, "mapq": 255},
        {"threads": 2.0, "mapq": 10.0},
    ],
)
def test_option_schema_and_validator_agree_on_valid_values(options):
    schema = build_hitrac_authoring_schema()
    assert Draft202012Validator(schema.option_schema).is_valid(options)
    result = validate_hitrac_inputs(_inputs(options=options))
    assert result.is_success
    assert type(result.value.options["threads"]) is int
    assert type(result.value.options["mapq"]) is int


@pytest.mark.parametrize(
    ("options", "path"),
    [
        ({"threads": 1}, "options.threads"),
        ({"threads": 9}, "options.threads"),
        ({"threads": True}, "options.threads"),
        ({"threads": 2.5}, "options.threads"),
        ({"threads": "2"}, "options.threads"),
        ({"mapq": -1}, "options.mapq"),
        ({"mapq": 256}, "options.mapq"),
        ({"mapq": False}, "options.mapq"),
        ({"mapq": 10.5}, "options.mapq"),
        ({"mapq": None}, "options.mapq"),
        ({"timeout": 60}, "options"),
        ({"executable": "/private marker/python"}, "options"),
    ],
)
def test_option_schema_and_validation_reject_bad_options(options, path):
    schema = build_hitrac_authoring_schema()
    assert not Draft202012Validator(schema.option_schema).is_valid(options)
    result = validate_hitrac_inputs(_inputs(options=options))
    assert result.is_failure
    assert result.issues[0].path == path
    assert result.issues[0].hint
    assert "/private marker" not in json.dumps(result.to_dict())


@pytest.mark.parametrize("value", [float("nan"), float("inf"), -float("inf")])
def test_nonfinite_options_never_enter_science(value):
    assert validate_hitrac_inputs(_inputs(options={"threads": value})).is_failure


@pytest.mark.parametrize("sample_id", ["A", "a_1.2-3", "sample name", "a" * 128])
def test_h2_display_id_alphabet_and_length_remain_supported(sample_id):
    row = _sample(sample_id)
    schema = build_hitrac_authoring_schema()
    assert Draft202012Validator(schema.sample_schema).is_valid([row])
    assert validate_hitrac_inputs(_inputs(samples=[row])).is_success


@pytest.mark.parametrize("sample_id", ["", "_leading", "a/b", "a;cmd", "a" * 129])
def test_invalid_sample_ids_have_fixed_field_location(sample_id):
    row = _sample(sample_id)
    schema = build_hitrac_authoring_schema()
    assert not Draft202012Validator(schema.sample_schema).is_valid([row])
    result = validate_hitrac_inputs(_inputs(samples=[row]))
    assert result.issues[0].code == "HITRAC_SAMPLE_ID_INVALID"
    assert result.issues[0].path == "samples[0].sample_id"


@pytest.mark.parametrize("missing", ["sample_id", "fastq_1", "fastq_2"])
def test_missing_mate_or_id_is_located_without_echoing_other_values(missing):
    row = _sample()
    del row[missing]
    result = validate_hitrac_inputs(_inputs(samples=[row]))
    assert result.is_failure
    assert result.issues[0].path == f"samples[0].{missing}"
    assert result.issues[0].code == "HITRAC_SAMPLE_FIELD_REQUIRED"


def test_duplicate_id_semantic_check_preserves_order_and_does_not_merge_lanes():
    rows = [_sample("same"), _sample("same", fastq_1="/lane2/one.gz")]
    # JSON Schema describes per-row structure; ID uniqueness is semantic validation.
    schema = build_hitrac_authoring_schema()
    assert Draft202012Validator(schema.sample_schema).is_valid(rows)
    result = validate_hitrac_inputs(_inputs(samples=rows))
    assert result.issues[0].code == "HITRAC_SAMPLE_ID_DUPLICATE"
    assert result.issues[0].path == "samples[1].sample_id"
    assert rows[1]["fastq_1"] == "/lane2/one.gz"


@pytest.mark.parametrize(
    "path",
    [
        "relative.gz",
        "/",
        "//server/file.gz",
        "/a//b.gz",
        "/a/../b.gz",
        "/a/./b.gz",
        "/a/file.gz/",
        "/a/back\\slash.gz",
    ],
)
def test_external_path_syntax_rejects_ambiguous_paths(path):
    row = _sample(fastq_1=path)
    schema = build_hitrac_authoring_schema()
    assert not Draft202012Validator(schema.sample_schema).is_valid([row])
    result = validate_hitrac_inputs(_inputs(samples=[row]))
    assert result.issues[0].path == "samples[0].fastq_1"
    assert result.issues[0].code == "HITRAC_FASTQ_PATH_INVALID"


@pytest.mark.parametrize("surface", ["config", "options", "samples"])
def test_unknown_fields_and_values_are_not_reflected(surface):
    secret = "PRIVATE_ENV_VALUE_and_path"
    kwargs = {}
    if surface == "samples":
        kwargs["samples"] = [_sample(**{secret: "/private/path"})]
    else:
        kwargs[surface] = {secret: "/private/path"}
    result = validate_hitrac_inputs(_inputs(**kwargs))
    assert result.is_failure
    assert result.issues[0].path == ("samples[0]" if surface == "samples" else surface)
    rendered = json.dumps(result.to_dict())
    assert secret not in rendered
    assert "/private/path" not in rendered
    assert result.issues[0].technical_message is None
    assert result.issues[0].context == {}


@pytest.mark.parametrize(
    "samples", [None, "/server/samples.tsv", Path("/server/samples.tsv")]
)
def test_other_sample_representations_are_not_silently_loaded(samples):
    result = validate_hitrac_inputs(WorkflowInputs(config={}, samples=samples))
    assert result.issues[0].code == "HITRAC_SAMPLES_INVALID"
    assert result.issues[0].path == "samples"
