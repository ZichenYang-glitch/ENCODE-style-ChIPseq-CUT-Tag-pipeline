"""Direct-API characterization tests for MNase config validation behavior.

These tests pin the current behavior of `_validate_mnase_config` so that
PR66 can safely extract MNase helpers without changing behavior.
"""

import pytest

from encode_pipeline.config.validate import ValidationError, validate_config


SAMPLES_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
)
SAMPLES_ROW = "M1\tR1.fq\tR2.fq\tPE\tmnase\tH3\tnucleosome\ths\tidx\n"


MNASE_DEFAULTS = {
    "mono_range": [140, 200],
    "fragments": {
        "sub": [1, 139],
        "mono": [140, 200],
        "di": [300, 400],
    },
    "dyad_range": [130, 200],
    "callers": {"danpos3": False, "inps": False, "sem": False},
}


def _make_config(tmp_path, **overrides):
    """Return a minimal valid config dict pointing at a temporary sample TSV."""
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(SAMPLES_HEADER + SAMPLES_ROW, encoding="utf-8")
    config = {"samples": str(samples_path), "use_control": False}
    config.update(overrides)
    return config


# ---------------------------------------------------------------------------
# Defaults and structure
# ---------------------------------------------------------------------------


def test_mnase_omitted_expands_to_defaults(tmp_path):
    validated = validate_config(_make_config(tmp_path))
    assert validated["mnase"] == MNASE_DEFAULTS


def test_mnase_empty_mapping_expands_to_defaults(tmp_path):
    validated = validate_config(_make_config(tmp_path, mnase={}))
    assert validated["mnase"] == MNASE_DEFAULTS


def test_mnase_ranges_return_lists_not_tuples(tmp_path):
    validated = validate_config(_make_config(tmp_path))
    assert isinstance(validated["mnase"]["mono_range"], list)
    assert isinstance(validated["mnase"]["dyad_range"], list)
    for value in validated["mnase"]["fragments"].values():
        assert isinstance(value, list)


@pytest.mark.parametrize("raw", ["yes", ["x"], True])
def test_mnase_must_be_mapping(tmp_path, raw):
    with pytest.raises(ValidationError, match="mnase must be a mapping"):
        validate_config(_make_config(tmp_path, mnase=raw))


def test_mnase_unknown_top_level_key_rejected(tmp_path):
    with pytest.raises(ValidationError, match="mnase: unknown key"):
        validate_config(_make_config(tmp_path, mnase={"unknown_key": True}))


# ---------------------------------------------------------------------------
# mono_range
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw",
    [
        [140, 200],
        (140, 200),
        ["140", "200"],
    ],
)
def test_mnase_mono_range_accepts_valid_inputs(tmp_path, raw):
    validated = validate_config(_make_config(tmp_path, mnase={"mono_range": raw}))
    assert validated["mnase"]["mono_range"] == [140, 200]
    assert isinstance(validated["mnase"]["mono_range"], list)


@pytest.mark.parametrize(
    "raw,label",
    [
        ([True, 200], "bool endpoint"),
        ([140.5, 200], "float endpoint"),
        ("140,200", "non list/tuple"),
    ],
)
def test_mnase_mono_range_rejects_invalid_element_types(
    tmp_path, raw, label
):
    if isinstance(raw, str):
        pattern = "mnase.mono_range must be a list of 2 positive ints"
    else:
        pattern = "mnase.mono_range elements must be integers"
    with pytest.raises(ValidationError, match=pattern):
        validate_config(_make_config(tmp_path, mnase={"mono_range": raw}))


@pytest.mark.parametrize("raw", [[140], [140, 200, 250]])
def test_mnase_mono_range_rejects_wrong_length(tmp_path, raw):
    with pytest.raises(
        ValidationError, match="mnase.mono_range must have exactly 2 elements"
    ):
        validate_config(_make_config(tmp_path, mnase={"mono_range": raw}))


@pytest.mark.parametrize("raw", [[0, 200], [140, 0], [-1, 200]])
def test_mnase_mono_range_rejects_non_positive_values(tmp_path, raw):
    with pytest.raises(
        ValidationError, match="mnase.mono_range values must be positive"
    ):
        validate_config(_make_config(tmp_path, mnase={"mono_range": raw}))


@pytest.mark.parametrize("raw", [[200, 200], [200, 100]])
def test_mnase_mono_range_rejects_min_not_less_than_max(tmp_path, raw):
    with pytest.raises(
        ValidationError, match="mnase.mono_range: min must be < max"
    ):
        validate_config(_make_config(tmp_path, mnase={"mono_range": raw}))


# ---------------------------------------------------------------------------
# fragments
# ---------------------------------------------------------------------------


def test_mnase_fragments_empty_mapping_uses_defaults(tmp_path):
    validated = validate_config(_make_config(tmp_path, mnase={"fragments": {}}))
    assert validated["mnase"]["fragments"] == MNASE_DEFAULTS["fragments"]


@pytest.mark.parametrize("raw", ["yes", ["x"], True])
def test_mnase_fragments_must_be_mapping(tmp_path, raw):
    with pytest.raises(
        ValidationError, match="mnase.fragments must be a mapping"
    ):
        validate_config(_make_config(tmp_path, mnase={"fragments": raw}))


def test_mnase_fragments_unknown_key_rejected(tmp_path):
    with pytest.raises(ValidationError, match="mnase.fragments: unknown key"):
        validate_config(
            _make_config(tmp_path, mnase={"fragments": {"unknown": [1, 2]}})
        )


def test_mnase_fragments_mono_overrides_only_fragments_mono(tmp_path):
    validated = validate_config(
        _make_config(
            tmp_path,
            mnase={
                "mono_range": [100, 180],
                "fragments": {"mono": [140, 200]},
            },
        )
    )
    assert validated["mnase"]["mono_range"] == [100, 180]
    assert validated["mnase"]["fragments"]["mono"] == [140, 200]


def test_mnase_fragments_partial_falls_back_to_defaults(tmp_path):
    validated = validate_config(
        _make_config(
            tmp_path,
            mnase={
                "mono_range": [150, 210],
                "fragments": {"sub": [10, 80]},
            },
        )
    )
    fragments = validated["mnase"]["fragments"]
    assert fragments["sub"] == [10, 80]
    assert fragments["mono"] == [150, 210]
    assert fragments["di"] == [300, 400]


@pytest.mark.parametrize(
    "raw",
    [
        [1, 139],
        (1, 139),
        ["1", "139"],
    ],
)
def test_mnase_fragments_sub_accepts_valid_inputs(tmp_path, raw):
    validated = validate_config(
        _make_config(tmp_path, mnase={"fragments": {"sub": raw}})
    )
    assert validated["mnase"]["fragments"]["sub"] == [1, 139]
    assert isinstance(validated["mnase"]["fragments"]["sub"], list)


@pytest.mark.parametrize(
    "raw,expected",
    [
        ([200, 100], "min must be < max"),
        ([0, 200], "values must be positive"),
        ([True, 200], "elements must be integers"),
        ([140], "exactly 2 elements"),
    ],
)
def test_mnase_fragments_sub_rejects_invalid_ranges(
    tmp_path, raw, expected
):
    with pytest.raises(
        ValidationError, match=f"mnase.fragments.sub.*{expected}"
    ):
        validate_config(
            _make_config(tmp_path, mnase={"fragments": {"sub": raw}})
        )


# ---------------------------------------------------------------------------
# dyad_range
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "raw",
    [
        [130, 200],
        (130, 200),
        ["130", "200"],
    ],
)
def test_mnase_dyad_range_accepts_valid_inputs(tmp_path, raw):
    validated = validate_config(_make_config(tmp_path, mnase={"dyad_range": raw}))
    assert validated["mnase"]["dyad_range"] == [130, 200]
    assert isinstance(validated["mnase"]["dyad_range"], list)


@pytest.mark.parametrize(
    "raw,expected",
    [
        ([130], "exactly 2 elements"),
        (["a", 200], "elements must be integers"),
        ([130.5, 200], "elements must be integers"),
        ([200, 100], "min must be < max"),
    ],
)
def test_mnase_dyad_range_rejects_invalid_ranges(tmp_path, raw, expected):
    with pytest.raises(
        ValidationError, match=f"mnase.dyad_range.*{expected}"
    ):
        validate_config(_make_config(tmp_path, mnase={"dyad_range": raw}))


# ---------------------------------------------------------------------------
# callers
# ---------------------------------------------------------------------------


def test_mnase_callers_empty_mapping_expands_to_defaults(tmp_path):
    validated = validate_config(_make_config(tmp_path, mnase={"callers": {}}))
    assert validated["mnase"]["callers"] == MNASE_DEFAULTS["callers"]


@pytest.mark.parametrize("raw", ["yes", ["x"], True])
def test_mnase_callers_must_be_mapping(tmp_path, raw):
    with pytest.raises(
        ValidationError, match="mnase.callers must be a mapping"
    ):
        validate_config(_make_config(tmp_path, mnase={"callers": raw}))


def test_mnase_callers_unknown_key_rejected(tmp_path):
    with pytest.raises(ValidationError, match="mnase.callers: unknown key"):
        validate_config(
            _make_config(tmp_path, mnase={"callers": {"unknown": False}})
        )


def test_mnase_callers_string_value_rejected(tmp_path):
    with pytest.raises(
        ValidationError, match="mnase.callers.danpos3 must be boolean"
    ):
        validate_config(
            _make_config(tmp_path, mnase={"callers": {"danpos3": "false"}})
        )


@pytest.mark.parametrize("caller", ["danpos3", "inps", "sem"])
def test_mnase_callers_true_rejected_as_not_implemented(tmp_path, caller):
    with pytest.raises(ValidationError, match="caller execution is not implemented"):
        validate_config(
            _make_config(tmp_path, mnase={"callers": {caller: True}})
        )


def test_mnase_callers_all_false_accepted(tmp_path):
    validated = validate_config(
        _make_config(
            tmp_path,
            mnase={
                "callers": {"danpos3": False, "inps": False, "sem": False},
            },
        )
    )
    assert validated["mnase"]["callers"] == {
        "danpos3": False,
        "inps": False,
        "sem": False,
    }


def test_mnase_callers_partial_false_fills_defaults(tmp_path):
    validated = validate_config(
        _make_config(tmp_path, mnase={"callers": {"danpos3": False}})
    )
    assert validated["mnase"]["callers"] == {
        "danpos3": False,
        "inps": False,
        "sem": False,
    }


# ---------------------------------------------------------------------------
# Explicit overrides
# ---------------------------------------------------------------------------


def test_mnase_explicit_values_override_defaults(tmp_path):
    validated = validate_config(
        _make_config(
            tmp_path,
            mnase={
                "mono_range": [100, 180],
                "fragments": {
                    "sub": [1, 120],
                    "mono": [130, 190],
                    "di": [280, 380],
                },
                "dyad_range": [120, 190],
                "callers": {"danpos3": False, "inps": False, "sem": False},
            },
        )
    )
    assert validated["mnase"]["mono_range"] == [100, 180]
    assert validated["mnase"]["fragments"] == {
        "sub": [1, 120],
        "mono": [130, 190],
        "di": [280, 380],
    }
    assert validated["mnase"]["dyad_range"] == [120, 190]
    assert validated["mnase"]["callers"] == {
        "danpos3": False,
        "inps": False,
        "sem": False,
    }
