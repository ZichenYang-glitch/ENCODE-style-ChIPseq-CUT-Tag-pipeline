"""Direct-API characterization tests for tool_parameters validation."""

import pytest

from encode_pipeline.config.validate import ValidationError, validate_config


SAMPLES_HEADER = "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
SAMPLES_ROW = "S1\tR1.fq\tR2.fq\tPE\tchipseq\tT\tnarrow\ths\tidx\n"


def _make_config(tmp_path, *, tool_parameters=None):
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(SAMPLES_HEADER + SAMPLES_ROW, encoding="utf-8")
    config = {"samples": str(samples_path), "use_control": False}
    if tool_parameters is not None:
        config["tool_parameters"] = tool_parameters
    return config


def _validate(tmp_path, tool_parameters):
    config = _make_config(tmp_path, tool_parameters=tool_parameters)
    return validate_config(config)["tool_parameters"]


# ---------------------------------------------------------------------------
# Structure and unknown keys
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("raw", [True, False, ["macs3"]])
def test_tool_parameters_must_be_mapping(tmp_path, raw):
    with pytest.raises(ValidationError, match="tool_parameters must be a mapping"):
        validate_config(_make_config(tmp_path, tool_parameters=raw))


def test_unknown_tool_block_is_rejected(tmp_path):
    with pytest.raises(
        ValidationError,
        match="tool_parameters: unknown tool block 'unknown_tool'",
    ):
        _validate(tmp_path, {"unknown_tool": {}})


def test_tool_block_must_be_mapping(tmp_path):
    with pytest.raises(
        ValidationError,
        match="tool_parameters.fastqc must be a mapping",
    ):
        _validate(tmp_path, {"fastqc": "bad"})


def test_unknown_tool_key_is_rejected(tmp_path):
    with pytest.raises(
        ValidationError,
        match="tool_parameters.macs3: unknown key 'bad_key'",
    ):
        _validate(tmp_path, {"macs3": {"bad_key": 1}})


# ---------------------------------------------------------------------------
# Per-tool normalization
# ---------------------------------------------------------------------------


def test_missing_tool_parameters_expands_known_tool_defaults(tmp_path):
    validated = validate_config(_make_config(tmp_path))["tool_parameters"]
    assert validated["fastqc"]["extra_args"] == ""
    assert validated["bowtie2"]["mode"] == ""
    assert validated["bamcoverage"]["normalize_using"] == "CPM"
    assert validated["macs3"]["qvalue"] == 0.01
    assert validated["idr_macs3"]["pvalue"] == 0.1


def test_fastqc_extra_args_must_be_string(tmp_path):
    with pytest.raises(
        ValidationError,
        match="tool_parameters.fastqc.extra_args must be a string",
    ):
        _validate(tmp_path, {"fastqc": {"extra_args": ["--quiet"]}})


def test_trim_galore_integer_fields_normalize_numeric_strings_and_empty(tmp_path):
    validated = _validate(
        tmp_path,
        {
            "trim_galore": {
                "quality": "20",
                "length": "",
                "stringency": 3,
                "extra_args": "--paired",
            }
        },
    )
    assert validated["trim_galore"] == {
        "quality": 20,
        "length": "",
        "stringency": 3,
        "extra_args": "--paired",
    }


@pytest.mark.parametrize("raw", [True, "3.5", -1])
def test_trim_galore_integer_fields_reject_invalid_values(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match=(
            "tool_parameters.trim_galore.quality must be an integer|"
            "tool_parameters.trim_galore.quality must be non-negative"
        ),
    ):
        _validate(tmp_path, {"trim_galore": {"quality": raw}})


def test_bowtie2_mode_and_bool_flags_are_normalized(tmp_path):
    validated = _validate(
        tmp_path,
        {
            "bowtie2": {
                "mode": "very-sensitive",
                "dovetail": "true",
                "no_mixed": False,
                "no_discordant": "false",
                "extra_args": "--seed 1",
            }
        },
    )
    assert validated["bowtie2"] == {
        "mode": "very-sensitive",
        "dovetail": True,
        "no_mixed": False,
        "no_discordant": False,
        "extra_args": "--seed 1",
    }


def test_bowtie2_rejects_invalid_mode(tmp_path):
    with pytest.raises(
        ValidationError,
        match="tool_parameters.bowtie2.mode must be one of",
    ):
        _validate(tmp_path, {"bowtie2": {"mode": "ultra"}})


def test_bowtie2_rejects_invalid_bool_flag(tmp_path):
    with pytest.raises(
        ValidationError,
        match="tool_parameters.bowtie2.dovetail must be true or false",
    ):
        _validate(tmp_path, {"bowtie2": {"dovetail": "yes"}})


@pytest.mark.parametrize(
    "raw,expected",
    [
        ("", ""),
        (2308, 2308),
        ("2308", 2308),
        ("0x904", 2308),
    ],
)
def test_samtools_filter_flags_normalization(tmp_path, raw, expected):
    validated = _validate(tmp_path, {"samtools_filter": {"filter_flags": raw}})
    assert validated["samtools_filter"]["filter_flags"] == expected


@pytest.mark.parametrize("raw", [True, 0, "0x0", "0xZZ", "abc"])
def test_samtools_filter_flags_reject_invalid_values(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match=(
            "tool_parameters.samtools_filter.filter_flags must be|"
            "tool_parameters.samtools_filter.filter_flags: invalid hex value"
        ),
    ):
        _validate(tmp_path, {"samtools_filter": {"filter_flags": raw}})


def test_picard_markduplicates_integer_field_normalizes_string(tmp_path):
    validated = _validate(
        tmp_path,
        {
            "picard_markduplicates": {
                "optical_duplicate_pixel_distance": "2500",
                "extra_args": "--VALIDATION_STRINGENCY SILENT",
            }
        },
    )
    assert (
        validated["picard_markduplicates"]["optical_duplicate_pixel_distance"] == 2500
    )
    assert (
        validated["picard_markduplicates"]["extra_args"]
        == "--VALIDATION_STRINGENCY SILENT"
    )


def test_bamcoverage_normalize_using_and_smooth_length(tmp_path):
    validated = _validate(
        tmp_path,
        {
            "bamcoverage": {
                "normalize_using": "None",
                "smooth_length": "50",
                "extra_args": "--extendReads",
            }
        },
    )
    assert validated["bamcoverage"] == {
        "normalize_using": "None",
        "smooth_length": 50,
        "extra_args": "--extendReads",
    }


def test_bamcoverage_rejects_unsupported_normalize_using(tmp_path):
    with pytest.raises(
        ValidationError,
        match="RPGC is not supported in this slice",
    ):
        _validate(tmp_path, {"bamcoverage": {"normalize_using": "RPGC"}})


def test_bamcoverage_smooth_length_must_be_positive(tmp_path):
    with pytest.raises(
        ValidationError,
        match="tool_parameters.bamcoverage.smooth_length must be non-negative",
    ):
        _validate(tmp_path, {"bamcoverage": {"smooth_length": 0}})


def test_macs3_positive_float_fields_normalize_strings(tmp_path):
    validated = _validate(
        tmp_path,
        {
            "macs3": {
                "qvalue": "0.05",
                "broad_cutoff": 0.2,
                "extra_args": "--nomodel",
            }
        },
    )
    assert validated["macs3"] == {
        "qvalue": 0.05,
        "broad_cutoff": 0.2,
        "extra_args": "--nomodel",
    }


@pytest.mark.parametrize("raw", [True, 0, -0.1, "bad"])
def test_macs3_positive_float_fields_reject_invalid_values(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match=(
            "tool_parameters.macs3.qvalue must be a number > 0|"
            "tool_parameters.macs3.qvalue must be positive"
        ),
    ):
        _validate(tmp_path, {"macs3": {"qvalue": raw}})


def test_multiqc_title_and_extra_args_must_be_strings(tmp_path):
    validated = _validate(
        tmp_path,
        {"multiqc": {"title": "ENCODE run", "extra_args": "--force"}},
    )
    assert validated["multiqc"] == {
        "title": "ENCODE run",
        "extra_args": "--force",
    }

    with pytest.raises(
        ValidationError,
        match="tool_parameters.multiqc.title must be a string",
    ):
        _validate(tmp_path, {"multiqc": {"title": 123}})


def test_idr_macs3_pvalue_normalizes_positive_float_string(tmp_path):
    validated = _validate(
        tmp_path,
        {"idr_macs3": {"pvalue": "0.1", "extra_args": "--keep-dup all"}},
    )
    assert validated["idr_macs3"] == {
        "pvalue": 0.1,
        "extra_args": "--keep-dup all",
    }


def test_idr_macs3_rejects_non_positive_pvalue(tmp_path):
    with pytest.raises(
        ValidationError,
        match="tool_parameters.idr_macs3.pvalue must be positive",
    ):
        _validate(tmp_path, {"idr_macs3": {"pvalue": 0}})
