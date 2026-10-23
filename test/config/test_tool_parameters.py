"""Direct-API characterization tests for tool_parameters validation."""

import json
import os
import shlex
import subprocess
import sys

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


@pytest.mark.parametrize(
    "raw",
    [
        True,
        0,
        -0.1,
        "bad",
        1.5,
        "1.0000001",
        float("nan"),
        float("inf"),
        float("-inf"),
        "nan",
        "inf",
        "-inf",
    ],
)
def test_macs3_probability_fields_reject_out_of_domain_values(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match=r"tool_parameters\.macs3\.qvalue must be a probability in \(0, 1\]",
    ):
        _validate(tmp_path, {"macs3": {"qvalue": raw}})


@pytest.mark.parametrize("raw", [0.01, 0.5, 1, "0.05", "1"])
def test_macs3_probability_fields_accept_in_domain_values(tmp_path, raw):
    validated = _validate(tmp_path, {"macs3": {"qvalue": raw}})
    assert 0 < validated["macs3"]["qvalue"] <= 1


@pytest.mark.parametrize("raw", [0, 2, float("nan"), float("inf"), "nan"])
def test_macs3_broad_cutoff_rejects_out_of_domain_values(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match=(
            r"tool_parameters\.macs3\.broad_cutoff must be a probability "
            r"in \(0, 1\]"
        ),
    ):
        _validate(tmp_path, {"macs3": {"broad_cutoff": raw}})


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
    # --keep-dup all is synonymous with the canonical flag the workflow
    # renders itself, so it is normalized out of extra_args (SCI-2).
    assert validated["idr_macs3"] == {
        "pvalue": 0.1,
        "extra_args": "",
    }


@pytest.mark.parametrize("raw", [0, 2, float("nan"), float("-inf"), "inf"])
def test_idr_macs3_rejects_out_of_domain_pvalue(tmp_path, raw):
    with pytest.raises(
        ValidationError,
        match=(
            r"tool_parameters\.idr_macs3\.pvalue must be a probability "
            r"in \(0, 1\]"
        ),
    ):
        _validate(tmp_path, {"idr_macs3": {"pvalue": raw}})


# ---------------------------------------------------------------------------
# --keep-dup migration (SCI-2)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
@pytest.mark.parametrize(
    "extra,expected",
    [
        ("--keep-dup all", ""),
        ("--keep-dup=all", ""),
        ("--keep-dup all --nomodel", "--nomodel"),
        ("--nomodel --keep-dup=all -B", "--nomodel -B"),
        ("", ""),
    ],
)
def test_keep_dup_all_extra_args_is_normalized(tmp_path, block, extra, expected):
    validated = _validate(tmp_path, {block: {"extra_args": extra}})
    assert validated[block]["extra_args"] == expected


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
@pytest.mark.parametrize(
    "extra",
    [
        "--keep-dup 1",
        "--keep-dup=auto",
        "--keep-dup 5",
        "--keep-dup",
    ],
)
def test_conflicting_keep_dup_extra_args_is_rejected(tmp_path, block, extra):
    with pytest.raises(ValidationError, match="conflicts with the BAM-layer"):
        _validate(tmp_path, {block: {"extra_args": extra}})


def test_keep_dup_migration_does_not_touch_unrelated_tokens(tmp_path):
    validated = _validate(
        tmp_path,
        {"macs3": {"extra_args": "--keep-duplicates --keep-dup all -q 0.05"}},
    )
    assert validated["macs3"]["extra_args"] == "--keep-duplicates -q 0.05"


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
def test_extra_args_without_keep_dup_is_returned_byte_identical(tmp_path, block):
    # No managed flag: quoting of unrelated arguments must be preserved
    # exactly instead of being split and re-joined.
    extra = "--tempdir '/tmp/temp with spaces' --nomodel"
    validated = _validate(tmp_path, {block: {"extra_args": extra}})
    assert validated[block]["extra_args"] == extra


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
def test_keep_dup_migration_preserves_shell_argument_boundaries(tmp_path, block):
    extra = "--tempdir '/tmp/temp with spaces' --keep-dup all"
    normalized = _validate(tmp_path, {block: {"extra_args": extra}})[block][
        "extra_args"
    ]
    # The migrated string must still parse to the same shell arguments.
    assert shlex.split(normalized) == ["--tempdir", "/tmp/temp with spaces"]


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
def test_keep_dup_migration_is_stable_across_two_validations(tmp_path, block):
    extra = '--tempdir "/tmp/it\'s here" --keep-dup all'
    first = _validate(tmp_path, {block: {"extra_args": extra}})[block]["extra_args"]
    assert shlex.split(first) == ["--tempdir", "/tmp/it's here"]
    # A second validate_config pass over the normalized value is a no-op.
    second = _validate(tmp_path, {block: {"extra_args": first}})[block]["extra_args"]
    assert second == first


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
@pytest.mark.parametrize(
    "extra,expected",
    [
        ('--tempdir "$TMPDIR" --keep-dup all', '--tempdir "$TMPDIR"'),
        ("--tempdir '$TMPDIR' --keep-dup all", "--tempdir '$TMPDIR'"),
        ('--keep-dup all --tempdir "$TMPDIR"', '--tempdir "$TMPDIR"'),
        (
            '--keep-d=all --tempdir "$TMPDIR" --nomodel',
            '--tempdir "$TMPDIR" --nomodel',
        ),
    ],
)
def test_keep_dup_migration_preserves_variable_quoting_verbatim(
    tmp_path, block, extra, expected
):
    # Only the managed tokens' original text is removed; double-quoted
    # expansions must stay double-quoted so the shell still expands them.
    validated = _validate(tmp_path, {block: {"extra_args": extra}})
    assert validated[block]["extra_args"] == expected
    # Second pass over the normalized value is a stable no-op.
    again = _validate(tmp_path, {block: {"extra_args": expected}})[block]["extra_args"]
    assert again == expected


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
def test_extra_args_invalid_quoting_raises_validation_error(tmp_path, block):
    with pytest.raises(
        ValidationError,
        match=rf"tool_parameters\.{block}\.extra_args: invalid shell quoting",
    ):
        _validate(tmp_path, {block: {"extra_args": "--tempdir '/unterminated"}})


# argparse in the pinned MACS3 accepts any unambiguous prefix of a long
# option; --keep-dup is the only callpeak option starting with "--k", so
# abbreviations like --keep-d are the same flag and must be managed too.


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
@pytest.mark.parametrize(
    "extra",
    ["--keep-d all", "--keep-du=all", "--k all"],
)
def test_keep_dup_abbreviation_all_is_normalized(tmp_path, block, extra):
    validated = _validate(tmp_path, {block: {"extra_args": extra}})
    assert validated[block]["extra_args"] == ""


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
@pytest.mark.parametrize(
    "extra",
    [
        "--keep-d 1",
        "--keep-d=1",
        "--k 2",
        "--keep-du auto",
        "--keep-d",
        "--k=1 --nomodel",
    ],
)
def test_keep_dup_abbreviation_conflict_is_rejected(tmp_path, block, extra):
    with pytest.raises(ValidationError, match="conflicts with the BAM-layer"):
        _validate(tmp_path, {block: {"extra_args": extra}})


def _bash_argv(fragment: str, *, env: dict[str, str] | None = None) -> list[str]:
    """Return the actual argv array produced by /bin/bash for ``fragment``."""
    probe = shlex.join(
        [sys.executable, "-c", "import json, sys; print(json.dumps(sys.argv[1:]))"]
    )
    completed = subprocess.run(
        f"{probe} {fragment}",
        shell=True,
        executable="/bin/bash",
        capture_output=True,
        text=True,
        env=env,
        check=True,
    )
    return json.loads(completed.stdout)


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
@pytest.mark.parametrize(
    "extra,expected_bash_argv",
    [
        ("--keep-dup\\\n 1", ["--keep-dup", "1"]),
        ("--keep-\\\ndup 1", ["--keep-dup", "1"]),
        ("--keep-\\\ndup=auto", ["--keep-dup=auto"]),
        ("--keep-dup\\\n=1", ["--keep-dup=1"]),
        ("--keep-dup=\\\n1", ["--keep-dup=1"]),
        ("--keep-dup \\\n 1", ["--keep-dup", "1"]),
        ("--keep-d\\\n 1", ["--keep-d", "1"]),
        ("--k\\\n=2", ["--k=2"]),
        ('"--keep-\\\ndup" 1', ["--keep-dup", "1"]),
        ('"--keep-d\\\n"=auto', ["--keep-d=auto"]),
        ('--keep-dup "a\\\nuto"', ["--keep-dup", "auto"]),
    ],
)
def test_keep_dup_line_continuation_conflict_is_rejected_matching_bash(
    tmp_path, block, extra, expected_bash_argv
):
    assert _bash_argv(extra) == expected_bash_argv
    with pytest.raises(ValidationError, match="conflicts with the BAM-layer"):
        _validate(tmp_path, {block: {"extra_args": extra}})


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
@pytest.mark.parametrize(
    "extra,expected_normalized,expected_bash_argv",
    [
        ("--keep-dup a\\\nll", "", []),
        ("--keep-dup a\\\nll --nomodel", "--nomodel", ["--nomodel"]),
        ("--keep-\\\ndup all --nomodel", "--nomodel", ["--nomodel"]),
        ("--keep-dup\\\n all --nomodel", "--nomodel", ["--nomodel"]),
        ("--keep-dup \\\n all --nomodel", "--nomodel", ["--nomodel"]),
        ("--keep-dup=a\\\nll -B", "-B", ["-B"]),
        ("--keep-d\\\n=a\\\nll --nomodel", "--nomodel", ["--nomodel"]),
        ("--k\\\n a\\\nll -B", "-B", ["-B"]),
        ('"--keep-\\\ndup" "a\\\nll" --nomodel', "--nomodel", ["--nomodel"]),
    ],
)
def test_keep_dup_line_continuation_synonym_all_is_normalized_matching_bash(
    tmp_path, block, extra, expected_normalized, expected_bash_argv
):
    first = _validate(tmp_path, {block: {"extra_args": extra}})[block]["extra_args"]
    assert first == expected_normalized
    assert _bash_argv(first) == expected_bash_argv
    second = _validate(tmp_path, {block: {"extra_args": first}})[block]["extra_args"]
    assert second == first


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
def test_single_quoted_backslash_newline_remains_literal_matching_bash(tmp_path, block):
    extra = "'--keep-dup\\\n' 1"
    assert _bash_argv(extra) == ["--keep-dup\\\n", "1"]
    first = _validate(tmp_path, {block: {"extra_args": extra}})[block]["extra_args"]
    assert first == extra
    second = _validate(tmp_path, {block: {"extra_args": first}})[block]["extra_args"]
    assert second == first


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
@pytest.mark.parametrize(
    "extra,expected_normalized,expected_bash_argv",
    [
        (
            "--tempdir /tmp/review-temp\\  --keep-dup all",
            "--tempdir /tmp/review-temp\\ ",
            ["--tempdir", "/tmp/review-temp "],
        ),
        (
            "--tempdir /tmp/review\\ path --keep-dup a\\\nll",
            "--tempdir /tmp/review\\ path",
            ["--tempdir", "/tmp/review path"],
        ),
        (
            "--tempdir '/tmp/temp with spaces' --keep-\\\ndup all",
            "--tempdir '/tmp/temp with spaces'",
            ["--tempdir", "/tmp/temp with spaces"],
        ),
        (
            '--tempdir "/tmp/it\'s here" --keep-dup a\\\nll',
            '--tempdir "/tmp/it\'s here"',
            ["--tempdir", "/tmp/it's here"],
        ),
    ],
)
def test_keep_dup_migration_preserves_escaped_trailing_space_and_quotes_matching_bash(
    tmp_path, block, extra, expected_normalized, expected_bash_argv
):
    first = _validate(tmp_path, {block: {"extra_args": extra}})[block]["extra_args"]
    assert first == expected_normalized
    assert _bash_argv(first) == expected_bash_argv
    second = _validate(tmp_path, {block: {"extra_args": first}})[block]["extra_args"]
    assert second == first


@pytest.mark.parametrize("block", ["macs3", "idr_macs3"])
def test_keep_dup_migration_does_not_expand_env_vars_in_validator_but_bash_does(
    tmp_path, block, monkeypatch
):
    monkeypatch.setenv("TMPDIR", "/should/not/be/inlined/by/validator")
    double_quoted = '--tempdir "$TMPDIR/sub dir" --keep-dup a\\\nll'
    normalized_dq = _validate(tmp_path, {block: {"extra_args": double_quoted}})[block][
        "extra_args"
    ]
    assert normalized_dq == '--tempdir "$TMPDIR/sub dir"'
    assert (
        _validate(tmp_path, {block: {"extra_args": normalized_dq}})[block]["extra_args"]
        == normalized_dq
    )
    env = os.environ.copy()
    env["TMPDIR"] = "/runtime/expanded dir"
    assert _bash_argv(normalized_dq, env=env) == [
        "--tempdir",
        "/runtime/expanded dir/sub dir",
    ]

    single_quoted = "--tempdir '$TMPDIR/sub dir' --keep-dup a\\\nll"
    normalized_sq = _validate(tmp_path, {block: {"extra_args": single_quoted}})[block][
        "extra_args"
    ]
    assert normalized_sq == "--tempdir '$TMPDIR/sub dir'"
    assert _bash_argv(normalized_sq, env=env) == ["--tempdir", "$TMPDIR/sub dir"]
