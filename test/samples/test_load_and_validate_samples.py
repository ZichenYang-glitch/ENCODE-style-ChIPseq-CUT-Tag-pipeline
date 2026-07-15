"""Direct behavior tests for sample-sheet loading and row validation."""

import pytest

from encode_pipeline.config.validate import ValidationError
from encode_pipeline.samples.load import load_and_validate_samples


REQUIRED_COLUMNS = [
    "sample",
    "fastq_1",
    "layout",
    "assay",
    "target",
    "peak_mode",
    "genome",
    "bowtie2_index",
]

# Default column order used for most test sheets. fastq_2 is optional but
# included by default because the default row is PE.
DEFAULT_COLUMNS = REQUIRED_COLUMNS + ["fastq_2"]

COLUMN_DEFAULTS = {
    "sample": "S1",
    "fastq_1": "R1.fq",
    "fastq_2": "R2.fq",
    "layout": "PE",
    "assay": "chipseq",
    "target": "T",
    "peak_mode": "narrow",
    "genome": "hs",
    "bowtie2_index": "idx",
}


def _write_tsv(path, content):
    path.write_text(content, encoding="utf-8")


def _make_header(columns):
    return "\t".join(columns)


def _make_row(overrides=None, columns=None):
    values = dict(COLUMN_DEFAULTS)
    values.update(overrides or {})
    columns = columns or DEFAULT_COLUMNS
    return "\t".join(str(values[col]) for col in columns)


def _make_tsv(path, overrides=None, extra_columns=None, extra_values=None):
    columns = list(DEFAULT_COLUMNS)
    if extra_columns:
        columns.extend(extra_columns)
    row_values = dict(COLUMN_DEFAULTS)
    row_values.update(overrides or {})
    if extra_values:
        row_values.update(extra_values)
    row = "\t".join(str(row_values[col]) for col in columns)
    _write_tsv(path, _make_header(columns) + "\n" + row + "\n")


# ---------------------------------------------------------------------------
# Happy path and normalization
# ---------------------------------------------------------------------------


def test_happy_path_pe_returns_normalized_dict(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(samples_path)
    samples = load_and_validate_samples(str(samples_path))
    assert len(samples) == 1
    assert samples[0] == {
        "id": "S1",
        "fq1": "R1.fq",
        "fq2": "R2.fq",
        "layout": "PE",
        "assay": "chipseq",
        "target": "T",
        "peak_mode": "narrow",
        "genome": "hs",
        "bt2_idx": "idx",
        "control_bam": "",
        "role": "treatment",
        "control_sample": "",
        "experiment": "S1",
        "condition": "T",
        "replicate": 1,
        "biological_replicate": 1,
        "technical_replicate": 1,
    }


def test_happy_path_se_accepts_missing_fastq_2(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        overrides={"fastq_2": "", "layout": "SE"},
    )
    samples = load_and_validate_samples(str(samples_path))
    assert len(samples) == 1
    assert samples[0]["layout"] == "SE"
    assert samples[0]["fq2"] == ""


def test_case_and_defaults_are_normalized(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _write_tsv(
        samples_path,
        "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
        "S1\tR1.fq\tR2.fq\tpe\tChIPseq\tT\tNARROW\ths\tidx\n",
    )
    samples = load_and_validate_samples(str(samples_path))
    s = samples[0]
    assert s["layout"] == "PE"
    assert s["assay"] == "chipseq"
    assert s["peak_mode"] == "narrow"
    assert s["role"] == "treatment"
    assert s["experiment"] == "S1"
    assert s["condition"] == "T"
    assert s["replicate"] == 1
    assert s["biological_replicate"] == 1
    assert s["technical_replicate"] == 1


def test_optional_replicate_columns_default_sensibly(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        overrides={"replicate": "3"},
        extra_columns=["replicate"],
    )
    samples = load_and_validate_samples(str(samples_path))
    s = samples[0]
    assert s["replicate"] == 3
    assert s["biological_replicate"] == 3
    assert s["technical_replicate"] == 1


# ---------------------------------------------------------------------------
# Required columns and basic structure
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("missing_col", REQUIRED_COLUMNS)
def test_missing_required_column_rejected(tmp_path, missing_col):
    samples_path = tmp_path / "samples.tsv"
    columns = [c for c in DEFAULT_COLUMNS if c != missing_col]
    header = _make_header(columns)
    row_values = {c: COLUMN_DEFAULTS[c] for c in columns}
    row = "\t".join(str(row_values[c]) for c in columns)
    _write_tsv(samples_path, header + "\n" + row + "\n")
    with pytest.raises(
        ValidationError, match=f"missing required column: {missing_col!r}"
    ):
        load_and_validate_samples(str(samples_path))


def test_sample_sheet_not_found(tmp_path):
    missing_path = tmp_path / "does_not_exist.tsv"
    with pytest.raises(ValidationError, match="Sample sheet not found"):
        load_and_validate_samples(str(missing_path))


# ---------------------------------------------------------------------------
# Per-row field validation
# ---------------------------------------------------------------------------


def test_empty_sample_id_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(samples_path, overrides={"sample": ""})
    with pytest.raises(ValidationError, match="empty 'sample'"):
        load_and_validate_samples(str(samples_path))


def test_invalid_sample_id_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(samples_path, overrides={"sample": "S1 Bad"})
    with pytest.raises(ValidationError, match="invalid characters"):
        load_and_validate_samples(str(samples_path))


def test_duplicate_sample_id_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    header = _make_header(DEFAULT_COLUMNS)
    row = _make_row()
    _write_tsv(samples_path, header + "\n" + row + "\n" + row + "\n")
    with pytest.raises(
        ValidationError, match="Duplicate sample ID in sample sheet: 'S1'"
    ):
        load_and_validate_samples(str(samples_path))


def test_pe_layout_requires_fastq_2(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(samples_path, overrides={"fastq_2": ""})
    with pytest.raises(ValidationError, match="PE layout requires 'fastq_2'"):
        load_and_validate_samples(str(samples_path))


@pytest.mark.parametrize("field", ["target", "genome", "bowtie2_index"])
def test_empty_required_string_field_rejected(tmp_path, field):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(samples_path, overrides={field: ""})
    with pytest.raises(ValidationError, match=f"empty '{field}'"):
        load_and_validate_samples(str(samples_path))


def test_invalid_layout_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(samples_path, overrides={"layout": "SEPE"})
    with pytest.raises(ValidationError, match="layout must be PE or SE"):
        load_and_validate_samples(str(samples_path))


def test_invalid_assay_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(samples_path, overrides={"assay": "rnaseq"})
    with pytest.raises(
        ValidationError, match="assay must be chipseq, cuttag, atac, or mnase"
    ):
        load_and_validate_samples(str(samples_path))


def test_invalid_peak_mode_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(samples_path, overrides={"peak_mode": "mixed"})
    with pytest.raises(
        ValidationError, match="peak_mode must be narrow, broad, or nucleosome"
    ):
        load_and_validate_samples(str(samples_path))


def test_invalid_role_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        extra_columns=["role"],
        extra_values={"role": "donor"},
    )
    with pytest.raises(ValidationError, match="role must be treatment or control"):
        load_and_validate_samples(str(samples_path))


# ---------------------------------------------------------------------------
# Assay / peak_mode coupling
# ---------------------------------------------------------------------------


def test_mnase_requires_pe_and_nucleosome_valid(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        overrides={"assay": "mnase", "peak_mode": "nucleosome"},
    )
    samples = load_and_validate_samples(str(samples_path))
    assert samples[0]["assay"] == "mnase"
    assert samples[0]["peak_mode"] == "nucleosome"


def test_mnase_se_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        overrides={
            "assay": "mnase",
            "peak_mode": "nucleosome",
            "layout": "SE",
            "fastq_2": "",
        },
    )
    with pytest.raises(ValidationError, match="assay=mnase requires paired-end layout"):
        load_and_validate_samples(str(samples_path))


def test_mnase_narrow_peak_mode_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        overrides={"assay": "mnase", "peak_mode": "narrow"},
    )
    with pytest.raises(
        ValidationError, match="assay=mnase requires peak_mode=nucleosome"
    ):
        load_and_validate_samples(str(samples_path))


def test_atac_broad_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        overrides={
            "assay": "atac",
            "peak_mode": "broad",
            "fastq_2": "",
            "layout": "SE",
        },
    )
    with pytest.raises(
        ValidationError, match="assay=atac currently supports peak_mode=narrow only"
    ):
        load_and_validate_samples(str(samples_path))


def test_chipseq_nucleosome_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        overrides={"peak_mode": "nucleosome"},
    )
    with pytest.raises(
        ValidationError, match="peak_mode=nucleosome is only allowed for assay=mnase"
    ):
        load_and_validate_samples(str(samples_path))


# ---------------------------------------------------------------------------
# Identifier sanitization
# ---------------------------------------------------------------------------


def test_experiment_and_condition_are_sanitized(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        extra_columns=["experiment", "condition"],
        extra_values={"experiment": "EXP 1", "condition": "Pol II"},
    )
    samples = load_and_validate_samples(str(samples_path))
    assert samples[0]["experiment"] == "EXP_1"
    assert samples[0]["condition"] == "Pol_II"


# ---------------------------------------------------------------------------
# Control gating
# ---------------------------------------------------------------------------


def test_use_control_false_ignores_invalid_control_sample(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        extra_columns=["control_sample"],
        extra_values={"control_sample": "MissingCtrl"},
    )
    samples = load_and_validate_samples(str(samples_path), use_control=False)
    assert samples[0]["control_sample"] == "MissingCtrl"


def test_use_control_false_ignores_missing_control_bam(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(
        samples_path,
        extra_columns=["control_bam"],
        extra_values={"control_bam": "/does/not/exist.bam"},
    )
    samples = load_and_validate_samples(str(samples_path), use_control=False)
    assert samples[0]["control_bam"] == "/does/not/exist.bam"


def test_use_control_true_missing_control_sample_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\trole\tcontrol_sample\n"
    )
    row = "S1\tR1.fq\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tCtrl1\n"
    _write_tsv(samples_path, header + row)
    with pytest.raises(ValidationError, match="not found in sample sheet"):
        load_and_validate_samples(str(samples_path), use_control=True)


def test_use_control_true_control_sample_wrong_role_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\trole\tcontrol_sample\n"
    )
    rows = (
        "S1\tR1.fq\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tCtrl1\n"
        "Ctrl1\tR1.fq\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\t\n"
    )
    _write_tsv(samples_path, header + rows)
    with pytest.raises(ValidationError, match="expected role=control"):
        load_and_validate_samples(str(samples_path), use_control=True)


def test_use_control_true_self_reference_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\trole\tcontrol_sample\n"
    )
    row = "S1\tR1.fq\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tS1\n"
    _write_tsv(samples_path, header + row)
    with pytest.raises(ValidationError, match="control_sample cannot reference itself"):
        load_and_validate_samples(str(samples_path), use_control=True)


def test_use_control_true_mutually_exclusive_controls_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\trole\tcontrol_sample\tcontrol_bam\n"
    )
    rows = (
        "S1\tR1.fq\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tCtrl1\t/tmp/a.bam\n"
        "Ctrl1\tR1.fq\tSE\tchipseq\tT\tnarrow\ths\tidx\tcontrol\t\t\n"
    )
    _write_tsv(samples_path, header + rows)
    with pytest.raises(
        ValidationError, match="cannot set both control_sample and control_bam"
    ):
        load_and_validate_samples(str(samples_path), use_control=True)


def test_use_control_true_missing_control_bam_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\trole\tcontrol_bam\n"
    )
    row = "S1\tR1.fq\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\t/does/not/exist.bam\n"
    _write_tsv(samples_path, header + row)
    with pytest.raises(ValidationError, match="control_bam file not found"):
        load_and_validate_samples(str(samples_path), use_control=True)


def test_use_control_true_valid_control_sample_passes(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\trole\tcontrol_sample\n"
    )
    rows = (
        "S1\tR1.fq\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\tCtrl1\n"
        "Ctrl1\tR1.fq\tSE\tchipseq\tIgG\tnarrow\ths\tidx\tcontrol\t\n"
    )
    _write_tsv(samples_path, header + rows)
    samples = load_and_validate_samples(str(samples_path), use_control=True)
    assert len(samples) == 2
    s1 = next(s for s in samples if s["id"] == "S1")
    ctrl = next(s for s in samples if s["id"] == "Ctrl1")
    assert s1["control_sample"] == "Ctrl1"
    assert ctrl["role"] == "control"
    assert ctrl["control_sample"] == ""


def test_use_control_true_valid_control_bam_passes(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    bam_path = tmp_path / "control.bam"
    bam_path.write_text("", encoding="utf-8")
    header = (
        "sample\tfastq_1\tlayout\tassay\ttarget\tpeak_mode\tgenome\t"
        "bowtie2_index\trole\tcontrol_bam\n"
    )
    row = f"S1\tR1.fq\tSE\tchipseq\tT\tnarrow\ths\tidx\ttreatment\t{bam_path}\n"
    _write_tsv(samples_path, header + row)
    samples = load_and_validate_samples(str(samples_path), use_control=True)
    assert len(samples) == 1
    assert samples[0]["control_bam"] == str(bam_path)


# ---------------------------------------------------------------------------
# Strict input files and index formats
# ---------------------------------------------------------------------------


def test_strict_inputs_missing_fastq_1_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    fq1 = tmp_path / "missing_R1.fq"
    _make_tsv(samples_path, overrides={"fastq_1": str(fq1)})
    with pytest.raises(ValidationError, match="fastq_1 file not found"):
        load_and_validate_samples(str(samples_path), strict_inputs=True)


def test_strict_inputs_pe_missing_fastq_2_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    fq1 = tmp_path / "R1.fq"
    fq1.write_text("")
    _make_tsv(samples_path, overrides={"fastq_1": str(fq1), "fastq_2": ""})
    with pytest.raises(ValidationError, match="PE layout requires 'fastq_2'"):
        load_and_validate_samples(str(samples_path), strict_inputs=True)


def test_strict_inputs_complete_bt2_set_accepted(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    fq1 = tmp_path / "R1.fq"
    fq2 = tmp_path / "R2.fq"
    fq1.write_text("")
    fq2.write_text("")
    prefix = tmp_path / "idx"
    for ext in (".1.bt2", ".2.bt2", ".3.bt2", ".4.bt2", ".rev.1.bt2", ".rev.2.bt2"):
        (tmp_path / (prefix.name + ext)).write_text("")
    _make_tsv(
        samples_path,
        overrides={
            "fastq_1": str(fq1),
            "fastq_2": str(fq2),
            "bowtie2_index": str(prefix),
        },
    )
    samples = load_and_validate_samples(str(samples_path), strict_inputs=True)
    assert len(samples) == 1
    assert samples[0]["bt2_idx"] == str(prefix)


def test_strict_inputs_missing_bt2_set_rejected(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    fq1 = tmp_path / "R1.fq"
    fq1.write_text("")
    prefix = tmp_path / "missing_idx"
    _make_tsv(
        samples_path,
        overrides={
            "fastq_1": str(fq1),
            "fastq_2": "",
            "layout": "SE",
            "bowtie2_index": str(prefix),
        },
    )
    with pytest.raises(ValidationError, match="Bowtie2 index not found"):
        load_and_validate_samples(str(samples_path), strict_inputs=True)


def test_strict_inputs_complete_bt2l_set_accepted(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    fq1 = tmp_path / "R1.fq"
    fq2 = tmp_path / "R2.fq"
    fq1.write_text("")
    fq2.write_text("")
    prefix = tmp_path / "idx"
    for ext in (
        ".1.bt2l",
        ".2.bt2l",
        ".3.bt2l",
        ".4.bt2l",
        ".rev.1.bt2l",
        ".rev.2.bt2l",
    ):
        (tmp_path / (prefix.name + ext)).write_text("")
    _make_tsv(
        samples_path,
        overrides={
            "fastq_1": str(fq1),
            "fastq_2": str(fq2),
            "bowtie2_index": str(prefix),
        },
    )

    samples = load_and_validate_samples(str(samples_path), strict_inputs=True)

    assert samples[0]["bt2_idx"] == str(prefix)


def test_non_strict_inputs_accept_placeholder_paths(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    _make_tsv(samples_path)

    samples = load_and_validate_samples(str(samples_path), strict_inputs=False)

    assert samples[0]["fq1"] == "R1.fq"
    assert samples[0]["bt2_idx"] == "idx"


def test_strict_inputs_accept_gzipped_fastq_names(tmp_path):
    samples_path = tmp_path / "samples.tsv"
    fq1 = tmp_path / "R1.fq.gz"
    fq2 = tmp_path / "R2.fq.gz"
    fq1.write_text("")
    fq2.write_text("")
    prefix = tmp_path / "idx"
    for ext in (
        ".1.bt2",
        ".2.bt2",
        ".3.bt2",
        ".4.bt2",
        ".rev.1.bt2",
        ".rev.2.bt2",
    ):
        (tmp_path / (prefix.name + ext)).write_text("")
    _make_tsv(
        samples_path,
        overrides={
            "fastq_1": str(fq1),
            "fastq_2": str(fq2),
            "bowtie2_index": str(prefix),
        },
    )

    samples = load_and_validate_samples(str(samples_path), strict_inputs=True)

    assert samples[0]["fq1"] == str(fq1)
    assert samples[0]["fq2"] == str(fq2)
