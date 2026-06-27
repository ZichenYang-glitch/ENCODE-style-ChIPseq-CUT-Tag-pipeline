"""Direct-API characterization tests for genome resource validation."""

import pytest

from encode_pipeline.config.validate import (
    ValidationError,
    validate_config,
    validate_picard_reference_resources,
    validate_tss_annotation_resources,
)


SAMPLES_HEADER = (
    "sample\tfastq_1\tfastq_2\tlayout\tassay\ttarget\tpeak_mode\tgenome\tbowtie2_index\n"
)
SAMPLES_ROW = "S1\tR1.fq\tR2.fq\tPE\tchipseq\tT\tnarrow\ths\tidx\n"


def _make_config(tmp_path, *, genome_resources=None, qc=None):
    samples_path = tmp_path / "samples.tsv"
    samples_path.write_text(SAMPLES_HEADER + SAMPLES_ROW, encoding="utf-8")
    config = {"samples": str(samples_path), "use_control": False}
    if genome_resources is not None:
        config["genome_resources"] = genome_resources
    if qc is not None:
        config["qc"] = qc
    return config


def _treatment_sample(genome="hs"):
    return {
        "id": "S1",
        "role": "treatment",
        "genome": genome,
    }


def _control_sample(genome="mm"):
    return {
        "id": "C1",
        "role": "control",
        "genome": genome,
    }


# ---------------------------------------------------------------------------
# effective_genome_size
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("raw", ["hs", "mm", "2913022398", 2913022398])
def test_effective_genome_size_accepts_shortcuts_and_positive_ints(tmp_path, raw):
    resources = {"hs": {"effective_genome_size": raw}}
    validated = validate_config(_make_config(tmp_path, genome_resources=resources))
    assert validated["genome_resources"]["hs"]["effective_genome_size"] == raw


@pytest.mark.parametrize("raw", [True, False, 0, -1, "0", "-1", "hg38", "3.5"])
def test_effective_genome_size_rejects_invalid_values(tmp_path, raw):
    resources = {"hs": {"effective_genome_size": raw}}
    with pytest.raises(
        ValidationError,
        match=(
            "genome_resources.hs: effective_genome_size must be "
            "'hs', 'mm', or a positive integer"
        ),
    ):
        validate_config(_make_config(tmp_path, genome_resources=resources))


def test_genome_resources_must_be_mapping(tmp_path):
    with pytest.raises(ValidationError, match="genome_resources must be a mapping"):
        validate_config(_make_config(tmp_path, genome_resources=["hs"]))


def test_genome_resource_entry_must_be_mapping(tmp_path):
    with pytest.raises(
        ValidationError,
        match="genome_resources.hs must be a mapping",
    ):
        validate_config(_make_config(tmp_path, genome_resources={"hs": "bad"}))


def test_effective_genome_size_is_required(tmp_path):
    with pytest.raises(
        ValidationError,
        match="genome_resources.hs: effective_genome_size is required",
    ):
        validate_config(_make_config(tmp_path, genome_resources={"hs": {}}))


# ---------------------------------------------------------------------------
# Optional path fields
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("field", ["chrom_sizes", "blacklist", "gtf", "reference_fasta"])
def test_optional_genome_resource_path_fields_allow_empty_strings(tmp_path, field):
    resources = {"hs": {"effective_genome_size": "hs", field: ""}}
    validated = validate_config(_make_config(tmp_path, genome_resources=resources))
    assert validated["genome_resources"]["hs"][field] == ""


@pytest.mark.parametrize("field", ["chrom_sizes", "blacklist", "gtf", "reference_fasta"])
def test_optional_genome_resource_path_fields_require_existing_files(tmp_path, field):
    missing = tmp_path / f"missing-{field}.txt"
    resources = {"hs": {"effective_genome_size": "hs", field: str(missing)}}
    with pytest.raises(
        ValidationError,
        match=f"genome_resources.hs.{field}: file not found",
    ):
        validate_config(_make_config(tmp_path, genome_resources=resources))


@pytest.mark.parametrize("field", ["chrom_sizes", "blacklist", "gtf", "reference_fasta"])
def test_optional_genome_resource_path_fields_accept_existing_files(tmp_path, field):
    resource_path = tmp_path / f"{field}.txt"
    resource_path.write_text("placeholder\n", encoding="utf-8")
    resources = {"hs": {"effective_genome_size": "hs", field: str(resource_path)}}
    validated = validate_config(_make_config(tmp_path, genome_resources=resources))
    assert validated["genome_resources"]["hs"][field] == str(resource_path)


# ---------------------------------------------------------------------------
# Picard / TSS resource gates
# ---------------------------------------------------------------------------


def test_picard_reference_gate_is_disabled_when_qc_false():
    validate_picard_reference_resources(
        {"qc": {"picard_metrics": False}, "genome_resources": {}},
        [_treatment_sample("hs")],
    )


def test_picard_reference_gate_requires_treatment_genome_reference():
    with pytest.raises(
        ValidationError,
        match="qc.picard_metrics is true but reference_fasta is missing for genome",
    ):
        validate_picard_reference_resources(
            {
                "qc": {"picard_metrics": True},
                "genome_resources": {"hs": {"reference_fasta": ""}},
            },
            [_treatment_sample("hs")],
        )


def test_picard_reference_gate_ignores_control_only_genomes():
    validate_picard_reference_resources(
        {
            "qc": {"picard_metrics": True},
            "genome_resources": {"hs": {"reference_fasta": "ref.fa"}},
        },
        [_treatment_sample("hs"), _control_sample("mm")],
    )


def test_tss_annotation_gate_is_disabled_when_qc_false():
    validate_tss_annotation_resources(
        {"qc": {"tss_enrichment": False}, "genome_resources": {}},
        [_treatment_sample("hs")],
    )


def test_tss_annotation_gate_requires_treatment_genome_gtf():
    with pytest.raises(
        ValidationError,
        match="qc.tss_enrichment is true but gtf is missing for genome",
    ):
        validate_tss_annotation_resources(
            {
                "qc": {"tss_enrichment": True},
                "genome_resources": {"hs": {"gtf": ""}},
            },
            [_treatment_sample("hs")],
        )


def test_tss_annotation_gate_ignores_control_only_genomes():
    validate_tss_annotation_resources(
        {
            "qc": {"tss_enrichment": True},
            "genome_resources": {"hs": {"gtf": "annotation.gtf"}},
        },
        [_treatment_sample("hs"), _control_sample("mm")],
    )
