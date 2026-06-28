"""Unit tests for genome resource validation helpers."""

import pytest

from encode_pipeline.config.genome import (
    validate_effective_genome_size,
    validate_genome_resources,
    validate_picard_reference_resources,
    validate_tss_annotation_resources,
)


class CustomError(Exception):
    """Sentinel exception for testing error injection."""


def test_validate_effective_genome_size_accepts_shortcuts_and_ints():
    validate_effective_genome_size("hs", "hs")
    validate_effective_genome_size("mm", "mm")
    validate_effective_genome_size("custom", 12345)


def test_validate_effective_genome_size_uses_custom_error_class():
    with pytest.raises(CustomError, match="effective_genome_size"):
        validate_effective_genome_size("hs", "hg38", error_cls=CustomError)


def test_validate_genome_resources_uses_custom_error_class():
    with pytest.raises(CustomError, match="genome_resources must be a mapping"):
        validate_genome_resources([], error_cls=CustomError)


def test_validate_picard_reference_resources_uses_custom_error_class():
    with pytest.raises(CustomError, match="reference_fasta is missing"):
        validate_picard_reference_resources(
            {"qc": {"picard_metrics": True}, "genome_resources": {}},
            [{"role": "treatment", "genome": "hs"}],
            error_cls=CustomError,
        )


def test_validate_tss_annotation_resources_uses_custom_error_class():
    with pytest.raises(CustomError, match="gtf is missing"):
        validate_tss_annotation_resources(
            {"qc": {"tss_enrichment": True}, "genome_resources": {}},
            [{"role": "treatment", "genome": "hs"}],
            error_cls=CustomError,
        )
