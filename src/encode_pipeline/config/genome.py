"""Genome resource validation helpers."""

import os

from encode_pipeline.config import defaults

__all__ = [
    "validate_effective_genome_size",
    "validate_genome_resources",
    "validate_picard_reference_resources",
    "validate_tss_annotation_resources",
]


def validate_effective_genome_size(genome: str, value, error_cls=ValueError) -> None:
    """Validate MACS3 effective genome size shortcut or positive integer."""
    if isinstance(value, bool):
        valid = False
    elif isinstance(value, int):
        valid = value > 0
    elif isinstance(value, str):
        text = value.strip()
        valid = text in ("hs", "mm") or (text.isdigit() and int(text) > 0)
    else:
        valid = False

    if not valid:
        raise error_cls(
            f"genome_resources.{genome}: effective_genome_size must be "
            f"'hs', 'mm', or a positive integer, got {value!r}"
        )


def validate_genome_resources(resources: dict, error_cls=ValueError) -> dict:
    """Validate the genome_resources config block."""
    if not isinstance(resources, dict):
        raise error_cls(
            f"genome_resources must be a mapping, got {type(resources).__name__}"
        )

    for genome, entry in resources.items():
        if not isinstance(entry, dict):
            raise error_cls(
                f"genome_resources.{genome} must be a mapping, "
                f"got {type(entry).__name__}"
            )

        egs = entry.get("effective_genome_size")
        if egs is None or egs == "":
            raise error_cls(
                f"genome_resources.{genome}: effective_genome_size is required"
            )

        validate_effective_genome_size(genome, egs, error_cls=error_cls)

        for field in defaults.GENOME_RESOURCE_PATH_FIELDS:
            path = entry.get(field, "")
            if path and not os.path.isfile(path):
                raise error_cls(
                    f"genome_resources.{genome}.{field}: file not found: {path}"
                )

    return resources


def validate_picard_reference_resources(
    validated_config: dict,
    samples: list[dict],
    error_cls=ValueError,
) -> None:
    """Validate reference_fasta when qc.picard_metrics is enabled."""
    qc = validated_config.get("qc", {})
    if not qc.get("picard_metrics", False):
        return
    genome_resources = validated_config.get("genome_resources", {})
    treatment_genomes = {s["genome"] for s in samples if s["role"] == "treatment"}
    missing = []
    for genome in sorted(treatment_genomes):
        ref = genome_resources.get(genome, {}).get("reference_fasta", "")
        if not ref:
            missing.append(genome)
    if missing:
        raise error_cls(
            f"qc.picard_metrics is true but reference_fasta is missing "
            f"for genome(s): {', '.join(missing)}. "
            f"Set genome_resources[{missing[0]}].reference_fasta "
            f"or set qc.picard_metrics: false."
        )


def validate_tss_annotation_resources(
    validated_config: dict,
    samples: list[dict],
    error_cls=ValueError,
) -> None:
    """Validate GTF annotation when qc.tss_enrichment is enabled."""
    qc = validated_config.get("qc", {})
    if not qc.get("tss_enrichment", False):
        return
    genome_resources = validated_config.get("genome_resources", {})
    treatment_genomes = {s["genome"] for s in samples if s["role"] == "treatment"}
    missing = []
    for genome in sorted(treatment_genomes):
        gtf = genome_resources.get(genome, {}).get("gtf", "")
        if not gtf:
            missing.append(genome)
    if missing:
        raise error_cls(
            f"qc.tss_enrichment is true but gtf is missing "
            f"for genome(s): {', '.join(missing)}. "
            f"Set genome_resources[{missing[0]}].gtf "
            f"or set qc.tss_enrichment: false."
        )
