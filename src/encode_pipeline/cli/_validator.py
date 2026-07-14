"""Core validator CLI implementation (no structured logging wrapper)."""

import argparse
import os
import sys

from encode_pipeline.config import validator as validator_module
from encode_pipeline.config.yaml_loader import load_yaml


def main(argv=None):
    """Validate ChIP-seq/CUT&Tag/ATAC-seq pipeline config and sample sheet."""
    parser = argparse.ArgumentParser(
        description="Validate ChIP-seq/CUT&Tag/ATAC-seq pipeline config and "
        "sample sheet."
    )
    parser.add_argument(
        "--config", required=True, help="Path to config YAML (e.g. config/config.yaml)"
    )
    parser.add_argument(
        "--strict-inputs",
        action="store_true",
        default=False,
        help="Validate FASTQ and Bowtie2 index file existence",
    )
    args = parser.parse_args(argv)

    config_path = args.config
    if not os.path.isfile(config_path):
        print(f"ERROR: Config file not found: {config_path}", file=sys.stderr)
        sys.exit(1)

    config = load_yaml(config_path)

    try:
        validated = validator_module.validate_config(config)
        repro = validated.get("reproducibility", {})
        atac_idr_enabled = repro.get("enabled", False) and repro.get("idr", {}).get(
            "atac_narrow", False
        )
        cuttag_idr_enabled = repro.get("enabled", False) and repro.get("idr", {}).get(
            "cuttag_narrow", False
        )
        broad_chipseq_idr_enabled = repro.get("enabled", False) and repro.get(
            "idr", {}
        ).get("chipseq_broad_experimental", False)
        broad_cuttag_idr_enabled = repro.get("enabled", False) and repro.get(
            "idr", {}
        ).get("cuttag_broad_experimental", False)
        # Lazy import to break circular dependency with samples.load.
        from encode_pipeline.samples.load import load_and_validate_samples

        samples = load_and_validate_samples(
            validated["samples"],
            use_control=validated["use_control"],
            stage5_enabled=validated.get("stage5", False),
            strict_inputs=args.strict_inputs,
            reproducibility_idr_atac_narrow=atac_idr_enabled,
            reproducibility_idr_cuttag_narrow=cuttag_idr_enabled,
            reproducibility_idr_chipseq_broad=broad_chipseq_idr_enabled,
            reproducibility_idr_cuttag_broad=broad_cuttag_idr_enabled,
        )
        validator_module.validate_picard_reference_resources(validated, samples)
        validator_module.validate_tss_annotation_resources(validated, samples)
    except validator_module.ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    n_treatment = sum(1 for s in samples if s["role"] == "treatment")
    n_control = sum(1 for s in samples if s["role"] == "control")
    print(
        f"OK: {len(samples)} sample(s) validated "
        f"({n_treatment} treatment, {n_control} control)"
    )
    print(f"     use_control: {validated['use_control']}")
    if validated.get("genome_resources"):
        genomes = ", ".join(validated["genome_resources"])
        print(f"     genome resources: {genomes}")
