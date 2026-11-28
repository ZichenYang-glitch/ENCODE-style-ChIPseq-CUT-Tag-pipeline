"""Core validator CLI implementation (no structured logging wrapper)."""

import argparse
import os
import sys

from encode_pipeline.config import validator as validator_module
from encode_pipeline.config.reproducibility import idr_sample_flags_from_validated
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
        validator_module.validate_execution_support(validated)
        idr_flags = idr_sample_flags_from_validated(validated)
        # Lazy import to break circular dependency with samples.load.
        from encode_pipeline.samples.load import load_and_validate_samples
        from encode_pipeline.samples.mixed_replicates import (
            mixed_replicate_notice,
            resolve_mixed_replicate_routes,
        )

        samples = load_and_validate_samples(
            validated["samples"],
            use_control=validated["use_control"],
            chipseq_idr_enabled=idr_flags["chipseq_narrow"],
            strict_inputs=args.strict_inputs,
            mixed_replicate_config=validated,
            reproducibility_idr_atac_narrow=idr_flags["atac_narrow"],
            reproducibility_idr_cuttag_narrow=idr_flags["cuttag_narrow"],
            reproducibility_idr_chipseq_broad=idr_flags["chipseq_broad"],
            reproducibility_idr_cuttag_broad=idr_flags["cuttag_broad"],
        )
        validator_module.validate_bowtie2_insert_bounds(validated, samples)
        validator_module.validate_picard_reference_resources(validated, samples)
        validator_module.validate_tss_annotation_resources(validated, samples)
        from encode_pipeline.samples.admission import (
            validate_reference_control_admission,
        )

        validate_reference_control_admission(
            samples,
            validated,
            reference_mode=validated["reference_mode"],
            external_control_policy=validated["external_control_policy"],
            error_cls=validator_module.ValidationError,
        )
    except validator_module.ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    routes = resolve_mixed_replicate_routes(samples, validated, idr_flags)
    if routes:
        print(f"WARNING: {mixed_replicate_notice(routes)}", file=sys.stderr)

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
