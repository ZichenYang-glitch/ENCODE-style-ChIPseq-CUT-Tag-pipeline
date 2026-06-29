#!/usr/bin/env python3
"""Sample sheet and config validation for the ChIP-seq / CUT&Tag / ATAC-seq / MNase-seq pipeline.

Importable by workflow/Snakefile and runnable as a standalone CLI.

Usage:
    python3 scripts/validate_samples.py --config config/config.yaml
"""

import csv
import os
import sys

from encode_pipeline.config import cuttag as cuttag_validation
from encode_pipeline.config import defaults
from encode_pipeline.config import errors
from encode_pipeline.config import genome as genome_validation
from encode_pipeline.config import mnase as mnase_validation
from encode_pipeline.config import qc as qc_validation
from encode_pipeline.config import reproducibility as reproducibility_validation
from encode_pipeline.config import tools as tools_validation
from encode_pipeline.config.coercion import coerce_int


ValidationError = errors.ValidationError

# Keep module-level aliases for backward compatibility with any code that may
# have imported these private names. These aliases are deprecated; prefer
# encode_pipeline.config.defaults directly.
_SAMPLE_ID_RE = defaults.SAMPLE_ID_RE
_SANITIZE_RE = defaults.SANITIZE_RE
_BT2_STANDARD = defaults.BT2_STANDARD
_BT2_LARGE = defaults.BT2_LARGE


# ---------------------------------------------------------------------------
# Sample loading compatibility surface (lazy to avoid circular imports)
# ---------------------------------------------------------------------------

_load_and_validate_samples = None
_validate_replicate_groups = None


def __getattr__(name: str):
    """Lazy re-exports for sample-loading functions.

    Importing config.validator during config/__init__ must not trigger a
    partially initialized samples.load module. These aliases are only resolved
    on first attribute access.
    """
    global _load_and_validate_samples, _validate_replicate_groups
    if name == "load_and_validate_samples":
        if _load_and_validate_samples is None:
            from encode_pipeline.samples.load import load_and_validate_samples
            _load_and_validate_samples = load_and_validate_samples
        return _load_and_validate_samples
    if name == "validate_replicate_groups":
        if _validate_replicate_groups is None:
            from encode_pipeline.samples.replicates import validate_replicate_groups
            _validate_replicate_groups = validate_replicate_groups
        return _validate_replicate_groups
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------

def _coerce_int(value, *, name: str, minimum: int) -> int:
    """Return a strictly parsed integer, rejecting bools and floats."""
    return coerce_int(
        value,
        name=name,
        minimum=minimum,
        error_cls=ValidationError,
    )


def _validate_effective_genome_size(genome: str, value) -> None:
    """Validate MACS3 effective genome size shortcut or positive integer."""
    return genome_validation.validate_effective_genome_size(
        genome,
        value,
        error_cls=ValidationError,
    )


def validate_config(config: dict) -> dict:
    """Validate and normalize workflow configuration.

    Returns a normalized dict with validated values.
    Raises ValidationError on invalid input.
    """
    validated: dict = {}

    # samples — required, path must exist
    samples_path = config.get("samples", "")
    if not samples_path:
        raise ValidationError(
            "config.yaml must set 'samples' to the TSV path."
        )
    if not os.path.isfile(samples_path):
        raise ValidationError(
            f"Sample sheet not found: {samples_path}"
        )
    validated["samples"] = samples_path

    # outdir — string, default "results"
    validated["outdir"] = str(config.get("outdir", "results"))

    # threads — positive integer
    threads = _coerce_int(
        config.get("threads", 8),
        name="threads",
        minimum=1,
    )
    validated["threads"] = threads

    # mapq — non-negative integer
    mapq = _coerce_int(
        config.get("mapq", 30),
        name="mapq",
        minimum=0,
    )
    validated["mapq"] = mapq

    # binsize — positive integer
    binsize = _coerce_int(
        config.get("binsize", 10),
        name="binsize",
        minimum=1,
    )
    validated["binsize"] = binsize

    # remove_dup — "auto", "yes", or "no"
    remove_dup = str(config.get("remove_dup", "auto"))
    if remove_dup not in defaults.REMOVE_DUP_KEYWORDS:
        raise ValidationError(
            f"config remove_dup must be auto, yes, or no, "
            f"got {remove_dup!r}"
        )
    validated["remove_dup"] = remove_dup

    # trim — boolean or string boolean, normalize to string
    trim_raw = config.get("trim", True)
    if isinstance(trim_raw, bool):
        validated["trim"] = str(trim_raw).lower()
    elif str(trim_raw).lower() in ("true", "false"):
        validated["trim"] = str(trim_raw).lower()
    else:
        raise ValidationError(
            f"config trim must be true or false, got {trim_raw!r}"
        )

    # extend_reads — "auto", "yes", "no", or positive integer string
    ext_raw = str(config.get("extend_reads", "auto"))
    if ext_raw not in defaults.EXTEND_READS_KEYWORDS and not (
        ext_raw.isdigit() and int(ext_raw) > 0
    ):
        raise ValidationError(
            f"config extend_reads must be auto, yes, no, or a "
            f"positive integer, got {ext_raw!r}"
        )
    validated["extend_reads"] = ext_raw

    # use_control — boolean or string boolean, normalize to bool
    use_raw = config.get("use_control", False)
    if isinstance(use_raw, bool):
        validated["use_control"] = use_raw
    else:
        _use = str(use_raw).lower()
        if _use in ("true", "yes", "1"):
            validated["use_control"] = True
        elif _use in ("false", "no", "0"):
            validated["use_control"] = False
        else:
            raise ValidationError(
                f"config use_control must be true or false, "
                f"got {use_raw!r}"
            )

    # multiqc — boolean or string boolean, normalize to bool
    mqc_raw = config.get("multiqc", True)
    if isinstance(mqc_raw, bool):
        validated["multiqc"] = mqc_raw
    elif str(mqc_raw).lower() in ("true", "false"):
        validated["multiqc"] = str(mqc_raw).lower() == "true"
    else:
        raise ValidationError(
            f"config multiqc must be true or false, got {mqc_raw!r}"
        )

    # genome_resources — optional, validate if present
    validated["genome_resources"] = _validate_genome_resources(
        config.get("genome_resources", {})
    )

    # qc — optional QC switches (legacy Stage 3 defaults true; heavier modules false)
    validated["qc"] = _validate_qc_config(config.get("qc", {}))

    # stage4b — optional Stage 4b replicate-aware outputs, default true
    stage4b_raw = config.get("stage4b", True)
    if isinstance(stage4b_raw, bool):
        validated["stage4b"] = stage4b_raw
    elif str(stage4b_raw).lower() in ("true", "false"):
        validated["stage4b"] = str(stage4b_raw).lower() == "true"
    else:
        raise ValidationError(
            f"config stage4b must be true or false, got {stage4b_raw!r}"
        )

    # tool_parameters — optional Stage 4c structured tool config, default empty
    validated["tool_parameters"] = _validate_tool_params(
        config.get("tool_parameters", {})
    )

    # stage5 — optional Stage 5 IDR, default false
    stage5_raw = config.get("stage5", False)
    if isinstance(stage5_raw, bool):
        validated["stage5"] = stage5_raw
    elif str(stage5_raw).lower() in ("true", "false"):
        validated["stage5"] = str(stage5_raw).lower() == "true"
    else:
        raise ValidationError(
            f"config stage5 must be true or false, got {stage5_raw!r}"
        )

    # stage5 requires stage4b
    if validated["stage5"] and not validated.get("stage4b", True):
        raise ValidationError(
            "config: stage5=true requires stage4b=true. "
            "Stage 5a depends on Stage 4b biorep BAMs and pooled control BAMs."
        )

    # reproducibility — Stage 53+ replicate-validated peak outputs
    validated["reproducibility"] = _validate_reproducibility(
        config.get("reproducibility", {}), validated
    )

    # Stage 55: determine whether ATAC IDR is enabled from validated reproducibility
    repro = validated["reproducibility"]
    atac_idr_enabled = (
        repro.get("enabled", False)
        and repro.get("idr", {}).get("atac_narrow", False)
    )

    # ATAC IDR requires stage4b (Stage 55)
    if atac_idr_enabled and not validated.get("stage4b", True):
        raise ValidationError(
            "config: reproducibility.idr.atac_narrow=true requires "
            "stage4b=true."
        )

    # Stage 64: determine whether CUT&Tag IDR is enabled
    cuttag_idr_enabled = (
        repro.get("enabled", False)
        and repro.get("idr", {}).get("cuttag_narrow", False)
    )

    # CUT&Tag IDR requires stage4b (Stage 64)
    if cuttag_idr_enabled and not validated.get("stage4b", True):
        raise ValidationError(
            "config: reproducibility.idr.cuttag_narrow=true requires "
            "stage4b=true."
        )

    # Stage 65: determine whether broad IDR is enabled
    broad_chipseq_idr_enabled = (
        repro.get("enabled", False)
        and repro.get("idr", {}).get("chipseq_broad_experimental", False)
    )
    broad_cuttag_idr_enabled = (
        repro.get("enabled", False)
        and repro.get("idr", {}).get("cuttag_broad_experimental", False)
    )

    # Broad IDR requires stage4b (Stage 65)
    if (broad_chipseq_idr_enabled or broad_cuttag_idr_enabled) \
       and not validated.get("stage4b", True):
        raise ValidationError(
            "config: reproducibility.idr.chipseq_broad_experimental=true "
            "or cuttag_broad_experimental=true requires stage4b=true."
        )

    # idr settings validated when stage5 or any IDR mode is enabled
    if (validated["stage5"] or atac_idr_enabled or cuttag_idr_enabled
        or broad_chipseq_idr_enabled or broad_cuttag_idr_enabled):
        validated["idr"] = _validate_idr_settings(config.get("idr", {}))
    else:
        validated["idr"] = {"threshold": 0.05, "rank": "p.value", "seed": 42}

    # cuttag — optional CUT&Tag-specific config (Stage 7b)
    validated["cuttag"] = _validate_cuttag_config(config.get("cuttag", {}))

    # mnase — optional MNase-seq config (Stage 39)
    validated["mnase"] = _validate_mnase_config(config.get("mnase", {}))

    return validated


def _validate_cuttag_config(cuttag: dict) -> dict:
    """Validate the cuttag config block. Returns normalized dict.

    Absent block → all defaults. Only validates keys for Stage 7b.
    """
    return cuttag_validation.validate_cuttag_config(
        cuttag,
        error_cls=ValidationError,
    )


def _validate_mnase_config(mnase: dict) -> dict:
    """Validate the mnase config block. Returns normalized dict.

    Absent block -> all defaults. Stage 39 keys + Stage 40 fragments,
    dyad_range, and callers.

    Fragment range precedence: fragments.mono > mono_range > [140, 200].
    """
    return mnase_validation.validate_mnase_config(
        mnase,
        error_cls=ValidationError,
    )


def _validate_qc_config(qc: dict) -> dict:
    """Validate and normalize the qc config block.

    Thin wrapper around encode_pipeline.config.qc.validate_qc_config.
    """
    return qc_validation.validate_qc_config(qc, error_cls=ValidationError)


def validate_picard_reference_resources(
    validated_config: dict, samples: list[dict]
):
    """Validate reference_fasta when qc.picard_metrics is enabled.

    Scans treatment sample genomes and checks that each has a non-empty
    reference_fasta in genome_resources. Raises ValidationError if
    qc.picard_metrics is true and any treatment genome is missing
    reference_fasta. Does nothing when qc.picard_metrics is false.
    """
    return genome_validation.validate_picard_reference_resources(
        validated_config,
        samples,
        error_cls=ValidationError,
    )


def validate_tss_annotation_resources(validated_config: dict, samples: list[dict]):
    """Validate GTF annotation when qc.tss_enrichment is enabled.

    Scans treatment sample genomes and checks that each has a non-empty gtf in
    genome_resources. Raises ValidationError if qc.tss_enrichment is true and
    any treatment genome is missing a GTF annotation. Does nothing when
    qc.tss_enrichment is false.
    """
    return genome_validation.validate_tss_annotation_resources(
        validated_config,
        samples,
        error_cls=ValidationError,
    )


def _validate_genome_resources(resources: dict) -> dict:
    """Validate genome_resources config block.

    Returns the resources dict unchanged if valid.
    Raises ValidationError on invalid entries.
    """
    return genome_validation.validate_genome_resources(
        resources,
        error_cls=ValidationError,
    )


def _validate_tool_params(tool_params) -> dict:
    """Validate and normalize the Stage 4c tool_parameters config block.

    Returns a normalized dict keyed by tool name.
    Raises ValidationError on unknown tools, unknown keys, or invalid types.
    Missing tool_parameters or empty blocks are valid and return an empty dict.
    """
    return tools_validation.validate_tool_params(
        tool_params,
        error_cls=ValidationError,
    )


def _validate_idr_settings(idr):
    """Validate the idr config block. Returns normalized dict.

    Only called when stage5 is true.
    """
    return reproducibility_validation.validate_idr_settings(
        idr,
        error_cls=ValidationError,
    )


def _validate_reproducibility(raw, validated_config):
    """Validate the reproducibility config block (Stage 53+).

    Returns a dict with validated reproducibility settings.
    When enabled is false/absent, returns {'enabled': False} without
    validating sub-keys.
    """
    return reproducibility_validation.validate_reproducibility(
        raw,
        validated_config,
        error_cls=ValidationError,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _load_yaml(path: str) -> dict:
    """Load YAML config, preferring PyYAML if available, with stdlib fallback.

    PyYAML is available in the chipseq Conda environment (transitive
    dependency of Snakemake). The stdlib fallback handles bare-minimum
    YAML parsing for standalone use outside Conda.
    """
    try:
        import yaml
    except ImportError:
        return _parse_config_minimal(path)

    with open(path) as fh:
        return yaml.safe_load(fh)


def _parse_config_minimal(path: str) -> dict:
    """Minimal YAML config parser using stdlib only.

    Extracts flat keys plus the genome_resources and qc nested blocks.
    This covers the config structure used by validate_config().
    """
    config: dict = {}
    genome_resources: dict = {}
    qc: dict = {}
    tool_parameters: dict = {}
    cuttag: dict = {}
    seacr_sub: dict = {}
    section: str | None = None
    current_genome: str | None = None
    current_entry: dict | None = None
    current_tool: str | None = None
    current_tool_entry: dict | None = None

    def parse_scalar(value: str):
        value = value.strip().strip('"').strip("'")
        if value == "":
            return ""
        if value.lower() in ("true", "false"):
            return value.lower() == "true"
        if value.isdigit():
            return int(value)
        return value

    def save_current_genome() -> None:
        nonlocal current_genome, current_entry
        if current_genome is not None and current_entry is not None:
            genome_resources[current_genome] = current_entry
        current_genome = None
        current_entry = None

    def save_current_tool() -> None:
        nonlocal current_tool, current_tool_entry
        if current_tool is not None and current_tool_entry is not None:
            tool_parameters[current_tool] = current_tool_entry
        current_tool = None
        current_tool_entry = None

    with open(path) as fh:
        for line in fh:
            stripped = line.rstrip()

            # Skip comments and empty lines
            if not stripped or stripped.lstrip().startswith("#"):
                continue

            indent = len(line) - len(line.lstrip(" "))

            # Top-level keys start or end nested blocks.
            if indent == 0:
                save_current_genome()
                save_current_tool()

                if stripped.startswith("genome_resources:"):
                    section = "genome_resources"
                    continue
                if stripped.startswith("qc:"):
                    section = "qc"
                    continue
                if stripped.startswith("tool_parameters:"):
                    section = "tool_parameters"
                    continue
                if stripped.startswith("cuttag:"):
                    section = "cuttag"
                    continue

                section = None
                if ":" in stripped:
                    k, v = stripped.split(":", 1)
                    config[k.strip()] = parse_scalar(v)
                continue

            if section == "genome_resources" and indent == 2:
                # Genome key: "  hs:", "  mm10:", etc.
                key = stripped.strip().rstrip(":")
                if key and not key.startswith("#"):
                    save_current_genome()
                    current_genome = key
                    current_entry = {}
                continue

            if section == "genome_resources" and indent == 4:
                # Field within a genome entry
                field = stripped.strip()
                if ":" in field and current_genome is not None and current_entry is not None:
                    k, v = field.split(":", 1)
                    current_entry[k.strip()] = parse_scalar(v)
                continue

            if section == "qc" and indent == 2:
                field = stripped.strip()
                if ":" in field:
                    k, v = field.split(":", 1)
                    qc[k.strip()] = parse_scalar(v)
                continue

            if section == "cuttag" and indent == 2:
                key = stripped.strip().rstrip(":")
                if key == "seacr":
                    seacr_sub = {}
                elif ":" in stripped:
                    k, v = stripped.split(":", 1)
                    cuttag[k.strip()] = parse_scalar(v)
                continue

            if section == "cuttag" and indent == 4:
                field = stripped.strip()
                if ":" in field:
                    k, v = field.split(":", 1)
                    seacr_sub[k.strip()] = parse_scalar(v)
                continue

            if section == "tool_parameters" and indent == 2:
                # Tool block key: "  fastqc:", "  trim_galore:", etc.
                key = stripped.strip().rstrip(":")
                if key and not key.startswith("#"):
                    save_current_tool()
                    current_tool = key
                    current_tool_entry = {}
                continue

            if section == "tool_parameters" and indent == 4:
                # Field within a tool block
                field = stripped.strip()
                if ":" in field and current_tool is not None and current_tool_entry is not None:
                    k, v = field.split(":", 1)
                    current_tool_entry[k.strip()] = parse_scalar(v)
                continue

    # Save last genome entry
    save_current_genome()
    save_current_tool()

    if genome_resources:
        config["genome_resources"] = genome_resources
    if qc:
        config["qc"] = qc
    if tool_parameters:
        config["tool_parameters"] = tool_parameters
    if seacr_sub:
        cuttag["seacr"] = seacr_sub
    if cuttag:
        config["cuttag"] = cuttag

    return config


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate ChIP-seq/CUT&Tag/ATAC-seq pipeline config and "
                    "sample sheet."
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to config YAML (e.g. config/config.yaml)"
    )
    parser.add_argument(
        "--strict-inputs", action="store_true", default=False,
        help="Validate FASTQ and Bowtie2 index file existence"
    )
    args = parser.parse_args(argv)

    config_path = args.config
    if not os.path.isfile(config_path):
        print(f"ERROR: Config file not found: {config_path}",
              file=sys.stderr)
        sys.exit(1)

    config = _load_yaml(config_path)

    try:
        validated = validate_config(config)
        repro = validated.get("reproducibility", {})
        atac_idr_enabled = (
            repro.get("enabled", False)
            and repro.get("idr", {}).get("atac_narrow", False)
        )
        cuttag_idr_enabled = (
            repro.get("enabled", False)
            and repro.get("idr", {}).get("cuttag_narrow", False)
        )
        broad_chipseq_idr_enabled = (
            repro.get("enabled", False)
            and repro.get("idr", {}).get(
                "chipseq_broad_experimental", False)
        )
        broad_cuttag_idr_enabled = (
            repro.get("enabled", False)
            and repro.get("idr", {}).get(
                "cuttag_broad_experimental", False)
        )
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
        validate_picard_reference_resources(validated, samples)
        validate_tss_annotation_resources(validated, samples)
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    n_treatment = sum(1 for s in samples if s["role"] == "treatment")
    n_control = sum(1 for s in samples if s["role"] == "control")
    print(f"OK: {len(samples)} sample(s) validated "
          f"({n_treatment} treatment, {n_control} control)")
    print(f"     use_control: {validated['use_control']}")
    if validated.get("genome_resources"):
        genomes = ", ".join(validated["genome_resources"])
        print(f"     genome resources: {genomes}")


if __name__ == "__main__":
    main()
