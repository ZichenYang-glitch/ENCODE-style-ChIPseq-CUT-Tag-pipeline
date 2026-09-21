#!/usr/bin/env python3
"""Sample sheet and config validation for the ChIP-seq / CUT&Tag / ATAC-seq / MNase-seq pipeline.

Importable by workflow/Snakefile and runnable as a standalone CLI.

Usage:
    python3 scripts/validate_samples.py --config config/config.yaml
"""

import os

from encode_pipeline.config import cuttag as cuttag_validation
from encode_pipeline.config import defaults
from encode_pipeline.config import errors
from encode_pipeline.config import genome as genome_validation
from encode_pipeline.config import mnase as mnase_validation
from encode_pipeline.config import qc as qc_validation
from encode_pipeline.config import reproducibility as reproducibility_validation
from encode_pipeline.config import tools as tools_validation
from encode_pipeline.config.coercion import coerce_int
from encode_pipeline.config.yaml_loader import load_yaml as _load_yaml  # noqa: F401
from encode_pipeline.config.yaml_loader import (  # noqa: F401
    parse_config_minimal as _parse_config_minimal,
)


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


def validate_config(config: dict) -> dict:
    """Validate and normalize workflow configuration.

    Returns a normalized dict with validated values.
    Raises ValidationError on invalid input.
    """
    validated: dict = {}

    # samples — required, path must exist
    samples_path = config.get("samples", "")
    if not samples_path:
        raise ValidationError("config.yaml must set 'samples' to the TSV path.")
    if not os.path.isfile(samples_path):
        raise ValidationError(f"Sample sheet not found: {samples_path}")
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
            f"config remove_dup must be auto, yes, or no, got {remove_dup!r}"
        )
    validated["remove_dup"] = remove_dup

    # trim — boolean or string boolean, normalize to string
    trim_raw = config.get("trim", True)
    if isinstance(trim_raw, bool):
        validated["trim"] = str(trim_raw).lower()
    elif str(trim_raw).lower() in ("true", "false"):
        validated["trim"] = str(trim_raw).lower()
    else:
        raise ValidationError(f"config trim must be true or false, got {trim_raw!r}")

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
                f"config use_control must be true or false, got {use_raw!r}"
            )

    # multiqc — boolean or string boolean, normalize to bool
    mqc_raw = config.get("multiqc", True)
    if isinstance(mqc_raw, bool):
        validated["multiqc"] = mqc_raw
    elif str(mqc_raw).lower() in ("true", "false"):
        validated["multiqc"] = str(mqc_raw).lower() == "true"
    else:
        raise ValidationError(f"config multiqc must be true or false, got {mqc_raw!r}")

    # genome_resources — optional, validate if present
    validated["genome_resources"] = _validate_genome_resources(
        config.get("genome_resources", {})
    )

    # qc — optional QC switches (standard defaults true; heavier modules false)
    validated["qc"] = _validate_qc_config(config.get("qc", {}))

    # replicate_analysis — optional replicate-aware outputs, default true
    replicate_analysis_raw = config.get("replicate_analysis", True)
    if isinstance(replicate_analysis_raw, bool):
        validated["replicate_analysis"] = replicate_analysis_raw
    elif str(replicate_analysis_raw).lower() in ("true", "false"):
        validated["replicate_analysis"] = str(replicate_analysis_raw).lower() == "true"
    else:
        raise ValidationError(
            f"config replicate_analysis must be true or false, got {replicate_analysis_raw!r}"
        )

    # tool_parameters — optional Stage 4c structured tool config, default empty
    validated["tool_parameters"] = _validate_tool_params(
        config.get("tool_parameters", {})
    )

    # chipseq_idr — optional ChIP-seq IDR, default false
    chipseq_idr_raw = config.get("chipseq_idr", False)
    if isinstance(chipseq_idr_raw, bool):
        validated["chipseq_idr"] = chipseq_idr_raw
    elif str(chipseq_idr_raw).lower() in ("true", "false"):
        validated["chipseq_idr"] = str(chipseq_idr_raw).lower() == "true"
    else:
        raise ValidationError(
            f"config chipseq_idr must be true or false, got {chipseq_idr_raw!r}"
        )

    # chipseq_idr requires replicate_analysis
    if validated["chipseq_idr"] and not validated.get("replicate_analysis", True):
        raise ValidationError(
            "config: chipseq_idr=true requires replicate_analysis=true. "
            "ChIP-seq IDR depends on replicate analysis biorep BAMs and pooled control BAMs."
        )

    # reproducibility — replicate-validated peak outputs
    validated["reproducibility"] = _validate_reproducibility(
        config.get("reproducibility", {}), validated
    )

    # Determine whether ATAC IDR is enabled from validated reproducibility
    repro = validated["reproducibility"]
    atac_idr_enabled = repro.get("enabled", False) and repro.get("idr", {}).get(
        "atac_narrow", False
    )

    # ATAC IDR requires replicate_analysis (ATAC narrow)
    if atac_idr_enabled and not validated.get("replicate_analysis", True):
        raise ValidationError(
            "config: reproducibility.idr.atac_narrow=true requires replicate_analysis=true."
        )

    # Stage 64: determine whether CUT&Tag IDR is enabled
    cuttag_idr_enabled = repro.get("enabled", False) and repro.get("idr", {}).get(
        "cuttag_narrow", False
    )

    # CUT&Tag IDR requires replicate_analysis (Stage 64)
    if cuttag_idr_enabled and not validated.get("replicate_analysis", True):
        raise ValidationError(
            "config: reproducibility.idr.cuttag_narrow=true requires replicate_analysis=true."
        )

    # Stage 65: determine whether broad IDR is enabled
    broad_chipseq_idr_enabled = repro.get("enabled", False) and repro.get(
        "idr", {}
    ).get("chipseq_broad_experimental", False)
    broad_cuttag_idr_enabled = repro.get("enabled", False) and repro.get("idr", {}).get(
        "cuttag_broad_experimental", False
    )

    # Broad IDR requires replicate_analysis (Stage 65)
    if (broad_chipseq_idr_enabled or broad_cuttag_idr_enabled) and not validated.get(
        "replicate_analysis", True
    ):
        raise ValidationError(
            "config: reproducibility.idr.chipseq_broad_experimental=true "
            "or cuttag_broad_experimental=true requires replicate_analysis=true."
        )

    # idr settings validated when chipseq_idr or any IDR mode is enabled
    if (
        validated["chipseq_idr"]
        or atac_idr_enabled
        or cuttag_idr_enabled
        or broad_chipseq_idr_enabled
        or broad_cuttag_idr_enabled
    ):
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


def validate_picard_reference_resources(validated_config: dict, samples: list[dict]):
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

    Only called when chipseq_idr is true.
    """
    return reproducibility_validation.validate_idr_settings(
        idr,
        error_cls=ValidationError,
    )


def _validate_reproducibility(raw, validated_config):
    """Validate the reproducibility config block.

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
# YAML loading compatibility surface
# ---------------------------------------------------------------------------

# The imports above retain deprecated aliases for code that may still import
# these private names.


# ---------------------------------------------------------------------------
# CLI compatibility surface (lazy to avoid circular imports)
# ---------------------------------------------------------------------------


def main(argv=None):
    """Compatibility wrapper around ``encode_pipeline.cli._validator.main``."""
    from encode_pipeline.cli._validator import main as _main

    return _main(argv)


# _load_yaml and _parse_config_minimal were moved to
# encode_pipeline.config.yaml_loader. The aliases above preserve compatibility.
# main() is lazily delegated to encode_pipeline.cli._validator.main to avoid
# a config.validator -> cli._validator -> config.validator import cycle at
# module load time.
