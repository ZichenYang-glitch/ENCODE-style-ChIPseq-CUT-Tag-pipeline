"""Reproducibility and IDR validation helpers."""

import warnings

from encode_pipeline.config import defaults
from encode_pipeline.config.coercion import coerce_bool

__all__ = ["validate_idr_settings", "validate_reproducibility"]


def validate_idr_settings(idr, error_cls=ValueError):
    """Validate the idr config block and return a normalized dict."""
    if not isinstance(idr, dict):
        raise error_cls(f"idr must be a mapping, got {type(idr).__name__}")

    threshold = idr.get("threshold", 0.05)
    if isinstance(threshold, bool):
        raise error_cls(f"idr.threshold must be a float in (0, 1), got {threshold!r}")
    try:
        threshold = float(threshold)
    except (ValueError, TypeError):
        raise error_cls(f"idr.threshold must be a float in (0, 1), got {threshold!r}")
    if not (0 < threshold < 1):
        raise error_cls(f"idr.threshold must be in (0, 1), got {threshold}")

    rank = str(idr.get("rank", "p.value"))
    if rank not in defaults.IDR_RANKS:
        raise error_cls(f"idr.rank must be 'p.value' or 'signal.value', got {rank!r}")

    seed = idr.get("seed", 42)
    if isinstance(seed, bool):
        raise error_cls(f"idr.seed must be a positive integer, got {seed!r}")
    if isinstance(seed, int):
        if seed <= 0:
            raise error_cls(f"idr.seed must be positive, got {seed}")
    elif isinstance(seed, str):
        text = seed.strip()
        if not text.isdigit() or int(text) <= 0:
            raise error_cls(f"idr.seed must be a positive integer, got {seed!r}")
        seed = int(text)
    else:
        raise error_cls(f"idr.seed must be a positive integer, got {seed!r}")

    known = defaults.IDR_KEYS
    for key in idr:
        if key not in known:
            raise error_cls(f"idr: unknown key {key!r}. Known: {sorted(known)}")

    return {"threshold": threshold, "rank": rank, "seed": seed}


def validate_reproducibility(raw, validated_config, error_cls=ValueError):
    """Validate the reproducibility config block and return normalized settings."""
    if not isinstance(raw, dict):
        raise error_cls(f"reproducibility must be a mapping, got {type(raw).__name__}")

    enabled_raw = raw.get("enabled", False)
    enabled = coerce_bool(
        enabled_raw,
        "reproducibility.enabled",
        error_cls=error_cls,
    )
    if not enabled:
        return {"enabled": False}

    result = {"enabled": True}

    consensus_raw = raw.get("consensus", {})
    if not isinstance(consensus_raw, dict):
        raise error_cls(
            f"reproducibility.consensus must be a mapping, "
            f"got {type(consensus_raw).__name__}"
        )
    consensus = {}
    consensus["enabled"] = coerce_bool(
        consensus_raw.get("enabled", True),
        "reproducibility.consensus.enabled",
        error_cls=error_cls,
    )

    min_reps = consensus_raw.get("min_replicates", 2)
    if isinstance(min_reps, bool):
        raise error_cls(
            f"reproducibility.consensus.min_replicates must be an integer "
            f">= 2, got {min_reps!r}"
        )
    try:
        min_reps = int(min_reps)
    except (ValueError, TypeError):
        raise error_cls(
            f"reproducibility.consensus.min_replicates must be an integer "
            f">= 2, got {min_reps!r}"
        )
    if min_reps < 2:
        raise error_cls(
            f"reproducibility.consensus.min_replicates must be >= 2, got {min_reps}"
        )
    consensus["min_replicates"] = min_reps

    overlap = consensus_raw.get("reciprocal_overlap", 0.5)
    if isinstance(overlap, bool):
        raise error_cls(
            f"reproducibility.consensus.reciprocal_overlap must be a float "
            f"in (0, 1], got {overlap!r}"
        )
    try:
        overlap = float(overlap)
    except (ValueError, TypeError):
        raise error_cls(
            f"reproducibility.consensus.reciprocal_overlap must be a float "
            f"in (0, 1], got {overlap!r}"
        )
    if not (0 < overlap <= 1):
        raise error_cls(
            f"reproducibility.consensus.reciprocal_overlap must be in (0, 1], "
            f"got {overlap}"
        )
    consensus["reciprocal_overlap"] = overlap

    known_consensus = defaults.REPRODUCIBILITY_CONSENSUS_KEYS
    for key in consensus_raw:
        if key not in known_consensus:
            raise error_cls(
                f"reproducibility.consensus: unknown key {key!r}. "
                f"Known: {sorted(known_consensus)}"
            )
    result["consensus"] = consensus

    idr_raw = raw.get("idr", {})
    if not isinstance(idr_raw, dict):
        raise error_cls(
            f"reproducibility.idr must be a mapping, got {type(idr_raw).__name__}"
        )
    idr_result = {}

    csn = idr_raw.get("chipseq_narrow", None)
    if csn is None:
        idr_result["chipseq_narrow"] = None
    else:
        idr_result["chipseq_narrow"] = coerce_bool(
            csn,
            "reproducibility.idr.chipseq_narrow",
            error_cls=error_cls,
        )
        if not idr_result["chipseq_narrow"] and validated_config.get("stage5", False):
            warnings.warn(
                "Config contradiction: reproducibility.idr.chipseq_narrow is "
                "explicitly false but stage5 is true. Legacy Stage 5 IDR will "
                "still run. Set reproducibility.idr.chipseq_narrow to null "
                "(omitted) to infer from stage5, or set stage5 to false to "
                "disable legacy IDR."
            )

    for flag in ("atac_narrow", "cuttag_narrow"):
        idr_result[flag] = coerce_bool(
            idr_raw.get(flag, False),
            f"reproducibility.idr.{flag}",
            error_cls=error_cls,
        )

    for flag in ("chipseq_broad_experimental", "cuttag_broad_experimental"):
        val = coerce_bool(
            idr_raw.get(flag, False),
            f"reproducibility.idr.{flag}",
            error_cls=error_cls,
        )
        if val:
            warnings.warn(
                f"Experimental IDR flag enabled: reproducibility.idr.{flag}=true. "
                f"Consensus remains the primary final reproducibility method for "
                f"broad-peak modes. Experimental IDR outputs appear under idr/ only "
                f"and do not replace consensus in final/."
            )
        idr_result[flag] = val

    known_idr = defaults.REPRODUCIBILITY_IDR_KEYS
    for key in idr_raw:
        if key not in known_idr:
            raise error_cls(
                f"reproducibility.idr: unknown key {key!r}. Known: {sorted(known_idr)}"
            )
    result["idr"] = idr_result

    known_top = defaults.REPRODUCIBILITY_TOP_KEYS
    for key in raw:
        if key not in known_top:
            raise error_cls(
                f"reproducibility: unknown key {key!r}. Known: {sorted(known_top)}"
            )

    return result
