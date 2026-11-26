"""Reproducibility and IDR validation helpers."""

import warnings

from encode_pipeline.config import defaults
from encode_pipeline.config.coercion import coerce_bool, coerce_int

__all__ = [
    "idr_sample_flags_from_validated",
    "resolve_idr_sample_flags",
    "validate_idr_config",
    "validate_idr_settings",
    "validate_reproducibility",
]


def idr_sample_flags_from_validated(config):
    """Select sample-validation switches from normalized IDR settings."""
    idr = config["reproducibility"].get("idr", {})
    return {
        "chipseq_narrow": config["chipseq_idr"],
        "atac_narrow": idr.get("atac_narrow", False),
        "cuttag_narrow": idr.get("cuttag_narrow", False),
        "chipseq_broad": idr.get("chipseq_broad_experimental", False),
        "cuttag_broad": idr.get("cuttag_broad_experimental", False),
    }


def resolve_idr_sample_flags(config, *, error_cls=ValueError):
    """Resolve sample admission for callers that receive a raw config or run record."""
    if not isinstance(config, dict):
        raise error_cls(f"config must be a mapping, got {type(config).__name__}")
    replicate_analysis = coerce_bool(
        config.get("replicate_analysis", True), "config replicate_analysis", error_cls
    )
    normalized = validate_idr_config(
        config, replicate_analysis, error_cls=error_cls, emit_warnings=False
    )
    return idr_sample_flags_from_validated(normalized)


def validate_idr_config(
    config, replicate_analysis, *, error_cls=ValueError, emit_warnings=True
):
    """Normalize the IDR and consensus settings that govern replicate admission."""
    chipseq_idr = coerce_bool(
        config.get("chipseq_idr", False), "config chipseq_idr", error_cls
    )
    if chipseq_idr and not replicate_analysis:
        raise error_cls(
            "config: chipseq_idr=true requires replicate_analysis=true. "
            "ChIP-seq IDR depends on replicate analysis biorep BAMs and pooled control BAMs."
        )

    repro = validate_reproducibility(
        config.get("reproducibility", {}),
        {"chipseq_idr": chipseq_idr},
        error_cls=error_cls,
        emit_warnings=emit_warnings,
    )
    idr = repro.get("idr", {})
    if repro.get("mixed_replicate_policy") and not (
        chipseq_idr
        or any(
            idr.get(key) for key in ("chipseq_narrow", "cuttag_narrow", "atac_narrow")
        )
    ):
        raise error_cls(
            "per_group_fallback requires an enabled narrow-peak MACS3 IDR method."
        )
    if not replicate_analysis:
        for key in ("chipseq_narrow", "atac_narrow", "cuttag_narrow"):
            if idr.get(key):
                raise error_cls(
                    f"config: reproducibility.idr.{key}=true requires replicate_analysis=true."
                )
        if idr.get("chipseq_broad_experimental") or idr.get(
            "cuttag_broad_experimental"
        ):
            raise error_cls(
                "config: reproducibility.idr.chipseq_broad_experimental=true "
                "or cuttag_broad_experimental=true requires replicate_analysis=true."
            )

    active = chipseq_idr or any(idr.values())
    settings = validate_idr_settings(config.get("idr", {}) if active else {}, error_cls)
    return {
        "chipseq_idr": chipseq_idr or bool(idr.get("chipseq_narrow")),
        "reproducibility": repro,
        "idr": settings,
    }


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


def validate_reproducibility(
    raw, validated_config, error_cls=ValueError, *, emit_warnings=True
):
    """Validate the reproducibility config block and return normalized settings."""
    if not isinstance(raw, dict):
        raise error_cls(f"reproducibility must be a mapping, got {type(raw).__name__}")

    enabled_raw = raw.get("enabled", False)
    enabled = coerce_bool(
        enabled_raw,
        "reproducibility.enabled",
        error_cls=error_cls,
    )
    policy = raw.get("mixed_replicate_policy", "strict")
    if policy not in ("strict", "per_group_fallback"):
        raise error_cls(
            "reproducibility.mixed_replicate_policy must be 'strict' or 'per_group_fallback'."
        )
    mixed = policy == "per_group_fallback"
    if mixed and not enabled:
        raise error_cls("per_group_fallback requires reproducibility.enabled=true.")
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
    if mixed:
        if not consensus["enabled"]:
            raise error_cls(
                "per_group_fallback requires reproducibility.consensus.enabled=true."
            )
        min_reps = coerce_int(
            min_reps,
            name="reproducibility.consensus.min_replicates",
            minimum=1,
            error_cls=error_cls,
        )
        if min_reps not in (2, 3):
            raise error_cls(
                "consensus.min_replicates must be 2 or 3 in per_group_fallback."
            )
        # Omitting strict keeps both the scientific identity and the serialized
        # effective configuration identical to those of historical runs.
        result["mixed_replicate_policy"] = policy
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
        if (
            emit_warnings
            and not idr_result["chipseq_narrow"]
            and validated_config.get("chipseq_idr", False)
        ):
            warnings.warn(
                "Config contradiction: reproducibility.idr.chipseq_narrow is "
                "explicitly false but chipseq_idr is true. ChIP-seq narrow IDR will "
                "still run. Set reproducibility.idr.chipseq_narrow to null "
                "(omitted) to infer from chipseq_idr, or set chipseq_idr to false to "
                "disable ChIP-seq narrow IDR."
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
        if val and emit_warnings:
            warnings.warn(
                f"Experimental IDR flag enabled: reproducibility.idr.{flag}=true. "
                "For eligible broad-peak experiments, IDR supplies the primary "
                "result in final/. Consensus remains available when enabled."
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
