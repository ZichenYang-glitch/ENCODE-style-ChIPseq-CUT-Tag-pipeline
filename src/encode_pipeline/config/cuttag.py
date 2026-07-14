"""CUT&Tag config validation helpers."""

from encode_pipeline.config import defaults

__all__ = ["validate_cuttag_config"]


_CUTTAG_DEFAULTS = {
    "peak_caller": "macs3",
    "seacr": {
        "enabled": False,
        "mode": "stringent",
        "normalization": "non",
        "threshold": 0.01,
    },
}


def validate_cuttag_config(cuttag: dict, error_cls=ValueError) -> dict:
    """Validate the cuttag config block. Returns normalized dict.

    Absent block → all defaults. Only validates keys for Stage 7b.
    """
    if not isinstance(cuttag, dict):
        raise error_cls(f"cuttag must be a mapping, got {type(cuttag).__name__}")

    known = defaults.CUTTAG_TOP_KEYS
    for key in cuttag:
        if key not in known:
            raise error_cls(f"cuttag: unknown key {key!r}. Known: {sorted(known)}")

    peak_caller = str(cuttag.get("peak_caller", _CUTTAG_DEFAULTS["peak_caller"]))
    if peak_caller != "macs3":
        raise error_cls(
            f"cuttag.peak_caller must be 'macs3' in Stage 7b, got {peak_caller!r}"
        )

    seacr_raw = cuttag.get("seacr", {})
    if isinstance(seacr_raw, bool):
        raise error_cls("cuttag.seacr must be a mapping, got boolean")
    if not isinstance(seacr_raw, dict):
        raise error_cls(
            f"cuttag.seacr must be a mapping, got {type(seacr_raw).__name__}"
        )

    seacr_known = defaults.CUTTAG_SEACR_KEYS
    for key in seacr_raw:
        if key not in seacr_known:
            raise error_cls(
                f"cuttag.seacr: unknown key {key!r}. Known: {sorted(seacr_known)}"
            )

    enabled_raw = seacr_raw.get("enabled", _CUTTAG_DEFAULTS["seacr"]["enabled"])
    if isinstance(enabled_raw, bool):
        seacr_enabled = enabled_raw
    elif str(enabled_raw).lower() in ("true", "false"):
        seacr_enabled = str(enabled_raw).lower() == "true"
    else:
        raise error_cls(
            f"cuttag.seacr.enabled must be true or false, got {enabled_raw!r}"
        )

    mode = str(seacr_raw.get("mode", _CUTTAG_DEFAULTS["seacr"]["mode"]))
    if mode not in defaults.CUTTAG_SEACR_MODES:
        raise error_cls(
            f"cuttag.seacr.mode must be 'stringent' or 'relaxed', got {mode!r}"
        )

    normalization = str(
        seacr_raw.get("normalization", _CUTTAG_DEFAULTS["seacr"]["normalization"])
    )
    if normalization != "non":
        raise error_cls(
            f"cuttag.seacr.normalization must be 'non' in Stage 7b, "
            f"got {normalization!r}"
        )

    threshold_raw = seacr_raw.get("threshold", _CUTTAG_DEFAULTS["seacr"]["threshold"])
    if isinstance(threshold_raw, bool):
        raise error_cls(
            f"cuttag.seacr.threshold must be a float in (0, 1), got {threshold_raw!r}"
        )
    try:
        threshold = float(threshold_raw)
    except (ValueError, TypeError):
        raise error_cls(
            f"cuttag.seacr.threshold must be a float in (0, 1), got {threshold_raw!r}"
        )
    if not (0 < threshold < 1):
        raise error_cls(f"cuttag.seacr.threshold must be in (0, 1), got {threshold}")

    return {
        "peak_caller": peak_caller,
        "seacr": {
            "enabled": seacr_enabled,
            "mode": mode,
            "normalization": normalization,
            "threshold": threshold,
        },
    }
