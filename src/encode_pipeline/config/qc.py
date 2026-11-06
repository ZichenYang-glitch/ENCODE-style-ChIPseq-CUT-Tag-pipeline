"""QC config validation helpers."""

from encode_pipeline.config.coercion import coerce_bool

__all__ = ["validate_qc_config"]


# Known QC config keys and their default values. These are intentionally pinned
# here rather than in defaults.py because they are tied to the normalization
# behavior of validate_qc_config.
_QC_DEFAULTS = {
    "blacklist_filter": False,
    "frip": True,
    "library_complexity": True,
    "nrf_pbc": True,
    "signal_tracks": True,
    "summary": True,
    "cuttag_fragment_size": True,
    "cross_correlation": False,
    "cross_correlation_target_type": "auto",
    "preseq_complexity": False,
    "picard_metrics": False,
    "tss_enrichment": False,
}

_QC_STRING_DEFAULTS = {
    "cross_correlation_target_type": frozenset({"auto", "tf", "histone"}),
}


def validate_qc_config(qc: dict, error_cls=ValueError) -> dict:
    """Validate and normalize the qc config block.

    Missing keys default to the values in _QC_DEFAULTS. Accepts boolean or
    string boolean. Returns a normalized dict with boolean values.
    Unknown keys are rejected so typos cannot silently disable QC switches.
    """
    if not isinstance(qc, dict):
        raise error_cls(f"qc must be a mapping, got {type(qc).__name__}")

    for key in qc:
        if key not in _QC_DEFAULTS:
            raise error_cls(f"qc: unknown key {key!r}. Known: {sorted(_QC_DEFAULTS)}")

    normalized = {}
    for key, default in _QC_DEFAULTS.items():
        if key in _QC_STRING_DEFAULTS:
            value = qc.get(key, default)
            if not isinstance(value, str) or value not in _QC_STRING_DEFAULTS[key]:
                raise error_cls(
                    f"qc.{key} must be one of "
                    f"{sorted(_QC_STRING_DEFAULTS[key])}, got {value!r}"
                )
            normalized[key] = value
        else:
            normalized[key] = coerce_bool(
                qc.get(key, default),
                name=f"qc.{key}",
                error_cls=error_cls,
            )
    return normalized
