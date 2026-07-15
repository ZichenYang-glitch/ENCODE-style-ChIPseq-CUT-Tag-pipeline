"""QC config validation helpers."""

from encode_pipeline.config.coercion import coerce_bool

__all__ = ["validate_qc_config"]


# Known QC config keys and their default values. These are intentionally pinned
# here rather than in defaults.py because they are tied to the normalization
# behavior of validate_qc_config.
_QC_DEFAULTS = {
    "blacklist_filter": True,
    "frip": True,
    "library_complexity": True,
    "nrf_pbc": True,
    "signal_tracks": True,
    "summary": True,
    "cuttag_fragment_size": True,
    "cross_correlation": False,
    "preseq_complexity": False,
    "picard_metrics": False,
    "tss_enrichment": False,
}


def validate_qc_config(qc: dict, error_cls=ValueError) -> dict:
    """Validate and normalize the qc config block.

    Missing keys default to the values in _QC_DEFAULTS. Accepts boolean or
    string boolean. Returns a normalized dict with boolean values.
    Unknown keys are silently ignored.
    """
    if not isinstance(qc, dict):
        raise error_cls(f"qc must be a mapping, got {type(qc).__name__}")

    return {
        key: coerce_bool(
            qc.get(key, default),
            name=f"qc.{key}",
            error_cls=error_cls,
        )
        for key, default in _QC_DEFAULTS.items()
    }
