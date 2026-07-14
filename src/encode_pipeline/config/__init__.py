"""Lighten encode_pipeline.config package-level imports."""

__all__ = [
    "validate_config",
    "validate_picard_reference_resources",
    "validate_tss_annotation_resources",
]

_LAZY_EXPORTS = {
    "validate_config",
    "validate_picard_reference_resources",
    "validate_tss_annotation_resources",
}


def __getattr__(name: str):
    """Lazy re-export of public config validation functions."""
    if name not in _LAZY_EXPORTS:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")
    from encode_pipeline.config.validate import (
        validate_config as _validate_config,
        validate_picard_reference_resources as _validate_picard_reference_resources,
        validate_tss_annotation_resources as _validate_tss_annotation_resources,
    )

    exports = {
        "validate_config": _validate_config,
        "validate_picard_reference_resources": _validate_picard_reference_resources,
        "validate_tss_annotation_resources": _validate_tss_annotation_resources,
    }
    return exports[name]
