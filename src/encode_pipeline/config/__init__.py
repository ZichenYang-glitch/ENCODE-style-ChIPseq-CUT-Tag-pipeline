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
        validate_config,
        validate_picard_reference_resources,
        validate_tss_annotation_resources,
    )
    return locals()[name]
