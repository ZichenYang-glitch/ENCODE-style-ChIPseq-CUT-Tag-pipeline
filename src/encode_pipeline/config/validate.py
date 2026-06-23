"""Config-level validation.

Phase 0 re-exports the legacy validator implementation. The long-term
plan is to gradually move logic from encode_pipeline.config.validator
into this module and break it into smaller helpers.
"""

from encode_pipeline.config.validator import (
    ValidationError,
    validate_config,
    validate_picard_reference_resources,
    validate_tss_annotation_resources,
)

__all__ = [
    "ValidationError",
    "validate_config",
    "validate_picard_reference_resources",
    "validate_tss_annotation_resources",
]
