"""ENCODE resource-field projection for the restricted authoring assistant.

The policy names the only fields the model may see. Live schema constraints are
intersected with this policy. Descriptions, defaults, and examples stay on the
server.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from encode_pipeline.adapters.encode_authoring import MAX_CORES
from encode_pipeline.platform.authoring import (
    AuthoringField,
    FieldConstraint,
    IssueLabel,
    SAMPLE_SELECTION_REQUIRED,
    WORKFLOW_SAMPLE_INPUT_UNAVAILABLE,
    WORKFLOW_SAMPLE_PLAN_INVALID,
    WORKFLOW_SAMPLE_PLAN_REQUIRED,
)
from encode_pipeline.platform.registry import ENCODE_WORKFLOW_ID


_FIELDS = (
    ("config.threads", "Threads", "config", "config_schema", "threads", 1, None),
    ("options.cores", "Cores", "options", "option_schema", "cores", 1, MAX_CORES),
)

_ISSUES = (
    ("ENCODE_CONFIG_INVALID", "The configuration does not match the workflow schema."),
    (
        "ENCODE_CONFIG_SEMANTIC_INVALID",
        "The configuration uses an unsupported combination.",
    ),
    ("ENCODE_OPTIONS_INVALID", "The options do not match the workflow schema."),
    ("ENCODE_SAMPLES_INVALID", "The sample table does not match the workflow schema."),
    ("ENCODE_SAMPLE_WARNING", "A sample row has a non-blocking warning."),
    (
        "REFERENCE_PROFILE_REQUIRED",
        "A reference profile revision is required before validation.",
    ),
    ("REFERENCE_PROFILE_DISABLED", "The selected reference profile is disabled."),
    ("REFERENCE_PROFILE_STALE", "The selected reference profile is stale."),
    (
        "REFERENCE_PROFILE_INCOMPATIBLE",
        "The selected reference profile does not match this workflow.",
    ),
    (
        "REFERENCE_PROFILE_UNAVAILABLE",
        "The selected reference profile is unavailable.",
    ),
    (
        "INPUT_BINDING_SELECTION_INVALID",
        "The project and sample revision selection is incomplete.",
    ),
    (
        "DATA_BINDING_SELECTION_INVALID",
        "The project and sample revision selection is incomplete.",
    ),
    (
        "INPUT_USE_CAPABILITY_UNAVAILABLE",
        "Managed input selection is not available for this workflow.",
    ),
    (
        SAMPLE_SELECTION_REQUIRED,
        "Select at least one sample revision before validation.",
    ),
    (
        WORKFLOW_SAMPLE_PLAN_REQUIRED,
        "Each selected sample revision needs an explicit workflow sample plan.",
    ),
    (
        WORKFLOW_SAMPLE_PLAN_INVALID,
        "The workflow sample plan does not match the selected revisions.",
    ),
    (
        WORKFLOW_SAMPLE_INPUT_UNAVAILABLE,
        "A selected input revision could not be verified.",
    ),
)


def explainable_issues() -> tuple[IssueLabel, ...]:
    return tuple(IssueLabel(code, label) for code, label in _ISSUES)


def explainable_codes() -> frozenset[str]:
    return frozenset(code for code, _label in _ISSUES)


def project_encode_fields(
    workflow_id: str,
    schema_payload: Mapping[str, Any],
) -> tuple[AuthoringField, ...] | None:
    """Return the approved fields, or None when this workflow has no policy."""

    if workflow_id != ENCODE_WORKFLOW_ID:
        return None
    projected: list[AuthoringField] = []
    for token, label, document, schema_key, key, policy_min, policy_max in _FIELDS:
        schema = schema_payload.get(schema_key)
        constraint = _intersect_integer(schema, key, policy_min, policy_max)
        if constraint is None:
            return None
        projected.append(
            AuthoringField(
                token=token,
                label=label,
                document=document,
                key=key,
                constraint=constraint,
            )
        )
    return tuple(projected)


def _intersect_integer(
    schema: object,
    key: str,
    policy_min: int,
    policy_max: int | None,
) -> FieldConstraint | None:
    if not isinstance(schema, Mapping):
        return None
    properties = schema.get("properties")
    if not isinstance(properties, Mapping):
        return None
    prop = properties.get(key)
    if not isinstance(prop, Mapping) or prop.get("type") != "integer":
        return None
    minimum = policy_min
    schema_min = prop.get("minimum")
    if isinstance(schema_min, bool) or not isinstance(schema_min, int):
        if schema_min is not None:
            return None
    else:
        minimum = max(policy_min, schema_min)
    maximum = policy_max
    schema_max = prop.get("maximum")
    if schema_max is not None:
        if isinstance(schema_max, bool) or not isinstance(schema_max, int):
            return None
        maximum = schema_max if maximum is None else min(maximum, schema_max)
    if maximum is not None and minimum > maximum:
        return None
    return FieldConstraint(value_type="integer", minimum=minimum, maximum=maximum)
