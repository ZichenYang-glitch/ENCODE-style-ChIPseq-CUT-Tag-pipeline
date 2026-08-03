"""Public Reference Profile API projections stay path-free and exact."""

from datetime import datetime, timezone

import pytest
from pydantic import ValidationError

from encode_pipeline.api.models import (
    ReferenceProfileRevisionResponse,
    ValidationRequest,
)
from encode_pipeline.platform.reference_profiles import (
    REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME,
    ReferenceProfileRevisionSummary,
)


def _summary() -> ReferenceProfileRevisionSummary:
    return ReferenceProfileRevisionSummary(
        profile_id="refp_" + "1" * 32,
        revision_id="refpr_" + "2" * 32,
        safe_key="operator-private-grch38",
        revision_number=3,
        display_name="GRCh38 laboratory reference",
        organism="Homo sapiens",
        assembly="GRCh38",
        public_identity_scheme=REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME,
        public_identity_sha256="a" * 64,
        compatible_workflow_ids=("bulk-rnaseq",),
        enabled=True,
        created_at=datetime(2026, 8, 3, tzinfo=timezone.utc),
    )


def test_reference_profile_projection_is_exact_and_omits_operator_coordinates():
    projected = ReferenceProfileRevisionResponse.from_summary(_summary())

    assert projected.model_dump() == {
        "profile_id": "refp_" + "1" * 32,
        "revision_id": "refpr_" + "2" * 32,
        "revision_number": 3,
        "display_name": "GRCh38 laboratory reference",
        "organism": "Homo sapiens",
        "assembly": "GRCh38",
        "identity_sha256": "a" * 64,
    }
    rendered = projected.model_dump_json()
    assert "safe_key" not in rendered
    assert "config_key" not in rendered
    assert "/" not in rendered


def test_validation_request_accepts_only_one_exact_revision_selection():
    revision_id = "refpr_" + "2" * 32

    request = ValidationRequest(
        config={},
        reference_profile_revision_id=revision_id,
    )

    assert request.reference_profile_revision_id == revision_id

    with pytest.raises(ValidationError):
        ValidationRequest(config={}, reference_profile_revision_id="latest")
    with pytest.raises(ValidationError):
        ValidationRequest(config={}, reference_profile_config_key="private")
