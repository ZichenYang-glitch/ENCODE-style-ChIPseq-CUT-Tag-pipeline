"""Reference Profile domain contract tests."""

from __future__ import annotations

from dataclasses import replace
from datetime import datetime, timezone

import pytest

from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.reference_profiles import (
    ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME,
    REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME,
    AdapterReferenceBindingIdentity,
    BoundWorkflowReference,
    ReferenceProfile,
    ReferenceProfileWorkflowBinding,
    build_reference_profile_revision,
    build_reference_profile_revision_binding,
)


NOW = datetime(2026, 8, 3, 8, 0, tzinfo=timezone.utc)
PROFILE_ID = "refp_" + "1" * 32
REVISION_ID = "refpr_" + "2" * 32


def _adapter_identity() -> AdapterReferenceBindingIdentity:
    return AdapterReferenceBindingIdentity(
        workflow_id="bulk-rnaseq",
        contract_version="bulk-rnaseq-reference-v1",
        identity_sha256="3" * 64,
    )


def test_profile_revision_public_identity_is_path_free_and_binding_sensitive() -> None:
    binding = ReferenceProfileWorkflowBinding.from_adapter_identity(_adapter_identity())
    revision = build_reference_profile_revision(
        revision_id=REVISION_ID,
        profile_id=PROFILE_ID,
        revision_number=1,
        display_name="GRCh38 primary",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-primary",
        workflow_bindings=(binding,),
        created_at=NOW,
    )

    assert revision.public_identity_scheme == (
        REFERENCE_PROFILE_REVISION_IDENTITY_SCHEME
    )
    assert len(revision.public_identity_sha256) == 64
    assert revision.workflow_ids == ("bulk-rnaseq",)
    assert "grch38-primary" not in repr(revision)
    assert (
        "config_key"
        not in revision.public_summary(safe_key="grch38", enabled=True).to_public_dict()
    )

    changed = build_reference_profile_revision(
        revision_id="refpr_" + "4" * 32,
        profile_id=PROFILE_ID,
        revision_number=2,
        display_name="GRCh38 primary",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-primary",
        workflow_bindings=(replace(binding, identity_sha256="5" * 64),),
        created_at=NOW,
    )
    assert changed.public_identity_sha256 != revision.public_identity_sha256


def test_profile_and_revision_reject_unsafe_or_ambiguous_identity() -> None:
    with pytest.raises(ValueError, match="safe_key"):
        ReferenceProfile(
            profile_id=PROFILE_ID,
            safe_key="../secret",
            created_at=NOW,
        )

    with pytest.raises(ValueError, match="duplicate workflow"):
        build_reference_profile_revision(
            revision_id=REVISION_ID,
            profile_id=PROFILE_ID,
            revision_number=1,
            display_name="GRCh38",
            organism="Homo sapiens",
            assembly="GRCh38",
            config_key="grch38",
            workflow_bindings=(
                ReferenceProfileWorkflowBinding.from_adapter_identity(
                    _adapter_identity()
                ),
                ReferenceProfileWorkflowBinding.from_adapter_identity(
                    _adapter_identity()
                ),
            ),
            created_at=NOW,
        )


def test_bound_reference_and_snapshot_run_evidence_pin_exact_identity() -> None:
    identity = _adapter_identity()
    inputs = WorkflowInputs(config={}, samples=None, options={})
    adapter = object()
    bound = BoundWorkflowReference(
        inputs=inputs,
        adapter=adapter,
        identity=identity,
    )
    evidence = build_reference_profile_revision_binding(
        profile_id=PROFILE_ID,
        revision_id=REVISION_ID,
        workflow_id="bulk-rnaseq",
        revision_public_identity_sha256="6" * 64,
        adapter_identity=identity,
        bound_at=NOW,
    )

    assert bound.inputs is inputs
    assert bound.adapter is adapter
    assert bound.identity.identity_scheme == (ADAPTER_REFERENCE_BINDING_IDENTITY_SCHEME)
    assert evidence.revision_id == REVISION_ID
    assert len(evidence.binding_digest) == 64
    with pytest.raises(ValueError, match="binding_digest"):
        replace(evidence, binding_digest="0" * 64)
