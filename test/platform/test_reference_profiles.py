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
    ResolvedReferenceProfile,
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


@pytest.mark.parametrize(
    ("case", "expected_exception", "message"),
    (
        (
            "adapter-digest-shape",
            ValueError,
            "identity_sha256 must be a lowercase SHA-256 digest",
        ),
        ("adapter-scheme", ValueError, "identity_scheme is unsupported"),
        (
            "summary-workflow-order",
            ValueError,
            "compatible_workflow_ids must be sorted and unique",
        ),
        ("summary-enabled-type", ValueError, "enabled must be boolean"),
        (
            "revision-empty-bindings",
            ValueError,
            "workflow_bindings must contain at least one binding",
        ),
        (
            "revision-scheme",
            ValueError,
            "public_identity_scheme is unsupported",
        ),
        (
            "revision-digest",
            ValueError,
            "public_identity_sha256 does not match revision",
        ),
        (
            "evidence-binding-scheme",
            ValueError,
            "binding_digest_scheme is unsupported",
        ),
        (
            "builder-adapter-type",
            ValueError,
            "adapter_identity must be an AdapterReferenceBindingIdentity",
        ),
        (
            "builder-workflow-mismatch",
            ValueError,
            "adapter identity workflow does not match evidence",
        ),
        (
            "resolved-identity-mismatch",
            ValueError,
            "resolved reference identities do not match",
        ),
    ),
)
def test_reference_profile_integrity_guards_fail_closed(
    case: str,
    expected_exception: type[Exception],
    message: str,
) -> None:
    identity = _adapter_identity()
    workflow_binding = ReferenceProfileWorkflowBinding.from_adapter_identity(identity)
    revision = build_reference_profile_revision(
        revision_id=REVISION_ID,
        profile_id=PROFILE_ID,
        revision_number=1,
        display_name="GRCh38 primary",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-primary",
        workflow_bindings=(workflow_binding,),
        created_at=NOW,
    )
    summary = revision.public_summary(safe_key="grch38", enabled=True)
    evidence = build_reference_profile_revision_binding(
        profile_id=PROFILE_ID,
        revision_id=REVISION_ID,
        workflow_id="bulk-rnaseq",
        revision_public_identity_sha256=revision.public_identity_sha256,
        adapter_identity=identity,
        bound_at=NOW,
    )
    bound = BoundWorkflowReference(
        inputs=WorkflowInputs(config={}, samples=None, options={}),
        adapter=object(),
        identity=identity,
    )
    operations = {
        "adapter-digest-shape": lambda: replace(
            identity,
            identity_sha256="3" * 63,
        ),
        "adapter-scheme": lambda: replace(
            identity,
            identity_scheme="unsupported",
        ),
        "summary-workflow-order": lambda: replace(
            summary,
            compatible_workflow_ids=("encode-chipseq", "bulk-rnaseq"),
        ),
        "summary-enabled-type": lambda: replace(summary, enabled=1),
        "revision-empty-bindings": lambda: replace(
            revision,
            workflow_bindings=(),
        ),
        "revision-scheme": lambda: replace(
            revision,
            public_identity_scheme="unsupported",
        ),
        "revision-digest": lambda: replace(
            revision,
            public_identity_sha256="0" * 64,
        ),
        "evidence-binding-scheme": lambda: replace(
            evidence,
            binding_digest_scheme="unsupported",
        ),
        "builder-adapter-type": lambda: build_reference_profile_revision_binding(
            profile_id=PROFILE_ID,
            revision_id=REVISION_ID,
            workflow_id="bulk-rnaseq",
            revision_public_identity_sha256=revision.public_identity_sha256,
            adapter_identity=object(),
            bound_at=NOW,
        ),
        "builder-workflow-mismatch": lambda: build_reference_profile_revision_binding(
            profile_id=PROFILE_ID,
            revision_id=REVISION_ID,
            workflow_id="encode-chipseq",
            revision_public_identity_sha256=revision.public_identity_sha256,
            adapter_identity=identity,
            bound_at=NOW,
        ),
        "resolved-identity-mismatch": lambda: ResolvedReferenceProfile(
            summary=replace(summary, revision_id="refpr_" + "9" * 32),
            evidence=evidence,
            bound_reference=bound,
        ),
    }

    with pytest.raises(expected_exception, match=message):
        operations[case]()


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
