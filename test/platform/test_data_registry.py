"""Tests for workflow-neutral Project and Sample registry primitives."""

from __future__ import annotations

from dataclasses import FrozenInstanceError
from datetime import datetime, timedelta, timezone

import pytest

from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_DISPLAY_NAME,
    LEGACY_PROJECT_ID,
    PROJECT_SAMPLE_BINDING_DIGEST_SCHEME,
    SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
    BindingMode,
    BindingProvenance,
    Project,
    ProjectKind,
    ProjectSampleSelection,
    Sample,
    SampleRevision,
    SampleRevisionBindingRef,
    build_legacy_project,
    build_legacy_project_sample_binding,
    build_project_sample_binding,
    build_sample_revision_payload_digest,
    canonical_project_sample_binding,
    canonical_sample_revision_payload,
)


NOW = datetime(2026, 7, 26, 8, 0, tzinfo=timezone.utc)
PROJECT_ID = "prj_" + "1" * 32
OTHER_PROJECT_ID = "prj_" + "2" * 32
SAMPLE_ID = "smp_" + "3" * 32
REVISION_A_ID = "smpr_" + "4" * 32
REVISION_B_ID = "smpr_" + "5" * 32
WORKFLOW_INPUTS_DIGEST = "6" * 64


def _project(**changes: object) -> Project:
    values: dict[str, object] = {
        "project_id": PROJECT_ID,
        "display_name": "Pilot cohort",
        "kind": ProjectKind.USER,
        "created_at": NOW,
        "archived_at": None,
    }
    values.update(changes)
    return Project(**values)  # type: ignore[arg-type]


def _sample(**changes: object) -> Sample:
    values: dict[str, object] = {
        "sample_id": SAMPLE_ID,
        "project_id": PROJECT_ID,
        "stable_key": "donor-01",
        "created_at": NOW,
    }
    values.update(changes)
    return Sample(**values)  # type: ignore[arg-type]


def _revision(
    *,
    sample_revision_id: str = REVISION_A_ID,
    revision_number: int = 1,
    payload: dict[str, object] | None = None,
) -> SampleRevision:
    canonical_payload = canonical_sample_revision_payload(
        display_name="Donor 01",
        attributes=payload or {"organism": "human", "replicate": 1},
    )
    return SampleRevision(
        sample_revision_id=sample_revision_id,
        project_id=PROJECT_ID,
        sample_id=SAMPLE_ID,
        revision_number=revision_number,
        canonical_payload=canonical_payload,
        payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
        payload_digest=build_sample_revision_payload_digest(canonical_payload),
        created_at=NOW,
    )


def test_project_and_sample_records_normalize_aware_timestamps_to_utc() -> None:
    offset = timezone(timedelta(hours=8))

    project = _project(
        created_at=datetime(2026, 7, 26, 16, 0, tzinfo=offset),
        archived_at=datetime(2026, 7, 27, 16, 0, tzinfo=offset),
    )
    sample = _sample(created_at=datetime(2026, 7, 26, 16, 0, tzinfo=offset))
    revision = _revision()

    assert project.created_at == NOW
    assert project.created_at.tzinfo is timezone.utc
    assert project.archived_at == NOW + timedelta(days=1)
    assert sample.created_at == NOW
    assert sample.created_at.tzinfo is timezone.utc
    assert revision.created_at.tzinfo is timezone.utc


@pytest.mark.parametrize(
    ("factory", "changes"),
    [
        (_project, {"project_id": "project-user"}),
        (_project, {"display_name": " "}),
        (_project, {"created_at": datetime(2026, 7, 26, 8, 0)}),
        (_project, {"archived_at": NOW - timedelta(seconds=1)}),
        (_sample, {"sample_id": "sample-1"}),
        (_sample, {"project_id": "project-1"}),
        (_sample, {"stable_key": "../unsafe"}),
        (_sample, {"created_at": datetime(2026, 7, 26, 8, 0)}),
    ],
)
def test_project_and_sample_records_reject_invalid_identity_or_time(
    factory: object,
    changes: dict[str, object],
) -> None:
    with pytest.raises(ValueError):
        factory(**changes)  # type: ignore[operator]


def test_reserved_legacy_project_is_deterministic_and_immutable() -> None:
    first = build_legacy_project(NOW)
    second = build_legacy_project(NOW + timedelta(hours=1))

    assert first.project_id == second.project_id == LEGACY_PROJECT_ID
    assert first.display_name == LEGACY_PROJECT_DISPLAY_NAME
    assert first.kind is ProjectKind.SYSTEM
    assert first.archived_at is None
    with pytest.raises(ValueError, match="reserved Legacy Project"):
        _project(project_id=LEGACY_PROJECT_ID)
    with pytest.raises(ValueError, match="reserved Legacy Project"):
        Project(
            project_id=LEGACY_PROJECT_ID,
            display_name="Renamed legacy",
            kind=ProjectKind.SYSTEM,
            created_at=NOW,
        )
    with pytest.raises(ValueError, match="Legacy Project"):
        _sample(project_id=LEGACY_PROJECT_ID)


def test_sample_payload_is_canonical_json_and_framed_digest_is_stable() -> None:
    left = canonical_sample_revision_payload(
        display_name="Donor α",
        attributes={"z": ["α", True], "a": {"count": 2}},
    )
    right = canonical_sample_revision_payload(
        display_name="Donor α",
        attributes={"a": {"count": 2}, "z": ["α", True]},
    )

    assert (
        left
        == right
        == ('{"attributes":{"a":{"count":2},"z":["α",true]},"display_name":"Donor α"}')
    )
    assert build_sample_revision_payload_digest(left) == (
        build_sample_revision_payload_digest(right)
    )
    assert len(build_sample_revision_payload_digest(left)) == 64


@pytest.mark.parametrize(
    "payload",
    [
        {"bad": float("nan")},
        {"bad": 2**53},
        {"bad": object()},
        {1: "non-string key"},
    ],
)
def test_sample_payload_rejects_non_json_safe_values(payload: object) -> None:
    with pytest.raises(ValueError, match="JSON-safe"):
        canonical_sample_revision_payload(
            display_name="Sample",
            attributes=payload,  # type: ignore[arg-type]
        )


def test_sample_revision_is_immutable_and_verifies_canonical_payload_digest() -> None:
    revision = _revision()

    assert revision.to_payload() == {
        "attributes": {"organism": "human", "replicate": 1},
        "display_name": "Donor 01",
    }
    assert revision.display_name == "Donor 01"
    assert revision.attributes == {"organism": "human", "replicate": 1}
    with pytest.raises(FrozenInstanceError):
        revision.revision_number = 2  # type: ignore[misc]
    with pytest.raises(ValueError, match="canonical"):
        SampleRevision(
            sample_revision_id=REVISION_A_ID,
            project_id=PROJECT_ID,
            sample_id=SAMPLE_ID,
            revision_number=1,
            canonical_payload='{ "organism": "human" }',
            payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
            payload_digest="0" * 64,
            created_at=NOW,
        )
    with pytest.raises(ValueError, match="digest"):
        SampleRevision(
            sample_revision_id=REVISION_A_ID,
            project_id=PROJECT_ID,
            sample_id=SAMPLE_ID,
            revision_number=1,
            canonical_payload=(
                '{"attributes":{"organism":"human"},"display_name":"Donor 01"}'
            ),
            payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
            payload_digest="0" * 64,
            created_at=NOW,
        )


def test_project_sample_selection_preserves_order_and_rejects_duplicates() -> None:
    selection = ProjectSampleSelection(
        project_id=PROJECT_ID,
        sample_revision_ids=(REVISION_B_ID, REVISION_A_ID),
    )

    assert selection.sample_revision_ids == (REVISION_B_ID, REVISION_A_ID)
    with pytest.raises(ValueError, match="duplicate"):
        ProjectSampleSelection(
            project_id=PROJECT_ID,
            sample_revision_ids=(REVISION_A_ID, REVISION_A_ID),
        )
    with pytest.raises(ValueError, match="at least one"):
        ProjectSampleSelection(project_id=PROJECT_ID, sample_revision_ids=())
    with pytest.raises(ValueError, match="Legacy"):
        ProjectSampleSelection(
            project_id=LEGACY_PROJECT_ID,
            sample_revision_ids=(REVISION_A_ID,),
        )


def test_bound_project_sample_binding_is_ordered_and_stage2_scoped() -> None:
    revision_a = _revision()
    revision_b = _revision(
        sample_revision_id=REVISION_B_ID,
        revision_number=2,
        payload={"organism": "human", "replicate": 2},
    )
    refs = tuple(
        SampleRevisionBindingRef(
            sample_revision_id=revision.sample_revision_id,
            payload_digest=revision.payload_digest,
        )
        for revision in (revision_a, revision_b)
    )

    binding = build_project_sample_binding(
        project_id=PROJECT_ID,
        binding_mode=BindingMode.BOUND_V1,
        provenance=BindingProvenance.RESOLVED,
        workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
        sample_revisions=refs,
    )
    reversed_binding = build_project_sample_binding(
        project_id=PROJECT_ID,
        binding_mode=BindingMode.BOUND_V1,
        provenance=BindingProvenance.RESOLVED,
        workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
        sample_revisions=tuple(reversed(refs)),
    )

    assert (
        binding.digest_scheme
        == "sha256-framed-project-sample-binding-v1"
        == PROJECT_SAMPLE_BINDING_DIGEST_SCHEME
    )
    assert binding.sample_revision_ids == (REVISION_A_ID, REVISION_B_ID)
    assert binding.digest != reversed_binding.digest
    canonical = canonical_project_sample_binding(
        project_id=PROJECT_ID,
        binding_mode=BindingMode.BOUND_V1,
        provenance=BindingProvenance.RESOLVED,
        workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
        sample_revisions=refs,
    )
    assert "input_revisions" not in canonical


def test_legacy_binding_is_unresolved_and_cannot_claim_sample_revisions() -> None:
    binding = build_legacy_project_sample_binding(WORKFLOW_INPUTS_DIGEST)

    assert binding.project_id == LEGACY_PROJECT_ID
    assert binding.binding_mode is BindingMode.LEGACY_V1
    assert binding.provenance is BindingProvenance.UNRESOLVED
    assert binding.sample_revisions == ()
    with pytest.raises(ValueError):
        build_project_sample_binding(
            project_id=LEGACY_PROJECT_ID,
            binding_mode=BindingMode.LEGACY_V1,
            provenance=BindingProvenance.UNRESOLVED,
            workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
            sample_revisions=(
                SampleRevisionBindingRef(
                    sample_revision_id=REVISION_A_ID,
                    payload_digest="9" * 64,
                ),
            ),
        )
