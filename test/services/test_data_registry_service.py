"""Tests for Project/Sample registry orchestration and in-memory storage."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from itertools import count

import pytest

from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    BindingMode,
    BindingProvenance,
    ProjectKind,
    ProjectSampleSelection,
)
from encode_pipeline.services.data_registry import (
    DataRegistryService,
    SampleImportSpec,
)
from encode_pipeline.services.data_registry_repositories import (
    DataRegistryConflictError,
    InMemoryDataRegistryRepository,
)


NOW = datetime(2026, 7, 26, 8, 0, tzinfo=timezone.utc)
WORKFLOW_INPUTS_DIGEST = "a" * 64


def _ids(prefix: str):
    serial = count(1)
    return lambda: f"{prefix}_{next(serial):032x}"


def _service() -> DataRegistryService:
    return DataRegistryService(
        repository=InMemoryDataRegistryRepository(legacy_created_at=NOW),
        project_id_factory=_ids("prj"),
        sample_id_factory=_ids("smp"),
        sample_revision_id_factory=_ids("smpr"),
        now_factory=lambda: NOW,
    )


def test_project_create_list_archive_preserves_reserved_legacy_project() -> None:
    service = _service()

    created = service.create_project("  Pilot project  ")

    assert created.project_id == "prj_" + "1".zfill(32)
    assert created.display_name == "Pilot project"
    assert created.kind is ProjectKind.USER
    assert created.created_at == NOW
    assert service.list_projects() == (
        service.get_project(LEGACY_PROJECT_ID),
        created,
    )

    archived = service.archive_project(created.project_id)

    assert archived.archived_at == NOW
    assert service.list_projects() == (service.get_project(LEGACY_PROJECT_ID),)
    assert service.list_projects(include_archived=True) == (
        service.get_project(LEGACY_PROJECT_ID),
        archived,
    )
    with pytest.raises(ValueError, match="reserved Legacy Project"):
        service.archive_project(LEGACY_PROJECT_ID)


def test_sample_import_creates_stable_identity_and_initial_immutable_revision() -> None:
    service = _service()
    project = service.create_project("Pilot")

    result = service.import_sample(
        project.project_id,
        stable_key="donor-01",
        display_name="Donor 01",
        attributes={"organism": "human", "sex": "female"},
    )

    assert result.sample.sample_id == "smp_" + "1".zfill(32)
    assert result.sample.project_id == project.project_id
    assert result.sample.stable_key == "donor-01"
    assert result.revision.sample_revision_id == "smpr_" + "1".zfill(32)
    assert result.revision.project_id == project.project_id
    assert result.revision.sample_id == result.sample.sample_id
    assert result.revision.revision_number == 1
    assert result.revision.to_payload() == {
        "attributes": {"organism": "human", "sex": "female"},
        "display_name": "Donor 01",
    }
    assert service.list_samples(project.project_id)[0].sample == result.sample
    assert service.list_samples(project.project_id)[0].revision == result.revision
    assert service.list_samples(project.project_id)[0].display_name == "Donor 01"
    assert service.list_samples(project.project_id)[0].attributes == {
        "organism": "human",
        "sex": "female",
    }
    assert service.list_sample_revisions(result.sample.sample_id) == (result.revision,)


def test_sample_import_rejects_duplicate_key_legacy_or_archived_project() -> None:
    service = _service()
    project = service.create_project("Pilot")
    service.import_sample(
        project.project_id,
        stable_key="donor-01",
        display_name="Donor 01",
        attributes={},
    )

    with pytest.raises(DataRegistryConflictError, match="stable_key"):
        service.import_sample(
            project.project_id,
            stable_key="donor-01",
            display_name="Another name",
            attributes={},
        )
    with pytest.raises(ValueError, match="Legacy"):
        service.import_sample(
            LEGACY_PROJECT_ID,
            stable_key="legacy",
            display_name="Legacy",
            attributes={},
        )

    archived_project = service.create_project("Archived")
    service.archive_project(archived_project.project_id)
    with pytest.raises(ValueError, match="archived"):
        service.import_sample(
            archived_project.project_id,
            stable_key="new",
            display_name="New",
            attributes={},
        )


def test_batch_sample_import_is_atomic_when_a_later_row_conflicts() -> None:
    service = _service()
    project = service.create_project("Pilot")

    with pytest.raises(DataRegistryConflictError, match="stable_key"):
        service.import_samples(
            project.project_id,
            rows=(
                SampleImportSpec(
                    stable_key="duplicate",
                    display_name="First",
                    attributes={"row": 1},
                ),
                SampleImportSpec(
                    stable_key="duplicate",
                    display_name="Second",
                    attributes={"row": 2},
                ),
            ),
        )

    assert service.list_samples(project.project_id) == ()


def test_sample_revise_appends_without_overwriting_history() -> None:
    revision_ids = _ids("smpr")
    times = iter((NOW, NOW, NOW + timedelta(hours=1)))
    service = DataRegistryService(
        repository=InMemoryDataRegistryRepository(legacy_created_at=NOW),
        project_id_factory=_ids("prj"),
        sample_id_factory=_ids("smp"),
        sample_revision_id_factory=revision_ids,
        now_factory=lambda: next(times),
    )
    project = service.create_project("Pilot")
    imported = service.import_sample(
        project.project_id,
        stable_key="donor-01",
        display_name="Donor 01",
        attributes={"state": "initial"},
    )

    revised = service.revise_sample(
        imported.sample.sample_id,
        display_name="Donor 01 reviewed",
        attributes={"state": "reviewed"},
    )

    assert revised.revision_number == 2
    assert revised.sample_revision_id != imported.revision.sample_revision_id
    assert revised.to_payload() == {
        "attributes": {"state": "reviewed"},
        "display_name": "Donor 01 reviewed",
    }
    assert imported.revision.to_payload() == {
        "attributes": {"state": "initial"},
        "display_name": "Donor 01",
    }
    current = service.list_samples(project.project_id)
    assert current[0].display_name == "Donor 01 reviewed"
    assert current[0].revision == revised
    assert service.list_sample_revisions(imported.sample.sample_id) == (
        imported.revision,
        revised,
    )


def test_sample_revise_rejects_archived_project() -> None:
    service = _service()
    project = service.create_project("Pilot")
    imported = service.import_sample(
        project.project_id,
        stable_key="donor-01",
        display_name="Donor 01",
        attributes={},
    )
    service.archive_project(project.project_id)

    with pytest.raises(ValueError, match="archived"):
        service.revise_sample(
            imported.sample.sample_id,
            display_name="Donor 01",
            attributes={"new": True},
        )


def test_selection_resolution_is_exact_ordered_and_builds_binding() -> None:
    service = _service()
    project = service.create_project("Pilot")
    first = service.import_sample(
        project.project_id,
        stable_key="first",
        display_name="First",
        attributes={"sample": 1},
    )
    second = service.import_sample(
        project.project_id,
        stable_key="second",
        display_name="Second",
        attributes={"sample": 2},
    )
    selection = ProjectSampleSelection(
        project_id=project.project_id,
        sample_revision_ids=(
            second.revision.sample_revision_id,
            first.revision.sample_revision_id,
        ),
    )

    binding = service.resolve_project_sample_selection(
        selection,
        workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
    )

    assert binding.project_id == project.project_id
    assert binding.binding_mode is BindingMode.BOUND_V1
    assert binding.provenance is BindingProvenance.RESOLVED
    assert binding.sample_revision_ids == selection.sample_revision_ids
    assert binding.input_revisions == ()


def test_selection_resolution_rejects_cross_project_and_archived_project() -> None:
    service = _service()
    first_project = service.create_project("First")
    second_project = service.create_project("Second")
    imported = service.import_sample(
        first_project.project_id,
        stable_key="sample",
        display_name="Sample",
        attributes={},
    )
    wrong_project_selection = ProjectSampleSelection(
        project_id=second_project.project_id,
        sample_revision_ids=(imported.revision.sample_revision_id,),
    )

    with pytest.raises(ValueError, match="same Project"):
        service.resolve_project_sample_selection(
            wrong_project_selection,
            workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
        )

    service.archive_project(first_project.project_id)
    archived_selection = ProjectSampleSelection(
        project_id=first_project.project_id,
        sample_revision_ids=(imported.revision.sample_revision_id,),
    )
    with pytest.raises(ValueError, match="archived"):
        service.resolve_project_sample_selection(
            archived_selection,
            workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST,
        )


def test_legacy_binding_is_explicit_and_never_resolves_samples() -> None:
    service = _service()

    binding = service.build_legacy_binding(
        workflow_inputs_digest=WORKFLOW_INPUTS_DIGEST
    )

    assert binding.project_id == LEGACY_PROJECT_ID
    assert binding.binding_mode is BindingMode.LEGACY_V1
    assert binding.provenance is BindingProvenance.UNRESOLVED
    assert binding.sample_revision_ids == ()
