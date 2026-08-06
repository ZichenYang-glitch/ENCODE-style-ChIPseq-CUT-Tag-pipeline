"""Public request-model coverage for exact Project/SampleRevision selection."""

from __future__ import annotations

import pytest
from pydantic import ValidationError

from encode_pipeline.api.models import ValidationRequest
from encode_pipeline.api.routes.workflows import validate_workflow
from encode_pipeline.platform.results import Result
from encode_pipeline.services.data_registry import DataRegistryService
from api_test_client import ApiTestClient


PROJECT_ID = "prj_11111111111111111111111111111111"
REVISION_IDS = (
    "smpr_22222222222222222222222222222222",
    "smpr_33333333333333333333333333333333",
)
REFERENCE_PROFILE_REVISION_ID = "refpr_44444444444444444444444444444444"


def test_validation_request_accepts_ordered_exact_sample_revision_selection() -> None:
    request = ValidationRequest(
        config={},
        project_id=PROJECT_ID,
        sample_revision_ids=list(REVISION_IDS),
    )

    assert request.project_id == PROJECT_ID
    assert request.sample_revision_ids == list(REVISION_IDS)


def test_validation_request_keeps_omitted_selection_for_transitional_legacy() -> None:
    request = ValidationRequest(config={})

    assert request.project_id is None
    assert request.sample_revision_ids == []


@pytest.mark.parametrize(
    ("project_id", "revision_ids"),
    [
        (PROJECT_ID, []),
        (None, [REVISION_IDS[0]]),
        (PROJECT_ID, [REVISION_IDS[0], REVISION_IDS[0]]),
        ("prj_00000000000000000000000000000000", [REVISION_IDS[0]]),
        ("project-readable-name", [REVISION_IDS[0]]),
        (PROJECT_ID, ["sample-readable-name"]),
    ],
)
def test_validation_request_rejects_incomplete_or_non_exact_selection(
    project_id: str | None,
    revision_ids: list[str],
) -> None:
    with pytest.raises(ValidationError):
        ValidationRequest(
            config={},
            project_id=project_id,
            sample_revision_ids=revision_ids,
        )


def test_validate_route_passes_selection_outside_adapter_owned_inputs() -> None:
    observed: dict[str, object] = {}

    class RecordingService:
        def validate(
            self,
            workflow_id: str,
            inputs: object,
            *,
            project_sample_selection: object,
            reference_profile_revision_id: object,
        ) -> Result[None]:
            observed.update(
                workflow_id=workflow_id,
                inputs=inputs,
                selection=project_sample_selection,
                reference_profile_revision_id=reference_profile_revision_id,
            )
            return Result.success(None)

    request = ValidationRequest(
        config={"threads": 1},
        project_id=PROJECT_ID,
        sample_revision_ids=list(REVISION_IDS),
        reference_profile_revision_id=REFERENCE_PROFILE_REVISION_ID,
    )

    response = validate_workflow(
        "workflow-a",
        request,
        validation_service=RecordingService(),  # type: ignore[arg-type]
    )

    selection = observed["selection"]
    assert getattr(selection, "project_id") == PROJECT_ID
    assert getattr(selection, "sample_revision_ids") == REVISION_IDS
    assert observed["workflow_id"] == "workflow-a"
    assert observed["reference_profile_revision_id"] == REFERENCE_PROFILE_REVISION_ID
    assert response.ok is True  # type: ignore[union-attr]


def test_validate_api_freezes_and_projects_safe_bound_revision_ids(
    tmp_path,
    reference_ready_app,
) -> None:
    app = reference_ready_app
    project_ids = iter((PROJECT_ID, "prj_44444444444444444444444444444444"))
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository,
        project_id_factory=lambda: next(project_ids),
        sample_id_factory=lambda: "smp_55555555555555555555555555555555",
        sample_revision_id_factory=lambda: REVISION_IDS[0],
    )
    project = registry.create_project("Private cohort display")
    imported = registry.import_sample(
        project.project_id,
        stable_key="sample-a",
        display_name="Private sample display",
        attributes={"private_note": "must-not-leak"},
    )
    row = {
        "sample": "adapter-row",
        "fastq_1": str((tmp_path / "private.fastq.gz").resolve()),
        "layout": "SE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
    }

    with ApiTestClient(app) as client:
        response = client.post(
            "/api/v1/workflows/encode-style-chipseq-cuttag-atac-mnase/validate",
            json={
                "config": {},
                "samples": [row],
                "project_id": project.project_id,
                "sample_revision_ids": [imported.revision.sample_revision_id],
                "reference_profile_revision_id": (
                    app.state.test_reference_profile.revision_id
                ),
            },
        )

    assert response.status_code == 200
    body = response.json()
    assert body["snapshot"]["project_id"] == project.project_id
    assert body["snapshot"]["binding_mode"] == "bound_v1"
    assert body["snapshot"]["provenance"] == "resolved"
    assert body["snapshot"]["sample_revision_ids"] == [
        imported.revision.sample_revision_id
    ]
    assert len(body["snapshot"]["binding_digest"]) == 64
    assert "Private cohort display" not in response.text
    assert "Private sample display" not in response.text
    assert "private_note" not in response.text
    assert str(tmp_path) not in response.text


def test_validate_api_rejects_cross_project_selection_without_snapshot(
    tmp_path,
    reference_ready_app,
) -> None:
    app = reference_ready_app
    project_ids = iter((PROJECT_ID, "prj_44444444444444444444444444444444"))
    registry = DataRegistryService(
        repository=app.state.persistence.data_registry_repository,
        project_id_factory=lambda: next(project_ids),
        sample_id_factory=lambda: "smp_55555555555555555555555555555555",
        sample_revision_id_factory=lambda: REVISION_IDS[0],
    )
    owner = registry.create_project("Owner")
    imported = registry.import_sample(
        owner.project_id,
        stable_key="sample-a",
        display_name="Sample A",
        attributes={},
    )
    other = registry.create_project("Other")
    row = {
        "sample": "adapter-row",
        "fastq_1": str((tmp_path / "private.fastq.gz").resolve()),
        "layout": "SE",
        "assay": "chipseq",
        "target": "CTCF",
        "peak_mode": "narrow",
    }

    with ApiTestClient(app) as client:
        response = client.post(
            "/api/v1/workflows/encode-style-chipseq-cuttag-atac-mnase/validate",
            json={
                "config": {},
                "samples": [row],
                "project_id": other.project_id,
                "sample_revision_ids": [imported.revision.sample_revision_id],
                "reference_profile_revision_id": (
                    app.state.test_reference_profile.revision_id
                ),
            },
        )

    assert response.status_code == 409
    body = response.json()
    assert body["snapshot"] is None
    assert body["issues"][0]["code"] == "DATA_BINDING_SELECTION_INVALID"
    with app.state.persistence.engine.connect() as connection:
        assert (
            connection.exec_driver_sql(
                "SELECT count(*) FROM validated_input_snapshots"
            ).scalar_one()
            == 0
        )
