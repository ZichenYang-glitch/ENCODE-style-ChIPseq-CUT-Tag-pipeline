"""Public request boundary for opaque InputFile revision selection."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
import json

import pytest
from pydantic import ValidationError

from api_test_client import ApiTestClient
from encode_pipeline.api.main import create_app
from encode_pipeline.api.models import ValidatedInputSnapshotResponse, ValidationRequest
from encode_pipeline.api.routes.workflows import validate_workflow
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.data_registry import (
    BindingMode,
    BindingProvenance,
    SampleRevisionBindingRef,
    build_project_sample_binding,
)
from encode_pipeline.platform.input_registry import (
    InputFileRevisionBindingRef,
    InputProvenanceMode,
    PlannedInputUse,
    build_input_use_binding,
    build_input_use_binding_envelope,
)
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.snapshots import (
    PAYLOAD_DIGEST_SCHEME,
    VALIDATION_EVIDENCE_OUTCOME,
    ValidatedInputSnapshot,
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)


PROJECT_ID = "prj_11111111111111111111111111111111"
SAMPLE_REVISION_ID = "smpr_22222222222222222222222222222222"
INPUT_FILE_REVISION_ID = "inpfr_33333333333333333333333333333333"
NOW = datetime(2026, 7, 26, 10, 0, tzinfo=timezone.utc)


def _selection() -> dict[str, object]:
    return {
        "input_use_key": "primary_reads",
        "occurrence": 0,
        "input_file_revision_ids": [INPUT_FILE_REVISION_ID],
    }


def test_validation_request_accepts_only_opaque_input_revision_selection() -> None:
    request = ValidationRequest(
        config={},
        project_id=PROJECT_ID,
        sample_revision_ids=[SAMPLE_REVISION_ID],
        input_selections=[_selection()],
    )

    assert request.input_selections[0].input_use_key == "primary_reads"
    assert request.input_selections[0].occurrence == 0
    assert request.input_selections[0].input_file_revision_ids == [
        INPUT_FILE_REVISION_ID
    ]
    dumped = request.model_dump(mode="json")
    assert "path" not in str(dumped).lower()
    assert "config_key" not in str(dumped)


def test_input_selection_requires_trusted_project_sample_binding() -> None:
    with pytest.raises(ValidationError):
        ValidationRequest(config={}, input_selections=[_selection()])


def test_input_selection_coordinates_must_be_unique() -> None:
    with pytest.raises(ValidationError):
        ValidationRequest(
            config={},
            project_id=PROJECT_ID,
            sample_revision_ids=[SAMPLE_REVISION_ID],
            input_selections=[_selection(), _selection()],
        )


def test_input_selection_revision_ids_must_be_unique() -> None:
    selection = _selection()
    selection["input_file_revision_ids"] = [
        INPUT_FILE_REVISION_ID,
        INPUT_FILE_REVISION_ID,
    ]

    with pytest.raises(ValidationError):
        ValidationRequest(
            config={},
            project_id=PROJECT_ID,
            sample_revision_ids=[SAMPLE_REVISION_ID],
            input_selections=[selection],
        )


def test_duplicate_input_revision_ids_return_request_error_not_server_error(
    tmp_path,
) -> None:
    app = create_app(
        database_url=f"sqlite:///{tmp_path / 'platform.db'}",
        workspace_root=tmp_path / "workspaces",
    )
    selection = _selection()
    selection["input_file_revision_ids"] = [
        INPUT_FILE_REVISION_ID,
        INPUT_FILE_REVISION_ID,
    ]

    with ApiTestClient(app) as client:
        response = client.post(
            "/api/v1/workflows/encode-style-chipseq-cuttag-atac-mnase/validate",
            json={
                "config": {},
                "project_id": PROJECT_ID,
                "sample_revision_ids": [SAMPLE_REVISION_ID],
                "input_selections": [selection],
            },
        )

    assert response.status_code == 400
    assert response.json()["issues"][0]["code"] == "API_REQUEST_INVALID"


def test_validation_request_rejects_more_than_public_input_use_limit() -> None:
    selections = [
        {
            **_selection(),
            "input_use_key": f"use_{index}",
        }
        for index in range(257)
    ]

    with pytest.raises(ValidationError):
        ValidationRequest(
            config={},
            project_id=PROJECT_ID,
            sample_revision_ids=[SAMPLE_REVISION_ID],
            input_selections=selections,
        )


def test_validate_route_passes_opaque_selection_outside_workflow_inputs() -> None:
    observed: dict[str, object] = {}

    class RecordingService:
        def validate(
            self,
            workflow_id: str,
            inputs: object,
            *,
            project_sample_selection: object,
            input_file_revision_selections: object,
        ) -> Result[None]:
            observed.update(
                workflow_id=workflow_id,
                inputs=inputs,
                project_sample_selection=project_sample_selection,
                input_file_revision_selections=input_file_revision_selections,
            )
            return Result.success(None)

    response = validate_workflow(
        "workflow-a",
        ValidationRequest(
            config={"threads": 1},
            project_id=PROJECT_ID,
            sample_revision_ids=[SAMPLE_REVISION_ID],
            input_selections=[_selection()],
        ),
        validation_service=RecordingService(),  # type: ignore[arg-type]
    )

    selections = observed["input_file_revision_selections"]
    assert len(selections) == 1  # type: ignore[arg-type]
    selection = selections[0]  # type: ignore[index]
    assert getattr(selection, "input_use_key") == "primary_reads"
    assert getattr(selection, "occurrence") == 0
    assert getattr(selection, "input_file_revision_ids") == (INPUT_FILE_REVISION_ID,)
    assert response.ok is True  # type: ignore[union-attr]


def test_input_use_capability_unavailable_maps_to_stable_conflict() -> None:
    class UnavailableService:
        def validate(self, workflow_id: str, inputs: object, **kwargs) -> Result[None]:
            return Result.failure(
                [
                    Issue(
                        code="INPUT_USE_CAPABILITY_UNAVAILABLE",
                        message="Managed provenance is not qualified.",
                        source="input_registry",
                        path="input_selections",
                    )
                ]
            )

    response = validate_workflow(
        "workflow-a",
        ValidationRequest(config={}),
        validation_service=UnavailableService(),  # type: ignore[arg-type]
    )

    assert response.status_code == 409  # type: ignore[union-attr]
    body = json.loads(response.body)  # type: ignore[union-attr]
    assert body["issues"][0]["code"] == "INPUT_USE_CAPABILITY_UNAVAILABLE"


def test_snapshot_response_allowlists_path_free_input_provenance() -> None:
    canonical_payload = canonical_workflow_inputs_json(WorkflowInputs(config={}))
    payload_digest = build_workflow_inputs_digest(canonical_payload)
    snapshot = ValidatedInputSnapshot(
        snapshot_id="vsnap_0123456789abcdef0123456789abcdef",
        workflow_id="workflow-a",
        adapter_version="1.0.0",
        schema_version="1.0.0",
        schema_dialect="https://json-schema.org/draft/2020-12/schema",
        workflow_build_identity=WorkflowBuildIdentity(
            workflow_id="workflow-a",
            adapter_version="1.0.0",
            scheme="sha256-tree-v1",
            logical_entrypoint="workflow/Snakefile",
            digest="a" * 64,
            captured_at=NOW,
        ),
        canonical_payload=canonical_payload,
        payload_digest_scheme=PAYLOAD_DIGEST_SCHEME,
        payload_digest=payload_digest,
        validation_outcome=VALIDATION_EVIDENCE_OUTCOME,
        validation_issue_codes=(),
        validated_at=NOW,
        expires_at=NOW + timedelta(minutes=30),
    )
    project_samples = build_project_sample_binding(
        project_id=PROJECT_ID,
        binding_mode=BindingMode.BOUND_V1,
        provenance=BindingProvenance.RESOLVED,
        workflow_inputs_digest=payload_digest,
        sample_revisions=(
            SampleRevisionBindingRef(
                sample_revision_id=SAMPLE_REVISION_ID,
                payload_digest="b" * 64,
            ),
        ),
    )
    planned = PlannedInputUse(
        key="primary_reads",
        occurrence=0,
        capability_version="regular-file-v1",
        closure_contract_version="regular_file_v1",
        provenance_mode=InputProvenanceMode.MANAGED_REVISION_V1,
        input_file_revision_ids=(INPUT_FILE_REVISION_ID,),
    )
    input_use = build_input_use_binding(
        planned,
        members=(
            InputFileRevisionBindingRef(
                logical_member_key="file",
                input_file_id="inpf_44444444444444444444444444444444",
                input_file_revision_id=INPUT_FILE_REVISION_ID,
                revision_digest="c" * 64,
                size_bytes=123,
                content_sha256="d" * 64,
            ),
        ),
    )
    input_binding = build_input_use_binding_envelope(
        project_id=PROJECT_ID,
        project_sample_binding_digest=project_samples.digest,
        workflow_id="workflow-a",
        adapter_contract_version="workflow-a-inputs-v1",
        workflow_inputs_digest=payload_digest,
        input_uses=(input_use,),
    )

    dumped = ValidatedInputSnapshotResponse.from_snapshot(
        snapshot,
        project_samples,
        input_binding,
    ).model_dump(mode="json")

    assert dumped["input_binding"] == {
        "mode": "declared_input_uses_v1",
        "adapter_contract_version": "workflow-a-inputs-v1",
        "digest": input_binding.digest,
        "fully_managed": True,
        "input_uses": [
            {
                "input_use_key": "primary_reads",
                "occurrence": 0,
                "capability_version": "regular-file-v1",
                "closure_contract_version": "regular_file_v1",
                "provenance_mode": "managed_revision_v1",
                "input_file_revision_ids": [INPUT_FILE_REVISION_ID],
            }
        ],
    }
    serialized = str(dumped).lower()
    for private_name in (
        "relative_path",
        "config_key",
        "content_sha256",
        "checksum",
        "input_file_id",
        "logical_member_key",
        "revision_digest",
        "size_bytes",
        "closure_digest",
        "workflow_inputs_digest",
        "project_sample_binding_digest",
    ):
        assert private_name not in serialized


@pytest.mark.parametrize(
    "private_field",
    (
        {"physical_path": "/srv/private/reads.fastq.gz"},
        {"pool_relative_path": "reads/private.fastq.gz"},
        {"storage_pool_config_key": "operator-primary"},
        {"provenance_mode": "managed_revision_v1"},
        {"capability_version": "managed-regular-file-v1"},
        {"closure_digest": "a" * 64},
    ),
)
def test_input_selection_rejects_private_or_server_owned_fields(
    private_field: dict[str, str],
) -> None:
    with pytest.raises(ValidationError):
        ValidationRequest(
            config={},
            project_id=PROJECT_ID,
            sample_revision_ids=[SAMPLE_REVISION_ID],
            input_selections=[{**_selection(), **private_field}],
        )
