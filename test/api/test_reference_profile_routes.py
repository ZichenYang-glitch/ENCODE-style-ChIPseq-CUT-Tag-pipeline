"""Public API contract tests for Reference Profile selection and projection."""

from __future__ import annotations

from datetime import datetime, timedelta, timezone
from pathlib import Path

import pytest

from api_test_client import ApiTestClient
from encode_pipeline.api.main import create_app
from encode_pipeline.platform.adapters import (
    CommandSpec,
    DagPreview,
    WorkflowAvailability,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
    WorkflowSchema,
    WorkspacePlan,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.reference_profiles import (
    AdapterReferenceBindingIdentity,
    BoundWorkflowReference,
)
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Result
from encode_pipeline.services.private_reference_profiles import (
    PrivateReferenceProfileConfig,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService
from encode_pipeline.services.reference_profile_runtime import (
    ReferenceProfileBindingService,
)
from encode_pipeline.services.runs import RunService
from encode_pipeline.services.validated_inputs import (
    ValidatedInputService,
    ValidatedRunCreationService,
)
from encode_pipeline.services.validation import ValidationService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider


NOW = datetime(2026, 8, 3, 11, 0, tzinfo=timezone.utc)
WORKFLOW_ID = "reference-api-workflow"
PROFILE_ID = "refp_" + "a" * 32
REVISION_ID = "refpr_" + "b" * 32
PRIVATE_REFERENCE_PATH = "/operator/private/references/mm10.fa"


class _ApiReferenceAdapter:
    metadata = WorkflowMetadata(
        workflow_id=WORKFLOW_ID,
        name="Reference API workflow",
        version="1.0.0",
    )
    capabilities = WorkflowCapabilities(
        supports=("validation", "workspace_plan", "command")
    )

    def schema(self) -> WorkflowSchema:
        return WorkflowSchema(schema_version="1.0.0")

    def validate(self, inputs: WorkflowInputs) -> Result[object]:
        if inputs.config.get("resolved_reference_path") != PRIVATE_REFERENCE_PATH:
            raise AssertionError("API validation bypassed exact reference resolution")
        return Result.success({"accepted": True})

    def preview_dag(self, inputs: WorkflowInputs) -> Result[DagPreview]:
        return Result.success(DagPreview())

    def plan_workspace(
        self,
        inputs: WorkflowInputs,
        workspace: str | Path,
    ) -> Result[WorkspacePlan]:
        return Result.success(WorkspacePlan(directories=[str(workspace)]))

    def build_command(
        self,
        plan: WorkspacePlan,
        workspace: str | Path,
    ) -> Result[CommandSpec]:
        return Result.success(CommandSpec(argv=["reference-api-test"]))

    def extract_artifacts(self, inputs, workspace):
        return Result.success(())

    def execution_availability(self) -> WorkflowAvailability:
        return WorkflowAvailability(
            execution="available",
            reason_code="WORKFLOW_EXECUTION_READY",
        )

    def capture_build_identity(self) -> Result[WorkflowBuildIdentity]:
        return Result.success(
            WorkflowBuildIdentity(
                workflow_id=WORKFLOW_ID,
                adapter_version="1.0.0",
                scheme="sha256-tree-v1",
                logical_entrypoint="tests/reference-api-workflow",
                digest="d" * 64,
                captured_at=NOW,
            )
        )

    def verify_reference_profile_binding(self, payload):
        return Result.success(
            AdapterReferenceBindingIdentity(
                workflow_id=WORKFLOW_ID,
                contract_version="reference-api-binding-v1",
                identity_sha256=str(payload["asset_sha256"]),
            )
        )

    def bind_reference_profile(self, inputs: WorkflowInputs, payload):
        identity = self.verify_reference_profile_binding(payload).value
        assert identity is not None
        return Result.success(
            BoundWorkflowReference(
                inputs=WorkflowInputs(
                    config={
                        **inputs.config,
                        "resolved_reference_path": str(payload["private_path"]),
                    },
                    samples=inputs.samples,
                    options=inputs.options,
                ),
                adapter=self,
                identity=identity,
            )
        )


class _ApiHarness:
    def __init__(self, tmp_path: Path) -> None:
        self.app = create_app(
            database_url=f"sqlite:///{tmp_path / 'reference-api.db'}",
            workspace_root=(tmp_path / "workspaces").resolve(),
            project_root=tmp_path.resolve(),
        )
        self.adapter = _ApiReferenceAdapter()
        self.registry = WorkflowRegistry([self.adapter])
        self.builds = WorkflowBuildIdentityProvider(
            self.registry,
            project_root=tmp_path.resolve(),
        )
        config = PrivateReferenceProfileConfig(
            {
                "mm10-private": {
                    WORKFLOW_ID: {
                        "private_path": PRIVATE_REFERENCE_PATH,
                        "asset_sha256": "e" * 64,
                    }
                }
            }
        )
        self.profiles = ReferenceProfileService(
            repository=self.app.state.persistence.reference_profile_repository,
            private_config_provider=lambda: config,
            adapter_provider=self.registry.get,
            profile_id_factory=lambda: PROFILE_ID,
            revision_id_factory=lambda: REVISION_ID,
            now_factory=lambda: NOW,
        )
        self.bindings = ReferenceProfileBindingService(
            repository=self.app.state.persistence.reference_profile_repository,
            private_config_provider=lambda: config,
            adapter_provider=self.registry.get,
            now_factory=lambda: NOW,
        )
        self.validation = ValidatedInputService(
            registry=self.registry,
            validation_service=ValidationService(self.registry),
            build_identity_provider=self.builds,
            repository=self.app.state.persistence.repository,
            reference_profile_binding_service=self.bindings,
            reference_profile_catalog=self.profiles,
            snapshot_id_factory=lambda: "vsnap_" + "c" * 32,
            clock=lambda: NOW,
            snapshot_ttl=timedelta(hours=1),
        )
        self.runs = RunService(
            self.registry,
            id_factory=lambda: "run-reference-api",
            repository=self.app.state.persistence.repository,
        )
        self.creation = ValidatedRunCreationService(
            run_service=self.runs,
            build_identity_provider=self.builds,
            reference_profile_binding_service=self.bindings,
            clock=lambda: NOW + timedelta(minutes=1),
        )
        self.app.state.registry = self.registry
        self.app.state.validation_service = ValidationService(self.registry)
        self.app.state.validated_input_service = self.validation
        self.app.state.reference_profile_service = self.profiles
        self.app.state.reference_profile_binding_service = self.bindings
        self.app.state.validated_run_creation_service = self.creation
        self.app.state.run_service = self.runs

    def register_and_enable(self):
        summary = self.profiles.register(
            safe_key="mouse-reference",
            display_name="mm10 primary",
            organism="Mus musculus",
            assembly="mm10",
            config_key="mm10-private",
        )
        self.profiles.enable(PROFILE_ID, revision_id=summary.revision_id)
        return summary

    def close(self) -> None:
        self.app.state.run_queue.close()
        self.app.state.persistence.close()


@pytest.fixture
def harness(tmp_path: Path):
    value = _ApiHarness(tmp_path)
    try:
        yield value
    finally:
        value.close()


def test_enabled_compatible_list_is_exact_path_free_and_handles_empty_and_unknown(
    harness: _ApiHarness,
) -> None:
    with ApiTestClient(harness.app) as client:
        empty = client.get(f"/api/v1/workflows/{WORKFLOW_ID}/reference-profiles")
        assert empty.status_code == 200
        assert empty.json() == {
            "ok": True,
            "workflow_id": WORKFLOW_ID,
            "profiles": [],
            "issues": [],
        }

        unknown = client.get("/api/v1/workflows/missing/reference-profiles")
        assert unknown.status_code == 404
        assert unknown.json()["profiles"] == []
        assert unknown.json()["issues"][0]["code"] == "WORKFLOW_NOT_FOUND"

        summary = harness.register_and_enable()
        response = client.get(f"/api/v1/workflows/{WORKFLOW_ID}/reference-profiles")

        assert response.status_code == 200
        body = response.json()
        assert body["ok"] is True
        assert len(body["profiles"]) == 1
        profile = body["profiles"][0]
        assert profile == {
            "profile_id": PROFILE_ID,
            "revision_id": REVISION_ID,
            "revision_number": 1,
            "display_name": "mm10 primary",
            "organism": "Mus musculus",
            "assembly": "mm10",
            "identity_sha256": summary.public_identity_sha256,
        }
        rendered = response.text
        for private_value in (
            PRIVATE_REFERENCE_PATH,
            "mm10-private",
            "mouse-reference",
            "config_key",
            "private_path",
        ):
            assert private_value not in rendered


def test_reference_profile_http_surface_is_read_only(harness: _ApiHarness) -> None:
    path = "/api/v1/workflows/{workflow_id}/reference-profiles"
    operations = harness.app.openapi()["paths"][path]

    assert set(operations) == {"get"}


def test_reference_profile_catalog_failure_is_stable_and_redacted(
    harness: _ApiHarness,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    def fail_closed(_workflow_id: str):
        raise RuntimeError(PRIVATE_REFERENCE_PATH)

    monkeypatch.setattr(harness.profiles, "list_enabled", fail_closed)
    with ApiTestClient(harness.app) as client:
        response = client.get(f"/api/v1/workflows/{WORKFLOW_ID}/reference-profiles")

    assert response.status_code == 503
    assert response.json() == {
        "ok": False,
        "workflow_id": WORKFLOW_ID,
        "profiles": [],
        "issues": [
            {
                "code": "REFERENCE_PROFILE_UNAVAILABLE",
                "message": "Reference Profiles are unavailable.",
                "severity": "error",
                "path": "reference_profile_revision_id",
                "source": "reference_profile",
                "technical_message": None,
                "hint": None,
                "context": {},
            }
        ],
    }
    assert PRIVATE_REFERENCE_PATH not in response.text


def test_validation_pins_exact_revision_and_run_detail_keeps_historical_summary(
    harness: _ApiHarness,
) -> None:
    summary = harness.register_and_enable()
    with ApiTestClient(harness.app) as client:
        validation = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/validate",
            json={
                "config": {"user_setting": "visible"},
                "samples": [{"sample": "S1"}],
                "options": {},
                "reference_profile_revision_id": summary.revision_id,
            },
        )

        assert validation.status_code == 200
        validation_body = validation.json()
        assert validation_body["ok"] is True
        snapshot = validation_body["snapshot"]
        assert snapshot["reference_profile"] == {
            "profile_id": PROFILE_ID,
            "revision_id": REVISION_ID,
            "revision_number": 1,
            "display_name": "mm10 primary",
            "organism": "Mus musculus",
            "assembly": "mm10",
            "identity_sha256": summary.public_identity_sha256,
        }
        frozen = harness.validation.get_validated_reference_binding(
            snapshot["snapshot_id"]
        )
        assert frozen is not None
        assert frozen.revision_id == summary.revision_id
        assert frozen.revision_public_identity_sha256 == (
            summary.public_identity_sha256
        )
        assert PRIVATE_REFERENCE_PATH not in validation.text

        created = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json={"snapshot_id": snapshot["snapshot_id"], "tags": {}},
        )
        assert created.status_code == 201
        run_id = created.json()["run"]["run_id"]
        harness.profiles.disable(PROFILE_ID)

        detail = client.get(f"/api/v1/runs/{run_id}")

        assert detail.status_code == 200
        assert (
            detail.json()["run"]["reference_profile"] == (snapshot["reference_profile"])
        )
        assert harness.runs.get_run_reference_binding(run_id) == frozen
        for private_value in (
            PRIVATE_REFERENCE_PATH,
            "mm10-private",
            "mouse-reference",
            "config_key",
            "private_path",
        ):
            assert private_value not in detail.text


def test_validation_route_forwards_the_selected_revision_without_fallback(
    harness: _ApiHarness,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    summary = harness.register_and_enable()
    observed: list[str | None] = []
    original = harness.validation.validate

    def record_selection(workflow_id, inputs, **kwargs):
        observed.append(kwargs.get("reference_profile_revision_id"))
        return original(workflow_id, inputs, **kwargs)

    monkeypatch.setattr(harness.validation, "validate", record_selection)
    with ApiTestClient(harness.app) as client:
        response = client.post(
            f"/api/v1/workflows/{WORKFLOW_ID}/validate",
            json={
                "config": {"user_setting": "visible"},
                "reference_profile_revision_id": summary.revision_id,
            },
        )

    assert response.status_code == 200
    assert observed == [summary.revision_id]
