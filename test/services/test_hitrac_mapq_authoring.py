"""Hi-TrAC mapq suggestions through the real authoring service.

The bridge returns a fixed proposal. It is not a model and does not stand in
for an online provider. Validation uses the registered Hi-TrAC adapter.
"""

from __future__ import annotations

from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
from types import SimpleNamespace

import pytest

from encode_pipeline.adapters.bulk_rnaseq.adapter import WORKFLOW_ID as BULK_WORKFLOW
from encode_pipeline.adapters.encode_authoring import build_encode_authoring_schema
from encode_pipeline.adapters.encode_authoring_projection import (
    explainable_issues,
    project_encode_fields,
)
from encode_pipeline.adapters.hitrac_preprocess.adapter import (
    HiTracPreprocessAdapter,
)
from encode_pipeline.adapters.hitrac_preprocess.admission import INDEX_SUFFIXES
from encode_pipeline.adapters.hitrac_preprocess.execution import (
    WORKFLOW_ID as HITRAC_WORKFLOW,
)
from encode_pipeline.adapters.hitrac_preprocess.results import (
    HiTracPreprocessResultsAdapter,
)
from encode_pipeline.persistence import (
    create_database_engine,
    create_session_factory,
    upgrade_database,
)
from encode_pipeline.persistence.authentication import (
    SqlAlchemyAuthenticationRepository,
)
from encode_pipeline.persistence.authoring import SqlAlchemyAuthoringRepository
from encode_pipeline.persistence.data_registry import SqlAlchemyDataRegistryRepository
from encode_pipeline.platform.adapters import WorkflowInputs
from encode_pipeline.platform.authoring import (
    CAPABILITY_UNAVAILABLE,
    INVALID_REQUEST,
    MODEL_OUTPUT_REJECTED,
    SAMPLE_SELECTION_REQUIRED,
    AuthoringActor,
    AuthoringRejected,
    projection_identity_digest,
)
from encode_pipeline.platform.registry import ENCODE_WORKFLOW_ID
from encode_pipeline.services.authentication import (
    Argon2idParameters,
    PasswordManager,
    digest_session_token,
)
from encode_pipeline.services.authentication_service import (
    AccountAdministrationService,
    AuthenticationService,
)
from encode_pipeline.services.authoring import AuthoringService
from encode_pipeline.services.authoring_bridge import BridgeAdmission, BridgeStatus
from encode_pipeline.services.data_registry import DataRegistryService
from encode_pipeline.services.defaults import (
    create_default_validation_service,
    create_default_workflow_registry,
)
from encode_pipeline.services.private_reference_profiles import (
    PrivateReferenceProfileConfig,
)
from encode_pipeline.services.reference_profile_runtime import (
    ReferenceProfileBindingService,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService
from encode_pipeline.services.run_repositories import InMemoryRunRepository
from encode_pipeline.services.validated_inputs import ValidatedInputService
from encode_pipeline.services.workflow_builds import WorkflowBuildIdentityProvider

NOW = datetime(2026, 11, 18, 12, 0, tzinfo=timezone.utc)
ENCODE_PROJECTION_DIGEST = (
    "b3ea5bb01c2457fc788508187d967445c1eca0f1c95f26d0a6560b6a497ac71b"
)
SAMPLE_REVISION = "smpr_" + ("ab" * 16)
SAMPLE = {
    "sample_id": "sample 1",
    "fastq_1": "/data/reads_1.fastq.gz",
    "fastq_2": "/data/reads_2.fastq.gz",
}


class Clock:
    def __init__(self) -> None:
        self.now = NOW

    def __call__(self) -> datetime:
        return self.now


class ControlledBridge:
    """Returns one authored proposal. This is not a model client."""

    def __init__(self) -> None:
        self.epoch = "owner-epoch-mapq"
        self.value = 20
        self.field_token = "options.mapq"

    def status(self) -> BridgeStatus:
        return BridgeStatus(available=True)

    def owner_epoch(self) -> str:
        return self.epoch

    def admit(self, request) -> BridgeAdmission:
        del request
        return BridgeAdmission(instance_id="instance01", state="ready")

    def execute_turn(self, request) -> str:
        body = {
            "base_revision": request.base_revision,
            "changes": [{"field_token": self.field_token, "value": self.value}],
            "kind": "proposal",
            "reason_codes": [],
            "schema_digest": request.schema_digest,
        }
        return json.dumps(body, separators=(",", ":"), sort_keys=True)

    def close(self, request) -> None:
        del request


def _digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _reference_payload(root: Path) -> dict[str, str]:
    fasta = root / "reference.fa"
    fasta.write_text(">chrA\nACGTACGT\n")
    prefix = root / "reference"
    files = {"fasta": _digest(fasta)}
    for suffix in INDEX_SUFFIXES:
        index = Path(f"{prefix}{suffix}")
        index.write_text(f"synthetic-index:{suffix}\n")
        files[suffix] = _digest(index)
    binding = root / "reference-binding.json"
    binding.write_text(
        json.dumps(
            {
                "schema_version": "hitrac-reference-binding-v1",
                "fasta": str(fasta),
                "prefix": str(prefix),
                "files": files,
                "contigs": {"chrA": 8},
            }
        )
    )
    return {
        "schema_version": "hitrac-reference-profile-v1",
        "binding": str(binding),
        "sha256": _digest(binding),
    }


def _bound_case_fixture():
    """Reuse the admitted-runtime stub from the Hi-TrAC adapter contract tests."""

    path = Path(__file__).resolve().parents[1] / "adapters" / "test_hitrac_adapter.py"
    spec = importlib.util.spec_from_file_location("hitrac_adapter_bound_case", path)
    if spec is None or spec.loader is None:
        raise RuntimeError("bound_case fixture is missing")
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module.bound_case


bound_case = _bound_case_fixture()


def _authoring_stack(tmp_path: Path, environ: dict[str, str]) -> SimpleNamespace:
    tmp_path.mkdir(parents=True, exist_ok=True)
    database_url = f"sqlite:///{tmp_path / 'authoring.db'}"
    upgrade_database(database_url)
    engine = create_database_engine(database_url)
    clock = Clock()
    auth_repository = SqlAlchemyAuthenticationRepository(create_session_factory(engine))
    passwords = PasswordManager(
        Argon2idParameters(memory_cost_kib=19456, time_cost=2, parallelism=1)
    )
    AccountAdministrationService(
        repository=auth_repository,
        password_manager=passwords,
        now_factory=clock,
    ).bootstrap_initial_administrator("owner", "correct horse battery staple")
    login = AuthenticationService(
        repository=auth_repository,
        password_manager=passwords,
        now_factory=clock,
    ).login("owner", "correct horse battery staple", client_identity="127.0.0.1")
    actor = AuthoringActor(
        user_id=login.principal.user_id,
        session_digest=digest_session_token(login.secrets.session_token),
    )
    data = DataRegistryService(
        repository=SqlAlchemyDataRegistryRepository(create_session_factory(engine)),
        now_factory=clock,
    )
    project = data.create_project("Mapq")
    registry = create_default_workflow_registry(environ=environ)
    payload = _reference_payload(tmp_path)
    private = PrivateReferenceProfileConfig({"tiny": {HITRAC_WORKFLOW: payload}})
    catalog = ReferenceProfileService(
        private_config_provider=lambda: private,
        adapter_provider=registry.get,
    )
    reference = catalog.register(
        safe_key="tiny",
        display_name="Tiny synthetic contract",
        organism="synthetic",
        assembly="tiny",
        config_key="tiny",
    )
    catalog.enable(reference.profile_id, revision_id=reference.revision_id)
    binder = ReferenceProfileBindingService(
        repository=catalog.repository,
        private_config_provider=lambda: private,
        adapter_provider=registry.get,
    )
    validated = ValidatedInputService(
        registry=registry,
        validation_service=create_default_validation_service(registry=registry),
        build_identity_provider=WorkflowBuildIdentityProvider(registry),
        repository=InMemoryRunRepository(),
        reference_profile_binding_service=binder,
        clock=clock,
    )
    service = AuthoringService(
        repository=SqlAlchemyAuthoringRepository(
            create_session_factory(engine),
            clock=clock,
        ),
        registry=registry,
        data_registry=data,
        validated_inputs=validated,
        snapshots=SimpleNamespace(),
        bridge=ControlledBridge(),
        clock=clock,
    )
    return SimpleNamespace(
        actor=actor,
        bridge=service._bridge,
        engine=engine,
        project_id=project.project_id,
        reference_id=reference.revision_id,
        registry=registry,
        service=service,
        validation=validated._validation_service,
    )


@pytest.fixture
def stack(tmp_path: Path):
    built = _authoring_stack(tmp_path, {})
    yield built
    built.engine.dispose()


def _open(stack, workflow_id: str):
    created = stack.service.create_session(
        stack.actor,
        project_id=stack.project_id,
        workflow_id=workflow_id,
    )
    opened = stack.service.open_connection(stack.actor, created.sid)
    assert opened.connection_state == "ready"
    return opened


def _suggest(stack, sid: str, revision: int, turn_id: str, *, preference: int | None):
    preferences = []
    if preference is not None:
        preferences.append({"field_token": "options.mapq", "value": preference})
    return stack.service.start_turn(
        stack.actor,
        sid,
        client_turn_id=turn_id,
        expected_revision=revision,
        question={
            "intent": "suggest",
            "field_tokens": ["options.mapq"],
            "preferences": preferences,
        },
    )


def test_hitrac_mapq_can_be_adopted_rejected_and_validated(stack) -> None:
    opened = _open(stack, HITRAC_WORKFLOW)
    capability = opened.capability
    assert capability["available"] is True
    assert [field["token"] for field in capability["fields"]] == ["options.mapq"]
    constraints = capability["fields"][0]["constraints"]
    assert constraints == {"type": "integer", "minimum": 0, "maximum": 255}
    assert [issue["code"] for issue in capability["allowed_issue_codes"]] == [
        "HITRAC_OPTION_INVALID",
        "HITRAC_OPTIONS_INVALID",
    ]
    assert "ENCODE_" not in json.dumps(capability["allowed_issue_codes"])

    stack.bridge.value = 20
    suggested = _suggest(stack, opened.sid, 0, "legal", preference=20)
    assert suggested.status == "completed"
    assert suggested.adoptable is True
    assert suggested.proposal["changes"][0]["value"] == 20
    adopted = stack.service.adopt(
        stack.actor,
        opened.sid,
        suggested.proposal_id,
        expected_revision=0,
        proposal_digest_value=suggested.proposal_digest,
    )
    current = stack.service.get_session(stack.actor, opened.sid)
    assert current.draft["options"]["mapq"] == 20
    assert adopted.draft_revision == current.draft_revision

    stack.bridge.value = 256
    rejected = _suggest(
        stack,
        opened.sid,
        current.draft_revision,
        "out-of-bounds",
        preference=None,
    )
    assert rejected.status == "failed"
    assert rejected.failure_reason == MODEL_OUTPUT_REJECTED
    unchanged = stack.service.get_session(stack.actor, opened.sid)
    assert unchanged.draft["options"]["mapq"] == 20
    assert unchanged.draft_revision == current.draft_revision

    with pytest.raises(AuthoringRejected) as preference:
        _suggest(
            stack,
            opened.sid,
            unchanged.draft_revision,
            "preference-bounds",
            preference=256,
        )
    assert preference.value.reason_code == INVALID_REQUEST
    assert (
        stack.service.get_session(stack.actor, opened.sid).draft["options"]["mapq"]
        == 20
    )

    edited = dict(unchanged.draft)
    edited["options"] = {"mapq": 30}
    edited["samples"] = [dict(SAMPLE)]
    edited["sample_revision_ids"] = [SAMPLE_REVISION]
    edited["reference_profile_revision_id"] = stack.reference_id
    saved = stack.service.save_draft(
        stack.actor,
        opened.sid,
        expected_revision=unchanged.draft_revision,
        content=edited,
    )
    assert saved.draft["options"]["mapq"] == 30
    checked = stack.service.validate(
        stack.actor,
        opened.sid,
        expected_revision=saved.draft_revision,
    )
    assert [issue.code for issue in checked.result.issues] == [
        "REFERENCE_PROFILE_BINDING_INVALID"
    ]
    direct = stack.validation.validate(
        HITRAC_WORKFLOW,
        WorkflowInputs(config={}, samples=[dict(SAMPLE)], options={"mapq": 30}),
    )
    assert direct.is_success
    assert direct.value.options == {"threads": 2, "mapq": 30}

    edited["options"] = {"mapq": 256}
    saved_again = stack.service.save_draft(
        stack.actor,
        opened.sid,
        expected_revision=saved.draft_revision,
        content=edited,
    )
    rejected_value = stack.service.validate(
        stack.actor,
        opened.sid,
        expected_revision=saved_again.draft_revision,
    )
    assert rejected_value.result.issues[0].code == "HITRAC_OPTION_INVALID"
    assert rejected_value.result.issues[0].path == "options.mapq"
    direct_rejected = stack.validation.validate(
        HITRAC_WORKFLOW,
        WorkflowInputs(config={}, samples=[dict(SAMPLE)], options={"mapq": 256}),
    )
    assert direct_rejected.issues[0].code == "HITRAC_OPTION_INVALID"


def test_existing_encode_projection_still_opens_and_validates(stack) -> None:
    schema = build_encode_authoring_schema().to_dict()
    fields = project_encode_fields(ENCODE_WORKFLOW_ID, schema)
    assert fields is not None
    historical = projection_identity_digest(
        ENCODE_WORKFLOW_ID,
        fields,
        explainable_issues(),
    )
    assert historical == ENCODE_PROJECTION_DIGEST

    opened = _open(stack, ENCODE_WORKFLOW_ID)
    assert opened.projection_digest == ENCODE_PROJECTION_DIGEST
    assert {field["token"] for field in opened.capability["fields"]} == {
        "config.threads",
        "options.cores",
    }
    empty = stack.service.validate(
        stack.actor,
        opened.sid,
        expected_revision=opened.draft_revision,
    )
    assert empty.result.issues[0].code == SAMPLE_SELECTION_REQUIRED

    edited = dict(opened.draft)
    edited["config"] = {"threads": 4}
    saved = stack.service.save_draft(
        stack.actor,
        opened.sid,
        expected_revision=opened.draft_revision,
        content=edited,
    )
    assert saved.projection_digest == ENCODE_PROJECTION_DIGEST
    assert saved.draft["config"]["threads"] == 4
    assert saved.connection_state == "ready"
    again = stack.service.validate(
        stack.actor,
        opened.sid,
        expected_revision=saved.draft_revision,
    )
    assert again.result.issues[0].code == SAMPLE_SELECTION_REQUIRED


def test_bulk_rnaseq_still_has_no_suggestion_catalog(stack) -> None:
    created = stack.service.create_session(
        stack.actor,
        project_id=stack.project_id,
        workflow_id=BULK_WORKFLOW,
    )
    assert created.capability["available"] is False
    assert created.capability["reason_code"] == CAPABILITY_UNAVAILABLE
    assert created.capability["fields"] == []
    with pytest.raises(AuthoringRejected) as denied:
        stack.service.open_connection(stack.actor, created.sid)
    assert denied.value.reason_code == CAPABILITY_UNAVAILABLE


def test_results_adapter_adopts_the_same_mapq_projection(tmp_path: Path, bound_case):
    base = _authoring_stack(tmp_path / "authoring-only", {})
    ready = _authoring_stack(
        tmp_path / "results-ready",
        {
            "HELIXWEAVE_HITRAC_RUNTIME_BINDING": str(bound_case["runtime_path"]),
            "HELIXWEAVE_HITRAC_RUNTIME_SHA256": bound_case["runtime_sha"],
        },
    )
    try:
        base_adapter = base.registry.get(HITRAC_WORKFLOW)
        ready_adapter = ready.registry.get(HITRAC_WORKFLOW)
        assert type(base_adapter) is HiTracPreprocessAdapter
        assert type(ready_adapter) is HiTracPreprocessResultsAdapter
        assert ready_adapter.execution_availability().execution == "available"

        opened = _open(base, HITRAC_WORKFLOW)
        ready_opened = _open(ready, HITRAC_WORKFLOW)
        assert ready_opened.connection_state == "ready"
        assert ready_opened.projection_digest == opened.projection_digest
        assert ready_opened.capability["fields"] == opened.capability["fields"]
        assert ready_opened.capability["fields"] == [
            {
                "constraints": {"type": "integer", "minimum": 0, "maximum": 255},
                "label": "MAPQ threshold",
                "token": "options.mapq",
                "type": "integer",
            }
        ]

        ready.bridge.value = 20
        suggested = _suggest(ready, ready_opened.sid, 0, "legal", preference=20)
        assert suggested.status == "completed"
        assert suggested.adoptable is True
        adopted = ready.service.adopt(
            ready.actor,
            ready_opened.sid,
            suggested.proposal_id,
            expected_revision=0,
            proposal_digest_value=suggested.proposal_digest,
        )
        current = ready.service.get_session(ready.actor, ready_opened.sid)
        assert current.draft["options"]["mapq"] == 20
        assert current.draft_revision == opened.draft_revision + 1
        assert adopted.draft_revision == current.draft_revision
    finally:
        base.engine.dispose()
        ready.engine.dispose()
