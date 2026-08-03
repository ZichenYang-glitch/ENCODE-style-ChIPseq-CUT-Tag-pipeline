"""Reference Profile service and in-memory persistence contract tests."""

from __future__ import annotations

from datetime import datetime, timezone
import json

from encode_pipeline.platform.adapters import (
    VALIDATION_CAPABILITY,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
)
from encode_pipeline.platform.reference_profiles import (
    AdapterReferenceBindingIdentity,
    BoundWorkflowReference,
)
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.private_reference_profiles import (
    PrivateReferenceProfileConfig,
    load_private_reference_profile_config,
)
from encode_pipeline.services.reference_profile_repositories import (
    InMemoryReferenceProfileRepository,
)
from encode_pipeline.services.reference_profile_runtime import (
    ReferenceProfileBindingService,
    ReferenceProfileRuntimeResolver,
)
from encode_pipeline.services.reference_profiles import ReferenceProfileService


NOW = datetime(2026, 8, 3, 9, 0, tzinfo=timezone.utc)


class _Adapter:
    def __init__(self, workflow_id: str, identity: str) -> None:
        self.workflow_id = workflow_id
        self.identity = identity
        self.metadata = WorkflowMetadata(
            workflow_id=workflow_id,
            name="Test workflow",
            version="test-adapter-v1",
        )
        self.capabilities = WorkflowCapabilities(supports=(VALIDATION_CAPABILITY,))

    def schema(self):
        return object()

    def validate(self, inputs):
        return Result.success(inputs)

    def preview_dag(self, inputs):
        return Result.success(object())

    def plan_workspace(self, inputs, workspace):
        return Result.success(object())

    def build_command(self, plan, workspace):
        return Result.success(object())

    def extract_artifacts(self, inputs, workspace):
        return Result.success(())

    def verify_reference_profile_binding(self, payload):
        return Result.success(
            AdapterReferenceBindingIdentity(
                workflow_id=self.workflow_id,
                contract_version=f"{self.workflow_id}-reference-v1",
                identity_sha256=self.identity,
            )
        )

    def bind_reference_profile(self, inputs, payload):
        identity = self.verify_reference_profile_binding(payload).value
        assert identity is not None
        return Result.success(
            BoundWorkflowReference(
                inputs=inputs,
                adapter=type(self)(self.workflow_id, self.identity),
                identity=identity,
            )
        )


def _service(
    tmp_path,
) -> tuple[
    ReferenceProfileService,
    ReferenceProfileBindingService,
    _Adapter,
    object,
]:
    config_path = tmp_path / "references.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": "helixweave-reference-profiles-v1",
                "profiles": {
                    "grch38-private": {
                        "bindings": {"bulk-rnaseq": {"fasta": "/operator/GRCh38.fa"}}
                    }
                },
            }
        ),
        encoding="utf-8",
    )
    repository = InMemoryReferenceProfileRepository()
    adapter = _Adapter("bulk-rnaseq", "a" * 64)
    service = ReferenceProfileService(
        repository=repository,
        private_config_provider=lambda: load_private_reference_profile_config(
            config_path
        ),
        adapter_provider=lambda workflow_id: (
            adapter if workflow_id == "bulk-rnaseq" else None
        ),
        profile_id_factory=lambda: "refp_" + "1" * 32,
        revision_id_factory=lambda: "refpr_" + "2" * 32,
        now_factory=lambda: NOW,
    )
    binding_service = ReferenceProfileBindingService(
        repository=repository,
        private_config_provider=lambda: load_private_reference_profile_config(
            config_path
        ),
        adapter_provider=lambda workflow_id: (
            adapter if workflow_id == "bulk-rnaseq" else None
        ),
        now_factory=lambda: NOW,
    )
    return service, binding_service, adapter, repository


def test_register_is_stable_profile_append_only_revision_and_enable_is_verified(
    tmp_path,
) -> None:
    service, _binding_service, _adapter, repository = _service(tmp_path)

    first = service.register(
        safe_key="grch38",
        display_name="GRCh38",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )
    service._revision_id_factory = lambda: "refpr_" + "3" * 32
    second = service.register(
        safe_key="grch38",
        display_name="GRCh38 patched metadata",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )

    assert first.profile_id == second.profile_id
    assert (first.revision_number, second.revision_number) == (1, 2)
    assert repository.get_profile(first.profile_id).enabled_revision_id is None
    enabled = service.enable(first.profile_id, revision_id=second.revision_id)
    assert enabled.enabled is True
    listed = service.list()
    assert [item.revision_id for item in listed] == [
        first.revision_id,
        second.revision_id,
    ]
    assert [item.enabled for item in listed] == [False, True]
    assert service.list_enabled("bulk-rnaseq") == (enabled,)
    service.disable(first.profile_id)
    assert service.list_enabled("bulk-rnaseq") == ()


def test_catalog_service_does_not_own_execution_binding(tmp_path) -> None:
    catalog, binding_service, _adapter, _repository = _service(tmp_path)

    assert not hasattr(catalog, "resolve_selection")
    assert not hasattr(catalog, "resolve_evidence")
    assert callable(binding_service.resolve_selection)
    assert callable(binding_service.resolve_evidence)


def test_multiple_enabled_profiles_coexist_for_one_workflow() -> None:
    repository = InMemoryReferenceProfileRepository()
    adapter = _Adapter("bulk-rnaseq", "a" * 64)
    profile_ids = iter(("refp_" + "1" * 32, "refp_" + "2" * 32))
    revision_ids = iter(("refpr_" + "3" * 32, "refpr_" + "4" * 32))
    config = PrivateReferenceProfileConfig(
        {
            "grch38-private": {"bulk-rnaseq": {}},
            "mm10-private": {"bulk-rnaseq": {}},
        }
    )
    service = ReferenceProfileService(
        repository=repository,
        private_config_provider=lambda: config,
        adapter_provider=lambda _workflow_id: adapter,
        profile_id_factory=lambda: next(profile_ids),
        revision_id_factory=lambda: next(revision_ids),
        now_factory=lambda: NOW,
    )

    grch38 = service.register(
        safe_key="grch38",
        display_name="GRCh38 primary",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )
    mm10 = service.register(
        safe_key="mm10",
        display_name="mm10 primary",
        organism="Mus musculus",
        assembly="mm10",
        config_key="mm10-private",
    )
    service.enable(grch38.profile_id, revision_id=grch38.revision_id)
    service.enable(mm10.profile_id, revision_id=mm10.revision_id)

    enabled = service.list_enabled("bulk-rnaseq")
    assert [item.assembly for item in enabled] == ["GRCh38", "mm10"]
    assert all(item.enabled for item in enabled)


def test_one_profile_can_bind_two_workflows_and_filters_compatibility() -> None:
    repository = InMemoryReferenceProfileRepository()
    adapters = {
        "encode-workflow": _Adapter("encode-workflow", "a" * 64),
        "bulk-rnaseq": _Adapter("bulk-rnaseq", "b" * 64),
    }
    config = PrivateReferenceProfileConfig(
        {"shared-grch38": {workflow_id: {} for workflow_id in adapters}}
    )
    service = ReferenceProfileService(
        repository=repository,
        private_config_provider=lambda: config,
        adapter_provider=adapters.__getitem__,
        profile_id_factory=lambda: "refp_" + "5" * 32,
        revision_id_factory=lambda: "refpr_" + "6" * 32,
        now_factory=lambda: NOW,
    )

    revision = service.register(
        safe_key="shared-grch38",
        display_name="GRCh38 shared",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="shared-grch38",
    )
    service.enable(revision.profile_id, revision_id=revision.revision_id)

    assert revision.compatible_workflow_ids == ("bulk-rnaseq", "encode-workflow")
    assert service.list_enabled("bulk-rnaseq")[0].revision_id == revision.revision_id
    assert (
        service.list_enabled("encode-workflow")[0].revision_id == revision.revision_id
    )
    assert service.list_enabled("unrelated-workflow") == ()


def test_resolve_selection_reloads_private_config_and_detects_identity_drift(
    tmp_path,
) -> None:
    service, binding_service, adapter, _repository = _service(tmp_path)
    summary = service.register(
        safe_key="grch38",
        display_name="GRCh38",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )
    service.enable(summary.profile_id, revision_id=summary.revision_id)
    inputs = WorkflowInputs(config={}, samples=None, options={})

    resolved = binding_service.resolve_selection(
        "bulk-rnaseq", summary.revision_id, inputs, require_enabled=True
    )
    assert resolved.is_success
    assert resolved.value is not None
    assert resolved.value.evidence.revision_id == summary.revision_id
    assert resolved.value.bound_reference.adapter is not adapter

    adapter.identity = "b" * 64
    stale = binding_service.resolve_selection(
        "bulk-rnaseq", summary.revision_id, inputs, require_enabled=True
    )
    assert stale.is_failure
    assert [issue.code for issue in stale.issues] == [
        "REFERENCE_PROFILE_IDENTITY_MISMATCH"
    ]


def test_binding_rejects_repository_revision_misdirection(tmp_path) -> None:
    catalog, binding_service, _adapter, repository = _service(tmp_path)
    first = catalog.register(
        safe_key="grch38",
        display_name="GRCh38",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )
    catalog._revision_id_factory = lambda: "refpr_" + "3" * 32
    second = catalog.register(
        safe_key="grch38",
        display_name="GRCh38 v2",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )
    second_revision = repository.get_revision(second.revision_id)
    repository.get_revision = lambda _revision_id: second_revision  # type: ignore[method-assign]

    result = binding_service.resolve_selection(
        "bulk-rnaseq",
        first.revision_id,
        WorkflowInputs(config={}, samples=None, options={}),
        require_enabled=False,
    )

    assert [issue.code for issue in result.issues] == [
        "REFERENCE_PROFILE_IDENTITY_MISMATCH"
    ]


def test_disabled_or_stale_revision_fails_closed_with_stable_issue(tmp_path) -> None:
    service, binding_service, _adapter, _repository = _service(tmp_path)
    summary = service.register(
        safe_key="grch38",
        display_name="GRCh38",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )
    inputs = WorkflowInputs(config={}, samples=None, options={})

    disabled = binding_service.resolve_selection(
        "bulk-rnaseq", summary.revision_id, inputs, require_enabled=True
    )
    assert [issue.code for issue in disabled.issues] == ["REFERENCE_PROFILE_DISABLED"]

    service.enable(summary.profile_id, revision_id=summary.revision_id)
    service._revision_id_factory = lambda: "refpr_" + "3" * 32
    second = service.register(
        safe_key="grch38",
        display_name="GRCh38 v2",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )
    service.enable(summary.profile_id, revision_id=second.revision_id)
    stale = binding_service.resolve_selection(
        "bulk-rnaseq", summary.revision_id, inputs, require_enabled=True
    )
    assert [issue.code for issue in stale.issues] == ["REFERENCE_PROFILE_STALE"]


def test_historical_summary_is_path_free_and_does_not_require_private_verify(
    tmp_path,
) -> None:
    service, _binding_service, adapter, _repository = _service(tmp_path)
    summary = service.register(
        safe_key="grch38",
        display_name="GRCh38",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )
    adapter.identity = "b" * 64

    observed = service.get_revision_summary(summary.revision_id)

    assert observed == summary
    assert "config_key" not in observed.to_public_dict()


def test_verify_is_read_only_for_catalog_state(tmp_path) -> None:
    service, _binding_service, _adapter, repository = _service(tmp_path)
    summary = service.register(
        safe_key="grch38",
        display_name="GRCh38",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )
    profiles_before = repository.list_profiles()
    revisions_before = repository.list_revisions(summary.profile_id)

    verified = service.verify(summary.revision_id)

    assert verified == summary
    assert repository.list_profiles() == profiles_before
    assert repository.list_revisions(summary.profile_id) == revisions_before


def test_resolve_selection_preserves_safe_user_binding_conflict(tmp_path) -> None:
    service, binding_service, adapter, _repository = _service(tmp_path)
    summary = service.register(
        safe_key="grch38",
        display_name="GRCh38",
        organism="Homo sapiens",
        assembly="GRCh38",
        config_key="grch38-private",
    )
    service.enable(summary.profile_id, revision_id=summary.revision_id)
    adapter.bind_reference_profile = lambda _inputs, _payload: Result.failure(  # type: ignore[method-assign]
        [
            Issue(
                code="TEST_REFERENCE_BINDING_CONFLICT",
                message="Submitted inputs conflict with the selected reference.",
                source="adapter",
                path="samples",
            )
        ]
    )

    result = binding_service.resolve_selection(
        "bulk-rnaseq",
        summary.revision_id,
        WorkflowInputs(config={}, samples=None, options={}),
    )

    assert [issue.code for issue in result.issues] == [
        "TEST_REFERENCE_BINDING_CONFLICT"
    ]


def test_runtime_resolver_uses_explicit_workflow_and_requires_capable_binding(
    tmp_path,
) -> None:
    _service_catalog, binding_service, adapter, _repository = _service(tmp_path)
    inputs = WorkflowInputs(config={}, samples=None, options={})

    class Runs:
        evidence = None

        def get_run_reference_binding(self, run_id):
            return self.evidence

    class Registry:
        def __init__(self, selected):
            self.selected = selected

        def get(self, workflow_id):
            return self.selected

    runs = Runs()
    capable = ReferenceProfileRuntimeResolver(runs, Registry(adapter), binding_service)
    missing = capable.resolve_run("run-1", "bulk-rnaseq", inputs, require_enabled=True)
    assert [issue.code for issue in missing.issues] == ["REFERENCE_PROFILE_REQUIRED"]

    legacy = ReferenceProfileRuntimeResolver(runs, Registry(object()), binding_service)
    assert (
        legacy.resolve_run(
            "run-2", "legacy-workflow", inputs, require_enabled=True
        ).value
        is None
    )
