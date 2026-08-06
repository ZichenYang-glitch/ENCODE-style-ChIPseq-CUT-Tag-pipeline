"""Execution-owned binding and resolution of immutable reference evidence."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from datetime import datetime, timezone
from typing import Protocol, cast

from encode_pipeline.platform.adapters import (
    ReferenceProfileBindingAdapter,
    WorkflowAdapter,
    WorkflowCapabilities,
    WorkflowInputs,
    WorkflowMetadata,
)
from encode_pipeline.platform.reference_profiles import (
    AdapterReferenceBindingIdentity,
    BoundWorkflowReference,
    ReferenceProfileRevisionBinding,
    ResolvedReferenceProfile,
    build_reference_profile_revision_binding,
)
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.private_reference_profiles import (
    PrivateReferenceProfileConfig,
    PrivateReferenceProfileConfigError,
)
from encode_pipeline.services.reference_profile_repositories import (
    ReferenceProfileRepository,
)


class _ReferenceBindingAdapter(Protocol):
    def verify_reference_profile_binding(
        self,
        payload: Mapping[str, object],
    ) -> Result[AdapterReferenceBindingIdentity]: ...

    def bind_reference_profile(
        self,
        inputs: WorkflowInputs,
        payload: Mapping[str, object],
    ) -> Result[BoundWorkflowReference]: ...


class _RunReferenceBindingRepository(Protocol):
    def get_run_reference_binding(
        self,
        run_id: str,
    ) -> ReferenceProfileRevisionBinding | None: ...


class _Registry(Protocol):
    def get(self, workflow_id: str) -> object: ...


PrivateConfigProvider = Callable[[], PrivateReferenceProfileConfig]
AdapterProvider = Callable[[str], object]


class ReferenceProfileBindingService:
    """Bind exact catalog revisions to adapter-owned execution inputs."""

    def __init__(
        self,
        *,
        repository: ReferenceProfileRepository,
        private_config_provider: PrivateConfigProvider,
        adapter_provider: AdapterProvider,
        now_factory: Callable[[], datetime] | None = None,
    ) -> None:
        self._repository = repository
        self._private_config_provider = private_config_provider
        self._adapter_provider = adapter_provider
        self._now_factory = now_factory or (lambda: datetime.now(timezone.utc))

    def resolve_selection(
        self,
        workflow_id: str,
        revision_id: str,
        inputs: WorkflowInputs,
        *,
        require_enabled: bool = True,
    ) -> Result[ResolvedReferenceProfile]:
        """Re-read private config and bind one selected exact revision."""
        try:
            revision = self._repository.get_revision(revision_id)
            profile = self._repository.get_profile(revision.profile_id)
        except (KeyError, ValueError):
            return _failure("REFERENCE_PROFILE_UNAVAILABLE")
        if (
            revision.revision_id != revision_id
            or profile.profile_id != revision.profile_id
        ):
            return _failure("REFERENCE_PROFILE_IDENTITY_MISMATCH")
        try:
            stored_identity = revision.binding_for(workflow_id).to_adapter_identity()
        except (KeyError, ValueError):
            return _failure("REFERENCE_PROFILE_INCOMPATIBLE")
        if require_enabled:
            if profile.enabled_revision_id is None:
                return _failure("REFERENCE_PROFILE_DISABLED")
            if profile.enabled_revision_id != revision.revision_id:
                return _failure("REFERENCE_PROFILE_STALE")
        try:
            config = self._private_config_provider()
            payload = config.binding_for(revision.config_key, workflow_id)
            adapter = self._adapter_for(workflow_id)
        except (
            LookupError,
            PrivateReferenceProfileConfigError,
            RuntimeError,
            ValueError,
        ):
            return _failure("REFERENCE_PROFILE_CONFIG_INVALID")
        verified = _call_verify(adapter, payload)
        if verified.is_failure or verified.value is None:
            return _failure("REFERENCE_PROFILE_BINDING_INVALID")
        if verified.value != stored_identity:
            return _failure("REFERENCE_PROFILE_IDENTITY_MISMATCH")
        bound = _call_bind(adapter, inputs, payload)
        if bound.is_failure or bound.value is None:
            operator_failure_suffixes = (
                "_REFERENCE_BINDING_INVALID",
                "_REFERENCE_BINDING_UNAVAILABLE",
            )
            has_private_operator_failure = any(
                issue.code.endswith(operator_failure_suffixes) for issue in bound.issues
            )
            if bound.issues and not has_private_operator_failure:
                return Result.failure(bound.issues)
            return _failure("REFERENCE_PROFILE_BINDING_INVALID")
        if (
            not _is_compatible_bound_adapter(
                source=adapter,
                bound=bound.value.adapter,
                workflow_id=workflow_id,
            )
            or bound.value.identity != stored_identity
            or bound.value.identity != verified.value
        ):
            return _failure("REFERENCE_PROFILE_IDENTITY_MISMATCH")
        evidence = build_reference_profile_revision_binding(
            profile_id=profile.profile_id,
            revision_id=revision.revision_id,
            workflow_id=workflow_id,
            revision_public_identity_sha256=revision.public_identity_sha256,
            adapter_identity=stored_identity,
            bound_at=self._now_factory(),
        )
        summary = revision.public_summary(
            safe_key=profile.safe_key,
            enabled=profile.enabled_revision_id == revision.revision_id,
        )
        return Result.success(
            ResolvedReferenceProfile(
                summary=summary,
                evidence=evidence,
                bound_reference=bound.value,
            )
        )

    def resolve_evidence(
        self,
        evidence: ReferenceProfileRevisionBinding,
        inputs: WorkflowInputs,
        *,
        require_enabled: bool,
    ) -> Result[ResolvedReferenceProfile]:
        """Resolve stored snapshot/run evidence and reject any identity drift."""
        if not isinstance(evidence, ReferenceProfileRevisionBinding):
            return _failure("REFERENCE_PROFILE_UNAVAILABLE")
        resolved = self.resolve_selection(
            evidence.workflow_id,
            evidence.revision_id,
            inputs,
            require_enabled=require_enabled,
        )
        if resolved.is_failure or resolved.value is None:
            return resolved
        current = resolved.value.evidence
        comparable = (
            "profile_id",
            "revision_id",
            "workflow_id",
            "revision_public_identity_scheme",
            "revision_public_identity_sha256",
            "adapter_contract_version",
            "adapter_identity_scheme",
            "adapter_identity_sha256",
            "binding_digest_scheme",
            "binding_digest",
        )
        if any(
            getattr(current, name) != getattr(evidence, name) for name in comparable
        ):
            return _failure("REFERENCE_PROFILE_IDENTITY_MISMATCH")
        return resolved

    def _adapter_for(self, workflow_id: str) -> _ReferenceBindingAdapter:
        adapter = self._adapter_provider(workflow_id)
        if (
            adapter is None
            or not callable(getattr(adapter, "verify_reference_profile_binding", None))
            or not callable(getattr(adapter, "bind_reference_profile", None))
        ):
            raise LookupError("adapter is not reference-capable")
        return cast(_ReferenceBindingAdapter, adapter)


class ReferenceProfileRuntimeResolver:
    """Resolve immutable run evidence to private adapter inputs at use time."""

    def __init__(
        self,
        run_repository: _RunReferenceBindingRepository,
        registry: _Registry,
        binding_service: ReferenceProfileBindingService,
    ) -> None:
        self._run_repository = run_repository
        self._registry = registry
        self._binding_service = binding_service

    def resolve_run(
        self,
        run_id: str,
        workflow_id: str,
        inputs: WorkflowInputs,
        *,
        require_enabled: bool,
    ) -> Result[BoundWorkflowReference | None]:
        try:
            evidence = self._run_repository.get_run_reference_binding(run_id)
        except (LookupError, RuntimeError, ValueError):
            return _failure("REFERENCE_PROFILE_UNAVAILABLE")
        try:
            adapter = self._registry.get(workflow_id)
        except (KeyError, ValueError):
            return _failure("REFERENCE_PROFILE_UNAVAILABLE")
        reference_capable = isinstance(adapter, ReferenceProfileBindingAdapter)
        if evidence is None:
            if reference_capable:
                return _failure("REFERENCE_PROFILE_REQUIRED")
            return Result.success(None)
        if evidence.workflow_id != workflow_id or not reference_capable:
            return _failure("REFERENCE_PROFILE_IDENTITY_MISMATCH")
        resolved = self._binding_service.resolve_evidence(
            evidence,
            inputs,
            require_enabled=require_enabled,
        )
        if resolved.is_failure or resolved.value is None:
            return Result.failure(resolved.issues)
        return Result.success(resolved.value.bound_reference)


def _call_verify(
    adapter: _ReferenceBindingAdapter,
    payload: Mapping[str, object],
) -> Result[AdapterReferenceBindingIdentity]:
    try:
        result = adapter.verify_reference_profile_binding(payload)
    except Exception:
        return _failure("REFERENCE_PROFILE_BINDING_INVALID")
    if not isinstance(result, Result):
        return _failure("REFERENCE_PROFILE_BINDING_INVALID")
    if result.is_success and not isinstance(
        result.value,
        AdapterReferenceBindingIdentity,
    ):
        return _failure("REFERENCE_PROFILE_BINDING_INVALID")
    return result


def _call_bind(
    adapter: _ReferenceBindingAdapter,
    inputs: WorkflowInputs,
    payload: Mapping[str, object],
) -> Result[BoundWorkflowReference]:
    try:
        result = adapter.bind_reference_profile(inputs, payload)
    except Exception:
        return _failure("REFERENCE_PROFILE_BINDING_INVALID")
    if not isinstance(result, Result):
        return _failure("REFERENCE_PROFILE_BINDING_INVALID")
    if result.is_success and not isinstance(result.value, BoundWorkflowReference):
        return _failure("REFERENCE_PROFILE_BINDING_INVALID")
    return result


def _is_compatible_bound_adapter(
    *,
    source: object,
    bound: object,
    workflow_id: str,
) -> bool:
    if not isinstance(source, WorkflowAdapter) or not isinstance(
        bound, WorkflowAdapter
    ):
        return False
    source_metadata = source.metadata
    bound_metadata = bound.metadata
    source_capabilities = source.capabilities
    bound_capabilities = bound.capabilities
    if (
        not isinstance(source_metadata, WorkflowMetadata)
        or not isinstance(bound_metadata, WorkflowMetadata)
        or not isinstance(source_capabilities, WorkflowCapabilities)
        or not isinstance(bound_capabilities, WorkflowCapabilities)
    ):
        return False
    if (
        source_metadata.workflow_id != workflow_id
        or bound_metadata.workflow_id != workflow_id
        or bound_metadata.version != source_metadata.version
    ):
        return False
    return set(source_capabilities.supports).issubset(bound_capabilities.supports)


def _failure(code: str) -> Result:
    messages = {
        "REFERENCE_PROFILE_REQUIRED": "A Reference Profile revision is required.",
        "REFERENCE_PROFILE_UNAVAILABLE": (
            "The selected Reference Profile is unavailable."
        ),
        "REFERENCE_PROFILE_DISABLED": "The selected Reference Profile is disabled.",
        "REFERENCE_PROFILE_STALE": (
            "The selected Reference Profile revision is stale."
        ),
        "REFERENCE_PROFILE_INCOMPATIBLE": (
            "The selected Reference Profile is incompatible with this workflow."
        ),
        "REFERENCE_PROFILE_CONFIG_INVALID": (
            "The private Reference Profile configuration is unavailable."
        ),
        "REFERENCE_PROFILE_BINDING_INVALID": (
            "The selected Reference Profile binding could not be verified."
        ),
        "REFERENCE_PROFILE_IDENTITY_MISMATCH": (
            "The selected Reference Profile identity has changed."
        ),
    }
    return Result.failure((Issue(code=code, message=messages[code]),))
