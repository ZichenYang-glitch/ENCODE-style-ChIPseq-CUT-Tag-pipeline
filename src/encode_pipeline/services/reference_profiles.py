"""Reference Profile catalog and exact revision verification."""

from __future__ import annotations

from collections.abc import Callable, Mapping
from datetime import datetime, timezone
from typing import Protocol, cast
from uuid import uuid4

from encode_pipeline.platform.adapters import WorkflowAdapter
from encode_pipeline.platform.reference_profiles import (
    AdapterReferenceBindingIdentity,
    ReferenceProfile,
    ReferenceProfileRevision,
    ReferenceProfileRevisionSummary,
    ReferenceProfileWorkflowBinding,
    build_reference_profile_revision,
)
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.services.private_reference_profiles import (
    PrivateReferenceProfileConfig,
    PrivateReferenceProfileConfigError,
)
from encode_pipeline.services.reference_profile_repositories import (
    InMemoryReferenceProfileRepository,
    ReferenceProfileConflictError,
    ReferenceProfileRepository,
)


REFERENCE_PROFILE_ISSUE_MESSAGES = {
    "REFERENCE_PROFILE_REQUIRED": "A Reference Profile revision is required.",
    "REFERENCE_PROFILE_UNAVAILABLE": "The selected Reference Profile is unavailable.",
    "REFERENCE_PROFILE_DISABLED": "The selected Reference Profile is disabled.",
    "REFERENCE_PROFILE_STALE": "The selected Reference Profile revision is stale.",
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


class ReferenceProfileAdminError(RuntimeError):
    """Stable redacted administrator command error."""

    def __init__(self, reason_code: str) -> None:
        self.reason_code = reason_code
        super().__init__(REFERENCE_PROFILE_ISSUE_MESSAGES[reason_code])


class _ReferenceBindingAdapter(Protocol):
    def verify_reference_profile_binding(
        self,
        payload: Mapping[str, object],
    ) -> Result[AdapterReferenceBindingIdentity]: ...


PrivateConfigProvider = Callable[[], PrivateReferenceProfileConfig]
AdapterProvider = Callable[[str], object]


class ReferenceProfileService:
    """Application service for administrator catalog mutations and public reads."""

    def __init__(
        self,
        *,
        repository: ReferenceProfileRepository | None = None,
        private_config_provider: PrivateConfigProvider,
        adapter_provider: AdapterProvider,
        profile_id_factory: Callable[[], str] | None = None,
        revision_id_factory: Callable[[], str] | None = None,
        now_factory: Callable[[], datetime] | None = None,
    ) -> None:
        self._repository = (
            repository
            if repository is not None
            else InMemoryReferenceProfileRepository()
        )
        self._private_config_provider = private_config_provider
        self._adapter_provider = adapter_provider
        self._profile_id_factory = profile_id_factory or (lambda: f"refp_{uuid4().hex}")
        self._revision_id_factory = revision_id_factory or (
            lambda: f"refpr_{uuid4().hex}"
        )
        self._now_factory = now_factory or (lambda: datetime.now(timezone.utc))

    @property
    def repository(self) -> ReferenceProfileRepository:
        return self._repository

    def register(
        self,
        *,
        safe_key: str,
        display_name: str,
        organism: str,
        assembly: str,
        config_key: str,
    ) -> ReferenceProfileRevisionSummary:
        """Create a stable profile or append a verified immutable revision."""
        config = self._load_config_for_admin()
        identities = self._verify_configured_bindings_for_admin(config, config_key)
        try:
            profile = self._repository.get_profile_by_safe_key(safe_key)
        except KeyError:
            now = self._now_factory()
            profile = ReferenceProfile(
                profile_id=self._profile_id_factory(),
                safe_key=safe_key,
                created_at=now,
            )
            revision = build_reference_profile_revision(
                revision_id=self._revision_id_factory(),
                profile_id=profile.profile_id,
                revision_number=1,
                display_name=display_name,
                organism=organism,
                assembly=assembly,
                config_key=config_key,
                workflow_bindings=tuple(
                    ReferenceProfileWorkflowBinding.from_adapter_identity(identity)
                    for identity in identities
                ),
                created_at=now,
            )
            try:
                self._repository.create_profile(profile, revision)
            except ReferenceProfileConflictError as exc:
                raise ReferenceProfileAdminError(
                    "REFERENCE_PROFILE_UNAVAILABLE"
                ) from exc
        else:
            revisions = self._repository.list_revisions(profile.profile_id)
            if not revisions:
                raise ReferenceProfileAdminError("REFERENCE_PROFILE_UNAVAILABLE")
            current = revisions[-1]
            revision = build_reference_profile_revision(
                revision_id=self._revision_id_factory(),
                profile_id=profile.profile_id,
                revision_number=current.revision_number + 1,
                display_name=display_name,
                organism=organism,
                assembly=assembly,
                config_key=config_key,
                workflow_bindings=tuple(
                    ReferenceProfileWorkflowBinding.from_adapter_identity(identity)
                    for identity in identities
                ),
                created_at=self._now_factory(),
            )
            try:
                self._repository.append_revision(
                    revision,
                    expected_previous_revision_number=current.revision_number,
                )
            except ReferenceProfileConflictError as exc:
                raise ReferenceProfileAdminError(
                    "REFERENCE_PROFILE_UNAVAILABLE"
                ) from exc
        return revision.public_summary(
            safe_key=profile.safe_key,
            enabled=profile.enabled_revision_id == revision.revision_id,
        )

    def verify(self, revision_id: str) -> ReferenceProfileRevisionSummary:
        """Verify one revision without mutating catalog state."""
        try:
            revision = self._repository.get_revision(revision_id)
            profile = self._repository.get_profile(revision.profile_id)
        except KeyError as exc:
            raise ReferenceProfileAdminError("REFERENCE_PROFILE_UNAVAILABLE") from exc
        result = self._verify_revision(revision)
        if result.is_failure:
            raise ReferenceProfileAdminError(result.issues[0].code)
        return revision.public_summary(
            safe_key=profile.safe_key,
            enabled=profile.enabled_revision_id == revision.revision_id,
        )

    def list(self) -> tuple[ReferenceProfileRevisionSummary, ...]:
        """List every append-only revision without private configuration fields."""
        summaries: list[ReferenceProfileRevisionSummary] = []
        for profile in self._repository.list_profiles():
            summaries.extend(
                revision.public_summary(
                    safe_key=profile.safe_key,
                    enabled=(profile.enabled_revision_id == revision.revision_id),
                )
                for revision in self._repository.list_revisions(profile.profile_id)
            )
        return tuple(summaries)

    def get_revision_summary(
        self,
        revision_id: str,
    ) -> ReferenceProfileRevisionSummary:
        """Return path-free historical metadata without private re-verification."""
        revision = self._repository.get_revision(revision_id)
        profile = self._repository.get_profile(revision.profile_id)
        return revision.public_summary(
            safe_key=profile.safe_key,
            enabled=profile.enabled_revision_id == revision.revision_id,
        )

    def enable(
        self,
        profile_id: str,
        *,
        revision_id: str | None = None,
    ) -> ReferenceProfileRevisionSummary:
        """Verify and enable an exact revision; never mutate verification state."""
        try:
            profile = self._repository.get_profile(profile_id)
            revisions = self._repository.list_revisions(profile_id)
            target = (
                revisions[-1]
                if revision_id is None
                else self._repository.get_revision(revision_id)
            )
        except (IndexError, KeyError) as exc:
            raise ReferenceProfileAdminError("REFERENCE_PROFILE_UNAVAILABLE") from exc
        if target.profile_id != profile.profile_id:
            raise ReferenceProfileAdminError("REFERENCE_PROFILE_UNAVAILABLE")
        verified = self._verify_revision(target)
        if verified.is_failure:
            raise ReferenceProfileAdminError(verified.issues[0].code)
        profile = self._repository.set_enabled_revision(profile_id, target.revision_id)
        return target.public_summary(safe_key=profile.safe_key, enabled=True)

    def disable(self, profile_id: str) -> ReferenceProfileRevisionSummary:
        """Disable new validation/create/start while preserving revision history."""
        try:
            profile = self._repository.set_enabled_revision(profile_id, None)
            revisions = self._repository.list_revisions(profile_id)
        except (IndexError, KeyError) as exc:
            raise ReferenceProfileAdminError("REFERENCE_PROFILE_UNAVAILABLE") from exc
        return revisions[-1].public_summary(safe_key=profile.safe_key, enabled=False)

    def list_enabled(
        self,
        workflow_id: str,
    ) -> tuple[ReferenceProfileRevisionSummary, ...]:
        """Return enabled exact revisions compatible with one workflow."""
        return tuple(
            revision.public_summary(safe_key=profile.safe_key, enabled=True)
            for profile, revision in self._repository.list_enabled_for_workflow(
                workflow_id
            )
        )

    def _verify_revision(
        self,
        revision: ReferenceProfileRevision,
    ) -> Result[tuple[AdapterReferenceBindingIdentity, ...]]:
        try:
            config = self._private_config_provider()
        except (PrivateReferenceProfileConfigError, RuntimeError, ValueError):
            return _failure("REFERENCE_PROFILE_CONFIG_INVALID")
        identities = []
        for stored in revision.workflow_bindings:
            try:
                payload = config.binding_for(revision.config_key, stored.workflow_id)
                adapter = self._adapter_for(stored.workflow_id)
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
            if verified.value != stored.to_adapter_identity():
                return _failure("REFERENCE_PROFILE_IDENTITY_MISMATCH")
            identities.append(verified.value)
        return Result.success(tuple(identities))

    def _load_config_for_admin(self) -> PrivateReferenceProfileConfig:
        try:
            return self._private_config_provider()
        except (PrivateReferenceProfileConfigError, RuntimeError, ValueError) as exc:
            raise ReferenceProfileAdminError(
                "REFERENCE_PROFILE_CONFIG_INVALID"
            ) from exc

    def _verify_configured_bindings_for_admin(
        self,
        config: PrivateReferenceProfileConfig,
        config_key: str,
    ) -> tuple[AdapterReferenceBindingIdentity, ...]:
        try:
            workflow_ids = config.workflow_ids_for(config_key)
        except PrivateReferenceProfileConfigError as exc:
            raise ReferenceProfileAdminError(
                "REFERENCE_PROFILE_CONFIG_INVALID"
            ) from exc
        identities = []
        for workflow_id in workflow_ids:
            try:
                adapter = self._adapter_for(workflow_id)
                payload = config.binding_for(config_key, workflow_id)
            except (
                LookupError,
                PrivateReferenceProfileConfigError,
                RuntimeError,
                ValueError,
            ) as exc:
                raise ReferenceProfileAdminError(
                    "REFERENCE_PROFILE_CONFIG_INVALID"
                ) from exc
            result = _call_verify(adapter, payload)
            if result.is_failure or result.value is None:
                raise ReferenceProfileAdminError("REFERENCE_PROFILE_BINDING_INVALID")
            if result.value.workflow_id != workflow_id:
                raise ReferenceProfileAdminError("REFERENCE_PROFILE_IDENTITY_MISMATCH")
            identities.append(result.value)
        return tuple(identities)

    def _adapter_for(self, workflow_id: str) -> _ReferenceBindingAdapter:
        adapter = self._adapter_provider(workflow_id)
        if not isinstance(adapter, WorkflowAdapter) or not callable(
            getattr(adapter, "verify_reference_profile_binding", None)
        ):
            raise LookupError("adapter is not reference-capable")
        return cast(_ReferenceBindingAdapter, adapter)


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
    return result


def _failure(code: str) -> Result:
    return Result.failure(
        (Issue(code=code, message=REFERENCE_PROFILE_ISSUE_MESSAGES[code]),)
    )
