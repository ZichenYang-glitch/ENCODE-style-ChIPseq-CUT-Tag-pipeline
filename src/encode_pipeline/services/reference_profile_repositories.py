"""Persistence boundary for stable profiles and append-only reference revisions."""

from __future__ import annotations

from dataclasses import replace
from threading import RLock
from typing import Protocol

from encode_pipeline.platform.security_audit import SecurityAuditEvent
from encode_pipeline.platform.reference_profiles import (
    ReferenceProfile,
    ReferenceProfileRevision,
)


class ReferenceProfileConflictError(RuntimeError):
    """A stable identity or append-only Reference Profile invariant conflicted."""


class ReferenceProfileRepository(Protocol):
    """Atomic storage contract consumed by ReferenceProfileService."""

    def create_profile(
        self,
        profile: ReferenceProfile,
        initial_revision: ReferenceProfileRevision,
        *,
        security_audit: SecurityAuditEvent | None = None,
    ) -> tuple[ReferenceProfile, ReferenceProfileRevision]: ...

    def get_profile(self, profile_id: str) -> ReferenceProfile: ...

    def get_profile_by_safe_key(self, safe_key: str) -> ReferenceProfile: ...

    def get_profile_for_revision(self, revision_id: str) -> ReferenceProfile: ...

    def list_profiles(self) -> tuple[ReferenceProfile, ...]: ...

    def get_revision(self, revision_id: str) -> ReferenceProfileRevision: ...

    def list_revisions(
        self,
        profile_id: str,
    ) -> tuple[ReferenceProfileRevision, ...]: ...

    def append_revision(
        self,
        revision: ReferenceProfileRevision,
        *,
        expected_previous_revision_number: int,
        security_audit: SecurityAuditEvent | None = None,
    ) -> ReferenceProfileRevision: ...

    def set_enabled_revision(
        self,
        profile_id: str,
        revision_id: str | None,
        *,
        security_audit: SecurityAuditEvent | None = None,
    ) -> ReferenceProfile: ...

    def list_enabled_for_workflow(
        self,
        workflow_id: str,
    ) -> tuple[tuple[ReferenceProfile, ReferenceProfileRevision], ...]: ...


class InMemoryReferenceProfileRepository:
    """Lock-protected implementation with SQL-equivalent invariants."""

    def __init__(self) -> None:
        self._profiles: dict[str, ReferenceProfile] = {}
        self._profile_id_by_safe_key: dict[str, str] = {}
        self._revisions: dict[str, ReferenceProfileRevision] = {}
        self._revision_ids_by_profile: dict[str, list[str]] = {}
        self._security_audits: list[SecurityAuditEvent] = []
        self._lock = RLock()

    def create_profile(
        self,
        profile: ReferenceProfile,
        initial_revision: ReferenceProfileRevision,
        *,
        security_audit: SecurityAuditEvent | None = None,
    ) -> tuple[ReferenceProfile, ReferenceProfileRevision]:
        if not isinstance(profile, ReferenceProfile):
            raise ValueError("profile must be a ReferenceProfile")
        if not isinstance(initial_revision, ReferenceProfileRevision):
            raise ValueError("initial_revision must be a ReferenceProfileRevision")
        if profile.enabled_revision_id is not None:
            raise ValueError("a new ReferenceProfile must begin disabled")
        if (
            initial_revision.profile_id != profile.profile_id
            or initial_revision.revision_number != 1
            or initial_revision.created_at < profile.created_at
        ):
            raise ValueError("initial revision does not match ReferenceProfile")
        with self._lock:
            if profile.profile_id in self._profiles:
                raise ReferenceProfileConflictError("duplicate profile identity")
            if profile.safe_key in self._profile_id_by_safe_key:
                raise ReferenceProfileConflictError("profile safe_key already exists")
            if initial_revision.revision_id in self._revisions:
                raise ReferenceProfileConflictError("duplicate revision identity")
            self._profiles[profile.profile_id] = profile
            self._profile_id_by_safe_key[profile.safe_key] = profile.profile_id
            self._revisions[initial_revision.revision_id] = initial_revision
            self._revision_ids_by_profile[profile.profile_id] = [
                initial_revision.revision_id
            ]
            if security_audit is not None:
                self._security_audits.append(security_audit)
            return profile, initial_revision

    def get_profile(self, profile_id: str) -> ReferenceProfile:
        with self._lock:
            try:
                return self._profiles[profile_id]
            except KeyError:
                raise KeyError("unknown reference profile") from None

    def get_profile_by_safe_key(self, safe_key: str) -> ReferenceProfile:
        with self._lock:
            try:
                return self.get_profile(self._profile_id_by_safe_key[safe_key])
            except KeyError:
                raise KeyError("unknown reference profile") from None

    def get_profile_for_revision(self, revision_id: str) -> ReferenceProfile:
        revision = self.get_revision(revision_id)
        return self.get_profile(revision.profile_id)

    def list_profiles(self) -> tuple[ReferenceProfile, ...]:
        with self._lock:
            return tuple(
                sorted(
                    self._profiles.values(),
                    key=lambda profile: (profile.created_at, profile.profile_id),
                )
            )

    def get_revision(self, revision_id: str) -> ReferenceProfileRevision:
        with self._lock:
            try:
                return self._revisions[revision_id]
            except KeyError:
                raise KeyError("unknown reference profile revision") from None

    def list_revisions(
        self,
        profile_id: str,
    ) -> tuple[ReferenceProfileRevision, ...]:
        with self._lock:
            self.get_profile(profile_id)
            return tuple(
                self._revisions[revision_id]
                for revision_id in self._revision_ids_by_profile[profile_id]
            )

    def append_revision(
        self,
        revision: ReferenceProfileRevision,
        *,
        expected_previous_revision_number: int,
        security_audit: SecurityAuditEvent | None = None,
    ) -> ReferenceProfileRevision:
        if not isinstance(revision, ReferenceProfileRevision):
            raise ValueError("revision must be a ReferenceProfileRevision")
        if (
            not isinstance(expected_previous_revision_number, int)
            or isinstance(expected_previous_revision_number, bool)
            or expected_previous_revision_number < 1
        ):
            raise ValueError("expected_previous_revision_number must be positive")
        with self._lock:
            profile = self.get_profile(revision.profile_id)
            revisions = self.list_revisions(profile.profile_id)
            current = revisions[-1]
            if current.revision_number != expected_previous_revision_number:
                raise ReferenceProfileConflictError(
                    "reference revision changed concurrently"
                )
            if revision.revision_number != current.revision_number + 1:
                raise ValueError("revision number must be append-only")
            if revision.created_at < current.created_at:
                raise ValueError("revision cannot precede previous revision")
            if revision.revision_id in self._revisions:
                raise ReferenceProfileConflictError("duplicate revision identity")
            self._revisions[revision.revision_id] = revision
            self._revision_ids_by_profile[profile.profile_id].append(
                revision.revision_id
            )
            if security_audit is not None:
                self._security_audits.append(security_audit)
            return revision

    def set_enabled_revision(
        self,
        profile_id: str,
        revision_id: str | None,
        *,
        security_audit: SecurityAuditEvent | None = None,
    ) -> ReferenceProfile:
        with self._lock:
            profile = self.get_profile(profile_id)
            if revision_id is not None:
                revision = self.get_revision(revision_id)
                if revision.profile_id != profile.profile_id:
                    raise ValueError("enabled revision does not belong to profile")
            updated = replace(profile, enabled_revision_id=revision_id)
            self._profiles[profile_id] = updated
            if security_audit is not None:
                self._security_audits.append(security_audit)
            return updated

    def list_enabled_for_workflow(
        self,
        workflow_id: str,
    ) -> tuple[tuple[ReferenceProfile, ReferenceProfileRevision], ...]:
        with self._lock:
            result = []
            for profile in self.list_profiles():
                if profile.enabled_revision_id is None:
                    continue
                revision = self.get_revision(profile.enabled_revision_id)
                if workflow_id in revision.workflow_ids:
                    result.append((profile, revision))
            return tuple(result)
