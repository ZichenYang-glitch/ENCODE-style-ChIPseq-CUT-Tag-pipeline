"""Persistence boundary for the workflow-neutral Project/Sample registry."""

from __future__ import annotations

from datetime import datetime, timezone
from dataclasses import replace
from threading import RLock
from typing import Protocol

from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    BindingMode,
    BindingProvenance,
    Project,
    ProjectSampleBinding,
    ProjectSampleSelection,
    Sample,
    SampleRevision,
    SampleRevisionBindingRef,
    build_legacy_project,
    build_project_sample_binding,
)


class DataRegistryConflictError(RuntimeError):
    """A registry identity or append-only lifecycle invariant conflicted."""


class DataRegistryRepository(Protocol):
    """Atomic storage contract consumed by :class:`DataRegistryService`."""

    def create_project(self, project: Project) -> Project: ...

    def get_project(self, project_id: str) -> Project: ...

    def list_projects(
        self,
        *,
        include_archived: bool = False,
    ) -> tuple[Project, ...]: ...

    def archive_project(
        self,
        project_id: str,
        *,
        archived_at: datetime,
    ) -> Project: ...

    def create_sample(
        self,
        sample: Sample,
        initial_revision: SampleRevision,
    ) -> tuple[Sample, SampleRevision]: ...

    def create_samples(
        self,
        entries: tuple[tuple[Sample, SampleRevision], ...],
    ) -> tuple[tuple[Sample, SampleRevision], ...]: ...

    def get_sample(self, sample_id: str) -> Sample: ...

    def list_samples(self, project_id: str) -> tuple[Sample, ...]: ...

    def get_sample_revision(self, sample_revision_id: str) -> SampleRevision: ...

    def list_sample_revisions(
        self,
        sample_id: str,
    ) -> tuple[SampleRevision, ...]: ...

    def append_sample_revision(
        self,
        revision: SampleRevision,
        *,
        expected_previous_revision_number: int,
    ) -> SampleRevision: ...

    def resolve_project_sample_selection(
        self,
        selection: ProjectSampleSelection,
        *,
        workflow_inputs_digest: str,
    ) -> ProjectSampleBinding: ...


class InMemoryDataRegistryRepository:
    """Lock-protected in-memory implementation with SQL-equivalent invariants."""

    def __init__(self, *, legacy_created_at: datetime | None = None) -> None:
        if legacy_created_at is None:
            legacy_created_at = datetime.now(timezone.utc)
        legacy = build_legacy_project(legacy_created_at)
        self._projects: dict[str, Project] = {legacy.project_id: legacy}
        self._samples: dict[str, Sample] = {}
        self._sample_ids_by_project: dict[str, list[str]] = {}
        self._sample_id_by_stable_key: dict[tuple[str, str], str] = {}
        self._revisions: dict[str, SampleRevision] = {}
        self._revision_ids_by_sample: dict[str, list[str]] = {}
        self._lock = RLock()

    def create_project(self, project: Project) -> Project:
        if not isinstance(project, Project):
            raise ValueError("project must be a Project")
        with self._lock:
            if project.project_id in self._projects:
                raise DataRegistryConflictError(
                    f"Duplicate project_id: {project.project_id!r}"
                )
            if project.project_id == LEGACY_PROJECT_ID:
                raise DataRegistryConflictError(
                    "reserved Legacy Project already exists"
                )
            self._projects[project.project_id] = project
            return project

    def get_project(self, project_id: str) -> Project:
        with self._lock:
            try:
                return self._projects[project_id]
            except KeyError:
                raise KeyError(f"Unknown project_id: {project_id!r}") from None

    def list_projects(
        self,
        *,
        include_archived: bool = False,
    ) -> tuple[Project, ...]:
        with self._lock:
            projects = tuple(self._projects.values())
            if include_archived:
                return projects
            return tuple(project for project in projects if project.archived_at is None)

    def archive_project(
        self,
        project_id: str,
        *,
        archived_at: datetime,
    ) -> Project:
        with self._lock:
            project = self.get_project(project_id)
            if project_id == LEGACY_PROJECT_ID:
                raise ValueError("reserved Legacy Project cannot be archived")
            if project.archived_at is not None:
                raise DataRegistryConflictError("Project is already archived")
            archived = replace(project, archived_at=archived_at)
            self._projects[project_id] = archived
            return archived

    def create_sample(
        self,
        sample: Sample,
        initial_revision: SampleRevision,
    ) -> tuple[Sample, SampleRevision]:
        return self.create_samples(((sample, initial_revision),))[0]

    def create_samples(
        self,
        entries: tuple[tuple[Sample, SampleRevision], ...],
    ) -> tuple[tuple[Sample, SampleRevision], ...]:
        entries = tuple(entries)
        with self._lock:
            pending_sample_ids: set[str] = set()
            pending_revision_ids: set[str] = set()
            pending_stable_identities: set[tuple[str, str]] = set()
            for sample, initial_revision in entries:
                if not isinstance(sample, Sample):
                    raise ValueError("sample must be a Sample")
                if not isinstance(initial_revision, SampleRevision):
                    raise ValueError("initial_revision must be a SampleRevision")
                project = self.get_project(sample.project_id)
                self._require_active_user_project(project)
                if (
                    sample.sample_id in self._samples
                    or sample.sample_id in pending_sample_ids
                ):
                    raise DataRegistryConflictError(
                        f"Duplicate sample_id: {sample.sample_id!r}"
                    )
                stable_identity = (sample.project_id, sample.stable_key)
                if (
                    stable_identity in self._sample_id_by_stable_key
                    or stable_identity in pending_stable_identities
                ):
                    raise DataRegistryConflictError(
                        "Sample stable_key already exists in this Project"
                    )
                if (
                    initial_revision.sample_revision_id in self._revisions
                    or initial_revision.sample_revision_id in pending_revision_ids
                ):
                    raise DataRegistryConflictError("Duplicate sample_revision_id")
                if (
                    initial_revision.project_id != sample.project_id
                    or initial_revision.sample_id != sample.sample_id
                    or initial_revision.revision_number != 1
                ):
                    raise ValueError(
                        "initial revision must be revision 1 of the new Sample"
                    )
                if initial_revision.created_at < sample.created_at:
                    raise ValueError("initial revision cannot precede Sample creation")
                pending_sample_ids.add(sample.sample_id)
                pending_revision_ids.add(initial_revision.sample_revision_id)
                pending_stable_identities.add(stable_identity)

            for sample, initial_revision in entries:
                self._samples[sample.sample_id] = sample
                self._sample_ids_by_project.setdefault(
                    sample.project_id,
                    [],
                ).append(sample.sample_id)
                self._sample_id_by_stable_key[
                    (sample.project_id, sample.stable_key)
                ] = sample.sample_id
                self._revisions[initial_revision.sample_revision_id] = initial_revision
                self._revision_ids_by_sample[sample.sample_id] = [
                    initial_revision.sample_revision_id
                ]
            return entries

    def get_sample(self, sample_id: str) -> Sample:
        with self._lock:
            try:
                return self._samples[sample_id]
            except KeyError:
                raise KeyError(f"Unknown sample_id: {sample_id!r}") from None

    def list_samples(self, project_id: str) -> tuple[Sample, ...]:
        with self._lock:
            self.get_project(project_id)
            return tuple(
                self._samples[sample_id]
                for sample_id in self._sample_ids_by_project.get(project_id, ())
            )

    def get_sample_revision(self, sample_revision_id: str) -> SampleRevision:
        with self._lock:
            try:
                return self._revisions[sample_revision_id]
            except KeyError:
                raise KeyError(
                    f"Unknown sample_revision_id: {sample_revision_id!r}"
                ) from None

    def list_sample_revisions(
        self,
        sample_id: str,
    ) -> tuple[SampleRevision, ...]:
        with self._lock:
            self.get_sample(sample_id)
            return tuple(
                self._revisions[revision_id]
                for revision_id in self._revision_ids_by_sample[sample_id]
            )

    def append_sample_revision(
        self,
        revision: SampleRevision,
        *,
        expected_previous_revision_number: int,
    ) -> SampleRevision:
        if not isinstance(revision, SampleRevision):
            raise ValueError("revision must be a SampleRevision")
        if (
            not isinstance(expected_previous_revision_number, int)
            or isinstance(expected_previous_revision_number, bool)
            or expected_previous_revision_number < 1
        ):
            raise ValueError("expected_previous_revision_number must be positive")
        with self._lock:
            sample = self.get_sample(revision.sample_id)
            project = self.get_project(sample.project_id)
            self._require_active_user_project(project)
            if (
                revision.project_id != sample.project_id
                or revision.sample_id != sample.sample_id
            ):
                raise ValueError("revision must belong to the existing Sample")
            if revision.sample_revision_id in self._revisions:
                raise DataRegistryConflictError("Duplicate sample_revision_id")
            revision_ids = self._revision_ids_by_sample[sample.sample_id]
            previous = self._revisions[revision_ids[-1]]
            if previous.revision_number != expected_previous_revision_number:
                raise DataRegistryConflictError("Sample revision changed concurrently")
            if revision.revision_number != previous.revision_number + 1:
                raise ValueError("Sample revision number must append exactly once")
            if revision.created_at < previous.created_at:
                raise ValueError("new Sample revision cannot precede history")

            self._revisions[revision.sample_revision_id] = revision
            revision_ids.append(revision.sample_revision_id)
            return revision

    def resolve_project_sample_selection(
        self,
        selection: ProjectSampleSelection,
        *,
        workflow_inputs_digest: str,
    ) -> ProjectSampleBinding:
        if not isinstance(selection, ProjectSampleSelection):
            raise ValueError("selection must be a ProjectSampleSelection")
        with self._lock:
            project = self.get_project(selection.project_id)
            self._require_active_user_project(project)
            revisions: list[SampleRevision] = []
            for revision_id in selection.sample_revision_ids:
                revision = self.get_sample_revision(revision_id)
                sample = self.get_sample(revision.sample_id)
                if (
                    revision.project_id != selection.project_id
                    or sample.project_id != selection.project_id
                ):
                    raise ValueError(
                        "all SampleRevisions must belong to the same Project"
                    )
                revisions.append(revision)
            return build_project_sample_binding(
                project_id=selection.project_id,
                binding_mode=BindingMode.BOUND_V1,
                provenance=BindingProvenance.RESOLVED,
                workflow_inputs_digest=workflow_inputs_digest,
                sample_revisions=tuple(
                    SampleRevisionBindingRef(
                        sample_revision_id=revision.sample_revision_id,
                        payload_digest=revision.payload_digest,
                    )
                    for revision in revisions
                ),
            )

    @staticmethod
    def _require_active_user_project(project: Project) -> None:
        if project.project_id == LEGACY_PROJECT_ID:
            raise ValueError("Legacy Project cannot receive trusted registry data")
        if project.archived_at is not None:
            raise ValueError("Project is archived")
