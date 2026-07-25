"""Workflow-neutral orchestration for Project and Sample registry operations."""

from __future__ import annotations

from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass
from datetime import datetime, timezone
import json
from threading import RLock
from types import MappingProxyType
from typing import Any
from uuid import uuid4

from encode_pipeline.platform.data_registry import (
    LEGACY_PROJECT_ID,
    SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
    Project,
    ProjectKind,
    ProjectSampleBinding,
    ProjectSampleSelection,
    Sample,
    SampleRevision,
    build_legacy_project_sample_binding,
    build_sample_revision_payload_digest,
    canonical_sample_revision_payload,
)
from encode_pipeline.services.data_registry_repositories import (
    DataRegistryRepository,
    InMemoryDataRegistryRepository,
)


def _freeze_json(value: object) -> object:
    if isinstance(value, dict):
        return MappingProxyType(
            {key: _freeze_json(item) for key, item in value.items()}
        )
    if isinstance(value, list):
        return tuple(_freeze_json(item) for item in value)
    return value


@dataclass(frozen=True)
class SampleImportSpec:
    """One validated operator row for an atomic Sample import batch."""

    stable_key: str
    display_name: str
    attributes: Mapping[str, object]

    def __post_init__(self) -> None:
        canonical_payload = canonical_sample_revision_payload(
            display_name=self.display_name,
            attributes=self.attributes,
        )
        decoded = json.loads(canonical_payload)
        object.__setattr__(self, "display_name", decoded["display_name"])
        object.__setattr__(
            self,
            "attributes",
            _freeze_json(decoded["attributes"]),
        )


@dataclass(frozen=True)
class SampleImportResult:
    """Stable Sample identity and the initial immutable revision."""

    sample: Sample
    revision: SampleRevision

    def __post_init__(self) -> None:
        if not isinstance(self.sample, Sample):
            raise ValueError("sample must be a Sample")
        if not isinstance(self.revision, SampleRevision):
            raise ValueError("revision must be a SampleRevision")
        if (
            self.revision.sample_id != self.sample.sample_id
            or self.revision.project_id != self.sample.project_id
            or self.revision.revision_number != 1
        ):
            raise ValueError("initial revision does not match Sample")


@dataclass(frozen=True)
class SampleCurrentView:
    """Read-only projection of a Sample and its current immutable revision."""

    sample: Sample
    revision: SampleRevision

    def __post_init__(self) -> None:
        if not isinstance(self.sample, Sample):
            raise ValueError("sample must be a Sample")
        if not isinstance(self.revision, SampleRevision):
            raise ValueError("revision must be a SampleRevision")
        if (
            self.revision.sample_id != self.sample.sample_id
            or self.revision.project_id != self.sample.project_id
        ):
            raise ValueError("current revision does not match Sample")

    @property
    def display_name(self) -> str:
        """Return display metadata from the current immutable revision."""
        return self.revision.display_name

    @property
    def attributes(self) -> dict[str, Any]:
        """Return a fresh copy of current workflow-neutral attributes."""
        return self.revision.attributes


class DataRegistryService:
    """Application service for local administrator registry operations."""

    def __init__(
        self,
        *,
        repository: DataRegistryRepository | None = None,
        project_id_factory: Callable[[], str] | None = None,
        sample_id_factory: Callable[[], str] | None = None,
        sample_revision_id_factory: Callable[[], str] | None = None,
        now_factory: Callable[[], datetime] | None = None,
    ) -> None:
        self._repository = (
            repository if repository is not None else InMemoryDataRegistryRepository()
        )
        self._project_id_factory = (
            project_id_factory
            if project_id_factory is not None
            else lambda: f"prj_{uuid4().hex}"
        )
        self._sample_id_factory = (
            sample_id_factory
            if sample_id_factory is not None
            else lambda: f"smp_{uuid4().hex}"
        )
        self._sample_revision_id_factory = (
            sample_revision_id_factory
            if sample_revision_id_factory is not None
            else lambda: f"smpr_{uuid4().hex}"
        )
        self._now_factory = (
            now_factory
            if now_factory is not None
            else lambda: datetime.now(timezone.utc)
        )
        self._lock = RLock()

    def create_project(self, display_name: str) -> Project:
        """Create one active user Project with an opaque stable identity."""
        with self._lock:
            project = Project(
                project_id=self._project_id_factory(),
                display_name=display_name,
                kind=ProjectKind.USER,
                created_at=self._now_factory(),
            )
            return self._repository.create_project(project)

    def get_project(self, project_id: str) -> Project:
        """Return one Project without exposing persistence rows."""
        with self._lock:
            return self._repository.get_project(project_id)

    def list_projects(
        self,
        *,
        include_archived: bool = False,
    ) -> tuple[Project, ...]:
        """Return Projects in stable creation order."""
        with self._lock:
            return self._repository.list_projects(include_archived=include_archived)

    def archive_project(self, project_id: str) -> Project:
        """Archive a user Project without deleting historical evidence."""
        with self._lock:
            if project_id == LEGACY_PROJECT_ID:
                raise ValueError("reserved Legacy Project cannot be archived")
            return self._repository.archive_project(
                project_id,
                archived_at=self._now_factory(),
            )

    def import_sample(
        self,
        project_id: str,
        *,
        stable_key: str,
        display_name: str,
        attributes: Mapping[str, object],
    ) -> SampleImportResult:
        """Import one Sample through the same atomic batch boundary as TSV."""
        return self.import_samples(
            project_id,
            rows=(
                SampleImportSpec(
                    stable_key=stable_key,
                    display_name=display_name,
                    attributes=attributes,
                ),
            ),
        )[0]

    def import_samples(
        self,
        project_id: str,
        *,
        rows: Iterable[SampleImportSpec],
    ) -> tuple[SampleImportResult, ...]:
        """Atomically import a non-empty batch of Samples and revision 1."""
        rows = tuple(rows)
        if not rows:
            raise ValueError("sample import batch must not be empty")
        if any(not isinstance(row, SampleImportSpec) for row in rows):
            raise ValueError("rows must contain SampleImportSpec values")
        with self._lock:
            project = self._repository.get_project(project_id)
            self._require_active_user_project(project)
            created_at = self._now_factory()
            entries: list[tuple[Sample, SampleRevision]] = []
            for row in rows:
                sample = Sample(
                    sample_id=self._sample_id_factory(),
                    project_id=project.project_id,
                    stable_key=row.stable_key,
                    created_at=created_at,
                )
                canonical_payload = canonical_sample_revision_payload(
                    display_name=row.display_name,
                    attributes=row.attributes,
                )
                revision = SampleRevision(
                    sample_revision_id=self._sample_revision_id_factory(),
                    project_id=project.project_id,
                    sample_id=sample.sample_id,
                    revision_number=1,
                    canonical_payload=canonical_payload,
                    payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
                    payload_digest=build_sample_revision_payload_digest(
                        canonical_payload
                    ),
                    created_at=created_at,
                )
                entries.append((sample, revision))

            created = self._repository.create_samples(tuple(entries))
            return tuple(
                SampleImportResult(sample=sample, revision=revision)
                for sample, revision in created
            )

    def get_sample(self, sample_id: str) -> Sample:
        """Return one stable Sample identity."""
        with self._lock:
            return self._repository.get_sample(sample_id)

    def list_samples(self, project_id: str) -> tuple[SampleCurrentView, ...]:
        """Return each Sample with its latest committed immutable revision."""
        with self._lock:
            samples = self._repository.list_samples(project_id)
            views: list[SampleCurrentView] = []
            for sample in samples:
                revisions = self._repository.list_sample_revisions(sample.sample_id)
                if not revisions:
                    raise RuntimeError("Sample has no immutable revision")
                views.append(SampleCurrentView(sample=sample, revision=revisions[-1]))
            return tuple(views)

    def list_sample_revisions(
        self,
        sample_id: str,
    ) -> tuple[SampleRevision, ...]:
        """Return complete append-only revision history in revision order."""
        with self._lock:
            return self._repository.list_sample_revisions(sample_id)

    def revise_sample(
        self,
        sample_id: str,
        *,
        display_name: str,
        attributes: Mapping[str, object],
    ) -> SampleRevision:
        """Append a complete new revision without mutating Sample identity."""
        with self._lock:
            sample = self._repository.get_sample(sample_id)
            project = self._repository.get_project(sample.project_id)
            self._require_active_user_project(project)
            revisions = self._repository.list_sample_revisions(sample_id)
            if not revisions:
                raise RuntimeError("Sample has no immutable revision")
            previous = revisions[-1]
            canonical_payload = canonical_sample_revision_payload(
                display_name=display_name,
                attributes=attributes,
            )
            revision = SampleRevision(
                sample_revision_id=self._sample_revision_id_factory(),
                project_id=sample.project_id,
                sample_id=sample.sample_id,
                revision_number=previous.revision_number + 1,
                canonical_payload=canonical_payload,
                payload_digest_scheme=SAMPLE_REVISION_PAYLOAD_DIGEST_SCHEME,
                payload_digest=build_sample_revision_payload_digest(canonical_payload),
                created_at=self._now_factory(),
            )
            return self._repository.append_sample_revision(
                revision,
                expected_previous_revision_number=previous.revision_number,
            )

    def resolve_project_sample_selection(
        self,
        selection: ProjectSampleSelection,
        *,
        workflow_inputs_digest: str,
    ) -> ProjectSampleBinding:
        """Atomically resolve exact active revisions and freeze their digests."""
        with self._lock:
            return self._repository.resolve_project_sample_selection(
                selection,
                workflow_inputs_digest=workflow_inputs_digest,
            )

    def build_legacy_binding(
        self,
        *,
        workflow_inputs_digest: str,
    ) -> ProjectSampleBinding:
        """Return conservative unresolved binding evidence for v1 inputs."""
        return build_legacy_project_sample_binding(workflow_inputs_digest)

    @staticmethod
    def _require_active_user_project(project: Project) -> None:
        if project.project_id == LEGACY_PROJECT_ID:
            raise ValueError("Legacy Project cannot receive trusted Samples")
        if project.archived_at is not None:
            raise ValueError("Project is archived")
