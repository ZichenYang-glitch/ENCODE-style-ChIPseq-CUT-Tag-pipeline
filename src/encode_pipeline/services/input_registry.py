"""Workflow-neutral orchestration for approved local Input registry metadata."""

from __future__ import annotations

from collections.abc import Callable, Iterator
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from threading import RLock
from uuid import uuid4

from encode_pipeline.platform.input_registry import (
    InputFile,
    InputFileRevision,
    InputUseBindingEnvelope,
    InputUseBindingPlan,
    ProjectStoragePoolBinding,
    StoragePool,
    build_input_file_revision,
    validate_input_file_stable_key,
)
from encode_pipeline.platform.security_audit import AuditAction
from encode_pipeline.services.admin_security_audit import (
    build_storage_action_event,
)
from encode_pipeline.services.authentication_service import AuthenticationActor
from encode_pipeline.services.input_file_access import FileObservation
from encode_pipeline.services.input_file_access import InputFileAccess
from encode_pipeline.services.input_registry_repositories import (
    InMemoryInputRegistryRepository,
    InputRegistryConflictError,
    InputRegistryRepository,
)
from encode_pipeline.services.private_storage_pools import (
    PrivateStoragePoolConfig,
    load_private_storage_pool_config,
)


@dataclass(frozen=True)
class InputFileCurrentView:
    """Stable InputFile identity paired with its current immutable revision."""

    input_file: InputFile
    revision: InputFileRevision

    def __post_init__(self) -> None:
        if not isinstance(self.input_file, InputFile):
            raise ValueError("input_file must be an InputFile")
        if not isinstance(self.revision, InputFileRevision):
            raise ValueError("revision must be an InputFileRevision")
        if (
            self.revision.input_file_id != self.input_file.input_file_id
            or self.revision.project_id != self.input_file.project_id
            or self.revision.storage_pool_id != self.input_file.storage_pool_id
        ):
            raise ValueError("current revision does not match InputFile")


class InputRegistryService:
    """Local-administrator service for logical pool and regular-file metadata."""

    _MAX_REGISTRATION_ATTEMPTS = 8

    def __init__(
        self,
        *,
        repository: InputRegistryRepository | None = None,
        storage_pool_id_factory: Callable[[], str] | None = None,
        input_file_id_factory: Callable[[], str] | None = None,
        input_file_revision_id_factory: Callable[[], str] | None = None,
        now_factory: Callable[[], datetime] | None = None,
    ) -> None:
        self._repository = (
            repository if repository is not None else InMemoryInputRegistryRepository()
        )
        self._storage_pool_id_factory = (
            storage_pool_id_factory
            if storage_pool_id_factory is not None
            else lambda: f"stgp_{uuid4().hex}"
        )
        self._input_file_id_factory = (
            input_file_id_factory
            if input_file_id_factory is not None
            else lambda: f"inpf_{uuid4().hex}"
        )
        self._input_file_revision_id_factory = (
            input_file_revision_id_factory
            if input_file_revision_id_factory is not None
            else lambda: f"inpfr_{uuid4().hex}"
        )
        self._now_factory = (
            now_factory
            if now_factory is not None
            else lambda: datetime.now(timezone.utc)
        )
        self._lock = RLock()

    def create_storage_pool(
        self,
        *,
        display_name: str,
        config_key: str,
        security_audit_actor: AuthenticationActor | None = None,
    ) -> StoragePool:
        """Create one approved logical pool without accepting a physical root."""
        with self._lock:
            pool = StoragePool(
                storage_pool_id=self._storage_pool_id_factory(),
                config_key=config_key,
                display_name=display_name,
                created_at=self._now_factory(),
            )
            return self._repository.create_storage_pool(
                pool,
                security_audit=(
                    None
                    if security_audit_actor is None
                    else build_storage_action_event(
                        AuditAction.STORAGE_REGISTER,
                        security_audit_actor,
                        pool.storage_pool_id,
                        occurred_at=pool.created_at,
                    )
                ),
            )

    def get_storage_pool(self, storage_pool_id: str) -> StoragePool:
        """Return allowlisted logical pool metadata."""
        with self._lock:
            return self._repository.get_storage_pool(storage_pool_id)

    def list_storage_pools(
        self,
        *,
        include_archived: bool = False,
    ) -> tuple[StoragePool, ...]:
        """Return pools in stable creation order."""
        with self._lock:
            return self._repository.list_storage_pools(
                include_archived=include_archived
            )

    def archive_storage_pool(self, storage_pool_id: str) -> StoragePool:
        """Archive a pool without deleting input or snapshot evidence."""
        with self._lock:
            return self._repository.archive_storage_pool(
                storage_pool_id,
                archived_at=self._now_factory(),
            )

    def bind_project_storage_pool(
        self,
        project_id: str,
        storage_pool_id: str,
    ) -> ProjectStoragePoolBinding:
        """Bind one active Project once to an active approved pool."""
        with self._lock:
            binding = ProjectStoragePoolBinding(
                project_id=project_id,
                storage_pool_id=storage_pool_id,
                bound_at=self._now_factory(),
            )
            return self._repository.create_project_storage_pool_binding(binding)

    def get_project_storage_pool_binding(
        self,
        project_id: str,
    ) -> ProjectStoragePoolBinding:
        """Return the immutable Project-to-pool authorization binding."""
        with self._lock:
            return self._repository.get_project_storage_pool_binding(project_id)

    def get_active_project_storage_pool_binding(
        self,
        project_id: str,
    ) -> ProjectStoragePoolBinding:
        """Return a binding only while its Project and pool remain active."""
        with self._lock:
            return self._repository.get_active_project_storage_pool_binding(project_id)

    def register_input_file(
        self,
        project_id: str,
        *,
        stable_key: str,
        observation: FileObservation,
    ) -> InputFileCurrentView:
        """Create or idempotently append a descriptor-observed file revision."""
        if not isinstance(observation, FileObservation):
            raise ValueError("observation must be a FileObservation")
        with self._lock:
            for _attempt in range(self._MAX_REGISTRATION_ATTEMPTS):
                binding = self._repository.get_active_project_storage_pool_binding(
                    project_id
                )
                try:
                    input_file = self._repository.get_input_file_by_stable_key(
                        project_id,
                        stable_key,
                    )
                except KeyError:
                    input_file = InputFile(
                        input_file_id=self._input_file_id_factory(),
                        project_id=project_id,
                        storage_pool_id=binding.storage_pool_id,
                        stable_key=stable_key,
                        created_at=self._now_factory(),
                    )
                    revision = self._build_revision(
                        input_file,
                        revision_number=1,
                        observation=observation,
                    )
                    try:
                        created = self._repository.create_input_file(
                            input_file,
                            revision,
                        )
                    except InputRegistryConflictError:
                        try:
                            self._repository.get_input_file_by_stable_key(
                                project_id,
                                stable_key,
                            )
                        except KeyError:
                            raise
                        continue
                    return InputFileCurrentView(
                        input_file=created[0],
                        revision=created[1],
                    )

                if input_file.archived_at is not None:
                    raise ValueError("InputFile is archived")
                if (
                    input_file.project_id != project_id
                    or input_file.storage_pool_id != binding.storage_pool_id
                ):
                    raise ValueError(
                        "InputFile does not match the active Project pool binding"
                    )
                revisions = self._repository.list_input_file_revisions(
                    input_file.input_file_id
                )
                if not revisions:
                    raise RuntimeError("InputFile has no immutable revision")
                current = revisions[-1]
                if self._matches_observation(current, observation):
                    return InputFileCurrentView(
                        input_file=input_file,
                        revision=current,
                    )
                revision = self._build_revision(
                    input_file,
                    revision_number=current.revision_number + 1,
                    observation=observation,
                )
                try:
                    appended = self._repository.append_input_file_revision(
                        revision,
                        expected_previous_revision_number=current.revision_number,
                    )
                except InputRegistryConflictError:
                    continue
                return InputFileCurrentView(
                    input_file=input_file,
                    revision=appended,
                )
        raise InputRegistryConflictError(
            "InputFile registration changed concurrently too many times"
        )

    def get_input_file(self, input_file_id: str) -> InputFile:
        """Return one stable InputFile identity."""
        with self._lock:
            return self._repository.get_input_file(input_file_id)

    def list_input_files(
        self,
        project_id: str,
        *,
        include_archived: bool = False,
    ) -> tuple[InputFileCurrentView, ...]:
        """Return InputFiles with their latest committed immutable revisions."""
        with self._lock:
            views: list[InputFileCurrentView] = []
            for input_file in self._repository.list_input_files(
                project_id,
                include_archived=include_archived,
            ):
                revisions = self._repository.list_input_file_revisions(
                    input_file.input_file_id
                )
                if not revisions:
                    raise RuntimeError("InputFile has no immutable revision")
                views.append(
                    InputFileCurrentView(
                        input_file=input_file,
                        revision=revisions[-1],
                    )
                )
            return tuple(views)

    def list_input_file_revisions(
        self,
        input_file_id: str,
    ) -> tuple[InputFileRevision, ...]:
        """Return complete append-only revision history in revision order."""
        with self._lock:
            return self._repository.list_input_file_revisions(input_file_id)

    def archive_input_file(self, input_file_id: str) -> InputFile:
        """Archive a logical file without deleting historical revisions."""
        with self._lock:
            return self._repository.archive_input_file(
                input_file_id,
                archived_at=self._now_factory(),
            )

    def resolve_input_use_binding_plan(
        self,
        plan: InputUseBindingPlan,
        *,
        project_sample_binding_digest: str,
        workflow_inputs_digest: str,
    ) -> InputUseBindingEnvelope:
        """Authorize exact active revisions and return path-free evidence."""
        with self._lock:
            return self._repository.resolve_input_use_binding_plan(
                plan,
                project_sample_binding_digest=project_sample_binding_digest,
                workflow_inputs_digest=workflow_inputs_digest,
            )

    def _build_revision(
        self,
        input_file: InputFile,
        *,
        revision_number: int,
        observation: FileObservation,
    ) -> InputFileRevision:
        return build_input_file_revision(
            input_file_revision_id=self._input_file_revision_id_factory(),
            input_file_id=input_file.input_file_id,
            project_id=input_file.project_id,
            storage_pool_id=input_file.storage_pool_id,
            revision_number=revision_number,
            relative_path=observation.relative_path,
            size_bytes=observation.size_bytes,
            content_sha256=observation.content_sha256,
            created_at=self._now_factory(),
        )

    @staticmethod
    def _matches_observation(
        revision: InputFileRevision,
        observation: FileObservation,
    ) -> bool:
        return (
            revision.relative_path == observation.relative_path
            and revision.size_bytes == observation.size_bytes
            and revision.content_sha256 == observation.content_sha256
        )


class InputRegistryAdminService:
    """Private local CLI facade that observes bytes without returning paths."""

    def __init__(
        self,
        service: InputRegistryService,
        *,
        private_config: PrivateStoragePoolConfig | None,
    ) -> None:
        if not isinstance(service, InputRegistryService):
            raise ValueError("service must be an InputRegistryService")
        if private_config is not None and not isinstance(
            private_config,
            PrivateStoragePoolConfig,
        ):
            raise ValueError("private_config must be a PrivateStoragePoolConfig")
        self._service = service
        self._private_config = private_config

    def register_storage_pool(
        self,
        *,
        display_name: str,
        config_key: str,
        security_audit_actor: AuthenticationActor | None = None,
    ) -> StoragePool:
        """Approve a key only when the private mapping contains its root."""
        config = self._require_private_config()
        config.root_for(config_key)
        return self._service.create_storage_pool(
            display_name=display_name,
            config_key=config_key,
            security_audit_actor=security_audit_actor,
        )

    def bind_project_storage_pool(
        self,
        *,
        project_id: str,
        storage_pool_id: str,
    ) -> ProjectStoragePoolBinding:
        """Bind logical identities without needing to resolve a root."""
        return self._service.bind_project_storage_pool(
            project_id,
            storage_pool_id,
        )

    def register_input_file(
        self,
        *,
        project_id: str,
        stable_key: str,
        pool_relative_path: str,
    ) -> dict[str, object]:
        """Observe a private path and return only path-free logical evidence."""
        validate_input_file_stable_key(stable_key)
        binding = self._service.get_active_project_storage_pool_binding(project_id)
        pool = self._service.get_storage_pool(binding.storage_pool_id)
        access = InputFileAccess.from_config(
            self._require_private_config(),
            pool.config_key,
        )
        observation = access.observe(pool_relative_path)
        registered = self._service.register_input_file(
            project_id,
            stable_key=stable_key,
            observation=observation,
        )
        return {
            "input_file_id": registered.input_file.input_file_id,
            "input_file_revision_id": (registered.revision.input_file_revision_id),
            "revision_number": registered.revision.revision_number,
            "size_bytes": registered.revision.size_bytes,
            "content_sha256": registered.revision.content_sha256,
        }

    def _require_private_config(self) -> PrivateStoragePoolConfig:
        if self._private_config is None:
            raise ValueError("private storage pool configuration is required")
        return self._private_config


@contextmanager
def open_input_registry_admin(
    *,
    database_url: str,
    storage_pool_config: Path | None,
) -> Iterator[InputRegistryAdminService]:
    """Migrate and compose the operator-only Input registry CLI boundary."""
    from encode_pipeline.persistence.database import (
        create_database_engine,
        create_session_factory,
    )
    from encode_pipeline.persistence.input_registry import (
        SqlAlchemyInputRegistryRepository,
    )
    from encode_pipeline.persistence.migrations import upgrade_database
    from encode_pipeline.persistence.runtime import resolve_database_url

    resolved_url = resolve_database_url(database_url)
    upgrade_database(resolved_url)
    engine = create_database_engine(resolved_url)
    try:
        private_config = (
            None
            if storage_pool_config is None
            else load_private_storage_pool_config(storage_pool_config)
        )
        repository = SqlAlchemyInputRegistryRepository(create_session_factory(engine))
        yield InputRegistryAdminService(
            InputRegistryService(repository=repository),
            private_config=private_config,
        )
    finally:
        engine.dispose()
