"""Persistence boundary for approved StoragePools and immutable InputFiles."""

from __future__ import annotations

from collections.abc import Iterator
from contextlib import AbstractContextManager, contextmanager
from dataclasses import replace
from datetime import datetime
from threading import RLock
from typing import Protocol

from encode_pipeline.platform.data_registry import LEGACY_PROJECT_ID, Project
from encode_pipeline.platform.input_registry import (
    InputFile,
    InputFileRevision,
    InputFileRevisionBindingRef,
    InputProvenanceMode,
    InputUseBindingEnvelope,
    InputUseBindingPlan,
    ProjectStoragePoolBinding,
    StoragePool,
    build_input_use_binding,
    build_input_use_binding_envelope,
)
from encode_pipeline.services.data_registry_repositories import (
    DataRegistryRepository,
    InMemoryDataRegistryRepository,
)


class InputRegistryConflictError(RuntimeError):
    """A stable identity or append-only Input registry invariant conflicted."""


class InputRegistryRepository(Protocol):
    """Atomic persistence contract consumed by :class:`InputRegistryService`."""

    def create_storage_pool(self, storage_pool: StoragePool) -> StoragePool: ...

    def get_storage_pool(self, storage_pool_id: str) -> StoragePool: ...

    def list_storage_pools(
        self,
        *,
        include_archived: bool = False,
    ) -> tuple[StoragePool, ...]: ...

    def archive_storage_pool(
        self,
        storage_pool_id: str,
        *,
        archived_at: datetime,
    ) -> StoragePool: ...

    def create_project_storage_pool_binding(
        self,
        binding: ProjectStoragePoolBinding,
    ) -> ProjectStoragePoolBinding: ...

    def get_project_storage_pool_binding(
        self,
        project_id: str,
    ) -> ProjectStoragePoolBinding: ...

    def get_active_project_storage_pool_binding(
        self,
        project_id: str,
    ) -> ProjectStoragePoolBinding: ...

    def create_input_file(
        self,
        input_file: InputFile,
        initial_revision: InputFileRevision,
    ) -> tuple[InputFile, InputFileRevision]: ...

    def get_input_file(self, input_file_id: str) -> InputFile: ...

    def get_input_file_by_stable_key(
        self,
        project_id: str,
        stable_key: str,
    ) -> InputFile: ...

    def list_input_files(
        self,
        project_id: str,
        *,
        include_archived: bool = False,
    ) -> tuple[InputFile, ...]: ...

    def archive_input_file(
        self,
        input_file_id: str,
        *,
        archived_at: datetime,
    ) -> InputFile: ...

    def get_input_file_revision(
        self,
        input_file_revision_id: str,
    ) -> InputFileRevision: ...

    def list_input_file_revisions(
        self,
        input_file_id: str,
    ) -> tuple[InputFileRevision, ...]: ...

    def append_input_file_revision(
        self,
        revision: InputFileRevision,
        *,
        expected_previous_revision_number: int,
    ) -> InputFileRevision: ...

    def resolve_input_use_binding_plan(
        self,
        plan: InputUseBindingPlan,
        *,
        project_sample_binding_digest: str,
        workflow_inputs_digest: str,
    ) -> InputUseBindingEnvelope: ...

    def locked_input_use_binding_plan(
        self,
        plan: InputUseBindingPlan,
        *,
        project_sample_binding_digest: str,
        workflow_inputs_digest: str,
    ) -> AbstractContextManager[InputUseBindingEnvelope]: ...


class InMemoryInputRegistryRepository:
    """Lock-protected Input registry with SQL-equivalent lifecycle checks."""

    def __init__(
        self,
        *,
        project_repository: DataRegistryRepository | None = None,
    ) -> None:
        self._project_repository = (
            project_repository
            if project_repository is not None
            else InMemoryDataRegistryRepository()
        )
        self._storage_pools: dict[str, StoragePool] = {}
        self._storage_pool_id_by_config_key: dict[str, str] = {}
        self._pool_bindings: dict[str, ProjectStoragePoolBinding] = {}
        self._input_files: dict[str, InputFile] = {}
        self._input_file_ids_by_project: dict[str, list[str]] = {}
        self._input_file_id_by_stable_key: dict[tuple[str, str], str] = {}
        self._revisions: dict[str, InputFileRevision] = {}
        self._revision_ids_by_input_file: dict[str, list[str]] = {}
        self._lock = RLock()

    def create_storage_pool(self, storage_pool: StoragePool) -> StoragePool:
        if not isinstance(storage_pool, StoragePool):
            raise ValueError("storage_pool must be a StoragePool")
        with self._lock:
            if storage_pool.storage_pool_id in self._storage_pools:
                raise InputRegistryConflictError(
                    f"Duplicate storage_pool_id: {storage_pool.storage_pool_id!r}"
                )
            if storage_pool.config_key in self._storage_pool_id_by_config_key:
                raise InputRegistryConflictError(
                    "StoragePool config_key already exists"
                )
            self._storage_pools[storage_pool.storage_pool_id] = storage_pool
            self._storage_pool_id_by_config_key[storage_pool.config_key] = (
                storage_pool.storage_pool_id
            )
            return storage_pool

    def get_storage_pool(self, storage_pool_id: str) -> StoragePool:
        with self._lock:
            try:
                return self._storage_pools[storage_pool_id]
            except KeyError:
                raise KeyError(
                    f"Unknown storage_pool_id: {storage_pool_id!r}"
                ) from None

    def list_storage_pools(
        self,
        *,
        include_archived: bool = False,
    ) -> tuple[StoragePool, ...]:
        with self._lock:
            pools = tuple(
                sorted(
                    self._storage_pools.values(),
                    key=lambda pool: (pool.created_at, pool.storage_pool_id),
                )
            )
            if include_archived:
                return pools
            return tuple(pool for pool in pools if pool.archived_at is None)

    def archive_storage_pool(
        self,
        storage_pool_id: str,
        *,
        archived_at: datetime,
    ) -> StoragePool:
        with self._lock:
            current = self.get_storage_pool(storage_pool_id)
            if current.archived_at is not None:
                raise InputRegistryConflictError("StoragePool is already archived")
            archived = replace(current, archived_at=archived_at)
            self._storage_pools[storage_pool_id] = archived
            return archived

    def create_project_storage_pool_binding(
        self,
        binding: ProjectStoragePoolBinding,
    ) -> ProjectStoragePoolBinding:
        if not isinstance(binding, ProjectStoragePoolBinding):
            raise ValueError("binding must be a ProjectStoragePoolBinding")
        with self._lock:
            project = self._project_repository.get_project(binding.project_id)
            self._require_active_user_project(project)
            pool = self.get_storage_pool(binding.storage_pool_id)
            self._require_active_storage_pool(pool)
            if (
                binding.bound_at < project.created_at
                or binding.bound_at < pool.created_at
            ):
                raise ValueError("Project StoragePool binding cannot precede its scope")
            if binding.project_id in self._pool_bindings:
                raise InputRegistryConflictError(
                    "Project is already bound to a StoragePool"
                )
            self._pool_bindings[binding.project_id] = binding
            return binding

    def get_project_storage_pool_binding(
        self,
        project_id: str,
    ) -> ProjectStoragePoolBinding:
        with self._lock:
            try:
                return self._pool_bindings[project_id]
            except KeyError:
                raise KeyError(
                    f"Unknown Project StoragePool binding: {project_id!r}"
                ) from None

    def get_active_project_storage_pool_binding(
        self,
        project_id: str,
    ) -> ProjectStoragePoolBinding:
        """Return one binding only while its Project and pool are active."""
        with self._lock:
            project = self._project_repository.get_project(project_id)
            self._require_active_user_project(project)
            binding = self.get_project_storage_pool_binding(project_id)
            pool = self.get_storage_pool(binding.storage_pool_id)
            self._require_active_storage_pool(pool)
            return binding

    def create_input_file(
        self,
        input_file: InputFile,
        initial_revision: InputFileRevision,
    ) -> tuple[InputFile, InputFileRevision]:
        if not isinstance(input_file, InputFile):
            raise ValueError("input_file must be an InputFile")
        if not isinstance(initial_revision, InputFileRevision):
            raise ValueError("initial_revision must be an InputFileRevision")
        self._validate_initial_revision(input_file, initial_revision)
        with self._lock:
            self._require_active_registration_scope(input_file)
            if input_file.input_file_id in self._input_files:
                raise InputRegistryConflictError(
                    f"Duplicate input_file_id: {input_file.input_file_id!r}"
                )
            stable_identity = (input_file.project_id, input_file.stable_key)
            if stable_identity in self._input_file_id_by_stable_key:
                raise InputRegistryConflictError(
                    "InputFile stable_key already exists in this Project"
                )
            if initial_revision.input_file_revision_id in self._revisions:
                raise InputRegistryConflictError("Duplicate input_file_revision_id")

            self._input_files[input_file.input_file_id] = input_file
            self._input_file_ids_by_project.setdefault(
                input_file.project_id,
                [],
            ).append(input_file.input_file_id)
            self._input_file_id_by_stable_key[stable_identity] = (
                input_file.input_file_id
            )
            self._revisions[initial_revision.input_file_revision_id] = initial_revision
            self._revision_ids_by_input_file[input_file.input_file_id] = [
                initial_revision.input_file_revision_id
            ]
            return input_file, initial_revision

    def get_input_file(self, input_file_id: str) -> InputFile:
        with self._lock:
            try:
                return self._input_files[input_file_id]
            except KeyError:
                raise KeyError(f"Unknown input_file_id: {input_file_id!r}") from None

    def get_input_file_by_stable_key(
        self,
        project_id: str,
        stable_key: str,
    ) -> InputFile:
        with self._lock:
            try:
                input_file_id = self._input_file_id_by_stable_key[
                    (project_id, stable_key)
                ]
            except KeyError:
                raise KeyError(
                    "Unknown InputFile stable_key for Project: "
                    f"{project_id!r}, {stable_key!r}"
                ) from None
            return self._input_files[input_file_id]

    def list_input_files(
        self,
        project_id: str,
        *,
        include_archived: bool = False,
    ) -> tuple[InputFile, ...]:
        with self._lock:
            self._project_repository.get_project(project_id)
            input_files = tuple(
                sorted(
                    (
                        self._input_files[input_file_id]
                        for input_file_id in self._input_file_ids_by_project.get(
                            project_id,
                            (),
                        )
                    ),
                    key=lambda item: (item.created_at, item.input_file_id),
                )
            )
            if include_archived:
                return input_files
            return tuple(
                input_file
                for input_file in input_files
                if input_file.archived_at is None
            )

    def archive_input_file(
        self,
        input_file_id: str,
        *,
        archived_at: datetime,
    ) -> InputFile:
        with self._lock:
            current = self.get_input_file(input_file_id)
            project = self._project_repository.get_project(current.project_id)
            self._require_active_user_project(project)
            if current.archived_at is not None:
                raise InputRegistryConflictError("InputFile is already archived")
            archived = replace(current, archived_at=archived_at)
            self._input_files[input_file_id] = archived
            return archived

    def get_input_file_revision(
        self,
        input_file_revision_id: str,
    ) -> InputFileRevision:
        with self._lock:
            try:
                return self._revisions[input_file_revision_id]
            except KeyError:
                raise KeyError(
                    f"Unknown input_file_revision_id: {input_file_revision_id!r}"
                ) from None

    def list_input_file_revisions(
        self,
        input_file_id: str,
    ) -> tuple[InputFileRevision, ...]:
        with self._lock:
            self.get_input_file(input_file_id)
            return tuple(
                self._revisions[revision_id]
                for revision_id in self._revision_ids_by_input_file[input_file_id]
            )

    def append_input_file_revision(
        self,
        revision: InputFileRevision,
        *,
        expected_previous_revision_number: int,
    ) -> InputFileRevision:
        if not isinstance(revision, InputFileRevision):
            raise ValueError("revision must be an InputFileRevision")
        if (
            not isinstance(expected_previous_revision_number, int)
            or isinstance(expected_previous_revision_number, bool)
            or expected_previous_revision_number < 1
        ):
            raise ValueError("expected_previous_revision_number must be positive")
        with self._lock:
            input_file = self.get_input_file(revision.input_file_id)
            self._require_active_registration_scope(input_file)
            if input_file.archived_at is not None:
                raise ValueError("InputFile is archived")
            if revision.input_file_revision_id in self._revisions:
                raise InputRegistryConflictError("Duplicate input_file_revision_id")
            if (
                revision.project_id != input_file.project_id
                or revision.storage_pool_id != input_file.storage_pool_id
            ):
                raise ValueError("revision must belong to the existing InputFile scope")
            revision_ids = self._revision_ids_by_input_file[input_file.input_file_id]
            previous = self._revisions[revision_ids[-1]]
            if previous.revision_number != expected_previous_revision_number:
                raise InputRegistryConflictError(
                    "InputFile revision changed concurrently"
                )
            if revision.revision_number != previous.revision_number + 1:
                raise ValueError("InputFile revision number must append exactly once")
            if revision.created_at < previous.created_at:
                raise ValueError("new InputFile revision cannot precede history")
            self._revisions[revision.input_file_revision_id] = revision
            revision_ids.append(revision.input_file_revision_id)
            return revision

    def resolve_input_use_binding_plan(
        self,
        plan: InputUseBindingPlan,
        *,
        project_sample_binding_digest: str,
        workflow_inputs_digest: str,
    ) -> InputUseBindingEnvelope:
        """Resolve exact active revisions to path-free immutable evidence."""
        if not isinstance(plan, InputUseBindingPlan):
            raise ValueError("plan must be an InputUseBindingPlan")
        with self._lock:
            project = self._project_repository.get_project(plan.project_id)
            self._require_active_user_project(project)
            binding: ProjectStoragePoolBinding | None = None
            if any(
                planned_use.provenance_mode is InputProvenanceMode.MANAGED_REVISION_V1
                for planned_use in plan.input_uses
            ):
                binding = self.get_project_storage_pool_binding(plan.project_id)
                pool = self.get_storage_pool(binding.storage_pool_id)
                self._require_active_storage_pool(pool)
            resolved_uses = []
            for planned_use in plan.input_uses:
                if (
                    planned_use.provenance_mode
                    is InputProvenanceMode.TRANSITIONAL_UNMANAGED_V1
                ):
                    resolved_uses.append(
                        build_input_use_binding(planned_use, members=())
                    )
                    continue
                assert binding is not None
                if planned_use.closure_contract_version != "regular_file_v1":
                    raise ValueError(
                        "managed closure contract is unsupported by this registry"
                    )
                members = []
                for revision_id in planned_use.input_file_revision_ids:
                    revision = self.get_input_file_revision(revision_id)
                    input_file = self.get_input_file(revision.input_file_id)
                    self._require_active_binding_member(
                        plan_project_id=plan.project_id,
                        binding=binding,
                        input_file=input_file,
                        revision=revision,
                    )
                    members.append(_input_file_revision_binding_ref(revision))
                resolved_uses.append(
                    build_input_use_binding(
                        planned_use,
                        members=tuple(members),
                    )
                )
            return build_input_use_binding_envelope(
                project_id=plan.project_id,
                project_sample_binding_digest=project_sample_binding_digest,
                workflow_id=plan.workflow_id,
                adapter_contract_version=plan.adapter_contract_version,
                workflow_inputs_digest=workflow_inputs_digest,
                input_uses=tuple(resolved_uses),
            )

    @contextmanager
    def locked_input_use_binding_plan(
        self,
        plan: InputUseBindingPlan,
        *,
        project_sample_binding_digest: str,
        workflow_inputs_digest: str,
    ) -> Iterator[InputUseBindingEnvelope]:
        """Hold the registry lock until the caller stores snapshot evidence."""
        with self._lock:
            yield self.resolve_input_use_binding_plan(
                plan,
                project_sample_binding_digest=project_sample_binding_digest,
                workflow_inputs_digest=workflow_inputs_digest,
            )

    def _require_active_registration_scope(self, input_file: InputFile) -> None:
        project = self._project_repository.get_project(input_file.project_id)
        self._require_active_user_project(project)
        pool = self.get_storage_pool(input_file.storage_pool_id)
        self._require_active_storage_pool(pool)
        binding = self.get_project_storage_pool_binding(input_file.project_id)
        if binding.storage_pool_id != input_file.storage_pool_id:
            raise ValueError("InputFile StoragePool does not match the Project binding")
        if (
            input_file.created_at < project.created_at
            or input_file.created_at < pool.created_at
            or input_file.created_at < binding.bound_at
        ):
            raise ValueError("InputFile cannot precede its authorization scope")

    def _require_active_binding_member(
        self,
        *,
        plan_project_id: str,
        binding: ProjectStoragePoolBinding,
        input_file: InputFile,
        revision: InputFileRevision,
    ) -> None:
        if input_file.archived_at is not None:
            raise ValueError("InputFile is archived")
        if (
            input_file.project_id != plan_project_id
            or input_file.storage_pool_id != binding.storage_pool_id
            or revision.project_id != plan_project_id
            or revision.storage_pool_id != binding.storage_pool_id
            or revision.input_file_id != input_file.input_file_id
        ):
            raise ValueError(
                "InputFile revision is not authorized by the Project pool binding"
            )

    @staticmethod
    def _validate_initial_revision(
        input_file: InputFile,
        revision: InputFileRevision,
    ) -> None:
        if (
            revision.input_file_id != input_file.input_file_id
            or revision.project_id != input_file.project_id
            or revision.storage_pool_id != input_file.storage_pool_id
            or revision.revision_number != 1
        ):
            raise ValueError("initial revision must be revision 1 of the new InputFile")
        if revision.created_at < input_file.created_at:
            raise ValueError("initial revision cannot precede InputFile creation")

    @staticmethod
    def _require_active_user_project(project: Project) -> None:
        if project.project_id == LEGACY_PROJECT_ID:
            raise ValueError("Legacy Project cannot receive trusted InputFiles")
        if project.archived_at is not None:
            raise ValueError("Project is archived")

    @staticmethod
    def _require_active_storage_pool(storage_pool: StoragePool) -> None:
        if storage_pool.archived_at is not None:
            raise ValueError("StoragePool is archived")


def _input_file_revision_binding_ref(
    revision: InputFileRevision,
) -> InputFileRevisionBindingRef:
    """Map one regular file to a fixed generic member coordinate."""
    return InputFileRevisionBindingRef(
        logical_member_key="file",
        input_file_id=revision.input_file_id,
        input_file_revision_id=revision.input_file_revision_id,
        revision_digest=revision.digest,
        size_bytes=revision.size_bytes,
        content_sha256=revision.content_sha256,
    )
