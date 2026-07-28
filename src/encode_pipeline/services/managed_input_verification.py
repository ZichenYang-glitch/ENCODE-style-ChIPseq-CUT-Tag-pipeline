"""Fail-closed verification of exact managed input revisions before use."""

from __future__ import annotations

from pathlib import Path

from encode_pipeline.platform.adapters import (
    ManagedInputUseValidatingAdapter,
    WorkflowInputs,
)
from encode_pipeline.platform.input_registry import (
    InputFile,
    InputFileRevision,
    InputProvenanceMode,
    InputUseBinding,
    InputUseBindingEnvelope,
    ProjectStoragePoolBinding,
    StoragePool,
)
from encode_pipeline.platform.results import Issue, Result
from encode_pipeline.platform.snapshots import (
    build_workflow_inputs_digest,
    canonical_workflow_inputs_json,
)
from encode_pipeline.services.input_file_access import InputFileAccess
from encode_pipeline.services.input_registry_repositories import (
    InputRegistryRepository,
)
from encode_pipeline.services.private_storage_pools import (
    load_private_storage_pool_config,
)
from encode_pipeline.services.run_repositories import RunRepository


_FAILURE_CODE = "MANAGED_INPUT_VERIFICATION_FAILED"
_FAILURE_MESSAGE = "Managed input evidence could not be verified."
_EXECUTION_UNAVAILABLE_CODE = "MANAGED_INPUT_EXECUTION_UNAVAILABLE"
_EXECUTION_UNAVAILABLE_MESSAGE = "Managed input execution handoff is not qualified."


class ManagedInputVerificationService:
    """Reverify catalog evidence, private bytes, and adapter-owned semantics."""

    __slots__ = (
        "_input_registry_repository",
        "_run_repository",
        "_storage_pool_config_path",
    )

    def __init__(
        self,
        *,
        run_repository: RunRepository,
        input_registry_repository: InputRegistryRepository,
        storage_pool_config_path: Path | None,
    ) -> None:
        if storage_pool_config_path is not None and not isinstance(
            storage_pool_config_path,
            Path,
        ):
            raise ValueError("storage_pool_config_path must be a Path or None")
        self._run_repository = run_repository
        self._input_registry_repository = input_registry_repository
        self._storage_pool_config_path = storage_pool_config_path

    def verify_run(self, run_id: str) -> Result[None]:
        """Admit only input modes whose execution handoff is qualified."""
        try:
            record = self._run_repository.get_run(run_id)
            binding = self._run_repository.get_run_input_use_binding(run_id)
            if (
                not isinstance(binding, InputUseBindingEnvelope)
                or record.workflow_id != binding.workflow_id
            ):
                return _failure()
            if not _managed_uses(binding):
                return Result.success(None)
        except Exception:
            return _failure()
        return _execution_unavailable()

    def verify_binding(
        self,
        inputs: WorkflowInputs,
        binding: InputUseBindingEnvelope,
        adapter: object,
    ) -> Result[None]:
        """Verify one resolved envelope inside snapshot or worker boundaries."""
        try:
            if not isinstance(binding, InputUseBindingEnvelope):
                return _failure()
            managed_uses = _managed_uses(binding)
            if not managed_uses:
                return Result.success(None)
            if not isinstance(inputs, WorkflowInputs):
                return _failure()
            if (
                build_workflow_inputs_digest(canonical_workflow_inputs_json(inputs))
                != binding.workflow_inputs_digest
            ):
                return _failure()
            if not isinstance(adapter, ManagedInputUseValidatingAdapter):
                return _failure()
            metadata = getattr(adapter, "metadata", None)
            if getattr(metadata, "workflow_id", None) != binding.workflow_id:
                return _failure()

            self._verify_catalog_and_bytes(binding, managed_uses)
            adapter_result = adapter.validate_managed_input_uses(inputs, binding)
            if (
                not isinstance(adapter_result, Result)
                or adapter_result.is_failure
                or adapter_result.value is None
            ):
                return _failure()
        except Exception:
            return _failure()
        return Result.success(None)

    def _verify_catalog_and_bytes(
        self,
        binding: InputUseBindingEnvelope,
        managed_uses: tuple[InputUseBinding, ...],
    ) -> None:
        repository = self._input_registry_repository
        authorization = repository.get_project_storage_pool_binding(binding.project_id)
        if (
            not isinstance(authorization, ProjectStoragePoolBinding)
            or authorization.project_id != binding.project_id
        ):
            raise ValueError("invalid authorization")

        resolved: list[tuple[InputFileRevision, StoragePool]] = []
        for input_use in managed_uses:
            for member in input_use.members:
                revision = repository.get_input_file_revision(
                    member.input_file_revision_id
                )
                if not isinstance(revision, InputFileRevision):
                    raise ValueError("invalid revision")
                input_file = repository.get_input_file(revision.input_file_id)
                if not isinstance(input_file, InputFile):
                    raise ValueError("invalid input file")
                pool = repository.get_storage_pool(revision.storage_pool_id)
                if not isinstance(pool, StoragePool):
                    raise ValueError("invalid storage pool")
                if not _catalog_evidence_matches(
                    binding=binding,
                    authorization=authorization,
                    member=member,
                    revision=revision,
                    input_file=input_file,
                    pool=pool,
                ):
                    raise ValueError("catalog evidence differs")
                resolved.append((revision, pool))

        config_path = self._storage_pool_config_path
        if config_path is None:
            raise ValueError("private configuration is unavailable")
        config = load_private_storage_pool_config(config_path)
        for revision, pool in resolved:
            InputFileAccess.from_config(config, pool.config_key).reverify(revision)

    def __repr__(self) -> str:
        return f"{type(self).__name__}(storage_pool_config_path=<redacted>)"


def _managed_uses(
    binding: InputUseBindingEnvelope,
) -> tuple[InputUseBinding, ...]:
    return tuple(
        input_use
        for input_use in binding.input_uses
        if input_use.provenance_mode is InputProvenanceMode.MANAGED_REVISION_V1
    )


def _catalog_evidence_matches(
    *,
    binding: InputUseBindingEnvelope,
    authorization: ProjectStoragePoolBinding,
    member: object,
    revision: InputFileRevision,
    input_file: InputFile,
    pool: StoragePool,
) -> bool:
    return (
        getattr(member, "input_file_revision_id", None)
        == revision.input_file_revision_id
        and getattr(member, "input_file_id", None) == revision.input_file_id
        and getattr(member, "revision_digest", None) == revision.digest
        and getattr(member, "size_bytes", None) == revision.size_bytes
        and getattr(member, "content_sha256", None) == revision.content_sha256
        and revision.input_file_id == input_file.input_file_id
        and revision.project_id == input_file.project_id == binding.project_id
        and revision.storage_pool_id
        == input_file.storage_pool_id
        == pool.storage_pool_id
        == authorization.storage_pool_id
    )


def _failure() -> Result[None]:
    return Result.failure(
        [
            Issue(
                code=_FAILURE_CODE,
                message=_FAILURE_MESSAGE,
                source="input_registry",
                path="input_binding",
            )
        ]
    )


def _execution_unavailable() -> Result[None]:
    return Result.failure(
        [
            Issue(
                code=_EXECUTION_UNAVAILABLE_CODE,
                message=_EXECUTION_UNAVAILABLE_MESSAGE,
                source="input_registry",
                path="input_binding",
            )
        ]
    )
