"""Seal one exact track into the operator archive.

Source paths come from a registered revision or a published run artifact.
Callers name those identities. They do not name a filename, and a later
revision is not substituted for the one they named. The copy is what lets
the original scratch or pool file be removed under a later retention step.
"""

from __future__ import annotations

from collections.abc import Mapping
from contextlib import suppress
import os
from pathlib import Path
import re
import stat
import tempfile
from types import SimpleNamespace

from encode_pipeline.adapters.browser.scaffolds import viewer_scaffold
from encode_pipeline.persistence.models import (
    ExperimentRevisionRow,
    ProjectRow,
    ReferenceProfileRevisionRow,
    RunReferenceBindingRow,
)
from encode_pipeline.platform.experiment_catalog import CatalogError
from encode_pipeline.platform.registered_results import (
    ParentEdge,
    RegisteredResultPayloadError,
    build_revision_digest,
    metadata_from_canonical_json,
    validate_project_id,
    validate_result_id,
    validate_result_revision_id,
)
from encode_pipeline.platform.result_generations import (
    validate_artifact_generation,
    validate_artifact_revision,
)
from encode_pipeline.services.artifact_downloads import ArtifactDownloadStreamError
from encode_pipeline.services.private_storage_pools import PrivateStoragePoolConfig
from encode_pipeline.services.registered_result_bytes import (
    RegisteredResultBytesError,
    discard_private_copy,
    seal_private_copy,
)
from encode_pipeline.services.registered_results import open_verified_copy


_CODES = frozenset(
    {
        "archive_version_not_found",
        "archive_version_conflict",
        "archive_offline",
        "archive_integrity",
        "archive_invalid",
    }
)
_TRACK_TYPES = frozenset({"bigwig", "bed", "longrange"})
_TRACK_KEY = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:-]{0,127}$")
_EXPERIMENT_ID = re.compile(r"^exp_[0-9a-f]{32}$")
_EXPERIMENT_REVISION_ID = re.compile(r"^expr_[0-9a-f]{32}$")
_SAFE_RUN_ID = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$")
_SAFE_ARTIFACT_ID = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]{0,127}$")
_SHA256 = re.compile(r"^[0-9a-f]{64}$")
_RECEIPT_FIELDS = frozenset(
    {
        "archive_version_id",
        "project_id",
        "experiment_id",
        "experiment_revision_id",
        "track_key",
        "version_number",
        "display_name",
        "assembly",
        "track_type",
        "source",
        "has_index",
        "members",
    }
)
_REGISTERED_SOURCE = frozenset(
    {
        "kind",
        "result_id",
        "result_revision_id",
        "index_result_id",
        "index_result_revision_id",
    }
)
_RUN_SOURCE = frozenset(
    {
        "kind",
        "run_id",
        "artifact_id",
        "generation",
        "revision",
        "index_artifact_id",
        "index_generation",
        "index_revision",
    }
)
_STORE_MESSAGES = {
    "archive_version_not_found": "The sealed version or its source was not found.",
    "archive_version_conflict": "The source no longer matches the requested identity.",
    "archive_offline": "The archive copy is unavailable.",
    "archive_integrity": "The sealed bytes do not match the source.",
    "archive_invalid": "The seal request is not valid.",
}


class ArchiveSealError(Exception):
    """Closed seal failure. The message names no path and no digest."""

    def __init__(self, code: str, message: str) -> None:
        if code not in _CODES:
            raise ValueError("archive error code is invalid")
        if not isinstance(message, str) or not message or "/" in message:
            raise ValueError("archive error message is invalid")
        super().__init__(message)
        self.code = code
        self.message = message


class InputRegistryPoolLookup:
    """Resolve a storage pool id to the operator config key, never a path."""

    def __init__(self, registry: object) -> None:
        self._registry = registry

    def config_key_for(self, storage_pool_id: str) -> str:
        try:
            pool = self._registry.get_storage_pool(storage_pool_id)  # type: ignore[attr-defined]
        except KeyError:
            raise ArchiveSealError(
                "archive_offline",
                "The registered storage pool is unavailable.",
            ) from None
        if getattr(pool, "archived_at", None) is not None:
            raise ArchiveSealError(
                "archive_offline",
                "The registered storage pool is unavailable.",
            )
        key = getattr(pool, "config_key", None)
        if not isinstance(key, str) or not key:
            raise ArchiveSealError(
                "archive_offline",
                "The registered storage pool is unavailable.",
            )
        return key


class ReferenceAssemblyLookup:
    """Read the assembly already pinned to a reference-profile revision."""

    def __init__(self, session_factory: object) -> None:
        self._sessions = session_factory

    def registered_assembly(self, reference: Mapping[str, object]) -> str | None:
        if reference.get("kind") != "reference_profile_revision":
            return None
        revision_id = reference.get("reference_profile_revision_id")
        expected = reference.get("public_identity_sha256")
        if not isinstance(revision_id, str) or not isinstance(expected, str):
            raise ArchiveSealError(
                "archive_integrity",
                "The bound reference is not intact.",
            )
        row = self._revision(revision_id)
        if (
            row is None
            or row.public_identity_sha256 != expected
            or not isinstance(row.assembly, str)
            or not row.assembly.strip()
        ):
            raise ArchiveSealError(
                "archive_integrity",
                "The bound reference is not intact.",
            )
        return row.assembly

    def run_assembly(self, run_id: str) -> str | None:
        try:
            with self._sessions() as session:  # type: ignore[operator]
                binding = session.get(RunReferenceBindingRow, run_id)
                if binding is None:
                    return None
                row = session.get(ReferenceProfileRevisionRow, binding.revision_id)
                if (
                    row is None
                    or row.profile_id != binding.profile_id
                    or row.public_identity_sha256
                    != binding.revision_public_identity_sha256
                    or not isinstance(row.assembly, str)
                    or not row.assembly.strip()
                ):
                    raise ArchiveSealError(
                        "archive_integrity",
                        "The bound reference is not intact.",
                    )
                return row.assembly
        except ArchiveSealError:
            raise
        except Exception:
            raise ArchiveSealError(
                "archive_offline",
                "The platform database is unavailable.",
            ) from None

    def _revision(self, revision_id: str) -> SimpleNamespace | None:
        try:
            with self._sessions() as session:  # type: ignore[operator]
                row = session.get(ReferenceProfileRevisionRow, revision_id)
                if row is None:
                    return None
                return SimpleNamespace(
                    public_identity_sha256=row.public_identity_sha256,
                    assembly=row.assembly,
                )
        except ArchiveSealError:
            raise
        except Exception:
            raise ArchiveSealError(
                "archive_offline",
                "The platform database is unavailable.",
            ) from None


class CatalogExperimentGate:
    """Confirm the named experiment revision, not the current head."""

    def __init__(self, catalog: object, session_factory: object) -> None:
        self._catalog = catalog
        self._session_factory = session_factory

    def require_selected_revision(
        self,
        *,
        project_id: str,
        experiment_id: str,
        experiment_revision_id: str,
    ) -> None:
        self._require_row(project_id, experiment_id, experiment_revision_id)
        try:
            detail = self._catalog.get_experiment_revision(  # type: ignore[attr-defined]
                experiment_id,
                experiment_revision_id,
            )
        except CatalogError as exc:
            raise _catalog_error(
                exc,
                missing="The experiment revision was not found.",
            ) from None
        if (
            not isinstance(detail, Mapping)
            or detail.get("experiment_id") != experiment_id
            or detail.get("experiment_revision_id") != experiment_revision_id
        ):
            raise ArchiveSealError(
                "archive_invalid",
                "The experiment revision cannot be sealed.",
            )

    def require_run_publication(
        self,
        *,
        project_id: str,
        experiment_id: str,
        experiment_revision_id: str,
        run_id: str,
        artifact_id: str,
        generation: str,
    ) -> None:
        self.require_selected_revision(
            project_id=project_id,
            experiment_id=experiment_id,
            experiment_revision_id=experiment_revision_id,
        )
        try:
            trace = self._catalog.get_publication_experiment_trace(  # type: ignore[attr-defined]
                project_id=project_id,
                run_id=run_id,
                artifact_id=artifact_id,
                artifact_generation=generation,
            )
        except CatalogError as exc:
            raise _catalog_error(
                exc,
                missing="The run publication was not found.",
            ) from None
        if not _publication_matches(
            trace,
            project_id=project_id,
            experiment_revision_id=experiment_revision_id,
            run_id=run_id,
            artifact_id=artifact_id,
            generation=generation,
        ):
            raise ArchiveSealError(
                "archive_invalid",
                "The run publication is not bound to this experiment revision.",
            )

    def _require_row(
        self,
        project_id: str,
        experiment_id: str,
        experiment_revision_id: str,
    ) -> None:
        try:
            with self._session_factory() as session:  # type: ignore[operator]
                row = session.get(ExperimentRevisionRow, experiment_revision_id)
                project = session.get(ProjectRow, project_id)
        except Exception:
            raise ArchiveSealError(
                "archive_offline",
                "The platform database is unavailable.",
            ) from None
        if (
            row is None
            or row.experiment_id != experiment_id
            or row.project_id != project_id
            or project is None
            or project.archived_at is not None
            or project.project_id != project_id
        ):
            raise ArchiveSealError(
                "archive_version_not_found",
                "The experiment revision was not found.",
            )


class ArchiveSealingService:
    """Copy exact source bytes once and publish the sealed version."""

    def __init__(
        self,
        *,
        archive_root: Path | None,
        objects: object,
        repository: object,
        experiments: object,
        downloads: object | None = None,
        registered_revisions: object | None = None,
        pool_lookup: object | None = None,
        pool_config_path: Path | None = None,
        scratch_root: Path | None = None,
        assemblies: object | None = None,
    ) -> None:
        self._archive_root = archive_root
        self._objects = objects
        self._repository = repository
        self._experiments = experiments
        self._downloads = downloads
        self._registered = registered_revisions
        self._pool_lookup = pool_lookup
        self._pool_config_path = pool_config_path
        self._scratch_root = scratch_root
        self._assemblies = assemblies

    def seal_registered(
        self,
        *,
        project_id: str,
        experiment_id: str,
        experiment_revision_id: str,
        track_key: str,
        display_name: str,
        assembly: str,
        track_type: str,
        result_id: str,
        result_revision_id: str,
        index_result_id: str | None = None,
        index_result_revision_id: str | None = None,
    ) -> dict[str, object]:
        self._require_root()
        project_id = _project_id(project_id)
        experiment_id = _experiment_id(experiment_id)
        experiment_revision_id = _experiment_revision_id(experiment_revision_id)
        track_key, display_name, assembly, track_type = _labels(
            track_key=track_key,
            display_name=display_name,
            assembly=assembly,
            track_type=track_type,
        )
        result_id = _result_id(result_id)
        result_revision_id = _result_revision_id(result_revision_id)
        _complete_pair(
            index_result_id,
            index_result_revision_id,
            "The index source is incomplete.",
        )
        record = self._verified_record(
            project_id=project_id,
            result_id=result_id,
            result_revision_id=result_revision_id,
        )
        index_record = None
        if index_result_id is not None and index_result_revision_id is not None:
            index_result_id = _result_id(index_result_id)
            index_result_revision_id = _result_revision_id(index_result_revision_id)
            if index_result_revision_id == result_revision_id:
                raise ArchiveSealError(
                    "archive_invalid",
                    "The index source must be a different member.",
                )
            index_record = self._verified_record(
                project_id=project_id,
                result_id=index_result_id,
                result_revision_id=index_result_revision_id,
            )
        self._require_registered_formats(record, index_record, track_type)
        self._require_source_assembly(
            assembly,
            self._registered_bound_assembly(record),
            self._registered_bound_assembly(index_record),
        )
        self._experiments.require_selected_revision(  # type: ignore[attr-defined]
            project_id=project_id,
            experiment_id=experiment_id,
            experiment_revision_id=experiment_revision_id,
        )
        source = {
            "kind": "registered",
            "result_id": result_id,
            "result_revision_id": result_revision_id,
        }
        if index_record is not None:
            source["index_result_id"] = index_result_id
            source["index_result_revision_id"] = index_result_revision_id
        members = {"data": self._seal_registered_member(record)}
        if index_record is not None:
            members["index"] = self._seal_registered_member(index_record)
        return self._publish(
            project_id=project_id,
            experiment_id=experiment_id,
            experiment_revision_id=experiment_revision_id,
            track_key=track_key,
            display_name=display_name,
            assembly=assembly,
            track_type=track_type,
            source=source,
            members=members,
        )

    def seal_run(
        self,
        *,
        project_id: str,
        experiment_id: str,
        experiment_revision_id: str,
        track_key: str,
        display_name: str,
        assembly: str,
        track_type: str,
        run_id: str,
        artifact_id: str,
        generation: str,
        revision: str,
        index_artifact_id: str | None = None,
        index_generation: str | None = None,
        index_revision: str | None = None,
    ) -> dict[str, object]:
        self._require_root()
        project_id = _project_id(project_id)
        experiment_id = _experiment_id(experiment_id)
        experiment_revision_id = _experiment_revision_id(experiment_revision_id)
        track_key, display_name, assembly, track_type = _labels(
            track_key=track_key,
            display_name=display_name,
            assembly=assembly,
            track_type=track_type,
        )
        run_id = _run_id(run_id)
        artifact_id = _artifact_id(artifact_id)
        generation = _generation(generation)
        revision = _revision(revision)
        index = _run_index(
            artifact_id,
            generation,
            revision,
            index_artifact_id,
            index_generation,
            index_revision,
        )
        self._experiments.require_run_publication(  # type: ignore[attr-defined]
            project_id=project_id,
            experiment_id=experiment_id,
            experiment_revision_id=experiment_revision_id,
            run_id=run_id,
            artifact_id=artifact_id,
            generation=generation,
        )
        if index is not None:
            self._experiments.require_run_publication(  # type: ignore[attr-defined]
                project_id=project_id,
                experiment_id=experiment_id,
                experiment_revision_id=experiment_revision_id,
                run_id=run_id,
                artifact_id=index[0],
                generation=index[1],
            )
        if self._assemblies is not None:
            self._require_source_assembly(
                assembly,
                self._assemblies.run_assembly(run_id),  # type: ignore[attr-defined]
            )
        source = {
            "kind": "run",
            "run_id": run_id,
            "artifact_id": artifact_id,
            "generation": generation,
            "revision": revision,
        }
        if index is not None:
            source["index_artifact_id"] = index[0]
            source["index_generation"] = index[1]
            source["index_revision"] = index[2]
        data_plan = self._open_run(run_id, artifact_id, generation, revision)
        index_plan = None
        try:
            if index is not None:
                index_plan = self._open_run(run_id, index[0], index[1], index[2])
            _require_declared_format(
                getattr(data_plan, "filename", ""),
                track_type,
                index_present=index_plan is not None,
                role="data",
            )
            if index_plan is not None:
                _require_declared_format(
                    getattr(index_plan, "filename", ""),
                    track_type,
                    index_present=True,
                    role="index",
                )
            members = {"data": self._seal_download(data_plan)}
            if index_plan is not None:
                members["index"] = self._seal_download(index_plan)
            return self._publish(
                project_id=project_id,
                experiment_id=experiment_id,
                experiment_revision_id=experiment_revision_id,
                track_key=track_key,
                display_name=display_name,
                assembly=assembly,
                track_type=track_type,
                source=source,
                members=members,
            )
        finally:
            if index_plan is not None:
                index_plan.close()
            data_plan.close()

    def _require_registered_formats(
        self,
        record: object,
        index_record: object | None,
        track_type: str,
    ) -> None:
        data_format = metadata_from_canonical_json(
            record.canonical_payload  # type: ignore[attr-defined]
        ).format
        _require_declared_format(
            data_format,
            track_type,
            index_present=index_record is not None,
            role="data",
        )
        if index_record is None:
            return
        index_format = metadata_from_canonical_json(
            index_record.canonical_payload  # type: ignore[attr-defined]
        ).format
        _require_declared_format(
            index_format,
            track_type,
            index_present=True,
            role="index",
        )

    def _registered_bound_assembly(self, record: object | None) -> str | None:
        if record is None or self._assemblies is None:
            return None
        metadata = metadata_from_canonical_json(
            record.canonical_payload  # type: ignore[attr-defined]
        )
        return self._assemblies.registered_assembly(  # type: ignore[attr-defined]
            metadata.reference
        )

    def _require_source_assembly(self, operator: str, *bound: str | None) -> None:
        """Keep a verifiable reference assembly and refuse a conflicting label."""

        for value in bound:
            if value is None:
                continue
            if viewer_scaffold(value) != operator:
                raise ArchiveSealError(
                    "archive_invalid",
                    "The assembly does not match the bound reference.",
                )

    def _require_root(self) -> None:
        root = self._archive_root
        if root is None:
            raise ArchiveSealError(
                "archive_invalid",
                "ENCODE_PIPELINE_ARCHIVE_ROOT is not configured.",
            )
        if not isinstance(root, Path) or root.is_symlink() or not root.is_dir():
            raise ArchiveSealError(
                "archive_offline",
                "The archive root is unavailable.",
            )

    def _verified_record(
        self,
        *,
        project_id: str,
        result_id: str,
        result_revision_id: str,
    ) -> object:
        if self._registered is None:
            raise ArchiveSealError(
                "archive_offline",
                "Registered results are unavailable.",
            )
        try:
            record = self._registered.get_revision(result_revision_id)  # type: ignore[attr-defined]
        except Exception:
            raise ArchiveSealError(
                "archive_offline",
                "The platform database is unavailable.",
            ) from None
        if (
            record is None
            or record.project_id != project_id
            or record.result_id != result_id
            or record.result_revision_id != result_revision_id
            or record.project_archived_at is not None
            or record.result_archived_at is not None
        ):
            raise ArchiveSealError(
                "archive_version_not_found",
                "The registered revision was not found.",
            )
        try:
            metadata = metadata_from_canonical_json(record.canonical_payload)
            parents = tuple(
                ParentEdge(
                    result_revision_id=parent.result_revision_id,
                    role=parent.role,
                    revision_digest=parent.revision_digest,
                )
                for parent in record.parents
            )
            if len(parents) != len(metadata.parent_refs):
                raise RegisteredResultPayloadError("parents")
            digest = build_revision_digest(
                project_id=record.project_id,
                result_id=record.result_id,
                result_revision_id=record.result_revision_id,
                revision_number=record.revision_number,
                canonical_payload=record.canonical_payload,
                storage_pool_id=record.storage_pool_id,
                relative_path=record.relative_path,
                size_bytes=record.size_bytes,
                content_sha256=record.content_sha256,
                parents=parents,
            )
        except RegisteredResultPayloadError:
            raise ArchiveSealError(
                "archive_integrity",
                "The registered revision record is not intact.",
            ) from None
        if digest != record.revision_digest:
            raise ArchiveSealError(
                "archive_integrity",
                "The registered revision record is not intact.",
            )
        return record

    def _pool_root(self, storage_pool_id: str) -> Path:
        if self._pool_config_path is None:
            raise ArchiveSealError(
                "archive_invalid",
                "ENCODE_PIPELINE_STORAGE_POOL_CONFIG is not configured.",
            )
        if self._pool_lookup is None:
            raise ArchiveSealError(
                "archive_offline",
                "The registered storage pool is unavailable.",
            )
        try:
            key = self._pool_lookup.config_key_for(storage_pool_id)  # type: ignore[attr-defined]
            config = PrivateStoragePoolConfig.from_file(self._pool_config_path)
            return config.root_for(key)
        except ArchiveSealError:
            raise
        except Exception:
            raise ArchiveSealError(
                "archive_offline",
                "The registered storage pool is unavailable.",
            ) from None

    def _seal_registered_member(self, record: object) -> object:
        root = self._pool_root(record.storage_pool_id)  # type: ignore[attr-defined]
        response_root = self._scratch_root
        if response_root is None:
            response_root = Path(tempfile.gettempdir())
        copy_path: Path | None = None
        directory: Path | None = None
        descriptor: int | None = None
        try:
            copy_path, directory = seal_private_copy(
                root=root,
                relative_path=record.relative_path,  # type: ignore[attr-defined]
                expected_size=record.size_bytes,  # type: ignore[attr-defined]
                expected_sha256=record.content_sha256,  # type: ignore[attr-defined]
                response_root=response_root,
            )
            descriptor = open_verified_copy(copy_path)
            info = os.fstat(descriptor)
            if (
                not stat.S_ISREG(info.st_mode)
                or info.st_size != record.size_bytes  # type: ignore[attr-defined]
            ):
                raise ArchiveSealError(
                    "archive_integrity",
                    "The registered bytes do not match the stored digest.",
                )
            return self._seal_descriptor(
                descriptor,
                record.content_sha256,  # type: ignore[attr-defined]
            )
        except ArchiveSealError:
            raise
        except RegisteredResultBytesError:
            raise ArchiveSealError(
                "archive_integrity",
                "The registered bytes do not match the stored digest.",
            ) from None
        except OSError:
            raise ArchiveSealError(
                "archive_offline",
                "The registered bytes are unavailable.",
            ) from None
        finally:
            if descriptor is not None:
                with suppress(OSError):
                    os.close(descriptor)
            if copy_path is not None and directory is not None:
                discard_private_copy(copy_path, directory)

    def _open_run(
        self,
        run_id: str,
        artifact_id: str,
        generation: str,
        revision: str,
    ) -> object:
        if self._downloads is None:
            raise ArchiveSealError(
                "archive_offline",
                "The run artifact is unavailable.",
            )
        result = self._downloads.prepare(  # type: ignore[attr-defined]
            run_id,
            artifact_id,
            expected_generation=generation,
            expected_revision=revision,
        )
        if result.is_failure or result.value is None:
            raise _download_failure(result)
        return result.value

    def _seal_download(self, plan: object) -> object:
        handle = tempfile.TemporaryFile()
        try:
            written = 0
            try:
                for chunk in plan.iter_bytes():  # type: ignore[attr-defined]
                    if not isinstance(chunk, bytes):
                        raise ArchiveSealError(
                            "archive_integrity",
                            "The run artifact bytes do not match the indexed revision.",
                        )
                    handle.write(chunk)
                    written += len(chunk)
            except ArchiveSealError:
                raise
            except ArtifactDownloadStreamError:
                raise ArchiveSealError(
                    "archive_integrity",
                    "The run artifact bytes do not match the indexed revision.",
                ) from None
            except OSError:
                raise ArchiveSealError(
                    "archive_offline",
                    "The run artifact is unavailable.",
                ) from None
            if written != plan.size_bytes:  # type: ignore[attr-defined]
                raise ArchiveSealError(
                    "archive_integrity",
                    "The run artifact bytes do not match the indexed revision.",
                )
            handle.flush()
            os.lseek(handle.fileno(), 0, os.SEEK_SET)
            member = self._seal_descriptor(handle.fileno(), None)
            if getattr(member, "size_bytes", None) != plan.size_bytes:  # type: ignore[attr-defined]
                raise ArchiveSealError(
                    "archive_integrity",
                    "The sealed member is incomplete.",
                )
            return member
        finally:
            handle.close()

    def _seal_descriptor(
        self,
        descriptor: int,
        expected_sha256: str | None,
    ) -> object:
        try:
            member = self._objects.seal_fd(  # type: ignore[attr-defined]
                descriptor,
                expected_sha256=expected_sha256,
            )
        except ArchiveSealError:
            raise
        except Exception as exc:
            raise _store_failure(exc) from None
        _require_member(member, expected_sha256)
        return member

    def _publish(self, **kwargs: object) -> dict[str, object]:
        try:
            version = self._repository.publish(**kwargs)  # type: ignore[attr-defined]
        except ArchiveSealError:
            raise
        except Exception as exc:
            raise _store_failure(exc) from None
        _require_version(version, **kwargs)
        return require_operator_receipt(_receipt(version))


def require_operator_receipt(receipt: object) -> dict[str, object]:
    """Return the operator JSON object, rejecting paths and extra fields."""

    if not isinstance(receipt, Mapping) or set(receipt) != _RECEIPT_FIELDS:
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    version_id = receipt["archive_version_id"]
    if not _public_token(version_id):
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    for field in (
        "project_id",
        "experiment_id",
        "experiment_revision_id",
        "track_key",
        "display_name",
        "assembly",
        "track_type",
    ):
        if not isinstance(receipt[field], str) or not receipt[field]:
            raise ArchiveSealError(
                "archive_integrity",
                "The archive record does not match the seal.",
            )
    number = receipt["version_number"]
    if isinstance(number, bool) or not isinstance(number, int) or number < 1:
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    source = receipt["source"]
    members = receipt["members"]
    has_index = receipt["has_index"]
    if not isinstance(source, Mapping) or not isinstance(members, Mapping):
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    kind = source.get("kind")
    if kind == "registered":
        allowed = _REGISTERED_SOURCE
        required = {"kind", "result_id", "result_revision_id"}
    elif kind == "run":
        allowed = _RUN_SOURCE
        required = {"kind", "run_id", "artifact_id", "generation", "revision"}
    else:
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    if not required <= set(source) or not set(source) <= allowed:
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    if any(not isinstance(key, str) or not isinstance(value, str) for key, value in source.items()):
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    index_keys = (
        ("index_result_id", "index_result_revision_id")
        if kind == "registered"
        else ("index_artifact_id", "index_generation", "index_revision")
    )
    indexed = any(key in source for key in index_keys)
    if indexed and any(key not in source for key in index_keys):
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    if has_index is not indexed or set(members) != (
        {"data", "index"} if indexed else {"data"}
    ):
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    rendered_members: dict[str, dict[str, object]] = {}
    for role in ("data", "index"):
        if role not in members:
            continue
        body = members[role]
        if not isinstance(body, Mapping) or set(body) != {"sha256", "size_bytes"}:
            raise ArchiveSealError(
                "archive_integrity",
                "The archive record does not match the seal.",
            )
        size_bytes = body["size_bytes"]
        digest = body["sha256"]
        if (
            isinstance(size_bytes, bool)
            or not isinstance(size_bytes, int)
            or size_bytes < 0
            or not isinstance(digest, str)
            or _SHA256.fullmatch(digest) is None
        ):
            raise ArchiveSealError(
                "archive_integrity",
                "The archive record does not match the seal.",
            )
        rendered_members[role] = {"sha256": digest, "size_bytes": size_bytes}
    return {
        "archive_version_id": version_id,
        "project_id": receipt["project_id"],
        "experiment_id": receipt["experiment_id"],
        "experiment_revision_id": receipt["experiment_revision_id"],
        "track_key": receipt["track_key"],
        "version_number": number,
        "display_name": receipt["display_name"],
        "assembly": receipt["assembly"],
        "track_type": receipt["track_type"],
        "source": {str(key): str(source[key]) for key in sorted(source)},
        "has_index": indexed,
        "members": rendered_members,
    }


def _receipt(version: object) -> dict[str, object]:
    source = getattr(version, "source", None)
    members = getattr(version, "members", None)
    if not isinstance(source, Mapping) or not isinstance(members, Mapping):
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    rendered_members = {}
    for role, member in members.items():
        rendered_members[str(role)] = {
            "sha256": getattr(member, "sha256", None),
            "size_bytes": getattr(member, "size_bytes", None),
        }
    return {
        "archive_version_id": getattr(version, "archive_version_id", None),
        "project_id": getattr(version, "project_id", None),
        "experiment_id": getattr(version, "experiment_id", None),
        "experiment_revision_id": getattr(version, "experiment_revision_id", None),
        "track_key": getattr(version, "track_key", None),
        "version_number": getattr(version, "version_number", None),
        "display_name": getattr(version, "display_name", None),
        "assembly": getattr(version, "assembly", None),
        "track_type": getattr(version, "track_type", None),
        "source": dict(source),
        "has_index": getattr(version, "has_index", None),
        "members": rendered_members,
    }


def _require_version(version: object, **kwargs: object) -> None:
    source = kwargs["source"]
    members = kwargs["members"]
    observed_source = getattr(version, "source", None)
    observed_members = getattr(version, "members", None)
    if (
        getattr(version, "project_id", None) != kwargs["project_id"]
        or getattr(version, "experiment_id", None) != kwargs["experiment_id"]
        or getattr(version, "experiment_revision_id", None)
        != kwargs["experiment_revision_id"]
        or getattr(version, "track_key", None) != kwargs["track_key"]
        or getattr(version, "display_name", None) != kwargs["display_name"]
        or getattr(version, "assembly", None) != kwargs["assembly"]
        or getattr(version, "track_type", None) != kwargs["track_type"]
        or not isinstance(observed_source, Mapping)
        or dict(observed_source) != source
        or not isinstance(observed_members, Mapping)
        or set(observed_members) != set(members)  # type: ignore[arg-type]
    ):
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )
    assert isinstance(members, Mapping)
    for role, member in members.items():
        stored = observed_members[role]
        if (
            getattr(stored, "size_bytes", None) != getattr(member, "size_bytes", None)
            or getattr(stored, "sha256", None) != getattr(member, "sha256", None)
        ):
            raise ArchiveSealError(
                "archive_integrity",
                "The archive record does not match the seal.",
            )
    if getattr(version, "has_index", None) is not ("index" in members):
        raise ArchiveSealError(
            "archive_integrity",
            "The archive record does not match the seal.",
        )


def _require_member(member: object, expected_sha256: str | None) -> None:
    size_bytes = getattr(member, "size_bytes", None)
    digest = getattr(member, "sha256", None)
    if (
        isinstance(size_bytes, bool)
        or not isinstance(size_bytes, int)
        or size_bytes < 0
        or not isinstance(digest, str)
        or _SHA256.fullmatch(digest) is None
        or (expected_sha256 is not None and digest != expected_sha256)
    ):
        raise ArchiveSealError(
            "archive_integrity",
            "The sealed member is incomplete.",
        )


def _store_failure(exc: BaseException) -> ArchiveSealError:
    code = getattr(exc, "code", None)
    if isinstance(code, str) and code in _STORE_MESSAGES:
        return ArchiveSealError(code, _STORE_MESSAGES[code])
    return ArchiveSealError(
        "archive_offline",
        "The archive store is unavailable.",
    )


def _catalog_error(exc: CatalogError, *, missing: str) -> ArchiveSealError:
    if exc.status_code == 404:
        return ArchiveSealError("archive_version_not_found", missing)
    return ArchiveSealError(
        "archive_invalid",
        "The experiment revision cannot be sealed.",
    )


def _publication_matches(
    trace: object,
    *,
    project_id: str,
    experiment_revision_id: str,
    run_id: str,
    artifact_id: str,
    generation: str,
) -> bool:
    return (
        isinstance(trace, Mapping)
        and trace.get("trace_status") == "experiment_snapshot"
        and trace.get("project_id") == project_id
        and trace.get("experiment_revision_id") == experiment_revision_id
        and trace.get("run_id") == run_id
        and trace.get("artifact_id") == artifact_id
        and trace.get("artifact_generation") == generation
    )


def _download_failure(result: object) -> ArchiveSealError:
    errors = getattr(result, "errors", ())
    issue = errors[0] if errors else None
    if issue is None:
        return ArchiveSealError(
            "archive_offline",
            "The run artifact is unavailable.",
        )
    reason = ""
    context = getattr(issue, "context", None)
    if isinstance(context, Mapping):
        raw = context.get("reason_code")
        if isinstance(raw, str):
            reason = raw
    code = getattr(issue, "code", "")
    if code in {"RUN_NOT_FOUND", "RUN_ARTIFACT_NOT_FOUND"}:
        return ArchiveSealError(
            "archive_version_not_found",
            "The run artifact was not found.",
        )
    if reason in {
        "ARTIFACT_DOWNLOAD_GENERATION_CHANGED",
        "ARTIFACT_DOWNLOAD_REVISION_CHANGED",
    }:
        return ArchiveSealError(
            "archive_version_conflict",
            "The artifact generation or revision does not match.",
        )
    if code == "RUN_ARTIFACT_DOWNLOAD_UNAVAILABLE":
        return ArchiveSealError(
            "archive_offline",
            "The run artifact is unavailable.",
        )
    if (
        code == "RUN_ARTIFACT_DOWNLOAD_DATA_INVALID"
        or reason == "ARTIFACT_DOWNLOAD_EXPECTATION_INVALID"
    ):
        return ArchiveSealError(
            "archive_invalid",
            "The run artifact record is not sealable.",
        )
    return ArchiveSealError(
        "archive_integrity",
        "The run artifact bytes do not match the indexed revision.",
    )


def _explicit_track(value: str) -> str | None:
    """Return a track type only when the source name states one."""

    lower = value.lower()
    if lower.endswith((".tbi", ".idx", ".bai")):
        return None
    if lower.endswith((".bw", ".bigwig")) or lower in {"bigwig", "bw"}:
        return "bigwig"
    if lower.endswith(".bed") or lower == "bed":
        return "bed"
    token = "".join(character for character in lower if character.isalnum())
    if "longrange" in token or lower.endswith((".pairs", ".hic", ".bedpe")):
        return "longrange"
    return None


def _format_requires_index(value: str) -> bool:
    token = "".join(character for character in value.lower() if character.isalnum())
    return "longrange" in token and "index" in token


def _require_declared_format(
    declared: object,
    track_type: str,
    *,
    index_present: bool,
    role: str,
) -> None:
    if not isinstance(declared, str) or not declared:
        raise ArchiveSealError(
            "archive_invalid",
            "The source format does not match the track type.",
        )
    explicit = _explicit_track(declared)
    if role == "index":
        if explicit is not None:
            raise ArchiveSealError(
                "archive_invalid",
                "The index source does not match the track.",
            )
        return
    if explicit is not None and explicit != track_type:
        raise ArchiveSealError(
            "archive_invalid",
            "The source format does not match the track type.",
        )
    if _format_requires_index(declared) and not index_present:
        raise ArchiveSealError(
            "archive_invalid",
            "The index source is incomplete.",
        )


def _labels(
    *,
    track_key: object,
    display_name: object,
    assembly: object,
    track_type: object,
) -> tuple[str, str, str, str]:
    if not isinstance(track_key, str) or _TRACK_KEY.fullmatch(track_key) is None:
        raise ArchiveSealError("archive_invalid", "The seal request is not valid.")
    if not isinstance(display_name, str):
        raise ArchiveSealError("archive_invalid", "The seal request is not valid.")
    normalized = display_name.strip()
    if (
        not normalized
        or len(normalized) > 255
        or any(ord(character) < 32 or ord(character) == 127 for character in normalized)
        or "/" in normalized
        or "\\" in normalized
    ):
        raise ArchiveSealError("archive_invalid", "The seal request is not valid.")
    if not isinstance(assembly, str):
        raise ArchiveSealError("archive_invalid", "The seal request is not valid.")
    canonical = viewer_scaffold(assembly)
    if canonical is None:
        raise ArchiveSealError(
            "archive_invalid",
            "The assembly is not a supported WashU scaffold.",
        )
    if track_type not in _TRACK_TYPES:
        raise ArchiveSealError("archive_invalid", "The seal request is not valid.")
    return track_key, normalized, canonical, str(track_type)


def _complete_pair(left: object, right: object, message: str) -> None:
    if (left is None) != (right is None):
        raise ArchiveSealError("archive_invalid", message)


def _run_index(
    artifact_id: str,
    generation: str,
    revision: str,
    index_artifact_id: str | None,
    index_generation: str | None,
    index_revision: str | None,
) -> tuple[str, str, str] | None:
    present = (
        index_artifact_id is not None,
        index_generation is not None,
        index_revision is not None,
    )
    if not any(present):
        return None
    if not all(present):
        raise ArchiveSealError("archive_invalid", "The index source is incomplete.")
    assert index_artifact_id is not None
    assert index_generation is not None
    assert index_revision is not None
    index_artifact_id = _artifact_id(index_artifact_id)
    index_generation = _generation(index_generation)
    index_revision = _revision(index_revision)
    if (
        index_artifact_id == artifact_id
        and index_generation == generation
        and index_revision == revision
    ):
        raise ArchiveSealError(
            "archive_invalid",
            "The index source must be a different member.",
        )
    return index_artifact_id, index_generation, index_revision


def _invalid() -> ArchiveSealError:
    return ArchiveSealError("archive_invalid", "The seal request is not valid.")


def _project_id(value: object) -> str:
    try:
        return validate_project_id(value)
    except RegisteredResultPayloadError:
        raise _invalid() from None


def _result_id(value: object) -> str:
    try:
        return validate_result_id(value)
    except RegisteredResultPayloadError:
        raise _invalid() from None


def _result_revision_id(value: object) -> str:
    try:
        return validate_result_revision_id(value)
    except RegisteredResultPayloadError:
        raise _invalid() from None


def _experiment_id(value: object) -> str:
    if not isinstance(value, str) or _EXPERIMENT_ID.fullmatch(value) is None:
        raise _invalid()
    return value


def _experiment_revision_id(value: object) -> str:
    if (
        not isinstance(value, str)
        or _EXPERIMENT_REVISION_ID.fullmatch(value) is None
    ):
        raise _invalid()
    return value


def _run_id(value: object) -> str:
    if not isinstance(value, str) or _SAFE_RUN_ID.fullmatch(value) is None:
        raise _invalid()
    return value


def _artifact_id(value: object) -> str:
    if not isinstance(value, str) or _SAFE_ARTIFACT_ID.fullmatch(value) is None:
        raise _invalid()
    return value


def _generation(value: object) -> str:
    try:
        return validate_artifact_generation(value)
    except (TypeError, ValueError):
        raise _invalid() from None


def _revision(value: object) -> str:
    try:
        return validate_artifact_revision(value)
    except (TypeError, ValueError):
        raise _invalid() from None


def _public_token(value: object) -> bool:
    return (
        isinstance(value, str)
        and bool(value)
        and len(value) <= 200
        and all(
            not character.isspace() and 31 < ord(character) != 127
            for character in value
        )
    )
