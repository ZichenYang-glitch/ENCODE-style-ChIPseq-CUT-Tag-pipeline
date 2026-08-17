"""Content-addressed identities for controlled workflow source bundles."""

from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path
import re

from encode_pipeline.platform.adapters import (
    WorkflowAdapter,
    WorkflowBuildIdentityProvidingAdapter,
)
from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result


BUILD_IDENTITY_SCHEME = "sha256-tree-v1"
ENCODE_LOGICAL_ENTRYPOINT = "workflow/Snakefile"
_CACHE_DIRECTORY_NAMES = frozenset({"__pycache__", ".mypy_cache", ".pytest_cache"})
_SCRIPT_REFERENCE = re.compile(rb"(?<![0-9A-Za-z_.-])scripts/([^\x00\s'\"`\\|;&()<>]+)")
_SCRIPT_NAME = re.compile(r"^[0-9A-Za-z][0-9A-Za-z_.-]{0,254}$")


class WorkflowBuildIdentityProvider:
    """Fingerprint the source bundle selected by the registered adapter."""

    def __init__(
        self,
        registry: WorkflowRegistry,
        *,
        project_root: Path | None = None,
    ) -> None:
        if not isinstance(registry, WorkflowRegistry):
            raise ValueError("registry must be a WorkflowRegistry")
        root = (
            Path(__file__).resolve().parents[3]
            if project_root is None
            else project_root
        )
        if not isinstance(root, Path) or not root.is_absolute():
            raise ValueError("project_root must be an absolute pathlib.Path")
        self._registry = registry
        self._project_root = root

    @property
    def project_root(self) -> Path:
        """Return the local source root used to read controlled files."""
        return self._project_root

    @property
    def registry(self) -> WorkflowRegistry:
        """Return the exact registry whose adapters provide build identities."""
        return self._registry

    def source_manifest(self) -> tuple[tuple[str, bytes], ...]:
        """Return the controlled ENCODE runtime bytes in canonical path order."""
        return self._source_manifest()

    def capture(self, workflow_id: str) -> Result[WorkflowBuildIdentity]:
        """Return the current build identity without leaking local paths."""
        try:
            adapter = self._registry.get(workflow_id)
        except (KeyError, ValueError):
            return Result.failure(
                [
                    Issue(
                        code="WORKFLOW_BUILD_WORKFLOW_NOT_FOUND",
                        message="Workflow build identity could not be resolved.",
                        severity="error",
                        path="workflow_id",
                        source="workflow_build_identity",
                    )
                ]
            )

        return self._capture_resolved_adapter(
            adapter,
            allow_encode_fallback=self._registry.uses_encode_execution_fallback(
                adapter
            ),
        )

    def _capture_resolved_adapter(
        self,
        adapter: WorkflowAdapter,
        *,
        allow_encode_fallback: bool,
    ) -> Result[WorkflowBuildIdentity]:
        if isinstance(adapter, WorkflowBuildIdentityProvidingAdapter):
            return self._capture_adapter_identity(adapter)
        if not allow_encode_fallback:
            return _source_unavailable()

        try:
            manifest = self._source_manifest()
            digest = self._digest_manifest(
                workflow_id=adapter.metadata.workflow_id,
                adapter_version=adapter.metadata.version,
                manifest=manifest,
            )
            identity = WorkflowBuildIdentity(
                workflow_id=adapter.metadata.workflow_id,
                adapter_version=adapter.metadata.version,
                scheme=BUILD_IDENTITY_SCHEME,
                logical_entrypoint=ENCODE_LOGICAL_ENTRYPOINT,
                digest=digest,
                captured_at=datetime.now(timezone.utc),
            )
        except (OSError, ValueError):
            return Result.failure(
                [
                    Issue(
                        code="WORKFLOW_BUILD_SOURCE_UNAVAILABLE",
                        message="Controlled workflow source could not be fingerprinted.",
                        severity="error",
                        path="workflow",
                        source="workflow_build_identity",
                    )
                ]
            )
        return Result.success(identity)

    def capture_executable(
        self,
        workflow_id: str,
    ) -> Result[WorkflowBuildIdentity]:
        """Capture identity only while the registered execution is admitted."""
        try:
            adapter = self._registry.get(workflow_id)
        except (KeyError, ValueError):
            return _source_unavailable()
        from encode_pipeline.services.workflow_info import (
            resolve_workflow_availability,
        )

        if (
            resolve_workflow_availability(adapter, registry=self._registry).execution
            != "available"
        ):
            return _source_unavailable()
        return self.capture(workflow_id)

    def capture_resolved_executable(
        self,
        adapter: WorkflowAdapter,
    ) -> Result[WorkflowBuildIdentity]:
        """Capture one verified run-scoped adapter against registry coordinates."""
        if not isinstance(adapter, WorkflowAdapter):
            return _source_unavailable()
        try:
            registered = self._registry.get(adapter.metadata.workflow_id)
            compatible = (
                adapter.metadata.workflow_id == registered.metadata.workflow_id
                and adapter.metadata.version == registered.metadata.version
                and adapter.metadata.engines == registered.metadata.engines
                and set(registered.capabilities.supports).issubset(
                    adapter.capabilities.supports
                )
            )
        except (AttributeError, KeyError, TypeError, ValueError):
            return _source_unavailable()
        if not compatible:
            return _source_unavailable()

        allow_encode_fallback = self._registry.uses_encode_execution_fallback(
            registered
        )
        if not allow_encode_fallback:
            from encode_pipeline.services.workflow_info import (
                resolve_workflow_availability,
            )

            if resolve_workflow_availability(adapter).execution != "available":
                return _source_unavailable()
        return self._capture_resolved_adapter(
            adapter,
            allow_encode_fallback=allow_encode_fallback,
        )

    @staticmethod
    def _capture_adapter_identity(
        adapter: WorkflowBuildIdentityProvidingAdapter,
    ) -> Result[WorkflowBuildIdentity]:
        """Accept only a matching identity while hiding adapter internals."""
        try:
            adapter_result = adapter.capture_build_identity()
            if not isinstance(adapter_result, Result) or adapter_result.is_failure:
                return _source_unavailable()
            identity = adapter_result.value
            if not isinstance(identity, WorkflowBuildIdentity):
                return _source_unavailable()
            metadata = adapter.metadata
            if (
                identity.workflow_id != metadata.workflow_id
                or identity.adapter_version != metadata.version
            ):
                return _source_unavailable()
        except Exception:
            return _source_unavailable()
        return Result.success(identity)

    def _source_manifest(self) -> tuple[tuple[str, bytes], ...]:
        """Capture only the independently switchable ENCODE runtime closure."""
        root = self._project_root
        if not root.is_dir() or root.is_symlink():
            raise ValueError("project root is unavailable")

        entries: dict[str, bytes] = {}
        self._add_required_file(
            entries,
            root / "docs" / "architecture" / "artifact-inventory.yaml",
        )
        self._add_tree(entries, root / "workflow", suffixes=None)
        self._add_tree(
            entries,
            root / "profiles" / "default",
            suffixes=None,
        )
        for script_name in self._referenced_runtime_scripts(entries):
            self._add_required_file(entries, root / "scripts" / script_name)
        entrypoint = ENCODE_LOGICAL_ENTRYPOINT
        if entrypoint not in entries:
            raise ValueError("workflow entrypoint is missing")
        return tuple(sorted(entries.items()))

    @staticmethod
    def _referenced_runtime_scripts(entries: dict[str, bytes]) -> tuple[str, ...]:
        """Return the exact top-level scripts named by workflow/profile bytes."""
        names: set[str] = set()
        for logical_path, content in entries.items():
            if not (
                logical_path.startswith("workflow/")
                or logical_path.startswith("profiles/default/")
            ):
                continue
            for match in _SCRIPT_REFERENCE.finditer(content):
                try:
                    name = match.group(1).decode("ascii")
                except UnicodeDecodeError as exc:
                    raise ValueError("runtime script reference is invalid") from exc
                if _SCRIPT_NAME.fullmatch(name) is None:
                    raise ValueError("runtime script reference is invalid")
                names.add(name)
        return tuple(sorted(names))

    def _add_tree(
        self,
        entries: dict[str, bytes],
        directory: Path,
        *,
        suffixes: frozenset[str] | None,
        recursive: bool = True,
    ) -> None:
        if not directory.is_dir() or directory.is_symlink():
            raise ValueError("required source directory is unavailable")
        added = 0
        candidates = directory.rglob("*") if recursive else directory.iterdir()
        for path in sorted(candidates):
            relative_to_tree = path.relative_to(directory)
            if any(part in _CACHE_DIRECTORY_NAMES for part in relative_to_tree.parts):
                continue
            if path.suffix == ".pyc":
                continue
            if path.is_symlink():
                raise ValueError("controlled source must not contain symlinks")
            if path.is_dir():
                continue
            if not path.is_file():
                raise ValueError("controlled source must contain regular files")
            if (
                suffixes is not None
                and path.name != "Snakefile"
                and path.suffix not in suffixes
            ):
                continue
            self._add_required_file(entries, path)
            added += 1
        if added == 0:
            raise ValueError("required source directory has no controlled files")

    def _add_required_file(self, entries: dict[str, bytes], path: Path) -> None:
        if path.is_symlink() or not path.is_file():
            raise ValueError("required source file is unavailable")
        try:
            logical_path = path.relative_to(self._project_root).as_posix()
        except ValueError as exc:
            raise ValueError("controlled source escaped the project root") from exc
        if logical_path in entries:
            raise ValueError("controlled source path was duplicated")
        entries[logical_path] = path.read_bytes()

    @staticmethod
    def _digest_manifest(
        *,
        workflow_id: str,
        adapter_version: str,
        manifest: tuple[tuple[str, bytes], ...],
    ) -> str:
        digest = sha256()
        for value in (
            b"encode-pipeline-workflow-build",
            BUILD_IDENTITY_SCHEME.encode("utf-8"),
            workflow_id.encode("utf-8"),
            adapter_version.encode("utf-8"),
            b"snakemake",
            ENCODE_LOGICAL_ENTRYPOINT.encode("utf-8"),
        ):
            _update_framed(digest, value)
        for logical_path, contents in manifest:
            _update_framed(digest, logical_path.encode("utf-8"))
            _update_framed(digest, sha256(contents).digest())
        return digest.hexdigest()


def _update_framed(digest, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, byteorder="big", signed=False))
    digest.update(value)


def _source_unavailable() -> Result[WorkflowBuildIdentity]:
    return Result.failure(
        [
            Issue(
                code="WORKFLOW_BUILD_SOURCE_UNAVAILABLE",
                message="Controlled workflow source could not be fingerprinted.",
                severity="error",
                path="workflow",
                source="workflow_build_identity",
            )
        ]
    )
