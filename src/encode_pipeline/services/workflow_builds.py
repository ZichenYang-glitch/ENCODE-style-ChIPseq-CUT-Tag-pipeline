"""Content-addressed identities for controlled workflow source bundles."""

from __future__ import annotations

from datetime import datetime, timezone
from hashlib import sha256
from pathlib import Path

from encode_pipeline.platform.builds import WorkflowBuildIdentity
from encode_pipeline.platform.registry import WorkflowRegistry
from encode_pipeline.platform.results import Issue, Result


BUILD_IDENTITY_SCHEME = "sha256-tree-v1"
ENCODE_LOGICAL_ENTRYPOINT = "workflow/Snakefile"
_CACHE_DIRECTORY_NAMES = frozenset({"__pycache__", ".mypy_cache", ".pytest_cache"})


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

    def _source_manifest(self) -> tuple[tuple[str, bytes], ...]:
        root = self._project_root
        if not root.is_dir() or root.is_symlink():
            raise ValueError("project root is unavailable")

        entries: dict[str, bytes] = {}
        self._add_required_file(entries, root / "pyproject.toml")
        self._add_tree(
            entries,
            root / "src" / "encode_pipeline",
            suffixes=frozenset({".py"}),
        )
        self._add_tree(entries, root / "workflow", suffixes=None)
        self._add_tree(
            entries,
            root / "profiles" / "default",
            suffixes=None,
        )
        self._add_tree(
            entries,
            root / "scripts",
            suffixes=None,
            recursive=False,
        )
        entrypoint = ENCODE_LOGICAL_ENTRYPOINT
        if entrypoint not in entries:
            raise ValueError("workflow entrypoint is missing")
        return tuple(sorted(entries.items()))

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
