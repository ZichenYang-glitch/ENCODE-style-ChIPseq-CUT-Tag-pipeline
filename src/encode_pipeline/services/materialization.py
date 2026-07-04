"""Filesystem materialization boundary for workspace plans."""

from __future__ import annotations

from pathlib import Path

from encode_pipeline.platform.adapters import WorkspacePlan
from encode_pipeline.platform.planning import WorkspacePathError, WorkspacePathPolicy
from encode_pipeline.platform.results import Issue, Result

_WORKSPACE_PATH_ERROR_MESSAGES = {
    "WORKSPACE_BASE_DIR_RELATIVE": "Workspace path policy rejected a planned path.",
    "WORKSPACE_PATH_ABSOLUTE": "Planned path must not be absolute.",
    "WORKSPACE_PATH_TRAVERSAL": "Planned path must not contain '.' or '..' components.",
    "WORKSPACE_PATH_INVALID": "Planned path is invalid.",
    "WORKSPACE_PATH_ESCAPE": "Planned path escapes base_dir.",
}


class WorkspaceMaterializer:
    """Pure filesystem materialization boundary for workspace plans."""

    def materialize(
        self,
        plan: WorkspacePlan,
        base_dir: Path,
    ) -> Result[None]:
        """Create directories and files described by ``plan`` under ``base_dir``."""
        if not isinstance(plan, WorkspacePlan):
            return Result.failure(
                [
                    Issue(
                        code="WORKSPACE_MATERIALIZATION_INVALID_PLAN",
                        message="plan must be a WorkspacePlan.",
                        severity="error",
                        path="plan",
                        source="workspace_materializer",
                    )
                ]
            )

        if not isinstance(base_dir, Path) or not base_dir.is_absolute():
            return Result.failure(
                [
                    Issue(
                        code="WORKSPACE_MATERIALIZATION_BASE_DIR_RELATIVE",
                        message="base_dir must be an absolute Path.",
                        severity="error",
                        path="base_dir",
                        source="workspace_materializer",
                    )
                ]
            )

        # Reject any symlink in the logical path of base_dir, including
        # symlinked existing ancestors and a broken-symlink base_dir itself.
        for component in (*base_dir.parents[::-1], base_dir):
            try:
                if component.is_symlink():
                    return Result.failure(
                        [
                            Issue(
                                code="WORKSPACE_MATERIALIZATION_SYMLINK",
                                message="base_dir path contains a symlink.",
                                severity="error",
                                path="base_dir",
                                source="workspace_materializer",
                            )
                        ]
                    )
            except OSError:
                continue

        try:
            base_dir.mkdir(parents=True, exist_ok=True)
        except OSError:
            return Result.failure(
                [
                    Issue(
                        code="WORKSPACE_MATERIALIZATION_BASE_DIR_CREATE_ERROR",
                        message="base_dir could not be created.",
                        severity="error",
                        path="base_dir",
                        source="workspace_materializer",
                    )
                ]
            )

        policy = WorkspacePathPolicy(base_dir=base_dir)

        resolved_files: list[tuple[Path, bytes, str]] = []
        for index, (file_path, contents) in enumerate(plan.files):
            path_locator = f"workspace_plan.files[{index}]"
            try:
                resolved = policy.resolve(file_path)
            except WorkspacePathError as exc:
                return Result.failure(
                    [
                        Issue(
                            code=f"WORKSPACE_MATERIALIZATION_PATH_POLICY_{exc.code.removeprefix('WORKSPACE_')}",
                            message=_WORKSPACE_PATH_ERROR_MESSAGES.get(
                                exc.code, "Planned path violates workspace path policy."
                            ),
                            severity="error",
                            path=path_locator,
                            source="workspace_materializer",
                        )
                    ]
                )

            for parent in resolved.relative_to(base_dir).parents:
                parent_path = base_dir / parent
                if parent_path.is_symlink():
                    return Result.failure(
                        [
                            Issue(
                                code="WORKSPACE_MATERIALIZATION_SYMLINK",
                                message="Parent directory under base_dir is a symlink.",
                                severity="error",
                                path=path_locator,
                                source="workspace_materializer",
                            )
                        ]
                    )

            if resolved.is_symlink():
                return Result.failure(
                    [
                        Issue(
                            code="WORKSPACE_MATERIALIZATION_SYMLINK",
                            message="Planned file path is a symlink.",
                            severity="error",
                            path=path_locator,
                            source="workspace_materializer",
                        )
                    ]
                )

            if resolved.exists():
                if not resolved.is_file():
                    return Result.failure(
                        [
                            Issue(
                                code="WORKSPACE_MATERIALIZATION_WRONG_TYPE",
                                message="Planned file path exists but is not a regular file.",
                                severity="error",
                                path=path_locator,
                                source="workspace_materializer",
                            )
                        ]
                    )
                return Result.failure(
                    [
                        Issue(
                            code="WORKSPACE_MATERIALIZATION_ALREADY_EXISTS",
                            message="Planned file already exists.",
                            severity="error",
                            path=path_locator,
                            source="workspace_materializer",
                        )
                    ]
                )

            resolved_files.append((resolved, contents, path_locator))

        resolved_directories: list[Path] = []
        for index, directory in enumerate(plan.directories):
            path_locator = f"workspace_plan.directories[{index}]"
            try:
                resolved = policy.resolve(directory)
            except WorkspacePathError as exc:
                return Result.failure(
                    [
                        Issue(
                            code=f"WORKSPACE_MATERIALIZATION_PATH_POLICY_{exc.code.removeprefix('WORKSPACE_')}",
                            message=_WORKSPACE_PATH_ERROR_MESSAGES.get(
                                exc.code, "Planned path violates workspace path policy."
                            ),
                            severity="error",
                            path=path_locator,
                            source="workspace_materializer",
                        )
                    ]
                )

            if resolved.is_symlink():
                return Result.failure(
                    [
                        Issue(
                            code="WORKSPACE_MATERIALIZATION_SYMLINK",
                            message="Planned directory path is a symlink.",
                            severity="error",
                            path=path_locator,
                            source="workspace_materializer",
                        )
                    ]
                )

            for parent in resolved.parents:
                try:
                    parent.relative_to(base_dir)
                except ValueError:
                    continue
                if parent == base_dir:
                    continue
                try:
                    if parent.is_symlink():
                        return Result.failure(
                            [
                                Issue(
                                    code="WORKSPACE_MATERIALIZATION_SYMLINK",
                                    message="Planned directory path has a symlinked parent under base_dir.",
                                    severity="error",
                                    path=path_locator,
                                    source="workspace_materializer",
                                )
                            ]
                        )
                except OSError:
                    continue

            if resolved.exists() and not resolved.is_dir():
                return Result.failure(
                    [
                        Issue(
                            code="WORKSPACE_MATERIALIZATION_WRONG_TYPE",
                            message="Planned directory path exists but is not a directory.",
                            severity="error",
                            path=path_locator,
                            source="workspace_materializer",
                        )
                    ]
                )

            resolved_directories.append(resolved)

        for resolved in resolved_directories:
            try:
                resolved.mkdir(parents=True, exist_ok=True)
            except OSError:
                return Result.failure(
                    [
                        Issue(
                            code="WORKSPACE_MATERIALIZATION_WRITE_ERROR",
                            message="Failed to create directory.",
                            severity="error",
                            path="workspace_plan.directories",
                            source="workspace_materializer",
                        )
                    ]
                )

        for resolved, contents, path_locator in resolved_files:
            try:
                resolved.parent.mkdir(parents=True, exist_ok=True)
                with open(resolved, "xb") as handle:
                    handle.write(contents)
            except FileExistsError:
                return Result.failure(
                    [
                        Issue(
                            code="WORKSPACE_MATERIALIZATION_ALREADY_EXISTS",
                            message="Planned file already exists.",
                            severity="error",
                            path=path_locator,
                            source="workspace_materializer",
                        )
                    ]
                )
            except OSError:
                return Result.failure(
                    [
                        Issue(
                            code="WORKSPACE_MATERIALIZATION_WRITE_ERROR",
                            message="Failed to write file.",
                            severity="error",
                            path=path_locator,
                            source="workspace_materializer",
                        )
                    ]
                )

        return Result.success(None)
