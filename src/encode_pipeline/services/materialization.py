"""Filesystem materialization boundary for workspace plans."""

from __future__ import annotations

from pathlib import Path

from encode_pipeline.platform.adapters import WorkspacePlan
from encode_pipeline.platform.results import Issue, Result


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

        return Result.success(None)
