"""ENCODE-owned Snakemake commands and private workspace execution contract."""

from __future__ import annotations

from dataclasses import dataclass
import hashlib
import json
import os
from pathlib import Path
from typing import Any

from encode_pipeline.adapters.encode_authoring import MAX_CORES
from encode_pipeline.platform.adapters import CommandSpec, WorkspacePlan
from encode_pipeline.platform.planning import WorkspacePathPolicy
from encode_pipeline.platform.results import Issue, Result


EXECUTION_CONFIG_PATH = "config/encode-execution.json"
EXECUTION_CONFIG_SCHEMA_VERSION = "1.0.0"


@dataclass(frozen=True, init=False)
class EncodeExecutionBinding:
    """Operator-owned source/runtime coordinates; never supplied by run inputs."""

    project_root: Path
    snakemake_executable: Path | None
    conda_prefix: Path | None

    def __init__(
        self,
        *,
        project_root: Path | None = None,
        snakemake_executable: Path | None = None,
        conda_prefix: Path | None = None,
    ) -> None:
        root = (
            Path(__file__).resolve().parents[3]
            if project_root is None
            else project_root
        )
        if not isinstance(root, Path) or not root.is_absolute():
            raise ValueError("project_root must be an absolute pathlib.Path")
        if (snakemake_executable is None) != (conda_prefix is None):
            raise ValueError(
                "snakemake_executable and conda_prefix must be configured together"
            )
        if snakemake_executable is not None:
            if (
                not isinstance(snakemake_executable, Path)
                or not snakemake_executable.is_absolute()
                or snakemake_executable.name != "snakemake"
                or any(
                    character in str(snakemake_executable)
                    for character in ("\x00", "\n", "\r")
                )
            ):
                raise ValueError("snakemake_executable is invalid")
            if (
                not isinstance(conda_prefix, Path)
                or not conda_prefix.is_absolute()
                or any(
                    character in str(conda_prefix) for character in ("\x00", "\n", "\r")
                )
            ):
                raise ValueError("conda_prefix is invalid")
        object.__setattr__(self, "project_root", root)
        object.__setattr__(self, "snakemake_executable", snakemake_executable)
        object.__setattr__(self, "conda_prefix", conda_prefix)

    def _scientific_runtime(
        self,
        workspace: Path,
    ) -> tuple[str, tuple[str, ...], dict[str, str]] | Result[CommandSpec]:
        executable = self.snakemake_executable
        prefix = self.conda_prefix
        if executable is None or prefix is None:
            return "snakemake", (), {}
        runtime_root = executable.parent.parent.parent
        mamba_root = runtime_root / "mamba-root"
        conda_executable = executable.parent / "conda"
        activate = mamba_root / "bin" / "activate"
        micromamba = runtime_root / "runner" / "libexec" / "micromamba"
        try:
            observed_files = tuple(
                (path, path.lstat())
                for path in (executable, conda_executable, activate, micromamba)
            )
            observed_directories = tuple(
                (path, path.lstat()) for path in (runtime_root, mamba_root, prefix)
            )
            if (
                executable.parent.parent != runtime_root / "runner"
                or prefix != runtime_root / "conda-envs"
                or any(
                    path.is_symlink()
                    or not path.is_file()
                    or not os.access(path, os.X_OK)
                    or witness.st_mode & 0o022
                    for path, witness in observed_files
                )
                or any(
                    path.is_symlink() or not path.is_dir() or witness.st_mode & 0o022
                    for path, witness in observed_directories
                )
            ):
                raise OSError
        except OSError:
            return Result.failure(
                [
                    Issue(
                        code="COMMAND_BUILD_SCIENTIFIC_RUNTIME_UNAVAILABLE",
                        message="The admitted scientific runtime is unavailable.",
                        severity="error",
                        path="workflow",
                        source="command_builder",
                    )
                ]
            )
        path = ":".join(
            (
                str(executable.parent),
                "/usr/sbin",
                "/usr/bin",
                "/sbin",
                "/bin",
            )
        )
        return (
            str(executable),
            (
                "--use-conda",
                "--conda-prefix",
                str(prefix),
                "--conda-base-path",
                str(mamba_root),
                "--conda-frontend",
                "conda",
            ),
            {
                "PATH": path,
                "CONDA_DEFAULT_ENV": "",
                "CONDA_EXE": str(conda_executable),
                "CONDA_PREFIX": "",
                "CONDA_SHLVL": "0",
                "HOME": str(workspace),
                "MAMBA_ROOT_PREFIX": str(mamba_root),
                "PYTHONDONTWRITEBYTECODE": "1",
                "PYTHONNOUSERSITE": "1",
                "TMPDIR": str(workspace),
                "XDG_CACHE_HOME": str(workspace / ".snakemake"),
                "_CONDA_EXE": str(conda_executable),
                "_CONDA_ROOT": str(mamba_root),
            },
        )


def _resolve_config_path(
    base_dir: Path,
    workspace_plan: WorkspacePlan,
) -> Path | Result[Any]:
    """Find config/config.yaml in the workspace plan and resolve it safely."""
    for index, (file_path, _) in enumerate(workspace_plan.files):
        if file_path == "config/config.yaml":
            policy = WorkspacePathPolicy(base_dir=base_dir)
            try:
                return policy.resolve(file_path)
            except Exception:
                return Result.failure(
                    [
                        Issue(
                            code="COMMAND_BUILD_MISSING_CONFIG",
                            message="Workspace config file could not be resolved.",
                            severity="error",
                            path=f"workspace_plan.files[{index}]",
                            source="command_builder",
                        )
                    ]
                )

    return Result.failure(
        [
            Issue(
                code="COMMAND_BUILD_MISSING_CONFIG",
                message="Workspace plan must include config/config.yaml.",
                severity="error",
                path="workspace_plan.files",
                source="command_builder",
            )
        ]
    )


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _workspace_digest(plan: WorkspacePlan, cores: int) -> str:
    return hashlib.sha256(
        _json_bytes(
            {
                "cores": cores,
                "directories": plan.directories,
                "files": [
                    [path, hashlib.sha256(content).hexdigest()]
                    for path, content in plan.files
                    if path != EXECUTION_CONFIG_PATH
                ],
            }
        )
    ).hexdigest()


def execution_config_bytes(plan: WorkspacePlan, cores: int) -> bytes:
    """Seal the ENCODE execution settings and the corresponding planned files."""
    return _json_bytes(
        {
            "schema_version": EXECUTION_CONFIG_SCHEMA_VERSION,
            "cores": cores,
            "workspace_contract_sha256": _workspace_digest(plan, cores),
        }
    )


def _unique_json_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
    value: dict[str, Any] = {}
    for key, item in pairs:
        if key in value:
            raise ValueError("duplicate execution configuration key")
        value[key] = item
    return value


class _InvalidCores(ValueError):
    pass


def _resolve_cores(plan: WorkspacePlan, workspace: Path) -> int:
    matches = [content for path, content in plan.files if path == EXECUTION_CONFIG_PATH]
    if len(matches) != 1 or len(matches[0]) > 1024:
        raise ValueError("execution configuration is missing or invalid")
    raw = matches[0]
    value = json.loads(raw, object_pairs_hook=_unique_json_object)
    if not isinstance(value, dict) or set(value) != {
        "schema_version",
        "cores",
        "workspace_contract_sha256",
    }:
        raise ValueError("execution configuration fields are invalid")
    cores = value["cores"]
    if type(cores) is not int or not 1 <= cores <= MAX_CORES:
        raise _InvalidCores("execution configuration cores are invalid")
    if value["schema_version"] != EXECUTION_CONFIG_SCHEMA_VERSION or not isinstance(
        value["workspace_contract_sha256"], str
    ):
        raise ValueError("execution configuration identity is invalid")
    if raw != execution_config_bytes(plan, cores):
        raise ValueError("execution configuration identity mismatch")
    path = workspace / EXECUTION_CONFIG_PATH
    if path.is_symlink():
        raise ValueError("execution configuration must not be a symlink")
    path = WorkspacePathPolicy(base_dir=workspace).resolve(EXECUTION_CONFIG_PATH)
    if path.exists():
        if (
            not path.is_file()
            or path.stat().st_size != len(raw)
            or path.read_bytes() != raw
        ):
            raise ValueError("materialized execution configuration mismatch")
    return cores


def build_encode_command(
    plan: WorkspacePlan,
    workspace: str | Path,
    *,
    binding: EncodeExecutionBinding,
) -> Result[CommandSpec]:
    """Build the ENCODE command from a revalidated, adapter-owned workspace plan."""
    if not isinstance(plan, WorkspacePlan):
        return _execution_config_failure()
    try:
        base_dir = Path(workspace)
    except (TypeError, ValueError):
        return _execution_config_failure()
    if not base_dir.is_absolute():
        return _execution_config_failure()
    snakefile = binding.project_root / "workflow" / "Snakefile"
    if not snakefile.is_file():
        return Result.failure(
            [
                Issue(
                    code="COMMAND_BUILD_SNAKEFILE_NOT_FOUND",
                    message="Bundled Snakefile was not found.",
                    severity="error",
                    path="workflow",
                    source="command_builder",
                )
            ]
        )
    config_path = _resolve_config_path(base_dir, plan)
    if isinstance(config_path, Result):
        return config_path
    try:
        cores = _resolve_cores(plan, base_dir)
    except _InvalidCores:
        return Result.failure(
            [
                Issue(
                    code="COMMAND_BUILD_INVALID_CORES",
                    message="cores must be a positive integer.",
                    severity="error",
                    path="plan.inputs_snapshot.options.cores",
                    source="command_builder",
                )
            ]
        )
    except (OSError, TypeError, ValueError, UnicodeError):
        return _execution_config_failure()
    runtime = binding._scientific_runtime(base_dir)
    if isinstance(runtime, Result):
        return runtime
    executable, runtime_arguments, environment = runtime
    argv = (
        executable,
        "--snakefile",
        str(snakefile),
        "--directory",
        str(base_dir),
        "--configfile",
        str(config_path),
        "--cores",
        str(cores),
        *runtime_arguments,
    )
    return Result.success(
        CommandSpec(
            argv=argv,
            cwd=None,
            env=environment,
            preflight_argv=argv + ("-n",),
        )
    )


def _execution_config_failure() -> Result[CommandSpec]:
    return Result.failure(
        [
            Issue(
                code="ENCODE_EXECUTION_CONFIG_INVALID",
                message="The planned execution configuration could not be verified.",
                severity="error",
                path="workspace_plan",
                source="adapter",
            )
        ]
    )
