"""Immutable, offline-first runtime assets for pinned nf-core/rnaseq.

The platform does not download, unpack, or repair any asset at run time. An
operator stages the exact source tree, Nextflow distribution, plugin archive
and expanded plugin tree, plus raw distribution manifests and content-addressed
Docker archives. This module verifies that closed set, then checks exact images
in one server-owned local Docker daemon. It never pulls or loads an image.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping
from dataclasses import dataclass
import hashlib
from importlib import resources
import io
import json
import os
from pathlib import Path, PurePosixPath
import selectors
import stat
import subprocess
import tarfile
import time
from typing import Any
import zipfile

from jsonschema import Draft202012Validator

from encode_pipeline.platform.results import Issue, Result


SOURCE_MANIFEST_FILE = "source-manifest-3.26.0.json"
SOURCE_MANIFEST_SHA256 = (
    "dc75d105ad26b381197268ef67c44da3107b694e788a3c693237924a86ead774"
)
SOURCE_TREE_SHA256 = "4f779bd8934d41896fbf137ff31158b02daacd7efbb40ed8cee55b9c8f757722"
SOURCE_FILE_COUNT = 829

RUNTIME_IDENTITY_FILE = "runtime-identity-3.26.0.json"
RUNTIME_IDENTITY_SHA256 = (
    "61209ee930b7df7819e8a63e58a87f5fed12477deff3451a77981af1dd4aad8b"
)
CONTAINER_LOCK_SCHEMA_FILE = "container-availability-lock-1.0.0.schema.json"
CONTAINER_LOCK_SCHEMA_SHA256 = (
    "1774c6dfbdefbb8d03ca99a76ad88e893a07bc913d2534855911353d210fa636"
)
CONTAINER_INVENTORY_FILE = "container-inventory-3.26.0.json"
CONTAINER_INVENTORY_SHA256 = (
    "3cae9f36c2b872958dd82b2233cc178aa53fe7ab6111f7bdb0346f8e2bbbe9cf"
)
CONTAINER_INVENTORY_ENTRIES_SHA256 = (
    "684ad35e9b3af3c0b825f6976a9568e34ad13ba62e1e1a0e316853c263f51365"
)
CONTAINER_PROCESS_COUNT = 56
CONTAINER_PROCESS_AUDIT_FILE = "container-process-audit-3.26.0.json"
CONTAINER_PROCESS_AUDIT_SHA256 = (
    "45a2cabd6ec74de359f8d53faa2ca9423ac97cdd89cfa04e2c79c005418fd1c2"
)
CONTAINER_PROCESS_AUDIT_PROCESSES_SHA256 = (
    "4acefb6a517a805845bc034943108ec8e40a6db9e62026b17c62b27a020cd71c"
)
CONTAINER_PROCESS_UNIVERSE_COUNT = 78
CONTAINER_EXCLUDED_PROCESS_COUNT = 22
CONTAINER_CONFIG_ASSIGNMENT_COUNT = 78
CONTAINER_CONFIG_ASSIGNMENTS_SHA256 = (
    "3cf6992db452554a0e93296a811cfdc6d8d7977ffc298f9fbc26fbce35b384d6"
)
CONTAINER_DEFAULT_CONFIG_ASSIGNMENT_COUNT = 1
CONTAINER_DEFAULT_CONFIG_ASSIGNMENTS_SHA256 = (
    "c407972e3fd7a293da2a1ebee45e0614ec6684a1140e58c5db41a57e08873c34"
)
CONTAINER_DEFAULT_CONFIG_DENIED_SELECTORS = (".*PARABRICKS_STARGENOMEGENERATE",)
CONTAINER_RESERVED_DEFAULT_DENY_LABEL = "helixweave_verified_container_v1"

NFCORE_RNASEQ_COMMIT = "e7ca46272c8f9d5ceee3f71759f4ba551d3217a4"
NEXTFLOW_VERSION = "25.04.3"
NEXTFLOW_SHA256 = "53c232cdd8a9419d2c205dc7c6c4dd2646182c997300e6439a453099e28aa21a"
NF_SCHEMA_VERSION = "2.5.1"
NF_SCHEMA_ARCHIVE_SHA256 = (
    "f3833d0c29a51dc5e3759e00a6af87fcb1e989c94d1c1e7995d06e3f7090c461"
)

# The execution overlay must preserve these exact controls. They are exported
# for command composition; the availability lock independently asserts them.
NEXTFLOW_OFFLINE_ENV = {"NXF_OFFLINE": "true"}
DOCKER_REQUIRED_RUN_OPTIONS = ("--pull=never", "--network=none")

_CONTRACT_PACKAGE = "encode_pipeline.contracts.nfcore_rnaseq"
_MAX_CONTRACT_BYTES = 2 * 1024 * 1024
_MAX_SOURCE_FILE_BYTES = 128 * 1024 * 1024
_MAX_CONTAINER_LOCK_BYTES = 2 * 1024 * 1024
_MAX_CONTAINER_ARCHIVE_BYTES = 100 * 1024 * 1024 * 1024
_MAX_DISTRIBUTION_MANIFEST_BYTES = 16 * 1024 * 1024
_MAX_IMAGE_CONFIG_BYTES = 16 * 1024 * 1024
_MAX_DOCKER_ARCHIVE_ENTRIES = 65_536
_MAX_DOCKER_INSPECT_BYTES = 8 * 1024 * 1024
_DOCKER_INSPECT_TIMEOUT_SECONDS = 15.0
_READ_CHUNK_BYTES = 1024 * 1024
_DOCKER_MANIFEST_MEDIA_TYPES = frozenset(
    {
        "application/vnd.docker.distribution.manifest.v2+json",
        "application/vnd.oci.image.manifest.v1+json",
    }
)


@dataclass(frozen=True)
class RuntimeAssetBinding:
    """Operator-owned locations for one closed, offline runtime asset set."""

    root: Path
    source_tree: str = "source/nf-core-rnaseq-e7ca46272c8f9d5ceee3f71759f4ba551d3217a4"
    nextflow_executable: str = "nextflow/nextflow-25.04.3-dist"
    plugin_archive: str = "plugins/nf-schema-2.5.1.zip"
    plugin_meta: str = "plugins/nf-schema-2.5.1-meta.json"
    plugin_tree: str = "plugins/nf-schema-2.5.1"
    container_lock: str = "containers/availability-lock.json"
    container_root: str = "containers/assets"
    docker_executable: Path = Path("/usr/bin/docker")
    docker_socket: Path = Path("/var/run/docker.sock")

    def __post_init__(self) -> None:
        root = Path(self.root)
        if not root.is_absolute():
            raise ValueError("Runtime asset root must be absolute")
        object.__setattr__(self, "root", root)
        for name in (
            "source_tree",
            "nextflow_executable",
            "plugin_archive",
            "plugin_meta",
            "plugin_tree",
            "container_lock",
            "container_root",
        ):
            _validate_relative_path(getattr(self, name))
        for name in ("docker_executable", "docker_socket"):
            value = getattr(self, name)
            if not isinstance(value, Path) or not value.is_absolute():
                raise ValueError(f"{name} must be an absolute pathlib.Path")
            rendered = str(value)
            if (
                rendered != str(Path(rendered))
                or any(character in rendered for character in ("\x00", "\n", "\r"))
                or any(part == ".." for part in value.parts)
            ):
                raise ValueError(f"{name} must be a canonical absolute path")
        if self.docker_executable.name != "docker" or ":" in str(
            self.docker_executable.parent
        ):
            raise ValueError("docker_executable must be a fixed Docker CLI path")


@dataclass(frozen=True)
class VerifiedContainerAsset:
    """One verified process-to-image-to-local-archive binding."""

    process: str
    image: str
    oci_digest: str
    local_asset: Path
    local_sha256: str
    size_bytes: int
    distribution_manifest: Path
    distribution_manifest_sha256: str
    config_digest: str
    runtime_image: str
    rootfs_diff_ids: tuple[str, ...]


@dataclass(frozen=True)
class VerifiedRuntimeAssets:
    """Paths and identities safe to compose into a server-owned command."""

    root: Path
    source_tree: Path
    nextflow_executable: Path
    plugin_root: Path
    plugin_archive: Path
    plugin_meta: Path
    plugin_tree: Path
    container_lock: Path
    containers: tuple[VerifiedContainerAsset, ...]
    source_tree_sha256: str
    runtime_identity_sha256: str
    nextflow_sha256: str
    plugin_archive_sha256: str
    plugin_tree_sha256: str
    container_inventory_sha256: str
    container_lock_sha256: str
    container_process_audit_sha256: str = CONTAINER_PROCESS_AUDIT_SHA256


@dataclass(frozen=True)
class RuntimeAssetDoctorReport:
    """Redacted operator diagnostic; it intentionally contains no identities."""

    ready: bool
    issues: tuple[Issue, ...]


class _AssetFault(Exception):
    def __init__(self, component: str, reason: str) -> None:
        self.component = component
        self.reason = reason
        super().__init__(reason)


@dataclass(frozen=True)
class _RuntimeAssetContract:
    """Private injection seam for tiny filesystem-boundary tests."""

    identity: Mapping[str, Any]
    source_manifest: Mapping[str, Any]
    container_schema: Mapping[str, Any]
    container_inventory: Mapping[str, Any]
    container_process_audit: Mapping[str, Any]


def verify_runtime_assets(
    binding: RuntimeAssetBinding,
    *,
    _contract: _RuntimeAssetContract | None = None,
    _docker_probe: Callable[[RuntimeAssetBinding, tuple[str, ...]], bytes]
    | None = None,
) -> Result[VerifiedRuntimeAssets]:
    """Verify every required runtime asset without fetching or repairing it."""
    if _contract is None:
        try:
            contract = _load_runtime_contract()
            _validate_embedded_contracts(
                contract.identity,
                contract.source_manifest,
                contract.container_inventory,
                contract.container_process_audit,
            )
        except (OSError, ValueError, json.JSONDecodeError) as exc:
            del exc
            return Result.failure((_issue("contract", "invalid"),))
    else:
        contract = _contract

    identity = contract.identity
    source_manifest = contract.source_manifest
    container_schema = contract.container_schema
    container_inventory = contract.container_inventory
    container_process_audit = contract.container_process_audit
    reserved_label = container_process_audit.get("reserved_default_deny_label")
    if not _valid_reserved_label(reserved_label):
        return Result.failure((_issue("contract", "invalid"),))
    assert isinstance(reserved_label, str)

    try:
        root_fd = _open_root(binding.root)
    except _AssetFault as fault:
        return Result.failure((_issue(fault.component, fault.reason),))

    issues: list[Issue] = []
    source_tree_sha256: str | None = None
    plugin_tree_sha256: str | None = None
    container_lock_sha256: str | None = None
    containers: tuple[VerifiedContainerAsset, ...] = ()
    try:
        try:
            source_tree_sha256 = _verify_source_tree(
                root_fd,
                binding.source_tree,
                source_manifest,
                forbidden_label=reserved_label,
            )
        except _AssetFault as fault:
            issues.append(_issue(fault.component, fault.reason))

        try:
            _verify_nextflow(root_fd, binding.nextflow_executable, identity)
        except _AssetFault as fault:
            issues.append(_issue(fault.component, fault.reason))

        try:
            plugin_tree_sha256 = _verify_plugin(
                root_fd,
                binding,
                identity,
            )
        except _AssetFault as fault:
            issues.append(_issue(fault.component, fault.reason))

        try:
            container_lock_sha256, containers = _verify_container_lock(
                root_fd,
                binding,
                container_schema,
                container_inventory,
            )
        except _AssetFault as fault:
            issues.append(_issue(fault.component, fault.reason))
    finally:
        os.close(root_fd)

    if issues:
        return Result.failure(issues)
    if (
        source_tree_sha256 is None
        or plugin_tree_sha256 is None
        or container_lock_sha256 is None
    ):
        return Result.failure((_issue("contract", "invalid"),))
    try:
        _verify_docker_availability(
            binding,
            containers,
            probe=_docker_probe,
        )
    except _AssetFault as fault:
        return Result.failure((_issue(fault.component, fault.reason),))

    return Result.success(
        VerifiedRuntimeAssets(
            root=binding.root,
            source_tree=binding.root / binding.source_tree,
            nextflow_executable=binding.root / binding.nextflow_executable,
            plugin_root=(binding.root / binding.plugin_tree).parent,
            plugin_archive=binding.root / binding.plugin_archive,
            plugin_meta=binding.root / binding.plugin_meta,
            plugin_tree=binding.root / binding.plugin_tree,
            container_lock=binding.root / binding.container_lock,
            containers=containers,
            source_tree_sha256=source_tree_sha256,
            runtime_identity_sha256=RUNTIME_IDENTITY_SHA256,
            nextflow_sha256=identity["nextflow"]["sha256"],
            plugin_archive_sha256=identity["plugins"][0]["archive_sha256"],
            plugin_tree_sha256=plugin_tree_sha256,
            container_inventory_sha256=_canonical_sha256(
                container_inventory["entries"],
                include_executable=False,
            ),
            container_lock_sha256=container_lock_sha256,
            container_process_audit_sha256=CONTAINER_PROCESS_AUDIT_SHA256,
        )
    )


def doctor_runtime_assets(
    binding: RuntimeAssetBinding,
    *,
    _contract: _RuntimeAssetContract | None = None,
    _docker_probe: Callable[[RuntimeAssetBinding, tuple[str, ...]], bytes]
    | None = None,
) -> RuntimeAssetDoctorReport:
    """Return a path- and digest-free readiness report for operators."""
    result = verify_runtime_assets(
        binding,
        _contract=_contract,
        _docker_probe=_docker_probe,
    )
    return RuntimeAssetDoctorReport(ready=result.is_success, issues=result.issues)


def _load_runtime_contract() -> _RuntimeAssetContract:
    identity = _load_embedded_json(
        RUNTIME_IDENTITY_FILE,
        RUNTIME_IDENTITY_SHA256,
    )
    source_manifest = _load_embedded_json(
        SOURCE_MANIFEST_FILE,
        SOURCE_MANIFEST_SHA256,
    )
    container_schema = _load_embedded_json(
        CONTAINER_LOCK_SCHEMA_FILE,
        CONTAINER_LOCK_SCHEMA_SHA256,
    )
    container_inventory = _load_embedded_json(
        CONTAINER_INVENTORY_FILE,
        CONTAINER_INVENTORY_SHA256,
    )
    container_process_audit = _load_embedded_json(
        CONTAINER_PROCESS_AUDIT_FILE,
        CONTAINER_PROCESS_AUDIT_SHA256,
    )
    if not all(
        isinstance(value, Mapping)
        for value in (
            identity,
            source_manifest,
            container_schema,
            container_inventory,
            container_process_audit,
        )
    ):
        raise ValueError("embedded runtime contracts must be objects")
    return _RuntimeAssetContract(
        identity=identity,
        source_manifest=source_manifest,
        container_schema=container_schema,
        container_inventory=container_inventory,
        container_process_audit=container_process_audit,
    )


def _load_embedded_json(filename: str, expected_sha256: str) -> Any:
    content = resources.files(_CONTRACT_PACKAGE).joinpath(filename).read_bytes()
    if len(content) > _MAX_CONTRACT_BYTES:
        raise ValueError("embedded contract is too large")
    if hashlib.sha256(content).hexdigest() != expected_sha256:
        raise ValueError("embedded contract identity differs")
    return _strict_json_loads(content)


def _validate_embedded_contracts(
    identity: Any,
    source_manifest: Any,
    container_inventory: Any,
    container_process_audit: Any,
) -> None:
    if (
        not isinstance(identity, Mapping)
        or not isinstance(source_manifest, Mapping)
        or not isinstance(container_inventory, Mapping)
        or not isinstance(container_process_audit, Mapping)
    ):
        raise ValueError("embedded contract must be an object")
    source = identity.get("source")
    nextflow = identity.get("nextflow")
    plugins = identity.get("plugins")
    network = identity.get("network_policy")
    if not isinstance(source, Mapping) or not isinstance(nextflow, Mapping):
        raise ValueError("embedded identity is incomplete")
    if not isinstance(plugins, list) or len(plugins) != 1:
        raise ValueError("embedded plugin identity is incomplete")
    if not isinstance(network, Mapping) or any(network.values()):
        raise ValueError("runtime network policy must deny every fetch")
    if source.get("commit") != NFCORE_RNASEQ_COMMIT:
        raise ValueError("source commit is invalid")
    if source.get("tree_sha256") != SOURCE_TREE_SHA256:
        raise ValueError("source tree is invalid")
    if source_manifest.get("tree_sha256") != SOURCE_TREE_SHA256:
        raise ValueError("source manifest is invalid")
    if source_manifest.get("file_count") != SOURCE_FILE_COUNT:
        raise ValueError("source manifest count is invalid")
    if nextflow.get("version") != NEXTFLOW_VERSION:
        raise ValueError("Nextflow version is invalid")
    if nextflow.get("sha256") != NEXTFLOW_SHA256:
        raise ValueError("Nextflow identity is invalid")
    plugin = plugins[0]
    if not isinstance(plugin, Mapping):
        raise ValueError("plugin identity is invalid")
    if plugin.get("version") != NF_SCHEMA_VERSION:
        raise ValueError("plugin version is invalid")
    if plugin.get("archive_sha256") != NF_SCHEMA_ARCHIVE_SHA256:
        raise ValueError("plugin identity is invalid")
    execution_policy = identity.get("execution_policy")
    if not isinstance(execution_policy, Mapping):
        raise ValueError("execution policy is invalid")
    if execution_policy.get("nxf_offline") is not True:
        raise ValueError("execution policy is invalid")
    if tuple(execution_policy.get("required_docker_run_options", ())) != (
        DOCKER_REQUIRED_RUN_OPTIONS
    ):
        raise ValueError("execution policy is invalid")
    containers = identity.get("containers")
    if (
        not isinstance(containers, Mapping)
        or containers.get("availability_lock_schema_sha256")
        != CONTAINER_LOCK_SCHEMA_SHA256
        or containers.get("inventory_sha256") != CONTAINER_INVENTORY_SHA256
        or containers.get("inventory_entries_sha256")
        != CONTAINER_INVENTORY_ENTRIES_SHA256
        or containers.get("required_process_count") != CONTAINER_PROCESS_COUNT
        or containers.get("process_audit") != CONTAINER_PROCESS_AUDIT_FILE
        or containers.get("process_audit_sha256") != CONTAINER_PROCESS_AUDIT_SHA256
        or containers.get("process_audit_processes_sha256")
        != CONTAINER_PROCESS_AUDIT_PROCESSES_SHA256
        or containers.get("process_universe_count") != CONTAINER_PROCESS_UNIVERSE_COUNT
        or containers.get("excluded_process_count") != CONTAINER_EXCLUDED_PROCESS_COUNT
        or containers.get("config_container_assignment_count")
        != CONTAINER_CONFIG_ASSIGNMENT_COUNT
        or containers.get("config_container_assignments_sha256")
        != CONTAINER_CONFIG_ASSIGNMENTS_SHA256
        or containers.get("default_config_container_assignment_count")
        != CONTAINER_DEFAULT_CONFIG_ASSIGNMENT_COUNT
        or containers.get("default_config_container_assignments_sha256")
        != CONTAINER_DEFAULT_CONFIG_ASSIGNMENTS_SHA256
        or containers.get("reserved_default_deny_label")
        != CONTAINER_RESERVED_DEFAULT_DENY_LABEL
        or containers.get("asset_format") != "docker-archive+distribution-manifest"
        or containers.get("runtime_pull_allowed") is not False
    ):
        raise ValueError("container identity is invalid")
    inventory_entries = container_inventory.get("entries")
    if (
        container_inventory.get("commit") != NFCORE_RNASEQ_COMMIT
        or container_inventory.get("process_count") != CONTAINER_PROCESS_COUNT
        or container_inventory.get("entries_sha256")
        != CONTAINER_INVENTORY_ENTRIES_SHA256
        or not isinstance(inventory_entries, list)
        or len(inventory_entries) != CONTAINER_PROCESS_COUNT
    ):
        raise ValueError("container inventory is invalid")
    source_files = {
        entry["path"]: entry["sha256"] for entry in source_manifest["files"]
    }
    seen_processes: set[str] = set()
    for entry in inventory_entries:
        if not isinstance(entry, Mapping):
            raise ValueError("container inventory is invalid")
        process = entry.get("process")
        source_file = entry.get("source_file")
        if (
            not _valid_process_name(process)
            or process in seen_processes
            or not isinstance(source_file, str)
            or source_files.get(source_file) != entry.get("source_file_sha256")
            or not _valid_image_coordinate(entry.get("image_coordinate"))
        ):
            raise ValueError("container inventory is invalid")
        seen_processes.add(process)
    canonical_entries = _canonical_sha256(inventory_entries, include_executable=False)
    if canonical_entries != CONTAINER_INVENTORY_ENTRIES_SHA256:
        raise ValueError("container inventory is invalid")
    if "RIBODETECTOR" in seen_processes:
        raise ValueError("unsupported container is present")
    _validate_container_process_audit(
        container_process_audit,
        inventory_entries=inventory_entries,
        source_files=source_files,
    )


def _validate_container_process_audit(
    audit: Mapping[str, Any],
    *,
    inventory_entries: list[Any],
    source_files: Mapping[str, str],
) -> None:
    required_top_level = {
        "schema_version",
        "project",
        "release",
        "commit",
        "source_tree_sha256",
        "reserved_default_deny_label",
        "process_universe_count",
        "included_process_count",
        "excluded_process_count",
        "config_container_assignment_count",
        "config_container_assignments_sha256",
        "default_config_container_assignment_count",
        "default_config_container_assignments_sha256",
        "default_config_container_assignments",
        "processes_sha256",
        "processes",
    }
    if set(audit) != required_top_level or (
        audit.get("schema_version") != "1.1.0"
        or audit.get("project") != "nf-core/rnaseq"
        or audit.get("release") != "3.26.0"
        or audit.get("commit") != NFCORE_RNASEQ_COMMIT
        or audit.get("source_tree_sha256") != SOURCE_TREE_SHA256
        or audit.get("reserved_default_deny_label")
        != CONTAINER_RESERVED_DEFAULT_DENY_LABEL
        or audit.get("process_universe_count") != CONTAINER_PROCESS_UNIVERSE_COUNT
        or audit.get("included_process_count") != CONTAINER_PROCESS_COUNT
        or audit.get("excluded_process_count") != CONTAINER_EXCLUDED_PROCESS_COUNT
        or audit.get("config_container_assignment_count")
        != CONTAINER_CONFIG_ASSIGNMENT_COUNT
        or audit.get("config_container_assignments_sha256")
        != CONTAINER_CONFIG_ASSIGNMENTS_SHA256
        or audit.get("default_config_container_assignment_count")
        != CONTAINER_DEFAULT_CONFIG_ASSIGNMENT_COUNT
        or audit.get("default_config_container_assignments_sha256")
        != CONTAINER_DEFAULT_CONFIG_ASSIGNMENTS_SHA256
        or audit.get("processes_sha256") != CONTAINER_PROCESS_AUDIT_PROCESSES_SHA256
    ):
        raise ValueError("container process audit is invalid")
    default_config_assignments = audit.get("default_config_container_assignments")
    if (
        not isinstance(default_config_assignments, list)
        or len(default_config_assignments) != CONTAINER_DEFAULT_CONFIG_ASSIGNMENT_COUNT
        or _canonical_value_sha256(default_config_assignments)
        != CONTAINER_DEFAULT_CONFIG_ASSIGNMENTS_SHA256
    ):
        raise ValueError("container process audit is invalid")
    processes = audit.get("processes")
    if (
        not isinstance(processes, list)
        or len(processes) != CONTAINER_PROCESS_UNIVERSE_COUNT
        or _canonical_value_sha256(processes)
        != CONTAINER_PROCESS_AUDIT_PROCESSES_SHA256
    ):
        raise ValueError("container process audit is invalid")

    inventory_by_process = {
        entry["process"]: entry
        for entry in inventory_entries
        if isinstance(entry, Mapping)
    }
    denied_selectors: list[str] = []
    for entry in default_config_assignments:
        if not isinstance(entry, Mapping) or set(entry) != {
            "source_file",
            "source_file_sha256",
            "selector",
            "image_coordinates",
            "matched_aliases",
            "load_scope",
            "platform_override",
            "override_reason",
        }:
            raise ValueError("container process audit is invalid")
        selector = entry.get("selector")
        source_file = entry.get("source_file")
        coordinates = entry.get("image_coordinates")
        matched_aliases = entry.get("matched_aliases")
        if (
            selector not in CONTAINER_DEFAULT_CONFIG_DENIED_SELECTORS
            or not isinstance(source_file, str)
            or source_file != "conf/modules/prepare_genome.config"
            or source_files.get(source_file) != entry.get("source_file_sha256")
            or not isinstance(coordinates, list)
            or not all(isinstance(value, str) for value in coordinates)
            or coordinates != sorted(set(coordinates))
            or not coordinates
            or not all(_valid_image_coordinate(value) for value in coordinates)
            or not isinstance(matched_aliases, list)
            or not matched_aliases
            or entry.get("load_scope") != "default"
            or entry.get("platform_override") != "deny"
            or entry.get("override_reason") != "unsupported_parabricks_alias"
        ):
            raise ValueError("container process audit is invalid")
        assert isinstance(selector, str)
        for alias in matched_aliases:
            if not isinstance(alias, Mapping) or set(alias) != {
                "alias",
                "source_process",
                "source_file",
                "source_file_sha256",
            }:
                raise ValueError("container process audit is invalid")
            alias_name = alias.get("alias")
            alias_source = alias.get("source_file")
            if (
                alias_name != "PARABRICKS_STARGENOMEGENERATE"
                or alias.get("source_process") != "STAR_GENOMEGENERATE"
                or "STAR_GENOMEGENERATE" not in inventory_by_process
                or not isinstance(alias_source, str)
                or alias_source != "subworkflows/local/prepare_genome/main.nf"
                or source_files.get(alias_source) != alias.get("source_file_sha256")
            ):
                raise ValueError("container process audit is invalid")
        denied_selectors.append(selector)
    if tuple(denied_selectors) != CONTAINER_DEFAULT_CONFIG_DENIED_SELECTORS:
        raise ValueError("container process audit is invalid")
    seen: set[str] = set()
    included: set[str] = set()
    excluded: set[str] = set()
    expected_order: list[str] = []
    for entry in processes:
        if not isinstance(entry, Mapping):
            raise ValueError("container process audit is invalid")
        disposition = entry.get("disposition")
        required_keys = {
            "process",
            "source_file",
            "source_file_sha256",
            "image_coordinates",
            "disposition",
        }
        if disposition == "excluded":
            required_keys.add("exclusion_reason")
        if set(entry) != required_keys:
            raise ValueError("container process audit is invalid")
        process = entry.get("process")
        source_file = entry.get("source_file")
        coordinates = entry.get("image_coordinates")
        try:
            _validate_relative_path(source_file)
        except ValueError as exc:
            raise ValueError("container process audit is invalid") from exc
        if (
            not _valid_process_name(process)
            or process in seen
            or not isinstance(source_file, str)
            or source_files.get(source_file) != entry.get("source_file_sha256")
            or not isinstance(coordinates, list)
            or not coordinates
            or len(coordinates) > 4
            or coordinates != sorted(set(coordinates))
            or not all(_valid_image_coordinate(value) for value in coordinates)
        ):
            raise ValueError("container process audit is invalid")
        assert isinstance(process, str)
        seen.add(process)
        expected_order.append(process)
        if disposition == "included":
            inventory_entry = inventory_by_process.get(process)
            if inventory_entry is None or (
                source_file != inventory_entry.get("source_file")
                or entry.get("source_file_sha256")
                != inventory_entry.get("source_file_sha256")
                or coordinates != [inventory_entry.get("image_coordinate")]
            ):
                raise ValueError("container process audit is invalid")
            included.add(process)
        elif disposition == "excluded":
            reason = entry.get("exclusion_reason")
            if (
                process in inventory_by_process
                or not isinstance(reason, str)
                or not reason.startswith("unsupported_")
                or not all(
                    char.islower() or char.isdigit() or char == "_" for char in reason
                )
            ):
                raise ValueError("container process audit is invalid")
            excluded.add(process)
        else:
            raise ValueError("container process audit is invalid")
    if (
        expected_order != sorted(expected_order)
        or included != set(inventory_by_process)
        or len(included) != CONTAINER_PROCESS_COUNT
        or len(excluded) != CONTAINER_EXCLUDED_PROCESS_COUNT
        or included & excluded
        or seen != included | excluded
        or "RIBODETECTOR" not in excluded
    ):
        raise ValueError("container process audit is invalid")


def _verify_source_tree(
    root_fd: int,
    relative_tree: str,
    manifest: Mapping[str, Any],
    *,
    forbidden_label: str,
) -> str:
    files = manifest.get("files")
    if (
        not isinstance(files, list)
        or isinstance(manifest.get("file_count"), bool)
        or not isinstance(manifest.get("file_count"), int)
        or len(files) != manifest.get("file_count")
    ):
        raise _AssetFault("source", "contract")
    expected: dict[str, _ExpectedFile] = {}
    for item in files:
        if not isinstance(item, Mapping):
            raise _AssetFault("source", "contract")
        try:
            path = item["path"]
            size_bytes = item["size_bytes"]
            sha256 = item["sha256"]
            executable = item["executable"]
        except KeyError as exc:
            raise _AssetFault("source", "contract") from exc
        _validate_relative_path(path)
        if (
            not isinstance(size_bytes, int)
            or isinstance(size_bytes, bool)
            or size_bytes < 0
            or size_bytes > _MAX_SOURCE_FILE_BYTES
            or not _valid_sha256(sha256)
            or not isinstance(executable, bool)
            or path in expected
        ):
            raise _AssetFault("source", "contract")
        expected[path] = _ExpectedFile(size_bytes, sha256, executable)

    tree_fd = _open_directory_at(root_fd, relative_tree, "source")
    try:
        verified = _verify_exact_tree(
            tree_fd,
            expected,
            "source",
            forbidden_bytes=forbidden_label.encode("ascii"),
        )
    finally:
        os.close(tree_fd)
    digest = _canonical_sha256(verified)
    if digest != manifest.get("tree_sha256"):
        raise _AssetFault("source", "identity")
    return digest


def _verify_nextflow(
    root_fd: int, relative_path: str, identity: Mapping[str, Any]
) -> None:
    nextflow = identity["nextflow"]
    expected = _ExpectedFile(
        size_bytes=nextflow["size_bytes"],
        sha256=nextflow["sha256"],
        executable=True,
    )
    _verify_file_at(root_fd, relative_path, expected, "nextflow")


def _verify_plugin(
    root_fd: int,
    binding: RuntimeAssetBinding,
    identity: Mapping[str, Any],
) -> str:
    plugin = identity["plugins"][0]
    archive = _read_file_at(
        root_fd,
        binding.plugin_archive,
        _ExpectedFile(plugin["archive_size_bytes"], plugin["archive_sha256"]),
        "plugin",
        _MAX_CONTRACT_BYTES,
    )
    meta = _read_file_at(
        root_fd,
        binding.plugin_meta,
        _ExpectedFile(plugin["meta_size_bytes"], plugin["meta_sha256"]),
        "plugin",
        _MAX_CONTRACT_BYTES,
    )

    expected: dict[str, _ExpectedFile] = {}
    try:
        with zipfile.ZipFile(io.BytesIO(archive)) as package:
            for entry in package.infolist():
                if entry.is_dir():
                    continue
                _validate_relative_path(entry.filename)
                mode = (entry.external_attr >> 16) & 0o170000
                if mode == stat.S_IFLNK or entry.flag_bits & 0x1:
                    raise _AssetFault("plugin", "contract")
                if entry.filename in expected:
                    raise _AssetFault("plugin", "contract")
                content = package.read(entry)
                expected[entry.filename] = _ExpectedFile(
                    len(content),
                    hashlib.sha256(content).hexdigest(),
                )
    except (OSError, ValueError, zipfile.BadZipFile, RuntimeError) as exc:
        raise _AssetFault("plugin", "contract") from exc
    meta_name = f"nf-schema-{NF_SCHEMA_VERSION}-meta.json"
    if meta_name in expected:
        raise _AssetFault("plugin", "contract")
    expected[meta_name] = _ExpectedFile(
        len(meta),
        hashlib.sha256(meta).hexdigest(),
    )

    tree_fd = _open_directory_at(root_fd, binding.plugin_tree, "plugin")
    try:
        verified = _verify_exact_tree(tree_fd, expected, "plugin")
    finally:
        os.close(tree_fd)
    digest = _canonical_sha256(verified, include_executable=False)
    if digest != plugin.get("tree_sha256"):
        raise _AssetFault("plugin", "identity")
    return digest


def _verify_container_lock(
    root_fd: int,
    binding: RuntimeAssetBinding,
    schema: Mapping[str, Any],
    inventory: Mapping[str, Any],
) -> tuple[str, tuple[VerifiedContainerAsset, ...]]:
    content = _read_file_at(
        root_fd,
        binding.container_lock,
        expected=None,
        component="containers",
        max_bytes=_MAX_CONTAINER_LOCK_BYTES,
    )
    try:
        lock = _strict_json_loads(content)
        errors = tuple(Draft202012Validator(schema).iter_errors(lock))
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise _AssetFault("containers", "contract") from exc
    if errors or not isinstance(lock, Mapping):
        raise _AssetFault("containers", "contract")

    entries = lock["entries"]
    processes = [entry["process"] for entry in entries]
    if len(processes) != len(set(processes)):
        raise _AssetFault("containers", "contract")
    expected_by_process = {entry["process"]: entry for entry in inventory["entries"]}
    if set(processes) != set(expected_by_process):
        raise _AssetFault("containers", "process_set")

    verified: list[VerifiedContainerAsset] = []
    closure_cache: dict[
        tuple[str, int, str, str, int, str],
        _ContainerArchiveClosure,
    ] = {}
    for entry in sorted(entries, key=lambda value: value["process"]):
        coordinate = expected_by_process[entry["process"]]["image_coordinate"]
        if entry["image"] != f"{coordinate}@{entry['oci_digest']}":
            raise _AssetFault("containers", "identity")
        try:
            _validate_relative_path(entry["local_asset"])
            _validate_relative_path(entry["distribution_manifest_asset"])
        except ValueError as exc:
            raise _AssetFault("containers", "contract") from exc
        relative_asset = _join_relative(
            binding.container_root,
            entry["local_asset"],
        )
        relative_manifest = _join_relative(
            binding.container_root,
            entry["distribution_manifest_asset"],
        )
        cache_key = (
            relative_asset,
            entry["size_bytes"],
            entry["sha256"],
            relative_manifest,
            entry["distribution_manifest_size_bytes"],
            entry["oci_digest"],
        )
        closure = closure_cache.get(cache_key)
        if closure is None:
            distribution = _verify_distribution_manifest(
                root_fd,
                relative_manifest,
                size_bytes=entry["distribution_manifest_size_bytes"],
                oci_digest=entry["oci_digest"],
            )
            closure = _verify_docker_archive(
                root_fd,
                relative_asset,
                expected=_ExpectedFile(
                    entry["size_bytes"],
                    entry["sha256"],
                ),
                distribution=distribution,
            )
            closure_cache[cache_key] = closure
        verified.append(
            VerifiedContainerAsset(
                process=entry["process"],
                image=entry["image"],
                oci_digest=entry["oci_digest"],
                local_asset=binding.root / relative_asset,
                local_sha256=entry["sha256"],
                size_bytes=entry["size_bytes"],
                distribution_manifest=binding.root / relative_manifest,
                distribution_manifest_sha256=entry["oci_digest"].removeprefix(
                    "sha256:"
                ),
                config_digest=closure.config_digest,
                runtime_image=closure.config_digest,
                rootfs_diff_ids=closure.layer_diff_ids,
            )
        )
    return hashlib.sha256(content).hexdigest(), tuple(verified)


@dataclass(frozen=True)
class _DistributionManifest:
    config_digest: str
    config_size_bytes: int
    layer_count: int


@dataclass(frozen=True)
class _ContainerArchiveClosure:
    config_digest: str
    layer_diff_ids: tuple[str, ...]


def _verify_distribution_manifest(
    root_fd: int,
    relative_path: str,
    *,
    size_bytes: int,
    oci_digest: str,
) -> _DistributionManifest:
    digest = _digest_hex(oci_digest)
    content = _read_file_at(
        root_fd,
        relative_path,
        _ExpectedFile(size_bytes=size_bytes, sha256=digest),
        "containers",
        _MAX_DISTRIBUTION_MANIFEST_BYTES,
    )
    try:
        value = _strict_json_loads(content)
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError) as exc:
        raise _AssetFault("containers", "contract") from exc
    if not isinstance(value, Mapping) or value.get("schemaVersion") != 2:
        raise _AssetFault("containers", "contract")
    media_type = value.get("mediaType")
    if media_type is not None and media_type not in _DOCKER_MANIFEST_MEDIA_TYPES:
        raise _AssetFault("containers", "contract")
    config = value.get("config")
    layers = value.get("layers")
    if not isinstance(config, Mapping) or not isinstance(layers, list):
        raise _AssetFault("containers", "contract")
    config_digest = config.get("digest")
    config_size = config.get("size")
    if (
        not _valid_digest(config_digest)
        or isinstance(config_size, bool)
        or not isinstance(config_size, int)
        or config_size <= 0
        or config_size > _MAX_IMAGE_CONFIG_BYTES
        or len(layers) > _MAX_DOCKER_ARCHIVE_ENTRIES
    ):
        raise _AssetFault("containers", "contract")
    for layer in layers:
        if not isinstance(layer, Mapping):
            raise _AssetFault("containers", "contract")
        layer_size = layer.get("size")
        if (
            not _valid_digest(layer.get("digest"))
            or isinstance(layer_size, bool)
            or not isinstance(layer_size, int)
            or layer_size < 0
            or layer_size > _MAX_CONTAINER_ARCHIVE_BYTES
        ):
            raise _AssetFault("containers", "contract")
    return _DistributionManifest(
        config_digest=config_digest,
        config_size_bytes=config_size,
        layer_count=len(layers),
    )


def _verify_docker_archive(
    root_fd: int,
    relative_path: str,
    *,
    expected: "_ExpectedFile",
    distribution: _DistributionManifest,
) -> _ContainerArchiveClosure:
    parent_fd, name = _open_parent(root_fd, relative_path, "containers")
    descriptor = -1
    try:
        descriptor, before = _open_regular_entry(parent_fd, name, "containers")
        if (
            before.st_size != expected.size_bytes
            or before.st_size > _MAX_CONTAINER_ARCHIVE_BYTES
        ):
            raise _AssetFault("containers", "identity")
        digest = _hash_open_descriptor(
            descriptor,
            maximum_bytes=_MAX_CONTAINER_ARCHIVE_BYTES,
            component="containers",
        )
        if digest != expected.sha256:
            raise _AssetFault("containers", "identity")
        os.lseek(descriptor, 0, os.SEEK_SET)
        with os.fdopen(os.dup(descriptor), "rb") as handle:
            closure = _parse_docker_archive(handle, distribution=distribution)
        after = os.fstat(descriptor)
        _check_stable(before, after, "containers")
        return closure
    except _AssetFault:
        raise
    except (OSError, ValueError, json.JSONDecodeError, tarfile.TarError) as exc:
        raise _AssetFault("containers", "contract") from exc
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        os.close(parent_fd)


def _parse_docker_archive(
    handle: Any,
    *,
    distribution: _DistributionManifest,
) -> _ContainerArchiveClosure:
    members: dict[str, tarfile.TarInfo] = {}
    total_declared_bytes = 0
    try:
        with tarfile.open(fileobj=handle, mode="r:") as package:
            for member in package:
                if len(members) >= _MAX_DOCKER_ARCHIVE_ENTRIES:
                    raise _AssetFault("containers", "bounds")
                normalized = member.name[:-1] if member.isdir() else member.name
                _validate_relative_path(normalized)
                if normalized in members:
                    raise _AssetFault("containers", "contract")
                if not member.isdir() and not member.isreg():
                    raise _AssetFault("containers", "file_type")
                if member.size < 0 or member.size > _MAX_CONTAINER_ARCHIVE_BYTES:
                    raise _AssetFault("containers", "bounds")
                total_declared_bytes += member.size
                if total_declared_bytes > _MAX_CONTAINER_ARCHIVE_BYTES:
                    raise _AssetFault("containers", "bounds")
                members[normalized] = member

            manifest_member = _required_regular_tar_member(members, "manifest.json")
            manifest_content = _read_tar_member(
                package,
                manifest_member,
                maximum_bytes=_MAX_DISTRIBUTION_MANIFEST_BYTES,
            )
            manifest = _strict_json_loads(manifest_content)
            if not isinstance(manifest, list) or len(manifest) != 1:
                raise _AssetFault("containers", "contract")
            image = manifest[0]
            if not isinstance(image, Mapping):
                raise _AssetFault("containers", "contract")
            config_path = image.get("Config")
            layer_paths = image.get("Layers")
            repo_tags = image.get("RepoTags")
            if (
                not isinstance(config_path, str)
                or not isinstance(layer_paths, list)
                or not all(isinstance(path, str) for path in layer_paths)
                or repo_tags is not None
                and (
                    not isinstance(repo_tags, list)
                    or not all(isinstance(tag, str) for tag in repo_tags)
                )
            ):
                raise _AssetFault("containers", "contract")
            _validate_relative_path(config_path)
            for layer_path in layer_paths:
                _validate_relative_path(layer_path)
            if len(layer_paths) != distribution.layer_count:
                raise _AssetFault("containers", "identity")

            config_member = _required_regular_tar_member(members, config_path)
            config_content = _read_tar_member(
                package,
                config_member,
                maximum_bytes=_MAX_IMAGE_CONFIG_BYTES,
            )
            if (
                len(config_content) != distribution.config_size_bytes
                or f"sha256:{hashlib.sha256(config_content).hexdigest()}"
                != distribution.config_digest
            ):
                raise _AssetFault("containers", "identity")
            config = _strict_json_loads(config_content)
            if not isinstance(config, Mapping):
                raise _AssetFault("containers", "contract")
            rootfs = config.get("rootfs")
            if not isinstance(rootfs, Mapping) or rootfs.get("type") != "layers":
                raise _AssetFault("containers", "contract")
            diff_ids = rootfs.get("diff_ids")
            if (
                not isinstance(diff_ids, list)
                or len(diff_ids) != len(layer_paths)
                or not all(_valid_digest(value) for value in diff_ids)
            ):
                raise _AssetFault("containers", "contract")
            for layer_path, expected_diff_id in zip(
                layer_paths,
                diff_ids,
                strict=True,
            ):
                layer_member = _required_regular_tar_member(members, layer_path)
                actual_diff_id = f"sha256:{_hash_tar_member(package, layer_member)}"
                if actual_diff_id != expected_diff_id:
                    raise _AssetFault("containers", "identity")
    except _AssetFault:
        raise
    except (OSError, ValueError, json.JSONDecodeError, tarfile.TarError) as exc:
        raise _AssetFault("containers", "contract") from exc
    return _ContainerArchiveClosure(
        config_digest=distribution.config_digest,
        layer_diff_ids=tuple(diff_ids),
    )


def _required_regular_tar_member(
    members: Mapping[str, tarfile.TarInfo],
    name: str,
) -> tarfile.TarInfo:
    member = members.get(name)
    if member is None or not member.isreg():
        raise _AssetFault("containers", "contract")
    return member


def _read_tar_member(
    package: tarfile.TarFile,
    member: tarfile.TarInfo,
    *,
    maximum_bytes: int,
) -> bytes:
    if member.size > maximum_bytes:
        raise _AssetFault("containers", "bounds")
    extracted = package.extractfile(member)
    if extracted is None:
        raise _AssetFault("containers", "contract")
    content = extracted.read(maximum_bytes + 1)
    if len(content) != member.size or len(content) > maximum_bytes:
        raise _AssetFault("containers", "bounds")
    return content


def _hash_tar_member(package: tarfile.TarFile, member: tarfile.TarInfo) -> str:
    extracted = package.extractfile(member)
    if extracted is None:
        raise _AssetFault("containers", "contract")
    digest = hashlib.sha256()
    total = 0
    while True:
        chunk = extracted.read(_READ_CHUNK_BYTES)
        if not chunk:
            break
        total += len(chunk)
        if total > member.size or total > _MAX_CONTAINER_ARCHIVE_BYTES:
            raise _AssetFault("containers", "bounds")
        digest.update(chunk)
    if total != member.size:
        raise _AssetFault("containers", "contract")
    return digest.hexdigest()


def _hash_open_descriptor(
    descriptor: int,
    *,
    maximum_bytes: int,
    component: str,
) -> str:
    digest = hashlib.sha256()
    total = 0
    while True:
        chunk = os.read(descriptor, _READ_CHUNK_BYTES)
        if not chunk:
            break
        total += len(chunk)
        if total > maximum_bytes:
            raise _AssetFault(component, "bounds")
        digest.update(chunk)
    return digest.hexdigest()


def _verify_docker_availability(
    binding: RuntimeAssetBinding,
    containers: tuple[VerifiedContainerAsset, ...],
    *,
    probe: Callable[[RuntimeAssetBinding, tuple[str, ...]], bytes] | None,
) -> None:
    expected: dict[str, tuple[str, ...]] = {}
    for container in containers:
        if (
            not _valid_digest(container.config_digest)
            or container.runtime_image != container.config_digest
            or not isinstance(container.rootfs_diff_ids, tuple)
            or not all(_valid_digest(value) for value in container.rootfs_diff_ids)
        ):
            raise _AssetFault("docker", "contract")
        assert isinstance(container.config_digest, str)
        previous = expected.setdefault(
            container.config_digest,
            container.rootfs_diff_ids,
        )
        if previous != container.rootfs_diff_ids:
            raise _AssetFault("docker", "identity")
    if not expected:
        raise _AssetFault("docker", "contract")
    images = tuple(sorted(expected))
    try:
        output = (probe or _run_docker_inspect)(binding, images)
    except _AssetFault:
        raise
    except Exception as exc:
        raise _AssetFault("docker", "unavailable") from exc
    if not isinstance(output, bytes) or len(output) > _MAX_DOCKER_INSPECT_BYTES:
        raise _AssetFault("docker", "unavailable")
    try:
        inspected = _strict_json_loads(output)
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError) as exc:
        raise _AssetFault("docker", "unavailable") from exc
    if not isinstance(inspected, list) or len(inspected) != len(images):
        raise _AssetFault("docker", "unavailable")
    observed: set[str] = set()
    for value in inspected:
        if not isinstance(value, Mapping):
            raise _AssetFault("docker", "unavailable")
        image_id = value.get("Id")
        rootfs = value.get("RootFS")
        if (
            not isinstance(image_id, str)
            or image_id not in expected
            or image_id in observed
            or not isinstance(rootfs, Mapping)
            or rootfs.get("Type") != "layers"
            or rootfs.get("Layers") != list(expected[image_id])
        ):
            raise _AssetFault("docker", "identity")
        observed.add(image_id)
    if observed != set(expected):
        raise _AssetFault("docker", "identity")


def _run_docker_inspect(
    binding: RuntimeAssetBinding,
    images: tuple[str, ...],
) -> bytes:
    endpoint_identity = _verify_local_docker_endpoint(binding)
    argv = (
        str(binding.docker_executable),
        "--host",
        f"unix://{binding.docker_socket}",
        "image",
        "inspect",
        *images,
    )
    try:
        process = subprocess.Popen(
            argv,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE,
            stderr=subprocess.DEVNULL,
            shell=False,
            env={"HOME": "/", "LANG": "C", "LC_ALL": "C"},
        )
    except OSError as exc:
        raise _AssetFault("docker", "unavailable") from exc
    try:
        output = _read_process_stdout_bounded(process)
    except _AssetFault:
        _kill_probe_process(process)
        raise
    try:
        final_endpoint_identity = _verify_local_docker_endpoint(binding)
    except _AssetFault as exc:
        raise _AssetFault("docker", "race") from exc
    if final_endpoint_identity != endpoint_identity:
        raise _AssetFault("docker", "race")
    return output


def _read_process_stdout_bounded(process: subprocess.Popen[bytes]) -> bytes:
    stdout = process.stdout
    if stdout is None:
        raise _AssetFault("docker", "unavailable")
    selector = selectors.DefaultSelector()
    content = bytearray()
    deadline = time.monotonic() + _DOCKER_INSPECT_TIMEOUT_SECONDS
    try:
        os.set_blocking(stdout.fileno(), False)
        selector.register(stdout, selectors.EVENT_READ)
        reading = True
        while reading:
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                raise _AssetFault("docker", "unavailable")
            events = selector.select(timeout=remaining)
            if not events:
                raise _AssetFault("docker", "unavailable")
            for key, _mask in events:
                try:
                    chunk = os.read(
                        key.fd,
                        min(
                            _READ_CHUNK_BYTES,
                            _MAX_DOCKER_INSPECT_BYTES + 1 - len(content),
                        ),
                    )
                except BlockingIOError:
                    continue
                if not chunk:
                    selector.unregister(stdout)
                    reading = False
                    break
                content.extend(chunk)
                if len(content) > _MAX_DOCKER_INSPECT_BYTES:
                    raise _AssetFault("docker", "unavailable")
        remaining = deadline - time.monotonic()
        if remaining <= 0:
            raise _AssetFault("docker", "unavailable")
        try:
            returncode = process.wait(timeout=remaining)
        except subprocess.TimeoutExpired as exc:
            raise _AssetFault("docker", "unavailable") from exc
        if returncode != 0:
            raise _AssetFault("docker", "unavailable")
        return bytes(content)
    except OSError as exc:
        raise _AssetFault("docker", "unavailable") from exc
    finally:
        selector.close()
        stdout.close()


def _kill_probe_process(process: subprocess.Popen[bytes]) -> None:
    try:
        if process.poll() is None:
            process.kill()
        process.wait(timeout=1.0)
    except (OSError, subprocess.TimeoutExpired):
        pass


def _verify_local_docker_endpoint(
    binding: RuntimeAssetBinding,
) -> tuple[tuple[int, ...], tuple[int, ...]]:
    try:
        executable = os.lstat(binding.docker_executable)
        socket_info = os.lstat(binding.docker_socket)
    except OSError as exc:
        raise _AssetFault("docker", "missing") from exc
    if (
        not stat.S_ISREG(executable.st_mode)
        or executable.st_nlink != 1
        or executable.st_mode & 0o111 == 0
        or executable.st_mode & 0o022 != 0
        or executable.st_uid not in {0, os.geteuid()}
        or not stat.S_ISSOCK(socket_info.st_mode)
        or socket_info.st_nlink != 1
    ):
        raise _AssetFault("docker", "file_type")
    return (
        (
            executable.st_dev,
            executable.st_ino,
            executable.st_mode,
            executable.st_size,
            executable.st_mtime_ns,
            executable.st_ctime_ns,
        ),
        (socket_info.st_dev, socket_info.st_ino, socket_info.st_mode),
    )


def _digest_hex(value: object) -> str:
    if not _valid_digest(value):
        raise _AssetFault("containers", "contract")
    assert isinstance(value, str)
    return value.removeprefix("sha256:")


@dataclass(frozen=True)
class _ExpectedFile:
    size_bytes: int
    sha256: str
    executable: bool | None = None


def _verify_exact_tree(
    tree_fd: int,
    expected: Mapping[str, _ExpectedFile],
    component: str,
    *,
    forbidden_bytes: bytes | None = None,
) -> list[dict[str, object]]:
    expected_directories = {""}
    for path in expected:
        parent = PurePosixPath(path).parent
        while str(parent) != ".":
            expected_directories.add(parent.as_posix())
            parent = parent.parent

    found: set[str] = set()
    verified: list[dict[str, object]] = []

    def walk(directory_fd: int, prefix: str) -> None:
        try:
            names = sorted(os.listdir(directory_fd))
        except OSError as exc:
            raise _AssetFault(component, "unreadable") from exc
        for name in names:
            relative = f"{prefix}/{name}" if prefix else name
            try:
                info = os.stat(name, dir_fd=directory_fd, follow_symlinks=False)
            except OSError as exc:
                raise _AssetFault(component, "race") from exc
            if stat.S_ISDIR(info.st_mode):
                if relative not in expected_directories:
                    raise _AssetFault(component, "file_set")
                child_fd = _open_single_directory(directory_fd, name, component)
                try:
                    walk(child_fd, relative)
                finally:
                    os.close(child_fd)
                continue
            if not stat.S_ISREG(info.st_mode):
                raise _AssetFault(component, "file_type")
            item = expected.get(relative)
            if item is None:
                raise _AssetFault(component, "file_set")
            digest, final_info = _hash_regular_entry(
                directory_fd,
                name,
                component,
                max_bytes=max(item.size_bytes, 1),
                forbidden_bytes=forbidden_bytes,
            )
            if final_info.st_size != item.size_bytes or digest != item.sha256:
                raise _AssetFault(component, "identity")
            executable = bool(final_info.st_mode & 0o111)
            if item.executable is not None and executable != item.executable:
                raise _AssetFault(component, "mode")
            found.add(relative)
            verified.append(
                {
                    "path": relative,
                    "size_bytes": item.size_bytes,
                    "sha256": digest,
                    "executable": executable,
                }
            )

    walk(tree_fd, "")
    if found != set(expected):
        raise _AssetFault(component, "file_set")
    return sorted(verified, key=lambda entry: str(entry["path"]))


def _verify_file_at(
    root_fd: int,
    relative_path: str,
    expected: _ExpectedFile,
    component: str,
    *,
    max_bytes: int | None = None,
) -> None:
    parent_fd, name = _open_parent(root_fd, relative_path, component)
    try:
        digest, info = _hash_regular_entry(
            parent_fd,
            name,
            component,
            max_bytes=max_bytes or max(expected.size_bytes, 1),
        )
    finally:
        os.close(parent_fd)
    if info.st_size != expected.size_bytes or digest != expected.sha256:
        raise _AssetFault(component, "identity")
    if expected.executable is not None:
        if bool(info.st_mode & 0o111) != expected.executable:
            raise _AssetFault(component, "mode")


def _read_file_at(
    root_fd: int,
    relative_path: str,
    expected: _ExpectedFile | None,
    component: str,
    max_bytes: int,
) -> bytes:
    parent_fd, name = _open_parent(root_fd, relative_path, component)
    try:
        file_fd, before = _open_regular_entry(parent_fd, name, component)
        try:
            if before.st_size > max_bytes:
                raise _AssetFault(component, "bounds")
            content = _read_bounded(file_fd, max_bytes, component)
            after = os.fstat(file_fd)
        finally:
            os.close(file_fd)
    finally:
        os.close(parent_fd)
    _check_stable(before, after, component)
    digest = hashlib.sha256(content).hexdigest()
    if expected is not None and (
        after.st_size != expected.size_bytes or digest != expected.sha256
    ):
        raise _AssetFault(component, "identity")
    return content


def _hash_regular_entry(
    parent_fd: int,
    name: str,
    component: str,
    *,
    max_bytes: int,
    forbidden_bytes: bytes | None = None,
) -> tuple[str, os.stat_result]:
    file_fd, before = _open_regular_entry(parent_fd, name, component)
    try:
        if before.st_size > max_bytes:
            raise _AssetFault(component, "bounds")
        digest = hashlib.sha256()
        total = 0
        scan_tail = b""
        while True:
            chunk = os.read(file_fd, _READ_CHUNK_BYTES)
            if not chunk:
                break
            total += len(chunk)
            if total > max_bytes:
                raise _AssetFault(component, "bounds")
            if forbidden_bytes is not None:
                scan_window = scan_tail + chunk
                if forbidden_bytes in scan_window:
                    raise _AssetFault(component, "reserved_label")
                scan_tail = scan_window[-max(len(forbidden_bytes) - 1, 0) :]
            digest.update(chunk)
        after = os.fstat(file_fd)
    finally:
        os.close(file_fd)
    _check_stable(before, after, component)
    if total != after.st_size:
        raise _AssetFault(component, "race")
    return digest.hexdigest(), after


def _open_root(root: Path) -> int:
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        if not root.is_absolute():
            raise _AssetFault("root", "contract")
        descriptor = os.open("/", flags)
    except OSError as exc:
        raise _AssetFault("root", "unreadable") from exc
    try:
        for component in root.parts[1:]:
            child = _open_single_directory(descriptor, component, "root")
            os.close(descriptor)
            descriptor = child
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def _open_directory_at(root_fd: int, relative_path: str, component: str) -> int:
    descriptor = os.dup(root_fd)
    try:
        for part in PurePosixPath(relative_path).parts:
            child = _open_single_directory(descriptor, part, component)
            os.close(descriptor)
            descriptor = child
        return descriptor
    except Exception:
        os.close(descriptor)
        raise


def _open_single_directory(parent_fd: int, name: str, component: str) -> int:
    flags = os.O_RDONLY | os.O_DIRECTORY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    try:
        info = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
        if not stat.S_ISDIR(info.st_mode) or stat.S_ISLNK(info.st_mode):
            raise _AssetFault(component, "file_type")
        descriptor = os.open(name, flags, dir_fd=parent_fd)
        opened = os.fstat(descriptor)
    except _AssetFault:
        raise
    except OSError as exc:
        raise _AssetFault(component, "missing") from exc
    if (opened.st_dev, opened.st_ino) != (info.st_dev, info.st_ino):
        os.close(descriptor)
        raise _AssetFault(component, "race")
    return descriptor


def _open_parent(root_fd: int, relative_path: str, component: str) -> tuple[int, str]:
    parts = PurePosixPath(relative_path).parts
    descriptor = os.dup(root_fd)
    try:
        for part in parts[:-1]:
            child = _open_single_directory(descriptor, part, component)
            os.close(descriptor)
            descriptor = child
        return descriptor, parts[-1]
    except Exception:
        os.close(descriptor)
        raise


def _open_regular_entry(
    parent_fd: int,
    name: str,
    component: str,
) -> tuple[int, os.stat_result]:
    try:
        before = os.stat(name, dir_fd=parent_fd, follow_symlinks=False)
    except OSError as exc:
        raise _AssetFault(component, "missing") from exc
    if not stat.S_ISREG(before.st_mode) or stat.S_ISLNK(before.st_mode):
        raise _AssetFault(component, "file_type")
    flags = os.O_RDONLY | os.O_CLOEXEC
    if hasattr(os, "O_NOFOLLOW"):
        flags |= os.O_NOFOLLOW
    if hasattr(os, "O_NONBLOCK"):
        flags |= os.O_NONBLOCK
    try:
        descriptor = os.open(name, flags, dir_fd=parent_fd)
        opened = os.fstat(descriptor)
    except OSError as exc:
        raise _AssetFault(component, "race") from exc
    if not stat.S_ISREG(opened.st_mode) or (
        opened.st_dev,
        opened.st_ino,
    ) != (before.st_dev, before.st_ino):
        os.close(descriptor)
        raise _AssetFault(component, "race")
    return descriptor, opened


def _check_stable(
    before: os.stat_result,
    after: os.stat_result,
    component: str,
) -> None:
    before_identity = (
        before.st_dev,
        before.st_ino,
        before.st_size,
        before.st_mtime_ns,
        before.st_ctime_ns,
    )
    after_identity = (
        after.st_dev,
        after.st_ino,
        after.st_size,
        after.st_mtime_ns,
        after.st_ctime_ns,
    )
    if before_identity != after_identity:
        raise _AssetFault(component, "race")


def _read_bounded(file_fd: int, max_bytes: int, component: str) -> bytes:
    content = bytearray()
    while True:
        chunk = os.read(file_fd, min(_READ_CHUNK_BYTES, max_bytes + 1 - len(content)))
        if not chunk:
            return bytes(content)
        content.extend(chunk)
        if len(content) > max_bytes:
            raise _AssetFault(component, "bounds")


def _strict_json_loads(content: bytes) -> Any:
    def reject_duplicates(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        result: dict[str, Any] = {}
        for key, value in pairs:
            if key in result:
                raise ValueError("duplicate JSON object key")
            result[key] = value
        return result

    return json.loads(content, object_pairs_hook=reject_duplicates)


def _canonical_sha256(
    entries: list[dict[str, object]],
    *,
    include_executable: bool = True,
) -> str:
    if all(isinstance(entry.get("path"), str) for entry in entries):
        entries = sorted(entries, key=lambda entry: PurePosixPath(str(entry["path"])))
    if not include_executable:
        entries = [
            {key: value for key, value in entry.items() if key != "executable"}
            for entry in entries
        ]
    payload = json.dumps(
        entries,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _canonical_value_sha256(value: object) -> str:
    payload = json.dumps(
        value,
        ensure_ascii=False,
        separators=(",", ":"),
        sort_keys=True,
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _validate_relative_path(value: object) -> None:
    if not isinstance(value, str) or not value or len(value) > 512:
        raise ValueError("Runtime asset relative path is invalid")
    if "\\" in value or "\x00" in value or any(ord(char) < 32 for char in value):
        raise ValueError("Runtime asset relative path is invalid")
    path = PurePosixPath(value)
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise ValueError("Runtime asset relative path is invalid")
    if path.as_posix() != value:
        raise ValueError("Runtime asset relative path is invalid")


def _join_relative(parent: str, child: str) -> str:
    joined = (PurePosixPath(parent) / PurePosixPath(child)).as_posix()
    _validate_relative_path(joined)
    return joined


def _valid_process_name(value: object) -> bool:
    if not isinstance(value, str) or not value:
        return False
    return all(
        segment
        and segment[0].isalpha()
        and segment[0].isupper()
        and all(char.isupper() or char.isdigit() or char == "_" for char in segment)
        for segment in value.split(":")
    )


def _valid_reserved_label(value: object) -> bool:
    return (
        isinstance(value, str)
        and 1 <= len(value) <= 64
        and value[0] in "abcdefghijklmnopqrstuvwxyz"
        and all(char in "abcdefghijklmnopqrstuvwxyz0123456789_" for char in value)
    )


def _valid_image_coordinate(value: object) -> bool:
    if not isinstance(value, str) or not value or len(value) > 512:
        return False
    if "@" in value or any(char.isspace() or ord(char) < 32 for char in value):
        return False
    final_segment = value.rsplit("/", 1)[-1]
    return "/" in value and ":" in final_segment


def _valid_sha256(value: object) -> bool:
    return (
        isinstance(value, str)
        and len(value) == 64
        and all(char in "0123456789abcdef" for char in value)
    )


def _valid_digest(value: object) -> bool:
    return (
        isinstance(value, str)
        and value.startswith("sha256:")
        and _valid_sha256(value.removeprefix("sha256:"))
    )


def _issue(component: str, reason: str) -> Issue:
    return Issue(
        code=f"BULK_RNASEQ_RUNTIME_ASSET_{reason.upper()}",
        message="A required offline bulk RNA-seq runtime asset is not ready.",
        severity="error",
        source="bulk_rnaseq_runtime_doctor",
        hint="Stage the pinned runtime asset set and run the doctor again.",
        context={"component": component},
    )
