#!/usr/bin/env python3
"""Explicitly stage the pinned nf-core/rnaseq offline runtime asset closure.

The default ``plan`` phase is read-only.  The ``stage`` phase requires separate
network and mutation acknowledgements, writes into a private sibling directory,
admits the complete candidate, atomically publishes it, and admits the final
path again.  Runtime execution remains offline and never invokes this tool.
"""

from __future__ import annotations

import argparse
from collections.abc import Callable, Mapping
import ctypes
from dataclasses import dataclass
import errno
import hashlib
import json
import os
from pathlib import Path, PurePosixPath
import re
import shutil
import stat
import subprocess
import sys
import tarfile
import tempfile
from typing import Any
import zipfile

from jsonschema import Draft202012Validator


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (  # noqa: E402
    NFCORE_RNASEQ_COMMIT,
    NF_SCHEMA_VERSION,
    RuntimeAssetBinding,
    _load_runtime_contract,
    _validate_embedded_contracts,
    verify_runtime_assets,
)


CURL_EXECUTABLE = "/usr/bin/curl"
SOURCE_ARCHIVE_URL = (
    f"https://github.com/nf-core/rnaseq/archive/{NFCORE_RNASEQ_COMMIT}.tar.gz"
)
NEXTFLOW_DISTRIBUTION_URL = (
    "https://github.com/nextflow-io/nextflow/releases/download/"
    "v25.04.3/nextflow-25.04.3-dist"
)
CORRETTO_ARCHIVE_URL = (
    "https://corretto.aws/downloads/resources/21.0.7.6.1/"
    "amazon-corretto-21.0.7.6.1-linux-x64.tar.gz"
)
NF_SCHEMA_ARCHIVE_URL = (
    "https://github.com/nextflow-io/nf-schema/releases/download/"
    f"{NF_SCHEMA_VERSION}/nf-schema-{NF_SCHEMA_VERSION}.zip"
)
NF_SCHEMA_META_URL = (
    "https://github.com/nextflow-io/nf-schema/releases/download/"
    f"{NF_SCHEMA_VERSION}/nf-schema-{NF_SCHEMA_VERSION}-meta.json"
)

_PHASES = frozenset({"plan", "stage", "verify"})
_MAX_SOURCE_ARCHIVE_BYTES = 512 * 1024 * 1024
_MAX_TREE_FILE_BYTES = 256 * 1024 * 1024
_MAX_PLUGIN_TREE_BYTES = 64 * 1024 * 1024
_MAX_RAW_MANIFEST_BYTES = 16 * 1024 * 1024
_MAX_DOCKER_ARCHIVE_BYTES = 100 * 1024 * 1024 * 1024
_MAX_COMMAND_OUTPUT_BYTES = 16 * 1024 * 1024
_COMMAND_TIMEOUT_SECONDS = 7_200
_AT_FDCWD = -100
_RENAME_NOREPLACE = 1
_DIGEST_RE = re.compile(r"^sha256:[a-f0-9]{64}$")
_IMAGE_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._:/-]*:[A-Za-z0-9._-]+$")
_MANIFEST_MEDIA_TYPES = frozenset(
    {
        "application/vnd.docker.distribution.manifest.v2+json",
        "application/vnd.oci.image.manifest.v1+json",
    }
)
_INDEX_MEDIA_TYPES = frozenset(
    {
        "application/vnd.docker.distribution.manifest.list.v2+json",
        "application/vnd.oci.image.index.v1+json",
    }
)


@dataclass(frozen=True)
class StagingContract:
    """Committed coordinates required to construct one runtime asset root."""

    identity: Mapping[str, Any]
    source_manifest: Mapping[str, Any]
    container_inventory: Mapping[str, Any]
    container_schema: Mapping[str, Any]


@dataclass(frozen=True)
class StageOptions:
    """Operator-owned locations and explicit staging authority."""

    asset_root: Path
    docker_executable: Path
    docker_socket: Path
    phase: str = "plan"
    allow_network: bool = False
    allow_mutation: bool = False


@dataclass(frozen=True)
class StagingReport:
    """Path-free staging evidence suitable for logs and CI summaries."""

    phase: str
    container_process_count: int
    unique_image_count: int
    source_tree_sha256: str
    container_lock_sha256: str | None


@dataclass(frozen=True)
class _ImageClosure:
    coordinate: str
    oci_digest: str
    manifest_name: str
    manifest_size_bytes: int
    archive_name: str
    archive_size_bytes: int
    archive_sha256: str


class StagingError(RuntimeError):
    """A stable, path-free failure raised by the staging boundary."""

    def __init__(self, code: str) -> None:
        self.code = code
        super().__init__(code)


CommandRunner = Callable[..., bytes]
AdmissionVerifier = Callable[[RuntimeAssetBinding], None]


def load_committed_contract() -> StagingContract:
    """Load and authenticate the embedded runtime staging coordinates."""

    try:
        contract = _load_runtime_contract()
        _validate_embedded_contracts(
            contract.identity,
            contract.source_manifest,
            contract.container_inventory,
            contract.container_process_audit,
        )
    except (OSError, TypeError, ValueError, json.JSONDecodeError) as exc:
        raise StagingError("committed_contract_invalid") from exc
    return StagingContract(
        identity=contract.identity,
        source_manifest=contract.source_manifest,
        container_inventory=contract.container_inventory,
        container_schema=contract.container_schema,
    )


def stage_runtime_assets(
    options: StageOptions,
    *,
    contract: StagingContract | None = None,
    command_runner: CommandRunner | None = None,
    admission_verifier: AdmissionVerifier | None = None,
) -> StagingReport:
    """Plan, stage, or verify the exact pinned runtime asset closure."""

    if not isinstance(options, StageOptions):
        raise StagingError("options_invalid")
    selected = contract if contract is not None else load_committed_contract()
    _validate_contract(selected)
    _validate_options(options, require_existing_root=options.phase == "verify")
    inventory_entries = _inventory_entries(selected)
    unique_images = {str(entry["image_coordinate"]) for entry in inventory_entries}
    source_tree_sha256 = _required_digest(
        selected.source_manifest.get("tree_sha256"),
        "committed_contract_invalid",
        prefixed=False,
    )

    if options.phase == "plan":
        if options.allow_network or options.allow_mutation:
            raise StagingError("plan_authority_invalid")
        return StagingReport(
            phase="plan",
            container_process_count=len(inventory_entries),
            unique_image_count=len(unique_images),
            source_tree_sha256=source_tree_sha256,
            container_lock_sha256=None,
        )
    if options.phase == "verify":
        if options.allow_network or options.allow_mutation:
            raise StagingError("verify_authority_invalid")
        _verify_docker_endpoint(options)
        _verify_candidate_files(options.asset_root, selected)
        verifier = admission_verifier or _verify_production_admission
        verifier(_binding(options.asset_root, options, selected.identity))
        lock_digest = _sha256_file(
            options.asset_root / "containers/availability-lock.json",
            2 * 1024 * 1024,
        )
        return StagingReport(
            phase="verify",
            container_process_count=len(inventory_entries),
            unique_image_count=len(unique_images),
            source_tree_sha256=source_tree_sha256,
            container_lock_sha256=lock_digest,
        )
    if not options.allow_network or not options.allow_mutation:
        raise StagingError("stage_authority_required")

    _verify_docker_endpoint(options)
    if options.asset_root.exists() or options.asset_root.is_symlink():
        raise StagingError("asset_root_exists")

    work_path: Path | None = None
    published_identity: tuple[int, int] | None = None
    try:
        work_path = Path(
            tempfile.mkdtemp(
                prefix=f".{options.asset_root.name}.stage-",
                dir=options.asset_root.parent,
            )
        )
        work_path.chmod(0o700)
        candidate = work_path / "assets"
        candidate.mkdir(mode=0o700)
        if command_runner is None:
            command_home = work_path / "command-home"
            command_home.mkdir(mode=0o700)
            (command_home / "docker").mkdir(mode=0o700)

            def runner(
                argv: tuple[str, ...],
                *,
                capture_stdout: bool,
            ) -> bytes:
                return _run_external(
                    argv,
                    capture_stdout=capture_stdout,
                    command_home=command_home,
                )

        else:
            runner = command_runner

        _stage_static_assets(candidate, work_path, selected, runner)
        closures = _stage_container_assets(candidate, options, selected, runner)
        lock_content = _availability_lock_bytes(selected, closures)
        _atomic_write(
            candidate / "containers/availability-lock.json",
            lock_content,
            mode=0o600,
        )

        verifier = admission_verifier or _verify_production_admission
        _verify_candidate_files(candidate, selected)
        verifier(_binding(candidate, options, selected.identity))

        _publish_candidate(candidate, options.asset_root)
        published = os.lstat(options.asset_root)
        published_identity = (published.st_dev, published.st_ino)

        _verify_candidate_files(options.asset_root, selected)
        verifier(_binding(options.asset_root, options, selected.identity))
        return StagingReport(
            phase="stage",
            container_process_count=len(inventory_entries),
            unique_image_count=len(unique_images),
            source_tree_sha256=source_tree_sha256,
            container_lock_sha256=hashlib.sha256(lock_content).hexdigest(),
        )
    except StagingError:
        if published_identity is not None:
            _remove_published_root(options.asset_root, published_identity)
        raise
    except (OSError, ValueError, json.JSONDecodeError, tarfile.TarError) as exc:
        if published_identity is not None:
            _remove_published_root(options.asset_root, published_identity)
        raise StagingError("staging_failed") from exc
    finally:
        if work_path is not None:
            shutil.rmtree(work_path, ignore_errors=True)


def _validate_contract(contract: StagingContract) -> None:
    if not isinstance(contract, StagingContract):
        raise StagingError("committed_contract_invalid")
    identity = contract.identity
    source = identity.get("source") if isinstance(identity, Mapping) else None
    nextflow = identity.get("nextflow") if isinstance(identity, Mapping) else None
    jdk = identity.get("jdk") if isinstance(identity, Mapping) else None
    plugins = identity.get("plugins") if isinstance(identity, Mapping) else None
    if (
        not isinstance(source, Mapping)
        or source.get("project") != "nf-core/rnaseq"
        or source.get("release") != "3.26.0"
        or source.get("commit") != NFCORE_RNASEQ_COMMIT
        or not isinstance(nextflow, Mapping)
        or nextflow.get("version") != "25.04.3"
        or not isinstance(jdk, Mapping)
        or jdk.get("distribution") != "amazon-corretto"
        or jdk.get("version") != "21.0.7.6.1"
        or jdk.get("runtime_version") != "21.0.7+6-LTS"
        or not isinstance(plugins, list)
        or len(plugins) != 1
        or not isinstance(plugins[0], Mapping)
        or plugins[0].get("id") != "nf-schema"
        or plugins[0].get("version") != NF_SCHEMA_VERSION
    ):
        raise StagingError("committed_contract_invalid")
    for value in (
        source.get("default_relative_path"),
        nextflow.get("default_relative_path"),
        jdk.get("default_archive_relative_path"),
        jdk.get("default_tree_relative_path"),
        plugins[0].get("default_archive_relative_path"),
        plugins[0].get("default_meta_relative_path"),
        plugins[0].get("default_tree_relative_path"),
    ):
        _relative_path(value, "committed_contract_invalid")
    _required_digest(
        contract.source_manifest.get("source_archive_sha256"),
        "committed_contract_invalid",
        prefixed=False,
    )
    _required_digest(
        contract.source_manifest.get("tree_sha256"),
        "committed_contract_invalid",
        prefixed=False,
    )
    _inventory_entries(contract)


def _inventory_entries(contract: StagingContract) -> tuple[Mapping[str, Any], ...]:
    inventory = contract.container_inventory
    entries = inventory.get("entries") if isinstance(inventory, Mapping) else None
    if not isinstance(entries, list) or not entries:
        raise StagingError("committed_contract_invalid")
    processes: set[str] = set()
    normalized: list[Mapping[str, Any]] = []
    for entry in entries:
        if not isinstance(entry, Mapping):
            raise StagingError("committed_contract_invalid")
        process = entry.get("process")
        coordinate = entry.get("image_coordinate")
        if (
            not isinstance(process, str)
            or not process
            or process in processes
            or not isinstance(coordinate, str)
            or _IMAGE_RE.fullmatch(coordinate) is None
        ):
            raise StagingError("committed_contract_invalid")
        processes.add(process)
        normalized.append(entry)
    if normalized != sorted(normalized, key=lambda value: str(value["process"])):
        raise StagingError("committed_contract_invalid")
    return tuple(normalized)


def _validate_options(options: StageOptions, *, require_existing_root: bool) -> None:
    if options.phase not in _PHASES:
        raise StagingError("phase_invalid")
    for name in ("allow_network", "allow_mutation"):
        if not isinstance(getattr(options, name), bool):
            raise StagingError("options_invalid")
    _absolute_path(options.asset_root, "asset_root_invalid", allow_missing=True)
    _absolute_path(
        options.docker_executable,
        "docker_endpoint_invalid",
        allow_missing=True,
    )
    _absolute_path(
        options.docker_socket,
        "docker_endpoint_invalid",
        allow_missing=True,
    )
    if not options.asset_root.parent.is_dir():
        raise StagingError("asset_root_invalid")
    _reject_symlink_ancestors(options.asset_root.parent, "asset_root_invalid")
    if require_existing_root:
        try:
            info = os.lstat(options.asset_root)
        except OSError as exc:
            raise StagingError("asset_root_missing") from exc
        if not stat.S_ISDIR(info.st_mode):
            raise StagingError("asset_root_invalid")


def _verify_docker_endpoint(options: StageOptions) -> None:
    try:
        executable = os.lstat(options.docker_executable)
        endpoint = os.lstat(options.docker_socket)
    except OSError as exc:
        raise StagingError("docker_endpoint_invalid") from exc
    if (
        options.docker_executable.name != "docker"
        or not stat.S_ISREG(executable.st_mode)
        or executable.st_nlink != 1
        or executable.st_mode & 0o111 == 0
        or executable.st_mode & 0o022
        or executable.st_uid not in {0, os.geteuid()}
        or not stat.S_ISSOCK(endpoint.st_mode)
        or endpoint.st_nlink != 1
    ):
        raise StagingError("docker_endpoint_invalid")


def _stage_static_assets(
    candidate: Path,
    work_path: Path,
    contract: StagingContract,
    runner: CommandRunner,
) -> None:
    identity = contract.identity
    source = _mapping(identity.get("source"))
    nextflow = _mapping(identity.get("nextflow"))
    jdk = _mapping(identity.get("jdk"))
    plugin = _mapping(_sequence(identity.get("plugins"))[0])

    downloads = work_path / "downloads"
    downloads.mkdir(mode=0o700)
    source_archive = downloads / "source.tar.gz"
    _download(
        SOURCE_ARCHIVE_URL,
        source_archive,
        maximum_bytes=_MAX_SOURCE_ARCHIVE_BYTES,
        runner=runner,
    )
    _verify_file(
        source_archive,
        expected_size=None,
        expected_sha256=_required_digest(
            contract.source_manifest.get("source_archive_sha256"),
            "committed_contract_invalid",
            prefixed=False,
        ),
        executable=None,
        maximum_bytes=_MAX_SOURCE_ARCHIVE_BYTES,
    )
    _extract_source_archive(
        source_archive,
        candidate / str(source["default_relative_path"]),
        contract.source_manifest,
    )

    _download_fixed_file(
        NEXTFLOW_DISTRIBUTION_URL,
        candidate / str(nextflow["default_relative_path"]),
        expected_size=_positive_int(nextflow.get("size_bytes")),
        expected_sha256=_required_digest(
            nextflow.get("sha256"),
            "committed_contract_invalid",
            prefixed=False,
        ),
        mode=0o755,
        runner=runner,
    )
    jdk_archive = candidate / str(jdk["default_archive_relative_path"])
    _download_fixed_file(
        CORRETTO_ARCHIVE_URL,
        jdk_archive,
        expected_size=_positive_int(jdk.get("archive_size_bytes")),
        expected_sha256=_required_digest(
            jdk.get("archive_sha256"),
            "committed_contract_invalid",
            prefixed=False,
        ),
        mode=0o644,
        runner=runner,
    )
    _extract_jdk_archive(
        jdk_archive,
        candidate / str(jdk["default_tree_relative_path"]),
        jdk,
    )

    plugin_archive = candidate / str(plugin["default_archive_relative_path"])
    plugin_meta = candidate / str(plugin["default_meta_relative_path"])
    _download_fixed_file(
        NF_SCHEMA_ARCHIVE_URL,
        plugin_archive,
        expected_size=_positive_int(plugin.get("archive_size_bytes")),
        expected_sha256=_required_digest(
            plugin.get("archive_sha256"),
            "committed_contract_invalid",
            prefixed=False,
        ),
        mode=0o644,
        runner=runner,
    )
    _download_fixed_file(
        NF_SCHEMA_META_URL,
        plugin_meta,
        expected_size=_positive_int(plugin.get("meta_size_bytes")),
        expected_sha256=_required_digest(
            plugin.get("meta_sha256"),
            "committed_contract_invalid",
            prefixed=False,
        ),
        mode=0o644,
        runner=runner,
    )
    _extract_plugin_archive(
        plugin_archive,
        plugin_meta,
        candidate / str(plugin["default_tree_relative_path"]),
        plugin,
    )


def _stage_container_assets(
    candidate: Path,
    options: StageOptions,
    contract: StagingContract,
    runner: CommandRunner,
) -> dict[str, _ImageClosure]:
    assets = candidate / "containers/assets"
    assets.mkdir(parents=True, mode=0o700)
    coordinates = sorted(
        {str(entry["image_coordinate"]) for entry in _inventory_entries(contract)}
    )
    host = f"unix://{options.docker_socket}"
    closures: dict[str, _ImageClosure] = {}
    by_digest: dict[str, _ImageClosure] = {}
    for coordinate in coordinates:
        index_or_manifest = _run_capture(
            runner,
            (
                str(options.docker_executable),
                "--host",
                host,
                "buildx",
                "imagetools",
                "inspect",
                "--raw",
                coordinate,
            ),
        )
        digest, expected_size = _select_linux_amd64_manifest(index_or_manifest)
        digest_reference = f"{coordinate}@{digest}"
        selected_output = _run_capture(
            runner,
            (
                str(options.docker_executable),
                "--host",
                host,
                "buildx",
                "imagetools",
                "inspect",
                "--raw",
                digest_reference,
            ),
        )
        selected_manifest, config_digest = _validated_selected_manifest(
            selected_output,
            expected_digest=digest,
            expected_size=expected_size,
        )

        _run_no_capture(
            runner,
            (
                str(options.docker_executable),
                "--host",
                host,
                "image",
                "pull",
                "--platform",
                "linux/amd64",
                digest_reference,
            ),
        )
        inspected = _run_capture(
            runner,
            (
                str(options.docker_executable),
                "--host",
                host,
                "image",
                "inspect",
                digest_reference,
            ),
        )
        _verify_pulled_image(inspected, config_digest=config_digest)

        manifest_name = f"{digest.removeprefix('sha256:')}.manifest.json"
        archive_name = f"{digest.removeprefix('sha256:')}.docker.tar"
        manifest_path = assets / manifest_name
        archive_path = assets / archive_name
        previous = by_digest.get(digest)
        if previous is None:
            _atomic_write(manifest_path, selected_manifest, mode=0o600)
            temporary_archive = assets / f".{archive_name}.part"
            _run_no_capture(
                runner,
                (
                    str(options.docker_executable),
                    "--host",
                    host,
                    "image",
                    "save",
                    "--output",
                    str(temporary_archive),
                    digest_reference,
                ),
            )
            archive_size = _regular_file_size(
                temporary_archive,
                maximum_bytes=_MAX_DOCKER_ARCHIVE_BYTES,
            )
            if archive_size <= 0:
                raise StagingError("container_archive_invalid")
            archive_sha256 = _sha256_file(
                temporary_archive,
                _MAX_DOCKER_ARCHIVE_BYTES,
            )
            temporary_archive.chmod(0o600)
            os.replace(temporary_archive, archive_path)
            closure = _ImageClosure(
                coordinate=coordinate,
                oci_digest=digest,
                manifest_name=manifest_name,
                manifest_size_bytes=len(selected_manifest),
                archive_name=archive_name,
                archive_size_bytes=archive_size,
                archive_sha256=archive_sha256,
            )
            by_digest[digest] = closure
        else:
            if manifest_path.read_bytes() != selected_manifest:
                raise StagingError("container_manifest_invalid")
            closure = _ImageClosure(
                coordinate=coordinate,
                oci_digest=previous.oci_digest,
                manifest_name=previous.manifest_name,
                manifest_size_bytes=previous.manifest_size_bytes,
                archive_name=previous.archive_name,
                archive_size_bytes=previous.archive_size_bytes,
                archive_sha256=previous.archive_sha256,
            )
        closures[coordinate] = closure
    return closures


def _availability_lock_bytes(
    contract: StagingContract,
    closures: Mapping[str, _ImageClosure],
) -> bytes:
    identity = contract.identity
    source = _mapping(identity.get("source"))
    inventory = contract.container_inventory
    inventory_digest = _required_digest(
        inventory.get("entries_sha256"),
        "committed_contract_invalid",
        prefixed=False,
    )
    entries = []
    for item in _inventory_entries(contract):
        coordinate = str(item["image_coordinate"])
        closure = closures.get(coordinate)
        if closure is None:
            raise StagingError("container_closure_incomplete")
        entries.append(
            {
                "process": item["process"],
                "image": f"{coordinate}@{closure.oci_digest}",
                "oci_digest": closure.oci_digest,
                "distribution_manifest_asset": closure.manifest_name,
                "distribution_manifest_size_bytes": closure.manifest_size_bytes,
                "local_asset": closure.archive_name,
                "size_bytes": closure.archive_size_bytes,
                "sha256": closure.archive_sha256,
            }
        )
    value = {
        "schema_version": "1.0.0",
        "pipeline": {
            "project": "nf-core/rnaseq",
            "release": "3.26.0",
            "commit": NFCORE_RNASEQ_COMMIT,
            "source_tree_sha256": source["tree_sha256"],
            "container_inventory_sha256": inventory_digest,
        },
        "asset_format": "docker-archive+distribution-manifest",
        "network_policy": {
            "container_pull": False,
            "wave": False,
            "tower": False,
        },
        "execution_policy": {
            "nxf_offline": True,
            "docker_pull": "never",
            "docker_network": "none",
        },
        "entries": entries,
    }
    try:
        errors = tuple(
            Draft202012Validator(contract.container_schema).iter_errors(value)
        )
    except (TypeError, ValueError) as exc:
        raise StagingError("availability_lock_invalid") from exc
    if errors:
        raise StagingError("availability_lock_invalid")
    return _json_bytes(value)


def _verify_candidate_files(root: Path, contract: StagingContract) -> None:
    identity = contract.identity
    source = _mapping(identity.get("source"))
    nextflow = _mapping(identity.get("nextflow"))
    jdk = _mapping(identity.get("jdk"))
    plugin = _mapping(_sequence(identity.get("plugins"))[0])
    _verify_source_tree(
        root / str(source["default_relative_path"]),
        contract.source_manifest,
    )
    _verify_file(
        root / str(nextflow["default_relative_path"]),
        expected_size=_positive_int(nextflow.get("size_bytes")),
        expected_sha256=_required_digest(
            nextflow.get("sha256"),
            "committed_contract_invalid",
            prefixed=False,
        ),
        executable=True,
        maximum_bytes=_positive_int(nextflow.get("size_bytes")),
    )
    _verify_file(
        root / str(jdk["default_archive_relative_path"]),
        expected_size=_positive_int(jdk.get("archive_size_bytes")),
        expected_sha256=_required_digest(
            jdk.get("archive_sha256"),
            "committed_contract_invalid",
            prefixed=False,
        ),
        executable=False,
        maximum_bytes=_positive_int(jdk.get("archive_size_bytes")),
    )
    _verify_jdk_tree(root / str(jdk["default_tree_relative_path"]), jdk)
    _verify_file(
        root / str(plugin["default_archive_relative_path"]),
        expected_size=_positive_int(plugin.get("archive_size_bytes")),
        expected_sha256=_required_digest(
            plugin.get("archive_sha256"),
            "committed_contract_invalid",
            prefixed=False,
        ),
        executable=False,
        maximum_bytes=_positive_int(plugin.get("archive_size_bytes")),
    )
    _verify_file(
        root / str(plugin["default_meta_relative_path"]),
        expected_size=_positive_int(plugin.get("meta_size_bytes")),
        expected_sha256=_required_digest(
            plugin.get("meta_sha256"),
            "committed_contract_invalid",
            prefixed=False,
        ),
        executable=False,
        maximum_bytes=_positive_int(plugin.get("meta_size_bytes")),
    )
    _verify_plugin_tree(
        root / str(plugin["default_archive_relative_path"]),
        root / str(plugin["default_meta_relative_path"]),
        root / str(plugin["default_tree_relative_path"]),
        plugin,
    )
    _verify_availability_lock(root, contract)


def _verify_availability_lock(root: Path, contract: StagingContract) -> None:
    lock_path = root / "containers/availability-lock.json"
    content = _read_regular_file(lock_path, 2 * 1024 * 1024)
    try:
        value = _strict_json(content)
        errors = tuple(
            Draft202012Validator(contract.container_schema).iter_errors(value)
        )
    except (TypeError, ValueError, json.JSONDecodeError) as exc:
        raise StagingError("availability_lock_invalid") from exc
    if errors or not isinstance(value, Mapping):
        raise StagingError("availability_lock_invalid")
    entries = value.get("entries")
    if not isinstance(entries, list):
        raise StagingError("availability_lock_invalid")
    inventory_by_process = {
        str(entry["process"]): entry for entry in _inventory_entries(contract)
    }
    if {entry.get("process") for entry in entries if isinstance(entry, Mapping)} != set(
        inventory_by_process
    ):
        raise StagingError("availability_lock_invalid")
    assets = root / "containers/assets"
    expected_asset_names: set[str] = set()
    coordinate_closures: dict[str, tuple[object, ...]] = {}
    for entry in entries:
        if not isinstance(entry, Mapping):
            raise StagingError("availability_lock_invalid")
        process = str(entry["process"])
        coordinate = str(inventory_by_process[process]["image_coordinate"])
        digest = _required_digest(
            entry.get("oci_digest"),
            "availability_lock_invalid",
            prefixed=True,
        )
        if entry.get("image") != f"{coordinate}@{digest}":
            raise StagingError("availability_lock_invalid")
        manifest_name = _single_component(
            entry.get("distribution_manifest_asset"),
            "availability_lock_invalid",
        )
        archive_name = _single_component(
            entry.get("local_asset"),
            "availability_lock_invalid",
        )
        expected_asset_names.update({manifest_name, archive_name})
        closure = (
            digest,
            manifest_name,
            entry.get("distribution_manifest_size_bytes"),
            archive_name,
            entry.get("size_bytes"),
            entry.get("sha256"),
        )
        if (
            coordinate in coordinate_closures
            and coordinate_closures[coordinate] != closure
        ):
            raise StagingError("availability_lock_invalid")
        coordinate_closures[coordinate] = closure
        manifest = _read_regular_file(
            assets / manifest_name,
            _MAX_RAW_MANIFEST_BYTES,
        )
        normalized, _config = _validated_selected_manifest(
            manifest,
            expected_digest=digest,
            expected_size=_positive_int(entry.get("distribution_manifest_size_bytes")),
        )
        if normalized != manifest:
            raise StagingError("container_manifest_invalid")
        archive_size = _positive_int(entry.get("size_bytes"))
        archive_sha = _required_digest(
            entry.get("sha256"),
            "availability_lock_invalid",
            prefixed=False,
        )
        _verify_file(
            assets / archive_name,
            expected_size=archive_size,
            expected_sha256=archive_sha,
            executable=False,
            maximum_bytes=_MAX_DOCKER_ARCHIVE_BYTES,
        )
    try:
        actual_names = {entry.name for entry in os.scandir(assets)}
    except OSError as exc:
        raise StagingError("availability_lock_invalid") from exc
    if actual_names != expected_asset_names:
        raise StagingError("availability_lock_invalid")


def _extract_source_archive(
    archive: Path,
    destination: Path,
    manifest: Mapping[str, Any],
) -> None:
    files = manifest.get("files")
    if not isinstance(files, list):
        raise StagingError("source_manifest_invalid")
    expected: dict[str, Mapping[str, Any]] = {}
    for entry in files:
        if not isinstance(entry, Mapping):
            raise StagingError("source_manifest_invalid")
        path = _relative_path(entry.get("path"), "source_manifest_invalid")
        if path in expected:
            raise StagingError("source_manifest_invalid")
        expected[path] = entry
    destination.mkdir(parents=True, mode=0o700)
    found: set[str] = set()
    prefix: str | None = None
    try:
        with tarfile.open(archive, mode="r:gz") as package:
            for member in package:
                parts = _archive_parts(member.name, "source_archive_invalid")
                if prefix is None:
                    prefix = parts[0]
                if parts[0] != prefix:
                    raise StagingError("source_archive_invalid")
                if member.isdir():
                    continue
                if not member.isreg() or len(parts) < 2:
                    raise StagingError("source_archive_invalid")
                relative = PurePosixPath(*parts[1:]).as_posix()
                entry = expected.get(relative)
                if entry is None or relative in found:
                    raise StagingError("source_archive_invalid")
                content = _read_tar_member(
                    package,
                    member,
                    _positive_nonnegative_int(entry.get("size_bytes")),
                    "source_archive_invalid",
                )
                expected_digest = _required_digest(
                    entry.get("sha256"),
                    "source_manifest_invalid",
                    prefixed=False,
                )
                if _sha256_bytes(content) != expected_digest:
                    raise StagingError("asset_identity_invalid")
                mode = 0o755 if entry.get("executable") is True else 0o644
                _atomic_write(destination / relative, content, mode=mode)
                found.add(relative)
    except (OSError, tarfile.TarError) as exc:
        raise StagingError("source_archive_invalid") from exc
    if found != set(expected):
        raise StagingError("source_archive_invalid")
    _verify_source_tree(destination, manifest)


def _extract_jdk_archive(
    archive: Path,
    destination: Path,
    identity: Mapping[str, Any],
) -> None:
    destination.mkdir(parents=True, mode=0o700)
    prefix: str | None = None
    found: set[str] = set()
    maximum_total = _positive_int(identity.get("tree_size_bytes"))
    total = 0
    try:
        with tarfile.open(archive, mode="r:gz") as package:
            for member in package:
                parts = _archive_parts(member.name, "jdk_archive_invalid")
                if prefix is None:
                    prefix = parts[0]
                if parts[0] != prefix:
                    raise StagingError("jdk_archive_invalid")
                if member.isdir():
                    continue
                if not member.isreg() or len(parts) < 2:
                    raise StagingError("jdk_archive_invalid")
                relative = PurePosixPath(*parts[1:]).as_posix()
                if relative in found or member.size > _MAX_TREE_FILE_BYTES:
                    raise StagingError("jdk_archive_invalid")
                total += member.size
                if total > maximum_total:
                    raise StagingError("jdk_archive_invalid")
                content = _read_tar_member(
                    package,
                    member,
                    member.size,
                    "jdk_archive_invalid",
                )
                mode = 0o755 if member.mode & 0o111 else 0o644
                _atomic_write(destination / relative, content, mode=mode)
                found.add(relative)
    except (OSError, tarfile.TarError) as exc:
        raise StagingError("jdk_archive_invalid") from exc
    _verify_jdk_tree(destination, identity)


def _extract_plugin_archive(
    archive: Path,
    meta: Path,
    destination: Path,
    identity: Mapping[str, Any],
) -> None:
    destination.mkdir(parents=True, mode=0o700)
    seen: set[str] = set()
    total = 0
    try:
        with zipfile.ZipFile(archive) as package:
            for entry in package.infolist():
                if entry.is_dir():
                    continue
                relative = _relative_path(entry.filename, "plugin_archive_invalid")
                file_type = (entry.external_attr >> 16) & 0o170000
                if (
                    relative in seen
                    or entry.flag_bits & 0x1
                    or file_type == stat.S_IFLNK
                ):
                    raise StagingError("plugin_archive_invalid")
                content = package.read(entry)
                total += len(content)
                if total > _MAX_PLUGIN_TREE_BYTES:
                    raise StagingError("plugin_archive_invalid")
                _atomic_write(destination / relative, content, mode=0o644)
                seen.add(relative)
    except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
        raise StagingError("plugin_archive_invalid") from exc
    meta_name = f"nf-schema-{NF_SCHEMA_VERSION}-meta.json"
    if meta_name in seen:
        raise StagingError("plugin_archive_invalid")
    _atomic_write(
        destination / meta_name,
        _read_regular_file(meta, 2 * 1024 * 1024),
        mode=0o644,
    )
    _verify_plugin_tree(archive, meta, destination, identity)


def _verify_source_tree(tree: Path, manifest: Mapping[str, Any]) -> None:
    files = manifest.get("files")
    if not isinstance(files, list):
        raise StagingError("source_manifest_invalid")
    expected = {
        str(entry.get("path")): entry for entry in files if isinstance(entry, Mapping)
    }
    if len(expected) != len(files):
        raise StagingError("source_manifest_invalid")
    observed = _observed_tree(
        tree,
        maximum_files=len(expected),
        maximum_total_bytes=sum(
            _positive_nonnegative_int(entry.get("size_bytes"))
            for entry in expected.values()
        ),
    )
    if {str(entry["path"]) for entry in observed} != set(expected):
        raise StagingError("asset_identity_invalid")
    for item in observed:
        configured = expected[str(item["path"])]
        if (
            item["size_bytes"] != configured.get("size_bytes")
            or item["sha256"] != configured.get("sha256")
            or item["executable"] is not configured.get("executable")
        ):
            raise StagingError("asset_identity_invalid")
    expected_tree_digest = _required_digest(
        manifest.get("tree_sha256"),
        "source_manifest_invalid",
        prefixed=False,
    )
    if _canonical_tree_digest(observed) != expected_tree_digest:
        raise StagingError("asset_identity_invalid")


def _verify_jdk_tree(tree: Path, identity: Mapping[str, Any]) -> None:
    maximum_files = _positive_int(identity.get("tree_file_count"))
    maximum_total = _positive_int(identity.get("tree_size_bytes"))
    observed = _observed_tree(
        tree,
        maximum_files=maximum_files,
        maximum_total_bytes=maximum_total,
    )
    if (
        len(observed) != maximum_files
        or sum(int(entry["size_bytes"]) for entry in observed) != maximum_total
        or _canonical_tree_digest(observed)
        != _required_digest(
            identity.get("tree_sha256"),
            "committed_contract_invalid",
            prefixed=False,
        )
    ):
        raise StagingError("asset_identity_invalid")
    by_path = {str(entry["path"]): entry for entry in observed}
    for path_key, size_key, digest_key, executable in (
        ("java_relative_path", "java_size_bytes", "java_sha256", True),
        ("license_file", None, "license_file_sha256", False),
        (
            "assembly_exception_file",
            None,
            "assembly_exception_file_sha256",
            False,
        ),
    ):
        item = by_path.get(str(identity.get(path_key)))
        if (
            item is None
            or item["sha256"] != identity.get(digest_key)
            or size_key is not None
            and item["size_bytes"] != identity.get(size_key)
            or executable
            and item["executable"] is not True
        ):
            raise StagingError("asset_identity_invalid")


def _verify_plugin_tree(
    archive: Path,
    meta: Path,
    tree: Path,
    identity: Mapping[str, Any],
) -> None:
    expected: dict[str, dict[str, object]] = {}
    try:
        with zipfile.ZipFile(archive) as package:
            for entry in package.infolist():
                if entry.is_dir():
                    continue
                relative = _relative_path(entry.filename, "plugin_archive_invalid")
                if relative in expected:
                    raise StagingError("plugin_archive_invalid")
                content = package.read(entry)
                expected[relative] = {
                    "path": relative,
                    "size_bytes": len(content),
                    "sha256": _sha256_bytes(content),
                }
    except (OSError, RuntimeError, zipfile.BadZipFile) as exc:
        raise StagingError("plugin_archive_invalid") from exc
    meta_content = _read_regular_file(meta, 2 * 1024 * 1024)
    meta_name = f"nf-schema-{NF_SCHEMA_VERSION}-meta.json"
    expected[meta_name] = {
        "path": meta_name,
        "size_bytes": len(meta_content),
        "sha256": _sha256_bytes(meta_content),
    }
    observed = _observed_tree(
        tree,
        maximum_files=_positive_int(identity.get("tree_file_count")),
        maximum_total_bytes=_MAX_PLUGIN_TREE_BYTES,
    )
    without_modes = [
        {
            "path": entry["path"],
            "size_bytes": entry["size_bytes"],
            "sha256": entry["sha256"],
        }
        for entry in observed
    ]
    if (
        without_modes != [expected[path] for path in sorted(expected)]
        or len(observed) != _positive_int(identity.get("tree_file_count"))
        or _canonical_tree_digest(without_modes)
        != _required_digest(
            identity.get("tree_sha256"),
            "committed_contract_invalid",
            prefixed=False,
        )
    ):
        raise StagingError("asset_identity_invalid")


def _select_linux_amd64_manifest(content: bytes) -> tuple[str, int]:
    if len(content) == 0 or len(content) > _MAX_RAW_MANIFEST_BYTES:
        raise StagingError("container_manifest_invalid")
    try:
        value = _strict_json(content)
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError) as exc:
        raise StagingError("container_manifest_invalid") from exc
    if not isinstance(value, Mapping) or value.get("schemaVersion") != 2:
        raise StagingError("container_manifest_invalid")
    manifests = value.get("manifests")
    if manifests is None:
        _validate_manifest_shape(value)
        return f"sha256:{_sha256_bytes(content)}", len(content)
    if value.get("mediaType") not in _INDEX_MEDIA_TYPES or not isinstance(
        manifests, list
    ):
        raise StagingError("container_manifest_invalid")
    selected: list[tuple[str, int]] = []
    for descriptor in manifests:
        if not isinstance(descriptor, Mapping):
            raise StagingError("container_manifest_invalid")
        platform = descriptor.get("platform")
        if not isinstance(platform, Mapping):
            continue
        if (
            platform.get("os") == "linux"
            and platform.get("architecture") == "amd64"
            and platform.get("variant") in {None, ""}
        ):
            digest = _required_digest(
                descriptor.get("digest"),
                "container_manifest_invalid",
                prefixed=True,
            )
            selected.append((digest, _positive_int(descriptor.get("size"))))
    if len(selected) != 1:
        raise StagingError("container_platform_invalid")
    return selected[0]


def _validated_selected_manifest(
    content: bytes,
    *,
    expected_digest: str,
    expected_size: int,
) -> tuple[bytes, str]:
    candidates = [content]
    if content.endswith(b"\n"):
        candidates.append(content[:-1])
    normalized = next(
        (
            candidate
            for candidate in candidates
            if f"sha256:{_sha256_bytes(candidate)}" == expected_digest
            and len(candidate) == expected_size
        ),
        None,
    )
    if normalized is None or len(normalized) > _MAX_RAW_MANIFEST_BYTES:
        raise StagingError("container_manifest_invalid")
    try:
        value = _strict_json(normalized)
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError) as exc:
        raise StagingError("container_manifest_invalid") from exc
    _validate_manifest_shape(value)
    config = _mapping(value.get("config"))
    return normalized, _required_digest(
        config.get("digest"),
        "container_manifest_invalid",
        prefixed=True,
    )


def _validate_manifest_shape(value: object) -> None:
    if not isinstance(value, Mapping) or value.get("schemaVersion") != 2:
        raise StagingError("container_manifest_invalid")
    if value.get("mediaType") not in _MANIFEST_MEDIA_TYPES:
        raise StagingError("container_manifest_invalid")
    config = value.get("config")
    layers = value.get("layers")
    if not isinstance(config, Mapping) or not isinstance(layers, list):
        raise StagingError("container_manifest_invalid")
    _required_digest(
        config.get("digest"),
        "container_manifest_invalid",
        prefixed=True,
    )
    _positive_int(config.get("size"))
    for layer in layers:
        if not isinstance(layer, Mapping):
            raise StagingError("container_manifest_invalid")
        _required_digest(
            layer.get("digest"),
            "container_manifest_invalid",
            prefixed=True,
        )
        size = layer.get("size")
        if isinstance(size, bool) or not isinstance(size, int) or size < 0:
            raise StagingError("container_manifest_invalid")


def _verify_pulled_image(content: bytes, *, config_digest: str) -> None:
    try:
        value = _strict_json(content)
    except (UnicodeDecodeError, ValueError, json.JSONDecodeError) as exc:
        raise StagingError("docker_image_invalid") from exc
    if (
        not isinstance(value, list)
        or len(value) != 1
        or not isinstance(value[0], Mapping)
    ):
        raise StagingError("docker_image_invalid")
    image = value[0]
    rootfs = image.get("RootFS")
    if (
        image.get("Id") != config_digest
        or image.get("Os") != "linux"
        or image.get("Architecture") != "amd64"
        or not isinstance(rootfs, Mapping)
        or rootfs.get("Type") != "layers"
        or not isinstance(rootfs.get("Layers"), list)
        or not all(_DIGEST_RE.fullmatch(item) for item in rootfs["Layers"])
    ):
        raise StagingError("docker_image_invalid")


def _download_fixed_file(
    url: str,
    destination: Path,
    *,
    expected_size: int,
    expected_sha256: str,
    mode: int,
    runner: CommandRunner,
) -> None:
    destination.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    temporary = destination.with_name(f".{destination.name}.part")
    _download(url, temporary, maximum_bytes=expected_size, runner=runner)
    _verify_file(
        temporary,
        expected_size=expected_size,
        expected_sha256=expected_sha256,
        executable=None,
        maximum_bytes=expected_size,
    )
    temporary.chmod(mode)
    os.replace(temporary, destination)


def _download(
    url: str,
    destination: Path,
    *,
    maximum_bytes: int,
    runner: CommandRunner,
) -> None:
    if destination.exists() or destination.is_symlink():
        raise StagingError("staging_collision")
    fixed_url = _fixed_https_url(url)
    _run_no_capture(
        runner,
        (
            CURL_EXECUTABLE,
            "--fail",
            "--location",
            "--silent",
            "--show-error",
            "--proto",
            "=https",
            "--max-filesize",
            str(maximum_bytes),
            "--output",
            str(destination),
            fixed_url,
        ),
    )
    _regular_file_size(destination, maximum_bytes=maximum_bytes)


def _run_capture(runner: CommandRunner, argv: tuple[str, ...]) -> bytes:
    try:
        content = runner(argv, capture_stdout=True)
    except StagingError:
        raise
    except Exception as exc:
        raise StagingError("external_command_failed") from exc
    if not isinstance(content, bytes) or len(content) > _MAX_COMMAND_OUTPUT_BYTES:
        raise StagingError("external_command_failed")
    return content


def _run_no_capture(runner: CommandRunner, argv: tuple[str, ...]) -> None:
    try:
        content = runner(argv, capture_stdout=False)
    except StagingError:
        raise
    except Exception as exc:
        raise StagingError("external_command_failed") from exc
    if content not in {b"", None}:
        raise StagingError("external_command_failed")


def _run_external(
    argv: tuple[str, ...],
    *,
    capture_stdout: bool,
    command_home: Path,
) -> bytes:
    if (
        not isinstance(argv, tuple)
        or not argv
        or not all(
            isinstance(value, str)
            and value
            and not any(character in value for character in ("\x00", "\n", "\r"))
            for value in argv
        )
    ):
        raise StagingError("external_command_invalid")
    environment = {
        "DOCKER_CONFIG": str(command_home / "docker"),
        "HOME": str(command_home),
        "LANG": "C.UTF-8",
        "LC_ALL": "C.UTF-8",
        "PATH": "/usr/bin:/bin",
    }
    try:
        completed = subprocess.run(
            argv,
            check=False,
            env=environment,
            stdin=subprocess.DEVNULL,
            stdout=subprocess.PIPE if capture_stdout else subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            timeout=_COMMAND_TIMEOUT_SECONDS,
            shell=False,
        )
    except (OSError, subprocess.TimeoutExpired) as exc:
        raise StagingError("external_command_failed") from exc
    if completed.returncode != 0:
        raise StagingError("external_command_failed")
    output = completed.stdout if capture_stdout else b""
    if not isinstance(output, bytes) or len(output) > _MAX_COMMAND_OUTPUT_BYTES:
        raise StagingError("external_command_failed")
    return output


def _binding(
    root: Path,
    options: StageOptions,
    identity: Mapping[str, Any],
) -> RuntimeAssetBinding:
    source = _mapping(identity.get("source"))
    nextflow = _mapping(identity.get("nextflow"))
    jdk = _mapping(identity.get("jdk"))
    plugin = _mapping(_sequence(identity.get("plugins"))[0])
    try:
        return RuntimeAssetBinding(
            root=root,
            source_tree=str(source["default_relative_path"]),
            nextflow_executable=str(nextflow["default_relative_path"]),
            jdk_archive=str(jdk["default_archive_relative_path"]),
            jdk_tree=str(jdk["default_tree_relative_path"]),
            plugin_archive=str(plugin["default_archive_relative_path"]),
            plugin_meta=str(plugin["default_meta_relative_path"]),
            plugin_tree=str(plugin["default_tree_relative_path"]),
            docker_executable=options.docker_executable,
            docker_socket=options.docker_socket,
        )
    except (TypeError, ValueError) as exc:
        raise StagingError("asset_binding_invalid") from exc


def _verify_production_admission(binding: RuntimeAssetBinding) -> None:
    result = verify_runtime_assets(binding)
    if result.is_failure:
        raise StagingError("runtime_admission_failed")


def _observed_tree(
    root: Path,
    *,
    maximum_files: int,
    maximum_total_bytes: int,
) -> list[dict[str, object]]:
    try:
        root_info = os.lstat(root)
    except OSError as exc:
        raise StagingError("asset_identity_invalid") from exc
    if not stat.S_ISDIR(root_info.st_mode):
        raise StagingError("asset_identity_invalid")
    observed: list[dict[str, object]] = []
    total = 0

    def walk(directory: Path, prefix: str) -> None:
        nonlocal total
        try:
            entries = sorted(os.scandir(directory), key=lambda value: value.name)
        except OSError as exc:
            raise StagingError("asset_identity_invalid") from exc
        if not entries:
            raise StagingError("asset_identity_invalid")
        for entry in entries:
            relative = f"{prefix}/{entry.name}" if prefix else entry.name
            info = entry.stat(follow_symlinks=False)
            if stat.S_ISDIR(info.st_mode):
                walk(Path(entry.path), relative)
                continue
            if (
                not stat.S_ISREG(info.st_mode)
                or info.st_nlink != 1
                or info.st_mode & 0o022
                or len(observed) >= maximum_files
                or info.st_size > _MAX_TREE_FILE_BYTES
            ):
                raise StagingError("asset_identity_invalid")
            total += info.st_size
            if total > maximum_total_bytes:
                raise StagingError("asset_identity_invalid")
            observed.append(
                {
                    "path": relative,
                    "size_bytes": info.st_size,
                    "sha256": _sha256_file(
                        Path(entry.path),
                        _MAX_TREE_FILE_BYTES,
                    ),
                    "executable": bool(info.st_mode & 0o111),
                }
            )

    walk(root, "")
    return sorted(observed, key=lambda value: str(value["path"]))


def _verify_file(
    path: Path,
    *,
    expected_size: int | None,
    expected_sha256: str,
    executable: bool | None,
    maximum_bytes: int,
) -> None:
    size = _regular_file_size(path, maximum_bytes=maximum_bytes)
    if expected_size is not None and size != expected_size:
        raise StagingError("asset_identity_invalid")
    if _sha256_file(path, maximum_bytes) != expected_sha256:
        raise StagingError("asset_identity_invalid")
    if executable is not None:
        info = os.lstat(path)
        if bool(info.st_mode & 0o111) is not executable:
            raise StagingError("asset_identity_invalid")


def _regular_file_size(path: Path, *, maximum_bytes: int) -> int:
    try:
        info = os.lstat(path)
    except OSError as exc:
        raise StagingError("asset_identity_invalid") from exc
    if (
        not stat.S_ISREG(info.st_mode)
        or info.st_nlink != 1
        or info.st_size < 0
        or info.st_size > maximum_bytes
    ):
        raise StagingError("asset_identity_invalid")
    return info.st_size


def _read_regular_file(path: Path, maximum_bytes: int) -> bytes:
    size = _regular_file_size(path, maximum_bytes=maximum_bytes)
    try:
        content = path.read_bytes()
    except OSError as exc:
        raise StagingError("asset_identity_invalid") from exc
    if len(content) != size:
        raise StagingError("asset_identity_invalid")
    return content


def _sha256_file(path: Path, maximum_bytes: int) -> str:
    _regular_file_size(path, maximum_bytes=maximum_bytes)
    digest = hashlib.sha256()
    total = 0
    try:
        with path.open("rb") as handle:
            while True:
                chunk = handle.read(1024 * 1024)
                if not chunk:
                    break
                total += len(chunk)
                if total > maximum_bytes:
                    raise StagingError("asset_identity_invalid")
                digest.update(chunk)
    except OSError as exc:
        raise StagingError("asset_identity_invalid") from exc
    return digest.hexdigest()


def _atomic_write(path: Path, content: bytes, *, mode: int) -> None:
    if not isinstance(content, bytes):
        raise StagingError("staging_write_invalid")
    path.parent.mkdir(parents=True, exist_ok=True, mode=0o700)
    descriptor = -1
    temporary: Path | None = None
    try:
        descriptor, rendered = tempfile.mkstemp(
            prefix=f".{path.name}.part-",
            dir=path.parent,
        )
        temporary = Path(rendered)
        os.fchmod(descriptor, mode)
        with os.fdopen(descriptor, "wb", closefd=False) as handle:
            handle.write(content)
            handle.flush()
            os.fsync(handle.fileno())
        os.close(descriptor)
        descriptor = -1
        os.replace(temporary, path)
        temporary = None
    except OSError as exc:
        raise StagingError("staging_write_failed") from exc
    finally:
        if descriptor >= 0:
            os.close(descriptor)
        if temporary is not None:
            try:
                temporary.unlink()
            except OSError:
                pass


def _publish_candidate(candidate: Path, asset_root: Path) -> None:
    """Atomically publish without ever replacing a concurrently created root."""

    try:
        library = ctypes.CDLL(None, use_errno=True)
        renameat2 = library.renameat2
    except (AttributeError, OSError) as exc:
        raise StagingError("asset_publish_unsupported") from exc
    renameat2.argtypes = (
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_int,
        ctypes.c_char_p,
        ctypes.c_uint,
    )
    renameat2.restype = ctypes.c_int
    result = renameat2(
        _AT_FDCWD,
        os.fsencode(candidate),
        _AT_FDCWD,
        os.fsencode(asset_root),
        _RENAME_NOREPLACE,
    )
    if result == 0:
        return
    error_number = ctypes.get_errno()
    if error_number in {errno.EEXIST, errno.ENOTEMPTY}:
        raise StagingError("asset_root_exists")
    if error_number in {errno.EINVAL, errno.ENOSYS, errno.EOPNOTSUPP}:
        raise StagingError("asset_publish_unsupported")
    raise StagingError("asset_publish_failed")


def _remove_published_root(root: Path, expected: tuple[int, int]) -> None:
    try:
        current = os.lstat(root)
    except OSError:
        return
    if (
        stat.S_ISDIR(current.st_mode)
        and (current.st_dev, current.st_ino) == expected
        and not root.is_symlink()
    ):
        shutil.rmtree(root, ignore_errors=True)


def _read_tar_member(
    package: tarfile.TarFile,
    member: tarfile.TarInfo,
    expected_size: int,
    code: str,
) -> bytes:
    if member.size != expected_size or member.size > _MAX_TREE_FILE_BYTES:
        raise StagingError(code)
    extracted = package.extractfile(member)
    if extracted is None:
        raise StagingError(code)
    content = extracted.read(expected_size + 1)
    if len(content) != expected_size:
        raise StagingError(code)
    return content


def _archive_parts(value: object, code: str) -> tuple[str, ...]:
    if not isinstance(value, str) or "\\" in value:
        raise StagingError(code)
    path = PurePosixPath(value.rstrip("/"))
    if path.is_absolute() or any(part in {"", ".", ".."} for part in path.parts):
        raise StagingError(code)
    if any(any(ord(character) < 32 for character in part) for part in path.parts):
        raise StagingError(code)
    return path.parts


def _relative_path(value: object, code: str) -> str:
    if not isinstance(value, str) or not value or "\\" in value:
        raise StagingError(code)
    path = PurePosixPath(value)
    if (
        path.is_absolute()
        or any(part in {"", ".", ".."} for part in path.parts)
        or any(ord(character) < 32 for character in value)
    ):
        raise StagingError(code)
    return path.as_posix()


def _single_component(value: object, code: str) -> str:
    normalized = _relative_path(value, code)
    if len(PurePosixPath(normalized).parts) != 1:
        raise StagingError(code)
    return normalized


def _absolute_path(path: Path, code: str, *, allow_missing: bool) -> None:
    if not isinstance(path, Path):
        raise StagingError(code)
    rendered = str(path)
    if (
        not path.is_absolute()
        or rendered != str(Path(rendered))
        or any(part == ".." for part in path.parts)
        or any(
            character in rendered
            for character in ("\x00", "\n", "\r", "'", "\\", "$", ":")
        )
    ):
        raise StagingError(code)
    if not allow_missing and not path.exists():
        raise StagingError(code)


def _reject_symlink_ancestors(path: Path, code: str) -> None:
    current = Path(path.anchor)
    for part in path.parts[1:]:
        current /= part
        try:
            info = os.lstat(current)
        except OSError as exc:
            raise StagingError(code) from exc
        if stat.S_ISLNK(info.st_mode) or not stat.S_ISDIR(info.st_mode):
            raise StagingError(code)


def _strict_json(content: bytes) -> Any:
    def reject_constant(_value: str) -> None:
        raise ValueError("non-finite JSON number")

    def unique_object(pairs: list[tuple[str, Any]]) -> dict[str, Any]:
        value: dict[str, Any] = {}
        for key, item in pairs:
            if key in value:
                raise ValueError("duplicate JSON key")
            value[key] = item
        return value

    return json.loads(
        content,
        object_pairs_hook=unique_object,
        parse_constant=reject_constant,
    )


def _canonical_tree_digest(entries: list[dict[str, object]]) -> str:
    ordered = sorted(
        entries,
        key=lambda entry: PurePosixPath(str(entry["path"])),
    )
    return _sha256_bytes(
        json.dumps(
            ordered,
            ensure_ascii=False,
            sort_keys=True,
            separators=(",", ":"),
        ).encode("utf-8")
    )


def _sha256_bytes(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _json_bytes(value: object) -> bytes:
    return (
        json.dumps(value, ensure_ascii=False, sort_keys=True, separators=(",", ":"))
        + "\n"
    ).encode()


def _required_digest(value: object, code: str, *, prefixed: bool) -> str:
    pattern = _DIGEST_RE if prefixed else re.compile(r"^[a-f0-9]{64}$")
    if not isinstance(value, str) or pattern.fullmatch(value) is None:
        raise StagingError(code)
    return value


def _positive_int(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value <= 0:
        raise StagingError("committed_contract_invalid")
    return value


def _positive_nonnegative_int(value: object) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise StagingError("committed_contract_invalid")
    return value


def _mapping(value: object) -> Mapping[str, Any]:
    if not isinstance(value, Mapping):
        raise StagingError("committed_contract_invalid")
    return value


def _sequence(value: object) -> list[Any]:
    if not isinstance(value, list) or not value:
        raise StagingError("committed_contract_invalid")
    return value


def _fixed_https_url(value: object) -> str:
    if (
        not isinstance(value, str)
        or not value.startswith("https://")
        or any(character in value for character in ("\x00", "\n", "\r", " "))
    ):
        raise StagingError("download_url_invalid")
    return value


def _plan_json(report: StagingReport) -> str:
    return json.dumps(
        {
            "phase": report.phase,
            "container_process_count": report.container_process_count,
            "unique_image_count": report.unique_image_count,
            "source_tree_sha256": report.source_tree_sha256,
            "container_lock_sha256": report.container_lock_sha256,
            "network_used": report.phase == "stage",
            "mutation_used": report.phase == "stage",
        },
        sort_keys=True,
        separators=(",", ":"),
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description="Stage the pinned nf-core/rnaseq runtime asset closure.",
    )
    parser.add_argument("--asset-root", required=True, type=Path)
    parser.add_argument("--docker-executable", required=True, type=Path)
    parser.add_argument("--docker-socket", required=True, type=Path)
    parser.add_argument("--phase", choices=sorted(_PHASES), default="plan")
    parser.add_argument("--allow-network", action="store_true")
    parser.add_argument("--allow-mutation", action="store_true")
    args = parser.parse_args(argv)
    try:
        report = stage_runtime_assets(
            StageOptions(
                asset_root=args.asset_root,
                docker_executable=args.docker_executable,
                docker_socket=args.docker_socket,
                phase=args.phase,
                allow_network=args.allow_network,
                allow_mutation=args.allow_mutation,
            )
        )
    except StagingError as exc:
        print(f"staging failed: {exc.code}", file=sys.stderr)
        return 1
    print(_plan_json(report))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
