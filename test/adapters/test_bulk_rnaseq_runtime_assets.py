"""Closed-set and no-follow tests for the bulk RNA-seq runtime assets."""

from __future__ import annotations

from collections.abc import Callable
from copy import deepcopy
from dataclasses import dataclass, replace
import hashlib
import io
import json
import os
from pathlib import Path
import stat
import subprocess
import tarfile
from types import SimpleNamespace
import zipfile

import pytest

import encode_pipeline.adapters.bulk_rnaseq.runtime_assets as runtime_assets_module
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    CONTAINER_CONFIG_ASSIGNMENTS_SHA256,
    CONTAINER_CONFIG_ASSIGNMENT_COUNT,
    CONTAINER_DEFAULT_CONFIG_ASSIGNMENTS_SHA256,
    CONTAINER_DEFAULT_CONFIG_ASSIGNMENT_COUNT,
    CONTAINER_DEFAULT_CONFIG_DENIED_SELECTORS,
    CONTAINER_EXCLUDED_PROCESS_COUNT,
    CONTAINER_INVENTORY_ENTRIES_SHA256,
    CONTAINER_PROCESS_AUDIT_PROCESSES_SHA256,
    CONTAINER_PROCESS_COUNT,
    CONTAINER_PROCESS_UNIVERSE_COUNT,
    CONTAINER_RESERVED_DEFAULT_DENY_LABEL,
    DOCKER_REQUIRED_RUN_OPTIONS,
    JAVA_EXECUTABLE_SHA256,
    JDK_ARCHITECTURE,
    JDK_ARCHIVE_SHA256,
    JDK_RUNTIME_VERSION,
    JDK_TREE_FILE_COUNT,
    JDK_TREE_SHA256,
    JDK_VENDOR,
    JDK_VERSION,
    NEXTFLOW_OFFLINE_ENV,
    NFCORE_RNASEQ_COMMIT,
    NF_SCHEMA_ARCHIVE_SHA256,
    NF_SCHEMA_VERSION,
    NEXTFLOW_SHA256,
    NEXTFLOW_VERSION,
    SOURCE_FILE_COUNT,
    SOURCE_TREE_SHA256,
    RuntimeAssetBinding,
    RuntimeAssetAdmission,
    _RuntimeAssetContract,
    _load_runtime_contract,
    _validate_embedded_contracts,
    _verify_docker_availability,
    doctor_runtime_assets,
    verify_runtime_assets,
)


@dataclass
class _Fixture:
    binding: RuntimeAssetBinding
    contract: _RuntimeAssetContract
    root: Path
    source: Path
    nextflow: Path
    jdk_archive: Path
    jdk_tree: Path
    java_executable: Path
    plugin_tree: Path
    container_lock: Path
    container_asset: Path
    distribution_manifest: Path
    config_digest: str
    rootfs_diff_ids: tuple[str, ...]
    docker_probe: Callable[[RuntimeAssetBinding, tuple[str, ...]], bytes]
    lock_payload: dict[str, object]


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _canonical_sha256(value: object) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _write_json(path: Path, value: object) -> None:
    path.write_text(json.dumps(value, indent=2) + "\n", encoding="utf-8")


def _json_bytes(value: object) -> bytes:
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def _tar_bytes(files: dict[str, bytes]) -> bytes:
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w", format=tarfile.USTAR_FORMAT) as package:
        for name, content in sorted(files.items()):
            member = tarfile.TarInfo(name)
            member.size = len(content)
            member.mode = 0o644
            member.mtime = 0
            member.uid = 0
            member.gid = 0
            member.uname = ""
            member.gname = ""
            package.addfile(member, io.BytesIO(content))
    return output.getvalue()


@dataclass(frozen=True)
class _DockerImageFixture:
    archive: bytes
    distribution_manifest: bytes
    config_digest: str
    layer_tar: bytes
    rootfs_diff_ids: tuple[str, ...]


def _docker_image_fixture(
    *,
    layer_tar: bytes | None = None,
    diff_id: str | None = None,
) -> _DockerImageFixture:
    layer = layer_tar if layer_tar is not None else _tar_bytes({"tiny.txt": b"tiny\n"})
    configured_diff_id = diff_id or f"sha256:{_sha256(layer)}"
    config = _json_bytes(
        {
            "architecture": "amd64",
            "os": "linux",
            "rootfs": {"type": "layers", "diff_ids": [configured_diff_id]},
        }
    )
    config_digest = f"sha256:{_sha256(config)}"
    archive_manifest = _json_bytes(
        [
            {
                "Config": f"{config_digest.removeprefix('sha256:')}.json",
                "RepoTags": ["quay.io/helixweave/star:test"],
                "Layers": ["layer/layer.tar"],
            }
        ]
    )
    archive = _tar_bytes(
        {
            f"{config_digest.removeprefix('sha256:')}.json": config,
            "layer/layer.tar": layer,
            "manifest.json": archive_manifest,
        }
    )
    distribution_manifest = _json_bytes(
        {
            "schemaVersion": 2,
            "mediaType": "application/vnd.docker.distribution.manifest.v2+json",
            "config": {
                "mediaType": "application/vnd.docker.container.image.v1+json",
                "size": len(config),
                "digest": config_digest,
            },
            "layers": [
                {
                    "mediaType": "application/vnd.docker.image.rootfs.diff.tar.gzip",
                    "size": len(layer),
                    "digest": f"sha256:{_sha256(layer)}",
                }
            ],
        }
    )
    return _DockerImageFixture(
        archive=archive,
        distribution_manifest=distribution_manifest,
        config_digest=config_digest,
        layer_tar=layer,
        rootfs_diff_ids=(configured_diff_id,),
    )


def _replace_archive(fixture: _Fixture, content: bytes) -> None:
    fixture.container_asset.write_bytes(content)
    entries = fixture.lock_payload["entries"]
    assert isinstance(entries, list)
    entries[0]["size_bytes"] = len(content)
    entries[0]["sha256"] = _sha256(content)
    _write_json(fixture.container_lock, fixture.lock_payload)


def _unsafe_tar(kind: str) -> bytes:
    output = io.BytesIO()
    with tarfile.open(fileobj=output, mode="w", format=tarfile.USTAR_FORMAT) as package:
        first = tarfile.TarInfo("unsafe")
        first.mtime = 0
        if kind == "traversal":
            first.name = "../escape"
            first.size = 1
            package.addfile(first, io.BytesIO(b"x"))
        elif kind == "symlink":
            first.type = tarfile.SYMTYPE
            first.linkname = "target"
            package.addfile(first)
        elif kind == "device":
            first.type = tarfile.CHRTYPE
            package.addfile(first)
        else:
            first.size = 1
            package.addfile(first, io.BytesIO(b"x"))
            duplicate = tarfile.TarInfo("unsafe")
            duplicate.size = 1
            duplicate.mtime = 0
            package.addfile(duplicate, io.BytesIO(b"y"))
    return output.getvalue()


def _tar_file_contents(content: bytes) -> dict[str, bytes]:
    files: dict[str, bytes] = {}
    with tarfile.open(fileobj=io.BytesIO(content), mode="r:") as package:
        for member in package:
            if not member.isreg():
                continue
            extracted = package.extractfile(member)
            assert extracted is not None
            files[member.name] = extracted.read()
    return files


def _tiny_assets(tmp_path: Path) -> _Fixture:
    root = tmp_path / "offline-assets"
    source = root / ("source/nf-core-rnaseq-e7ca46272c8f9d5ceee3f71759f4ba551d3217a4")
    source_files = {
        "bin/tool.sh": b"#!/bin/sh\nexit 0\n",
        "conf/settings.config": b"params.enabled = true\n",
        "main.nf": b"nextflow.enable.dsl=2\n",
    }
    for relative, content in source_files.items():
        path = source / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        path.chmod(0o755 if relative == "bin/tool.sh" else 0o644)
    source_entries = [
        {
            "path": relative,
            "size_bytes": len(content),
            "sha256": _sha256(content),
            "executable": relative == "bin/tool.sh",
        }
        for relative, content in sorted(source_files.items())
    ]
    source_tree_sha256 = _canonical_sha256(source_entries)
    source_manifest = {
        "file_count": len(source_entries),
        "tree_sha256": source_tree_sha256,
        "files": source_entries,
    }

    nextflow_content = b"#!/bin/sh\necho 25.04.3\n"
    nextflow = root / "nextflow/nextflow-25.04.3-dist"
    nextflow.parent.mkdir(parents=True)
    nextflow.write_bytes(nextflow_content)
    nextflow.chmod(0o755)

    jdk_archive_content = b"synthetic Corretto archive"
    jdk_archive = root / "jdk/amazon-corretto-21.0.7.6.1-linux-x64.tar.gz"
    jdk_archive.parent.mkdir(parents=True)
    jdk_archive.write_bytes(jdk_archive_content)
    jdk_tree = root / "jdk/amazon-corretto-21.0.7.6.1-linux-x64"
    jdk_files = {
        "ASSEMBLY_EXCEPTION": b"synthetic classpath exception\n",
        "LICENSE": b"synthetic GPLv2 license\n",
        "bin/java": b"#!/bin/sh\necho synthetic java >&2\n",
    }
    for relative, content in jdk_files.items():
        path = jdk_tree / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        path.chmod(0o755 if relative == "bin/java" else 0o644)
    jdk_entries = [
        {
            "path": relative,
            "size_bytes": len(content),
            "sha256": _sha256(content),
            "executable": relative == "bin/java",
        }
        for relative, content in sorted(jdk_files.items())
    ]

    plugin_file = b"synthetic plugin bytecode"
    archive_buffer = io.BytesIO()
    with zipfile.ZipFile(archive_buffer, "w") as package:
        package.writestr("classes/Test.class", plugin_file)
    archive = archive_buffer.getvalue()
    meta = json.dumps(
        {
            "version": "2.5.1",
            "requires": ">=25.04.0",
            "sha512sum": hashlib.sha512(archive).hexdigest(),
        },
        separators=(",", ":"),
    ).encode("utf-8")
    plugin_archive = root / "plugins/nf-schema-2.5.1.zip"
    plugin_archive.parent.mkdir(parents=True)
    plugin_archive.write_bytes(archive)
    plugin_meta = root / "plugins/nf-schema-2.5.1-meta.json"
    plugin_meta.write_bytes(meta)
    plugin_tree = root / "plugins/nf-schema-2.5.1"
    plugin_tree.joinpath("classes").mkdir(parents=True)
    plugin_tree.joinpath("classes/Test.class").write_bytes(plugin_file)
    plugin_tree.joinpath("nf-schema-2.5.1-meta.json").write_bytes(meta)
    plugin_entries = [
        {
            "path": "classes/Test.class",
            "size_bytes": len(plugin_file),
            "sha256": _sha256(plugin_file),
        },
        {
            "path": "nf-schema-2.5.1-meta.json",
            "size_bytes": len(meta),
            "sha256": _sha256(meta),
        },
    ]

    inventory_entries = [
        {
            "process": "STAR_ALIGN",
            "image_coordinate": "quay.io/helixweave/star:test",
            "source_file": "main.nf",
            "source_file_sha256": _sha256(source_files["main.nf"]),
        }
    ]
    inventory_sha256 = _canonical_sha256(inventory_entries)
    inventory = {
        "entries": inventory_entries,
        "entries_sha256": inventory_sha256,
    }

    docker_image = _docker_image_fixture()
    container_content = docker_image.archive
    container_asset = root / "containers/assets/star.tar"
    container_asset.parent.mkdir(parents=True)
    container_asset.write_bytes(container_content)
    distribution_manifest = root / "containers/assets/star.manifest.json"
    distribution_manifest.write_bytes(docker_image.distribution_manifest)
    oci_digest = f"sha256:{_sha256(docker_image.distribution_manifest)}"
    lock_payload: dict[str, object] = {
        "schema_version": "1.0.0",
        "pipeline": {
            "project": "nf-core/rnaseq",
            "release": "3.26.0",
            "commit": NFCORE_RNASEQ_COMMIT,
            "source_tree_sha256": source_tree_sha256,
            "container_inventory_sha256": inventory_sha256,
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
        "entries": [
            {
                "process": "STAR_ALIGN",
                "image": f"quay.io/helixweave/star:test@{oci_digest}",
                "oci_digest": oci_digest,
                "distribution_manifest_asset": "star.manifest.json",
                "distribution_manifest_size_bytes": len(
                    docker_image.distribution_manifest
                ),
                "local_asset": "star.tar",
                "size_bytes": len(container_content),
                "sha256": _sha256(container_content),
            }
        ],
    }
    container_lock = root / "containers/availability-lock.json"
    _write_json(container_lock, lock_payload)

    production = _load_runtime_contract()
    container_schema = deepcopy(production.container_schema)
    pipeline = container_schema["properties"]["pipeline"]["properties"]
    pipeline["source_tree_sha256"]["const"] = source_tree_sha256
    pipeline["container_inventory_sha256"]["const"] = inventory_sha256
    identity = {
        "nextflow": {
            "build": "5949",
            "size_bytes": len(nextflow_content),
            "sha256": _sha256(nextflow_content),
        },
        "jdk": {
            "archive_size_bytes": len(jdk_archive_content),
            "archive_sha256": _sha256(jdk_archive_content),
            "tree_file_count": len(jdk_entries),
            "tree_size_bytes": sum(len(content) for content in jdk_files.values()),
            "tree_sha256": _canonical_sha256(jdk_entries),
            "java_relative_path": "bin/java",
            "java_size_bytes": len(jdk_files["bin/java"]),
            "java_sha256": _sha256(jdk_files["bin/java"]),
            "license_file": "LICENSE",
            "license_file_sha256": _sha256(jdk_files["LICENSE"]),
            "assembly_exception_file": "ASSEMBLY_EXCEPTION",
            "assembly_exception_file_sha256": _sha256(jdk_files["ASSEMBLY_EXCEPTION"]),
            "java_version_output_sha256": _sha256(b"synthetic java version\n"),
        },
        "plugins": [
            {
                "archive_size_bytes": len(archive),
                "archive_sha256": _sha256(archive),
                "meta_size_bytes": len(meta),
                "meta_sha256": _sha256(meta),
                "tree_sha256": _canonical_sha256(plugin_entries),
            }
        ],
    }
    contract = _RuntimeAssetContract(
        identity=identity,
        source_manifest=source_manifest,
        container_schema=container_schema,
        container_inventory=inventory,
        container_process_audit={
            "reserved_default_deny_label": CONTAINER_RESERVED_DEFAULT_DENY_LABEL
        },
    )

    def exact_docker_probe(
        _binding: RuntimeAssetBinding,
        images: tuple[str, ...],
    ) -> bytes:
        assert images == (docker_image.config_digest,)
        return _json_bytes(
            [
                {
                    "Id": docker_image.config_digest,
                    "RepoDigests": [],
                    "RootFS": {
                        "Type": "layers",
                        "Layers": [*docker_image.rootfs_diff_ids],
                    },
                }
                for _image in images
            ]
        )

    return _Fixture(
        binding=RuntimeAssetBinding(root=root),
        contract=contract,
        root=root,
        source=source,
        nextflow=nextflow,
        jdk_archive=jdk_archive,
        jdk_tree=jdk_tree,
        java_executable=jdk_tree / "bin/java",
        plugin_tree=plugin_tree,
        container_lock=container_lock,
        container_asset=container_asset,
        distribution_manifest=distribution_manifest,
        config_digest=docker_image.config_digest,
        rootfs_diff_ids=docker_image.rootfs_diff_ids,
        docker_probe=exact_docker_probe,
        lock_payload=lock_payload,
    )


def _add_process_sharing_container_coordinate(fixture: _Fixture) -> None:
    inventory_entries = fixture.contract.container_inventory["entries"]
    assert isinstance(inventory_entries, list)
    shared_inventory = deepcopy(inventory_entries[0])
    shared_inventory["process"] = "ZZZ_ALIGN"
    inventory_entries.append(shared_inventory)
    inventory_sha256 = _canonical_sha256(inventory_entries)
    fixture.contract.container_inventory["entries_sha256"] = inventory_sha256
    pipeline_schema = fixture.contract.container_schema["properties"]["pipeline"][
        "properties"
    ]
    pipeline_schema["container_inventory_sha256"]["const"] = inventory_sha256
    pipeline_lock = fixture.lock_payload["pipeline"]
    assert isinstance(pipeline_lock, dict)
    pipeline_lock["container_inventory_sha256"] = inventory_sha256
    lock_entries = fixture.lock_payload["entries"]
    assert isinstance(lock_entries, list)
    shared_lock = deepcopy(lock_entries[0])
    shared_lock["process"] = "ZZZ_ALIGN"
    lock_entries.append(shared_lock)
    _write_json(fixture.container_lock, fixture.lock_payload)


def test_committed_runtime_contracts_are_immutable_and_offline() -> None:
    contract = _load_runtime_contract()
    _validate_embedded_contracts(
        contract.identity,
        contract.source_manifest,
        contract.container_inventory,
        contract.container_process_audit,
    )

    assert contract.source_manifest["file_count"] == SOURCE_FILE_COUNT == 829
    assert contract.source_manifest["tree_sha256"] == SOURCE_TREE_SHA256
    assert contract.identity["nextflow"] == {
        "version": NEXTFLOW_VERSION,
        "build": "5949",
        "default_relative_path": "nextflow/nextflow-25.04.3-dist",
        "size_bytes": 31_307_907,
        "sha256": NEXTFLOW_SHA256,
    }
    jdk = contract.identity["jdk"]
    assert jdk["vendor"] == JDK_VENDOR == "Amazon Corretto"
    assert jdk["version"] == JDK_VERSION == "21.0.7.6.1"
    assert jdk["runtime_version"] == JDK_RUNTIME_VERSION == "21.0.7+6-LTS"
    assert jdk["operating_system"] == "linux"
    assert jdk["architecture"] == JDK_ARCHITECTURE == "x64"
    assert jdk["archive_size_bytes"] == 208_603_382
    assert jdk["archive_sha256"] == JDK_ARCHIVE_SHA256
    assert jdk["archive_url"].startswith(
        "https://github.com/corretto/corretto-21/releases/download/21.0.7.6.1/"
    )
    assert jdk["tree_file_count"] == JDK_TREE_FILE_COUNT == 457
    assert jdk["tree_size_bytes"] == 362_486_940
    assert jdk["tree_sha256"] == JDK_TREE_SHA256
    assert jdk["java_size_bytes"] == 12_944
    assert jdk["java_sha256"] == JAVA_EXECUTABLE_SHA256
    assert jdk["license"] == "GPL-2.0-only WITH Classpath-exception-2.0"
    assert jdk["license_url"].endswith("/blob/21.0.7.6.1/LICENSE")
    plugin = contract.identity["plugins"][0]
    assert plugin["version"] == NF_SCHEMA_VERSION
    assert plugin["archive_sha256"] == NF_SCHEMA_ARCHIVE_SHA256
    assert contract.container_inventory["process_count"] == CONTAINER_PROCESS_COUNT
    assert (
        contract.container_inventory["entries_sha256"]
        == CONTAINER_INVENTORY_ENTRIES_SHA256
    )
    assert "RIBODETECTOR" not in {
        entry["process"] for entry in contract.container_inventory["entries"]
    }
    processes = contract.container_process_audit["processes"]
    included = {
        entry["process"] for entry in processes if entry["disposition"] == "included"
    }
    excluded = {
        entry["process"] for entry in processes if entry["disposition"] == "excluded"
    }
    assert len(processes) == CONTAINER_PROCESS_UNIVERSE_COUNT == 78
    assert len(included) == CONTAINER_PROCESS_COUNT == 56
    assert len(excluded) == CONTAINER_EXCLUDED_PROCESS_COUNT == 22
    assert not included & excluded
    assert included | excluded == {entry["process"] for entry in processes}
    assert included == {
        entry["process"] for entry in contract.container_inventory["entries"]
    }
    assert "RIBODETECTOR" in excluded
    assert (
        contract.container_process_audit["processes_sha256"]
        == CONTAINER_PROCESS_AUDIT_PROCESSES_SHA256
    )
    assert (
        contract.container_process_audit["reserved_default_deny_label"]
        == CONTAINER_RESERVED_DEFAULT_DENY_LABEL
    )
    default_config_assignments = contract.container_process_audit[
        "default_config_container_assignments"
    ]
    assert (
        contract.container_process_audit["config_container_assignment_count"]
        == CONTAINER_CONFIG_ASSIGNMENT_COUNT
        == 78
    )
    assert (
        contract.container_process_audit["config_container_assignments_sha256"]
        == CONTAINER_CONFIG_ASSIGNMENTS_SHA256
    )
    assert (
        len(default_config_assignments)
        == CONTAINER_DEFAULT_CONFIG_ASSIGNMENT_COUNT
        == 1
    )
    assert (
        contract.container_process_audit["default_config_container_assignments_sha256"]
        == CONTAINER_DEFAULT_CONFIG_ASSIGNMENTS_SHA256
    )
    assert tuple(entry["selector"] for entry in default_config_assignments) == (
        CONTAINER_DEFAULT_CONFIG_DENIED_SELECTORS
    )
    assert default_config_assignments[0]["platform_override"] == "deny"
    assert default_config_assignments[0]["matched_aliases"] == [
        {
            "alias": "PARABRICKS_STARGENOMEGENERATE",
            "source_process": "STAR_GENOMEGENERATE",
            "source_file": "subworkflows/local/prepare_genome/main.nf",
            "source_file_sha256": (
                "e565b2e67e460cc353361b927067db81ae14169307003a261fa27e6900c9645f"
            ),
        }
    ]
    container_identity = contract.identity["containers"]
    assert (
        container_identity["config_container_assignment_count"]
        == CONTAINER_CONFIG_ASSIGNMENT_COUNT
    )
    assert (
        container_identity["config_container_assignments_sha256"]
        == CONTAINER_CONFIG_ASSIGNMENTS_SHA256
    )
    assert (
        container_identity["default_config_container_assignment_count"]
        == CONTAINER_DEFAULT_CONFIG_ASSIGNMENT_COUNT
    )
    assert (
        container_identity["default_config_container_assignments_sha256"]
        == CONTAINER_DEFAULT_CONFIG_ASSIGNMENTS_SHA256
    )
    assert not any(contract.identity["network_policy"].values())
    assert contract.identity["execution_policy"]["required_docker_run_options"] == [
        *DOCKER_REQUIRED_RUN_OPTIONS
    ]
    assert (
        contract.identity["containers"]["asset_format"]
        == "docker-archive+distribution-manifest"
    )
    assert NEXTFLOW_OFFLINE_ENV == {"NXF_OFFLINE": "true"}


def test_committed_contract_records_real_scale_admission_cost() -> None:
    contract = _load_runtime_contract()
    unique_images = {
        entry["image_coordinate"] for entry in contract.container_inventory["entries"]
    }
    plugin_files = contract.identity["plugins"][0]["tree_file_count"]

    assert SOURCE_FILE_COUNT == 829
    assert JDK_TREE_FILE_COUNT == 457
    assert plugin_files == 106
    assert len(unique_images) == 34
    # One heavyweight admission performs these top-level content passes:
    # source files, Nextflow, JDK archive/tree, plugin archive/meta/tree,
    # container lock, and one manifest/hash/closure pass per unique OCI asset.
    maximum_top_level_content_passes = (
        SOURCE_FILE_COUNT
        + 1
        + 1
        + JDK_TREE_FILE_COUNT
        + 2
        + plugin_files
        + 1
        + len(unique_images) * 3
    )
    assert maximum_top_level_content_passes == 1_499
    assert 4 * len(unique_images) == 136  # old preflight archive hash passes
    assert 3 * len(unique_images) == 102  # old worker archive hash passes


def test_embedded_contract_rejects_tampered_default_config_container_alias() -> None:
    contract = _load_runtime_contract()
    process_audit = deepcopy(contract.container_process_audit)
    process_audit["default_config_container_assignments"][0]["selector"] = ".*FASTQC"

    with pytest.raises(ValueError, match="container process audit"):
        _validate_embedded_contracts(
            contract.identity,
            contract.source_manifest,
            contract.container_inventory,
            process_audit,
        )


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("vendor", "host-openjdk"),
        ("version", "21.0.8"),
        ("architecture", "aarch64"),
        ("java_sha256", "0" * 64),
        ("java_version_output_sha256", "0" * 64),
    ],
)
def test_embedded_contract_rejects_changed_jdk_identity(
    field: str,
    value: str,
) -> None:
    contract = _load_runtime_contract()
    identity = deepcopy(contract.identity)
    identity["jdk"][field] = value

    with pytest.raises(ValueError, match="JDK identity"):
        _validate_embedded_contracts(
            identity,
            contract.source_manifest,
            contract.container_inventory,
            contract.container_process_audit,
        )


def test_tiny_closed_asset_set_verifies_deterministically(tmp_path: Path) -> None:
    fixture = _tiny_assets(tmp_path)

    first = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )
    second = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert first.is_success
    assert second.is_success
    assert first.value == second.value
    assert first.value is not None
    assert first.value.source_tree == fixture.source
    assert first.value.nextflow_executable == fixture.nextflow
    assert tuple(item.process for item in first.value.containers) == ("STAR_ALIGN",)
    assert first.value.containers[0].runtime_image == fixture.config_digest
    assert first.value.containers[0].rootfs_diff_ids == fixture.rootfs_diff_ids
    assert doctor_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    ).ready


@pytest.mark.parametrize(
    "mutation",
    ["missing", "archive", "extra", "symlink", "fifo", "java_mode"],
)
def test_jdk_archive_tree_and_java_identity_fail_closed(
    tmp_path: Path,
    mutation: str,
) -> None:
    fixture = _tiny_assets(tmp_path)
    if mutation == "missing":
        fixture.jdk_archive.unlink()
    elif mutation == "archive":
        fixture.jdk_archive.write_bytes(b"replaced archive")
    elif mutation == "extra":
        fixture.jdk_tree.joinpath("unexpected").write_bytes(b"extra")
    elif mutation == "symlink":
        fixture.java_executable.unlink()
        fixture.java_executable.symlink_to("/bin/false")
    elif mutation == "fifo":
        fixture.java_executable.unlink()
        os.mkfifo(fixture.java_executable)
    else:
        fixture.java_executable.chmod(0o644)

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert any(issue.context == {"component": "jdk"} for issue in result.issues)


def test_runtime_version_probe_is_mandatory_when_composed(tmp_path: Path) -> None:
    fixture = _tiny_assets(tmp_path)
    calls = 0

    def reject_version(
        _binding: RuntimeAssetBinding,
        _identity: object,
        _containers: object,
    ) -> None:
        nonlocal calls
        calls += 1
        raise runtime_assets_module._AssetFault("jdk", "version")

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
        _runtime_probe=reject_version,
    )

    assert result.is_failure
    assert calls == 1
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_VERSION"
    assert result.issues[0].context == {"component": "jdk"}


def test_runtime_canary_uses_exact_jdk_and_hard_config_with_poison_sources(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    verified = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    ).value
    assert verified is not None
    calls: list[dict[str, object]] = []

    def run(argv, **kwargs):
        argv_tuple = tuple(argv)
        record: dict[str, object] = {"argv": argv_tuple, **kwargs}
        if "-C" in argv_tuple:
            config_path = Path(argv_tuple[argv_tuple.index("-C") + 1])
            record["hard_config"] = config_path.read_text(encoding="utf-8")
            record["launch_poison"] = (
                Path(kwargs["cwd"])
                .joinpath("nextflow.config")
                .read_text(encoding="utf-8")
            )
            environment = kwargs["env"]
            record["home_poison"] = (
                Path(environment["NXF_HOME"])
                .joinpath("config")
                .read_text(encoding="utf-8")
            )
            output = "\n".join(
                (
                    "params.helixweave_runtime_admission_marker = 'hard-config-only'",
                    "process.executor = 'local'",
                    "docker.enabled = true",
                    "docker.registry = ''",
                    "wave.enabled = false",
                    "tower.enabled = false",
                    "fusion.enabled = false",
                    "manifest.name = 'nf-core/rnaseq'",
                    "manifest.version = '3.26.0'",
                    "manifest.nextflowVersion = '!>=25.04.3'",
                    "plugins = ['nf-schema@2.5.1']",
                    "shifter.enabled = false",
                    fixture.config_digest,
                )
            ).encode()
        elif argv_tuple[-1] == "-version" and argv_tuple[0].endswith("bin/java"):
            output = b"synthetic java version\n"
        else:
            output = b"version 25.04.3 build 5949\n"
        calls.append(record)
        return SimpleNamespace(returncode=0, stdout=output)

    monkeypatch.setattr(runtime_assets_module.subprocess, "run", run)

    runtime_assets_module._verify_runtime_canary(
        fixture.binding,
        fixture.contract.identity,
        verified.containers,
    )

    assert len(calls) == 3
    java_call, version_call, config_call = calls
    java_home = fixture.jdk_tree
    assert java_call["argv"] == (str(fixture.java_executable), "-version")
    assert version_call["argv"] == (str(fixture.nextflow), "-version")
    argv = config_call["argv"]
    assert argv.count("-C") == 1
    assert "-c" not in argv
    assert argv.index("-C") < argv.index("config")
    hard_config = config_call["hard_config"]
    assert hard_config.count("includeConfig ") == 1
    assert hard_config.splitlines()[0] == (
        f"includeConfig '{fixture.source}/nextflow.config'"
    )
    assert "helixweave_launch_poison" in config_call["launch_poison"]
    assert "helixweave_home_poison" in config_call["home_poison"]
    environment = config_call["env"]
    assert environment["JAVA_HOME"] == str(java_home)
    assert environment["NXF_JAVA_HOME"] == str(java_home)
    assert environment["JAVA_CMD"] == str(fixture.java_executable)
    assert environment["PATH"] == f"{java_home / 'bin'}:/usr/bin:/bin"
    assert environment["LD_LIBRARY_PATH"] == ""
    assert "JAVA_TOOL_OPTIONS" not in environment
    assert "JDK_JAVA_OPTIONS" not in environment
    assert "_JAVA_OPTIONS" not in environment
    assert config_call["shell"] is False


def test_admission_hashes_heavy_assets_once_and_reuses_verified_object(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    heavy_calls = 0
    docker_calls = 0
    original_verify = runtime_assets_module.verify_runtime_assets

    def counted_verify(*args, **kwargs):
        nonlocal heavy_calls
        heavy_calls += 1
        return original_verify(*args, **kwargs)

    def counted_docker(
        binding: RuntimeAssetBinding,
        images: tuple[str, ...],
    ) -> bytes:
        nonlocal docker_calls
        docker_calls += 1
        return fixture.docker_probe(binding, images)

    monkeypatch.setattr(runtime_assets_module, "verify_runtime_assets", counted_verify)
    admission = RuntimeAssetAdmission(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=counted_docker,
        _runtime_probe=lambda *_args: None,
    )

    first = admission.acquire()
    hash_calls = 0
    original_hash = runtime_assets_module._hash_regular_entry

    def counted_hash(*args, **kwargs):
        nonlocal hash_calls
        hash_calls += 1
        return original_hash(*args, **kwargs)

    monkeypatch.setattr(runtime_assets_module, "_hash_regular_entry", counted_hash)
    second = admission.acquire()
    third = admission.acquire()

    assert first.is_success and second.is_success and third.is_success
    assert first.value is second.value is third.value
    assert heavy_calls == 1
    assert docker_calls == 3
    assert hash_calls == 0


def test_admission_metadata_change_invalidates_and_readmits(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    heavy_calls = 0
    original_verify = runtime_assets_module.verify_runtime_assets

    def counted_verify(*args, **kwargs):
        nonlocal heavy_calls
        heavy_calls += 1
        return original_verify(*args, **kwargs)

    monkeypatch.setattr(runtime_assets_module, "verify_runtime_assets", counted_verify)
    admission = RuntimeAssetAdmission(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
        _runtime_probe=lambda *_args: None,
    )

    first = admission.acquire()
    current = fixture.nextflow.stat()
    os.utime(
        fixture.nextflow,
        ns=(current.st_atime_ns, current.st_mtime_ns + 1_000_000),
    )
    second = admission.acquire()

    assert first.is_success and second.is_success
    assert first.value is not second.value
    assert heavy_calls == 2


def test_admission_never_uses_stale_evidence_after_asset_tamper(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    heavy_calls = 0
    original_verify = runtime_assets_module.verify_runtime_assets

    def counted_verify(*args, **kwargs):
        nonlocal heavy_calls
        heavy_calls += 1
        return original_verify(*args, **kwargs)

    monkeypatch.setattr(runtime_assets_module, "verify_runtime_assets", counted_verify)
    admission = RuntimeAssetAdmission(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
        _runtime_probe=lambda *_args: None,
    )
    assert admission.acquire().is_success
    content = bytearray(fixture.container_asset.read_bytes())
    content[-1] ^= 1
    fixture.container_asset.write_bytes(content)
    changed_stat = fixture.container_asset.stat()
    os.utime(
        fixture.container_asset,
        ns=(changed_stat.st_atime_ns, changed_stat.st_mtime_ns + 1_000_000),
    )

    changed = admission.acquire()

    assert changed.is_failure
    assert changed.value is None
    assert heavy_calls == 2


def test_admission_same_size_replacement_and_pid_change_force_readmission(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    heavy_calls = 0
    original_verify = runtime_assets_module.verify_runtime_assets

    def counted_verify(*args, **kwargs):
        nonlocal heavy_calls
        heavy_calls += 1
        return original_verify(*args, **kwargs)

    monkeypatch.setattr(runtime_assets_module, "verify_runtime_assets", counted_verify)
    admission = RuntimeAssetAdmission(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
        _runtime_probe=lambda *_args: None,
    )
    assert admission.acquire().is_success
    replacement = fixture.nextflow.with_name("replacement")
    replacement.write_bytes(fixture.nextflow.read_bytes())
    replacement.chmod(0o755)
    replacement.replace(fixture.nextflow)
    assert admission.acquire().is_success
    initial_pid = os.getpid()
    monkeypatch.setattr(runtime_assets_module.os, "getpid", lambda: initial_pid + 1)
    assert admission.acquire().is_success

    assert heavy_calls == 3


def test_admission_rechecks_docker_liveness_without_rehashing(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    heavy_calls = 0
    docker_calls = 0
    original_verify = runtime_assets_module.verify_runtime_assets

    def counted_verify(*args, **kwargs):
        nonlocal heavy_calls
        heavy_calls += 1
        return original_verify(*args, **kwargs)

    def changing_docker(
        binding: RuntimeAssetBinding,
        images: tuple[str, ...],
    ) -> bytes:
        nonlocal docker_calls
        docker_calls += 1
        if docker_calls == 1:
            return fixture.docker_probe(binding, images)
        return b"[]"

    monkeypatch.setattr(runtime_assets_module, "verify_runtime_assets", counted_verify)
    admission = RuntimeAssetAdmission(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=changing_docker,
        _runtime_probe=lambda *_args: None,
    )

    assert admission.acquire().is_success
    second = admission.acquire()

    assert second.is_failure
    assert heavy_calls == 1
    assert docker_calls == 2
    assert second.issues[0].context == {"component": "docker"}


def test_docker_endpoint_inode_change_invalidates_admission(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    docker = tmp_path / "server-tools/docker"
    docker.parent.mkdir()
    docker.write_bytes(b"synthetic docker client")
    docker.chmod(0o755)
    socket_path = tmp_path / "server-run/docker.sock"
    socket_path.parent.mkdir()
    binding = RuntimeAssetBinding(
        root=fixture.root,
        docker_executable=docker,
        docker_socket=socket_path,
    )
    heavy_calls = 0
    original_verify = runtime_assets_module.verify_runtime_assets

    def counted_verify(*args, **kwargs):
        nonlocal heavy_calls
        heavy_calls += 1
        return original_verify(*args, **kwargs)

    monkeypatch.setattr(runtime_assets_module, "verify_runtime_assets", counted_verify)
    monkeypatch.setattr(
        runtime_assets_module,
        "_verify_docker_availability",
        lambda *_args, **_kwargs: None,
    )
    real_lstat = runtime_assets_module.os.lstat
    endpoint_inode = 12

    def lstat(path):
        if Path(path) == socket_path:
            return SimpleNamespace(
                st_dev=11,
                st_ino=endpoint_inode,
                st_mode=stat.S_IFSOCK | 0o660,
                st_nlink=1,
                st_uid=os.geteuid(),
                st_gid=os.getegid(),
                st_size=0,
                st_mtime_ns=1,
                st_ctime_ns=1,
            )
        return real_lstat(path)

    monkeypatch.setattr(runtime_assets_module.os, "lstat", lstat)
    admission = RuntimeAssetAdmission(
        binding,
        _contract=fixture.contract,
        _runtime_probe=lambda *_args: None,
    )
    assert admission.acquire().is_success
    endpoint_inode = 13
    assert admission.acquire().is_success

    assert heavy_calls == 2


def test_fresh_worker_admission_does_not_trust_another_process_object(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    heavy_calls = 0
    original_verify = runtime_assets_module.verify_runtime_assets

    def counted_verify(*args, **kwargs):
        nonlocal heavy_calls
        heavy_calls += 1
        return original_verify(*args, **kwargs)

    monkeypatch.setattr(runtime_assets_module, "verify_runtime_assets", counted_verify)
    first = RuntimeAssetAdmission(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
        _runtime_probe=lambda *_args: None,
    )
    worker = RuntimeAssetAdmission(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
        _runtime_probe=lambda *_args: None,
    )

    assert first.acquire().is_success
    assert worker.acquire().is_success
    assert heavy_calls == 2


@pytest.mark.parametrize("case", ["missing", "extra"])
def test_container_process_set_is_exact(
    tmp_path: Path,
    case: str,
) -> None:
    fixture = _tiny_assets(tmp_path)
    entries = fixture.lock_payload["entries"]
    assert isinstance(entries, list)
    if case == "missing":
        entries.clear()
    else:
        extra = deepcopy(entries[0])
        extra["process"] = "FASTQC"
        entries.append(extra)
    _write_json(fixture.container_lock, fixture.lock_payload)

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert {issue.code for issue in result.issues} == {
        "BULK_RNASEQ_RUNTIME_ASSET_PROCESS_SET"
    }


def test_shared_image_coordinate_verifies_one_canonical_archive_closure(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    _add_process_sharing_container_coordinate(fixture)
    archive_calls = 0
    original = runtime_assets_module._verify_docker_archive

    def counted(*args, **kwargs):
        nonlocal archive_calls
        archive_calls += 1
        return original(*args, **kwargs)

    monkeypatch.setattr(runtime_assets_module, "_verify_docker_archive", counted)

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_success
    assert tuple(item.process for item in result.value.containers) == (
        "STAR_ALIGN",
        "ZZZ_ALIGN",
    )
    assert archive_calls == 1


@pytest.mark.parametrize("divergence", ["digest", "archive_path", "manifest_path"])
def test_shared_image_coordinate_rejects_divergent_operator_closure(
    tmp_path: Path,
    divergence: str,
) -> None:
    fixture = _tiny_assets(tmp_path)
    _add_process_sharing_container_coordinate(fixture)
    entries = fixture.lock_payload["entries"]
    assert isinstance(entries, list)
    shared = next(entry for entry in entries if entry["process"] == "ZZZ_ALIGN")
    if divergence == "digest":
        shared["oci_digest"] = "sha256:" + "0" * 64
        coordinate = fixture.contract.container_inventory["entries"][0][
            "image_coordinate"
        ]
        shared["image"] = f"{coordinate}@{shared['oci_digest']}"
    elif divergence == "archive_path":
        shared["local_asset"] = "star-copy.tar"
    else:
        shared["distribution_manifest_asset"] = "star-copy.manifest.json"
    _write_json(fixture.container_lock, fixture.lock_payload)

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_IDENTITY"


def test_container_image_coordinate_and_digest_are_both_fixed(tmp_path: Path) -> None:
    fixture = _tiny_assets(tmp_path)
    entries = fixture.lock_payload["entries"]
    assert isinstance(entries, list)
    entries[0]["image"] = f"quay.io/helixweave/other:test@{entries[0]['oci_digest']}"
    _write_json(fixture.container_lock, fixture.lock_payload)

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_IDENTITY"


@pytest.mark.parametrize("mutation", ["tamper", "symlink", "fifo"])
def test_container_archive_rejects_tamper_and_non_regular_files(
    tmp_path: Path,
    mutation: str,
) -> None:
    fixture = _tiny_assets(tmp_path)
    fixture.container_asset.unlink()
    if mutation == "tamper":
        fixture.container_asset.write_bytes(b"tampered")
    elif mutation == "symlink":
        target = fixture.root / "outside.tar"
        target.write_bytes(_docker_image_fixture().archive)
        fixture.container_asset.symlink_to(target)
    else:
        os.mkfifo(fixture.container_asset)

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code in {
        "BULK_RNASEQ_RUNTIME_ASSET_FILE_TYPE",
        "BULK_RNASEQ_RUNTIME_ASSET_IDENTITY",
    }


def test_distribution_manifest_bytes_must_equal_declared_oci_digest(
    tmp_path: Path,
) -> None:
    fixture = _tiny_assets(tmp_path)
    fixture.distribution_manifest.write_bytes(
        fixture.distribution_manifest.read_bytes() + b"\n"
    )

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_IDENTITY"


def test_archive_config_blob_must_match_distribution_manifest(
    tmp_path: Path,
) -> None:
    fixture = _tiny_assets(tmp_path)
    replacement = _docker_image_fixture(
        layer_tar=_tar_bytes({"different.txt": b"different\n"})
    )
    _replace_archive(fixture, replacement.archive)

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_IDENTITY"


def test_archive_layer_hash_must_match_config_rootfs_diff_id(tmp_path: Path) -> None:
    fixture = _tiny_assets(tmp_path)
    original = _docker_image_fixture()
    replacement = _docker_image_fixture(
        layer_tar=_tar_bytes({"different.txt": b"different\n"}),
        diff_id=f"sha256:{_sha256(original.layer_tar)}",
    )
    assert replacement.config_digest == fixture.config_digest
    _replace_archive(fixture, replacement.archive)

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_IDENTITY"


@pytest.mark.parametrize("kind", ["traversal", "symlink", "device", "duplicate"])
def test_docker_archive_rejects_unsafe_tar_members(
    tmp_path: Path,
    kind: str,
) -> None:
    fixture = _tiny_assets(tmp_path)
    _replace_archive(fixture, _unsafe_tar(kind))

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code in {
        "BULK_RNASEQ_RUNTIME_ASSET_CONTRACT",
        "BULK_RNASEQ_RUNTIME_ASSET_FILE_TYPE",
    }


def test_docker_archive_rejects_multiple_images(tmp_path: Path) -> None:
    fixture = _tiny_assets(tmp_path)
    files = _tar_file_contents(fixture.container_asset.read_bytes())
    manifest = json.loads(files["manifest.json"])
    files["manifest.json"] = _json_bytes([manifest[0], manifest[0]])
    _replace_archive(fixture, _tar_bytes(files))

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_CONTRACT"


@pytest.mark.parametrize(
    ("section", "field", "value"),
    [
        ("network_policy", "container_pull", True),
        ("network_policy", "wave", True),
        ("network_policy", "tower", True),
        ("execution_policy", "nxf_offline", False),
        ("execution_policy", "docker_pull", "missing"),
        ("execution_policy", "docker_network", "default"),
    ],
)
def test_container_lock_cannot_enable_network_or_pulls(
    tmp_path: Path,
    section: str,
    field: str,
    value: object,
) -> None:
    fixture = _tiny_assets(tmp_path)
    policy = fixture.lock_payload[section]
    assert isinstance(policy, dict)
    policy[field] = value
    _write_json(fixture.container_lock, fixture.lock_payload)

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_CONTRACT"


@pytest.mark.parametrize(
    "case",
    ["missing", "wrong_config_digest", "wrong_rootfs", "no_repo_digests", "exact"],
)
def test_doctor_requires_exact_image_in_local_daemon(
    tmp_path: Path,
    case: str,
) -> None:
    fixture = _tiny_assets(tmp_path)

    def probe(
        _binding: RuntimeAssetBinding,
        images: tuple[str, ...],
    ) -> bytes:
        if case == "missing":
            raise OSError("daemon unavailable")
        config_digest = fixture.config_digest
        if case == "wrong_config_digest":
            config_digest = f"sha256:{'f' * 64}"
        rootfs_diff_ids = fixture.rootfs_diff_ids
        if case == "wrong_rootfs":
            rootfs_diff_ids = (f"sha256:{'e' * 64}",)
        return _json_bytes(
            [
                {
                    "Id": config_digest,
                    "RepoDigests": [] if case == "no_repo_digests" else ["ignored"],
                    "RootFS": {
                        "Type": "layers",
                        "Layers": [*rootfs_diff_ids],
                    },
                }
                for _image in images
            ]
        )

    report = doctor_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=probe,
    )
    rendered = json.dumps([issue.to_dict() for issue in report.issues])

    assert report.ready is (case in {"exact", "no_repo_digests"})
    assert str(fixture.root) not in rendered
    assert fixture.config_digest not in rendered
    if case in {"wrong_config_digest", "wrong_rootfs"}:
        assert report.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_IDENTITY"
    if case == "missing":
        assert report.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_UNAVAILABLE"


def test_doctor_probes_a_shared_config_image_id_only_once(tmp_path: Path) -> None:
    fixture = _tiny_assets(tmp_path)
    verified = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )
    assert verified.value is not None
    first = verified.value.containers[0]
    containers = (first, replace(first, process="FASTQC"))
    observed: list[tuple[str, ...]] = []

    def probe(
        _binding: RuntimeAssetBinding,
        images: tuple[str, ...],
    ) -> bytes:
        observed.append(images)
        return _json_bytes(
            [
                {
                    "Id": fixture.config_digest,
                    "RepoDigests": [],
                    "RootFS": {
                        "Type": "layers",
                        "Layers": [*fixture.rootfs_diff_ids],
                    },
                }
            ]
        )

    _verify_docker_availability(fixture.binding, containers, probe=probe)

    assert observed == [(fixture.config_digest,)]


def test_production_docker_probe_is_one_fixed_local_inspect_without_shell_or_fetch(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    docker = tmp_path / "server-tools/docker"
    docker.parent.mkdir()
    docker.write_bytes(b"synthetic docker client")
    docker.chmod(0o755)
    socket_path = tmp_path / "server-run/docker.sock"
    socket_path.parent.mkdir()
    binding = RuntimeAssetBinding(
        root=fixture.root,
        docker_executable=docker,
        docker_socket=socket_path,
    )
    calls: list[tuple[tuple[str, ...], dict[str, object]]] = []

    class FakePopen:
        def __init__(self, argv, **kwargs):
            read_fd, write_fd = os.pipe()
            output = _json_bytes(
                [
                    {
                        "Id": fixture.config_digest,
                        "RepoDigests": [],
                        "RootFS": {
                            "Type": "layers",
                            "Layers": [*fixture.rootfs_diff_ids],
                        },
                    }
                ]
            )
            os.write(write_fd, output)
            os.close(write_fd)
            self.stdout = os.fdopen(read_fd, "rb", buffering=0)
            self.returncode = 0
            calls.append((argv, kwargs))

        def wait(self, timeout=None):
            return self.returncode

        def poll(self):
            return self.returncode

        def kill(self):
            self.returncode = -9

    def popen(argv, **kwargs):
        return FakePopen(argv, **kwargs)

    monkeypatch.setattr(runtime_assets_module.subprocess, "Popen", popen)
    real_lstat = runtime_assets_module.os.lstat

    def lstat(path):
        if Path(path) == socket_path:
            return SimpleNamespace(
                st_dev=11,
                st_ino=12,
                st_mode=stat.S_IFSOCK | 0o660,
                st_nlink=1,
                st_uid=os.geteuid(),
                st_gid=os.getegid(),
                st_size=0,
                st_mtime_ns=1,
                st_ctime_ns=1,
            )
        return real_lstat(path)

    monkeypatch.setattr(runtime_assets_module.os, "lstat", lstat)
    result = verify_runtime_assets(binding, _contract=fixture.contract)

    assert result.is_success
    assert len(calls) == 1
    argv, kwargs = calls[0]
    assert argv == (
        str(docker),
        "--host",
        f"unix://{socket_path}",
        "image",
        "inspect",
        fixture.config_digest,
    )
    assert kwargs["shell"] is False
    assert kwargs["stdin"] is subprocess.DEVNULL
    assert kwargs["stdout"] is subprocess.PIPE
    assert kwargs["stderr"] is subprocess.DEVNULL
    assert kwargs["env"] == {"HOME": "/", "LANG": "C", "LC_ALL": "C"}
    assert not {"pull", "load", "login", "run"} & set(argv)
    assert not any(
        token.startswith(("http://", "https://", "tcp://")) for token in argv
    )


def test_docker_inspect_command_failure_is_doctor_not_ready_and_redacted(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    monkeypatch.setattr(
        runtime_assets_module,
        "_verify_local_docker_endpoint",
        lambda _binding: None,
    )

    class FailingPopen:
        def __init__(self, _argv, **_kwargs):
            read_fd, write_fd = os.pipe()
            os.close(write_fd)
            self.stdout = os.fdopen(read_fd, "rb", buffering=0)
            self.returncode = 1

        def wait(self, timeout=None):
            return self.returncode

        def poll(self):
            return self.returncode

        def kill(self):
            self.returncode = -9

    monkeypatch.setattr(runtime_assets_module.subprocess, "Popen", FailingPopen)

    report = doctor_runtime_assets(fixture.binding, _contract=fixture.contract)
    rendered = json.dumps([issue.to_dict() for issue in report.issues])

    assert not report.ready
    assert report.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_UNAVAILABLE"
    assert str(fixture.root) not in rendered
    assert fixture.config_digest not in rendered


def test_doctor_rejects_a_group_or_world_writable_docker_cli(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    docker = tmp_path / "server-tools/docker"
    docker.parent.mkdir()
    docker.write_bytes(b"synthetic docker client")
    docker.chmod(0o777)
    socket_path = tmp_path / "server-run/docker.sock"
    binding = RuntimeAssetBinding(
        root=fixture.root,
        docker_executable=docker,
        docker_socket=socket_path,
    )
    real_lstat = runtime_assets_module.os.lstat

    def lstat(path):
        if Path(path) == socket_path:
            return SimpleNamespace(
                st_dev=11,
                st_ino=12,
                st_mode=stat.S_IFSOCK | 0o660,
                st_nlink=1,
                st_uid=os.geteuid(),
                st_gid=os.getegid(),
                st_size=0,
                st_mtime_ns=1,
                st_ctime_ns=1,
            )
        return real_lstat(path)

    monkeypatch.setattr(runtime_assets_module.os, "lstat", lstat)
    result = verify_runtime_assets(binding, _contract=fixture.contract)

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_FILE_TYPE"


def test_doctor_rechecks_docker_endpoint_identity_after_inspect(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    fixture = _tiny_assets(tmp_path)
    docker = tmp_path / "server-tools/docker"
    docker.parent.mkdir()
    docker.write_bytes(b"synthetic docker client")
    docker.chmod(0o755)
    socket_path = tmp_path / "server-run/docker.sock"
    binding = RuntimeAssetBinding(
        root=fixture.root,
        docker_executable=docker,
        docker_socket=socket_path,
    )
    output = _json_bytes(
        [
            {
                "Id": fixture.config_digest,
                "RootFS": {
                    "Type": "layers",
                    "Layers": [*fixture.rootfs_diff_ids],
                },
            }
        ]
    )

    class FakePopen:
        def __init__(self, _argv, **_kwargs):
            read_fd, write_fd = os.pipe()
            os.write(write_fd, output)
            os.close(write_fd)
            self.stdout = os.fdopen(read_fd, "rb", buffering=0)
            self.returncode = 0

        def wait(self, timeout=None):
            return self.returncode

        def poll(self):
            return self.returncode

        def kill(self):
            self.returncode = -9

    monkeypatch.setattr(runtime_assets_module.subprocess, "Popen", FakePopen)
    real_lstat = runtime_assets_module.os.lstat
    socket_reads = 0

    def lstat(path):
        nonlocal socket_reads
        if Path(path) == socket_path:
            socket_reads += 1
            return SimpleNamespace(
                st_dev=11,
                st_ino=11 + socket_reads,
                st_mode=stat.S_IFSOCK | 0o660,
                st_nlink=1,
                st_uid=os.geteuid(),
                st_gid=os.getegid(),
                st_size=0,
                st_mtime_ns=1,
                st_ctime_ns=1,
            )
        return real_lstat(path)

    monkeypatch.setattr(runtime_assets_module.os, "lstat", lstat)
    result = verify_runtime_assets(binding, _contract=fixture.contract)

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_RACE"
    assert socket_reads == 2


def test_source_tree_rejects_extra_and_symlink_entries(tmp_path: Path) -> None:
    extra_fixture = _tiny_assets(tmp_path / "extra")
    extra_fixture.source.joinpath("unexpected.nf").write_text("process X {}\n")
    extra = verify_runtime_assets(
        extra_fixture.binding,
        _contract=extra_fixture.contract,
        _docker_probe=extra_fixture.docker_probe,
    )

    symlink_fixture = _tiny_assets(tmp_path / "symlink")
    main = symlink_fixture.source / "main.nf"
    main.unlink()
    main.symlink_to(symlink_fixture.source / "conf/settings.config")
    symlink = verify_runtime_assets(
        symlink_fixture.binding,
        _contract=symlink_fixture.contract,
        _docker_probe=symlink_fixture.docker_probe,
    )

    assert extra.is_failure
    assert extra.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_FILE_SET"
    assert symlink.is_failure
    assert symlink.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_FILE_TYPE"


def test_runtime_root_rejects_symlinked_ancestor(tmp_path: Path) -> None:
    real_parent = tmp_path / "real-parent"
    fixture = _tiny_assets(real_parent / "runtime")
    linked_parent = tmp_path / "linked-parent"
    linked_parent.symlink_to(real_parent, target_is_directory=True)
    binding = replace(fixture.binding, root=linked_parent / "runtime")

    result = verify_runtime_assets(
        binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_FILE_TYPE"


def test_source_tree_rejects_reserved_default_deny_label_even_when_hashed(
    tmp_path: Path,
) -> None:
    fixture = _tiny_assets(tmp_path)
    content = (
        fixture.source.joinpath("main.nf").read_bytes()
        + f"\nlabel '{CONTAINER_RESERVED_DEFAULT_DENY_LABEL}'\n".encode()
    )
    fixture.source.joinpath("main.nf").write_bytes(content)
    files = fixture.contract.source_manifest["files"]
    assert isinstance(files, list)
    for entry in files:
        if entry["path"] == "main.nf":
            entry["size_bytes"] = len(content)
            entry["sha256"] = _sha256(content)
    fixture.contract.source_manifest["tree_sha256"] = _canonical_sha256(files)

    result = verify_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )

    assert result.is_failure
    assert result.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_RESERVED_LABEL"


def test_plugin_tree_and_nextflow_mode_fail_closed(tmp_path: Path) -> None:
    plugin_fixture = _tiny_assets(tmp_path / "plugin")
    plugin_fixture.plugin_tree.joinpath("unexpected.jar").write_bytes(b"extra")
    plugin = verify_runtime_assets(
        plugin_fixture.binding,
        _contract=plugin_fixture.contract,
        _docker_probe=plugin_fixture.docker_probe,
    )

    nextflow_fixture = _tiny_assets(tmp_path / "nextflow")
    nextflow_fixture.nextflow.chmod(0o644)
    nextflow = verify_runtime_assets(
        nextflow_fixture.binding,
        _contract=nextflow_fixture.contract,
        _docker_probe=nextflow_fixture.docker_probe,
    )

    assert plugin.is_failure
    assert plugin.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_FILE_SET"
    assert nextflow.is_failure
    assert nextflow.issues[0].code == "BULK_RNASEQ_RUNTIME_ASSET_MODE"


def test_doctor_issues_never_disclose_paths_or_digests(tmp_path: Path) -> None:
    fixture = _tiny_assets(tmp_path)
    fixture.nextflow.unlink()

    report = doctor_runtime_assets(
        fixture.binding,
        _contract=fixture.contract,
        _docker_probe=fixture.docker_probe,
    )
    rendered = json.dumps([issue.to_dict() for issue in report.issues])

    assert not report.ready
    assert str(fixture.root) not in rendered
    assert NEXTFLOW_SHA256 not in rendered
    assert fixture.contract.identity["nextflow"]["sha256"] not in rendered
    assert all(issue.path is None for issue in report.issues)
    assert all(issue.technical_message is None for issue in report.issues)


@pytest.mark.parametrize(
    "relative_path",
    ["../escape", "/absolute/path", "plugins/../escape", "windows\\path"],
)
def test_binding_rejects_paths_outside_asset_root(
    tmp_path: Path,
    relative_path: str,
) -> None:
    with pytest.raises(ValueError, match="relative path is invalid"):
        RuntimeAssetBinding(root=tmp_path, container_lock=relative_path)


@pytest.mark.parametrize("suffix", ["quoted'root", "dollar$root", "back\\slash", "a:b"])
def test_binding_rejects_runtime_roots_unsafe_for_hard_config_or_path(
    tmp_path: Path,
    suffix: str,
) -> None:
    with pytest.raises(ValueError, match="safe canonical path"):
        RuntimeAssetBinding(root=(tmp_path / suffix).resolve())


@pytest.mark.parametrize("source_tree", ["source/quoted'tree", "source/dollar$tree"])
def test_binding_rejects_source_path_unsafe_for_fixed_config(
    tmp_path: Path,
    source_tree: str,
) -> None:
    with pytest.raises(ValueError, match="fixed configuration"):
        RuntimeAssetBinding(root=tmp_path.resolve(), source_tree=source_tree)


def test_binding_rejects_jdk_path_with_path_separator(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="controlled PATH"):
        RuntimeAssetBinding(root=tmp_path.resolve(), jdk_tree="jdk/corretto:host")


@pytest.mark.parametrize("field", ["docker_executable", "docker_socket"])
def test_binding_requires_server_owned_absolute_docker_endpoint(
    tmp_path: Path,
    field: str,
) -> None:
    with pytest.raises(ValueError, match="absolute pathlib.Path"):
        RuntimeAssetBinding(root=tmp_path, **{field: Path("relative")})


def test_binding_requires_a_fixed_docker_cli_basename(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="fixed Docker CLI path"):
        RuntimeAssetBinding(
            root=tmp_path,
            docker_executable=Path("/usr/bin/docker-wrapper"),
        )


def test_binding_rejects_string_docker_endpoints(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="absolute pathlib.Path"):
        RuntimeAssetBinding(root=tmp_path, docker_executable="/usr/bin/docker")
