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
    NEXTFLOW_OFFLINE_ENV,
    NFCORE_RNASEQ_COMMIT,
    NF_SCHEMA_ARCHIVE_SHA256,
    NF_SCHEMA_VERSION,
    NEXTFLOW_SHA256,
    NEXTFLOW_VERSION,
    SOURCE_FILE_COUNT,
    SOURCE_TREE_SHA256,
    RuntimeAssetBinding,
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
            "size_bytes": len(nextflow_content),
            "sha256": _sha256(nextflow_content),
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
        plugin_tree=plugin_tree,
        container_lock=container_lock,
        container_asset=container_asset,
        distribution_manifest=distribution_manifest,
        config_digest=docker_image.config_digest,
        rootfs_diff_ids=docker_image.rootfs_diff_ids,
        docker_probe=exact_docker_probe,
        lock_payload=lock_payload,
    )


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
