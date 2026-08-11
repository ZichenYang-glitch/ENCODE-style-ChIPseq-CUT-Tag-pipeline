from __future__ import annotations

from dataclasses import replace
from hashlib import sha256
import hashlib
import os
from pathlib import Path
from types import SimpleNamespace

import pytest

import encode_pipeline.deployment.bundle_builder as builder_module
import encode_pipeline.deployment.native_contracts as native_module
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    RUNTIME_IDENTITY_SHA256,
    VerifiedRuntimeAssets,
)
from encode_pipeline.deployment.bundle import BundleStore
from encode_pipeline.deployment.bundle_builder import (
    build_bulk_rnaseq_runtime_bundle,
    build_encode_runtime_bundle,
    build_platform_bundle,
)
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    ENCODE_RUNTIME,
    PLATFORM,
)
from encode_pipeline.deployment.native_contracts import (
    ENCODE_MICROMAMBA_PATH,
    ENCODE_PACKAGE_ARCHIVE_ROOT,
    ENCODE_RUNTIME_INDEX_PATH,
    PLATFORM_BULK_REFERENCES_PATH,
    PLATFORM_ENCODE_REFERENCES_PATH,
    PLATFORM_FRONTEND_PATH,
    PLATFORM_METADATA_PATH,
    PLATFORM_MIGRATIONS_PATH,
    PLATFORM_WHEEL_PATH,
    ProductionNativeContractResolver,
)
from encode_pipeline.deployment.platform_runtime import (
    PLATFORM_RUNTIME_LAUNCHER_PATH,
    PLATFORM_RUNTIME_LOCK_PATH,
    PLATFORM_RUNTIME_SITE_PACKAGES,
    PLATFORM_RUNTIME_WHEELHOUSE_ROOT,
    LockedWheel,
    PlatformRuntimeClosure,
    PlatformRuntimeFile,
    PlatformWheelLock,
)
from encode_pipeline.platform.results import Issue, Result


def _platform_inputs(
    tmp_path: Path,
    monkeypatch,
    wheel: bytes,
) -> tuple[Path, Path, Path, Path]:
    filename = "helixweave-1.4.0-py3-none-any.whl"
    wheel_path = tmp_path / "helixweave.whl"
    wheel_path.write_bytes(wheel)
    wheelhouse = tmp_path / "wheelhouse"
    wheelhouse.mkdir()
    (wheelhouse / filename).write_bytes(wheel)
    lock = PlatformWheelLock((LockedWheel(filename, sha256(wheel).hexdigest()),))
    lock_path = tmp_path / "python-runtime-wheel-lock.json"
    lock_path.write_bytes(lock.to_bytes())
    scratch = tmp_path / "scratch"
    scratch.mkdir(mode=0o700)

    def collect(
        selected_wheelhouse: Path,
        selected_lock: Path,
        destination: Path,
    ) -> PlatformRuntimeClosure:
        assert selected_wheelhouse == wheelhouse
        assert selected_lock == lock_path
        files = {
            PLATFORM_RUNTIME_LOCK_PATH: lock.to_bytes(),
            f"{PLATFORM_RUNTIME_WHEELHOUSE_ROOT}/{filename}": wheel,
            f"{PLATFORM_RUNTIME_SITE_PACKAGES}/encode_pipeline/__init__.py": (
                b"expanded platform\n"
            ),
            PLATFORM_RUNTIME_LAUNCHER_PATH: b"#!/bin/sh\n",
        }
        records: list[PlatformRuntimeFile] = []
        for logical_path, content in files.items():
            path = destination.joinpath(*Path(logical_path).parts)
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(content)
            mode = 0o555 if logical_path == PLATFORM_RUNTIME_LAUNCHER_PATH else 0o444
            path.chmod(mode)
            records.append(
                PlatformRuntimeFile(
                    logical_path,
                    path,
                    mode,
                    len(content),
                    sha256(content).hexdigest(),
                )
            )
        return PlatformRuntimeClosure(lock.identity, tuple(sorted(records)))

    monkeypatch.setattr(builder_module, "collect_platform_runtime_closure", collect)
    return wheel_path, wheelhouse, lock_path, scratch


def _static_elf() -> bytes:
    content = bytearray(120)
    content[:7] = b"\x7fELF\x02\x01\x01"
    content[18:20] = (62).to_bytes(2, "little")
    content[32:40] = (64).to_bytes(8, "little")
    content[54:56] = (56).to_bytes(2, "little")
    content[56:58] = (1).to_bytes(2, "little")
    content[64:68] = (1).to_bytes(4, "little")
    return bytes(content)


def _bulk_inputs(tmp_path: Path) -> tuple[Path, VerifiedRuntimeAssets]:
    root = tmp_path / "bulk-runtime"
    files = {
        "source/rnaseq/main.nf": b"workflow rnaseq {}\n",
        "nextflow/nextflow": b"#!/bin/sh\n",
        "jdk/runtime.tar.gz": b"jdk archive",
        "plugins/plugin.zip": b"plugin archive",
        "containers/image.tar": b"container archive",
    }
    for relative, content in files.items():
        path = root.joinpath(*Path(relative).parts)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
        path.chmod(0o555 if relative == "nextflow/nextflow" else 0o444)
    assets = VerifiedRuntimeAssets(
        root=root,
        source_tree=root / "source/rnaseq",
        nextflow_executable=root / "nextflow/nextflow",
        jdk_archive=root / "jdk/runtime.tar.gz",
        jdk_tree=root / "jdk/runtime",
        java_executable=root / "jdk/runtime/bin/java",
        plugin_root=root / "plugins",
        plugin_archive=root / "plugins/plugin.zip",
        plugin_meta=root / "plugins/plugin.json",
        plugin_tree=root / "plugins/plugin",
        container_lock=root / "containers/availability-lock.json",
        containers=(),
        source_tree_sha256="1" * 64,
        runtime_identity_sha256=RUNTIME_IDENTITY_SHA256,
        nextflow_sha256="2" * 64,
        jdk_archive_sha256="3" * 64,
        jdk_tree_sha256="4" * 64,
        java_executable_sha256="5" * 64,
        plugin_archive_sha256="6" * 64,
        plugin_tree_sha256="7" * 64,
        container_inventory_sha256="8" * 64,
        container_lock_sha256="9" * 64,
    )
    return root, assets


def _explicit_lock(
    filename: str,
    content: bytes,
    *,
    platform: str,
    channel: str = "conda-forge",
) -> bytes:
    digest = hashlib.md5(content, usedforsecurity=False).hexdigest()
    return (
        "# Generated by conda-lock.\n"
        f"# platform: {platform}\n"
        "@EXPLICIT\n"
        f"https://conda.anaconda.org/{channel}/{platform}/{filename}#{digest}\n"
    ).encode()


def _encode_inputs(
    tmp_path: Path,
    monkeypatch,
    *,
    lock_platform: str = "linux-64",
    basename_collision: bool = False,
    runner_filename: str = "runner-1.conda",
) -> tuple[Path, Path, Path, tuple[tuple[str, bytes], ...]]:
    project = tmp_path / "project"
    inventory = (
        Path(__file__).resolve().parents[2]
        / "docs/architecture/artifact-inventory.yaml"
    )
    catalog = project / "docs/architecture/artifact-inventory.yaml"
    catalog.parent.mkdir(parents=True)
    catalog.write_bytes(inventory.read_bytes())
    runner_archive = b"runner archive"
    tool_archive = runner_archive if basename_collision else b"tool archive"
    tool_filename = "runner-1.conda" if basename_collision else "tool-1.conda"
    source_manifest = (
        (
            "docs/architecture/artifact-inventory.yaml",
            inventory.read_bytes(),
        ),
        ("workflow/Snakefile", b'include: "rules/test.smk"\n'),
        (
            "workflow/rules/test.smk",
            b'rule test:\n    conda:\n        "../envs/tool.yml"\n',
        ),
        ("workflow/envs/runner.yml", b"name: runner\n"),
        (
            "workflow/envs/runner.lock",
            _explicit_lock(runner_filename, runner_archive, platform=lock_platform),
        ),
        ("workflow/envs/tool.yml", b"name: tool\n"),
        (
            "workflow/envs/tool.lock",
            _explicit_lock(
                tool_filename,
                tool_archive,
                platform=lock_platform,
                channel="bioconda" if basename_collision else "conda-forge",
            ),
        ),
    )
    for logical_path, content in source_manifest:
        path = project.joinpath(*Path(logical_path).parts)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(content)
    micromamba = tmp_path / "micromamba"
    micromamba.write_bytes(_static_elf())
    micromamba.chmod(0o755)
    cache = tmp_path / "archives"
    cache.mkdir()
    (cache / runner_filename).write_bytes(runner_archive)
    if not basename_collision:
        (cache / "tool-1.conda").write_bytes(tool_archive)

    class Provider:
        def __init__(self, _registry, *, project_root):
            assert project_root.is_absolute()

        def source_manifest(self):
            return source_manifest

        def capture(self, _workflow_id):
            return SimpleNamespace(
                is_failure=False,
                value=SimpleNamespace(digest="d" * 64, adapter_version="1.0.0"),
            )

    monkeypatch.setattr(builder_module, "WorkflowBuildIdentityProvider", Provider)
    monkeypatch.setattr(native_module, "WorkflowBuildIdentityProvider", Provider)
    return project, micromamba, cache, source_manifest


def test_encode_builder_accepts_a_real_conda_leading_underscore_basename(
    tmp_path: Path,
    monkeypatch,
) -> None:
    project, micromamba, cache, _source = _encode_inputs(
        tmp_path,
        monkeypatch,
        runner_filename="_openmp_mutex-4.5-20_gnu.conda",
    )

    manifest = build_encode_runtime_bundle(
        project,
        micromamba,
        cache,
        tmp_path / "encode-underscore.tar",
    )

    assert any(
        record.path.endswith("/_openmp_mutex-4.5-20_gnu.conda")
        for record in manifest.files
    )


def test_bulk_builder_is_canonical_and_indexes_only_verified_runtime_bytes(
    tmp_path: Path,
) -> None:
    root, assets = _bulk_inputs(tmp_path)
    observed_roots: list[Path] = []

    def verify(binding):
        observed_roots.append(binding.root)
        return Result.success(assets)

    first = tmp_path / "bulk-a.tar"
    second = tmp_path / "bulk-b.tar"
    manifest_a = build_bulk_rnaseq_runtime_bundle(
        root,
        first,
        runtime_verifier=verify,
    )
    manifest_b = build_bulk_rnaseq_runtime_bundle(
        root,
        second,
        runtime_verifier=verify,
    )

    assert manifest_a == manifest_b
    assert first.read_bytes() == second.read_bytes()
    assert manifest_a.component == BULK_RNASEQ_RUNTIME
    assert observed_roots == [root, root]
    assert {record.path for record in manifest_a.files} == {
        native_module.BULK_RUNTIME_IDENTITY_PATH,
        "payload/runtime/source/rnaseq/main.nf",
        "payload/runtime/nextflow/nextflow",
        "payload/runtime/jdk/runtime.tar.gz",
        "payload/runtime/plugins/plugin.zip",
        "payload/runtime/containers/image.tar",
    }
    assert (
        BundleStore(DeploymentLayout.isolated(tmp_path / "bulk-host")).inspect(first)
        == manifest_a
    )


@pytest.mark.parametrize("fault", ("exception", "failure", "wrong-identity"))
def test_bulk_builder_rejects_unverified_runtime_closure(
    tmp_path: Path,
    fault: str,
) -> None:
    root, assets = _bulk_inputs(tmp_path)

    def verify(_binding):
        if fault == "exception":
            raise RuntimeError(str(tmp_path))
        if fault == "failure":
            return Result.failure((Issue("RUNTIME_INVALID", "Runtime is invalid."),))
        return Result.success(replace(assets, runtime_identity_sha256="f" * 64))

    output = tmp_path / f"bulk-{fault}.tar"
    with pytest.raises(DeploymentError) as caught:
        build_bulk_rnaseq_runtime_bundle(root, output, runtime_verifier=verify)

    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_SOURCE_INVALID"
    assert str(tmp_path) not in str(caught.value)
    assert not output.exists()


@pytest.mark.parametrize("fault", ("extra-root", "symlink-root", "source-drift"))
def test_bulk_builder_rejects_unverified_tree_changes(
    tmp_path: Path,
    fault: str,
) -> None:
    root, assets = _bulk_inputs(tmp_path)
    if fault == "extra-root":
        (root / "unexpected").mkdir()
    elif fault == "symlink-root":
        source = root / "source"
        relocated = tmp_path / "source-real"
        source.rename(relocated)
        source.symlink_to(relocated, target_is_directory=True)

    def verify(_binding):
        if fault == "source-drift":
            (root / "containers/image.tar").chmod(0o666)
        return Result.success(assets)

    output = tmp_path / f"bulk-tree-{fault}.tar"
    with pytest.raises(DeploymentError) as caught:
        build_bulk_rnaseq_runtime_bundle(root, output, runtime_verifier=verify)

    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_SOURCE_INVALID"
    assert not output.exists()


def test_bulk_builder_never_overwrites_an_existing_output(
    tmp_path: Path,
) -> None:
    root, assets = _bulk_inputs(tmp_path)
    output = tmp_path / "bulk.tar"
    build_bulk_rnaseq_runtime_bundle(
        root,
        output,
        runtime_verifier=lambda _binding: Result.success(assets),
    )
    original = output.read_bytes()

    with pytest.raises(DeploymentError) as caught:
        build_bulk_rnaseq_runtime_bundle(
            root,
            output,
            runtime_verifier=lambda _binding: Result.success(assets),
        )

    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_OUTPUT_EXISTS"
    assert output.read_bytes() == original
    assert not tuple(tmp_path.glob("*.partial"))


def test_platform_builder_indexes_only_authoritative_candidate_bytes(
    tmp_path: Path,
    monkeypatch,
) -> None:
    wheel = b"already-built-wheel"
    wheel_path, wheelhouse, lock_path, scratch = _platform_inputs(
        tmp_path,
        monkeypatch,
        wheel,
    )
    facts = SimpleNamespace(
        metadata=b"Name: helixweave\nVersion: 1.4.0\n",
        frontend_manifest=b'{"frontend":"native"}\n',
        migration_inventory=b'{"migration":"native"}\n',
        encode_reference_source=b'ENCODE_REFERENCE_BINDING_CONTRACT = "1.0.0"\n',
        bulk_reference_source=b'BULK_RNASEQ_REFERENCE_BINDING_CONTRACT = "1.0.0"\n',
        frontend_identity=f"sha256-{'1' * 64}",
        reference_identity=f"sha256-{'2' * 64}",
    )
    monkeypatch.setattr(
        builder_module,
        "verify_platform_wheel",
        lambda content: facts if content == wheel else None,
    )

    output_a = tmp_path / "platform-a.tar"
    output_b = tmp_path / "platform-b.tar"
    manifest_a = build_platform_bundle(
        wheel_path,
        wheelhouse,
        lock_path,
        output_a,
        scratch_root=scratch,
    )
    manifest_b = build_platform_bundle(
        wheel_path,
        wheelhouse,
        lock_path,
        output_b,
        scratch_root=scratch,
    )

    assert manifest_a == manifest_b
    assert output_a.read_bytes() == output_b.read_bytes()
    assert manifest_a.component == PLATFORM
    assert {record.path for record in manifest_a.files} == {
        PLATFORM_METADATA_PATH,
        PLATFORM_WHEEL_PATH,
        PLATFORM_FRONTEND_PATH,
        PLATFORM_MIGRATIONS_PATH,
        PLATFORM_ENCODE_REFERENCES_PATH,
        PLATFORM_BULK_REFERENCES_PATH,
        PLATFORM_RUNTIME_LOCK_PATH,
        f"{PLATFORM_RUNTIME_WHEELHOUSE_ROOT}/helixweave-1.4.0-py3-none-any.whl",
        f"{PLATFORM_RUNTIME_SITE_PACKAGES}/encode_pipeline/__init__.py",
        PLATFORM_RUNTIME_LAUNCHER_PATH,
    }
    assert all("version.json" not in record.path for record in manifest_a.files)
    assert (
        BundleStore(DeploymentLayout.isolated(tmp_path / "host")).inspect(output_a)
        == manifest_a
    )


def test_encode_builder_is_canonical_and_reuses_native_build_identity(
    tmp_path: Path,
    monkeypatch,
) -> None:
    project_root, micromamba, cache, source_manifest = _encode_inputs(
        tmp_path,
        monkeypatch,
    )
    output_a = tmp_path / "encode-a.tar"
    output_b = tmp_path / "encode-b.tar"

    manifest_a = build_encode_runtime_bundle(
        project_root,
        micromamba,
        cache,
        output_a,
    )
    manifest_b = build_encode_runtime_bundle(
        project_root,
        micromamba,
        cache,
        output_b,
    )

    assert manifest_a == manifest_b
    assert (
        sha256(output_a.read_bytes()).digest() == sha256(output_b.read_bytes()).digest()
    )
    assert manifest_a.component == ENCODE_RUNTIME
    paths = {record.path for record in manifest_a.files}
    assert ENCODE_RUNTIME_INDEX_PATH in paths
    assert ENCODE_MICROMAMBA_PATH in paths
    assert {
        f"payload/contracts/encode-runtime/{logical_path}"
        for logical_path, _content in source_manifest
    } < paths
    assert {
        record.path
        for record in manifest_a.files
        if record.path.startswith(f"{ENCODE_PACKAGE_ARCHIVE_ROOT}/")
    }
    assert all("/config/" not in record.path for record in manifest_a.files)
    layout = DeploymentLayout.isolated(tmp_path / "installed")
    store = BundleStore(layout)
    staged = store.stage(
        output_a,
        installed_owner_uid=os.getuid(),
        installed_owner_gid=os.getgid(),
    )
    facts = ProductionNativeContractResolver().resolve(
        layout.component_store(ENCODE_RUNTIME) / staged.identity,
        staged,
    )

    assert facts.component == ENCODE_RUNTIME
    assert facts.deployment_identity == manifest_a.identity
    assert facts.contracts == manifest_a.provides
    assert facts.database_heads == ()


@pytest.mark.parametrize("fault", ("missing", "extra", "symlink"))
def test_encode_builder_refuses_an_inexact_or_unsafe_archive_cache(
    tmp_path: Path,
    monkeypatch,
    fault: str,
) -> None:
    project, micromamba, cache, _source = _encode_inputs(tmp_path, monkeypatch)
    if fault == "missing":
        (cache / "tool-1.conda").unlink()
    elif fault == "extra":
        (cache / "unrequested-1.conda").write_bytes(b"extra")
    else:
        archive = cache / "tool-1.conda"
        archive.unlink()
        archive.symlink_to(cache / "runner-1.conda")

    with pytest.raises(DeploymentError) as caught:
        build_encode_runtime_bundle(
            project,
            micromamba,
            cache,
            tmp_path / f"{fault}.tar",
        )

    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_SOURCE_INVALID"


def test_encode_builder_refuses_a_non_linux_explicit_lock(
    tmp_path: Path,
    monkeypatch,
) -> None:
    project, micromamba, cache, _source = _encode_inputs(
        tmp_path,
        monkeypatch,
        lock_platform="osx-64",
    )

    with pytest.raises(DeploymentError) as caught:
        build_encode_runtime_bundle(
            project,
            micromamba,
            cache,
            tmp_path / "osx.tar",
        )

    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_SOURCE_INVALID"


def test_encode_builder_refuses_ambiguous_archive_basenames(
    tmp_path: Path,
    monkeypatch,
) -> None:
    project, micromamba, cache, _source = _encode_inputs(
        tmp_path,
        monkeypatch,
        basename_collision=True,
    )

    with pytest.raises(DeploymentError) as caught:
        build_encode_runtime_bundle(
            project,
            micromamba,
            cache,
            tmp_path / "collision.tar",
        )

    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_SOURCE_INVALID"


def test_builder_refuses_a_world_writable_output_boundary(
    tmp_path: Path,
    monkeypatch,
) -> None:
    wheel_path, wheelhouse, lock_path, scratch = _platform_inputs(
        tmp_path,
        monkeypatch,
        b"wheel",
    )
    facts = SimpleNamespace(
        metadata=b"metadata",
        frontend_manifest=b"frontend",
        migration_inventory=b"migrations",
        encode_reference_source=b"encode",
        bulk_reference_source=b"bulk",
        frontend_identity=f"sha256-{'1' * 64}",
        reference_identity=f"sha256-{'2' * 64}",
    )
    monkeypatch.setattr(builder_module, "verify_platform_wheel", lambda _raw: facts)
    unsafe = tmp_path / "unsafe"
    unsafe.mkdir()
    unsafe.chmod(0o777)

    with pytest.raises(DeploymentError) as caught:
        build_platform_bundle(
            wheel_path,
            wheelhouse,
            lock_path,
            unsafe / "platform.tar",
            scratch_root=scratch,
        )

    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_OUTPUT_INVALID"
    assert not tuple(unsafe.iterdir())


@pytest.mark.parametrize(
    ("limit_name", "limit_value"),
    (
        ("MAX_BUNDLE_FILES", 1),
        ("MAX_MANIFEST_BYTES", 1),
        ("MAX_BUNDLE_BYTES", 1),
    ),
)
def test_builder_never_writes_a_bundle_the_consumer_bounds_reject(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
    limit_name: str,
    limit_value: int,
) -> None:
    wheel_path, wheelhouse, lock_path, scratch = _platform_inputs(
        tmp_path,
        monkeypatch,
        b"wheel",
    )
    facts = SimpleNamespace(
        metadata=b"metadata",
        frontend_manifest=b"frontend",
        migration_inventory=b"migrations",
        encode_reference_source=b"encode",
        bulk_reference_source=b"bulk",
        frontend_identity=f"sha256-{'1' * 64}",
        reference_identity=f"sha256-{'2' * 64}",
    )
    monkeypatch.setattr(builder_module, "verify_platform_wheel", lambda _raw: facts)
    monkeypatch.setattr(builder_module, limit_name, limit_value)
    monkeypatch.setattr(
        builder_module,
        "_write_canonical_tar",
        lambda *_args, **_kwargs: pytest.fail("tar writer reached"),
    )
    output = tmp_path / "oversized.tar"

    with pytest.raises(DeploymentError) as caught:
        build_platform_bundle(
            wheel_path,
            wheelhouse,
            lock_path,
            output,
            scratch_root=scratch,
        )

    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_SOURCE_INVALID"
    assert not output.exists()


def test_builder_preflights_output_capacity_before_tar_creation(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    wheel_path, wheelhouse, lock_path, scratch = _platform_inputs(
        tmp_path,
        monkeypatch,
        b"wheel",
    )
    facts = SimpleNamespace(
        metadata=b"metadata",
        frontend_manifest=b"frontend",
        migration_inventory=b"migrations",
        encode_reference_source=b"encode",
        bulk_reference_source=b"bulk",
        frontend_identity=f"sha256-{'1' * 64}",
        reference_identity=f"sha256-{'2' * 64}",
    )
    monkeypatch.setattr(builder_module, "verify_platform_wheel", lambda _raw: facts)
    monkeypatch.setattr(
        builder_module.shutil,
        "disk_usage",
        lambda _path: SimpleNamespace(free=0),
    )
    monkeypatch.setattr(
        builder_module,
        "_write_canonical_tar",
        lambda *_args, **_kwargs: pytest.fail("tar writer reached"),
    )
    output = tmp_path / "capacity.tar"

    with pytest.raises(DeploymentError) as caught:
        build_platform_bundle(
            wheel_path,
            wheelhouse,
            lock_path,
            output,
            scratch_root=scratch,
        )

    assert caught.value.issue.code == "DEPLOYMENT_CAPACITY_INSUFFICIENT"
    assert caught.value.issue.recoverable is True
    assert not output.exists()


def test_directory_fsync_failure_does_not_publish_a_bundle(
    tmp_path: Path,
    monkeypatch,
) -> None:
    wheel_path, wheelhouse, lock_path, scratch = _platform_inputs(
        tmp_path,
        monkeypatch,
        b"wheel",
    )
    facts = SimpleNamespace(
        metadata=b"metadata",
        frontend_manifest=b"frontend",
        migration_inventory=b"migrations",
        encode_reference_source=b"encode",
        bulk_reference_source=b"bulk",
        frontend_identity=f"sha256-{'1' * 64}",
        reference_identity=f"sha256-{'2' * 64}",
    )
    monkeypatch.setattr(builder_module, "verify_platform_wheel", lambda _raw: facts)
    original_fsync = builder_module.os.fsync
    calls = 0

    def fail_directory_fsync(descriptor: int) -> None:
        nonlocal calls
        calls += 1
        if calls == 2:
            raise OSError
        original_fsync(descriptor)

    monkeypatch.setattr(builder_module.os, "fsync", fail_directory_fsync)
    output = tmp_path / "platform.tar"

    with pytest.raises(DeploymentError) as caught:
        build_platform_bundle(
            wheel_path,
            wheelhouse,
            lock_path,
            output,
            scratch_root=scratch,
        )

    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_BUILD_FAILED"
    assert not output.exists()
    assert not tuple(tmp_path.glob("*.partial"))
