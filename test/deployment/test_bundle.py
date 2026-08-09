from __future__ import annotations

from dataclasses import replace
import hashlib
import inspect
import os
import json
from pathlib import Path
from types import SimpleNamespace
import tarfile

import pytest

from encode_pipeline.deployment.admission import DatabaseSchemaObservation
from encode_pipeline.deployment.bundle import BundleStore
import encode_pipeline.deployment.bundle as bundle_module
from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.manager import DeploymentManager, DeploymentOwnership
from encode_pipeline.deployment.models import (
    ENCODE_RUNTIME,
    PLATFORM,
    BundleManifest,
    ContractRequirement,
    FileRecord,
)
from .support import (
    FixtureNativeContractResolver,
    manager_for,
    manifest_for,
    write_bundle,
)


def test_bundle_stages_and_verifies_immutable_content(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manifest, payload = manifest_for(ENCODE_RUNTIME)
    bundle = write_bundle(tmp_path / "encode.tar", manifest, payload)
    store = BundleStore(layout)

    staged = store.stage(bundle)

    assert staged == manifest
    assert store.verify_installed(ENCODE_RUNTIME, manifest.identity) == manifest
    release = layout.encode_runtimes / manifest.identity
    assert release.stat().st_mode & 0o777 == 0o555
    assert (release / manifest.files[0].path).stat().st_mode & 0o777 == 0o444


def test_bundle_rejects_symlink_member_and_unexpected_member(tmp_path: Path) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manifest, payload = manifest_for(ENCODE_RUNTIME)

    def make_symlink(_index, member: tarfile.TarInfo) -> None:
        member.type = tarfile.SYMTYPE
        member.linkname = "/etc/passwd"
        member.size = 0

    symlink_bundle = write_bundle(
        tmp_path / "symlink.tar",
        manifest,
        payload,
        mutate_member=make_symlink,
    )
    with pytest.raises(DeploymentError) as captured:
        BundleStore(layout).stage(symlink_bundle)
    assert captured.value.issue.code == "DEPLOYMENT_BUNDLE_INVALID"

    extra_bundle = write_bundle(
        tmp_path / "extra.tar",
        manifest,
        payload,
        extra_member=("payload/extra", b"extra"),
    )
    with pytest.raises(DeploymentError):
        BundleStore(layout).stage(extra_bundle)


@pytest.mark.parametrize("field", ("mtime", "owner-name", "trailing"))
def test_bundle_rejects_noncanonical_header_or_padding(
    tmp_path: Path,
    field: str,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manifest, payload = manifest_for(ENCODE_RUNTIME)

    def mutate(_index, member: tarfile.TarInfo) -> None:
        if field == "mtime":
            member.mtime = 1
        elif field == "owner-name":
            member.uname = "root"

    bundle = write_bundle(
        tmp_path / f"{field}.tar",
        manifest,
        payload,
        mutate_member=mutate,
    )
    if field == "trailing":
        with bundle.open("r+b") as stream:
            stream.seek(-1, os.SEEK_END)
            stream.write(b"X")

    with pytest.raises(DeploymentError) as caught:
        BundleStore(layout).stage(bundle)
    assert caught.value.issue.code == "DEPLOYMENT_BUNDLE_INVALID"


def test_bundle_rejects_group_writable_ingress_and_tampered_install(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manifest, payload = manifest_for(ENCODE_RUNTIME)
    bundle = write_bundle(tmp_path / "encode.tar", manifest, payload)
    bundle.chmod(0o664)

    with pytest.raises(DeploymentError) as captured:
        BundleStore(layout).stage(bundle)
    assert captured.value.issue.code == "DEPLOYMENT_BUNDLE_INVALID"

    bundle.chmod(0o644)
    hardlink = tmp_path / "hardlinked-encode.tar"
    os.link(bundle, hardlink)
    with pytest.raises(DeploymentError) as hardlink_error:
        BundleStore(layout).stage(bundle)
    assert hardlink_error.value.issue.code == "DEPLOYMENT_BUNDLE_INVALID"
    hardlink.unlink()

    store = BundleStore(layout)
    store.stage(bundle)
    with pytest.raises(DeploymentError) as owner_error:
        store.verify_installed(
            ENCODE_RUNTIME,
            manifest.identity,
            expected_owner_uid=os.getuid(),
            expected_owner_gid=os.getgid() + 1,
        )
    assert owner_error.value.issue.code == "DEPLOYMENT_RELEASE_INVALID"

    installed = layout.encode_runtimes / manifest.identity / manifest.files[0].path
    installed.chmod(0o644)
    installed.write_bytes(b"tampered")
    with pytest.raises(DeploymentError) as captured:
        store.verify_installed(ENCODE_RUNTIME, manifest.identity)
    assert captured.value.issue.code == "DEPLOYMENT_RELEASE_INVALID"


def test_committed_release_without_state_update_is_reported_as_interrupted(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manifest, payload = manifest_for(ENCODE_RUNTIME)
    bundle = write_bundle(tmp_path / "encode.tar", manifest, payload)
    manager = manager_for(layout)

    def crash(point: str) -> None:
        if point == "bundle:release-committed":
            raise RuntimeError("injected")

    with pytest.raises(RuntimeError, match="injected"):
        manager.stage(bundle, fault=crash)

    status = manager.verify()
    assert status.interrupted is True
    assert status.orphaned_deployments == ((ENCODE_RUNTIME, manifest.identity),)


def test_platform_and_runtime_activate_and_rollback_independently(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    runtime, runtime_payload = manifest_for(
        ENCODE_RUNTIME,
        provider_identity="encode-v1",
    )
    runtime_bundle = write_bundle(tmp_path / "encode-v1.tar", runtime, runtime_payload)
    manager.stage(runtime_bundle)
    manager.activate(
        ENCODE_RUNTIME,
        expected_staged_identity=runtime.identity,
    )

    runtime_contract_identity = next(
        item.identity
        for item in runtime.provides
        if item.contract == "helixweave.runtime.encode"
    )
    requirement = ContractRequirement(
        "helixweave.runtime.encode", (runtime_contract_identity,)
    )
    platform_v1, platform_v1_payload = manifest_for(
        PLATFORM,
        requirements=(requirement,),
    )
    manager.stage(
        write_bundle(tmp_path / "platform-v1.tar", platform_v1, platform_v1_payload)
    )
    manager.activate(
        PLATFORM,
        expected_staged_identity=platform_v1.identity,
    )

    platform_v2, platform_v2_payload = manifest_for(
        PLATFORM,
        version="1.1.0",
        requirements=(requirement,),
    )
    manager.stage(
        write_bundle(tmp_path / "platform-v2.tar", platform_v2, platform_v2_payload)
    )
    upgraded = manager.activate(
        PLATFORM,
        expected_staged_identity=platform_v2.identity,
    )
    assert upgraded.components[PLATFORM].active == platform_v2.identity
    assert upgraded.components[PLATFORM].previous == platform_v1.identity

    rolled_back = manager.rollback(
        PLATFORM,
        expected_previous_identity=platform_v1.identity,
    )
    assert rolled_back.components[PLATFORM].active == platform_v1.identity
    assert rolled_back.components[PLATFORM].previous == platform_v2.identity


def test_incompatible_runtime_upgrade_leaves_active_state_unchanged(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    runtime_v1, payload_v1 = manifest_for(
        ENCODE_RUNTIME,
        provider_identity="encode-v1",
    )
    manager.stage(write_bundle(tmp_path / "runtime-v1.tar", runtime_v1, payload_v1))
    manager.activate(
        ENCODE_RUNTIME,
        expected_staged_identity=runtime_v1.identity,
    )
    runtime_v1_identity = next(
        item.identity
        for item in runtime_v1.provides
        if item.contract == "helixweave.runtime.encode"
    )
    platform, platform_payload = manifest_for(
        PLATFORM,
        requirements=(
            ContractRequirement("helixweave.runtime.encode", (runtime_v1_identity,)),
        ),
    )
    manager.stage(write_bundle(tmp_path / "platform.tar", platform, platform_payload))
    manager.activate(
        PLATFORM,
        expected_staged_identity=platform.identity,
    )
    before = manager.states.read()

    runtime_v2, payload_v2 = manifest_for(
        ENCODE_RUNTIME,
        version="2.0.0",
        provider_identity="encode-v2",
    )
    manager.stage(write_bundle(tmp_path / "runtime-v2.tar", runtime_v2, payload_v2))
    with pytest.raises(DeploymentError) as captured:
        manager.activate(
            ENCODE_RUNTIME,
            expected_staged_identity=runtime_v2.identity,
        )
    assert captured.value.issue.code == "DEPLOYMENT_COMPATIBILITY_FAILED"
    assert manager.states.read().components[ENCODE_RUNTIME].active == (
        before.components[ENCODE_RUNTIME].active
    )


def test_verify_rechecks_active_contract_and_schema_compatibility(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    runtime, runtime_payload = manifest_for(ENCODE_RUNTIME)
    manager.stage(write_bundle(tmp_path / "runtime.tar", runtime, runtime_payload))
    manager.activate(ENCODE_RUNTIME, expected_staged_identity=runtime.identity)
    platform, platform_payload = manifest_for(PLATFORM)
    manager.stage(write_bundle(tmp_path / "platform.tar", platform, platform_payload))
    manager.activate(PLATFORM, expected_staged_identity=platform.identity)

    native = FixtureNativeContractResolver()

    class IncompatibleResolver:
        def resolve(self, root, manifest):
            facts = native.resolve(root, manifest)
            if manifest.component != PLATFORM:
                return facts
            return replace(
                facts,
                requirements=(
                    ContractRequirement(
                        "helixweave.runtime.encode",
                        ("sha256-" + "f" * 64,),
                    ),
                ),
            )

    manager.contract_resolver = IncompatibleResolver()
    with pytest.raises(DeploymentError) as contract_error:
        manager.verify()
    assert contract_error.value.issue.code == "DEPLOYMENT_COMPATIBILITY_FAILED"

    manager.contract_resolver = native

    class IncompatibleSchemaObserver:
        def observe(self, state):
            return DatabaseSchemaObservation.create(
                provider_identity="sha256-" + "d" * 64,
                state_identity=state.identity,
                database_identity="sha256-" + "e" * 64,
                heads=("schema-v2",),
            )

    manager.schema_observer = IncompatibleSchemaObserver()
    with pytest.raises(DeploymentError) as schema_error:
        manager.verify()
    assert schema_error.value.issue.code == "DEPLOYMENT_SCHEMA_INCOMPATIBLE"


def test_activation_and_rollback_require_the_explicit_slot_identity(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    runtime, payload = manifest_for(ENCODE_RUNTIME)
    manager.stage(write_bundle(tmp_path / "runtime.tar", runtime, payload))

    with pytest.raises(DeploymentError) as activation_error:
        manager.activate(
            ENCODE_RUNTIME,
            expected_staged_identity="sha256-" + "f" * 64,
        )
    assert activation_error.value.issue.code == ("DEPLOYMENT_STAGED_IDENTITY_MISMATCH")
    assert manager.states.read().components[ENCODE_RUNTIME].active is None

    manager.activate(
        ENCODE_RUNTIME,
        expected_staged_identity=runtime.identity,
    )
    with pytest.raises(DeploymentError) as rollback_error:
        manager.rollback(
            ENCODE_RUNTIME,
            expected_previous_identity="sha256-" + "e" * 64,
        )
    assert rollback_error.value.issue.code == ("DEPLOYMENT_PREVIOUS_IDENTITY_MISMATCH")


def test_native_contract_resolver_rejects_a_reindexed_false_provider_claim(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifest, payload = manifest_for(ENCODE_RUNTIME)
    binding = manifest.contracts[0]
    document = json.loads(payload[binding.path])
    document["native_revision"] = "tampered-v2"
    changed = canonical_json_bytes(document)
    payload[binding.path] = changed
    files = tuple(
        FileRecord(
            path=item.path,
            size_bytes=len(changed),
            sha256=hashlib.sha256(changed).hexdigest(),
            mode=item.mode,
        )
        if item.path == binding.path
        else item
        for item in manifest.files
    )
    tampered = BundleManifest.create(
        component=manifest.component,
        contracts=manifest.contracts,
        files=files,
    )
    manager.stage(write_bundle(tmp_path / "tampered.tar", tampered, payload))

    with pytest.raises(DeploymentError) as caught:
        manager.activate(
            ENCODE_RUNTIME,
            expected_staged_identity=tampered.identity,
        )

    assert caught.value.issue.code == "DEPLOYMENT_CONTRACT_ADMISSION_FAILED"


def test_phase_a_manager_cannot_activate_without_native_contract_admission(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = DeploymentManager(
        layout,
        ownership=DeploymentOwnership(os.getuid(), os.getgid()),
    )
    manifest, payload = manifest_for(ENCODE_RUNTIME)
    manager.stage(write_bundle(tmp_path / "runtime.tar", manifest, payload))

    with pytest.raises(DeploymentError) as caught:
        manager.activate(
            ENCODE_RUNTIME,
            expected_staged_identity=manifest.identity,
        )

    assert caught.value.issue.code == "DEPLOYMENT_CONTRACT_ADMISSION_DEFERRED"


def test_activation_uses_a_state_bound_schema_observer_not_caller_heads(
    tmp_path: Path,
) -> None:
    assert (
        "observed_database_heads"
        not in inspect.signature(DeploymentManager.activate).parameters
    )
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    platform, payload = manifest_for(PLATFORM)
    manager.stage(write_bundle(tmp_path / "platform.tar", platform, payload))

    class WrongStateObserver:
        def observe(self, _state):
            return DatabaseSchemaObservation.create(
                provider_identity="sha256-" + "d" * 64,
                state_identity="sha256-" + "0" * 64,
                database_identity="sha256-" + "e" * 64,
                heads=("schema-v1",),
            )

    manager.schema_observer = WrongStateObserver()
    with pytest.raises(DeploymentError) as caught:
        manager.activate(
            PLATFORM,
            expected_staged_identity=platform.identity,
        )

    assert caught.value.issue.code == "DEPLOYMENT_SCHEMA_OBSERVATION_FAILED"
    assert manager.states.read().components[PLATFORM].active is None


def test_completed_historical_releases_are_not_reported_as_interrupted(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manager = manager_for(layout)
    manifests: list[BundleManifest] = []

    for index in range(3):
        manifest, payload = manifest_for(
            ENCODE_RUNTIME,
            version=f"1.{index}.0",
            provider_identity=f"encode-v{index}",
        )
        manifests.append(manifest)
        manager.stage(
            write_bundle(tmp_path / f"runtime-{index}.tar", manifest, payload)
        )
        manager.activate(
            ENCODE_RUNTIME,
            expected_staged_identity=manifest.identity,
        )

    status = manager.verify()
    assert status.state.components[ENCODE_RUNTIME].active == manifests[-1].identity
    assert status.state.components[ENCODE_RUNTIME].previous == manifests[-2].identity
    assert status.orphaned_deployments == ()
    assert status.interrupted is False


def test_stage_fails_before_extraction_when_capacity_is_insufficient(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    manifest, payload = manifest_for(ENCODE_RUNTIME)
    bundle = write_bundle(tmp_path / "encode.tar", manifest, payload)
    monkeypatch.setattr(
        bundle_module.shutil,
        "disk_usage",
        lambda _path: SimpleNamespace(free=0),
    )

    with pytest.raises(DeploymentError) as caught:
        BundleStore(layout).stage(bundle)

    assert caught.value.issue.code == "DEPLOYMENT_CAPACITY_INSUFFICIENT"
    assert not (layout.encode_runtimes / manifest.identity).exists()
    assert BundleStore(layout).partial_staging() == ()
