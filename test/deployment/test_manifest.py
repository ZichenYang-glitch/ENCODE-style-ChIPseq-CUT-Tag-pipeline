from __future__ import annotations

import copy

import pytest

from encode_pipeline.deployment.admission import (
    ResolvedContractFacts,
    compatible_resolved_facts,
)
from encode_pipeline.deployment.canonical import canonical_identity
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.models import (
    BUNDLE_IDENTITY_SCHEME,
    ENCODE_RUNTIME,
    PLATFORM,
    BundleManifest,
    ContractDocument,
    ContractRequirement,
    FileRecord,
)
from .support import manifest_for


def test_manifest_identity_is_canonical_and_round_trips() -> None:
    manifest, _payload = manifest_for(PLATFORM)

    assert BundleManifest.from_dict(manifest.to_dict()) == manifest
    assert manifest.identity.startswith("sha256-")
    assert manifest.to_dict()["files"] == sorted(
        manifest.to_dict()["files"], key=lambda item: item["path"]
    )


def test_manifest_rejects_identity_tampering_and_unknown_fields() -> None:
    manifest, _payload = manifest_for(PLATFORM)
    document = manifest.to_dict()
    document["version"] = "1.0.1"

    with pytest.raises(DeploymentError) as captured:
        BundleManifest.from_dict(document)
    assert captured.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"

    document = manifest.to_dict()
    document["private_path"] = "/secret/reference"
    document["identity"] = canonical_identity(
        {key: value for key, value in document.items() if key != "identity"},
        scheme=BUNDLE_IDENTITY_SCHEME,
    )
    with pytest.raises(DeploymentError):
        BundleManifest.from_dict(document)

    document = manifest.to_dict()
    document["requires"] = []
    document["identity"] = canonical_identity(
        {key: value for key, value in document.items() if key != "identity"},
        scheme=BUNDLE_IDENTITY_SCHEME,
    )
    with pytest.raises(DeploymentError):
        BundleManifest.from_dict(document)


def test_manifest_rejects_traversal_and_case_colliding_payloads() -> None:
    with pytest.raises(DeploymentError):
        BundleManifest.create(
            component=PLATFORM,
            contracts=manifest_for(PLATFORM)[0].contracts,
            files=(FileRecord("payload/../escape", 0, "0" * 64, 0o444),),
        )

    original, _payload = manifest_for(PLATFORM)
    record = original.files[0]
    with pytest.raises(DeploymentError):
        BundleManifest.create(
            component=PLATFORM,
            contracts=original.contracts,
            files=(
                *original.files,
                FileRecord(
                    record.path.upper(),
                    record.size_bytes,
                    record.sha256,
                    record.mode,
                ),
            ),
        )


def test_manifest_requires_component_owned_fact_contracts() -> None:
    manifest, _payload = manifest_for(PLATFORM)
    document = copy.deepcopy(manifest.to_dict())
    document["contracts"] = document["contracts"][:-1]
    document["identity"] = canonical_identity(
        {key: value for key, value in document.items() if key != "identity"},
        scheme=BUNDLE_IDENTITY_SCHEME,
    )

    with pytest.raises(DeploymentError):
        BundleManifest.from_dict(document)


def test_runtime_manifest_cannot_compete_with_platform_fact_sources() -> None:
    runtime, _payload = manifest_for(ENCODE_RUNTIME)

    with pytest.raises(DeploymentError) as caught:
        BundleManifest.create(
            component=ENCODE_RUNTIME,
            contracts=(
                *runtime.contracts,
                ContractDocument(
                    contract="helixweave.platform.distribution",
                    identity="sha256-" + "a" * 64,
                    path=runtime.contracts[0].path,
                ),
            ),
            files=runtime.files,
        )

    assert caught.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"


def test_active_contract_compatibility_is_fail_closed() -> None:
    runtime, _payload = manifest_for(
        ENCODE_RUNTIME,
        provider_identity="encode-runtime-v1",
    )
    runtime_identity = next(
        item.identity
        for item in runtime.provides
        if item.contract == "helixweave.runtime.encode"
    )
    platform, _payload = manifest_for(
        PLATFORM,
        requirements=(
            ContractRequirement("helixweave.runtime.encode", (runtime_identity,)),
        ),
    )
    incompatible, _payload = manifest_for(
        ENCODE_RUNTIME,
        version="2.0.0",
        provider_identity="encode-runtime-v2",
    )

    requirement = ContractRequirement("helixweave.runtime.encode", (runtime_identity,))
    platform_facts = ResolvedContractFacts(
        component=PLATFORM,
        deployment_identity=platform.identity,
        version="1.0.0",
        contracts=platform.provides,
        requirements=(requirement,),
        database_heads=("schema-v1",),
    )
    runtime_facts = ResolvedContractFacts(
        component=ENCODE_RUNTIME,
        deployment_identity=runtime.identity,
        version="1.0.0",
        contracts=runtime.provides,
    )
    incompatible_facts = ResolvedContractFacts(
        component=ENCODE_RUNTIME,
        deployment_identity=incompatible.identity,
        version="2.0.0",
        contracts=incompatible.provides,
    )

    assert compatible_resolved_facts((platform_facts, runtime_facts)) is True
    assert compatible_resolved_facts((platform_facts, incompatible_facts)) is False


def test_manifest_rejects_duplicate_file_paths_with_different_records() -> None:
    original, _payload = manifest_for(ENCODE_RUNTIME)
    record = original.files[0]

    with pytest.raises(DeploymentError) as caught:
        BundleManifest.create(
            component=ENCODE_RUNTIME,
            contracts=original.contracts,
            files=(
                *original.files,
                FileRecord(
                    record.path,
                    record.size_bytes + 1,
                    "f" * 64,
                    record.mode,
                ),
            ),
        )

    assert caught.value.issue.code == "DEPLOYMENT_MANIFEST_INVALID"
