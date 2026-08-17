from __future__ import annotations

import hashlib
import io
import json
import os
from pathlib import Path
import tarfile

from encode_pipeline.deployment.admission import (
    DatabaseSchemaObservation,
    ResolvedContractFacts,
)
from encode_pipeline.deployment.canonical import canonical_json_bytes
from encode_pipeline.deployment.canonical import canonical_identity
from encode_pipeline.deployment.manager import DeploymentManager, DeploymentOwnership
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    ENCODE_RUNTIME,
    PLATFORM,
    BundleManifest,
    ContractDocument,
    ContractIdentity,
    ContractRequirement,
    FileRecord,
    REQUIRED_PROVIDERS,
)
from encode_pipeline.deployment.state import render_platform_environment


FIXTURE_CONTRACT_SCHEME = "helixweave-test-native-contract-v1"


def sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def component_contract(component: str) -> str:
    return {
        PLATFORM: "helixweave.platform.distribution",
        ENCODE_RUNTIME: "helixweave.runtime.encode",
        BULK_RNASEQ_RUNTIME: "helixweave.runtime.bulk-rnaseq",
    }[component]


def manifest_for(
    component: str,
    *,
    version: str = "1.0.0",
    provider_identity: str | None = None,
    requirements: tuple[ContractRequirement, ...] = (),
    database_heads: tuple[str, ...] | None = None,
    extra_payload: dict[str, bytes] | None = None,
) -> tuple[BundleManifest, dict[str, bytes]]:
    primary = component_contract(component)
    heads = (
        ("schema-v1",)
        if database_heads is None and component == PLATFORM
        else (database_heads or ())
    )
    payload: dict[str, bytes] = {}
    contracts: list[ContractDocument] = []
    for contract in sorted(REQUIRED_PROVIDERS[component]):
        path = f"payload/contracts/{contract.replace('.', '-')}.json"
        document = {
            "schema_version": "helixweave-test-native-contract-v1",
            "component": component,
            "contract": contract,
            "native_revision": (
                provider_identity or version if contract == primary else version
            ),
            "database_heads": list(heads)
            if contract == "helixweave.platform.database-migrations"
            else [],
            "requires": [item.to_dict() for item in sorted(requirements)]
            if contract == primary
            else [],
        }
        content = canonical_json_bytes(document)
        payload[path] = content
        contracts.append(
            ContractDocument(
                contract=contract,
                identity=canonical_identity(
                    document,
                    scheme=FIXTURE_CONTRACT_SCHEME,
                ),
                path=path,
            )
        )
    payload.update(extra_payload or {})
    records = tuple(
        FileRecord(
            path=path,
            size_bytes=len(content),
            sha256=sha256(content),
            mode=0o555 if path.endswith(".sh") else 0o444,
        )
        for path, content in sorted(payload.items())
    )
    manifest = BundleManifest.create(
        component=component,
        contracts=contracts,
        files=records,
    )
    return manifest, payload


class FixtureNativeContractResolver:
    """Finite test resolver proving manifest claims come from indexed bytes."""

    def resolve(
        self,
        root: Path,
        manifest: BundleManifest,
    ) -> ResolvedContractFacts:
        contracts: list[ContractIdentity] = []
        version: str | None = None
        heads: tuple[str, ...] = ()
        requirements: tuple[ContractRequirement, ...] = ()
        primary = component_contract(manifest.component)
        for binding in manifest.contracts:
            raw = (root / binding.path).read_bytes()
            document = json.loads(raw)
            if raw != canonical_json_bytes(document) or set(document) != {
                "schema_version",
                "component",
                "contract",
                "native_revision",
                "database_heads",
                "requires",
            }:
                raise ValueError("invalid native test contract")
            if (
                document["schema_version"] != "helixweave-test-native-contract-v1"
                or document["component"] != manifest.component
                or document["contract"] != binding.contract
                or canonical_identity(
                    document,
                    scheme=FIXTURE_CONTRACT_SCHEME,
                )
                != binding.identity
            ):
                raise ValueError("invalid native test contract")
            contracts.append(ContractIdentity(binding.contract, binding.identity))
            if binding.contract == primary:
                version = document["native_revision"]
                requirements = tuple(
                    ContractRequirement.from_dict(item) for item in document["requires"]
                )
            elif document["requires"]:
                raise ValueError("non-primary test contract has requirements")
            if binding.contract == "helixweave.platform.database-migrations":
                heads = tuple(document["database_heads"])
        if version is None:
            raise ValueError("missing primary native test contract")
        return ResolvedContractFacts(
            component=manifest.component,
            deployment_identity=manifest.identity,
            version=version,
            contracts=tuple(sorted(contracts)),
            requirements=requirements,
            database_heads=heads,
        )


class FixtureDatabaseSchemaObserver:
    def observe(self, state) -> DatabaseSchemaObservation:
        return DatabaseSchemaObservation.create(
            provider_identity="sha256-" + "d" * 64,
            state_identity=state.identity,
            database_identity="sha256-" + "e" * 64,
            heads=("schema-v1",),
        )


def manager_for(layout) -> DeploymentManager:
    return DeploymentManager(
        layout,
        ownership=DeploymentOwnership(os.getuid(), os.getgid()),
        contract_resolver=FixtureNativeContractResolver(),
        schema_observer=FixtureDatabaseSchemaObserver(),
        platform_environment_renderer=lambda state: render_platform_environment(
            layout,
            state,
            api_contract_sha256="a" * 64,
        ),
    )


def write_bundle(
    path: Path,
    manifest: BundleManifest,
    payload: dict[str, bytes],
    *,
    mutate_member=None,
    extra_member: tuple[str, bytes] | None = None,
    archive_format: int = tarfile.USTAR_FORMAT,
) -> Path:
    with tarfile.open(path, mode="w", format=archive_format) as archive:
        manifest_bytes = canonical_json_bytes(manifest.to_dict())
        header = tarfile.TarInfo("manifest.json")
        header.size = len(manifest_bytes)
        header.mode = 0o444
        header.uid = 0
        header.gid = 0
        header.mtime = 0
        archive.addfile(header, io.BytesIO(manifest_bytes))
        records = {record.path: record for record in manifest.files}
        for index, name in enumerate(sorted(payload)):
            content = payload[name]
            member = tarfile.TarInfo(name)
            member.size = len(content)
            member.mode = records[name].mode
            member.uid = 0
            member.gid = 0
            member.mtime = 0
            if mutate_member is not None:
                mutate_member(index, member)
            archive.addfile(member, io.BytesIO(content) if member.isreg() else None)
        if extra_member is not None:
            name, content = extra_member
            member = tarfile.TarInfo(name)
            member.size = len(content)
            member.mode = 0o444
            member.uid = 0
            member.gid = 0
            member.mtime = 0
            archive.addfile(member, io.BytesIO(content))
    path.chmod(0o644)
    return path
