from __future__ import annotations

import base64
from copy import deepcopy
import csv
from hashlib import sha256
from importlib import resources
import io
import json
import sqlite3
from pathlib import Path
from types import SimpleNamespace
import zipfile

import pytest

import encode_pipeline.deployment.native_contracts as native_module
from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
    RUNTIME_IDENTITY_FILE,
    RUNTIME_IDENTITY_SHA256,
    VerifiedRuntimeAssets,
)
from encode_pipeline.deployment.canonical import (
    canonical_identity,
    canonical_json_bytes,
    without_key,
)
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import (
    BULK_RNASEQ_RUNTIME,
    PLATFORM,
    BundleManifest,
    ContractDocument,
    DeploymentState,
    FileRecord,
)
from encode_pipeline.deployment.native_contracts import (
    InventoryDatabaseSchemaObserver,
    ProductionNativeContractResolver,
    encode_environment_paths,
    encode_runtime_index_bytes,
    parse_encode_explicit_lock,
    parse_encode_runtime_index,
    verify_static_micromamba,
)
from encode_pipeline.persistence.migration_admission import (
    MigrationInventoryTrustAnchor,
    verify_migration_execution_inventory,
)
from encode_pipeline.platform.results import Result


def _write_database(path: Path, head: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    connection = sqlite3.connect(path)
    try:
        connection.execute(
            "CREATE TABLE alembic_version (version_num VARCHAR(128) NOT NULL)"
        )
        connection.execute(
            "INSERT INTO alembic_version (version_num) VALUES (?)",
            (head,),
        )
        connection.commit()
    finally:
        connection.close()
    path.chmod(0o440)


def _explicit_lock(*urls: str) -> bytes:
    return (
        "# platform: linux-64\n@EXPLICIT\n" + "".join(f"{url}\n" for url in urls)
    ).encode()


def _package(*, suffix: str = "tool") -> dict[str, object]:
    digest = ("1" if suffix == "tool" else "2") * 64
    filename = f"{suffix}-1.0-0.conda"
    url = f"https://conda.anaconda.org/conda-forge/linux-64/{filename}#{'a' * 32}"
    return {
        "url": url,
        "platform": "linux-64",
        "filename": filename,
        "md5": "a" * 32,
        "archive_path": (
            f"{native_module.ENCODE_PACKAGE_ARCHIVE_ROOT}/{digest}/{filename}"
        ),
        "size_bytes": 7,
        "sha256": digest,
    }


def _runtime_index_document() -> dict[str, object]:
    package = _package()
    content = encode_runtime_index_bytes(
        workflow_build_identity=f"sha256-{'3' * 64}",
        micromamba={
            "path": native_module.ENCODE_MICROMAMBA_PATH,
            "size_bytes": 120,
            "sha256": "4" * 64,
        },
        packages=(package,),
        environments=(
            {
                "environment_path": "workflow/envs/runner.yml",
                "environment_sha256": "5" * 64,
                "lock_path": "workflow/envs/runner.lock",
                "lock_sha256": "6" * 64,
                "packages": [package["url"]],
            },
        ),
    )
    return json.loads(content)


def _canonical_index(document: dict[str, object]) -> bytes:
    document["identity"] = canonical_identity(
        without_key(document, "identity"),
        scheme=native_module.ENCODE_RUNTIME_INDEX_IDENTITY_SCHEME,
    )
    return canonical_json_bytes(document)


def _verified_bulk_assets(root: Path) -> VerifiedRuntimeAssets:
    return VerifiedRuntimeAssets(
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


def _bulk_manifest(root: Path, content: bytes) -> BundleManifest:
    path = native_module.BULK_RUNTIME_IDENTITY_PATH
    target = root.joinpath(*Path(path).parts)
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_bytes(content)
    target.chmod(0o444)
    digest = sha256(content).hexdigest()
    return BundleManifest.create(
        component=BULK_RNASEQ_RUNTIME,
        contracts=(
            ContractDocument(
                native_module.BULK_RUNTIME_CONTRACT,
                f"sha256-{digest}",
                path,
            ),
        ),
        files=(FileRecord(path, len(content), digest, 0o444),),
    )


def test_checked_in_encode_locks_use_the_supported_package_filename_grammar() -> None:
    repository = Path(__file__).resolve().parents[2]
    locks = tuple(sorted((repository / "workflow" / "envs").glob("*.lock")))

    coordinates = tuple(
        coordinate
        for lock in locks
        for coordinate in parse_encode_explicit_lock(lock.read_bytes())
    )

    assert locks
    assert coordinates
    assert any(item.filename.startswith("_") for item in coordinates)
    assert all("/" not in item.filename for item in coordinates)


@pytest.mark.parametrize(
    "content",
    (
        b"\xff",
        b"# platform: linux-64\n@EXPLICIT",
        b"unexpected\n# platform: linux-64\n@EXPLICIT\n",
        b"# platform: linux-64\n@EXPLICIT\n\n",
        _explicit_lock("http://attacker.invalid/tool.conda#" + "a" * 32),
        b"# platform: linux-64\n@EXPLICIT\n",
        _explicit_lock(
            "https://conda.anaconda.org/conda-forge/linux-64/tool.conda#" + "a" * 32,
            "https://conda.anaconda.org/conda-forge/linux-64/tool.conda#" + "a" * 32,
        ),
    ),
)
def test_explicit_lock_rejects_noncanonical_or_ambiguous_coordinates(
    content: bytes,
) -> None:
    with pytest.raises(native_module._NativeContractFault):
        parse_encode_explicit_lock(content)


@pytest.mark.parametrize(
    "manifest",
    (
        (("workflow/envs/runner.yml", b"one"), ("workflow/envs/runner.yml", b"two")),
        (
            ("workflow/envs/runner.yml", b"name: runner\n"),
            ("workflow/envs/runner.lock", b"lock\n"),
            ("workflow/rules/bad.smk", b"\xff"),
        ),
        (
            ("workflow/envs/runner.yml", b"name: runner\n"),
            ("workflow/envs/runner.lock", b"lock\n"),
            ("workflow/rules/bad.smk", b"rule bad:\n    conda: dynamic_name\n"),
        ),
        (
            ("workflow/envs/runner.yml", b"name: runner\n"),
            ("workflow/envs/runner.lock", b"lock\n"),
            (
                "workflow/rules/bad.smk",
                b'rule bad:\n    conda: "../../../outside.yml"\n',
            ),
        ),
        (("workflow/envs/runner.yml", b"name: runner\n"),),
    ),
)
def test_environment_discovery_fails_closed_on_untrusted_source_manifest(
    manifest: tuple[tuple[str, bytes], ...],
) -> None:
    with pytest.raises(native_module._NativeContractFault):
        encode_environment_paths(manifest)


@pytest.mark.parametrize("fault", ("header", "table", "dynamic"))
def test_micromamba_requires_a_bounded_supported_linux_elf(fault: str) -> None:
    content = bytearray(120)
    content[:7] = b"\x7fELF\x02\x01\x01"
    content[18:20] = (62).to_bytes(2, "little")
    content[32:40] = (64).to_bytes(8, "little")
    content[54:56] = (56).to_bytes(2, "little")
    content[56:58] = (1).to_bytes(2, "little")
    content[64:68] = (1).to_bytes(4, "little")
    if fault == "header":
        content[:7] = b"invalid"
    elif fault == "table":
        content[32:40] = (119).to_bytes(8, "little")
    else:
        content[64:68] = (3).to_bytes(4, "little")

    with pytest.raises(native_module._NativeContractFault):
        verify_static_micromamba(bytes(content))


def test_micromamba_accepts_only_the_supported_glibc_loader() -> None:
    interpreter = b"/lib64/ld-linux-x86-64.so.2\x00"
    content = bytearray(192)
    content[:7] = b"\x7fELF\x02\x01\x01"
    content[18:20] = (62).to_bytes(2, "little")
    content[32:40] = (64).to_bytes(8, "little")
    content[54:56] = (56).to_bytes(2, "little")
    content[56:58] = (2).to_bytes(2, "little")
    content[64:68] = (1).to_bytes(4, "little")
    content[120:124] = (3).to_bytes(4, "little")
    content[128:136] = (160).to_bytes(8, "little")
    content[152:160] = len(interpreter).to_bytes(8, "little")
    content[160 : 160 + len(interpreter)] = interpreter

    verify_static_micromamba(bytes(content))

    content[160 : 160 + len(interpreter)] = b"/tmp/ld-linux-x86-64.so.2".ljust(
        len(interpreter), b"\x00"
    )
    with pytest.raises(native_module._NativeContractFault):
        verify_static_micromamba(bytes(content))


@pytest.mark.parametrize(
    "fault",
    (
        "schema",
        "micromamba",
        "no-packages",
        "package-shape",
        "package-binding",
        "package-order",
        "no-environments",
        "environment-shape",
        "environment-binding",
        "runner-missing",
    ),
)
def test_encode_runtime_index_rejects_each_ambiguous_native_binding(
    fault: str,
) -> None:
    document = deepcopy(_runtime_index_document())
    packages = document["packages"]
    environments = document["environments"]
    assert isinstance(packages, list)
    assert isinstance(environments, list)
    if fault == "schema":
        document["schema_version"] = "future"
    elif fault == "micromamba":
        document["micromamba"]["size_bytes"] = 0
    elif fault == "no-packages":
        packages.clear()
    elif fault == "package-shape":
        packages[0].pop("archive_path")
    elif fault == "package-binding":
        packages[0]["archive_path"] = "payload/native/encode/packages/wrong"
    elif fault == "package-order":
        second = _package(suffix="ztool")
        packages.append(second)
        packages.reverse()
        environments[0]["packages"] = sorted((packages[0]["url"], packages[1]["url"]))
    elif fault == "no-environments":
        environments.clear()
    elif fault == "environment-shape":
        environments[0].pop("lock_path")
    elif fault == "environment-binding":
        environments[0]["lock_path"] = "workflow/envs/other.lock"
    else:
        environments[0]["environment_path"] = "workflow/envs/tool.yml"
        environments[0]["lock_path"] = "workflow/envs/tool.lock"

    with pytest.raises(native_module._NativeContractFault):
        parse_encode_runtime_index(_canonical_index(document))


def test_bulk_native_admission_defaults_to_the_static_runtime_closure() -> None:
    from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
        verify_runtime_asset_closure,
    )

    resolver = ProductionNativeContractResolver()

    assert resolver._bulk_runtime_verifier is verify_runtime_asset_closure


def test_bulk_native_resolver_admits_only_the_committed_runtime_identity(
    tmp_path: Path,
) -> None:
    content = (
        resources.files("encode_pipeline.contracts.nfcore_rnaseq")
        .joinpath(RUNTIME_IDENTITY_FILE)
        .read_bytes()
    )
    root = tmp_path / "bulk-release"
    manifest = _bulk_manifest(root, content)
    runtime_root = root.joinpath(*Path(native_module.BULK_RUNTIME_ROOT_PATH).parts)
    observed_bindings = []

    def verify(binding):
        observed_bindings.append(binding)
        return Result.success(_verified_bulk_assets(binding.root))

    resolver = ProductionNativeContractResolver(bulk_runtime_verifier=verify)
    facts = resolver.resolve(root.resolve(), manifest)

    assert facts.component == BULK_RNASEQ_RUNTIME
    assert facts.deployment_identity == manifest.identity
    assert facts.contracts == manifest.provides
    assert facts.version
    assert [binding.root for binding in observed_bindings] == [runtime_root]

    tampered_root = tmp_path / "tampered-bulk-release"
    tampered = _bulk_manifest(tampered_root, b'{"source":{"release":"3.0.0"}}\n')
    with pytest.raises(DeploymentError) as caught:
        resolver.resolve(tampered_root.resolve(), tampered)
    assert caught.value.issue.code == "DEPLOYMENT_CONTRACT_ADMISSION_FAILED"
    assert len(observed_bindings) == 1


def test_platform_native_resolver_composes_verified_candidate_facts(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    identity = f"sha256-{'a' * 64}"
    wheel = b"candidate wheel"
    payload = {
        native_module.PLATFORM_METADATA_PATH: b"metadata",
        native_module.PLATFORM_WHEEL_PATH: wheel,
        native_module.PLATFORM_RUNTIME_LOCK_PATH: b"runtime lock",
        native_module.PLATFORM_FRONTEND_PATH: b"frontend",
        native_module.PLATFORM_MIGRATIONS_PATH: b"migrations",
        native_module.PLATFORM_ENCODE_REFERENCES_PATH: b"encode refs",
        native_module.PLATFORM_BULK_REFERENCES_PATH: b"bulk refs",
    }
    facts = SimpleNamespace(
        version="2.0.0",
        metadata=payload[native_module.PLATFORM_METADATA_PATH],
        frontend_manifest=payload[native_module.PLATFORM_FRONTEND_PATH],
        migration_inventory=payload[native_module.PLATFORM_MIGRATIONS_PATH],
        encode_reference_source=payload[native_module.PLATFORM_ENCODE_REFERENCES_PATH],
        bulk_reference_source=payload[native_module.PLATFORM_BULK_REFERENCES_PATH],
        frontend_identity=f"sha256-{'b' * 64}",
        migration=SimpleNamespace(heads=("head",)),
        reference_identity=f"sha256-{'c' * 64}",
        bulk_runtime_identity=f"sha256-{'d' * 64}",
    )
    locked = SimpleNamespace(
        identity=identity,
        wheels=(SimpleNamespace(sha256=sha256(wheel).hexdigest()),),
    )
    closure = SimpleNamespace(
        identity=identity, lock_identity=identity, files=(object(),)
    )
    requested: list[tuple[str, str]] = []

    def contract_bytes(_root, _manifest, contract, path, *, maximum_bytes):
        assert maximum_bytes > 0
        requested.append((contract, path))
        return payload[path]

    monkeypatch.setattr(native_module, "_contract_bytes", contract_bytes)
    monkeypatch.setattr(
        native_module,
        "_indexed_bytes",
        lambda _root, _manifest, path, **_kwargs: payload[path],
    )
    monkeypatch.setattr(native_module, "verify_platform_wheel", lambda value: facts)
    monkeypatch.setattr(
        native_module.PlatformWheelLock,
        "from_bytes",
        lambda value: locked,
    )
    monkeypatch.setattr(
        native_module,
        "inspect_platform_runtime_closure",
        lambda *_args, **_kwargs: closure,
    )
    monkeypatch.setattr(
        native_module, "_verify_platform_runtime_records", lambda *_: None
    )
    manifest = SimpleNamespace(component=PLATFORM, identity=identity, files=())
    resolver = ProductionNativeContractResolver()

    resolved = resolver.resolve(tmp_path.resolve(), manifest)

    assert resolved.component == PLATFORM
    assert resolved.version == "2.0.0"
    assert resolved.database_heads == ("head",)
    assert resolved.requirements[0].accepted_identities == (
        facts.bulk_runtime_identity,
    )
    assert len(requested) == 6

    locked.wheels = ()
    with pytest.raises(DeploymentError) as caught:
        resolver.resolve(tmp_path.resolve(), manifest)
    assert caught.value.issue.code == "DEPLOYMENT_CONTRACT_ADMISSION_FAILED"


def test_native_resolver_rejects_invalid_coordinates_and_sanitizes_faults(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    resolver = ProductionNativeContractResolver()
    manifest = SimpleNamespace(component=PLATFORM, identity=f"sha256-{'a' * 64}")

    with pytest.raises(DeploymentError) as relative:
        resolver.resolve(Path("relative"), manifest)
    assert relative.value.issue.code == "DEPLOYMENT_CONTRACT_ADMISSION_FAILED"

    unknown = SimpleNamespace(component="unknown", identity=manifest.identity)
    with pytest.raises(DeploymentError) as unsupported:
        resolver.resolve(tmp_path.resolve(), unknown)
    assert unsupported.value.issue.code == "DEPLOYMENT_CONTRACT_ADMISSION_FAILED"

    monkeypatch.setattr(
        resolver,
        "_resolve_platform",
        lambda *_args: (_ for _ in ()).throw(RuntimeError(str(tmp_path))),
    )
    with pytest.raises(DeploymentError) as failed:
        resolver.resolve(tmp_path.resolve(), manifest)
    assert failed.value.issue.code == "DEPLOYMENT_CONTRACT_ADMISSION_FAILED"
    assert str(tmp_path) not in str(failed.value)


def test_database_observer_uses_the_verified_native_migration_inventory(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    inventory = verify_migration_execution_inventory()
    assert len(inventory.heads) == 1
    _write_database(layout.database, inventory.heads[0])
    state = DeploymentState.initial()

    observation = InventoryDatabaseSchemaObserver(layout).observe(state)

    assert observation.provider_identity == f"sha256-{inventory.contract_sha256}"
    assert observation.state_identity == state.identity
    assert observation.heads == inventory.heads
    assert observation.database_identity.startswith("sha256-")


def test_database_observer_rejects_a_symlink_without_disclosing_its_target(
    tmp_path: Path,
) -> None:
    layout = DeploymentLayout.isolated(tmp_path / "host")
    inventory = verify_migration_execution_inventory()
    private = tmp_path / "private.sqlite"
    _write_database(private, inventory.heads[0])
    layout.database.parent.mkdir(parents=True, exist_ok=True)
    layout.database.symlink_to(private)

    with pytest.raises(DeploymentError) as caught:
        InventoryDatabaseSchemaObserver(layout).observe(DeploymentState.initial())

    assert caught.value.issue.code == "DEPLOYMENT_SCHEMA_OBSERVATION_FAILED"
    assert str(private) not in str(caught.value)


def _platform_wheel_members() -> dict[str, bytes]:
    return {
        "helixweave-2.0.0.dist-info/METADATA": (
            b"Name: helixweave\nVersion: 2.0.0\n\n"
        ),
        native_module.PLATFORM_WHEEL_MEMBER_FRONTEND: b"frontend\n",
        native_module.PLATFORM_WHEEL_MEMBER_MIGRATIONS: b'{"migrations":true}\n',
        native_module.PLATFORM_WHEEL_MEMBER_MIGRATION_ADMISSION: b"anchor source\n",
        native_module.PLATFORM_WHEEL_MEMBER_ENCODE_REFERENCES: (
            b'ENCODE_REFERENCE_BINDING_CONTRACT = "1.0.0"\n'
        ),
        native_module.PLATFORM_WHEEL_MEMBER_BULK_REFERENCES: (
            b'BULK_RNASEQ_REFERENCE_BINDING_CONTRACT = "1.0.0"\n'
        ),
        native_module.PLATFORM_WHEEL_MEMBER_IMPLEMENTATION: b"{}\n",
        native_module.PLATFORM_WHEEL_MEMBER_QUALIFICATION: (
            b'{"runtime":{"identity_sha256":"' + b"a" * 64 + b'"}}\n'
        ),
    }


def _raw_wheel(
    members: dict[str, bytes],
    *,
    record_content: bytes | None,
    extra_records: tuple[tuple[str, bytes], ...] = (),
    modes: dict[str, int] | None = None,
) -> bytes:
    archive_bytes = io.BytesIO()
    with zipfile.ZipFile(archive_bytes, mode="w") as archive:
        for name, content in members.items():
            info = zipfile.ZipInfo(name)
            info.external_attr = ((modes or {}).get(name, 0o100644)) << 16
            archive.writestr(info, content)
        if record_content is not None:
            archive.writestr(
                "helixweave-2.0.0.dist-info/RECORD",
                record_content,
            )
        for name, content in extra_records:
            archive.writestr(name, content)
    return archive_bytes.getvalue()


def _record_row(name: str, content: bytes) -> str:
    digest = base64.urlsafe_b64encode(sha256(content).digest()).rstrip(b"=")
    return f"{name},sha256={digest.decode('ascii')},{len(content)}\n"


@pytest.mark.parametrize(
    "fault",
    (
        "empty",
        "unsafe-path",
        "symlink",
        "missing-record",
        "multiple-records",
        "record-encoding",
        "record-width",
        "unknown-member",
        "record-self-hash",
        "unsupported-hash",
        "digest-mismatch",
        "undeclared-member",
    ),
)
def test_platform_wheel_rejects_untrusted_record_inventory(fault: str) -> None:
    member_name = "encode_pipeline/__init__.py"
    members = {member_name: b"VALUE = 1\n"}
    record_name = "helixweave-2.0.0.dist-info/RECORD"
    record = _record_row(member_name, members[member_name]) + f"{record_name},,\n"
    extra_records: tuple[tuple[str, bytes], ...] = ()
    modes: dict[str, int] = {}
    if fault == "empty":
        members = {}
        record_content = None
    elif fault == "unsafe-path":
        members = {"../escape.py": b"escape\n"}
        record_content = b"placeholder"
    elif fault == "symlink":
        modes[member_name] = 0o120777
        record_content = record.encode()
    elif fault == "missing-record":
        record_content = None
    elif fault == "multiple-records":
        record_content = record.encode()
        extra_records = (("other-1.0.dist-info/RECORD", b""),)
    elif fault == "record-encoding":
        record_content = b"\xff"
    elif fault == "record-width":
        record_content = f"{member_name},sha256=x\n".encode()
    elif fault == "unknown-member":
        record_content = b"missing.py,sha256=x,1\n"
    elif fault == "record-self-hash":
        record_content = f"{record_name},sha256=x,1\n".encode()
    elif fault == "unsupported-hash":
        record_content = f"{member_name},md5=x,{len(members[member_name])}\n".encode()
    elif fault == "digest-mismatch":
        record_content = (
            f"{member_name},sha256={'a' * 43},{len(members[member_name])}\n"
        ).encode()
    else:
        record_content = f"{record_name},,\n".encode()

    wheel = _raw_wheel(
        members,
        record_content=record_content,
        extra_records=extra_records,
        modes=modes,
    )
    with pytest.raises(native_module._NativeContractFault) as caught:
        native_module.verify_platform_wheel(wheel)
    assert str(caught.value) == ""


@pytest.mark.parametrize(
    "fault",
    (
        "missing-frontend",
        "duplicate-metadata",
        "metadata-name",
        "implementation",
        "qualification",
        "runtime-identity",
        "reference-source",
    ),
)
def test_platform_wheel_rejects_semantically_invalid_native_contracts(
    monkeypatch: pytest.MonkeyPatch,
    fault: str,
) -> None:
    members = _platform_wheel_members()
    if fault == "missing-frontend":
        members.pop(native_module.PLATFORM_WHEEL_MEMBER_FRONTEND)
    elif fault == "duplicate-metadata":
        members["other-2.0.0.dist-info/METADATA"] = b"Name: other\nVersion: 2.0.0\n\n"
    elif fault == "metadata-name":
        members["helixweave-2.0.0.dist-info/METADATA"] = (
            b"Name: attacker\nVersion: 2.0.0\n\n"
        )
    elif fault == "runtime-identity":
        members[native_module.PLATFORM_WHEEL_MEMBER_QUALIFICATION] = (
            b'{"runtime":{"identity_sha256":"short"}}\n'
        )
    elif fault == "reference-source":
        members[native_module.PLATFORM_WHEEL_MEMBER_ENCODE_REFERENCES] = (
            b'ENCODE_REFERENCE_BINDING_CONTRACT = str("1.0.0")\n'
        )

    anchor = MigrationInventoryTrustAnchor(123, "b" * 64)
    monkeypatch.setattr(
        native_module,
        "migration_inventory_trust_anchor_from_source",
        lambda _content: anchor,
    )
    monkeypatch.setattr(
        native_module,
        "_verify_wheel_migrations",
        lambda *_args, **_kwargs: SimpleNamespace(heads=("head",)),
    )
    monkeypatch.setattr(
        native_module,
        "verify_execution_implementation",
        lambda **_kwargs: SimpleNamespace(
            is_failure=fault == "implementation",
            value=object(),
        ),
    )
    monkeypatch.setattr(
        native_module,
        "load_default_execution_qualification",
        lambda *_args, **_kwargs: SimpleNamespace(is_failure=fault == "qualification"),
    )
    monkeypatch.setattr(
        native_module,
        "parse_manifest_bytes",
        lambda _content: SimpleNamespace(identity=f"sha256-{'c' * 64}"),
    )

    with pytest.raises(native_module._NativeContractFault):
        native_module.verify_platform_wheel(_wheel_bytes(members))


def test_platform_wheel_uses_its_record_owned_migration_trust_anchor(
    monkeypatch,
) -> None:
    admission = b"candidate migration admission source\n"
    inventory = b'{"candidate":"migration inventory"}\n'
    members = {
        "helixweave-2.0.0.dist-info/METADATA": (
            b"Name: helixweave\nVersion: 2.0.0\n\n"
        ),
        native_module.PLATFORM_WHEEL_MEMBER_FRONTEND: b"frontend\n",
        native_module.PLATFORM_WHEEL_MEMBER_MIGRATIONS: inventory,
        native_module.PLATFORM_WHEEL_MEMBER_MIGRATION_ADMISSION: admission,
        native_module.PLATFORM_WHEEL_MEMBER_ENCODE_REFERENCES: (
            b'ENCODE_REFERENCE_BINDING_CONTRACT = "1.0.0"\n'
        ),
        native_module.PLATFORM_WHEEL_MEMBER_BULK_REFERENCES: (
            b'BULK_RNASEQ_REFERENCE_BINDING_CONTRACT = "1.0.0"\n'
        ),
        native_module.PLATFORM_WHEEL_MEMBER_IMPLEMENTATION: b"{}\n",
        native_module.PLATFORM_WHEEL_MEMBER_QUALIFICATION: (
            b'{"runtime":{"identity_sha256":"' + b"a" * 64 + b'"}}\n'
        ),
    }
    anchor = MigrationInventoryTrustAnchor(123, "b" * 64)
    captured: dict[str, object] = {}

    def parse_anchor(content: bytes) -> MigrationInventoryTrustAnchor:
        captured["admission"] = content
        return anchor

    def verify_migrations(
        archive: zipfile.ZipFile,
        content: bytes,
        *,
        trust_anchor: MigrationInventoryTrustAnchor,
    ) -> object:
        assert archive.getinfo(native_module.PLATFORM_WHEEL_MEMBER_MIGRATIONS)
        captured["inventory"] = content
        captured["anchor"] = trust_anchor
        return SimpleNamespace(heads=("candidate-head",))

    monkeypatch.setattr(
        native_module,
        "migration_inventory_trust_anchor_from_source",
        parse_anchor,
    )
    monkeypatch.setattr(native_module, "_verify_wheel_migrations", verify_migrations)
    monkeypatch.setattr(
        native_module,
        "verify_execution_implementation",
        lambda **_kwargs: SimpleNamespace(is_failure=False, value=object()),
    )
    monkeypatch.setattr(
        native_module,
        "load_default_execution_qualification",
        lambda *_args, **_kwargs: SimpleNamespace(is_failure=False),
    )
    monkeypatch.setattr(
        native_module,
        "parse_manifest_bytes",
        lambda _content: SimpleNamespace(identity=f"sha256-{'c' * 64}"),
    )

    verified = native_module.verify_platform_wheel(_wheel_bytes(members))

    assert verified.migration.heads == ("candidate-head",)
    assert captured == {
        "admission": admission,
        "inventory": inventory,
        "anchor": anchor,
    }


def _wheel_bytes(members: dict[str, bytes]) -> bytes:
    record = "helixweave-2.0.0.dist-info/RECORD"
    rows = []
    for name, content in members.items():
        digest = base64.urlsafe_b64encode(sha256(content).digest()).rstrip(b"=")
        rows.append((name, f"sha256={digest.decode('ascii')}", str(len(content))))
    rows.append((record, "", ""))
    stream = io.StringIO(newline="")
    csv.writer(stream, lineterminator="\n").writerows(rows)
    archive_bytes = io.BytesIO()
    with zipfile.ZipFile(archive_bytes, mode="w") as archive:
        for name, content in members.items():
            archive.writestr(name, content)
        archive.writestr(record, stream.getvalue().encode("utf-8"))
    return archive_bytes.getvalue()
