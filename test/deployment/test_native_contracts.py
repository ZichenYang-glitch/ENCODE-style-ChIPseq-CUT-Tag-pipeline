from __future__ import annotations

import base64
import csv
from hashlib import sha256
import io
import sqlite3
from pathlib import Path
from types import SimpleNamespace
import zipfile

import pytest

import encode_pipeline.deployment.native_contracts as native_module
from encode_pipeline.deployment.errors import DeploymentError
from encode_pipeline.deployment.layout import DeploymentLayout
from encode_pipeline.deployment.models import DeploymentState
from encode_pipeline.deployment.native_contracts import (
    InventoryDatabaseSchemaObserver,
    ProductionNativeContractResolver,
    parse_encode_explicit_lock,
)
from encode_pipeline.persistence.migration_admission import (
    MigrationInventoryTrustAnchor,
    verify_migration_execution_inventory,
)


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


def test_bulk_native_admission_defaults_to_the_static_runtime_closure() -> None:
    from encode_pipeline.adapters.bulk_rnaseq.runtime_assets import (
        verify_runtime_asset_closure,
    )

    resolver = ProductionNativeContractResolver()

    assert resolver._bulk_runtime_verifier is verify_runtime_asset_closure


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
