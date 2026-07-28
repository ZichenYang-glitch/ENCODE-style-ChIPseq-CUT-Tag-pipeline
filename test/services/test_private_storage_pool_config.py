"""Tests for operator-private StoragePool root configuration."""

from __future__ import annotations

from collections import Counter
import json
import os
from pathlib import Path

import pytest

from encode_pipeline.services import private_storage_pools as storage_pool_module
from encode_pipeline.services.private_storage_pools import (
    PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
    PrivateStoragePoolConfig,
    PrivateStoragePoolConfigError,
    load_private_storage_pool_config,
)


def test_private_config_loads_bounded_opaque_key_to_absolute_directory(
    tmp_path,
) -> None:
    root = tmp_path / "private-inputs"
    root.mkdir()
    config_path = tmp_path / "operator-storage.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"primary-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )

    config = load_private_storage_pool_config(config_path)

    assert config.config_keys == ("primary-inputs",)
    assert config.root_for("primary-inputs") == root
    assert str(root) not in repr(config)
    assert str(config_path) not in repr(config)


@pytest.mark.parametrize(
    "document",
    (
        '{"schema_version":"storage-pool-roots-v1",'
        '"schema_version":"storage-pool-roots-v1","storage_pool_roots":{}}',
        '{"schema_version":"storage-pool-roots-v1","storage_pool_roots":'
        '{"primary":"/private/a","primary":"/private/b"}}',
        '{"schema_version":"storage-pool-roots-v1",'
        '"storage_pool_roots":{},"extra":true}',
        '{"schema_version":"future","storage_pool_roots":{}}',
    ),
)
def test_private_config_rejects_duplicates_extra_fields_and_wrong_version(
    tmp_path,
    document: str,
) -> None:
    config_path = tmp_path / "operator-storage.json"
    config_path.write_text(document, encoding="utf-8")

    with pytest.raises(PrivateStoragePoolConfigError):
        load_private_storage_pool_config(config_path)


def test_private_config_rejects_relative_or_missing_roots_without_disclosure(
    tmp_path,
) -> None:
    secret_root = tmp_path / "sensitive" / "missing"
    config_path = tmp_path / "secret-name.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {
                    "relative": "../escape",
                    "missing": str(secret_root),
                },
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(PrivateStoragePoolConfigError) as captured:
        PrivateStoragePoolConfig.from_file(config_path)

    rendered = f"{captured.value!s} {captured.value!r}"
    assert str(secret_root) not in rendered
    assert str(config_path) not in rendered
    assert "../escape" not in rendered


def test_private_config_unknown_key_is_stable_and_does_not_list_roots(
    tmp_path,
) -> None:
    root = tmp_path / "private-inputs"
    root.mkdir()
    config_path = tmp_path / "operator-storage.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"primary-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )
    config = PrivateStoragePoolConfig.from_file(config_path)

    with pytest.raises(PrivateStoragePoolConfigError) as captured:
        config.root_for("not-configured")

    rendered = f"{captured.value!s} {captured.value!r}"
    assert str(root) not in rendered
    assert "primary-inputs" not in rendered


def test_private_config_rejects_symlink_config_file(tmp_path) -> None:
    root = tmp_path / "private-inputs"
    root.mkdir()
    real_config = tmp_path / "real-storage.json"
    real_config.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"primary-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )
    linked_config = tmp_path / "linked-storage.json"
    linked_config.symlink_to(real_config)

    with pytest.raises(PrivateStoragePoolConfigError):
        load_private_storage_pool_config(linked_config)


def test_private_config_rejects_fifo_without_blocking(tmp_path) -> None:
    fifo = tmp_path / "operator-storage.fifo"
    os.mkfifo(fifo)

    with pytest.raises(PrivateStoragePoolConfigError):
        load_private_storage_pool_config(fifo)


def test_private_config_rejects_symlink_storage_root(tmp_path) -> None:
    real_root = tmp_path / "real-inputs"
    real_root.mkdir()
    linked_root = tmp_path / "linked-inputs"
    linked_root.symlink_to(real_root, target_is_directory=True)
    config_path = tmp_path / "operator-storage.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"primary-inputs": str(linked_root)},
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(PrivateStoragePoolConfigError):
        load_private_storage_pool_config(config_path)


def test_private_config_rejects_double_slash_config_path(tmp_path) -> None:
    root = tmp_path / "private-inputs"
    root.mkdir()
    config_path = tmp_path / "operator-storage.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"primary-inputs": str(root)},
            }
        ),
        encoding="utf-8",
    )
    ambiguous_path = Path(f"/{config_path}")
    assert os.fspath(ambiguous_path).startswith("//")

    with pytest.raises(PrivateStoragePoolConfigError):
        load_private_storage_pool_config(ambiguous_path)


@pytest.mark.parametrize("ambiguous_kind", ("double-prefix", "empty-component"))
def test_private_config_rejects_noncanonical_absolute_root(
    tmp_path,
    ambiguous_kind: str,
) -> None:
    root = tmp_path / "private-inputs"
    root.mkdir()
    raw_root = str(root)
    if ambiguous_kind == "double-prefix":
        raw_root = f"/{raw_root}"
    else:
        parent, leaf = raw_root.rsplit("/", 1)
        raw_root = f"{parent}//{leaf}"
    config_path = tmp_path / "operator-storage.json"
    config_path.write_text(
        json.dumps(
            {
                "schema_version": PRIVATE_STORAGE_POOL_SCHEMA_VERSION,
                "storage_pool_roots": {"primary-inputs": raw_root},
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(PrivateStoragePoolConfigError):
        load_private_storage_pool_config(config_path)


def test_private_config_closes_leaf_when_initial_fstat_fails(
    tmp_path,
    monkeypatch,
) -> None:
    config_path = tmp_path / "operator-storage.json"
    config_path.write_text("{}", encoding="utf-8")
    real_open = os.open
    real_close = os.close
    opened: list[int] = []
    closed: list[int] = []

    def tracking_open(*args, **kwargs):
        descriptor = real_open(*args, **kwargs)
        opened.append(descriptor)
        return descriptor

    def tracking_close(descriptor: int) -> None:
        closed.append(descriptor)
        real_close(descriptor)

    def failing_fstat(_descriptor: int):
        raise OSError("simulated fstat failure")

    monkeypatch.setattr(storage_pool_module.os, "open", tracking_open)
    monkeypatch.setattr(storage_pool_module.os, "close", tracking_close)
    monkeypatch.setattr(storage_pool_module.os, "fstat", failing_fstat)

    with pytest.raises(PrivateStoragePoolConfigError):
        load_private_storage_pool_config(config_path)

    assert opened
    assert Counter(opened) == Counter(closed)
