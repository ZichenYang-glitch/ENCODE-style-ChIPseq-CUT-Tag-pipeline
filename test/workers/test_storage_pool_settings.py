"""Private StoragePool configuration path shared by API and workers."""

from __future__ import annotations

from encode_pipeline.workers.settings import (
    STORAGE_POOL_CONFIG_ENV,
    load_worker_settings,
)


def test_worker_settings_reads_optional_private_pool_config_without_repr_leak(
    tmp_path,
) -> None:
    private_config = tmp_path / "operator-storage-pools.json"

    configured = load_worker_settings({STORAGE_POOL_CONFIG_ENV: str(private_config)})

    assert configured.storage_pool_config == private_config
    assert str(private_config) not in repr(configured)


def test_worker_settings_leave_pool_config_unset_for_legacy_compatibility() -> None:
    configured = load_worker_settings({})

    assert configured.storage_pool_config is None
