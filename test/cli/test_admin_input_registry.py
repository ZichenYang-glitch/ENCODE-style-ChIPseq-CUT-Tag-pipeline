"""Administrator-only StoragePool and InputFile CLI boundaries."""

from __future__ import annotations

from contextlib import contextmanager
import json
from pathlib import Path

from encode_pipeline.cli import admin


DATABASE_URL = "sqlite:////tmp/helixweave-input-registry-test.db"


class RecordingInputRegistry:
    def __init__(self) -> None:
        self.calls: list[tuple[str, dict[str, object]]] = []

    def register_storage_pool(
        self,
        *,
        display_name: str,
        config_key: str,
    ) -> dict[str, str]:
        self.calls.append(
            (
                "register_storage_pool",
                {"display_name": display_name, "config_key": config_key},
            )
        )
        return {"storage_pool_id": "stgp_" + "1" * 32}

    def bind_project_storage_pool(
        self,
        *,
        project_id: str,
        storage_pool_id: str,
    ) -> dict[str, str]:
        self.calls.append(
            (
                "bind_project_storage_pool",
                {
                    "project_id": project_id,
                    "storage_pool_id": storage_pool_id,
                },
            )
        )
        return {
            "project_id": project_id,
            "storage_pool_id": storage_pool_id,
        }

    def register_input_file(
        self,
        *,
        project_id: str,
        stable_key: str,
        pool_relative_path: str,
    ) -> dict[str, str]:
        self.calls.append(
            (
                "register_input_file",
                {
                    "project_id": project_id,
                    "stable_key": stable_key,
                    "pool_relative_path": pool_relative_path,
                },
            )
        )
        return {"input_file_revision_id": "inpfr_" + "2" * 32}


def _factory(
    registry: RecordingInputRegistry,
    observed: list[tuple[str, Path | None]],
):
    @contextmanager
    def factory(database_url: str, storage_pool_config: Path | None):
        observed.append((database_url, storage_pool_config))
        yield registry

    return factory


def test_storage_pool_registration_requires_explicit_private_config(
    tmp_path,
    capsys,
) -> None:
    config_path = tmp_path / "operator-pools.json"
    registry = RecordingInputRegistry()
    observed: list[tuple[str, Path | None]] = []

    exit_code = admin.main(
        [
            "--database-url",
            DATABASE_URL,
            "--storage-pool-config",
            str(config_path),
            "storage-pool",
            "register",
            "--display-name",
            "Approved ingress",
            "--config-key",
            "ingress-primary",
        ],
        input_registry_factory=_factory(registry, observed),
    )

    assert exit_code == 0
    assert observed == [(DATABASE_URL, config_path)]
    assert registry.calls == [
        (
            "register_storage_pool",
            {
                "display_name": "Approved ingress",
                "config_key": "ingress-primary",
            },
        )
    ]
    assert json.loads(capsys.readouterr().out) == {
        "storage_pool_id": "stgp_" + "1" * 32
    }


def test_project_pool_binding_never_accepts_a_root_or_path(
    tmp_path,
    capsys,
) -> None:
    registry = RecordingInputRegistry()
    observed: list[tuple[str, Path | None]] = []
    project_id = "prj_" + "1" * 32
    pool_id = "stgp_" + "2" * 32

    exit_code = admin.main(
        [
            "--database-url",
            DATABASE_URL,
            "project",
            "bind-storage-pool",
            "--project-id",
            project_id,
            "--storage-pool-id",
            pool_id,
        ],
        input_registry_factory=_factory(registry, observed),
    )

    assert exit_code == 0
    assert observed == [(DATABASE_URL, None)]
    assert registry.calls[-1] == (
        "bind_project_storage_pool",
        {"project_id": project_id, "storage_pool_id": pool_id},
    )
    output = capsys.readouterr().out
    assert str(tmp_path) not in output
    assert "config_key" not in output


def test_input_file_registration_accepts_only_pool_relative_path(
    tmp_path,
    capsys,
) -> None:
    config_path = tmp_path / "operator-pools.json"
    registry = RecordingInputRegistry()
    observed: list[tuple[str, Path | None]] = []
    project_id = "prj_" + "1" * 32

    exit_code = admin.main(
        [
            "--database-url",
            DATABASE_URL,
            "--storage-pool-config",
            str(config_path),
            "input-file",
            "register",
            "--project-id",
            project_id,
            "--stable-key",
            "donor-01-r1",
            "--pool-relative-path",
            "reads/donor-01.fastq.gz",
        ],
        input_registry_factory=_factory(registry, observed),
    )

    assert exit_code == 0
    assert registry.calls[-1] == (
        "register_input_file",
        {
            "project_id": project_id,
            "stable_key": "donor-01-r1",
            "pool_relative_path": "reads/donor-01.fastq.gz",
        },
    )
    output = capsys.readouterr().out
    assert str(tmp_path) not in output
    assert "reads/donor-01.fastq.gz" not in output
