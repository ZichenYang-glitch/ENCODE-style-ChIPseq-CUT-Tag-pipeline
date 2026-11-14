#!/usr/bin/env python3
"""Seed one user project and sample into the running browser runtime database."""

from __future__ import annotations

import json
import os
from pathlib import Path

from encode_pipeline.persistence import create_session_factory
from encode_pipeline.persistence.data_registry import SqlAlchemyDataRegistryRepository
from encode_pipeline.persistence.database import create_database_engine
from encode_pipeline.persistence.input_registry import SqlAlchemyInputRegistryRepository
from encode_pipeline.services.data_registry import DataRegistryService, SampleImportSpec
from encode_pipeline.services.input_registry import (
    InputRegistryAdminService,
    InputRegistryService,
)
from encode_pipeline.services.private_storage_pools import (
    load_private_storage_pool_config,
)

PROJECT_NAME = "Browser authoring"
SAMPLE_NAME = "Browser sample"
POOL_KEY = "authoring-inputs"
READ_KEY = "browser-sample-r1"
READ_RELATIVE = "reads/s1.fastq.gz"


def main() -> None:
    root = Path(os.environ["ENCODE_PIPELINE_E2E_ROOT"])
    engine = create_database_engine(f"sqlite:///{root / 'platform.db'}")
    try:
        registry = DataRegistryService(
            repository=SqlAlchemyDataRegistryRepository(create_session_factory(engine))
        )
        project = next(
            (
                item
                for item in registry.list_projects()
                if item.display_name == PROJECT_NAME
            ),
            None,
        )
        if project is None:
            project = registry.create_project(PROJECT_NAME)
        if not registry.list_samples(project.project_id):
            registry.import_samples(
                project.project_id,
                rows=(
                    SampleImportSpec(
                        stable_key="browser-sample",
                        display_name=SAMPLE_NAME,
                        attributes={},
                    ),
                ),
            )
        _ensure_admitted_read(engine, project.project_id, root)
    finally:
        engine.dispose()
    print(json.dumps({"project": PROJECT_NAME, "sample": SAMPLE_NAME}))


def _ensure_admitted_read(engine, project_id: str, root: Path) -> None:
    manifest = json.loads((root / "runtime.json").read_text(encoding="utf-8"))
    config_path = Path(manifest["storagePoolConfig"])
    session_factory = create_session_factory(engine)
    repository = SqlAlchemyInputRegistryRepository(session_factory)
    service = InputRegistryService(repository=repository)
    admin = InputRegistryAdminService(
        service,
        private_config=load_private_storage_pool_config(config_path),
    )
    pool = next(
        (item for item in service.list_storage_pools() if item.config_key == POOL_KEY),
        None,
    )
    if pool is None:
        pool = admin.register_storage_pool(
            display_name="Authoring inputs",
            config_key=POOL_KEY,
        )
    try:
        service.get_active_project_storage_pool_binding(project_id)
    except KeyError:
        admin.bind_project_storage_pool(
            project_id=project_id,
            storage_pool_id=pool.storage_pool_id,
        )
    files = repository.list_input_files(project_id)
    if any(item.stable_key == READ_KEY for item in files):
        return
    admin.register_input_file(
        project_id=project_id,
        stable_key=READ_KEY,
        pool_relative_path=READ_RELATIVE,
    )


if __name__ == "__main__":
    main()
