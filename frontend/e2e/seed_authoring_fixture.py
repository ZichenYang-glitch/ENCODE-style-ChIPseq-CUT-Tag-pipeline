#!/usr/bin/env python3
"""Seed one user project and sample into the running browser runtime database."""

from __future__ import annotations

import json
import os
from pathlib import Path

from encode_pipeline.persistence import create_session_factory
from encode_pipeline.persistence.data_registry import SqlAlchemyDataRegistryRepository
from encode_pipeline.persistence.database import create_database_engine
from encode_pipeline.services.data_registry import DataRegistryService, SampleImportSpec

PROJECT_NAME = "Browser authoring"
SAMPLE_NAME = "Browser sample"


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
    finally:
        engine.dispose()
    print(json.dumps({"project": PROJECT_NAME, "sample": SAMPLE_NAME}))


if __name__ == "__main__":
    main()
