"""Export the FastAPI OpenAPI schema to a JSON file."""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from unittest.mock import patch

from encode_pipeline.api.main import create_app
from encode_pipeline.workers.settings import (
    DEFAULT_JOB_TIMEOUT_SECONDS,
    DEFAULT_REDIS_API_READ_TIMEOUT_SECONDS,
    DEFAULT_REDIS_CONNECT_TIMEOUT_SECONDS,
    JOB_TIMEOUT_SECONDS_ENV,
    QUEUE_NAME_ENV,
    REDIS_API_READ_TIMEOUT_SECONDS_ENV,
    REDIS_CONNECT_TIMEOUT_SECONDS_ENV,
    REDIS_URL_ENV,
)


OPENAPI_REDIS_URL = "redis://127.0.0.1:6379/0"
OPENAPI_QUEUE_NAME = "encode-pipeline-openapi-export"


def export_openapi(output_path: Path) -> None:
    """Write the OpenAPI schema without opening the configured runtime state."""
    with TemporaryDirectory(prefix="encode-pipeline-openapi-") as temporary_dir:
        runtime_root = Path(temporary_dir)
        isolated_worker_environment = {
            REDIS_URL_ENV: OPENAPI_REDIS_URL,
            QUEUE_NAME_ENV: OPENAPI_QUEUE_NAME,
            JOB_TIMEOUT_SECONDS_ENV: str(DEFAULT_JOB_TIMEOUT_SECONDS),
            REDIS_CONNECT_TIMEOUT_SECONDS_ENV: str(
                DEFAULT_REDIS_CONNECT_TIMEOUT_SECONDS
            ),
            REDIS_API_READ_TIMEOUT_SECONDS_ENV: str(
                DEFAULT_REDIS_API_READ_TIMEOUT_SECONDS
            ),
        }
        with patch.dict(os.environ, isolated_worker_environment):
            app = create_app(
                database_url=f"sqlite:///{runtime_root / 'platform.db'}",
                workspace_root=runtime_root / "workspaces",
            )
        try:
            schema = app.openapi()
            output_path.write_text(
                json.dumps(schema, indent=2) + "\n",
                encoding="utf-8",
            )
        finally:
            try:
                app.state.run_queue.close()
            finally:
                app.state.persistence.close()


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Export FastAPI OpenAPI schema to JSON."
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("frontend/openapi.json"),
        help="Output JSON path (default: frontend/openapi.json)",
    )
    args = parser.parse_args()
    export_openapi(args.output)


if __name__ == "__main__":
    main()
