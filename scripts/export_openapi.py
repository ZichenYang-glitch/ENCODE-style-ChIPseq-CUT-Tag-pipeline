"""Export the FastAPI OpenAPI schema to a JSON file."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from encode_pipeline.api.main import create_app


def export_openapi(output_path: Path) -> None:
    """Write the current FastAPI OpenAPI schema to *output_path*."""
    app = create_app()
    schema = app.openapi()
    output_path.write_text(json.dumps(schema, indent=2) + "\n", encoding="utf-8")


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
