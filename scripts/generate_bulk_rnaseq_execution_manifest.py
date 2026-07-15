#!/usr/bin/env python3
"""Regenerate the committed bulk RNA-seq execution implementation manifest."""

from __future__ import annotations

import hashlib
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (  # noqa: E402
    EXECUTION_IMPLEMENTATION_MANIFEST_FILE,
    build_execution_implementation_manifest,
    canonical_execution_manifest_bytes,
)


def main() -> int:
    """Write canonical bytes and print non-sensitive identity evidence."""
    manifest = build_execution_implementation_manifest(PROJECT_ROOT)
    content = canonical_execution_manifest_bytes(manifest)
    destination = (
        PROJECT_ROOT
        / "src/encode_pipeline/contracts/nfcore_rnaseq"
        / EXECUTION_IMPLEMENTATION_MANIFEST_FILE
    )
    destination.write_bytes(content)
    print(f"files={manifest['file_count']}")
    print(f"aggregate_sha256={manifest['aggregate_sha256']}")
    print(f"manifest_sha256={hashlib.sha256(content).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
