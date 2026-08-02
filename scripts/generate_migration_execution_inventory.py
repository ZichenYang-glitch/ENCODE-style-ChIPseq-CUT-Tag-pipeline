#!/usr/bin/env python3
"""Regenerate the source-owned Alembic migration execution inventory."""

from __future__ import annotations

import hashlib
from pathlib import Path
import sys


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from encode_pipeline.persistence.migration_admission import (  # noqa: E402
    MIGRATION_EXECUTION_INVENTORY_FILE,
    build_migration_execution_inventory,
    canonical_migration_execution_inventory_bytes,
)


def main() -> int:
    """Write canonical bytes and print coordinates for source-anchor review."""

    persistence_root = PROJECT_ROOT / "src/encode_pipeline/persistence"
    document = build_migration_execution_inventory(persistence_root)
    content = canonical_migration_execution_inventory_bytes(document)
    destination = persistence_root / MIGRATION_EXECUTION_INVENTORY_FILE
    destination.write_bytes(content)
    print(f"revisions={document['revision_count']}")
    print(f"inventory_identity_sha256={document['inventory_sha256']}")
    print(f"inventory_size_bytes={len(content)}")
    print(f"inventory_file_sha256={hashlib.sha256(content).hexdigest()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
