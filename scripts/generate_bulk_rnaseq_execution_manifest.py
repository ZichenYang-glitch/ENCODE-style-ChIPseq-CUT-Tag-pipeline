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
    verify_execution_implementation,
)
from encode_pipeline.adapters.bulk_rnaseq.qualification import (  # noqa: E402
    DEFAULT_EXECUTION_QUALIFICATION_FILE,
    canonical_qualification_bytes,
    qualification_document,
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
    implementation = verify_execution_implementation(
        manifest_bytes=content,
        package_root=PROJECT_ROOT / "src/encode_pipeline",
    )
    if implementation.is_failure:
        raise ValueError("generated execution implementation could not be verified")
    qualification_content = canonical_qualification_bytes(
        qualification_document(implementation.value)
    )
    qualification_destination = destination.with_name(
        DEFAULT_EXECUTION_QUALIFICATION_FILE
    )
    qualification_destination.write_bytes(qualification_content)
    print(f"files={manifest['file_count']}")
    print(f"aggregate_sha256={manifest['aggregate_sha256']}")
    print(f"manifest_sha256={hashlib.sha256(content).hexdigest()}")
    print(
        "persistence_contract_version="
        f"{implementation.value.persistence_contract.contract_version}"
    )
    print(
        "persistence_contract_sha256="
        f"{implementation.value.persistence_contract.sha256}"
    )
    print(f"qualification_sha256={hashlib.sha256(qualification_content).hexdigest()}")
    print(
        "qualification_record_sha256="
        f"{qualification_document(implementation.value)['record_sha256']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
