#!/usr/bin/env python3
"""Regenerate the committed bulk RNA-seq execution implementation manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
import sys
import tempfile


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (  # noqa: E402
    EXECUTION_IMPLEMENTATION_MANIFEST_FILE,
    EXECUTION_PERSISTENCE_CONTRACT_PATH,
    build_execution_implementation_manifest,
    canonical_execution_manifest_bytes,
    verify_execution_implementation,
)
from encode_pipeline.adapters.bulk_rnaseq.qualification import (  # noqa: E402
    DEFAULT_EXECUTION_QUALIFICATION_FILE,
    canonical_qualification_bytes,
    qualification_document,
)
from encode_pipeline.adapters.bulk_rnaseq.persistence_projection import (  # noqa: E402
    parse_schema_projection_spec,
    schema_projection_sha256,
    schema_projection_spec_document,
    with_schema_projection_sha256,
)
from encode_pipeline.persistence import create_database_engine  # noqa: E402
from encode_pipeline.persistence.migrations import upgrade_database  # noqa: E402


def main() -> int:
    """Write canonical bytes and print non-sensitive identity evidence."""
    projection_sha256 = _refresh_persistence_projection()
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
    print(f"schema_projection_sha256={projection_sha256}")
    print(f"qualification_sha256={hashlib.sha256(qualification_content).hexdigest()}")
    print(
        "qualification_record_sha256="
        f"{qualification_document(implementation.value)['record_sha256']}"
    )
    return 0


def _refresh_persistence_projection() -> str:
    contract_path = PROJECT_ROOT.joinpath(
        *EXECUTION_PERSISTENCE_CONTRACT_PATH.split("/")
    )
    document = json.loads(
        contract_path.read_text(encoding="utf-8"),
        object_pairs_hook=_unique_json_object,
    )
    if not isinstance(document, dict) or "schema_projection" not in document:
        raise ValueError("execution persistence contract is invalid")
    spec = parse_schema_projection_spec(document["schema_projection"])
    with tempfile.TemporaryDirectory(
        prefix="helixweave-bulk-persistence-projection-"
    ) as temporary:
        database_url = f"sqlite:///{Path(temporary) / 'projection.db'}"
        upgrade_database(database_url)
        engine = create_database_engine(database_url)
        try:
            digest = schema_projection_sha256(engine, spec)
        finally:
            engine.dispose()
    document["schema_projection"] = schema_projection_spec_document(
        with_schema_projection_sha256(spec, digest)
    )
    contract_path.write_bytes(
        (
            json.dumps(
                document,
                ensure_ascii=False,
                sort_keys=True,
                separators=(",", ":"),
            )
            + "\n"
        ).encode("utf-8")
    )
    return digest


def _unique_json_object(pairs: list[tuple[str, object]]) -> dict[str, object]:
    result: dict[str, object] = {}
    for key, value in pairs:
        if key in result:
            raise ValueError("execution persistence contract contains duplicate fields")
        result[key] = value
    return result


if __name__ == "__main__":
    raise SystemExit(main())
