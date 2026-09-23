"""Validation-only real API fixture; no queue worker or admitted scientific runtime.

The tiny reference supports the real private binding verifier. It is not a
scientific fixture. Callers own the temporary root and must close the app.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


def prepare_reference(root: Path) -> Path:
    root.mkdir(parents=True, exist_ok=True)

    def resource(name: str, content: bytes):
        path = root / name
        path.write_bytes(content)
        return str(path), hashlib.sha256(content).hexdigest()

    fasta, fasta_sha = resource("genome.fa", b">chr1\nACGT\n")
    gtf, gtf_sha = resource(
        "genes.gtf", b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n'
    )
    transcript, transcript_sha = resource("transcripts.fa", b">tx1\nACGT\n")
    binding = {
        "schema_version": "bulk-rnaseq-reference-binding-v1",
        "reference": {
            "reference_id": "tiny-ref",
            "fasta": fasta,
            "fasta_sha256": fasta_sha,
            "gtf": gtf,
            "gtf_sha256": gtf_sha,
            "annotation_style": "ensembl",
        },
        "transcriptome": {
            "reference_id": "tiny-ref",
            "fasta_sha256": fasta_sha,
            "gtf_sha256": gtf_sha,
            "transcript_fasta": transcript,
            "transcript_fasta_sha256": transcript_sha,
        },
    }
    path = root / "reference-profiles.json"
    path.write_text(
        json.dumps(
            {
                "schema_version": "helixweave-reference-profiles-v1",
                "profiles": {"tiny": {"bindings": {"bulk-rnaseq": binding}}},
            }
        ),
        encoding="utf-8",
    )
    return path


def create_validation_app(root: Path):
    # The caller explicitly configures the private reference file before composition.
    from encode_pipeline.adapters.bulk_rnaseq import (
        BulkRnaSeqExecutionBinding,
        BulkRnaSeqWorkflowAdapter,
        RuntimeAssetBinding,
    )
    from encode_pipeline.adapters.bulk_rnaseq.execution_identity import (
        ExecutionImplementationQualification,
        verify_execution_implementation,
    )
    from encode_pipeline.api.main import create_app
    from encode_pipeline.platform.registry import WorkflowRegistry

    verified = verify_execution_implementation()
    assert verified.is_success, verified.issues
    binding = BulkRnaSeqExecutionBinding(
        assets=RuntimeAssetBinding(root=root / "absent-runtime"),
        implementation_qualification=ExecutionImplementationQualification.from_verified(
            verified.value
        ),
    )
    app = create_app(
        database_url=f"sqlite:///{root / 'platform.db'}",
        workspace_root=root / "workspaces",
        project_root=root,
        registry=WorkflowRegistry([BulkRnaSeqWorkflowAdapter(execution=binding)]),
    )
    try:
        summary = app.state.reference_profile_service.register(
            safe_key="tiny",
            display_name="Tiny validation reference",
            organism="Synthetic",
            assembly="tiny",
            config_key="tiny",
        )
        app.state.test_reference_profile = app.state.reference_profile_service.enable(
            summary.profile_id,
            revision_id=summary.revision_id,
        )
        app.state.account_administration_service.bootstrap_initial_administrator(
            "e2e-admin",
            "e2e playwright admin password",
        )
        return app
    except BaseException:
        app.state.run_queue.close()
        app.state.persistence.close()
        raise
