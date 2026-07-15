"""Contract-only Omics Intake Bundle 0.2 fixtures for consumer tests."""

from __future__ import annotations

import csv
import hashlib
import io
import json
from pathlib import Path

import yaml


SAMPLE_COLUMNS = (
    "sample",
    "fastq_1",
    "fastq_2",
    "layout",
    "assay",
    "target",
    "peak_mode",
    "genome",
    "bowtie2_index",
    "control_bam",
    "role",
    "control_sample",
    "experiment",
    "condition",
    "replicate",
    "biological_replicate",
    "technical_replicate",
)

BT2_STANDARD = (
    "{prefix}.1.bt2",
    "{prefix}.2.bt2",
    "{prefix}.3.bt2",
    "{prefix}.4.bt2",
    "{prefix}.rev.1.bt2",
    "{prefix}.rev.2.bt2",
)
BT2_LARGE = tuple(template.replace(".bt2", ".bt2l") for template in BT2_STANDARD)


def bowtie2_shard_paths(root: Path, *, large: bool = False) -> tuple[Path, ...]:
    """Return one deterministic logical Bowtie2 candidate set."""
    templates = BT2_LARGE if large else BT2_STANDARD
    prefix = "resources/mm10/index/genome"
    return tuple(root / template.format(prefix=prefix) for template in templates)


def configure_bowtie2_shards(
    bundle_path: Path,
    *,
    standard_count: int,
    large_count: int,
) -> None:
    """Replace tiny test shards without importing producer or consumer code."""
    if not 0 <= standard_count <= len(BT2_STANDARD):
        raise ValueError("standard_count is outside the test candidate set")
    if not 0 <= large_count <= len(BT2_LARGE):
        raise ValueError("large_count is outside the test candidate set")
    root = bundle_path.parent
    standard = bowtie2_shard_paths(root)
    large = bowtie2_shard_paths(root, large=True)
    for path in (*standard, *large):
        if path.is_symlink() or path.exists():
            path.unlink()
    for path in (*standard[:standard_count], *large[:large_count]):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"tiny-index\n")


def create_runnable_bundle(root: Path) -> Path:
    """Create a tiny schema-valid Bundle without importing producer code."""
    root.mkdir()
    fastq_relative = "downloads/ena/SRR000001.fastq.gz"
    fastq = root / fastq_relative
    fastq.parent.mkdir(parents=True)
    fastq.write_bytes(b"@read\nACGT\n+\n!!!!\n")

    index_prefix = "resources/mm10/index/genome"
    for template in BT2_STANDARD:
        path = root / template.format(prefix=index_prefix)
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"tiny-index\n")

    sample_row = {
        "sample": "SRR000001.b1.t1",
        "fastq_1": fastq_relative,
        "fastq_2": "",
        "layout": "SE",
        "assay": "chipseq",
        "target": "H3K4me1",
        "peak_mode": "narrow",
        "genome": "mm10",
        "bowtie2_index": index_prefix,
        "control_bam": "",
        "role": "treatment",
        "control_sample": "",
        "experiment": "SRX000001",
        "condition": "H3K4me1",
        "replicate": "1",
        "biological_replicate": "1",
        "technical_replicate": "1",
    }
    samples_buffer = io.StringIO(newline="")
    writer = csv.DictWriter(
        samples_buffer,
        fieldnames=SAMPLE_COLUMNS,
        delimiter="\t",
        lineterminator="\n",
    )
    writer.writeheader()
    writer.writerow(sample_row)
    (root / "samples.encode.tsv").write_text(
        samples_buffer.getvalue(), encoding="utf-8"
    )

    config = {
        "samples": "./samples.encode.tsv",
        "outdir": "results",
        "threads": 4,
        "mapq": 30,
        "binsize": 10,
        "use_control": False,
        "genome_resources": {"mm10": {"effective_genome_size": "mm"}},
    }
    (root / "config.encode.yaml").write_text(
        yaml.safe_dump(config, allow_unicode=True, sort_keys=False),
        encoding="utf-8",
    )
    (root / "intake.json").write_text(
        json.dumps({"schema_version": "0.2"}, indent=2) + "\n",
        encoding="utf-8",
    )

    canonical = _artifact(root, "intake.json", "canonical_project", "application/json")
    sample_artifact = _artifact(
        root,
        "samples.encode.tsv",
        "sample_sheet",
        "text/tab-separated-values",
    )
    config_artifact = _artifact(
        root,
        "config.encode.yaml",
        "workflow_config",
        "text/yaml",
    )
    fastq_bytes = fastq.read_bytes()
    fastq_sha256 = hashlib.sha256(fastq_bytes).hexdigest()
    fingerprint = "sha256:" + "1" * 64
    artifact_fingerprint = "sha256:" + "2" * 64
    payload = {
        "bundle_schema_version": "0.2",
        "producer": {"name": "omics-intake", "version": "0.2.0"},
        "workflow": {
            "name": "encode-epigenomics",
            "genome": "mm10",
            "render_contract": "encode-render-v3",
        },
        "canonical_project": {
            "source_schema_version": "0.2",
            "identity": f"sha256:{canonical['sha256']}",
            "record": canonical,
        },
        "artifacts": [sample_artifact, config_artifact],
        "files": [
            {
                "file_id": "insdc:SRR000001:1",
                "kind": "local",
                "scope": "project",
                "path": fastq_relative,
                "file_format": "fastq",
                "role": "treatment",
                "read_number": 1,
                "size_bytes": len(fastq_bytes),
                "checksum": {
                    "status": "verified",
                    "algorithm": "sha256",
                    "digest": fastq_sha256,
                },
            },
            {
                "file_id": "insdc:SRR000001:1",
                "kind": "planned",
                "scope": "project",
                "path": fastq_relative,
                "file_format": "fastq",
                "role": "treatment",
                "read_number": 1,
                "size_bytes": len(fastq_bytes),
                "checksum": {
                    "status": "declared",
                    "algorithm": "md5",
                    "digest": hashlib.md5(
                        fastq_bytes, usedforsecurity=False
                    ).hexdigest(),
                },
            },
        ],
        "readiness": {
            "export_status": "ready",
            "execution_status": "runnable",
            "validator_status": "passed",
            "strict_validator_status": "passed",
            "unresolved_requirements": [],
            "source_basis_fingerprint": fingerprint,
        },
        "validations": [
            _validation("non_strict", artifact_fingerprint, None),
            _validation("strict", artifact_fingerprint, "sha256:" + "3" * 64),
        ],
        "issues": [],
    }
    bundle_path = root / "intake-bundle.json"
    write_bundle(bundle_path, payload)
    return bundle_path


def read_bundle(bundle_path: Path) -> dict:
    return json.loads(bundle_path.read_text(encoding="utf-8"))


def write_bundle(bundle_path: Path, payload: dict) -> None:
    bundle_path.write_text(
        json.dumps(payload, indent=2, sort_keys=False) + "\n",
        encoding="utf-8",
    )


def refresh_artifact(bundle_path: Path, role: str) -> None:
    payload = read_bundle(bundle_path)
    artifact = next(item for item in payload["artifacts"] if item["role"] == role)
    body = (bundle_path.parent / artifact["path"]).read_bytes()
    artifact["size_bytes"] = len(body)
    artifact["sha256"] = hashlib.sha256(body).hexdigest()
    write_bundle(bundle_path, payload)


def tree_digest(root: Path) -> dict[str, str]:
    return {
        path.relative_to(root).as_posix(): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def _artifact(
    root: Path,
    path: str,
    role: str,
    media_type: str,
) -> dict[str, object]:
    body = (root / path).read_bytes()
    return {
        "role": role,
        "path": path,
        "media_type": media_type,
        "size_bytes": len(body),
        "sha256": hashlib.sha256(body).hexdigest(),
    }


def _validation(
    mode: str,
    artifact_fingerprint: str,
    strict_inputs_fingerprint: str | None,
) -> dict[str, object]:
    return {
        "mode": mode,
        "freshness": "passed",
        "recorded": True,
        "invoked": True,
        "outcome": "passed",
        "reason_code": None,
        "return_code": 0,
        "completed_at": "2026-07-15T00:00:00Z",
        "identity_state": "recorded",
        "identity_fingerprint": "sha256:" + "4" * 64,
        "artifact_fingerprint": artifact_fingerprint,
        "strict_inputs_fingerprint": strict_inputs_fingerprint,
    }
