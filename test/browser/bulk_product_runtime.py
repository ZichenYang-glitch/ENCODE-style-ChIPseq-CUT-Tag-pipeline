"""Prepare verified product inputs for the protected Bulk RNA-seq browser Gate."""

from __future__ import annotations

from collections.abc import Callable, Mapping
import csv
from dataclasses import dataclass
import json
import os
from pathlib import Path
import sys

from encode_pipeline.adapters.bulk_rnaseq.deployment import (
    TRANSCRIPTOME_BINDING_MANIFEST_ENV,
    TRANSCRIPTOME_BINDING_SCHEMA_VERSION,
)
from encode_pipeline.services.defaults import create_default_workflow_registry
from encode_pipeline.services.workflow_info import WorkflowInfoService


TEST_ROOT = Path(__file__).resolve().parents[1]
if str(TEST_ROOT) not in sys.path:
    sys.path.insert(0, str(TEST_ROOT))

from bulk_rnaseq_real_execution.support import (  # noqa: E402
    AcceptanceFixture,
    GateSettings,
    load_acceptance_fixture,
    require_gate_settings,
)


_SAMPLE_FIELDS = (
    "sample",
    "library",
    "lane",
    "layout",
    "fastq_1",
    "fastq_2",
    "strandedness",
    "platform",
)


@dataclass(frozen=True)
class BulkProductBrowserRuntime:
    """Browser-visible fixture fields and one private deployment coordinate."""

    manifest_fields: Mapping[str, object]
    deployment_environment: Mapping[str, str]


def prepare_bulk_product_browser_runtime(
    runtime_root: Path,
    environ: Mapping[str, str] | None = None,
    *,
    settings_loader: Callable[[Mapping[str, str] | None], GateSettings] = (
        require_gate_settings
    ),
    fixture_loader: Callable[[Path], AcceptanceFixture] = load_acceptance_fixture,
    admission_probe: Callable[[Mapping[str, str]], None] | None = None,
) -> BulkProductBrowserRuntime:
    """Derive public UI inputs from the already-verified synthetic fixture."""
    if (
        not isinstance(runtime_root, Path)
        or not runtime_root.is_absolute()
        or runtime_root.resolve(strict=False) != runtime_root
        or not runtime_root.is_dir()
        or runtime_root.is_symlink()
    ):
        raise ValueError("browser runtime root must be a canonical directory")
    source = dict(os.environ if environ is None else environ)
    settings = settings_loader(source)
    if not isinstance(settings, GateSettings):
        raise ValueError("gate settings loader returned an invalid value")
    fixture = fixture_loader(settings.fixture_manifest)
    if not isinstance(fixture, AcceptanceFixture):
        raise ValueError("fixture loader returned an invalid value")

    product_root = runtime_root / "bulk-product"
    product_root.mkdir(mode=0o700)
    samples_path = product_root / "samples.tsv"
    _write_samples(samples_path, fixture)
    binding_manifest = product_root / "transcriptome-binding.json"
    _write_binding_manifest(binding_manifest, fixture)

    inputs = fixture.workflow_inputs.to_dict()
    manifest_fields = {
        "bulkWorkflowId": "bulk-rnaseq",
        "bulkSamplesPath": str(samples_path),
        "bulkConfig": inputs["config"],
        "bulkOptions": inputs["options"],
        "bulkExpectedExecution": "available",
        "bulkRequiredArtifactOutputTypes": list(fixture.required_artifact_output_types),
        "bulkRequiredQcMetricKeys": list(fixture.required_qc_metric_keys),
    }
    deployment_environment = {TRANSCRIPTOME_BINDING_MANIFEST_ENV: str(binding_manifest)}
    combined_environment = {**source, **deployment_environment}
    (admission_probe or _require_product_available)(combined_environment)
    return BulkProductBrowserRuntime(
        manifest_fields=manifest_fields,
        deployment_environment=deployment_environment,
    )


def _require_product_available(environ: Mapping[str, str]) -> None:
    registry = create_default_workflow_registry(environ=environ)
    descriptor = WorkflowInfoService(registry).get_descriptor("bulk-rnaseq")
    if descriptor.is_failure or descriptor.value is None:
        raise AssertionError("default Bulk RNA-seq product descriptor is unavailable")
    value = descriptor.value
    upstream = value.upstream_identity
    required_capabilities = {
        "workspace_plan",
        "command",
        "artifact_extract",
        "qc_summary_extract",
    }
    if (
        value.availability.execution != "available"
        or value.availability.reason_code != "WORKFLOW_EXECUTION_READY"
        or upstream is None
        or upstream.name != "nf-core/rnaseq"
        or upstream.version != "3.26.0"
        or upstream.revision != "e7ca46272c8f9d5ceee3f71759f4ba551d3217a4"
        or not required_capabilities <= set(value.capabilities.supports)
    ):
        raise AssertionError("default Bulk RNA-seq product admission failed")


def _write_samples(path: Path, fixture: AcceptanceFixture) -> None:
    samples = fixture.workflow_inputs.samples
    if not isinstance(samples, list) or not samples:
        raise ValueError("acceptance samples are unavailable")
    allowed = set(_SAMPLE_FIELDS)
    for row in samples:
        if not isinstance(row, dict) or not set(row) <= allowed:
            raise ValueError("acceptance sample fields are invalid")
        if set(row) - {"fastq_2"} != allowed - {"fastq_2"}:
            raise ValueError("acceptance sample fields are incomplete")
        if any(
            not isinstance(value, str)
            or any(character in value for character in ("\x00", "\t", "\n", "\r"))
            for value in row.values()
        ):
            raise ValueError("acceptance sample values are invalid")
    with path.open("x", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=_SAMPLE_FIELDS,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="raise",
        )
        writer.writeheader()
        writer.writerows(samples)


def _write_binding_manifest(path: Path, fixture: AcceptanceFixture) -> None:
    transcriptome = fixture.transcriptome
    payload = {
        "schema_version": TRANSCRIPTOME_BINDING_SCHEMA_VERSION,
        "reference_id": transcriptome.reference_id,
        "fasta_sha256": transcriptome.fasta_sha256,
        "gtf_sha256": transcriptome.gtf_sha256,
        "transcript_fasta": str(transcriptome.transcript_fasta),
        "transcript_fasta_sha256": transcriptome.transcript_fasta_sha256,
    }
    with path.open("x", encoding="utf-8") as handle:
        json.dump(payload, handle, ensure_ascii=True, sort_keys=True)
        handle.write("\n")
