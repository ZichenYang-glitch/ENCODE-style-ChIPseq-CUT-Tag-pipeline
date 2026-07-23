"""Deterministic contracts for the bulk RNA-seq acceptance harness."""

from __future__ import annotations

from dataclasses import replace
from hashlib import sha256
import json
import os
from pathlib import Path
import subprocess

import pytest

from encode_pipeline.adapters.bulk_rnaseq import BulkRnaSeqResultsWorkflowAdapter
from encode_pipeline.adapters.bulk_rnaseq.reference_closure import (
    REFERENCE_INDEX_MANIFEST,
    verify_reference_index,
)
from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    SORTMERNA_INDEX_MANIFEST_FILENAME,
    SORTMERNA_INDEX_MANIFEST_SCHEMA_VERSION,
    SORTMERNA_VERSION,
    verify_ribo_database_manifest,
    verify_sortmerna_index,
)
from encode_pipeline.services.defaults import create_default_workflow_registry
from scripts.generate_bulk_rnaseq_tiny_fixture import _reference_index_sidecar

from . import support as support_module
from .support import (
    ACCEPTANCE_EVIDENCE_SCHEMA_VERSION,
    FIXTURE_MANIFEST_SCHEMA_VERSION,
    AcceptanceEvidence,
    AcceptanceEvidenceStaleError,
    AcceptanceEvidenceValues,
    GateSettings,
    assert_no_managed_containers,
    build_results_composition,
    load_acceptance_fixture,
    read_tested_head,
    require_gate_settings,
)


_HEX_A = "a" * 64
_HEX_B = "b" * 64
_HEX_C = "c" * 64
_HEX_D = "d" * 64
_HEX_E = "e" * 64
_HEX_F = "f" * 64
_GIT_SHA1 = "a" * 40


def _gate_environment(tmp_path: Path) -> dict[str, str]:
    runtime_root = tmp_path / "runtime"
    runtime_root.mkdir()
    fixture_manifest = tmp_path / "fixture.json"
    fixture_manifest.write_text("{}\n", encoding="utf-8")
    docker_executable = tmp_path / "bin" / "docker"
    docker_executable.parent.mkdir()
    docker_executable.write_text("#!/bin/sh\nexit 0\n", encoding="utf-8")
    docker_executable.chmod(0o755)
    socket_path = tmp_path / "docker.sock"
    socket_path.write_text("test-only socket coordinate\n", encoding="utf-8")
    environment = {
        "HELIXWEAVE_REQUIRE_BULK_RNASEQ_REAL_EXECUTION": "1",
        "HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT": str(runtime_root),
        "HELIXWEAVE_BULK_RNASEQ_FIXTURE_MANIFEST": str(fixture_manifest),
        "ENCODE_PIPELINE_TEST_REDIS_URL": "redis://127.0.0.1:6379/15",
        "ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE": str(docker_executable),
        "ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET": str(socket_path),
    }
    return environment


def _fixture_document(fixture_root: Path) -> dict[str, object]:
    fixture_root.mkdir(parents=True, exist_ok=True)
    repeated_pe = {
        "sample": "PE1",
        "library": "lib1",
        "lane": "L001",
        "layout": "PE",
        "fastq_1": str(fixture_root / "reads/PE1_L001_R1.fastq.gz"),
        "fastq_2": str(fixture_root / "reads/PE1_L001_R2.fastq.gz"),
        "strandedness": "auto",
        "platform": "ILLUMINA",
    }
    artifact_coordinates = sorted(
        {
            *(
                ("PE1", f"bulk_rnaseq.fastqc.raw.{role}.zip")
                for role in ("read1", "read2")
            ),
            *(
                ("PE1", f"bulk_rnaseq.rrna.sortmerna.filtered.{role}")
                for role in ("read1", "read2")
            ),
            *(
                ("PE1", value)
                for value in (
                    "bulk_rnaseq.salmon.meta_info",
                    "bulk_rnaseq.salmon.quant_gene",
                    "bulk_rnaseq.star.bam",
                    "bulk_rnaseq.star.log_final",
                )
            ),
            ("SE1", "bulk_rnaseq.fastqc.raw.single.zip"),
            ("SE1", "bulk_rnaseq.rrna.sortmerna.filtered.single"),
            *(
                ("SE1", value)
                for value in (
                    "bulk_rnaseq.salmon.meta_info",
                    "bulk_rnaseq.salmon.quant_gene",
                    "bulk_rnaseq.star.bam",
                    "bulk_rnaseq.star.log_final",
                )
            ),
        }
    )
    metric_coordinates = sorted(
        {
            *(
                ("PE1", f"fastqc.raw.{role}.total_sequences")
                for role in ("read1", "read2")
            ),
            *(
                ("PE1", f"trimming.{role}.retained_reads")
                for role in ("read1", "read2")
            ),
            *(
                ("PE1", value)
                for value in (
                    "salmon.mapping_fraction",
                    "salmon.processed_fragments",
                    "star.input_templates",
                    "star.uniquely_mapped_template_fraction",
                )
            ),
            ("SE1", "fastqc.raw.single.total_sequences"),
            ("SE1", "trimming.single.retained_reads"),
            *(
                ("SE1", value)
                for value in (
                    "salmon.mapping_fraction",
                    "salmon.processed_fragments",
                    "star.input_templates",
                    "star.uniquely_mapped_template_fraction",
                )
            ),
        }
    )
    document: dict[str, object] = {
        "schema_version": FIXTURE_MANIFEST_SCHEMA_VERSION,
        "transcriptome_binding": {
            "reference_id": "tiny-reference",
            "fasta_sha256": _HEX_A,
            "gtf_sha256": _HEX_B,
            "transcript_fasta": str(fixture_root / "reference/transcripts.fa"),
            "transcript_fasta_sha256": "0" * 64,
        },
        "workflow_inputs": {
            "config": {
                "advanced": {
                    "min_mapped_reads": 0,
                    "min_trimmed_reads": 1,
                },
                "standard": {
                    "reference": {
                        "reference_id": "tiny-reference",
                        "fasta": str(fixture_root / "reference/genome.fa"),
                        "fasta_sha256": _HEX_A,
                        "gtf": str(fixture_root / "reference/genes.gtf"),
                        "gtf_sha256": _HEX_B,
                        "annotation_style": "ensembl",
                        "star_index": {
                            "path": str(fixture_root / "indexes/star"),
                            "identity_sha256": _HEX_C,
                        },
                        "salmon_index": {
                            "path": str(fixture_root / "indexes/salmon"),
                            "identity_sha256": _HEX_D,
                        },
                    },
                    "analysis": {
                        "alignment": "star",
                        "quantification": "salmon",
                    },
                    "trimming": {"enabled": True, "tool": "trimgalore"},
                    "ribosomal_rna_removal": {
                        "enabled": True,
                        "tool": "sortmerna",
                        "save_filtered_reads": True,
                        "database_manifest": {
                            "path": str(fixture_root / "rrna/database-manifest.txt"),
                            "identity_sha256": _HEX_E,
                        },
                        "sortmerna_index": {
                            "path": str(fixture_root / "indexes/sortmerna"),
                            "identity_sha256": _HEX_F,
                        },
                    },
                },
            },
            "samples": [
                repeated_pe,
                {**repeated_pe, "lane": "L002"},
                {
                    "sample": "SE1",
                    "library": "lib1",
                    "lane": "L001",
                    "layout": "SE",
                    "fastq_1": str(fixture_root / "reads/SE1_L001_R1.fastq.gz"),
                    "strandedness": "auto",
                    "platform": "ILLUMINA",
                },
            ],
            "options": {},
        },
        "required_artifact_output_types": sorted(
            {output_type for _sample, output_type in artifact_coordinates}
        ),
        "required_artifact_sample_output_types": [
            list(value) for value in artifact_coordinates
        ],
        "required_qc_metric_keys": sorted(
            {metric_key for _sample, metric_key in metric_coordinates}
        ),
        "required_qc_sample_metric_keys": [list(value) for value in metric_coordinates],
        "required_qc_sample_metric_values": [
            ["PE1", "fastqc.raw.read1.total_sequences", "768"],
            ["PE1", "fastqc.raw.read2.total_sequences", "768"],
            ["SE1", "fastqc.raw.single.total_sequences", "384"],
        ],
        "required_sample_ids": ["PE1", "SE1"],
    }
    _materialize_fixture_identity_closure(fixture_root, document)
    return document


def _canonical_json_bytes(value: object) -> bytes:
    return (
        json.dumps(
            value,
            allow_nan=False,
            ensure_ascii=True,
            separators=(",", ":"),
            sort_keys=True,
        )
        + "\n"
    ).encode("ascii")


def _framed_identity(scheme: str, payload: dict[str, object]) -> str:
    digest = sha256()
    for value in (scheme.encode("ascii"), _canonical_json_bytes(payload)):
        digest.update(len(value).to_bytes(8, "big"))
        digest.update(value)
    return digest.hexdigest()


def _materialize_fixture_identity_closure(
    fixture_root: Path,
    acceptance: dict[str, object],
) -> None:
    inputs = acceptance["workflow_inputs"]
    assert isinstance(inputs, dict)
    samples = inputs["samples"]
    assert isinstance(samples, list)
    source_samples: list[dict[str, object]] = []
    relative_files: set[str] = {
        "reference/genome.fa",
        "reference/genes.gtf",
        "reference/transcripts.fa",
        "rrna/database-manifest.txt",
        "rrna/tiny-rrna.fa",
    }
    for value in samples:
        assert isinstance(value, dict)
        row = dict(value)
        for key in ("fastq_1", "fastq_2"):
            if key in row:
                relative = Path(row[key]).relative_to(fixture_root).as_posix()
                row[key] = relative
                relative_files.add(relative)
        source_samples.append(row)
    contents = {
        "reference/genome.fa": b">chrTiny\nACGTACGT\n",
        "reference/genes.gtf": b'chrTiny\ttiny\tgene\t1\t8\t.\t+\t.\tgene_id "g1";\n',
        "reference/transcripts.fa": b">t1\nACGTACGT\n",
        "rrna/database-manifest.txt": b"tiny-rrna.fa\n",
        "rrna/tiny-rrna.fa": b">tiny-rRNA\nACGT\n",
    }
    for relative in sorted(relative_files):
        path = fixture_root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        content = contents.get(relative, b"controlled-fastq\n")
        path.write_bytes(content)
    file_entries = []
    for relative in sorted(relative_files):
        content = (fixture_root / relative).read_bytes()
        file_entries.append(
            {
                "path": relative,
                "size_bytes": len(content),
                "sha256": sha256(content).hexdigest(),
            }
        )
    source: dict[str, object] = {
        "schema_version": support_module.SOURCE_MANIFEST_SCHEMA_VERSION,
        "fixture_id": support_module.FIXTURE_ID,
        "biological_validity": False,
        "limitations": "controlled unit fixture; no biological validity",
        "generator": {"seed": 1},
        "reference": {"contig": "chrTiny"},
        "samples": source_samples,
        "files": file_entries,
        "index_build_contracts": {},
    }
    source["source_identity_sha256"] = _framed_identity(
        support_module.SOURCE_IDENTITY_SCHEME,
        source,
    )
    source_raw = _canonical_json_bytes(source)
    (fixture_root / support_module.SOURCE_MANIFEST_FILENAME).write_bytes(source_raw)

    config = inputs["config"]
    assert isinstance(config, dict)
    standard = config["standard"]
    assert isinstance(standard, dict)
    reference = standard["reference"]
    rrna = standard["ribosomal_rna_removal"]
    assert isinstance(reference, dict) and isinstance(rrna, dict)
    file_digests = {item["path"]: item["sha256"] for item in file_entries}
    reference["fasta_sha256"] = file_digests["reference/genome.fa"]
    reference["gtf_sha256"] = file_digests["reference/genes.gtf"]
    transcriptome = acceptance["transcriptome_binding"]
    assert isinstance(transcriptome, dict)
    transcriptome["fasta_sha256"] = reference["fasta_sha256"]
    transcriptome["gtf_sha256"] = reference["gtf_sha256"]
    transcriptome["transcript_fasta_sha256"] = file_digests["reference/transcripts.fa"]
    database = rrna["database_manifest"]
    assert isinstance(database, dict)
    database["identity_sha256"] = file_digests["rrna/database-manifest.txt"]

    index_closures: dict[str, str] = {}
    for kind in ("star", "salmon"):
        binding = reference[f"{kind}_index"]
        assert isinstance(binding, dict)
        root = Path(binding["path"])
        root.mkdir(parents=True, exist_ok=True)
        content = f"controlled-{kind}-index\n".encode()
        index_file = root / f"{kind}.bin"
        index_file.write_bytes(content)
        files = [
            {
                "path": index_file.name,
                "size_bytes": len(content),
                "sha256": sha256(content).hexdigest(),
            }
        ]
        sidecar = _reference_index_sidecar(
            kind=kind,
            container_image=support_module._EXPECTED_INDEX_PRODUCER_IMAGES[kind],
            fasta_sha256=reference["fasta_sha256"],
            gtf_sha256=reference["gtf_sha256"],
            transcript_fasta_sha256=transcriptome["transcript_fasta_sha256"],
            files=files,
            genome_sa_index_nbases=1,
        )
        sidecar_path = root / REFERENCE_INDEX_MANIFEST
        sidecar_path.write_bytes(sidecar)
        binding["identity_sha256"] = sha256(sidecar).hexdigest()
        result = verify_reference_index(
            root,
            kind=kind,
            expected_manifest_sha256=binding["identity_sha256"],
            fasta_sha256=reference["fasta_sha256"],
            gtf_sha256=reference["gtf_sha256"],
            transcript_fasta_sha256=(
                transcriptome["transcript_fasta_sha256"] if kind == "salmon" else None
            ),
            annotation_style=reference["annotation_style"],
            expected_container_image=support_module._EXPECTED_INDEX_PRODUCER_IMAGES[
                kind
            ],
        )
        assert result.is_success and result.value is not None
        index_closures[kind] = result.value.identity_sha256

    database_result = verify_ribo_database_manifest(
        database["path"],
        expected_manifest_sha256=database["identity_sha256"],
    )
    assert database_result.is_success and database_result.value is not None
    sortmerna_binding = rrna["sortmerna_index"]
    assert isinstance(sortmerna_binding, dict)
    sortmerna_root = Path(sortmerna_binding["path"])
    sortmerna_root.mkdir(parents=True, exist_ok=True)
    sortmerna_content = b"controlled-sortmerna-index\n"
    sortmerna_file = sortmerna_root / "tiny.idx"
    sortmerna_file.write_bytes(sortmerna_content)
    sortmerna_sidecar = _canonical_json_bytes(
        {
            "schema_version": SORTMERNA_INDEX_MANIFEST_SCHEMA_VERSION,
            "database_closure_sha256": database_result.value.identity_sha256,
            "sortmerna_version": SORTMERNA_VERSION,
            "files": [
                {
                    "path": sortmerna_file.name,
                    "size_bytes": len(sortmerna_content),
                    "sha256": sha256(sortmerna_content).hexdigest(),
                }
            ],
        }
    )
    (sortmerna_root / SORTMERNA_INDEX_MANIFEST_FILENAME).write_bytes(sortmerna_sidecar)
    sortmerna_binding["identity_sha256"] = sha256(sortmerna_sidecar).hexdigest()
    sortmerna_result = verify_sortmerna_index(
        sortmerna_root,
        expected_index_sha256=sortmerna_binding["identity_sha256"],
        database_closure=database_result.value,
    )
    assert sortmerna_result.is_success and sortmerna_result.value is not None
    index_closures["sortmerna"] = sortmerna_result.value.identity_sha256

    indexes: dict[str, object] = {}
    route_values = {
        "star": ("STAR_GENOMEGENERATE", "2.7.11b", reference["star_index"]),
        "salmon": ("SALMON_INDEX", "1.10.3", reference["salmon_index"]),
        "sortmerna": ("SORTMERNA", "4.3.7", rrna["sortmerna_index"]),
    }
    for kind, (process, version, binding) in route_values.items():
        assert isinstance(binding, dict)
        argv = [kind, "--controlled"]
        indexes[kind] = {
            "relative_index_root": Path(binding["path"])
            .relative_to(fixture_root)
            .as_posix(),
            "producer_process": process,
            "tool": kind,
            "tool_version": version,
            "container_image": support_module._EXPECTED_INDEX_PRODUCER_IMAGES[kind],
            "immutable_config_sha256": "1" * 64,
            "command_argv": argv,
            "command_argv_sha256": sha256(_canonical_json_bytes(argv)).hexdigest(),
            "sidecar_manifest_sha256": binding["identity_sha256"],
            "index_closure_sha256": index_closures[kind],
        }
    provenance: dict[str, object] = {
        "schema_version": support_module.INDEX_PROVENANCE_SCHEMA_VERSION,
        "fixture_id": support_module.FIXTURE_ID,
        "biological_validity": False,
        "source_manifest_sha256": sha256(source_raw).hexdigest(),
        "source_identity_sha256": source["source_identity_sha256"],
        "indexes": indexes,
    }
    provenance["provenance_identity_sha256"] = _framed_identity(
        support_module.PROVENANCE_IDENTITY_SCHEME,
        provenance,
    )
    provenance_raw = _canonical_json_bytes(provenance)
    (fixture_root / support_module.INDEX_PROVENANCE_FILENAME).write_bytes(
        provenance_raw
    )
    acceptance.update(
        {
            "source_manifest_sha256": sha256(source_raw).hexdigest(),
            "source_identity_sha256": source["source_identity_sha256"],
            "index_provenance_manifest_sha256": sha256(provenance_raw).hexdigest(),
            "index_provenance_identity_sha256": provenance[
                "provenance_identity_sha256"
            ],
        }
    )


def _evidence_values() -> AcceptanceEvidenceValues:
    return AcceptanceEvidenceValues(
        tested_head=_GIT_SHA1,
        workflow_build_digest=_HEX_B,
        fixture_acceptance_manifest_sha256="d" * 64,
        fixture_source_manifest_sha256="e" * 64,
        fixture_source_identity_sha256="f" * 64,
        fixture_index_provenance_manifest_sha256="0" * 64,
        fixture_index_provenance_identity_sha256="1" * 64,
        validated_snapshot_id=f"vsnap_{'0' * 32}",
        validated_payload_digest=_HEX_C,
        cache_identity_sha256=_HEX_D,
        input_closure_sha256=_HEX_E,
        ribo_database_closure_sha256=_HEX_F,
        workspace_identity_sha256="1" * 64,
        workspace_contract_sha256="2" * 64,
        execution_implementation_manifest_sha256="3" * 64,
        execution_implementation_aggregate_sha256="4" * 64,
        container_process_audit_sha256="5" * 64,
        run_id="run-acceptance-1",
        job_id="run-execution-acceptance-1",
        lifecycle_status="succeeded",
        artifact_attempt_id=f"resultattempt-{'6' * 64}",
        artifact_revision=1,
        artifact_generation=f"artifactgen-{'7' * 64}",
        artifact_manifest_digest="8" * 64,
        artifact_content_sha256=(("artifact-a", "9" * 64),),
        qc_attempt_id=f"resultattempt-{'a' * 64}",
        qc_revision=1,
        qc_generation=f"qcgen-{'b' * 64}",
        qc_manifest_digest="c" * 64,
        qc_artifact_generation=f"artifactgen-{'7' * 64}",
        artifact_output_types=("bulk_rnaseq.star.bam",),
        qc_metric_keys=("star.input_templates",),
        qc_sample_ids=("S1",),
        artifact_sample_output_types=(("S1", "bulk_rnaseq.star.bam"),),
        qc_sample_metric_keys=(("S1", "star.input_templates"),),
        qc_sample_metric_values=(("S1", "star.input_templates", "384"),),
        cleanup_confirmed=True,
    )


def test_gate_environment_is_explicit_and_missing_values_fail_without_skip(
    tmp_path: Path,
) -> None:
    with pytest.raises(AssertionError, match="HELIXWEAVE_REQUIRE"):
        require_gate_settings({})

    environment = _gate_environment(tmp_path)
    settings = require_gate_settings(environment, _socket_probe=lambda _path: True)

    assert settings.runtime_root == tmp_path / "runtime"
    assert settings.fixture_manifest == tmp_path / "fixture.json"
    assert settings.docker_executable.name == "docker"
    assert settings.docker_socket == tmp_path / "docker.sock"
    assert "redis" not in repr(settings).lower()


def test_gate_requires_every_private_runtime_coordinate(tmp_path: Path) -> None:
    environment = _gate_environment(tmp_path)
    for name in tuple(environment):
        incomplete = dict(environment)
        del incomplete[name]
        with pytest.raises(AssertionError, match=name):
            require_gate_settings(incomplete, _socket_probe=lambda _path: True)


def test_container_residual_audit_fails_before_any_repair(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    calls: list[tuple[str, ...]] = []

    class FakeCleaner:
        executable = Path("/usr/bin/docker")
        local_docker_host = "unix:///tmp/test-docker.sock"

        @staticmethod
        def verify_endpoint():
            return type("EndpointResult", (), {"is_failure": False})()

    def run(argv, **_kwargs):
        calls.append(tuple(argv))
        return subprocess.CompletedProcess(argv, 0, stdout=f"{'a' * 64}\n")

    monkeypatch.setattr(support_module.subprocess, "run", run)

    with pytest.raises(AssertionError, match="left a managed container residual"):
        assert_no_managed_containers(FakeCleaner(), "b" * 64)  # type: ignore[arg-type]

    assert len(calls) == 1
    assert "ps" in calls[0]
    assert not set(calls[0]) & {"stop", "kill", "rm"}


def test_tested_head_rejects_untracked_execution_sources(tmp_path: Path) -> None:
    repository = tmp_path / "repository"
    repository.mkdir()
    environment = {
        **os.environ,
        "GIT_AUTHOR_NAME": "HelixWeave test",
        "GIT_AUTHOR_EMAIL": "helixweave@example.invalid",
        "GIT_COMMITTER_NAME": "HelixWeave test",
        "GIT_COMMITTER_EMAIL": "helixweave@example.invalid",
    }

    def git(*arguments: str) -> subprocess.CompletedProcess[str]:
        return subprocess.run(
            ("git", *arguments),
            cwd=repository,
            env=environment,
            check=True,
            capture_output=True,
            text=True,
        )

    git("init", "--initial-branch=main")
    (repository / "tracked.py").write_text("VALUE = 1\n", encoding="utf-8")
    git("add", "tracked.py")
    git("commit", "-m", "test fixture")

    assert read_tested_head(repository) == git("rev-parse", "HEAD").stdout.strip()

    (repository / "untracked_execution.py").write_text(
        "VALUE = 2\n",
        encoding="utf-8",
    )
    with pytest.raises(AssertionError, match="clean exact Git HEAD"):
        read_tested_head(repository)


def test_fixture_manifest_requires_full_star_salmon_sortmerna_representation(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "fixture.json"
    manifest.write_text(
        json.dumps(_fixture_document(tmp_path), sort_keys=True),
        encoding="utf-8",
    )

    fixture = load_acceptance_fixture(manifest)

    assert {row["layout"] for row in fixture.workflow_inputs.samples} == {"SE", "PE"}
    assert fixture.required_sample_ids == ("PE1", "SE1")
    assert "bulk_rnaseq.star.bam" in fixture.required_artifact_output_types
    assert "salmon.mapping_fraction" in fixture.required_qc_metric_keys

    missing_lane = _fixture_document(tmp_path)
    missing_lane["workflow_inputs"]["samples"] = missing_lane["workflow_inputs"][
        "samples"
    ][1:]
    manifest.write_text(json.dumps(missing_lane), encoding="utf-8")
    with pytest.raises(AssertionError, match="repeated lane"):
        load_acceptance_fixture(manifest)


def test_fixture_manifest_rejects_non_string_evidence_tokens(tmp_path: Path) -> None:
    manifest = tmp_path / "fixture.json"
    document = _fixture_document(tmp_path)
    document["required_artifact_output_types"] = [{}]
    manifest.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(AssertionError, match="artifact output types"):
        load_acceptance_fixture(manifest)


def test_fixture_manifest_rejects_a_tampered_fastqc_oracle(tmp_path: Path) -> None:
    manifest = tmp_path / "fixture.json"
    document = _fixture_document(tmp_path)
    document["required_qc_sample_metric_values"][0][2] = "767"
    manifest.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(AssertionError, match="fixed source"):
        load_acceptance_fixture(manifest)


@pytest.mark.parametrize(
    "advanced",
    (
        {},
        {"min_mapped_reads": 0, "min_trimmed_reads": 10_000},
        {"min_mapped_reads": 0, "min_trimmed_reads": True},
    ),
)
def test_fixture_manifest_requires_exact_tiny_execution_thresholds(
    tmp_path: Path,
    advanced: dict[str, object],
) -> None:
    manifest = tmp_path / "fixture.json"
    document = _fixture_document(tmp_path)
    document["workflow_inputs"]["config"]["advanced"] = advanced
    manifest.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(AssertionError, match="tiny execution thresholds"):
        load_acceptance_fixture(manifest)


def test_fixture_manifest_rejects_resealed_absolute_provenance_argv(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "fixture.json"
    document = _fixture_document(tmp_path)
    provenance_path = tmp_path / support_module.INDEX_PROVENANCE_FILENAME
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    argv = ["star", "--genomeDir", "/private/runner/index"]
    provenance["indexes"]["star"]["command_argv"] = argv
    provenance["indexes"]["star"]["command_argv_sha256"] = sha256(
        _canonical_json_bytes(argv)
    ).hexdigest()
    provenance.pop("provenance_identity_sha256")
    provenance["provenance_identity_sha256"] = _framed_identity(
        support_module.PROVENANCE_IDENTITY_SCHEME,
        provenance,
    )
    provenance_raw = _canonical_json_bytes(provenance)
    provenance_path.write_bytes(provenance_raw)
    document["index_provenance_manifest_sha256"] = sha256(provenance_raw).hexdigest()
    document["index_provenance_identity_sha256"] = provenance[
        "provenance_identity_sha256"
    ]
    manifest.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(AssertionError, match="path-free"):
        load_acceptance_fixture(manifest)


def test_fixture_manifest_rejects_a_resealed_false_index_closure(
    tmp_path: Path,
) -> None:
    manifest = tmp_path / "fixture.json"
    document = _fixture_document(tmp_path)
    provenance_path = tmp_path / support_module.INDEX_PROVENANCE_FILENAME
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    provenance["indexes"]["star"]["index_closure_sha256"] = "9" * 64
    provenance.pop("provenance_identity_sha256")
    provenance["provenance_identity_sha256"] = _framed_identity(
        support_module.PROVENANCE_IDENTITY_SCHEME,
        provenance,
    )
    provenance_raw = _canonical_json_bytes(provenance)
    provenance_path.write_bytes(provenance_raw)
    document["index_provenance_manifest_sha256"] = sha256(provenance_raw).hexdigest()
    document["index_provenance_identity_sha256"] = provenance[
        "provenance_identity_sha256"
    ]
    manifest.write_text(json.dumps(document), encoding="utf-8")

    with pytest.raises(AssertionError, match="reference index closure changed"):
        load_acceptance_fixture(manifest)


def test_results_registry_is_explicit_and_default_remains_authoring_only_without_env(
    tmp_path: Path,
) -> None:
    settings = GateSettings(
        runtime_root=(tmp_path / "runtime").resolve(),
        fixture_manifest=(tmp_path / "fixture.json").resolve(),
        redis_url="redis://127.0.0.1:6379/15",
        docker_executable=Path("/usr/bin/docker"),
        docker_socket=Path("/var/run/docker.sock"),
    )
    settings.fixture_manifest.write_text(
        json.dumps(_fixture_document(tmp_path)),
        encoding="utf-8",
    )

    composition = build_results_composition(settings)

    default_bulk = create_default_workflow_registry(environ={}).get("bulk-rnaseq")
    assert not isinstance(default_bulk, BulkRnaSeqResultsWorkflowAdapter)
    assert default_bulk.capabilities.supports == ("validation", "input_authoring")
    adapter = composition.registry.get("bulk-rnaseq")
    assert isinstance(adapter, BulkRnaSeqResultsWorkflowAdapter)
    assert composition.build_identity_provider.registry is composition.registry
    assert composition.binding.assets.root == settings.runtime_root
    assert composition.binding.assets.docker_executable == settings.docker_executable
    assert composition.binding.assets.docker_socket == settings.docker_socket


def test_platform_harness_rejects_an_unbounded_worker_timeout(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    from .platform_harness import PlatformAcceptanceHarness

    settings = GateSettings(
        runtime_root=(tmp_path / "runtime").resolve(),
        fixture_manifest=(tmp_path / "fixture.json").resolve(),
        redis_url="redis://127.0.0.1:6379/15",
        docker_executable=Path("/usr/bin/docker"),
        docker_socket=Path("/var/run/docker.sock"),
    )
    monkeypatch.setenv("ENCODE_PIPELINE_JOB_TIMEOUT_SECONDS", "14401")

    with pytest.raises(ValueError, match="timeout exceeds"):
        PlatformAcceptanceHarness(
            gate_settings=settings,
            repository_root=tmp_path.resolve(),
            temporary_root=(tmp_path / "acceptance").resolve(),
        )


@pytest.mark.parametrize("tested_head", ("a" * 39, "a" * 41, "A" * 40))
def test_acceptance_evidence_rejects_an_invalid_git_commit_identity(
    tested_head: str,
) -> None:
    with pytest.raises(ValueError, match="tested_head"):
        replace(_evidence_values(), tested_head=tested_head)


@pytest.mark.parametrize(
    ("field", "replacement_value"),
    [
        ("tested_head", "d" * 40),
        ("workflow_build_digest", "d" * 64),
        ("fixture_acceptance_manifest_sha256", "e" * 64),
        ("fixture_source_identity_sha256", "e" * 64),
        ("validated_payload_digest", "d" * 64),
        ("cache_identity_sha256", "0" * 64),
        ("input_closure_sha256", "d" * 64),
        ("ribo_database_closure_sha256", "d" * 64),
        ("workspace_contract_sha256", "d" * 64),
        ("execution_implementation_manifest_sha256", "d" * 64),
        ("execution_implementation_aggregate_sha256", "d" * 64),
        ("container_process_audit_sha256", "d" * 64),
        ("artifact_generation", f"artifactgen-{'d' * 64}"),
        ("artifact_content_sha256", (("artifact-a", "d" * 64),)),
        ("qc_generation", f"qcgen-{'d' * 64}"),
    ],
)
def test_evidence_becomes_stale_after_source_input_runtime_or_result_change(
    field: str,
    replacement_value: object,
) -> None:
    original = AcceptanceEvidence.create(_evidence_values())
    changes = {field: replacement_value}
    if field == "artifact_generation":
        changes["qc_artifact_generation"] = replacement_value
    changed_values = replace(_evidence_values(), **changes)
    changed = AcceptanceEvidence.create(changed_values)

    with pytest.raises(AcceptanceEvidenceStaleError):
        original.assert_matches(changed)


def test_evidence_round_trip_is_canonical_and_rejects_tampering() -> None:
    evidence = AcceptanceEvidence.create(_evidence_values())
    payload = evidence.to_dict()

    assert payload["schema_version"] == ACCEPTANCE_EVIDENCE_SCHEMA_VERSION
    assert AcceptanceEvidence.from_dict(payload) == evidence
    assert "/" not in json.dumps(payload, sort_keys=True)

    payload["input_closure_sha256"] = "d" * 64
    with pytest.raises(AcceptanceEvidenceStaleError):
        AcceptanceEvidence.from_dict(payload)

    malformed = evidence.to_dict()
    malformed["artifact_content_sha256"] = [["artifact-a"]]
    with pytest.raises(ValueError, match="entries are invalid"):
        AcceptanceEvidence.from_dict(malformed)


def test_evidence_requires_distinct_artifact_and_qc_attempts() -> None:
    values = _evidence_values()

    with pytest.raises(ValueError, match="attempt identities must be distinct"):
        replace(values, qc_attempt_id=values.artifact_attempt_id)


def _workspace_identity_documents(
    workspace: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    cache_coordinates: dict[str, object] = {
        "schema_version": support_module.WORKSPACE_SCHEMA_VERSION,
        "adapter_version": "1.0.0",
        "execution_mode": "standard-v1",
        "workflow_build_sha256": _HEX_B,
        "execution_implementation_manifest_sha256": _HEX_C,
        "execution_implementation_aggregate_sha256": _HEX_D,
        "input_closure_sha256": _HEX_E,
        "ribo_database_closure_sha256": _HEX_F,
        "sortmerna_index_build_strategy": None,
        "normalized_inputs_sha256": "1" * 64,
        "nextflow_version": support_module.NEXTFLOW_VERSION,
        "resume_scope": "single-run-workspace",
    }
    cache = {
        **cache_coordinates,
        "identity_sha256": sha256(
            support_module._workspace_json_bytes(cache_coordinates)
        ).hexdigest(),
    }
    execution = {
        "schema_version": support_module.WORKSPACE_SCHEMA_VERSION,
        "workflow_id": "bulk-rnaseq",
        "adapter_version": cache["adapter_version"],
        "execution_mode": cache["execution_mode"],
        "build_identity_sha256": cache["workflow_build_sha256"],
        "input_identity_sha256": cache["input_closure_sha256"],
        "ribo_database_closure_sha256": cache["ribo_database_closure_sha256"],
        "sortmerna_index_build_strategy": None,
        "cache_identity_sha256": cache["identity_sha256"],
        "execution_implementation_manifest_sha256": cache[
            "execution_implementation_manifest_sha256"
        ],
        "execution_implementation_aggregate_sha256": cache[
            "execution_implementation_aggregate_sha256"
        ],
        "container_process_audit_sha256": "2" * 64,
        "workspace_identity_sha256": support_module.managed_container_scope(workspace),
        "workspace_contract_sha256": "3" * 64,
        "resume_enabled": False,
    }
    return cache, execution


def test_workspace_evidence_recomputes_cache_and_managed_scope(
    tmp_path: Path,
) -> None:
    workspace = (tmp_path / "workspace").resolve()
    cache, execution = _workspace_identity_documents(workspace)

    assert support_module._verified_workspace_identity(
        cache_identity=cache,
        execution_identity=execution,
        workspace=workspace,
        workflow_build_digest=_HEX_B,
    ) == support_module.managed_container_scope(workspace)

    changed_cache = dict(cache)
    changed_cache["input_closure_sha256"] = "4" * 64
    with pytest.raises(AssertionError, match="cache identity closure"):
        support_module._verified_workspace_identity(
            cache_identity=changed_cache,
            execution_identity=execution,
            workspace=workspace,
            workflow_build_digest=_HEX_B,
        )

    changed_execution = dict(execution)
    changed_execution["workspace_identity_sha256"] = "5" * 64
    with pytest.raises(AssertionError, match="execution identity closure"):
        support_module._verified_workspace_identity(
            cache_identity=cache,
            execution_identity=changed_execution,
            workspace=workspace,
            workflow_build_digest=_HEX_B,
        )
