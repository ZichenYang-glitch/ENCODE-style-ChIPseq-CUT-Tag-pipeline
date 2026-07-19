"""Contracts for the controlled bulk RNA-seq real-execution fixture tool."""

from __future__ import annotations

from dataclasses import replace
import gzip
import hashlib
import json
from pathlib import Path
import sys

import pytest


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from encode_pipeline.adapters.bulk_rnaseq.reference_closure import (
    REFERENCE_INDEX_MANIFEST,
    verify_reference_closure,
)
from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (
    SORTMERNA_INDEX_MANIFEST_FILENAME,
    verify_ribo_database_manifest,
    verify_sortmerna_index,
)
from scripts.generate_bulk_rnaseq_tiny_fixture import (
    ADAPTER_READS_PER_FASTQ,
    FIXTURE_MANIFEST_FILENAME,
    INDEX_PROVENANCE_FILENAME,
    READ_LENGTH,
    READS_PER_FASTQ,
    RRNA_READS_PER_FASTQ,
    FixtureAssembly,
    generate_fixture,
    main,
)
from bulk_rnaseq_real_execution.support import load_acceptance_fixture


STAR_IMAGE = (
    "community.wave.seqera.io/library/htslib_samtools_star_gawk:"
    "ae438e9a604351a4@sha256:"
    "4a468118dbd7491a69bf9813c68233afa8558d1f3380fd8cab03e0e3d3135190"
)
SALMON_IMAGE = (
    "quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2@sha256:"
    "f83ebb158845ee8138d793347f83b92c75e83c58dd8f4600c6fea2a2453ef08e"
)
SORTMERNA_IMAGE = (
    "community.wave.seqera.io/library/sortmerna:4.3.7--b730cad73fc42b8e"
    "@sha256:"
    "3c873f2a4c007c17b3b30aedab6b1d0d0670a9c629033686508c2d80d780a4af"
)


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _fastq_records(path: Path) -> list[tuple[str, str, str]]:
    lines = gzip.decompress(path.read_bytes()).decode("ascii").splitlines()
    assert len(lines) % 4 == 0
    records = []
    for offset in range(0, len(lines), 4):
        name, sequence, separator, quality = lines[offset : offset + 4]
        assert name.startswith("@")
        assert separator == "+"
        records.append((name, sequence, quality))
    return records


def _write_index_files(root: Path, kind: str) -> None:
    root.mkdir(parents=True)
    (root / "nested").mkdir()
    (root / f"{kind}.bin").write_bytes(f"real-{kind}-index\n".encode())
    (root / "nested" / "metadata.json").write_text(
        json.dumps({"kind": kind}, sort_keys=True) + "\n",
        encoding="utf-8",
    )


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _assembly(root: Path, *, output: Path | None = None) -> FixtureAssembly:
    return FixtureAssembly(
        fixture_root=root,
        star_index_root=(root / "indexes/star").resolve(),
        salmon_index_root=(root / "indexes/salmon").resolve(),
        sortmerna_index_root=(root / "indexes/sortmerna").resolve(),
        output_manifest=(output or root / FIXTURE_MANIFEST_FILENAME).resolve(),
        star_container_image=STAR_IMAGE,
        salmon_container_image=SALMON_IMAGE,
        sortmerna_container_image=SORTMERNA_IMAGE,
        star_tool_version="2.7.11b",
        salmon_tool_version="1.10.3",
        sortmerna_tool_version="4.3.7",
        star_config_sha256="d" * 64,
        salmon_config_sha256="e" * 64,
        sortmerna_config_sha256="f" * 64,
        star_command_argv=(
            "docker",
            "run",
            "--network=none",
            "{star_container_image}",
            "STAR",
            "--runMode",
            "genomeGenerate",
            "--genomeDir",
            "{star_index_root}",
        ),
        salmon_command_argv=(
            "docker",
            "run",
            "--network=none",
            "{salmon_container_image}",
            "salmon",
            "index",
            "-i",
            "{salmon_index_root}",
        ),
        sortmerna_command_argv=(
            "docker",
            "run",
            "--network=none",
            "{sortmerna_container_image}",
            "sortmerna",
            "--idx-dir",
            "{sortmerna_index_root}",
        ),
        genome_sa_index_nbases=7,
    )


def _prepare_index_trees(root: Path) -> None:
    for kind in ("star", "salmon", "sortmerna"):
        _write_index_files(root / f"indexes/{kind}", kind)


def test_generate_is_byte_deterministic_path_free_and_contains_real_read_routes(
    tmp_path: Path,
) -> None:
    first_root = (tmp_path / "first").resolve()
    second_root = (tmp_path / "second").resolve()

    first = generate_fixture(first_root)
    second = generate_fixture(second_root)

    assert first.source_identity_sha256 == second.source_identity_sha256
    assert first.manifest_sha256 == second.manifest_sha256
    assert first.manifest.read_bytes() == second.manifest.read_bytes()
    document = json.loads(first.manifest.read_bytes())
    assert document["schema_version"] == "1.0.0"
    assert document["biological_validity"] is False
    assert "execution and contract behavior" in document["limitations"]
    assert str(first_root) not in first.manifest.read_text(encoding="utf-8")
    assert document["generator"] == {
        "adapter_reads_per_fastq": ADAPTER_READS_PER_FASTQ,
        "read_length": READ_LENGTH,
        "reads_per_fastq": READS_PER_FASTQ,
        "rrna_reads_per_fastq": RRNA_READS_PER_FASTQ,
        "seed": 3260153,
    }
    assert [item["path"] for item in document["files"]] == sorted(
        item["path"] for item in document["files"]
    )
    for item in document["files"]:
        left = first_root / item["path"]
        right = second_root / item["path"]
        assert left.read_bytes() == right.read_bytes()
        assert item["size_bytes"] == left.stat().st_size
        assert item["sha256"] == _sha256(left)

    rrna = (first_root / "rrna/tiny-rrna.fa").read_text(encoding="ascii")
    genome = (first_root / "reference/genome.fa").read_text(encoding="ascii")
    adapter_prefix = "AGATCGGAAGAGC"
    fastqs = sorted((first_root / "reads").glob("*.fastq.gz"))
    assert len(fastqs) == 5
    for fastq in fastqs:
        records = _fastq_records(fastq)
        assert len(records) == READS_PER_FASTQ
        assert all(len(sequence) == READ_LENGTH for _, sequence, _ in records)
        assert all(quality == "I" * READ_LENGTH for _, _, quality in records)
        assert all(
            sequence in rrna or _reverse_complement(sequence) in rrna
            for _, sequence, _ in records[:RRNA_READS_PER_FASTQ]
        )
        assert all(
            adapter_prefix in sequence
            for _, sequence, _ in records[
                RRNA_READS_PER_FASTQ : RRNA_READS_PER_FASTQ + ADAPTER_READS_PER_FASTQ
            ]
        )
        assert any(
            sequence in genome or _reverse_complement(sequence) in genome
            for _, sequence, _ in records[
                RRNA_READS_PER_FASTQ + ADAPTER_READS_PER_FASTQ :
            ]
        )

    assert (first_root / "rrna/database-manifest.txt").read_text() == ("tiny-rrna.fa\n")
    gtf = (first_root / "reference/genes.gtf").read_text(encoding="ascii")
    transcript_headers = [
        line
        for line in (first_root / "reference/transcripts.fa")
        .read_text(encoding="ascii")
        .splitlines()
        if line.startswith(">")
    ]
    assert gtf.count("\texon\t") == len(transcript_headers) == 3


def test_generate_refuses_existing_or_non_absolute_destination(tmp_path: Path) -> None:
    existing = tmp_path / "existing"
    existing.mkdir()
    with pytest.raises(ValueError, match="must not already exist"):
        generate_fixture(existing.resolve())
    with pytest.raises(ValueError, match="absolute"):
        generate_fixture(Path("relative"))


def test_finalize_writes_and_verifies_production_sidecars_and_harness_manifest(
    tmp_path: Path,
) -> None:
    root = (tmp_path / "fixture").resolve()
    generated = generate_fixture(root)
    _prepare_index_trees(root)

    finalized = generated.finalize(_assembly(root))

    assert finalized.acceptance_manifest == root / FIXTURE_MANIFEST_FILENAME
    assert finalized.provenance_manifest == root / INDEX_PROVENANCE_FILENAME
    acceptance = json.loads(finalized.acceptance_manifest.read_bytes())
    source = json.loads(generated.manifest.read_bytes())
    provenance_bytes = finalized.provenance_manifest.read_bytes()
    provenance = json.loads(provenance_bytes)
    assert acceptance["schema_version"] == "1.1.0"
    assert (
        acceptance["source_manifest_sha256"]
        == generated.manifest_sha256
        == _sha256(generated.manifest)
    )
    assert (
        acceptance["source_identity_sha256"]
        == generated.source_identity_sha256
        == source["source_identity_sha256"]
    )
    assert (
        acceptance["index_provenance_manifest_sha256"]
        == finalized.provenance_manifest_sha256
        == _sha256(finalized.provenance_manifest)
    )
    assert (
        acceptance["index_provenance_identity_sha256"]
        == finalized.provenance_identity_sha256
        == provenance["provenance_identity_sha256"]
    )
    artifact_pairs = acceptance["required_artifact_sample_output_types"]
    qc_pairs = acceptance["required_qc_sample_metric_keys"]
    qc_values = acceptance["required_qc_sample_metric_values"]
    artifact_pair_tuples = [tuple(pair) for pair in artifact_pairs]
    qc_pair_tuples = [tuple(pair) for pair in qc_pairs]
    qc_value_tuples = [tuple(value) for value in qc_values]
    assert artifact_pair_tuples == sorted(set(artifact_pair_tuples))
    assert qc_pair_tuples == sorted(set(qc_pair_tuples))
    assert qc_value_tuples == sorted(set(qc_value_tuples))
    assert {pair for pair in artifact_pair_tuples if pair[0] == "PE1"} == {
        ("PE1", "bulk_rnaseq.fastqc.raw.read1.zip"),
        ("PE1", "bulk_rnaseq.fastqc.raw.read2.zip"),
        ("PE1", "bulk_rnaseq.rrna.sortmerna.filtered.read1"),
        ("PE1", "bulk_rnaseq.rrna.sortmerna.filtered.read2"),
        ("PE1", "bulk_rnaseq.salmon.meta_info"),
        ("PE1", "bulk_rnaseq.salmon.quant_gene"),
        ("PE1", "bulk_rnaseq.star.bam"),
        ("PE1", "bulk_rnaseq.star.log_final"),
    }
    assert {pair for pair in qc_pair_tuples if pair[0] == "PE1"} == {
        ("PE1", "fastqc.raw.read1.total_sequences"),
        ("PE1", "fastqc.raw.read2.total_sequences"),
        ("PE1", "salmon.mapping_fraction"),
        ("PE1", "salmon.processed_fragments"),
        ("PE1", "star.input_templates"),
        ("PE1", "star.uniquely_mapped_template_fraction"),
        ("PE1", "trimming.read1.retained_reads"),
        ("PE1", "trimming.read2.retained_reads"),
    }
    lane_counts = {
        sample_id: sum(
            1
            for row in acceptance["workflow_inputs"]["samples"]
            if row["sample"] == sample_id
        )
        for sample_id in ("PE1", "SE1")
    }
    assert set(qc_value_tuples) == {
        (
            "PE1",
            "fastqc.raw.read1.total_sequences",
            str(READS_PER_FASTQ * lane_counts["PE1"]),
        ),
        (
            "PE1",
            "fastqc.raw.read2.total_sequences",
            str(READS_PER_FASTQ * lane_counts["PE1"]),
        ),
        (
            "SE1",
            "fastqc.raw.single.total_sequences",
            str(READS_PER_FASTQ * lane_counts["SE1"]),
        ),
    }
    assert lane_counts == {"PE1": 2, "SE1": 1}
    assert {
        (sample_id, metric_key) for sample_id, metric_key, _value in qc_value_tuples
    }.issubset(qc_pair_tuples)
    assert acceptance["required_artifact_output_types"] == sorted(
        {output_type for _, output_type in artifact_pairs}
    )
    assert acceptance["required_qc_metric_keys"] == sorted(
        {metric_key for _, metric_key in qc_pairs}
    )
    fixture = load_acceptance_fixture(finalized.acceptance_manifest)
    assert {row["layout"] for row in fixture.workflow_inputs.samples} == {"SE", "PE"}
    assert fixture.required_sample_ids == ("PE1", "SE1")
    assert fixture.transcriptome.transcript_fasta == root / "reference/transcripts.fa"
    assert all(
        Path(row["fastq_1"]).is_absolute() for row in fixture.workflow_inputs.samples
    )

    workflow_inputs = json.loads(finalized.acceptance_manifest.read_bytes())[
        "workflow_inputs"
    ]
    reference = workflow_inputs["config"]["standard"]["reference"]
    reference_result = verify_reference_closure(
        reference,
        producer_images={
            "STAR_GENOMEGENERATE": STAR_IMAGE,
            "SALMON_INDEX": SALMON_IMAGE,
        },
        transcript_fasta_sha256=fixture.transcriptome.transcript_fasta_sha256,
    )
    assert reference_result.is_success
    database_config = workflow_inputs["config"]["standard"]["ribosomal_rna_removal"]
    database_result = verify_ribo_database_manifest(
        Path(database_config["database_manifest"]["path"]),
        expected_manifest_sha256=database_config["database_manifest"][
            "identity_sha256"
        ],
    )
    assert database_result.is_success
    sortmerna_result = verify_sortmerna_index(
        Path(database_config["sortmerna_index"]["path"]),
        expected_index_sha256=database_config["sortmerna_index"]["identity_sha256"],
        database_closure=database_result.value,
    )
    assert sortmerna_result.is_success

    assert str(root) not in provenance_bytes.decode("utf-8")
    assert provenance["source_manifest_sha256"] == generated.manifest_sha256
    assert provenance["indexes"]["star"]["tool_version"] == "2.7.11b"
    assert provenance["indexes"]["salmon"]["tool_version"] == "1.10.3"
    assert provenance["indexes"]["sortmerna"]["tool_version"] == "4.3.7"
    assert provenance["indexes"]["star"]["immutable_config_sha256"] == "d" * 64
    assert provenance["indexes"]["star"]["command_argv"][0] == "docker"
    assert len(provenance["provenance_identity_sha256"]) == 64

    for root_path, sidecar in (
        (root / "indexes/star", REFERENCE_INDEX_MANIFEST),
        (root / "indexes/salmon", REFERENCE_INDEX_MANIFEST),
        (root / "indexes/sortmerna", SORTMERNA_INDEX_MANIFEST_FILENAME),
    ):
        assert (root_path / sidecar).is_file()


def test_finalize_fails_closed_on_source_tamper_and_absolute_command(
    tmp_path: Path,
) -> None:
    tampered_root = (tmp_path / "tampered").resolve()
    generated = generate_fixture(tampered_root)
    _prepare_index_trees(tampered_root)
    (tampered_root / "reference/genome.fa").write_text(
        ">chrTiny\nchanged\n", encoding="ascii"
    )
    with pytest.raises(ValueError, match="source fixture closure"):
        generated.finalize(_assembly(tampered_root))
    assert not (tampered_root / FIXTURE_MANIFEST_FILENAME).exists()

    unsafe_root = (tmp_path / "unsafe").resolve()
    unsafe = generate_fixture(unsafe_root)
    _prepare_index_trees(unsafe_root)
    bad = _assembly(unsafe_root)
    bad = FixtureAssembly(
        **{
            **vars(bad),
            "star_command_argv": ("STAR", "--genomeDir", str(bad.star_index_root)),
        }
    )
    with pytest.raises(ValueError, match="path-free"):
        unsafe.finalize(bad)
    assert not (unsafe_root / "indexes/star" / REFERENCE_INDEX_MANIFEST).exists()


def test_finalize_detects_index_change_instead_of_rebinding_existing_sidecar(
    tmp_path: Path,
) -> None:
    root = (tmp_path / "fixture").resolve()
    generated = generate_fixture(root)
    _prepare_index_trees(root)
    assembly = _assembly(root)
    generated.finalize(assembly)
    (root / "indexes/star/star.bin").write_bytes(b"changed-real-index\n")

    with pytest.raises(ValueError, match="existing sidecar"):
        generated.finalize(assembly)


def test_finalize_refuses_tampered_bound_provenance_manifest(tmp_path: Path) -> None:
    root = (tmp_path / "fixture").resolve()
    generated = generate_fixture(root)
    _prepare_index_trees(root)
    assembly = _assembly(root)
    finalized = generated.finalize(assembly)
    acceptance_before = finalized.acceptance_manifest.read_bytes()
    finalized.provenance_manifest.write_bytes(
        finalized.provenance_manifest.read_bytes() + b" "
    )

    with pytest.raises(ValueError, match="existing output manifest"):
        generated.finalize(assembly)

    assert finalized.acceptance_manifest.read_bytes() == acceptance_before


def test_index_rebuild_changes_bound_provenance_but_not_source_identity(
    tmp_path: Path,
) -> None:
    first_root = (tmp_path / "first").resolve()
    second_root = (tmp_path / "second").resolve()
    first_generated = generate_fixture(first_root)
    second_generated = generate_fixture(second_root)
    _prepare_index_trees(first_root)
    _prepare_index_trees(second_root)
    (second_root / "indexes/star/star.bin").write_bytes(b"rebuilt-real-star-index\n")

    first_finalized = first_generated.finalize(_assembly(first_root))
    second_finalized = second_generated.finalize(_assembly(second_root))
    first_acceptance = json.loads(first_finalized.acceptance_manifest.read_bytes())
    second_acceptance = json.loads(second_finalized.acceptance_manifest.read_bytes())

    assert (
        first_acceptance["source_manifest_sha256"]
        == second_acceptance["source_manifest_sha256"]
    )
    assert (
        first_acceptance["source_identity_sha256"]
        == second_acceptance["source_identity_sha256"]
    )
    assert (
        first_acceptance["index_provenance_manifest_sha256"]
        != second_acceptance["index_provenance_manifest_sha256"]
    )
    assert (
        first_acceptance["index_provenance_identity_sha256"]
        != second_acceptance["index_provenance_identity_sha256"]
    )
    for acceptance, generated, finalized in (
        (first_acceptance, first_generated, first_finalized),
        (second_acceptance, second_generated, second_finalized),
    ):
        provenance = json.loads(finalized.provenance_manifest.read_bytes())
        assert acceptance["source_manifest_sha256"] == _sha256(generated.manifest)
        assert acceptance["source_identity_sha256"] == generated.source_identity_sha256
        assert acceptance["index_provenance_manifest_sha256"] == _sha256(
            finalized.provenance_manifest
        )
        assert (
            acceptance["index_provenance_identity_sha256"]
            == provenance["provenance_identity_sha256"]
        )


def test_finalize_rejects_star_index_value_that_is_not_the_pinned_nfcore_auto(
    tmp_path: Path,
) -> None:
    root = (tmp_path / "fixture").resolve()
    generated = generate_fixture(root)
    _prepare_index_trees(root)
    assembly = replace(_assembly(root), genome_sa_index_nbases=5)

    with pytest.raises(ValueError, match="nf-core auto"):
        generated.finalize(assembly)


def test_cli_requires_contract_tool_versions(tmp_path: Path) -> None:
    root = (tmp_path / "cli-fixture").resolve()
    assert main(["generate", "--output-root", str(root)]) == 0
    _prepare_index_trees(root)

    with pytest.raises(SystemExit):
        main(
            [
                "finalize",
                "--fixture-root",
                str(root),
                "--star-index-root",
                str((root / "indexes/star").resolve()),
                "--salmon-index-root",
                str((root / "indexes/salmon").resolve()),
                "--sortmerna-index-root",
                str((root / "indexes/sortmerna").resolve()),
                "--output-manifest",
                str((root / FIXTURE_MANIFEST_FILENAME).resolve()),
                "--star-container-image",
                STAR_IMAGE,
                "--salmon-container-image",
                SALMON_IMAGE,
                "--sortmerna-container-image",
                SORTMERNA_IMAGE,
                "--star-tool-version",
                "2.7.10b",
                "--salmon-tool-version",
                "1.10.3",
                "--sortmerna-tool-version",
                "4.3.7",
                "--star-config-sha256",
                "d" * 64,
                "--salmon-config-sha256",
                "e" * 64,
                "--sortmerna-config-sha256",
                "f" * 64,
                "--star-command-json",
                '["STAR"]',
                "--salmon-command-json",
                '["salmon","index"]',
                "--sortmerna-command-json",
                '["sortmerna"]',
                "--genome-sa-index-nbases",
                "6",
            ]
        )
