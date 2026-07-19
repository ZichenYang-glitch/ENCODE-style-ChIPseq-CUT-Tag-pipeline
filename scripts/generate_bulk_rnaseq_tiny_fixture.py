#!/usr/bin/env python3
"""Generate and close a controlled tiny bulk RNA-seq execution fixture.

The ``generate`` phase writes only deterministic synthetic input bytes.  It
does not invoke a container or access the network.  Operators then build STAR,
Salmon, and SortMeRNA indexes with the admitted pinned runtime.  The
``finalize`` phase records those already-built trees using the production
sidecar schemas, verifies the same production closures, and emits the exact
private acceptance-harness manifest.

This tiny fixture demonstrates execution and contract behavior only.  It is
not evidence of biological validity or production-scale scientific quality.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import hashlib
import json
import math
from pathlib import Path, PurePosixPath
import random
import re
import shutil
import stat
import struct
import sys
from typing import Any, Mapping, Sequence
import zlib


PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from encode_pipeline.adapters.bulk_rnaseq import (  # noqa: E402
    BulkRnaSeqWorkflowAdapter,
)
from encode_pipeline.adapters.bulk_rnaseq.reference_closure import (  # noqa: E402
    REFERENCE_INDEX_MANIFEST,
    REFERENCE_INDEX_SCHEMA_VERSION,
    verify_reference_closure,
)
from encode_pipeline.adapters.bulk_rnaseq.resource_closure import (  # noqa: E402
    SORTMERNA_INDEX_MANIFEST_FILENAME,
    SORTMERNA_INDEX_MANIFEST_SCHEMA_VERSION,
    SORTMERNA_VERSION,
    verify_ribo_database_manifest,
    verify_sortmerna_index,
)
from encode_pipeline.platform.adapters import WorkflowInputs  # noqa: E402


SOURCE_MANIFEST_FILENAME = "helixweave-bulk-rnaseq-tiny-source.json"
INDEX_PROVENANCE_FILENAME = "helixweave-bulk-rnaseq-index-provenance.json"
FIXTURE_MANIFEST_FILENAME = "helixweave-bulk-rnaseq-acceptance-fixture.json"
SOURCE_MANIFEST_SCHEMA_VERSION = "1.0.0"
INDEX_PROVENANCE_SCHEMA_VERSION = "1.0.0"
FIXTURE_MANIFEST_SCHEMA_VERSION = "1.1.0"
FIXTURE_ID = "helixweave-bulk-rnaseq-tiny-v1"
SOURCE_IDENTITY_SCHEME = "sha256-framed-bulk-rnaseq-tiny-source-v1"
PROVENANCE_IDENTITY_SCHEME = "sha256-framed-bulk-rnaseq-index-provenance-v1"
SOURCE_LIMITATIONS = (
    "This controlled synthetic tiny fixture proves execution and contract "
    "behavior only; it does not establish biological validity, reference "
    "quality, or production-scale performance."
)

SEED = 3_260_153
READ_LENGTH = 101
READS_PER_FASTQ = 384
RRNA_READS_PER_FASTQ = 24
ADAPTER_READS_PER_FASTQ = 24
GENOME_LENGTH = 70_000
RRNA_LENGTH = 5_000
ADAPTER_SEQUENCE = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"

STAR_TOOL_VERSION = "2.7.11b"
SALMON_TOOL_VERSION = "1.10.3"
SALMON_KMER_SIZE = 31
STAR_SJDB_OVERHANG = READ_LENGTH - 1

_HEX64 = re.compile(r"^[0-9a-f]{64}$")
_IMMUTABLE_IMAGE = re.compile(r"^[^@\s]+@sha256:[0-9a-f]{64}$")
_WINDOWS_ABSOLUTE = re.compile(r"^[A-Za-z]:[\\/]")
_EMBEDDED_ABSOLUTE = re.compile(r"(?:^|[=,\s])/(?!/)")
_SOURCE_MANIFEST_KEYS = {
    "schema_version",
    "fixture_id",
    "biological_validity",
    "limitations",
    "generator",
    "reference",
    "samples",
    "files",
    "index_build_contracts",
    "source_identity_sha256",
}
_SOURCE_RELATIVE_FILES = (
    "reads/PE1_L001_R1.fastq.gz",
    "reads/PE1_L001_R2.fastq.gz",
    "reads/PE1_L002_R1.fastq.gz",
    "reads/PE1_L002_R2.fastq.gz",
    "reads/SE1_L001_R1.fastq.gz",
    "reference/genes.gtf",
    "reference/genome.fa",
    "reference/transcripts.fa",
    "rrna/database-manifest.txt",
    "rrna/tiny-rrna.fa",
)


@dataclass(frozen=True)
class Gene:
    gene_id: str
    transcript_id: str
    start: int
    length: int
    strand: str = "+"

    @property
    def end(self) -> int:
        return self.start + self.length - 1


GENES = (
    Gene("GENE_SE", "TX_SE", 1_001, 5_000),
    Gene("GENE_PE_A", "TX_PE_A", 12_001, 5_000),
    Gene("GENE_PE_B", "TX_PE_B", 24_001, 5_000),
)


@dataclass(frozen=True)
class GeneratedFixture:
    """Identity and location of one generated source-only fixture."""

    root: Path
    manifest: Path
    manifest_sha256: str
    source_identity_sha256: str

    def finalize(self, assembly: "FixtureAssembly") -> "FinalizedFixture":
        """Close real index trees and assemble the private acceptance input."""
        if assembly.fixture_root != self.root:
            raise ValueError("assembly fixture_root does not match generated fixture")
        return finalize_fixture(assembly)


@dataclass(frozen=True)
class FixtureAssembly:
    """Explicit identities and path-free receipts for already-built indexes."""

    fixture_root: Path
    star_index_root: Path
    salmon_index_root: Path
    sortmerna_index_root: Path
    output_manifest: Path
    star_container_image: str
    salmon_container_image: str
    sortmerna_container_image: str
    star_tool_version: str
    salmon_tool_version: str
    sortmerna_tool_version: str
    star_config_sha256: str
    salmon_config_sha256: str
    sortmerna_config_sha256: str
    star_command_argv: tuple[str, ...]
    salmon_command_argv: tuple[str, ...]
    sortmerna_command_argv: tuple[str, ...]
    genome_sa_index_nbases: int


@dataclass(frozen=True)
class FinalizedFixture:
    """Final acceptance and path-free index-provenance manifests."""

    acceptance_manifest: Path
    acceptance_manifest_sha256: str
    provenance_manifest: Path
    provenance_manifest_sha256: str
    provenance_identity_sha256: str


def generate_fixture(output_root: Path) -> GeneratedFixture:
    """Create deterministic FASTA/GTF/FASTQ inputs under a new absolute root."""
    root = _new_output_root(output_root)
    created = False
    try:
        root.mkdir()
        created = True
        for relative in ("reads", "reference", "rrna"):
            (root / relative).mkdir()

        sequence_rng = random.Random(SEED)
        genome = _random_dna(sequence_rng, GENOME_LENGTH)
        rrna = _random_dna(sequence_rng, RRNA_LENGTH)
        transcripts = {
            gene.transcript_id: genome[gene.start - 1 : gene.end] for gene in GENES
        }

        _write_ascii(root / "reference/genome.fa", _fasta("chrTiny", genome))
        _write_ascii(root / "reference/genes.gtf", _gtf())
        _write_ascii(
            root / "reference/transcripts.fa",
            "".join(
                _fasta(transcript_id, sequence)
                for transcript_id, sequence in transcripts.items()
            ),
        )
        _write_ascii(root / "rrna/tiny-rrna.fa", _fasta("tiny_rRNA", rrna))
        _write_ascii(root / "rrna/database-manifest.txt", "tiny-rrna.fa\n")

        samples = _sample_rows()
        for row in samples:
            read1, read2 = _reads_for_row(row, transcripts=transcripts, rrna=rrna)
            _write_deterministic_gzip_fastq(root / row["fastq_1"], read1)
            if read2 is not None:
                _write_deterministic_gzip_fastq(root / row["fastq_2"], read2)

        files = [_file_entry(root, relative) for relative in _SOURCE_RELATIVE_FILES]
        payload: dict[str, object] = {
            "schema_version": SOURCE_MANIFEST_SCHEMA_VERSION,
            "fixture_id": FIXTURE_ID,
            "biological_validity": False,
            "limitations": SOURCE_LIMITATIONS,
            "generator": {
                "seed": SEED,
                "read_length": READ_LENGTH,
                "reads_per_fastq": READS_PER_FASTQ,
                "rrna_reads_per_fastq": RRNA_READS_PER_FASTQ,
                "adapter_reads_per_fastq": ADAPTER_READS_PER_FASTQ,
            },
            "reference": _source_reference(),
            "samples": samples,
            "files": files,
            "index_build_contracts": _index_build_contracts(),
        }
        source_identity = _identity(SOURCE_IDENTITY_SCHEME, payload)
        payload["source_identity_sha256"] = source_identity
        manifest = root / SOURCE_MANIFEST_FILENAME
        _write_new(manifest, _canonical_json(payload))
        return GeneratedFixture(
            root=root,
            manifest=manifest,
            manifest_sha256=_sha256_file(manifest),
            source_identity_sha256=source_identity,
        )
    except BaseException:
        if created:
            shutil.rmtree(root)
        raise


def finalize_fixture(assembly: FixtureAssembly) -> FinalizedFixture:
    """Bind and verify real index trees, then write exact acceptance inputs."""
    _validate_assembly(assembly)
    source = _verify_source_fixture(assembly.fixture_root)
    expected_genome_sa_index_nbases = _nfcore_star_genome_sa_index_nbases(
        source["reference"]["genome_length"]
    )
    if assembly.genome_sa_index_nbases != expected_genome_sa_index_nbases:
        raise ValueError(
            "genome_sa_index_nbases does not match the pinned nf-core auto strategy"
        )
    source_manifest = assembly.fixture_root / SOURCE_MANIFEST_FILENAME
    source_manifest_sha256 = _sha256_file(source_manifest)

    star_files = _index_files(
        assembly.star_index_root,
        excluded=REFERENCE_INDEX_MANIFEST,
    )
    salmon_files = _index_files(
        assembly.salmon_index_root,
        excluded=REFERENCE_INDEX_MANIFEST,
    )
    sortmerna_files = _index_files(
        assembly.sortmerna_index_root,
        excluded=SORTMERNA_INDEX_MANIFEST_FILENAME,
    )
    fasta_sha256 = _source_file_digest(source, "reference/genome.fa")
    gtf_sha256 = _source_file_digest(source, "reference/genes.gtf")
    transcript_fasta_sha256 = _source_file_digest(
        source,
        "reference/transcripts.fa",
    )

    star_sidecar = _reference_index_sidecar(
        kind="star",
        container_image=assembly.star_container_image,
        fasta_sha256=fasta_sha256,
        gtf_sha256=gtf_sha256,
        transcript_fasta_sha256=transcript_fasta_sha256,
        files=star_files,
        genome_sa_index_nbases=assembly.genome_sa_index_nbases,
    )
    salmon_sidecar = _reference_index_sidecar(
        kind="salmon",
        container_image=assembly.salmon_container_image,
        fasta_sha256=fasta_sha256,
        gtf_sha256=gtf_sha256,
        transcript_fasta_sha256=transcript_fasta_sha256,
        files=salmon_files,
        genome_sa_index_nbases=assembly.genome_sa_index_nbases,
    )

    database_manifest = assembly.fixture_root / "rrna/database-manifest.txt"
    database_result = verify_ribo_database_manifest(
        database_manifest,
        expected_manifest_sha256=_sha256_file(database_manifest),
    )
    if database_result.is_failure or database_result.value is None:
        raise ValueError("rRNA database closure does not satisfy production contract")
    sortmerna_sidecar = _canonical_json(
        {
            "schema_version": SORTMERNA_INDEX_MANIFEST_SCHEMA_VERSION,
            "database_closure_sha256": database_result.value.identity_sha256,
            "sortmerna_version": SORTMERNA_VERSION,
            "files": sortmerna_files,
        }
    )

    sidecars = (
        (assembly.star_index_root / REFERENCE_INDEX_MANIFEST, star_sidecar),
        (assembly.salmon_index_root / REFERENCE_INDEX_MANIFEST, salmon_sidecar),
        (
            assembly.sortmerna_index_root / SORTMERNA_INDEX_MANIFEST_FILENAME,
            sortmerna_sidecar,
        ),
    )
    _require_sidecars_absent_or_identical(sidecars)
    created_sidecars: list[Path] = []
    created_outputs: list[Path] = []
    try:
        for path, content in sidecars:
            if not path.exists():
                _write_new(path, content)
                created_sidecars.append(path)

        reference_document = {
            "reference_id": FIXTURE_ID,
            "fasta": str((assembly.fixture_root / "reference/genome.fa").resolve()),
            "fasta_sha256": fasta_sha256,
            "gtf": str((assembly.fixture_root / "reference/genes.gtf").resolve()),
            "gtf_sha256": gtf_sha256,
            "annotation_style": "ensembl",
            "star_index": {
                "path": str(assembly.star_index_root),
                "identity_sha256": _sha256_file(sidecars[0][0]),
            },
            "salmon_index": {
                "path": str(assembly.salmon_index_root),
                "identity_sha256": _sha256_file(sidecars[1][0]),
            },
        }
        producer_images = {
            "STAR_GENOMEGENERATE": assembly.star_container_image,
            "SALMON_INDEX": assembly.salmon_container_image,
        }
        reference_result = verify_reference_closure(
            reference_document,
            producer_images=producer_images,
            transcript_fasta_sha256=transcript_fasta_sha256,
        )
        if reference_result.is_failure or reference_result.value is None:
            raise ValueError(
                "reference index closures do not satisfy production contract"
            )
        sortmerna_result = verify_sortmerna_index(
            assembly.sortmerna_index_root,
            expected_index_sha256=_sha256_file(sidecars[2][0]),
            database_closure=database_result.value,
        )
        if sortmerna_result.is_failure or sortmerna_result.value is None:
            raise ValueError(
                "SortMeRNA index closure does not satisfy production contract"
            )

        provenance = _provenance_document(
            assembly=assembly,
            source_manifest_sha256=source_manifest_sha256,
            source_identity_sha256=source["source_identity_sha256"],
            reference_closure=reference_result.value,
            sortmerna_closure=sortmerna_result.value,
        )
        provenance_bytes = _canonical_json(provenance)
        acceptance = _acceptance_document(
            assembly=assembly,
            source=source,
            reference=reference_document,
            database_manifest_sha256=database_result.value.manifest_sha256,
            sortmerna_manifest_sha256=sortmerna_result.value.manifest_sha256,
            source_manifest_sha256=source_manifest_sha256,
            source_identity_sha256=source["source_identity_sha256"],
            provenance_manifest_sha256=hashlib.sha256(provenance_bytes).hexdigest(),
            provenance_identity_sha256=provenance["provenance_identity_sha256"],
        )
        validation = BulkRnaSeqWorkflowAdapter().validate(
            WorkflowInputs(
                config=acceptance["workflow_inputs"]["config"],
                samples=acceptance["workflow_inputs"]["samples"],
                options=acceptance["workflow_inputs"]["options"],
            )
        )
        if validation.is_failure:
            raise ValueError("assembled workflow inputs fail the production adapter")

        acceptance_bytes = _canonical_json(acceptance)
        provenance_path = assembly.fixture_root / INDEX_PROVENANCE_FILENAME
        outputs = (
            (assembly.output_manifest, acceptance_bytes),
            (provenance_path, provenance_bytes),
        )
        _require_sidecars_absent_or_identical(outputs, label="output manifest")
        for path, content in outputs:
            if not path.exists():
                _write_new(path, content)
                created_outputs.append(path)
        return FinalizedFixture(
            acceptance_manifest=assembly.output_manifest,
            acceptance_manifest_sha256=_sha256_file(assembly.output_manifest),
            provenance_manifest=provenance_path,
            provenance_manifest_sha256=_sha256_file(provenance_path),
            provenance_identity_sha256=provenance["provenance_identity_sha256"],
        )
    except BaseException:
        for path in reversed(created_outputs):
            path.unlink(missing_ok=True)
        for path in reversed(created_sidecars):
            path.unlink(missing_ok=True)
        raise


def _new_output_root(value: Path) -> Path:
    if not isinstance(value, Path) or not value.is_absolute():
        raise ValueError("output root must be an absolute Path")
    if value != value.resolve(strict=False):
        raise ValueError("output root must be canonical and contain no symlinks")
    if value.exists():
        raise ValueError("output root must not already exist")
    if not value.parent.is_dir() or value.parent.is_symlink():
        raise ValueError("output root parent must be an existing regular directory")
    if value in {Path("/"), Path("/tmp"), PROJECT_ROOT, Path.home().resolve()}:
        raise ValueError("output root is unsafe")
    return value


def _random_dna(rng: random.Random, length: int) -> str:
    return "".join(rng.choice("ACGT") for _ in range(length))


def _fasta(name: str, sequence: str) -> str:
    return f">{name}\n{sequence}\n"


def _gtf() -> str:
    lines = []
    for gene in GENES:
        common = (
            f'gene_id "{gene.gene_id}"; '
            f'transcript_id "{gene.transcript_id}"; '
            f'gene_name "{gene.gene_id}"; '
            'gene_biotype "protein_coding"; '
            'transcript_biotype "protein_coding";'
        )
        gene_attributes = (
            f'gene_id "{gene.gene_id}"; gene_name "{gene.gene_id}"; '
            'gene_biotype "protein_coding";'
        )
        lines.extend(
            (
                _gtf_line(gene, "gene", gene_attributes),
                _gtf_line(gene, "transcript", common),
                _gtf_line(gene, "exon", common + ' exon_number "1";'),
            )
        )
    return "".join(lines)


def _gtf_line(gene: Gene, feature: str, attributes: str) -> str:
    return (
        f"chrTiny\tHelixWeave\t{feature}\t{gene.start}\t{gene.end}\t.\t"
        f"{gene.strand}\t.\t{attributes}\n"
    )


def _sample_rows() -> list[dict[str, str]]:
    return [
        {
            "sample": "PE1",
            "library": "libPE",
            "lane": "L001",
            "layout": "PE",
            "fastq_1": "reads/PE1_L001_R1.fastq.gz",
            "fastq_2": "reads/PE1_L001_R2.fastq.gz",
            "strandedness": "auto",
            "platform": "ILLUMINA",
        },
        {
            "sample": "PE1",
            "library": "libPE",
            "lane": "L002",
            "layout": "PE",
            "fastq_1": "reads/PE1_L002_R1.fastq.gz",
            "fastq_2": "reads/PE1_L002_R2.fastq.gz",
            "strandedness": "auto",
            "platform": "ILLUMINA",
        },
        {
            "sample": "SE1",
            "library": "libSE",
            "lane": "L001",
            "layout": "SE",
            "fastq_1": "reads/SE1_L001_R1.fastq.gz",
            "strandedness": "auto",
            "platform": "ILLUMINA",
        },
    ]


def _reads_for_row(
    row: Mapping[str, str],
    *,
    transcripts: Mapping[str, str],
    rrna: str,
) -> tuple[list[tuple[str, str]], list[tuple[str, str]] | None]:
    sample = row["sample"]
    lane = row["lane"]
    layout = row["layout"]
    lane_offset = 0 if lane == "L001" else 173
    read1: list[tuple[str, str]] = []
    read2: list[tuple[str, str]] | None = [] if layout == "PE" else None
    targets = (
        (transcripts["TX_PE_A"], transcripts["TX_PE_B"])
        if layout == "PE"
        else (transcripts["TX_SE"],)
    )
    for index in range(READS_PER_FASTQ):
        name = f"{sample}.{lane}.{index + 1:06d}"
        if index < RRNA_READS_PER_FASTQ:
            if layout == "PE":
                fragment = _window(rrna, index * 31 + lane_offset, 260)
                first = fragment[:READ_LENGTH]
                second = _reverse_complement(fragment[-READ_LENGTH:])
            else:
                first = _window(rrna, index * 31 + lane_offset, READ_LENGTH)
                second = None
        elif index < RRNA_READS_PER_FASTQ + ADAPTER_READS_PER_FASTQ:
            target = targets[index % len(targets)]
            insert = _window(target, index * 37 + lane_offset, 70)
            first = (insert + ADAPTER_SEQUENCE)[:READ_LENGTH]
            second = (
                (_reverse_complement(insert) + ADAPTER_SEQUENCE)[:READ_LENGTH]
                if layout == "PE"
                else None
            )
        elif layout == "PE":
            target = targets[index % len(targets)]
            fragment = _window(target, index * 43 + lane_offset, 260)
            first = fragment[:READ_LENGTH]
            second = _reverse_complement(fragment[-READ_LENGTH:])
        else:
            target = targets[index % len(targets)]
            first = _window(target, index * 43 + lane_offset, READ_LENGTH)
            second = None
        read1.append((f"{name}/1", first))
        if read2 is not None:
            assert second is not None
            read2.append((f"{name}/2", second))
    return read1, read2


def _window(sequence: str, offset: int, length: int) -> str:
    start = offset % (len(sequence) - length + 1)
    return sequence[start : start + length]


def _reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGT", "TGCA"))[::-1]


def _write_deterministic_gzip_fastq(
    path: Path,
    records: Sequence[tuple[str, str]],
) -> None:
    raw = "".join(
        f"@{name}\n{sequence}\n+\n{'I' * len(sequence)}\n" for name, sequence in records
    ).encode("ascii")
    compressor = zlib.compressobj(level=9, wbits=-zlib.MAX_WBITS)
    compressed = compressor.compress(raw) + compressor.flush()
    header = b"\x1f\x8b\x08\x00\x00\x00\x00\x00\x02\xff"
    footer = struct.pack("<II", zlib.crc32(raw) & 0xFFFFFFFF, len(raw) & 0xFFFFFFFF)
    _write_new(path, header + compressed + footer)


def _write_ascii(path: Path, content: str) -> None:
    _write_new(path, content.encode("ascii"))


def _write_new(path: Path, content: bytes) -> None:
    created = False
    try:
        with path.open("xb") as stream:
            created = True
            stream.write(content)
    except BaseException:
        if created:
            path.unlink(missing_ok=True)
        raise


def _canonical_json(value: object) -> bytes:
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


def _identity(scheme: str, payload: Mapping[str, object]) -> str:
    digest = hashlib.sha256()
    _frame(digest, scheme.encode("ascii"))
    _frame(digest, _canonical_json(payload))
    return digest.hexdigest()


def _frame(digest: Any, value: bytes) -> None:
    digest.update(len(value).to_bytes(8, "big"))
    digest.update(value)


def _sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        while chunk := stream.read(1024 * 1024):
            digest.update(chunk)
    return digest.hexdigest()


def _file_entry(root: Path, relative: str) -> dict[str, object]:
    path = root / relative
    return {
        "path": relative,
        "size_bytes": path.stat().st_size,
        "sha256": _sha256_file(path),
    }


def _index_build_contracts() -> dict[str, object]:
    return {
        "star": {
            "producer_process": "STAR_GENOMEGENERATE",
            "tool": "star",
            "tool_version": STAR_TOOL_VERSION,
            "build_contract": "nfcore-rnaseq-3.26.0-star-genomegenerate-default-v1",
            "sjdb_overhang": STAR_SJDB_OVERHANG,
            "genome_sa_index_nbases_strategy": "nfcore-rnaseq-3.26.0-auto",
        },
        "salmon": {
            "producer_process": "SALMON_INDEX",
            "tool": "salmon",
            "tool_version": SALMON_TOOL_VERSION,
            "build_contract": "nfcore-rnaseq-3.26.0-salmon-index-default-v1",
            "decoy_mode": "gentrome",
            "kmer_size": SALMON_KMER_SIZE,
        },
        "sortmerna": {
            "tool": "sortmerna",
            "tool_version": SORTMERNA_VERSION,
            "database_entries": 1,
        },
    }


def _source_reference() -> dict[str, object]:
    return {
        "contig": "chrTiny",
        "genome_length": GENOME_LENGTH,
        "rrna_length": RRNA_LENGTH,
        "gene_ids": [gene.gene_id for gene in GENES],
        "transcript_ids": [gene.transcript_id for gene in GENES],
    }


def _validate_assembly(value: FixtureAssembly) -> None:
    if not isinstance(value, FixtureAssembly):
        raise ValueError("assembly must be a FixtureAssembly")
    for name in (
        "fixture_root",
        "star_index_root",
        "salmon_index_root",
        "sortmerna_index_root",
        "output_manifest",
    ):
        path = getattr(value, name)
        if (
            not isinstance(path, Path)
            or not path.is_absolute()
            or path != path.resolve(strict=False)
        ):
            raise ValueError(f"{name} must be a canonical absolute Path")
    if not value.fixture_root.is_dir() or value.fixture_root.is_symlink():
        raise ValueError("fixture_root must be a regular directory")
    index_roots = (
        value.star_index_root,
        value.salmon_index_root,
        value.sortmerna_index_root,
    )
    if len(set(index_roots)) != 3:
        raise ValueError("index roots must be distinct")
    for root in index_roots:
        if value.fixture_root not in root.parents:
            raise ValueError("index roots must be controlled fixture descendants")
        if not root.is_dir() or root.is_symlink():
            raise ValueError("every index root must be a regular directory")
    if value.output_manifest != value.fixture_root / FIXTURE_MANIFEST_FILENAME:
        raise ValueError("output_manifest must use the controlled fixture filename")
    expected_versions = {
        "star_tool_version": STAR_TOOL_VERSION,
        "salmon_tool_version": SALMON_TOOL_VERSION,
        "sortmerna_tool_version": SORTMERNA_VERSION,
    }
    for name, expected in expected_versions.items():
        if getattr(value, name) != expected:
            raise ValueError(f"{name} does not match the production contract")
    for name in (
        "star_container_image",
        "salmon_container_image",
        "sortmerna_container_image",
    ):
        if _IMMUTABLE_IMAGE.fullmatch(getattr(value, name)) is None:
            raise ValueError(f"{name} must be an immutable image coordinate")
    for name in (
        "star_config_sha256",
        "salmon_config_sha256",
        "sortmerna_config_sha256",
    ):
        if _HEX64.fullmatch(getattr(value, name)) is None:
            raise ValueError(f"{name} must be a SHA-256 digest")
    for name in (
        "star_command_argv",
        "salmon_command_argv",
        "sortmerna_command_argv",
    ):
        _require_path_free_argv(getattr(value, name), fixture_root=value.fixture_root)
    if (
        isinstance(value.genome_sa_index_nbases, bool)
        or not isinstance(value.genome_sa_index_nbases, int)
        or not 1 <= value.genome_sa_index_nbases <= 14
    ):
        raise ValueError("genome_sa_index_nbases must be between 1 and 14")


def _require_path_free_argv(argv: object, *, fixture_root: Path) -> None:
    if (
        isinstance(argv, str)
        or not isinstance(argv, Sequence)
        or not argv
        or len(argv) > 512
    ):
        raise ValueError("generation command argv must be a non-empty sequence")
    total = 0
    for token in argv:
        if not isinstance(token, str) or not token or len(token) > 4096:
            raise ValueError("generation command argv contains an invalid token")
        total += len(token.encode("utf-8"))
        if (
            total > 64 * 1024
            or str(fixture_root) in token
            or "\x00" in token
            or "\n" in token
            or "\r" in token
            or _WINDOWS_ABSOLUTE.match(token)
            or _EMBEDDED_ABSOLUTE.search(token)
            or token.startswith("file:")
        ):
            raise ValueError("generation command provenance must be path-free")


def _verify_source_fixture(root: Path) -> dict[str, Any]:
    manifest_path = root / SOURCE_MANIFEST_FILENAME
    try:
        manifest_status = manifest_path.lstat()
        if (
            not stat.S_ISREG(manifest_status.st_mode)
            or manifest_status.st_nlink != 1
            or manifest_status.st_size <= 0
            or manifest_status.st_size > 1024 * 1024
        ):
            raise ValueError
        raw = manifest_path.read_bytes()
        source = json.loads(raw)
    except (OSError, UnicodeError, ValueError, json.JSONDecodeError):
        raise ValueError("source fixture closure is invalid") from None
    if (
        not isinstance(source, dict)
        or set(source) != _SOURCE_MANIFEST_KEYS
        or source.get("schema_version") != SOURCE_MANIFEST_SCHEMA_VERSION
        or source.get("fixture_id") != FIXTURE_ID
        or source.get("biological_validity") is not False
        or source.get("limitations") != SOURCE_LIMITATIONS
        or source.get("reference") != _source_reference()
        or source.get("generator")
        != {
            "seed": SEED,
            "read_length": READ_LENGTH,
            "reads_per_fastq": READS_PER_FASTQ,
            "rrna_reads_per_fastq": RRNA_READS_PER_FASTQ,
            "adapter_reads_per_fastq": ADAPTER_READS_PER_FASTQ,
        }
        or source.get("samples") != _sample_rows()
        or source.get("index_build_contracts") != _index_build_contracts()
    ):
        raise ValueError("source fixture closure is invalid")
    identity_payload = dict(source)
    observed_identity = identity_payload.pop("source_identity_sha256", None)
    if not isinstance(observed_identity, str) or observed_identity != _identity(
        SOURCE_IDENTITY_SCHEME, identity_payload
    ):
        raise ValueError("source fixture closure is invalid")
    files = source.get("files")
    if not isinstance(files, list) or [
        item.get("path") if isinstance(item, dict) else None for item in files
    ] != list(_SOURCE_RELATIVE_FILES):
        raise ValueError("source fixture closure is invalid")
    observed_paths: list[str] = []
    for item in files:
        if not isinstance(item, dict) or set(item) != {"path", "size_bytes", "sha256"}:
            raise ValueError("source fixture closure is invalid")
        relative = item["path"]
        if not isinstance(relative, str):
            raise ValueError("source fixture closure is invalid")
        pure = PurePosixPath(relative)
        if (
            pure.is_absolute()
            or not pure.parts
            or any(part in {"", ".", ".."} for part in pure.parts)
        ):
            raise ValueError("source fixture closure is invalid")
        path = root.joinpath(*pure.parts)
        try:
            parent = root
            for part in pure.parts[:-1]:
                parent /= part
                parent_status = parent.lstat()
                if not stat.S_ISDIR(parent_status.st_mode):
                    raise ValueError
            status = path.lstat()
        except (OSError, ValueError):
            raise ValueError("source fixture closure is invalid") from None
        if (
            not stat.S_ISREG(status.st_mode)
            or status.st_nlink != 1
            or isinstance(item["size_bytes"], bool)
            or not isinstance(item["size_bytes"], int)
            or item["size_bytes"] != status.st_size
            or not isinstance(item["sha256"], str)
            or _HEX64.fullmatch(item["sha256"]) is None
            or item["sha256"] != _sha256_file(path)
        ):
            raise ValueError("source fixture closure is invalid")
        observed_paths.append(relative)
    if observed_paths != sorted(observed_paths) or len(observed_paths) != len(
        set(observed_paths)
    ):
        raise ValueError("source fixture closure is invalid")
    return source


def _source_file_digest(source: Mapping[str, Any], relative: str) -> str:
    for item in source["files"]:
        if item["path"] == relative:
            return item["sha256"]
    raise ValueError("source fixture closure is incomplete")


def _index_files(root: Path, *, excluded: str) -> list[dict[str, object]]:
    entries: list[dict[str, object]] = []
    directories: set[str] = set()
    for path in sorted(
        root.rglob("*"), key=lambda item: item.relative_to(root).as_posix()
    ):
        relative = path.relative_to(root).as_posix()
        status = path.lstat()
        if stat.S_ISLNK(status.st_mode):
            raise ValueError("index trees must not contain symlinks")
        if stat.S_ISDIR(status.st_mode):
            directories.add(relative)
            continue
        if not stat.S_ISREG(status.st_mode) or status.st_nlink != 1:
            raise ValueError("index trees must contain only regular private files")
        if relative == excluded:
            continue
        entries.append(
            {
                "path": relative,
                "size_bytes": status.st_size,
                "sha256": _sha256_file(path),
            }
        )
    if not entries:
        raise ValueError("index tree must contain real generated files")
    required_directories: set[str] = set()
    for entry in entries:
        parent = PurePosixPath(entry["path"]).parent
        while parent.as_posix() != ".":
            required_directories.add(parent.as_posix())
            parent = parent.parent
    if directories != required_directories:
        raise ValueError("index trees must not contain empty directories")
    return entries


def _reference_index_sidecar(
    *,
    kind: str,
    container_image: str,
    fasta_sha256: str,
    gtf_sha256: str,
    transcript_fasta_sha256: str,
    files: list[dict[str, object]],
    genome_sa_index_nbases: int,
) -> bytes:
    if kind == "star":
        producer = {
            "process": "STAR_GENOMEGENERATE",
            "tool": "star",
            "tool_version": STAR_TOOL_VERSION,
            "build_contract": "nfcore-rnaseq-3.26.0-star-genomegenerate-default-v1",
            "container_image": container_image,
        }
        build_parameters = {
            "extra_args": [],
            "genome_sa_index_nbases": genome_sa_index_nbases,
            "genome_sa_index_nbases_strategy": "nfcore-rnaseq-3.26.0-auto",
            "sjdb_gtf_feature": "exon",
            "sjdb_overhang": STAR_SJDB_OVERHANG,
            "skip_gtf_filter": False,
            "skip_gtf_transcript_filter": False,
        }
    elif kind == "salmon":
        producer = {
            "process": "SALMON_INDEX",
            "tool": "salmon",
            "tool_version": SALMON_TOOL_VERSION,
            "build_contract": "nfcore-rnaseq-3.26.0-salmon-index-default-v1",
            "container_image": container_image,
        }
        build_parameters = {
            "decoy_mode": "gentrome",
            "extra_args": [],
            "gencode": False,
            "gffread_transcript_fasta": False,
            "kmer_size": SALMON_KMER_SIZE,
            "skip_gtf_filter": False,
            "skip_gtf_transcript_filter": False,
        }
    else:
        raise ValueError("unsupported reference index kind")
    reference = {
        "fasta_sha256": fasta_sha256,
        "gtf_sha256": gtf_sha256,
    }
    if kind == "salmon":
        reference["transcript_fasta_sha256"] = transcript_fasta_sha256
    return _canonical_json(
        {
            "schema_version": REFERENCE_INDEX_SCHEMA_VERSION,
            "index_kind": kind,
            "producer": producer,
            "build_parameters": build_parameters,
            "reference": reference,
            "files": files,
        }
    )


def _nfcore_star_genome_sa_index_nbases(genome_length: int) -> int:
    """Mirror the pinned STAR_GENOMEGENERATE gawk rounding formula."""
    if isinstance(genome_length, bool) or not isinstance(genome_length, int):
        raise ValueError("genome length is invalid")
    if genome_length <= 0:
        raise ValueError("genome length is invalid")
    unbounded = (math.log(genome_length) / math.log(2)) / 2 - 1
    return min(14, int(math.floor(unbounded + 0.5)))


def _require_sidecars_absent_or_identical(
    documents: Sequence[tuple[Path, bytes]],
    *,
    label: str = "sidecar",
) -> None:
    for path, expected in documents:
        if not path.exists():
            continue
        try:
            status = path.lstat()
            observed = path.read_bytes()
        except OSError:
            raise ValueError(f"existing {label} is invalid") from None
        if (
            not stat.S_ISREG(status.st_mode)
            or status.st_nlink != 1
            or observed != expected
        ):
            raise ValueError(f"existing {label} conflicts with observed closure")


def _acceptance_document(
    *,
    assembly: FixtureAssembly,
    source: Mapping[str, Any],
    reference: Mapping[str, object],
    database_manifest_sha256: str,
    sortmerna_manifest_sha256: str,
    source_manifest_sha256: str,
    source_identity_sha256: str,
    provenance_manifest_sha256: str,
    provenance_identity_sha256: str,
) -> dict[str, object]:
    samples = []
    for source_row in source["samples"]:
        row = dict(source_row)
        row["fastq_1"] = str(
            assembly.fixture_root.joinpath(*PurePosixPath(row["fastq_1"]).parts)
        )
        if "fastq_2" in row:
            row["fastq_2"] = str(
                assembly.fixture_root.joinpath(*PurePosixPath(row["fastq_2"]).parts)
            )
        samples.append(row)
    required_artifact_sample_output_types = sorted(
        {
            ("PE1", "bulk_rnaseq.fastqc.raw.read1.zip"),
            ("PE1", "bulk_rnaseq.fastqc.raw.read2.zip"),
            ("PE1", "bulk_rnaseq.rrna.sortmerna.filtered.read1"),
            ("PE1", "bulk_rnaseq.rrna.sortmerna.filtered.read2"),
            ("PE1", "bulk_rnaseq.salmon.meta_info"),
            ("PE1", "bulk_rnaseq.salmon.quant_gene"),
            ("PE1", "bulk_rnaseq.star.bam"),
            ("PE1", "bulk_rnaseq.star.log_final"),
            ("SE1", "bulk_rnaseq.fastqc.raw.single.zip"),
            ("SE1", "bulk_rnaseq.rrna.sortmerna.filtered.single"),
            ("SE1", "bulk_rnaseq.salmon.meta_info"),
            ("SE1", "bulk_rnaseq.salmon.quant_gene"),
            ("SE1", "bulk_rnaseq.star.bam"),
            ("SE1", "bulk_rnaseq.star.log_final"),
        }
    )
    required_qc_sample_metric_keys = sorted(
        {
            ("PE1", "fastqc.raw.read1.total_sequences"),
            ("PE1", "fastqc.raw.read2.total_sequences"),
            ("PE1", "salmon.mapping_fraction"),
            ("PE1", "salmon.processed_fragments"),
            ("PE1", "star.input_templates"),
            ("PE1", "star.uniquely_mapped_template_fraction"),
            ("PE1", "trimming.read1.retained_reads"),
            ("PE1", "trimming.read2.retained_reads"),
            ("SE1", "fastqc.raw.single.total_sequences"),
            ("SE1", "salmon.mapping_fraction"),
            ("SE1", "salmon.processed_fragments"),
            ("SE1", "star.input_templates"),
            ("SE1", "star.uniquely_mapped_template_fraction"),
            ("SE1", "trimming.single.retained_reads"),
        }
    )
    sample_lane_counts = {
        sample_id: sum(1 for row in source["samples"] if row["sample"] == sample_id)
        for sample_id in ("PE1", "SE1")
    }
    required_qc_sample_metric_values = sorted(
        {
            (
                "PE1",
                "fastqc.raw.read1.total_sequences",
                str(READS_PER_FASTQ * sample_lane_counts["PE1"]),
            ),
            (
                "PE1",
                "fastqc.raw.read2.total_sequences",
                str(READS_PER_FASTQ * sample_lane_counts["PE1"]),
            ),
            (
                "SE1",
                "fastqc.raw.single.total_sequences",
                str(READS_PER_FASTQ * sample_lane_counts["SE1"]),
            ),
        }
    )
    return {
        "schema_version": FIXTURE_MANIFEST_SCHEMA_VERSION,
        "source_manifest_sha256": source_manifest_sha256,
        "source_identity_sha256": source_identity_sha256,
        "index_provenance_manifest_sha256": provenance_manifest_sha256,
        "index_provenance_identity_sha256": provenance_identity_sha256,
        "transcriptome_binding": {
            "reference_id": reference["reference_id"],
            "fasta_sha256": reference["fasta_sha256"],
            "gtf_sha256": reference["gtf_sha256"],
            "transcript_fasta": str(assembly.fixture_root / "reference/transcripts.fa"),
            "transcript_fasta_sha256": _source_file_digest(
                source,
                "reference/transcripts.fa",
            ),
        },
        "workflow_inputs": {
            "config": {
                "standard": {
                    "reference": dict(reference),
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
                            "path": str(
                                assembly.fixture_root / "rrna/database-manifest.txt"
                            ),
                            "identity_sha256": database_manifest_sha256,
                        },
                        "sortmerna_index": {
                            "path": str(assembly.sortmerna_index_root),
                            "identity_sha256": sortmerna_manifest_sha256,
                        },
                    },
                    "qc": {
                        "enabled": True,
                        "fastqc": True,
                        "multiqc": True,
                        "rseqc": False,
                        "qualimap": False,
                        "dupradar": False,
                        "biotype": False,
                        "deseq2_pca": False,
                        "preseq": False,
                        "mark_duplicates": False,
                    },
                    "outputs": {
                        "bigwig": False,
                        "stringtie": False,
                        "trimmed_reads": False,
                        "unaligned_reads": False,
                        "reference_files": False,
                        "alignment_intermediates": False,
                        "merged_fastq": False,
                        "umi_intermediates": False,
                    },
                }
            },
            "samples": samples,
            "options": {},
        },
        "required_artifact_output_types": sorted(
            {output_type for _, output_type in required_artifact_sample_output_types}
        ),
        "required_artifact_sample_output_types": [
            [sample_id, output_type]
            for sample_id, output_type in required_artifact_sample_output_types
        ],
        "required_qc_metric_keys": sorted(
            {metric_key for _, metric_key in required_qc_sample_metric_keys}
        ),
        "required_qc_sample_metric_keys": [
            [sample_id, metric_key]
            for sample_id, metric_key in required_qc_sample_metric_keys
        ],
        "required_qc_sample_metric_values": [
            [sample_id, metric_key, value]
            for sample_id, metric_key, value in required_qc_sample_metric_values
        ],
        "required_sample_ids": ["PE1", "SE1"],
    }


def _provenance_document(
    *,
    assembly: FixtureAssembly,
    source_manifest_sha256: str,
    source_identity_sha256: str,
    reference_closure: Any,
    sortmerna_closure: Any,
) -> dict[str, Any]:
    star = reference_closure.star_index
    salmon = reference_closure.salmon_index
    if star is None or salmon is None:
        raise ValueError("reference closure is incomplete")
    indexes = {
        "star": _index_provenance(
            root=assembly.star_index_root,
            fixture_root=assembly.fixture_root,
            process="STAR_GENOMEGENERATE",
            tool="star",
            tool_version=assembly.star_tool_version,
            image=assembly.star_container_image,
            config_sha256=assembly.star_config_sha256,
            argv=assembly.star_command_argv,
            manifest_sha256=star.manifest_sha256,
            closure_sha256=star.identity_sha256,
        ),
        "salmon": _index_provenance(
            root=assembly.salmon_index_root,
            fixture_root=assembly.fixture_root,
            process="SALMON_INDEX",
            tool="salmon",
            tool_version=assembly.salmon_tool_version,
            image=assembly.salmon_container_image,
            config_sha256=assembly.salmon_config_sha256,
            argv=assembly.salmon_command_argv,
            manifest_sha256=salmon.manifest_sha256,
            closure_sha256=salmon.identity_sha256,
        ),
        "sortmerna": _index_provenance(
            root=assembly.sortmerna_index_root,
            fixture_root=assembly.fixture_root,
            process="SORTMERNA",
            tool="sortmerna",
            tool_version=assembly.sortmerna_tool_version,
            image=assembly.sortmerna_container_image,
            config_sha256=assembly.sortmerna_config_sha256,
            argv=assembly.sortmerna_command_argv,
            manifest_sha256=sortmerna_closure.manifest_sha256,
            closure_sha256=sortmerna_closure.identity_sha256,
        ),
    }
    payload: dict[str, Any] = {
        "schema_version": INDEX_PROVENANCE_SCHEMA_VERSION,
        "fixture_id": FIXTURE_ID,
        "biological_validity": False,
        "source_manifest_sha256": source_manifest_sha256,
        "source_identity_sha256": source_identity_sha256,
        "indexes": indexes,
    }
    payload["provenance_identity_sha256"] = _identity(
        PROVENANCE_IDENTITY_SCHEME,
        payload,
    )
    return payload


def _index_provenance(
    *,
    root: Path,
    fixture_root: Path,
    process: str,
    tool: str,
    tool_version: str,
    image: str,
    config_sha256: str,
    argv: tuple[str, ...],
    manifest_sha256: str,
    closure_sha256: str,
) -> dict[str, object]:
    command = list(argv)
    return {
        "relative_index_root": root.relative_to(fixture_root).as_posix(),
        "producer_process": process,
        "tool": tool,
        "tool_version": tool_version,
        "container_image": image,
        "immutable_config_sha256": config_sha256,
        "command_argv": command,
        "command_argv_sha256": hashlib.sha256(_canonical_json(command)).hexdigest(),
        "sidecar_manifest_sha256": manifest_sha256,
        "index_closure_sha256": closure_sha256,
    }


def _json_argv(value: str) -> tuple[str, ...]:
    try:
        document = json.loads(value)
    except json.JSONDecodeError as error:
        raise argparse.ArgumentTypeError("must be a JSON argv array") from error
    if (
        not isinstance(document, list)
        or not document
        or any(not isinstance(item, str) for item in document)
    ):
        raise argparse.ArgumentTypeError("must be a non-empty JSON string array")
    return tuple(document)


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    generate = subparsers.add_parser("generate")
    generate.add_argument("--output-root", required=True, type=Path)

    finalize = subparsers.add_parser("finalize")
    for name in (
        "fixture-root",
        "star-index-root",
        "salmon-index-root",
        "sortmerna-index-root",
        "output-manifest",
    ):
        finalize.add_argument(f"--{name}", required=True, type=Path)
    for name in (
        "star-container-image",
        "salmon-container-image",
        "sortmerna-container-image",
        "star-tool-version",
        "salmon-tool-version",
        "sortmerna-tool-version",
        "star-config-sha256",
        "salmon-config-sha256",
        "sortmerna-config-sha256",
    ):
        finalize.add_argument(f"--{name}", required=True)
    for name in (
        "star-command-json",
        "salmon-command-json",
        "sortmerna-command-json",
    ):
        finalize.add_argument(f"--{name}", required=True, type=_json_argv)
    finalize.add_argument("--genome-sa-index-nbases", required=True, type=int)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    parser = _parser()
    arguments = parser.parse_args(argv)
    try:
        if arguments.command == "generate":
            generated = generate_fixture(arguments.output_root)
            print(
                json.dumps(
                    {
                        "fixture_root": str(generated.root),
                        "source_manifest": str(generated.manifest),
                        "source_manifest_sha256": generated.manifest_sha256,
                        "source_identity_sha256": generated.source_identity_sha256,
                    },
                    sort_keys=True,
                )
            )
            return 0
        assembly = FixtureAssembly(
            fixture_root=arguments.fixture_root,
            star_index_root=arguments.star_index_root,
            salmon_index_root=arguments.salmon_index_root,
            sortmerna_index_root=arguments.sortmerna_index_root,
            output_manifest=arguments.output_manifest,
            star_container_image=arguments.star_container_image,
            salmon_container_image=arguments.salmon_container_image,
            sortmerna_container_image=arguments.sortmerna_container_image,
            star_tool_version=arguments.star_tool_version,
            salmon_tool_version=arguments.salmon_tool_version,
            sortmerna_tool_version=arguments.sortmerna_tool_version,
            star_config_sha256=arguments.star_config_sha256,
            salmon_config_sha256=arguments.salmon_config_sha256,
            sortmerna_config_sha256=arguments.sortmerna_config_sha256,
            star_command_argv=arguments.star_command_json,
            salmon_command_argv=arguments.salmon_command_json,
            sortmerna_command_argv=arguments.sortmerna_command_json,
            genome_sa_index_nbases=arguments.genome_sa_index_nbases,
        )
        finalized = finalize_fixture(assembly)
        print(
            json.dumps(
                {
                    "acceptance_manifest": str(finalized.acceptance_manifest),
                    "acceptance_manifest_sha256": (
                        finalized.acceptance_manifest_sha256
                    ),
                    "provenance_manifest": str(finalized.provenance_manifest),
                    "provenance_manifest_sha256": (
                        finalized.provenance_manifest_sha256
                    ),
                    "provenance_identity_sha256": (
                        finalized.provenance_identity_sha256
                    ),
                },
                sort_keys=True,
            )
        )
        return 0
    except ValueError as error:
        parser.error(str(error))
    raise AssertionError("unreachable")


if __name__ == "__main__":
    raise SystemExit(main())
