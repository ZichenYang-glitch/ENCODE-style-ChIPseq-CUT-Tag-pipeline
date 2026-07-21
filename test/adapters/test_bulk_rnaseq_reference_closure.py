"""Reference and prebuilt-index closure tests for the bulk RNA-seq runtime."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import pytest

from encode_pipeline.adapters.bulk_rnaseq.reference_closure import (
    REFERENCE_INDEX_MANIFEST,
    verify_reference_closure,
    verify_reference_index,
)
from encode_pipeline.adapters.bulk_rnaseq.resource_closure import ResourceClosurePolicy


_PRODUCER_IMAGES = {
    "STAR_GENOMEGENERATE": (
        "community.wave.seqera.io/library/htslib_samtools_star_gawk:"
        + "ae438e9a604351a4@sha256:"
        + "a" * 64
    ),
    "SALMON_INDEX": (
        "quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2@sha256:" + "b" * 64
    ),
}
_TRANSCRIPT_FASTA_SHA256 = hashlib.sha256(b">tx1\nACGT\n").hexdigest()


def _sha256(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def _write_file(path: Path, value: bytes) -> str:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(value)
    return _sha256(value)


def _write_index(
    root: Path,
    *,
    kind: str,
    fasta_sha256: str,
    gtf_sha256: str,
    annotation_style: str = "ensembl",
) -> str:
    entries = []
    for relative, content in (("index/a.bin", b"a"), ("index/b.bin", b"bb")):
        digest = _write_file(root / relative, content)
        entries.append({"path": relative, "size_bytes": len(content), "sha256": digest})
    reference = {
        "fasta_sha256": fasta_sha256,
        "gtf_sha256": gtf_sha256,
    }
    if kind == "salmon":
        reference["transcript_fasta_sha256"] = _TRANSCRIPT_FASTA_SHA256
    manifest = {
        "schema_version": "1.1.0",
        "index_kind": kind,
        "producer": (
            {
                "process": "STAR_GENOMEGENERATE",
                "tool": "star",
                "tool_version": "2.7.11b",
                "build_contract": (
                    "nfcore-rnaseq-3.26.0-star-genomegenerate-default-v1"
                ),
                "container_image": _PRODUCER_IMAGES["STAR_GENOMEGENERATE"],
            }
            if kind == "star"
            else {
                "process": "SALMON_INDEX",
                "tool": "salmon",
                "tool_version": "1.10.3",
                "build_contract": "nfcore-rnaseq-3.26.0-salmon-index-default-v1",
                "container_image": _PRODUCER_IMAGES["SALMON_INDEX"],
            }
        ),
        "build_parameters": (
            {
                "extra_args": [],
                "genome_sa_index_nbases": 9,
                "genome_sa_index_nbases_strategy": "nfcore-rnaseq-3.26.0-auto",
                "sjdb_gtf_feature": "exon",
                "sjdb_overhang": 100,
                "skip_gtf_filter": False,
                "skip_gtf_transcript_filter": False,
            }
            if kind == "star"
            else {
                "decoy_mode": "gentrome",
                "extra_args": [],
                "gencode": annotation_style == "gencode",
                "gffread_transcript_fasta": False,
                "kmer_size": 31,
                "skip_gtf_filter": False,
                "skip_gtf_transcript_filter": False,
            }
        ),
        "reference": reference,
        "files": entries,
    }
    content = json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()
    return _write_file(root / REFERENCE_INDEX_MANIFEST, content)


def _reference(
    tmp_path: Path,
    *,
    indexes: bool = True,
    annotation_style: str = "ensembl",
) -> dict[str, object]:
    fasta = tmp_path / "reference/genome.fa"
    gtf = tmp_path / "reference/genes.gtf"
    fasta_sha256 = _write_file(fasta, b">chr1\nACGT\n")
    gtf_sha256 = _write_file(gtf, b'chr1\ttest\texon\t1\t4\t.\t+\t.\tgene_id "g1";\n')
    result: dict[str, object] = {
        "reference_id": "tiny-reference",
        "fasta": str(fasta),
        "fasta_sha256": fasta_sha256,
        "gtf": str(gtf),
        "gtf_sha256": gtf_sha256,
        "annotation_style": annotation_style,
    }
    if indexes:
        for kind in ("star", "salmon"):
            root = tmp_path / f"reference/{kind}"
            identity = _write_index(
                root,
                kind=kind,
                fasta_sha256=fasta_sha256,
                gtf_sha256=gtf_sha256,
                annotation_style=annotation_style,
            )
            result[f"{kind}_index"] = {
                "path": str(root),
                "identity_sha256": identity,
            }
    return result


def test_reference_closure_is_deterministic_and_binds_indexes(tmp_path: Path):
    reference = _reference(tmp_path)

    first = verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
    )
    second = verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
    )

    assert first.is_success
    assert second.is_success
    assert first.value.identity_sha256 == second.value.identity_sha256
    assert first.value.star_index is not None
    assert first.value.salmon_index is not None
    assert first.value.star_index.producer_tool_version == "2.7.11b"
    assert first.value.salmon_index.producer_tool_version == "1.10.3"
    assert (
        first.value.star_index.manifest_sha256
        == reference["star_index"]["identity_sha256"]
    )
    assert [item.relative_path for item in first.value.star_index.files] == [
        "index/a.bin",
        "index/b.bin",
    ]


def test_gencode_salmon_index_build_contract_is_accepted(tmp_path: Path):
    reference = _reference(tmp_path, annotation_style="gencode")

    result = verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
    )

    assert result.is_success
    assert result.value.salmon_index is not None
    assert ("gencode", True) in result.value.salmon_index.build_parameters


@pytest.mark.parametrize(
    "run_parameters",
    [
        {"skip_gtf_filter": True},
        {"skip_gtf_transcript_filter": True},
        {"gffread_transcript_fasta": True},
    ],
)
def test_reference_indexes_reject_run_shaping_parameter_mismatch(
    tmp_path: Path,
    run_parameters: dict[str, bool],
):
    reference = _reference(tmp_path)

    result = verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        index_build_parameters=run_parameters,
        transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_INVALID"


@pytest.mark.parametrize(
    ("section", "name", "value"),
    [
        ("producer", "tool_version", "2.7.10b"),
        (
            "producer",
            "container_image",
            "community.wave.seqera.io/library/star:other@sha256:" + "c" * 64,
        ),
        ("build_parameters", "sjdb_overhang", 99),
        ("build_parameters", "extra_args", ["--unsafe"]),
    ],
)
def test_reference_index_rejects_incompatible_producer_or_build_contract(
    tmp_path: Path,
    section: str,
    name: str,
    value: object,
):
    reference = _reference(tmp_path)
    star = reference["star_index"]
    manifest_path = Path(star["path"], REFERENCE_INDEX_MANIFEST)
    manifest = json.loads(manifest_path.read_bytes())
    manifest[section][name] = value
    content = json.dumps(manifest, sort_keys=True, separators=(",", ":")).encode()
    manifest_path.write_bytes(content)
    star["identity_sha256"] = _sha256(content)

    result = verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_INVALID"


def test_reference_primary_file_tamper_fails_closed_without_path(tmp_path: Path):
    reference = _reference(tmp_path, indexes=False)
    Path(reference["fasta"]).write_bytes(b"changed")

    result = verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_INVALID"
    assert str(tmp_path) not in json.dumps(result.to_dict())


def test_reference_index_rejects_reference_identity_mismatch(tmp_path: Path):
    reference = _reference(tmp_path)
    star = reference["star_index"]

    result = verify_reference_index(
        star["path"],
        kind="star",
        expected_manifest_sha256=star["identity_sha256"],
        fasta_sha256="f" * 64,
        gtf_sha256=reference["gtf_sha256"],
        annotation_style="ensembl",
        expected_container_image=_PRODUCER_IMAGES["STAR_GENOMEGENERATE"],
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_INDEX_INVALID"


def test_salmon_index_rejects_transcript_fasta_identity_mismatch(
    tmp_path: Path,
) -> None:
    reference = _reference(tmp_path)

    result = verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        transcript_fasta_sha256="f" * 64,
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_INVALID"


def test_reference_index_rejects_unlisted_extra_file(tmp_path: Path):
    reference = _reference(tmp_path)
    star = reference["star_index"]
    Path(star["path"], "extra.bin").write_bytes(b"extra")

    result = verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
    )

    assert result.is_failure


def test_reference_index_rejects_unlisted_empty_directory(tmp_path: Path):
    reference = _reference(tmp_path)
    star = reference["star_index"]
    Path(star["path"], "unlisted-empty").mkdir()

    assert verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
    ).is_failure


def test_reference_index_directory_depth_is_bounded(tmp_path: Path):
    reference = _reference(tmp_path)

    result = verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
        policy=ResourceClosurePolicy(maximum_index_depth=1),
    )

    assert result.is_failure
    assert result.errors[0].code == "BULK_RNASEQ_REFERENCE_INVALID"


def test_reference_index_rejects_symlink_and_fifo(tmp_path: Path):
    reference = _reference(tmp_path)
    star = reference["star_index"]
    index_root = Path(star["path"])
    external = tmp_path / "external.bin"
    external.write_bytes(b"external")
    (index_root / "linked.bin").symlink_to(external)

    assert verify_reference_closure(
        reference,
        producer_images=_PRODUCER_IMAGES,
        transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
    ).is_failure

    (index_root / "linked.bin").unlink()
    fifo = index_root / "named.pipe"
    os.mkfifo(fifo)
    try:
        assert verify_reference_closure(
            reference,
            producer_images=_PRODUCER_IMAGES,
            transcript_fasta_sha256=_TRANSCRIPT_FASTA_SHA256,
        ).is_failure
    finally:
        fifo.unlink()


def test_reference_index_manifest_rejects_traversal(tmp_path: Path):
    root = tmp_path / "index"
    root.mkdir()
    manifest = {
        "schema_version": "1.1.0",
        "index_kind": "star",
        "producer": {
            "process": "STAR_GENOMEGENERATE",
            "tool": "star",
            "tool_version": "2.7.11b",
            "build_contract": "nfcore-rnaseq-3.26.0-star-genomegenerate-default-v1",
            "container_image": _PRODUCER_IMAGES["STAR_GENOMEGENERATE"],
        },
        "build_parameters": {
            "extra_args": [],
            "genome_sa_index_nbases": 9,
            "genome_sa_index_nbases_strategy": "nfcore-rnaseq-3.26.0-auto",
            "sjdb_gtf_feature": "exon",
            "sjdb_overhang": 100,
            "skip_gtf_filter": False,
            "skip_gtf_transcript_filter": False,
        },
        "reference": {"fasta_sha256": "a" * 64, "gtf_sha256": "b" * 64},
        "files": [{"path": "../escape", "size_bytes": 1, "sha256": "c" * 64}],
    }
    content = json.dumps(manifest).encode()
    (root / REFERENCE_INDEX_MANIFEST).write_bytes(content)

    result = verify_reference_index(
        root,
        kind="star",
        expected_manifest_sha256=_sha256(content),
        fasta_sha256="a" * 64,
        gtf_sha256="b" * 64,
        annotation_style="ensembl",
        expected_container_image=_PRODUCER_IMAGES["STAR_GENOMEGENERATE"],
    )

    assert result.is_failure


def test_reference_index_rejects_relative_root(tmp_path: Path):
    with pytest.MonkeyPatch.context() as monkeypatch:
        monkeypatch.chdir(tmp_path)
        result = verify_reference_index(
            "index",
            kind="star",
            expected_manifest_sha256="a" * 64,
            fasta_sha256="b" * 64,
            gtf_sha256="c" * 64,
            annotation_style="ensembl",
            expected_container_image=_PRODUCER_IMAGES["STAR_GENOMEGENERATE"],
        )
    assert result.is_failure
