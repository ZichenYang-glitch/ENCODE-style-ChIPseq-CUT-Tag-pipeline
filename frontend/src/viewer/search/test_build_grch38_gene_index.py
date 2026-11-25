"""The gene-index generator keeps spans from the GTF and drops everything else."""

from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

SCRIPT = Path(__file__).with_name("build_grch38_gene_index.py")
README = """\
name: GRCh38
original name: hg38
assembly_accession: GCA_000001405.15
genome url: https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
annotation url: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.primary_assembly.annotation.gtf.gz
genomepy version: 0.16.3
date: 2025-12-19 16:30:21
"""
REPORT = """\
Sequence-Name\tUCSC-style-name
7\tchr7
X\tchrX
Y\tchrY
3\tchr3
9\tchr9
17\tchr17
MT\tchrM
"""
HGNC = """\
hgnc_id\tsymbol\tensembl_gene_id\talias_symbol\tprev_symbol
HGNC:1\tEGFR\tENSG00000000001\tERBB1|ERRP\tERBB
HGNC:2\tBHLHE40\tENSG00000000005\tDEC1\t
HGNC:3\tDELEC1\tENSG00000000006\tDEC1\t
HGNC:4\tGHOST\tENSG00000000999\tNOPE\t
HGNC:5\tEGFR\tENSG00000000001\tEGFR\t
HGNC:6\tTP53\tENSG00000000002\tp53\t
"""


def _gene(
    chrom: str,
    start: int,
    end: int,
    gene_id: str,
    name: str,
    *,
    feature: str = "gene",
    include_name: bool = True,
) -> str:
    name_attr = f' gene_name "{name}";' if include_name else ""
    attrs = f'gene_id "{gene_id}"; gene_type "protein_coding";{name_attr}'
    return f"{chrom}\tHAVANA\t{feature}\t{start}\t{end}\t.\t+\t.\t{attrs}\n"


def _fixture_gtf() -> str:
    return "".join(
        [
            "# comment\n",
            "\n",
            _gene("chr7", 10, 40, "ENSG00000000001.1", "EGFR"),
            _gene("chr7", 10, 20, "ENSG00000000011.1", "EXONONLY", feature="exon"),
            _gene("chrX", 5, 8, "ENSG00000000003.1", "IL9R"),
            _gene("chrY", 9, 12, "ENSG00000000004.1", "IL9R"),
            _gene("chr3", 1, 3, "ENSG00000000005.1", "BHLHE40"),
            _gene("chr9", 4, 7, "ENSG00000000006.1", "DELEC1"),
            _gene("chr17", 100, 110, "ENSG00000000002.1", "TP53"),
            _gene("chr9_KI270000v1_alt", 1, 2, "ENSG00000000007.1", "DROPME"),
            _gene("chr7", 1, 2, "ENSG00000000008.1", "NONAME", include_name=False),
            _gene("chr7", 10, 4, "ENSG00000000009.1", "BADEND"),
        ]
    )


def _write_inputs(
    root: Path, *, readme: str = README, hgnc: str = HGNC
) -> dict[str, Path]:
    gtf = root / "annotation.gtf"
    gtf.write_text(_fixture_gtf(), encoding="utf-8")
    report = root / "assembly_report.txt"
    report.write_text(REPORT, encoding="utf-8")
    readme_path = root / "README.txt"
    readme_path.write_text(readme, encoding="utf-8")
    hgnc_path = root / "hgnc_complete_set.txt"
    hgnc_path.write_text(hgnc, encoding="utf-8")
    return {
        "gtf": gtf,
        "report": report,
        "readme": readme_path,
        "hgnc": hgnc_path,
        "output": root / "hg38-gene-index.json",
        "provenance": root / "hg38-gene-index.provenance.json",
    }


def _run(
    paths: dict[str, Path], *, readme: bool = True
) -> subprocess.CompletedProcess[str]:
    command = [
        sys.executable,
        str(SCRIPT),
        "--gtf",
        str(paths["gtf"]),
        "--assembly-report",
        str(paths["report"]),
        "--hgnc",
        str(paths["hgnc"]),
        "--output",
        str(paths["output"]),
        "--provenance",
        str(paths["provenance"]),
    ]
    if readme:
        command.extend(["--readme", str(paths["readme"])])
    return subprocess.run(
        command,
        check=False,
        capture_output=True,
        text=True,
    )


def _loci(index: dict, symbol: str) -> list[tuple]:
    chromosomes = index["chromosomes"]
    return [
        (chromosomes[row[2]], row[3], row[4], row[1])
        for row in index["rows"]
        if row[0] == symbol
    ]


def test_index_uses_gtf_spans_and_hgnc_aliases(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    result = _run(paths)
    assert result.returncode == 0, result.stderr
    index = json.loads(paths["output"].read_text(encoding="utf-8"))
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    assert index["provenance"] == provenance
    assert index["coordinate_system"] == "1-based-closed"
    assert (
        index["provenance"]["annotation"]["sha256"]
        == hashlib.sha256(paths["gtf"].read_bytes()).hexdigest()
    )
    assert _loci(index, "EGFR") == [("chr7", 10, 40, "ENSG00000000001.1")]
    assert _loci(index, "IL9R") == [
        ("chrX", 5, 8, "ENSG00000000003.1"),
        ("chrY", 9, 12, "ENSG00000000004.1"),
    ]
    assert _loci(index, "DROPME") == []
    assert _loci(index, "EXONONLY") == []
    symbols = {row[0] for row in index["rows"]}
    assert symbols == {"BHLHE40", "DELEC1", "EGFR", "IL9R", "TP53"}
    alias_targets = {
        alias.upper(): [
            index["rows"][row_index][0]
            for name, row_index in index["aliases"]
            if name.upper() == alias.upper()
        ]
        for alias in ("DEC1", "ERBB1", "ERRP", "ERBB", "p53", "NOPE", "EGFR")
    }
    assert alias_targets["DEC1"] == ["BHLHE40", "DELEC1"]
    assert alias_targets["ERBB1"] == ["EGFR"]
    assert alias_targets["P53"] == ["TP53"]
    assert alias_targets["NOPE"] == []
    assert alias_targets["EGFR"] == []
    for alias, row_index in index["aliases"]:
        assert alias.upper() != index["rows"][row_index][0].upper()
    assert "chrM" not in index["chromosomes"]
    assert index["provenance"]["counts"]["genes"] == 6
    assert index["provenance"]["counts"]["dropped_unknown_sequence"] == 1
    assert index["provenance"]["counts"]["dropped_incomplete"] == 2


def test_readme_for_a_different_release_is_rejected(tmp_path: Path) -> None:
    paths = _write_inputs(
        tmp_path,
        readme=README.replace(
            "gencode.v49.primary_assembly.annotation.gtf",
            "gencode.v40.primary_assembly.annotation.gtf",
        ),
    )
    result = _run(paths)
    assert result.returncode == 1
    assert "release 49" in result.stderr
    assert not paths["output"].exists()


def test_gzip_annotation_uses_the_compressed_checksum(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    plain = paths["gtf"].read_bytes()
    gz_path = tmp_path / "annotation.gtf.gz"
    with gzip.open(gz_path, "wb") as handle:
        handle.write(plain)
    paths["gtf"].unlink()
    paths["gtf"] = gz_path
    result = _run(paths)
    assert result.returncode == 0, result.stderr
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    annotation = provenance["annotation"]
    assert annotation["sha256"] == hashlib.sha256(gz_path.read_bytes()).hexdigest()
    assert annotation["uncompressed_sha256"] == hashlib.sha256(plain).hexdigest()
    assert provenance["aliases"]["source"] == "HGNC complete set (local mirror)"
    assert "source_url" not in provenance["aliases"]
    assert provenance["viewer_scaffold"]["sequences"][0] == "chr1"
    index = json.loads(paths["output"].read_text(encoding="utf-8"))
    assert _loci(index, "EGFR") == [("chr7", 10, 40, "ENSG00000000001.1")]


def test_plain_gtf_without_readme_is_rejected(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    result = _run(paths, readme=False)
    assert result.returncode == 1
    assert "README" in result.stderr
    assert not paths["output"].exists()


def test_ncbi_preamble_and_genbank_contig_are_read(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    paths["report"].write_text(
        "# Assembly name: GRCh38\n# " + REPORT + "KI270442.1\tchrUn_KI270442v1\n"
    )
    with paths["gtf"].open("a") as handle:
        handle.write(_gene("KI270442.1", 1, 2, "ENSG00000000010.1", "SCAFF"))
    result = _run(paths)
    assert result.returncode == 0, result.stderr
    index = json.loads(paths["output"].read_bytes())
    assert _loci(index, "SCAFF") == [("chrUn_KI270442v1", 1, 2, "ENSG00000000010.1")]
    assert index["provenance"]["counts"]["outside_viewer_scaffold"] == 1


def test_hgnc_without_alias_columns_is_rejected(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path, hgnc="hgnc_id\tsymbol\tensembl_gene_id\n")
    result = _run(paths)
    assert result.returncode == 1
    assert "missing columns" in result.stderr
