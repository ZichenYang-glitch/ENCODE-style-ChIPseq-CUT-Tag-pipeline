"""The mm10 index keeps GTF spans and mouse MGI aliases, not HGNC."""

from __future__ import annotations

import gzip
import hashlib
import json
import subprocess
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))

from build_grcm38_gene_index import PINNED_MM10_SEQUENCES

SCRIPT = Path(__file__).with_name("build_grcm38_gene_index.py")
README = """\
name: GRCm38
original name: mm10
assembly_accession: GCA_000001635.2
genome url: https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/chromFa.tar.gz
annotation url: ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
genomepy version: 0.16.3
date: 2025-12-15 19:29:42
"""
EXTERNAL_DB = (
    "1400\tMGI\t\\N\tKNOWNXREF\t270\tMGI Symbol\tPRIMARY_DB_SYNONYM\n"
    "1100\tHGNC\t\\N\tKNOWNXREF\t1\tHGNC Symbol\tPRIMARY_DB_SYNONYM\n"
)
GENE = (
    "1321000\tprotein_coding\t1\t1\t1\t2\t1\t\\N\tensembl\tEgfr\t1\t\\N\t"
    "ENSMUSG00000020122\t16\t2020-03-24 00:00:00\t2020-03-24 00:00:00\n"
    "1335377\tprotein_coding\t1\t1\t1\t2\t1\t\\N\tensembl\tTrp53\t1\t\\N\t"
    "ENSMUSG00000059552\t13\t2020-03-24 00:00:00\t2020-03-24 00:00:00\n"
    "1400001\tprotein_coding\t1\t1\t1\t2\t1\t\\N\tensembl\tOther\t1\t\\N\t"
    "ENSMUSG00000000003\t1\t2020-03-24 00:00:00\t2020-03-24 00:00:00\n"
    "1400002\tprotein_coding\t1\t1\t1\t2\t1\t\\N\tensembl\tKeepScaffold\t1\t"
    "\\N\tENSMUSG00000000004\t1\t2020-03-24 00:00:00\t2020-03-24 00:00:00\n"
)
XREF = (
    "2192528\t1400\tMGI:95294\tEgfr\t0\t\\N\tDIRECT\t\n"
    "2194877\t1400\tMGI:98834\tTrp53\t0\t\\N\tDIRECT\t\n"
    "2195001\t1400\tMGI:1\tOther\t0\t\\N\tDIRECT\t\n"
    "2195002\t1400\tMGI:2\tKeepScaffold\t0\t\\N\tDIRECT\t\n"
    "999\t1100\tHGNC:3236\tEGFR\t0\t\\N\tDIRECT\t\n"
)
SYNONYM = (
    "2192528\tErbb\n"
    "2192528\tWa5\n"
    "2192528\tEgfr\n"
    "2192528\tavian erythroblastic leukemia viral (v-erb-b) oncogene homolog\n"
    "2194877\tp53\n"
    "2195001\tWa5\n"
    "999\tERBB1\n"
)
OBJECT_XREF = (
    "1\t1321000\tGene\t2192528\t\\N\t1\n"
    "2\t1321000\tGene\t999\t\\N\t1\n"
    "3\t1335377\tGene\t2194877\t\\N\t1\n"
    "4\t1400001\tGene\t2195001\t\\N\t1\n"
    "5\t1400002\tGene\t2195002\t\\N\t1\n"
)


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
            _gene("chr11", 10, 40, "ENSMUSG00000020122.16", "Egfr"),
            _gene("chr11", 10, 20, "ENSMUSG00000000011.1", "EXONONLY", feature="exon"),
            _gene("chr11", 100, 110, "ENSMUSG00000059552.13", "Trp53"),
            _gene("chr11", 200, 210, "ENSMUSG00000000003.1", "Other"),
            _gene("chr1_GL456210_random", 1, 5, "ENSMUSG00000000004.1", "KeepScaffold"),
            _gene("chrUn_GL456239", 1, 5, "ENSMUSG00000000005.1", "DropMe"),
            _gene("chr11", 1, 5, "ENSG00000000001.1", "HumanLike"),
            _gene("chr11", 1, 2, "ENSMUSG00000000008.1", "NoName", include_name=False),
            _gene("chr11", 10, 4, "ENSMUSG00000000009.1", "BadEnd"),
        ]
    )


def _report() -> str:
    rows = ["Sequence-Name\tSequence-Length\tUCSC-style-name"]
    for name, length in PINNED_MM10_SEQUENCES:
        rows.append(f"{name}\t{length}\t{name}")
    rows.append("GL456210.1\t169725\tchr1_GL456210_random")
    return "\n".join(rows) + "\n"


def _write_inputs(root: Path, *, readme: str = README, report: str | None = None):
    root.mkdir(parents=True, exist_ok=True)
    gtf = root / "annotation.gtf"
    gtf.write_text(_fixture_gtf(), encoding="utf-8")
    report_path = root / "assembly_report.txt"
    report_path.write_text(report if report is not None else _report(), encoding="utf-8")
    readme_path = root / "README.txt"
    readme_path.write_text(readme, encoding="utf-8")
    tables = {
        "external_db": EXTERNAL_DB,
        "gene": GENE,
        "xref": XREF,
        "synonym": SYNONYM,
        "object_xref": OBJECT_XREF,
    }
    paths = {
        "gtf": gtf,
        "report": report_path,
        "readme": readme_path,
        "output": root / "mm10-gene-index.json",
        "provenance": root / "mm10-gene-index.provenance.json",
    }
    for name, text in tables.items():
        path = root / f"{name}.txt"
        path.write_text(text, encoding="utf-8")
        paths[name] = path
    return paths


def _run(paths: dict[str, Path], *, readme: bool = True):
    command = [
        sys.executable,
        str(SCRIPT),
        "--gtf",
        str(paths["gtf"]),
        "--assembly-report",
        str(paths["report"]),
        "--ensembl-external-db",
        str(paths["external_db"]),
        "--ensembl-gene",
        str(paths["gene"]),
        "--ensembl-xref",
        str(paths["xref"]),
        "--ensembl-synonym",
        str(paths["synonym"]),
        "--ensembl-object-xref",
        str(paths["object_xref"]),
        "--output",
        str(paths["output"]),
        "--provenance",
        str(paths["provenance"]),
    ]
    if readme:
        command.extend(["--readme", str(paths["readme"])])
    return subprocess.run(command, check=False, capture_output=True, text=True)


def _loci(index: dict, symbol: str) -> list[tuple]:
    chromosomes = index["chromosomes"]
    return [
        (chromosomes[row[2]], row[3], row[4], row[1])
        for row in index["rows"]
        if row[0] == symbol
    ]


def _alias_targets(index: dict, alias: str) -> list[str]:
    return [
        index["rows"][row_index][0]
        for name, row_index in index["aliases"]
        if name.upper() == alias.upper()
    ]


def test_index_uses_gtf_spans_and_mgi_aliases(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    result = _run(paths)
    assert result.returncode == 0, result.stderr
    index = json.loads(paths["output"].read_text(encoding="utf-8"))
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    assert index["provenance"] == provenance
    assert index["format"] == "helixweave.gene-index.v1"
    assert index["assembly"] == "mm10"
    assert index["coordinate_system"] == "1-based-closed"
    assert "mm39" not in index["assembly_match"]
    assert "GCA_000001635.9" not in index["assembly_match"]
    assert "GRCm38.p6" in index["assembly_match"]
    assert (
        index["provenance"]["annotation"]["sha256"]
        == hashlib.sha256(paths["gtf"].read_bytes()).hexdigest()
    )
    assert _loci(index, "Egfr") == [("chr11", 10, 40, "ENSMUSG00000020122.16")]
    assert _loci(index, "Trp53") == [("chr11", 100, 110, "ENSMUSG00000059552.13")]
    assert _loci(index, "KeepScaffold") == [
        ("chr1_GL456210_random", 1, 5, "ENSMUSG00000000004.1")
    ]
    assert _loci(index, "DropMe") == []
    assert _loci(index, "HumanLike") == []
    assert _loci(index, "EXONONLY") == []
    symbols = {row[0] for row in index["rows"]}
    assert symbols == {"Egfr", "KeepScaffold", "Other", "Trp53"}
    assert _alias_targets(index, "Erbb") == ["Egfr"]
    assert _alias_targets(index, "Wa5") == ["Egfr", "Other"]
    assert _alias_targets(index, "p53") == ["Trp53"]
    assert _alias_targets(index, "ERBB1") == []
    assert _alias_targets(index, "Egfr") == []
    assert _alias_targets(index, "avian erythroblastic leukemia viral (v-erb-b) oncogene homolog") == []
    for alias, row_index in index["aliases"]:
        assert alias.upper() != index["rows"][row_index][0].upper()
        assert " " not in alias
    assert index["provenance"]["counts"]["genes"] == 4
    assert index["provenance"]["counts"]["dropped_unknown_sequence"] == 1
    assert index["provenance"]["counts"]["dropped_incomplete"] == 3
    assert index["provenance"]["aliases"]["db_name"] == "MGI"
    assert index["provenance"]["aliases"]["excluded_db_names"] == ["HGNC"]
    assert index["provenance"]["aliases"]["ignored_hgnc_gene_links"] == 1
    assert index["provenance"]["viewer_scaffold"]["sequences"][-1] == "chrM"
    assert index["provenance"]["viewer_scaffold"]["lengths"][0] == 195471971
    assert "upstream_crosscheck" not in provenance
    assert "HGNC" not in index["provenance"]["aliases"]["license"] or (
        "not used" in index["provenance"]["aliases"]["license"]
    )


def test_readme_for_a_different_release_is_rejected(tmp_path: Path) -> None:
    paths = _write_inputs(
        tmp_path,
        readme=README.replace(
            "gencode.vM25.primary_assembly.annotation.gtf",
            "gencode.vM26.primary_assembly.annotation.gtf",
        ),
    )
    result = _run(paths)
    assert result.returncode == 1
    assert "M25" in result.stderr
    assert not paths["output"].exists()


def test_mm39_name_and_accession_are_rejected(tmp_path: Path) -> None:
    named = _write_inputs(
        tmp_path / "name",
        readme=README.replace("original name: mm10", "original name: mm39"),
    )
    result = _run(named)
    assert result.returncode == 1
    assert "mm10" in result.stderr
    assert not named["output"].exists()
    accession = _write_inputs(
        tmp_path / "accession",
        readme=README.replace("GCA_000001635.2", "GCA_000001635.9"),
    )
    result = _run(accession)
    assert result.returncode == 1
    assert "GRCm39" in result.stderr
    assert not accession["output"].exists()


def test_grcm38_patch_accession_is_accepted(tmp_path: Path) -> None:
    paths = _write_inputs(
        tmp_path,
        readme=README.replace("GCA_000001635.2", "GCA_000001635.8"),
    )
    result = _run(paths)
    assert result.returncode == 0, result.stderr
    provenance = json.loads(paths["provenance"].read_text(encoding="utf-8"))
    assert provenance["annotation"]["assembly_accession"] == "GCA_000001635.8"


def test_pinned_mm10_length_is_required(tmp_path: Path) -> None:
    report = _report().replace("chr1\t195471971\tchr1", "chr1\t195154279\tchr1")
    paths = _write_inputs(tmp_path, report=report)
    result = _run(paths)
    assert result.returncode == 1
    assert "pinned mm10 length" in result.stderr
    assert "195471971" in result.stderr
    assert not paths["output"].exists()


def test_gzip_annotation_records_both_checksums(tmp_path: Path) -> None:
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
    assert "release-100" in provenance["aliases"]["source_url"]
    assert provenance["aliases"]["source_url"].endswith("external_synonym.txt.gz")
    assert provenance["viewer_scaffold"]["sequences"][0] == "chr1"
    index = json.loads(paths["output"].read_text(encoding="utf-8"))
    assert _loci(index, "Egfr") == [("chr11", 10, 40, "ENSMUSG00000020122.16")]


def test_plain_gtf_without_readme_is_rejected(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    result = _run(paths, readme=False)
    assert result.returncode == 1
    assert "README" in result.stderr
    assert not paths["output"].exists()


def test_shipped_mm10_index_matches_gencode_m25() -> None:
    directory = Path(__file__).parent
    index_path = directory / "mm10-gene-index.json"
    provenance_path = directory / "mm10-gene-index.provenance.json"
    index = json.loads(index_path.read_text(encoding="utf-8"))
    provenance = json.loads(provenance_path.read_text(encoding="utf-8"))
    assert index["provenance"] == provenance
    assert hashlib.sha256(index_path.read_bytes()).hexdigest() == (
        "7514b737170d40b906d9ea3572062a10e29247735380e449d60b352cb2bdbe9c"
    )
    assert index["assembly"] == "mm10"
    assert index["format"] == "helixweave.gene-index.v1"
    assert index["coordinate_system"] == "1-based-closed"
    assert index["chromosomes"][:22] == [name for name, _length in PINNED_MM10_SEQUENCES]
    assert _loci(index, "Egfr") == [
        ("chr11", 16752203, 16918158, "ENSMUSG00000020122.16")
    ]
    assert _loci(index, "Trp53") == [
        ("chr11", 69580359, 69591873, "ENSMUSG00000059552.13")
    ]
    assert _alias_targets(index, "Erbb") == ["Egfr"]
    assert _alias_targets(index, "Wa5") == ["Egfr"]
    assert _alias_targets(index, "p53") == ["Trp53"]
    assert _alias_targets(index, "ERBB1") == []
    egfr = next(row for row in index["rows"] if row[0] == "Egfr")
    trp53 = next(row for row in index["rows"] if row[0] == "Trp53")
    assert egfr[5:] == ["+", "protein_coding"]
    assert trp53[5:] == ["+", "protein_coding"]
    assert provenance["counts"]["genes"] == 55487
    assert provenance["counts"]["aliases"] == 45612
    assert provenance["counts"]["chromosomes"] == 45
    assert provenance["annotation"]["sha256"] == (
        "dbd5e0488a5994cc5e9379e13ecaceb0d2affcb0c6c9f7584def7f1be01769b6"
    )
    assert provenance["annotation"]["bytes"] == 888182176
    cross = provenance["upstream_crosscheck"]
    assert cross["official_gzip_md5"] == "c5125258a0a2c5250ddb4c192abbf4e8"
    assert cross["official_gzip_sha256"] == (
        "b3ceeef170048fc86c5dead33ba7af61d71e50c9a19ca5b4b9784b830802f5cd"
    )
    assert provenance["aliases"]["db_name"] == "MGI"
    assert provenance["aliases"]["ignored_hgnc_gene_links"] == 2
    assert provenance["aliases"]["sha256"] == (
        "ad3d47113baeca83c1d83166e360cdb3d58c13d48a52ad5f45ee350eb44dca70"
    )
    assert all(
        item["sha256"] == item["canonical_sha256"] for item in provenance["aliases"]["files"]
    )
    pinned = set(provenance["viewer_scaffold"]["sequences"])
    off_scaffold = sum(
        1 for row in index["rows"] if index["chromosomes"][row[2]] not in pinned
    )
    assert off_scaffold == 86
    assert provenance["viewer_scaffold"]["lengths"] == [
        length for _name, length in PINNED_MM10_SEQUENCES
    ]


def test_alias_table_without_mgi_is_rejected(tmp_path: Path) -> None:
    paths = _write_inputs(tmp_path)
    paths["external_db"].write_text(
        "1100\tHGNC\t\\N\tKNOWNXREF\t1\tHGNC Symbol\tPRIMARY_DB_SYNONYM\n",
        encoding="utf-8",
    )
    result = _run(paths)
    assert result.returncode == 1
    assert "no MGI" in result.stderr
    assert not paths["output"].exists()
