#!/usr/bin/env python3
"""Rebuild the hg38 search index from the three pinned public downloads.

Download each source_url in hg38-gene-index.provenance.json, then run

    python3 frontend/src/viewer/search/build_grch38_gene_index.py \\
        --gtf gencode.v49.primary_assembly.annotation.gtf.gz \\
        --assembly-report GCA_000001405.15_GRCh38_assembly_report.txt \\
        --hgnc hgnc_complete_set_2025-10-07.tsv

To inspect a genomepy mirror instead, pass its GTF, assembly report and
--readme README.txt. The provenance then identifies those local input bytes.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import sys
from collections.abc import Iterable
from dataclasses import dataclass
from pathlib import Path

RELEASE = "GENCODE human 49"
FORMAT = "helixweave.gene-index.v1"
GENCODE_V49_GTF = "gencode.v49.primary_assembly.annotation.gtf"
GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/"
    + GENCODE_V49_GTF
    + ".gz"
)
REPORT_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/"
    "GCA_000001405.15_GRCh38_assembly_report.txt"
)
HGNC_URL = (
    "https://storage.googleapis.com/public-download-files/hgnc/archive/archive/quarterly/tsv/"
    "hgnc_complete_set_2025-10-07.tsv"
)
GTF_SHA256 = "5cabe88e93064851dc6cb826efe0e8b7babf7985e73df280cfb05910d05cd023"
REPORT_SHA256 = "aa733ce92719f6339b3a4540137d4322894434363314c4b7704ba6691e9eab66"
HGNC_SHA256 = "30c2a235934e0f2ec579cb5a7158b8df0091a6674793341d08b8ea9e42e797e5"
# The chromosome list in wuepgg 62.0.2's hg38 navigation context.
VIEWER_CHROMOSOMES = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY", "chrM"}
ASSEMBLY_MATCH = ["hg38", "GRCh38", "GRCh38.p14", "GCA_000001405.15"]
SCRIPT = "frontend/src/viewer/search/build_grch38_gene_index.py"

_ATTRIBUTE = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]*)"')
_GENE_ID = re.compile(r"^ENSG\d+\.\d+$")

GENCODE_LICENSE = (
    "GENCODE project data is open access "
    "(https://www.gencodegenes.org/pages/data_access.html). "
    "Reuse follows the EMBL-EBI Terms of Use: EMBL-EBI expects attribution "
    "and does not add restrictions beyond the data owner's terms. "
    "This index redistributes gene id, symbol, type, strand, and span for "
    "browser search. It is not the GTF."
)
HGNC_LICENSE = (
    "HGNC alias_symbol and prev_symbol are taken from the HGNC complete set, "
    "which HGNC releases under CC0. Attribution is requested: Seal et al., "
    "Nucleic Acids Research (2023); HGNC RRID:SCR_002827."
)
COORDINATE_NOTE = (
    "Spans are GENCODE gene-feature coordinates: 1-based and closed, written "
    "as chrom:start-end. Typed coordinates are passed through unchanged. The "
    "same span in 0-based half-open BED form uses start-1 and the same end."
)
ATTESTATION = (
    "The GTF has no ##version header. The release is the annotation url in "
    "the local genomepy README, and the file identity is the SHA-256 of the "
    "GTF bytes. Sequence names in that GTF are already UCSC-style."
)


class IndexBuildError(Exception):
    """The inputs cannot produce a traceable index."""


@dataclass(frozen=True)
class GeneRow:
    symbol: str
    gene_id: str
    chrom: str
    start: int
    end: int
    strand: str
    gene_type: str
    hgnc_id: str


def sha256_file(path: Path) -> tuple[str, int]:
    digest = hashlib.sha256()
    size = 0
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            size += len(chunk)
            digest.update(chunk)
    return digest.hexdigest(), size


def read_readme(path: Path) -> dict[str, str]:
    found: dict[str, str] = {}
    for line in path.read_text(encoding="utf-8").splitlines():
        if ":" not in line:
            continue
        key, value = line.split(":", 1)
        key = key.strip().lower()
        value = value.strip()
        if key and value:
            found.setdefault(key, value)
    url = found.get("annotation url", "")
    if GENCODE_V49_GTF not in url:
        raise IndexBuildError(
            "README annotation url is not GENCODE human release 49 primary assembly"
        )
    if found.get("original name") != "hg38":
        raise IndexBuildError("README original name is not hg38")
    if not found.get("assembly_accession"):
        raise IndexBuildError("README is missing assembly_accession")
    return found


def ucsc_names(path: Path) -> dict[str, str]:
    with path.open(encoding="utf-8", newline="") as handle:
        for line in handle:
            header = line.removeprefix("# ").rstrip("\r\n").split("\t")
            if "UCSC-style-name" in header:
                column = header.index("UCSC-style-name")
                break
        else:
            raise IndexBuildError("assembly report has no UCSC-style-name column")
        alias_columns = [
            header.index(name)
            for name in ("Sequence-Name", "GenBank-Accn", "RefSeq-Accn")
            if name in header
        ]
        names: dict[str, str] = {}
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) <= column:
                continue
            name = parts[column].strip()
            if name and name != "na":
                names[name] = name
                for index in alias_columns:
                    alias = parts[index].strip()
                    if alias and alias != "na":
                        names[alias] = name
    if not names:
        raise IndexBuildError("assembly report lists no UCSC sequence names")
    return names


def _decode_line(raw: bytes) -> str:
    text = raw.decode("utf-8")
    if text.endswith("\n"):
        text = text[:-1]
    if text.endswith("\r"):
        text = text[:-1]
    return text


def _attributes(field: str) -> dict[str, str]:
    found: dict[str, str] = {}
    for key, value in _ATTRIBUTE.findall(field):
        found.setdefault(key, value)
    return found


def _acceptable_symbol(value: str) -> bool:
    if not value or len(value) > 80 or any(char.isspace() for char in value):
        return False
    return all(32 < ord(char) < 127 for char in value)


def read_genes(
    lines: Iterable[bytes], admitted: dict[str, str]
) -> tuple[list[GeneRow], dict[str, int]]:
    rows: list[GeneRow] = []
    dropped_incomplete = 0
    dropped_sequence = 0
    outside_viewer = 0
    for raw in lines:
        if b"\tgene\t" not in raw:
            continue
        fields = _decode_line(raw).split("\t")
        if len(fields) < 9 or fields[2] != "gene":
            continue
        try:
            start = int(fields[3])
            end = int(fields[4])
        except ValueError:
            dropped_incomplete += 1
            continue
        attrs = _attributes(fields[8])
        gene_id = attrs.get("gene_id", "")
        symbol = attrs.get("gene_name", "")
        if (
            start < 1
            or end < start
            or _GENE_ID.fullmatch(gene_id) is None
            or not _acceptable_symbol(symbol)
        ):
            dropped_incomplete += 1
            continue
        chrom = admitted.get(fields[0])
        if chrom is None:
            dropped_sequence += 1
            continue
        if chrom not in VIEWER_CHROMOSOMES:
            outside_viewer += 1
        strand = fields[6] if fields[6] in {"+", "-", "."} else "."
        rows.append(
            GeneRow(
                symbol=symbol,
                gene_id=gene_id,
                chrom=chrom,
                start=start,
                end=end,
                strand=strand,
                gene_type=attrs.get("gene_type") or "unknown",
                hgnc_id=attrs.get("hgnc_id", ""),
            )
        )
    if not rows:
        raise IndexBuildError("no gene features were admitted")
    counts = {
        "dropped_incomplete": dropped_incomplete,
        "dropped_unknown_sequence": dropped_sequence,
        "outside_viewer_scaffold": outside_viewer,
    }
    return rows, counts


def _chrom_key(name: str) -> tuple[int, int, str]:
    body = name[3:] if name.startswith("chr") else name
    if body.isdigit():
        return (0, int(body), "")
    special = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    if body in special:
        return (1, special[body], "")
    return (2, 0, name)


def _split_hgnc_list(value: str) -> list[str]:
    text = value.strip()
    if not text:
        return []
    return [part.strip() for part in text.split("|") if part.strip()]


def read_hgnc_aliases(path: Path) -> dict[str, set[str]]:
    aliases: dict[str, set[str]] = {}
    with path.open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = set(reader.fieldnames or [])
        required = {
            "hgnc_id",
            "ensembl_gene_id",
            "symbol",
            "alias_symbol",
            "prev_symbol",
        }
        missing = required - fieldnames
        if missing:
            raise IndexBuildError(
                "HGNC file is missing columns: " + ", ".join(sorted(missing))
            )
        for record in reader:
            hgnc_id = (record.get("hgnc_id") or "").strip()
            if not hgnc_id:
                continue
            spellings = {(record.get("symbol") or "").strip()}
            for column in ("alias_symbol", "prev_symbol"):
                spellings.update(_split_hgnc_list(record.get(column) or ""))
            aliases.setdefault(hgnc_id, set()).update(spellings)
            ensembl = (record.get("ensembl_gene_id") or "").strip().split(".", 1)[0]
            if ensembl:
                aliases.setdefault(ensembl, set()).update(spellings)
    return aliases


def _alias_links(rows: list[GeneRow], hgnc: dict[str, set[str]]) -> list[list[object]]:
    links: dict[tuple[str, int], str] = {}
    for index, row in enumerate(rows):
        # A pseudoautosomal gene can have one HGNC ID and distinct X/Y Ensembl IDs.
        join_id = row.hgnc_id or row.gene_id.split(".", 1)[0]
        for alias in sorted(hgnc.get(join_id, ())):
            if not alias or len(alias) > 80 or not alias.isprintable():
                continue
            if alias.upper() != row.symbol.upper():
                links.setdefault((alias.upper(), index), alias)
    return [
        [links[key], key[1]]
        for key in sorted(links, key=lambda item: (item[0], item[1]))
    ]


def build_index(
    *,
    gtf: Path,
    assembly_report: Path,
    readme: Path | None,
    hgnc: Path,
) -> dict[str, object]:
    meta = read_readme(readme) if readme else None
    gtf_sha256, gtf_size = sha256_file(gtf)
    report_sha, _ = sha256_file(assembly_report)
    hgnc_sha, hgnc_size = sha256_file(hgnc)
    if meta is None:
        for path, digest, expected in (
            (gtf, gtf_sha256, GTF_SHA256),
            (assembly_report, report_sha, REPORT_SHA256),
            (hgnc, hgnc_sha, HGNC_SHA256),
        ):
            if digest != expected:
                raise IndexBuildError(
                    f"{path.name}: SHA-256 differs from pinned source"
                )
    admitted = ucsc_names(assembly_report)
    with gtf.open("rb") as raw:
        compressed = raw.read(2) == b"\x1f\x8b"
        raw.seek(0)
        if compressed:
            with gzip.GzipFile(fileobj=raw) as lines:
                rows, drop_counts = read_genes(lines, admitted)
        else:
            rows, drop_counts = read_genes(raw, admitted)
    rows.sort(
        key=lambda row: (
            row.symbol.upper(),
            _chrom_key(row.chrom),
            row.start,
            row.end,
            row.gene_id,
        )
    )
    chromosomes = sorted({row.chrom for row in rows}, key=_chrom_key)
    chrom_index = {name: index for index, name in enumerate(chromosomes)}
    alias_links = _alias_links(rows, read_hgnc_aliases(hgnc))
    annotation = {
        "release": RELEASE,
        "source_url": GTF_URL,
        "filename": GTF_URL.rsplit("/", 1)[-1],
        "bytes": gtf_size,
        "sha256": gtf_sha256,
        "feature": "gene",
        "assembly_accession": "GCA_000001405.29",
        "license": GENCODE_LICENSE,
    }
    if readme and meta:
        annotation.update(
            {
                "source_url": meta["annotation url"],
                "filename": gtf.name,
                "genome_url": meta.get("genome url", ""),
                "genomepy_version": meta.get("genomepy version", ""),
                "readme_date": meta.get("date", ""),
                "readme_filename": readme.name,
                "readme_sha256": sha256_file(readme)[0],
                "attestation": ATTESTATION,
            }
        )
    provenance = {
        "annotation": annotation,
        "assembly_report": {
            **({"source_url": REPORT_URL} if report_sha == REPORT_SHA256 else {}),
            "filename": REPORT_URL.rsplit("/", 1)[-1]
            if not readme
            else assembly_report.name,
            "assembly_accession": "GCA_000001405.15",
            "sha256": report_sha,
            "admitted_sequence_names": len(set(admitted.values())),
        },
        "aliases": {
            "source": "HGNC complete set",
            **({"source_url": HGNC_URL} if hgnc_sha == HGNC_SHA256 else {}),
            "filename": HGNC_URL.rsplit("/", 1)[-1] if not readme else hgnc.name,
            "bytes": hgnc_size,
            "sha256": hgnc_sha,
            "fields": ["symbol", "alias_symbol", "prev_symbol"],
            "join": "GENCODE hgnc_id, with versionless Ensembl ID for rows without an HGNC ID",
            "license": HGNC_LICENSE,
        },
        "generator": SCRIPT,
        "coordinate_system": "1-based-closed",
        "region_format": "chrom:start-end",
        "coordinate_note": COORDINATE_NOTE,
        "assembly_note": (
            "GENCODE 49 annotates GRCh38.p14. The hg38 navigation context in "
            "wuepgg 62.0.2 has 25 chromosomes: chr1–chr22, chrX, chrY and chrM. "
            "Only gene features on these UCSC names from the GRCh38 assembly "
            "report are selectable."
        ),
        "counts": {
            "genes": len(rows),
            "aliases": len(alias_links),
            "chromosomes": len(chromosomes),
            **drop_counts,
        },
    }
    return {
        "format": FORMAT,
        "assembly": "hg38",
        "assembly_match": list(ASSEMBLY_MATCH),
        "coordinate_system": "1-based-closed",
        "release": RELEASE,
        "chromosomes": chromosomes,
        "viewer_chromosomes": sorted(VIEWER_CHROMOSOMES, key=_chrom_key),
        "rows": [
            [
                row.symbol,
                row.gene_id,
                chrom_index[row.chrom],
                row.start,
                row.end,
                row.strand,
                row.gene_type,
            ]
            for row in rows
        ],
        "aliases": alias_links,
        "provenance": provenance,
    }


def _publish(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(text, encoding="utf-8")
    temporary.replace(path)


def write_index(payload: dict[str, object], output: Path, provenance: Path) -> None:
    encoded = json.dumps(payload, ensure_ascii=True, separators=(",", ":"))
    loaded = json.loads(encoded)
    counts = payload["provenance"]["counts"]
    if len(loaded["rows"]) != counts["genes"]:
        raise IndexBuildError("wrote a gene index that does not reload")
    _publish(output, encoded + "\n")
    pretty = json.dumps(payload["provenance"], ensure_ascii=True, indent=2)
    _publish(provenance, pretty + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--assembly-report", type=Path, required=True)
    parser.add_argument(
        "--readme", type=Path, help="genomepy README for a local mirror"
    )
    parser.add_argument("--hgnc", type=Path, required=True)
    parser.add_argument(
        "--output", type=Path, default=Path(__file__).with_name("hg38-gene-index.json")
    )
    parser.add_argument(
        "--provenance",
        type=Path,
        default=Path(__file__).with_name("hg38-gene-index.provenance.json"),
    )
    args = parser.parse_args(argv)
    try:
        payload = build_index(
            gtf=args.gtf,
            assembly_report=args.assembly_report,
            readme=args.readme,
            hgnc=args.hgnc,
        )
        write_index(payload, args.output, args.provenance)
    except (OSError, UnicodeError, csv.Error, IndexBuildError) as error:
        print(f"error: {error}", file=sys.stderr)
        return 1
    counts = payload["provenance"]["counts"]
    print(
        f"genes {counts['genes']} aliases {counts['aliases']} "
        f"dropped_sequence {counts['dropped_unknown_sequence']} "
        f"dropped_incomplete {counts['dropped_incomplete']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
