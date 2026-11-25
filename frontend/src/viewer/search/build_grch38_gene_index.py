#!/usr/bin/env python3
"""Build the bundled hg38 gene index from a local GENCODE GTF.

The browser ships the JSON this script writes. It does not ship the GTF.
Rebuilding reads the annotation, the assembly report, the genomepy README,
and an HGNC complete set, then records their checksums in the provenance.

Example:

    python3 frontend/src/viewer/search/build_grch38_gene_index.py \\
        --gtf /path/to/GRCh38.annotation.gtf \\
        --assembly-report /path/to/assembly_report.txt \\
        --readme /path/to/README.txt \\
        --hgnc /path/to/hgnc_complete_set.txt \\
        --output frontend/src/viewer/search/hg38-gene-index.json \\
        --provenance frontend/src/viewer/search/hg38-gene-index.provenance.json
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import sys
from dataclasses import dataclass
from pathlib import Path

RELEASE = "GENCODE human 49"
FORMAT = "helixweave.gene-index.v1"
GENCODE_V49_GTF = "gencode.v49.primary_assembly.annotation.gtf"
ASSEMBLY_MATCH = ["hg38", "GRCh38", "GRCh38.p14", "GCA_000001405.15"]
SCRIPT = "frontend/src/viewer/search/build_grch38_gene_index.py"

_ATTRIBUTE = re.compile(r'([A-Za-z0-9_]+)\s+"([^"]*)"')
_GENE_ID = re.compile(r"^ENSG\d+\.\d+$")
_ALIAS = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._+-]{1,63}$")

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


class _HashedLines:
    def __init__(self, path: Path) -> None:
        self.path = path
        self.digest = hashlib.sha256()
        self.size = 0

    def __iter__(self):
        with self.path.open("rb") as handle:
            for raw in handle:
                self.size += len(raw)
                self.digest.update(raw)
                yield raw


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
            "README annotation url is not GENCODE human release 49 "
            "primary assembly"
        )
    if found.get("original name") != "hg38":
        raise IndexBuildError("README original name is not hg38")
    if not found.get("assembly_accession"):
        raise IndexBuildError("README is missing assembly_accession")
    return found


def ucsc_names(path: Path) -> set[str]:
    with path.open(encoding="utf-8", newline="") as handle:
        header = handle.readline().rstrip("\r\n").split("\t")
        try:
            column = header.index("UCSC-style-name")
        except ValueError as error:
            raise IndexBuildError(
                "assembly report has no UCSC-style-name column"
            ) from error
        names: set[str] = set()
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\r\n").split("\t")
            if len(parts) <= column:
                continue
            name = parts[column].strip()
            if name and name != "na":
                names.add(name)
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
    path: Path, admitted: set[str]
) -> tuple[list[GeneRow], dict[str, int], str, int]:
    lines = _HashedLines(path)
    rows: list[GeneRow] = []
    dropped_incomplete = 0
    dropped_sequence = 0
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
        chrom = fields[0]
        if chrom not in admitted:
            dropped_sequence += 1
            continue
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
            )
        )
    if not rows:
        raise IndexBuildError("no gene features were admitted")
    if lines.size != path.stat().st_size:
        raise IndexBuildError("GTF changed while it was read")
    counts = {
        "dropped_incomplete": dropped_incomplete,
        "dropped_unknown_sequence": dropped_sequence,
    }
    return rows, counts, lines.digest.hexdigest(), lines.size


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
        required = {"ensembl_gene_id", "symbol", "alias_symbol", "prev_symbol"}
        missing = required - fieldnames
        if missing:
            raise IndexBuildError(
                "HGNC file is missing columns: " + ", ".join(sorted(missing))
            )
        for record in reader:
            ensembl = (record.get("ensembl_gene_id") or "").strip()
            base = ensembl.split(".", 1)[0]
            if not base.startswith("ENSG"):
                continue
            found = aliases.setdefault(base, set())
            symbol = (record.get("symbol") or "").strip()
            if symbol:
                found.add(symbol)
            for column in ("alias_symbol", "prev_symbol"):
                found.update(_split_hgnc_list(record.get(column) or ""))
    return aliases


def _alias_links(
    rows: list[GeneRow], hgnc: dict[str, set[str]]
) -> list[list[object]]:
    by_ensembl: dict[str, list[int]] = {}
    for index, row in enumerate(rows):
        by_ensembl.setdefault(row.gene_id.split(".", 1)[0], []).append(index)
    links: dict[tuple[str, int], str] = {}
    for base, spellings in hgnc.items():
        row_indexes = by_ensembl.get(base)
        if not row_indexes:
            continue
        for row_index in row_indexes:
            symbol = rows[row_index].symbol
            for alias in spellings:
                if _ALIAS.fullmatch(alias) is None:
                    continue
                if alias.upper() == symbol.upper():
                    continue
                links.setdefault((alias.upper(), row_index), alias)
    return [
        [links[key], key[1]]
        for key in sorted(links, key=lambda item: (item[0], item[1]))
    ]


def build_index(
    *,
    gtf: Path,
    assembly_report: Path,
    readme: Path,
    hgnc: Path,
) -> dict[str, object]:
    meta = read_readme(readme)
    admitted = ucsc_names(assembly_report)
    rows, drop_counts, gtf_sha256, gtf_size = read_genes(gtf, admitted)
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
    readme_sha, _readme_size = sha256_file(readme)
    report_sha, _report_size = sha256_file(assembly_report)
    hgnc_sha, hgnc_size = sha256_file(hgnc)
    provenance = {
        "annotation": {
            "release": RELEASE,
            "source_url": meta["annotation url"],
            "filename": gtf.name,
            "bytes": gtf_size,
            "sha256": gtf_sha256,
            "feature": "gene",
            "assembly_accession": meta["assembly_accession"],
            "genome_url": meta.get("genome url", ""),
            "genomepy_version": meta.get("genomepy version", ""),
            "readme_date": meta.get("date", ""),
            "readme_filename": readme.name,
            "readme_sha256": readme_sha,
            "attestation": ATTESTATION,
            "license": GENCODE_LICENSE,
        },
        "assembly_report": {
            "filename": assembly_report.name,
            "sha256": report_sha,
            "admitted_sequence_names": len(admitted),
        },
        "aliases": {
            "source": "HGNC complete set",
            "filename": hgnc.name,
            "bytes": hgnc_size,
            "sha256": hgnc_sha,
            "fields": ["symbol", "alias_symbol", "prev_symbol"],
            "join": "Ensembl gene id without version",
            "license": HGNC_LICENSE,
        },
        "generator": SCRIPT,
        "coordinate_system": "1-based-closed",
        "region_format": "chrom:start-end",
        "coordinate_note": COORDINATE_NOTE,
        "assembly_note": (
            "assembly_match lists browser scaffolds that use these primary "
            "chromosome coordinates. GENCODE 49 annotates GRCh38.p14. The "
            "local README identifies UCSC hg38 and its assembly accession. "
            "This primary-assembly GTF has no patch or alt loci. Sequences "
            "are kept only when the assembly report gives them a UCSC-style "
            "name."
        ),
        "counts": {
            "genes": len(rows),
            "aliases": len(alias_links),
            "chromosomes": len(chromosomes),
            "dropped_incomplete": drop_counts["dropped_incomplete"],
            "dropped_unknown_sequence": drop_counts["dropped_unknown_sequence"],
        },
    }
    return {
        "format": FORMAT,
        "assembly": "hg38",
        "assembly_match": list(ASSEMBLY_MATCH),
        "coordinate_system": "1-based-closed",
        "release": RELEASE,
        "chromosomes": chromosomes,
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
    parser.add_argument("--readme", type=Path, required=True)
    parser.add_argument("--hgnc", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
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
