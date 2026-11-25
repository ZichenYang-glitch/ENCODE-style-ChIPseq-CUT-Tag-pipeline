#!/usr/bin/env python3
"""Build the bundled mm10 gene index from GENCODE M25 and Ensembl 100 MGI.

The browser ships only the JSON this script writes. It does not download
the annotation or the alias tables at runtime. Aliases come from the mouse
MGI cross-references in Ensembl 100. HGNC rows in that database are ignored.

Pinned annotation:

* GENCODE mouse release M25 primary-assembly GTF gzip
  https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.primary_assembly.annotation.gtf.gz
  MD5 c5125258a0a2c5250ddb4c192abbf4e8 (the release MD5SUMS line)
  SHA-256 b3ceeef170048fc86c5dead33ba7af61d71e50c9a19ca5b4b9784b830802f5cd
  28556861 bytes. When ``--gtf`` is this gzip, ``--readme`` is unnecessary.

* The shipped index was built from the genomepy plain GTF
  ``GRCm38.annotation.gtf``.
  SHA-256 dbd5e0488a5994cc5e9379e13ecaceb0d2affcb0c6c9f7584def7f1be01769b6
  888182176 bytes. After the official gzip's five header lines, the first
  1872052 records match that plain file. The other 1019 records differ only
  because genomepy renamed scaffold accessions to UCSC-style names.

* Ensembl release 100 ``mus_musculus_core_100_38`` tables
  external_db, gene, xref, external_synonym, and object_xref.
  GENCODE M25 is that Ensembl release on GRCm38. Only ``db_name`` MGI
  (external_db_id 1400) supplies symbols and synonyms. Supply the release
  gzip files or their exact uncompressed contents.

* An assembly report for GCA_000001635.2 with UCSC-style-name and
  Sequence-Length. The GenBank-Accn column maps upstream scaffold names to
  the UCSC names in the local genomepy annotation. Both forms produce the
  same rows and aliases. Provenance records the actual inputs. The report is at
  https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/635/GCA_000001635.2_GRCm38/GCA_000001635.2_GRCm38_assembly_report.txt
  Every sequence in the pinned wuepgg 62.0.2 mm10 genome must be present
  at that genome's length.

Example:

    python3 frontend/src/viewer/search/build_grcm38_gene_index.py \\
        --gtf GRCm38.annotation.gtf \\
        --assembly-report assembly_report.txt \\
        --readme README.txt \\
        --ensembl-external-db external_db.txt.gz \\
        --ensembl-gene gene.txt.gz \\
        --ensembl-xref xref.txt.gz \\
        --ensembl-synonym external_synonym.txt.gz \\
        --ensembl-object-xref object_xref.txt.gz \\
        --output frontend/src/viewer/search/mm10-gene-index.json \\
        --provenance frontend/src/viewer/search/mm10-gene-index.provenance.json
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import re
import sys
from pathlib import Path

import build_grch38_gene_index as hg

RELEASE = "GENCODE mouse M25"
SCRIPT = "frontend/src/viewer/search/build_grcm38_gene_index.py"
ASSEMBLY_MATCH = [
    "mm10",
    "GRCm38",
    "GRCm38.p6",
    "GCA_000001635.2",
    "GCA_000001635.8",
]
# GCA_000001635.2 is GRCm38. GCA_000001635.8 is GRCm38.p6.
# GCA_000001635.9 is GRCm39 and must not label these coordinates.
GRCM38_PATCH_ACCESSIONS = {
    "GCA_000001635.2",
    "GCA_000001635.3",
    "GCA_000001635.4",
    "GCA_000001635.5",
    "GCA_000001635.6",
    "GCA_000001635.7",
    "GCA_000001635.8",
}
GRCM39_ACCESSION = "GCA_000001635.9"
OFFICIAL_GTF_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/"
    "gencode.vM25.primary_assembly.annotation.gtf.gz"
)
OFFICIAL_GTF_SHA = "b3ceeef170048fc86c5dead33ba7af61d71e50c9a19ca5b4b9784b830802f5cd"
OFFICIAL_GTF_BYTES = 28556861
OFFICIAL_GTF_MD5 = "c5125258a0a2c5250ddb4c192abbf4e8"
OFFICIAL_PLAIN_SHA = "009fe0274de483412832839740e92c959877f511833a2c6bcc41043477d1af9e"
OFFICIAL_PLAIN_BYTES = 888172458
KNOWN_LOCAL_SHA = "dbd5e0488a5994cc5e9379e13ecaceb0d2affcb0c6c9f7584def7f1be01769b6"
KNOWN_LOCAL_BYTES = 888182176
ENSEMBL_BASE = "https://ftp.ensembl.org/pub/release-100/mysql/mus_musculus_core_100_38/"
# SHA-256 and size of each release gzip, followed by SHA-256 of its contents.
ENSEMBL_CANONICAL = {
    "external_db.txt.gz": (
        "1627403c7f8ed443ad65087ca11d2e1133499874c322af2046f4e59448ce91d0",
        14232,
        "87f570b4e8eb8a298a876748487a87aa39ea81552914f01d491c679a428c48b3",
    ),
    "gene.txt.gz": (
        "7eb13f94fa7ce33b41d7b2d36ec973aa89e1b9741c6948f311bed14f201f4431",
        3146731,
        "15e6dcf7632553402e7c1dffb0ce8e0136b546798df72429d9a191c3e50d1382",
    ),
    "xref.txt.gz": (
        "6de286b59baaa4c05486284e5fdfec368f289691bd06916f34e22d5964946439",
        30340560,
        "cbed03b891c17c281c4aedc92ca843d4031fdc98be36974ad8f1d45dd070799c",
    ),
    "external_synonym.txt.gz": (
        "ad3d47113baeca83c1d83166e360cdb3d58c13d48a52ad5f45ee350eb44dca70",
        381823,
        "6354a830dca7acdbb23412844a3a5ff7d5870c7bc0bed8efe8dc9ac49c0b9428",
    ),
    "object_xref.txt.gz": (
        "44ee975e37c04f4cea75d69c2ae32d33f8de3f70b71b0b8babac49165e54ca87",
        27426363,
        "bebaecf0d79c00576edb60e3168f3b3db4f310cd8190ed230bc21773e610e4ea",
    ),
}
# Sequences constructed for genome "mm10" in the locked viewer 62.0.2.
# chrM is last. mm39 in that same viewer puts chrM before chrX and its
# chr1 length is 195154279. No random, unplaced, or alt names.
PINNED_MM10_SEQUENCES = (
    ("chr1", 195471971),
    ("chr2", 182113224),
    ("chr3", 160039680),
    ("chr4", 156508116),
    ("chr5", 151834684),
    ("chr6", 149736546),
    ("chr7", 145441459),
    ("chr8", 129401213),
    ("chr9", 124595110),
    ("chr10", 130694993),
    ("chr11", 122082543),
    ("chr12", 120129022),
    ("chr13", 120421639),
    ("chr14", 124902244),
    ("chr15", 104043685),
    ("chr16", 98207768),
    ("chr17", 94987271),
    ("chr18", 90702639),
    ("chr19", 61431566),
    ("chrX", 171031299),
    ("chrY", 91744698),
    ("chrM", 16299),
)
_ENSMUSG = re.compile(r"^ENSMUSG\d+\.\d+$")
MGI_LICENSE = (
    "Mouse aliases are the MGI symbol and external synonyms joined through "
    "Ensembl release 100 mus_musculus_core_100_38, external_db db_name MGI. "
    "HGNC rows in that database are not used. Ensembl redistributes the "
    "cross-reference under its public data policy. MGI asks for attribution: "
    "Mouse Genome Informatics, The Jackson Laboratory."
)
ATTESTATION = (
    "The uncompressed GTF matches a pinned GENCODE M25 annotation digest. "
    "The supplied genomepy README identifies the release and assembly."
)
OFFICIAL_ATTESTATION = (
    "The annotation file is the GENCODE mouse release M25 primary-assembly "
    "GTF gzip. Its MD5 matches the release MD5SUMS, and its SHA-256 is the "
    "checksum of that gzip. uncompressed_sha256 is the checksum of the "
    "decompressed bytes. Reference chromosome names are UCSC-style. Scaffold "
    "names in this gzip are GenBank accessions, so a scaffold is kept only "
    "when the assembly report lists that accession. A README is not required "
    "when the gzip checksum matches the pinned release."
)
LOCAL_CROSSCHECK_ATTESTATION = (
    "This plain GTF is the local genomepy copy of GENCODE mouse release M25. "
    "It has no ##description header. Its primary-chromosome records are "
    "byte-identical to gencode.vM25.primary_assembly.annotation.gtf.gz after "
    "that gzip's five-line header. The remaining records differ only by "
    "UCSC-style scaffold names; coordinates are unchanged. "
    "upstream_crosscheck records the official checksums and the comparison."
)
UPSTREAM_CROSSCHECK = {
    "official_gzip_url": OFFICIAL_GTF_URL,
    "official_gzip_md5": OFFICIAL_GTF_MD5,
    "official_gzip_sha256": OFFICIAL_GTF_SHA,
    "official_gzip_bytes": OFFICIAL_GTF_BYTES,
    "official_uncompressed_sha256": OFFICIAL_PLAIN_SHA,
    "official_uncompressed_bytes": OFFICIAL_PLAIN_BYTES,
    "header_lines": 5,
    "identical_records": 1872052,
    "scaffold_seqname_only_records": 1019,
    "note": (
        "Compared with the pinned gzip after removing its five header lines. "
        "The first 1872052 records match byte for byte. The remaining 1019 "
        "records differ only in the sequence name: genomepy rewrote GenBank "
        "scaffold accessions to UCSC-style names. Coordinates are unchanged."
    ),
}
ASSEMBLY_NOTE = (
    "assembly_match lists names that use these GRCm38 primary chromosome "
    "coordinates. GENCODE mouse M25 corresponds to Ensembl 100 on GRCm38. "
    "GCA_000001635.2 is GRCm38 and GCA_000001635.8 is GRCm38.p6. "
    "GCA_000001635.9 is GRCm39 and is rejected. The pinned viewer genome is "
    "the 22 sequences in viewer_scaffold, from wuepgg 62.0.2, with chrM last. "
    "Sequences are kept when the assembly report gives them a UCSC-style "
    "name; navigation still refuses a name the viewer genome lacks."
)


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
    if "gencode.vM25.primary_assembly.annotation.gtf" not in url:
        raise hg.IndexBuildError(
            "README annotation url is not GENCODE mouse release M25 primary assembly"
        )
    if found.get("original name") != "mm10":
        raise hg.IndexBuildError("README original name is not mm10")
    accession = found.get("assembly_accession", "")
    if not accession:
        raise hg.IndexBuildError("README is missing assembly_accession")
    if accession == GRCM39_ACCESSION or "GRCm39" in accession:
        raise hg.IndexBuildError(f"README assembly {accession} is GRCm39, not GRCm38")
    if accession not in GRCM38_PATCH_ACCESSIONS:
        raise hg.IndexBuildError(
            f"README assembly_accession {accession} is not a GRCm38 release"
        )
    return found


def read_assembly_report(path: Path) -> tuple[dict[str, int], dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        # NCBI includes a commented preamble and header. genomepy strips both
        # the preamble and the header's comment marker.
        headers = (line.lstrip("# ").rstrip("\r\n").split("\t") for line in handle)
        header = next((fields for fields in headers if "UCSC-style-name" in fields), [])
        if "Sequence-Length" not in header:
            raise hg.IndexBuildError(
                "assembly report needs UCSC-style-name and Sequence-Length"
            )
        lengths: dict[str, int] = {}
        sequence_aliases: dict[str, str] = {}
        data = (line for line in handle if not line.startswith("#"))
        for record in csv.DictReader(data, fieldnames=header, delimiter="\t"):
            name = (record["UCSC-style-name"] or "").strip()
            if not name or name == "na":
                continue
            try:
                length = int(record["Sequence-Length"] or "")
            except ValueError:
                continue
            if lengths.setdefault(name, length) != length:
                raise hg.IndexBuildError(
                    f"assembly report repeats {name} with different lengths"
                )
            accession = (record.get("GenBank-Accn") or "").strip()
            if accession and accession != "na":
                if sequence_aliases.setdefault(accession, name) != name:
                    raise hg.IndexBuildError(
                        f"assembly report maps {accession} to multiple UCSC names"
                    )
    if not lengths:
        raise hg.IndexBuildError("assembly report lists no UCSC sequence names")
    for name, expected in PINNED_MM10_SEQUENCES:
        got = lengths.get(name)
        if got != expected:
            raise hg.IndexBuildError(
                f"assembly report length for {name} is {got}, "
                f"not the pinned mm10 length {expected}"
            )
    return lengths, sequence_aliases


def _rows(path: Path):
    if path.suffix == ".gz":
        handle_cm = gzip.open(path, "rt", encoding="utf-8", newline="")
    else:
        handle_cm = path.open("r", encoding="utf-8", newline="")
    with handle_cm as handle:
        for line in handle:
            line = line.rstrip("\r\n")
            if line:
                yield line.split("\t")


def _keep_text(value: str) -> str:
    text = value.strip()
    if not text or text == r"\N":
        return ""
    return text


def read_mgi_aliases(
    external_db: Path,
    gene: Path,
    xref: Path,
    synonym: Path,
    object_xref: Path,
) -> tuple[dict[str, set[str]], dict[str, int]]:
    mgi_ids: set[str] = set()
    hgnc_ids: set[str] = set()
    for parts in _rows(external_db):
        if len(parts) < 2:
            continue
        if parts[1] == "MGI":
            mgi_ids.add(parts[0])
        elif parts[1] == "HGNC":
            hgnc_ids.add(parts[0])
    if not mgi_ids:
        raise hg.IndexBuildError("Ensembl external_db has no MGI database")
    labels: dict[str, set[str]] = {}
    hgnc_xrefs: set[str] = set()
    for parts in _rows(xref):
        if len(parts) < 4:
            continue
        if parts[1] in mgi_ids:
            spellings = labels.setdefault(parts[0], set())
            label = _keep_text(parts[3])
            if label:
                spellings.add(label)
        elif parts[1] in hgnc_ids:
            hgnc_xrefs.add(parts[0])
    if not labels:
        raise hg.IndexBuildError("Ensembl xref table has no MGI entries")
    for parts in _rows(synonym):
        if len(parts) < 2:
            continue
        bucket = labels.get(parts[0])
        if bucket is None:
            continue
        text = _keep_text(parts[1])
        if text:
            bucket.add(text)
    stable: dict[str, str] = {}
    for parts in _rows(gene):
        if len(parts) < 13:
            continue
        ensembl = parts[12]
        if ensembl.startswith("ENSMUSG") and ensembl[7:].isdigit():
            stable[parts[0]] = ensembl
    if not stable:
        raise hg.IndexBuildError("Ensembl gene table has no ENSMUSG ids")
    aliases: dict[str, set[str]] = {}
    linked = 0
    ignored_hgnc = 0
    for parts in _rows(object_xref):
        if len(parts) < 4 or parts[2] != "Gene":
            continue
        xref_id = parts[3]
        if xref_id in hgnc_xrefs:
            ignored_hgnc += 1
            continue
        spellings = labels.get(xref_id)
        if not spellings:
            continue
        base = stable.get(parts[1])
        if not base:
            continue
        aliases.setdefault(base, set()).update(spellings)
        linked += 1
    if linked == 0:
        raise hg.IndexBuildError("MGI xrefs did not join to mouse genes")
    return aliases, {
        "mgi_xrefs": len(labels),
        "mgi_gene_links": linked,
        "ignored_hgnc_gene_links": ignored_hgnc,
    }


def _file_record(path: Path, canonical_name: str) -> dict[str, object]:
    digest, size = hg.sha256_file(path)
    gzip_sha, gzip_size, plain_sha = ENSEMBL_CANONICAL[canonical_name]
    expected = gzip_sha if path.suffix == ".gz" else plain_sha
    if digest != expected:
        raise hg.IndexBuildError(
            f"{canonical_name} SHA-256 does not match Ensembl release 100"
        )
    return {
        "role": canonical_name,
        "filename": path.name,
        "bytes": size,
        "sha256": digest,
        "canonical_bytes": gzip_size,
        "canonical_sha256": gzip_sha,
        "url": ENSEMBL_BASE + canonical_name,
    }


def _mapping_sha(records: list[dict[str, object]]) -> str:
    digest = hashlib.sha256()
    for record in records:
        line = f"{record['role']} {record['sha256']}\n"
        digest.update(line.encode("ascii"))
    return digest.hexdigest()


def _annotation_meta(
    gtf: Path, file_sha: str, readme: Path | None
) -> tuple[dict[str, str], str, str, str, dict[str, object] | None]:
    official = gtf.suffix == ".gz" and file_sha == OFFICIAL_GTF_SHA
    if readme is None:
        if not official:
            raise hg.IndexBuildError(
                "README is required unless --gtf is the pinned GENCODE release gzip"
            )
        return (
            {
                "annotation url": OFFICIAL_GTF_URL,
                "assembly_accession": "GCA_000001635.2",
                "genome url": (
                    "https://hgdownload.soe.ucsc.edu/goldenPath/mm10/"
                    "bigZips/chromFa.tar.gz"
                ),
            },
            "",
            "",
            OFFICIAL_ATTESTATION,
            None,
        )
    meta = read_readme(readme)
    readme_sha, _readme_size = hg.sha256_file(readme)
    if official:
        meta["annotation url"] = OFFICIAL_GTF_URL
        return meta, readme.name, readme_sha, OFFICIAL_ATTESTATION, None
    if file_sha == KNOWN_LOCAL_SHA:
        return (
            meta,
            readme.name,
            readme_sha,
            LOCAL_CROSSCHECK_ATTESTATION,
            dict(UPSTREAM_CROSSCHECK),
        )
    return meta, readme.name, readme_sha, ATTESTATION, None


def build_index(
    *,
    gtf: Path,
    assembly_report: Path,
    readme: Path | None,
    external_db: Path,
    gene: Path,
    xref: Path,
    synonym: Path,
    object_xref: Path,
) -> dict[str, object]:
    lengths, sequence_aliases = read_assembly_report(assembly_report)
    (
        rows,
        drop_counts,
        gtf_sha256,
        gtf_size,
        plain_sha256,
        plain_size,
    ) = hg.read_genes(
        gtf, set(lengths), gene_id_pattern=_ENSMUSG, sequence_aliases=sequence_aliases
    )
    if plain_sha256 not in {OFFICIAL_PLAIN_SHA, KNOWN_LOCAL_SHA}:
        raise hg.IndexBuildError(
            "GTF SHA-256 does not match a pinned GENCODE M25 annotation"
        )
    if drop_counts["dropped_unknown_sequence"]:
        raise hg.IndexBuildError(
            "assembly report does not identify all sequences with M25 gene records. "
            "The upstream GTF needs the GenBank-Accn to UCSC-style-name mapping."
        )
    if (
        gtf.suffix == ".gz"
        and gtf_sha256 == OFFICIAL_GTF_SHA
        and (
            gtf_size != OFFICIAL_GTF_BYTES
            or plain_sha256 != OFFICIAL_PLAIN_SHA
            or plain_size != OFFICIAL_PLAIN_BYTES
        )
    ):
        raise hg.IndexBuildError(
            "pinned GENCODE gzip did not decompress to the pinned annotation"
        )
    if gtf_sha256 == KNOWN_LOCAL_SHA and gtf_size != KNOWN_LOCAL_BYTES:
        raise hg.IndexBuildError("known local GTF length does not match")
    meta, readme_name, readme_sha, attestation, upstream = _annotation_meta(
        gtf, gtf_sha256, readme
    )
    rows.sort(
        key=lambda row: (
            row.symbol.upper(),
            hg._chrom_key(row.chrom),
            row.start,
            row.end,
            row.gene_id,
        )
    )
    chromosomes = sorted({row.chrom for row in rows}, key=hg._chrom_key)
    chrom_index = {name: index for index, name in enumerate(chromosomes)}
    alias_map, alias_counts = read_mgi_aliases(
        external_db, gene, xref, synonym, object_xref
    )
    alias_links = hg._alias_links(rows, alias_map)
    report_sha, _report_size = hg.sha256_file(assembly_report)
    alias_files = [
        _file_record(path, name)
        for path, name in (
            (external_db, "external_db.txt.gz"),
            (gene, "gene.txt.gz"),
            (xref, "xref.txt.gz"),
            (synonym, "external_synonym.txt.gz"),
            (object_xref, "object_xref.txt.gz"),
        )
    ]
    synonym_record = alias_files[3]
    provenance: dict[str, object] = {
        "annotation": {
            "release": RELEASE,
            "source_url": meta["annotation url"],
            "filename": gtf.name,
            "bytes": gtf_size,
            "sha256": gtf_sha256,
            "uncompressed_bytes": plain_size,
            "uncompressed_sha256": plain_sha256,
            "feature": "gene",
            "assembly_accession": meta["assembly_accession"],
            "genome_url": meta.get("genome url", ""),
            "genomepy_version": meta.get("genomepy version", ""),
            "readme_date": meta.get("date", ""),
            "readme_filename": readme_name,
            "readme_sha256": readme_sha,
            "attestation": attestation,
            "license": hg.GENCODE_LICENSE,
        },
        "assembly_report": {
            "filename": assembly_report.name,
            "sha256": report_sha,
            "admitted_sequence_names": len(lengths),
            "source_note": (
                "UCSC-style-name and Sequence-Length columns. Every pinned "
                "wuepgg 62.0.2 mm10 sequence is present at that viewer "
                "length. This record identifies the report by sha256."
            ),
        },
        "aliases": {
            "source": "Ensembl 100 mus_musculus_core_100_38 MGI synonyms",
            "source_url": ENSEMBL_BASE + "external_synonym.txt.gz",
            "release": "Ensembl 100",
            "db_name": "MGI",
            "external_db_id": "1400",
            "excluded_db_names": ["HGNC"],
            "canonical_sha256": synonym_record["canonical_sha256"],
            "canonical_bytes": synonym_record["canonical_bytes"],
            "filename": synonym_record["filename"],
            "bytes": synonym_record["bytes"],
            "sha256": synonym_record["sha256"],
            "mapping_sha256": _mapping_sha(alias_files),
            "mapping_note": (
                "mapping_sha256 is SHA-256 over lines of 'role sha256' for "
                "external_db, gene, xref, external_synonym, and object_xref, "
                "in that order. sha256 itself is the synonym file."
            ),
            "files": alias_files,
            "fields": ["display_label", "synonym"],
            "join": "Ensembl gene id without version",
            "mgi_xrefs": alias_counts["mgi_xrefs"],
            "mgi_gene_links": alias_counts["mgi_gene_links"],
            "ignored_hgnc_gene_links": alias_counts["ignored_hgnc_gene_links"],
            "license": MGI_LICENSE,
        },
        "viewer_scaffold": {
            "id": "mm10",
            "version": "62.0.2",
            "sequences": [name for name, _length in PINNED_MM10_SEQUENCES],
            "lengths": [length for _name, length in PINNED_MM10_SEQUENCES],
            "note": (
                "These are the only sequences in the pinned viewer genome. "
                "chrM is last, which matches wuepgg 62.0.2 mm10 and not "
                "mm39. A region on any other sequence is not sent to the "
                "viewer, because its region parser throws when the sequence "
                "is missing and the viewer unmounts."
            ),
        },
        "generator": SCRIPT,
        "coordinate_system": "1-based-closed",
        "region_format": "chrom:start-end",
        "coordinate_note": hg.COORDINATE_NOTE,
        "assembly_note": ASSEMBLY_NOTE,
        "counts": {
            "genes": len(rows),
            "aliases": len(alias_links),
            "chromosomes": len(chromosomes),
            "dropped_incomplete": drop_counts["dropped_incomplete"],
            "dropped_unknown_sequence": drop_counts["dropped_unknown_sequence"],
        },
    }
    if upstream is not None:
        provenance["upstream_crosscheck"] = upstream
    return {
        "format": hg.FORMAT,
        "assembly": "mm10",
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


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--assembly-report", type=Path, required=True)
    parser.add_argument("--readme", type=Path)
    parser.add_argument("--ensembl-external-db", type=Path, required=True)
    parser.add_argument("--ensembl-gene", type=Path, required=True)
    parser.add_argument("--ensembl-xref", type=Path, required=True)
    parser.add_argument("--ensembl-synonym", type=Path, required=True)
    parser.add_argument("--ensembl-object-xref", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    args = parser.parse_args(argv)
    try:
        payload = build_index(
            gtf=args.gtf,
            assembly_report=args.assembly_report,
            readme=args.readme,
            external_db=args.ensembl_external_db,
            gene=args.ensembl_gene,
            xref=args.ensembl_xref,
            synonym=args.ensembl_synonym,
            object_xref=args.ensembl_object_xref,
        )
        hg.write_index(payload, args.output, args.provenance)
    except (OSError, UnicodeError, hg.IndexBuildError) as error:
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
