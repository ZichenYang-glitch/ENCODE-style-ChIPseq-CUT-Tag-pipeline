#!/usr/bin/env python3
"""Stage 27a: Public data validation input inventory (stdlib-only, no downloads).

Prints a planned dataset inventory as JSON or TSV.  Does NOT download data.
Intended for human review before manual data preparation or future automation.

Usage:
    python3 scripts/prepare_public_validation_inputs.py          # TSV to stdout
    python3 scripts/prepare_public_validation_inputs.py --json   # JSON to stdout
    python3 scripts/prepare_public_validation_inputs.py --dry-run
"""

import argparse
import json
import sys

_VALIDATION_DATASETS = [
    {
        "queue": "tf_chip_cebpb",
        "assay": "chipseq",
        "target": "CEBPB",
        "accession": "ENCSR000DYI",
        "source_url": "https://www.encodeproject.org/experiments/ENCSR000DYI/",
        "genome": "hg38",
        "layout": "PE",
        "peak_mode": "narrow",
        "n_treatment_bioreps": 2,
        "n_control_samples": 1,
        "has_control": True,
        "stage5_idr": True,
        "notes": "TF ChIP-seq with control — validates pooled, IDR, and control paths.",
    },
    {
        "queue": "broad_histone_h3k27me3",
        "assay": "chipseq",
        "target": "H3K27me3",
        "accession": "ENCSR000AKB",
        "source_url": "https://www.encodeproject.org/experiments/ENCSR000AKB/",
        "genome": "hg38",
        "layout": "PE",
        "peak_mode": "broad",
        "n_treatment_bioreps": 2,
        "n_control_samples": 1,
        "has_control": True,
        "stage5_idr": False,
        "notes": "Broad histone mark — validates broad-peak MACS3, histone pooled QC, no IDR.",
    },
    {
        "queue": "atac_keratinocyte",
        "assay": "atac",
        "target": "ATAC",
        "accession": "ENCSR254KDA",
        "source_url": "https://www.encodeproject.org/experiments/ENCSR254KDA/",
        "genome": "hg38",
        "layout": "PE",
        "peak_mode": "narrow",
        "n_treatment_bioreps": 2,
        "n_control_samples": 0,
        "has_control": False,
        "stage5_idr": False,
        "notes": "ATAC-seq baseline — validates ATAC dispatch, Tn5-aware MACS3, TSS/FRiP.",
    },
    {
        "queue": "cuttag_h3k27me3",
        "assay": "cuttag",
        "target": "H3K27me3",
        "accession": "GSE145187",
        "source_url": "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE145187",
        "genome": "hg38",
        "layout": "PE",
        "peak_mode": "narrow",
        "n_treatment_bioreps": 2,
        "n_control_samples": 1,
        "has_control": True,
        "stage5_idr": False,
        "seacr_enabled": True,
        "notes": "CUT&Tag with IgG control — validates SEACR sidecar, fragment-size QC, no CUT&Tag IDR.",
    },
]

_MANIFEST_COLUMNS = [
    "queue", "assay", "target", "accession", "source_url",
    "genome", "layout", "peak_mode", "n_treatment_bioreps",
    "n_control_samples", "stage5_idr", "seacr_enabled", "notes",
]


def _collect_columns(dataset):
    row = {}
    for col in _MANIFEST_COLUMNS:
        val = dataset.get(col, "")
        if isinstance(val, bool):
            val = str(val).lower()
        row[col] = str(val)
    return row


def _print_tsv(datasets):
    print("\t".join(_MANIFEST_COLUMNS))
    for ds in datasets:
        row = _collect_columns(ds)
        print("\t".join(row[col] for col in _MANIFEST_COLUMNS))


def _print_json(datasets):
    print(json.dumps(datasets, indent=2))


_METADATA_CHECK_FIELDS = [
    "assay", "target", "genome", "layout", "n_treatment_bioreps",
    "n_control_samples", "has_control", "peak_mode", "stage5_idr",
]


def _print_report_stubs(datasets):
    """Print queue name, accession, and expected report stub path."""
    print("=== Public Data Execution Report Stubs ===\n")
    for ds in datasets:
        stub_path = (
            f"docs/release-checks/public-data-runs/{ds['queue']}.md"
        )
        print(f"{ds['queue']:30s} {ds['accession']:15s}  {stub_path}")
    print()
    print("Template: docs/release-checks/public-data-execution-report-template.md")
    print("No downloads performed.")


def _print_metadata_checklist(datasets):
    """Print metadata fields to verify for each dataset (no network)."""
    print("=== Public Validation Dataset Metadata Checklist ===\n")
    print("Verify the following fields before downloading data.\n")
    print("ENCODE datasets: check at https://www.encodeproject.org/experiments/<accession>/")
    print("GEO datasets:    check at https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<accession>\n")

    for ds in datasets:
        print(f"--- {ds['queue']} ({ds['accession']}) ---")
        print(f"  Source URL:       {ds['source_url']}")
        for field in _METADATA_CHECK_FIELDS:
            val = ds.get(field, "")
            if isinstance(val, bool):
                val = str(val).lower()
            print(f"  {field:25s} {val}")
        print(f"  notes:             {ds['notes']}")
        print()

    print("=== Acceptance Criteria ===")
    print("- All metadata fields above must match the dataset's official record.")
    print("- If assay or target differ, the dataset is rejected for that queue.")
    print("- If replicate count differs, re-evaluate sample sheet structure.")
    print("- If genome build differs from hg38, accept only if mm10 alternative is suitable.")
    print("- If data is not publicly accessible or requires dbGaP login, reject.")
    print()
    print("=== Manual Verification Steps ===")
    print("1. Navigate to the source URL for each accession.")
    print("2. Confirm the organism, cell type/tissue, and treatment.")
    print("3. Confirm FASTQ or BAM availability (not just processed peaks/signals).")
    print("4. Record the exact file URLs or GEO/SRA run accessions for downstream fetch.")
    print("5. Note library layout (PE/SE) and read length.")
    print()
    print("No network requests were made by this script.")


def main():
    parser = argparse.ArgumentParser(
        description="Stage 27a: public data validation input inventory"
    )
    parser.add_argument(
        "--json", action="store_true", default=False,
        help="Output JSON instead of TSV",
    )
    parser.add_argument(
        "--dry-run", action="store_true", default=False,
        help="Print planned actions without executing",
    )
    parser.add_argument(
        "--check-metadata", action="store_true", default=False,
        help="Print metadata fields to verify for each dataset (no network)",
    )
    parser.add_argument(
        "--report-stubs", action="store_true", default=False,
        help="Print queue names and expected report stub paths (no download)",
    )
    args = parser.parse_args()

    if args.report_stubs:
        _print_report_stubs(_VALIDATION_DATASETS)
        return

    if args.check_metadata:
        _print_metadata_checklist(_VALIDATION_DATASETS)
        return

    if args.dry_run:
        print("[dry-run] Would print validation dataset inventory")
        print(f"[dry-run] {len(_VALIDATION_DATASETS)} dataset(s) defined")
        for ds in _VALIDATION_DATASETS:
            print(f"  {ds['queue']}: {ds['assay']} / {ds['target']} "
                  f"({ds['accession']}) — {ds['notes']}")
        print("[dry-run] No downloads performed")
        return

    if args.json:
        _print_json(_VALIDATION_DATASETS)
    else:
        _print_tsv(_VALIDATION_DATASETS)


if __name__ == "__main__":
    main()
