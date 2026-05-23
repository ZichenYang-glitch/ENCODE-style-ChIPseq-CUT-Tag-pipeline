#!/usr/bin/env python3
"""Generate a project-level result manifest recording core output existence.

Stage 25: Minimal manifest — existence-only, TSV format, stdlib-only.
Uses validate_samples to normalize config and sample defaults matching the DAG.

Usage (standalone):
    python3 scripts/make_manifest.py \\
        --config config/config.yaml \\
        --output results/multiqc/result_manifest.tsv \\
        [--strict]

Usage (from Snakemake rule):
    python3 scripts/make_manifest.py \\
        --config-json '<json_dump_of_validated_config>' \\
        --output results/multiqc/result_manifest.tsv

Output columns:
    sample_id  experiment_id  assay  target  genome  output_type  method  path  status  qc_flag
"""

import argparse
import csv
import json
import os
import sys

# Make scripts/ importable when run standalone.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_samples import (
    validate_config,
    load_and_validate_samples,
)


_NA = "NA"
_MANIFEST_COLUMNS = [
    "sample_id", "experiment_id", "assay", "target", "genome",
    "output_type", "method", "path", "status", "qc_flag",
]


def _resolve_path(outdir, *parts):
    return os.path.join(outdir, *parts)


def _check_exists(filepath):
    return "present" if os.path.exists(filepath) else "missing"


def _has_chrom_sizes(genome, genome_resources):
    entry = genome_resources.get(genome, {})
    return bool(entry.get("chrom_sizes", ""))


def _add_row(rows, sample_id, experiment_id, assay, target, genome,
             output_type, method, path, check_exists=True):
    """Append a manifest row dict."""
    status = _check_exists(path) if check_exists else "not_applicable"
    rows.append({
        "sample_id": sample_id or "",
        "experiment_id": experiment_id or "",
        "assay": assay or "",
        "target": target or "",
        "genome": genome or "",
        "output_type": output_type,
        "method": method,
        "path": path,
        "status": status,
        "qc_flag": _NA,
    })


def _build_sample_rows(samples, outdir, signal_tracks, genomic_resources):
    """Per-sample rows using validated sample dicts (key='id')."""
    rows = []
    for s in sorted(samples, key=lambda x: x["id"]):
        sid = s["id"]
        assay = s.get("assay", "")
        target = s.get("target", "")
        genome = s.get("genome", "")

        # Always
        _add_row(rows, sid, "", assay, target, genome,
                 "final_bam", "bowtie2+samtools",
                 _resolve_path(outdir, sid, "02_align", f"{sid}.final.bam"))
        _add_row(rows, sid, "", assay, target, genome,
                 "final_bai", "samtools index",
                 _resolve_path(outdir, sid, "02_align", f"{sid}.final.bam.bai"))
        _add_row(rows, sid, "", assay, target, genome,
                 "cpm_bigwig", "bamCoverage",
                 _resolve_path(outdir, sid, "03_bigwig", f"{sid}.CPM.bw"))
        _add_row(rows, sid, "", assay, target, genome,
                 "macs3_peak", "macs3_callpeak",
                 _resolve_path(outdir, sid, "04_peaks", sid))
        _add_row(rows, sid, "", assay, target, genome,
                 "qc_summary", "assemble_qc_summary",
                 _resolve_path(outdir, sid, "01_qc", f"{sid}.qc_summary.tsv"))

        # Gated: signal_tracks
        if signal_tracks:
            _add_row(rows, sid, "", assay, target, genome,
                     "macs3_fe_bdg", "macs3_bdgcmp",
                     _resolve_path(outdir, sid, "03_signal", f"{sid}.FE.bdg"))
            _add_row(rows, sid, "", assay, target, genome,
                     "macs3_ppois_bdg", "macs3_bdgcmp",
                     _resolve_path(outdir, sid, "03_signal", f"{sid}.ppois.bdg"))
        else:
            _add_row(rows, sid, "", assay, target, genome,
                     "macs3_fe_bdg", "macs3_bdgcmp", "", check_exists=False)
            _add_row(rows, sid, "", assay, target, genome,
                     "macs3_ppois_bdg", "macs3_bdgcmp", "", check_exists=False)

        # Gated: chrom_sizes + signal_tracks
        if signal_tracks and _has_chrom_sizes(genome, genomic_resources):
            _add_row(rows, sid, "", assay, target, genome,
                     "macs3_fe_bw", "macs3_bdgcmp+bedGraphToBigWig",
                     _resolve_path(outdir, sid, "03_signal", f"{sid}.FE.bw"))
            _add_row(rows, sid, "", assay, target, genome,
                     "macs3_ppois_bw", "macs3_bdgcmp+bedGraphToBigWig",
                     _resolve_path(outdir, sid, "03_signal", f"{sid}.ppois.bw"))
        else:
            _add_row(rows, sid, "", assay, target, genome,
                     "macs3_fe_bw", "macs3_bdgcmp+bedGraphToBigWig",
                     "", check_exists=False)
            _add_row(rows, sid, "", assay, target, genome,
                     "macs3_ppois_bw", "macs3_bdgcmp+bedGraphToBigWig",
                     "", check_exists=False)

    return rows


def _build_experiment_rows(samples, outdir, signal_tracks, genomic_resources,
                           stage4b):
    """Per-experiment rows using validated sample dicts.

    Matches Snakefile Stage 4b gating:
    - pooled outputs only for multi-biorep experiments (>=2 unique bio_rep values)
    - biorep rows follow _biorep_expand_pairs(): always for multi-biorep,
      or for single-biorep when that biorep has >=2 technical replicates
    """
    rows = []
    if not stage4b:
        return rows

    # Group treatment samples by experiment + count unique biological_replicates
    exp_bioreps: dict[str, set[int]] = {}
    exp_samples: dict[str, list[dict]] = {}
    for s in samples:
        exp = s.get("experiment", "")
        if not exp:
            continue
        exp_bioreps.setdefault(exp, set()).add(int(s.get("biological_replicate", 1)))
        exp_samples.setdefault(exp, []).append(s)

    for exp in sorted(exp_samples):
        all_samples = exp_samples[exp]
        unique_bioreps = exp_bioreps.get(exp, set())
        is_multi = len(unique_bioreps) >= 2

        first = all_samples[0]
        assay = first.get("assay", "")
        target = first.get("target", "")
        genome = first.get("genome", "")

        # Pooled outputs only for multi-biorep experiments
        if is_multi:
            _add_row(rows, "", exp, assay, target, genome,
                     "pooled_final_bam", "samtools merge",
                     _resolve_path(outdir, "experiments", exp, "02_align",
                                   f"{exp}.pooled.final.bam"))
            _add_row(rows, "", exp, assay, target, genome,
                     "pooled_final_bai", "samtools index",
                     _resolve_path(outdir, "experiments", exp, "02_align",
                                   f"{exp}.pooled.final.bam.bai"))
            _add_row(rows, "", exp, assay, target, genome,
                     "pooled_macs3_peak", "macs3_callpeak",
                     _resolve_path(outdir, "experiments", exp, "04_peaks",
                                   "pooled", f"{exp}_pooled_peaks"))
            _add_row(rows, "", exp, assay, target, genome,
                     "pooled_qc_summary", "pooled_qc_summary",
                     _resolve_path(outdir, "experiments", exp, "01_qc",
                                   f"{exp}.pooled_qc_summary.tsv"))

        # Biorep rows: match _biorep_expand_pairs()
        # - multi-biorep: all bioreps always included
        # - single-biorep: include only if that biorep has >=2 technical replicates
        for br in sorted(unique_bioreps):
            n_techreps = sum(1 for s in all_samples
                           if int(s.get("biological_replicate", 1)) == br)
            if is_multi or n_techreps >= 2:
                _add_row(rows, "", exp, assay, target, genome,
                         f"biorep{br}_final_bam", "samtools merge/symlink",
                         _resolve_path(outdir, "experiments", exp, "02_align",
                                       f"biorep{br}.final.bam"))
                _add_row(rows, "", exp, assay, target, genome,
                         f"biorep{br}_final_bai", "samtools index",
                         _resolve_path(outdir, "experiments", exp, "02_align",
                                       f"biorep{br}.final.bam.bai"))

        # Gated signal tracks for pooled (multi-biorep only)
        if is_multi:
            if signal_tracks:
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_fe_bdg", "macs3_bdgcmp",
                         _resolve_path(outdir, "experiments", exp, "03_signal",
                                       f"{exp}.pooled.FE.bdg"))
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_ppois_bdg", "macs3_bdgcmp",
                         _resolve_path(outdir, "experiments", exp, "03_signal",
                                       f"{exp}.pooled.ppois.bdg"))
            else:
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_fe_bdg", "macs3_bdgcmp", "", check_exists=False)
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_ppois_bdg", "macs3_bdgcmp", "", check_exists=False)

            if signal_tracks and _has_chrom_sizes(genome, genomic_resources):
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_fe_bw", "macs3_bdgcmp+bedGraphToBigWig",
                         _resolve_path(outdir, "experiments", exp, "03_signal",
                                       f"{exp}.pooled.FE.bw"))
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_ppois_bw", "macs3_bdgcmp+bedGraphToBigWig",
                         _resolve_path(outdir, "experiments", exp, "03_signal",
                                       f"{exp}.pooled.ppois.bw"))
            else:
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_fe_bw", "macs3_bdgcmp+bedGraphToBigWig",
                         "", check_exists=False)
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_ppois_bw", "macs3_bdgcmp+bedGraphToBigWig",
                         "", check_exists=False)

    return rows


def _build_idr_rows(samples, outdir, stage5):
    """IDR rows using validated sample dicts + Stage 5 DAG gating logic.

    Matches Snakefile Stage 5 semantics:
    - stage5: true
    - exactly 2 unique biological_replicate values (not sample rows)
    - assay: chipseq, peak_mode: narrow
    """
    rows = []
    if not stage5:
        return rows

    exp_bioreps: dict[str, set[int]] = {}
    exp_samples: dict[str, list[dict]] = {}
    for s in samples:
        exp = s.get("experiment", "")
        if not exp:
            continue
        exp_bioreps.setdefault(exp, set()).add(int(s.get("biological_replicate", 1)))
        exp_samples.setdefault(exp, []).append(s)

    for exp in sorted(exp_samples):
        unique_bioreps = exp_bioreps.get(exp, set())
        if len(unique_bioreps) != 2:
            continue
        first = exp_samples[exp][0]
        if first.get("assay") != "chipseq" or first.get("peak_mode") != "narrow":
            continue

        assay = first.get("assay", "")
        target = first.get("target", "")
        genome = first.get("genome", "")

        _add_row(rows, "", exp, assay, target, genome,
                 "idr_conservative", "idr",
                 _resolve_path(outdir, "experiments", exp, "06_idr",
                               "final", "conservative.narrowPeak"))
        _add_row(rows, "", exp, assay, target, genome,
                 "idr_optimal", "idr",
                 _resolve_path(outdir, "experiments", exp, "06_idr",
                               "final", "optimal.narrowPeak"))
        _add_row(rows, "", exp, assay, target, genome,
                 "idr_reproducibility_summary", "idr",
                 _resolve_path(outdir, "experiments", exp, "06_idr",
                               "final", "reproducibility_summary.tsv"))
    return rows


def _build_project_rows(outdir, multiqc_enabled):
    rows = []
    _add_row(rows, "", "", "", "", "",
             "stage3_qc_summary", "aggregate_qc_summary",
             _resolve_path(outdir, "multiqc", "stage3_qc_summary.tsv"))
    if multiqc_enabled:
        _add_row(rows, "", "", "", "", "",
                 "multiqc_report", "multiqc",
                 _resolve_path(outdir, "multiqc", "multiqc_report.html"))
    else:
        _add_row(rows, "", "", "", "", "",
                 "multiqc_report", "multiqc", "", check_exists=False)
    return rows


def main():
    parser = argparse.ArgumentParser(
        description="Generate minimal project-level result manifest (Stage 25)"
    )
    parser.add_argument("--config", default=None,
                        help="Path to config.yaml (standalone use)")
    parser.add_argument("--config-json", default=None,
                        help="JSON dump of validated config (from Snakemake)")
    parser.add_argument("--output", required=True, help="Path to output TSV")
    parser.add_argument("--strict", action="store_true", default=False,
                        help="Exit non-zero if any row has status=missing")
    args = parser.parse_args()

    # Load config — prefer JSON from Snakemake, fall back to YAML for standalone
    if args.config_json:
        cfg = json.loads(args.config_json)
    elif args.config:
        cfg = validate_config(
            # Use validate_samples._load_yaml for stdlib-safe YAML loading
            __import__("validate_samples")._load_yaml(args.config)
        )
    else:
        print("ERROR: --config or --config-json is required", file=sys.stderr)
        sys.exit(1)

    outdir = cfg.get("outdir", "results")
    qc = cfg.get("qc", {})
    signal_tracks = qc.get("signal_tracks", True)
    stage4b = cfg.get("stage4b", True)
    stage5 = cfg.get("stage5", False)
    multiqc_enabled = cfg.get("multiqc", True)
    genomic_resources = cfg.get("genome_resources", {})

    # Load and validate samples using the shared normalizer
    samples_path = cfg.get("samples", "config/samples.tsv")
    samples = load_and_validate_samples(
        samples_path,
        use_control=cfg.get("use_control", False),
        stage5_enabled=stage5,
    )
    # Filter to treatment only (manifest only records treatment outputs)
    treatment_samples = [s for s in samples if s.get("role") == "treatment"]

    all_rows = []
    all_rows.extend(_build_sample_rows(
        treatment_samples, outdir, signal_tracks, genomic_resources))
    all_rows.extend(_build_experiment_rows(
        treatment_samples, outdir, signal_tracks, genomic_resources, stage4b))
    all_rows.extend(_build_idr_rows(treatment_samples, outdir, stage5))
    all_rows.extend(_build_project_rows(outdir, multiqc_enabled))

    os.makedirs(os.path.dirname(args.output) or ".", exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_MANIFEST_COLUMNS, delimiter="\t",
            lineterminator="\n", extrasaction="ignore",
        )
        writer.writeheader()
        for row in all_rows:
            writer.writerow(row)

    missing_count = sum(1 for r in all_rows if r["status"] == "missing")
    na_count = sum(1 for r in all_rows if r["status"] == "not_applicable")
    print(f"Manifest written: {args.output}")
    print(f"  {len(all_rows)} rows: "
          f"{len(all_rows) - missing_count - na_count} present, "
          f"{missing_count} missing, "
          f"{na_count} not_applicable")

    if args.strict and missing_count > 0:
        print(f"ERROR: --strict mode, {missing_count} missing output(s)",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
