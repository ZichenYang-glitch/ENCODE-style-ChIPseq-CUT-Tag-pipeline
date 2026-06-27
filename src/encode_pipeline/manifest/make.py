"""Generate a project-level result manifest recording core output existence.

Stage 25: Minimal manifest — existence-only, TSV format, stdlib-only.
Uses validate_samples to normalize config and sample defaults matching the DAG.

Output columns:
    sample_id  experiment_id  assay  target  genome  output_type  method  path  status  qc_flag
"""

import csv
import json
import os

import yaml

from encode_pipeline.config.validate import validate_config
from encode_pipeline.samples.load import load_and_validate_samples


_NA = "NA"
_MANIFEST_COLUMNS = [
    "sample_id", "experiment_id", "assay", "target", "genome",
    "output_type", "method", "path", "status", "qc_flag",
]

MANIFEST_COLUMNS = list(_MANIFEST_COLUMNS)


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


def _is_mnase(sample_dict):
    """Return True if the sample has assay=mnase."""
    return sample_dict.get("assay", "") == "mnase"


def _build_sample_rows(samples, outdir, signal_tracks, genomic_resources):
    """Per-sample rows using validated sample dicts (key='id')."""
    rows = []
    for s in sorted(samples, key=lambda x: x["id"]):
        sid = s["id"]
        assay = s.get("assay", "")
        target = s.get("target", "")
        genome = s.get("genome", "")
        is_mn = _is_mnase(s)

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

        # Peak-centric outputs — gated on non-MNase
        if not is_mn:
            _add_row(rows, sid, "", assay, target, genome,
                     "macs3_peak", "macs3_callpeak",
                     _resolve_path(outdir, sid, "04_peaks", sid))
            _add_row(rows, sid, "", assay, target, genome,
                     "qc_summary", "assemble_qc_summary",
                     _resolve_path(outdir, sid, "01_qc", f"{sid}.qc_summary.tsv"))
        else:
            _add_row(rows, sid, "", assay, target, genome,
                     "macs3_peak", "macs3_callpeak",
                     "", check_exists=False)
            _add_row(rows, sid, "", assay, target, genome,
                     "qc_summary", "assemble_qc_summary",
                     "", check_exists=False)

        # Gated: signal_tracks (peak-centric only)
        if signal_tracks and not is_mn:
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

        # Gated: chrom_sizes + signal_tracks (peak-centric only)
        if signal_tracks and not is_mn and _has_chrom_sizes(genome, genomic_resources):
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

        # MNase-specific outputs (Stage 39-40)
        if is_mn:
            # Stage 39: mono BAM, dyad BW, mono occupancy BW
            _add_row(rows, sid, "", assay, target, genome,
                     "mnase_mono_bam", "alignmentSieve",
                     _resolve_path(outdir, sid, "03_fragments", f"{sid}.mono.bam"))
            _add_row(rows, sid, "", assay, target, genome,
                     "mnase_mono_bai", "samtools index",
                     _resolve_path(outdir, sid, "03_fragments", f"{sid}.mono.bam.bai"))
            _add_row(rows, sid, "", assay, target, genome,
                     "mnase_dyad_bigwig", "bamCoverage --MNase",
                     _resolve_path(outdir, sid, "04_signal", f"{sid}.dyad.CPM.bw"))
            _add_row(rows, sid, "", assay, target, genome,
                     "mnase_mono_bigwig", "bamCoverage",
                     _resolve_path(outdir, sid, "04_signal", f"{sid}.mono.CPM.bw"))
            # Stage 40: sub-nucleosome BAM, di-nucleosome BAM, QC summary
            _add_row(rows, sid, "", assay, target, genome,
                     "mnase_sub_bam", "alignmentSieve",
                     _resolve_path(outdir, sid, "03_fragments", f"{sid}.sub.bam"))
            _add_row(rows, sid, "", assay, target, genome,
                     "mnase_sub_bai", "samtools index",
                     _resolve_path(outdir, sid, "03_fragments", f"{sid}.sub.bam.bai"))
            _add_row(rows, sid, "", assay, target, genome,
                     "mnase_di_bam", "alignmentSieve",
                     _resolve_path(outdir, sid, "03_fragments", f"{sid}.di.bam"))
            _add_row(rows, sid, "", assay, target, genome,
                     "mnase_di_bai", "samtools index",
                     _resolve_path(outdir, sid, "03_fragments", f"{sid}.di.bam.bai"))
            _add_row(rows, sid, "", assay, target, genome,
                     "mnase_qc_summary", "mnase_qc_summary.py",
                     _resolve_path(outdir, sid, "01_qc", f"{sid}.mnase_qc_summary.tsv"))

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
        is_mn = _is_mnase(first)

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

            # Peak-centric pooled outputs — gated on non-MNase
            if not is_mn:
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_macs3_peak", "macs3_callpeak",
                         _resolve_path(outdir, "experiments", exp, "04_peaks",
                                       "pooled", f"{exp}_pooled_peaks"))
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_qc_summary", "pooled_qc_summary",
                         _resolve_path(outdir, "experiments", exp, "01_qc",
                                       f"{exp}.pooled_qc_summary.tsv"))
            else:
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_macs3_peak", "macs3_callpeak",
                         "", check_exists=False)
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_qc_summary", "pooled_qc_summary",
                         "", check_exists=False)

            # MNase-specific pooled outputs (Stage 39)
            if is_mn:
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_mnase_mono_bam", "alignmentSieve",
                         _resolve_path(outdir, "experiments", exp,
                                       "03_fragments", f"{exp}.pooled.mono.bam"))
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_mnase_mono_bai", "samtools index",
                         _resolve_path(outdir, "experiments", exp,
                                       "03_fragments", f"{exp}.pooled.mono.bam.bai"))
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_mnase_dyad_bigwig", "bamCoverage --MNase",
                         _resolve_path(outdir, "experiments", exp,
                                       "04_signal", f"{exp}.pooled.dyad.CPM.bw"))
                _add_row(rows, "", exp, assay, target, genome,
                         "pooled_mnase_mono_bigwig", "bamCoverage",
                         _resolve_path(outdir, "experiments", exp,
                                       "04_signal", f"{exp}.pooled.mono.CPM.bw"))

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

        # Gated signal tracks for pooled (multi-biorep only, peak-centric only)
        if is_multi:
            if signal_tracks and not is_mn:
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

            if signal_tracks and not is_mn and _has_chrom_sizes(genome, genomic_resources):
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


def _compute_reproducibility_eligibility(config, treatment_samples):
    """Recompute reproducibility eligibility from config + sample sheet.

    Must match validate_samples.py and metadata.smk logic. Returns a dict
    with per-mode experiment lists. Does NOT rely on Snakemake globals.
    """
    repro = config.get("reproducibility", {})
    repro_enabled = repro.get("enabled", False)
    consensus_cfg = repro.get("consensus", {})
    consensus_enabled = (
        repro_enabled and consensus_cfg.get("enabled", True)
    )
    idr_cfg = repro.get("idr", {})
    stage4b = config.get("stage4b", True)
    cuttag_cfg = config.get("cuttag", {})
    seacr_cfg = cuttag_cfg.get("seacr", {})
    seacr_enabled = (
        seacr_cfg.get("enabled", False)
        if isinstance(seacr_cfg, dict) else False
    )
    seacr_mode = str(seacr_cfg.get("mode", "stringent"))

    atac_narrow_idr = repro_enabled and idr_cfg.get("atac_narrow", False)
    cuttag_narrow_idr = repro_enabled and idr_cfg.get("cuttag_narrow", False)
    chipseq_broad_idr = repro_enabled and idr_cfg.get(
        "chipseq_broad_experimental", False)
    cuttag_broad_idr = repro_enabled and idr_cfg.get(
        "cuttag_broad_experimental", False)

    # Group samples by experiment
    exp_samples = {}
    for s in treatment_samples:
        exp = s.get("experiment", "")
        if exp:
            exp_samples.setdefault(exp, []).append(s)

    # Discover eligible experiments per mode
    consensus_exps = {}
    atac_idr_exps = []
    cuttag_idr_exps = []
    chipseq_broad_idr_exps = []
    cuttag_broad_idr_exps = []
    seacr_consensus_exps = []

    for exp, rows in exp_samples.items():
        first = rows[0]
        assay = first.get("assay", "")
        peak_mode = first.get("peak_mode", "")
        layout = first.get("layout", "PE")
        bioreps = sorted({r.get("biological_replicate", 1) for r in rows})

        if len(bioreps) < 2:
            continue

        # Consensus eligibility (≥2 bioreps needed, stage4b downstream)
        if (assay, peak_mode) in [
            ("chipseq", "narrow"), ("chipseq", "broad"),
            ("cuttag", "narrow"), ("cuttag", "broad"),
            ("atac", "narrow"),
        ]:
            consensus_exps.setdefault((assay, peak_mode), []).append(exp)

        # SEACR consensus: cuttag + PE + ≥2 bioreps
        if assay == "cuttag" and layout == "PE":
            seacr_consensus_exps.append(exp)

        # IDR eligibility (exactly 2 bioreps)
        if len(bioreps) != 2:
            continue

        if atac_narrow_idr and assay == "atac" and peak_mode == "narrow":
            atac_idr_exps.append(exp)
        if cuttag_narrow_idr and assay == "cuttag" and peak_mode == "narrow":
            cuttag_idr_exps.append(exp)
        if chipseq_broad_idr and assay == "chipseq" and peak_mode == "broad":
            chipseq_broad_idr_exps.append(exp)
        if cuttag_broad_idr and assay == "cuttag" and peak_mode == "broad":
            cuttag_broad_idr_exps.append(exp)

    return {
        "stage4b_enabled": stage4b,
        "repro_enabled": repro_enabled,
        "consensus_enabled": consensus_enabled,
        "atac_narrow_idr": atac_narrow_idr,
        "cuttag_narrow_idr": cuttag_narrow_idr,
        "chipseq_broad_idr": chipseq_broad_idr,
        "cuttag_broad_idr": cuttag_broad_idr,
        "seacr_enabled": seacr_enabled,
        "seacr_mode": seacr_mode,
        "consensus_experiments": consensus_exps,
        "atac_idr_experiments": atac_idr_exps,
        "cuttag_idr_experiments": cuttag_idr_exps,
        "chipseq_broad_idr_experiments": chipseq_broad_idr_exps,
        "cuttag_broad_idr_experiments": cuttag_broad_idr_exps,
        "seacr_consensus_experiments": seacr_consensus_exps,
    }


def _build_reproducibility_rows(samples, config, outdir):
    """Stage 66: Emit mode-specific reproducibility manifest rows.

    Only user-facing catalog outputs (consensus peak, consensus summary,
    final validated peak, IDR summary). No intermediate per-biorep files.
    Omitted entirely for disabled modes and when stage4b is false.

    Uses _add_row() with literal output_type strings so Stage 49 AST
    contract can discover the full vocabulary.
    """
    rows = []
    treatment = [s for s in samples if s.get("role") == "treatment"]
    e = _compute_reproducibility_eligibility(config, treatment)

    if not e["stage4b_enabled"]:
        return rows

    def _exp_meta(exp_id):
        exp_samples = [s for s in treatment
                       if s.get("experiment") == exp_id]
        if not exp_samples:
            return None
        first = exp_samples[0]
        return {
            "assay": first.get("assay", ""),
            "target": first.get("target", ""),
            "genome": first.get("genome", ""),
        }

    # Consensus peak + summary for every eligible mode
    if e["repro_enabled"] and e["consensus_enabled"]:
        for (assay, peak_mode), exps in sorted(e["consensus_experiments"].items()):
            for exp in sorted(exps):
                meta = _exp_meta(exp)
                if meta is None:
                    continue
                suffix = "narrowPeak" if peak_mode == "narrow" else "broadPeak"
                stem = f"{outdir}/experiments/{exp}/06_reproducibility/consensus"

                if assay == "chipseq" and peak_mode == "narrow":
                    _add_row(rows, "", exp, meta["assay"], meta["target"],
                             meta["genome"], "chipseq_macs3_narrow_consensus_peak",
                             "compute_consensus.py",
                             f"{stem}/{exp}.{assay}.macs3.{peak_mode}.consensus.{suffix}")
                    _add_row(rows, "", exp, meta["assay"], meta["target"],
                             meta["genome"], "chipseq_macs3_narrow_consensus_summary",
                             "compute_consensus.py",
                             f"{stem}/{exp}.{assay}.macs3.{peak_mode}.consensus.summary.tsv")
                elif assay == "chipseq" and peak_mode == "broad":
                    _add_row(rows, "", exp, meta["assay"], meta["target"],
                             meta["genome"], "chipseq_macs3_broad_consensus_peak",
                             "compute_consensus.py",
                             f"{stem}/{exp}.{assay}.macs3.{peak_mode}.consensus.{suffix}")
                    _add_row(rows, "", exp, meta["assay"], meta["target"],
                             meta["genome"], "chipseq_macs3_broad_consensus_summary",
                             "compute_consensus.py",
                             f"{stem}/{exp}.{assay}.macs3.{peak_mode}.consensus.summary.tsv")
                elif assay == "cuttag" and peak_mode == "narrow":
                    _add_row(rows, "", exp, meta["assay"], meta["target"],
                             meta["genome"], "cuttag_macs3_narrow_consensus_peak",
                             "compute_consensus.py",
                             f"{stem}/{exp}.{assay}.macs3.{peak_mode}.consensus.{suffix}")
                    _add_row(rows, "", exp, meta["assay"], meta["target"],
                             meta["genome"], "cuttag_macs3_narrow_consensus_summary",
                             "compute_consensus.py",
                             f"{stem}/{exp}.{assay}.macs3.{peak_mode}.consensus.summary.tsv")
                elif assay == "cuttag" and peak_mode == "broad":
                    _add_row(rows, "", exp, meta["assay"], meta["target"],
                             meta["genome"], "cuttag_macs3_broad_consensus_peak",
                             "compute_consensus.py",
                             f"{stem}/{exp}.{assay}.macs3.{peak_mode}.consensus.{suffix}")
                    _add_row(rows, "", exp, meta["assay"], meta["target"],
                             meta["genome"], "cuttag_macs3_broad_consensus_summary",
                             "compute_consensus.py",
                             f"{stem}/{exp}.{assay}.macs3.{peak_mode}.consensus.summary.tsv")
                elif assay == "atac" and peak_mode == "narrow":
                    _add_row(rows, "", exp, meta["assay"], meta["target"],
                             meta["genome"], "atac_macs3_narrow_consensus_peak",
                             "compute_consensus.py",
                             f"{stem}/{exp}.{assay}.macs3.{peak_mode}.consensus.{suffix}")
                    _add_row(rows, "", exp, meta["assay"], meta["target"],
                             meta["genome"], "atac_macs3_narrow_consensus_summary",
                             "compute_consensus.py",
                             f"{stem}/{exp}.{assay}.macs3.{peak_mode}.consensus.summary.tsv")

        # SEACR consensus
        if e["seacr_enabled"]:
            for exp in sorted(e["seacr_consensus_experiments"]):
                meta = _exp_meta(exp)
                if meta is None:
                    continue
                stem = f"{outdir}/experiments/{exp}/06_reproducibility/consensus"
                mode = e["seacr_mode"]
                _add_row(rows, "", exp, meta["assay"], meta["target"],
                         meta["genome"], "cuttag_seacr_consensus_peak",
                         "compute_consensus.py",
                         f"{stem}/{exp}.cuttag.seacr.{mode}.consensus.bed")
                _add_row(rows, "", exp, meta["assay"], meta["target"],
                         meta["genome"], "cuttag_seacr_consensus_summary",
                         "compute_consensus.py",
                         f"{stem}/{exp}.cuttag.seacr.{mode}.consensus.summary.tsv")

    # IDR final peaks + summaries
    step = f"{outdir}/experiments"

    if e["atac_narrow_idr"]:
        for exp in sorted(e["atac_idr_experiments"]):
            meta = _exp_meta(exp)
            if meta is None:
                continue
            stem = f"{step}/{exp}/06_reproducibility/final"
            _add_row(rows, "", exp, meta["assay"], meta["target"],
                     meta["genome"], "atac_macs3_narrow_idr_final_peak", "idr",
                     f"{stem}/{exp}.atac.macs3.narrow.replicate_validated.idr.narrowPeak")
            _add_row(rows, "", exp, meta["assay"], meta["target"],
                     meta["genome"], "atac_macs3_narrow_idr_summary", "idr",
                     f"{stem}/reproducibility_summary.tsv")

    if e["cuttag_narrow_idr"]:
        for exp in sorted(e["cuttag_idr_experiments"]):
            meta = _exp_meta(exp)
            if meta is None:
                continue
            stem = f"{step}/{exp}/06_reproducibility/final"
            _add_row(rows, "", exp, meta["assay"], meta["target"],
                     meta["genome"], "cuttag_macs3_narrow_idr_final_peak", "idr",
                     f"{stem}/{exp}.cuttag.macs3.narrow.replicate_validated.idr.narrowPeak")
            _add_row(rows, "", exp, meta["assay"], meta["target"],
                     meta["genome"], "cuttag_macs3_narrow_idr_summary", "idr",
                     f"{stem}/reproducibility_summary.tsv")

    if e["chipseq_broad_idr"]:
        for exp in sorted(e["chipseq_broad_idr_experiments"]):
            meta = _exp_meta(exp)
            if meta is None:
                continue
            stem = f"{step}/{exp}/06_reproducibility/final"
            _add_row(rows, "", exp, meta["assay"], meta["target"],
                     meta["genome"], "chipseq_macs3_broad_idr_final_peak", "idr",
                     f"{stem}/{exp}.chipseq.macs3.broad.replicate_validated.idr.broadPeak")
            _add_row(rows, "", exp, meta["assay"], meta["target"],
                     meta["genome"], "chipseq_macs3_broad_idr_summary", "idr",
                     f"{stem}/reproducibility_summary.tsv")

    if e["cuttag_broad_idr"]:
        for exp in sorted(e["cuttag_broad_idr_experiments"]):
            meta = _exp_meta(exp)
            if meta is None:
                continue
            stem = f"{step}/{exp}/06_reproducibility/final"
            _add_row(rows, "", exp, meta["assay"], meta["target"],
                     meta["genome"], "cuttag_macs3_broad_idr_final_peak", "idr",
                     f"{stem}/{exp}.cuttag.macs3.broad.replicate_validated.idr.broadPeak")
            _add_row(rows, "", exp, meta["assay"], meta["target"],
                     meta["genome"], "cuttag_macs3_broad_idr_summary", "idr",
                     f"{stem}/reproducibility_summary.tsv")

    # Consensus final peaks (modes where IDR is NOT enabled)
    if e["repro_enabled"] and e["consensus_enabled"]:
        stem = f"{step}"

        # chipseq broad: consensus final when broad IDR not enabled
        if not e["chipseq_broad_idr"]:
            for exp in sorted(e["consensus_experiments"].get(
                    ("chipseq", "broad"), [])):
                meta = _exp_meta(exp)
                if meta is None:
                    continue
                p = f"{stem}/{exp}/06_reproducibility/final"
                _add_row(rows, "", exp, meta["assay"], meta["target"],
                         meta["genome"], "chipseq_macs3_broad_consensus_final_peak",
                         "cp",
                         f"{p}/{exp}.chipseq.macs3.broad.replicate_validated.consensus.broadPeak")

        # cuttag narrow: consensus final when cuttag narrow IDR not enabled
        if not e["cuttag_narrow_idr"]:
            for exp in sorted(e["consensus_experiments"].get(
                    ("cuttag", "narrow"), [])):
                meta = _exp_meta(exp)
                if meta is None:
                    continue
                p = f"{stem}/{exp}/06_reproducibility/final"
                _add_row(rows, "", exp, meta["assay"], meta["target"],
                         meta["genome"], "cuttag_macs3_narrow_consensus_final_peak",
                         "cp",
                         f"{p}/{exp}.cuttag.macs3.narrow.replicate_validated.consensus.narrowPeak")

        # cuttag broad: consensus final when broad IDR not enabled
        if not e["cuttag_broad_idr"]:
            for exp in sorted(e["consensus_experiments"].get(
                    ("cuttag", "broad"), [])):
                meta = _exp_meta(exp)
                if meta is None:
                    continue
                p = f"{stem}/{exp}/06_reproducibility/final"
                _add_row(rows, "", exp, meta["assay"], meta["target"],
                         meta["genome"], "cuttag_macs3_broad_consensus_final_peak",
                         "cp",
                         f"{p}/{exp}.cuttag.macs3.broad.replicate_validated.consensus.broadPeak")

        # SEACR consensus final
        if e["seacr_enabled"]:
            for exp in sorted(e["seacr_consensus_experiments"]):
                meta = _exp_meta(exp)
                if meta is None:
                    continue
                mode = e["seacr_mode"]
                p = f"{stem}/{exp}/06_reproducibility/final"
                _add_row(rows, "", exp, meta["assay"], meta["target"],
                         meta["genome"], "cuttag_seacr_consensus_final_peak",
                         "cp",
                         f"{p}/{exp}.cuttag.seacr.{mode}.replicate_validated.consensus.bed")

    return rows


def _build_project_rows(outdir, multiqc_enabled, has_peak_samples):
    rows = []
    if has_peak_samples:
        _add_row(rows, "", "", "", "", "",
                 "stage3_qc_summary", "aggregate_qc_summary",
                 _resolve_path(outdir, "multiqc", "stage3_qc_summary.tsv"))
    else:
        _add_row(rows, "", "", "", "", "",
                 "stage3_qc_summary", "aggregate_qc_summary",
                 "", check_exists=False)
    if multiqc_enabled:
        _add_row(rows, "", "", "", "", "",
                 "multiqc_report", "multiqc",
                 _resolve_path(outdir, "multiqc", "multiqc_report.html"))
    else:
        _add_row(rows, "", "", "", "", "",
                 "multiqc_report", "multiqc", "", check_exists=False)
    return rows


def build_manifest_rows(config, samples_path=None):
    """Build manifest row dicts from a validated config dict.

    Args:
        config: validated config dict.
        samples_path: optional override for samples TSV path.

    Returns:
        tuple (rows, missing_count, not_applicable_count)
    """
    outdir = config.get("outdir", "results")
    qc = config.get("qc", {})
    signal_tracks = qc.get("signal_tracks", True)
    stage4b = config.get("stage4b", True)
    stage5 = config.get("stage5", False)
    multiqc_enabled = config.get("multiqc", True)
    genomic_resources = config.get("genome_resources", {})

    samples_path = samples_path or config.get("samples", "config/samples.tsv")
    samples = load_and_validate_samples(
        samples_path,
        use_control=config.get("use_control", False),
        stage5_enabled=stage5,
    )
    treatment_samples = [s for s in samples if s.get("role") == "treatment"]
    has_peak_samples = any(_is_mnase(s) is False for s in treatment_samples)

    rows = []
    rows.extend(_build_sample_rows(
        treatment_samples, outdir, signal_tracks, genomic_resources))
    rows.extend(_build_experiment_rows(
        treatment_samples, outdir, signal_tracks, genomic_resources, stage4b))
    rows.extend(_build_idr_rows(treatment_samples, outdir, stage5))
    rows.extend(_build_reproducibility_rows(treatment_samples, config, outdir))
    rows.extend(_build_project_rows(outdir, multiqc_enabled, has_peak_samples))

    missing_count = sum(1 for r in rows if r["status"] == "missing")
    na_count = sum(1 for r in rows if r["status"] == "not_applicable")
    return rows, missing_count, na_count


def write_manifest_tsv(rows, output_path):
    """Write manifest rows to a TSV file."""
    os.makedirs(os.path.dirname(output_path) or ".", exist_ok=True)
    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh, fieldnames=_MANIFEST_COLUMNS, delimiter="\t",
            lineterminator="\n", extrasaction="ignore",
        )
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def main():
    # Keep argparse import local to main() so importing this module does not
    # depend on CLI machinery for library-style callers.
    import argparse
    import sys

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
        with open(args.config) as fh:
            raw_cfg = yaml.safe_load(fh)
        cfg = validate_config(raw_cfg)
    else:
        print("ERROR: --config or --config-json is required", file=sys.stderr)
        sys.exit(1)

    rows, missing_count, na_count = build_manifest_rows(cfg)
    write_manifest_tsv(rows, args.output)

    print(f"Manifest written: {args.output}")
    print(f"  {len(rows)} rows: "
          f"{len(rows) - missing_count - na_count} present, "
          f"{missing_count} missing, "
          f"{na_count} not_applicable")

    if args.strict and missing_count > 0:
        print(f"ERROR: --strict mode, {missing_count} missing output(s)",
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
