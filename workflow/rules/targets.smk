
# ---------------------------------------------------------------------------
# 7. Target builder helpers — return flat lists for rule all input
# ---------------------------------------------------------------------------


def _base_targets():
    """Core pipeline outputs: done markers, CPM BigWig, peaks, MultiQC."""
    targets = []
    targets += expand(
        "{outdir}/{sample}/logs/{sample}.pipeline.done",
        outdir=OUTDIR,
        sample=TREATMENT_SAMPLE_IDS,
    )
    targets += expand(
        "{outdir}/{sample}/03_bigwig/{sample}.CPM.bw",
        outdir=OUTDIR,
        sample=ACTIVE_SAMPLE_IDS,
    )
    targets += expand(
        "{outdir}/{sample}/04_peaks/{sample}",
        outdir=OUTDIR,
        sample=PEAK_SAMPLE_IDS,
    )
    if MULTIQC:
        targets += [f"{OUTDIR}/multiqc/multiqc_report.html"]
    return targets


def _manifest_dependency_targets():
    """Return all targets the manifest must wait for (all rule-all minus manifest)."""
    targets = []
    targets += _base_targets()
    targets += _blacklist_targets()
    targets += _single_sample_qc_targets()
    targets += _signal_targets()
    targets += _cuttag_targets()
    targets += _mnase_targets()
    targets += _advanced_qc_targets()
    targets += _tss_targets()
    targets += _replicate_targets()
    targets += _idr_targets()
    targets += _atac_idr_targets()
    # Remove the manifest itself (circular dependency)
    manifest_path = f"{OUTDIR}/multiqc/result_manifest.tsv"
    targets = [t for t in targets if t != manifest_path]
    return targets


import json as _json
_MANIFEST_CONFIG_JSON = _json.dumps(VALIDATED_CONFIG)


def _blacklist_targets():
    """Blacklist-filtered BAM + peaks (resource-gated per genome)."""
    targets = []
    if QC_CONFIG.get("blacklist_filter", True) and BLACKLIST_SAMPLES:
        targets += expand(
            "{outdir}/{sample}/02_align/{sample}.blacklist_filtered.bam",
            outdir=OUTDIR,
            sample=BLACKLIST_SAMPLES,
        )
        targets += expand(
            "{outdir}/{sample}/02_align/{sample}.blacklist_filtered.bam.bai",
            outdir=OUTDIR,
            sample=BLACKLIST_SAMPLES,
        )
        narrow_bl = [s for s in BLACKLIST_SAMPLES if SAMPLE_MAP[s]["peak_mode"] == "narrow"]
        broad_bl = [s for s in BLACKLIST_SAMPLES if SAMPLE_MAP[s]["peak_mode"] == "broad"]
        if narrow_bl:
            targets += expand(
                "{outdir}/{sample}/04_peaks/{sample}_blacklist_filtered/{sample}_peaks.blacklist_filtered.narrowPeak",
                outdir=OUTDIR,
                sample=narrow_bl,
            )
        if broad_bl:
            targets += expand(
                "{outdir}/{sample}/04_peaks/{sample}_blacklist_filtered/{sample}_peaks.blacklist_filtered.broadPeak",
                outdir=OUTDIR,
                sample=broad_bl,
            )
    return targets


def _single_sample_qc_targets():
    """Per-sample QC: peak counts, library complexity, NRF/PBC, FRiP, summaries.

    Peak-dependent QC (peak_counts, frip, qc_summary) is gated on
    PEAK_SAMPLE_IDS to exclude MNase samples. Library complexity and
    NRF/PBC are BAM-derived and apply to all treatment samples.
    """
    targets = []
    if QC_CONFIG.get("summary", True):
        targets += expand(
            "{outdir}/{sample}/01_qc/{sample}.peak_counts.tsv",
            outdir=OUTDIR,
            sample=PEAK_SAMPLE_IDS,
        )
    if QC_CONFIG.get("library_complexity", True) or QC_CONFIG.get("summary", True):
        targets += expand(
            "{outdir}/{sample}/01_qc/{sample}.library_complexity.tsv",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        )
    if QC_CONFIG.get("nrf_pbc", True) or QC_CONFIG.get("summary", True):
        targets += expand(
            "{outdir}/{sample}/01_qc/{sample}.nrf_pbc.tsv",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        )
    if QC_CONFIG.get("frip", True) or QC_CONFIG.get("summary", True):
        targets += expand(
            "{outdir}/{sample}/01_qc/{sample}.frip.tsv",
            outdir=OUTDIR,
            sample=PEAK_SAMPLE_IDS,
        )
    if QC_CONFIG.get("summary", True):
        targets += expand(
            "{outdir}/{sample}/01_qc/{sample}.qc_summary.tsv",
            outdir=OUTDIR,
            sample=PEAK_SAMPLE_IDS,
        )
    if QC_CONFIG.get("summary", True) and PEAK_SAMPLE_IDS:
        targets += [f"{OUTDIR}/multiqc/stage3_qc_summary.tsv"]
    return targets


def _signal_targets():
    """MACS3 FE/ppois bedGraph + BigWig for single-sample and pooled experiments.

    Gated on PEAK_SAMPLE_IDS / PEAK_MULTI_BIOREP_EXPERIMENTS to exclude MNase.
    """
    targets = []
    if QC_CONFIG.get("signal_tracks", True):
        # Single-sample bedGraph
        targets += expand(
            "{outdir}/{sample}/03_signal/{sample}.FE.bdg",
            outdir=OUTDIR,
            sample=PEAK_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/03_signal/{sample}.ppois.bdg",
            outdir=OUTDIR,
            sample=PEAK_SAMPLE_IDS,
        )
    # Stage 22: single-sample BigWig (gated on signal_tracks + chrom_sizes)
    if SIGNAL_BW_SAMPLE_IDS:
        targets += expand(
            "{outdir}/{sample}/03_signal/{sample}.FE.bw",
            outdir=OUTDIR,
            sample=SIGNAL_BW_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/03_signal/{sample}.ppois.bw",
            outdir=OUTDIR,
            sample=SIGNAL_BW_SAMPLE_IDS,
        )
    # Stage 6a: pooled bedGraph
    if STAGE4B and PEAK_MULTI_BIOREP_EXPERIMENTS and QC_CONFIG.get("signal_tracks", True):
        targets += expand(
            "{outdir}/experiments/{experiment}/03_signal/{experiment}.pooled.FE.bdg",
            outdir=OUTDIR,
            experiment=PEAK_MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/03_signal/{experiment}.pooled.ppois.bdg",
            outdir=OUTDIR,
            experiment=PEAK_MULTI_BIOREP_EXPERIMENTS,
        )
    # Stage 22: pooled BigWig (gated on signal_tracks + chrom_sizes)
    if STAGE4B and SIGNAL_BW_EXPERIMENTS:
        targets += expand(
            "{outdir}/experiments/{experiment}/03_signal/{experiment}.pooled.FE.bw",
            outdir=OUTDIR,
            experiment=SIGNAL_BW_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/03_signal/{experiment}.pooled.ppois.bw",
            outdir=OUTDIR,
            experiment=SIGNAL_BW_EXPERIMENTS,
        )
    return targets


def _cuttag_targets():
    """SEACR sidecar + CUT&Tag fragment-size QC."""
    targets = []
    if SEACR_ENABLED and SEACR_SAMPLE_IDS:
        targets += expand(
            "{outdir}/{sample}/04_peaks_seacr/{sample}.bedgraph",
            outdir=OUTDIR,
            sample=SEACR_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/04_peaks_seacr/{sample}/{sample}.seacr.{mode}.bed",
            outdir=OUTDIR,
            sample=SEACR_SAMPLE_IDS,
            mode=[SEACR_MODE],
        )
    if QC_CONFIG.get("cuttag_fragment_size", True) and CUTTAG_ACTIVE_SAMPLE_IDS:
        targets += expand(
            "{outdir}/{sample}/01_qc/{sample}.cuttag_fragment_size.tsv",
            outdir=OUTDIR,
            sample=CUTTAG_ACTIVE_SAMPLE_IDS,
        )
    return targets


def _advanced_qc_targets():
    """Opt-in QC: cross-correlation, preseq, Picard metrics."""
    targets = []
    if QC_CONFIG.get("cross_correlation", False) and TREATMENT_SAMPLE_IDS:
        targets += expand(
            "{outdir}/{sample}/05_qc/cross_correlation/{sample}.cc.qc",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/05_qc/cross_correlation/{sample}.cc.plot.pdf",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        )
        targets += [f"{OUTDIR}/multiqc/cross_correlation_summary.tsv"]
    if QC_CONFIG.get("preseq_complexity", False) and TREATMENT_SAMPLE_IDS:
        targets += expand(
            "{outdir}/{sample}/05_qc/preseq/{sample}.preseq.txt",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        )
    if QC_CONFIG.get("picard_metrics", False) and TREATMENT_SAMPLE_IDS:
        targets += expand(
            "{outdir}/{sample}/05_qc/picard/{sample}.alignment_summary_metrics",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/05_qc/picard/{sample}.insert_size_metrics",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/05_qc/picard/{sample}.quality_distribution_metrics",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/05_qc/picard/{sample}.insert_size_histogram.pdf",
            outdir=OUTDIR,
            sample=TREATMENT_SAMPLE_IDS,
        )
    return targets


def _tss_targets():
    """TSS enrichment profiles (requires genome_resources.<genome>.gtf)."""
    targets = []
    if QC_CONFIG.get("tss_enrichment", False) and TSS_GENOMES:
        targets += expand(
            "{outdir}/reference/{genome}.tss.bed",
            outdir=OUTDIR,
            genome=TSS_GENOMES,
        )
    if QC_CONFIG.get("tss_enrichment", False) and TSS_SAMPLE_IDS:
        targets += expand(
            "{outdir}/{sample}/05_qc/tss/{sample}.tss_matrix.gz",
            outdir=OUTDIR,
            sample=TSS_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/05_qc/tss/{sample}.tss_profile.tsv",
            outdir=OUTDIR,
            sample=TSS_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/05_qc/tss/{sample}.tss_profile.pdf",
            outdir=OUTDIR,
            sample=TSS_SAMPLE_IDS,
        )
    return targets


def _replicate_targets():
    """Stage 4b/6b: pooled BAMs, biorep BAMs, pooled peaks, pooled QC summary."""
    targets = []
    if STAGE4B and MULTI_BIOREP_EXPERIMENTS:
        targets += expand(
            "{outdir}/experiments/{experiment}/02_align/{experiment}.pooled.final.bam",
            outdir=OUTDIR,
            experiment=MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/02_align/{experiment}.pooled.final.bam.bai",
            outdir=OUTDIR,
            experiment=MULTI_BIOREP_EXPERIMENTS,
        )
    if STAGE4B and _EXP_LIST:
        targets += expand(
            "{outdir}/experiments/{experiment}/02_align/biorep{bio_rep}.final.bam",
            zip,
            outdir=OUTDIR,
            experiment=_EXP_LIST, bio_rep=_BR_LIST,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/02_align/biorep{bio_rep}.final.bam.bai",
            zip,
            outdir=OUTDIR,
            experiment=_EXP_LIST, bio_rep=_BR_LIST,
        )
    if STAGE4B and POOLED_CONTROL_EXPERIMENTS:
        targets += expand(
            "{outdir}/experiments/{experiment}/02_align/{experiment}.pooled.control.final.bam",
            outdir=OUTDIR,
            experiment=POOLED_CONTROL_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/02_align/{experiment}.pooled.control.final.bam.bai",
            outdir=OUTDIR,
            experiment=POOLED_CONTROL_EXPERIMENTS,
        )
    if STAGE4B and PEAK_MULTI_BIOREP_EXPERIMENTS:
        targets += expand(
            "{outdir}/experiments/{experiment}/04_peaks/pooled/{experiment}_pooled_peaks",
            outdir=OUTDIR,
            experiment=PEAK_MULTI_BIOREP_EXPERIMENTS,
        )
    if STAGE4B and PEAK_MULTI_BIOREP_EXPERIMENTS:
        targets += expand(
            "{outdir}/experiments/{experiment}/01_qc/{experiment}.pooled_qc_summary.tsv",
            outdir=OUTDIR,
            experiment=PEAK_MULTI_BIOREP_EXPERIMENTS,
        )
    return targets


def _mnase_targets():
    """Stage 39-40: MNase nucleosome-centric outputs (sample and pooled)."""
    targets = []
    # Sample-level MNase outputs
    if MNASE_SAMPLE_IDS:
        # Stage 39: mono BAM, dyad BW, mono occupancy BW
        targets += expand(
            mnase_fragment_bam("{sample}", "mono"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_fragment_bai("{sample}", "mono"),
            sample=MNASE_SAMPLE_IDS,
        )
        # Stage 40: sub-nucleosome and di-nucleosome BAMs
        targets += expand(
            mnase_fragment_bam("{sample}", "sub"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_fragment_bai("{sample}", "sub"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_fragment_bam("{sample}", "di"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_fragment_bai("{sample}", "di"),
            sample=MNASE_SAMPLE_IDS,
        )
        # Stage 40: MNase QC summary
        targets += expand(
            mnase_qc_summary_tsv("{sample}"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_signal_bw("{sample}", "dyad"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_signal_bw("{sample}", "mono"),
            sample=MNASE_SAMPLE_IDS,
        )
    # Pooled MNase outputs (>=2 biorep MNase experiments)
    if STAGE4B and MNASE_MULTI_BIOREP_EXPERIMENTS:
        targets += expand(
            mnase_pooled_fragment_bam("{experiment}", "mono"),
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            mnase_pooled_fragment_bai("{experiment}", "mono"),
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            mnase_pooled_signal_bw("{experiment}", "dyad"),
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            mnase_pooled_signal_bw("{experiment}", "mono"),
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
    return targets


def _idr_targets():
    """Stage 5a/5b: IDR-ready peaks, true-replicate IDR, pseudoreplicate IDR, final outputs."""
    targets = []
    if STAGE5 and IDR_EXPERIMENTS:
        targets += expand(
            "{outdir}/experiments/{experiment}/04_peaks/idr/{experiment}_biorep{bio_rep}_idr_peaks.narrowPeak",
            zip,
            outdir=OUTDIR,
            experiment=IDR_BIOREP_EXP_LIST, bio_rep=IDR_BIOREP_LIST,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_idr/true_replicates/idr.txt",
            outdir=OUTDIR,
            experiment=IDR_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_idr/true_replicates/idr.thresholded.narrowPeak",
            outdir=OUTDIR,
            experiment=IDR_EXPERIMENTS,
        )
    if STAGE5 and IDR_EXPERIMENTS:
        targets += expand(
            "{outdir}/experiments/{experiment}/04_peaks/idr/{experiment}_{source}_pr{pr}_idr_peaks.narrowPeak",
            zip,
            outdir=OUTDIR,
            experiment=IDR_PR_PEAK_EXP, source=IDR_PR_PEAK_SRC, pr=IDR_PR_PEAK_PR,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_idr/self_pseudoreplicates/biorep{bio_rep}.idr.txt",
            zip,
            outdir=OUTDIR,
            experiment=IDR_SELF_EXP, bio_rep=IDR_SELF_BR,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_idr/self_pseudoreplicates/biorep{bio_rep}.idr.thresholded.narrowPeak",
            zip,
            outdir=OUTDIR,
            experiment=IDR_SELF_EXP, bio_rep=IDR_SELF_BR,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_idr/pooled_pseudoreplicates/idr.txt",
            outdir=OUTDIR,
            experiment=IDR_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_idr/pooled_pseudoreplicates/idr.thresholded.narrowPeak",
            outdir=OUTDIR,
            experiment=IDR_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_idr/final/reproducibility_summary.tsv",
            outdir=OUTDIR,
            experiment=IDR_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_idr/final/conservative.narrowPeak",
            outdir=OUTDIR,
            experiment=IDR_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_idr/final/optimal.narrowPeak",
            outdir=OUTDIR,
            experiment=IDR_EXPERIMENTS,
        )
    return targets


def _atac_idr_targets():
    """Stage 55: ATAC narrow IDR — biorep peaks, true-rep IDR, pseudorep IDR, final outputs."""
    targets = []
    if ATAC_IDR_ENABLED and ATAC_IDR_EXPERIMENTS:
        # IDR-ready per-biorep MACS3 peaks
        targets += expand(
            "{outdir}/experiments/{experiment}/06_reproducibility/idr/"
            "idr_peaks/{experiment}_atac_biorep{bio_rep}_idr.narrowPeak",
            zip,
            outdir=OUTDIR,
            experiment=ATAC_IDR_BIOREP_EXP_LIST, bio_rep=ATAC_IDR_BIOREP_LIST,
        )
        # True-replicate IDR
        targets += expand(
            "{outdir}/experiments/{experiment}/06_reproducibility/idr/"
            "true_replicates/{experiment}_atac_idr.txt",
            outdir=OUTDIR,
            experiment=ATAC_IDR_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_reproducibility/idr/"
            "true_replicates/{experiment}_atac_idr.thresholded.narrowPeak",
            outdir=OUTDIR,
            experiment=ATAC_IDR_EXPERIMENTS,
        )
    if ATAC_IDR_ENABLED and ATAC_IDR_EXPERIMENTS:
        # IDR-ready pseudorep MACS3 peaks
        targets += expand(
            "{outdir}/experiments/{experiment}/06_reproducibility/idr/"
            "idr_peaks/{experiment}_atac_{source}_pr{pr}_idr.narrowPeak",
            zip,
            outdir=OUTDIR,
            experiment=ATAC_IDR_PR_PEAK_EXP,
            source=ATAC_IDR_PR_PEAK_SRC, pr=ATAC_IDR_PR_PEAK_PR,
        )
        # Self-pseudorep IDR
        targets += expand(
            "{outdir}/experiments/{experiment}/06_reproducibility/idr/"
            "self_pseudoreplicates/{experiment}_atac_biorep{bio_rep}_idr.txt",
            zip,
            outdir=OUTDIR,
            experiment=ATAC_IDR_SELF_EXP, bio_rep=ATAC_IDR_SELF_BR,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_reproducibility/idr/"
            "self_pseudoreplicates/{experiment}_atac_biorep{bio_rep}_idr.thresholded.narrowPeak",
            zip,
            outdir=OUTDIR,
            experiment=ATAC_IDR_SELF_EXP, bio_rep=ATAC_IDR_SELF_BR,
        )
        # Pooled-pseudorep IDR
        targets += expand(
            "{outdir}/experiments/{experiment}/06_reproducibility/idr/"
            "pooled_pseudoreplicates/{experiment}_atac_idr.txt",
            outdir=OUTDIR,
            experiment=ATAC_IDR_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_reproducibility/idr/"
            "pooled_pseudoreplicates/{experiment}_atac_idr.thresholded.narrowPeak",
            outdir=OUTDIR,
            experiment=ATAC_IDR_EXPERIMENTS,
        )
        # Final outputs under 06_reproducibility/final/
        targets += expand(
            "{outdir}/experiments/{experiment}/06_reproducibility/final/"
            "{experiment}.atac.macs3.narrow.replicate_validated.idr.narrowPeak",
            outdir=OUTDIR,
            experiment=ATAC_IDR_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/06_reproducibility/final/"
            "reproducibility_summary.tsv",
            outdir=OUTDIR,
            experiment=ATAC_IDR_EXPERIMENTS,
        )
    return targets
