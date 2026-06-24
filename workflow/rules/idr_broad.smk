# idr_broad.smk — Stage 65 broad-peak IDR (experimental opt-in)
# ===============================================================
# Shared rule file for chipseq broad and cuttag broad IDR.
# All helpers use _broad_idr_ prefix. Rules use broad_idr_ prefix.
#
# Supported modes:
#   reproducibility.idr.chipseq_broad_experimental
#   reproducibility.idr.cuttag_broad_experimental
#
# Rules:
#   broad_idr_macs3_biorep       — IDR-ready MACS3 on biorep BAM
#   broad_idr_true_replicates    — true-replicate IDR
#   broad_idr_split_pseudoreps   — BAM pseudorep split
#   broad_idr_macs3_pseudorep    — MACS3 on pseudorep BAMs
#   broad_idr_self_pseudoreps    — self-pseudorep IDR
#   broad_idr_pooled_pseudoreps  — pooled pseudorep IDR
#   broad_idr_summary            — reproducibility QC + final peak


# ---------------------------------------------------------------------------
# Helper: MACS3 args for broad IDR-ready calls
# ---------------------------------------------------------------------------

def _broad_idr_macs3_args(wildcards):
    """Return MACS3 args for broad-peak IDR calls."""
    return idr_macs3_args(wildcards.experiment, wildcards.assay, "broad")


# ---------------------------------------------------------------------------
# Helper: inputs for broad IDR biorep MACS3
# ---------------------------------------------------------------------------

def _broad_idr_biorep_inputs(wildcards):
    """Return inputs for broad IDR per-biorep MACS3."""
    return idr_biorep_peaks_inputs(
        wildcards.experiment, int(wildcards.bio_rep)
    )


def _broad_idr_assay(experiment):
    """Return 'chipseq' or 'cuttag' for a broad IDR experiment."""
    if experiment in BROAD_CHIPSEQ_IDR_EXPERIMENTS:
        return "chipseq"
    if experiment in BROAD_CUTTAG_IDR_EXPERIMENTS:
        return "cuttag"
    raise ValueError(
        f"Experiment {experiment} not in any broad IDR experiment list"
    )


# ============================================================================
# 1. Per-biorep IDR-ready MACS3 (broad)
# ============================================================================

rule broad_idr_macs3_biorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
        f"idr_peaks/{{experiment}}_broad_{{assay}}_biorep{{bio_rep}}"
        f"_idr.broadPeak",
    input:
        lambda wc: _broad_idr_biorep_inputs(wc),
    params:
        macs3_args = lambda wc: _broad_idr_macs3_args(wc),
        sample     = lambda wc: (
            f"{wc.experiment}_broad_{wc.assay}_biorep{wc.bio_rep}_idr"
        ),
    wildcard_constraints:
        assay   = r"chipseq|cuttag",
        bio_rep = r"\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_broad_{{assay}}_biorep{{bio_rep}}_idr.macs3.log",
    threads: THREADS,
    conda:
        "../envs/macs3.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"

        set -- {input:q}
        TREATMENT="$1"
        if [[ $# -ge 3 ]]; then
            macs3 callpeak \
                -t "$TREATMENT" \
                -c "$3" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} \
                2>&1 | tee {log:q}
        else
            macs3 callpeak \
                -t "$TREATMENT" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} \
                2>&1 | tee {log:q}
        fi

        # MACS3 writes {name}_peaks.broadPeak. Move to canonical path.
        MACS3_OUT="$(dirname {output:q})/{params.sample}_peaks.broadPeak"
        if [[ -f "$MACS3_OUT" ]]; then
            mv "$MACS3_OUT" {output:q}
        else
            echo "ERROR: Expected MACS3 output not found: $MACS3_OUT" >&2
            exit 1
        fi
        """


# ============================================================================
# 2. True-replicate IDR (broadPeak)
# ============================================================================

rule broad_idr_true_replicates:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/true_replicates/{{experiment}}_broad_{{assay}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/true_replicates/"
                  f"{{experiment}}_broad_{{assay}}_idr.thresholded.broadPeak",
    input:
        peaks1 = lambda wc: idr_repro_peak_input(wc.experiment, 0, wc.assay, "broadPeak"),
        peaks2 = lambda wc: idr_repro_peak_input(wc.experiment, 1, wc.assay, "broadPeak"),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"true_replicates/{wc.experiment}_broad_{wc.assay}_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_idr_true.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_idr_true_thr.log"
        ),
    wildcard_constraints:
        assay = r"chipseq|cuttag",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.idr_out:q})" "$(dirname {output.thr_out:q})"

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type broadPeak \
            --rank {params.rank:q} \
            --output-file {params.output_prefix:q}.txt \
            --plot \
            2>&1 | tee {params.raw_log_file:q}

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type broadPeak \
            --rank {params.rank:q} \
            --idr-threshold {params.threshold} \
            --output-file {params.output_prefix:q}.thresholded.broadPeak \
            2>&1 | tee {params.thr_log_file:q}
        """


# ============================================================================
# 3. BAM pseudorep split
# ============================================================================

rule broad_idr_split_pseudoreps:
    output:
        pr1 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_broad_{{assay}}_{{source}}.pr1.bam",
        p1i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_broad_{{assay}}_{{source}}.pr1.bam.bai",
        pr2 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_broad_{{assay}}_{{source}}.pr2.bam",
        p2i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_broad_{{assay}}_{{source}}.pr2.bam.bai",
    input:
        lambda wc: _broad_split_input(wildcards=wc),
    wildcard_constraints:
        assay  = r"chipseq|cuttag",
        source = r"pooled|biorep\d+",
    threads: THREADS,
    conda:
        "../envs/samtools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.pr1:q})"
        python3 scripts/split_pseudoreps.py \
            --input {input:q} \
            --out1 {output.pr1:q} \
            --out2 {output.pr2:q} \
            --seed {IDR_SEED} \
            --threads {threads}
        samtools index {output.pr1:q}
        samtools index {output.pr2:q}
        """


def _broad_split_input(wildcards):
    """Return the BAM to split for broad pseudoreps."""
    return idr_split_input(wildcards.experiment, wildcards.source)


# ============================================================================
# 4. MACS3 on pseudorep BAMs (broad)
# ============================================================================

rule broad_idr_macs3_pseudorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
        f"idr_peaks/{{experiment}}_broad_{{assay}}_{{source}}_pr{{pr}}"
        f"_idr.broadPeak",
    input:
        lambda wc: _broad_idr_pseudorep_inputs(wildcards=wc),
    params:
        macs3_args = lambda wc: _broad_idr_macs3_args(wc),
        sample     = lambda wc: (
            f"{wc.experiment}_broad_{wc.assay}_{wc.source}_pr{wc.pr}_idr"
        ),
    wildcard_constraints:
        assay  = r"chipseq|cuttag",
        source = r"pooled|biorep\d+",
        pr     = r"[12]",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_broad_{{assay}}_{{source}}_pr{{pr}}_idr.macs3.log",
    threads: THREADS,
    conda:
        "../envs/macs3.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"

        set -- {input:q}
        TREATMENT="$1"
        CONTROL=""
        if [[ $# -ge 3 ]]; then
            CONTROL="$3"
        fi

        if [[ -n "$CONTROL" ]]; then
            macs3 callpeak \
                -t "$TREATMENT" \
                -c "$CONTROL" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} \
                2>&1 | tee {log:q}
        else
            macs3 callpeak \
                -t "$TREATMENT" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} \
                2>&1 | tee {log:q}
        fi

        MACS3_OUT="$(dirname {output:q})/{params.sample}_peaks.broadPeak"
        if [[ -f "$MACS3_OUT" ]]; then
            mv "$MACS3_OUT" {output:q}
        else
            echo "ERROR: Expected MACS3 output not found: $MACS3_OUT" >&2
            exit 1
        fi
        """


def _broad_idr_pseudorep_inputs(wildcards):
    """Return inputs for broad IDR pseudorep MACS3."""
    return idr_pseudorep_peaks_inputs(
        wildcards.experiment,
        wildcards.source,
        wildcards.pr,
        source_prefix=f"broad_{wildcards.assay}_",
    )


# ============================================================================
# 5. Self-pseudorep IDR
# ============================================================================

rule broad_idr_self_pseudoreps:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/self_pseudoreplicates/"
                  f"{{experiment}}_broad_{{assay}}_biorep{{bio_rep}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/self_pseudoreplicates/"
                  f"{{experiment}}_broad_{{assay}}_biorep{{bio_rep}}_idr."
                  f"thresholded.broadPeak",
    input:
        peaks1 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"idr_peaks/{wc.experiment}_broad_{wc.assay}_"
            f"biorep{wc.bio_rep}_pr1_idr.broadPeak"
        ),
        peaks2 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"idr_peaks/{wc.experiment}_broad_{wc.assay}_"
            f"biorep{wc.bio_rep}_pr2_idr.broadPeak"
        ),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"self_pseudoreplicates/"
            f"{wc.experiment}_broad_{wc.assay}_biorep{wc.bio_rep}_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_"
            f"biorep{wc.bio_rep}_idr_self.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_"
            f"biorep{wc.bio_rep}_idr_self_thr.log"
        ),
    wildcard_constraints:
        assay   = r"chipseq|cuttag",
        bio_rep = r"\d+",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.idr_out:q})" "$(dirname {output.thr_out:q})"

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type broadPeak \
            --rank {params.rank:q} \
            --output-file {output.idr_out:q} \
            --plot \
            2>&1 | tee {params.raw_log_file:q}

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type broadPeak \
            --rank {params.rank:q} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            2>&1 | tee {params.thr_log_file:q}
        """


# ============================================================================
# 6. Pooled-pseudorep IDR
# ============================================================================

rule broad_idr_pooled_pseudoreps:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/pooled_pseudoreplicates/"
                  f"{{experiment}}_broad_{{assay}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/pooled_pseudoreplicates/"
                  f"{{experiment}}_broad_{{assay}}_idr.thresholded.broadPeak",
    input:
        peaks1 = lambda wc: idr_pooled_peak_input(wc.experiment, 1, wc.assay, "broadPeak"),
        peaks2 = lambda wc: idr_pooled_peak_input(wc.experiment, 2, wc.assay, "broadPeak"),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"pooled_pseudoreplicates/"
            f"{wc.experiment}_broad_{wc.assay}_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_idr_pooled.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_idr_pooled_thr.log"
        ),
    wildcard_constraints:
        assay = r"chipseq|cuttag",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.idr_out:q})"

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type broadPeak \
            --rank {params.rank:q} \
            --output-file {output.idr_out:q} \
            --plot \
            2>&1 | tee {params.raw_log_file:q}

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type broadPeak \
            --rank {params.rank:q} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            2>&1 | tee {params.thr_log_file:q}
        """


# ============================================================================
# 7. Reproducibility summary + final peak
# ============================================================================
#
# The summary path is intentionally generic:
#   final/reproducibility_summary.tsv
#
# Snakemake requires every output in one rule to contain the same wildcard set,
# so broad IDR uses two assay-specific summary rules inside this shared file.
# This keeps final peak paths assay-specific without making summary paths
# ambiguous across ATAC, CUT&Tag narrow, and broad IDR rules.

rule broad_idr_chipseq_summary:
    output:
        summary    = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/reproducibility_summary.tsv",
        final_peak = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/{{experiment}}.chipseq.macs3.broad."
                     f"replicate_validated.idr.broadPeak",
    input:
        true_thresh = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                      f"idr/true_replicates/"
                      f"{{experiment}}_broad_chipseq_idr.thresholded.broadPeak",
        pool_thresh = lambda wc: idr_pooled_thresh_path(wc.experiment, "chipseq", "broadPeak"),
        self1_thresh = lambda wc: idr_self_thresh_path(wc.experiment, 0, "chipseq", "broadPeak"),
        self2_thresh = lambda wc: idr_self_thresh_path(wc.experiment, 1, "chipseq", "broadPeak"),
    params:
        experiment    = lambda wc: wc.experiment,
        assay         = "chipseq",
        caller        = "macs3",
        peak_mode     = "broad",
        bio_rep_a     = lambda wc: idr_biorep_labels(wc.experiment)[0],
        bio_rep_b     = lambda wc: idr_biorep_labels(wc.experiment)[1],
        final_method  = "idr",
        final_output  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/final/"
            f"{wc.experiment}.chipseq.macs3.broad."
            f"replicate_validated.idr.broadPeak"
        ),
    wildcard_constraints:
        experiment = wildcard_choices(BROAD_CHIPSEQ_IDR_EXPERIMENTS),
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.summary:q})" "$(dirname {output.final_peak:q})"

        python3 scripts/broad_idr_summary.py \
            --true-peaks {input.true_thresh:q} \
            --pooled-peaks {input.pool_thresh:q} \
            --self1-peaks {input.self1_thresh:q} \
            --self2-peaks {input.self2_thresh:q} \
            --experiment {params.experiment:q} \
            --assay chipseq \
            --caller macs3 \
            --peak-mode broad \
            --bio-rep-a {params.bio_rep_a} \
            --bio-rep-b {params.bio_rep_b} \
            --final-method idr \
            --final-output {params.final_output:q} \
            --output-tsv {output.summary:q} \
            --output-peak {output.final_peak:q}
        """


rule broad_idr_cuttag_summary:
    output:
        summary    = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/reproducibility_summary.tsv",
        final_peak = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/{{experiment}}.cuttag.macs3.broad."
                     f"replicate_validated.idr.broadPeak",
    input:
        true_thresh = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                      f"idr/true_replicates/"
                      f"{{experiment}}_broad_cuttag_idr.thresholded.broadPeak",
        pool_thresh = lambda wc: idr_pooled_thresh_path(wc.experiment, "cuttag", "broadPeak"),
        self1_thresh = lambda wc: idr_self_thresh_path(wc.experiment, 0, "cuttag", "broadPeak"),
        self2_thresh = lambda wc: idr_self_thresh_path(wc.experiment, 1, "cuttag", "broadPeak"),
    params:
        experiment    = lambda wc: wc.experiment,
        assay         = "cuttag",
        caller        = "macs3",
        peak_mode     = "broad",
        bio_rep_a     = lambda wc: idr_biorep_labels(wc.experiment)[0],
        bio_rep_b     = lambda wc: idr_biorep_labels(wc.experiment)[1],
        final_method  = "idr",
        final_output  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/final/"
            f"{wc.experiment}.cuttag.macs3.broad."
            f"replicate_validated.idr.broadPeak"
        ),
    wildcard_constraints:
        experiment = wildcard_choices(BROAD_CUTTAG_IDR_EXPERIMENTS),
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.summary:q})" "$(dirname {output.final_peak:q})"

        python3 scripts/broad_idr_summary.py \
            --true-peaks {input.true_thresh:q} \
            --pooled-peaks {input.pool_thresh:q} \
            --self1-peaks {input.self1_thresh:q} \
            --self2-peaks {input.self2_thresh:q} \
            --experiment {params.experiment:q} \
            --assay cuttag \
            --caller macs3 \
            --peak-mode broad \
            --bio-rep-a {params.bio_rep_a} \
            --bio-rep-b {params.bio_rep_b} \
            --final-method idr \
            --final-output {params.final_output:q} \
            --output-tsv {output.summary:q} \
            --output-peak {output.final_peak:q}
        """
