# idr_reproducibility.smk — Unified ATAC/CUT&Tag/broad IDR reproducibility rules
# =============================================================================
# Replaces the assay-specific rule files idr_atac.smk, idr_cuttag.smk, and
# idr_broad.smk with a single parameterized file.
#
# The legacy Stage 5 ChIP-seq IDR rules in idr.smk are intentionally kept
# separate because they live in the 06_idr/ namespace and use a different
# filename convention.
#
# Supported profiles:
#   - ATAC narrow      (assay=atac)
#   - CUT&Tag narrow   (assay=cuttag)
#   - ChIP-seq broad   (assay=chipseq)
#   - CUT&Tag broad    (assay=cuttag)
#
# Narrow and broad modes are split into separate rules because the raw IDR
# text output does not contain a peak-mode suffix; mixing modes in one rule
# would violate Snakemake's requirement that all output files share the same
# wildcards.


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _idr_repro_biorep_inputs(wildcards):
    """Inputs for the per-biorep IDR-ready MACS3 call."""
    return idr_biorep_peaks_inputs(
        wildcards.experiment, int(wildcards.bio_rep)
    )


def _idr_repro_split_input(wildcards):
    """Input BAM for pseudoreplicate splitting."""
    return idr_split_input(wildcards.experiment, wildcards.source)


def _idr_repro_pseudorep_inputs(wildcards, source_prefix):
    """Inputs for the pseudorep IDR-ready MACS3 call."""
    return idr_pseudorep_peaks_inputs(
        wildcards.experiment,
        wildcards.source,
        wildcards.pr,
        source_prefix=source_prefix,
    )


# ===========================================================================
# Narrow mode (ATAC and CUT&Tag)
# ===========================================================================

# ---------------------------------------------------------------------------
# 1. Per-biorep IDR-ready MACS3 call
# ---------------------------------------------------------------------------

rule idr_macs3_biorep_narrow:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
        f"idr_peaks/{{experiment}}_{{assay}}_biorep{{bio_rep}}_idr.narrowPeak",
    input:
        lambda wc: _idr_repro_biorep_inputs(wc),
    params:
        macs3_args = lambda wc: idr_macs3_args(wc.experiment, wc.assay, "narrow"),
        sample     = lambda wc: (
            f"{wc.experiment}_{wc.assay}_biorep{wc.bio_rep}_idr"
        ),
    wildcard_constraints:
        assay   = r"atac|cuttag",
        bio_rep = r"\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_{{assay}}_biorep{{bio_rep}}_idr.narrow.macs3.log",
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

        MACS3_OUT="$(dirname {output:q})/{params.sample}_peaks.narrowPeak"
        if [[ -f "$MACS3_OUT" ]]; then
            mv "$MACS3_OUT" {output:q}
        elif [[ ! -f "{output}" ]]; then
            echo "ERROR: Expected MACS3 output not found: {output}" >&2
            exit 1
        fi
        """


# ---------------------------------------------------------------------------
# 2. True-replicate IDR
# ---------------------------------------------------------------------------

rule idr_true_replicates_narrow:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"true_replicates/{{experiment}}_{{assay}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"true_replicates/{{experiment}}_{{assay}}_idr."
                  f"thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: idr_repro_peak_input(
            wc.experiment, 0, wc.assay, "narrowPeak"
        ),
        peaks2 = lambda wc: idr_repro_peak_input(
            wc.experiment, 1, wc.assay, "narrowPeak"
        ),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"true_replicates/{wc.experiment}_{wc.assay}_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_{wc.assay}_idr_true.narrow.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_{wc.assay}_idr_true_thr.narrow.log"
        ),
    wildcard_constraints:
        assay = r"atac|cuttag",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_{{assay}}_idr_true_thr.narrow.log",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.idr_out:q})" \
                 "$(dirname {params.raw_log_file:q})" \
                 "$(dirname {params.thr_log_file:q})"

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank:q} \
            --output-file {output.idr_out:q} \
            --plot \
            2>&1 | tee {params.raw_log_file:q}

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank:q} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 3. BAM pseudoreplicate split
# ---------------------------------------------------------------------------

rule idr_split_pseudoreps_narrow:
    output:
        pr1 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_{{assay}}_{{source}}.pr1.bam",
        p1i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_{{assay}}_{{source}}.pr1.bam.bai",
        pr2 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_{{assay}}_{{source}}.pr2.bam",
        p2i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_{{assay}}_{{source}}.pr2.bam.bai",
    input:
        lambda wc: _idr_repro_split_input(wc),
    params:
        seed = IDR_SEED,
    wildcard_constraints:
        assay  = r"atac|cuttag",
        source = r"pooled|biorep\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_{{assay}}_{{source}}.narrow.split.log",
    threads: THREADS,
    conda:
        "../envs/samtools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.pr1:q})" "$(dirname {log:q})"
        python3 scripts/split_pseudoreps.py \
            --input {input:q} \
            --out1 {output.pr1:q} \
            --out2 {output.pr2:q} \
            --seed {params.seed} \
            --threads {threads} \
            2>&1 | tee {log:q}
        samtools index {output.pr1:q}
        samtools index {output.pr2:q}
        """


# ---------------------------------------------------------------------------
# 4. IDR-ready MACS3 on pseudoreplicate BAMs
# ---------------------------------------------------------------------------

rule idr_macs3_pseudorep_narrow:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
        f"idr_peaks/{{experiment}}_{{assay}}_{{source}}_pr{{pr}}_idr.narrowPeak",
    input:
        lambda wc: _idr_repro_pseudorep_inputs(wc, source_prefix=f"{wc.assay}_"),
    params:
        macs3_args = lambda wc: idr_macs3_args(wc.experiment, wc.assay, "narrow"),
        sample     = lambda wc: (
            f"{wc.experiment}_{wc.assay}_{wc.source}_pr{wc.pr}_idr"
        ),
    wildcard_constraints:
        assay  = r"atac|cuttag",
        source = r"pooled|biorep\d+",
        pr     = r"[12]",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_{{assay}}_{{source}}_pr{{pr}}_idr.narrow.macs3.log",
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

        MACS3_OUT="$(dirname {output:q})/{params.sample}_peaks.narrowPeak"
        if [[ -f "$MACS3_OUT" ]]; then
            mv "$MACS3_OUT" {output:q}
        elif [[ ! -f "{output}" ]]; then
            echo "ERROR: Expected MACS3 output not found: {output}" >&2
            exit 1
        fi
        """


# ---------------------------------------------------------------------------
# 5. Self-pseudoreplicate IDR
# ---------------------------------------------------------------------------

rule idr_self_pseudoreps_narrow:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"self_pseudoreplicates/"
                  f"{{experiment}}_{{assay}}_biorep{{bio_rep}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"self_pseudoreplicates/"
                  f"{{experiment}}_{{assay}}_biorep{{bio_rep}}_idr."
                  f"thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"idr_peaks/{wc.experiment}_{wc.assay}_"
            f"biorep{wc.bio_rep}_pr1_idr.narrowPeak"
        ),
        peaks2 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"idr_peaks/{wc.experiment}_{wc.assay}_"
            f"biorep{wc.bio_rep}_pr2_idr.narrowPeak"
        ),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"self_pseudoreplicates/{wc.experiment}_{wc.assay}_"
            f"biorep{wc.bio_rep}_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_{wc.assay}_"
            f"biorep{wc.bio_rep}_idr_self.narrow.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_{wc.assay}_"
            f"biorep{wc.bio_rep}_idr_self_thr.narrow.log"
        ),
    wildcard_constraints:
        assay   = r"atac|cuttag",
        bio_rep = r"\d+",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.idr_out:q})" \
                 "$(dirname {params.raw_log_file:q})" \
                 "$(dirname {params.thr_log_file:q})"

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank:q} \
            --output-file {output.idr_out:q} \
            --plot \
            2>&1 | tee {params.raw_log_file:q}

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank:q} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            2>&1 | tee {params.thr_log_file:q}
        """


# ---------------------------------------------------------------------------
# 6. Pooled-pseudoreplicate IDR
# ---------------------------------------------------------------------------

rule idr_pooled_pseudoreps_narrow:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"pooled_pseudoreplicates/{{experiment}}_{{assay}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"pooled_pseudoreplicates/{{experiment}}_{{assay}}_idr."
                  f"thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: idr_pooled_peak_input(
            wc.experiment, 1, wc.assay, "narrowPeak"
        ),
        peaks2 = lambda wc: idr_pooled_peak_input(
            wc.experiment, 2, wc.assay, "narrowPeak"
        ),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"pooled_pseudoreplicates/{wc.experiment}_{wc.assay}_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_{wc.assay}_idr_pooled.narrow.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_{wc.assay}_idr_pooled_thr.narrow.log"
        ),
    wildcard_constraints:
        assay = r"atac|cuttag",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.idr_out:q})"

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank:q} \
            --output-file {output.idr_out:q} \
            --plot \
            2>&1 | tee {params.raw_log_file:q}

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank:q} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            2>&1 | tee {params.thr_log_file:q}
        """


# ---------------------------------------------------------------------------
# 7. Reproducibility summary + final peak
# ---------------------------------------------------------------------------

rule idr_summary_atac_narrow:
    output:
        summary    = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/reproducibility_summary.tsv",
        final_peak = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/{{experiment}}.atac.macs3.narrow."
                     f"replicate_validated.idr.narrowPeak",
    input:
        true_thresh = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                      f"true_replicates/{{experiment}}_atac_idr.thresholded.narrowPeak",
        pool_thresh = lambda wc: idr_pooled_thresh_path(wc.experiment, "atac", "narrowPeak"),
        self1_thresh = lambda wc: idr_self_thresh_path(wc.experiment, 0, "atac", "narrowPeak"),
        self2_thresh = lambda wc: idr_self_thresh_path(wc.experiment, 1, "atac", "narrowPeak"),
    params:
        experiment    = lambda wc: wc.experiment,
        assay         = "atac",
        caller        = "macs3",
        peak_mode     = "narrow",
        bio_rep_a     = lambda wc: idr_biorep_labels(wc.experiment)[0],
        bio_rep_b     = lambda wc: idr_biorep_labels(wc.experiment)[1],
        final_method  = "idr",
        final_output  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/final/"
            f"{wc.experiment}.atac.macs3.narrow."
            f"replicate_validated.idr.narrowPeak"
        ),
    wildcard_constraints:
        experiment = wildcard_choices(ATAC_IDR_EXPERIMENTS),
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}.atac.macs3.narrow.idr.summary.log",
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.summary:q})" \
                 "$(dirname {output.final_peak:q})" \
                 "$(dirname {log:q})"

        python3 scripts/idr_reproducibility_summary.py \
            --true-peaks {input.true_thresh:q} \
            --pooled-peaks {input.pool_thresh:q} \
            --self1-peaks {input.self1_thresh:q} \
            --self2-peaks {input.self2_thresh:q} \
            --experiment {params.experiment:q} \
            --assay {params.assay} \
            --caller {params.caller} \
            --peak-mode {params.peak_mode} \
            --bio-rep-a {params.bio_rep_a} \
            --bio-rep-b {params.bio_rep_b} \
            --final-method {params.final_method} \
            --final-output {params.final_output:q} \
            --output-tsv {output.summary:q} \
            --output-peak {output.final_peak:q} \
            2>&1 | tee {log:q}
        """


rule idr_summary_cuttag_narrow:
    output:
        summary    = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/reproducibility_summary.tsv",
        final_peak = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/{{experiment}}.cuttag.macs3.narrow."
                     f"replicate_validated.idr.narrowPeak",
    input:
        true_thresh = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                      f"true_replicates/{{experiment}}_cuttag_idr.thresholded.narrowPeak",
        pool_thresh = lambda wc: idr_pooled_thresh_path(wc.experiment, "cuttag", "narrowPeak"),
        self1_thresh = lambda wc: idr_self_thresh_path(wc.experiment, 0, "cuttag", "narrowPeak"),
        self2_thresh = lambda wc: idr_self_thresh_path(wc.experiment, 1, "cuttag", "narrowPeak"),
    params:
        experiment    = lambda wc: wc.experiment,
        assay         = "cuttag",
        caller        = "macs3",
        peak_mode     = "narrow",
        bio_rep_a     = lambda wc: idr_biorep_labels(wc.experiment)[0],
        bio_rep_b     = lambda wc: idr_biorep_labels(wc.experiment)[1],
        final_method  = "idr",
        final_output  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/final/"
            f"{wc.experiment}.cuttag.macs3.narrow."
            f"replicate_validated.idr.narrowPeak"
        ),
    wildcard_constraints:
        experiment = wildcard_choices(CUTTAG_IDR_EXPERIMENTS),
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}.cuttag.macs3.narrow.idr.summary.log",
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.summary:q})" \
                 "$(dirname {output.final_peak:q})" \
                 "$(dirname {log:q})"

        python3 scripts/idr_reproducibility_summary.py \
            --true-peaks {input.true_thresh:q} \
            --pooled-peaks {input.pool_thresh:q} \
            --self1-peaks {input.self1_thresh:q} \
            --self2-peaks {input.self2_thresh:q} \
            --experiment {params.experiment:q} \
            --assay {params.assay} \
            --caller {params.caller} \
            --peak-mode {params.peak_mode} \
            --bio-rep-a {params.bio_rep_a} \
            --bio-rep-b {params.bio_rep_b} \
            --final-method {params.final_method} \
            --final-output {params.final_output:q} \
            --output-tsv {output.summary:q} \
            --output-peak {output.final_peak:q} \
            2>&1 | tee {log:q}
        """


# ===========================================================================
# Broad mode (ChIP-seq and CUT&Tag)
# ===========================================================================

# ---------------------------------------------------------------------------
# 1. Per-biorep IDR-ready MACS3 call
# ---------------------------------------------------------------------------

rule idr_macs3_biorep_broad:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
        f"idr_peaks/{{experiment}}_broad_{{assay}}_biorep{{bio_rep}}_idr.broadPeak",
    input:
        lambda wc: _idr_repro_biorep_inputs(wc),
    params:
        macs3_args = lambda wc: idr_macs3_args(wc.experiment, wc.assay, "broad"),
        sample     = lambda wc: (
            f"{wc.experiment}_broad_{wc.assay}_biorep{wc.bio_rep}_idr"
        ),
    wildcard_constraints:
        assay   = r"chipseq|cuttag",
        bio_rep = r"\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_broad_{{assay}}_biorep{{bio_rep}}_idr.broad.macs3.log",
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
        elif [[ ! -f "{output}" ]]; then
            echo "ERROR: Expected MACS3 output not found: {output}" >&2
            exit 1
        fi
        """


# ---------------------------------------------------------------------------
# 2. True-replicate IDR
# ---------------------------------------------------------------------------

rule idr_true_replicates_broad:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"true_replicates/{{experiment}}_broad_{{assay}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"true_replicates/{{experiment}}_broad_{{assay}}_idr."
                  f"thresholded.broadPeak",
    input:
        peaks1 = lambda wc: idr_repro_peak_input(
            wc.experiment, 0, wc.assay, "broadPeak"
        ),
        peaks2 = lambda wc: idr_repro_peak_input(
            wc.experiment, 1, wc.assay, "broadPeak"
        ),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"true_replicates/{wc.experiment}_broad_{wc.assay}_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_idr_true.broad.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_idr_true_thr.broad.log"
        ),
    wildcard_constraints:
        assay = r"chipseq|cuttag",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_broad_{{assay}}_idr_true_thr.broad.log",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.idr_out:q})" \
                 "$(dirname {params.raw_log_file:q})" \
                 "$(dirname {params.thr_log_file:q})"

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
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 3. BAM pseudoreplicate split
# ---------------------------------------------------------------------------

rule idr_split_pseudoreps_broad:
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
        lambda wc: _idr_repro_split_input(wc),
    params:
        seed = IDR_SEED,
    wildcard_constraints:
        assay  = r"chipseq|cuttag",
        source = r"pooled|biorep\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_broad_{{assay}}_{{source}}.broad.split.log",
    threads: THREADS,
    conda:
        "../envs/samtools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.pr1:q})" "$(dirname {log:q})"
        python3 scripts/split_pseudoreps.py \
            --input {input:q} \
            --out1 {output.pr1:q} \
            --out2 {output.pr2:q} \
            --seed {params.seed} \
            --threads {threads} \
            2>&1 | tee {log:q}
        samtools index {output.pr1:q}
        samtools index {output.pr2:q}
        """


# ---------------------------------------------------------------------------
# 4. IDR-ready MACS3 on pseudoreplicate BAMs
# ---------------------------------------------------------------------------

rule idr_macs3_pseudorep_broad:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
        f"idr_peaks/{{experiment}}_broad_{{assay}}_{{source}}_pr{{pr}}_idr.broadPeak",
    input:
        lambda wc: _idr_repro_pseudorep_inputs(
            wc, source_prefix=f"broad_{wc.assay}_"
        ),
    params:
        macs3_args = lambda wc: idr_macs3_args(wc.experiment, wc.assay, "broad"),
        sample     = lambda wc: (
            f"{wc.experiment}_broad_{wc.assay}_{wc.source}_pr{wc.pr}_idr"
        ),
    wildcard_constraints:
        assay  = r"chipseq|cuttag",
        source = r"pooled|biorep\d+",
        pr     = r"[12]",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_broad_{{assay}}_{{source}}_pr{{pr}}_idr.broad.macs3.log",
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
        elif [[ ! -f "{output}" ]]; then
            echo "ERROR: Expected MACS3 output not found: {output}" >&2
            exit 1
        fi
        """


# ---------------------------------------------------------------------------
# 5. Self-pseudoreplicate IDR
# ---------------------------------------------------------------------------

rule idr_self_pseudoreps_broad:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"self_pseudoreplicates/"
                  f"{{experiment}}_broad_{{assay}}_biorep{{bio_rep}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"self_pseudoreplicates/"
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
            f"self_pseudoreplicates/{wc.experiment}_broad_{wc.assay}_"
            f"biorep{wc.bio_rep}_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_"
            f"biorep{wc.bio_rep}_idr_self.broad.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_"
            f"biorep{wc.bio_rep}_idr_self_thr.broad.log"
        ),
    wildcard_constraints:
        assay   = r"chipseq|cuttag",
        bio_rep = r"\d+",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.idr_out:q})" \
                 "$(dirname {params.raw_log_file:q})" \
                 "$(dirname {params.thr_log_file:q})"

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


# ---------------------------------------------------------------------------
# 6. Pooled-pseudoreplicate IDR
# ---------------------------------------------------------------------------

rule idr_pooled_pseudoreps_broad:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"pooled_pseudoreplicates/{{experiment}}_broad_{{assay}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"pooled_pseudoreplicates/{{experiment}}_broad_{{assay}}_idr."
                  f"thresholded.broadPeak",
    input:
        peaks1 = lambda wc: idr_pooled_peak_input(
            wc.experiment, 1, wc.assay, "broadPeak"
        ),
        peaks2 = lambda wc: idr_pooled_peak_input(
            wc.experiment, 2, wc.assay, "broadPeak"
        ),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"pooled_pseudoreplicates/{wc.experiment}_broad_{wc.assay}_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_idr_pooled.broad.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_broad_{wc.assay}_idr_pooled_thr.broad.log"
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


# ---------------------------------------------------------------------------
# 7. Reproducibility summary + final peak
# ---------------------------------------------------------------------------

rule idr_summary_chipseq_broad:
    output:
        summary    = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/reproducibility_summary.tsv",
        final_peak = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/{{experiment}}.chipseq.macs3.broad."
                     f"replicate_validated.idr.broadPeak",
    input:
        true_thresh = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                      f"true_replicates/{{experiment}}_broad_chipseq_idr.thresholded.broadPeak",
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
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}.chipseq.macs3.broad.idr.summary.log",
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.summary:q})" \
                 "$(dirname {output.final_peak:q})" \
                 "$(dirname {log:q})"

        python3 scripts/idr_reproducibility_summary.py \
            --true-peaks {input.true_thresh:q} \
            --pooled-peaks {input.pool_thresh:q} \
            --self1-peaks {input.self1_thresh:q} \
            --self2-peaks {input.self2_thresh:q} \
            --experiment {params.experiment:q} \
            --assay {params.assay} \
            --caller {params.caller} \
            --peak-mode {params.peak_mode} \
            --bio-rep-a {params.bio_rep_a} \
            --bio-rep-b {params.bio_rep_b} \
            --final-method {params.final_method} \
            --final-output {params.final_output:q} \
            --output-tsv {output.summary:q} \
            --output-peak {output.final_peak:q} \
            2>&1 | tee {log:q}
        """


rule idr_summary_cuttag_broad:
    output:
        summary    = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/reproducibility_summary.tsv",
        final_peak = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/{{experiment}}.cuttag.macs3.broad."
                     f"replicate_validated.idr.broadPeak",
    input:
        true_thresh = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                      f"true_replicates/{{experiment}}_broad_cuttag_idr.thresholded.broadPeak",
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
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}.cuttag.macs3.broad.idr.summary.log",
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.summary:q})" \
                 "$(dirname {output.final_peak:q})" \
                 "$(dirname {log:q})"

        python3 scripts/idr_reproducibility_summary.py \
            --true-peaks {input.true_thresh:q} \
            --pooled-peaks {input.pool_thresh:q} \
            --self1-peaks {input.self1_thresh:q} \
            --self2-peaks {input.self2_thresh:q} \
            --experiment {params.experiment:q} \
            --assay {params.assay} \
            --caller {params.caller} \
            --peak-mode {params.peak_mode} \
            --bio-rep-a {params.bio_rep_a} \
            --bio-rep-b {params.bio_rep_b} \
            --final-method {params.final_method} \
            --final-output {params.final_output:q} \
            --output-tsv {output.summary:q} \
            --output-peak {output.final_peak:q} \
            2>&1 | tee {log:q}
        """
