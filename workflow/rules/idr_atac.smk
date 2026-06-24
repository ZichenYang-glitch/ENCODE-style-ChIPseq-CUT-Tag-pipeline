# idr_atac.smk — Stage 55 ATAC-seq narrowPeak IDR
# ===================================================
# Parallel path to legacy Stage 5 ChIP-seq IDR (idr.smk).
# Uses Stage 53 namespace: 06_reproducibility/idr/ and 06_reproducibility/final/.
# All helper functions use _atac_ prefix to avoid Snakemake shared-namespace
# collisions with legacy idr.smk helpers.
#
# Rules:
#   atac_macs3_idr_biorep      — IDR-ready MACS3 call on ATAC biorep BAM
#   atac_idr_true_replicates    — IDR between the two ATAC biorep peak sets
#   atac_split_pseudoreps       — deterministic BAM pseudorep split
#   atac_macs3_idr_pseudorep    — IDR-ready MACS3 on ATAC pseudorep BAM
#   atac_idr_self_pseudoreps    — self-IDR per ATAC biorep
#   atac_idr_pooled_pseudoreps  — pooled pseudorep IDR for ATAC
#   atac_idr_summary            — reproducibility QC + final peak (Stage 55)


# ---------------------------------------------------------------------------
# Helper: inputs for atac_macs3_idr_biorep
# ---------------------------------------------------------------------------

def _atac_idr_biorep_peaks_inputs(wildcards):
    """Return inputs for MACS3 IDR peak call on a single ATAC biorep BAM."""
    return idr_biorep_peaks_inputs(
        wildcards.experiment, int(wildcards.bio_rep)
    )


# ---------------------------------------------------------------------------
# Helper: MACS3 args for ATAC IDR-ready peak calls
# ---------------------------------------------------------------------------

def _atac_idr_macs3_args(wildcards):
    """Return MACS3 args for IDR-ready ATAC narrow peak calls."""
    return idr_macs3_args(wildcards.experiment, "atac", "narrow")


# ---------------------------------------------------------------------------
# 1. MACS3 IDR-ready peak call on a single ATAC biorep BAM
# ---------------------------------------------------------------------------

rule atac_macs3_idr_biorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
        f"idr_peaks/{{experiment}}_atac_biorep{{bio_rep}}_idr.narrowPeak",
    input:
        lambda wc: _atac_idr_biorep_peaks_inputs(wc),
    params:
        macs3_args = lambda wc: _atac_idr_macs3_args(wc),
        sample     = lambda wc: f"{wc.experiment}_atac_biorep{wc.bio_rep}_idr",
    wildcard_constraints:
        bio_rep = r"\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_atac_biorep{{bio_rep}}.idr.macs3.log",
    threads: THREADS,
    conda:
        "../envs/macs3.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output})" "$(dirname {log})"

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

        EXPECTED="{output}"
        if [[ ! -f "$EXPECTED" ]]; then
            echo "ERROR: Expected peak file not found: $EXPECTED" >&2
            exit 1
        fi
        """



rule atac_idr_true_replicates:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"true_replicates/{{experiment}}_atac_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"true_replicates/{{experiment}}_atac_idr.thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: idr_repro_peak_input(wc.experiment, 0, "atac", "narrowPeak"),
        peaks2 = lambda wc: idr_repro_peak_input(wc.experiment, 1, "atac", "narrowPeak"),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"true_replicates/{wc.experiment}_atac_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.atac_idr.raw.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.atac_idr.thresholded.log"
        ),
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        command -v idr >/dev/null 2>&1 || {{
            echo "ERROR: idr is required for ATAC IDR analysis but was not found in PATH." >&2
            exit 1
        }}

        mkdir -p "$(dirname {output.idr_out})" \
                 "$(dirname {params.raw_log_file})" \
                 "$(dirname {params.thr_log_file})"

        # Run 1: raw IDR output
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --output-file {output.idr_out:q} \
            --log-output-file {params.raw_log_file:q}

        # Run 2: thresholded narrowPeak output
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            --log-output-file {params.thr_log_file:q}
        """


# ============================================================================
# Stage 55 ATAC pseudoreplicate rules
# ============================================================================

# ---------------------------------------------------------------------------
# Helper: input BAM for ATAC pseudorep splitting
# ---------------------------------------------------------------------------

def _atac_split_input(wildcards):
    """Return the input BAM path for an ATAC pseudorep split."""
    return idr_split_input(wildcards.experiment, wildcards.source)


# ---------------------------------------------------------------------------
# Helper: inputs for ATAC pseudorep MACS3 peak call
# ---------------------------------------------------------------------------

def _atac_idr_pseudorep_inputs(wildcards):
    """Return inputs for MACS3 IDR peak call on an ATAC pseudorep BAM."""
    return idr_pseudorep_peaks_inputs(
        wildcards.experiment, wildcards.source, wildcards.pr, source_prefix="atac_"
    )



rule atac_split_pseudoreps:
    output:
        pr1 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_atac_{{source}}.pr1.bam",
        p1i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_atac_{{source}}.pr1.bam.bai",
        pr2 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_atac_{{source}}.pr2.bam",
        p2i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_atac_{{source}}.pr2.bam.bai",
    input:
        lambda wc: _atac_split_input(wc),
    params:
        seed = IDR_SEED,
    wildcard_constraints:
        source = r"pooled|biorep\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_atac_{{source}}.split.log",
    threads: THREADS,
    conda:
        "../envs/samtools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.pr1})" "$(dirname {log})"
        python3 scripts/split_pseudoreps.py \
            --input {input:q} \
            --out1 {output.pr1:q} \
            --out2 {output.pr2:q} \
            --seed {params.seed} \
            --threads {threads} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 4. IDR-ready MACS3 on ATAC pseudorep BAM
# ---------------------------------------------------------------------------

rule atac_macs3_idr_pseudorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
        f"idr_peaks/{{experiment}}_atac_{{source}}_pr{{pr}}_idr.narrowPeak",
    input:
        lambda wc: _atac_idr_pseudorep_inputs(wc),
    params:
        macs3_args = lambda wc: _atac_idr_macs3_args(wc),
        sample     = lambda wc: f"{wc.experiment}_atac_{wc.source}_pr{wc.pr}_idr",
    wildcard_constraints:
        source = r"pooled|biorep\d+",
        pr     = r"[12]",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_atac_{{source}}_pr{{pr}}.idr.macs3.log",
    threads: THREADS,
    conda:
        "../envs/macs3.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output})" "$(dirname {log})"

        set -- {input:q}
        TREATMENT="$1"
        if [[ $# -ge 3 ]]; then
            macs3 callpeak -t "$TREATMENT" -c "$3" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} 2>&1 | tee {log:q}
        else
            macs3 callpeak -t "$TREATMENT" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} 2>&1 | tee {log:q}
        fi

        [[ -f "{output}" ]] || {{ echo "ERROR: Expected {output}" >&2; exit 1; }}
        """


# ---------------------------------------------------------------------------
# 5. Self-pseudoreplicate IDR per ATAC biological replicate
# ---------------------------------------------------------------------------

rule atac_idr_self_pseudoreps:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"self_pseudoreplicates/{{experiment}}_atac_biorep{{bio_rep}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"self_pseudoreplicates/{{experiment}}_atac_biorep{{bio_rep}}_idr.thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"idr_peaks/{wc.experiment}_atac_biorep{wc.bio_rep}_pr1_idr.narrowPeak"
        ),
        peaks2 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"idr_peaks/{wc.experiment}_atac_biorep{wc.bio_rep}_pr2_idr.narrowPeak"
        ),
    params:
        threshold    = IDR_THRESHOLD,
        rank         = IDR_RANK,
        raw_log_file = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.atac_idr.self.biorep{wc.bio_rep}.raw.log"
        ),
        thr_log_file = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.atac_idr.self.biorep{wc.bio_rep}.thr.log"
        ),
    wildcard_constraints:
        bio_rep = r"\d+",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        command -v idr >/dev/null 2>&1 || {{
            echo "ERROR: idr is required for ATAC IDR analysis but was not found in PATH." >&2
            exit 1
        }}

        mkdir -p "$(dirname {output.idr_out})" \
                 "$(dirname {params.raw_log_file})" \
                 "$(dirname {params.thr_log_file})"

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --output-file {output.idr_out:q} \
            --log-output-file {params.raw_log_file:q}

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            --log-output-file {params.thr_log_file:q}
        """


# ---------------------------------------------------------------------------
# 6. Pooled-pseudoreplicate IDR for ATAC
# ---------------------------------------------------------------------------

rule atac_idr_pooled_pseudoreps:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"pooled_pseudoreplicates/{{experiment}}_atac_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                  f"pooled_pseudoreplicates/{{experiment}}_atac_idr.thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: idr_pooled_peak_input(wc.experiment, 1, "atac", "narrowPeak"),
        peaks2 = lambda wc: idr_pooled_peak_input(wc.experiment, 2, "atac", "narrowPeak"),
    params:
        threshold    = IDR_THRESHOLD,
        rank         = IDR_RANK,
        raw_log_file = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.atac_idr.pooled.raw.log"
        ),
        thr_log_file = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.atac_idr.pooled.thr.log"
        ),
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        command -v idr >/dev/null 2>&1 || {{
            echo "ERROR: idr is required for ATAC IDR analysis but was not found in PATH." >&2
            exit 1
        }}

        mkdir -p "$(dirname {output.idr_out})" \
                 "$(dirname {params.raw_log_file})" \
                 "$(dirname {params.thr_log_file})"

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --output-file {output.idr_out:q} \
            --log-output-file {params.raw_log_file:q}

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            --log-output-file {params.thr_log_file:q}
        """


# ---------------------------------------------------------------------------
# 7. ATAC IDR reproducibility QC summary and final peak set
# ---------------------------------------------------------------------------

rule atac_idr_summary:
    output:
        summary      = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/final/"
                       f"reproducibility_summary.tsv",
        final_peak   = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/final/"
                       f"{{experiment}}.atac.macs3.narrow.replicate_validated.idr.narrowPeak",
    input:
        true_thresh  = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
                       f"true_replicates/{{experiment}}_atac_idr.thresholded.narrowPeak",
        pool_thresh  = lambda wc: idr_pooled_thresh_path(wc.experiment, "atac", "narrowPeak"),
        self1_thresh = lambda wc: idr_self_thresh_path(wc.experiment, 0, "atac", "narrowPeak"),
        self2_thresh = lambda wc: idr_self_thresh_path(wc.experiment, 1, "atac", "narrowPeak"),
    params:
        experiment   = lambda wc: wc.experiment,
        assay        = "atac",
        caller       = "macs3",
        peak_mode    = "narrow",
        bio_rep_a    = lambda wc: idr_biorep_labels(wc.experiment)[0],
        bio_rep_b    = lambda wc: idr_biorep_labels(wc.experiment)[1],
        final_method = "idr",
        final_output = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/final/"
            f"{wc.experiment}.atac.macs3.narrow.replicate_validated.idr.narrowPeak"
        ),
    wildcard_constraints:
        experiment = wildcard_choices(ATAC_IDR_EXPERIMENTS),
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}.atac_idr.summary.log",
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.summary})" "$(dirname {log})"
        python3 scripts/atac_idr_summary.py \
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
