# idr_cuttag.smk — Stage 64 CUT&Tag narrow IDR
# ===============================================
# Mirrors workflow/rules/idr_atac.smk for CUT&Tag narrow MACS3 IDR.
# All helpers use _cuttag_ prefix to avoid namespace collisions.
#
# Rules:
#   cuttag_macs3_idr_biorep      — IDR-ready MACS3 on biorep BAM
#   cuttag_idr_true_replicates    — true-replicate IDR
#   cuttag_split_pseudoreps       — BAM pseudorep split
#   cuttag_macs3_idr_pseudorep    — MACS3 on pseudorep BAMs
#   cuttag_idr_self_pseudoreps    — self-pseudorep IDR
#   cuttag_idr_pooled_pseudoreps  — pooled pseudorep IDR
#   cuttag_idr_summary            — reproducibility QC + final peak


# ---------------------------------------------------------------------------
# Helper: MACS3 args for CUT&Tag IDR-ready calls
# ---------------------------------------------------------------------------

def _cuttag_idr_macs3_args(wildcards):
    """Return MACS3 args for CUT&Tag narrow IDR biorep peak call.

    Uses layout-aware format, Tn5 narrow shift, and RELAXED p-value (-p),
    not standard q-value (-q). Never emits both -p and -q.
    """
    experiment = wildcards.experiment
    treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    if not treatment_ids:
        return ""
    first = SAMPLE_MAP[treatment_ids[0]]
    layout = first.get("layout", "PE")
    fmt = "BAMPE" if layout == "PE" else "BAM"
    genome = _normalize_genome(first["genome"])
    pvalue = _tool_param("idr_macs3", "pvalue", 0.1)
    extra = _tool_param("idr_macs3", "extra_args", "")
    return (
        f"-f {fmt} -g {genome} "
        f"-p {pvalue} "
        f"--nomodel --shift -100 --extsize 200 "
        f"{extra}"
    ).strip()


# ---------------------------------------------------------------------------
# Helper: inputs for CUT&Tag IDR biorep MACS3
# ---------------------------------------------------------------------------

def _cuttag_idr_biorep_peaks_inputs(wildcards):
    """Return inputs for CUT&Tag IDR per-biorep MACS3.

    Same pooled-control policy as ATAC IDR helpers.
    """
    exp = wildcards.experiment
    br = int(wildcards.bio_rep)
    inputs = [
        f"{OUTDIR}/experiments/{exp}/02_align/biorep{br}.final.bam",
        f"{OUTDIR}/experiments/{exp}/02_align/biorep{br}.final.bam.bai",
    ]
    if exp in POOLED_CONTROL_EXPERIMENTS:
        inputs.append(
            f"{OUTDIR}/experiments/{exp}/02_align/"
            f"{exp}.pooled.control.final.bam"
        )
    return inputs


# ---------------------------------------------------------------------------
# Helper: resolve CUT&Tag IDR per-biorep peak files
# ---------------------------------------------------------------------------

def _cuttag_idr_peak_input(experiment, index):
    """Return the IDR-ready biorep peak file at the given 0-based index."""
    bioreps = sorted(_bioreps_for(experiment, "treatment"))
    br = bioreps[index]
    return (
        f"{OUTDIR}/experiments/{experiment}/06_reproducibility/idr/"
        f"idr_peaks/{experiment}_cuttag_biorep{br}_idr.narrowPeak"
    )


# ---------------------------------------------------------------------------
# Helper: self-IDR thresholded path
# ---------------------------------------------------------------------------

def _cuttag_self_thresh_path(experiment, index):
    """Return self-IDR thresholded path for the given 0-based biorep index."""
    bioreps = sorted(_bioreps_for(experiment, "treatment"))
    br = bioreps[index]
    return (
        f"{OUTDIR}/experiments/{experiment}/06_reproducibility/idr/"
        f"self_pseudoreplicates/"
        f"{experiment}_cuttag_biorep{br}_idr.thresholded.narrowPeak"
    )


# ============================================================================
# 1. Per-biorep IDR-ready MACS3
# ============================================================================

rule cuttag_macs3_idr_biorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
        f"idr_peaks/{{experiment}}_cuttag_biorep{{bio_rep}}_idr.narrowPeak",
    input:
        lambda wc: _cuttag_idr_biorep_peaks_inputs(wc),
    params:
        macs3_args = lambda wc: _cuttag_idr_macs3_args(wc),
        sample     = lambda wc: (
            f"{wc.experiment}_cuttag_biorep{wc.bio_rep}_idr"
        ),
    wildcard_constraints:
        bio_rep = r"\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_cuttag_biorep{{bio_rep}}_idr.macs3.log",
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

        # MACS3 writes {name}_peaks.narrowPeak. Move to canonical path.
        MACS3_OUT="$(dirname {output:q})/{params.sample}_peaks.narrowPeak"
        if [[ -f "$MACS3_OUT" ]]; then
            mv "$MACS3_OUT" {output:q}
        else
            echo "ERROR: Expected MACS3 output not found: $MACS3_OUT" >&2
            exit 1
        fi
        """


# ============================================================================
# 2. True-replicate IDR
# ============================================================================

rule cuttag_idr_true_replicates:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/true_replicates/{{experiment}}_cuttag_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/true_replicates/"
                  f"{{experiment}}_cuttag_idr.thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: _cuttag_idr_peak_input(wc.experiment, 0),
        peaks2 = lambda wc: _cuttag_idr_peak_input(wc.experiment, 1),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"true_replicates/{wc.experiment}_cuttag_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_cuttag_idr_true.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_cuttag_idr_true_thresholded.log"
        ),
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_cuttag_idr_true_thresholded.log",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.idr_out:q})" "$(dirname {log:q})"

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank:q} \
            --output-file {params.output_prefix:q}.txt \
            --plot \
            2>&1 | tee {params.raw_log_file:q}

        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank:q} \
            --idr-threshold {params.threshold} \
            --output-file {params.output_prefix:q}.thresholded.narrowPeak \
            2>&1 | tee {log:q}
        """


# ============================================================================
# 3. BAM pseudorep split
# ============================================================================

rule cuttag_split_pseudoreps:
    output:
        pr1 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_cuttag_{{source}}.pr1.bam",
        p1i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_cuttag_{{source}}.pr1.bam.bai",
        pr2 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_cuttag_{{source}}.pr2.bam",
        p2i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/"
              f"{{experiment}}_cuttag_{{source}}.pr2.bam.bai",
    input:
        lambda wc: _cuttag_split_input(wildcards=wc),
    wildcard_constraints:
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


def _cuttag_split_input(wildcards):
    """Return the BAM to split for pseudoreps."""
    exp = wildcards.experiment
    src = wildcards.source
    if src == "pooled":
        return (
            f"{OUTDIR}/experiments/{exp}/02_align/"
            f"{exp}.pooled.final.bam"
        )
    # src is "biorepN" — extract N and return the biorep BAM
    br = int(src.replace("biorep", ""))
    return (
        f"{OUTDIR}/experiments/{exp}/02_align/"
        f"biorep{br}.final.bam"
    )


# ============================================================================
# 4. MACS3 on pseudorep BAMs
# ============================================================================

rule cuttag_macs3_idr_pseudorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/idr/"
        f"idr_peaks/{{experiment}}_cuttag_{{source}}_pr{{pr}}_idr.narrowPeak",
    input:
        lambda wc: _cuttag_idr_pseudorep_inputs(wildcards=wc),
    params:
        macs3_args = lambda wc: _cuttag_idr_macs3_args(wc),
        sample     = lambda wc: (
            f"{wc.experiment}_cuttag_{wc.source}_pr{wc.pr}_idr"
        ),
    wildcard_constraints:
        source = r"pooled|biorep\d+",
        pr     = r"[12]",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}_cuttag_{{source}}_pr{{pr}}_idr.macs3.log",
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
        else
            echo "ERROR: Expected MACS3 output not found: $MACS3_OUT" >&2
            exit 1
        fi
        """


def _cuttag_idr_pseudorep_inputs(wildcards):
    """Return inputs for CUT&Tag IDR pseudorep MACS3."""
    exp = wildcards.experiment
    src = wildcards.source
    pr = wildcards.pr
    inputs = [
        f"{OUTDIR}/experiments/{exp}/05_pseudorep/"
        f"{exp}_cuttag_{src}.pr{pr}.bam",
        f"{OUTDIR}/experiments/{exp}/05_pseudorep/"
        f"{exp}_cuttag_{src}.pr{pr}.bam.bai",
    ]
    if exp in POOLED_CONTROL_EXPERIMENTS:
        inputs.append(
            f"{OUTDIR}/experiments/{exp}/02_align/"
            f"{exp}.pooled.control.final.bam"
        )
    return inputs


# ============================================================================
# 5. Self-pseudorep IDR
# ============================================================================

rule cuttag_idr_self_pseudoreps:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/self_pseudoreplicates/"
                  f"{{experiment}}_cuttag_biorep{{bio_rep}}_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/self_pseudoreplicates/"
                  f"{{experiment}}_cuttag_biorep{{bio_rep}}_idr."
                  f"thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"idr_peaks/{wc.experiment}_cuttag_biorep{wc.bio_rep}_pr1_idr.narrowPeak"
        ),
        peaks2 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"idr_peaks/{wc.experiment}_cuttag_biorep{wc.bio_rep}_pr2_idr.narrowPeak"
        ),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"self_pseudoreplicates/{wc.experiment}_cuttag_biorep{wc.bio_rep}_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_cuttag_biorep{wc.bio_rep}_idr_self.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_cuttag_biorep{wc.bio_rep}_idr_self_thr.log"
        ),
    wildcard_constraints:
        bio_rep = r"\d+",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.idr_out:q})" "$(dirname {output.thr_out:q})"

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


# ============================================================================
# 6. Pooled-pseudorep IDR
# ============================================================================

rule cuttag_idr_pooled_pseudoreps:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/pooled_pseudoreplicates/{{experiment}}_cuttag_idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"idr/pooled_pseudoreplicates/"
                  f"{{experiment}}_cuttag_idr.thresholded.narrowPeak",
    input:
        peaks1 = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                 f"idr/idr_peaks/{{experiment}}_cuttag_pooled_pr1_idr.narrowPeak",
        peaks2 = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                 f"idr/idr_peaks/{{experiment}}_cuttag_pooled_pr2_idr.narrowPeak",
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/idr/"
            f"pooled_pseudoreplicates/{wc.experiment}_cuttag_idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_cuttag_idr_pooled.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}_cuttag_idr_pooled_thr.log"
        ),
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


# ============================================================================
# 7. Reproducibility summary + final peak
# ============================================================================

rule cuttag_idr_summary:
    output:
        summary    = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/reproducibility_summary.tsv",
        final_peak = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                     f"final/{{experiment}}.cuttag.macs3.narrow."
                     f"replicate_validated.idr.narrowPeak",
    input:
        true_thresh = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                      f"idr/true_replicates/"
                      f"{{experiment}}_cuttag_idr.thresholded.narrowPeak",
        pool_thresh = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                      f"idr/pooled_pseudoreplicates/"
                      f"{{experiment}}_cuttag_idr.thresholded.narrowPeak",
        self1_thresh = lambda wc: _cuttag_self_thresh_path(wc.experiment, 0),
        self2_thresh = lambda wc: _cuttag_self_thresh_path(wc.experiment, 1),
    params:
        experiment    = lambda wc: wc.experiment,
        assay         = "cuttag",
        caller        = "macs3",
        peak_mode     = "narrow",
        bio_rep_a     = lambda wc: _idr_biorep_labels(wc.experiment)[0],
        bio_rep_b     = lambda wc: _idr_biorep_labels(wc.experiment)[1],
        final_method  = "idr",
        final_output  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_reproducibility/final/"
            f"{wc.experiment}.cuttag.macs3.narrow."
            f"replicate_validated.idr.narrowPeak"
        ),
    wildcard_constraints:
        experiment = wildcard_choices(CUTTAG_IDR_EXPERIMENTS),
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.summary:q})" "$(dirname {output.final_peak:q})"

        python3 scripts/cuttag_idr_summary.py \
            --true-peaks {input.true_thresh:q} \
            --pooled-peaks {input.pool_thresh:q} \
            --self1-peaks {input.self1_thresh:q} \
            --self2-peaks {input.self2_thresh:q} \
            --experiment {params.experiment:q} \
            --assay cuttag \
            --caller macs3 \
            --peak-mode narrow \
            --bio-rep-a {params.bio_rep_a} \
            --bio-rep-b {params.bio_rep_b} \
            --final-method idr \
            --final-output {params.final_output:q} \
            --output-tsv {output.summary:q} \
            --output-peak {output.final_peak:q}
        """


# ---------------------------------------------------------------------------
# Helper: resolve biorep labels for IDR summary
# ---------------------------------------------------------------------------

def _idr_biorep_labels(experiment):
    """Return (br_a, br_b) as strings for the experiment's 2 bioreps."""
    bioreps = sorted(_bioreps_for(experiment, "treatment"))
    return str(bioreps[0]), str(bioreps[1])
