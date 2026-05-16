# idr.smk — Stage 5a TF ChIP-seq True Replicate IDR
# ===================================================
# Produces IDR-ready MACS3 peak calls per biological replicate and
# runs true-replicate IDR between them. Pseudoreplicates, self-IDR,
# pooled-IDR, and conservative/optimal peak sets are deferred to Stage 5b.
#
# Rules:
#   macs3_idr_biorep    — IDR-ready MACS3 call on a single biorep BAM
#   idr_true_replicates — IDR between the two biorep IDR peak sets


# ---------------------------------------------------------------------------
# Helper: inputs for macs3_idr_biorep
# ---------------------------------------------------------------------------

def _idr_biorep_peaks_inputs(wildcards):
    """Return inputs for MACS3 IDR peak call on a single biorep BAM.

    Input order: [biorep_bam, biorep_bam.bai, ...optional_pooled_control_bam]
    The pooled control BAM is included as an explicit dependency so Snakemake
    schedules pool_control_bam before this rule.
    """
    exp = wildcards.experiment
    br = int(wildcards.bio_rep)
    inputs = [
        f"{OUTDIR}/experiments/{exp}/02_align/biorep{br}.final.bam",
        f"{OUTDIR}/experiments/{exp}/02_align/biorep{br}.final.bam.bai",
    ]
    if exp in POOLED_CONTROL_EXPERIMENTS:
        inputs.append(
            f"{OUTDIR}/experiments/{exp}/02_align/{exp}.pooled.control.final.bam"
        )
    return inputs


# ---------------------------------------------------------------------------
# Helper: MACS3 args for IDR-ready peak calls
# ---------------------------------------------------------------------------

def _idr_macs3_args(wildcards):
    """Return MACS3 args for IDR-ready peak calls on a biorep BAM.

    Uses the same layout/genome as per-sample MACS3, but replaces
    -q (q-value) with -p (p-value from idr_macs3 config).
    Never emits both -q and -p.
    """
    experiment = wildcards.experiment
    treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    if not treatment_ids:
        return ""
    s = SAMPLE_MAP[treatment_ids[0]]
    fmt = "BAMPE" if s["layout"] == "PE" else "BAM"
    genome = _normalize_genome(s["genome"])
    pvalue = _tool_param("idr_macs3", "pvalue", 0.1)
    extra = _tool_param("idr_macs3", "extra_args", "")
    return f"-f {fmt} -g {genome} -p {pvalue} {extra}".strip()


# ---------------------------------------------------------------------------
# 1. MACS3 IDR-ready peak call on a single biorep BAM
# ---------------------------------------------------------------------------

rule macs3_idr_biorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/04_peaks/idr/{{experiment}}_biorep{{bio_rep}}_idr_peaks.narrowPeak",
    input:
        lambda wc: _idr_biorep_peaks_inputs(wc),
    params:
        macs3_args = lambda wc: _idr_macs3_args(wc),
        sample     = lambda wc: f"{wc.experiment}_biorep{wc.bio_rep}_idr",
    wildcard_constraints:
        bio_rep = r"\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}_biorep{{bio_rep}}.idr.macs3.log",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output})" "$(dirname {log})"

        set -- {input:q}
        TREATMENT="$1"
        # Input order: treatment.bam ($1)[, treatment.bam.bai ($2)][, pooled.control.bam ($3)]
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


# ---------------------------------------------------------------------------
# Helper: resolve IDR peak input path for a given bio_rep index
# ---------------------------------------------------------------------------

def _idr_peak_input(experiment, index):
    """Return the IDR peak file for a bio_rep by 0-based index.

    index 0 -> first bio_rep (smallest label), index 1 -> second.
    The actual bio_rep numbers are derived from the sample sheet.
    """
    bioreps = _bioreps_for(experiment, "treatment")
    br = bioreps[index]
    return (
        f"{OUTDIR}/experiments/{experiment}/04_peaks/idr/"
        f"{experiment}_biorep{br}_idr_peaks.narrowPeak"
    )


# ---------------------------------------------------------------------------
# 2. True replicate IDR
# ---------------------------------------------------------------------------

rule idr_true_replicates:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/true_replicates/idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/true_replicates/idr.thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: _idr_peak_input(wc.experiment, 0),
        peaks2 = lambda wc: _idr_peak_input(wc.experiment, 1),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_idr/true_replicates/idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/{wc.experiment}.idr.raw.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/{wc.experiment}.idr.thresholded.log"
        ),
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        command -v idr >/dev/null 2>&1 || {{
            echo "ERROR: idr is required for Stage 5 IDR analysis but was not found in PATH." >&2
            exit 1
        }}

        mkdir -p "$(dirname {output.idr_out})" \
                 "$(dirname {params.raw_log_file})" \
                 "$(dirname {params.thr_log_file})"

        # Run 1: raw IDR output (no --idr-threshold).
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --output-file {output.idr_out:q} \
            --log-output-file {params.raw_log_file:q}

        # Run 2: thresholded narrowPeak output.
        # --idr-threshold N keeps peaks whose accumulated global IDR
        # score falls at or below N (smaller = more stringent).
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            --log-output-file {params.thr_log_file:q}
        """
