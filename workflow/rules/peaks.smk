# peaks.smk — MACS3 peak calling rule
# ====================================
# Runs only for treatment samples (expansion controlled by rule all in Snakefile).
# Assay-specific MACS3 parameters dispatched through get_macs3_args().
# Control input (if any) is tracked as a dependency via _macs3_inputs().

rule macs3_callpeak:
    output:
        peaks = directory(f"{OUTDIR}/{{sample}}/04_peaks/{{sample}}"),
    input:
        lambda wc: _macs3_inputs(wc),
    params:
        macs3_args = lambda wc: get_macs3_args(wc),
        bdg_args   = lambda wc: "-B" if QC_CONFIG.get("signal_tracks", True) else "",
        sample     = "{sample}",
        peak_mode  = lambda wc: SAMPLE_MAP[wc.sample]["peak_mode"],
        extra      = _tool_param("macs3", "extra_args", ""),
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.macs3.log",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p {output.peaks:q}

        set -- {input:q}
        TREATMENT="$1"
        # Input order: treatment.bam ($1), treatment.bam.bai ($2)[, control.bam ($3)]
        if [[ $# -ge 3 ]]; then
            macs3 callpeak \
                -t "$TREATMENT" \
                -c "$3" \
                -n {params.sample:q} \
                --outdir {output.peaks:q} \
                {params.macs3_args} \
                {params.bdg_args} \
                {params.extra} \
                2>&1 | tee {log:q}
        else
            macs3 callpeak \
                -t "$TREATMENT" \
                -n {params.sample:q} \
                --outdir {output.peaks:q} \
                {params.macs3_args} \
                {params.bdg_args} \
                {params.extra} \
                2>&1 | tee {log:q}
        fi

        # Verify expected peak file exists before Snakemake marks directory complete
        SUFFIX="narrowPeak"
        if [[ "{params.peak_mode}" == "broad" ]]; then
            SUFFIX="broadPeak"
        fi
        EXPECTED="{output.peaks}/{params.sample}_peaks.$SUFFIX"
        if [[ ! -f "$EXPECTED" ]]; then
            echo "ERROR: Expected peak file not found: $EXPECTED" >&2
            exit 1
        fi
        """


# ---------------------------------------------------------------------------
# Pooled MACS3 peak calling (Stage 4b)
# ---------------------------------------------------------------------------

def _pooled_macs3_args(experiment):
    """Return MACS3 args string for a pooled experiment.

    Delegates to get_macs3_args() so pooled peak calling uses the exact
    same ChIP-seq / CUT&Tag MACS3 policy as per-sample peak calling.
    """
    treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    if not treatment_ids:
        return ""
    class _Wildcards:
        sample = treatment_ids[0]
    return get_macs3_args(_Wildcards)


def _pooled_peak_suffix(experiment):
    """Return 'narrowPeak' or 'broadPeak' for the experiment."""
    treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    if not treatment_ids:
        return "narrowPeak"
    return "broadPeak" if SAMPLE_MAP[treatment_ids[0]]["peak_mode"] == "broad" else "narrowPeak"


def _pooled_peaks_inputs(wildcards):
    """Return inputs for macs3_pooled_peaks.

    Elements: [pooled.bam, pooled.bam.bai, ...optional_pooled_control]
    """
    exp = wildcards.experiment
    inputs = [
        f"{OUTDIR}/experiments/{exp}/02_align/{exp}.pooled.final.bam",
        f"{OUTDIR}/experiments/{exp}/02_align/{exp}.pooled.final.bam.bai",
    ]
    if exp in POOLED_CONTROL_EXPERIMENTS:
        inputs.append(
            f"{OUTDIR}/experiments/{exp}/02_align/{exp}.pooled.control.final.bam"
        )
    return inputs


rule macs3_pooled_peaks:
    output:
        peaks = directory(
            f"{OUTDIR}/experiments/{{experiment}}/04_peaks/pooled/{{experiment}}_pooled_peaks"
        ),
    input:
        lambda wc: _pooled_peaks_inputs(wc),
    params:
        macs3_args = lambda wc: _pooled_macs3_args(wc.experiment),
        bdg_args   = lambda wc: "-B" if QC_CONFIG.get("signal_tracks", True) else "",
        sample     = "{experiment}",
        peak_mode  = lambda wc: SAMPLE_MAP[
            TREATMENT_SAMPLES_BY_EXPERIMENT[wc.experiment][0]
        ]["peak_mode"],
        extra      = _tool_param("macs3", "extra_args", ""),
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.pooled.macs3.log",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p {output.peaks:q}

        set -- {input:q}
        TREATMENT="$1"
        # Input order: pooled.bam ($1), pooled.bam.bai ($2)[, pooled.control.bam ($3)]
        if [[ $# -ge 3 ]]; then
            macs3 callpeak \
                -t "$TREATMENT" \
                -c "$3" \
                -n {params.sample:q}_pooled \
                --outdir {output.peaks:q} \
                {params.macs3_args} \
                {params.bdg_args} \
                {params.extra} \
                2>&1 | tee {log:q}
        else
            macs3 callpeak \
                -t "$TREATMENT" \
                -n {params.sample:q}_pooled \
                --outdir {output.peaks:q} \
                {params.macs3_args} \
                {params.bdg_args} \
                {params.extra} \
                2>&1 | tee {log:q}
        fi

        # Verify expected peak file exists
        SUFFIX="narrowPeak"
        if [[ "{params.peak_mode}" == "broad" ]]; then
            SUFFIX="broadPeak"
        fi
        EXPECTED="{output.peaks}/{params.sample}_pooled_peaks.$SUFFIX"
        if [[ ! -f "$EXPECTED" ]]; then
            echo "ERROR: Expected peak file not found: $EXPECTED" >&2
            exit 1
        fi
        """
