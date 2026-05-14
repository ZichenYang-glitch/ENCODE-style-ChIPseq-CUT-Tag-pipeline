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
        sample     = "{sample}",
        peak_mode  = lambda wc: SAMPLE_MAP[wc.sample]["peak_mode"],
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.macs3.log",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
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
                2>&1 | tee {log:q}
        else
            macs3 callpeak \
                -t "$TREATMENT" \
                -n {params.sample:q} \
                --outdir {output.peaks:q} \
                {params.macs3_args} \
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
