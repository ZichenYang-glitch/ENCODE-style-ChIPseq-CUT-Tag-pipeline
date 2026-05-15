# report.smk — Final aggregation and sentinel rules
# ==================================================

# ---------------------------------------------------------------------------
# Pipeline done sentinel — treatment samples only
# ---------------------------------------------------------------------------

rule pipeline_done:
    output:
        done = f"{OUTDIR}/{{sample}}/logs/{{sample}}.pipeline.done",
    input:
        bw       = f"{OUTDIR}/{{sample}}/03_bigwig/{{sample}}.CPM.bw",
        peaks    = f"{OUTDIR}/{{sample}}/04_peaks/{{sample}}",
        fastqc   = f"{OUTDIR}/{{sample}}/logs/{{sample}}.fastqc.done",
        trim     = f"{OUTDIR}/{{sample}}/logs/{{sample}}.trim.done",
        bai      = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.sorted.bam.bai",
        flagstat = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.flagstat.txt",
        final_fs = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.final.flagstat.txt",
        idxstats = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.idxstats.txt",
        dup_met  = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.dup_metrics.txt",
        qc_summ  = (
            f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.qc_summary.tsv"
            if QC_CONFIG.get("summary", True) else []
        ),
    shell:
        "touch {output.done}"


# ---------------------------------------------------------------------------
# MultiQC aggregation
# ---------------------------------------------------------------------------

if MULTIQC:

    rule multiqc:
        input:
            # Treatment samples: pipeline.done ensures all QC + peaks + bigWig
            # are complete before MultiQC runs.
            [f"{OUTDIR}/{sid}/logs/{sid}.pipeline.done"
             for sid in TREATMENT_SAMPLE_IDS],
            # Control samples: no pipeline.done. Pull every QC artifact
            # explicitly so MultiQC does not run before side-branch QC exists.
            # Only active when use_control is true.
            [f"{OUTDIR}/{sid}/03_bigwig/{sid}.CPM.bw"
             for sid in CONTROL_SAMPLE_IDS] if USE_CONTROL else [],
            [f"{OUTDIR}/{sid}/logs/{sid}.fastqc.done"
             for sid in CONTROL_SAMPLE_IDS] if USE_CONTROL else [],
            [f"{OUTDIR}/{sid}/01_qc/{sid}.flagstat.txt"
             for sid in CONTROL_SAMPLE_IDS] if USE_CONTROL else [],
            [f"{OUTDIR}/{sid}/01_qc/{sid}.final.flagstat.txt"
             for sid in CONTROL_SAMPLE_IDS] if USE_CONTROL else [],
            [f"{OUTDIR}/{sid}/01_qc/{sid}.idxstats.txt"
             for sid in CONTROL_SAMPLE_IDS] if USE_CONTROL else [],
            [f"{OUTDIR}/{sid}/01_qc/{sid}.dup_metrics.txt"
             for sid in CONTROL_SAMPLE_IDS] if USE_CONTROL else [],
            [f"{OUTDIR}/{sid}/logs/{sid}.trim.done"
             for sid in CONTROL_SAMPLE_IDS] if USE_CONTROL else [],
        output:
            f"{OUTDIR}/multiqc/multiqc_report.html",
        log:
            f"{OUTDIR}/logs/multiqc.log",
        params:
            logdir      = f"{OUTDIR}/logs",
            multiqc_dir = f"{OUTDIR}/multiqc",
            sample_dirs = [f"{OUTDIR}/{sid}" for sid in ACTIVE_SAMPLE_IDS],
        conda:
            "../envs/chipseq.yml",
        shell:
            """
            command -v multiqc >/dev/null 2>&1 || {{
                echo "ERROR: multiqc not found. Install it or set multiqc: false in config.yaml." >&2
                exit 1
            }}
            mkdir -p {params.logdir:q} {params.multiqc_dir:q}
            multiqc {params.sample_dirs:q} -o {params.multiqc_dir:q} 2>&1 | tee {log:q}
            """
