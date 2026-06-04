# mnase.smk — MNase-seq nucleosome-centric rules and policy functions
# =====================================================================
# Policy functions define MNase-specific duplicate removal, read extension,
# and MACS3 (error-only) behavior. Executable rules produce mono-nucleosome
# BAM, dyad BigWig, and mono occupancy BigWig at sample and pooled levels.
#
# MNase samples reuse the shared preprocessing pipeline (FastQC, Trim Galore,
# Bowtie2, MAPQ filter, duplicate handling, final.bam) and skip the entire
# peak-centric path (MACS3, FRiP, peak counts, FE/ppois signal, IDR).
# Sub-nucleosome, di-nucleosome, DANPOS3, and iNPS are deferred to Stage 40.


# ---------------------------------------------------------------------------
# MNase config helper
# ---------------------------------------------------------------------------

def get_mono_range():
    """Return (min, max) mono-nucleosome fragment range from config.

    Defaults to [140, 200] when the mnase block is absent.
    """
    mnase_cfg = config.get("mnase", {})
    if isinstance(mnase_cfg, dict):
        mr = mnase_cfg.get("mono_range", [140, 200])
        if isinstance(mr, (list, tuple)) and len(mr) == 2:
            return int(mr[0]), int(mr[1])
    return 140, 200


# ---------------------------------------------------------------------------
# MNase policy functions — consumed by dispatch wrappers in workflow/Snakefile
# ---------------------------------------------------------------------------

def get_remove_dup_mnase(wildcards):
    """MNase duplicate removal policy.

    MNase is PE-only. auto -> yes (remove duplicates by default).
    Broad-peak exception does not apply.
    """
    policy = str(config.get("remove_dup", "auto"))
    if policy == "auto":
        return "yes"
    return policy


def get_extend_reads_mnase(wildcards):
    """MNase read extension policy for bamCoverage.

    PE MNase: no read extension. Fragment pairs are the biological unit.
    MNase is PE-only — always return empty string, regardless of the
    global extend_reads config value.
    """
    return ""


def get_macs3_args_mnase(wildcards):
    """MNase MACS3 guard — MNase samples must never reach MACS3.

    Raises ValueError because MNase samples should never be scheduled for
    macs3_callpeak. The target builder gates 04_peaks/ on PEAK_SAMPLE_IDS
    (which excludes MNase). If this function is reached, the DAG is
    misconfigured.
    """
    raise ValueError(
        f"Sample {wildcards.sample}: assay=mnase does not use MACS3 "
        f"peak calling. MNase samples should not reach get_macs3_args."
    )


# ---------------------------------------------------------------------------
# 1. Mono-nucleosome BAM — alignmentSieve from final.bam
# ---------------------------------------------------------------------------

rule mnase_split_mono:
    output:
        bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam",
        bai = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam.bai",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
        bai = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam.bai",
    params:
        mono_min = get_mono_range()[0],
        mono_max = get_mono_range()[1],
    threads: THREADS,
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p $(dirname {output.bam:q})

        alignmentSieve \
            -b {input.bam:q} \
            --minFragmentLength {params.mono_min} \
            --maxFragmentLength {params.mono_max} \
            -p {threads} \
            -o {output.bam:q}

        samtools index -@ {threads} {output.bam:q}
        """


# ---------------------------------------------------------------------------
# 2. Dyad BigWig — bamCoverage --MNase --binSize 1 from final.bam
# ---------------------------------------------------------------------------

rule mnase_dyad_bigwig:
    output:
        bw = f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.dyad.CPM.bw",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
        bai = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam.bai",
    params:
        normalize_using = f"--normalizeUsing {v}" if (v := _tool_param("bamcoverage", "normalize_using", "CPM")) != "" else "--normalizeUsing CPM",
        extra = _tool_param("bamcoverage", "extra_args", ""),
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.mnase_dyad_bw.log",
    threads: THREADS,
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p $(dirname {output.bw:q}) $(dirname {log:q})

        bamCoverage \
            -b {input.bam:q} \
            -o {output.bw:q} \
            --MNase \
            --binSize 1 \
            {params.normalize_using} \
            --numberOfProcessors {threads} \
            {params.extra} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 3. Mono occupancy BigWig — bamCoverage from mono.bam
# ---------------------------------------------------------------------------

rule mnase_mono_bigwig:
    output:
        bw = f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.mono.CPM.bw",
    input:
        bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam",
        bai = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam.bai",
    params:
        binsize = BINSIZE,
        normalize_using = f"--normalizeUsing {v}" if (v := _tool_param("bamcoverage", "normalize_using", "CPM")) != "" else "--normalizeUsing CPM",
        smooth_length = f"--smoothLength {v}" if (v := _tool_param("bamcoverage", "smooth_length", "")) != "" else "",
        extra = _tool_param("bamcoverage", "extra_args", ""),
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.mnase_mono_bw.log",
    threads: THREADS,
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p $(dirname {output.bw:q}) $(dirname {log:q})

        bamCoverage \
            -b {input.bam:q} \
            -o {output.bw:q} \
            {params.normalize_using} \
            --binSize {params.binsize} \
            {params.smooth_length} \
            --numberOfProcessors {threads} \
            {params.extra} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 4. Pooled mono-nucleosome BAM — alignmentSieve from pooled final.bam
# ---------------------------------------------------------------------------

rule mnase_pooled_mono:
    output:
        bam = f"{OUTDIR}/experiments/{{experiment}}/03_fragments/{{experiment}}.pooled.mono.bam",
        bai = f"{OUTDIR}/experiments/{{experiment}}/03_fragments/{{experiment}}.pooled.mono.bam.bai",
    input:
        bam = f"{OUTDIR}/experiments/{{experiment}}/02_align/{{experiment}}.pooled.final.bam",
        bai = f"{OUTDIR}/experiments/{{experiment}}/02_align/{{experiment}}.pooled.final.bam.bai",
    params:
        mono_min = get_mono_range()[0],
        mono_max = get_mono_range()[1],
    threads: THREADS,
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p $(dirname {output.bam:q})

        alignmentSieve \
            -b {input.bam:q} \
            --minFragmentLength {params.mono_min} \
            --maxFragmentLength {params.mono_max} \
            -p {threads} \
            -o {output.bam:q}

        samtools index -@ {threads} {output.bam:q}
        """


# ---------------------------------------------------------------------------
# 5. Pooled dyad BigWig — bamCoverage --MNase from pooled final.bam
# ---------------------------------------------------------------------------

rule mnase_pooled_dyad_bigwig:
    output:
        bw = f"{OUTDIR}/experiments/{{experiment}}/04_signal/{{experiment}}.pooled.dyad.CPM.bw",
    input:
        bam = f"{OUTDIR}/experiments/{{experiment}}/02_align/{{experiment}}.pooled.final.bam",
        bai = f"{OUTDIR}/experiments/{{experiment}}/02_align/{{experiment}}.pooled.final.bam.bai",
    params:
        normalize_using = f"--normalizeUsing {v}" if (v := _tool_param("bamcoverage", "normalize_using", "CPM")) != "" else "--normalizeUsing CPM",
        extra = _tool_param("bamcoverage", "extra_args", ""),
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.pooled.mnase_dyad_bw.log",
    threads: THREADS,
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p $(dirname {output.bw:q}) $(dirname {log:q})

        bamCoverage \
            -b {input.bam:q} \
            -o {output.bw:q} \
            --MNase \
            --binSize 1 \
            {params.normalize_using} \
            --numberOfProcessors {threads} \
            {params.extra} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 6. Pooled mono occupancy BigWig — bamCoverage from pooled mono.bam
# ---------------------------------------------------------------------------

rule mnase_pooled_mono_bigwig:
    output:
        bw = f"{OUTDIR}/experiments/{{experiment}}/04_signal/{{experiment}}.pooled.mono.CPM.bw",
    input:
        bam = f"{OUTDIR}/experiments/{{experiment}}/03_fragments/{{experiment}}.pooled.mono.bam",
        bai = f"{OUTDIR}/experiments/{{experiment}}/03_fragments/{{experiment}}.pooled.mono.bam.bai",
    params:
        binsize = BINSIZE,
        normalize_using = f"--normalizeUsing {v}" if (v := _tool_param("bamcoverage", "normalize_using", "CPM")) != "" else "--normalizeUsing CPM",
        smooth_length = f"--smoothLength {v}" if (v := _tool_param("bamcoverage", "smooth_length", "")) != "" else "",
        extra = _tool_param("bamcoverage", "extra_args", ""),
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.pooled.mnase_mono_bw.log",
    threads: THREADS,
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p $(dirname {output.bw:q}) $(dirname {log:q})

        bamCoverage \
            -b {input.bam:q} \
            -o {output.bw:q} \
            {params.normalize_using} \
            --binSize {params.binsize} \
            {params.smooth_length} \
            --numberOfProcessors {threads} \
            {params.extra} \
            2>&1 | tee {log:q}
        """
