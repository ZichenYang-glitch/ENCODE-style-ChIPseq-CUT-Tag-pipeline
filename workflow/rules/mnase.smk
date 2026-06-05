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
# MNase config helpers (Stage 39 + Stage 40)
# ---------------------------------------------------------------------------

def get_mono_range():
    """Return (min, max) mono-nucleosome fragment range from config.

    Priority: mnase.fragments.mono > mnase.mono_range > [140, 200].
    """
    mnase_cfg = config.get("mnase", {})
    if isinstance(mnase_cfg, dict):
        fragments = mnase_cfg.get("fragments", {})
        if isinstance(fragments, dict) and "mono" in fragments:
            r = fragments["mono"]
            if isinstance(r, (list, tuple)) and len(r) == 2:
                return int(r[0]), int(r[1])
        mr = mnase_cfg.get("mono_range", [140, 200])
        if isinstance(mr, (list, tuple)) and len(mr) == 2:
            return int(mr[0]), int(mr[1])
    return 140, 200


def get_fragment_range(class_name):
    """Return (min, max) for sub/mono/di. Never returns None.

    class_name: "sub" | "mono" | "di"
    Defaults: sub [1,139], mono [140,200], di [300,400]
    """
    DEFAULTS = {"sub": (1, 139), "mono": (140, 200), "di": (300, 400)}
    if class_name == "mono":
        return get_mono_range()
    mnase_cfg = config.get("mnase", {})
    if isinstance(mnase_cfg, dict):
        fragments = mnase_cfg.get("fragments", {})
        if isinstance(fragments, dict) and class_name in fragments:
            r = fragments[class_name]
            if isinstance(r, (list, tuple)) and len(r) == 2:
                return int(r[0]), int(r[1])
    return DEFAULTS[class_name]


def get_dyad_range():
    """Return (min, max) for dyad BigWig --MNase fragment length.

    Default [130, 200] preserves Stage 39 / deepTools --MNase default.
    """
    mnase_cfg = config.get("mnase", {})
    if isinstance(mnase_cfg, dict):
        dr = mnase_cfg.get("dyad_range", [130, 200])
        if isinstance(dr, (list, tuple)) and len(dr) == 2:
            return int(dr[0]), int(dr[1])
    return 130, 200


def _caller_enabled(name):
    """Return True if caller *name* is enabled in mnase.callers.

    Validation requires YAML boolean values (true/false), so only
    isinstance(val, bool) is checked here. String values are rejected
    by validate_samples.py; this helper is defense-in-depth.
    """
    mnase_cfg = config.get("mnase", {})
    if isinstance(mnase_cfg, dict):
        callers = mnase_cfg.get("callers", {})
        if isinstance(callers, dict):
            val = callers.get(name, False)
            if isinstance(val, bool):
                return val
    return False


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
# 2. Sub-nucleosome BAM — alignmentSieve from final.bam (Stage 40)
# ---------------------------------------------------------------------------

rule mnase_split_sub:
    output:
        bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.sub.bam",
        bai = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.sub.bam.bai",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
        bai = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam.bai",
    params:
        sub_min = get_fragment_range("sub")[0],
        sub_max = get_fragment_range("sub")[1],
    threads: THREADS,
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p $(dirname {output.bam:q})

        alignmentSieve \
            -b {input.bam:q} \
            --minFragmentLength {params.sub_min} \
            --maxFragmentLength {params.sub_max} \
            -p {threads} \
            -o {output.bam:q}

        samtools index -@ {threads} {output.bam:q}
        """


# ---------------------------------------------------------------------------
# 3. Di-nucleosome BAM — alignmentSieve from final.bam (Stage 40)
# ---------------------------------------------------------------------------

rule mnase_split_di:
    output:
        bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.di.bam",
        bai = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.di.bam.bai",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
        bai = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam.bai",
    params:
        di_min = get_fragment_range("di")[0],
        di_max = get_fragment_range("di")[1],
    threads: THREADS,
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p $(dirname {output.bam:q})

        alignmentSieve \
            -b {input.bam:q} \
            --minFragmentLength {params.di_min} \
            --maxFragmentLength {params.di_max} \
            -p {threads} \
            -o {output.bam:q}

        samtools index -@ {threads} {output.bam:q}
        """


# ---------------------------------------------------------------------------
# 4. Dyad BigWig — bamCoverage --MNase --binSize 1 from final.bam
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
        dyad_min = get_dyad_range()[0],
        dyad_max = get_dyad_range()[1],
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.mnase_dyad_bw.log",
    threads: THREADS,
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p $(dirname {output.bw:q}) $(dirname {log:q})

        # tool_parameters.bamcoverage.extra_args must NOT override
        # --MNase, --minFragmentLength, or --maxFragmentLength;
        # these are workflow-managed for MNase dyad BigWig.

        bamCoverage \
            -b {input.bam:q} \
            -o {output.bw:q} \
            --MNase \
            --binSize 1 \
            --minFragmentLength {params.dyad_min} \
            --maxFragmentLength {params.dyad_max} \
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
        dyad_min = get_dyad_range()[0],
        dyad_max = get_dyad_range()[1],
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.pooled.mnase_dyad_bw.log",
    threads: THREADS,
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p $(dirname {output.bw:q}) $(dirname {log:q})

        # tool_parameters.bamcoverage.extra_args must NOT override
        # --MNase, --minFragmentLength, or --maxFragmentLength;
        # these are workflow-managed for MNase dyad BigWig.

        bamCoverage \
            -b {input.bam:q} \
            -o {output.bw:q} \
            --MNase \
            --binSize 1 \
            --minFragmentLength {params.dyad_min} \
            --maxFragmentLength {params.dyad_max} \
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


# ---------------------------------------------------------------------------
# 7. MNase QC summary (Stage 40)
# ---------------------------------------------------------------------------

rule mnase_qc_summary:
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.mnase_qc_summary.tsv",
    input:
        sub_bam  = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.sub.bam",
        mono_bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam",
        di_bam   = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.di.bam",
        dyad_bw  = f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.dyad.CPM.bw",
        mono_bw  = f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.mono.CPM.bw",
    params:
        sample  = "{sample}",
        assay   = lambda wc: SAMPLE_MAP[wc.sample]["assay"],
        peak_mode = lambda wc: SAMPLE_MAP[wc.sample]["peak_mode"],
        sub_min  = get_fragment_range("sub")[0],
        sub_max  = get_fragment_range("sub")[1],
        mono_min = get_fragment_range("mono")[0],
        mono_max = get_fragment_range("mono")[1],
        di_min   = get_fragment_range("di")[0],
        di_max   = get_fragment_range("di")[1],
        dyad_min = get_dyad_range()[0],
        dyad_max = get_dyad_range()[1],
        insert_size_metrics = f"{OUTDIR}/{{sample}}/05_qc/picard/{{sample}}.insert_size_metrics",
        danpos3_enabled = str(_caller_enabled("danpos3")).lower(),
        inps_enabled    = str(_caller_enabled("inps")).lower(),
        sem_enabled     = str(_caller_enabled("sem")).lower(),
    conda:
        "../envs/deeptools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})"
        python3 scripts/mnase_qc_summary.py \
            --sample {params.sample:q} \
            --assay {params.assay:q} \
            --peak-mode {params.peak_mode:q} \
            --sub-min {params.sub_min} --sub-max {params.sub_max} \
            --mono-min {params.mono_min} --mono-max {params.mono_max} \
            --di-min {params.di_min} --di-max {params.di_max} \
            --dyad-min {params.dyad_min} --dyad-max {params.dyad_max} \
            --sub-bam {input.sub_bam:q} \
            --mono-bam {input.mono_bam:q} \
            --di-bam {input.di_bam:q} \
            --dyad-bw {input.dyad_bw:q} \
            --mono-bw {input.mono_bw:q} \
            --insert-size-metrics {params.insert_size_metrics:q} \
            --danpos3-enabled {params.danpos3_enabled} \
            --inps-enabled {params.inps_enabled} \
            --sem-enabled {params.sem_enabled} \
            --output {output:q}
        """
