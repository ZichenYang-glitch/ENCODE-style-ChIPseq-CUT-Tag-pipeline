# common.smk — Shared workflow mechanics for ChIP-seq and CUT&Tag
# ================================================================
# These rules run for ALL samples (treatment + control).
# Assay-specific policy is dispatched through params functions
# (get_remove_dup, get_extend_reads) defined in the Snakefile
# dispatch layer, which delegates to chipseq.smk / cuttag.smk.
# MAPQ tag is derived from config (MAPQ_TAG = f"mapq{MAPQ}").


# ---------------------------------------------------------------------------
# 1. FastQC — sentinel-based output, layout-aware input list
# ---------------------------------------------------------------------------

def _fastqc_inputs(wildcards):
    """Return [fq1] for SE, [fq1, fq2] for PE.  Never returns empty string."""
    s = SAMPLE_MAP[wildcards.sample]
    if s["layout"] == "PE" and s.get("fq2"):
        return [s["fq1"], s["fq2"]]
    return [s["fq1"]]


rule fastqc:
    output:
        done = f"{OUTDIR}/{{sample}}/logs/{{sample}}.fastqc.done",
    input:
        lambda wc: _fastqc_inputs(wc),
    params:
        qcdir = lambda wc: f"{OUTDIR}/{wc.sample}/01_qc",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        mkdir -p {params.qcdir}
        set -- {input:q}
        fastqc -t {threads} -o {params.qcdir:q} "$@"
        touch {output.done}
        """


# ---------------------------------------------------------------------------
# 2. Trim Galore — layout-aware input list + report preservation
# ---------------------------------------------------------------------------

def _trim_galore_inputs(wildcards):
    """Return [fq1] for SE, [fq1, fq2] for PE.  Never returns empty string."""
    s = SAMPLE_MAP[wildcards.sample]
    if s["layout"] == "PE" and s.get("fq2"):
        return [s["fq1"], s["fq2"]]
    return [s["fq1"]]


rule trim_galore:
    output:
        r1   = f"{OUTDIR}/{{sample}}/00_raw/{{sample}}_R1_val_1.fq.gz",
        r2   = f"{OUTDIR}/{{sample}}/00_raw/{{sample}}_R2_val_2.fq.gz",
        done = f"{OUTDIR}/{{sample}}/logs/{{sample}}.trim.done",
    input:
        lambda wc: _trim_galore_inputs(wc),
    params:
        layout  = lambda wc: SAMPLE_MAP[wc.sample]["layout"],
        do_trim = TRIM,
        trimdir = lambda wc: f"{OUTDIR}/{wc.sample}/01_qc/trim_galore",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -- {input:q}
        # $1 = R1, $2 = R2 (PE only)

        if [[ "{params.do_trim}" == "true" ]]; then
            TMPD=$(mktemp -d)
            trap 'rm -rf "$TMPD"' EXIT
            mkdir -p {params.trimdir:q}

            if [[ "{params.layout}" == "PE" ]]; then
                trim_galore --paired --cores {threads} -o "$TMPD" "$1" "$2"
                mv "$TMPD"/*_val_1.fq.gz {output.r1:q}
                mv "$TMPD"/*_val_2.fq.gz {output.r2:q}
            else
                trim_galore --cores {threads} -o "$TMPD" "$1"
                mv "$TMPD"/*_trimmed.fq.gz {output.r1:q}
                touch {output.r2:q}
            fi

            # Preserve trimming reports for MultiQC
            mv "$TMPD"/*_trimming_report.txt {params.trimdir:q}/ 2>/dev/null || true
            mv "$TMPD"/*_trimming_report.json {params.trimdir:q}/ 2>/dev/null || true

            trap - EXIT
            rm -rf "$TMPD"
        else
            # No trimming: symlink raw FASTQs to expected trimmed paths
            ln -sf "$(readlink -f "$1")" {output.r1:q}
            if [[ "{params.layout}" == "PE" ]]; then
                ln -sf "$(readlink -f "$2")" {output.r2:q}
            else
                touch {output.r2:q}
            fi
            mkdir -p {params.trimdir:q}
        fi
        touch {output.done:q}
        """


# ---------------------------------------------------------------------------
# 3. Bowtie2 alignment + samtools sort (SE/PE dispatch in shell)
# ---------------------------------------------------------------------------

rule bowtie2_align:
    output:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.sorted.bam",
    input:
        r1 = f"{OUTDIR}/{{sample}}/00_raw/{{sample}}_R1_val_1.fq.gz",
        r2 = f"{OUTDIR}/{{sample}}/00_raw/{{sample}}_R2_val_2.fq.gz",
    params:
        layout  = lambda wc: SAMPLE_MAP[wc.sample]["layout"],
        bt2_idx = lambda wc: SAMPLE_MAP[wc.sample]["bt2_idx"],
        sample  = "{sample}",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.bowtie2.log",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        RG="--rg-id {params.sample} --rg SM:{params.sample} --rg LB:{params.sample} --rg PL:ILLUMINA"

        if [[ "{params.layout}" == "PE" ]]; then
            bowtie2 $RG -x {params.bt2_idx:q} -1 {input.r1:q} -2 {input.r2:q} -p {threads} \
                2> {log:q} | samtools view -b | samtools sort -@ {threads} -o {output.bam:q}
        else
            bowtie2 $RG -x {params.bt2_idx:q} -U {input.r1:q} -p {threads} \
                2> {log:q} | samtools view -b | samtools sort -@ {threads} -o {output.bam:q}
        fi
        """


# ---------------------------------------------------------------------------
# 4. Samtools index — sorted BAM
# ---------------------------------------------------------------------------

rule samtools_index_sorted:
    output:
        bai = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.sorted.bam.bai",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.sorted.bam",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        "samtools index -@ {threads} {input.bam:q}"


# ---------------------------------------------------------------------------
# 5. MAPQ filter — remove low-quality, secondary, supplementary reads
# ---------------------------------------------------------------------------

rule samtools_filter:
    output:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.{MAPQ_TAG}.bam",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.sorted.bam",
    params:
        mapq = MAPQ,
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        samtools view -@ {threads} -b -q {params.mapq} -F 0x904 {input.bam:q} > {output.bam:q}
        """


# ---------------------------------------------------------------------------
# 6. Samtools index — filtered BAM
# ---------------------------------------------------------------------------

rule samtools_index_filt:
    output:
        bai = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.{MAPQ_TAG}.bam.bai",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.{MAPQ_TAG}.bam",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        "samtools index -@ {threads} {input.bam:q}"


# ---------------------------------------------------------------------------
# 7. Duplicate handling → final.bam contract
# ---------------------------------------------------------------------------

rule duplicate_handling:
    output:
        bam     = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
        bai     = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam.bai",
        metrics = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.dup_metrics.txt",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.{MAPQ_TAG}.bam",
        bai = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.{MAPQ_TAG}.bam.bai",
    params:
        remove_dup = lambda wc: get_remove_dup(wc),
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        if [[ "{params.remove_dup}" == "yes" ]]; then
            if command -v picard >/dev/null 2>&1; then
                picard MarkDuplicates \
                    I={input.bam:q} \
                    O={output.bam:q} \
                    M={output.metrics:q} \
                    REMOVE_DUPLICATES=true
            else
                # samtools markdup fallback
                OUT_BAM_DIR=$(dirname {output.bam:q})
                TMP_NAME="$OUT_BAM_DIR/{wildcards.sample}.name.bam"
                TMP_FIX="$OUT_BAM_DIR/{wildcards.sample}.fixmate.bam"
                TMP_POS="$OUT_BAM_DIR/{wildcards.sample}.pos.bam"
                samtools sort -n -@ {threads} -o "$TMP_NAME" {input.bam:q}
                samtools fixmate -m -@ {threads} "$TMP_NAME" "$TMP_FIX"
                samtools sort -@ {threads} -o "$TMP_POS" "$TMP_FIX"
                samtools markdup -r -@ {threads} "$TMP_POS" {output.bam:q}
                rm -f "$TMP_NAME" "$TMP_FIX" "$TMP_POS"
                echo "samtools markdup used (picard not found)." > {output.metrics}
            fi
            samtools index -@ {threads} {output.bam}
        else
            if command -v picard >/dev/null 2>&1; then
                picard MarkDuplicates \
                    I={input.bam:q} \
                    O=/dev/null \
                    M={output.metrics:q} \
                    REMOVE_DUPLICATES=false \
                    >/dev/null 2>&1 || true
            else
                echo "picard not found; duplicate metrics unavailable." > {output.metrics}
            fi
            ln -sf $(basename {input.bam}) {output.bam}
            ln -sf $(basename {input.bai}) {output.bai}
        fi
        """


# ---------------------------------------------------------------------------
# 8. Flagstat — sorted BAM (pre-filter QC)
# ---------------------------------------------------------------------------

rule samtools_flagstat:
    output:
        flagstat = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.flagstat.txt",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.sorted.bam",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        "samtools flagstat -@ {threads} {input.bam:q} > {output.flagstat:q}"


# ---------------------------------------------------------------------------
# 9. Flagstat — final BAM (post-filter, post-dedup QC)
# ---------------------------------------------------------------------------

rule samtools_flagstat_final:
    output:
        flagstat = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.final.flagstat.txt",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        "samtools flagstat -@ {threads} {input.bam:q} > {output.flagstat:q}"


# ---------------------------------------------------------------------------
# 10. Idxstats — chromosomal read distribution
# ---------------------------------------------------------------------------

rule samtools_idxstats:
    output:
        idxstats = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.idxstats.txt",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.sorted.bam",
        bai = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.sorted.bam.bai",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        "samtools idxstats {input.bam:q} > {output.idxstats:q}"


# ---------------------------------------------------------------------------
# 11. bigWig generation (deepTools bamCoverage)
# ---------------------------------------------------------------------------

def _bamcoverage_inputs(wildcards):
    """Return input list for bamcoverage: [final.bam, final.bam.bai, ...optional peaks].

    For SE ChIP-seq with extend_reads=auto/yes, the MACS3 peaks output is
    added as a dependency so the MACS3 log exists for fragment-size extraction.
    """
    s = SAMPLE_MAP[wildcards.sample]
    inputs = [
        f"{OUTDIR}/{wildcards.sample}/02_align/{wildcards.sample}.final.bam",
        f"{OUTDIR}/{wildcards.sample}/02_align/{wildcards.sample}.final.bam.bai",
    ]
    ext = str(config.get("extend_reads", "auto"))
    if (
        s["role"] == "treatment"
        and s["layout"] == "SE"
        and s["assay"] == "chipseq"
        and ext in ("auto", "yes")
    ):
        inputs.append(
            f"{OUTDIR}/{wildcards.sample}/04_peaks/{wildcards.sample}"
        )
    return inputs


rule bamcoverage:
    output:
        bw = f"{OUTDIR}/{{sample}}/03_bigwig/{{sample}}.CPM.bw",
    input:
        lambda wc: _bamcoverage_inputs(wc),
    params:
        ext_args = lambda wc: get_extend_reads(wc),
        binsize  = BINSIZE,
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.bamCoverage.log",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        set -- {input:q}
        BAM="$1"
        bamCoverage \
            -b "$BAM" \
            -o {output.bw:q} \
            --normalizeUsing CPM \
            --binSize {params.binsize} \
            {params.ext_args} \
            --numberOfProcessors {threads} \
            2>&1 | tee {log:q}
        """
