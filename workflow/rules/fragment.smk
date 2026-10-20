# fragment.smk — additive fragment-level CPM tracks (PR-E1)
# ========================================================
# Legacy/read-level .CPM.bw paths are unchanged. Updated runs can add
# .fragment.CPM.bw tracks for ATAC, CUT&Tag, and MNase mono-nucleosome data.

_NON_MNASE_FRAGMENT_SAMPLES = [sid for sid in FRAGMENT_CPM_SAMPLE_IDS if SAMPLE_MAP[sid]["assay"] != "mnase"]
_MNASE_FRAGMENT_SAMPLES = [sid for sid in FRAGMENT_CPM_SAMPLE_IDS if SAMPLE_MAP[sid]["assay"] == "mnase"]


def _fragment_extend_args(wildcards):
    """One fragment per pair in PE; explicit SE length with 200 bp fallback."""
    if SAMPLE_MAP[wildcards.sample]["layout"] == "PE":
        return "--extendReads --samFlagInclude 64"
    ext = str(config.get("extend_reads", "auto"))
    if ext.isdigit():
        return f"--extendReads {ext}"
    legacy_args = get_extend_reads(wildcards)
    if legacy_args:
        return legacy_args
    return "--extendReads 200"


def _fragment_model_input(wildcards):
    """MACS3 model dependency for SE broad CUT&Tag auto/yes."""
    sample = SAMPLE_MAP[wildcards.sample]
    ext = str(config.get("extend_reads", "auto"))
    if (
        sample["layout"] == "SE"
        and sample["assay"] == "cuttag"
        and sample["peak_mode"] == "broad"
        and ext in ("auto", "yes")
    ):
        return [f"{OUTDIR}/{wildcards.sample}/04_peaks/{wildcards.sample}"]
    return []


rule fragment_cpm:
    input:
        bam=lambda wc: analysis_bam(wc.sample),
        bai=lambda wc: analysis_bai(wc.sample),
        model=lambda wc: _fragment_model_input(wc),
    output:
        bw=f"{OUTDIR}/{{sample}}/03_bigwig/{{sample}}.fragment.CPM.bw",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.fragment_cpm.log",
    wildcard_constraints:
        sample=wildcard_choices(_NON_MNASE_FRAGMENT_SAMPLES),
    conda:
        "../envs/deeptools.yml"
    threads: THREADS
    params:
        extend_args=_fragment_extend_args,
        binsize=BINSIZE,
        extra=_tool_param("bamcoverage", "extra_args", ""),
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.bw:q})" "$(dirname {log:q})"
        bamCoverage \
            --bam {input.bam:q} \
            --outFileName {output.bw:q} \
            --outFileFormat bigwig \
            --normalizeUsing CPM \
            --binSize {params.binsize} \
            {params.extend_args} {params.extra} \
            2>&1 | tee {log:q}
        """


rule mnase_fragment_cpm:
    input:
        bam=lambda wc: mnase_fragment_bam(wc.sample, "mono"),
        bai=lambda wc: mnase_fragment_bai(wc.sample, "mono"),
    output:
        bw=f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.mono.fragment.CPM.bw",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.mnase_fragment_cpm.log",
    wildcard_constraints:
        sample=wildcard_choices(_MNASE_FRAGMENT_SAMPLES),
    conda:
        "../envs/deeptools.yml"
    threads: THREADS
    params:
        extend_args=_fragment_extend_args,
        binsize=BINSIZE,
        extra=_tool_param("bamcoverage", "extra_args", ""),
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.bw:q})" "$(dirname {log:q})"
        bamCoverage \
            --bam {input.bam:q} \
            --outFileName {output.bw:q} \
            --outFileFormat bigwig \
            --normalizeUsing CPM \
            --binSize {params.binsize} \
            {params.extend_args} {params.extra} \
            2>&1 | tee {log:q}
        """
