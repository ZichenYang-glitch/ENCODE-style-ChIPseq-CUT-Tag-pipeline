# analysis.smk — updated-profile derived analysis BAM
# ===================================================
# Legacy runs keep using final.bam directly. Updated runs keep a marked final
# BAM and derive analysis.bam for duplicate-excluded analysis consumers.

def _analysis_input(wildcards, kind):
    """Reject legacy execution even when the analysis target is requested."""
    if not is_updated_profile():
        raise ValueError(
            "analysis.bam is only available when analysis_profile=updated"
        )
    suffix = "final.bam" if kind == "bam" else "final.bam.bai"
    return f"{OUTDIR}/{wildcards.sample}/02_align/{wildcards.sample}.{suffix}"


rule prepare_analysis_bam:
    input:
        bam=lambda wc: _analysis_input(wc, "bam"),
        bai=lambda wc: _analysis_input(wc, "bai"),
    output:
        bam=f"{OUTDIR}/{{sample}}/02_align/{{sample}}.analysis.bam",
        bai=f"{OUTDIR}/{{sample}}/02_align/{{sample}}.analysis.bam.bai",
    conda:
        "../envs/samtools.yml"
    threads: THREADS
    params:
        exclude_dup=lambda wc: VALIDATED_CONFIG.get("duplicate_analysis") == "exclude",
    shell:
        """
        set -e -o pipefail
        if [[ "{params.exclude_dup}" == "True" ]]; then
            samtools view -@ {threads} -b -F 0x400 {input.bam:q} > {output.bam:q}
        else
            ln -sf "$(basename {input.bam:q})" {output.bam:q}
        fi
        samtools index -@ {threads} {output.bam:q}
        """
