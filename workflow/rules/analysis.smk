# analysis.smk — updated-profile derived analysis BAM
# ===================================================
# Legacy runs keep using final.bam directly. Updated runs keep a marked final
# BAM and derive analysis.bam for duplicate/organelle-filtered consumers.
# The template filter imports scripts/blacklist_filter.py


def _analysis_input(wildcards, kind):
    """Reject legacy execution even when an updated target is requested."""
    if not is_updated_profile():
        raise ValueError("analysis outputs are only available when analysis_profile=updated")
    suffix = "final.bam" if kind == "bam" else "final.bam.bai"
    return f"{RAW_OUTDIR}/{wildcards.sample}/02_align/{wildcards.sample}.{suffix}"


def _require_updated_input(wildcards):
    if not is_updated_profile():
        raise ValueError("analysis outputs are only available when analysis_profile=updated")
    return []


rule analysis_identity_marker:
    input:
        lambda wc: _require_updated_input(wc),
    output:
        f"{RAW_OUTDIR}/.helixweave/analysis_identity.json",
    params:
        identity=lambda wc: _analysis_identity_hash(),
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})"
        printf '{{"identity":"%s","analysis_profile":"updated"}}\n' \
            "{params.identity}" >{output:q}
        """


rule prepare_analysis_bam:
    input:
        identity=f"{RAW_OUTDIR}/.helixweave/analysis_identity.json",
        bam=lambda wc: _analysis_input(wc, "bam"),
        bai=lambda wc: _analysis_input(wc, "bai"),
        blacklist=lambda wc: (get_genome_resource(wc.sample, "blacklist") if ANALYSIS_FILTERED else []),
    output:
        bam=f"{OUTDIR}/{{sample}}/02_align/{{sample}}.analysis.bam",
        bai=f"{OUTDIR}/{{sample}}/02_align/{{sample}}.analysis.bam.bai",
        dropped=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.blacklist_dropped.tsv",
    conda:
        "../envs/samtools.yml"
    threads: THREADS
    params:
        exclude_dup=lambda wc: VALIDATED_CONFIG.get("duplicate_analysis") == "exclude",
        organelle_mode=lambda wc: get_organelle_mode(wc),
        organelle_chroms=lambda wc: " ".join(get_organelle_chroms(wc)),
        layout=lambda wc: SAMPLE_MAP[wc.sample]["layout"],
        blacklist_mode=BLACKLIST_MODE,
        script=f"{workflow.basedir}/../scripts/filter_blacklist_templates.py",
    shell:
        """
        set -e -o pipefail
        TMPD=$(mktemp -d)
        trap 'rm -rf "$TMPD"' EXIT
        TMP_BAM="$TMPD/analysis.bam"
        TMP_FILTERED="$TMPD/filtered.bam"
        FILTER_ARGS=()
        if [[ "{params.exclude_dup}" == "True" ]]; then
            FILTER_ARGS+=(-F 0x400)
        fi
        if [[ "{params.organelle_mode}" == "exclude" && -n "{params.organelle_chroms}" ]]; then
            samtools view -h -@ {threads} "${{FILTER_ARGS[@]}}" {input.bam:q} \
                | awk -v omit="{params.organelle_chroms}" \
                    'BEGIN {{ n=split(omit,a," "); for (i=1;i<=n;i++) drop[a[i]]=1 }}
                     /^@/ {{ print; next }}
                     !($3 in drop) {{ print }}' \
                | samtools view -b -@ {threads} -o "$TMP_BAM" -
        else
            samtools view -@ {threads} -b "${{FILTER_ARGS[@]}}" -o "$TMP_BAM" {input.bam:q}
        fi
        mkdir -p "$(dirname {output.dropped:q})"
        if [[ "{params.blacklist_mode}" == "analysis" ]]; then
            python3 {params.script:q} \
                --bam "$TMP_BAM" \
                --blacklist {input.blacklist:q} \
                --layout {params.layout:q} \
                --output "$TMP_FILTERED" \
                --dropped-count {output.dropped:q} \
                --sample {wildcards.sample:q}
        else
            mv "$TMP_BAM" "$TMP_FILTERED"
            printf 'sample\\tlayout\\tdropped_templates\\tdropped_alignments\\tkept_alignments\\ttotal_alignments\\n' \
                >{output.dropped:q}
            printf '%s\\t%s\\t0\\t0\\tNA\\tNA\\n' \
                {wildcards.sample:q} {params.layout:q} >>{output.dropped:q}
        fi
        mv "$TMP_FILTERED" {output.bam:q}
        samtools index -@ {threads} {output.bam:q}
        """


rule mito_fraction:
    input:
        identity=f"{RAW_OUTDIR}/.helixweave/analysis_identity.json",
        bam=lambda wc: _analysis_input(wc, "bam"),
    output:
        tsv=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.mito_fraction.tsv",
    conda:
        "../envs/samtools.yml"
    params:
        chroms=lambda wc: " ".join(get_organelle_chroms(wc)),
        mode=lambda wc: get_organelle_mode(wc),
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.tsv:q})"
        samtools idxstats {input.bam:q} \
            | awk -v omit="{params.chroms}" -v sample="{wildcards.sample}" -v mode="{params.mode}" \
                'BEGIN {{ n=split(omit,a," "); for (i=1;i<=n;i++) drop[a[i]]=1; total=0; organelle=0; print "sample\torganelle_reads\ttotal_mapped_reads\tfraction\tmode" }}
                 {{ total += $3; if ($1 in drop) organelle += $3 }}
                 END {{
                     frac = (total > 0) ? organelle / total : 0;
                     print sample "\t" organelle "\t" total "\t" sprintf("%.6f", frac) "\t" mode;
                 }}' \
                >{output.tsv:q}
        """
