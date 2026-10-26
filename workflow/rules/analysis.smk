# analysis.smk — run identity marker and canonical analysis BAM view
# =============================================================================
# The identity marker/effective-config rule runs in both profiles. Both profiles
# produce a canonical analysis.bam derived from final.bam to enforce resolved
# duplicate, organelle, and optional blacklist filtering policies.


def _analysis_input(wildcards, kind):
    """Return the input final BAM/BAI path for analysis view derivation."""
    suffix = "final.bam" if kind == "bam" else "final.bam.bai"
    return f"{RAW_OUTDIR}/{wildcards.sample}/02_align/{wildcards.sample}.{suffix}"


def _write_run_records(metadata_dir, records):
    """Write run records atomically, record first, marker last.

    Each record is written to an exclusively created temporary regular file
    in the verified metadata directory and atomically replaced into place,
    so pre-placed symlinks can never redirect the write and a valid marker
    never publishes without its record.
    """
    import os
    import tempfile

    os.makedirs(metadata_dir, exist_ok=True)
    for target, payload in records:
        fd, tmp_path = tempfile.mkstemp(
            dir=metadata_dir,
            prefix=".tmp-",
            suffix=f"-{os.path.basename(target)}",
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as handle:
                handle.write(f"{payload}\n")
            os.replace(tmp_path, target)
        except BaseException:
            if os.path.lexists(tmp_path):
                os.remove(tmp_path)
            raise


rule analysis_identity_marker:
    """Record the versioned run identity and the merged effective config.

    Scheduled in both profiles as a production-rule prerequisite; dry-runs
    never execute it, so dry-runs leave no run records behind.
    """
    output:
        marker=f"{RAW_OUTDIR}/.helixweave/analysis_identity.json",
        effective_config=f"{RAW_OUTDIR}/.helixweave/effective_config.json",
    log:
        f"{RAW_OUTDIR}/.helixweave/analysis_identity.log",
    params:
        marker_json=lambda wc: json.dumps(_analysis_identity_marker_payload(), indent=2, sort_keys=True),
        config_json=lambda wc: json.dumps(_effective_config_record_payload(), indent=2, sort_keys=True),
    run:
        _write_run_records(
            os.path.dirname(output.marker),
            ((output.effective_config, params.config_json), (output.marker, params.marker_json)),
        )


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
        stage_metrics=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.blacklist_dropped.stage_metrics.json",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.analysis_bam.log",
    conda:
        "../envs/samtools.yml"
    threads: THREADS
    params:
        exclude_dup=lambda wc: get_duplicate_policy(wc.sample) == "exclude",
        organelle_mode=lambda wc: get_organelle_mode(wc),
        organelle_chroms=lambda wc: " ".join(get_organelle_chroms(wc)),
        organelle_chroms_explicit=lambda wc: "true" if get_organelle_chroms_explicit(wc.sample) else "false",
        chrom_sizes=lambda wc: get_genome_resource(wc.sample, "chrom_sizes"),
        fai=lambda wc: (
            get_genome_resource(wc.sample, "reference_fasta") + ".fai"
            if (
                has_genome_resource(wc.sample, "reference_fasta")
                and os.path.isfile(get_genome_resource(wc.sample, "reference_fasta") + ".fai")
            )
            else ""
        ),
        layout=lambda wc: SAMPLE_MAP[wc.sample]["layout"],
        blacklist_mode=BLACKLIST_MODE,
        script=f"{workflow.basedir}/../scripts/filter_blacklist_templates.py",
        blacklist_module=f"{workflow.basedir}/../scripts/blacklist_filter.py",
    shell:
        """
        set -e -o pipefail
        EXTRA_ARGS=()
        if [[ "{params.exclude_dup}" == "True" ]]; then
            EXTRA_ARGS+=(--exclude-dup)
        fi
        if [[ "{params.organelle_chroms_explicit}" == "true" ]]; then
            EXTRA_ARGS+=(--organelle-chroms-explicit)
        fi
        if [[ -n "{params.organelle_chroms}" ]]; then
            EXTRA_ARGS+=(--organelle-chroms {params.organelle_chroms})
        fi
        if [[ -n "{params.chrom_sizes}" ]]; then
            EXTRA_ARGS+=(--chrom-sizes {params.chrom_sizes:q})
        fi
        if [[ -n "{params.fai}" ]]; then
            EXTRA_ARGS+=(--fai {params.fai:q})
        fi
        if [[ "{params.blacklist_mode}" == "analysis" && -n "{input.blacklist}" ]]; then
            EXTRA_ARGS+=(--blacklist {input.blacklist:q})
        fi
        mkdir -p "$(dirname {output.bam:q})" "$(dirname {output.dropped:q})" "$(dirname {log:q})"
        python3 {params.script:q} \
            --bam {input.bam:q} \
            --layout {params.layout:q} \
            --output {output.bam:q} \
            --dropped-count {output.dropped:q} \
            --stage-metrics {output.stage_metrics:q} \
            --sample {wildcards.sample:q} \
            --organelle-mode {params.organelle_mode:q} \
            --threads {threads} \
            "${{EXTRA_ARGS[@]}}" \
            2>&1 | tee {log:q}
        samtools index -@ {threads} {output.bam:q}
        """


rule analysis_template_metrics:
    """Count complete PE templates in the canonical analysis BAM.

    Runs for every sample independent of CUT&Tag optional QC.  SE samples are
    written explicitly as not applicable with NA counts.
    """
    input:
        bam=lambda wc: analysis_bam(wc.sample),
        stage_metrics=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.blacklist_dropped.stage_metrics.json",
    output:
        tsv=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.analysis_template_metrics.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.analysis_template_metrics.log",
    conda:
        "../envs/samtools.yml"
    params:
        script=f"{workflow.basedir}/../scripts/calc_analysis_template_metrics.py",
        layout=lambda wc: SAMPLE_MAP[wc.sample]["layout"],
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.tsv:q})" "$(dirname {log:q})"
        python3 {params.script:q} \
            --sample {wildcards.sample:q} \
            --layout {params.layout:q} \
            --bam {input.bam:q} \
            --stage-metrics {input.stage_metrics:q} \
            --output {output.tsv:q} \
            2>&1 | tee {log:q}
        """


rule mito_fraction:
    input:
        identity=f"{RAW_OUTDIR}/.helixweave/analysis_identity.json",
        bam=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.{MAPQ_TAG}.bam",
        bai=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.{MAPQ_TAG}.bam.bai",
    output:
        tsv=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.mito_fraction.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.mito_fraction.log",
    conda:
        "../envs/samtools.yml"
    threads: 1
    params:
        chroms=lambda wc: " ".join(get_organelle_chroms(wc)),
        is_explicit=lambda wc: "true" if get_organelle_chroms_explicit(wc.sample) else "false",
        mode=lambda wc: get_organelle_mode(wc),
        script=f"{workflow.basedir}/../scripts/calc_organelle_metrics.py",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.tsv:q})"
        EXTRA_ARGS=()
        if [[ "{params.is_explicit}" == "true" ]]; then
            EXTRA_ARGS+=(--is-explicit)
        fi
        if [[ -n "{params.chroms}" ]]; then
            EXTRA_ARGS+=(--chroms {params.chroms})
        fi
        python3 {params.script:q} \
            --bam {input.bam:q} \
            --sample {wildcards.sample:q} \
            --mode {params.mode:q} \
            --output {output.tsv:q} \
            "${{EXTRA_ARGS[@]}}"
        """


rule analysis_policy:
    """Record the resolved policy, input stages, and stage counts."""
    input:
        identity=f"{RAW_OUTDIR}/.helixweave/analysis_identity.json",
        stage_metrics=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.blacklist_dropped.stage_metrics.json",
        mito=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.mito_fraction.tsv",
        analysis_template_metrics=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.analysis_template_metrics.tsv",
        fragment_qc=lambda wc: (
            [f"{OUTDIR}/{wc.sample}/01_qc/{wc.sample}.cuttag_fragment_size.tsv"]
            if (QC_CONFIG.get("cuttag_fragment_size", True) and SAMPLE_MAP[wc.sample]["assay"] == "cuttag")
            else []
        ),
        mnase_qc=lambda wc: ([mnase_qc_summary_tsv(wc.sample)] if SAMPLE_MAP[wc.sample]["assay"] == "mnase" else []),
    output:
        f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.analysis_policy.json",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.analysis_policy.log",
    conda:
        "../envs/python.yml"
    params:
        script=f"{workflow.basedir}/../scripts/write_analysis_policy.py",
        sample="{sample}",
        assay=lambda wc: SAMPLE_MAP[wc.sample]["assay"],
        genome=lambda wc: SAMPLE_MAP[wc.sample]["genome"],
        profile=lambda wc: VALIDATED_CONFIG.get("analysis_profile"),
        policy_revision=SCIENTIFIC_POLICY_REVISION,
        duplicate_analysis=lambda wc: get_duplicate_policy(wc.sample),
        max_insert=lambda wc: get_max_insert(wc.sample),
        organelle_mode=lambda wc: get_organelle_mode(wc),
        blacklist_mode=BLACKLIST_MODE,
        blacklist_applied=lambda wc: (
            "true" if (ANALYSIS_FILTERED and has_genome_resource(wc.sample, "blacklist")) else "false"
        ),
        fragment_qc=lambda wc: (
            f"{OUTDIR}/{wc.sample}/01_qc/{wc.sample}.cuttag_fragment_size.tsv"
            if (QC_CONFIG.get("cuttag_fragment_size", True) and SAMPLE_MAP[wc.sample]["assay"] == "cuttag")
            else ""
        ),
        mnase_qc=lambda wc: (mnase_qc_summary_tsv(wc.sample) if SAMPLE_MAP[wc.sample]["assay"] == "mnase" else ""),
        analysis_template_metrics=lambda wc: (
            f"{OUTDIR}/{wc.sample}/05_qc/{wc.sample}.analysis_template_metrics.tsv"
        ),
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        FRAGMENT_QC_PATH={params.fragment_qc:q}
        FRAGMENT_QC_ARGS=()
        if [[ -n "$FRAGMENT_QC_PATH" ]]; then
            FRAGMENT_QC_ARGS+=(--fragment-qc "$FRAGMENT_QC_PATH")
        fi
        MNASE_QC_PATH={params.mnase_qc:q}
        MNASE_QC_ARGS=()
        if [[ -n "$MNASE_QC_PATH" ]]; then
            MNASE_QC_ARGS+=(--mnase-qc "$MNASE_QC_PATH")
        fi
        python3 {params.script:q} \
            --sample {params.sample:q} \
            --assay {params.assay:q} \
            --genome {params.genome:q} \
            --profile {params.profile:q} \
            --policy-revision {params.policy_revision} \
            --duplicate-analysis {params.duplicate_analysis:q} \
            --max-insert {params.max_insert:q} \
            --organelle-mode {params.organelle_mode:q} \
            --blacklist-mode {params.blacklist_mode:q} \
            --blacklist-applied {params.blacklist_applied:q} \
            --stage-metrics {input.stage_metrics:q} \
            --mito {input.mito:q} \
            --analysis-template-metrics {params.analysis_template_metrics:q} \
            "${{FRAGMENT_QC_ARGS[@]}}" \
            "${{MNASE_QC_ARGS[@]}}" \
            --output {output:q} \
            2>&1 | tee {log:q}
        """
