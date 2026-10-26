# qc.smk — Single-Sample QC rules
# ========================================
# Shared versioned-QC runtime helpers: scripts/qc_common.py and
# scripts/qc_template.py are imported by the QC writers below and are part
# of the controlled runtime closure.
# Blacklist filtering (BAM + peaks), peak counts, FRiP, library complexity,
# MACS3 signal tracks, and QC summary.
# All QC features are resource-gated — blacklist operations only run when
# a blacklist BED is configured for the sample's genome in genome_resources.
#
# Rules:
#   blacklist_filter_bam    → bedtools intersect -v + samtools index
#   blacklist_filter_peaks  → bedtools intersect -v on peak file
#   peak_counts             → count lines in peak files
#   frip                    → calc_frip.py (read-record based)
#   library_complexity      → parse duplicate metrics when available
#   signal_track_fe         → MACS3 bdgcmp FE bedGraph
#   signal_track_ppois      → MACS3 bdgcmp ppois bedGraph
#   qc_summary              → per-sample TSV assembly
#   project_qc_summary       → project-level aggregated summary


# Helper: peak file suffix for a sample
def _peak_suffix(sample_id):
    return "broadPeak" if SAMPLE_MAP[sample_id]["peak_mode"] == "broad" else "narrowPeak"


def _peak_file(sample_id):
    """Return the MACS3 peak file path for a sample."""
    suffix = _peak_suffix(sample_id)
    return f"{OUTDIR}/{sample_id}/04_peaks/{sample_id}/" f"{sample_id}_peaks.{suffix}"


def _blacklist_output_root():
    """Output root for filtered peak copies: filtered branch or raw copies."""
    return OUTDIR if ANALYSIS_FILTERED else RAW_OUTDIR


def _has_blacklist_qc(sample_id):
    """Return True when blacklist filtering is enabled and resource-backed."""
    return sample_id in BLACKLIST_SAMPLES and QC_CONFIG.get("blacklist_filter", True)


def _blacklist_peak_file(sample_id):
    """Return the blacklist-filtered peak file path for a sample."""
    suffix = _peak_suffix(sample_id)
    return (
        f"{OUTDIR}/{sample_id}/04_peaks/{sample_id}_blacklist_filtered/"
        f"{sample_id}_peaks.blacklist_filtered.{suffix}"
    )


def _analysis_filtered_peak_file(sample_id):
    """Return the final filtered peak file in the analysis_filtered branch."""
    suffix = _peak_suffix(sample_id)
    return (
        f"{OUTDIR}/{sample_id}/04_peaks/{sample_id}_blacklist_filtered/"
        f"{sample_id}_peaks.blacklist_filtered.{suffix}"
    )


def _blacklist_bam_file(sample_id):
    """Return the blacklist-filtered BAM file path for a sample."""
    return f"{OUTDIR}/{sample_id}/02_align/{sample_id}.blacklist_filtered.bam"


def _signal_track_file(sample_id, signal):
    """Return a MACS3-derived signal bedGraph path for a sample."""
    return f"{OUTDIR}/{sample_id}/03_signal/{sample_id}.{signal}.bdg"


def _pooled_signal_track_file(experiment, signal):
    """Return a pooled experiment MACS3-derived signal bedGraph path."""
    return f"{OUTDIR}/experiments/{experiment}/03_signal/" f"{experiment}.pooled.{signal}.bdg"


# ---------------------------------------------------------------------------
# 1. Blacklist filter — BAM
# ---------------------------------------------------------------------------


rule blacklist_filter_bam:
    input:
        bam=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
        bai=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam.bai",
    output:
        bam=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.blacklist_filtered.bam",
        bai=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.blacklist_filtered.bam.bai",
    conda:
        "../envs/samtools.yml"
    threads: THREADS
    params:
        blacklist=lambda wc: get_genome_resource(wc.sample, "blacklist"),
        validate_script=f"{workflow.basedir}/../scripts/validate_blacklist.py",
    shell:
        """
        set -e -o pipefail
        python3 {params.validate_script:q} \
            --bam {input.bam:q} \
            --blacklist {params.blacklist:q} \
            --sample {wildcards.sample:q}
        bedtools intersect -v -abam {input.bam:q} -b {params.blacklist:q} \
            >{output.bam:q}
        samtools index -@ {threads} {output.bam:q}
        """


# ---------------------------------------------------------------------------
# 2. Blacklist filter — peaks
# ---------------------------------------------------------------------------


rule blacklist_filter_peaks:
    input:
        # MACS3 peaks directory (Snakemake-tracked) — peak file assembled in shell
        peaks_dir=f"{OUTDIR}/{{sample}}/04_peaks/{{sample}}",
        bam=lambda wc: analysis_bam(wc.sample),
    output:
        f"{_blacklist_output_root()}/{{sample}}/04_peaks/{{sample}}_blacklist_filtered/{{sample}}_peaks.blacklist_filtered.{{suffix}}",
    conda:
        "../envs/samtools.yml"
    params:
        blacklist=lambda wc: get_genome_resource(wc.sample, "blacklist"),
        validate_script=f"{workflow.basedir}/../scripts/validate_blacklist.py",
    shell:
        """
        set -e -o pipefail
        PEAKS_FILE="{input.peaks_dir:q}/{wildcards.sample}_peaks.{wildcards.suffix}"
        mkdir -p $(dirname {output:q})
        python3 {params.validate_script:q} \
            --bam {input.bam:q} \
            --blacklist {params.blacklist:q} \
            --sample {wildcards.sample:q}
        bedtools intersect -v -a "$PEAKS_FILE" -b {params.blacklist:q} \
            >{output:q}
        """


rule blacklist_filter_pooled_peaks:
    input:
        peaks_dir=f"{OUTDIR}/experiments/{{experiment}}/04_peaks/pooled/{{experiment}}_pooled_peaks",
        bam=f"{OUTDIR}/experiments/{{experiment}}/02_align/{{experiment}}.pooled.final.bam",
    output:
        f"{_blacklist_output_root()}/experiments/{{experiment}}/04_peaks/pooled/{{experiment}}_pooled_peaks_blacklist_filtered/{{experiment}}_pooled_peaks.blacklist_filtered.{{suffix}}",
    conda:
        "../envs/samtools.yml"
    params:
        blacklist=lambda wc: get_genome_resource(TREATMENT_SAMPLES_BY_EXPERIMENT[wc.experiment][0], "blacklist"),
        validate_script=f"{workflow.basedir}/../scripts/validate_blacklist.py",
    shell:
        """
        set -e -o pipefail
        PEAKS_FILE="{input.peaks_dir:q}/{wildcards.experiment}_pooled_peaks.{wildcards.suffix}"
        mkdir -p $(dirname {output:q})
        python3 {params.validate_script:q} \
            --bam {input.bam:q} \
            --blacklist {params.blacklist:q} \
            --sample {wildcards.experiment:q}
        bedtools intersect -v -a "$PEAKS_FILE" -b {params.blacklist:q} \
            >{output:q}
        """


# ---------------------------------------------------------------------------
# 3. Peak counts
# ---------------------------------------------------------------------------


rule peak_counts:
    input:
        lambda wc: _peak_counts_inputs(wc),
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.peak_counts.tsv",
    conda:
        "../envs/python.yml"
    params:
        peak_mode=lambda wc: SAMPLE_MAP[wc.sample]["peak_mode"],
        sample="{sample}",
        analysis_mode="yes" if ANALYSIS_FILTERED else "no",
        bl_peak=lambda wc: (
            _analysis_filtered_peak_file(wc.sample)
            if ANALYSIS_FILTERED
            else (_blacklist_peak_file(wc.sample) if _has_blacklist_qc(wc.sample) else "")
        ),
    shell:
        """
        set -e -o pipefail
        PEAKS_DIR="{input[0]:q}"
        SUFFIX="{params.peak_mode}"
        [[ "$SUFFIX" == "broad" ]] && SUFFIX="broadPeak" || SUFFIX="narrowPeak"
        PEAKS_FILE="$PEAKS_DIR/{params.sample}_peaks.$SUFFIX"
        PEAK_COUNT=$(wc -l <"$PEAKS_FILE")

        BL_PEAK={params.bl_peak:q}
        if [[ -n "$BL_PEAK" && -f "$BL_PEAK" ]]; then
            BL_COUNT=$(wc -l <"$BL_PEAK")
            if [[ "{params.analysis_mode}" == "yes" ]]; then
                PEAK_COUNT="$BL_COUNT"
            fi
        else
            BL_COUNT="NA"
        fi

        printf "sample\\tpeak_mode\\tpeaks\\tblacklist_filtered_peaks\\n" \
            >{output:q}
        printf "%s\\t%s\\t%s\\t%s\\n" \
            {params.sample:q} {params.peak_mode:q} \
            "$PEAK_COUNT" "$BL_COUNT" \
            >>{output:q}
        """


def _peak_counts_inputs(wildcards):
    """Return [peaks_dir] or [peaks_dir, blacklist_filtered_peak_file]."""
    sid = wildcards.sample
    peaks_dir = f"{OUTDIR}/{sid}/04_peaks/{sid}"
    inputs = [peaks_dir]
    if ANALYSIS_FILTERED or _has_blacklist_qc(sid):
        inputs.append(_analysis_filtered_peak_file(sid) if ANALYSIS_FILTERED else _blacklist_peak_file(sid))
    return inputs


# ---------------------------------------------------------------------------
# 4. FRiP (Fraction of Reads in Peaks)
# ---------------------------------------------------------------------------


rule frip:
    input:
        lambda wc: _frip_inputs(wc),
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.frip.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.frip.log",
    conda:
        "../envs/samtools.yml"
    params:
        script=f"{workflow.basedir}/../scripts/calc_frip.py",
        peak_mode=lambda wc: SAMPLE_MAP[wc.sample]["peak_mode"],
    shell:
        """
        set -e -o pipefail
        set -- {input:q}
        BAM="$1"
        SUFFIX="narrowPeak"
        [[ "{params.peak_mode}" == "broad" ]] && SUFFIX="broadPeak"

        # $2 is either a peaks_dir (directory) or a blacklist-filtered peak file
        if [[ -d "$2" ]]; then
            PEAKS="$2/{wildcards.sample}_peaks.$SUFFIX"
        else
            PEAKS="$2"
        fi

        python3 {params.script:q} \
            --sample {wildcards.sample:q} \
            --bam "$BAM" \
            --peaks "$PEAKS" \
            --output {output:q} \
            2>&1 | tee {log:q}
        """


def _frip_inputs(wildcards):
    """Return [bam, peaks_dir_or_file] for primary FRiP calculation.

    Always consumes the canonical analysis BAM. In the analysis-filtered branch,
    pairs with the secondary blacklist-filtered peaks; otherwise pairs with the
    primary peaks.
    """
    sid = wildcards.sample
    if ANALYSIS_FILTERED:
        return [analysis_bam(sid), _analysis_filtered_peak_file(sid)]
    return [analysis_bam(sid), f"{OUTDIR}/{sid}/04_peaks/{sid}"]


# ---------------------------------------------------------------------------
# 5. Library complexity (duplication metrics)
# ---------------------------------------------------------------------------


rule library_complexity:
    input:
        f"{RAW_OUTDIR}/{{sample}}/01_qc/{{sample}}.dup_metrics.txt",
    output:
        f"{RAW_OUTDIR}/{{sample}}/01_qc/{{sample}}.library_complexity.tsv",
    log:
        f"{RAW_OUTDIR}/{{sample}}/logs/{{sample}}.library_complexity.log",
    conda:
        "../envs/python.yml"
    params:
        script=f"{workflow.basedir}/../scripts/parse_dup_metrics.py",
    shell:
        """
        set -e -o pipefail
        python3 {params.script:q} \
            --sample {wildcards.sample:q} \
            --metrics {input:q} \
            --output {output:q} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 5a. NRF/PBC library complexity (BAM-derived metrics)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# 5c. SEACR sidecar peak calling (Stage 7b)
# ---------------------------------------------------------------------------


rule seacr_bedgraph:
    input:
        lambda wc: analysis_bam(wc.sample),
    output:
        f"{OUTDIR}/{{sample}}/04_peaks_seacr/{{sample}}.bedgraph",
    conda:
        "../envs/seacr.yml"
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})"
        bedtools genomecov -ibam {input:q} -bg -pc > {output:q}
        """


rule seacr_call:
    input:
        f"{OUTDIR}/{{sample}}/04_peaks_seacr/{{sample}}.bedgraph",
    output:
        f"{OUTDIR}/{{sample}}/04_peaks_seacr/{{sample}}/{{sample}}.seacr.{{mode}}.bed",
    conda:
        "../envs/seacr.yml"
    params:
        mode=lambda wc: wc.mode,
        threshold=SEACR_THRESHOLD,
        norm=SEACR_NORMALIZATION,
    shell:
        """
        set -e -o pipefail
        SEACR_SCRIPT=$(command -v SEACR_1.3.sh) || {{ \\
                    echo "ERROR: SEACR_1.3.sh not found in PATH" >&2; exit 1; }}
        mkdir -p "$(dirname {output:q})"
        bash "$SEACR_SCRIPT" {input:q} {params.threshold} \\
        {params.norm:q} {params.mode:q} \\
        "$(dirname {output:q})/{wildcards.sample}.seacr"
        """


# 5b. CUT&Tag fragment-size QC (Stage 7a)
# ---------------------------------------------------------------------------


rule cuttag_fragment_size:
    input:
        bam=lambda wc: analysis_bam(wc.sample),
        stage_metrics=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.blacklist_dropped.stage_metrics.json",
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.cuttag_fragment_size.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.cuttag_fragment_size.log",
    conda:
        "../envs/samtools.yml"
    params:
        script=f"{workflow.basedir}/../scripts/calc_cuttag_fragment_size.py",
        layout=lambda wc: SAMPLE_MAP[wc.sample]["layout"],
        duplicate_analysis=lambda wc: get_duplicate_policy(wc.sample),
        blacklist_applied=lambda wc: (
            "true" if (ANALYSIS_FILTERED and has_genome_resource(wc.sample, "blacklist")) else "false"
        ),
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        python3 {params.script:q} \
            --sample {wildcards.sample:q} \
            --bam {input.bam:q} \
            --layout {params.layout:q} \
            --duplicate-analysis {params.duplicate_analysis:q} \
            --blacklist-applied {params.blacklist_applied:q} \
            --stage-metrics {input.stage_metrics:q} \
            --output {output:q} \
            2>&1 | tee {log:q}
        """


rule nrf_pbc:
    input:
        f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.{MAPQ_TAG}.bam",
    output:
        f"{RAW_OUTDIR}/{{sample}}/01_qc/{{sample}}.nrf_pbc.tsv",
    log:
        f"{RAW_OUTDIR}/{{sample}}/logs/{{sample}}.nrf_pbc.log",
    conda:
        "../envs/samtools.yml"
    params:
        script=f"{workflow.basedir}/../scripts/calc_nrf_pbc.py",
    shell:
        """
        set -e -o pipefail
        python3 {params.script:q} \
            --sample {wildcards.sample:q} \
            --bam {input:q} \
            --output {output:q} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 6. MACS3 signal tracks (fold enrichment and p-value)
# ---------------------------------------------------------------------------


rule signal_track_fe:
    input:
        peaks_dir=f"{OUTDIR}/{{sample}}/04_peaks/{{sample}}",
    output:
        f"{OUTDIR}/{{sample}}/03_signal/{{sample}}.FE.bdg",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.bdgcmp.FE.log",
    conda:
        "../envs/macs3.yml"
    shell:
        """
        set -e -o pipefail
        TREAT={input.peaks_dir:q}/{wildcards.sample}_treat_pileup.bdg
        LAMBDA={input.peaks_dir:q}/{wildcards.sample}_control_lambda.bdg
        if [[ ! -s "$TREAT" || ! -s "$LAMBDA" ]]; then
            echo "ERROR: MACS3 bedGraph inputs missing for FE signal track." >&2
            echo "Expected: $TREAT" >&2
            echo "Expected: $LAMBDA" >&2
            echo "Ensure qc.signal_tracks is true so macs3 callpeak runs with -B." >&2
            exit 1
        fi

        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        macs3 bdgcmp \
            -t "$TREAT" \
            -c "$LAMBDA" \
            -m FE \
            -p 1 \
            -o {output:q} \
            2>&1 | tee {log:q}
        """


rule signal_track_ppois:
    input:
        peaks_dir=f"{OUTDIR}/{{sample}}/04_peaks/{{sample}}",
    output:
        f"{OUTDIR}/{{sample}}/03_signal/{{sample}}.ppois.bdg",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.bdgcmp.ppois.log",
    conda:
        "../envs/macs3.yml"
    shell:
        """
        set -e -o pipefail
        TREAT={input.peaks_dir:q}/{wildcards.sample}_treat_pileup.bdg
        LAMBDA={input.peaks_dir:q}/{wildcards.sample}_control_lambda.bdg
        if [[ ! -s "$TREAT" || ! -s "$LAMBDA" ]]; then
            echo "ERROR: MACS3 bedGraph inputs missing for ppois signal track." >&2
            echo "Expected: $TREAT" >&2
            echo "Expected: $LAMBDA" >&2
            echo "Ensure qc.signal_tracks is true so macs3 callpeak runs with -B." >&2
            exit 1
        fi

        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        macs3 bdgcmp \
            -t "$TREAT" \
            -c "$LAMBDA" \
            -m ppois \
            -o {output:q} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 6b. Single-sample FE/ppois BigWig conversion (Stage 22)
# ---------------------------------------------------------------------------


rule signal_track_fe_bw:
    input:
        f"{OUTDIR}/{{sample}}/03_signal/{{sample}}.FE.bdg",
    output:
        f"{OUTDIR}/{{sample}}/03_signal/{{sample}}.FE.bw",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.bdg_to_bw.FE.log",
    conda:
        "../envs/ucsc.yml"
    resources:
        bigwig_convert=1,
    params:
        chrom_sizes=lambda wc: get_genome_resource(wc.sample, "chrom_sizes"),
        sort_tmpdir="/tmp",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        SORTED="$(dirname {output:q})/.tmp_{wildcards.sample}.FE.sorted.bdg"
        trap 'rm -f "$SORTED"' EXIT
        LC_ALL=C sort --temporary-directory={params.sort_tmpdir:q} -k1,1 -k2,2n {input:q} >"$SORTED"
        bedGraphToBigWig "$SORTED" {params.chrom_sizes:q} {output:q} 2>&1 | tee {log:q}
        """


rule signal_track_ppois_bw:
    input:
        f"{OUTDIR}/{{sample}}/03_signal/{{sample}}.ppois.bdg",
    output:
        f"{OUTDIR}/{{sample}}/03_signal/{{sample}}.ppois.bw",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.bdg_to_bw.ppois.log",
    conda:
        "../envs/ucsc.yml"
    resources:
        bigwig_convert=1,
    params:
        chrom_sizes=lambda wc: get_genome_resource(wc.sample, "chrom_sizes"),
        sort_tmpdir="/tmp",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        SORTED="$(dirname {output:q})/.tmp_{wildcards.sample}.ppois.sorted.bdg"
        trap 'rm -f "$SORTED"' EXIT
        LC_ALL=C sort --temporary-directory={params.sort_tmpdir:q} -k1,1 -k2,2n {input:q} >"$SORTED"
        bedGraphToBigWig "$SORTED" {params.chrom_sizes:q} {output:q} 2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 6a. Pooled MACS3 signal tracks (Stage 6a)
# ---------------------------------------------------------------------------


rule pooled_signal_track_fe:
    input:
        peaks_dir=f"{OUTDIR}/experiments/{{experiment}}/04_peaks/pooled/{{experiment}}_pooled_peaks",
    output:
        f"{OUTDIR}/experiments/{{experiment}}/03_signal/{{experiment}}.pooled.FE.bdg",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.pooled.bdgcmp.FE.log",
    conda:
        "../envs/macs3.yml"
    shell:
        """
        set -e -o pipefail
        TREAT={input.peaks_dir:q}/{wildcards.experiment}_pooled_treat_pileup.bdg
        LAMBDA={input.peaks_dir:q}/{wildcards.experiment}_pooled_control_lambda.bdg
        if [[ ! -s "$TREAT" || ! -s "$LAMBDA" ]]; then
            echo "ERROR: pooled MACS3 bedGraph inputs missing for FE signal track." >&2
            echo "Expected: $TREAT" >&2
            echo "Expected: $LAMBDA" >&2
            echo "Ensure qc.signal_tracks is true so pooled MACS3 runs with -B." >&2
            exit 1
        fi

        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        macs3 bdgcmp \
            -t "$TREAT" \
            -c "$LAMBDA" \
            -m FE \
            -p 1 \
            -o {output:q} \
            2>&1 | tee {log:q}
        """


rule pooled_signal_track_ppois:
    input:
        peaks_dir=f"{OUTDIR}/experiments/{{experiment}}/04_peaks/pooled/{{experiment}}_pooled_peaks",
    output:
        f"{OUTDIR}/experiments/{{experiment}}/03_signal/{{experiment}}.pooled.ppois.bdg",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.pooled.bdgcmp.ppois.log",
    conda:
        "../envs/macs3.yml"
    shell:
        """
        set -e -o pipefail
        TREAT={input.peaks_dir:q}/{wildcards.experiment}_pooled_treat_pileup.bdg
        LAMBDA={input.peaks_dir:q}/{wildcards.experiment}_pooled_control_lambda.bdg
        if [[ ! -s "$TREAT" || ! -s "$LAMBDA" ]]; then
            echo "ERROR: pooled MACS3 bedGraph inputs missing for ppois signal track." >&2
            echo "Expected: $TREAT" >&2
            echo "Expected: $LAMBDA" >&2
            echo "Ensure qc.signal_tracks is true so pooled MACS3 runs with -B." >&2
            exit 1
        fi

        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        macs3 bdgcmp \
            -t "$TREAT" \
            -c "$LAMBDA" \
            -m ppois \
            -o {output:q} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 6c. Pooled FE/ppois BigWig conversion (Stage 22)
# ---------------------------------------------------------------------------


rule pooled_signal_track_fe_bw:
    input:
        f"{OUTDIR}/experiments/{{experiment}}/03_signal/{{experiment}}.pooled.FE.bdg",
    output:
        f"{OUTDIR}/experiments/{{experiment}}/03_signal/{{experiment}}.pooled.FE.bw",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.pooled.bdg_to_bw.FE.log",
    conda:
        "../envs/ucsc.yml"
    resources:
        bigwig_convert=1,
    params:
        chrom_sizes=lambda wc: _pooled_chrom_sizes(wc.experiment),
        sort_tmpdir="/tmp",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        SORTED="$(dirname {output:q})/.tmp_{wildcards.experiment}.pooled.FE.sorted.bdg"
        trap 'rm -f "$SORTED"' EXIT
        LC_ALL=C sort --temporary-directory={params.sort_tmpdir:q} -k1,1 -k2,2n {input:q} >"$SORTED"
        bedGraphToBigWig "$SORTED" {params.chrom_sizes:q} {output:q} 2>&1 | tee {log:q}
        """


rule pooled_signal_track_ppois_bw:
    input:
        f"{OUTDIR}/experiments/{{experiment}}/03_signal/{{experiment}}.pooled.ppois.bdg",
    output:
        f"{OUTDIR}/experiments/{{experiment}}/03_signal/{{experiment}}.pooled.ppois.bw",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.pooled.bdg_to_bw.ppois.log",
    conda:
        "../envs/ucsc.yml"
    resources:
        bigwig_convert=1,
    params:
        chrom_sizes=lambda wc: _pooled_chrom_sizes(wc.experiment),
        sort_tmpdir="/tmp",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        SORTED="$(dirname {output:q})/.tmp_{wildcards.experiment}.pooled.ppois.sorted.bdg"
        trap 'rm -f "$SORTED"' EXIT
        LC_ALL=C sort --temporary-directory={params.sort_tmpdir:q} -k1,1 -k2,2n {input:q} >"$SORTED"
        bedGraphToBigWig "$SORTED" {params.chrom_sizes:q} {output:q} 2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 6b. Pooled experiment QC summary (Stage 6b)
# ---------------------------------------------------------------------------


def _classify(experiment):
    """Return histone classification dict for an experiment's first treatment."""
    sid = TREATMENT_SAMPLES_BY_EXPERIMENT[experiment][0]
    s = SAMPLE_MAP[sid]
    return classify_histone_target(s["target"], s["peak_mode"])


def _pooled_qc_peaks_input(experiment):
    """Return deliverable pooled peaks input for QC summary."""
    suffix = _pooled_peak_suffix(experiment)
    if ANALYSIS_FILTERED:
        return (
            f"{OUTDIR}/experiments/{experiment}/04_peaks/pooled/"
            f"{experiment}_pooled_peaks_blacklist_filtered/"
            f"{experiment}_pooled_peaks.blacklist_filtered.{suffix}"
        )
    return f"{OUTDIR}/experiments/{experiment}/04_peaks/pooled/{experiment}_pooled_peaks"


rule pooled_experiment_qc_summary:
    input:
        pooled_bam=f"{OUTDIR}/experiments/{{experiment}}/02_align/{{experiment}}.pooled.final.bam",
        pooled_peaks=lambda wc: _pooled_qc_peaks_input(wc.experiment),
        pooled_fe=(
            f"{OUTDIR}/experiments/{{experiment}}/03_signal/{{experiment}}.pooled.FE.bdg"
            if QC_CONFIG.get("signal_tracks", True)
            else []
        ),
        pooled_ppois=(
            f"{OUTDIR}/experiments/{{experiment}}/03_signal/{{experiment}}.pooled.ppois.bdg"
            if QC_CONFIG.get("signal_tracks", True)
            else []
        ),
    output:
        f"{OUTDIR}/experiments/{{experiment}}/01_qc/{{experiment}}.pooled_qc_summary.tsv",
    conda:
        "../envs/python.yml"
    params:
        script=f"{workflow.basedir}/../scripts/pooled_qc_summary.py",
        experiment="{experiment}",
        target=lambda wc: SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[wc.experiment][0]]["target"],
        assay=lambda wc: SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[wc.experiment][0]]["assay"],
        peak_mode=lambda wc: SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[wc.experiment][0]]["peak_mode"],
        n_bioreps=lambda wc: len(_bioreps_for(wc.experiment, "treatment")),
        bio_rep_labels=lambda wc: ",".join(str(b) for b in sorted(_bioreps_for(wc.experiment, "treatment"))),
        inferred_histone_class=lambda wc: _classify(wc.experiment)["inferred_histone_class"],
        expected_peak_mode=lambda wc: _classify(wc.experiment)["expected_peak_mode"],
        peak_mode_status=lambda wc: _classify(wc.experiment)["peak_mode_status"],
        pooled_fe_path=(
            f"{OUTDIR}/experiments/{{experiment}}/03_signal/{{experiment}}.pooled.FE.bdg"
            if QC_CONFIG.get("signal_tracks", True)
            else "NA"
        ),
        pooled_ppois_path=(
            f"{OUTDIR}/experiments/{{experiment}}/03_signal/{{experiment}}.pooled.ppois.bdg"
            if QC_CONFIG.get("signal_tracks", True)
            else "NA"
        ),
        signal_enabled=QC_CONFIG.get("signal_tracks", True),
    shell:
        """
        set -e -o pipefail
        if [[ "{params.signal_enabled}" == "True" ]]; then
            if [[ -f "{input.pooled_fe}" && -f "{input.pooled_ppois}" ]]; then
                SIG_STATUS="enabled_present"
            else
                SIG_STATUS="enabled_missing"
            fi
        else
            SIG_STATUS="disabled"
        fi

        mkdir -p "$(dirname {output:q})"
        python3 {params.script:q} \
            --experiment {params.experiment:q} \
            --assay {params.assay:q} \
            --target {params.target:q} \
            --peak-mode {params.peak_mode:q} \
            --inferred-histone-class {params.inferred_histone_class:q} \
            --expected-peak-mode {params.expected_peak_mode:q} \
            --peak-mode-status {params.peak_mode_status:q} \
            --n-bioreps {params.n_bioreps} \
            --bio-rep-labels {params.bio_rep_labels:q} \
            --pooled-bam {input.pooled_bam:q} \
            --pooled-peaks-dir {input.pooled_peaks:q} \
            --pooled-fe-bdg {params.pooled_fe_path:q} \
            --pooled-ppois-bdg {params.pooled_ppois_path:q} \
            --signal-tracks-status "$SIG_STATUS" \
            --output {output:q}
        """


# ---------------------------------------------------------------------------
# 7. Per-sample QC summary
# ---------------------------------------------------------------------------


rule qc_summary:
    input:
        peak_counts=f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.peak_counts.tsv",
        frip=f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.frip.tsv",
        library_complexity=f"{RAW_OUTDIR}/{{sample}}/01_qc/{{sample}}.library_complexity.tsv",
        nrf_pbc=f"{RAW_OUTDIR}/{{sample}}/01_qc/{{sample}}.nrf_pbc.tsv",
        final_bam=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
        analysis_bam=lambda wc: analysis_bam(wc.sample),
        mito=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.mito_fraction.tsv",
        analysis_policy=f"{OUTDIR}/{{sample}}/05_qc/{{sample}}.analysis_policy.json",
        fragment_qc=lambda wc: (
            [f"{OUTDIR}/{wc.sample}/01_qc/{wc.sample}.cuttag_fragment_size.tsv"]
            if (
                QC_CONFIG.get("cuttag_fragment_size", True)
                and SAMPLE_MAP[wc.sample]["assay"] == "cuttag"
            )
            else []
        ),
        cross_correlation_qc=lambda wc: (
            [f"{OUTDIR}/{wc.sample}/05_qc/cross_correlation/{wc.sample}.cc.qc"]
            if QC_CONFIG.get("cross_correlation", False)
            else []
        ),
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.qc_summary.tsv",
    conda:
        "../envs/python.yml"
    params:
        script=f"{workflow.basedir}/../scripts/assemble_qc_summary.py",
        sample="{sample}",
        assay=lambda wc: SAMPLE_MAP[wc.sample]["assay"],
        target=lambda wc: SAMPLE_MAP[wc.sample]["target"],
        genome=lambda wc: SAMPLE_MAP[wc.sample]["genome"],
        layout=lambda wc: SAMPLE_MAP[wc.sample]["layout"],
        peak_mode=lambda wc: SAMPLE_MAP[wc.sample]["peak_mode"],
        has_blacklist=lambda wc: ("yes" if (ANALYSIS_FILTERED or _has_blacklist_qc(wc.sample)) else "no"),
        blacklist=lambda wc: (
            get_genome_resource(wc.sample, "blacklist")
            if (ANALYSIS_FILTERED or _has_blacklist_qc(wc.sample))
            else "NA"
        ),
        peaks_file=lambda wc: (_analysis_filtered_peak_file(wc.sample) if ANALYSIS_FILTERED else _peak_file(wc.sample)),
        bl_bam=lambda wc: (
            analysis_bam(wc.sample)
            if ANALYSIS_FILTERED
            else (_blacklist_bam_file(wc.sample) if _has_blacklist_qc(wc.sample) else "NA")
        ),
        bl_peaks=lambda wc: (
            _analysis_filtered_peak_file(wc.sample)
            if ANALYSIS_FILTERED
            else (_blacklist_peak_file(wc.sample) if _has_blacklist_qc(wc.sample) else "NA")
        ),
        ctrl_type=lambda wc: _control_type(wc.sample),
        fragment_qc=lambda wc: (
            f"{OUTDIR}/{wc.sample}/01_qc/{wc.sample}.cuttag_fragment_size.tsv"
            if (
                QC_CONFIG.get("cuttag_fragment_size", True)
                and SAMPLE_MAP[wc.sample]["assay"] == "cuttag"
            )
            else ""
        ),
        cross_correlation_qc=lambda wc: (
            f"{OUTDIR}/{wc.sample}/05_qc/cross_correlation/{wc.sample}.cc.qc"
            if QC_CONFIG.get("cross_correlation", False)
            else ""
        ),
        cross_correlation_enabled=str(QC_CONFIG.get("cross_correlation", False)).lower(),
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})"
        python3 {params.script:q} \\
        --sample {params.sample:q} \\
        --assay {params.assay:q} \\
        --target {params.target:q} \\
        --genome {params.genome:q} \\
        --layout {params.layout:q} \\
        --peak-mode {params.peak_mode:q} \\
        --use-control "{USE_CONTROL}" \\
        --control-type {params.ctrl_type:q} \\
        --final-bam {input.final_bam:q} \\
        --analysis-bam {input.analysis_bam:q} \\
        --peaks-file {params.peaks_file:q} \\
        --has-blacklist {params.has_blacklist:q} \\
        --blacklist {params.blacklist:q} \\
        --bl-bam {params.bl_bam:q} \\
        --bl-peaks {params.bl_peaks:q} \\
        --peak-counts {input.peak_counts:q} \\
        --frip {input.frip:q} \\
        --library-complexity {input.library_complexity:q} \\
        --nrf-pbc {input.nrf_pbc:q} \\
        --mito-tsv {input.mito:q} \\
        --analysis-policy {input.analysis_policy:q} \\
        --fragment-qc {params.fragment_qc:q} \\
        --cross-correlation-qc {params.cross_correlation_qc:q} \\
        --cross-correlation-enabled {params.cross_correlation_enabled} \\
        --output {output:q}
        """


def _control_type(sample_id):
    """Return control type string for a sample: none, control_bam, or control_sample."""
    s = SAMPLE_MAP[sample_id]
    if not USE_CONTROL:
        return "none"
    if s.get("control_bam"):
        return "control_bam"
    if s.get("control_sample"):
        return "control_sample"
    return "none"


# ---------------------------------------------------------------------------
# 9. Stage 12: phantompeakqualtools cross-correlation
# ---------------------------------------------------------------------------


rule cross_correlation:
    input:
        bam=lambda wc: analysis_bam(wc.sample),
    output:
        qc=f"{OUTDIR}/{{sample}}/05_qc/cross_correlation/{{sample}}.cc.qc",
        pdf=f"{OUTDIR}/{{sample}}/05_qc/cross_correlation/{{sample}}.cc.plot.pdf",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.cross_correlation.log",
    conda:
        "../envs/phantompeakqualtools.yml"
    threads: THREADS
    shell:
        """
        set -e -o pipefail
        command -v run_spp.R >/dev/null 2>&1 || {{
                    echo "ERROR: run_spp.R not found. Install phantompeakqualtools conda env." >&2
                    exit 1
                }}
        mkdir -p "$(dirname {output.qc:q})" "$(dirname {log:q})"
        rm -f {output.qc:q} {output.pdf:q}
        run_spp.R -c={input.bam:q} -savp={output.pdf:q} -out={output.qc:q} \
            -x=-500:15 -rf -p={threads} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 10. Stage 12: preseq library complexity
# ---------------------------------------------------------------------------


rule preseq_complexity:
    input:
        bam=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.{MAPQ_TAG}.bam",
    output:
        f"{RAW_OUTDIR}/{{sample}}/05_qc/preseq/{{sample}}.preseq.txt",
    log:
        f"{RAW_OUTDIR}/{{sample}}/logs/{{sample}}.preseq.log",
    conda:
        "../envs/preseq.yml"
    params:
        paired=lambda wc: "-P" if SAMPLE_MAP[wc.sample]["layout"] == "PE" else "",
    shell:
        """
        set -e -o pipefail
        command -v preseq >/dev/null 2>&1 || {{
                    echo "ERROR: preseq not found. Install preseq conda env." >&2
                    exit 1
                }}
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"
        preseq lc_extrap -B {params.paired} -o {output:q} {input.bam:q} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 11. Stage 12: Picard CollectMultipleMetrics
# ---------------------------------------------------------------------------

_PICARD_PE_SAMPLES = [sid for sid in TREATMENT_SAMPLE_IDS if SAMPLE_MAP[sid]["layout"] == "PE"]
_PICARD_SE_SAMPLES = [sid for sid in TREATMENT_SAMPLE_IDS if SAMPLE_MAP[sid]["layout"] == "SE"]


rule picard_collect_multiple_metrics:
    input:
        bam=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
        bai=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam.bai",
    output:
        alignment=f"{RAW_OUTDIR}/{{sample}}/05_qc/picard/{{sample}}.alignment_summary_metrics",
        insert_size=f"{RAW_OUTDIR}/{{sample}}/05_qc/picard/{{sample}}.insert_size_metrics",
        quality=f"{RAW_OUTDIR}/{{sample}}/05_qc/picard/{{sample}}.quality_distribution_metrics",
        ins_hist=f"{RAW_OUTDIR}/{{sample}}/05_qc/picard/{{sample}}.insert_size_histogram.pdf",
    log:
        f"{RAW_OUTDIR}/{{sample}}/logs/{{sample}}.picard_metrics.log",
    wildcard_constraints:
        sample=wildcard_choices(_PICARD_PE_SAMPLES),
    conda:
        "../envs/picard.yml"
    params:
        ref=lambda wc: get_genome_resource(wc.sample, "reference_fasta"),
    shell:
        """
        set -e -o pipefail
        command -v picard >/dev/null 2>&1 || {{
                    echo "ERROR: picard not found. Install picard conda env." >&2
                    exit 1
                }}
        mkdir -p "$(dirname {output.alignment:q})" "$(dirname {log:q})"
        picard CollectMultipleMetrics \
            I={input.bam:q} \
            R={params.ref:q} \
            O={RAW_OUTDIR:q}/{wildcards.sample}/05_qc/picard/{wildcards.sample} \
            PROGRAM=CollectAlignmentSummaryMetrics \
            PROGRAM=CollectInsertSizeMetrics \
            PROGRAM=QualityScoreDistribution \
            VALIDATION_STRINGENCY=LENIENT \
            2>&1 | tee {log:q}
        """


rule picard_collect_multiple_metrics_se:
    input:
        bam=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
        bai=f"{RAW_OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam.bai",
    output:
        alignment=f"{RAW_OUTDIR}/{{sample}}/05_qc/picard/{{sample}}.alignment_summary_metrics",
        quality=f"{RAW_OUTDIR}/{{sample}}/05_qc/picard/{{sample}}.quality_distribution_metrics",
    log:
        f"{RAW_OUTDIR}/{{sample}}/logs/{{sample}}.picard_metrics.log",
    wildcard_constraints:
        sample=wildcard_choices(_PICARD_SE_SAMPLES),
    conda:
        "../envs/picard.yml"
    params:
        ref=lambda wc: get_genome_resource(wc.sample, "reference_fasta"),
    shell:
        """
        set -e -o pipefail
        command -v picard >/dev/null 2>&1 || {{
                    echo "ERROR: picard not found. Install picard conda env." >&2
                    exit 1
                }}
        mkdir -p "$(dirname {output.alignment:q})" "$(dirname {log:q})"
        picard CollectMultipleMetrics \
            I={input.bam:q} \
            R={params.ref:q} \
            O={RAW_OUTDIR:q}/{wildcards.sample}/05_qc/picard/{wildcards.sample} \
            PROGRAM=null \
            PROGRAM=CollectAlignmentSummaryMetrics \
            PROGRAM=QualityScoreDistribution \
            VALIDATION_STRINGENCY=LENIENT \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 12. Stage 18: TSS enrichment profile
# ---------------------------------------------------------------------------


rule tss_bed_from_gtf:
    input:
        gtf=lambda wc: GENOME_RESOURCES.get(wc.genome, {}).get("gtf", ""),
        identity=lambda wc: _analysis_identity_inputs(wc),
    output:
        f"{OUTDIR}/reference/{{genome}}.tss.bed",
    conda:
        "../envs/python.yml"
    params:
        script=f"{workflow.basedir}/../scripts/gtf_to_tss_bed.py",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})"
        python3 {params.script:q} \
            --gtf {input.gtf:q} \
            --output {output:q}
        """


rule tss_enrichment_profile:
    input:
        bw=f"{OUTDIR}/{{sample}}/03_bigwig/{{sample}}.CPM.bw",
        tss=lambda wc: f"{OUTDIR}/reference/{SAMPLE_MAP[wc.sample]['genome']}.tss.bed",
    output:
        matrix=f"{OUTDIR}/{{sample}}/05_qc/tss/{{sample}}.tss_matrix.gz",
        profile_tsv=f"{OUTDIR}/{{sample}}/05_qc/tss/{{sample}}.tss_profile.tsv",
        profile_pdf=f"{OUTDIR}/{{sample}}/05_qc/tss/{{sample}}.tss_profile.pdf",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.tss_enrichment.log",
    conda:
        "../envs/deeptools.yml"
    threads: THREADS
    shell:
        """
        set -e -o pipefail
        command -v computeMatrix >/dev/null 2>&1 || {{
                    echo "ERROR: computeMatrix not found. Install deeptools conda env." >&2
                    exit 1
                }}
        command -v plotProfile >/dev/null 2>&1 || {{
                    echo "ERROR: plotProfile not found. Install deeptools conda env." >&2
                    exit 1
                }}
        mkdir -p "$(dirname {output.matrix:q})" "$(dirname {log:q})"
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 3000 -a 3000 \
            -R {input.tss:q} \
            -S {input.bw:q} \
            --skipZeros \
            -o {output.matrix:q} \
            -p {threads} \
            2>&1 | tee {log:q}
        plotProfile \
            -m {output.matrix:q} \
            -out {output.profile_pdf:q} \
            --outFileNameData {output.profile_tsv:q} \
            --plotTitle {wildcards.sample:q} \
            2>&1 | tee -a {log:q}
        """


rule tss_fragment_profile:
    input:
        bw=lambda wc: fragment_cpm_bw(wc.sample),
        tss=lambda wc: f"{OUTDIR}/reference/{SAMPLE_MAP[wc.sample]['genome']}.tss.bed",
    output:
        matrix=f"{OUTDIR}/{{sample}}/05_qc/tss/{{sample}}.fragment.tss_matrix.gz",
        profile_tsv=f"{OUTDIR}/{{sample}}/05_qc/tss/{{sample}}.fragment.tss_profile.tsv",
        profile_pdf=f"{OUTDIR}/{{sample}}/05_qc/tss/{{sample}}.fragment.tss_profile.pdf",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.tss_fragment_profile.log",
    wildcard_constraints:
        sample=wildcard_choices(TSS_FRAGMENT_SAMPLE_IDS),
    conda:
        "../envs/deeptools.yml"
    threads: THREADS
    shell:
        """
        set -e -o pipefail
        command -v computeMatrix >/dev/null 2>&1 || {{
                    echo "ERROR: computeMatrix not found. Install deeptools conda env." >&2
                    exit 1
                }}
        command -v plotProfile >/dev/null 2>&1 || {{
                    echo "ERROR: plotProfile not found. Install deeptools conda env." >&2
                    exit 1
                }}
        mkdir -p "$(dirname {output.matrix:q})" "$(dirname {log:q})"
        computeMatrix reference-point \
            --referencePoint TSS \
            -b 3000 -a 3000 \
            -R {input.tss:q} \
            -S {input.bw:q} \
            --skipZeros \
            -o {output.matrix:q} \
            -p {threads} \
            2>&1 | tee {log:q}
        plotProfile \
            -m {output.matrix:q} \
            -out {output.profile_pdf:q} \
            --outFileNameData {output.profile_tsv:q} \
            --plotTitle "{wildcards.sample} fragment occupancy" \
            2>&1 | tee -a {log:q}
        """


# ---------------------------------------------------------------------------
# 8. Project-level QC summary
# ---------------------------------------------------------------------------


rule project_qc_summary:
    input:
        [f"{OUTDIR}/{sid}/01_qc/{sid}.qc_summary.tsv" for sid in PEAK_SAMPLE_IDS],
    output:
        f"{OUTDIR}/multiqc/project_qc_summary.tsv",
    conda:
        "../envs/python.yml"
    params:
        script=f"{workflow.basedir}/../scripts/aggregate_qc_summary.py",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})"
        python3 {params.script:q} \\
        --output {output:q} \\
        {input:q}
        """
