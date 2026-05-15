# qc.smk — Stage 3 Single-Sample QC rules
# ========================================
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
#   library_complexity      → parse Picard MarkDuplicates metrics
#   signal_track_fe         → MACS3 bdgcmp FE bedGraph
#   signal_track_ppois      → MACS3 bdgcmp ppois bedGraph
#   qc_summary              → per-sample TSV assembly
#   stage3_qc_summary       → project-level aggregated summary

# Helper: peak file suffix for a sample
def _peak_suffix(sample_id):
    return "broadPeak" if SAMPLE_MAP[sample_id]["peak_mode"] == "broad" else "narrowPeak"


def _peak_file(sample_id):
    """Return the MACS3 peak file path for a sample."""
    suffix = _peak_suffix(sample_id)
    return (
        f"{OUTDIR}/{sample_id}/04_peaks/{sample_id}/"
        f"{sample_id}_peaks.{suffix}"
    )


def _has_blacklist_qc(sample_id):
    """Return True when blacklist filtering is enabled and resource-backed."""
    return (
        sample_id in BLACKLIST_SAMPLES
        and QC_CONFIG.get("blacklist_filter", True)
    )


def _blacklist_peak_file(sample_id):
    """Return the blacklist-filtered peak file path for a sample."""
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


# ---------------------------------------------------------------------------
# 1. Blacklist filter — BAM
# ---------------------------------------------------------------------------

rule blacklist_filter_bam:
    output:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.blacklist_filtered.bam",
        bai = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.blacklist_filtered.bam.bai",
    input:
        bam = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
        bai = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam.bai",
    params:
        blacklist = lambda wc: get_genome_resource(wc.sample, "blacklist"),
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        bedtools intersect -v -abam {input.bam:q} -b {params.blacklist:q} \
            > {output.bam:q}
        samtools index -@ {threads} {output.bam:q}
        """


# ---------------------------------------------------------------------------
# 2. Blacklist filter — peaks
# ---------------------------------------------------------------------------

rule blacklist_filter_peaks:
    output:
        f"{OUTDIR}/{{sample}}/04_peaks/{{sample}}_blacklist_filtered/{{sample}}_peaks.blacklist_filtered.{{suffix}}",
    input:
        # MACS3 peaks directory (Snakemake-tracked) — peak file assembled in shell
        peaks_dir = f"{OUTDIR}/{{sample}}/04_peaks/{{sample}}",
    params:
        blacklist = lambda wc: get_genome_resource(wc.sample, "blacklist"),
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        PEAKS_FILE="{input.peaks_dir:q}/{wildcards.sample}_peaks.{wildcards.suffix}"
        mkdir -p $(dirname {output:q})
        bedtools intersect -v -a "$PEAKS_FILE" -b {params.blacklist:q} \
            > {output:q}
        """


# ---------------------------------------------------------------------------
# 3. Peak counts
# ---------------------------------------------------------------------------

rule peak_counts:
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.peak_counts.tsv",
    input:
        lambda wc: _peak_counts_inputs(wc),
    params:
        peak_mode = lambda wc: SAMPLE_MAP[wc.sample]["peak_mode"],
        sample    = "{sample}",
        bl_peak   = lambda wc: (
            _blacklist_peak_file(wc.sample)
            if _has_blacklist_qc(wc.sample) else ""
        ),
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        PEAKS_DIR="{input[0]:q}"
        SUFFIX="{params.peak_mode}"
        [[ "$SUFFIX" == "broad" ]] && SUFFIX="broadPeak" || SUFFIX="narrowPeak"
        PEAKS_FILE="$PEAKS_DIR/{params.sample}_peaks.$SUFFIX"
        PEAK_COUNT=$(wc -l < "$PEAKS_FILE")

        BL_PEAK={params.bl_peak:q}
        if [[ -n "$BL_PEAK" && -f "$BL_PEAK" ]]; then
            BL_COUNT=$(wc -l < "$BL_PEAK")
        else
            BL_COUNT="NA"
        fi

        printf "sample\\tpeak_mode\\tpeaks\\tblacklist_filtered_peaks\\n" \
            > {output:q}
        printf "%s\\t%s\\t%s\\t%s\\n" \
            {params.sample:q} {params.peak_mode:q} \
            "$PEAK_COUNT" "$BL_COUNT" \
            >> {output:q}
        """


def _peak_counts_inputs(wildcards):
    """Return [peaks_dir] or [peaks_dir, blacklist_filtered_peak_file]."""
    sid = wildcards.sample
    peaks_dir = f"{OUTDIR}/{sid}/04_peaks/{sid}"
    inputs = [peaks_dir]
    if _has_blacklist_qc(sid):
        inputs.append(_blacklist_peak_file(sid))
    return inputs


# ---------------------------------------------------------------------------
# 4. FRiP (Fraction of Reads in Peaks)
# ---------------------------------------------------------------------------

rule frip:
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.frip.tsv",
    input:
        lambda wc: _frip_inputs(wc),
    params:
        peak_mode = lambda wc: SAMPLE_MAP[wc.sample]["peak_mode"],
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.frip.log",
    conda:
        "../envs/chipseq.yml",
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

        python3 scripts/calc_frip.py \
            --sample {wildcards.sample:q} \
            --bam "$BAM" \
            --peaks "$PEAKS" \
            --output {output:q} \
            2>&1 | tee {log:q}
        """


def _frip_inputs(wildcards):
    """Return [bam, peaks_dir_or_file] for FRiP calculation.

    Prefers blacklist-filtered BAM and peaks only when BOTH are
    scheduled/available for this sample. Otherwise falls back to
    unfiltered BAM and peaks directory.
    """
    sid = wildcards.sample
    if _has_blacklist_qc(sid):
        return [_blacklist_bam_file(sid), _blacklist_peak_file(sid)]
    else:
        bam = f"{OUTDIR}/{sid}/02_align/{sid}.final.bam"
        peaks_dir = f"{OUTDIR}/{sid}/04_peaks/{sid}"
        return [bam, peaks_dir]


# ---------------------------------------------------------------------------
# 5. Library complexity (Stage 3b-1)
# ---------------------------------------------------------------------------

rule library_complexity:
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.library_complexity.tsv",
    input:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.dup_metrics.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.library_complexity.log",
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        python3 scripts/parse_dup_metrics.py \
            --sample {wildcards.sample:q} \
            --metrics {input:q} \
            --output {output:q} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 5a. NRF/PBC library complexity (Stage 3c-1)
# ---------------------------------------------------------------------------

rule nrf_pbc:
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.nrf_pbc.tsv",
    input:
        f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.nrf_pbc.log",
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        python3 scripts/calc_nrf_pbc.py \
            --sample {wildcards.sample:q} \
            --bam {input:q} \
            --output {output:q} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 6. MACS3 signal tracks (Stage 3b-2)
# ---------------------------------------------------------------------------

rule signal_track_fe:
    output:
        f"{OUTDIR}/{{sample}}/03_signal/{{sample}}.FE.bdg",
    input:
        peaks_dir = f"{OUTDIR}/{{sample}}/04_peaks/{{sample}}",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.bdgcmp.FE.log",
    conda:
        "../envs/chipseq.yml",
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
    output:
        f"{OUTDIR}/{{sample}}/03_signal/{{sample}}.ppois.bdg",
    input:
        peaks_dir = f"{OUTDIR}/{{sample}}/04_peaks/{{sample}}",
    log:
        f"{OUTDIR}/{{sample}}/logs/{{sample}}.bdgcmp.ppois.log",
    conda:
        "../envs/chipseq.yml",
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
# 7. Per-sample QC summary
# ---------------------------------------------------------------------------

rule qc_summary:
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.qc_summary.tsv",
    input:
        peak_counts         = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.peak_counts.tsv",
        frip                = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.frip.tsv",
        library_complexity  = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.library_complexity.tsv",
        nrf_pbc             = f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.nrf_pbc.tsv",
        final_bam           = f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
    params:
        sample     = "{sample}",
        assay      = lambda wc: SAMPLE_MAP[wc.sample]["assay"],
        target     = lambda wc: SAMPLE_MAP[wc.sample]["target"],
        genome     = lambda wc: SAMPLE_MAP[wc.sample]["genome"],
        layout     = lambda wc: SAMPLE_MAP[wc.sample]["layout"],
        peak_mode  = lambda wc: SAMPLE_MAP[wc.sample]["peak_mode"],
        has_blacklist = lambda wc: (
            "yes" if wc.sample in BLACKLIST_SAMPLES
            and QC_CONFIG.get("blacklist_filter", True) else "no"
        ),
        blacklist = lambda wc: (
            get_genome_resource(wc.sample, "blacklist")
            if _has_blacklist_qc(wc.sample) else "NA"
        ),
        peaks_file = lambda wc: _peak_file(wc.sample),
        bl_bam = lambda wc: (
            _blacklist_bam_file(wc.sample)
            if _has_blacklist_qc(wc.sample) else "NA"
        ),
        bl_peaks = lambda wc: (
            _blacklist_peak_file(wc.sample)
            if _has_blacklist_qc(wc.sample) else "NA"
        ),
        ctrl_type = lambda wc: _control_type(wc.sample),
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        # Read peak counts (skip header)
        PEAK_COUNT=$(tail -n +2 {input.peak_counts:q} | cut -f3)
        BL_PEAK_COUNT=$(tail -n +2 {input.peak_counts:q} | cut -f4)

        # Read FRiP values (skip header)
        TOTAL_READS=$(tail -n +2 {input.frip:q} | cut -f2)
        READS_IN_PEAKS=$(tail -n +2 {input.frip:q} | cut -f3)
        FRIP=$(tail -n +2 {input.frip:q} | cut -f4)

        # Read library complexity values (skip header and sample column)
        LC_DATA=$(tail -n +2 {input.library_complexity:q})
        METRICS_SOURCE=$(echo "$LC_DATA" | cut -f2)
        LC_UP_EXAM=$(echo "$LC_DATA" | cut -f3)
        LC_RP_EXAM=$(echo "$LC_DATA" | cut -f4)
        LC_SEC_SUP=$(echo "$LC_DATA" | cut -f5)
        LC_UNMAPPED=$(echo "$LC_DATA" | cut -f6)
        LC_UP_DUP=$(echo "$LC_DATA" | cut -f7)
        LC_RP_DUP=$(echo "$LC_DATA" | cut -f8)
        LC_RP_OPT=$(echo "$LC_DATA" | cut -f9)
        LC_PCT_DUP=$(echo "$LC_DATA" | cut -f10)
        LC_EST_SIZE=$(echo "$LC_DATA" | cut -f11)
        LC_TOT_EXAM=$(echo "$LC_DATA" | cut -f12)
        LC_EST_DUP=$(echo "$LC_DATA" | cut -f13)

        # Read NRF/PBC values (skip header)
        NRF_PBC_DATA=$(tail -n +2 {input.nrf_pbc:q})
        NRF_TOTAL=$(echo "$NRF_PBC_DATA" | cut -f2)
        NRF_DISTINCT=$(echo "$NRF_PBC_DATA" | cut -f3)
        NRF_ONE=$(echo "$NRF_PBC_DATA" | cut -f4)
        NRF_TWO=$(echo "$NRF_PBC_DATA" | cut -f5)
        NRF_NRF=$(echo "$NRF_PBC_DATA" | cut -f6)
        NRF_PBC1=$(echo "$NRF_PBC_DATA" | cut -f7)
        NRF_PBC2=$(echo "$NRF_PBC_DATA" | cut -f8)

        if [[ "{params.has_blacklist}" == "no" ]]; then
            if [[ "$BL_PEAK_COUNT" != "NA" ]]; then
                BL_PEAK_COUNT="NA"
            fi
        fi

        # Header
        printf "sample\\tassay\\ttarget\\tgenome\\tlayout\\tpeak_mode\\t" \
            > {output:q}
        printf "use_control\\tcontrol_type\\tfinal_bam\\tpeaks\\t" \
            >> {output:q}
        printf "blacklist\\tblacklist_filtered_bam\\t" \
            >> {output:q}
        printf "blacklist_filtered_peaks\\ttotal_reads\\t" \
            >> {output:q}
        printf "reads_in_peaks\\tfrip\\tpeak_count\\t" \
            >> {output:q}
        printf "blacklist_filtered_peak_count\\t" \
            >> {output:q}
        printf "metrics_source\\tunpaired_reads_examined\\t" \
            >> {output:q}
        printf "read_pairs_examined\\t" \
            >> {output:q}
        printf "secondary_or_supplementary_reads\\t" \
            >> {output:q}
        printf "unmapped_reads\\tunpaired_read_duplicates\\t" \
            >> {output:q}
        printf "read_pair_duplicates\\t" \
            >> {output:q}
        printf "read_pair_optical_duplicates\\t" \
            >> {output:q}
        printf "percent_duplication\\testimated_library_size\\t" \
            >> {output:q}
        printf "total_reads_examined\\t" \
            >> {output:q}
        printf "duplicate_reads_estimate\\t" \
            >> {output:q}
        printf "total_fragments\\tdistinct_fragments\\t" \
            >> {output:q}
        printf "one_read_fragments\\ttwo_read_fragments\\t" \
            >> {output:q}
        printf "nrf\\tpbc1\\tpbc2\\n" \
            >> {output:q}

        # Data
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t" \
            {params.sample:q} \
            {params.assay:q} \
            {params.target:q} \
            {params.genome:q} \
            {params.layout:q} \
            {params.peak_mode:q} \
            >> {output:q}
        printf "%s\\t%s\\t%s\\t%s\\t" \
            "{USE_CONTROL}" \
            {params.ctrl_type:q} \
            {input.final_bam:q} \
            {params.peaks_file:q} \
            >> {output:q}
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t" \
            {params.blacklist:q} \
            {params.bl_bam:q} \
            {params.bl_peaks:q} \
            "$TOTAL_READS" \
            "$READS_IN_PEAKS" \
            "$FRIP" \
            "$PEAK_COUNT" \
            "$BL_PEAK_COUNT" \
            >> {output:q}
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t" \
            "$METRICS_SOURCE" \
            "$LC_UP_EXAM" \
            "$LC_RP_EXAM" \
            "$LC_SEC_SUP" \
            "$LC_UNMAPPED" \
            "$LC_UP_DUP" \
            "$LC_RP_DUP" \
            "$LC_RP_OPT" \
            "$LC_PCT_DUP" \
            "$LC_EST_SIZE" \
            "$LC_TOT_EXAM" \
            "$LC_EST_DUP" \
            >> {output:q}
        printf "%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n" \
            "$NRF_TOTAL" \
            "$NRF_DISTINCT" \
            "$NRF_ONE" \
            "$NRF_TWO" \
            "$NRF_NRF" \
            "$NRF_PBC1" \
            "$NRF_PBC2" \
            >> {output:q}
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
# 8. Project-level QC summary
# ---------------------------------------------------------------------------

rule stage3_qc_summary:
    output:
        f"{OUTDIR}/multiqc/stage3_qc_summary.tsv",
    input:
        [f"{OUTDIR}/{sid}/01_qc/{sid}.qc_summary.tsv"
         for sid in TREATMENT_SAMPLE_IDS],
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -- {input:q}
        if [[ "$#" -eq 0 ]]; then
            printf "sample\\tassay\\ttarget\\tgenome\\tlayout\\tpeak_mode\\t" \
                > {output:q}
            printf "use_control\\tcontrol_type\\tfinal_bam\\tpeaks\\t" \
                >> {output:q}
            printf "blacklist\\tblacklist_filtered_bam\\t" \
                >> {output:q}
            printf "blacklist_filtered_peaks\\ttotal_reads\\t" \
                >> {output:q}
            printf "reads_in_peaks\\tfrip\\tpeak_count\\t" \
                >> {output:q}
            printf "blacklist_filtered_peak_count\\t" \
                >> {output:q}
            printf "metrics_source\\tunpaired_reads_examined\\t" \
                >> {output:q}
            printf "read_pairs_examined\\t" \
                >> {output:q}
            printf "secondary_or_supplementary_reads\\t" \
                >> {output:q}
            printf "unmapped_reads\\tunpaired_read_duplicates\\t" \
                >> {output:q}
            printf "read_pair_duplicates\\t" \
                >> {output:q}
            printf "read_pair_optical_duplicates\\t" \
                >> {output:q}
            printf "percent_duplication\\testimated_library_size\\t" \
                >> {output:q}
            printf "total_reads_examined\\t" \
                >> {output:q}
            printf "duplicate_reads_estimate\\t" \
                >> {output:q}
            printf "total_fragments\\tdistinct_fragments\\t" \
                >> {output:q}
            printf "one_read_fragments\\ttwo_read_fragments\\t" \
                >> {output:q}
            printf "nrf\\tpbc1\\tpbc2\\n" \
                >> {output:q}
            exit 0
        fi

        # Emit header from first file, then concatenate data rows
        head -n 1 "$1" > {output:q}
        for f in "$@"; do
            tail -n +2 "$f" >> {output:q}
        done
        """
