# idr.smk — Stage 5a + 5b TF ChIP-seq IDR
# =========================================
# Stage 5a: IDR-ready MACS3 per biorep, true-replicate IDR.
# Stage 5b: pseudorep BAM splitting, pseudorep MACS3, self-IDR,
# pooled-IDR, reproducibility summary.
#
# Rules:
#   macs3_idr_biorep     — IDR-ready MACS3 call on a single biorep BAM (5a)
#   idr_true_replicates  — IDR between the two biorep IDR peak sets (5a)
#   split_pseudoreps     — deterministic BAM pseudorep split (5b)
#   macs3_idr_pseudorep  — IDR-ready MACS3 on pseudorep BAM (5b)
#   idr_self_pseudoreps  — self-IDR per biorep (5b)
#   idr_pooled_pseudoreps— pooled pseudorep IDR (5b)
#   stage5b_summary      — reproducibility QC + final peaks (5b)


# ---------------------------------------------------------------------------
# Helper: inputs for macs3_idr_biorep
# ---------------------------------------------------------------------------

def _idr_biorep_peaks_inputs(wildcards):
    """Return inputs for MACS3 IDR peak call on a single biorep BAM.

    Input order: [biorep_bam, biorep_bam.bai, ...optional_pooled_control_bam]
    The pooled control BAM is included as an explicit dependency so Snakemake
    schedules pool_control_bam before this rule.
    """
    exp = wildcards.experiment
    br = int(wildcards.bio_rep)
    inputs = [
        idr_biorep_bam(exp, br),
        idr_biorep_bai(exp, br),
    ]
    if exp in POOLED_CONTROL_EXPERIMENTS:
        inputs.append(idr_pooled_control_bam(exp))
    return inputs


# ---------------------------------------------------------------------------
# Helper: MACS3 args for IDR-ready peak calls
# ---------------------------------------------------------------------------

def _idr_macs3_args(wildcards):
    """Return MACS3 args for IDR-ready peak calls on a biorep BAM.

    Uses the same layout/genome as per-sample MACS3, but replaces
    -q (q-value) with -p (p-value from idr_macs3 config).
    Never emits both -q and -p.
    """
    experiment = wildcards.experiment
    treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    if not treatment_ids:
        return ""
    s = SAMPLE_MAP[treatment_ids[0]]
    fmt = "BAMPE" if s["layout"] == "PE" else "BAM"
    genome = _normalize_genome(s["genome"])
    pvalue = _tool_param("idr_macs3", "pvalue", 0.1)
    extra = _tool_param("idr_macs3", "extra_args", "")
    return f"-f {fmt} -g {genome} -p {pvalue} {extra}".strip()


# ---------------------------------------------------------------------------
# 1. MACS3 IDR-ready peak call on a single biorep BAM
# ---------------------------------------------------------------------------

rule macs3_idr_biorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/04_peaks/idr/{{experiment}}_biorep{{bio_rep}}_idr_peaks.narrowPeak",
    input:
        lambda wc: _idr_biorep_peaks_inputs(wc),
    params:
        macs3_args = lambda wc: _idr_macs3_args(wc),
        sample     = lambda wc: f"{wc.experiment}_biorep{wc.bio_rep}_idr",
    wildcard_constraints:
        bio_rep = r"\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}_biorep{{bio_rep}}.idr.macs3.log",
    threads: THREADS,
    conda:
        "../envs/macs3.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output})" "$(dirname {log})"

        set -- {input:q}
        TREATMENT="$1"
        # Input order: treatment.bam ($1)[, treatment.bam.bai ($2)][, pooled.control.bam ($3)]
        if [[ $# -ge 3 ]]; then
            macs3 callpeak \
                -t "$TREATMENT" \
                -c "$3" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} \
                2>&1 | tee {log:q}
        else
            macs3 callpeak \
                -t "$TREATMENT" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} \
                2>&1 | tee {log:q}
        fi

        EXPECTED="{output}"
        if [[ ! -f "$EXPECTED" ]]; then
            echo "ERROR: Expected peak file not found: $EXPECTED" >&2
            exit 1
        fi
        """


# ---------------------------------------------------------------------------
# Helper: resolve IDR peak input path for a given bio_rep index
# ---------------------------------------------------------------------------

def _idr_peak_input(experiment, index):
    """Return the IDR peak file for a bio_rep by 0-based index.

    index 0 -> first bio_rep (smallest label), index 1 -> second.
    The actual bio_rep numbers are derived from the sample sheet.
    """
    bioreps = _bioreps_for(experiment, "treatment")
    br = bioreps[index]
    return (
        f"{OUTDIR}/experiments/{experiment}/04_peaks/idr/"
        f"{experiment}_biorep{br}_idr_peaks.narrowPeak"
    )


# ---------------------------------------------------------------------------
# 2. True replicate IDR
# ---------------------------------------------------------------------------

rule idr_true_replicates:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/true_replicates/idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/true_replicates/idr.thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: _idr_peak_input(wc.experiment, 0),
        peaks2 = lambda wc: _idr_peak_input(wc.experiment, 1),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_idr/true_replicates/idr"
        ),
        raw_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/{wc.experiment}.idr.raw.log"
        ),
        thr_log_file  = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/{wc.experiment}.idr.thresholded.log"
        ),
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        command -v idr >/dev/null 2>&1 || {{
            echo "ERROR: idr is required for Stage 5 IDR analysis but was not found in PATH." >&2
            exit 1
        }}

        mkdir -p "$(dirname {output.idr_out})" \
                 "$(dirname {params.raw_log_file})" \
                 "$(dirname {params.thr_log_file})"

        # Run 1: raw IDR output (no --idr-threshold).
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --output-file {output.idr_out:q} \
            --log-output-file {params.raw_log_file:q}

        # Run 2: thresholded narrowPeak output.
        # --idr-threshold N keeps peaks whose accumulated global IDR
        # score falls at or below N (smaller = more stringent).
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            --log-output-file {params.thr_log_file:q}
        """


# ============================================================================
# Stage 5b rules
# ============================================================================

# ---------------------------------------------------------------------------
# Helper: input BAM for pseudorep splitting
# ---------------------------------------------------------------------------

def _split_input(wildcards):
    """Return the input BAM path for a pseudorep split.

    source == "pooled" -> pooled treatment BAM
    source starts with "biorep" -> parse bio_rep label, return that biorep BAM
    """
    exp = wildcards.experiment
    src = wildcards.source
    if src == "pooled":
        return idr_pooled_treatment_bam(exp)
    # source format: "biorep<label>"
    br_label = src[len("biorep"):]
    return idr_biorep_bam(exp, br_label)


# ---------------------------------------------------------------------------
# Helper: inputs for pseudorep MACS3 peak call
# ---------------------------------------------------------------------------

def _idr_pseudorep_inputs(wildcards):
    """Return inputs for MACS3 IDR peak call on a pseudorep BAM.

    Input order: [pseudorep_bam, pseudorep_bam.bai, ...optional_pooled_control]
    """
    exp = wildcards.experiment
    src = wildcards.source
    pr = wildcards.pr
    inputs = [
        idr_pseudorep_bam(exp, src, pr),
        idr_pseudorep_bai(exp, src, pr),
    ]
    if exp in POOLED_CONTROL_EXPERIMENTS:
        inputs.append(idr_pooled_control_bam(exp))
    return inputs


# ---------------------------------------------------------------------------
# Helper: self-IDR thresholded path for a bio_rep index
# ---------------------------------------------------------------------------

def _self_thresh_path(experiment, index):
    """Return the thresholded self-IDR narrowPeak for a bio_rep at index."""
    bioreps = _bioreps_for(experiment, "treatment")
    br = bioreps[index]
    return (
        f"{OUTDIR}/experiments/{experiment}/06_idr/"
        f"self_pseudoreplicates/biorep{br}.idr.thresholded.narrowPeak"
    )


# ---------------------------------------------------------------------------
# 3. Pseudoreplicate BAM splitting — deterministic hash-based
# ---------------------------------------------------------------------------

rule split_pseudoreps:
    output:
        pr1 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/{{experiment}}_{{source}}.pr1.bam",
        p1i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/{{experiment}}_{{source}}.pr1.bam.bai",
        pr2 = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/{{experiment}}_{{source}}.pr2.bam",
        p2i = f"{OUTDIR}/experiments/{{experiment}}/05_pseudorep/{{experiment}}_{{source}}.pr2.bam.bai",
    input:
        lambda wc: _split_input(wc),
    params:
        seed = IDR_SEED,
    wildcard_constraints:
        source = r"pooled|biorep\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}_{{source}}.split.log",
    threads: THREADS,
    conda:
        "../envs/samtools.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.pr1})" "$(dirname {log})"
        python3 scripts/split_pseudoreps.py \
            --input {input:q} \
            --out1 {output.pr1:q} \
            --out2 {output.pr2:q} \
            --seed {params.seed} \
            --threads {threads} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 4. IDR-ready MACS3 on pseudorep BAM
# ---------------------------------------------------------------------------

rule macs3_idr_pseudorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/04_peaks/idr/{{experiment}}_{{source}}_pr{{pr}}_idr_peaks.narrowPeak",
    input:
        lambda wc: _idr_pseudorep_inputs(wc),
    params:
        macs3_args = lambda wc: _idr_macs3_args(wc),
        sample     = lambda wc: f"{wc.experiment}_{wc.source}_pr{wc.pr}_idr",
    wildcard_constraints:
        source = r"pooled|biorep\d+",
        pr     = r"[12]",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}_{{source}}_pr{{pr}}.idr.macs3.log",
    threads: THREADS,
    conda:
        "../envs/macs3.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output})" "$(dirname {log})"

        set -- {input:q}
        TREATMENT="$1"
        if [[ $# -ge 3 ]]; then
            macs3 callpeak -t "$TREATMENT" -c "$3" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} 2>&1 | tee {log:q}
        else
            macs3 callpeak -t "$TREATMENT" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} 2>&1 | tee {log:q}
        fi

        [[ -f "{output}" ]] || {{ echo "ERROR: Expected {output}" >&2; exit 1; }}
        """


# ---------------------------------------------------------------------------
# 5. Self-pseudoreplicate IDR (per biological replicate)
# ---------------------------------------------------------------------------

rule idr_self_pseudoreps:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/self_pseudoreplicates/biorep{{bio_rep}}.idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/self_pseudoreplicates/biorep{{bio_rep}}.idr.thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/04_peaks/idr/"
            f"{wc.experiment}_biorep{wc.bio_rep}_pr1_idr_peaks.narrowPeak"
        ),
        peaks2 = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/04_peaks/idr/"
            f"{wc.experiment}_biorep{wc.bio_rep}_pr2_idr_peaks.narrowPeak"
        ),
    params:
        threshold    = IDR_THRESHOLD,
        rank         = IDR_RANK,
        raw_log_file = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.idr.self.biorep{wc.bio_rep}.raw.log"
        ),
        thr_log_file = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.idr.self.biorep{wc.bio_rep}.thr.log"
        ),
    wildcard_constraints:
        bio_rep = r"\d+",
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        command -v idr >/dev/null 2>&1 || {{
            echo "ERROR: idr is required for Stage 5 IDR analysis but was not found in PATH." >&2
            exit 1
        }}

        mkdir -p "$(dirname {output.idr_out})" \
                 "$(dirname {params.raw_log_file})" \
                 "$(dirname {params.thr_log_file})"

        # Run 1: raw IDR output
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --output-file {output.idr_out:q} \
            --log-output-file {params.raw_log_file:q}

        # Run 2: thresholded narrowPeak output
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            --log-output-file {params.thr_log_file:q}
        """


# ---------------------------------------------------------------------------
# 6. Pooled-pseudoreplicate IDR
# ---------------------------------------------------------------------------

rule idr_pooled_pseudoreps:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/pooled_pseudoreplicates/idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/pooled_pseudoreplicates/idr.thresholded.narrowPeak",
    input:
        peaks1 = f"{OUTDIR}/experiments/{{experiment}}/04_peaks/idr/{{experiment}}_pooled_pr1_idr_peaks.narrowPeak",
        peaks2 = f"{OUTDIR}/experiments/{{experiment}}/04_peaks/idr/{{experiment}}_pooled_pr2_idr_peaks.narrowPeak",
    params:
        threshold    = IDR_THRESHOLD,
        rank         = IDR_RANK,
        raw_log_file = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.idr.pooled.raw.log"
        ),
        thr_log_file = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/logs/"
            f"{wc.experiment}.idr.pooled.thr.log"
        ),
    conda:
        "../envs/idr.yml",
    shell:
        """
        set -e -o pipefail
        command -v idr >/dev/null 2>&1 || {{
            echo "ERROR: idr is required for Stage 5 IDR analysis but was not found in PATH." >&2
            exit 1
        }}

        mkdir -p "$(dirname {output.idr_out})" \
                 "$(dirname {params.raw_log_file})" \
                 "$(dirname {params.thr_log_file})"

        # Run 1: raw IDR output
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --output-file {output.idr_out:q} \
            --log-output-file {params.raw_log_file:q}

        # Run 2: thresholded narrowPeak output
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            --log-output-file {params.thr_log_file:q}
        """


# ---------------------------------------------------------------------------
# 7. Reproducibility QC summary and final peak sets
# ---------------------------------------------------------------------------

rule stage5b_summary:
    output:
        summary      = f"{OUTDIR}/experiments/{{experiment}}/06_idr/final/reproducibility_summary.tsv",
        conservative = f"{OUTDIR}/experiments/{{experiment}}/06_idr/final/conservative.narrowPeak",
        optimal      = f"{OUTDIR}/experiments/{{experiment}}/06_idr/final/optimal.narrowPeak",
    input:
        true_thresh  = f"{OUTDIR}/experiments/{{experiment}}/06_idr/true_replicates/idr.thresholded.narrowPeak",
        pool_thresh  = f"{OUTDIR}/experiments/{{experiment}}/06_idr/pooled_pseudoreplicates/idr.thresholded.narrowPeak",
        self1_thresh = lambda wc: _self_thresh_path(wc.experiment, 0),
        self2_thresh = lambda wc: _self_thresh_path(wc.experiment, 1),
    params:
        experiment = lambda wc: wc.experiment,
        bio_rep_a  = lambda wc: str(sorted(_bioreps_for(wc.experiment, "treatment"))[0]),
        bio_rep_b  = lambda wc: str(sorted(_bioreps_for(wc.experiment, "treatment"))[1]),
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.stage5b.summary.log",
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.summary})" "$(dirname {log})"
        python3 scripts/stage5b_summary.py \
            --true-peaks {input.true_thresh:q} \
            --pooled-peaks {input.pool_thresh:q} \
            --self1-peaks {input.self1_thresh:q} \
            --self2-peaks {input.self2_thresh:q} \
            --experiment {params.experiment:q} \
            --bio-rep-a {params.bio_rep_a} \
            --bio-rep-b {params.bio_rep_b} \
            --output-tsv {output.summary:q} \
            --output-cons {output.conservative:q} \
            --output-opt {output.optimal:q} \
            2>&1 | tee {log:q}
        """
