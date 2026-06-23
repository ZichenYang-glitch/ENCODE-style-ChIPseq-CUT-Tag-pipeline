# consensus.smk — Stage 62 consensus DAG integration
# ======================================================
# Per-biorep standard MACS3 peak calls + compute_consensus.py + final promotion.
# All helpers use _consensus_ prefix to avoid Snakemake shared-namespace
# collisions.
#
# Modes in scope:
#   chipseq MACS3 narrow, chipseq MACS3 broad,
#   cuttag MACS3 narrow, cuttag MACS3 broad,
#   atac MACS3 narrow
#
# Stage 63 adds SEACR consensus for CUT&Tag PE experiments.
#
# Rules:
#   consensus_macs3_biorep_peaks  — standard MACS3 on biorep BAM
#   consensus_compute_narrow      — compute_consensus.py on narrowPeak inputs
#   consensus_compute_broad       — compute_consensus.py on broadPeak inputs
#   consensus_final               — promote to final/ for primary-consensus modes


# ---------------------------------------------------------------------------
# Helper: MACS3 args for consensus per-biorep peak calls
# ---------------------------------------------------------------------------

def _consensus_macs3_args(wildcards):
    """Return standard assay-specific MACS3 args for consensus biorep peak call.

    Uses the same layout/genome as per-sample MACS3 via get_macs3_args().
    This produces STANDARD q-value calls (-q 0.01), NOT IDR relaxed p-value.
    """
    experiment = wildcards.experiment
    treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    if not treatment_ids:
        return ""
    # Dispatch through the standard assay-specific args function
    class _WC:
        sample = treatment_ids[0]
    return get_macs3_args(_WC)


# ---------------------------------------------------------------------------
# Helper: validate consensus wildcards against experiment metadata
# ---------------------------------------------------------------------------

def _consensus_validate_mode(wildcards):
    """Validate assay/peak_mode/suffix wildcards for a consensus experiment."""
    experiment = wildcards.experiment
    treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    if not treatment_ids:
        raise ValueError(f"Consensus experiment has no treatment samples: {experiment}")

    first = SAMPLE_MAP[treatment_ids[0]]
    actual_assay = first["assay"]
    actual_peak_mode = first["peak_mode"]
    expected_suffix = "narrowPeak" if actual_peak_mode == "narrow" else "broadPeak"

    if (actual_assay, actual_peak_mode) not in CONSENSUS_MODES:
        raise ValueError(
            f"Experiment {experiment} is not eligible for MACS3 consensus: "
            f"{actual_assay}/{actual_peak_mode}"
        )
    if wildcards.assay != actual_assay:
        raise ValueError(
            f"Consensus assay wildcard mismatch for {experiment}: "
            f"target has {wildcards.assay}, sample sheet has {actual_assay}"
        )
    if wildcards.peak_mode != actual_peak_mode:
        raise ValueError(
            f"Consensus peak_mode wildcard mismatch for {experiment}: "
            f"target has {wildcards.peak_mode}, sample sheet has {actual_peak_mode}"
        )
    if wildcards.suffix != expected_suffix:
        raise ValueError(
            f"Consensus suffix wildcard mismatch for {experiment}: "
            f"target has {wildcards.suffix}, expected {expected_suffix}"
        )


# ---------------------------------------------------------------------------
# Helper: inputs for consensus per-biorep MACS3
# ---------------------------------------------------------------------------

def _consensus_macs3_inputs(wildcards):
    """Return inputs for consensus per-biorep MACS3 on a biorep BAM.

    Same pooled-control policy as IDR biorep helpers: if the experiment has
    controls, the pooled control BAM is added as -c.
    """
    _consensus_validate_mode(wildcards)
    exp = wildcards.experiment
    br = int(wildcards.bio_rep)
    inputs = [
        f"{OUTDIR}/experiments/{exp}/02_align/biorep{br}.final.bam",
        f"{OUTDIR}/experiments/{exp}/02_align/biorep{br}.final.bam.bai",
    ]
    if exp in POOLED_CONTROL_EXPERIMENTS:
        inputs.append(
            f"{OUTDIR}/experiments/{exp}/02_align/"
            f"{exp}.pooled.control.final.bam"
        )
    return inputs


# ---------------------------------------------------------------------------
# Helper: resolve per-biorep consensus peak files for compute_consensus.py
# ---------------------------------------------------------------------------

def _consensus_biorep_peaks(experiment, assay, peak_mode):
    """Return sorted list of per-biorep consensus peak files for an experiment."""
    bioreps = sorted(_bioreps_for(experiment, "treatment"))
    suffix = "narrowPeak" if peak_mode == "narrow" else "broadPeak"
    return [
        f"{OUTDIR}/experiments/{experiment}/06_reproducibility/consensus/"
        f"biorep_peaks/{experiment}.biorep{br}.{assay}.macs3."
        f"{peak_mode}_peaks.{suffix}"
        for br in bioreps
    ]


# ---------------------------------------------------------------------------
# Helper: consensus final output for primary-consensus modes
# ---------------------------------------------------------------------------

def _consensus_final_output(experiment, assay, peak_mode):
    """Return the final consensus output path, or empty string if not primary."""
    # Modes where consensus is primary final:
    #   chipseq broad (when IDR not enabled), cuttag narrow (when IDR not enabled),
    #   cuttag broad (when IDR not enabled)
    # Modes where consensus is NOT primary final:
    #   chipseq narrow (legacy stage5 IDR), atac narrow (policy: wait for IDR),
    #   cuttag narrow (when IDR enabled),
    #   chipseq broad / cuttag broad (when broad IDR enabled)
    if (assay == "chipseq" and peak_mode == "broad"
        and BROAD_CHIPSEQ_IDR_ENABLED
        and experiment in BROAD_CHIPSEQ_IDR_EXPERIMENTS):
        return (
            f"{OUTDIR}/experiments/{experiment}/06_reproducibility/final/"
            f"{experiment}.chipseq.macs3.broad."
            f"replicate_validated.idr.broadPeak"
        )
    if (assay == "cuttag" and peak_mode == "broad"
        and BROAD_CUTTAG_IDR_ENABLED
        and experiment in BROAD_CUTTAG_IDR_EXPERIMENTS):
        return (
            f"{OUTDIR}/experiments/{experiment}/06_reproducibility/final/"
            f"{experiment}.cuttag.macs3.broad."
            f"replicate_validated.idr.broadPeak"
        )
    if (assay == "cuttag" and peak_mode == "narrow"
        and CUTTAG_IDR_ENABLED
        and experiment in CUTTAG_IDR_EXPERIMENTS):
        return (
            f"{OUTDIR}/experiments/{experiment}/06_reproducibility/final/"
            f"{experiment}.cuttag.macs3.narrow."
            f"replicate_validated.idr.narrowPeak"
        )
    if (assay == "chipseq" and peak_mode == "broad") or \
       (assay == "cuttag" and peak_mode == "narrow") or \
       (assay == "cuttag" and peak_mode == "broad"):
        suffix = "narrowPeak" if peak_mode == "narrow" else "broadPeak"
        return (
            f"{OUTDIR}/experiments/{experiment}/06_reproducibility/final/"
            f"{experiment}.{assay}.macs3.{peak_mode}."
            f"replicate_validated.consensus.{suffix}"
        )
    return ""


# ============================================================================
# 1. Per-biorep standard MACS3 peak call
# ============================================================================

rule consensus_macs3_biorep_peaks:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/consensus/"
        f"biorep_peaks/{{experiment}}.biorep{{bio_rep}}.{{assay}}."
        f"macs3.{{peak_mode}}_peaks.{{suffix}}",
    input:
        lambda wc: _consensus_macs3_inputs(wc),
    params:
        macs3_args = lambda wc: _consensus_macs3_args(wc),
        sample     = lambda wc: (
            f"{wc.experiment}.biorep{wc.bio_rep}."
            f"{wc.assay}.macs3.{wc.peak_mode}"
        ),
    wildcard_constraints:
        assay     = r"chipseq|cuttag|atac",
        peak_mode = r"narrow|broad",
        bio_rep   = r"\d+",
        suffix    = r"narrowPeak|broadPeak",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}.biorep{{bio_rep}}.{{assay}}."
        f"macs3.{{peak_mode}}.{{suffix}}.consensus.macs3.log",
    threads: THREADS,
    conda:
        "../envs/macs3.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})" "$(dirname {log:q})"

        set -- {input:q}
        TREATMENT="$1"
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

        # MACS3 writes the name_peaks.suffix file. The canonical output path
        # uses that same name, so leave it in place when it already exists.
        MACS3_OUT="$(dirname {output:q})/{params.sample}_peaks.{wildcards.suffix}"
        if [[ -f {output:q} ]]; then
            :
        elif [[ -f "$MACS3_OUT" ]]; then
            mv "$MACS3_OUT" {output:q}
        else
            # Try the non-_peaks variant (MACS3 naming differs by version)
            ALT_OUT="$(dirname {output:q})/{params.sample}.{wildcards.suffix}"
            if [[ -f "$ALT_OUT" ]]; then
                mv "$ALT_OUT" {output:q}
            else
                echo "ERROR: Expected MACS3 output not found: $MACS3_OUT" >&2
                exit 1
            fi
        fi
        """


# ============================================================================
# 2. Consensus computation via compute_consensus.py
# ============================================================================

rule consensus_compute_narrow:
    output:
        peak    = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"consensus/{{experiment}}.{{assay}}.macs3.{{peak_mode}}."
                  f"consensus.narrowPeak",
        summary = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"consensus/{{experiment}}.{{assay}}.macs3.{{peak_mode}}."
                  f"consensus.summary.tsv",
    input:
        lambda wc: _consensus_biorep_peaks(
            wc.experiment, wc.assay, wc.peak_mode
        ),
    params:
        min_replicates     = lambda wc: REPRODUCIBILITY_CONFIG.get(
            "consensus", {}).get("min_replicates", 2),
        reciprocal_overlap = lambda wc: REPRODUCIBILITY_CONFIG.get(
            "consensus", {}).get("reciprocal_overlap", 0.5),
        bioreps_args       = lambda wc: " ".join(
            str(br) for br in sorted(_bioreps_for(wc.experiment, "treatment"))
        ),
        n_bioreps          = lambda wc: len(
            _bioreps_for(wc.experiment, "treatment")),
        final_method       = lambda wc: _consensus_final_method(
            wc.experiment, wc.assay, wc.peak_mode),
        final_output       = lambda wc: _consensus_final_output(
            wc.experiment, wc.assay, wc.peak_mode),
    wildcard_constraints:
        assay     = r"chipseq|cuttag|atac",
        peak_mode = r"narrow",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}.{{assay}}.macs3.{{peak_mode}}."
        f"narrowPeak.consensus.log",
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.peak:q})" "$(dirname {log:q})"

        python3 scripts/compute_consensus.py \
            --peaks {input:q} \
            --bioreps {params.bioreps_args} \
            --format narrowPeak \
            --min-replicates {params.min_replicates} \
            --reciprocal-overlap {params.reciprocal_overlap} \
            --output {output.peak:q} \
            --summary {output.summary:q} \
            --experiment {wildcards.experiment:q} \
            --assay {wildcards.assay:q} \
            --caller macs3 \
            --peak-mode {wildcards.peak_mode:q} \
            --final-method {params.final_method:q} \
            --final-output {params.final_output:q} \
            2>&1 | tee {log:q}
        """


rule consensus_compute_broad:
    output:
        peak    = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"consensus/{{experiment}}.{{assay}}.macs3.{{peak_mode}}."
                  f"consensus.broadPeak",
        summary = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"consensus/{{experiment}}.{{assay}}.macs3.{{peak_mode}}."
                  f"consensus.summary.tsv",
    input:
        lambda wc: _consensus_biorep_peaks(
            wc.experiment, wc.assay, wc.peak_mode
        ),
    params:
        min_replicates     = lambda wc: REPRODUCIBILITY_CONFIG.get(
            "consensus", {}).get("min_replicates", 2),
        reciprocal_overlap = lambda wc: REPRODUCIBILITY_CONFIG.get(
            "consensus", {}).get("reciprocal_overlap", 0.5),
        bioreps_args       = lambda wc: " ".join(
            str(br) for br in sorted(_bioreps_for(wc.experiment, "treatment"))
        ),
        n_bioreps          = lambda wc: len(
            _bioreps_for(wc.experiment, "treatment")),
        final_method       = lambda wc: _consensus_final_method(
            wc.experiment, wc.assay, wc.peak_mode),
        final_output       = lambda wc: _consensus_final_output(
            wc.experiment, wc.assay, wc.peak_mode),
    wildcard_constraints:
        assay     = r"chipseq|cuttag",
        peak_mode = r"broad",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}.{{assay}}.macs3.{{peak_mode}}."
        f"broadPeak.consensus.log",
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.peak:q})" "$(dirname {log:q})"

        python3 scripts/compute_consensus.py \
            --peaks {input:q} \
            --bioreps {params.bioreps_args} \
            --format broadPeak \
            --min-replicates {params.min_replicates} \
            --reciprocal-overlap {params.reciprocal_overlap} \
            --output {output.peak:q} \
            --summary {output.summary:q} \
            --experiment {wildcards.experiment:q} \
            --assay {wildcards.assay:q} \
            --caller macs3 \
            --peak-mode {wildcards.peak_mode:q} \
            --final-method {params.final_method:q} \
            --final-output {params.final_output:q} \
            2>&1 | tee {log:q}
        """


# ============================================================================
# 3. Final consensus promotion (primary-consensus modes only)
# ============================================================================

rule consensus_final:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/final/"
        f"{{experiment}}.{{assay}}.macs3.{{peak_mode}}."
        f"replicate_validated.consensus.{{suffix}}",
    input:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/consensus/"
        f"{{experiment}}.{{assay}}.macs3.{{peak_mode}}.consensus.{{suffix}}",
    wildcard_constraints:
        assay     = r"chipseq|cuttag|atac",
        peak_mode = r"narrow|broad",
        suffix    = r"narrowPeak|broadPeak",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})"
        cp {input:q} {output:q}
        """


# ---------------------------------------------------------------------------
# Helper: determine final_method for consensus summary
# ---------------------------------------------------------------------------

def _consensus_final_method(experiment, assay, peak_mode):
    """Return the final_method string for compute_consensus.py --final-method."""
    if assay == "chipseq" and peak_mode == "narrow":
        if STAGE5 and experiment in IDR_EXPERIMENTS:
            return "idr"
        return "none"
    if assay == "chipseq" and peak_mode == "broad":
        if (BROAD_CHIPSEQ_IDR_ENABLED
            and experiment in BROAD_CHIPSEQ_IDR_EXPERIMENTS):
            return "idr"
        return "consensus"
    if assay == "cuttag" and peak_mode == "broad":
        if (BROAD_CUTTAG_IDR_ENABLED
            and experiment in BROAD_CUTTAG_IDR_EXPERIMENTS):
            return "idr"
        return "consensus"
    if assay == "cuttag" and peak_mode == "narrow":
        if CUTTAG_IDR_ENABLED and experiment in CUTTAG_IDR_EXPERIMENTS:
            return "idr"
        return "consensus"
    if assay == "atac" and peak_mode == "narrow":
        if ATAC_IDR_ENABLED and experiment in ATAC_IDR_EXPERIMENTS:
            return "idr"
        return "none"
    # chipseq narrow → legacy stage5 or consensus
    return "consensus"


# ============================================================================
# Stage 63: SEACR consensus
# ============================================================================
# SEACR consensus is enabled only when ALL of:
#   - REPRO_ENABLED and CONSENSUS_ENABLED
#   - SEACR_ENABLED (cuttag.seacr.enabled: true)
#   - STAGE4B
#   - assay == cuttag, layout == PE, >= 2 biological replicates
#
# Helpers use _consensus_seacr_ prefix.


# ---------------------------------------------------------------------------
# Helper: SEACR consensus validation
# ---------------------------------------------------------------------------

def _consensus_seacr_validate(wildcards):
    """Validate SEACR consensus wildcards against config and metadata."""
    if wildcards.experiment not in SEACR_CONSENSUS_EXPERIMENTS:
        raise ValueError(
            f"Experiment {wildcards.experiment} is not eligible for SEACR consensus"
        )
    if hasattr(wildcards, "mode") and wildcards.mode != SEACR_MODE:
        raise ValueError(
            f"SEACR consensus mode mismatch for {wildcards.experiment}: "
            f"target has {wildcards.mode}, config has {SEACR_MODE}"
        )


def _consensus_seacr_bedgraph_input(wildcards):
    """Return per-biorep BAM input for SEACR consensus bedGraph."""
    _consensus_seacr_validate(wildcards)
    return (
        f"{OUTDIR}/experiments/{wildcards.experiment}/02_align/"
        f"biorep{wildcards.bio_rep}.final.bam"
    )


def _consensus_seacr_peak_inputs(wildcards):
    """Return sorted per-biorep SEACR BED inputs for compute_consensus.py."""
    _consensus_seacr_validate(wildcards)
    return [
        f"{OUTDIR}/experiments/{wildcards.experiment}/06_reproducibility/"
        f"consensus/biorep_peaks/{wildcards.experiment}.biorep{br}."
        f"cuttag.seacr.{wildcards.mode}.bed"
        for br in sorted(_bioreps_for(wildcards.experiment, "treatment"))
    ]


# ---------------------------------------------------------------------------
# Helper: SEACR consensus final output path
# ---------------------------------------------------------------------------

def _consensus_seacr_final_output(experiment):
    """Return the SEACR consensus final BED path."""
    return (
        f"{OUTDIR}/experiments/{experiment}/06_reproducibility/final/"
        f"{experiment}.cuttag.seacr.{SEACR_MODE}."
        f"replicate_validated.consensus.bed"
    )


# ---------------------------------------------------------------------------
# 5. Per-biorep SEACR bedGraph
# ---------------------------------------------------------------------------

rule consensus_seacr_bedgraph:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/consensus/"
        f"biorep_peaks/{{experiment}}.biorep{{bio_rep}}.cuttag.seacr.bedgraph",
    input:
        lambda wc: _consensus_seacr_bedgraph_input(wc),
    wildcard_constraints:
        bio_rep = r"\d+",
    conda:
        "../envs/seacr.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})"
        bedtools genomecov -ibam {input:q} -bg -pc > {output:q}
        """


# ---------------------------------------------------------------------------
# 6. Per-biorep SEACR peak call
# ---------------------------------------------------------------------------

rule consensus_seacr_biorep_peaks:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/consensus/"
        f"biorep_peaks/{{experiment}}.biorep{{bio_rep}}."
        f"cuttag.seacr.{{mode}}.bed",
    input:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/consensus/"
        f"biorep_peaks/{{experiment}}.biorep{{bio_rep}}.cuttag.seacr.bedgraph",
    params:
        threshold = SEACR_THRESHOLD,
        norm      = SEACR_NORMALIZATION,
    wildcard_constraints:
        bio_rep = r"\d+",
        mode    = r"stringent|relaxed",
    conda:
        "../envs/seacr.yml",
    shell:
        """
        set -e -o pipefail
        SEACR_SCRIPT=$(command -v SEACR_1.3.sh) || {{
            echo "ERROR: SEACR_1.3.sh not found in PATH" >&2; exit 1; }}
        mkdir -p "$(dirname {output:q})"

        # SEACR writes prefix.mode.bed.
        TMP_OUTDIR="$(dirname {output:q})"
        TMP_PREFIX="$TMP_OUTDIR/{wildcards.experiment}.biorep{wildcards.bio_rep}.cuttag.seacr"

        bash "$SEACR_SCRIPT" {input:q} {params.threshold} \
            {params.norm:q} {wildcards.mode:q} "$TMP_PREFIX"

        # SEACR_1.3.sh writes prefix.mode.bed, which is also the canonical path.
        SEACR_BED="$TMP_PREFIX.{wildcards.mode}.bed"
        if [[ -f {output:q} ]]; then
            :
        elif [[ -f "$SEACR_BED" ]]; then
            mv "$SEACR_BED" {output:q}
        else
            # Fallback: SEACR may append .bed to the prefix (resulting in .bed.bed)
            ALT_BED="$TMP_PREFIX.bed"
            if [[ -f "$ALT_BED" ]]; then
                mv "$ALT_BED" {output:q}
            else
                echo "ERROR: SEACR output not found: $SEACR_BED or $ALT_BED" >&2
                exit 1
            fi
        fi
        """


# ---------------------------------------------------------------------------
# 7. SEACR consensus computation
# ---------------------------------------------------------------------------

rule consensus_compute_seacr:
    output:
        peak    = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"consensus/{{experiment}}.cuttag.seacr.{{mode}}."
                  f"consensus.bed",
        summary = f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/"
                  f"consensus/{{experiment}}.cuttag.seacr.{{mode}}."
                  f"consensus.summary.tsv",
    input:
        lambda wc: _consensus_seacr_peak_inputs(wc),
    params:
        min_replicates     = lambda wc: REPRODUCIBILITY_CONFIG.get(
            "consensus", {}).get("min_replicates", 2),
        reciprocal_overlap = lambda wc: REPRODUCIBILITY_CONFIG.get(
            "consensus", {}).get("reciprocal_overlap", 0.5),
        bioreps_args       = lambda wc: " ".join(
            str(br) for br in sorted(_bioreps_for(wc.experiment, "treatment"))
        ),
        final_output       = lambda wc: _consensus_seacr_final_output(
            wc.experiment),
    wildcard_constraints:
        mode = r"stringent|relaxed",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/"
        f"{{experiment}}.cuttag.seacr.{{mode}}."
        f"consensus.log",
    conda:
        "../envs/python.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output.peak:q})" "$(dirname {log:q})"

        python3 scripts/compute_consensus.py \
            --peaks {input:q} \
            --bioreps {params.bioreps_args} \
            --format bed \
            --min-replicates {params.min_replicates} \
            --reciprocal-overlap {params.reciprocal_overlap} \
            --output {output.peak:q} \
            --summary {output.summary:q} \
            --experiment {wildcards.experiment:q} \
            --assay cuttag \
            --caller seacr \
            --peak-mode {wildcards.mode:q} \
            --final-method consensus \
            --final-output {params.final_output:q} \
            2>&1 | tee {log:q}
        """


# ---------------------------------------------------------------------------
# 8. SEACR consensus final promotion
# ---------------------------------------------------------------------------

rule consensus_final_seacr:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/final/"
        f"{{experiment}}.cuttag.seacr.{{mode}}."
        f"replicate_validated.consensus.bed",
    input:
        f"{OUTDIR}/experiments/{{experiment}}/06_reproducibility/consensus/"
        f"{{experiment}}.cuttag.seacr.{{mode}}.consensus.bed",
    wildcard_constraints:
        mode = r"stringent|relaxed",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})"
        cp {input:q} {output:q}
        """
