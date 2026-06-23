# Stage 65: IDR path helpers (zero-behavior-change pilot).
# Included into Snakefile namespace; references global OUTDIR.
# No rule declarations. No imports. No side effects.


def idr_biorep_bam(experiment, bio_rep):
    """Treatment biorep merged BAM path used by IDR rules."""
    return f"{OUTDIR}/experiments/{experiment}/02_align/biorep{bio_rep}.final.bam"


def idr_biorep_bai(experiment, bio_rep):
    """Treatment biorep merged BAM index path used by IDR rules."""
    return f"{OUTDIR}/experiments/{experiment}/02_align/biorep{bio_rep}.final.bam.bai"


def idr_pooled_treatment_bam(experiment):
    """Pooled treatment BAM path used by IDR pseudorep rules."""
    return f"{OUTDIR}/experiments/{experiment}/02_align/{experiment}.pooled.final.bam"


def idr_pooled_control_bam(experiment):
    """Pooled control BAM path used when controls are present."""
    return f"{OUTDIR}/experiments/{experiment}/02_align/{experiment}.pooled.control.final.bam"


def idr_pseudorep_bam(experiment, source, pr):
    """Pseudoreplicate BAM path used by IDR rules."""
    return f"{OUTDIR}/experiments/{experiment}/05_pseudorep/{experiment}_{source}.pr{pr}.bam"


def idr_pseudorep_bai(experiment, source, pr):
    """Pseudoreplicate BAM index path used by IDR rules."""
    return f"{OUTDIR}/experiments/{experiment}/05_pseudorep/{experiment}_{source}.pr{pr}.bam.bai"


def idr_biorep_peaks_inputs(experiment, bio_rep):
    """Input BAM paths for a single biorep IDR-ready MACS3 call.

    Order: [treatment_bam, treatment_bai, optional_pooled_control_bam].
    The pooled control is appended only when the experiment has controls,
    ensuring Snakemake schedules pool_control_bam before this rule.
    """
    inputs = [
        idr_biorep_bam(experiment, bio_rep),
        idr_biorep_bai(experiment, bio_rep),
    ]
    if experiment in POOLED_CONTROL_EXPERIMENTS:
        inputs.append(idr_pooled_control_bam(experiment))
    return inputs


def idr_split_input(experiment, source):
    """Input BAM path for pseudorep splitting.

    source == "pooled" -> pooled treatment BAM
    source starts with "biorep" -> the matching biorep treatment BAM
    """
    if source == "pooled":
        return idr_pooled_treatment_bam(experiment)
    br_label = source[len("biorep"):]
    return idr_biorep_bam(experiment, br_label)


def idr_pseudorep_peaks_inputs(experiment, source, pr, source_prefix=""):
    """Input BAM paths for a pseudorep IDR-ready MACS3 call.

    source_prefix lets ATAC/CUT&Tag/broad variants prepend their assay
    infix to the pseudorep source name (e.g. "atac_", "cuttag_",
    "broad_chipseq_") without changing the path helper contract.
    """
    prefixed = f"{source_prefix}{source}" if source_prefix else source
    inputs = [
        idr_pseudorep_bam(experiment, prefixed, pr),
        idr_pseudorep_bai(experiment, prefixed, pr),
    ]
    if experiment in POOLED_CONTROL_EXPERIMENTS:
        inputs.append(idr_pooled_control_bam(experiment))
    return inputs
