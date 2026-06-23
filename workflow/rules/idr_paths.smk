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
