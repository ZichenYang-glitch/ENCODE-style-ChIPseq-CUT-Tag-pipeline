# Stage 44: MNase path helpers — zero-behavior-change pilot.
# Included into Snakefile namespace; references global OUTDIR.
# No rule declarations. No imports. No side effects.


def mnase_fragment_bam(sample, fragment_class):
    """Sample-level MNase fragment BAM path.

    Args:
        sample: Snakemake wildcard string (e.g. "{sample}") or concrete id.
        fragment_class: "mono", "sub", or "di".
    """
    return f"{OUTDIR}/{sample}/03_fragments/{sample}.{fragment_class}.bam"


def mnase_fragment_bai(sample, fragment_class):
    """Sample-level MNase fragment BAM index path."""
    return f"{OUTDIR}/{sample}/03_fragments/{sample}.{fragment_class}.bam.bai"


def mnase_signal_bw(sample, kind):
    """Sample-level MNase BigWig path.

    Args:
        sample: Snakemake wildcard string or concrete id.
        kind: "dyad" or "mono".
    """
    return f"{OUTDIR}/{sample}/04_signal/{sample}.{kind}.CPM.bw"


def mnase_qc_summary_tsv(sample):
    """Sample-level MNase QC summary TSV path."""
    return f"{OUTDIR}/{sample}/01_qc/{sample}.mnase_qc_summary.tsv"


def mnase_pooled_fragment_bam(experiment, fragment_class):
    """Experiment-level pooled MNase fragment BAM path.

    Args:
        experiment: Snakemake wildcard string (e.g. "{experiment}") or concrete id.
        fragment_class: "mono", "sub", or "di".
    """
    return f"{OUTDIR}/experiments/{experiment}/03_fragments/{experiment}.pooled.{fragment_class}.bam"


def mnase_pooled_fragment_bai(experiment, fragment_class):
    """Experiment-level pooled MNase fragment BAM index path."""
    return f"{OUTDIR}/experiments/{experiment}/03_fragments/{experiment}.pooled.{fragment_class}.bam.bai"


def mnase_pooled_signal_bw(experiment, kind):
    """Experiment-level pooled MNase BigWig path.

    Args:
        experiment: Snakemake wildcard string or concrete id.
        kind: "dyad" or "mono".
    """
    return f"{OUTDIR}/experiments/{experiment}/04_signal/{experiment}.pooled.{kind}.CPM.bw"
