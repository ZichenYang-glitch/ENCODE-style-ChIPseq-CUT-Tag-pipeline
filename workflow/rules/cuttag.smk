# cuttag.smk — CUT&Tag biological policy functions
# ==================================================
# No rules in this file. These functions define assay-specific MACS3
# parameters, duplicate removal policy, and read extension policy.
# They are consumed by dispatch wrappers in workflow/Snakefile and
# called from rules in common.smk / peaks.smk.
#
# Stage 1: duplicate and extension policies are identical to ChIP-seq.
# Future stages may diverge (e.g., keep duplicates for CUT&Tag,
# shorter extension lengths, SEACR/GoPeaks hooks).


def get_macs3_args_cuttag(wildcards):
    """CUT&Tag MACS3 parameters with Tn5-aware shift (narrow peaks only).

    Narrow: --nomodel --shift -100 --extsize 200
    Broad:  --broad --broad-cutoff 0.1 (no Tn5 shift per ENCODE guidelines)
    The -c (control) flag is built from tracked input in peaks.smk,
    not from this function.
    """
    s = SAMPLE_MAP[wildcards.sample]
    fmt = "BAMPE" if s["layout"] == "PE" else "BAM"
    genome = _normalize_genome(s["genome"])
    if s["peak_mode"] == "narrow":
        return (
            f"-f {fmt} -g {genome} -q 0.01 "
            "--nomodel --shift -100 --extsize 200"
        )
    else:
        return f"-f {fmt} -g {genome} -q 0.01 --broad --broad-cutoff 0.1"


def get_remove_dup_cuttag(wildcards):
    """CUT&Tag duplicate removal policy.

    Stage 1: identical to ChIP-seq (auto → peak_mode-based).
    Future: CUT&Tag typically has few PCR duplicates; keeping them
    may improve sensitivity for TF marks.
    """
    return get_remove_dup_chipseq(wildcards)


def get_extend_reads_cuttag(wildcards):
    """CUT&Tag read extension policy for bamCoverage.

    Stage 1: identical to ChIP-seq.
    Future: CUT&Tag fragments are shorter (~150-300 bp for nucleosomes,
    ~50-80 bp for TFs). Extension defaults may be adjusted.
    """
    return get_extend_reads_chipseq(wildcards)
