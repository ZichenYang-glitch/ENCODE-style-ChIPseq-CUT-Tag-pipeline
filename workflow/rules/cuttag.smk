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

    Narrow SE (-f BAM):   --nomodel --shift -100 --extsize 200
    Narrow PE (-f BAMPE): no shift/extension flags. This is a deliberate
                          decision, not a parameter fix: under BAMPE, MACS3
                          neutralizes --shift to 0 and, with --nomodel,
                          ignores --extsize (using the real pair distance),
                          so the flags would be silently inert.
    Broad:  --broad --broad-cutoff 0.1 (no Tn5 shift per ENCODE guidelines)
    The -c (control) flag is built from tracked input in peaks.smk,
    not from this function.
    """
    s = SAMPLE_MAP[wildcards.sample]
    fmt = "BAMPE" if s["layout"] == "PE" else "BAM"
    genome = _normalize_genome(s["genome"])
    qvalue = _tool_param("macs3", "qvalue", 0.01)
    if s["peak_mode"] == "narrow":
        args = f"-f {fmt} -g {genome} -q {qvalue}"
        if s["layout"] != "PE":
            args += " --nomodel --shift -100 --extsize 200"
        return args
    else:
        broad_cutoff = _tool_param("macs3", "broad_cutoff", 0.1)
        return f"-f {fmt} -g {genome} -q {qvalue} --broad --broad-cutoff {broad_cutoff}"


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
