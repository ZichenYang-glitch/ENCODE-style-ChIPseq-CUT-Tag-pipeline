# atac.smk — ATAC-seq biological policy functions
# =================================================
# No rules in this file. These functions define baseline ATAC-seq MACS3
# parameters, duplicate removal policy, and read extension policy.
# They are consumed by dispatch wrappers in workflow/Snakefile and called
# from generic rules in common.smk / peaks.smk.


def get_macs3_args_atac(wildcards):
    """ATAC-seq MACS3 parameters.

    Stage 19 supports narrow peaks only. The Tn5-aware shift/extension policy
    mirrors common ATAC-seq callpeak practice while keeping inputs and outputs
    in the shared MACS3 rule.
    """
    s = SAMPLE_MAP[wildcards.sample]
    fmt = "BAMPE" if s["layout"] == "PE" else "BAM"
    genome = _normalize_genome(s["genome"])
    qvalue = _tool_param("macs3", "qvalue", 0.01)
    return (
        f"-f {fmt} -g {genome} -q {qvalue} "
        "--nomodel --shift -100 --extsize 200"
    )


def get_remove_dup_atac(wildcards):
    """ATAC-seq duplicate removal policy.

    config remove_dup=auto  → yes
    config remove_dup=yes/no → returned directly.
    """
    policy = config.get("remove_dup", "auto")
    if policy == "auto":
        return "yes"
    return policy


def get_extend_reads_atac(wildcards):
    """ATAC-seq read extension policy for bamCoverage.

    For paired-end data, bamCoverage infers fragments from pairs in auto mode.
    For single-end data, auto/yes use a conservative 200 bp extension.
    Explicit integer and no policies match the shared workflow semantics.
    """
    s = SAMPLE_MAP[wildcards.sample]
    ext = str(config.get("extend_reads", "auto"))
    is_pe = s["layout"] == "PE"

    if ext.isdigit():
        return f"--extendReads {ext}"
    if ext == "no":
        return ""
    if is_pe:
        return "" if ext == "auto" else "--extendReads"
    if ext in ("auto", "yes"):
        return "--extendReads 200"
    return ""
