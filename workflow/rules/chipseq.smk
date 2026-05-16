# chipseq.smk — ChIP-seq biological policy functions
# =====================================================
# No rules in this file. These functions define assay-specific MACS3
# parameters, duplicate removal policy, and read extension policy.
# They are consumed by dispatch wrappers in workflow/Snakefile and
# called from rules in common.smk / peaks.smk.


def get_macs3_args_chipseq(wildcards):
    """Standard ChIP-seq MACS3 parameters.

    Narrow: standard model-based peak calling.
    Broad:  --broad --broad-cutoff 0.1
    The -c (control) flag is built from tracked input in peaks.smk,
    not from this function.
    """
    s = SAMPLE_MAP[wildcards.sample]
    fmt = "BAMPE" if s["layout"] == "PE" else "BAM"
    genome = _normalize_genome(s["genome"])
    qvalue = _tool_param("macs3", "qvalue", 0.01)
    if s["peak_mode"] == "broad":
        broad_cutoff = _tool_param("macs3", "broad_cutoff", 0.1)
        return f"-f {fmt} -g {genome} -q {qvalue} --broad --broad-cutoff {broad_cutoff}"
    else:
        return f"-f {fmt} -g {genome} -q {qvalue}"


def get_remove_dup_chipseq(wildcards):
    """ChIP-seq duplicate removal policy.

    config remove_dup=auto  → yes for narrow peaks, no for broad peaks.
    config remove_dup=yes/no → returned directly.
    """
    s = SAMPLE_MAP[wildcards.sample]
    policy = config.get("remove_dup", "auto")
    if policy == "auto":
        return "yes" if s["peak_mode"] == "narrow" else "no"
    return policy


def get_extend_reads_chipseq(wildcards):
    """ChIP-seq read extension policy for bamCoverage.

    extend_reads=auto:  PE → no extension; SE → MACS3 frag size or 200.
    extend_reads=yes:   PE → --extendReads (pair inference); SE → same as auto.
    extend_reads=no:    no extension.
    extend_reads=<int>: --extendReads <int>.
    """
    s = SAMPLE_MAP[wildcards.sample]
    ext = str(config.get("extend_reads", "auto"))
    is_pe = s["layout"] == "PE"

    if ext.isdigit():
        return f"--extendReads {ext}"
    if ext == "no":
        return ""

    if is_pe:
        if ext == "auto":
            return ""
        else:  # yes
            return "--extendReads"
    else:
        # SE mode: try MACS3 predicted fragment size, fallback 200
        if ext in ("auto", "yes"):
            frag = _extract_macs3_fragment_size(wildcards)
            return f"--extendReads {frag if frag else '200'}"
        return ""


def _extract_macs3_fragment_size(wildcards):
    """Read predicted fragment length from MACS3 log file.

    Returns fragment size as a string, or None if unavailable.
    Only realistic sizes (>= 60 bp) are accepted.
    """
    log_path = f"{OUTDIR}/{wildcards.sample}/logs/{wildcards.sample}.macs3.log"
    try:
        if not os.path.exists(log_path):
            return None
        with open(log_path) as fh:
            for line in fh:
                m = re.search(r"predicted fragment length is (\d+)", line)
                if m:
                    size = int(m.group(1))
                    if size >= 60:
                        return str(size)
    except Exception:
        pass
    return None
