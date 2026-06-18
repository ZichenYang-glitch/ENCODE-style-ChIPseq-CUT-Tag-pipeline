
# ---------------------------------------------------------------------------
# 3a. Stage 4a replicate-aware metadata (derived structures, no scheduling changes)
# ---------------------------------------------------------------------------

EXPERIMENT_IDS = sorted(
    {s["experiment"] for s in SAMPLES}
)

SAMPLES_BY_EXPERIMENT: dict[str, list[str]] = {exp: [] for exp in EXPERIMENT_IDS}
for s in SAMPLES:
    SAMPLES_BY_EXPERIMENT[s["experiment"]].append(s["id"])

TREATMENT_SAMPLES_BY_EXPERIMENT: dict[str, list[str]] = {
    exp: [sid for sid in sids if SAMPLE_MAP[sid]["role"] == "treatment"]
    for exp, sids in SAMPLES_BY_EXPERIMENT.items()
}

CONTROL_SAMPLES_BY_EXPERIMENT: dict[str, list[str]] = {
    exp: [sid for sid in sids if SAMPLE_MAP[sid]["role"] == "control"]
    for exp, sids in SAMPLES_BY_EXPERIMENT.items()
}


# ---------------------------------------------------------------------------
# 3a2. Stage 4b replicate-aware grouped outputs
# ---------------------------------------------------------------------------

if STAGE4B:

    # Biological replicate groups: keyed by (experiment, role, bio_rep)
    # Role separation ensures treatment and control BAMs are never merged together.
    BIO_REP_GROUPS: dict[tuple[str, str, int], list[str]] = {}
    for s in SAMPLES:
        key = (s["experiment"], s["role"], s["biological_replicate"])
        BIO_REP_GROUPS.setdefault(key, []).append(s["id"])

    # Biological replicates for a specific experiment + role
    def _bioreps_for(experiment, role):
        return sorted(
            bio_rep for (exp, r, bio_rep) in BIO_REP_GROUPS
            if exp == experiment and r == role
        )

    # Pre-computed (experiment, bio_rep) pairs for treatment expand() calls.
    # Only includes bio-reps that need merging (multiple tech-reps) or are
    # part of a multi-biorep experiment (needed for pooling).
    def _biorep_expand_pairs():
        """Return (experiment_list, bio_rep_list) for expand(zip, ...).

        Gating:
        - bio-reps in multi-biorep experiments: always included (pooling needs them)
        - bio-reps in single-biorep experiments: only included if they have
          multiple tech-reps (need tech-rep merging)
        - single-sample bio-reps with 1 tech-rep: excluded (no-op, no DAG change)
        """
        exp_list = []
        br_list = []
        for exp in EXPERIMENT_IDS:
            for br in _bioreps_for(exp, "treatment"):
                key = (exp, "treatment", br)
                n_techreps = len(BIO_REP_GROUPS.get(key, []))
                n_bioreps = len(_bioreps_for(exp, "treatment"))
                if n_bioreps >= 2 or n_techreps >= 2:
                    exp_list.append(exp)
                    br_list.append(br)
        return exp_list, br_list

    _EXP_LIST, _BR_LIST = _biorep_expand_pairs()

    # Experiments with >=2 treatment biological replicates → pooled BAMs + peaks
    MULTI_BIOREP_EXPERIMENTS = sorted(
        exp for exp in EXPERIMENT_IDS
        if len(_bioreps_for(exp, "treatment")) >= 2
    )

    # Control BAM resolution for a given experiment.
    # Returns list of (source_type, path) tuples where source_type is
    # "control_sample" or "control_bam".
    def _resolve_experiment_controls(experiment):
        """Return unique non-empty control paths for an experiment's treatment rows."""
        treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
        if not treatment_ids:
            return []
        ctrl_samples = set()
        ctrl_bams = set()
        for sid in treatment_ids:
            s = SAMPLE_MAP[sid]
            cs = s.get("control_sample", "")
            cb = s.get("control_bam", "")
            if cs:
                ctrl_samples.add(cs)
            elif cb:
                ctrl_bams.add(cb)
        result = []
        for cs in sorted(ctrl_samples):
            result.append(("control_sample", f"{OUTDIR}/{cs}/02_align/{cs}.final.bam"))
        for cb in sorted(ctrl_bams):
            result.append(("control_bam", cb))
        return result

    # Experiments that have controls (for pooled control BAM/BAI)
    # Pooled controls are produced whenever a multi-biorep treatment
    # experiment has controls actually referenced by its treatment rows
    # (1 unique → symlink, >1 → merge).
    POOLED_CONTROL_EXPERIMENTS = sorted(
        exp for exp in MULTI_BIOREP_EXPERIMENTS
        if USE_CONTROL and _resolve_experiment_controls(exp)
    )

else:
    MULTI_BIOREP_EXPERIMENTS = []
    POOLED_CONTROL_EXPERIMENTS = []

# Stage 39: MNase and peak-centric experiment subsets (subset of MULTI_BIOREP_EXPERIMENTS)
# Used for assay-specific pooled target expansions (peaks vs nucleosome outputs).
PEAK_MULTI_BIOREP_EXPERIMENTS = [
    exp for exp in MULTI_BIOREP_EXPERIMENTS
    if SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[exp][0]]["assay"] != "mnase"
]

MNASE_MULTI_BIOREP_EXPERIMENTS = [
    exp for exp in MULTI_BIOREP_EXPERIMENTS
    if SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[exp][0]]["assay"] == "mnase"
]


# ---------------------------------------------------------------------------
# 3a3. Stage 5a IDR derived structures
# ---------------------------------------------------------------------------

if STAGE5:
    IDR_EXPERIMENTS = []
    IDR_BIOREP_EXP_LIST = []
    IDR_BIOREP_LIST = []

    for exp in MULTI_BIOREP_EXPERIMENTS:
        bioreps = _bioreps_for(exp, "treatment")
        if len(bioreps) == 2:
            first = SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[exp][0]]
            if first["assay"] == "chipseq" and first["peak_mode"] == "narrow":
                IDR_EXPERIMENTS.append(exp)
                for br in sorted(bioreps):
                    IDR_BIOREP_EXP_LIST.append(exp)
                    IDR_BIOREP_LIST.append(br)
else:
    IDR_EXPERIMENTS = []
    IDR_BIOREP_EXP_LIST = []
    IDR_BIOREP_LIST = []

# Stage 5b precomputed expansion lists (gated on STAGE5 and IDR_EXPERIMENTS)
if STAGE5 and IDR_EXPERIMENTS:
    IDR_SPLIT_SOURCE_EXP = []
    IDR_SPLIT_SOURCE_NAME = []
    IDR_PR_PEAK_EXP = []
    IDR_PR_PEAK_SRC = []
    IDR_PR_PEAK_PR = []
    IDR_SELF_EXP = []
    IDR_SELF_BR = []

    for exp in IDR_EXPERIMENTS:
        bioreps = sorted(_bioreps_for(exp, "treatment"))
        br_a = bioreps[0]
        br_b = bioreps[1]

        for br in (br_a, br_b):
            src = f"biorep{br}"
            IDR_SPLIT_SOURCE_EXP.append(exp)
            IDR_SPLIT_SOURCE_NAME.append(src)
            for pr in ("1", "2"):
                IDR_PR_PEAK_EXP.append(exp)
                IDR_PR_PEAK_SRC.append(src)
                IDR_PR_PEAK_PR.append(pr)
            IDR_SELF_EXP.append(exp)
            IDR_SELF_BR.append(br)

        IDR_SPLIT_SOURCE_EXP.append(exp)
        IDR_SPLIT_SOURCE_NAME.append("pooled")
        for pr in ("1", "2"):
            IDR_PR_PEAK_EXP.append(exp)
            IDR_PR_PEAK_SRC.append("pooled")
            IDR_PR_PEAK_PR.append(pr)
else:
    IDR_SPLIT_SOURCE_EXP = []
    IDR_SPLIT_SOURCE_NAME = []
    IDR_PR_PEAK_EXP = []
    IDR_PR_PEAK_SRC = []
    IDR_PR_PEAK_PR = []
    IDR_SELF_EXP = []
    IDR_SELF_BR = []


# ---------------------------------------------------------------------------
# 3a4. Stage 55 ATAC narrow IDR derived structures
# ---------------------------------------------------------------------------

if ATAC_IDR_ENABLED:
    ATAC_IDR_EXPERIMENTS = []
    ATAC_IDR_BIOREP_EXP_LIST = []
    ATAC_IDR_BIOREP_LIST = []

    for exp in MULTI_BIOREP_EXPERIMENTS:
        bioreps = _bioreps_for(exp, "treatment")
        if len(bioreps) == 2:
            first = SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[exp][0]]
            if first["assay"] == "atac" and first["peak_mode"] == "narrow":
                ATAC_IDR_EXPERIMENTS.append(exp)
                for br in sorted(bioreps):
                    ATAC_IDR_BIOREP_EXP_LIST.append(exp)
                    ATAC_IDR_BIOREP_LIST.append(br)
else:
    ATAC_IDR_EXPERIMENTS = []
    ATAC_IDR_BIOREP_EXP_LIST = []
    ATAC_IDR_BIOREP_LIST = []

# Stage 55 pseudorep expansion lists
if ATAC_IDR_ENABLED and ATAC_IDR_EXPERIMENTS:
    ATAC_IDR_SPLIT_SOURCE_EXP = []
    ATAC_IDR_SPLIT_SOURCE_NAME = []
    ATAC_IDR_PR_PEAK_EXP = []
    ATAC_IDR_PR_PEAK_SRC = []
    ATAC_IDR_PR_PEAK_PR = []
    ATAC_IDR_SELF_EXP = []
    ATAC_IDR_SELF_BR = []

    for exp in ATAC_IDR_EXPERIMENTS:
        bioreps = sorted(_bioreps_for(exp, "treatment"))
        br_a = bioreps[0]
        br_b = bioreps[1]

        for br in (br_a, br_b):
            src = f"biorep{br}"
            ATAC_IDR_SPLIT_SOURCE_EXP.append(exp)
            ATAC_IDR_SPLIT_SOURCE_NAME.append(src)
            for pr in ("1", "2"):
                ATAC_IDR_PR_PEAK_EXP.append(exp)
                ATAC_IDR_PR_PEAK_SRC.append(src)
                ATAC_IDR_PR_PEAK_PR.append(pr)
            ATAC_IDR_SELF_EXP.append(exp)
            ATAC_IDR_SELF_BR.append(br)

        ATAC_IDR_SPLIT_SOURCE_EXP.append(exp)
        ATAC_IDR_SPLIT_SOURCE_NAME.append("pooled")
        for pr in ("1", "2"):
            ATAC_IDR_PR_PEAK_EXP.append(exp)
            ATAC_IDR_PR_PEAK_SRC.append("pooled")
            ATAC_IDR_PR_PEAK_PR.append(pr)
else:
    ATAC_IDR_SPLIT_SOURCE_EXP = []
    ATAC_IDR_SPLIT_SOURCE_NAME = []
    ATAC_IDR_PR_PEAK_EXP = []
    ATAC_IDR_PR_PEAK_SRC = []
    ATAC_IDR_PR_PEAK_PR = []
    ATAC_IDR_SELF_EXP = []
    ATAC_IDR_SELF_BR = []


# ---------------------------------------------------------------------------
# 3b. Stage 3 QC configuration and genome resource helpers
# ---------------------------------------------------------------------------

QC_CONFIG = VALIDATED_CONFIG.get("qc", {
    "blacklist_filter": True,
    "frip": True,
    "library_complexity": True,
    "nrf_pbc": True,
    "signal_tracks": True,
    "summary": True,
    "cuttag_fragment_size": True,
    "cross_correlation": False,
    "preseq_complexity": False,
    "picard_metrics": False,
    "tss_enrichment": False,
})

# Stage 4c: structured tool parameters (absent → empty dict, defaults silently)
TOOL_PARAMS = VALIDATED_CONFIG.get("tool_parameters", {})


def _tool_param(tool, key, default):
    """Return a single tool parameter value, or default if not configured."""
    block = TOOL_PARAMS.get(tool, {})
    return block.get(key, default)


def _filter_flags_arg():
    """Build the -F argument for samtools_filter. Preserves hex for default."""
    block = TOOL_PARAMS.get("samtools_filter", {})
    raw = block.get("filter_flags", "")
    if raw == "" or raw is None:
        return "-F 0x904"
    if isinstance(raw, int):
        # Check if the original config had hex representation
        raw_config = (
            VALIDATED_CONFIG.get("tool_parameters", {})
            .get("samtools_filter", {})
            .get("filter_flags", "")
        )
        if isinstance(raw_config, str) and (
            raw_config.strip().startswith("0x")
            or raw_config.strip().startswith("0X")
        ):
            return f"-F {raw_config.strip()}"
        return f"-F {raw}"
    return f"-F {raw}"

# Samples with a configured blacklist path (resource-gated blacklist filtering)
BLACKLIST_SAMPLES = [
    s for s in TREATMENT_SAMPLE_IDS
    if GENOME_RESOURCES.get(SAMPLE_MAP[s]["genome"], {}).get("blacklist", "")
]


def get_genome_resource(sample, key, default=""):
    """Return a genome resource path for a sample, or *default* if unset."""
    genome = SAMPLE_MAP[sample]["genome"]
    return GENOME_RESOURCES.get(genome, {}).get(key, default)


def has_genome_resource(sample, key):
    """Return True if a non-empty genome resource is configured for *sample*."""
    return bool(get_genome_resource(sample, key, ""))


TSS_SAMPLE_IDS = [
    sid for sid in PEAK_SAMPLE_IDS
    if QC_CONFIG.get("tss_enrichment", False)
    and has_genome_resource(sid, "gtf")
]
TSS_GENOMES = sorted({SAMPLE_MAP[sid]["genome"] for sid in TSS_SAMPLE_IDS})


# Stage 22: FE/ppois BigWig gating — requires signal_tracks + chrom_sizes
# has_genome_resource() checks non-empty genome_resources.<genome>.chrom_sizes.
# Path existence for non-empty chrom_sizes is already enforced by validate_samples.py.
SIGNAL_BW_SAMPLE_IDS = [
    sid for sid in PEAK_SAMPLE_IDS
    if QC_CONFIG.get("signal_tracks", True)
    and has_genome_resource(sid, "chrom_sizes")
]

if STAGE4B:
    SIGNAL_BW_EXPERIMENTS = [
        exp for exp in PEAK_MULTI_BIOREP_EXPERIMENTS
        if QC_CONFIG.get("signal_tracks", True)
        and SIGNAL_BW_SAMPLE_IDS
        and has_genome_resource(
            TREATMENT_SAMPLES_BY_EXPERIMENT[exp][0], "chrom_sizes"
        )
    ]
else:
    SIGNAL_BW_EXPERIMENTS = []


def _pooled_chrom_sizes(experiment):
    """Return chrom_sizes path for experiment's first treatment sample genome.

    Assumption: all treatment samples in a pooled experiment share the same
    genome. Mixed-genome pooled experiments are not validated against in
    this slice.
    """
    tids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    if not tids:
        return ""
    return get_genome_resource(tids[0], "chrom_sizes")



# ---------------------------------------------------------------------------
# 5. Dispatch wrappers — wire assay-specific policy to generic names
# ---------------------------------------------------------------------------

def get_remove_dup(wildcards):
    """Return "yes" or "no" for duplicate removal, dispatched by assay."""
    s = SAMPLE_MAP[wildcards.sample]
    if s["assay"] == "chipseq":
        return get_remove_dup_chipseq(wildcards)
    if s["assay"] == "cuttag":
        return get_remove_dup_cuttag(wildcards)
    if s["assay"] == "atac":
        return get_remove_dup_atac(wildcards)
    if s["assay"] == "mnase":
        return get_remove_dup_mnase(wildcards)
    raise ValueError(f"Unsupported assay for duplicate policy: {s['assay']}")


def get_macs3_args(wildcards):
    """Return MACS3 arguments string, dispatched by assay."""
    s = SAMPLE_MAP[wildcards.sample]
    if s["assay"] == "chipseq":
        return get_macs3_args_chipseq(wildcards)
    if s["assay"] == "cuttag":
        return get_macs3_args_cuttag(wildcards)
    if s["assay"] == "atac":
        return get_macs3_args_atac(wildcards)
    if s["assay"] == "mnase":
        return get_macs3_args_mnase(wildcards)
    raise ValueError(f"Unsupported assay for MACS3 policy: {s['assay']}")


def get_extend_reads(wildcards):
    """Return bamCoverage extension arguments, dispatched by assay."""
    s = SAMPLE_MAP[wildcards.sample]
    if s["assay"] == "chipseq":
        return get_extend_reads_chipseq(wildcards)
    if s["assay"] == "cuttag":
        return get_extend_reads_cuttag(wildcards)
    if s["assay"] == "atac":
        return get_extend_reads_atac(wildcards)
    if s["assay"] == "mnase":
        return get_extend_reads_mnase(wildcards)
    raise ValueError(f"Unsupported assay for bamCoverage policy: {s['assay']}")


# ---------------------------------------------------------------------------
# 6. MACS3 input resolution — tracks control dependencies
# ---------------------------------------------------------------------------

def _macs3_inputs(wildcards):
    """Return input list for macs3_callpeak.

    Elements: [treatment_final.bam, treatment_final.bam.bai, ...optional_control]
    - treatment final.bam + final.bam.bai are always included.
    - use_control gates both FASTQ-based control_sample and external control_bam.
      When false, control_sample and control_bam are silently ignored.
    """
    s = SAMPLE_MAP[wildcards.sample]
    sid = wildcards.sample
    inputs = [
        f"{OUTDIR}/{sid}/02_align/{sid}.final.bam",
        f"{OUTDIR}/{sid}/02_align/{sid}.final.bam.bai",
    ]
    if USE_CONTROL:
        if s.get("control_sample"):
            cs = s["control_sample"]
            inputs.append(
                f"{OUTDIR}/{cs}/02_align/{cs}.final.bam"
            )
        elif s.get("control_bam"):
            inputs.append(s["control_bam"])
    return inputs
