# ---------------------------------------------------------------------------
# 3a. Stage 4a replicate-aware metadata (derived structures, no scheduling changes)
# ---------------------------------------------------------------------------

EXPERIMENT_IDS = sorted({s["experiment"] for s in SAMPLES})

SAMPLES_BY_EXPERIMENT: dict[str, list[str]] = {exp: [] for exp in EXPERIMENT_IDS}
for s in SAMPLES:
    SAMPLES_BY_EXPERIMENT[s["experiment"]].append(s["id"])

TREATMENT_SAMPLES_BY_EXPERIMENT: dict[str, list[str]] = {
    exp: [sid for sid in sids if SAMPLE_MAP[sid]["role"] == "treatment"] for exp, sids in SAMPLES_BY_EXPERIMENT.items()
}

CONTROL_SAMPLES_BY_EXPERIMENT: dict[str, list[str]] = {
    exp: [sid for sid in sids if SAMPLE_MAP[sid]["role"] == "control"] for exp, sids in SAMPLES_BY_EXPERIMENT.items()
}


# ---------------------------------------------------------------------------
# 3a2. Replicate-aware grouped outputs
# ---------------------------------------------------------------------------

if REPLICATE_ANALYSIS:

    # Biological replicate groups: keyed by (experiment, role, bio_rep)
    # Role separation ensures treatment and control BAMs are never merged together.
    BIO_REP_GROUPS: dict[tuple[str, str, int], list[str]] = {}
    for s in SAMPLES:
        key = (s["experiment"], s["role"], s["biological_replicate"])
        BIO_REP_GROUPS.setdefault(key, []).append(s["id"])

    # Biological replicates for a specific experiment + role
    def _bioreps_for(experiment, role):
        return sorted(bio_rep for (exp, r, bio_rep) in BIO_REP_GROUPS if exp == experiment and r == role)

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
    MULTI_BIOREP_EXPERIMENTS = sorted(exp for exp in EXPERIMENT_IDS if len(_bioreps_for(exp, "treatment")) >= 2)

    # Control BAM resolution for a given experiment.
    # Returns list of (source_type, path) tuples where source_type is
    # "control_sample" or "control_bam".
    def _resolve_experiment_controls(experiment):
        """Return unique non-empty control paths for an experiment's treatment rows."""
        if not USE_CONTROL:
            return []
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
                if VALIDATED_CONFIG.get("analysis_profile") == "updated":
                    raise ValueError("updated profile does not yet support external control_bam")
                ctrl_bams.add(cb)
        result = []
        for cs in sorted(ctrl_samples):
            path = (
                f"{OUTDIR}/{cs}/02_align/{cs}.analysis.bam"
                if VALIDATED_CONFIG.get("analysis_profile") == "updated"
                else f"{OUTDIR}/{cs}/02_align/{cs}.final.bam"
            )
            result.append(("control_sample", path))
        for cb in sorted(ctrl_bams):
            result.append(("control_bam", cb))
        return result

    # Experiments that have controls (for pooled control BAM/BAI)
    # Pooled controls are produced whenever a multi-biorep treatment
    # experiment has controls actually referenced by its treatment rows
    # (1 unique → symlink, >1 → merge).
    POOLED_CONTROL_EXPERIMENTS = sorted(
        exp for exp in MULTI_BIOREP_EXPERIMENTS if USE_CONTROL and _resolve_experiment_controls(exp)
    )

else:
    MULTI_BIOREP_EXPERIMENTS = []
    POOLED_CONTROL_EXPERIMENTS = []

# Stage 39: MNase and peak-centric experiment subsets (subset of MULTI_BIOREP_EXPERIMENTS)
# Used for assay-specific pooled target expansions (peaks vs nucleosome outputs).
PEAK_MULTI_BIOREP_EXPERIMENTS = [
    exp for exp in MULTI_BIOREP_EXPERIMENTS if SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[exp][0]]["assay"] != "mnase"
]

MNASE_MULTI_BIOREP_EXPERIMENTS = [
    exp for exp in MULTI_BIOREP_EXPERIMENTS if SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[exp][0]]["assay"] == "mnase"
]


# ---------------------------------------------------------------------------
# 3a3. IDR derived structures
# ---------------------------------------------------------------------------


def _idr_exp_matches(exp, assay, peak_mode):
    """Return True if experiment's first treatment sample matches assay/peak_mode."""
    treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(exp, [])
    if not treatment_ids:
        return False
    first = SAMPLE_MAP[treatment_ids[0]]
    return first["assay"] == assay and first["peak_mode"] == peak_mode


def _build_idr_experiment_lists(enabled, assay, peak_mode, multi_experiments):
    """Build IDR experiment and biorep expansion lists for a single assay profile.

    Returns (experiments, biorep_exp_list, biorep_list).  The three lists are
    empty when disabled or when no experiment matches the profile.
    """
    experiments = []
    biorep_exp_list = []
    biorep_list = []
    if not enabled:
        return experiments, biorep_exp_list, biorep_list

    for exp in multi_experiments:
        if not _idr_exp_matches(exp, assay, peak_mode):
            continue
        bioreps = sorted(_bioreps_for(exp, "treatment"))
        if len(bioreps) != 2:
            continue
        experiments.append(exp)
        for br in bioreps:
            biorep_exp_list.append(exp)
            biorep_list.append(br)
    return experiments, biorep_exp_list, biorep_list


def _build_idr_pseudorep_lists(experiments):
    """Build pseudoreplicate expansion lists for a narrow IDR profile.

    Returns (split_source_exp, split_source_name, pr_peak_exp, pr_peak_src,
             pr_peak_pr, self_exp, self_br).
    """
    split_source_exp = []
    split_source_name = []
    pr_peak_exp = []
    pr_peak_src = []
    pr_peak_pr = []
    self_exp = []
    self_br = []

    for exp in experiments:
        bioreps = sorted(_bioreps_for(exp, "treatment"))
        br_a, br_b = bioreps[0], bioreps[1]

        for br in (br_a, br_b):
            src = f"biorep{br}"
            split_source_exp.append(exp)
            split_source_name.append(src)
            for pr in ("1", "2"):
                pr_peak_exp.append(exp)
                pr_peak_src.append(src)
                pr_peak_pr.append(pr)
            self_exp.append(exp)
            self_br.append(br)

        split_source_exp.append(exp)
        split_source_name.append("pooled")
        for pr in ("1", "2"):
            pr_peak_exp.append(exp)
            pr_peak_src.append("pooled")
            pr_peak_pr.append(pr)

    return (
        split_source_exp,
        split_source_name,
        pr_peak_exp,
        pr_peak_src,
        pr_peak_pr,
        self_exp,
        self_br,
    )


def _build_broad_idr_pseudorep_lists(experiments, assays):
    """Build pseudoreplicate expansion lists for broad IDR, with assay dimension.

    experiments and assays must be aligned (one assay per experiment).  Returns
    (split_source_exp, split_source_name, split_source_assay,
     pr_peak_exp, pr_peak_src, pr_peak_assay, pr_peak_pr,
     self_exp, self_br, self_assay,
     biorep_exp_list, biorep_list, biorep_assay).
    """
    split_source_exp = []
    split_source_name = []
    split_source_assay = []
    pr_peak_exp = []
    pr_peak_src = []
    pr_peak_assay = []
    pr_peak_pr = []
    self_exp = []
    self_br = []
    self_assay = []
    biorep_exp_list = []
    biorep_list = []
    biorep_assay = []

    for exp, assay in zip(experiments, assays):
        bioreps = sorted(_bioreps_for(exp, "treatment"))
        br_a, br_b = bioreps[0], bioreps[1]

        for br in (br_a, br_b):
            src = f"biorep{br}"
            split_source_exp.append(exp)
            split_source_name.append(src)
            split_source_assay.append(assay)
            for pr in ("1", "2"):
                pr_peak_exp.append(exp)
                pr_peak_src.append(src)
                pr_peak_assay.append(assay)
                pr_peak_pr.append(pr)
            self_exp.append(exp)
            self_br.append(br)
            self_assay.append(assay)
            biorep_exp_list.append(exp)
            biorep_list.append(br)
            biorep_assay.append(assay)

        split_source_exp.append(exp)
        split_source_name.append("pooled")
        split_source_assay.append(assay)
        for pr in ("1", "2"):
            pr_peak_exp.append(exp)
            pr_peak_src.append("pooled")
            pr_peak_assay.append(assay)
            pr_peak_pr.append(pr)

    return (
        split_source_exp,
        split_source_name,
        split_source_assay,
        pr_peak_exp,
        pr_peak_src,
        pr_peak_assay,
        pr_peak_pr,
        self_exp,
        self_br,
        self_assay,
        biorep_exp_list,
        biorep_list,
        biorep_assay,
    )


# ChIP-seq true-replicate IDR
IDR_EXPERIMENTS, IDR_BIOREP_EXP_LIST, IDR_BIOREP_LIST = _build_idr_experiment_lists(
    CHIPSEQ_IDR, "chipseq", "narrow", MULTI_BIOREP_EXPERIMENTS
)

# ChIP-seq pseudoreplicate expansion lists
(
    IDR_SPLIT_SOURCE_EXP,
    IDR_SPLIT_SOURCE_NAME,
    IDR_PR_PEAK_EXP,
    IDR_PR_PEAK_SRC,
    IDR_PR_PEAK_PR,
    IDR_SELF_EXP,
    IDR_SELF_BR,
) = _build_idr_pseudorep_lists(IDR_EXPERIMENTS)


# ATAC narrow IDR
ATAC_IDR_EXPERIMENTS, ATAC_IDR_BIOREP_EXP_LIST, ATAC_IDR_BIOREP_LIST = _build_idr_experiment_lists(
    ATAC_IDR_ENABLED, "atac", "narrow", MULTI_BIOREP_EXPERIMENTS
)

# ATAC pseudoreplicate expansion lists
(
    ATAC_IDR_SPLIT_SOURCE_EXP,
    ATAC_IDR_SPLIT_SOURCE_NAME,
    ATAC_IDR_PR_PEAK_EXP,
    ATAC_IDR_PR_PEAK_SRC,
    ATAC_IDR_PR_PEAK_PR,
    ATAC_IDR_SELF_EXP,
    ATAC_IDR_SELF_BR,
) = _build_idr_pseudorep_lists(ATAC_IDR_EXPERIMENTS)


# Stage 64 CUT&Tag narrow IDR
CUTTAG_IDR_EXPERIMENTS, CUTTAG_IDR_BIOREP_EXP_LIST, CUTTAG_IDR_BIOREP_LIST = _build_idr_experiment_lists(
    CUTTAG_IDR_ENABLED, "cuttag", "narrow", MULTI_BIOREP_EXPERIMENTS
)
CUTTAG_IDR_EXPERIMENTS = sorted(CUTTAG_IDR_EXPERIMENTS)

(
    CUTTAG_IDR_SPLIT_SOURCE_EXP,
    CUTTAG_IDR_SPLIT_SOURCE_NAME,
    CUTTAG_IDR_PR_PEAK_EXP,
    CUTTAG_IDR_PR_PEAK_SRC,
    CUTTAG_IDR_PR_PEAK_PR,
    CUTTAG_IDR_SELF_EXP,
    CUTTAG_IDR_SELF_BR,
) = _build_idr_pseudorep_lists(CUTTAG_IDR_EXPERIMENTS)


# Stage 65 broad-peak IDR experiment lists (experimental opt-in)
BROAD_CHIPSEQ_IDR_EXPERIMENTS, BROAD_CHIPSEQ_IDR_BIOREP_EXP_LIST, BROAD_CHIPSEQ_IDR_BIOREP_LIST = (
    _build_idr_experiment_lists(
        BROAD_IDR_ENABLED and REPLICATE_ANALYSIS and BROAD_CHIPSEQ_IDR_ENABLED,
        "chipseq",
        "broad",
        MULTI_BIOREP_EXPERIMENTS,
    )
)
BROAD_CUTTAG_IDR_EXPERIMENTS, BROAD_CUTTAG_IDR_BIOREP_EXP_LIST, BROAD_CUTTAG_IDR_BIOREP_LIST = (
    _build_idr_experiment_lists(
        BROAD_IDR_ENABLED and REPLICATE_ANALYSIS and BROAD_CUTTAG_IDR_ENABLED,
        "cuttag",
        "broad",
        MULTI_BIOREP_EXPERIMENTS,
    )
)

BROAD_CHIPSEQ_IDR_EXPERIMENTS = sorted(BROAD_CHIPSEQ_IDR_EXPERIMENTS)
BROAD_CUTTAG_IDR_EXPERIMENTS = sorted(BROAD_CUTTAG_IDR_EXPERIMENTS)

BROAD_IDR_EXPERIMENTS = []
BROAD_IDR_EXPERIMENT_ASSAY = []
for exp, assay in sorted(
    [(e, "chipseq") for e in BROAD_CHIPSEQ_IDR_EXPERIMENTS] + [(e, "cuttag") for e in BROAD_CUTTAG_IDR_EXPERIMENTS]
):
    BROAD_IDR_EXPERIMENTS.append(exp)
    BROAD_IDR_EXPERIMENT_ASSAY.append(assay)

(
    BROAD_IDR_SPLIT_SOURCE_EXP,
    BROAD_IDR_SPLIT_SOURCE_NAME,
    BROAD_IDR_SPLIT_SOURCE_ASSAY,
    BROAD_IDR_PR_PEAK_EXP,
    BROAD_IDR_PR_PEAK_SRC,
    BROAD_IDR_PR_PEAK_ASSAY,
    BROAD_IDR_PR_PEAK_PR,
    BROAD_IDR_SELF_EXP,
    BROAD_IDR_SELF_BR,
    BROAD_IDR_SELF_ASSAY,
    BROAD_IDR_BIOREP_EXP_LIST,
    BROAD_IDR_BIOREP_LIST,
    BROAD_IDR_BIOREP_ASSAY,
) = _build_broad_idr_pseudorep_lists(BROAD_IDR_EXPERIMENTS, BROAD_IDR_EXPERIMENT_ASSAY)


# ---------------------------------------------------------------------------
# 3a5. Stage 62 consensus experiment lists
# ---------------------------------------------------------------------------

# Consensus modes: (assay, peak_mode) pairs with >=2 biological replicates
CONSENSUS_MODES = [
    ("chipseq", "narrow"),
    ("chipseq", "broad"),
    ("cuttag", "narrow"),
    ("cuttag", "broad"),
    ("atac", "narrow"),
]

if CONSENSUS_ENABLED and REPLICATE_ANALYSIS:
    CONSENSUS_EXPERIMENTS = {}  # (assay, peak_mode) -> [experiment_ids]
    CONSENSUS_BIOREP_EXP_LIST = []  # for expand(zip, ...)
    CONSENSUS_BIOREP_LIST = []

    for assay, peak_mode in CONSENSUS_MODES:
        exps = []
        for exp in MULTI_BIOREP_EXPERIMENTS:
            bioreps = _bioreps_for(exp, "treatment")
            if len(bioreps) >= 2:
                first = SAMPLE_MAP[TREATMENT_SAMPLES_BY_EXPERIMENT[exp][0]]
                if first["assay"] == assay and first["peak_mode"] == peak_mode:
                    exps.append(exp)
                    for br in sorted(bioreps):
                        CONSENSUS_BIOREP_EXP_LIST.append(exp)
                        CONSENSUS_BIOREP_LIST.append(br)
        CONSENSUS_EXPERIMENTS[(assay, peak_mode)] = sorted(exps)
else:
    CONSENSUS_EXPERIMENTS = {}
    CONSENSUS_BIOREP_EXP_LIST = []
    CONSENSUS_BIOREP_LIST = []


# ---------------------------------------------------------------------------
# 3a6. Stage 63 SEACR consensus experiment list
# ---------------------------------------------------------------------------

if CONSENSUS_ENABLED and REPLICATE_ANALYSIS and SEACR_ENABLED:
    SEACR_CONSENSUS_EXPERIMENTS = []
    for exp in MULTI_BIOREP_EXPERIMENTS:
        treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(exp, [])
        if not treatment_ids:
            continue
        bioreps = _bioreps_for(exp, "treatment")
        if len(bioreps) < 2:
            continue
        if all(SAMPLE_MAP[sid]["assay"] == "cuttag" and SAMPLE_MAP[sid]["layout"] == "PE" for sid in treatment_ids):
            SEACR_CONSENSUS_EXPERIMENTS.append(exp)
    SEACR_CONSENSUS_EXPERIMENTS = sorted(SEACR_CONSENSUS_EXPERIMENTS)
else:
    SEACR_CONSENSUS_EXPERIMENTS = []


# ---------------------------------------------------------------------------
# 3b. QC configuration and genome resource helpers
# ---------------------------------------------------------------------------

QC_CONFIG = VALIDATED_CONFIG.get(
    "qc",
    {
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
    },
)

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
        raw_config = VALIDATED_CONFIG.get("tool_parameters", {}).get("samtools_filter", {}).get("filter_flags", "")
        if isinstance(raw_config, str) and (raw_config.strip().startswith("0x") or raw_config.strip().startswith("0X")):
            return f"-F {raw_config.strip()}"
        return f"-F {raw}"
    return f"-F {raw}"


# Samples with a configured blacklist path (resource-gated blacklist filtering)
BLACKLIST_SAMPLES = [
    s for s in TREATMENT_SAMPLE_IDS if GENOME_RESOURCES.get(SAMPLE_MAP[s]["genome"], {}).get("blacklist", "")
]


def get_genome_resource(sample, key, default=""):
    """Return a genome resource path for a sample, or *default* if unset."""
    genome = SAMPLE_MAP[sample]["genome"]
    return GENOME_RESOURCES.get(genome, {}).get(key, default)


def has_genome_resource(sample, key):
    """Return True if a non-empty genome resource is configured for *sample*."""
    return bool(get_genome_resource(sample, key, ""))


def experiment_blacklist(experiment):
    """Return the experiment blacklist BED in the analysis_filtered branch."""
    if not ANALYSIS_FILTERED:
        return ""
    sample_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    if not sample_ids:
        return ""
    genome = SAMPLE_MAP[sample_ids[0]]["genome"]
    return GENOME_RESOURCES.get(genome, {}).get("blacklist", "")


TSS_SAMPLE_IDS = [
    sid for sid in PEAK_SAMPLE_IDS if QC_CONFIG.get("tss_enrichment", False) and has_genome_resource(sid, "gtf")
]
TSS_GENOMES = sorted({SAMPLE_MAP[sid]["genome"] for sid in TSS_SAMPLE_IDS})


# Stage 22: FE/ppois BigWig gating — requires signal_tracks + chrom_sizes
# has_genome_resource() checks non-empty genome_resources.<genome>.chrom_sizes.
# Path existence for non-empty chrom_sizes is already enforced by validate_samples.py.
SIGNAL_BW_SAMPLE_IDS = [
    sid for sid in PEAK_SAMPLE_IDS if QC_CONFIG.get("signal_tracks", True) and has_genome_resource(sid, "chrom_sizes")
]

if REPLICATE_ANALYSIS:
    SIGNAL_BW_EXPERIMENTS = [
        exp
        for exp in PEAK_MULTI_BIOREP_EXPERIMENTS
        if QC_CONFIG.get("signal_tracks", True)
        and SIGNAL_BW_SAMPLE_IDS
        and has_genome_resource(TREATMENT_SAMPLES_BY_EXPERIMENT[exp][0], "chrom_sizes")
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


# SCI-2 duplicate policy: an explicit duplicate_analysis=exclude/retain
# applies to every sample; "auto" resolves per sample by assay. CUT&Tag
# retains duplicates (apparent duplicates at high-affinity sites are mostly
# real independent fragments); all other assays exclude them. The
# narrow/broad peak mode no longer proxies this decision.
def get_duplicate_policy_for_assay(assay):
    """The auto duplicate policy for an assay: CUT&Tag retains, others exclude."""
    return "retain" if assay == "cuttag" else "exclude"


def get_duplicate_policy(sample_id):
    """Resolve the effective duplicate policy for a sample."""
    setting = VALIDATED_CONFIG.get("duplicate_analysis", "auto")
    if setting != "auto":
        return setting
    return get_duplicate_policy_for_assay(SAMPLE_MAP[sample_id]["assay"])


def get_remove_dup(wildcards):
    """Return the legacy BAM-layer switch for the resolved duplicate policy."""
    return "yes" if get_duplicate_policy(wildcards.sample) == "exclude" else "no"


def is_updated_profile():
    """True when the run uses the B/C updated analysis chain."""
    return VALIDATED_CONFIG.get("analysis_profile") == "updated"


FRAGMENT_CPM_MODE = VALIDATED_CONFIG.get("fragment_cpm", "off")


def is_fragment_cpm_enabled():
    """Return True when the updated profile requests fragment tracks."""
    return is_updated_profile() and FRAGMENT_CPM_MODE != "off"


def fragment_cpm_eligible(sample_id):
    """Return True when a treatment sample should get a fragment CPM track.

    ATAC, MNase, and known histone CUT&Tag targets are eligible in auto mode.
    Explicit ``on`` also allows CUT&Tag targets that the classifier cannot
    recognize; ChIP-seq remains read-level in PR-E1.
    """
    if not is_fragment_cpm_enabled():
        return False
    sample = SAMPLE_MAP[sample_id]
    assay = sample["assay"]
    if assay in ("atac", "mnase"):
        return True
    if assay == "cuttag":
        classifier = globals().get("classify_histone_target")
        if not callable(classifier):
            return FRAGMENT_CPM_MODE == "on"
        info = classifier(sample["target"], sample["peak_mode"])
        return FRAGMENT_CPM_MODE == "on" or info["inferred_histone_class"] != "unknown"
    return False


FRAGMENT_CPM_SAMPLE_IDS = [sid for sid in TREATMENT_SAMPLE_IDS if fragment_cpm_eligible(sid)]

if FRAGMENT_CPM_MODE == "on" and not FRAGMENT_CPM_SAMPLE_IDS:
    raise ValueError("fragment_cpm=on requested but no treatment sample is eligible for a " "fragment CPM track")

TSS_FRAGMENT_ENABLED = bool(VALIDATED_CONFIG.get("tss_fragment_profile", False))
_fragment_tss_candidates = [sid for sid in FRAGMENT_CPM_SAMPLE_IDS if SAMPLE_MAP[sid]["assay"] != "mnase"]
TSS_FRAGMENT_SAMPLE_IDS = [sid for sid in _fragment_tss_candidates if has_genome_resource(sid, "gtf")]
if TSS_FRAGMENT_ENABLED:
    _missing_gtf = [sid for sid in _fragment_tss_candidates if not has_genome_resource(sid, "gtf")]
    if not TSS_FRAGMENT_SAMPLE_IDS or _missing_gtf:
        raise ValueError(
            "tss_fragment_profile=true requires eligible ATAC/CUT&Tag "
            "fragment samples with genome_resources.<genome>.gtf; "
            f"missing GTF for: {', '.join(_missing_gtf) or 'none'}"
        )
TSS_FRAGMENT_GENOMES = sorted({SAMPLE_MAP[sid]["genome"] for sid in TSS_FRAGMENT_SAMPLE_IDS})


def analysis_bam(sample):
    """Return the analysis BAM for a sample.

    Legacy keeps the existing final.bam contract. Updated uses a derived
    analysis.bam so the raw marked BAM remains available.
    """
    if is_updated_profile():
        return f"{OUTDIR}/{sample}/02_align/{sample}.analysis.bam"
    return f"{OUTDIR}/{sample}/02_align/{sample}.final.bam"


def analysis_bai(sample):
    return f"{analysis_bam(sample)}.bai"


def get_duplicate_mode(wildcards):
    """Return the duplicate handling mode for the current profile."""
    if is_updated_profile():
        return "mark"
    return get_remove_dup(wildcards)


def _organelle_context_assay(sample_id):
    """Return the assay whose organelle policy applies to a sample.

    Controls inherit the policy of the treatment they are called for; a
    control row's own assay is not authoritative. A control shared across
    assays with different policies cannot be represented by one analysis BAM,
    so fail closed instead of silently choosing one policy.
    """
    sample = SAMPLE_MAP[sample_id]
    if sample.get("role") != "control":
        return sample["assay"]
    contexts = sorted(
        {SAMPLE_MAP[sid]["assay"] for sid in TREATMENT_SAMPLE_IDS if SAMPLE_MAP[sid].get("control_sample") == sample_id}
    )
    if not contexts:
        return sample["assay"]
    if len(contexts) > 1:
        raise ValueError(
            "updated profile cannot derive one organelle policy for control "
            f"{sample_id!r}; referenced by treatment assays: "
            f"{', '.join(contexts)}"
        )
    return contexts[0]


def get_organelle_mode(wildcards):
    """Resolve the effective organelle policy for a sample."""
    config_mode = str(VALIDATED_CONFIG.get("organelle_mode", "auto"))
    if config_mode != "auto":
        return config_mode
    assay = _organelle_context_assay(wildcards.sample)
    return "exclude" if assay in ("atac", "cuttag") else "report"


def get_organelle_chroms(wildcards):
    """Return organelle contig names for a sample's genome."""
    sample = SAMPLE_MAP[wildcards.sample]
    resources = VALIDATED_CONFIG.get("genome_resources", {})
    resource = resources.get(sample["genome"], {}) if isinstance(resources, dict) else {}
    chroms = resource.get("organelle_chroms") if isinstance(resource, dict) else None
    if isinstance(chroms, (list, tuple)) and chroms:
        return [str(chrom) for chrom in chroms]
    return ["chrM", "MT", "chrMT"]


# SCI-2: the BAM layer owns the duplicate decision, so every MACS3 call —
# per-sample, pooled, consensus, and IDR-ready — consumes the chosen
# molecule set faithfully with an explicit --keep-dup all instead of a
# second implicit coordinate-based collapse. idr_macs3_args must use this
# same constant.
MACS3_KEEP_DUP_ARGS = "--keep-dup all"


def get_macs3_args(wildcards):
    """Return MACS3 arguments string, dispatched by assay."""
    s = SAMPLE_MAP[wildcards.sample]
    if s["assay"] == "chipseq":
        assay_args = get_macs3_args_chipseq(wildcards)
    elif s["assay"] == "cuttag":
        assay_args = get_macs3_args_cuttag(wildcards)
    elif s["assay"] == "atac":
        assay_args = get_macs3_args_atac(wildcards)
    elif s["assay"] == "mnase":
        assay_args = get_macs3_args_mnase(wildcards)
    else:
        raise ValueError(f"Unsupported assay for MACS3 policy: {s['assay']}")
    return f"{assay_args} {MACS3_KEEP_DUP_ARGS}"


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


def _biorep_control_paths(experiment, bio_rep):
    """Return unique control BAMs referenced by one treatment biorep."""
    if not USE_CONTROL or not REPLICATE_ANALYSIS:
        return []
    sids = sorted(BIO_REP_GROUPS.get((experiment, "treatment", bio_rep), []))
    control_samples = set()
    control_bams = set()
    for sid in sids:
        sample = SAMPLE_MAP[sid]
        cs = sample.get("control_sample", "")
        cb = sample.get("control_bam", "")
        if cs:
            control_samples.add(cs)
        elif cb:
            if is_updated_profile():
                raise ValueError("updated profile does not yet support external control_bam")
            control_bams.add(cb)
    paths = []
    for cs in sorted(control_samples):
        paths.append(analysis_bam(cs) if is_updated_profile() else f"{OUTDIR}/{cs}/02_align/{cs}.final.bam")
    paths.extend(sorted(control_bams))
    return paths


def _biorep_control_bam(experiment, bio_rep):
    return f"{OUTDIR}/experiments/{experiment}/02_align/" f"biorep{bio_rep}.control.final.bam"


def _biorep_control_bai(experiment, bio_rep):
    return f"{_biorep_control_bam(experiment, bio_rep)}.bai"


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
        analysis_bam(sid),
        analysis_bai(sid),
    ]
    if USE_CONTROL:
        if s.get("control_sample"):
            cs = s["control_sample"]
            inputs.append(analysis_bam(cs))
        elif s.get("control_bam"):
            if is_updated_profile():
                raise ValueError("updated profile does not yet support external control_bam")
            inputs.append(s["control_bam"])
    return inputs
