"""Pure constants, allowed-value sets, and static defaults for config validation.

This module contains no control flow, no I/O, and no imports from
encode_pipeline.config.validator.
"""

import re

# Bowtie2 index file templates
BT2_STANDARD = (
    "{prefix}.1.bt2",
    "{prefix}.2.bt2",
    "{prefix}.3.bt2",
    "{prefix}.4.bt2",
    "{prefix}.rev.1.bt2",
    "{prefix}.rev.2.bt2",
)
BT2_LARGE = tuple(f.replace(".bt2", ".bt2l") for f in BT2_STANDARD)

# Sample identifier validation
SAMPLE_ID_PATTERN = r"^[A-Za-z0-9_.-]+$"
SANITIZE_PATTERN = r"[^A-Za-z0-9_.-]"
SAMPLE_ID_RE = re.compile(SAMPLE_ID_PATTERN)
SANITIZE_RE = re.compile(SANITIZE_PATTERN)

# Config keyword allowed values
REMOVE_DUP_KEYWORDS = ("auto", "yes", "no")
EXTEND_READS_KEYWORDS = ("auto", "yes", "no")

# Tool / mode allowed values
CUTTAG_SEACR_MODES = ("stringent", "relaxed")
IDR_RANKS = ("p.value", "signal.value")
BOWTIE2_MODES = ("", "very-fast", "fast", "sensitive", "very-sensitive")
BAMCOVERAGE_NORMALIZE_USING = ("CPM", "RPKM", "BPM", "None")

# Sample sheet allowed values
ASSAYS = ("chipseq", "cuttag", "atac", "mnase")
LAYOUTS = ("PE", "SE")
ROLES = ("treatment", "control")
PEAK_MODES = ("narrow", "broad", "nucleosome")

# Sample sheet required columns
SAMPLE_REQUIRED_COLUMNS = (
    "sample",
    "fastq_1",
    "layout",
    "assay",
    "target",
    "peak_mode",
    "genome",
    "bowtie2_index",
)

# Known config-block keys
CUTTAG_TOP_KEYS = frozenset({"peak_caller", "seacr"})
CUTTAG_SEACR_KEYS = frozenset({"enabled", "mode", "normalization", "threshold"})

MNASE_TOP_KEYS = frozenset({"mono_range", "fragments", "dyad_range", "callers"})
MNASE_FRAGMENTS_KEYS = frozenset({"sub", "mono", "di"})
MNASE_CALLERS = frozenset({"danpos3", "inps", "sem"})

IDR_KEYS = frozenset({"seed", "threshold", "rank"})

REPRODUCIBILITY_TOP_KEYS = frozenset({"enabled", "consensus", "idr"})
REPRODUCIBILITY_CONSENSUS_KEYS = frozenset(
    {"enabled", "min_replicates", "reciprocal_overlap"}
)
REPRODUCIBILITY_IDR_KEYS = frozenset(
    {
        "chipseq_narrow",
        "atac_narrow",
        "cuttag_narrow",
        "chipseq_broad_experimental",
        "cuttag_broad_experimental",
    }
)

TOOL_PARAMETERS_TOOLS = frozenset({
    "fastqc",
    "trim_galore",
    "bowtie2",
    "samtools_filter",
    "picard_markduplicates",
    "bamcoverage",
    "macs3",
    "multiqc",
    "idr_macs3",
})
TOOL_PARAMETERS_KEYS = {
    "fastqc": frozenset({"extra_args"}),
    "trim_galore": frozenset({"quality", "length", "stringency", "extra_args"}),
    "bowtie2": frozenset({"mode", "dovetail", "no_mixed", "no_discordant", "extra_args"}),
    "samtools_filter": frozenset({"filter_flags", "extra_args"}),
    "picard_markduplicates": frozenset({"optical_duplicate_pixel_distance", "extra_args"}),
    "bamcoverage": frozenset({"normalize_using", "smooth_length", "extra_args"}),
    "macs3": frozenset({"qvalue", "broad_cutoff", "extra_args"}),
    "multiqc": frozenset({"title", "extra_args"}),
    "idr_macs3": frozenset({"pvalue", "extra_args"}),
}

# Genome resource path fields validated for file existence
GENOME_RESOURCE_PATH_FIELDS = ("chrom_sizes", "blacklist", "gtf", "reference_fasta")

# MNase static default ranges. Validator converts these to lists when building
# the validated config so downstream consumers see the same output types as
# before.
MNASE_MONO_RANGE_DEFAULT = (140, 200)
MNASE_FRAGMENT_DEFAULTS = {
    "sub": (1, 139),
    "mono": (140, 200),
    "di": (300, 400),
}
MNASE_DYAD_RANGE_DEFAULT = (130, 200)
