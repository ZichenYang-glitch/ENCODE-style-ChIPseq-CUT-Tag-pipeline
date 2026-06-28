#!/usr/bin/env python3
"""Sample sheet and config validation for the ChIP-seq / CUT&Tag / ATAC-seq / MNase-seq pipeline.

Importable by workflow/Snakefile and runnable as a standalone CLI.

Usage:
    python3 scripts/validate_samples.py --config config/config.yaml
"""

import csv
import os
import sys

from encode_pipeline.config import cuttag as cuttag_validation
from encode_pipeline.config import defaults
from encode_pipeline.config import genome as genome_validation
from encode_pipeline.config import mnase as mnase_validation
from encode_pipeline.config import qc as qc_validation
from encode_pipeline.config import reproducibility as reproducibility_validation
from encode_pipeline.config import tools as tools_validation
from encode_pipeline.config.coercion import coerce_int


class ValidationError(Exception):
    """Raised when config or sample sheet validation fails."""


# Keep module-level aliases for backward compatibility with any code that may
# have imported these private names. These aliases are deprecated; prefer
# encode_pipeline.config.defaults directly.
_SAMPLE_ID_RE = defaults.SAMPLE_ID_RE
_SANITIZE_RE = defaults.SANITIZE_RE
_BT2_STANDARD = defaults.BT2_STANDARD
_BT2_LARGE = defaults.BT2_LARGE


# ---------------------------------------------------------------------------
# Stage 30: strict input validation helpers
# ---------------------------------------------------------------------------


def _check_fastq_exists(path: str, sample_id: str, label: str) -> None:
    """Raise ValidationError if *path* does not exist as a regular file."""
    if not os.path.isfile(path):
        raise ValidationError(
            f"Sample {sample_id!r}: {label} file not found: {path}"
        )


def _check_bowtie2_index(prefix: str, sample_id: str) -> None:
    """Raise ValidationError if neither the complete .bt2 nor .bt2l index set exists."""
    standard_set = [f.format(prefix=prefix) for f in _BT2_STANDARD]
    large_set = [f.format(prefix=prefix) for f in _BT2_LARGE]

    standard_ok = all(os.path.isfile(f) for f in standard_set)
    large_ok = all(os.path.isfile(f) for f in large_set)

    if standard_ok or large_ok:
        return

    missing = [f for f in standard_set if not os.path.isfile(f)]
    raise ValidationError(
        f"Sample {sample_id!r}: Bowtie2 index not found at {prefix!r}. "
        f"Missing {len(missing)} of {len(standard_set)} expected .bt2 files. "
        f"Neither complete .bt2 nor .bt2l set exists. "
        f"First missing: {missing[0] if missing else 'n/a'}"
    )


def _validate_strict_inputs(samples: list, strict_inputs: bool) -> None:
    """If strict_inputs, validate FASTQ and Bowtie2 index file existence."""
    if not strict_inputs:
        return
    for s in samples:
        sid = s.get("id", s.get("sample", "?"))
        fq1 = s.get("fq1", "")
        fq2 = s.get("fq2", "")
        bt2 = s.get("bt2_idx", "")
        layout = s.get("layout", "SE")

        if fq1:
            _check_fastq_exists(fq1, sid, "fastq_1")
        else:
            raise ValidationError(
                f"Sample {sid!r}: fastq_1 is empty"
            )
        if layout == "PE":
            if fq2:
                _check_fastq_exists(fq2, sid, "fastq_2")
            else:
                raise ValidationError(
                    f"Sample {sid!r}: PE layout requires fastq_2"
                )
        if bt2:
            _check_bowtie2_index(bt2, sid)
        else:
            raise ValidationError(
                f"Sample {sid!r}: bowtie2_index is empty"
            )


def _coerce_int(value, *, name: str, minimum: int) -> int:
    """Return a strictly parsed integer, rejecting bools and floats."""
    return coerce_int(
        value,
        name=name,
        minimum=minimum,
        error_cls=ValidationError,
    )


def _validate_effective_genome_size(genome: str, value) -> None:
    """Validate MACS3 effective genome size shortcut or positive integer."""
    return genome_validation.validate_effective_genome_size(
        genome,
        value,
        error_cls=ValidationError,
    )


# ---------------------------------------------------------------------------
# Config validation
# ---------------------------------------------------------------------------

def validate_config(config: dict) -> dict:
    """Validate and normalize workflow configuration.

    Returns a normalized dict with validated values.
    Raises ValidationError on invalid input.
    """
    validated: dict = {}

    # samples — required, path must exist
    samples_path = config.get("samples", "")
    if not samples_path:
        raise ValidationError(
            "config.yaml must set 'samples' to the TSV path."
        )
    if not os.path.isfile(samples_path):
        raise ValidationError(
            f"Sample sheet not found: {samples_path}"
        )
    validated["samples"] = samples_path

    # outdir — string, default "results"
    validated["outdir"] = str(config.get("outdir", "results"))

    # threads — positive integer
    threads = _coerce_int(
        config.get("threads", 8),
        name="threads",
        minimum=1,
    )
    validated["threads"] = threads

    # mapq — non-negative integer
    mapq = _coerce_int(
        config.get("mapq", 30),
        name="mapq",
        minimum=0,
    )
    validated["mapq"] = mapq

    # binsize — positive integer
    binsize = _coerce_int(
        config.get("binsize", 10),
        name="binsize",
        minimum=1,
    )
    validated["binsize"] = binsize

    # remove_dup — "auto", "yes", or "no"
    remove_dup = str(config.get("remove_dup", "auto"))
    if remove_dup not in defaults.REMOVE_DUP_KEYWORDS:
        raise ValidationError(
            f"config remove_dup must be auto, yes, or no, "
            f"got {remove_dup!r}"
        )
    validated["remove_dup"] = remove_dup

    # trim — boolean or string boolean, normalize to string
    trim_raw = config.get("trim", True)
    if isinstance(trim_raw, bool):
        validated["trim"] = str(trim_raw).lower()
    elif str(trim_raw).lower() in ("true", "false"):
        validated["trim"] = str(trim_raw).lower()
    else:
        raise ValidationError(
            f"config trim must be true or false, got {trim_raw!r}"
        )

    # extend_reads — "auto", "yes", "no", or positive integer string
    ext_raw = str(config.get("extend_reads", "auto"))
    if ext_raw not in defaults.EXTEND_READS_KEYWORDS and not (
        ext_raw.isdigit() and int(ext_raw) > 0
    ):
        raise ValidationError(
            f"config extend_reads must be auto, yes, no, or a "
            f"positive integer, got {ext_raw!r}"
        )
    validated["extend_reads"] = ext_raw

    # use_control — boolean or string boolean, normalize to bool
    use_raw = config.get("use_control", False)
    if isinstance(use_raw, bool):
        validated["use_control"] = use_raw
    else:
        _use = str(use_raw).lower()
        if _use in ("true", "yes", "1"):
            validated["use_control"] = True
        elif _use in ("false", "no", "0"):
            validated["use_control"] = False
        else:
            raise ValidationError(
                f"config use_control must be true or false, "
                f"got {use_raw!r}"
            )

    # multiqc — boolean or string boolean, normalize to bool
    mqc_raw = config.get("multiqc", True)
    if isinstance(mqc_raw, bool):
        validated["multiqc"] = mqc_raw
    elif str(mqc_raw).lower() in ("true", "false"):
        validated["multiqc"] = str(mqc_raw).lower() == "true"
    else:
        raise ValidationError(
            f"config multiqc must be true or false, got {mqc_raw!r}"
        )

    # genome_resources — optional, validate if present
    validated["genome_resources"] = _validate_genome_resources(
        config.get("genome_resources", {})
    )

    # qc — optional QC switches (legacy Stage 3 defaults true; heavier modules false)
    validated["qc"] = _validate_qc_config(config.get("qc", {}))

    # stage4b — optional Stage 4b replicate-aware outputs, default true
    stage4b_raw = config.get("stage4b", True)
    if isinstance(stage4b_raw, bool):
        validated["stage4b"] = stage4b_raw
    elif str(stage4b_raw).lower() in ("true", "false"):
        validated["stage4b"] = str(stage4b_raw).lower() == "true"
    else:
        raise ValidationError(
            f"config stage4b must be true or false, got {stage4b_raw!r}"
        )

    # tool_parameters — optional Stage 4c structured tool config, default empty
    validated["tool_parameters"] = _validate_tool_params(
        config.get("tool_parameters", {})
    )

    # stage5 — optional Stage 5 IDR, default false
    stage5_raw = config.get("stage5", False)
    if isinstance(stage5_raw, bool):
        validated["stage5"] = stage5_raw
    elif str(stage5_raw).lower() in ("true", "false"):
        validated["stage5"] = str(stage5_raw).lower() == "true"
    else:
        raise ValidationError(
            f"config stage5 must be true or false, got {stage5_raw!r}"
        )

    # stage5 requires stage4b
    if validated["stage5"] and not validated.get("stage4b", True):
        raise ValidationError(
            "config: stage5=true requires stage4b=true. "
            "Stage 5a depends on Stage 4b biorep BAMs and pooled control BAMs."
        )

    # reproducibility — Stage 53+ replicate-validated peak outputs
    validated["reproducibility"] = _validate_reproducibility(
        config.get("reproducibility", {}), validated
    )

    # Stage 55: determine whether ATAC IDR is enabled from validated reproducibility
    repro = validated["reproducibility"]
    atac_idr_enabled = (
        repro.get("enabled", False)
        and repro.get("idr", {}).get("atac_narrow", False)
    )

    # ATAC IDR requires stage4b (Stage 55)
    if atac_idr_enabled and not validated.get("stage4b", True):
        raise ValidationError(
            "config: reproducibility.idr.atac_narrow=true requires "
            "stage4b=true."
        )

    # Stage 64: determine whether CUT&Tag IDR is enabled
    cuttag_idr_enabled = (
        repro.get("enabled", False)
        and repro.get("idr", {}).get("cuttag_narrow", False)
    )

    # CUT&Tag IDR requires stage4b (Stage 64)
    if cuttag_idr_enabled and not validated.get("stage4b", True):
        raise ValidationError(
            "config: reproducibility.idr.cuttag_narrow=true requires "
            "stage4b=true."
        )

    # Stage 65: determine whether broad IDR is enabled
    broad_chipseq_idr_enabled = (
        repro.get("enabled", False)
        and repro.get("idr", {}).get("chipseq_broad_experimental", False)
    )
    broad_cuttag_idr_enabled = (
        repro.get("enabled", False)
        and repro.get("idr", {}).get("cuttag_broad_experimental", False)
    )

    # Broad IDR requires stage4b (Stage 65)
    if (broad_chipseq_idr_enabled or broad_cuttag_idr_enabled) \
       and not validated.get("stage4b", True):
        raise ValidationError(
            "config: reproducibility.idr.chipseq_broad_experimental=true "
            "or cuttag_broad_experimental=true requires stage4b=true."
        )

    # idr settings validated when stage5 or any IDR mode is enabled
    if (validated["stage5"] or atac_idr_enabled or cuttag_idr_enabled
        or broad_chipseq_idr_enabled or broad_cuttag_idr_enabled):
        validated["idr"] = _validate_idr_settings(config.get("idr", {}))
    else:
        validated["idr"] = {"threshold": 0.05, "rank": "p.value", "seed": 42}

    # cuttag — optional CUT&Tag-specific config (Stage 7b)
    validated["cuttag"] = _validate_cuttag_config(config.get("cuttag", {}))

    # mnase — optional MNase-seq config (Stage 39)
    validated["mnase"] = _validate_mnase_config(config.get("mnase", {}))

    return validated


def _validate_cuttag_config(cuttag: dict) -> dict:
    """Validate the cuttag config block. Returns normalized dict.

    Absent block → all defaults. Only validates keys for Stage 7b.
    """
    return cuttag_validation.validate_cuttag_config(
        cuttag,
        error_cls=ValidationError,
    )


def _validate_mnase_config(mnase: dict) -> dict:
    """Validate the mnase config block. Returns normalized dict.

    Absent block -> all defaults. Stage 39 keys + Stage 40 fragments,
    dyad_range, and callers.

    Fragment range precedence: fragments.mono > mono_range > [140, 200].
    """
    return mnase_validation.validate_mnase_config(
        mnase,
        error_cls=ValidationError,
    )


def _validate_qc_config(qc: dict) -> dict:
    """Validate and normalize the qc config block.

    Thin wrapper around encode_pipeline.config.qc.validate_qc_config.
    """
    return qc_validation.validate_qc_config(qc, error_cls=ValidationError)


def validate_picard_reference_resources(
    validated_config: dict, samples: list[dict]
):
    """Validate reference_fasta when qc.picard_metrics is enabled.

    Scans treatment sample genomes and checks that each has a non-empty
    reference_fasta in genome_resources. Raises ValidationError if
    qc.picard_metrics is true and any treatment genome is missing
    reference_fasta. Does nothing when qc.picard_metrics is false.
    """
    return genome_validation.validate_picard_reference_resources(
        validated_config,
        samples,
        error_cls=ValidationError,
    )


def validate_tss_annotation_resources(validated_config: dict, samples: list[dict]):
    """Validate GTF annotation when qc.tss_enrichment is enabled.

    Scans treatment sample genomes and checks that each has a non-empty gtf in
    genome_resources. Raises ValidationError if qc.tss_enrichment is true and
    any treatment genome is missing a GTF annotation. Does nothing when
    qc.tss_enrichment is false.
    """
    return genome_validation.validate_tss_annotation_resources(
        validated_config,
        samples,
        error_cls=ValidationError,
    )


def _validate_genome_resources(resources: dict) -> dict:
    """Validate genome_resources config block.

    Returns the resources dict unchanged if valid.
    Raises ValidationError on invalid entries.
    """
    return genome_validation.validate_genome_resources(
        resources,
        error_cls=ValidationError,
    )


def _validate_tool_params(tool_params) -> dict:
    """Validate and normalize the Stage 4c tool_parameters config block.

    Returns a normalized dict keyed by tool name.
    Raises ValidationError on unknown tools, unknown keys, or invalid types.
    Missing tool_parameters or empty blocks are valid and return an empty dict.
    """
    return tools_validation.validate_tool_params(
        tool_params,
        error_cls=ValidationError,
    )


def _validate_idr_settings(idr):
    """Validate the idr config block. Returns normalized dict.

    Only called when stage5 is true.
    """
    return reproducibility_validation.validate_idr_settings(
        idr,
        error_cls=ValidationError,
    )


def _validate_reproducibility(raw, validated_config):
    """Validate the reproducibility config block (Stage 53+).

    Returns a dict with validated reproducibility settings.
    When enabled is false/absent, returns {'enabled': False} without
    validating sub-keys.
    """
    return reproducibility_validation.validate_reproducibility(
        raw,
        validated_config,
        error_cls=ValidationError,
    )


def _sanitize_identifier(value: str) -> str:
    """Replace any character not in the safe identifier set with ``_``.

    This lets ``target`` values like ``Pol II`` or ``p300/CBP`` survive
    the default-for-``condition`` path without needing an explicit safe
    ``condition`` column.
    """
    return _SANITIZE_RE.sub("_", value)


def _parse_positive_int(value, *, default, name, sample, row):
    """Parse an optional positive integer with a default.

    Blanks and missing values use *default*.  Strictly rejects
    zero, negative, and non-integer values.
    """
    raw = (value or "").strip()
    if not raw:
        return default
    try:
        parsed = int(raw)
    except (ValueError, TypeError):
        raise ValidationError(
            f"Row {row}: sample {sample!r}: {name} must be a positive "
            f"integer, got {raw!r}"
        )
    if parsed <= 0:
        raise ValidationError(
            f"Row {row}: sample {sample!r}: {name} must be a positive "
            f"integer, got {parsed}"
        )
    return parsed


# ---------------------------------------------------------------------------
# Sample sheet validation
# ---------------------------------------------------------------------------

def load_and_validate_samples(
    sample_tsv: str,
    *,
    use_control: bool = False,
    stage5_enabled: bool = False,
    strict_inputs: bool = False,
    reproducibility_idr_atac_narrow: bool = False,
    reproducibility_idr_cuttag_narrow: bool = False,
    reproducibility_idr_chipseq_broad: bool = False,
    reproducibility_idr_cuttag_broad: bool = False,
) -> list[dict]:
    """Load and validate a sample sheet TSV.

    Returns a list of sample dicts. Raises ValidationError on invalid input.

    Sample dict keys:
        id, fq1, fq2, layout, assay, target, peak_mode, genome,
        bt2_idx, control_bam, role, control_sample,
        experiment, condition, replicate,
        biological_replicate, technical_replicate
    """
    if not os.path.isfile(sample_tsv):
        raise ValidationError(
            f"Sample sheet not found: {sample_tsv}"
        )

    samples: list[dict] = []
    seen_ids: set[str] = set()

    with open(sample_tsv) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []

        # Required columns
        required = list(defaults.SAMPLE_REQUIRED_COLUMNS)
        for col in required:
            if col not in fieldnames:
                raise ValidationError(
                    f"Sample sheet missing required column: {col!r}"
                )

        for i, row in enumerate(reader, start=2):
            sid = (row.get("sample") or "").strip()
            fq1 = (row.get("fastq_1") or "").strip()
            fq2 = (row.get("fastq_2") or "").strip()
            lo = (row.get("layout") or "").strip().upper()
            assay = (row.get("assay") or "").strip().lower()
            tgt = (row.get("target") or "").strip()
            pmode = (row.get("peak_mode") or "").strip().lower()
            gnm = (row.get("genome") or "").strip()
            bt2 = (row.get("bowtie2_index") or "").strip()

            # Optional columns with defaults
            ctrl_bam = (row.get("control_bam") or "").strip()
            role = (row.get("role") or "treatment").strip().lower()
            ctrl_sample = (row.get("control_sample") or "").strip()

            # Stage 4a optional columns with defaults
            experiment = _sanitize_identifier(
                (row.get("experiment") or "").strip() or sid
            )
            condition  = _sanitize_identifier(
                (row.get("condition") or "").strip() or tgt
            )
            replicate = _parse_positive_int(
                row.get("replicate"), default=1, name="replicate",
                sample=sid, row=i,
            )
            bio_rep = _parse_positive_int(
                row.get("biological_replicate"), default=replicate,
                name="biological_replicate", sample=sid, row=i,
            )
            tech_rep = _parse_positive_int(
                row.get("technical_replicate"), default=1,
                name="technical_replicate", sample=sid, row=i,
            )

            # --- Pass 1: per-row validation ---
            if not sid:
                raise ValidationError(
                    f"Row {i} in sample sheet has empty 'sample'"
                )
            if not _SAMPLE_ID_RE.match(sid):
                raise ValidationError(
                    f"Row {i}: sample ID {sid!r} contains invalid "
                    f"characters. Allowed: A-Z, a-z, 0-9, underscore, "
                    f"period, hyphen."
                )
            if sid in seen_ids:
                raise ValidationError(
                    f"Duplicate sample ID in sample sheet: {sid!r}"
                )
            seen_ids.add(sid)

            if not fq1:
                raise ValidationError(
                    f"Sample {sid!r} has empty 'fastq_1'"
                )
            if lo not in defaults.LAYOUTS:
                raise ValidationError(
                    f"Sample {sid!r}: layout must be PE or SE, "
                    f"got {lo!r}"
                )
            if lo == "PE" and not fq2:
                raise ValidationError(
                    f"Sample {sid!r}: PE layout requires 'fastq_2'"
                )
            if assay == "mnase" and lo != "PE":
                raise ValidationError(
                    f"Sample {sid!r}: assay=mnase requires paired-end "
                    f"layout (PE), got {lo!r}"
                )
            if assay not in defaults.ASSAYS:
                raise ValidationError(
                    f"Sample {sid!r}: assay must be chipseq, cuttag, "
                    f"atac, or mnase, got {assay!r}"
                )
            if not tgt:
                raise ValidationError(
                    f"Sample {sid!r} has empty 'target'"
                )
            if pmode not in defaults.PEAK_MODES:
                raise ValidationError(
                    f"Sample {sid!r}: peak_mode must be narrow, broad, "
                    f"or nucleosome, got {pmode!r}"
                )
            if assay == "atac" and pmode != "narrow":
                raise ValidationError(
                    f"Sample {sid!r}: assay=atac currently supports "
                    f"peak_mode=narrow only, got {pmode!r}"
                )
            if assay == "mnase" and pmode != "nucleosome":
                raise ValidationError(
                    f"Sample {sid!r}: assay=mnase requires "
                    f"peak_mode=nucleosome, got {pmode!r}"
                )
            if assay != "mnase" and pmode == "nucleosome":
                raise ValidationError(
                    f"Sample {sid!r}: peak_mode=nucleosome is only "
                    f"allowed for assay=mnase, got assay={assay!r}"
                )
            if not gnm:
                raise ValidationError(
                    f"Sample {sid!r} has empty 'genome'"
                )
            if not bt2:
                raise ValidationError(
                    f"Sample {sid!r} has empty 'bowtie2_index'"
                )
            if not _SAMPLE_ID_RE.match(experiment):
                raise ValidationError(
                    f"Row {i}: sample {sid!r}: experiment {experiment!r} "
                    f"contains invalid characters. "
                    f"Allowed: A-Z, a-z, 0-9, underscore, period, hyphen."
                )
            if not _SAMPLE_ID_RE.match(condition):
                raise ValidationError(
                    f"Row {i}: sample {sid!r}: condition {condition!r} "
                    f"contains invalid characters. "
                    f"Allowed: A-Z, a-z, 0-9, underscore, period, hyphen."
                )
            if not experiment:
                raise ValidationError(
                    f"Sample {sid!r} has empty 'experiment' after defaulting"
                )
            if not condition:
                raise ValidationError(
                    f"Sample {sid!r} has empty 'condition' after defaulting"
                )
            if role not in defaults.ROLES:
                raise ValidationError(
                    f"Sample {sid!r}: role must be treatment or control, "
                    f"got {role!r}"
                )
            if use_control:
                if ctrl_sample == sid:
                    raise ValidationError(
                        f"Sample {sid!r}: control_sample cannot reference "
                        f"itself"
                    )
                if ctrl_sample and ctrl_bam:
                    raise ValidationError(
                        f"Sample {sid!r}: cannot set both control_sample "
                        f"and control_bam"
                    )
                if ctrl_bam and not os.path.isfile(ctrl_bam):
                    raise ValidationError(
                        f"Sample {sid!r}: control_bam file not found: "
                        f"{ctrl_bam}"
                    )

            samples.append({
                "id": sid,
                "fq1": fq1,
                "fq2": fq2,
                "layout": lo,
                "assay": assay,
                "target": tgt,
                "peak_mode": pmode,
                "genome": gnm,
                "bt2_idx": bt2,
                "control_bam": ctrl_bam,
                "role": role,
                "control_sample": ctrl_sample,
                "experiment": experiment,
                "condition": condition,
                "replicate": replicate,
                "biological_replicate": bio_rep,
                "technical_replicate": tech_rep,
            })

    # --- Pass 2: cross-reference validation ---
    # Control references are intentionally ignored when use_control is false.
    if use_control:
        all_ids = {s["id"] for s in samples}
        for s in samples:
            cs = s.get("control_sample", "")
            if cs and cs not in all_ids:
                raise ValidationError(
                    f"Sample {s['id']!r}: control_sample {cs!r} "
                    f"not found in sample sheet"
                )
            if cs:
                cs_row = next(r for r in samples if r["id"] == cs)
                if cs_row["role"] != "control":
                    raise ValidationError(
                        f"Sample {s['id']!r}: control_sample {cs!r} "
                        f"has role={cs_row['role']!r}, expected role=control"
                    )

    # --- Pass 3: replicate group validation (Stage 4b) ---
    validate_replicate_groups(
        samples, use_control, stage5_enabled,
        reproducibility_idr_atac_narrow=reproducibility_idr_atac_narrow,
        reproducibility_idr_cuttag_narrow=reproducibility_idr_cuttag_narrow,
        reproducibility_idr_chipseq_broad=reproducibility_idr_chipseq_broad,
        reproducibility_idr_cuttag_broad=reproducibility_idr_cuttag_broad,
    )

    # --- Pass 4: strict input validation (Stage 30, optional) ---
    _validate_strict_inputs(samples, strict_inputs)

    return samples


def validate_replicate_groups(samples, use_control, stage5_enabled=False,
                              reproducibility_idr_atac_narrow=False,
                              reproducibility_idr_cuttag_narrow=False,
                              reproducibility_idr_chipseq_broad=False,
                              reproducibility_idr_cuttag_broad=False):
    """Pass 3: cross-replicate validation for Stage 4b replicate-aware outputs.

    Raises ValidationError on:
    - Inconsistent assay/target/genome/peak_mode/layout within an experiment
      (treatment samples only)
    - Duplicate (biological_replicate, technical_replicate) combos within
      an experiment
    - Mixed control types (control_bam vs control_sample) within an experiment
    - Partial controls: some treatment rows have controls and others do not
    """
    # Group treatment samples by experiment
    exp_treatments: dict[str, list[dict]] = {}
    exp_controls: dict[str, list[dict]] = {}
    for s in samples:
        if s["role"] == "treatment":
            exp_treatments.setdefault(s["experiment"], []).append(s)
        else:
            exp_controls.setdefault(s["experiment"], []).append(s)

    for exp, rows in exp_treatments.items():
        # --- Consistency checks ---
        if len(rows) < 2:
            continue

        first = rows[0]
        for field in ("assay", "target", "condition", "genome", "peak_mode", "layout"):
            values = {r[field] for r in rows}
            if len(values) > 1:
                raise ValidationError(
                    f"Experiment {exp!r}: treatment samples disagree on "
                    f"{field}: {sorted(values)}"
                )

        # Duplicate (biological_replicate, technical_replicate) combos
        seen_bt: set[tuple[int, int]] = set()
        for r in rows:
            bt = (r["biological_replicate"], r["technical_replicate"])
            if bt in seen_bt:
                raise ValidationError(
                    f"Experiment {exp!r}: duplicate "
                    f"(biological_replicate={bt[0]}, "
                    f"technical_replicate={bt[1]}) combination"
                )
            seen_bt.add(bt)

        # --- Control consistency (only when use_control is true) ---
        if not use_control:
            continue

        has_bam = [r for r in rows if r.get("control_bam")]
        has_sample = [r for r in rows if r.get("control_sample")]

        # Partial controls: not all have controls and not all lack them
        n_with_ctrl = len(has_bam) + len(has_sample)
        if n_with_ctrl not in (0, len(rows)):
            raise ValidationError(
                f"Experiment {exp!r}: partial controls detected — "
                f"{n_with_ctrl} of {len(rows)} treatment rows have "
                f"controls. All treatment rows in an experiment must "
                f"either all have controls or all lack controls."
            )

        # Mixed control types
        if has_bam and has_sample:
            raise ValidationError(
                f"Experiment {exp!r}: mixed control types — some rows "
                f"use control_bam and others use control_sample. "
                f"All treatment rows in an experiment must use the "
                f"same control type."
            )

        # Control sample references must belong to the same experiment
        if has_sample:
            for r in has_sample:
                cs = r["control_sample"]
                cs_row = next(s for s in samples if s["id"] == cs)
                if cs_row["experiment"] != exp:
                    raise ValidationError(
                        f"Experiment {exp!r}: sample {r['id']!r} references "
                        f"control_sample {cs!r} which belongs to experiment "
                        f"{cs_row['experiment']!r}, not {exp!r}"
                    )

    # --- Stage 5 IDR eligibility (only when stage5_enabled) ---
    if stage5_enabled:
        chipseq_narrow_exps = []
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]

            if first["assay"] != "chipseq":
                continue       # skip non-chipseq silently

            if first["peak_mode"] != "narrow":
                continue       # skip chipseq broad

            bio_reps = sorted({r["biological_replicate"] for r in rows})
            if len(bio_reps) != 2:
                raise ValidationError(
                    f"Stage 5 IDR: experiment {exp!r} has "
                    f"{len(bio_reps)} biological replicate(s) ({bio_reps}). "
                    f"Stage 5 requires exactly 2."
                )
            chipseq_narrow_exps.append(exp)

        if not chipseq_narrow_exps:
            raise ValidationError(
                "stage5=true but no eligible ChIP-seq narrow experiments "
                "were found."
            )

    # --- reproducibility.idr.atac_narrow eligibility (Stage 55) ---
    if reproducibility_idr_atac_narrow:
        atac_narrow_exps = []
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]
            if first["assay"] != "atac":
                continue       # skip non-ATAC silently
            if first["peak_mode"] != "narrow":
                continue       # ATAC broad — no IDR
            bio_reps = sorted({r["biological_replicate"] for r in rows})
            if len(bio_reps) != 2:
                raise ValidationError(
                    f"reproducibility.idr.atac_narrow: ATAC narrow "
                    f"experiment {exp!r} has {len(bio_reps)} biological "
                    f"replicate(s). IDR requires exactly 2."
                )
            atac_narrow_exps.append(exp)

        if not atac_narrow_exps:
            raise ValidationError(
                "reproducibility.idr.atac_narrow is true but no eligible "
                "ATAC narrow experiments with exactly 2 biological "
                "replicates were found."
            )

    # --- reproducibility.idr.cuttag_narrow eligibility (Stage 64) ---
    if reproducibility_idr_cuttag_narrow:
        cuttag_narrow_exps = []
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]

            if first["assay"] != "cuttag":
                continue       # skip non-CUT&Tag silently

            if first["peak_mode"] != "narrow":
                continue       # skip CUT&Tag broad

            bio_reps = sorted({r["biological_replicate"] for r in rows})
            if len(bio_reps) != 2:
                raise ValidationError(
                    f"reproducibility.idr.cuttag_narrow: CUT&Tag narrow "
                    f"experiment {exp!r} has {len(bio_reps)} biological "
                    f"replicate(s). IDR requires exactly 2."
                )
            cuttag_narrow_exps.append(exp)

        if not cuttag_narrow_exps:
            raise ValidationError(
                "reproducibility.idr.cuttag_narrow is true but no eligible "
                "CUT&Tag narrow experiments with exactly 2 biological "
                "replicates were found."
            )

    # --- reproducibility.idr.chipseq_broad_experimental eligibility (Stage 65) ---
    if reproducibility_idr_chipseq_broad:
        chipseq_broad_exps = []
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]

            if first["assay"] != "chipseq":
                continue       # skip non-chipseq silently

            if first["peak_mode"] != "broad":
                continue       # skip chipseq narrow

            bio_reps = sorted({r["biological_replicate"] for r in rows})
            if len(bio_reps) != 2:
                raise ValidationError(
                    f"reproducibility.idr.chipseq_broad_experimental: "
                    f"ChIP-seq broad experiment {exp!r} has "
                    f"{len(bio_reps)} biological replicate(s). "
                    f"IDR requires exactly 2."
                )
            chipseq_broad_exps.append(exp)

        if not chipseq_broad_exps:
            raise ValidationError(
                "reproducibility.idr.chipseq_broad_experimental is true "
                "but no eligible ChIP-seq broad experiments with exactly "
                "2 biological replicates were found."
            )

    # --- reproducibility.idr.cuttag_broad_experimental eligibility (Stage 65) ---
    if reproducibility_idr_cuttag_broad:
        cuttag_broad_exps = []
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]

            if first["assay"] != "cuttag":
                continue       # skip non-CUT&Tag silently

            if first["peak_mode"] != "broad":
                continue       # skip CUT&Tag narrow

            bio_reps = sorted({r["biological_replicate"] for r in rows})
            if len(bio_reps) != 2:
                raise ValidationError(
                    f"reproducibility.idr.cuttag_broad_experimental: "
                    f"CUT&Tag broad experiment {exp!r} has "
                    f"{len(bio_reps)} biological replicate(s). "
                    f"IDR requires exactly 2."
                )
            cuttag_broad_exps.append(exp)

        if not cuttag_broad_exps:
            raise ValidationError(
                "reproducibility.idr.cuttag_broad_experimental is true "
                "but no eligible CUT&Tag broad experiments with exactly "
                "2 biological replicates were found."
            )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _load_yaml(path: str) -> dict:
    """Load YAML config, preferring PyYAML if available, with stdlib fallback.

    PyYAML is available in the chipseq Conda environment (transitive
    dependency of Snakemake). The stdlib fallback handles bare-minimum
    YAML parsing for standalone use outside Conda.
    """
    try:
        import yaml
    except ImportError:
        return _parse_config_minimal(path)

    with open(path) as fh:
        return yaml.safe_load(fh)


def _parse_config_minimal(path: str) -> dict:
    """Minimal YAML config parser using stdlib only.

    Extracts flat keys plus the genome_resources and qc nested blocks.
    This covers the config structure used by validate_config().
    """
    config: dict = {}
    genome_resources: dict = {}
    qc: dict = {}
    tool_parameters: dict = {}
    cuttag: dict = {}
    seacr_sub: dict = {}
    section: str | None = None
    current_genome: str | None = None
    current_entry: dict | None = None
    current_tool: str | None = None
    current_tool_entry: dict | None = None

    def parse_scalar(value: str):
        value = value.strip().strip('"').strip("'")
        if value == "":
            return ""
        if value.lower() in ("true", "false"):
            return value.lower() == "true"
        if value.isdigit():
            return int(value)
        return value

    def save_current_genome() -> None:
        nonlocal current_genome, current_entry
        if current_genome is not None and current_entry is not None:
            genome_resources[current_genome] = current_entry
        current_genome = None
        current_entry = None

    def save_current_tool() -> None:
        nonlocal current_tool, current_tool_entry
        if current_tool is not None and current_tool_entry is not None:
            tool_parameters[current_tool] = current_tool_entry
        current_tool = None
        current_tool_entry = None

    with open(path) as fh:
        for line in fh:
            stripped = line.rstrip()

            # Skip comments and empty lines
            if not stripped or stripped.lstrip().startswith("#"):
                continue

            indent = len(line) - len(line.lstrip(" "))

            # Top-level keys start or end nested blocks.
            if indent == 0:
                save_current_genome()
                save_current_tool()

                if stripped.startswith("genome_resources:"):
                    section = "genome_resources"
                    continue
                if stripped.startswith("qc:"):
                    section = "qc"
                    continue
                if stripped.startswith("tool_parameters:"):
                    section = "tool_parameters"
                    continue
                if stripped.startswith("cuttag:"):
                    section = "cuttag"
                    continue

                section = None
                if ":" in stripped:
                    k, v = stripped.split(":", 1)
                    config[k.strip()] = parse_scalar(v)
                continue

            if section == "genome_resources" and indent == 2:
                # Genome key: "  hs:", "  mm10:", etc.
                key = stripped.strip().rstrip(":")
                if key and not key.startswith("#"):
                    save_current_genome()
                    current_genome = key
                    current_entry = {}
                continue

            if section == "genome_resources" and indent == 4:
                # Field within a genome entry
                field = stripped.strip()
                if ":" in field and current_genome is not None and current_entry is not None:
                    k, v = field.split(":", 1)
                    current_entry[k.strip()] = parse_scalar(v)
                continue

            if section == "qc" and indent == 2:
                field = stripped.strip()
                if ":" in field:
                    k, v = field.split(":", 1)
                    qc[k.strip()] = parse_scalar(v)
                continue

            if section == "cuttag" and indent == 2:
                key = stripped.strip().rstrip(":")
                if key == "seacr":
                    seacr_sub = {}
                elif ":" in stripped:
                    k, v = stripped.split(":", 1)
                    cuttag[k.strip()] = parse_scalar(v)
                continue

            if section == "cuttag" and indent == 4:
                field = stripped.strip()
                if ":" in field:
                    k, v = field.split(":", 1)
                    seacr_sub[k.strip()] = parse_scalar(v)
                continue

            if section == "tool_parameters" and indent == 2:
                # Tool block key: "  fastqc:", "  trim_galore:", etc.
                key = stripped.strip().rstrip(":")
                if key and not key.startswith("#"):
                    save_current_tool()
                    current_tool = key
                    current_tool_entry = {}
                continue

            if section == "tool_parameters" and indent == 4:
                # Field within a tool block
                field = stripped.strip()
                if ":" in field and current_tool is not None and current_tool_entry is not None:
                    k, v = field.split(":", 1)
                    current_tool_entry[k.strip()] = parse_scalar(v)
                continue

    # Save last genome entry
    save_current_genome()
    save_current_tool()

    if genome_resources:
        config["genome_resources"] = genome_resources
    if qc:
        config["qc"] = qc
    if tool_parameters:
        config["tool_parameters"] = tool_parameters
    if seacr_sub:
        cuttag["seacr"] = seacr_sub
    if cuttag:
        config["cuttag"] = cuttag

    return config


def main(argv=None):
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate ChIP-seq/CUT&Tag/ATAC-seq pipeline config and "
                    "sample sheet."
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to config YAML (e.g. config/config.yaml)"
    )
    parser.add_argument(
        "--strict-inputs", action="store_true", default=False,
        help="Validate FASTQ and Bowtie2 index file existence"
    )
    args = parser.parse_args(argv)

    config_path = args.config
    if not os.path.isfile(config_path):
        print(f"ERROR: Config file not found: {config_path}",
              file=sys.stderr)
        sys.exit(1)

    config = _load_yaml(config_path)

    try:
        validated = validate_config(config)
        repro = validated.get("reproducibility", {})
        atac_idr_enabled = (
            repro.get("enabled", False)
            and repro.get("idr", {}).get("atac_narrow", False)
        )
        cuttag_idr_enabled = (
            repro.get("enabled", False)
            and repro.get("idr", {}).get("cuttag_narrow", False)
        )
        broad_chipseq_idr_enabled = (
            repro.get("enabled", False)
            and repro.get("idr", {}).get(
                "chipseq_broad_experimental", False)
        )
        broad_cuttag_idr_enabled = (
            repro.get("enabled", False)
            and repro.get("idr", {}).get(
                "cuttag_broad_experimental", False)
        )
        samples = load_and_validate_samples(
            validated["samples"],
            use_control=validated["use_control"],
            stage5_enabled=validated.get("stage5", False),
            strict_inputs=args.strict_inputs,
            reproducibility_idr_atac_narrow=atac_idr_enabled,
            reproducibility_idr_cuttag_narrow=cuttag_idr_enabled,
            reproducibility_idr_chipseq_broad=broad_chipseq_idr_enabled,
            reproducibility_idr_cuttag_broad=broad_cuttag_idr_enabled,
        )
        validate_picard_reference_resources(validated, samples)
        validate_tss_annotation_resources(validated, samples)
    except ValidationError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    n_treatment = sum(1 for s in samples if s["role"] == "treatment")
    n_control = sum(1 for s in samples if s["role"] == "control")
    print(f"OK: {len(samples)} sample(s) validated "
          f"({n_treatment} treatment, {n_control} control)")
    print(f"     use_control: {validated['use_control']}")
    if validated.get("genome_resources"):
        genomes = ", ".join(validated["genome_resources"])
        print(f"     genome resources: {genomes}")


if __name__ == "__main__":
    main()
