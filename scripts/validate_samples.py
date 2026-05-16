#!/usr/bin/env python3
"""Sample sheet and config validation for the ChIP-seq/CUT&Tag Snakemake pipeline.

Importable by workflow/Snakefile and runnable as a standalone CLI.

Usage:
    python3 scripts/validate_samples.py --config config/config.yaml
"""

import csv
import os
import re
import sys

_SAMPLE_ID_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
_SANITIZE_RE = re.compile(r"[^A-Za-z0-9_.-]")


class ValidationError(Exception):
    """Raised when config or sample sheet validation fails."""


def _coerce_int(value, *, name: str, minimum: int) -> int:
    """Return a strictly parsed integer, rejecting bools and floats."""
    if isinstance(value, bool):
        raise ValidationError(
            f"config {name} must be an integer, got {value!r}"
        )
    if isinstance(value, int):
        parsed = value
    elif isinstance(value, str):
        text = value.strip()
        if not text.isdigit():
            raise ValidationError(
                f"config {name} must be an integer, got {value!r}"
            )
        parsed = int(text)
    else:
        raise ValidationError(
            f"config {name} must be an integer, got {value!r}"
        )

    if parsed < minimum:
        if minimum == 1:
            raise ValidationError(
                f"config {name} must be positive, got {parsed}"
            )
        raise ValidationError(
            f"config {name} must be non-negative, got {parsed}"
        )
    return parsed


def _validate_effective_genome_size(genome: str, value) -> None:
    """Validate MACS3 effective genome size shortcut or positive integer."""
    if isinstance(value, bool):
        valid = False
    elif isinstance(value, int):
        valid = value > 0
    elif isinstance(value, str):
        text = value.strip()
        valid = text in ("hs", "mm") or (text.isdigit() and int(text) > 0)
    else:
        valid = False

    if not valid:
        raise ValidationError(
            f"genome_resources.{genome}: effective_genome_size must be "
            f"'hs', 'mm', or a positive integer, got {value!r}"
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
    if remove_dup not in ("auto", "yes", "no"):
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
    if ext_raw not in ("auto", "yes", "no") and not (
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

    # qc — optional Stage3 QC switches, default all true
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

    return validated


def _validate_qc_config(qc: dict) -> dict:
    """Validate and normalize the qc config block.

    Missing keys default to True. Accepts boolean or string boolean.
    Returns a normalized dict with boolean values.
    """
    if not isinstance(qc, dict):
        raise ValidationError(
            f"qc must be a mapping, got {type(qc).__name__}"
        )

    def _normalize_bool(key: str) -> bool:
        raw = qc.get(key, True)
        if isinstance(raw, bool):
            return raw
        val = str(raw).lower()
        if val in ("true", "false"):
            return val == "true"
        raise ValidationError(
            f"qc.{key} must be true or false, got {raw!r}"
        )

    return {
        "blacklist_filter": _normalize_bool("blacklist_filter"),
        "frip": _normalize_bool("frip"),
        "library_complexity": _normalize_bool("library_complexity"),
        "nrf_pbc": _normalize_bool("nrf_pbc"),
        "signal_tracks": _normalize_bool("signal_tracks"),
        "summary": _normalize_bool("summary"),
    }


def _validate_genome_resources(resources: dict) -> dict:
    """Validate genome_resources config block.

    Returns the resources dict unchanged if valid.
    Raises ValidationError on invalid entries.
    """
    if not isinstance(resources, dict):
        raise ValidationError(
            f"genome_resources must be a mapping, "
            f"got {type(resources).__name__}"
        )

    for genome, entry in resources.items():
        if not isinstance(entry, dict):
            raise ValidationError(
                f"genome_resources.{genome} must be a mapping, "
                f"got {type(entry).__name__}"
            )

        egs = entry.get("effective_genome_size")
        if egs is None or egs == "":
            raise ValidationError(
                f"genome_resources.{genome}: "
                f"effective_genome_size is required"
            )

        _validate_effective_genome_size(genome, egs)

        # Optional path fields: if non-empty, validate existence
        for field in ("chrom_sizes", "blacklist", "gtf", "reference_fasta"):
            path = entry.get(field, "")
            if path and not os.path.isfile(path):
                raise ValidationError(
                    f"genome_resources.{genome}.{field}: "
                    f"file not found: {path}"
                )

    return resources


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
        required = [
            "sample", "fastq_1", "layout", "assay",
            "target", "peak_mode", "genome", "bowtie2_index",
        ]
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
            if lo not in ("PE", "SE"):
                raise ValidationError(
                    f"Sample {sid!r}: layout must be PE or SE, "
                    f"got {lo!r}"
                )
            if lo == "PE" and not fq2:
                raise ValidationError(
                    f"Sample {sid!r}: PE layout requires 'fastq_2'"
                )
            if assay not in ("chipseq", "cuttag"):
                raise ValidationError(
                    f"Sample {sid!r}: assay must be chipseq or cuttag, "
                    f"got {assay!r}"
                )
            if not tgt:
                raise ValidationError(
                    f"Sample {sid!r} has empty 'target'"
                )
            if pmode not in ("narrow", "broad"):
                raise ValidationError(
                    f"Sample {sid!r}: peak_mode must be narrow or broad, "
                    f"got {pmode!r}"
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
            if role not in ("treatment", "control"):
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
    validate_replicate_groups(samples, use_control)

    return samples


def validate_replicate_groups(samples, use_control):
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
    section: str | None = None
    current_genome: str | None = None
    current_entry: dict | None = None

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

                if stripped.startswith("genome_resources:"):
                    section = "genome_resources"
                    continue
                if stripped.startswith("qc:"):
                    section = "qc"
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

    # Save last genome entry
    save_current_genome()

    if genome_resources:
        config["genome_resources"] = genome_resources
    if qc:
        config["qc"] = qc

    return config


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Validate ChIP-seq/CUT&Tag pipeline config and "
                    "sample sheet."
    )
    parser.add_argument(
        "--config", required=True,
        help="Path to config YAML (e.g. config/config.yaml)"
    )
    args = parser.parse_args()

    config_path = args.config
    if not os.path.isfile(config_path):
        print(f"ERROR: Config file not found: {config_path}",
              file=sys.stderr)
        sys.exit(1)

    config = _load_yaml(config_path)

    try:
        validated = validate_config(config)
        samples = load_and_validate_samples(
            validated["samples"],
            use_control=validated["use_control"],
        )
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
