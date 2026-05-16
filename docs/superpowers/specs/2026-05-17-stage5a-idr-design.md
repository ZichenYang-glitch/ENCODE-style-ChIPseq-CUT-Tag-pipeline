# Stage 5a: TF ChIP-seq True Replicate IDR Foundation

## Scope

Stage 5a delivers the core IDR result: **true replicate IDR** between two biological
replicate IDR-ready MACS3 peak calls without pseudoreplicates.

**Stage 5a includes:**
- Config schema and validation for `stage5`, `idr`, and `tool_parameters.idr_macs3`
- Per-biorep IDR-ready MACS3 peak calling with relaxed p-value threshold
- True replicate IDR (`idr --samples biorep1_peaks biorep2_peaks`)
- Raw `idr.txt` and `idr.thresholded.narrowPeak` output

**Deferred to Stage 5b:**
- Pseudoreplicate BAM generation (deterministic hash split script)
- Self-pseudoreplicate IDR and pooled-pseudoreplicate IDR
- Final conservative/optimal peak set assembly
- Reproducibility summary TSV
- `idr.seed` (not used in 5a — only needed for pseudoreplicate splitting)

**Deferred to Stage 6/7:**
- Histone broad-peak IDR
- CUT&Tag IDR
- MultiQC integration of IDR results

## Configuration

### config.yaml

```yaml
# Stage 5: TF ChIP-seq IDR and reproducibility (default off)
stage5: false

# Stage 5 IDR analysis settings
idr:
  threshold: 0.05       # IDR threshold for thresholded peak set (float in (0,1))
  rank: "p.value"       # ranking measure: p.value or signal.value

# In tool_parameters block (Stage 4c):
tool_parameters:
  idr_macs3:
    pvalue: 0.1         # relaxed -p for IDR peak calls (float >0)
    extra_args: ""
```

`idr.seed` is deferred to Stage 5b (only needed for pseudoreplicate splitting).

### Config Keys

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `stage5` | boolean/string boolean | `false` | Master on/off. When false, no Stage 5 outputs. |
| `idr.threshold` | float in (0, 1) | 0.05 | IDR threshold for thresholded narrowPeak output |
| `idr.rank` | `"p.value"` or `"signal.value"` | `"p.value"` | Column to rank by in IDR |
| `tool_parameters.idr_macs3.pvalue` | float >0 | 0.1 | Relaxed p-value for IDR-ready MACS3 calls |
| `tool_parameters.idr_macs3.extra_args` | string | `""` | Additional MACS3 args for IDR calls |

### Stage 5 requires Stage 4b

Stage 5a depends on Stage 4b biorep BAMs (`merge_biorep_bam`) and pooled control
BAMs (`pool_control_bam`). Enabling Stage 5 while Stage 4b is disabled is a
configuration error.

**Validation rule:** If `stage5` is true and `stage4b` is false, raise:

```
"config: stage5=true requires stage4b=true. Stage 5a depends on
Stage 4b biorep BAMs and pooled control BAMs."
```

This check runs in `validate_config()` after both `stage4b` and `stage5` are
normalized.

## Validation Semantics

### When `stage5: false` (or absent)

- **No Stage 5 eligibility validation.** The `load_and_validate_samples()` and
  `validate_replicate_groups()` functions do not check chipseq-only, narrow-only,
  or exact-2-biorep constraints. Any sample sheet is valid regardless of assay,
  peak_mode, or replicate count.
- **No Stage 5 outputs or targets.** `IDR_EXPERIMENTS` is empty, all IDR expand()
  calls return `[]`.
- **Config block validation is lightweight but runs.** If an `idr` block is
  present in config.yaml when `stage5: false`, `_validate_idr_settings()` is NOT
  called. Only `stage5` itself is validated (boolean/string boolean). The
  `tool_parameters.idr_macs3` block is validated by `_validate_tool_params`
  because `KNOWN_TOOLS` includes it; this is harmless — the block is accepted
  but unused when `stage5: false`.
- **`idr` block absent when `stage5: false`** → no validation, no error. Defaults
  are silently unused.

### When `stage5: true`

- **`stage4b: true` must also be set.** Validation error otherwise.
- **`idr` block validated.** `_validate_idr_settings(config.get("idr", {}))` is
  called. Missing `idr` block → defaults used. Invalid threshold/rank → error.
- **Stage 5 eligibility validated.** `validate_replicate_groups()` checks:
  chipseq assay only, narrow peak_mode only, exactly 2 treatment bio-reps per
  experiment.
- **`tool_parameters.idr_macs3` validated.** By existing `_validate_tool_params()`,
  same as any other tool block.

### Implementation wiring

```python
# In validate_config():
validated["stage5"] = _normalize_stage5_bool(config.get("stage5", False))

# Stage 5 requires Stage 4b
if validated["stage5"] and not validated["stage4b"]:
    raise ValidationError(
        "config: stage5=true requires stage4b=true. "
        "Stage 5a depends on Stage 4b biorep BAMs and pooled control BAMs."
    )

# Only validate idr settings block when stage5 is true
if validated["stage5"]:
    validated["idr"] = _validate_idr_settings(config.get("idr", {}))
else:
    validated["idr"] = {"threshold": 0.05, "rank": "p.value"}  # defaults, unused
```

```python
# In load_and_validate_samples():
stage5_enabled = validated_config.get("stage5", False)  # passed from caller
```

```python
# In validate_replicate_groups(samples, use_control, stage5_enabled=False):
# Stage 5 eligibility block runs only when stage5_enabled is True
```

## Validation — `validate_samples.py`

### `_validate_idr_settings`

```python
def _validate_idr_settings(idr: dict) -> dict:
    """Validate the idr config block. Returns normalized dict.

    Only called when stage5 is true.
    """
    if not isinstance(idr, dict):
        raise ValidationError(
            "idr must be a mapping, got {}".format(type(idr).__name__)
        )

    threshold = idr.get("threshold", 0.05)
    if isinstance(threshold, bool):
        raise ValidationError(
            "idr.threshold must be a float in (0, 1), got {!r}".format(threshold)
        )
    try:
        threshold = float(threshold)
    except (ValueError, TypeError):
        raise ValidationError(
            "idr.threshold must be a float in (0, 1), got {!r}".format(threshold)
        )
    if not (0 < threshold < 1):
        raise ValidationError(
            "idr.threshold must be in (0, 1), got {}".format(threshold)
        )

    rank = str(idr.get("rank", "p.value"))
    if rank not in ("p.value", "signal.value"):
        raise ValidationError(
            "idr.rank must be 'p.value' or 'signal.value', got {!r}".format(rank)
        )

    return {"threshold": threshold, "rank": rank}
```

### Stage 5 eligibility (in `validate_replicate_groups`)

Added after existing Pass 3 consistency checks, guarded by `stage5_enabled`:

```python
    # --- Stage 5 IDR eligibility (only when stage5_enabled) ---
    if stage5_enabled:
        for exp, rows in exp_treatments.items():
            if len(rows) == 0:
                continue
            first = rows[0]

            if first["assay"] != "chipseq":
                raise ValidationError(
                    "Stage 5 IDR: experiment {!r} has assay {!r}. "
                    "Stage 5 supports chipseq only.".format(exp, first["assay"])
                )

            if first["peak_mode"] != "narrow":
                raise ValidationError(
                    "Stage 5 IDR: experiment {!r} has peak_mode {!r}. "
                    "Stage 5 supports narrowPeak only.".format(exp, first["peak_mode"])
                )

            bio_reps = sorted({r["biological_replicate"] for r in rows})
            if len(bio_reps) != 2:
                raise ValidationError(
                    "Stage 5 IDR: experiment {!r} has {} "
                    "biological replicate(s) ({}). "
                    "Stage 5 requires exactly 2.".format(
                        exp, len(bio_reps), bio_reps)
                )
```

### `idr_macs3` in `KNOWN_TOOLS`

In `_validate_tool_params`, add:
```python
KNOWN_TOOLS = {..., "idr_macs3"}
KNOWN_KEYS["idr_macs3"] = {"pvalue", "extra_args"}
```

`pvalue` reuses `_normalize_positive_float`. Validated regardless of `stage5`
setting (same as all other tool blocks — harmless when unused).

## Output Namespace

```
results/experiments/<experiment>/
├── 04_peaks/
│   └── idr/                                       # Stage 5a NEW
│       ├── <experiment>_biorep<bio_rep1>_idr_peaks.narrowPeak
│       └── <experiment>_biorep<bio_rep2>_idr_peaks.narrowPeak
│       (bio-rep labels are the actual positive integers from the sample sheet)
├── 06_idr/                                        # Stage 5a NEW
│   └── true_replicates/
│       ├── idr.txt                                # raw IDR output
│       └── idr.thresholded.narrowPeak             # thresholded peaks
└── logs/
    ├── <experiment>.idr.log                       # raw IDR log
    └── <experiment>.idr.thresholded.log           # thresholded IDR log
```

## Snakefile Changes

### Config extraction

```python
STAGE5 = VALIDATED_CONFIG.get("stage5", False)
IDR_SETTINGS = VALIDATED_CONFIG.get("idr", {})
IDR_THRESHOLD = float(IDR_SETTINGS.get("threshold", 0.05))
IDR_RANK = str(IDR_SETTINGS.get("rank", "p.value"))
```

### Sample loading (updated signature)

```python
SAMPLES = load_and_validate_samples(
    VALIDATED_CONFIG["samples"],
    use_control=USE_CONTROL,
    stage5_enabled=STAGE5,
)
```

### Derived structures

```python
if STAGE5:
    # Pre-computed (experiment, bio_rep) pairs for IDR expansion.
    # Bio-rep labels are derived from the actual sample sheet data,
    # never hardcoded as 1 and 2.
    IDR_EXPERIMENTS = []
    IDR_BIOREP_EXP_LIST = []   # zip-compatible experiment column
    IDR_BIOREP_LIST = []       # zip-compatible bio_rep column

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
```

### Rule all targets

```python
# Stage 5a: IDR-ready per-biorep peaks
idr_biorep_peaks = (
    expand(
        "{outdir}/experiments/{experiment}/04_peaks/idr/{experiment}_biorep{bio_rep}_idr_peaks.narrowPeak",
        zip,
        outdir=OUTDIR,
        experiment=IDR_BIOREP_EXP_LIST, bio_rep=IDR_BIOREP_LIST,
    )
    if STAGE5 and IDR_EXPERIMENTS else []
),

# Stage 5a: true replicate IDR (raw + thresholded)
true_idr_raw = (
    expand(
        "{outdir}/experiments/{experiment}/06_idr/true_replicates/idr.txt",
        outdir=OUTDIR,
        experiment=IDR_EXPERIMENTS,
    )
    if STAGE5 and IDR_EXPERIMENTS else []
),
true_idr_thresh = (
    expand(
        "{outdir}/experiments/{experiment}/06_idr/true_replicates/idr.thresholded.narrowPeak",
        outdir=OUTDIR,
        experiment=IDR_EXPERIMENTS,
    )
    if STAGE5 and IDR_EXPERIMENTS else []
),
```

## New Rules

### `macs3_idr_biorep` (in `workflow/rules/idr.smk` — NEW FILE)

IDR-ready MACS3 call on a single biorep BAM.

```python
rule macs3_idr_biorep:
    output:
        f"{OUTDIR}/experiments/{{experiment}}/04_peaks/idr/{{experiment}}_biorep{{bio_rep}}_idr_peaks.narrowPeak",
    input:
        lambda wc: _idr_biorep_peaks_inputs(wc),
    params:
        macs3_args = lambda wc: _idr_macs3_args(wc),
        sample     = lambda wc: f"{wc.experiment}_biorep{wc.bio_rep}_idr",
    wildcard_constraints:
        bio_rep = r"\d+",
    log:
        f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}_biorep{{bio_rep}}.idr.macs3.log",
    threads: THREADS,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output})"

        set -- {input:q}
        TREATMENT="$1"
        # Input order: treatment.bam ($1)[, treatment.bam.bai ($2)][, pooled.control.bam ($3)]
        if [[ $# -ge 3 ]]; then
            macs3 callpeak \
                -t "$TREATMENT" \
                -c "$3" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} \
                2>&1 | tee {log:q}
        else
            macs3 callpeak \
                -t "$TREATMENT" \
                -n {params.sample:q} \
                --outdir "$(dirname {output:q})" \
                {params.macs3_args} \
                2>&1 | tee {log:q}
        fi

        EXPECTED="{output}"
        if [[ ! -f "$EXPECTED" ]]; then
            echo "ERROR: Expected peak file not found: $EXPECTED" >&2
            exit 1
        fi
        """
```

**`_idr_biorep_peaks_inputs(wildcards)`:**

```python
def _idr_biorep_peaks_inputs(wildcards):
    """Return inputs for MACS3 IDR peak call on a single biorep BAM.

    Input order: [biorep_bam, biorep_bam.bai, ...optional_pooled_control_bam]
    The pooled control BAM is included as an explicit dependency so Snakemake
    schedules pool_control_bam before this rule.
    """
    exp = wildcards.experiment
    br = int(wildcards.bio_rep)
    inputs = [
        f"{OUTDIR}/experiments/{exp}/02_align/biorep{br}.final.bam",
        f"{OUTDIR}/experiments/{exp}/02_align/biorep{br}.final.bam.bai",
    ]
    if exp in POOLED_CONTROL_EXPERIMENTS:
        inputs.append(
            f"{OUTDIR}/experiments/{exp}/02_align/{exp}.pooled.control.final.bam"
        )
    return inputs
```

**`_idr_macs3_args(wildcards)`:**

```python
def _idr_macs3_args(wildcards):
    """Return MACS3 args for IDR-ready peak calls on a biorep BAM.

    Uses the same layout/genome as per-sample MACS3, but replaces
    -q (q-value) with -p (p-value from idr_macs3 config). Never emits
    both -q and -p.
    """
    experiment = wildcards.experiment
    treatment_ids = TREATMENT_SAMPLES_BY_EXPERIMENT.get(experiment, [])
    s = SAMPLE_MAP[treatment_ids[0]]
    fmt = "BAMPE" if s["layout"] == "PE" else "BAM"
    genome = _normalize_genome(s["genome"])
    pvalue = _tool_param("idr_macs3", "pvalue", 0.1)
    extra = _tool_param("idr_macs3", "extra_args", "")
    return f"-f {fmt} -g {genome} -p {pvalue} {extra}".strip()
```

### `idr_true_replicates` (in `workflow/rules/idr.smk`)

Runs the `idr` tool twice: once for the raw IDR table, once for the thresholded
narrowPeak output. Each invocation writes to a separate `--log-output-file`.

```python
rule idr_true_replicates:
    output:
        idr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/true_replicates/idr.txt",
        thr_out = f"{OUTDIR}/experiments/{{experiment}}/06_idr/true_replicates/idr.thresholded.narrowPeak",
    input:
        peaks1 = lambda wc: _idr_peak_input(wc.experiment, 0),
        peaks2 = lambda wc: _idr_peak_input(wc.experiment, 1),
    params:
        threshold     = IDR_THRESHOLD,
        rank          = IDR_RANK,
        output_prefix = lambda wc: (
            f"{OUTDIR}/experiments/{wc.experiment}/06_idr/true_replicates/idr"
        ),
    log:
        raw_log  = f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.idr.log",
        thr_log  = f"{OUTDIR}/experiments/{{experiment}}/logs/{{experiment}}.idr.thresholded.log",
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        command -v idr >/dev/null 2>&1 || {
            echo "ERROR: idr is required for Stage 5 IDR analysis but was not found in PATH." >&2
            exit 1
        }

        mkdir -p "$(dirname {output.idr_out})"

        # Run 1: raw IDR output (no --idr-threshold).
        # Log goes to {log.raw_log}.
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --output-file {output.idr_out:q} \
            --log-output-file {log.raw_log:q}

        # Run 2: thresholded narrowPeak output.
        # Log goes to {log.thr_log}.
        # --idr-threshold N keeps peaks whose accumulated global IDR
        # score falls at or below N (more stringent = smaller value).
        idr --samples {input.peaks1:q} {input.peaks2:q} \
            --input-file-type narrowPeak \
            --rank {params.rank} \
            --idr-threshold {params.threshold} \
            --output-file {output.thr_out:q} \
            --log-output-file {log.thr_log:q}
        """
```

**`_idr_peak_input(experiment, index)` — index is 0-based into sorted bio_rep list:**

```python
def _idr_peak_input(experiment, index):
    """Return the IDR peak file for a bio_rep by 0-based index.

    index 0 → first bio_rep (smallest label), index 1 → second.
    The actual bio_rep numbers are derived from the sample sheet.
    """
    bioreps = _bioreps_for(experiment, "treatment")
    br = bioreps[index]
    return (
        f"{OUTDIR}/experiments/{experiment}/04_peaks/idr/"
        f"{experiment}_biorep{br}_idr_peaks.narrowPeak"
    )
```

Call sites: `_idr_peak_input(wc.experiment, 0)` and `_idr_peak_input(wc.experiment, 1)`.

## IDR Output Policy

**Two invocations, two output files, two log files:**

| Invocation | `--idr-threshold` | `--output-file` | `--log-output-file` | Produces |
|---|---|---|---|---|
| Run 1 | absent | `.../idr.txt` | `.../<exp>.idr.log` | Raw IDR table (all peaks with IDR scores) |
| Run 2 | `0.05` | `.../idr.thresholded.narrowPeak` | `.../<exp>.idr.thresholded.log` | Peaks passing IDR threshold |

**IDR threshold semantics:**

`idr --idr-threshold N` keeps peaks whose global IDR value is **at or below** N.
A threshold of 0.05 means only peaks with global IDR ≤ 0.05 are retained — this
is a *more stringent* filter than 0.10. Smaller threshold = more conservative
peak set.

IDR output table columns (from nboley/idr docs):

| Col | Field | Notes |
|-----|-------|-------|
| 1 | chrom | |
| 2 | start | |
| 3 | end | |
| 4 | name | |
| 5 | score | `-125 * log2(IDR)` |
| 6 | strand | |
| 7 | signalValue | |
| 8 | p-value | |
| 9 | q-value | |
| 10 | summit | |
| 11 | localIDR | |
| 12 | globalIDR | Used by `--idr-threshold` |

The `--idr-threshold` flag filters on the globalIDR column. We use the tool's
native filtering; no manual awk/post-processing is needed.

## Files Changed

| File | Change |
|------|--------|
| `config/config.yaml` | Add `stage5: false` and commented `idr` block |
| `scripts/validate_samples.py` | Add `stage5` boolean, `_validate_idr_settings()` (gated on `stage5: true`), `idr_macs3` in `KNOWN_TOOLS`, `stage5=true` requires `stage4b=true` check, extend `load_and_validate_samples` + `validate_replicate_groups` with `stage5_enabled`, Stage 5 eligibility checks (gated on `stage5_enabled`) |
| `workflow/Snakefile` | Extract `STAGE5`, `IDR_SETTINGS`, `IDR_EXPERIMENTS`, `IDR_BIOREP_EXP_LIST`, `IDR_BIOREP_LIST`; pass `stage5_enabled` to validation; add IDR targets to `rule all` |
| `workflow/rules/idr.smk` | **New** — `macs3_idr_biorep`, `_idr_biorep_peaks_inputs`, `_idr_macs3_args`, `idr_true_replicates`, `_idr_peak_input` |
| `workflow/envs/chipseq.yml` | Add `idr` Conda dependency |
| `workflow/schemas/config.schema.yaml` | Document `stage5`, `idr`, `idr_macs3` |
| `README.md` | Add Stage 5a section |
| `KNOWN_ISSUES.md` | Update Stage 5 items |

## Test Plan (`test/test_stage5a_stress.py`)

Reuses Stage 4b/4c harness: `SNAKEMAKE` env var with conda fallback, temp file
cleanup, `scripts/__pycache__` removal.

**Validation tests:**
1. No `stage5` key → validates, same DAG
2. `stage5: false` → validates, same DAG
3. `stage5: true` + non-chipseq assay → rejected
4. `stage5: true` + broad peak_mode → rejected
5. `stage5: true` + !=2 bio-reps (1 replicate) → rejected
6. `stage5: true` + !=2 bio-reps (3 replicates) → rejected
7. Invalid `idr.threshold` (0, 1.5, "high") → rejected
8. Invalid `idr.rank` ("score") → rejected
9. `stage5: true` + `stage4b: false` → rejected with clear message about dependency

**DAG tests:**
10. `stage5: true` + 2 chipseq narrow bio-reps → validates, IDR rules in DAG
11. `stage5: true` + arbitrary bio_rep labels (2 and 4) → correct IDR paths
    in DAG (e.g., `_biorep2_idr_peaks.narrowPeak` and `_biorep4_idr_peaks.narrowPeak`)
12. `idr_macs3.pvalue` in dry-run `-n -p` → `-p 0.1` visible
13. `idr --idr-threshold` in dry-run `-n -p` → threshold flag visible

## Verification

```bash
python3 scripts/validate_samples.py --config config/config.yaml
conda run -n chipseq snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet
python3 test/test_validation_stress.py
python3 test/test_stage4b_stress.py
python3 test/test_stage4c_stress.py
python3 test/test_stage5a_stress.py
```
