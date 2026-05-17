# Stage 7b: CUT&Tag Peak Caller Abstraction and Optional SEACR Support

## Scope

Stage 7b introduces a CUT&Tag peak caller abstraction layer and optional SEACR
sidecar outputs. MACS3 remains the canonical peak caller for all downstream
analyses. SEACR produces independent sidecar BED files that do not interact
with FRiP, peak_counts, blacklist filtering, pooled peaks, or signal tracks.

**Stage 7b includes:**
- `cuttag` top-level config block with `peak_caller` and `seacr` sub-blocks
- `cuttag.seacr.enabled: false` (default off) — sidecar SEACR bedGraph + BED
  production for CUT&Tag PE samples
- A SEACR input bedGraph rule using `bedtools genomecov -bg`
- A SEACR invocation rule
- Validation for new config keys

**Deferred to Stage 7c+:**
- Switching canonical peak caller from MACS3 to SEACR/GoPeaks
- GoPeaks support
- Downstream QC integration with SEACR outputs
- Spike-in normalization

## Configuration

### config.yaml

```yaml
cuttag:
  peak_caller: "macs3"      # canonical peak caller (future: seacr, gopeaks)
  seacr:
    enabled: false           # sidecar SEACR outputs (default off)
    mode: "stringent"        # stringent | relaxed
    normalization: "non"     # non only in 7b (no control bedGraph)
    threshold: 0.01         # numeric FDR threshold for treatment-only mode
```

### Config Keys

| Key | Type | Default | Description |
|-----|------|---------|-------------|
| `cuttag.peak_caller` | `"macs3"` | `"macs3"` | Canonical peak caller (future extension). Not honored in 7b. |
| `cuttag.seacr.enabled` | boolean | `false` | Enable sidecar SEACR outputs |
| `cuttag.seacr.mode` | `"stringent"` or `"relaxed"` | `"stringent"` | SEACR threshold mode |
| `cuttag.seacr.normalization` | `"non"` only in 7b | `"non"` | SEACR normalization (no control bedGraph) |
| `cuttag.seacr.threshold` | float in (0, 1) | `0.01` | Numeric FDR threshold for treatment-only mode |

`cuttag.peak_caller` is a forward-compatibility placeholder — Stage 7b always
uses MACS3 for canonical peaks. SEACR is only a sidecar when
`cuttag.seacr.enabled: true`.

## Validation

New function `_validate_cuttag_config()` in `validate_samples.py`:

- `cuttag` absent → defaults used (all keys default)
- `cuttag.peak_caller` must be `"macs3"` (only accepted value in 7b)
- `cuttag.seacr.enabled` must be boolean or string boolean
- `cuttag.seacr.mode` must be `"stringent"` or `"relaxed"`
- `cuttag.seacr.normalization` must be `"non"` (only accepted value in 7b)
- `cuttag.seacr.threshold` must be float in (0, 1), default 0.01
- Unknown keys rejected
- Unknown keys rejected

### Minimal YAML parser

`_parse_config_minimal()` must parse the new `cuttag:` top-level block.
Parser behavior for:

```yaml
cuttag:
  peak_caller: "macs3"
  seacr:
    enabled: false
```

- `indent=0`: `cuttag:` starts the section
- `indent=2`: keys like `peak_caller` are scalar; `seacr:` is detected as
  a subsection (starts a nested dict, same pattern as genome names under
  `genome_resources:`)
- `indent=4`: `seacr` subkeys (`enabled`, `mode`, `normalization`, `threshold`)
  are parsed as scalars

Called from `validate_config()` after other block validations.

### config.schema.yaml

Add `cuttag` section with `peak_caller`, `seacr`, and sub-keys documented.

## Derived Structures

### Snakefile additions

```python
CUTTAG_CONFIG = VALIDATED_CONFIG.get("cuttag", {})
SEACR_ENABLED = (
    CUTTAG_CONFIG.get("seacr", {}).get("enabled", False)
    if isinstance(CUTTAG_CONFIG.get("seacr", {}), dict)
    else False
)

# Precompute CUT&Tag PE treatment samples eligible for SEACR
if SEACR_ENABLED:
    SEACR_SAMPLE_IDS = [
        sid for sid in TREATMENT_SAMPLE_IDS
        if SAMPLE_MAP[sid]["assay"] == "cuttag"
        and SAMPLE_MAP[sid]["layout"] == "PE"
    ]
else:
    SEACR_SAMPLE_IDS = []
```

SEACR only supports PE samples (requires paired-end fragments for fragment-aware
bedGraph generation). SE samples are excluded with no error — the rule simply
isn't scheduled for them.

## New Rules

### `seacr_bedgraph` — generate bedGraph from final.bam (in qc.smk or new file)

```python
rule seacr_bedgraph:
    output:
        f"{OUTDIR}/{{sample}}/04_peaks_seacr/{{sample}}.bedgraph",
    input:
        f"{OUTDIR}/{{sample}}/02_align/{{sample}}.final.bam",
    params:
        # bedtools genomecov -bg -pc (pair coverage) for proper fragment handling
        extra = "-bg -pc",
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        mkdir -p "$(dirname {output:q})"
        bedtools genomecov -ibam {input:q} -bg -pc > {output:q}
        """
```

Uses `bedtools genomecov -ibam final.bam -bg -pc`:
- `-bg`: UCSC bedGraph output (zero-signal regions omitted, as SEACR expects)
- `-pc`: pair coverage — counts each properly paired fragment once rather than
  each read. For PE CUT&Tag, this produces a fragment-density bedGraph matching
  the SEACR recommendation. The `-fs` flag is NOT used (it applies a fixed
  fragment-size extension which is not appropriate for variable-length CUT&Tag
  fragments).
- No control bedGraph is generated — Stage 7b uses SEACR's numeric threshold
  mode (FDR 0.01) for treatment-only peak calling.

bedGraph generation tests should verify UCSC 4-column format and that the
output can be consumed by SEACR (at minimum, the first few lines conform to
`chr start end value` with non-zero values).

Output directory: `04_peaks_seacr/sample/` — intentionally separate from `04_peaks/sample/` to keep MACS3 and SEACR namespaces distinct.

### `seacr_call` — run SEACR (in qc.smk or new file)

```python
rule seacr_call:
    output:
        f"{OUTDIR}/{{sample}}/04_peaks_seacr/{{sample}}/{{sample}}.seacr.{{mode}}.bed",
    input:
        bedgraph = f"{OUTDIR}/{{sample}}/04_peaks_seacr/{{sample}}.bedgraph",
    params:
        mode      = lambda wc: wc.mode,
        norm      = SEACR_NORMALIZATION,
        threshold = SEACR_THRESHOLD,
    conda:
        "../envs/chipseq.yml",
    shell:
        """
        set -e -o pipefail
        SEACR_SCRIPT=$(command -v SEACR_1.3.sh) || {{ \\
            echo "ERROR: SEACR_1.3.sh not found in PATH" >&2; exit 1; }}
        mkdir -p "$(dirname {output:q})"
        bash "$SEACR_SCRIPT" {input.bedgraph:q} {params.threshold} \\
            {params.norm:q} {params.mode:q} \\
            "$(dirname {output:q})/{wildcards.sample}.seacr"
        """
```

**SEACR CLI mapping:**

| Config key | SEACR positional arg | Notes |
|---|---|---|
| `--input` (bedgraph) | arg 1: target bedgraph | Always present |
| `cuttag.seacr.threshold` | arg 2: FDR threshold | From config; default 0.01 (top 1% AUC) |
| `cuttag.seacr.normalization` | arg 3: `"non"` only in 7b | No control bedGraph support yet |
| `cuttag.seacr.mode` | arg 4: `"stringent"` or `"relaxed"` | |
| output prefix | arg 5 | `{sample}.seacr` → produces `{sample}.seacr.stringent.bed` or `.relaxed.bed` |

Treatment-only mode: When no control bedGraph is provided, SEACR accepts a
numeric FDR threshold (0.01 = top 1% of AUC) in position 2. This avoids the
complexity of generating a control bedGraph in Stage 7b.

Snakemake output path: `{OUTDIR}/{sample}/04_peaks_seacr/{sample}/{sample}.seacr.{mode}.bed`.
Examples: `s1.seacr.stringent.bed`, `s1.seacr.relaxed.bed`.

### Rule placement

Both rules go in `workflow/rules/qc.smk` for now (they produce QC sidecar outputs,
not canonical peak files). If the file grows unwieldy, a future slice can extract
CUT&Tag-specific rules into a new file.

### Rule-all targets

```python
seacr_bedgraph_targets = (
    expand(
        "{outdir}/{sample}/04_peaks_seacr/{sample}.bedgraph",
        outdir=OUTDIR, sample=SEACR_SAMPLE_IDS,
    )
    if SEACR_ENABLED and SEACR_SAMPLE_IDS else []
),
seacr_call_targets = (
    expand(
        "{outdir}/{sample}/04_peaks_seacr/{sample}/{sample}.seacr.{mode}.bed",
        outdir=OUTDIR, sample=SEACR_SAMPLE_IDS,
        mode=[SEACR_MODE],
    )
    if SEACR_ENABLED and SEACR_SAMPLE_IDS else []
),
```

## Output Namespace

```
results/<sample>/
├── 04_peaks/                    # MACS3 canonical peaks (unchanged)
│   └── <sample>/
│       └── <sample>_peaks.narrowPeak
└── 04_peaks_seacr/              # Stage 7b NEW — SEACR sidecar outputs
    ├── <sample>.bedgraph        # SEACR input bedGraph
    └── <sample>/                # SEACR output directory
        └── <sample>.seacr.stringent.bed   # or .relaxed.bed
```

## Dependencies

### Conda environment

Add SEACR to `workflow/envs/chipseq.yml`:

```yaml
  # --- Stage 7 CUT&Tag ---
  - seacr
```

`seacr` is available on bioconda as `bioconda::seacr` (v1.3-2, noarch).
It requires `r-base` and `bedtools` on PATH — both already in the workflow env.

### bedtools

`bedtools genomecov` is already in the Conda environment (used by blacklist
filtering rules).

## Backward Compatibility

| Scenario | Behavior |
|---|---|
| No `cuttag` block in config | All keys default; `seacr.enabled` = false; zero DAG change |
| `cuttag.seacr.enabled: false` | No SEACR targets; zero DAG change |
| No CUT&Tag PE samples | `SEACR_SAMPLE_IDS` is empty; zero DAG change |
| ChIP-seq samples | Excluded (assay != "cuttag"), no SEACR targets |
| Existing MACS3 outputs | Completely unchanged — `04_peaks/` untouched |
| Existing QC rules | Unchanged — all still consume MACS3 outputs |
| Stage 5/6 IDR | Unchanged |

## Test Plan

### `test/test_stage7b_stress.py` — 11 tests

1. **No `cuttag` block** — validates, no SEACR targets in DAG
2. **`cuttag.seacr.enabled: false`** — no SEACR targets
3. **`cuttag.seacr.enabled: true` + CUT&Tag PE** — `seacr_bedgraph` and
   `seacr_call` in DAG; output paths contain `04_peaks_seacr/`
4. **`cuttag.seacr.enabled: true` + ChIP-seq** — no SEACR targets
5. **`cuttag.seacr.enabled: true` + CUT&Tag SE** — no SEACR targets
   (PE-only gating)
6. **Invalid mode rejected** — `cuttag.seacr.mode: "extreme"` fails validation
7. **Invalid normalization rejected** — `cuttag.seacr.normalization: "nope"`
   fails validation
8. **`cuttag.seacr.threshold: 0`** — rejected (must be in (0,1))
9. **`cuttag.seacr.threshold: 1`** — rejected (must be < 1)
10. **`cuttag.seacr.threshold: "bad"`** — rejected (not a float)
11. **`--list-rules`** — shows `seacr_bedgraph` and `seacr_call` when enabled

### `test/test_seacr_integration.py` — 2 integration tests (optional, deferred)

- **bedGraph generation** — a tiny PE SAM → BAM fixture, run
  `bedtools genomecov -bg -pc`, verify UCSC bedGraph format
- **SEACR output** — if SEACR is installed, run the actual SEACR script on
  the generated bedGraph, verify `.stringent.bed` is produced with 6 columns

These are marked deferred if SEACR dependency resolution is complex in CI.

## Files Changed

| File | Change |
|------|--------|
| `config/config.yaml` | Add commented `cuttag` block |
| `scripts/validate_samples.py` | Add `_validate_cuttag_config()` |
| `workflow/Snakefile` | Extract CUTTAG_CONFIG, SEACR_ENABLED, SEACR_SAMPLE_IDS; rule-all targets |
| `workflow/rules/qc.smk` | Add `seacr_bedgraph` and `seacr_call` rules |
| `workflow/envs/chipseq.yml` | Add `seacr` dependency |
| `workflow/schemas/config.schema.yaml` | Document `cuttag` block |
| `README.md` | Add Stage 7b mention |
| `KNOWN_ISSUES.md` | Update Stage 7 items |
| `test/test_stage7b_stress.py` | **New** — 8 tests |

## Verification

```bash
python3 scripts/validate_samples.py --config config/config.yaml
conda run -n chipseq snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet
python3 test/test_validation_stress.py
python3 test/test_stage7a_stress.py
python3 test/test_stage7b_stress.py
```
