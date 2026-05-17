# Stage 6b: Histone ChIP-seq Policy and Pooled QC Summary

## Scope

Stage 6b adds a histone target classification helper and a pooled experiment QC
summary. It supports histone ChIP-seq expectations without forcing histone assays
through TF-style IDR.

**Stage 6b includes:**
- `scripts/histone_utils.py` — histone target classification helper
- `pooled_experiment_qc_summary` rule in `qc.smk`
- Per-experiment QC TSV output under `results/experiments/<exp>/01_qc/`

**Deferred:**
- Pooled FRiP, NRF/PBC, library complexity aggregation
- Per-sample QC aggregation into experiment summaries
- MultiQC integration of experiment-level QC
- Auto-flagging or automatic peak_mode adjustment

## Configuration

No new config switches. Stage 6b is gated on `STAGE4B` and
`MULTI_BIOREP_EXPERIMENTS` — the same infrastructure that gates Stage 4b
pooled BAMs and peaks. When these conditions are met, the pooled QC
summary is always produced.

## `scripts/histone_utils.py`

Pure stdlib module. No file I/O, no Snakemake dependency. Directly unit-tested.

### API

```python
def classify_histone_target(target, configured_peak_mode=None):
    """Classify a histone target for QC purposes.

    Returns a dict with keys:
      target_normalized       — uppercase, stripped, de-punctuated
      inferred_histone_class  — "broad_like" | "narrow_like" |
                                 "context_dependent" | "unknown"
      expected_peak_mode      — "broad" | "narrow" | "broad_or_narrow" |
                                 "unknown"
      compatible_peak_modes   — ["broad"], ["narrow"], ["broad","narrow"], []
      peak_mode_status        — "ok" | "mismatch" | "unknown"

    Does NOT raise exceptions. Mismatches are reported via status only.
    """
```

### Target Groups

**BROAD_LIKE** (expected_peak_mode: "broad", compatible: ["broad"]):

H3K27me3, H3K36me3, H3K9me3, H3K79me2, H3K79me3,
H4K20me1, H4K20me3

**NARROW_LIKE** (expected_peak_mode: "narrow", compatible: ["narrow"]):

H3K4me3, H3K9ac, H3K14ac, H3K18ac, H3K23ac, H3K56ac,
H4K5ac, H4K8ac, H4K12ac, H4K16ac, H2AZ, H2A.Z, H2AFZ

**CONTEXT_DEPENDENT** (expected_peak_mode: "broad_or_narrow",
compatible: ["broad", "narrow"]):

H3K27ac, H3K4me1, H3K4me2

**UNKNOWN** (expected_peak_mode: "unknown", compatible: []):

Anything not matching the above.

### peak_mode_status logic

| inferred_class | configured == expected? | status |
|---|---|---|
| broad_like, peak_mode=="broad" | yes | ok |
| broad_like, peak_mode=="narrow" | no | mismatch |
| narrow_like, peak_mode=="narrow" | yes | ok |
| narrow_like, peak_mode=="broad" | no | mismatch |
| context_dependent, any peak_mode | n/a (compatible) | ok |
| unknown | n/a | unknown |

### Matching

- Case-insensitive: uppercase input, strip whitespace.
- Remove hyphens and underscores before matching.
- `H2A.Z` and `H2AZ` both match via explicit entries for both forms.
- `h3k27me3`, `H3K27ME3`, `H3K27ac` all resolve correctly.

## `pooled_experiment_qc_summary` Rule (in `qc.smk`)

### Output

```
results/experiments/{experiment}/01_qc/{experiment}.pooled_qc_summary.tsv
```

### Inputs

```python
input:
    pooled_bam   = f"{OUTDIR}/experiments/{{exp}}/02_align/{{exp}}.pooled.final.bam",
    pooled_peaks = f"{OUTDIR}/experiments/{{exp}}/04_peaks/pooled/{{exp}}_pooled_peaks",
    pooled_fe    = (f"{OUTDIR}/experiments/{{exp}}/03_signal/{{exp}}.pooled.FE.bdg"
                    if QC_CONFIG.get("signal_tracks", True) else []),
    pooled_ppois = (f"{OUTDIR}/experiments/{{exp}}/03_signal/{{exp}}.pooled.ppois.bdg"
                    if QC_CONFIG.get("signal_tracks", True) else []),
```

When `qc.signal_tracks` is true, the FE and ppois bedGraphs are input dependencies
so Snakemake builds them before the summary. When false, they are not inputs.

### Params

```python
params:
    experiment              = "{experiment}",
    target                  = lambda wc: SAMPLE_MAP[
                               TREATMENT_SAMPLES_BY_EXPERIMENT[wc.experiment][0]]["target"],
    assay                   = lambda wc: SAMPLE_MAP[
                               TREATMENT_SAMPLES_BY_EXPERIMENT[wc.experiment][0]]["assay"],
    peak_mode               = lambda wc: SAMPLE_MAP[
                               TREATMENT_SAMPLES_BY_EXPERIMENT[wc.experiment][0]]["peak_mode"],
    n_bioreps               = lambda wc: len(_bioreps_for(wc.experiment, "treatment")),
    bio_rep_labels          = lambda wc: ",".join(
        str(b) for b in sorted(_bioreps_for(wc.experiment, "treatment"))
    ),
    inferred_histone_class  = lambda wc: _classify(wc.experiment)["inferred_histone_class"],
    expected_peak_mode      = lambda wc: _classify(wc.experiment)["expected_peak_mode"],
    peak_mode_status        = lambda wc: _classify(wc.experiment)["peak_mode_status"],
    signal_enabled          = QC_CONFIG.get("signal_tracks", True),
    pooled_fe_path          = _pooled_fe_path if QC_CONFIG.get(...) else "NA",
    pooled_ppois_path       = _pooled_ppois_path if QC_CONFIG.get(...) else "NA",
```

### Shell (inline peak counting)

```bash
set -e -o pipefail
mkdir -p "$(dirname {output})"

# Resolve peak file by peak_mode
if [[ "{params.peak_mode}" == "broad" ]]; then
    PK_FILE="{input.pooled_peaks}/{wildcards.experiment}_pooled_peaks.broadPeak"
else
    PK_FILE="{input.pooled_peaks}/{wildcards.experiment}_pooled_peaks.narrowPeak"
fi

if [[ ! -f "$PK_FILE" ]]; then
    echo "ERROR: pooled peak file not found: $PK_FILE" >&2
    exit 1
fi

PK_COUNT=$(wc -l < "$PK_FILE")

# Signal tracks status
if [[ "{params.signal_enabled}" == "True" ]]; then
    SIG_STATUS="enabled_present"
else
    SIG_STATUS="disabled"
fi

# Write TSV (14 columns)
printf "experiment\tassay\ttarget\tinferred_histone_class\t" \
       "expected_peak_mode\tconfigured_peak_mode\tpeak_mode_status\t" \
       "biological_replicates\tbiological_replicate_labels\t" \
       "pooled_bam\tpooled_peaks\tpooled_peak_count\t" \
       "pooled_FE_bdg\tpooled_ppois_bdg\tsignal_tracks_status\n" \
       > {output}

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    {params.experiment} \
    {params.assay} \
    {params.target} \
    {params.inferred_histone_class} \
    {params.expected_peak_mode} \
    {params.peak_mode} \
    {params.peak_mode_status} \
    {params.n_bioreps} \
    {params.bio_rep_labels} \
    {input.pooled_bam} \
    {input.pooled_peaks} \
    "$PK_COUNT" \
    {params.pooled_fe_path} \
    {params.pooled_ppois_path} \
    "$SIG_STATUS" \
    >> {output}
```

Histone classification fields are extracted into separate scalar params
(`inferred_histone_class`, `expected_peak_mode`, `peak_mode_status`) —
no dictionary key access needed in the shell block.

### signal_tracks_status

- `disabled`: `qc.signal_tracks=false` — no FE/ppois inputs, summary still produced
- `enabled_present`: `qc.signal_tracks=true` — FE/ppois are input dependencies,
  Snakemake builds them before the summary; they should be present at runtime
- `enabled_missing`: defensive state for manual invocation or edge cases where
  signal tracks are enabled but files are absent; normally Snakemake input
  dependencies prevent this state from occurring

### Why H3K27ac, H3K4me1, H3K4me2 Are `context_dependent`

These marks can produce either narrow or broad peaks depending on biological
context. H3K27ac is commonly analyzed as narrow peaks at promoters/enhancers
but can show broad domain-like enrichment. H3K4me1 and H3K4me2 mark enhancers
and promoters respectively but can also appear in broad domains. Stage 6b
accepts both `broad` and `narrow` peak_mode for these marks and reports
`peak_mode_status = ok` for either. No warning is emitted.

### TSV Columns (15 columns)

| # | Column | Source |
|---|--------|--------|
| 1 | `experiment` | wildcard |
| 2 | `assay` | sample sheet |
| 3 | `target` | sample sheet |
| 4 | `inferred_histone_class` | `classify_histone_target()` |
| 5 | `expected_peak_mode` | `classify_histone_target()` |
| 6 | `configured_peak_mode` | sample sheet |
| 7 | `peak_mode_status` | `classify_histone_target()` |
| 8 | `biological_replicates` | `len(_bioreps_for(...))` |
| 9 | `biological_replicate_labels` | `"2,4"` style |
| 10 | `pooled_bam` | input path |
| 11 | `pooled_peaks` | input directory |
| 12 | `pooled_peak_count` | `wc -l` on peak file |
| 13 | `pooled_FE_bdg` | param path or "NA" |
| 14 | `pooled_ppois_bdg` | param path or "NA" |
| 15 | `signal_tracks_status` | `enabled_present` / `enabled_missing` / `disabled` |

## Snakefile Changes

### Import

At the top of the Snakefile (after the existing `from validate_samples import`):

```python
from histone_utils import classify_histone_target
```

### Rule-all target

```python
pooled_qc_summary = (
    expand(
        "{outdir}/experiments/{experiment}/01_qc/{experiment}.pooled_qc_summary.tsv",
        outdir=OUTDIR,
        experiment=MULTI_BIOREP_EXPERIMENTS,
    )
    if STAGE4B and MULTI_BIOREP_EXPERIMENTS else []
),
```

## Test Plan

### `test/test_histone_utils.py` — unit tests (8 tests)

1. H3K27me3 + broad → broad_like, ok
2. H3K27me3 + narrow → broad_like, mismatch
3. H3K27ac + broad → context_dependent, ok
4. H3K27ac + narrow → context_dependent, ok
5. H3K4me3 + narrow → narrow_like, ok
6. CTCF + narrow → unknown, unknown
7. h3k27me3 (lowercase) → broad_like (case-insensitive)
8. H2A.Z → narrow_like (punctuation variant)

### `test/test_stage6b_stress.py` — integration tests (7 tests)

Reuses Stage 4b/5a/6a harness pattern.

9. Multi-biorep → `pooled_experiment_qc_summary` in DAG
10. Single-sample → no pooled QC summary in DAG
11. `qc.signal_tracks: false` → summary produced, status=disabled
12. Broad pooled experiment → correct broadPeak path in summary
13. Dry-run compatible with Stage 5 disabled (`stage5: false`)
14. `pooled_experiment_qc_summary` appears in `--list-rules`
15. Inline peak count using `wc -l` visible in dry-run -p

## Files Changed

| File | Change |
|------|--------|
| `scripts/histone_utils.py` | **New** — histone target classification helper |
| `workflow/Snakefile` | Import classify_histone_target, add rule-all target |
| `workflow/rules/qc.smk` | Add `pooled_experiment_qc_summary` rule |
| `README.md` | Add Stage 6b section |
| `KNOWN_ISSUES.md` | Mark Stage 6b items completed |
| `test/test_histone_utils.py` | **New** — 8 unit tests |
| `test/test_stage6b_stress.py` | **New** — 7 integration tests |

## Verification

```bash
python3 scripts/validate_samples.py --config config/config.yaml
conda run -n chipseq snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet
python3 test/test_validation_stress.py
python3 test/test_stage4b_stress.py
python3 test/test_stage6a_stress.py
python3 test/test_histone_utils.py
python3 test/test_stage6b_stress.py
```
