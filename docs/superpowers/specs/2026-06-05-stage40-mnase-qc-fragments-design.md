# Stage 40: MNase Fragment Stratification and QC Summary — Design Spec

**Date:** 2026-06-05
**Status:** in design
**Scope:** Add sample-level sub/di fragment BAMs, explicit dyad_range, MNase QC summary
**Excluded:** Pooled sub/di, DANPOS3/iNPS/SEM execution, NFR/TSS MNase QC, rotational periodicity

---

## 1. Goals

1. Add `mnase.fragments` config block for sub/mono/di fragment ranges, with defaults so Stage 39 profiles automatically schedule Stage 40 outputs.
2. Add sample-level sub-nucleosome and di-nucleosome BAM rules (`mnase_split_sub`, `mnase_split_di`) alongside existing `mnase_split_mono`.
3. Add explicit `mnase.dyad_range` config key (default [130, 200]), wired into sample-level and pooled dyad BigWig rules via `--minFragmentLength` / `--maxFragmentLength`.
4. Add sample-level MNase QC summary (`scripts/mnase_qc_summary.py`, stdlib-only) producing `results/<sample>/01_qc/<sample>.mnase_qc_summary.tsv`.
5. Reserve `mnase.callers` config surface with validation that rejects `true` ("execution deferred"), no caller execution.
6. Wire Stage 40 outputs into `_mnase_targets()`, `pipeline_done`, manifest, and MultiQC custom content.
7. Zero new Conda environments. Zero new tools beyond existing `alignmentSieve`, `bamCoverage`, `samtools`, and Python stdlib.

**Not in Stage 40:** pooled sub/di, DANPOS3/iNPS/SEM execution, NFR/TSS profiles, rotational periodicity QC, Artifact model, AssayPolicy YAML, global target resolver.

---

## 2. Config Design: Thin MNase Policy

### 2.1 Full config surface

```yaml
mnase:
  # --- Legacy (Stage 39, still accepted) ---
  mono_range: [140, 200]

  # --- Stage 40 fragment stratification ---
  fragments:
    sub: [1, 139]
    mono: [140, 200]
    di: [300, 400]

  # --- Stage 40 explicit dyad range ---
  dyad_range: [130, 200]

  # --- Stage 40 caller surface (execution deferred) ---
  callers:
    danpos3: false
    inps: false
    sem: false
```

### 2.2 Defaults and backward compatibility

All three fragment classes always have values. The resolution chain:

| Class | Priority 1 | Priority 2 | Priority 3 (hard default) |
|:---|:---|:---|:---|
| `sub` | `mnase.fragments.sub` | — | `[1, 139]` |
| `mono` | `mnase.fragments.mono` | `mnase.mono_range` | `[140, 200]` |
| `di` | `mnase.fragments.di` | — | `[300, 400]` |

- `dyad_range`: `mnase.dyad_range` → default `[130, 200]`.
- `mono_range`: still accepted (not an unknown-key error). Used only as fallback for mono fragment range.
- No `enabled` / `disabled` switches. If a user removes `fragments.sub` from config, it silently falls back to the default, not to "no sub BAM."

### 2.3 Caller validation

- `callers` block: optional. If present, known keys are `danpos3`, `inps`, `sem`.
- All must be boolean.
- If any is `true`: `ValidationError("caller execution is not implemented in v0.2; set <caller>: false")`.
- If all are `false` or the block is absent: validation passes.
- No placeholder targets, no dummy rules, no `rule all` entries for callers.

### 2.4 Validation rules (extend `_validate_mnase_config`)

| # | Condition | Error |
|---|-----------|-------|
| 1 | `mnase` not a mapping | `mnase must be a mapping` |
| 2 | Unknown key (not in `mono_range, fragments, dyad_range, callers`) | `unknown key` |
| 3 | `mono_range` wrong shape | Existing Stage 39 checks |
| 4 | `fragments` not a mapping (when present) | `mnase.fragments must be a mapping` |
| 5 | `fragments.<class>` not `[int, int]` with `min < max`, both > 0 | `invalid fragment range` |
| 6 | `dyad_range` not `[int, int]` with `min < max`, both > 0 | `invalid dyad_range` |
| 7 | `callers` not a mapping (when present) | `callers must be a mapping` |
| 8 | Unknown caller key | `unknown caller key` |
| 9 | Caller value is `true` | `caller execution is not implemented` |
| 10 | All defaults (mnase block absent) | All fragment ranges = defaults, dyad_range = [130,200], callers = all false |

---

## 3. Policy Helpers (in mnase.smk)

Three helpers replace/augment the existing `get_mono_range()`:

```python
def get_mono_range():
    """Existing helper — extended with fragments.mono precedence."""
    mnase_cfg = config.get("mnase", {})
    if isinstance(mnase_cfg, dict):
        fragments = mnase_cfg.get("fragments", {})
        if isinstance(fragments, dict) and "mono" in fragments:
            r = fragments["mono"]
            if isinstance(r, (list, tuple)) and len(r) == 2:
                return int(r[0]), int(r[1])
        mr = mnase_cfg.get("mono_range", [140, 200])
        if isinstance(mr, (list, tuple)) and len(mr) == 2:
            return int(mr[0]), int(mr[1])
    return 140, 200


def get_fragment_range(class_name):
    """Return (min, max) for sub/mono/di. Never returns None.

    class_name: "sub" | "mono" | "di"
    """
    DEFAULTS = {"sub": (1, 139), "mono": (140, 200), "di": (300, 400)}
    if class_name == "mono":
        return get_mono_range()
    mnase_cfg = config.get("mnase", {})
    if isinstance(mnase_cfg, dict):
        fragments = mnase_cfg.get("fragments", {})
        if isinstance(fragments, dict) and class_name in fragments:
            r = fragments[class_name]
            if isinstance(r, (list, tuple)) and len(r) == 2:
                return int(r[0]), int(r[1])
    return DEFAULTS[class_name]


def get_dyad_range():
    """Return (min, max) for dyad BigWig. Default [130, 200]."""
    mnase_cfg = config.get("mnase", {})
    if isinstance(mnase_cfg, dict):
        dr = mnase_cfg.get("dyad_range", [130, 200])
        if isinstance(dr, (list, tuple)) and len(dr) == 2:
            return int(dr[0]), int(dr[1])
    return 130, 200
```

---

## 4. Fragment Rules

### 4.1 New rules

Two new rules following the exact pattern of `mnase_split_mono`:

| Rule | Output | Fragment range |
|:---|:---|:---|
| `mnase_split_sub` | `<outdir>/<sample>/03_fragments/<sample>.sub.bam` + `.bai` | `get_fragment_range("sub")` |
| `mnase_split_di` | `<outdir>/<sample>/03_fragments/<sample>.di.bam` + `.bai` | `get_fragment_range("di")` |

Each rule:
- Input: `final.bam` + `final.bam.bai`
- Tool: `alignmentSieve --minFragmentLength {min} --maxFragmentLength {max}`
- Threads: `THREADS`
- Conda: `../envs/deeptools.yml`
- Shell: `set -e -o pipefail; mkdir -p ...; alignmentSieve ...; samtools index ...`

### 4.2 Existing rule unchanged

`mnase_split_mono` — output path, rule name, and shell unchanged. `get_mono_range()` now considers `fragments.mono` first, but the output contract is identical.

### 4.3 No wildcard rule

The three fragment rules are independent. They target distinct output paths (`sub.bam` vs `mono.bam` vs `di.bam`). No ambiguity.

### 4.4 No pooled sub/di

Pooled fragment outputs (`mnase_pooled_mono`, etc.) remain mono-only. Pooled sub/di is deferred beyond Stage 40.

---

## 5. Dyad BigWig Rules

### 5.1 Sample-level (`mnase_dyad_bigwig`)

Add params and CLI flags:

```python
params:
    dyad_min = get_dyad_range()[0],
    dyad_max = get_dyad_range()[1],
    ...
```

Shell change:
```bash
bamCoverage \
    -b {input.bam} -o {output.bw} \
    --MNase --binSize 1 \
    --minFragmentLength {params.dyad_min} \
    --maxFragmentLength {params.dyad_max} \
    {params.normalize_using} \
    --numberOfProcessors {threads} \
    {params.extra}
```

### 5.2 Pooled (`mnase_pooled_dyad_bigwig`)

Same pattern — add `dyad_min` / `dyad_max` params, pass `--minFragmentLength` / `--maxFragmentLength`.

### 5.3 Default behavior preservation

With default `dyad_range: [130, 200]`, the behavior is identical to Stage 39 (deepTools `--MNase` default is 130–200 bp). This is a no-op for existing runs, but makes the range explicit and configurable.

### 5.4 Documentation note

The spec and docs must state: `dyad_range` and `fragments.mono` can differ. A user may choose `fragments.mono: [140, 200]` but `dyad_range: [130, 200]` to capture a broader dyad signal than the mono occupancy window.

---

## 6. MNase QC Summary

### 6.1 Script: `scripts/mnase_qc_summary.py`

- Python stdlib only (`argparse`, `csv`, `subprocess`, `os.path`).
- Accepts CLI args for all column values (see below).
- Computes read counts via `samtools view -c <bam>` for sub/mono/di BAMs.
- Writes single-row TSV to `--output`.
- Does NOT depend on Picard metrics as a required input.

### 6.2 Output columns

```
sample  assay  peak_mode
sub_min  sub_max  mono_min  mono_max  di_min  di_max
dyad_min  dyad_max
sub_bam  mono_bam  di_bam  dyad_bigwig  mono_bigwig
insert_size_metrics
caller_danpos3_enabled  caller_inps_enabled  caller_sem_enabled
sub_reads  mono_reads  di_reads
```

- `sub_reads` / `mono_reads` / `di_reads`: integer from `samtools view -c`, or `NA` if BAM missing.
- `insert_size_metrics`: path string or `NA` if Picard not configured / file absent.
- Caller columns: `true` / `false` strings.
- `*_min` / `*_max`: integers from config.
- BAM/bw paths: resolved absolute paths (as written by Snakemake).

### 6.3 Snakemake rule: `mnase_qc_summary`

```python
rule mnase_qc_summary:
    output:
        f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.mnase_qc_summary.tsv"
    input:
        sub_bam  = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.sub.bam",
        mono_bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam",
        di_bam   = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.di.bam",
        dyad_bw  = f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.dyad.CPM.bw",
        mono_bw  = f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.mono.CPM.bw",
    params:
        sample = "{sample}",
        assay  = lambda wc: SAMPLE_MAP[wc.sample]["assay"],
        peak_mode = lambda wc: SAMPLE_MAP[wc.sample]["peak_mode"],
        sub_min  = get_fragment_range("sub")[0],
        sub_max  = get_fragment_range("sub")[1],
        mono_min = get_fragment_range("mono")[0],
        mono_max = get_fragment_range("mono")[1],
        di_min   = get_fragment_range("di")[0],
        di_max   = get_fragment_range("di")[1],
        dyad_min = get_dyad_range()[0],
        dyad_max = get_dyad_range()[1],
        insert_size_metrics = f"{OUTDIR}/{{sample}}/05_qc/picard/{{sample}}.insert_size_metrics",
        danpos3_enabled = ...
        inps_enabled = ...
        sem_enabled = ...
    conda:
        "../envs/python.yml"
    shell:
        """
        python3 scripts/mnase_qc_summary.py \
            --sample {params.sample:q} \
            --assay {params.assay:q} \
            --peak-mode {params.peak_mode:q} \
            --sub-min {params.sub_min} --sub-max {params.sub_max} \
            --mono-min {params.mono_min} --mono-max {params.mono_max} \
            --di-min {params.di_min} --di-max {params.di_max} \
            --dyad-min {params.dyad_min} --dyad-max {params.dyad_max} \
            --sub-bam {input.sub_bam:q} \
            --mono-bam {input.mono_bam:q} \
            --di-bam {input.di_bam:q} \
            --dyad-bw {input.dyad_bw:q} \
            --mono-bw {input.mono_bw:q} \
            --insert-size-metrics {params.insert_size_metrics:q} \
            --danpos3-enabled {params.danpos3_enabled} \
            --inps-enabled {params.inps_enabled} \
            --sem-enabled {params.sem_enabled} \
            --output {output:q}
        """
```

Key design points:
- All inputs are concrete paths (no conditional empty strings — sub/di always exist).
- `insert_size_metrics` is a params path, not an input. The script checks existence; writes `NA` if absent.
- No dependency on Picard being enabled — the rule always runs, column reflects availability.

---

## 7. DAG Wiring

### 7.1 `_mnase_targets()` in Snakefile

Add to the existing `_mnase_targets()`:
```python
targets += expand("{outdir}/{sample}/03_fragments/{sample}.sub.bam", ...)
targets += expand("{outdir}/{sample}/03_fragments/{sample}.sub.bam.bai", ...)
targets += expand("{outdir}/{sample}/03_fragments/{sample}.di.bam", ...)
targets += expand("{outdir}/{sample}/03_fragments/{sample}.di.bam.bai", ...)
targets += expand("{outdir}/{sample}/01_qc/{sample}.mnase_qc_summary.tsv", ...)
```

All gated on `MNASE_SAMPLE_IDS` only (no config gating — defaults guarantee these always apply).

### 7.2 `pipeline_done` in report.smk

Add to the MNase conditional lambda inputs:
```python
mnase_sub_bam = lambda wc: (
    f"{OUTDIR}/{wc.sample}/03_fragments/{wc.sample}.sub.bam"
    if _is_mnase(wc) else []
),
mnase_sub_bai = lambda wc: (...),
mnase_di_bam = lambda wc: (...),
mnase_di_bai = lambda wc: (...),
mnase_qc_summ = lambda wc: (
    f"{OUTDIR}/{wc.sample}/01_qc/{wc.sample}.mnase_qc_summary.tsv"
    if _is_mnase(wc) else []
),
```

This ensures `pipeline.done` depends on all Stage 40 MNase outputs.

### 7.3 MultiQC dependency

The `multiqc` rule already depends on `pipeline.done` for treatment samples. `pipeline.done` now depends on `mnase_qc_summary.tsv`. So `multiqc` automatically waits for the MNase QC summary via the existing DAG chain. No change needed in the multiqc rule's input list.

### 7.4 MultiQC custom content

Add to `workflow/multiqc_config.yaml` (alongside existing `cross_correlation_qc`):

```yaml
custom_data:
  mnase_qc:
    section_name: "MNase-seq QC"
    description: "Per-sample MNase-seq fragment stratification and QC summary."
    plot_type: "table"
    pconfig:
      id: "mnase_qc"
      title: "MNase-seq QC Summary"
      col1_header: "Sample"

sp:
  mnase_qc:
    fn: "*mnase_qc_summary.tsv"
```

The `sp.mnase_qc.fn` glob matches per-sample `mnase_qc_summary.tsv` files in the search paths (already includes each sample's `results/<sample>/` directory).

Do NOT modify or remove the existing `cross_correlation_qc` sections.

---

## 8. Manifest (`make_manifest.py`)

For every MNase treatment sample, always emit rows:

| output_type | method | path |
|:---|:---|:---|
| `mnase_sub_bam` | `alignmentSieve` | `results/<s>/03_fragments/<s>.sub.bam` |
| `mnase_sub_bai` | `samtools index` | `results/<s>/03_fragments/<s>.sub.bam.bai` |
| `mnase_di_bam` | `alignmentSieve` | `results/<s>/03_fragments/<s>.di.bam` |
| `mnase_di_bai` | `samtools index` | `results/<s>/03_fragments/<s>.di.bam.bai` |
| `mnase_qc_summary` | `mnase_qc_summary.py` | `results/<s>/01_qc/<s>.mnase_qc_summary.tsv` |

Existing MNase mono/dyad/monoBW rows unchanged.

---

## 9. Test Plan

### 9.1 Validation stress (extend `test_validation_stress.py`)

| # | Case | Expected |
|---|------|----------|
| 1 | Default config (no mnase block) → all defaults applied | PASS |
| 2 | `mnase.fragments.sub: [1, 139], mono: [140, 200], di: [300, 400]` | PASS |
| 3 | Old `mnase.mono_range: [130, 190]` only → mono = [130,190], sub/di = defaults | PASS |
| 4 | Both `fragments.mono` and `mono_range` → fragments.mono wins | PASS |
| 5 | `dyad_range: [140, 200]` | PASS |
| 6 | Missing dyad_range → default [130, 200] | PASS |
| 7 | `fragments.sub: [200, 100]` (min >= max) → error | FAIL |
| 8 | `fragments.mono: [0, 200]` (not positive) → error | FAIL |
| 9 | `dyad_range: [100]` (not length 2) → error | FAIL |
| 10 | `dyad_range: ["a", 200]` (not int) → error | FAIL |
| 11 | `callers.danpos3: true` → error with "not implemented" | FAIL |
| 12 | `callers.inps: true` → error with "not implemented" | FAIL |
| 13 | `callers.sem: true` → error with "not implemented" | FAIL |
| 14 | `callers: {danpos3: false, inps: false, sem: false}` → passes | PASS |
| 15 | Unknown caller key → error | FAIL |
| 16 | MNase peak_mode=narrow/broad still fails | FAIL |
| 17 | MNase SE layout still fails | FAIL |

### 9.2 Stage 39/40 stress (extend `test_stage39_mnase_stress.py`)

| # | Case | Expected |
|---|------|----------|
| 1 | MNase dry-run schedules `mnase_split_sub`, `mnase_split_di` | PASS |
| 2 | MNase dry-run schedules `mnase_qc_summary` | PASS |
| 3 | MNase dry-run does NOT schedule MACS3/FRiP/IDR | PASS |
| 4 | Direct `pipeline.done` schedules all Stage 40 MNase outputs | PASS |
| 5 | Dyad BigWig rule uses `--minFragmentLength` / `--maxFragmentLength` params | PASS |
| 6 | Default config (no mnase block) → all MNase targets scheduled with defaults | PASS |

### 9.3 Manifest stress (extend `test_stage25_manifest_stress.py`)

| # | Case | Expected |
|---|------|----------|
| 1 | MNase sample manifest includes sub_bam, sub_bai, di_bam, di_bai, qc_summary | PASS |
| 2 | Non-MNase sample manifest does NOT include MNase rows | PASS |

### 9.4 MultiQC config (extend `test_stage12_stress.py`)

| # | Case | Expected |
|---|------|----------|
| 1 | `multiqc_config.yaml` declares `mnase_qc` custom content | PASS |
| 2 | Existing `cross_correlation_qc` custom content still present | PASS |

### 9.5 Regression safety

All existing suites must pass unchanged:
- `test_validation_stress.py` (existing + new Stage 40 cases)
- `test_stage8_smoke_profiles.py` (MNase profile + 7 others)
- `test_stage22_bigwig_stress.py`
- `test_stage24_qc_summary_unit.py`
- `test_stage25_manifest_stress.py`
- `test_no_hardcoded_paths.py`
- `test_stage27_public_validation_plan.py`
- `test_stage27b_metadata_ci_plan.py`
- `test_stage27c_ci_workflow.py`
- `test_stage28_release_readiness.py`

---

## 10. Files Summary

| File | Action | Purpose |
|:---|:---|:---|
| `scripts/validate_samples.py` | Modify | Extend `_validate_mnase_config()` |
| `workflow/schemas/config.schema.yaml` | Modify | Document new mnase keys |
| `workflow/rules/mnase.smk` | Modify | Policy helpers, sub/di rules, dyad params, qc_summary rule |
| `workflow/Snakefile` | Modify | `_mnase_targets()` with Stage 40 targets |
| `workflow/rules/report.smk` | Modify | `pipeline_done` MNase inputs |
| `scripts/mnase_qc_summary.py` | **Create** | stdlib-only QC summary |
| `scripts/make_manifest.py` | Modify | Add sub/di/qc_summary manifest rows |
| `workflow/multiqc_config.yaml` | Modify | MNase QC custom content |
| `docs/assay-policy.md` | Modify | Update MNase section |
| `docs/configuration.md` | Modify | Document new config keys |
| `docs/output-contract.md` | Modify | Add new output types |
| `README.md` | Modify | Update MNase features / limitations |
| `CHANGELOG.md` | Modify | Stage 40 entry |
| `test/test_validation_stress.py` | Modify | Stage 40 validation cases |
| `test/test_stage39_mnase_stress.py` | Modify | Stage 40 DAG/target tests |
| `test/test_stage25_manifest_stress.py` | Modify | MNase sub/di/qc manifest rows |
| `test/test_stage12_stress.py` | Modify | MNase QC custom content assertions |

---

## 11. Risks

| Risk | Mitigation |
|:---|:---|
| `get_mono_range()` changed behavior breaks Stage 39 | New chain: fragments.mono → mono_range → [140,200]. With neither new keys set, result is identical to Stage 39. |
| sub/di rules create large BAMs | Same read count as mono rule (per-fragment, not duplication). Users already handle mono BAM size. |
| MNase QC summary script bugs | Stdlib-only, single-file, CLI-driven. Easy to test in isolation. `samtools view -c` failure writes `NA`. |
| MultiQC MNase section conflicts with cross_correlation | Separate `custom_data` keys and `sp` entries. Independent. |
| Test profile dry-run now schedules more rules | MNase test profile was already dry-run only. Extra rules in dry-run add negligible time. |
