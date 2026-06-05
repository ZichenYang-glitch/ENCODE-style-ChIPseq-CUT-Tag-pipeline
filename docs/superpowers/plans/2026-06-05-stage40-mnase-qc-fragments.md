# Stage 40: MNase Fragment Stratification and QC Summary — Implementation Plan

**Date:** 2026-06-05
**Status:** ready
**Spec:** `docs/superpowers/specs/2026-06-05-stage40-mnase-qc-fragments-design.md`

---

## Task Breakdown

### Task 1: Config validation — extend `_validate_mnase_config()`

**File:** `scripts/validate_samples.py`
**Lines:** ~408–456 (existing `_validate_mnase_config`)

**Changes:**
- Expand `known` set to `{"mono_range", "fragments", "dyad_range", "callers"}`.
- Add `fragments` validation: optional mapping. If present, validate each key (`sub`, `mono`, `di`) as `[int, int]` with `min < max`, both `> 0`.
- Add `dyad_range` validation: `[int, int]`, `min < max`, both `> 0`. Default `[130, 200]`.
- Add `callers` validation: optional mapping. Known keys: `danpos3`, `inps`, `sem`. All must be boolean. Any `true` → `ValidationError`.
- Return normalized dict: `{"mono_range": [...], "fragments": {...}, "dyad_range": [...], "callers": {...}}`.
- Backward compat: `mono_range` still accepted. If `fragments.mono` exists, it takes precedence in the returned dict (the helper `get_mono_range()` handles the fallback chain at read time).

**Test cases to add to `test/test_validation_stress.py`:**
- Default config (no mnase block) → all defaults applied
- `fragments.sub: [1, 139], mono: [140, 200], di: [300, 400]` → passes
- Old `mono_range: [130, 190]` only → passes
- `dyad_range: [140, 200]` → passes
- Missing dyad_range → passes (defaults)
- `fragments.sub: [200, 100]` → error
- `fragments.mono: [0, 200]` → error
- `dyad_range: [100]` → error
- `dyad_range: ["a", 200]` → error
- `callers.danpos3: true` → error with "not implemented"
- `callers.inps: true` → error with "not implemented"  
- `callers.sem: true` → error with "not implemented"
- `callers: {danpos3: false, inps: false, sem: false}` → passes
- Unknown caller key → error

---

### Task 2: Config schema — update `config.schema.yaml`

**File:** `workflow/schemas/config.schema.yaml`
**Lines:** ~304–321 (existing `mnase:` block)

**Changes:**
- Add `fragments` block documentation.
- Add `dyad_range` documentation.
- Add `callers` block documentation.
- Mark `mono_range` as "deprecated, use fragments.mono instead."
- Document default values.

---

### Task 3: Policy helpers — extend `workflow/rules/mnase.smk`

**File:** `workflow/rules/mnase.smk`
**Lines:** 1–68 (config helper + policy functions)

**Changes:**
- Refactor `get_mono_range()` to check `fragments.mono` first, then `mono_range`, then `[140, 200]`.
- Add `get_fragment_range(class_name)` — returns `(min, max)` for sub/mono/di. Never returns None.
- Add `get_dyad_range()` — returns `(min, max)`. Default `[130, 200]`.
- Add `get_caller_enabled(name)` — returns `True`/`False` from config (validation guarantees value is boolean and false, but helper reads it anyway).

---

### Task 4: Fragment rules — add `mnase_split_sub`, `mnase_split_di`

**File:** `workflow/rules/mnase.smk`
**Lines:** after existing `mnase_split_mono` (~line 100)

**Changes:**
- Add `mnase_split_sub` rule: same pattern as `mnase_split_mono`. Outputs: `<outdir>/<sample>/03_fragments/<sample>.sub.bam` + `.bai`. Uses `get_fragment_range("sub")`.
- Add `mnase_split_di` rule: same pattern. Outputs: `<outdir>/<sample>/03_fragments/<sample>.di.bam` + `.bai`. Uses `get_fragment_range("di")`.
- Existing `mnase_split_mono` unchanged in output paths, but `get_mono_range()` now respects `fragments.mono` precedence.

**Conda:** `../envs/deeptools.yml` (same as existing mono rule — provides alignmentSieve + samtools).

---

### Task 5: Dyad BigWig rules — add explicit fragment length flags

**File:** `workflow/rules/mnase.smk`
**Lines:** ~107–136 (`mnase_dyad_bigwig`), ~212–240 (`mnase_pooled_dyad_bigwig`)

**Changes:**
- Add `dyad_min` / `dyad_max` params calling `get_dyad_range()`.
- Add `--minFragmentLength {params.dyad_min} --maxFragmentLength {params.dyad_max}` to both rules' `bamCoverage` commands.
- Document in a comment that `tool_parameters.bamcoverage.extra_args` must NOT override these flags.

---

### Task 6: `bamCoverage` extra_args conflict guard

**File:** `scripts/validate_samples.py` or `workflow/rules/mnase.smk` comment

**Implementation note from user:**
When `bamCoverage --MNase` runs with explicit `--minFragmentLength/--maxFragmentLength`, a user might accidentally put conflicting flags in `tool_parameters.bamcoverage.extra_args`. 

**Decision:** Do not add a runtime validation guard for Stage 40 (no existing pattern for tool parameter conflict checking in the codebase). Instead, document in:
1. A comment in the dyad rules.
2. `docs/configuration.md` in the `bamcoverage` section.
3. The spec document.

The comment should read: "extra_args must not override --minFragmentLength, --maxFragmentLength, or --MNase. These are workflow-managed for MNase dyad BigWig."

---

### Task 7: MNase QC summary script

**File:** `scripts/mnase_qc_summary.py` (CREATE)

**Implementation:**
- stdlib-only: `argparse`, `csv`, `subprocess`, `os.path`.
- CLI flags: one per column (see spec Section 6.2).
- Read counts: `samtools view -c <bam>` for sub/mono/di BAMs. Writes `NA` if file missing or command fails.
- `insert_size_metrics` column: writes the path string if file exists, else `NA`.
- Writes single-row TSV with header to `--output`.

**Conda for the Snakemake rule:** `../envs/deeptools.yml` (provides Python 3 + samtools). `python.yml` does NOT include samtools, so `deeptools.yml` is required.

**Key constraint:** The script must NOT import any non-stdlib modules. `subprocess.run(["samtools", "view", "-c", bam_path])` is the only external tool call.

---

### Task 8: MNase QC summary Snakemake rule

**File:** `workflow/rules/mnase.smk`
**Lines:** after fragment rules

**Changes:**
- Add `mnase_qc_summary` rule (spec Section 6.3).
- Inputs: sub_bam, mono_bam, di_bam, dyad_bw, mono_bw (all concrete paths).
- Params: all column values from policy helpers + sample metadata + insert_size_metrics path + caller booleans.
- Conda: `../envs/deeptools.yml`.
- Shell: `python3 scripts/mnase_qc_summary.py ... --output {output}`.

---

### Task 9: DAG wiring — `_mnase_targets()` in Snakefile

**File:** `workflow/Snakefile`
**Lines:** ~860–907 (`_mnase_targets()`)

**Changes:**
- Add `expand()` calls for sub.bam, sub.bam.bai, di.bam, di.bam.bai, mnase_qc_summary.tsv.
- All gated on `MNASE_SAMPLE_IDS` only (no config gating needed — defaults guarantee targets are valid).

---

### Task 10: DAG wiring — `pipeline_done` in report.smk

**File:** `workflow/rules/report.smk`
**Lines:** ~20–95 (`pipeline_done` rule)

**Changes:**
- Add MNase lambda inputs: `mnase_sub_bam`, `mnase_sub_bai`, `mnase_di_bam`, `mnase_di_bai`, `mnase_qc_summ`.
- Each gated on `_is_mnase(wc)`.

---

### Task 11: Manifest — `make_manifest.py`

**File:** `scripts/make_manifest.py`
**Lines:** ~146–159 (existing MNase rows)

**Changes:**
- After existing MNase mono/dyad/monoBW rows, add 5 new rows for MNase samples (always, no config gating):
  - `mnase_sub_bam` → `alignmentSieve` → `03_fragments/<s>.sub.bam`
  - `mnase_sub_bai` → `samtools index` → `03_fragments/<s>.sub.bam.bai`
  - `mnase_di_bam` → `alignmentSieve` → `03_fragments/<s>.di.bam`
  - `mnase_di_bai` → `samtools index` → `03_fragments/<s>.di.bam.bai`
  - `mnase_qc_summary` → `mnase_qc_summary.py` → `01_qc/<s>.mnase_qc_summary.tsv`

---

### Task 12: MultiQC config — MNase QC custom content

**File:** `workflow/multiqc_config.yaml`

**Changes:**
- Add `custom_data.mnase_qc` section (spec Section 7.4).
- Add `sp.mnase_qc` section with `fn: "*mnase_qc_summary.tsv"`.
- Verify syntax against MultiQC 1.35 (the version in `workflow/envs/multiqc.yml`).
- Do NOT modify or remove existing `cross_correlation_qc` sections.

---

### Task 13: Tests — validation stress

**File:** `test/test_validation_stress.py`

**Changes:**
- Add Stage 40 validation test cases (see Task 1 list).
- Follow existing test patterns (tempdir, _write_yaml, _run_validate).

---

### Task 14: Tests — MNase stress (DAG/targets)

**File:** `test/test_stage39_mnase_stress.py`

**Changes:**
- Add assertions for Stage 40 rules in dry-run output:
  - `mnase_split_sub` scheduled
  - `mnase_split_di` scheduled
  - `mnase_qc_summary` scheduled
  - MACS3/FRiP/IDR still not scheduled
  - `pipeline.done` schedules sub/di/qc_summary
- Add test: default config (no mnase block) → all Stage 40 targets scheduled with defaults.

---

### Task 15: Tests — manifest stress

**File:** `test/test_stage25_manifest_stress.py`

**Changes:**
- Add assertions: MNase sample manifest includes `mnase_sub_bam`, `mnase_sub_bai`, `mnase_di_bam`, `mnase_di_bai`, `mnase_qc_summary`.

---

### Task 16: Tests — MultiQC config

**File:** `test/test_stage12_stress.py`

**Changes:**
- Add assertion: `multiqc_config.yaml` declares `mnase_qc` custom content.
- Add assertion: existing `cross_correlation_qc` still present.
- Add assertion: `sp.mnase_qc.fn` is `"*mnase_qc_summary.tsv"`.

---

### Task 17: Documentation — assay-policy.md

**File:** `docs/assay-policy.md`
**Lines:** ~156–229 (MNase section)

**Changes:**
- Update "Not implemented in v0.2" → move sub/di fragment classes and MNase QC summary to "implemented (Stage 40)."
- Update config example to show full `mnase:` block.
- Add note: `dyad_range` and `fragments.mono` can differ.
- Add note: caller config surface is reserved, execution deferred.

---

### Task 18: Documentation — configuration.md

**File:** `docs/configuration.md`

**Changes:**
- Document `mnase.fragments` block.
- Document `mnase.dyad_range`.
- Document `mnase.callers` surface.
- Document backward compatibility (mono_range still works).
- Document `extra_args` constraint for `bamcoverage` with MNase.

---

### Task 19: Documentation — output-contract.md

**File:** `docs/output-contract.md`

**Changes:**
- Add `mnase_sub_bam`, `mnase_sub_bai`, `mnase_di_bam`, `mnase_di_bai`, `mnase_qc_summary` to MNase single-sample outputs table.
- Mark all as "stable (Stage 40)."

---

### Task 20: Documentation — README.md

**File:** `README.md`

**Changes:**
- Update MNase limitations paragraph: remove "sub-nucleosome and di-nucleosome fragment classes" from deferred list.
- Add "MNase per-sample QC summary" to implemented features.
- Keep "Nucleosome calling (DANPOS3/iNPS/SEM)" in deferred list.

---

### Task 21: Documentation — CHANGELOG.md

**File:** `CHANGELOG.md`

**Changes:**
- Add Stage 40 entries under `[Unreleased]`:
  - `### Added`: sub-nucleosome and di-nucleosome fragment BAMs, MNase QC summary, explicit dyad_range config, caller config surface (execution deferred).
  - `### Changed`: `dyad_range` wired into sample and pooled dyad BigWig rules.

---

### Task 22: Documentation — Spec and plan files

**Files to create:**
- `docs/superpowers/specs/2026-06-05-stage40-mnase-qc-fragments-design.md` (already created)
- `docs/superpowers/plans/2026-06-05-stage40-mnase-qc-fragments.md` (this file)

---

## Execution Order

Tasks are ordered for incremental testability:

```
Task 1 (validate) → Task 2 (schema) → Task 3 (helpers)
    → Task 4 (sub/di rules) → Task 5 (dyad rules) → Task 6 (extra_args doc)
    → Task 7 (qc script) → Task 8 (qc rule)
    → Task 9 (targets) → Task 10 (pipeline_done)
    → Task 11 (manifest) → Task 12 (multiqc config)
    → Tasks 13–16 (tests) run after each code change
    → Tasks 17–22 (docs) last
```

## Test Commands

```bash
# 1. Validation
python3 test/test_validation_stress.py

# 2. Smoke profiles (requires SNAKEMAKE env var)
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake python3 test/test_stage8_smoke_profiles.py

# 3. MNase stress
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake python3 test/test_stage39_mnase_stress.py

# 4. Manifest stress
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake python3 test/test_stage25_manifest_stress.py

# 5. MultiQC/stage12
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake python3 test/test_stage12_stress.py

# 6. No hardcoded paths
python3 test/test_no_hardcoded_paths.py

# 7. Repo hygiene
git diff --check
```

## Stage 40 Config Example

```yaml
# Full mnase block (config/config.yaml)
mnase:
  fragments:
    sub: [1, 139]
    mono: [140, 200]
    di: [300, 400]
  dyad_range: [130, 200]
  callers:
    danpos3: false
    inps: false
    sem: false

# Minimal — all defaults apply
mnase: {}

# Legacy Stage 39 — still works
mnase:
  mono_range: [140, 200]
```
