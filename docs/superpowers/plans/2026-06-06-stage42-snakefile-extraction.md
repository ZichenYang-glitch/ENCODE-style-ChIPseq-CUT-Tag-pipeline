# Stage 42: Snakefile Responsibility Extraction — Implementation Plan

**Date:** 2026-06-06
**Status:** implemented
**Spec:** `docs/superpowers/specs/2026-06-06-stage42-snakefile-extraction-design.md`

---

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development or superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Extract derived metadata + dispatch functions into `metadata.smk` and target-builder functions into `targets.smk`. Reduce `workflow/Snakefile` from ~1037 to ~200 lines with zero behavior change.

**Architecture:** Snakemake `include:` literal text inclusion. The extraction is a pure move — no variable renaming, no logic changes, no path changes. `include:` order preserves the dependency chain: assay policy → metadata (dispatch) → targets → rule all → executable rules.

**Tech Stack:** Snakemake, Python (stdlib)

---

## Pre-requisite: Stage 41 merge + branch setup

Stage 42 MUST NOT start before Stage 41 (`stage41-artifact-readiness`) is merged to main.

```bash
git fetch origin
git switch main
git reset --hard origin/main
git switch -c stage42-snakefile-extraction
```

---

### Task 1: Create `workflow/rules/metadata.smk`

**Files:**
- Create: `workflow/rules/metadata.smk`

Cut the following line ranges from `workflow/Snakefile` and paste into `metadata.smk`:

| Lines | Content |
|:---|---:|
| 142–266 | Stage 4a/b experiment groups, biorep structures, pooled experiment lists |
| 268–331 | Stage 5a/b IDR precomputed lists |
| 333–442 | QC_CONFIG, TOOL_PARAMS, blacklist/signal/TSS gating, genome resource helpers |
| 452–525 | Dispatch wrappers + `_macs3_inputs()` |

Note: lines 444–450 (assay policy includes) are NOT moved — they stay in the Snakefile.

- [ ] **Step 1: Verify Snakefile line ranges are correct**

```bash
sed -n '142,266p' workflow/Snakefile | head -3
# Expected: first line of EXPERIMENT_IDS section
sed -n '452,525p' workflow/Snakefile | head -3
# Expected: def get_remove_dup
```

- [ ] **Step 2: Create metadata.smk from the extracted ranges**

```bash
# Extract the two blocks (gap at lines 444-450 is the assay policy includes)
sed -n '142,442p' workflow/Snakefile > workflow/rules/metadata.smk
echo "" >> workflow/rules/metadata.smk
sed -n '452,525p' workflow/Snakefile >> workflow/rules/metadata.smk
```

- [ ] **Step 3: Verify metadata.smk contains no rules**

```bash
grep -c "^rule " workflow/rules/metadata.smk
# Expected: 0
```

---

### Task 2: Create `workflow/rules/targets.smk`

**Files:**
- Create: `workflow/rules/targets.smk`

Cut lines 527–1002 from `workflow/Snakefile` and paste into `targets.smk`.

- [ ] **Step 1: Extract targets.smk**

```bash
sed -n '527,1002p' workflow/Snakefile > workflow/rules/targets.smk
```

- [ ] **Step 2: Verify targets.smk contains no rules**

```bash
grep -c "^rule " workflow/rules/targets.smk
# Expected: 0
```

- [ ] **Step 3: Verify function count**

```bash
grep -c "^def _.*_targets" workflow/rules/targets.smk
# Expected: 11 (10 rule-all target helpers + _manifest_dependency_targets)
```

---

### Task 3: Modify `workflow/Snakefile`

**Files:**
- Modify: `workflow/Snakefile`

Replace the extracted blocks with two `include:` statements. The assay policy includes (lines 448–450) stay and are placed BEFORE the metadata.smk include.

- [ ] **Step 1: Build the new Snakefile**

The new Snakefile structure (line by line):

```
Lines 1-140:   Keep (imports, config, sample loading)
               # End of line 140 is the PEAK_SAMPLE_IDS definition

Lines 448-450: Keep (assay policy includes — chipseq.smk, cuttag.smk, atac.smk)
               # These must precede metadata.smk because dispatch functions
               # in metadata.smk call assay policy functions.

New line:      include: "rules/metadata.smk"
               # Replaces lines 142-442 + 452-525 (gap at 444-450 is assay policy)

New line:      include: "rules/targets.smk"
               # Replaces lines 527-1002

Lines 1005-1036: Keep (rule all + executable rule includes)
```

- [ ] **Step 2: Construct by removing extracted blocks**

The simplest safe approach: comment out the extracted blocks, add the includes at the right positions.

Script:
```bash
# Build new Snakefile:
# 1. Lines 1-140 (unchanged)
sed -n '1,140p' workflow/Snakefile > /tmp/new_snakefile

# 2. Assay policy includes (lines 448-450)
echo "" >> /tmp/new_snakefile
sed -n '448,450p' workflow/Snakefile >> /tmp/new_snakefile

# 3. metadata.smk include
echo "" >> /tmp/new_snakefile
echo 'include: "rules/metadata.smk"' >> /tmp/new_snakefile

# 4. targets.smk include
echo "" >> /tmp/new_snakefile
echo 'include: "rules/targets.smk"' >> /tmp/new_snakefile

# 5. Lines 1004-1037 (rule all comment through end)
echo "" >> /tmp/new_snakefile
sed -n '1004,1037p' workflow/Snakefile >> /tmp/new_snakefile

# Replace
mv /tmp/new_snakefile workflow/Snakefile
```

- [ ] **Step 3: Verify Snakefile structure**

```bash
grep -n "^include:" workflow/Snakefile
# Expected order:
#   include: "rules/chipseq.smk"
#   include: "rules/cuttag.smk"
#   include: "rules/atac.smk"
#   include: "rules/metadata.smk"
#   include: "rules/targets.smk"
#   include: "rules/common.smk"
#   include: "rules/peaks.smk"
#   ...
```

- [ ] **Step 4: Verify rule all is first rule-bearing declaration**

```bash
grep -n "^rule " workflow/Snakefile | head -1
# Expected: "rule all:" (line number somewhere after the metadata/targets includes)
```

---

### Task 4: Run acceptance tests

- [ ] **Step 1: Validation tests**

```bash
python3 test/test_validation_stress.py
# Expected: 40/40 PASS
```

- [ ] **Step 2: No hardcoded paths**

```bash
python3 test/test_no_hardcoded_paths.py
# Expected: PASS
```

- [ ] **Step 3: Release readiness**

```bash
python3 test/test_stage28_release_readiness.py
# Expected: 11/11 PASS
```

- [ ] **Step 4: Smoke profiles (requires chipseq env)**

```bash
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage8_smoke_profiles.py
# Expected: 8/8 PASS
```

- [ ] **Step 5: MNase stress (requires chipseq env)**

```bash
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage39_mnase_stress.py
# Expected: all PASS
```

- [ ] **Step 6: Stage 12 / MultiQC stress (requires chipseq env)**

```bash
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage12_stress.py
# Expected: all PASS
```

- [ ] **Step 7: Manifest stress (requires chipseq env)**

```bash
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage25_manifest_stress.py
# Expected: all PASS
```

- [ ] **Step 8: Whitespace hygiene**

```bash
git diff --check
# Expected: clean (no output)
```

- [ ] **Step 9: Snakemake dry-run comparison**

```bash
/home/irenadler/miniconda3/envs/chipseq/bin/snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet
# Expected: identical output (no warnings, no errors with default config)
```

---

### Task 5: Commit

- [ ] **Step 1: Stage and commit**

```bash
git add workflow/rules/metadata.smk workflow/rules/targets.smk workflow/Snakefile
git add docs/superpowers/specs/2026-06-06-stage42-snakefile-extraction-design.md
git add docs/superpowers/plans/2026-06-06-stage42-snakefile-extraction.md
git commit -m "refactor: extract metadata.smk and targets.smk from Snakefile

Move derived metadata computation, dispatch functions, and QC gating
to workflow/rules/metadata.smk. Move ten target-builder functions to
workflow/rules/targets.smk. Snakefile reduced from ~1037 to ~200 lines.

Zero behavior change — Snakemake include: is literal text inclusion.
All variable names, function signatures, include order, and DAG
behavior are preserved identically."
```

---

### Task 6: Verify final file sizes

- [ ] **Step 1: Line counts**

```bash
wc -l workflow/Snakefile workflow/rules/metadata.smk workflow/rules/targets.smk
# Expected approximate:
#   ~200 workflow/Snakefile
#   ~380 workflow/rules/metadata.smk
#   ~475 workflow/rules/targets.smk
```
