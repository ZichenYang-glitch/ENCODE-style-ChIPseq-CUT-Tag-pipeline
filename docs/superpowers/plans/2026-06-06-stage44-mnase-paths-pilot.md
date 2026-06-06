# Stage 44: MNase paths.smk Pilot — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Create a zero-behavior-change MNase-only path helper file (`paths.smk`) and replace inline path f-strings in three consumer files, with no manifest/contract/artifact changes.

**Architecture:** `paths.smk` is an include-only file with 7 thin helper functions referencing global `OUTDIR`. It is included in `workflow/Snakefile` before its three consumers (`targets.smk`, `mnase.smk`, `report.smk`). Every helper substitution produces byte-identical path strings — the DAG does not change.

**Tech Stack:** Snakemake (include-based namespace), Python 3 stdlib f-strings. No packages, no imports, no external dependencies.

---

### Task 1: Create `workflow/rules/paths.smk`

**Files:**
- Create: `workflow/rules/paths.smk`

- [ ] **Step 1: Write the 7-helper file**

```python
# Stage 44: MNase path helpers — zero-behavior-change pilot.
# Included into Snakefile namespace; references global OUTDIR.
# No rule declarations. No imports. No side effects.


def mnase_fragment_bam(sample, fragment_class):
    """Sample-level MNase fragment BAM path.

    Args:
        sample: Snakemake wildcard string (e.g. "{sample}") or concrete id.
        fragment_class: "mono", "sub", or "di".
    """
    return f"{OUTDIR}/{sample}/03_fragments/{sample}.{fragment_class}.bam"


def mnase_fragment_bai(sample, fragment_class):
    """Sample-level MNase fragment BAM index path."""
    return f"{OUTDIR}/{sample}/03_fragments/{sample}.{fragment_class}.bam.bai"


def mnase_signal_bw(sample, kind):
    """Sample-level MNase BigWig path.

    Args:
        sample: Snakemake wildcard string or concrete id.
        kind: "dyad" or "mono".
    """
    return f"{OUTDIR}/{sample}/04_signal/{sample}.{kind}.CPM.bw"


def mnase_qc_summary_tsv(sample):
    """Sample-level MNase QC summary TSV path."""
    return f"{OUTDIR}/{sample}/01_qc/{sample}.mnase_qc_summary.tsv"


def mnase_pooled_fragment_bam(experiment, fragment_class):
    """Experiment-level pooled MNase fragment BAM path.

    Args:
        experiment: Snakemake wildcard string (e.g. "{experiment}") or concrete id.
        fragment_class: "mono", "sub", or "di".
    """
    return f"{OUTDIR}/experiments/{experiment}/03_fragments/{experiment}.pooled.{fragment_class}.bam"


def mnase_pooled_fragment_bai(experiment, fragment_class):
    """Experiment-level pooled MNase fragment BAM index path."""
    return f"{OUTDIR}/experiments/{experiment}/03_fragments/{experiment}.pooled.{fragment_class}.bam.bai"


def mnase_pooled_signal_bw(experiment, kind):
    """Experiment-level pooled MNase BigWig path.

    Args:
        experiment: Snakemake wildcard string or concrete id.
        kind: "dyad" or "mono".
    """
    return f"{OUTDIR}/experiments/{experiment}/04_signal/{experiment}.pooled.{kind}.CPM.bw"
```

- [ ] **Step 2: Verify no rules in paths.smk**

```bash
rg -n '^rule ' workflow/rules/paths.smk
```
Expected: no output (exit 1).

---

### Task 2: Include `paths.smk` in `workflow/Snakefile`

**Files:**
- Modify: `workflow/Snakefile`

- [ ] **Step 1: Insert include after assay policy includes, before metadata.smk**

Insert `include: "rules/paths.smk"` between lines 149 and 154 (after `atac.smk`, before the "Derived metadata" section header and `metadata.smk`).

Before:
```python
include: "rules/chipseq.smk"
include: "rules/cuttag.smk"
include: "rules/atac.smk"

# ---------------------------------------------------------------------------
# 5. Derived metadata, QC gating, dispatch wrappers
# ---------------------------------------------------------------------------
include: "rules/metadata.smk"
```

After:
```python
include: "rules/chipseq.smk"
include: "rules/cuttag.smk"
include: "rules/atac.smk"
include: "rules/paths.smk"

# ---------------------------------------------------------------------------
# 5. Derived metadata, QC gating, dispatch wrappers
# ---------------------------------------------------------------------------
include: "rules/metadata.smk"
```

- [ ] **Step 2: Add section header comment before the new include**

Insert these two comment lines above the new include line:
```python
# ---------------------------------------------------------------------------
# 4b. MNase path helpers (no rules — define-only file)
# ---------------------------------------------------------------------------
```

The final block should read:
```python
include: "rules/cuttag.smk"
include: "rules/atac.smk"

# ---------------------------------------------------------------------------
# 4b. MNase path helpers (no rules — define-only file)
# ---------------------------------------------------------------------------
include: "rules/paths.smk"

# ---------------------------------------------------------------------------
# 5. Derived metadata, QC gating, dispatch wrappers
# ---------------------------------------------------------------------------
include: "rules/metadata.smk"
```

- [ ] **Step 3: Verify rule all is still the first rule-bearing declaration**

```bash
rg -n '^rule ' workflow/Snakefile workflow/rules/paths.smk \
  workflow/rules/metadata.smk workflow/rules/targets.smk
```
Expected: First match is `workflow/Snakefile:<line>:rule all:`. No matches in `workflow/rules/paths.smk`.

---

### Task 3: Replace MNase paths in `workflow/rules/mnase.smk`

**Files:**
- Modify: `workflow/rules/mnase.smk`

Replacement table:

| Line | Current | Replacement |
|:---|:---|:---|
| 130 | `bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam",` | `bam = mnase_fragment_bam("{sample}", "mono"),` |
| 131 | `bai = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam.bai",` | `bai = mnase_fragment_bai("{sample}", "mono"),` |
| 163 | `bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.sub.bam",` | `bam = mnase_fragment_bam("{sample}", "sub"),` |
| 164 | `bai = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.sub.bam.bai",` | `bai = mnase_fragment_bai("{sample}", "sub"),` |
| 196 | `bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.di.bam",` | `bam = mnase_fragment_bam("{sample}", "di"),` |
| 197 | `bai = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.di.bam.bai",` | `bai = mnase_fragment_bai("{sample}", "di"),` |
| 229 | `bw = f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.dyad.CPM.bw",` | `bw = mnase_signal_bw("{sample}", "dyad"),` |
| 272 | `bw = f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.mono.CPM.bw",` | `bw = mnase_signal_bw("{sample}", "mono"),` |
| 274 | `bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam",` | `bam = mnase_fragment_bam("{sample}", "mono"),` |
| 275 | `bai = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam.bai",` | `bai = mnase_fragment_bai("{sample}", "mono"),` |
| 309 | `bam = f"{OUTDIR}/experiments/{{experiment}}/03_fragments/{{experiment}}.pooled.mono.bam",` | `bam = mnase_pooled_fragment_bam("{experiment}", "mono"),` |
| 310 | `bai = f"{OUTDIR}/experiments/{{experiment}}/03_fragments/{{experiment}}.pooled.mono.bam.bai",` | `bai = mnase_pooled_fragment_bai("{experiment}", "mono"),` |
| 342 | `bw = f"{OUTDIR}/experiments/{{experiment}}/04_signal/{{experiment}}.pooled.dyad.CPM.bw",` | `bw = mnase_pooled_signal_bw("{experiment}", "dyad"),` |
| 385 | `bw = f"{OUTDIR}/experiments/{{experiment}}/04_signal/{{experiment}}.pooled.mono.CPM.bw",` | `bw = mnase_pooled_signal_bw("{experiment}", "mono"),` |
| 387 | `bam = f"{OUTDIR}/experiments/{{experiment}}/03_fragments/{{experiment}}.pooled.mono.bam",` | `bam = mnase_pooled_fragment_bam("{experiment}", "mono"),` |
| 388 | `bai = f"{OUTDIR}/experiments/{{experiment}}/03_fragments/{{experiment}}.pooled.mono.bam.bai",` | `bai = mnase_pooled_fragment_bai("{experiment}", "mono"),` |
| 422 | `f"{OUTDIR}/{{sample}}/01_qc/{{sample}}.mnase_qc_summary.tsv",` | `mnase_qc_summary_tsv("{sample}"),` |
| 424 | `sub_bam  = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.sub.bam",` | `sub_bam  = mnase_fragment_bam("{sample}", "sub"),` |
| 425 | `mono_bam = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.mono.bam",` | `mono_bam = mnase_fragment_bam("{sample}", "mono"),` |
| 426 | `di_bam   = f"{OUTDIR}/{{sample}}/03_fragments/{{sample}}.di.bam",` | `di_bam   = mnase_fragment_bam("{sample}", "di"),` |
| 427 | `dyad_bw  = f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.dyad.CPM.bw",` | `dyad_bw  = mnase_signal_bw("{sample}", "dyad"),` |
| 428 | `mono_bw  = f"{OUTDIR}/{{sample}}/04_signal/{{sample}}.mono.CPM.bw",` | `mono_bw  = mnase_signal_bw("{sample}", "mono"),` |

**22 lines total across 9 rules.** Not changed: `final.bam`/`final.bai` inputs, log file paths, `insert_size_metrics` param path.

- [ ] **Step 1: Apply all 22 replacements**

The replacements are listed in the table above. Apply them one rule at a time.

---

### Task 4: Replace `expand()` calls in `workflow/rules/targets.smk`

**Files:**
- Modify: `workflow/rules/targets.smk` (lines 338-408, `_mnase_targets()`)

- [ ] **Step 1: Replace the 9 sample-level expand() calls**

Current (lines 340-386):
```python
    if MNASE_SAMPLE_IDS:
        # Stage 39: mono BAM, dyad BW, mono occupancy BW
        targets += expand(
            "{outdir}/{sample}/03_fragments/{sample}.mono.bam",
            outdir=OUTDIR,
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/03_fragments/{sample}.mono.bam.bai",
            outdir=OUTDIR,
            sample=MNASE_SAMPLE_IDS,
        )
        # Stage 40: sub-nucleosome and di-nucleosome BAMs
        targets += expand(
            "{outdir}/{sample}/03_fragments/{sample}.sub.bam",
            outdir=OUTDIR,
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/03_fragments/{sample}.sub.bam.bai",
            outdir=OUTDIR,
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/03_fragments/{sample}.di.bam",
            outdir=OUTDIR,
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/03_fragments/{sample}.di.bam.bai",
            outdir=OUTDIR,
            sample=MNASE_SAMPLE_IDS,
        )
        # Stage 40: MNase QC summary
        targets += expand(
            "{outdir}/{sample}/01_qc/{sample}.mnase_qc_summary.tsv",
            outdir=OUTDIR,
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/04_signal/{sample}.dyad.CPM.bw",
            outdir=OUTDIR,
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            "{outdir}/{sample}/04_signal/{sample}.mono.CPM.bw",
            outdir=OUTDIR,
            sample=MNASE_SAMPLE_IDS,
        )
```

Replacement:
```python
    if MNASE_SAMPLE_IDS:
        # Stage 39: mono BAM, dyad BW, mono occupancy BW
        targets += expand(
            mnase_fragment_bam("{sample}", "mono"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_fragment_bai("{sample}", "mono"),
            sample=MNASE_SAMPLE_IDS,
        )
        # Stage 40: sub-nucleosome and di-nucleosome BAMs
        targets += expand(
            mnase_fragment_bam("{sample}", "sub"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_fragment_bai("{sample}", "sub"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_fragment_bam("{sample}", "di"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_fragment_bai("{sample}", "di"),
            sample=MNASE_SAMPLE_IDS,
        )
        # Stage 40: MNase QC summary
        targets += expand(
            mnase_qc_summary_tsv("{sample}"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_signal_bw("{sample}", "dyad"),
            sample=MNASE_SAMPLE_IDS,
        )
        targets += expand(
            mnase_signal_bw("{sample}", "mono"),
            sample=MNASE_SAMPLE_IDS,
        )
```

- [ ] **Step 2: Replace the 4 experiment-level expand() calls**

Current (lines 388-408):
```python
    # Pooled MNase outputs (>=2 biorep MNase experiments)
    if STAGE4B and MNASE_MULTI_BIOREP_EXPERIMENTS:
        targets += expand(
            "{outdir}/experiments/{experiment}/03_fragments/{experiment}.pooled.mono.bam",
            outdir=OUTDIR,
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/03_fragments/{experiment}.pooled.mono.bam.bai",
            outdir=OUTDIR,
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/04_signal/{experiment}.pooled.dyad.CPM.bw",
            outdir=OUTDIR,
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            "{outdir}/experiments/{experiment}/04_signal/{experiment}.pooled.mono.CPM.bw",
            outdir=OUTDIR,
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
```

Replacement:
```python
    # Pooled MNase outputs (>=2 biorep MNase experiments)
    if STAGE4B and MNASE_MULTI_BIOREP_EXPERIMENTS:
        targets += expand(
            mnase_pooled_fragment_bam("{experiment}", "mono"),
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            mnase_pooled_fragment_bai("{experiment}", "mono"),
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            mnase_pooled_signal_bw("{experiment}", "dyad"),
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
        targets += expand(
            mnase_pooled_signal_bw("{experiment}", "mono"),
            experiment=MNASE_MULTI_BIOREP_EXPERIMENTS,
        )
```

Key change: `outdir=OUTDIR` is removed from every `expand()` call because the helper already embeds global `OUTDIR`.

---

### Task 5: Replace MNase lambda inputs in `workflow/rules/report.smk`

**Files:**
- Modify: `workflow/rules/report.smk` (lines 61-95, `pipeline_done` rule)

- [ ] **Step 1: Replace the 9 MNase conditional lambda inputs**

Current (lines 61-95):
```python
        mnase_mono_bam = lambda wc: (
            f"{OUTDIR}/{wc.sample}/03_fragments/{wc.sample}.mono.bam"
            if _is_mnase(wc) else []
        ),
        mnase_mono_bai = lambda wc: (
            f"{OUTDIR}/{wc.sample}/03_fragments/{wc.sample}.mono.bam.bai"
            if _is_mnase(wc) else []
        ),
        mnase_dyad_bw = lambda wc: (
            f"{OUTDIR}/{wc.sample}/04_signal/{wc.sample}.dyad.CPM.bw"
            if _is_mnase(wc) else []
        ),
        mnase_mono_bw = lambda wc: (
            f"{OUTDIR}/{wc.sample}/04_signal/{wc.sample}.mono.CPM.bw"
            if _is_mnase(wc) else []
        ),
        mnase_sub_bam = lambda wc: (
            f"{OUTDIR}/{wc.sample}/03_fragments/{wc.sample}.sub.bam"
            if _is_mnase(wc) else []
        ),
        mnase_sub_bai = lambda wc: (
            f"{OUTDIR}/{wc.sample}/03_fragments/{wc.sample}.sub.bam.bai"
            if _is_mnase(wc) else []
        ),
        mnase_di_bam = lambda wc: (
            f"{OUTDIR}/{wc.sample}/03_fragments/{wc.sample}.di.bam"
            if _is_mnase(wc) else []
        ),
        mnase_di_bai = lambda wc: (
            f"{OUTDIR}/{wc.sample}/03_fragments/{wc.sample}.di.bam.bai"
            if _is_mnase(wc) else []
        ),
        mnase_qc_summ = lambda wc: (
            f"{OUTDIR}/{wc.sample}/01_qc/{wc.sample}.mnase_qc_summary.tsv"
            if _is_mnase(wc) else []
        ),
```

Replacement:
```python
        mnase_mono_bam = lambda wc: (
            mnase_fragment_bam(wc.sample, "mono") if _is_mnase(wc) else []
        ),
        mnase_mono_bai = lambda wc: (
            mnase_fragment_bai(wc.sample, "mono") if _is_mnase(wc) else []
        ),
        mnase_dyad_bw = lambda wc: (
            mnase_signal_bw(wc.sample, "dyad") if _is_mnase(wc) else []
        ),
        mnase_mono_bw = lambda wc: (
            mnase_signal_bw(wc.sample, "mono") if _is_mnase(wc) else []
        ),
        mnase_sub_bam = lambda wc: (
            mnase_fragment_bam(wc.sample, "sub") if _is_mnase(wc) else []
        ),
        mnase_sub_bai = lambda wc: (
            mnase_fragment_bai(wc.sample, "sub") if _is_mnase(wc) else []
        ),
        mnase_di_bam = lambda wc: (
            mnase_fragment_bam(wc.sample, "di") if _is_mnase(wc) else []
        ),
        mnase_di_bai = lambda wc: (
            mnase_fragment_bai(wc.sample, "di") if _is_mnase(wc) else []
        ),
        mnase_qc_summ = lambda wc: (
            mnase_qc_summary_tsv(wc.sample) if _is_mnase(wc) else []
        ),
```

Key invariant: `helper(wc.sample, ...) if _is_mnase(wc) else []` — no list wrapping on the helper result, exactly matching the current lambda return shape.

---

### Task 6: Run hard acceptance tests

- [ ] **Step 1: Verify paths.smk has zero rule declarations**

```bash
rg -n '^rule ' workflow/rules/paths.smk
```
Expected: no output (exit 1).

- [ ] **Step 2: Verify rule all remains first rule-bearing declaration**

```bash
rg -n '^rule ' workflow/Snakefile workflow/rules/paths.smk \
  workflow/rules/metadata.smk workflow/rules/targets.smk
```
Expected: First match is `workflow/Snakefile:<line>:rule all:`. No matches in `workflow/rules/paths.smk`.

- [ ] **Step 3: Run existing test suite**

```bash
python3 test/test_stage43_artifact_inventory.py
```
Expected: 10/10 PASS.

```bash
python3 test/test_validation_stress.py
```
Expected: 40/40 PASS.

```bash
python3 test/test_stage28_release_readiness.py
```
Expected: 11/11 PASS.

```bash
python3 test/test_no_hardcoded_paths.py
```
Expected: PASS.

```bash
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage39_mnase_stress.py
```
Expected: all tests pass.

```bash
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage8_smoke_profiles.py
```
Expected: all tests pass.

```bash
SNAKEMAKE=/home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  python3 test/test_stage25_manifest_stress.py
```
Expected: all tests pass.

- [ ] **Step 4: Check git diff for whitespace issues**

```bash
git diff --check
```
Expected: no output (clean).

---

### Task 7: Commit

- [ ] **Step 1: Stage and commit**

```bash
git add workflow/rules/paths.smk workflow/Snakefile \
  workflow/rules/mnase.smk workflow/rules/targets.smk \
  workflow/rules/report.smk \
  docs/superpowers/specs/2026-06-06-stage44-mnase-paths-pilot-design.md \
  docs/superpowers/plans/2026-06-06-stage44-mnase-paths-pilot.md

git commit -m "$(cat <<'EOF'
feat: add MNase path helpers pilot (Stage 44)

Zero-behavior-change: extracts 7 thin MNase path helpers into
workflow/rules/paths.smk (no rule declarations), replacing inline
f-string paths in mnase.smk (9 rules), targets.smk (_mnase_targets),
and report.smk (pipeline_done MNase inputs).
EOF
)"
```

- [ ] **Step 2: Verify commit landed**

```bash
git log --oneline -1
git status
```

---

### Optional Manual Check: Dry-run DAG equivalence

This requires capturing a "before" snapshot before implementation, or running on a separate checkout.

```bash
# After implementation:
XDG_CACHE_HOME=/tmp/chipseq-snakemake-cache \
  /home/irenadler/miniconda3/envs/chipseq/bin/snakemake \
  -s workflow/Snakefile --configfile config/config.yaml -n --quiet \
  2>&1 | sort > /tmp/dag_after.txt

# Compare against saved "before" (if available)
diff /tmp/dag_before.txt /tmp/dag_after.txt
```
Expected: no differences.
