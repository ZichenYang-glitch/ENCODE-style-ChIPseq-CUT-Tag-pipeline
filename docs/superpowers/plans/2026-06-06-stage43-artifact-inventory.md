# Stage 43: Artifact Inventory — Implementation Plan

**Date:** 2026-06-06
**Status:** implemented
**Spec:** `docs/superpowers/specs/2026-06-06-stage43-artifact-inventory-design.md`

---

> **For agentic workers:** Use superpowers:executing-plans to implement.

**Goal:** Create a machine-readable YAML inventory of all 62 stable pipeline output types with metadata fields, plus a consistency test.

**Architecture:** Documentation-only. No code changes to Snakemake rules, targets, manifest, or DAG.

**Tech Stack:** YAML, Python stdlib, required PyYAML for parsing.

---

## Pre-requisite: Branch setup

```bash
git fetch origin
git switch -c stage43-artifact-inventory origin/main
```

---

### Task 1: Create the inventory file

**Files:**
- Create: `docs/architecture/artifact-inventory.yaml`

Structure: a single top-level list of artifact entries, ordered by scope (sample → experiment → IDR → project → reference). Each entry follows the schema from the spec.

Sections within the file:

1. **Sample-level outputs (36 entries):** shared preprocessing, peak-centric signal/QC, opt-in QC, TSS, CUT&Tag, and MNase sample-level outputs.

2. **Experiment-level outputs (21 entries):** standard pooled BAM/BAI, biological-replicate BAM/BAI, pooled signal/QC, MNase pooled outputs, and IDR outputs.

3. **Project-level outputs (4 entries):** multiqc_report, stage3_qc_summary, result_manifest, cross_correlation_summary.

4. **Reference outputs (1 entry):** tss_bed.

Stage 43 also updates `docs/output-contract.md` with two docs-only contract rows that were emitted by `make_manifest.py` but not documented: `pooled_final_bai` and `biorep_final_bai`.

File format:

```yaml
# Artifact Inventory — machine-readable reference for all pipeline outputs.
# Stage 43. Hand-maintained. Does NOT drive DAG behavior.
# See docs/architecture/artifact-roadmap.md for context.

artifacts:
  - id: final_bam
    description: Final filtered/dedup BAM for downstream use
    scope: sample
    level: per_sample
    assay_gate: all
    path_template: results/<sample>/02_align/<sample>.final.bam
    producing_rule: (common.smk)
    tool: bowtie2+samtools
    manifest_output_type: final_bam
    pipeline_done: false
    rule_all: false
    config_gate: null
    notes: downstream dependency; not a named pipeline_done input

  # ... 61 more entries ...
```

- [ ] **Step 1: Write the inventory with all 62 entries**

- [ ] **Step 2: Verify YAML syntax**

```bash
python3 -c "import yaml; yaml.safe_load(open('docs/architecture/artifact-inventory.yaml')); print('YAML OK')"
```

---

### Task 2: Create the consistency test

**Files:**
- Create: `test/test_stage43_artifact_inventory.py`

- [ ] **Step 1: Write the test file**

The test uses required `import yaml`. 10 checks:
1. YAML parses successfully; import/load failures fail the test
2. Inventory contains entries
3. All `id` values are unique
4. No `path_template` contains absolute workstation paths (`/home/`, `/data/`)
5. No entry mentions future-feature keywords as implemented in its `notes` field
6. Every manifest `output_type` has a matching inventory entry (AST extraction of `_add_row()` calls)
7. Every non-null inventory `manifest_output_type` is emitted by the manifest
8. All output-contract types are present in the inventory (section-based parsing)
9. All inventory IDs are present in the output-contract (bidirectional exact match)
10. Schema completeness: exactly 13 required fields per entry, `pipeline_done` and `rule_all` are booleans

Implementation details:
- Use `yaml.safe_load()` and fail if PyYAML is unavailable.
- Use Python `ast` to inspect `_add_row()` calls in `scripts/make_manifest.py`, including f-string `JoinedStr` output types such as `biorep<N>_final_bam`.
- Parse only the `docs/output-contract.md` section between `## Current Output Types` and `## Manifest Field Schema`.
- Assert both directions of output-contract coverage: output-contract to inventory and inventory to output-contract.

- [ ] **Step 2: Run the test**

```bash
python3 test/test_stage43_artifact_inventory.py
# Expected: 10/10 PASS
```

---

### Task 3: Run existing tests

- [ ] **Step 1: Release readiness**

```bash
python3 test/test_stage28_release_readiness.py
# Expected: 11/11 PASS
```

- [ ] **Step 2: No hardcoded paths**

```bash
python3 test/test_no_hardcoded_paths.py
# Expected: PASS
```

- [ ] **Step 3: Validation (no changes, regression check)**

```bash
python3 test/test_validation_stress.py
# Expected: 40/40 PASS
```

- [ ] **Step 4: Whitespace**

```bash
git diff --check
# Expected: clean
```

---

### Task 4: Commit

```bash
git add docs/architecture/artifact-inventory.yaml \
        test/test_stage43_artifact_inventory.py \
        docs/output-contract.md \
        test/test_no_hardcoded_paths.py \
        docs/superpowers/specs/2026-06-06-stage43-artifact-inventory-design.md \
        docs/superpowers/plans/2026-06-06-stage43-artifact-inventory.md

git commit -m "docs: add artifact inventory (Stage 43)

Create machine-readable YAML inventory of all 62 stable pipeline
output types with assay gating, path templates, and rule references.
Add consistency test cross-referencing make_manifest.py and
output-contract.md vocabularies. No DAG behavior changes."
```

---

### Task 5: Verify scope boundaries

- [ ] **Step 1: No accidental code changes**

```bash
git diff --stat origin/main
# Expected: docs/architecture/, test/test_stage43*.py, docs/superpowers/,
# docs/output-contract.md, and test/test_no_hardcoded_paths.py
# No workflow/*.smk, no workflow/Snakefile, no scripts/make_manifest.py
```
