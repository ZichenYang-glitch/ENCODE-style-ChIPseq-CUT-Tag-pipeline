# Stage 10a: Public Beta Release Packaging — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add release metadata (CITATION.cff, CHANGELOG.md), a beta-status sentence to README, and fix the LICENSE typo/copyright holder. Four files, metadata/docs only.

**Architecture:** Make all four file changes first, then run verification, show diff/status for review, and commit once after approval. No pipeline or schema changes.

**Tech Stack:** YAML (CITATION.cff), Markdown (CHANGELOG.md, README.md), plain text (LICENSE)

---

### Task 1: Fix LICENSE typo and copyright holder

**Files:**
- Modify: `LICENSE` (lines 1 and 3)

- [ ] **Step 1: Fix the first line typo `shencMIT License` → `MIT License`**

Current line 1:
```
shencMIT License
```

Replace with:
```
MIT License
```

- [ ] **Step 2: Update the copyright holder**

Current line 3:
```
Copyright (c) 2025 YangZiChen-glitch
```

Replace with:
```
Copyright (c) 2025 Zichen Yang
```

- [ ] **Step 3: Verify LICENSE after edit**

Expected file start:
```
MIT License

Copyright (c) 2025 Zichen Yang

Permission is hereby granted, free of charge, to any person obtaining a copy
...
```

---

### Task 2: Create CITATION.cff

**Files:**
- Create: `CITATION.cff`

- [ ] **Step 1: Write CITATION.cff**

```yaml
cff-version: 1.2.0
title: "ENCODE-style ChIP-seq and CUT&Tag Pipeline"
message: "If you use this software in your research, please cite it as below."
type: software
authors:
  - given-names: Zichen
    family-names: Yang
repository-code: "https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline"
abstract: >
  A Snakemake-based pipeline suite for ChIP-seq and CUT&Tag data analysis,
  supporting single-sample preprocessing, multi-replicate experiments with
  pooled outputs, single-sample QC, TF ChIP-seq IDR reproducibility analysis,
  histone pooled QC, and CUT&Tag-specific fragment-size QC with optional
  SEACR peak calling.
keywords:
  - ChIP-seq
  - CUT&Tag
  - Snakemake
  - peak-calling
  - IDR
  - bioinformatics
license: MIT
version: 0.1.0-beta
```

- [ ] **Step 2: Verify CITATION.cff structure**

Manual checks:
- `repository-code` is a single-line quoted YAML scalar (not split across two lines)
- No `date-released` field present
- `version: 0.1.0-beta` (no `v` prefix, string not float)
- No `orcid`, `affiliation`, `doi`, `identifiers`, or `url` fields
- `given-names: Zichen`, `family-names: Yang` matches LICENSE copyright holder
- All 6 keywords present

- [ ] **Step 3: Try cffconvert validation (optional, if available)**

```bash
cd /home/irenadler/workflow/chipseq && cffconvert --validate 2>&1 || echo "cffconvert not available — manual YAML check sufficient"
```

If `cffconvert` is not installed, skip — manual YAML structure review in Step 2 is sufficient for a valid CFF 1.2.0 document.

---

### Task 3: Create CHANGELOG.md

**Files:**
- Create: `CHANGELOG.md`

- [ ] **Step 1: Write CHANGELOG.md**

```markdown
# Changelog

All notable changes to this project will be documented in this file.
The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [v0.1.0-beta] - 2026-05-19

### Added
- Shared preprocessing pipeline: FastQC, Trim Galore, Bowtie2 alignment, MAPQ
  filtering, Picard duplicate handling, flagstat, idxstats, BigWig generation
- ChIP-seq / CUT&Tag assay-aware MACS3 peak calling with configurable parameters
- Optional control support via external BAM or FASTQ-based control sample rows
- Single-sample QC: FRiP, library complexity (Picard + NRF/PBC), MACS3 FE/ppois
  bedGraph signal tracks, per-sample QC summaries, and project-level aggregation
- Replicate-aware experiment model with technical replicate merging, biological
  replicate BAMs, and pooled treatment/control outputs
- TF ChIP-seq IDR: true-replicate IDR, pseudoreplicate-based self-IDR and
  pooled-IDR, reproducibility metrics, conservative/optimal peak sets
- Histone pooled QC summaries with target classification and peak-mode
  compatibility status
- CUT&Tag fragment-size QC and optional SEACR sidecar peak calling
- Configurable tool parameters for MACS3, bamCoverage, and IDR
- Schema-based config and sample sheet validation
- 7 test profiles covering major SE/PE/control/SEACR/IDR dispatch paths
- Tiny real-execution smoke test with synthetic data (no binary fixtures)
- GitHub Actions CI: fast validation + dry-run on PR/push, full execution on
  manual dispatch
- Streamlined environment setup with a full bioinformatics runtime environment
  and a lightweight CI environment
- Release checklist in README Developer Notes

### Documentation
- README with quick start, sample sheet reference, configuration guide,
  output structure, limitations, and developer notes
- KNOWN_ISSUES.md with full roadmap and current status
- Design specs and implementation plans under docs/superpowers/

### Known limitations
- Cross-correlation metrics (NSC/RSC) not yet implemented
- FE/ppois BigWig conversion not yet implemented
- IDR requires exactly 2 treatment biological replicates (no 3+ selection)
- CUT&Tag and broad-peak IDR not yet supported
- No bundled real public dataset
- ATAC-seq not supported

See README Limitations and KNOWN_ISSUES.md for full details.
```

- [ ] **Step 2: Verify CHANGELOG.md hygiene**

```bash
cd /home/irenadler/workflow/chipseq && grep -n 'Stage [0-9]' CHANGELOG.md
```

Expected: no output (no internal stage numbers in the changelog).

Also manually check:
- `[Unreleased]` section exists and is empty
- `[v0.1.0-beta] - 2026-05-19` date format correct
- Sections: `Added`, `Documentation`, `Known limitations`
- Environment bullet says: "Streamlined environment setup with a full bioinformatics runtime environment and a lightweight CI environment"
- No `Changed` / `Fixed` / `Removed` sections (first release, nothing to compare against)

---

### Task 4: Add beta status sentence to README

**Files:**
- Modify: `README.md` (insert after Overview paragraph, before `## Key Features` heading)

- [ ] **Step 1: Insert beta status sentence**

Insert the following paragraph after the line `Dependencies are managed with Conda.` and before `## Key Features`:

```markdown
This is a **v0.1.0-beta research workflow release**: functional and tested
with smoke and tiny-execution profiles, but not a fully ENCODE-compliant
production pipeline. See [Limitations](#limitations) for known gaps.
```

- [ ] **Step 2: Verify the README edit**

```bash
cd /home/irenadler/workflow/chipseq && grep -n 'v0.1.0-beta' README.md
```

Expected: one match, showing the beta status sentence with correct version string.

Also manually check:
- H1 title is unchanged: `# ENCODE-style ChIP-seq and CUT&Tag Pipeline`
- No new badge was added to the badge row
- The sentence appears between the Overview paragraph and `## Key Features`

---

### Task 5: Full verification (run after all four file changes are complete)

**Files:**
- Inspect: `CITATION.cff`, `CHANGELOG.md`, `README.md`, `LICENSE`

No new files created in this task — verification only.

- [ ] **Step 1: Confirm file scope — exactly 4 target files changed**

```bash
cd /home/irenadler/workflow/chipseq && git status --short
```

Expected: `LICENSE` and `README.md` modified; `CITATION.cff` and `CHANGELOG.md` untracked. No workflow, scripts, config, test, or schema files.

- [ ] **Step 2: Run validation stress tests — confirm zero regressions**

```bash
cd /home/irenadler/workflow/chipseq && python3 test/test_validation_stress.py
```

Expected: 15/15 PASS.

- [ ] **Step 3: Run smoke profiles — confirm zero regressions**

```bash
cd /home/irenadler/workflow/chipseq && python3 test/test_stage8_smoke_profiles.py
```

Expected: 7/7 PASS.

- [ ] **Step 4: Version consistency check**

```bash
cd /home/irenadler/workflow/chipseq && echo "=== CITATION.cff ===" && grep 'version' CITATION.cff && echo "=== CHANGELOG.md ===" && grep 'v0.1.0-beta' CHANGELOG.md && echo "=== README.md ===" && grep 'v0.1.0-beta' README.md
```

Expected:
- CITATION.cff: `version: 0.1.0-beta` (no `v` prefix)
- CHANGELOG.md: `## [v0.1.0-beta] - 2026-05-19`
- README.md: `This is a **v0.1.0-beta research workflow release**...`

- [ ] **Step 5: CITATION.cff field check**

```bash
cd /home/irenadler/workflow/chipseq && grep -c 'date-released' CITATION.cff
```

Expected: `0` (no `date-released` field present).

- [ ] **Step 6: LICENSE correctness**

```bash
cd /home/irenadler/workflow/chipseq && head -3 LICENSE
```

Expected output:
```
MIT License

Copyright (c) 2025 Zichen Yang
```

- [ ] **Step 7: No stray Stage numbers in CHANGELOG**

```bash
cd /home/irenadler/workflow/chipseq && grep -c 'Stage [0-9]' CHANGELOG.md
```

Expected: `0`.

- [ ] **Step 8: Show final diff for review**

```bash
cd /home/irenadler/workflow/chipseq && git diff && echo "=== Untracked files ===" && cat CITATION.cff && echo "---" && cat CHANGELOG.md
```

- [ ] **Step 9: Single commit (only after approval)**

```bash
cd /home/irenadler/workflow/chipseq && git add CITATION.cff CHANGELOG.md LICENSE README.md && git commit -m "$(cat <<'EOF'
feat: add v0.1.0-beta release packaging metadata

Add CITATION.cff (CFF 1.2.0), CHANGELOG.md (Keep a Changelog format),
beta status sentence in README Overview, and fix LICENSE header typo
plus copyright holder update.
EOF
)"
```

---

## Specification Reference

All design decisions are in `docs/superpowers/specs/2026-05-19-stage10a-release-packaging-design.md`. No deviation from the spec is planned.
