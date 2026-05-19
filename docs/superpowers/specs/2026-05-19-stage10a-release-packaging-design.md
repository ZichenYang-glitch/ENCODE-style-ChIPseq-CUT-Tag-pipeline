# Stage 10a: Public Beta Release Packaging — Design Spec

**Date:** 2026-05-19
**Status:** design review
**Scope:** metadata and documentation files only — CITATION.cff, CHANGELOG.md, README.md (1 sentence), LICENSE typo fix and copyright holder update
**Excluded:** new pipeline features, config/schema changes, env modifications, workflow/script/test changes, version tagging

---

## 1. Goals

1. Add `CITATION.cff` with project metadata for a first public beta release (v0.1.0-beta).
2. Create `CHANGELOG.md` in Keep a Changelog format with `[Unreleased]` and `[v0.1.0-beta]` entries.
3. Add a single beta-status sentence to the README Overview (no H1 change, no badge).
4. Fix the LICENSE typo (`shencMIT` → `MIT`) and update the copyright holder to `Zichen Yang`.
5. No pipeline logic, workflow, script, schema, config, or test changes.

---

## 2. CITATION.cff

New file at the repository root.

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

### Field decisions

- `cff-version: 1.2.0` — current CFF schema version; compatible with GitHub and Zenodo.
- No `date-released` — the release date will be added in Stage 10d when the tag is created.
- No `doi`, `identifiers`, `url` — no DOI has been minted yet.
- No `orcid` or `affiliation` — deferred until the author chooses to add them.
- `version: 0.1.0-beta` — string, not float, per CFF spec.
- `repository-code` — single-line quoted scalar pointing to the public GitHub repo.
- `abstract` — multi-line folded block scalar (`>`), readable for GitHub/Zenodo indexing.
- `license: MIT` — matches `LICENSE` file.
- `authors` — single entry, `given-names: Zichen`, `family-names: Yang`, consistent with the updated LICENSE copyright line.

---

## 3. CHANGELOG.md

New file at the repository root. Keep a Changelog format with `[Unreleased]` and `[v0.1.0-beta]`.

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

### Structure decisions

- `[Unreleased]` at top — empty; future PRs append entries here. At the next release, `[Unreleased]` content moves under a new version heading.
- `[v0.1.0-beta] - 2026-05-19` — release date matches the day the entry is committed.
- Sections: `Added`, `Documentation`, `Known limitations`. No `Changed`/`Fixed`/`Removed` yet (first release).
- `Known limitations` is a non-standard Keep a Changelog section but appropriate for a beta release. Honest and discoverable.
- No internal Stage numbers anywhere in the changelog text.
- Environment bullet uses the approved wording: "Streamlined environment setup with a full bioinformatics runtime environment and a lightweight CI environment."

---

## 4. README.md

### 4.1 Beta status sentence

After the Overview paragraph ending with "Dependencies are managed with Conda." (line 11 of current README), insert one paragraph:

```markdown
This is a **v0.1.0-beta research workflow release**: functional and tested
with smoke and tiny-execution profiles, but not a fully ENCODE-compliant
production pipeline. See [Limitations](#limitations) for known gaps.
```

The sentence is between the Overview paragraph and the `## Key Features` heading. No other changes to README.

### 4.2 What stays unchanged

- H1 title: `# ENCODE-style ChIP-seq and CUT&Tag Pipeline` — unchanged.
- Badge row — no beta badge added.
- All other sections (Key Features, Quick Start, Sample Sheet, Configuration, Workflow Capabilities, Limitations, Output Structure, Developer Notes, License) — untouched.

---

## 5. LICENSE

Current state (line 1):

```
shencMIT License
```

Replace with:

```
MIT License
```

Current copyright line:

```
Copyright (c) 2025 YangZiChen-glitch
```

Replace with:

```
Copyright (c) 2025 Zichen Yang
```

Only these two lines change. The rest of the MIT license text remains verbatim.

The updated copyright holder `Zichen Yang` matches the citation-facing author name in `CITATION.cff`.

---

## 6. Files changed

| File | Action | Lines changed |
|------|--------|---------------|
| `CITATION.cff` | Create | ~22 lines |
| `CHANGELOG.md` | Create | ~56 lines |
| `README.md` | Edit — insert 3 lines after Overview | ~3 lines inserted |
| `LICENSE` | Edit — fix 2 lines | 2 lines changed |

Total: 4 files. No workflow, scripts, config, test, schema, or docs/superpowers changes.

---

## 7. Verification plan

1. **CITATION.cff validity:** `cffconvert --validate` if `cffconvert` is available in the environment; otherwise manual YAML structure check.
2. **File scope:** `git diff --stat` — exactly 4 files: `CITATION.cff`, `CHANGELOG.md`, `README.md`, `LICENSE`.
3. **No pipeline regressions:**
   - `python3 test/test_validation_stress.py` — 15/15 PASS.
   - `python3 test/test_stage8_smoke_profiles.py` — 7/7 PASS.
4. **CHANGELOG hygiene:** Grep for internal stage numbers (`Stage 2`, `Stage 3`, etc.) — must be zero matches.
5. **Version consistency:**
   - CITATION.cff: `version: 0.1.0-beta` (no `date-released`).
   - CHANGELOG.md: `[v0.1.0-beta]` heading.
   - README.md: `v0.1.0-beta` in the Overview sentence.
6. **LICENSE:** First line is `MIT License` (not `shencMIT`). Copyright line is `Copyright (c) 2025 Zichen Yang`.

---

## 8. Self-review

- **Placeholder scan:** No TBDs, TODOs, or incomplete sections.
- **Internal consistency:** Version string `v0.1.0-beta` / `0.1.0-beta` is referenced consistently across all four files. README uses `v0.1.0-beta` with the `v` prefix (user-facing prose); CITATION.cff uses `0.1.0-beta` without the prefix (CFF spec convention).
- **Scope check:** Four files, metadata/docs only. No pipeline, config, schema, or test modifications. Appropriate for a single implementation session.
- **Ambiguity check:** CITATION.cff YAML shape is explicit — no optional fields that might be interpreted differently. CHANGELOG structure is prescribed down to section headings. README insertion point is precise (after the Overview paragraph, before `## Key Features`).
- **No feature creep:** Release notes file (`docs/release-notes/v0.1.0-beta.md`) is explicitly excluded. Beta badge is explicitly excluded. H1 title change is explicitly excluded. `date-released` is explicitly excluded from CITATION.cff. Tagging is deferred to Stage 10d.
