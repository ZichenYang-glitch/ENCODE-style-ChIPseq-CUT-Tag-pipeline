# Stage 9: Release Readiness and User-Facing Polish — Design Spec

**Date:** 2026-05-19
**Status:** design review
**Scope:** README Limitations restructure, KNOWN_ISSUES cleanup, release checklist
**Excluded:** new pipeline features, config/schema changes, env refactoring

---

## 1. Goals

1. Restructure README Limitations into grouped, user-facing sections covering
   assay support, IDR constraints, QC gaps, and data limitations.
2. Update KNOWN_ISSUES.md: refresh the Scope Gap section (remove items now
   completed), add a concise Current Status snapshot, mark High Priority #3
   completed.
3. Add a release checklist to README Developer Notes.
4. No new pipeline features, no config/samples changes, no workflow/schema
   modifications.

---

## 2. README changes

### 2.1 Limitations section — full replacement

Replace the current Limitations section (7 flat bullet points) with:

```markdown
## Limitations

### Assay support

- **ChIP-seq** and **CUT&Tag** are supported. ATAC-seq is not included.

### TF ChIP-seq IDR

- Requires `chipseq` assay, `narrowPeak` mode, and **exactly 2 treatment
  biological replicates** per experiment.
- 3+ replicate IDR is not yet supported (automatic pairwise selection
  among 3+ replicates is not implemented).
- CUT&Tag IDR is not yet supported.

### Histone broad marks

- Broad-peak IDR is not yet supported. Histone experiments benefit from
  pooled QC summaries (Stage 6b) but do not produce IDR peak sets.

### QC gaps

- Cross-correlation metrics (NSC/RSC) are not yet implemented.
- FE/ppois bedGraph tracks are produced; BigWig conversion for those
  tracks is not yet implemented.
- `bamCoverage` RPGC normalization requires `--effectiveGenomeSize` and
  is not yet wired up in `tool_parameters`.
- MultiQC integration for experiment-level IDR and pooled QC summaries
  is not yet complete.

### Data and reproducibility

- No real public dataset is bundled with the pipeline.
- When duplicate removal is disabled (`remove_dup: "no"` or
  `auto` + `broad` peak mode), `final.bam` is the workflow's final
  post-filter BAM and may be a symlink to the filtered BAM. It should
  not be assumed to be duplicate-removed in all modes.
```

### 2.2 Sample sheet example — add optional-columns note

After the README sample sheet example table (~line 142), append:

```markdown
This is a compact replicate example. Optional replicate/control columns
(`experiment`, `condition`, `replicate`, `biological_replicate`,
`technical_replicate`, `role`, `control_sample`, `control_bam`)
are documented in `workflow/schemas/samples.schema.yaml`. Minimal
single-sample sheets may omit optional replicate and control columns.
```

### 2.3 Release checklist

Add under Developer Notes, after the Local execution section:

```markdown
### Release checklist

Before tagging a release, run these checks locally:

# 1. Validation
python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py

# 2. Default DAG check
snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet

# 3. Dry-run smoke profiles (7 profiles, <30 s)
python3 test/test_stage8_smoke_profiles.py

# 4. Tiny real execution (preprocessing + signal, <60 s)
python3 test/test_stage8b_tiny_execution.py

# 5. Repo hygiene
git status --short --untracked-files=all
# Expect: no results/, .snakemake/, *.fq, *.fq.gz, *.bam, *.bai, *.bw

# 6. CI status
# Check GitHub Actions for green on main / current branch.
```

---

## 3. KNOWN_ISSUES.md changes

### 3.1 Scope Gap — refresh

Replace the current Scope Gap section (lines 31-48) with:

```markdown
## Scope Gap

The current workflow is ENCODE-inspired, not a full ENCODE-compliant
ChIP-seq pipeline. It now provides:

- Single-sample preprocessing, MACS3 peak calling, CPM BigWig generation
- Single-sample QC: FRiP, library complexity (Picard + NRF/PBC),
  MACS3 FE/ppois bedGraph signal tracks
- Replicate-aware experiment model (Stage 4)
- TF ChIP-seq IDR for narrowPeak experiments with exactly 2 treatment
  biological replicates (Stage 5)
- Histone pooled QC summaries (Stage 6)
- CUT&Tag fragment-size QC and optional SEACR sidecar peaks (Stage 7)
- Test profiles, CI, and execution documentation (Stage 8)

Still missing:

- Cross-correlation metrics (NSC/RSC)
- BigWig conversion for FE/ppois bedGraph signal tracks
- MultiQC custom config for improved report layout
- Full preseq-style library complexity
- Browser-ready genome resource management
```

### 3.2 Current Status — replace existing section

Replace the existing `## Current Status` section (lines 6-29) with:

```markdown
## Current Status

- Core workflow, validation, replicate model, TF ChIP-seq IDR,
  histone/CUT&Tag support, smoke tests, and CI are implemented
  (Stages 2-8 complete).
- Remaining roadmap: NSC/RSC cross-correlation, FE/ppois BigWig
  conversion, optional public real-data validation, MultiQC tuning,
  GoPeaks support, advanced IDR strategies.
```

### 3.3 High Priority #3 — mark complete + add pending follow-up

Replace:

```markdown
3. Add small, reproducible smoke-test data.
   - The current bundled FASTQs are useful but large.
   - Add tiny FASTQ subsets and a small index for quick validation.
   - Include at least one treatment/control sample pair.
```

With:

```markdown
3. ✅ Add synthetic tiny FASTQ/index smoke execution (one ChIP-seq PE
   no-control real execution). (Stage 8b)
   - Stage 8a covers control paths by dry-run.
   - Pending: optional small public real-data treatment/control
     validation set.
```

---

## 4. Files changed

| File | Action |
|------|--------|
| `README.md` | Edit — Limitations restructure, sample-sheet note, release checklist |
| `KNOWN_ISSUES.md` | Edit — Scope Gap refresh, Current Status, High Priority #3 |

No config, samples, schemas, or workflow files are modified.

---

## 5. Verification plan

1. Visual check: README Limitations section covers all known gaps.
2. Visual check: KNOWN_ISSUES Scope Gap no longer lists completed features
   as "missing."
3. `python3 test/test_validation_stress.py` — 15/15 PASS.
4. `python3 test/test_stage8_smoke_profiles.py` — 7/7 PASS.
5. `git status --short --untracked-files=all` — only expected files.

---

## 6. Self-review notes

- No new pipeline features — documentation polish only.
- README does not grow excessively; Limitations goes from 7 flat bullets to
  5 grouped subsections (same total content, better structure).
- Sample sheet example stays minimal with a discoverability sentence.
- IDR limitations are explicit: exactly 2 bioreps, narrowPeak, chipseq only,
  no 3+ replicate auto-pairwise, no CUT&Tag IDR.
- No config/samples/schema changes.
- FE/ppois BigWig wording clearly separates existing CPM BigWig from missing
  FE/ppois BigWig conversion.
- Release checklist includes default DAG check.
- No TBDs or TODOs.
