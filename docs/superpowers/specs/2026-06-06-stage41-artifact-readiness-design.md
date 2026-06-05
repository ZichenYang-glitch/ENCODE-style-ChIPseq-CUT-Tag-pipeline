# Stage 41: Artifact Readiness / Release Stabilization — Design Spec

**Date:** 2026-06-06
**Status:** implemented
**Scope:** Documentation-only — archive old report, create artifact roadmap, add developer checklist, output-contract consistency pass
**Excluded:** Artifact dataclass, AssayPolicy YAML, artifact_path(), global target resolver, DAG changes

---

## 1. Goals

1. Archive the shell-era reconstruction report with a clear header noting its outdated migration plan.
2. Create a staged artifact roadmap defining the path from current architecture toward artifact-oriented design.
3. Add a developer checklist for adding new outputs or assays.
4. Cross-reference `output-contract.md`, `make_manifest.py`, `Snakefile` target helpers, and `pipeline_done` for consistency.

---

## 2. Deliverables

### A. Archive reconstruction report

- **Move** `reconstruction-deep-research-report.md` → `docs/archive/reconstruction-deep-research-report-shell-era.md`
- Use `mv` (untracked file, not `git mv`).
- Prepend an archive header:

```
<!--
ARCHIVED — 2026-06-06

This document was written during the shell-pipeline era (scripts/chipseq.sh).
It remains useful as directional artifact-oriented research and the staged
roadmap concepts (artifact graph, assay policy, targets resolver) still inform
the project's long-term architecture thinking.

Its concrete migration plan is OUTDATED: the project is now a Snakemake
workflow with target-helper functions, assay-specific rule files, and
four supported assays (ChIP-seq, CUT&Tag, ATAC-seq, MNase-seq).

Do NOT treat this as a current implementation plan.
-->
```

- Do NOT touch `reconstruction-deep-research-report.md:Zone.Identifier`.

### B. Artifact roadmap

**Path:** `docs/architecture/artifact-roadmap.md`

Sections:
1. **Current baseline** — Snakemake + target-helper functions (`_base_targets()`, `_mnase_targets()`, etc.) + assay-specific rule files. Artifact-oriented design is future direction, not current implementation.
2. **Staged path** — Stages 41 through 50+ with goal, scope, non-goals, and trigger condition for each.
3. **Hard gates** — what must NOT be introduced before the trigger condition is met.

### C. Developer checklist

**Path:** `docs/developer/add-output-or-assay-checklist.md`

Flat numbered checklist (11 items) covering every file/subsystem that needs updating when a new output type or assay is added. Each item names the specific file and what to change.

### D. Output contract consistency pass

Cross-reference four sources:
- `docs/output-contract.md` tables
- `scripts/make_manifest.py` `output_type` vocabulary
- `workflow/Snakefile` `_*_targets()` expansions
- `workflow/rules/report.smk` `pipeline_done` inputs

Record result. Fix `docs/output-contract.md` only if a concrete mismatch is found. Do not modify `.smk`, `Snakefile`, or manifest code.

---

## 3. Files

| File | Action |
|:---|:---|
| `docs/archive/reconstruction-deep-research-report-shell-era.md` | Move from root + add header |
| `docs/architecture/artifact-roadmap.md` | Create |
| `docs/developer/add-output-or-assay-checklist.md` | Create |
| `docs/output-contract.md` | Possibly modify (consistency fixes only) |
| `docs/superpowers/specs/2026-06-06-stage41-artifact-readiness-design.md` | This file |
| `docs/superpowers/plans/2026-06-06-stage41-artifact-readiness.md` | Plan file |

---

## 4. Non-Goals (repeated for emphasis)

- No Artifact dataclass
- No AssayPolicy YAML
- No `artifact_path()`
- No global target resolver
- No `rule all` rewrite
- No DAG behavior changes
- No `.smk` file modifications (unless a concrete output contract inconsistency is found and approved)
- No test code changes (ditto)
