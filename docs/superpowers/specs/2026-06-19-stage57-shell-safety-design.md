# Stage 57: Shell Safety Hardening — Design Spec

**Date:** 2026-06-19
**Status:** Implemented
**Author:** YangZiChen-glitch (Kaslana)
**Scope:** P0 safety — add `set -e -o pipefail` to all Snakemake shell blocks

## 1. Purpose

Ensure every Snakemake shell block fails fast on command errors and pipe
failures. Without this guard, a failed intermediate command in a pipeline
or multi-step shell block can go undetected, producing incomplete or
corrupt output files that Snakemake treats as valid.

## 2. Audit Results

| File | Total shell blocks | Already guarded | Missing | Fixed |
|------|-------------------|-----------------|---------|-------|
| `common.smk` | 11 | 1 | 10 | +10 |
| `qc.smk` | 25 | 19 | 6 | +6 |
| `report.smk` | 4 | 2 | 2 | +2 |
| `peaks.smk` | 2 | 2 | 0 | — |
| `idr.smk` | 7 | 7 | 0 | — |
| `idr_atac.smk` | 7 | 7 | 0 | — |
| `replicates.smk` | 3 | 3 | 0 | — |
| `mnase.smk` | 9 | 9 | 0 | — |
| **Total** | **68** | **50** | **18** | **+18** |

## 3. Changes

### 3.1 Shell block fixes (3 files, 18 blocks)

Add `set -e -o pipefail` as the first executable line in every missing
shell block.

- **common.smk:** fastqc, trim_galore, bowtie2_align, samtools_index_sorted,
  samtools_filter, samtools_index_filt, duplicate_handling, samtools_flagstat,
  samtools_flagstat_final, samtools_idxstats
- **qc.smk:** blacklist_filter_bam, blacklist_filter_peaks, peak_counts,
  pooled_experiment_qc_summary, qc_summary, stage3_qc_summary
- **report.smk:** tss_profile_qc, result_manifest

### 3.2 Enforcement test

**File:** `test/test_stage57_shell_safety.py`

Parses all `workflow/Snakefile` and `workflow/rules/*.smk` files, extracts
`shell:` blocks, and asserts each block starts with `set -e -o pipefail` or
has an exemption marker `# no pipefail: <reason>` with a non-empty reason.

Exits non-zero if any block fails the check.

## 4. Design decisions

| Decision | Rationale |
|----------|----------|
| `set -e -o pipefail` not `set -euo pipefail` | `-u` would break scripts referencing unset Snakemake variables |
| Every shell block, including trivial `touch` | Consistency prevents drift and simplifies enforcement |
| `script:` blocks excluded | Snakemake's `script:` runs Python directly, no bash |
| Exemption marker requires a reason | Prevents silent exemptions; documents intent |
| Conservative parser fails loud | If a shell block can't be parsed, error rather than silently skip |

## 5. Files

| File | Action |
|------|--------|
| `workflow/rules/common.smk` | +10 guards |
| `workflow/rules/qc.smk` | +6 guards |
| `workflow/rules/report.smk` | +2 guards |
| `test/test_stage57_shell_safety.py` | Create |
| `docs/superpowers/specs/2026-06-19-stage57-shell-safety-design.md` | Create |
| `docs/superpowers/plans/2026-06-19-stage57-shell-safety.md` | Create |

## 6. Non-goals

- No rule logic changes
- No output path, target, config, manifest, or reproducibility changes
- No `set -euo pipefail`
- No exemption markers added (none needed — all blocks get the guard)

## 7. Verification

```bash
python3 test/test_stage57_shell_safety.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```
