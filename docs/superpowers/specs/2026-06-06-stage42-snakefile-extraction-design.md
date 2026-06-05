# Stage 42: Snakefile Responsibility Extraction — Design Spec

**Date:** 2026-06-06
**Status:** implemented
**Scope:** Create metadata.smk + targets.smk; reduce Snakefile from ~1037 to ~200 lines; zero behavior change
**Excluded:** Artifact dataclass, AssayPolicy YAML, artifact_path(), paths.smk, artifact inventory, global target resolver, path changes, DAG changes

---

## 1. Goals

1. Extract derived metadata computation and dispatch functions from `workflow/Snakefile` into `workflow/rules/metadata.smk`.
2. Extract ten target-builder functions from `workflow/Snakefile` into `workflow/rules/targets.smk`.
3. Replace extracted sections with `include:` statements.
4. Preserve every variable name, function signature, `include:` order, and DAG behavior identically.

---

## 2. Extraction Boundaries

### 2.1 Stay in workflow/Snakefile

| Section | Current lines | Content |
|:---|:---|---:|
| Header + imports | 1–24 | Docstring, `import os/sys/re`, `sys.path.insert`, imports from `scripts/` |
| Config loading | 26–83 | `configfile:`, `VALIDATED_CONFIG`, derived constants, `_normalize_genome()` |
| Sample loading | 85–140 | `SAMPLES`, validation calls, `SAMPLE_MAP`, all derived ID lists (`ALL_SAMPLE_IDS` through `PEAK_SAMPLE_IDS`) |
| Assay policy includes | 448–450 | `include: "rules/chipseq.smk"`, `cuttag.smk`, `atac.smk` |
| metadata.smk include | (new) | `include: "rules/metadata.smk"` |
| targets.smk include | (new) | `include: "rules/targets.smk"` |
| rule all | 1005–1023 | `rule all:` with ten target-helper calls + manifest |
| Executable rule includes | 1026–1036 | `include: "rules/common.smk"` through `report.smk` |

### 2.2 Move to workflow/rules/metadata.smk

| Section | Current lines | Content |
|:---|:---|---:|
| Stage 4a replicate metadata | 142–164 | `EXPERIMENT_IDS`, `SAMPLES_BY_EXPERIMENT`, `TREATMENT_SAMPLES_BY_EXPERIMENT`, `CONTROL_SAMPLES_BY_EXPERIMENT` |
| Stage 4b grouped outputs | 166–266 | `BIO_REP_GROUPS`, `_bioreps_for()`, `_biorep_expand_pairs()`, `MULTI_BIOREP_EXPERIMENTS`, `_resolve_experiment_controls()`, `POOLED_CONTROL_EXPERIMENTS`, `PEAK_MULTI_BIOREP_EXPERIMENTS`, `MNASE_MULTI_BIOREP_EXPERIMENTS` |
| Stage 5a/b IDR structures | 268–331 | `IDR_EXPERIMENTS`, IDR precomputed expansion lists |
| QC config + tool params | 333–381 | `QC_CONFIG`, `TOOL_PARAMS`, `_tool_param()`, `_filter_flags_arg()` |
| Blacklist/signal/TSS gating | 383–442 | `BLACKLIST_SAMPLES`, `get_genome_resource()`, `has_genome_resource()`, `TSS_SAMPLE_IDS`, `TSS_GENOMES`, `SIGNAL_BW_SAMPLE_IDS`, `SIGNAL_BW_EXPERIMENTS`, `_pooled_chrom_sizes()` |
| Dispatch wrappers | 452–496 | `get_remove_dup()`, `get_macs3_args()`, `get_extend_reads()` |
| MACS3 input resolution | 499–525 | `_macs3_inputs()` |

**Total:** ~380 lines. No Snakemake rules. Uses variables from Snakefile scope: `SAMPLE_MAP`, `OUTDIR`, `TREATMENT_SAMPLE_IDS`, etc.

### 2.3 Move to workflow/rules/targets.smk

| Section | Current lines | Content |
|:---|:---|---:|
| Target helpers | 528–1002 | `_base_targets()`, `_manifest_dependency_targets()`, `_MANIFEST_CONFIG_JSON` (lines 575–576), `_blacklist_targets()`, `_single_sample_qc_targets()`, `_signal_targets()`, `_cuttag_targets()`, `_advanced_qc_targets()`, `_tss_targets()`, `_replicate_targets()`, `_mnase_targets()`, `_idr_targets()` |

**Total:** ~475 lines. No Snakemake rules. Uses variables from Snakefile + metadata.smk scope: `OUTDIR`, `MULTIQC`, `TREATMENT_SAMPLE_IDS`, `PEAK_SAMPLE_IDS`, `MNASE_SAMPLE_IDS`, `QC_CONFIG`, etc.

---

## 3. Include Order

The critical dependency chain:

```
assay policy (chipseq.smk etc.) → dispatch functions (metadata.smk)
sample loading → metadata.smk → targets.smk → rule all
```

Final Snakefile include order (top to bottom):

```
1. imports (lines 1-24)
2. config loading (lines 26-83)
3. sample loading (lines 85-140)
4. include: "rules/chipseq.smk"      ← assay policy (must precede dispatch)
5. include: "rules/cuttag.smk"
6. include: "rules/atac.smk"
7. include: "rules/metadata.smk"     ← dispatch + derived structures
8. include: "rules/targets.smk"      ← target builders
9. rule all: (lines 1005-1023)
10. include: "rules/common.smk" ...   ← executable rules
```

Assay policy includes (items 4-6) are placed before metadata.smk because the dispatch functions inside metadata.smk call `get_remove_dup_chipseq()` etc. — functions defined in the assay policy files.

---

## 4. Files

| File | Action |
|:---|:---|
| `workflow/rules/metadata.smk` | **Create** — lines 142–525 of current Snakefile, minus assay policy includes (lines 448–450) |
| `workflow/rules/targets.smk` | **Create** — lines 527–1002 of current Snakefile |
| `workflow/Snakefile` | Modify — remove ~770 lines, add 2 `include:` statements |

---

## 5. Non-Goals (repeated for emphasis)

- No logic changes (identical byte-for-byte after Snakemake include resolution)
- No path changes
- No variable renaming
- No function signature changes
- No `rule all` rewrite
- No include order changes for executable rules
- No `paths.smk`
- No artifact inventory
- No Artifact dataclass
- No `artifact_path()`
- No AssayPolicy YAML
- No global target resolver
- No test assertion changes

---

## 6. Risks

| Risk | Likelihood | Mitigation |
|:---|:---|:---|
| Assay policy included after dispatch → NameError | Low | Include assay policy BEFORE metadata.smk (dispatch functions need the assay policy) |
| Variable used in targets.smk not visible | Low | All variables used by targets are defined in the Snakefile top or metadata.smk; both are included before targets.smk |
| `_MANIFEST_CONFIG_JSON` depends on `import json` at line 575 | Low | The `import json as _json` is inside the targets block (currently line 575); it moves WITH `_MANIFEST_CONFIG_JSON` to targets.smk |
| `rule all` no longer the first rule | Medium | `rule all:` stays in Snakefile at lines 1005-1023. The includes before it (metadata, targets, assay policy) define no rules. The first rule-bearing declaration is still `rule all`. Verified by Snakemake's "first rule = default target" behavior. |
| Cut at wrong line | Low | Exact line ranges verified against current `workflow/Snakefile`. Full test suite catches any missing variable references as NameError. |

---

## 7. Acceptance Tests

Every test must produce identical PASS/FAIL results before and after extraction:

1. `python3 test/test_validation_stress.py` — 40/40
2. `SNAKEMAKE=... python3 test/test_stage8_smoke_profiles.py` — 8/8 profiles
3. `SNAKEMAKE=... python3 test/test_stage39_mnase_stress.py` — all PASS
4. `SNAKEMAKE=... python3 test/test_stage12_stress.py` — all PASS
5. `SNAKEMAKE=... python3 test/test_stage25_manifest_stress.py` — all PASS
6. `python3 test/test_stage28_release_readiness.py` — 11/11
7. `python3 test/test_no_hardcoded_paths.py` — PASS
8. `git diff --check` — clean
9. Manual: `snakemake -n --quiet` with default config produces identical output
