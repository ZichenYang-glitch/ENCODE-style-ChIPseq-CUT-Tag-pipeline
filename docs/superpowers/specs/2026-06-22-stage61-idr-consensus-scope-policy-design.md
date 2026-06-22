# Stage 61: IDR/Consensus Scope Policy Alignment — Design Spec

**Date:** 2026-06-22
**Status:** Implemented
**Author:** YangZiChen-glitch (Kaslana)
**Scope:** Policy/documentation alignment — no runtime changes

## 1. Purpose

Stage 53 defined the reproducibility policy. After implementing ATAC narrow
IDR (Stage 55) and revisiting the Artifact boundary (Stage 56 ADR), the policy
docs had accumulated minor inconsistencies. CUT&Tag narrow IDR was classified
as "experimental" when it should be "supported opt-in," and the final-output
table was missing a CUT&Tag narrow IDR-conditional row that ATAC narrow
already had.

Stage 61 aligns the policy docs, config schema, and tests to the following
agreed matrix without modifying any Snakemake rules, targets, or runtime
behavior.

## 2. Policy Matrix (Canonical)

| # | Assay | Peak mode | Caller | IDR classification | Final when IDR enabled |
|---|-------|-----------|--------|-------------------|----------------------|
| 1 | chipseq | narrow | MACS3 | Production-supported (legacy stage5) | IDR final, consensus secondary |
| 2 | chipseq | broad | MACS3 | Experimental opt-in | Consensus final; experimental IDR never silently replaces |
| 3 | cuttag | narrow | MACS3 | Supported opt-in | IDR final when `cuttag_narrow: true`; consensus otherwise |
| 4 | cuttag | broad | MACS3 | Experimental opt-in | Consensus final; experimental IDR never silently replaces |
| 5 | cuttag | — | SEACR | No IDR | Consensus only |
| 6 | atac | narrow | MACS3 | Production-supported | IDR final when `atac_narrow: true`; consensus otherwise |
| — | mnase | — | — | No IDR | Outside peak IDR policy |

## 3. Files changed

| File | Change |
|------|--------|
| `docs/reproducibility-policy.md` | Update §2.3 IDR classification, §3 row 3, §4.4 final table |
| `docs/superpowers/specs/2026-06-18-stage53-reproducibility-policy-design.md` | Update §2.2 row 1/3/6, §4.4 final table |
| `workflow/schemas/config.schema.yaml` | Fix `cuttag_narrow` description: "opt-in, experimental" → "opt-in, supported" |
| `test/test_stage53_reproducibility_policy_contract.py` | Lock canonical IDR/consensus policy matrix |
| `docs/superpowers/specs/2026-06-22-stage61-idr-consensus-scope-policy-design.md` | Create |
| `docs/superpowers/plans/2026-06-22-stage61-idr-consensus-scope-policy.md` | Create |

## 4. Non-goals

- No new IDR rules, no consensus DAG integration
- No `.smk`, `Snakefile`, `validate_samples.py`, config, or Artifact changes
- No CUT&Tag narrow IDR implementation (deferred to future stage)
- No Co-Authored-By

## 5. Verification

```bash
python3 test/test_stage53_config_validation.py
python3 test/test_stage53_stage5_invariant.py
python3 test/test_stage53_pooled_not_validated.py
python3 test/test_stage53_experimental_warnings.py
python3 test/test_stage53_reproducibility_policy_contract.py
python3 test/test_stage53_output_path_templates.py
python3 test/test_stage55_config_validation.py
python3 test/test_stage55_stage5_invariant.py
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py
git diff --check
```
