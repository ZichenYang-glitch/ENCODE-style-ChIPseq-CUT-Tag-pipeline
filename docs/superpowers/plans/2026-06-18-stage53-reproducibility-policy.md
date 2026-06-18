# Stage 53 Implementation Plan: Reproducibility Policy

**Date:** 2026-06-18
**Status:** Ready for implementation
**Depends on:** None (policy/config layer only)
**Blocks:** Stage 54 (consensus engine), Stage 55 (ATAC narrow IDR)

## Context

The pipeline currently has one reproducibility path: Stage 5 IDR for ChIP-seq
narrow with exactly 2 biological replicates. All other peak-calling modes have
no replicate-validated outputs. Pooled peaks (`04_peaks/pooled/`) are
aggregate-signal, not validated.

Stage 53 is the **policy/config semantics foundation** only. It defines rules,
config surface, output naming conventions, and tests that verify the policy is
correctly encoded — without implementing any new runtime behavior.

## Files in scope

### Already written (spec phase)

| File | Purpose |
|------|---------|
| `docs/reproducibility-policy.md` | User-facing reproducibility policy document |
| `docs/superpowers/specs/2026-06-18-stage53-reproducibility-policy-design.md` | Full design specification |

### To modify

| File | Change |
|------|--------|
| `scripts/validate_samples.py` | Add `_validate_reproducibility()` helper; call from `validate_config()` |
| `config/config.yaml` | Add commented-out `reproducibility` block as reference |
| `workflow/schemas/config.schema.yaml` | Add `reproducibility` block schema documentation (human-readable contract) |

### To create

| File | What it tests |
|------|--------------|
| `test/test_stage53_config_validation.py` | Config parsing: defaults, bad values, null inference |
| `test/test_stage53_stage5_invariant.py` | `stage5` behavior not suppressed by `reproducibility` |
| `test/test_stage53_pooled_not_validated.py` | Policy doc + output-contract: pooled peaks ≠ validated |
| `test/test_stage53_experimental_warnings.py` | Experimental IDR flags emit informational warnings |
| `test/test_stage53_reproducibility_policy_contract.py` | Policy matrix completeness across all 6 modes |
| `test/test_stage53_output_path_templates.py` | Path template conformance (doc-level, no file existence) |

**Total: 3 files modified, 6 test files created, 2 spec files already written.**

### Explicitly NOT in scope

- ❌ No roadmap/changelog update (deferred to Stage 58 release hardening)
- ❌ No `CHANGELOG.md`, `KNOWN_ISSUES.md`, or `README.md` changes
- ❌ No `workflow/rules/*.smk` changes
- ❌ No `workflow/Snakefile` changes
- ❌ No new scripts under `scripts/`
- ❌ No manifest/output-contract runtime changes

## Implementation steps

### Step 1: Add `reproducibility` reference block to config

File: `config/config.yaml`

Insert after the `stage5` / `idr` section:

```yaml
# Stage 53+: Reproducibility policy for replicate-validated peak outputs.
# Added in Stage 53 — does not affect existing stage5 ChIP-seq narrow IDR.
# reproducibility:
#   enabled: false
#   consensus:
#     enabled: true
#     min_replicates: 2
#     reciprocal_overlap: 0.5
#   idr:
#     chipseq_narrow: null
#     atac_narrow: false
#     cuttag_narrow: false
#     chipseq_broad_experimental: false
#     cuttag_broad_experimental: false
```

### Step 2: Add `reproducibility` schema documentation

File: `workflow/schemas/config.schema.yaml`

Insert after the existing `idr:` block (after `# ---- Stage 5` section,
before `# ---- Stage 7b` section):

```yaml
# ---- Stage 53+: Reproducibility Policy (optional) ----
reproducibility:
  type: mapping
  required: false
  default:
    enabled: false
  description: >
    Stage 53+ reproducibility policy for replicate-validated peak outputs.
    When enabled is false, no new reproducibility outputs are generated.
    Does not affect existing stage5 ChIP-seq narrow IDR.
  schema:
    enabled:
      type: boolean | string
      default: false
      allowed: [true, false, "true", "false"]
      description: >
        Master switch for new reproducibility outputs. When false,
        all sub-keys are ignored for DAG purposes. Orthogonal to stage5.
    consensus:
      type: mapping
      required: false
      description: Consensus peak strategy (N-of-M replicate support).
      schema:
        enabled:
          type: boolean | string
          default: true
          description: Enable consensus peak generation.
        min_replicates:
          type: integer
          default: 2
          constraints: ">= 2"
          description: Minimum supporting biological replicates to retain a peak.
        reciprocal_overlap:
          type: float
          default: 0.5
          constraints: "> 0, <= 1"
          description: Reciprocal overlap threshold for pairwise peak overlap.
    idr:
      type: mapping
      required: false
      description: >
        IDR configuration for non-legacy modes. Does not affect stage5.
        No seacr_experimental key — SEACR IDR is not planned.
      schema:
        chipseq_narrow:
          type: boolean | null
          default: null
          description: >
            null/omitted → inferred from stage5. When explicitly false and
            stage5 is true, a config contradiction warning is emitted but
            legacy Stage 5 still runs.
        atac_narrow:
          type: boolean
          default: false
          description: ATAC-seq narrow-peak IDR (uses same narrowPeak machinery as ChIP-seq).
        cuttag_narrow:
          type: boolean
          default: false
          description: CUT&Tag narrow-peak IDR (opt-in, experimental).
        chipseq_broad_experimental:
          type: boolean
          default: false
          description: >
            Experimental ChIP-seq broad-peak IDR. Consensus remains primary
            final. When true, an informational warning is emitted.
        cuttag_broad_experimental:
          type: boolean
          default: false
          description: >
            Experimental CUT&Tag broad-peak IDR. Consensus remains primary
            final. When true, an informational warning is emitted.
```

### Step 3: Add `_validate_reproducibility()` helper

File: `scripts/validate_samples.py`

Add after the existing `stage5` validation block (after the `idr` settings
validation, approximately line 296).

Function signature and behavior:

```python
def _validate_reproducibility(raw, validated_config):
    """Validate the reproducibility config block (Stage 53+).

    Returns a dict with validated reproducibility settings.
    When enabled is false/absent, returns {'enabled': False}
    without validating sub-keys.
    """
    if not raw or not raw.get("enabled", False):
        return {"enabled": False}

    result = {"enabled": True}

    # --- consensus ---
    consensus = raw.get("consensus", {})
    result["consensus"] = {
        "enabled": _parse_bool(consensus.get("enabled", True), "consensus.enabled"),
        "min_replicates": _validate_min_replicates(consensus),
        "reciprocal_overlap": _validate_reciprocal_overlap(consensus),
    }

    # --- idr ---
    idr = raw.get("idr", {})
    result["idr"] = _validate_idr_block(idr, validated_config)

    return result
```

### Step 4: Wire into `validate_config()`

In `validate_config()`, after the `idr` settings block:

```python
# Stage 53: reproducibility policy
validated["reproducibility"] = _validate_reproducibility(
    config.get("reproducibility", {}), validated
)
```

### Step 5: Validation rules detail

1. `consensus.min_replicates`: default 2, must be int ≥ 2
2. `consensus.reciprocal_overlap`: default 0.5, must be float in (0, 1]
3. `idr.chipseq_narrow`: true / false / null. Default null. If explicitly false and `stage5: true`: emit warning; legacy Stage 5 still runs.
4. `idr.atac_narrow`: bool only, default false
5. `idr.cuttag_narrow`: bool only, default false
6. `idr.chipseq_broad_experimental`: bool only, default false. If true: emit informational warning.
7. `idr.cuttag_broad_experimental`: bool only, default false. If true: emit informational warning.
8. No `seacr_experimental` key — SEACR IDR is not planned.

### Step 6: Write tests

All 6 tests listed above. See test plan in the design spec for detailed
per-test expectations.

## Verification

```bash
# Run the Stage 53 test suite (standalone scripts)
cd /home/irenadler/workflow/chipseq
python3 test/test_stage53_config_validation.py
python3 test/test_stage53_stage5_invariant.py
python3 test/test_stage53_pooled_not_validated.py
python3 test/test_stage53_experimental_warnings.py
python3 test/test_stage53_reproducibility_policy_contract.py
python3 test/test_stage53_output_path_templates.py

# Run broader regression safety checks
python3 test/test_stage28_release_readiness.py
python3 test/test_no_hardcoded_paths.py

# Verify existing dry-run still passes (no regression)
conda run -n chipseq snakemake -s workflow/Snakefile \
  --configfile config/config.yaml -n --use-conda

# Git consistency check
git diff --check
```

## Non-goals (Stage 53)

- No `scripts/compute_consensus.py`
- No new Snakemake rules
- No DAG target expansion
- No `06_reproducibility/` directory creation
- No ATAC/CUT&Tag/broad IDR rules
- No manifest/output-contract runtime changes
- No artifact runtime adoption
- No roadmap/changelog/README updates
- Existing `stage5` and `06_idr/` behavior unchanged
