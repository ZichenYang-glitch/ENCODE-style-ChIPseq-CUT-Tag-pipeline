# Legacy Test Migration Guide

This document tracks the legacy `test_stage*.py` scripts and the plan for
migrating or retiring them. The single source of truth for classification is
`LEGACY_STAGE_CLASSIFICATION` in `test/test_stage_shim.py`. This document
summarizes that mapping; if the two drift, trust the Python mapping.

## Current counts

| Category | Count |
|---|---|
| `manual-integration` | 14 |
| `delete-candidate` | 6 |
| `real-execution-only` | 1 |
| `migrate-to-pytest` | 19 |
| **Remaining classified legacy scripts** | **40** |

The documentation reset retired all 10 `obsolete-plan-doc` scripts. They
validated historical plans, checklists, or reports rather than product or
scientific behavior and are no longer part of the classification.

`test_stage8_smoke_profiles.py` is a pytest-native module and is not part of
the legacy classification.

## Category definitions

- **`migrate-to-pytest`**: Logic is already covered by native pytest tests, or
the script should be converted into native pytest tests.
- **`delete-candidate`**: Superseded by native tests or no longer valuable.
- **`manual-integration`**: Requires real data, external tools, or a full
environment not available in fast CI.
- **`real-execution-only`**: Belongs in a real-execution or container-smoke
harness rather than the fast test suite.
- **`keep-quarantined-for-now`**: Not yet ready for action; needs further
review.

## Scripts by category

### `manual-integration` (14)

- `test_stage12_stress.py`
- `test_stage18_19_stress.py`
- `test_stage22_bigwig_stress.py`
- `test_stage30_strict_inputs_stress.py`
- `test_stage39_mnase_stress.py`
- `test_stage4b_stress.py`
- `test_stage4c_stress.py`
- `test_stage5a_stress.py`
- `test_stage5b_stress.py`
- `test_stage60_legacy_script.py`
- `test_stage6a_stress.py`
- `test_stage6b_stress.py`
- `test_stage7a_stress.py`
- `test_stage7b_stress.py`

### `delete-candidate` (6)

These assert old IDR rule names that no longer exist after Phase 1
consolidation:

- `test_stage55_atac_idr_dryrun.py`
- `test_stage55_atac_idr_summary.py`
- `test_stage64_cuttag_idr_dryrun.py`
- `test_stage64_cuttag_idr_summary.py`
- `test_stage65_broad_idr_dryrun.py`
- `test_stage65_broad_idr_summary.py`

### `real-execution-only` (1)

- `test_stage8b_tiny_execution.py`

### `migrate-to-pytest` (19)

Config validation superseded by `test/config/test_validation.py`:

- `test_stage53_config_validation.py`
- `test_stage53_experimental_warnings.py`
- `test_stage53_output_path_templates.py`
- `test_stage53_pooled_not_validated.py`
- `test_stage53_reproducibility_policy_contract.py`
- `test_stage53_stage5_invariant.py`
- `test_stage55_config_validation.py`
- `test_stage55_stage5_invariant.py`
- `test_stage64_cuttag_idr_config_validation.py`
- `test_stage65_broad_idr_config_validation.py`

Other contracts covered by native pytest harnesses:

- `test_stage24_qc_summary_unit.py`
- `test_stage47_mnase_path_contract.py`
- `test_stage49_manifest_artifact_contract.py`
- `test_stage50_output_contract_dry_run.py`
- `test_stage54_consensus.py`
- `test_stage58_mixed_idr_validation.py`
- `test_stage62_consensus_dryrun.py`
- `test_stage63_seacr_consensus_dryrun.py`
- `test_stage66_reproducibility_manifest.py`

## Completed Batch 1 (PR49)

PR49 migrated the following low-risk scripts to native pytest and deleted the
legacy originals:

1. `test_stage57_shell_safety.py` → `test/test_shell_safety.py`
2. `test_stage59_env_pinning.py` → `test/test_env_pinning.py`
3. `test_stage43_artifact_inventory.py` and `test_stage45_artifact_model.py` →
   extended `test/artifacts/test_catalog_contracts.py`

Future batches should continue with the remaining `migrate-to-pytest` scripts
using the same rules below.

## Rules for retiring or migrating a legacy script

1. Verify native pytest coverage before deleting a `migrate-to-pytest` script.
2. If migrating, rewrite the logic in `test/` using pytest conventions, shared
   fixtures, and assertions. Do not preserve hand-rolled pass/fail counters.
3. If deleting, remove the script and its entry from
   `LEGACY_STAGE_CLASSIFICATION` in the same PR.
4. If a script is behavior-critical and cannot be covered by native tests yet,
   change its category to `keep-quarantined-for-now` and document the blocker.
5. Never move a script to `delete-candidate` without explaining what covers its
   behavior.
6. Stress tests and real-execution tests should move to a dedicated harness,
   not into the fast pytest suite.

## Status

The structured classification remains temporary migration machinery. The
documentation reset removed all plan-document tests; the maintenance baseline
must give every remaining entry a final action and then delete this guide and
the shim.
