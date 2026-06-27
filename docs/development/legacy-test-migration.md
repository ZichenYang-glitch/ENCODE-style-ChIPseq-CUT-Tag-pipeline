# Legacy Test Migration Guide

This document tracks the legacy `test_stage*.py` scripts and the plan for
migrating or retiring them. The single source of truth for classification is
`LEGACY_STAGE_CLASSIFICATION` in `test/test_stage_shim.py`. This document
summarizes that mapping; if the two drift, trust the Python mapping.

## Current counts

| Category | Count |
|---|---|
| `obsolete-plan-doc` | 10 |
| `manual-integration` | 14 |
| `delete-candidate` | 6 |
| `real-execution-only` | 1 |
| `migrate-to-pytest` | 23 |
| **Total classified legacy scripts** | **54** |

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
- **`obsolete-plan-doc`**: Validates planning or checklist documents that are
no longer active.
- **`keep-quarantined-for-now`**: Not yet ready for action; needs further
review.

## Scripts by category

### `obsolete-plan-doc` (10)

- `test_stage27_public_validation_plan.py`
- `test_stage27b_metadata_ci_plan.py`
- `test_stage27c_ci_workflow.py`
- `test_stage28_release_readiness.py`
- `test_stage32_public_report_scaffold.py`
- `test_stage33_containerization_plan.py`
- `test_stage34_runner_container_files.py`
- `test_stage35_docker_smoke_report.py`
- `test_stage36_singularity_smoke_report.py`
- `test_stage37_container_ux.py`

### `manual-integration` (15)

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

### `migrate-to-pytest` (22)

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
- `test_stage43_artifact_inventory.py`
- `test_stage45_artifact_model.py`
- `test_stage47_mnase_path_contract.py`
- `test_stage49_manifest_artifact_contract.py`
- `test_stage50_output_contract_dry_run.py`
- `test_stage54_consensus.py`
- `test_stage57_shell_safety.py`
- `test_stage58_mixed_idr_validation.py`
- `test_stage59_env_pinning.py`
- `test_stage62_consensus_dryrun.py`
- `test_stage63_seacr_consensus_dryrun.py`
- `test_stage66_reproducibility_manifest.py`

## Recommended Batch 1 for PR49

Start with low-risk scripts that have no external tool dependency, no Snakemake
execution, and no working-tree side effects:

1. `test_stage57_shell_safety.py` — scans rule files for `set -e -o pipefail`.
2. `test_stage59_env_pinning.py` — parses conda env YAMLs for version pins.
3. `test_stage43_artifact_inventory.py` — reads artifact inventory via
   `lib.artifact`.
4. `test_stage45_artifact_model.py` — tests `lib.artifact` models and helpers.

Batch 1 should either be deleted (if native pytest coverage is complete) or
rewritten as native pytest tests. Do not add them to the shim allowlist.

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

PR48 introduced the structured classification only. PR49 will begin with Batch
1 migration or deletion.
