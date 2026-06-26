#!/usr/bin/env python3
"""Pytest shim for a curated allowlist of legacy test_stage*.py scripts.

Legacy scripts use hand-rolled pass/fail counters and expose one of three
entry-point conventions:

* a ``main()`` function
* a ``__main__`` guard block that calls test functions
* neither (pure pytest module, left to normal discovery)

Only scripts explicitly added to ``ALLOWLIST`` run through this shim. All
other ``test_stage*.py`` files are either collected normally by pytest (if
pytest-native) or ignored entirely (if legacy but not yet allowlisted).

Migration policy:
* A legacy script may be added to the allowlist only if it:
  - is legacy ``main()`` / ``__main__`` style, not pytest-native;
  - runs cleanly under the shim in the fast test environment;
  - has no tool-missing failures;
  - has no stale implementation-detail assertions (e.g. old IDR rule names);
  - does not mutate the repository working tree;
  - tests real pipeline behavior, not docs/CI/planning meta-structure.
* ``SystemExit(2)`` is treated as ``pytest.skip``, not failure, because some
  legacy scripts use that code to indicate "skipped due to missing tool".
* Stale or side-effectful legacy scripts remain quarantined here and in
  ``docs/development/harness.md`` until they are migrated to native pytest or
  retired.
"""

import importlib
import os
import runpy
import sys
from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
TEST_DIR = REPO_ROOT / "test"

# Curated allowlist of legacy stage scripts that are safe to run through the
# shim. Start empty until a script is proven stable against the criteria above.
ALLOWLIST: tuple[str, ...] = ()

# Known legacy scripts that are intentionally quarantined pending migration or
# retirement. Keep reasons explicit so future maintainers know why each one is
# not in the allowlist.
QUARANTINED: dict[str, str] = {
    # CI / planning meta-structure; not real pipeline behavior tests.
    "test_stage27_public_validation_plan.py": "validates planning docs, not pipeline behavior",
    "test_stage27b_metadata_ci_plan.py": "validates planning docs, not pipeline behavior",
    "test_stage27c_ci_workflow.py": "validates CI workflow structure; references removed stress tests",
    "test_stage28_release_readiness.py": "release checklist meta-test",
    "test_stage32_public_report_scaffold.py": "public report scaffold meta-test",
    "test_stage33_containerization_plan.py": "containerization plan meta-test; flags existing .sif artifact",
    "test_stage34_runner_container_files.py": "container file list meta-test",
    "test_stage35_docker_smoke_report.py": "container smoke report meta-test",
    "test_stage36_singularity_smoke_report.py": "container smoke report meta-test",
    "test_stage37_container_ux.py": "container UX meta-test",
    # Stress tests that may fail or mutate state in fast test environments.
    "test_stage12_stress.py": "legacy stress test; side effects in working tree",
    "test_stage18_19_stress.py": "legacy stress test",
    "test_stage22_bigwig_stress.py": "legacy stress test; may require bigWig tools",
    "test_stage30_strict_inputs_stress.py": "legacy stress test",
    "test_stage39_mnase_stress.py": "legacy stress test",
    "test_stage4b_stress.py": "legacy stress test",
    "test_stage4c_stress.py": "legacy stress test",
    "test_stage5a_stress.py": "legacy stress test; SystemExit(1) on failure",
    "test_stage5b_stress.py": "legacy stress test",
    "test_stage60_legacy_script.py": "legacy script smoke test",
    "test_stage6a_stress.py": "legacy stress test",
    "test_stage6b_stress.py": "legacy stress test",
    "test_stage7a_stress.py": "legacy stress test",
    "test_stage7b_stress.py": "legacy stress test",
    # Stale after Phase 1 rule consolidation; assert old IDR rule names.
    "test_stage55_atac_idr_dryrun.py": "asserts old ATAC IDR rule names after Phase 1 consolidation",
    "test_stage55_atac_idr_summary.py": "asserts old ATAC IDR rule names after Phase 1 consolidation",
    "test_stage64_cuttag_idr_dryrun.py": "asserts old CUT&Tag IDR rule names after Phase 1 consolidation",
    "test_stage64_cuttag_idr_summary.py": "asserts old CUT&Tag IDR rule names after Phase 1 consolidation",
    "test_stage65_broad_idr_dryrun.py": "asserts old broad IDR rule names after Phase 1 consolidation",
    "test_stage65_broad_idr_summary.py": "asserts old broad IDR rule names after Phase 1 consolidation",
    # Integration / tiny real execution; requires tools not present in fast env.
    "test_stage8b_tiny_execution.py": "tiny real execution; skips/fails when cutadapt is missing",
    # Config validation scripts largely superseded by native pytest harness.
    "test_stage53_config_validation.py": "superseded by test/config/test_validation.py",
    "test_stage53_experimental_warnings.py": "superseded by test/config/test_validation.py",
    "test_stage53_output_path_templates.py": "superseded by native pytest contracts",
    "test_stage53_pooled_not_validated.py": "superseded by native pytest contracts",
    "test_stage53_reproducibility_policy_contract.py": "superseded by native pytest contracts",
    "test_stage53_stage5_invariant.py": "superseded by native pytest contracts",
    "test_stage55_config_validation.py": "superseded by test/config/test_validation.py",
    "test_stage55_stage5_invariant.py": "superseded by native pytest contracts",
    "test_stage64_cuttag_idr_config_validation.py": "superseded by test/config/test_validation.py",
    "test_stage65_broad_idr_config_validation.py": "superseded by test/config/test_validation.py",
    # Other contracts now covered by native pytest harnesses.
    "test_stage24_qc_summary_unit.py": "covered by native pytest QC tests",
    "test_stage43_artifact_inventory.py": "covered by test/artifacts/test_catalog_contracts.py",
    "test_stage45_artifact_model.py": "covered by test/artifacts/test_catalog_contracts.py",
    "test_stage47_mnase_path_contract.py": "covered by native MNase pytest tests",
    "test_stage49_manifest_artifact_contract.py": "covered by test/artifacts/test_catalog_contracts.py",
    "test_stage50_output_contract_dry_run.py": "covered by manifest and DAG harnesses",
    "test_stage54_consensus.py": "covered by native consensus pytest tests",
    "test_stage57_shell_safety.py": "covered by native shell safety tests",
    "test_stage58_mixed_idr_validation.py": "covered by test/config/test_validation.py",
    "test_stage59_env_pinning.py": "covered by lockfile CI checks",
    "test_stage62_consensus_dryrun.py": "covered by native consensus pytest tests",
    "test_stage63_seacr_consensus_dryrun.py": "covered by native consensus pytest tests",
    "test_stage66_reproducibility_manifest.py": "covered by manifest harness",
}


def _allowlist_paths():
    """Return Paths for the explicitly allowlisted legacy scripts."""
    return [TEST_DIR / name for name in ALLOWLIST]


def _has_main(src):
    return "def main(" in src


def _has_main_guard(src):
    return 'if __name__ == "__main__":' in src


def _run_main(module_name):
    """Import the module and call its main() function."""
    module = importlib.import_module(module_name)
    try:
        result = module.main()
    except SystemExit as exc:
        code = exc.code if exc.code is not None else 0
        if code == 2:
            pytest.skip(f"{module_name}.main() exited with code 2 (skip)")
        if code != 0:
            pytest.fail(f"{module_name}.main() exited with code {code}")
        return
    if result is not None and result != 0:
        pytest.fail(f"{module_name}.main() returned {result}")


def _run_guard(path):
    """Execute the module as __main__ so its guard block runs."""
    try:
        runpy.run_path(str(path), run_name="__main__")
    except SystemExit as exc:
        code = exc.code if exc.code is not None else 0
        if code == 2:
            pytest.skip(f"{path.name} __main__ guard exited with code 2 (skip)")
        if code != 0:
            pytest.fail(f"{path.name} __main__ guard exited with code {code}")


def _module_name(path):
    """Return a dotted module name for an import."""
    return path.stem


@pytest.mark.parametrize(
    "stage_path",
    _allowlist_paths(),
    ids=lambda p: p.name,
)
def test_stage_shim(stage_path):
    """Run an allowlisted legacy test_stage*.py script under pytest."""
    if not stage_path.exists():
        pytest.fail(f"Allowlisted script missing: {stage_path}")

    src = stage_path.read_text(encoding="utf-8")
    module_name = _module_name(stage_path)

    if _has_main(src):
        _run_main(module_name)
    elif _has_main_guard(src):
        _run_guard(stage_path)
    else:
        pytest.skip(f"{stage_path.name} has no shim entry point; collected normally")
