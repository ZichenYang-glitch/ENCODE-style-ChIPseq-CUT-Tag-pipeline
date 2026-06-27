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

# Approved classification categories for legacy stage scripts.
_APPROVED_CATEGORIES: frozenset[str] = frozenset(
    {
        "migrate-to-pytest",
        "delete-candidate",
        "manual-integration",
        "real-execution-only",
        "obsolete-plan-doc",
        "keep-quarantined-for-now",
    }
)

# Structured classification of every legacy-style test_stage*.py script.
# This is the single source of truth; QUARANTINED is derived from it.
LEGACY_STAGE_CLASSIFICATION: dict[str, dict[str, str]] = {
    # CI / planning meta-structure; not real pipeline behavior tests.
    "test_stage27_public_validation_plan.py": {
        "category": "obsolete-plan-doc",
        "rationale": "validates planning docs, not pipeline behavior",
        "next_action": "delete when planning docs are archived or refreshed",
    },
    "test_stage27b_metadata_ci_plan.py": {
        "category": "obsolete-plan-doc",
        "rationale": "validates planning docs, not pipeline behavior",
        "next_action": "delete when planning docs are archived or refreshed",
    },
    "test_stage27c_ci_workflow.py": {
        "category": "obsolete-plan-doc",
        "rationale": "validates CI workflow structure; references removed stress tests",
        "next_action": "migrate useful assertions to native pytest CI contract or delete",
    },
    "test_stage28_release_readiness.py": {
        "category": "obsolete-plan-doc",
        "rationale": "release checklist meta-test",
        "next_action": "delete; release readiness is a manual process",
    },
    "test_stage32_public_report_scaffold.py": {
        "category": "obsolete-plan-doc",
        "rationale": "public report scaffold meta-test",
        "next_action": "delete when public report template is retired",
    },
    "test_stage33_containerization_plan.py": {
        "category": "obsolete-plan-doc",
        "rationale": "containerization plan meta-test; flags existing .sif artifact",
        "next_action": "delete; container support is now implemented",
    },
    "test_stage34_runner_container_files.py": {
        "category": "obsolete-plan-doc",
        "rationale": "container file list meta-test",
        "next_action": "delete; runner container files are now implemented",
    },
    "test_stage35_docker_smoke_report.py": {
        "category": "obsolete-plan-doc",
        "rationale": "container smoke report meta-test",
        "next_action": "delete; container smoke reports are superseded by CI",
    },
    "test_stage36_singularity_smoke_report.py": {
        "category": "obsolete-plan-doc",
        "rationale": "container smoke report meta-test",
        "next_action": "delete; container smoke reports are superseded by CI",
    },
    "test_stage37_container_ux.py": {
        "category": "obsolete-plan-doc",
        "rationale": "container UX meta-test",
        "next_action": "delete; container UX is now implemented",
    },
    # Stress tests that may fail or mutate state in fast test environments.
    "test_stage12_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test; side effects in working tree",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage18_19_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage22_bigwig_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test; may require bigWig tools",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage30_strict_inputs_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test",
        "next_action": "migrate coverage to test/config/test_validation.py or retire",
    },
    "test_stage39_mnase_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage4b_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage4c_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage5a_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test; SystemExit(1) on failure",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage5b_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage60_legacy_script.py": {
        "category": "manual-integration",
        "rationale": "legacy script smoke test",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage6a_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage6b_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage7a_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test",
        "next_action": "migrate to real-execution harness or retire",
    },
    "test_stage7b_stress.py": {
        "category": "manual-integration",
        "rationale": "legacy stress test",
        "next_action": "migrate to real-execution harness or retire",
    },
    # Stale after Phase 1 rule consolidation; assert old IDR rule names.
    "test_stage55_atac_idr_dryrun.py": {
        "category": "delete-candidate",
        "rationale": "asserts old ATAC IDR rule names after Phase 1 consolidation",
        "next_action": "delete; ATAC IDR behavior is covered by native pytest",
    },
    "test_stage55_atac_idr_summary.py": {
        "category": "delete-candidate",
        "rationale": "asserts old ATAC IDR rule names after Phase 1 consolidation",
        "next_action": "delete; ATAC IDR behavior is covered by native pytest",
    },
    "test_stage64_cuttag_idr_dryrun.py": {
        "category": "delete-candidate",
        "rationale": "asserts old CUT&Tag IDR rule names after Phase 1 consolidation",
        "next_action": "delete; CUT&Tag IDR behavior is covered by native pytest",
    },
    "test_stage64_cuttag_idr_summary.py": {
        "category": "delete-candidate",
        "rationale": "asserts old CUT&Tag IDR rule names after Phase 1 consolidation",
        "next_action": "delete; CUT&Tag IDR behavior is covered by native pytest",
    },
    "test_stage65_broad_idr_dryrun.py": {
        "category": "delete-candidate",
        "rationale": "asserts old broad IDR rule names after Phase 1 consolidation",
        "next_action": "delete; broad IDR behavior is covered by native pytest",
    },
    "test_stage65_broad_idr_summary.py": {
        "category": "delete-candidate",
        "rationale": "asserts old broad IDR rule names after Phase 1 consolidation",
        "next_action": "delete; broad IDR behavior is covered by native pytest",
    },
    # Integration / tiny real execution; requires tools not present in fast env.
    "test_stage8b_tiny_execution.py": {
        "category": "real-execution-only",
        "rationale": "tiny real execution; skips/fails when cutadapt is missing",
        "next_action": "migrate to container/real-execution smoke harness",
    },
    # Config validation scripts largely superseded by native pytest harness.
    "test_stage53_config_validation.py": {
        "category": "migrate-to-pytest",
        "rationale": "superseded by test/config/test_validation.py",
        "next_action": "migrate remaining coverage to test/config/test_validation.py then delete",
    },
    "test_stage53_experimental_warnings.py": {
        "category": "migrate-to-pytest",
        "rationale": "superseded by test/config/test_validation.py",
        "next_action": "migrate remaining coverage to test/config/test_validation.py then delete",
    },
    "test_stage53_output_path_templates.py": {
        "category": "migrate-to-pytest",
        "rationale": "superseded by native pytest contracts",
        "next_action": "migrate remaining assertions to native pytest contract tests then delete",
    },
    "test_stage53_pooled_not_validated.py": {
        "category": "migrate-to-pytest",
        "rationale": "superseded by native pytest contracts",
        "next_action": "migrate remaining assertions to native pytest contract tests then delete",
    },
    "test_stage53_reproducibility_policy_contract.py": {
        "category": "migrate-to-pytest",
        "rationale": "superseded by native pytest contracts",
        "next_action": "migrate remaining assertions to native pytest contract tests then delete",
    },
    "test_stage53_stage5_invariant.py": {
        "category": "migrate-to-pytest",
        "rationale": "superseded by native pytest contracts",
        "next_action": "migrate remaining assertions to native pytest contract tests then delete",
    },
    "test_stage55_config_validation.py": {
        "category": "migrate-to-pytest",
        "rationale": "superseded by test/config/test_validation.py",
        "next_action": "migrate remaining coverage to test/config/test_validation.py then delete",
    },
    "test_stage55_stage5_invariant.py": {
        "category": "migrate-to-pytest",
        "rationale": "superseded by native pytest contracts",
        "next_action": "migrate remaining assertions to native pytest contract tests then delete",
    },
    "test_stage64_cuttag_idr_config_validation.py": {
        "category": "migrate-to-pytest",
        "rationale": "superseded by test/config/test_validation.py",
        "next_action": "migrate remaining coverage to test/config/test_validation.py then delete",
    },
    "test_stage65_broad_idr_config_validation.py": {
        "category": "migrate-to-pytest",
        "rationale": "superseded by test/config/test_validation.py",
        "next_action": "migrate remaining coverage to test/config/test_validation.py then delete",
    },
    # Other contracts now covered by native pytest harnesses.
    "test_stage24_qc_summary_unit.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by native pytest QC tests",
        "next_action": "verify coverage and delete",
    },
    "test_stage43_artifact_inventory.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by test/artifacts/test_catalog_contracts.py",
        "next_action": "verify coverage and delete",
    },
    "test_stage45_artifact_model.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by test/artifacts/test_catalog_contracts.py",
        "next_action": "verify coverage and delete",
    },
    "test_stage47_mnase_path_contract.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by native MNase pytest tests",
        "next_action": "verify coverage and delete",
    },
    "test_stage49_manifest_artifact_contract.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by test/artifacts/test_catalog_contracts.py",
        "next_action": "verify coverage and delete",
    },
    "test_stage50_output_contract_dry_run.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by manifest and DAG harnesses",
        "next_action": "verify coverage and delete",
    },
    "test_stage54_consensus.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by native consensus pytest tests",
        "next_action": "verify coverage and delete",
    },
    "test_stage57_shell_safety.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by native shell safety tests",
        "next_action": "verify coverage and delete",
    },
    "test_stage58_mixed_idr_validation.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by test/config/test_validation.py",
        "next_action": "verify coverage and delete",
    },
    "test_stage59_env_pinning.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by lockfile CI checks",
        "next_action": "verify coverage and delete",
    },
    "test_stage62_consensus_dryrun.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by native consensus pytest tests",
        "next_action": "verify coverage and delete",
    },
    "test_stage63_seacr_consensus_dryrun.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by native consensus pytest tests",
        "next_action": "verify coverage and delete",
    },
    "test_stage66_reproducibility_manifest.py": {
        "category": "migrate-to-pytest",
        "rationale": "covered by manifest harness",
        "next_action": "verify coverage and delete",
    },
}

# Known legacy scripts that are intentionally quarantined pending migration or
# retirement. Derived from LEGACY_STAGE_CLASSIFICATION so it stays in sync.
QUARANTINED: dict[str, str] = {
    name: meta["rationale"]
    for name, meta in LEGACY_STAGE_CLASSIFICATION.items()
    if name not in ALLOWLIST
}


def _is_legacy_script(path):
    """Return True if the file looks like a legacy main()/__main__ script."""
    src = path.read_text(encoding="utf-8")
    return _has_main(src) or _has_main_guard(src)


def _discover_legacy_scripts():
    """Return all root-level legacy-style test_stage*.py files except this shim."""
    return sorted(
        p
        for p in TEST_DIR.glob("test_stage*.py")
        if p.name != "test_stage_shim.py" and _is_legacy_script(p)
    )


def test_legacy_scripts_are_classified():
    """Guard: every legacy-style stage script must be allowlisted or quarantined."""
    classified = set(ALLOWLIST) | set(QUARANTINED)
    legacy_scripts = _discover_legacy_scripts()
    unclassified = [p.name for p in legacy_scripts if p.name not in classified]

    if unclassified:
        pytest.fail(
            "New legacy stage script(s) must be allowlisted or quarantined: "
            + ", ".join(sorted(unclassified))
        )


def test_allowlist_and_quarantine_refer_to_existing_scripts():
    """Every allowlisted or quarantined entry must exist at the test root."""
    classified = set(ALLOWLIST) | set(QUARANTINED)
    missing = [name for name in classified if not (TEST_DIR / name).exists()]
    if missing:
        pytest.fail(
            "ALLOWLIST/QUARANTINED refer to missing scripts: "
            + ", ".join(sorted(missing))
        )


def test_allowlist_and_quarantine_are_disjoint():
    """A script must not be both allowlisted and quarantined."""
    overlap = set(ALLOWLIST) & set(QUARANTINED)
    if overlap:
        pytest.fail(
            "Script(s) appear in both ALLOWLIST and QUARANTINED: "
            + ", ".join(sorted(overlap))
        )


def test_classifications_use_approved_categories():
    """Every classification entry must use an approved category."""
    bad = [
        (name, meta.get("category"))
        for name, meta in LEGACY_STAGE_CLASSIFICATION.items()
        if meta.get("category") not in _APPROVED_CATEGORIES
    ]
    if bad:
        details = "; ".join(f"{n}: {c!r}" for n, c in bad)
        pytest.fail(f"Unapproved classification categories: {details}")


def test_classifications_have_required_fields():
    """Every classification entry must have non-empty rationale and next_action."""
    bad = []
    for name, meta in LEGACY_STAGE_CLASSIFICATION.items():
        for field in ("category", "rationale", "next_action"):
            if not meta.get(field):
                bad.append(f"{name}: missing or empty {field}")
    if bad:
        pytest.fail("Incomplete classification entries: " + "; ".join(bad))


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
