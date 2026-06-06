#!/usr/bin/env python3
"""Stage 47: MNase path-helper inventory contract tests.

Verifies that all 7 helpers in workflow/rules/paths.smk produce path
strings matching the 13 MNase entries in the artifact inventory.
"""

import os
import sys

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(REPO_ROOT, "workflow"))

from lib.artifact import load_artifacts, filter_artifacts, artifacts_by_id

_PASS = 0
_FAIL = 0


def _check(name, condition, detail=""):
    global _PASS, _FAIL
    if condition:
        print(f"PASS: {name}")
        _PASS += 1
    else:
        msg = f"FAIL: {name}"
        if detail:
            msg += f"  -- {detail}"
        print(msg)
        _FAIL += 1


def _load_paths_helpers():
    """exec() paths.smk in an isolated namespace with OUTDIR='results'."""
    paths_file = os.path.join(REPO_ROOT, "workflow", "rules", "paths.smk")
    namespace = {"OUTDIR": "results"}
    with open(paths_file) as fh:
        code = fh.read()
    exec(compile(code, paths_file, "exec"), namespace)
    return namespace


TEST_SAMPLE = "TEST_S1"
TEST_EXPERIMENT = "TEST_EXP"

# exact inventory IDs from artifact-inventory.yaml
MNASE_CONTRACT = {
    "mnase_mono_bam":       ("mnase_fragment_bam",       [TEST_SAMPLE, "mono"]),
    "mnase_mono_bai":       ("mnase_fragment_bai",       [TEST_SAMPLE, "mono"]),
    "mnase_sub_bam":        ("mnase_fragment_bam",       [TEST_SAMPLE, "sub"]),
    "mnase_sub_bai":        ("mnase_fragment_bai",       [TEST_SAMPLE, "sub"]),
    "mnase_di_bam":         ("mnase_fragment_bam",       [TEST_SAMPLE, "di"]),
    "mnase_di_bai":         ("mnase_fragment_bai",       [TEST_SAMPLE, "di"]),
    "mnase_dyad_bigwig":    ("mnase_signal_bw",          [TEST_SAMPLE, "dyad"]),
    "mnase_mono_bigwig":    ("mnase_signal_bw",          [TEST_SAMPLE, "mono"]),
    "mnase_qc_summary":     ("mnase_qc_summary_tsv",     [TEST_SAMPLE]),
    "pooled_mnase_mono_bam":      ("mnase_pooled_fragment_bam", [TEST_EXPERIMENT, "mono"]),
    "pooled_mnase_mono_bai":      ("mnase_pooled_fragment_bai", [TEST_EXPERIMENT, "mono"]),
    "pooled_mnase_dyad_bigwig":   ("mnase_pooled_signal_bw",    [TEST_EXPERIMENT, "dyad"]),
    "pooled_mnase_mono_bigwig":   ("mnase_pooled_signal_bw",    [TEST_EXPERIMENT, "mono"]),
}


def main():
    global _PASS, _FAIL

    inv_path = os.path.join(REPO_ROOT, "docs", "architecture",
                            "artifact-inventory.yaml")

    # 1. Load artifacts
    try:
        artifacts = load_artifacts(inv_path)
    except Exception as e:
        _check("1-load_artifacts_succeeds", False,
               f"load_artifacts raised {type(e).__name__}: {e}")
        print(f"\n{_PASS} passed, {_FAIL} failed, {_PASS + _FAIL} total")
        sys.exit(1)
    _check("1-load_artifacts_succeeds", True)

    # 2. Filter MNase entries
    mnase_arts = filter_artifacts(artifacts, assay_gate="mnase")
    _check("2-mnase_entries_13", len(mnase_arts) == 13,
           f"expected 13, got {len(mnase_arts)}")

    # 3. paths.smk has zero rule declarations
    paths_file = os.path.join(REPO_ROOT, "workflow", "rules", "paths.smk")
    with open(paths_file) as fh:
        paths_lines = fh.readlines()
    rule_lines = [ln.strip() for ln in paths_lines
                  if ln.strip().startswith("rule ")]
    _check("3-paths_smk_no_rules", len(rule_lines) == 0,
           f"found rule declarations: {rule_lines}")

    # 4. Load helpers and verify all 7 present
    helpers = _load_paths_helpers()
    expected_helpers = {
        "mnase_fragment_bam", "mnase_fragment_bai",
        "mnase_signal_bw", "mnase_qc_summary_tsv",
        "mnase_pooled_fragment_bam", "mnase_pooled_fragment_bai",
        "mnase_pooled_signal_bw",
    }
    missing_helpers = []
    not_callable = []
    for name in sorted(expected_helpers):
        if name not in helpers:
            missing_helpers.append(name)
        elif not callable(helpers[name]):
            not_callable.append(name)
    _check("4-helpers_all_7_present",
           len(missing_helpers) == 0 and len(not_callable) == 0,
           f"missing: {missing_helpers}, not callable: {not_callable}")

    # 5-17. Contract check per MNase entry
    inv_by_id = artifacts_by_id(mnase_arts)
    for inv_id in sorted(MNASE_CONTRACT):
        helper_name, args = MNASE_CONTRACT[inv_id]
        artifact = inv_by_id.get(inv_id)
        if artifact is None:
            _check(f"contract_{inv_id}", False,
                   f"inventory entry not found for id '{inv_id}'")
            continue

        func = helpers.get(helper_name)
        if not callable(func):
            _check(f"contract_{inv_id}", False,
                   f"helper '{helper_name}' missing or not callable")
            continue

        expected = artifact.path_template.replace(
            "<sample>", TEST_SAMPLE
        ).replace(
            "<experiment>", TEST_EXPERIMENT
        )
        actual = func(*args)
        _check(f"contract_{inv_id}", actual == expected,
               f"expected '{expected}', got '{actual}'")

    total = _PASS + _FAIL
    print(f"\n{_PASS} passed, {_FAIL} failed, {total} total")
    sys.exit(0 if _FAIL == 0 else 1)


if __name__ == "__main__":
    main()
