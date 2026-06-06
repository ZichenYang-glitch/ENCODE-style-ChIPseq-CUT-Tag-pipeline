# Stage 47: MNase Path-Helper Inventory Contract Tests — Design Spec

**Date:** 2026-06-07
**Branch:** stage47-mnase-path-contract-tests
**Status:** implemented 2026-06-07

## Purpose

Add a test that verifies the 7 MNase path helpers in `workflow/rules/paths.smk` produce path strings matching the 13 MNase entries in `docs/architecture/artifact-inventory.yaml`.

## Scope

### In scope

- Create `test/test_stage47_mnase_path_contract.py`.
- Load artifacts via `load_artifacts()` from `workflow/lib/artifact.py`.
- Load paths.smk helpers via `exec()` in an isolated namespace with `OUTDIR="results"`.
- One contract check per MNase inventory entry: helper output == expected path from inventory `path_template` with placeholders substituted.
- Structural checks: paths.smk has zero `rule` declarations, all 7 helpers present and callable.
- Update `docs/architecture/artifact-roadmap.md`.

### Out of scope

- No changes to `workflow/rules/paths.smk`.
- No changes to `workflow/Snakefile`, `.smk` files, `rule all`.
- No changes to `make_manifest.py`, `artifact.py`.
- No `artifact_path()`.
- No DAG/runtime changes.

## Design

### Helper loading

```python
def _load_paths_helpers():
    paths_file = os.path.join(REPO_ROOT, "workflow", "rules", "paths.smk")
    namespace = {"OUTDIR": "results"}
    with open(paths_file) as fh:
        code = fh.read()
    exec(compile(code, paths_file, "exec"), namespace)
    return namespace
```

### Contract table

Built from the exact inventory IDs:

```python
TEST_SAMPLE = "TEST_S1"
TEST_EXPERIMENT = "TEST_EXP"

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
```

13 entries, each mapping an inventory `id` → `(helper_name, args)`.

### Expected path computation

For each entry:
1. Look up inventory artifact by `id`: `a.path_template`
2. Substitute: `pt.replace("<sample>", TEST_SAMPLE).replace("<experiment>", TEST_EXPERIMENT)`
3. Call helper: `actual = helpers[helper_name](*args)`
4. Assert: `actual == expected_path`

### Test checks (17 total)

| # | Check name | What |
|:---|:---|:---|
| 1 | `load_artifacts_succeeds` | `load_artifacts()` returns list |
| 2 | `mnase_entries_13` | Exactly 13 MNase entries |
| 3 | `paths_smk_no_rules` | No line in paths.smk starts with `rule ` |
| 4 | `helpers_all_7_present` | All 7 helper names exist in namespace and are callable |
| 5-17 | One per MNase entry (13) | `contract_<id>`: helper output matches inventory path_template |

### Roadmap update

Minimal: add Stage 47 implementation note.

## Verification

```bash
python3 test/test_stage47_mnase_path_contract.py    # 17/17 PASS
python3 test/test_stage43_artifact_inventory.py     # 10/10 PASS
python3 test/test_stage45_artifact_model.py         # 22/22 PASS
python3 test/test_stage28_release_readiness.py      # 11/11 PASS
python3 test/test_no_hardcoded_paths.py             # PASS
git diff --check                                      # clean
```
