# Development Harness

This document describes the current development verification harness for the
pipeline. Run these commands locally before opening a PR and rely on the same
commands in CI.

## Lint

Two lint commands guard the workflow rules:

```bash
# Compare current snakemake --lint output against the accepted baseline.
python3 test/check_snakemake_lint.py

# Check workflow rule formatting.
snakefmt --check workflow/
```

`test/check_snakemake_lint.py` records an accepted `--lint` baseline in
`docs/operations/snakemake-lint-warnings.txt`. The test fails if new warnings
appear; intentional changes require updating the baseline.

## CI contract

`test/test_stage27c_ci_workflow.py` asserts that `.github/workflows/ci.yml`
contains the expected fast-check commands. Run it with:

```bash
python3 test/test_stage27c_ci_workflow.py
```

## Full pytest suite

Run the full test suite the same way CI does:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -q -p no:cacheprovider
```

## DAG snapshots and smoke profiles

- `test/test_dag_snapshots.py` validates DAG target snapshots.
- `test/test_stage8_smoke_profiles.py` smoke-tests profile execution without
  running heavy computation.

## Legacy stage shim

`test/test_stage_shim.py` provides a pytest-compatible wrapper for legacy
`test_stage*.py` scripts. Some legacy scripts are quarantined because they are
slow, depend on real data, or reference outputs that no longer exist.

The shim uses an allowlist/quarantine classification:

- **Allowlisted** scripts are imported and run as normal pytest tests.
- **Quarantined** scripts are skipped by default and only run when explicitly
  requested.

This lets CI keep fast, reliable tests while preserving the legacy scripts for
manual archaeology.

## Lock-check behavior

The `lock-check` workflow verifies that every `workflow/envs/*.yml` has a
matching `.lock` file and that PRs modifying an env YAML also modify the
corresponding lockfile. It does not install every environment or prove lockfile
freshness; CI separately installs the fast-check environment from
`workflow/envs/ci-fast.lock`.

## When to run focused vs full verification

Use focused verification while iterating:

```bash
python3 test/check_snakemake_lint.py
python3 test/test_stage27c_ci_workflow.py
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest test/test_stage_shim.py -q -p no:cacheprovider
```

Run the full suite before marking a PR ready:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest -q -p no:cacheprovider
```
