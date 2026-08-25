# Development Harness

This document describes the maintained verification layers. Run focused checks
while iterating, then use the layer that matches the change's risk. CI keeps
each layer visible instead of hiding several independent suites behind one
opaque command.

## CI layers

| Layer | Triggers | Selection and purpose |
| --- | --- | --- |
| PR fast | Pull requests | Deterministic unit, contract, validator, and scientific DAG smoke tests, executed as two mutually exclusive path-hash shards with one combined coverage artifact; excludes `full_main` and all three real-execution markers. Changed Python lines must be at least 80% covered. |
| Full main | Push to `main`, `workflow_dispatch`, nightly schedule, published release | The complete deterministic suite in the same two-shard form, including `full_main`; excludes all three real-execution markers. Enforces the repository and core-module floors in the [coverage policy](coverage-policy.md). |
| Platform real | `workflow_dispatch`, nightly schedule, published release | Real Redis/RQ, SIGALRM timeout, process-group cancellation, and tiny Snakemake execution under `platform_real_execution`. |
| Scientific real | `workflow_dispatch`, nightly schedule, published release | Complete `test/real_execution` suite with real scientific tools under `real_execution`. |
| Bulk RNA-seq real | Explicit `workflow_dispatch` opt-in on a protected self-hosted runner | Verifies the pre-staged nf-core/rnaseq 3.26.0 closure, then runs rapid quantification, full STAR+Salmon/SortMeRNA, and lifecycle gates under `bulk_rnaseq_real_execution`. |
| Container smoke | `workflow_dispatch`, nightly schedule, published release | Builds the runner image and exercises the default container profile without publishing an image. |

The four real-execution jobs are intentionally not PR-required checks. A
high-risk PR can be validated by manually dispatching its exact branch. A
nightly or release trigger uses the workflow on the default branch, so a
workflow change cannot claim those external results until it is present there.

### One logical deterministic suite, two physical shards

Each event runs one logical deterministic suite as two mutually exclusive and
exhaustive physical shards. `scripts/split_deterministic_shards.py` assigns
every collected `test_*.py` file to a shard by `SHA-256(relative POSIX path)
% 2`; the partition is a pure function of the checkout, so both shards are
reproducible anywhere and their union is the complete tree. Each shard job
checks out the exact event SHA, runs only its own files with the event's
marker expression, and produces a shard JUnit report plus an unrendered
`.coverage.shard-N` data file:

```bash
python3 scripts/split_deterministic_shards.py --root test --shard 1 --of 2
```

The pull-request shards select:

```bash
-m "not full_main and not platform_real_execution and not real_execution and not bulk_rnaseq_real_execution"
```

Full-main shards select:

```bash
-m "not platform_real_execution and not real_execution and not bulk_rnaseq_real_execution"
```

Both shard selections keep `--cov-config=pyproject.toml --cov-context=test`,
`--cov-fail-under=0`, and an empty `--cov-report=`; no shard renders XML, JSON,
or terminal coverage output. The `fast-checks` aggregate job waits for both
shards, requires the exact two-shard evidence set, reruns the zero-skip JUnit
gate over both reports, and then runs `coverage combine` once to produce the
canonical `.coverage`, `coverage.xml`, and `coverage.json`.

The `coverage` job consumes the canonical artifact from that aggregate. It
enforces changed lines on PRs and the global/core floors on complete-suite
events; it does not invoke pytest a second time. These reports are uploaded as
short-lived CI artifacts and are never committed.

All five pytest entry points—PR fast, full main, platform real, scientific
real, and bulk RNA-seq real—must report zero skips and zero xfails in CI.
Markers are registered and checked strictly, and XPASS is a failure. The bulk
RNA-seq tier is dispatched only when its protected self-hosted runner and
pre-staged immutable assets are available; missing coordinates fail instead
of becoming skips.

## Stable checks and budgets

The recommended required contexts are:

- `fast-checks`
- `frontend`
- `browser-e2e`
- `lint`
- `lock-check`
- `coverage`

This is a repository policy recommendation, not a repository-settings change.
`platform-real-execution`, `real-execution`, `bulk-rnaseq-real-execution`, and
`container-smoke` should remain non-required unless they are later redesigned
as stable always-reporting checks.

CI reports the slowest tests and records wall time in the job summary. Soft
budgets are review signals; a cold cache crossing one does not fail the job.
Hard timeouts only stop hangs.

| Check or tier | Soft budget | Hard timeout |
| --- | ---: | ---: |
| PR `fast-checks` shards | 5 min | 25 min |
| PR `fast-checks` aggregate | 1 min | 10 min |
| Coverage artifact gate | 1 min | 5 min |
| Full-main Python | 20 min | 25 min |
| Frontend | 2 min | 8 min |
| Browser E2E | 4 min | 10 min |
| Lint | 1 min | 5 min |
| Lock check | 1 min | 3 min |
| Platform real execution | 6 min | 12 min |
| Scientific real execution | 20 min | 30 min |
| Container smoke | 20 min | 30 min |

## Lint

The `lint` check runs Python lint/format gates plus the workflow-specific
checks:

```bash
python3 -m ruff check src scripts test containers workflow/lib
python3 -m ruff format --check src scripts test containers workflow/lib
python3 test/check_snakemake_lint.py
snakefmt --check workflow/
```

`test/check_snakemake_lint.py` compares `snakemake --lint` output with the
accepted baseline in `docs/operations/snakemake-lint-warnings.txt`. Intentional
changes require a reviewed baseline update.

## Documentation contract

The maintained documentation suite checks that relative links resolve:

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  pytest test/docs/test_internal_links.py -v
```

## DAG snapshots and smoke profiles

- `test/workflow/test_smoke_profiles.py` is the PR-fast scientific validator
  and dry-run smoke.
- `test/test_dag_snapshots.py` and the broader assay DAG contracts carry the
  `full_main` marker.

## Lock-check behavior

The `lock-check` workflow verifies that every `workflow/envs/*.yml` has a
matching `.lock` file and that an event modifying an environment specification
also modifies its lock. It handles PR and push comparison bases separately.
It does not prove that the lock input hash or every transitive package is
fresh; CI separately installs the environments it uses.

## Local iteration

Use focused checks while editing:

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  pytest test/docs/test_internal_links.py -v
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  pytest test/workflow/test_smoke_profiles.py -q
python3 test/check_snakemake_lint.py
```

The local default pytest selection runs the complete deterministic suite and
excludes only both real-execution markers:

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest test -ra
```

See the [coverage policy](coverage-policy.md) for complete coverage commands
and [real-execution harness](real-execution-harness.md) for the explicit real
tiers.
