# Real-Execution and Container Smoke Harness

This document defines where tests that need real services, processes,
scientific tools, or containers run. They stay outside the default
deterministic suite and are selected by explicit markers or entry points.

## Execution tiers

| Tier / CI job | Trigger | Environment | Purpose |
| --- | --- | --- | --- |
| `platform-real-execution` | Manual `workflow_dispatch`, nightly `schedule`, published `release` | `workflow/envs/ci-fast.lock` plus Redis 7 | Redis/RQ worker, SIGALRM timeout, process-group cancellation, and tiny real Snakemake execution. |
| `real-execution` | Manual `workflow_dispatch`, nightly `schedule`, published `release` | `workflow/envs/chipseq.lock` | Complete scientific real-tool suite, including focused samtools contracts and the tiny end-to-end preprocessing run. |
| `container-smoke` | Manual `workflow_dispatch`, nightly `schedule`, published `release` | Docker runner image plus `profiles/default` | Build and smoke-test the runner image without pushing it. |
| Local container smoke | Developer opt-in | Docker, Apptainer, or SingularityCE plus `profiles/default` | Check mounts, Conda prefix, runner tools, and a workflow dry-run. |
| Manual HPC smoke | Site-specific opt-in | `profiles/hpc` | Validate scheduler, shadow/storage behavior, and site integration. |

The three CI real jobs are deliberately non-required. They do not run on every
pull request. Dispatch the exact PR branch when a high-risk worker, scientific,
environment, profile, or container change needs them. Nightly and published
release events use the workflow from the default branch; before a workflow
change is merged, only its manual branch dispatch can validate those new
definitions.

The nightly cron is `17 19 * * *` UTC (03:17 the next day in Asia/Shanghai).

## Concrete entry points

- `.github/workflows/ci.yml` defines `platform-real-execution`,
  `real-execution`, and `container-smoke`.
- `test/workers/test_redis_process_integration.py` covers Redis/RQ and SIGALRM.
- `test/workers/test_tiny_execution_e2e.py` covers the API-to-worker tiny
  Snakemake lifecycle.
- `test/workers/test_cancellation_e2e.py` covers process-group cancellation.
- `test/real_execution/test_scientific_tiny_preprocessing.py` is the tiny
  scientific end-to-end run.
- `test/real_execution/test_pseudoreplicate_splitting.py` and
  `test/real_execution/test_cuttag_fragment_size.py` are focused real-tool
  contracts that require `samtools`.
- `scripts/smoke_container_runner.sh` is the container smoke runner.
- `profiles/default/config.yaml` is the local container profile.
- `profiles/hpc/config.yaml` is the manual HPC profile template.
- `docs/container-usage.md` documents the runner-only image.

## Platform real execution

The platform tier is marked `platform_real_execution`. It uses a real Redis
service, a separate RQ worker, file-backed SQLite, real POSIX processes, and
the real Snakemake executable:

```bash
ENCODE_PIPELINE_TEST_REDIS_URL=redis://127.0.0.1:6379/15 \
  python3 -m pytest -m platform_real_execution test/workers -ra -v
```

Without the dedicated Redis URL, an opt-in local run may skip tests whose
service was not requested. In CI the URL and Redis service are always present;
missing Snakemake, POSIX process control, or another required prerequisite is
a failure. The job emits JUnit and rejects any skip or xfail, so an empty or
partially skipped run cannot pass.

## Scientific real execution

The scientific tier is marked `real_execution` and runs generated or tiny
synthetic data through real tools:

```bash
python3 -m pytest -m real_execution test/real_execution -ra -v
```

The default pytest selection excludes this marker. The `chipseq` lock includes
pytest and every executable used by the tier. CI sets
`HELIXWEAVE_REQUIRE_REAL_EXECUTION=1`; a missing executable therefore fails
instead of becoming an all-skipped success. Its JUnit report must also contain
zero skips and zero xfails.

## Container smoke

The CI entry point builds the runner image locally and calls the maintained
smoke script:

```bash
docker build -f containers/Dockerfile.runner -t helixweave-runner:ci .
bash scripts/smoke_container_runner.sh docker helixweave-runner:ci
```

It validates startup, bind mounts, Conda prefix behavior, and the default
profile dry-run. It does not log in to a registry, push an image, or run heavy
scientific compute.

## Manual HPC smoke

HPC validation remains site-specific. It covers the `profiles/hpc` scheduler,
shadow directory, storage plugins, and cluster integration. It is not a GitHub
required check and requires access to the target cluster.

## Zero-skip CI contract

PR fast, full main, platform real, and scientific real are four separate
pytest entry points. Each writes a JUnit report, runs with strict registered
markers and strict xfail behavior, and must finish with zero skips and zero
xfails. Local missing-prerequisite skips are a development convenience only;
the two real CI jobs are fail-closed.
