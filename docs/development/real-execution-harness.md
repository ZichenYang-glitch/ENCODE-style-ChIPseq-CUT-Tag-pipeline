# Real-Execution and Container Smoke Harness

This document defines where tests that need real services, processes,
scientific tools, or containers run. They stay outside the default
deterministic suite and are selected by explicit markers or entry points.

## Execution tiers

| Tier / CI job | Trigger | Environment | Purpose |
| --- | --- | --- | --- |
| `platform-real-execution` | Manual `workflow_dispatch`, nightly `schedule`, published `release` | `workflow/envs/ci-fast.lock` plus Redis 7 | Redis/RQ worker, SIGALRM timeout, process-group cancellation, and tiny real Snakemake execution. |
| `real-execution` | Manual `workflow_dispatch`, nightly `schedule`, published `release` | `workflow/envs/chipseq.lock` | Complete scientific real-tool suite, including focused samtools contracts and the tiny end-to-end preprocessing run. |
| `bulk-rnaseq-real-execution` | Explicit manual `workflow_dispatch` opt-in | Protected self-hosted Linux/x64 runner with pre-staged immutable runtime, controlled fixture, Redis, and local Docker | Runtime admission, rapid quantification, full STAR+Salmon/SortMeRNA, cancellation, and timeout through canonical platform paths. |
| `container-smoke` | Manual `workflow_dispatch`, nightly `schedule`, published `release` | Docker runner image plus `profiles/default` | Build and smoke-test the runner image without pushing it. |
| Local container smoke | Developer opt-in | Docker, Apptainer, or SingularityCE plus `profiles/default` | Check mounts, Conda prefix, runner tools, and a workflow dry-run. |
| Manual HPC smoke | Site-specific opt-in | `profiles/hpc` | Validate scheduler, shadow/storage behavior, and site integration. |

The four CI real jobs are deliberately non-required. They do not run on every
pull request. Dispatch the exact PR branch when a high-risk worker, scientific,
environment, profile, or container change needs them. Nightly and published
release events use the workflow from the default branch; before a workflow
change is merged, only its manual branch dispatch can validate those new
definitions.

The nightly cron is `17 19 * * *` UTC (03:17 the next day in Asia/Shanghai).

## Concrete entry points

- `.github/workflows/ci.yml` defines `platform-real-execution`,
  `real-execution`, `bulk-rnaseq-real-execution`, and `container-smoke`.
- `test/workers/test_redis_process_integration.py` covers Redis/RQ and SIGALRM.
- `test/workers/test_tiny_execution_e2e.py` covers the API-to-worker tiny
  Snakemake lifecycle.
- `test/workers/test_cancellation_e2e.py` covers process-group cancellation.
- `test/real_execution/test_scientific_tiny_preprocessing.py` is the tiny
  scientific end-to-end run.
- `test/real_execution/test_pseudoreplicate_splitting.py` and
  `test/real_execution/test_cuttag_fragment_size.py` are focused real-tool
  contracts that require `samtools`.
- `test/bulk_rnaseq_real_execution` contains the private, explicit
  nf-core/rnaseq 3.26.0 real-execution gates.
- `scripts/stage_bulk_rnaseq_runtime_assets.py --phase verify` re-admits the
  complete pre-staged runtime without network or mutation.
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

## Bulk RNA-seq real execution

The `bulk_rnaseq_real_execution` marker is selected only when the manual
`bulk_rnaseq_real_execution` dispatch input is true. The protected runner must
provide absolute runtime, fixture, Redis, Docker executable, and Docker socket
coordinates. It must also provide the contract-pinned `/usr/bin/unshare`
2.39.3 executable and allow unprivileged user/network namespaces in the same
account that owns the worker. Admission hashes that launcher, verifies its
exact version output, and executes the JDK, Nextflow, and local Docker probe
inside a network namespace with no interface or route. Pipeline containers
remain separately constrained by `--pull=never --network=none`. The job first
runs the staging tool in read-only `verify` mode, then runs all rapid,
full-platform, cancellation, and timeout gates. Missing assets, a stale
closure, unavailable namespace isolation, a non-local Docker endpoint, any
skip/xfail, or an unclean exact checkout fails the job. OCI archives, JDKs,
references, indexes, and biological data remain outside Git.

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

PR fast, full main, platform real, scientific real, and bulk RNA-seq real are
five separate pytest entry points. Each writes a JUnit report, runs with strict
registered markers and strict xfail behavior, and must finish with zero skips
and zero xfails. Local missing-prerequisite skips are a development convenience
only; every selected real CI job is fail-closed.
