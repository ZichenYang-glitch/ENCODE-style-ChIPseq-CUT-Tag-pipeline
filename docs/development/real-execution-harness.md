# Real-Execution and Container Smoke Harness

This document defines where real-execution tests run so external tools and
heavier compute do not leak into the fast pytest suite by default.

## Execution tiers

| Tier | Trigger | Environment | Purpose |
|---|---|---|---|
| `fast-pytest` | Every PR / push | `workflow/envs/ci-fast.lock` via `.github/workflows/ci.yml` | Unit, contract, dry-run, and lint checks that finish in minutes. |
| `workflow-dispatch-real-execution` | Manual `workflow_dispatch` in GitHub Actions | `workflow/envs/chipseq.lock` | Complete `test/real_execution` suite with actual tools, including the end-to-end tiny scientific run. |
| `local-container-smoke` | Local developer, optional CI artifact | Docker / Apptainer / SingularityCE runner image + `profiles/default` | Lightweight container smoke tests, dry-runs, and environment sanity checks that need a real container runtime but not full compute. |
| `manual-hpc-smoke` | Site-specific HPC scheduler | `profiles/hpc` | Large-scale execution checks, multi-node profiles, and manual integration checks that require HPC resources. |

## Concrete entry points

- `.github/workflows/ci.yml` defines the `real-execution` job.
- `test/real_execution/test_scientific_tiny_preprocessing.py` is the current `workflow-dispatch-real-execution` test.
- `test/real_execution/test_pseudoreplicate_splitting.py` and
  `test/real_execution/test_cuttag_fragment_size.py` are focused real-tool
  contracts that require `samtools`.
- `scripts/smoke_container_runner.sh` is the local container smoke runner.
- `profiles/default/config.yaml` is the local container smoke profile.
- `profiles/hpc/config.yaml` is the manual HPC smoke profile template.
- `docs/container-usage.md` documents the runner-only container image.

## What belongs in each tier

### workflow-dispatch-real-execution
Runs real tools on generated or tiny public data. Must be triggered manually,
not on PR/push. Expected to take minutes to tens of minutes.

Run the complete local tier explicitly with:

```bash
python3 -m pytest -m real_execution test/real_execution -v
```

The default pytest selection excludes this marker.
The `chipseq` lock includes pytest as well as every executable used by this
tier. The CI job sets `HELIXWEAVE_REQUIRE_REAL_EXECUTION=1`, so a missing
executable fails the run instead of turning it into an all-skipped success.

### local-container-smoke
Validates that the runner image starts, that bind mounts and conda prefix work,
and that the workflow dry-runs inside the container. Does not run heavy compute
by default.

### manual-hpc-smoke
Site-specific validation of the HPC profile, shadow directory, storage plugins,
and scheduler integration. Requires access to a cluster.
