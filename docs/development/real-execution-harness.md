# Real-Execution and Container Smoke Harness

This document defines where legacy and future real-execution tests run, so they
do not leak into the fast pytest suite by default.

## Execution tiers

| Tier | Trigger | Environment | Purpose |
|---|---|---|---|
| `fast-pytest` | Every PR / push | `workflow/envs/ci-fast.lock` via `.github/workflows/ci.yml` | Unit, contract, dry-run, and lint checks that finish in minutes. |
| `workflow-dispatch-real-execution` | Manual `workflow_dispatch` in GitHub Actions | `workflow/envs/chipseq.lock` | End-to-end tiny real execution with actual tools. Currently runs `test/test_stage8b_tiny_execution.py`. |
| `local-container-smoke` | Local developer, optional CI artifact | Docker / Apptainer / SingularityCE runner image + `profiles/default` | Lightweight container smoke tests, dry-runs, and environment sanity checks that need a real container runtime but not full compute. |
| `manual-hpc-smoke` | Site-specific HPC scheduler | `profiles/hpc` | Larger stress tests, multi-node profiles, and manual integration checks that require HPC resources. |
| `legacy-quarantined` | None by default | N/A | Legacy `test_stage*.py` scripts classified as `obsolete-plan-doc`, `delete-candidate`, or not yet migrated; tracked in `test/test_stage_shim.py`. |

## Legacy category mapping

- `real-execution-only` → `workflow-dispatch-real-execution`
- `manual-integration` → `local-container-smoke` or `manual-hpc-smoke` until individually migrated or retired
- `migrate-to-pytest`, `delete-candidate`, `obsolete-plan-doc` → remain outside this harness

## Concrete entry points

- `.github/workflows/ci.yml` defines the `real-execution` job.
- `test/test_stage8b_tiny_execution.py` is the current `workflow-dispatch-real-execution` test.
- `scripts/smoke_container_runner.sh` is the local container smoke runner.
- `profiles/default/config.yaml` is the local container smoke profile.
- `profiles/hpc/config.yaml` is the manual HPC smoke profile template.
- `docs/container-usage.md` documents the runner-only container image.

## What belongs in each tier

### workflow-dispatch-real-execution
Runs real tools on generated or tiny public data. Must be triggered manually,
not on PR/push. Expected to take minutes to tens of minutes.

### local-container-smoke
Validates that the runner image starts, that bind mounts and conda prefix work,
and that the workflow dry-runs inside the container. Does not run heavy compute
by default.

### manual-hpc-smoke
Site-specific validation of the HPC profile, shadow directory, storage plugins,
and scheduler integration. Requires access to a cluster.

## manual-integration scripts (14)

There are currently 14 `manual-integration` legacy scripts. They are candidates
for `local-container-smoke` or `manual-hpc-smoke` migration, or for retirement
if they only validate meta-structure.

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

## Migration checklist for a `manual-integration` script

Before moving a `manual-integration` script into any tier:

1. Decide whether it tests real pipeline behavior or only meta-structure.
2. If meta-structure, delete or migrate to a docs/native-pytest contract.
3. If real behavior, choose the smallest tier that can run it:
   - Can it run on GitHub Actions with generated data? → `workflow-dispatch-real-execution`
   - Does it only need a container runtime and a dry-run? → `local-container-smoke`
   - Does it need cluster resources? → `manual-hpc-smoke`
4. Update this doc and any guard test before changing the script's category.
