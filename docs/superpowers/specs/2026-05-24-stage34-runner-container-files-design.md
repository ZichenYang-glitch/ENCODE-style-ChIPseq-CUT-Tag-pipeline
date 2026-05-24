# Stage 34: Runner-Only Container Files — Design Spec

**Date:** 2026-05-24
**Status:** implemented (static files complete; 18/18 tests pass; build deferred to Stage 35)
**Scope:** Static container definition files only — no build, no push, no data downloads
**Prerequisite:** Stage 33 containerization plan (complete)

---

## 1. Goals

1. Create `Dockerfile.runner` and `Apptainer.runner.def` that produce a runner-only Snakemake image.
2. Create `containers/README.md` documenting usage, bind mounts, `--conda-prefix`, and HPC guidance.
3. Create static-file tests that validate container file content and contract.
4. No Docker build, no Apptainer build, no image artifacts.
5. Defer actual build/smoke execution to Stage 35.

## 2. Design Decisions

### Base image: `condaforge/miniforge3:24.11.3-0`

- Verified real tag from Docker Hub (active, last pushed 2025-03-01).
- `miniforge3` ships `conda` and `mamba`; no `micromamba` binary.
- Both `Dockerfile.runner` and `Apptainer.runner.def` MUST use the identical `FROM` string.
- No `:latest`, no `<pinned-version>` placeholder.

### Env creation: `conda env create`

- Both files use `conda env create -f /opt/pipeline/workflow/envs/runner.yml`.
- No `|| mamba ...` fallback — masks real failures.
- `conda clean -afy` after env creation to keep image small.

### COPY scope: `workflow/envs/runner.yml` only

- Explicit `COPY workflow/envs/runner.yml /opt/pipeline/workflow/envs/runner.yml`.
- No `COPY . .` — defense-in-depth even with `.dockerignore`.
- The workflow repo is bind-mounted at `/workspace` at runtime.
- The runner image does NOT contain scripts, tests, configs, or reference data.

### Stage 34 is static-file-only

- No `docker build` or `apptainer build` execution.
- Tests validate file existence, content, and contract — not build output.
- Stage 35 will add build verification and smoke tests.

### Security: container executes mounted workflow code

The container executes the mounted workflow code. Users should only run trusted
workflow repositories and configs. This is equivalent to running Snakemake locally
on a trusted repo. Bind-mounting the repo rather than copying it into the image
means the user always runs the code they can inspect on their filesystem.

## 3. Files

| File | Purpose |
|------|---------|
| `containers/Dockerfile.runner` | Docker runner image definition |
| `containers/Apptainer.runner.def` | Apptainer/Singularity runner image definition |
| `containers/README.md` | Usage guide with Docker and Apptainer examples |
| `test/test_stage34_runner_container_files.py` | Static-file contract tests |

## 4. Contract Tests

| # | Assertion |
|---|-----------|
| 1 | `containers/Dockerfile.runner` exists |
| 2 | `containers/Apptainer.runner.def` exists |
| 3 | `containers/README.md` exists |
| 4 | Both files use identical `FROM condaforge/miniforge3:24.11.3-0` |
| 5 | Neither file contains `:latest` or `<pinned-version>` |
| 6 | Dockerfile does NOT contain `COPY . .` (defense-in-depth) |
| 7 | Dockerfile copies `workflow/envs/runner.yml` |
| 8 | Apptainer.def copies `workflow/envs/runner.yml` |
| 9 | README mentions `--conda-prefix` |
| 10 | README mentions `-u $(id -u):$(id -g)` (Docker user mapping) |
| 11 | README mentions `--pwd /workspace` (Apptainer) |
| 12 | README mentions bind-mount paths (`/workspace`, `/data`, `/reference`, `/conda_cache`) |
| 13 | No `wget`/`curl`/`fastq-dump` in any container file |
| 14 | `.dockerignore` excludes `.sif`, `.tar`, data files |

## 5. README Coverage

| Section | Content |
| :--- | :--- |
| What this image is | Snakemake runtime only; bioinformatics tools from `--use-conda` |
| What this image is NOT | Full bioinformatics; does not contain Bowtie2, samtools, MACS3 |
| Docker usage example | `-v` bind mounts, `-u $(id -u):$(id -g)`, `--conda-prefix` |
| Apptainer usage example | `--bind`, `--pwd /workspace`, `--conda-prefix` |
| HPC guidance | Set `--conda-prefix` to scratch/project directory, not home |
| Build instructions | Refer to Stage 35 (not part of Stage 34) |
| Verification | `snakemake --version`, `python -c "import yaml; print('ok')"` |

## 6. Out of Scope

- `docker build` or `apptainer build` execution
- Image registry push
- Full bioinformatics tool image
- Multi-architecture builds
- CI workflow changes
- `snakemake --containerize` integration
- Conda lockfile reproducibility (future work)
