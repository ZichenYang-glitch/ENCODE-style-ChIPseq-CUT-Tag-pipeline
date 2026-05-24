# Stage 33: Containerization Plan

**Date:** 2026-05-24
**Status:** planning (no images built, no Dockerfile yet)

## Summary

This document records the containerization strategy for v0.3-dev.
The primary goal is to provide Docker and Apptainer/Singularity options
without disrupting the existing Conda-first `--use-conda` workflow.

## Runner-Only Container Contract

**Runner-only containers are execution wrappers, not full bioinformatics
runtime images.** They require writable bind-mounted directories for
Snakemake conda environments via `--conda-prefix`.

Key constraints:
- The runner image does NOT contain the full workflow repository; users must bind-mount the repo into `/workspace`.
- Snakemake `--use-conda` creates rule-specific environments at runtime. These must land on a writable filesystem (not inside the read-only container image).
- HPC home directories are often small — always use `--conda-prefix` pointing to a larger filesystem (scratch, project, or data directory).
- The runner image provides only Python + PyYAML + Snakemake. All bioinformatics tools come from `workflow/envs/*.yml` via `--use-conda`.
- Full bioinformatics container image is deferred to v0.3+ or post-v0.3.

## Design Decisions

| Decision | Choice | Rationale |
| :--- | :--- | :--- |
| Primary target | Apptainer/Singularity | HPC compatibility, no daemon required |
| Secondary target | Docker | CI/CD and local portability |
| First deliverable | Runner-only wrapper image | Small (~200 MB), fast build, preserves `--use-conda` for rule envs |
| Full-tool image | Deferred (v0.3+ or post-v0.3) | Large (5-20 GB), complex to maintain, conflicts with env files |
| Conda env source of truth | `workflow/envs/*.yml` | Unchanged; container uses same env files via `--use-conda` |
| Image artifacts | None committed | `.sif`, `.tar`, OCI layers excluded via `.gitignore` |
| CI testing | Build dry-run on PR; smoke test on workflow_dispatch | Full pipeline runs too large for CI |
| `--conda-prefix` | Mandatory for container runs | Writable cache outside read-only image; avoids filling home dir |
| Base image tag | Pinned placeholder (Stage 34 must pin real tag) | Unpinned tags are unreproducible |
| Conda tool | Must match selected base image | `mambaforge` → `mamba`; `miniforge3` → `conda` or `mamba` |

## Apptainer/Singularity Guidance (HPC)

Planned `Apptainer.def` recipe structure (base image tag is a placeholder —
Stage 34 must replace with a real pinned tag):

```
Bootstrap: docker
From: condaforge/miniforge3:<pinned-version>

%files
    workflow/envs/runner.yml /opt/pipeline/workflow/envs/runner.yml

%post
    conda env create -f /opt/pipeline/workflow/envs/runner.yml -y || \
      mamba env create -f /opt/pipeline/workflow/envs/runner.yml -y
    conda clean -afy

%environment
    export PATH=/opt/conda/envs/chipseq-runner/bin:$PATH

%runscript
    exec snakemake "$@"
```

**Usage example (HPC):**

```bash
apptainer run \
    --pwd /workspace \
    --bind "$PWD":/workspace,/data:/data,/reference:/reference,/data/conda_cache:/conda_cache \
    chipseq-runner.sif \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --cores 16 \
    --use-conda \
    --conda-prefix /conda_cache/snakemake
```

## Docker Guidance (CI / Local)

Planned `Dockerfile` structure (base image tag is a placeholder —
Stage 34 must replace with a real pinned tag):

```dockerfile
FROM condaforge/miniforge3:<pinned-version>
COPY workflow/envs/runner.yml /opt/pipeline/workflow/envs/runner.yml
RUN conda env create -f /opt/pipeline/workflow/envs/runner.yml -y && \
    conda clean -afy
ENV PATH=/opt/conda/envs/chipseq-runner/bin:$PATH
WORKDIR /workspace
ENTRYPOINT ["snakemake"]
```

**Usage example (local Docker):**

```bash
docker run \
    -v "$PWD":/workspace \
    -v /data:/data \
    -v /reference:/reference \
    -v /data/conda_cache:/conda_cache \
    chipseq-runner \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --cores 16 \
    --use-conda \
    --conda-prefix /conda_cache/snakemake
```

## Build Smoke Checks (Stage 34)

After building the runner image, verify:

```bash
snakemake --version
python -c "import yaml; print('runner ok')"
```

## CI Integration Plan

| Check | Trigger | What it does |
| :--- | :--- | :--- |
| Dockerfile syntax check | PR | `docker build --check .` or equivalent lint |
| Apptainer recipe lint | PR | Validate `.def` syntax (no build) |
| Runner image build + smoke | `workflow_dispatch` | Build image, run `snakemake --version` and Python yaml check |
| Full pipeline in container | Manual / external | Not automated — too large for CI runners |

## What Is NOT in This Stage

- Docker/Apptainer image building or committing
- Container registry publishing
- Full bioinformatics-tool container image
- `snakemake --containerize` auto-generation
- Multi-architecture images
- Kubernetes or cloud batch execution
- Modification of `workflow/envs/*.yml`
