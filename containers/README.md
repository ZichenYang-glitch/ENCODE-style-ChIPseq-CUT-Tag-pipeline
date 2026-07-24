# Local ENCODE Snakemake Runner Container

This directory contains container definition files that produce a **runner-only**
Snakemake image for the bundled ENCODE-style epigenomics workflow. The image
includes Python, PyYAML, and Snakemake.

HelixWeave v0.3.0 does not publish this image or include it in the release
assets. These definitions are provided for optional local builds. They are not
the HelixWeave API/frontend application image and are not the pinned
Bulk RNA-seq OCI closure.

## What This Image Is

- A lightweight Snakemake runtime (~200 MB)
- Python 3.10+ with PyYAML and `snakemake-minimal >= 8.0`
- Built from `condaforge/miniforge3:24.11.3-0`

## What This Image Is NOT

- It does NOT contain bioinformatics tools (Bowtie2, samtools, MACS3, etc.).
  All tools come from `workflow/envs/*.yml` via `snakemake --use-conda`.
- It does NOT contain the workflow repository, reference data, or sample data.
  The repository must be bind-mounted at `/workspace`.
- It does NOT provision the Bulk RNA-seq Nextflow/JDK/plugin/container,
  reference, STAR/Salmon index, or SortMeRNA assets.

## Key Requirements

- **Writable `--conda-prefix` is mandatory.** The container image is read-only.
  Snakemake creates rule-specific Conda environments at runtime, which must land
  on a writable filesystem. Always use `--conda-prefix` pointing to a directory
  that exists and is writable.
- **Bind-mount the workflow repository** at `/workspace`.
- **Bind-mount data and reference directories** as needed (e.g., `/data`, `/reference`).

## Docker Usage

```bash
docker run \
    -v "$PWD":/workspace \
    -v /data:/data \
    -v /reference:/reference \
    -v /data/conda_cache:/conda_cache \
    -e HOME=/conda_cache/home \
    -e XDG_CACHE_HOME=/conda_cache/xdg-cache \
    -u $(id -u):$(id -g) \
    helixweave-encode-runner:0.3.0-local \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --cores 16 \
    --use-conda \
    --conda-prefix /conda_cache/snakemake
```

## Apptainer / Singularity Usage

```bash
apptainer exec \
    --pwd /workspace \
    --bind "$PWD":/workspace,/data:/data,/reference:/reference,/data/conda_cache:/conda_cache \
    --env XDG_CACHE_HOME=/conda_cache/xdg-cache \
    helixweave-encode-runner-0.3.0-local.sif \
    snakemake \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --cores 16 \
    --use-conda \
    --conda-prefix /conda_cache/snakemake
```

## HPC Guidance

- Set `--conda-prefix` to a scratch, project, or data directory — not your home
  directory, which is often too small for Conda environments.
- Create the cache directories before running:
  ```bash
  mkdir -p /data/conda_cache/{snakemake,home,xdg-cache}
  ```
- When using `-u` (Docker user mapping), the container has no home directory for
  the mapped UID. Snakemake/Python may try to write to `/.cache` or `~/.cache`,
  causing a `PermissionError`. Provide writable HOME and XDG_CACHE_HOME via
  environment variables:
  - `-e HOME=/conda_cache/home`
  - `-e XDG_CACHE_HOME=/conda_cache/xdg-cache`
- Apptainer runs as the invoking user by default; no `--user` flag needed.
- SingularityCE does not allow overriding `HOME` with `--env HOME=...`.
  Use `--env XDG_CACHE_HOME=/conda_cache/xdg-cache` and a writable
  `--conda-prefix`. If a custom home directory is required, use the runtime's
  `--home` option instead of `--env HOME=...`.
- Docker users should always include `-u $(id -u):$(id -g)` to match file
  ownership of bind-mounted directories.

## Full User Guide

See [`docs/container-usage.md`](../docs/container-usage.md) for the full local
ENCODE runner guide covering Docker and Apptainer/SingularityCE builds, bind
mounts, troubleshooting, and the non-publishing boundary.

## Build Instructions

Historical Docker build and smoke evidence is retained in
[`docs/release-checks/stage35-docker-runner-smoke.md`](../docs/release-checks/stage35-docker-runner-smoke.md)
for review.

Historical SingularityCE build and smoke evidence is retained in
[`docs/release-checks/stage36-singularity-runner-smoke.md`](../docs/release-checks/stage36-singularity-runner-smoke.md)
for review.

Build commands (for reference):

```bash
# Docker
docker build -f containers/Dockerfile.runner \
  -t helixweave-encode-runner:0.3.0-local .

# Apptainer
apptainer build helixweave-encode-runner-0.3.0-local.sif \
  containers/Apptainer.runner.def

# SingularityCE
singularity build helixweave-encode-runner-0.3.0-local.sif \
  containers/Apptainer.runner.def
```

## Verification

After building, verify the image works:

```bash
# Check Snakemake and PyYAML are available
snakemake --version
python -c "import yaml; print('runner ok')"
```
