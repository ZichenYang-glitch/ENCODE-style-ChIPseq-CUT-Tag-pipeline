# Runner-Only Container

This directory contains container definition files that produce a **runner-only**
Snakemake image. The image includes Python + PyYAML + Snakemake and can execute
any Snakemake workflow.

## What This Image Is

- A lightweight Snakemake runtime (~200 MB)
- Python 3.10+ with PyYAML and `snakemake-minimal >= 8.0`
- Built from `condaforge/miniforge3:24.11.3-0`

## What This Image Is NOT

- It does NOT contain bioinformatics tools (Bowtie2, samtools, MACS3, etc.).
  All tools come from `workflow/envs/*.yml` via `snakemake --use-conda`.
- It does NOT contain the workflow repository, reference data, or sample data.
  The repository must be bind-mounted at `/workspace`.

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
    -u $(id -u):$(id -g) \
    chipseq-runner \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --cores 16 \
    --use-conda \
    --conda-prefix /conda_cache/snakemake
```

## Apptainer / Singularity Usage

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

## HPC Guidance

- Set `--conda-prefix` to a scratch, project, or data directory — not your home
  directory, which is often too small for Conda environments.
- Create the cache directory before running: `mkdir -p /data/conda_cache`.
- Apptainer runs as the invoking user by default; no `--user` flag needed.
- Docker users should always include `-u $(id -u):$(id -g)` to match file
  ownership of bind-mounted directories.

## Build Instructions

Building and testing the image is deferred to Stage 35. The Dockerfile and
Apptainer recipe in this directory are static files that can be reviewed and
built manually.

Expected build commands (Stage 35+):

```bash
# Docker
docker build -f containers/Dockerfile.runner -t chipseq-runner .

# Apptainer
apptainer build chipseq-runner.sif containers/Apptainer.runner.def
```

## Verification

After building, verify the image works:

```bash
# Check Snakemake and PyYAML are available
snakemake --version
python -c "import yaml; print('runner ok')"
```
