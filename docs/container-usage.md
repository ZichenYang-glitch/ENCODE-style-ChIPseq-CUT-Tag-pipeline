# Local ENCODE Runner Container Guide

This document explains how to build and run the bundled ENCODE-style
epigenomics workflow with the repository's runner-only Docker and
Apptainer/SingularityCE definitions.

HelixWeave v0.3.0 does not publish a container image. These definitions are
optional local build inputs, not HelixWeave application release assets. They
do not provide the Bulk RNA-seq runtime, whose admitted nf-core/rnaseq 3.26.0
OCI closure remains operator-owned.

## Runner-Only Image Contract

The runner image is a **lightweight Snakemake runtime** (~320 MB), not a full
bioinformatics environment. It contains:

- Python 3.10+
- PyYAML
- Snakemake-minimal >= 8.0
- `conda` (for `--use-conda` rule environments)

It does **NOT** contain: Bowtie2, samtools, MACS3, FastQC, deepTools, Picard,
or any other bioinformatics tools. All tools come from `workflow/envs/*.yml`
via `snakemake --use-conda`.

The runner image does **NOT** contain the workflow repository, reference data,
or sample data. The repository must be bind-mounted at `/workspace`.
It also does **NOT** contain Nextflow, the Bulk RNA-seq JDK/plugin/runtime
closure, STAR/Salmon indexes, SortMeRNA assets, or Bulk RNA-seq references.

**`--conda-prefix` is mandatory.** The container image is read-only. Snakemake
creates rule-specific Conda environments at runtime, which must land on a
writable bind-mounted filesystem.

## Building the Image

### Docker

```bash
docker build -f containers/Dockerfile.runner \
  -t helixweave-encode-runner:0.3.0-local .
```

### Apptainer / SingularityCE

```bash
# Apptainer
apptainer build helixweave-encode-runner-0.3.0-local.sif \
  containers/Apptainer.runner.def

# SingularityCE
singularity build helixweave-encode-runner-0.3.0-local.sif \
  containers/Apptainer.runner.def
```

Historical local build verification is retained in the
[Docker smoke report](release-checks/stage35-docker-runner-smoke.md) and
[SingularityCE smoke report](release-checks/stage36-singularity-runner-smoke.md).

## Docker Usage

### Prepare cache directories

```bash
mkdir -p /path/to/conda_cache/{snakemake,home,xdg-cache}
```

### Run

```bash
docker run --rm \
    -v "$PWD":/workspace \
    -v /data:/data \
    -v /reference:/reference \
    -v /path/to/conda_cache:/conda_cache \
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

### User mapping

When using `-u $(id -u):$(id -g)`, the container maps your local UID/GID.
This ensures output file ownership matches your host user. However, the
container has no home directory for the mapped UID, so Snakemake/Python
may try to write to `/.cache` or `~/.cache`, causing a `PermissionError`.

The fix: provide writable `HOME` and `XDG_CACHE_HOME` via environment
variables pointing to bind-mounted cache directories.

## Apptainer / SingularityCE Usage

### Prepare cache directories

```bash
mkdir -p /path/to/conda_cache/{snakemake,xdg-cache}
```

### Apptainer

```bash
apptainer exec \
    --pwd /workspace \
    --bind "$PWD":/workspace,/data:/data,/reference:/reference,/path/to/conda_cache:/conda_cache \
    --env XDG_CACHE_HOME=/conda_cache/xdg-cache \
    helixweave-encode-runner-0.3.0-local.sif \
    snakemake \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --cores 16 \
    --use-conda \
    --conda-prefix /conda_cache/snakemake
```

### SingularityCE

```bash
singularity exec \
    --pwd /workspace \
    --bind "$PWD":/workspace,/data:/data,/reference:/reference,/path/to/conda_cache:/conda_cache \
    --env XDG_CACHE_HOME=/conda_cache/xdg-cache \
    helixweave-encode-runner-0.3.0-local.sif \
    snakemake \
    -s workflow/Snakefile \
    --configfile config/config.yaml \
    --cores 16 \
    --use-conda \
    --conda-prefix /conda_cache/snakemake
```

### SingularityCE HOME warning

SingularityCE does not allow overriding `HOME` with `--env HOME=...`.
Use `--env XDG_CACHE_HOME=/conda_cache/xdg-cache` and a writable
`--conda-prefix`. If a custom home directory is required, use the
runtime's `--home` option instead of `--env HOME=...`.

## HPC Guidance

- Bind-mount scratch/project directories for both data and conda cache —
  home directories are often too small.
- Create all cache directories before running.
- Apptainer/Singularity runs as the invoking user by default; no `--user`
  flag needed.
- Docker users should always include `-u $(id -u):$(id -g)`.

## Required Bind Mounts

| Mount | Purpose | Required |
| :--- | :--- | :--- |
| repo → `/workspace` | Workflow files (Snakefile, rules, configs) | Yes |
| data → `/data` | FASTQ files, output directory | For real runs |
| reference → `/reference` | Bowtie2 index, chrom_sizes, blacklist, GTF, FASTA | For real runs |
| conda_cache → `/conda_cache` | Snakemake `--conda-prefix`, cache directories | Yes |

## Troubleshooting

### Docker Hub pull timeout or proxy errors

1. Verify Docker daemon is running: `docker run --rm hello-world`
2. Check network connectivity to Docker Hub: `docker pull alpine:latest`
3. If behind a proxy, configure Docker daemon proxy settings.

### PermissionError: `/.cache` or `~/.cache`

Root cause: Docker `-u` maps your host UID but the container has no
home directory for that UID.

Fix: Add `-e HOME=/conda_cache/home -e XDG_CACHE_HOME=/conda_cache/xdg-cache`
and ensure those directories exist in the bind-mounted conda cache.

### Missing `samples.tsv` or config

The container's working directory is `/workspace`. Ensure:
1. The repo is correctly bind-mounted at `/workspace`.
2. `--configfile` paths are relative to `/workspace`, or use absolute
   paths like `--configfile /workspace/config/config.yaml`.
3. `samples.tsv` paths are accessible from inside the container — use
   absolute paths in the sample sheet or ensure relative paths resolve
   correctly from the workflow directory (`--directory`).

### Conda cache permissions

The `--conda-prefix` directory must be writable by the container user.
If using Docker with `-u`, ensure the bind-mounted cache directory is
writable by your host UID: `chmod 755 /path/to/conda_cache`.

## Smoke Testing

A smoke test script is available at [`scripts/smoke_container_runner.sh`](../scripts/smoke_container_runner.sh).
It creates a temporary workspace under `/tmp`, copies a test profile,
creates placeholder FASTQs, and runs a bind-mounted dry-run. Usage:

```bash
# Docker
bash scripts/smoke_container_runner.sh \
  docker helixweave-encode-runner:0.3.0-local

# SingularityCE
bash scripts/smoke_container_runner.sh \
  singularity helixweave-encode-runner-0.3.0-local.sif
```

The script does NOT download public data and does NOT write outputs into
the repository.

## v0.3.0 release boundary

No Docker, GHCR, Apptainer, Singularity, or other container artifact is part of
the HelixWeave v0.3.0 release contract. Do not interpret the local
`0.3.0-local` examples as published tags.

Build the optional ENCODE runner locally from `Dockerfile.runner` or
`Apptainer.runner.def`. Publishing an application image, this runner image, or
any Bulk RNA-seq runtime image requires a separate versioning, provenance,
security, and release decision.
