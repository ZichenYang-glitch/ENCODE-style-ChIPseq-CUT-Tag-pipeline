# Stage 36: Singularity Runner Build Smoke Report

**Date:** 2026-05-25
**Status:** completed (SingularityCE smoke)
**Image artifact:** `chipseq-runner-stage36.sif` (local only, ignored, not committed)

---

## 1. Runtime Verification

```bash
$ singularity --version
singularity-ce version 4.1.1
```

The Ubuntu 24.04 `singularity-container` package provides the
SingularityCE command-line runtime. This satisfies the Apptainer/Singularity
runner validation target for HPC-style execution.

---

## 2. Build

The image was built from `containers/Apptainer.runner.def`:

```bash
export http_proxy=http://127.0.0.1:10808
export https_proxy=http://127.0.0.1:10808
export HTTP_PROXY=http://127.0.0.1:10808
export HTTPS_PROXY=http://127.0.0.1:10808

sudo -E singularity build chipseq-runner-stage36.sif containers/Apptainer.runner.def
```

The recipe uses the pinned base image
`condaforge/miniforge3:24.11.3-0`, matching the Docker runner image.

---

## 3. Smoke Checks

### 3.1 Snakemake version

```bash
$ singularity exec chipseq-runner-stage36.sif snakemake --version
8.30.0
```

### 3.2 PyYAML import

```bash
$ singularity exec chipseq-runner-stage36.sif python -c "import yaml; print('runner ok')"
runner ok
```

Both smoke checks pass: Snakemake 8.30.0 and PyYAML are available.

---

## 4. Bind-Mounted Dry-Run

A dry-run was executed with the workflow repository mounted at `/workspace`, a
temporary smoke workspace mounted at `/smoke`, and a writable Conda cache mounted
at `/conda_cache`:

```bash
rm -rf /tmp/chipseq-singularity-smoke /tmp/chipseq-singularity-conda-cache
mkdir -p /tmp/chipseq-singularity-smoke
mkdir -p /tmp/chipseq-singularity-conda-cache/snakemake
mkdir -p /tmp/chipseq-singularity-conda-cache/xdg-cache
mkdir -p /tmp/chipseq-singularity-conda-cache/home

cp test/profiles/chipseq_pe_noctrl/config.yaml /tmp/chipseq-singularity-smoke/config.yaml
cp test/profiles/chipseq_pe_noctrl/samples.tsv /tmp/chipseq-singularity-smoke/samples.tsv
touch /tmp/chipseq-singularity-smoke/R1.fq /tmp/chipseq-singularity-smoke/R2.fq

singularity exec \
    --bind "$PWD":/workspace \
    --bind /tmp/chipseq-singularity-smoke:/smoke \
    --bind /tmp/chipseq-singularity-conda-cache:/conda_cache \
    --env XDG_CACHE_HOME=/conda_cache/xdg-cache \
    chipseq-runner-stage36.sif \
    snakemake \
        -s /workspace/workflow/Snakefile \
        --directory /smoke \
        --configfile /smoke/config.yaml \
        --cores 1 \
        --use-conda \
        --conda-prefix /conda_cache/snakemake \
        -n
```

Dry-run completed without errors. The DAG resolved to 14 jobs:

```text
job                        count
-----------------------  -------
all                            1
bamcoverage                    1
bowtie2_align                  1
duplicate_handling             1
fastqc                         1
macs3_callpeak                 1
pipeline_done                  1
samtools_filter                1
samtools_flagstat              1
samtools_flagstat_final        1
samtools_idxstats              1
samtools_index_filt            1
samtools_index_sorted          1
trim_galore                    1
total                         14
```

The command ended with:

```text
This was a dry-run (flag -n).
```

---

## 5. HOME Environment Note

An earlier command included:

```bash
--env HOME=/conda_cache/home
```

SingularityCE reported:

```text
WARNING: Overriding HOME environment variable with SINGULARITYENV_HOME is not permitted
```

This warning is non-fatal, and the dry-run still passed. Unlike the Docker
`-u $(id -u):$(id -g)` case, Singularity runs as the invoking user by default and
does not need the Docker-specific `HOME` override. For Singularity/Apptainer,
use a writable `--conda-prefix` and, if needed, `XDG_CACHE_HOME`. If a custom home
directory is required, use the runtime's `--home` option rather than
`--env HOME=...`.

---

## 6. Observations

| Item | Result |
| :--- | :--- |
| Runtime | PASS - `singularity-ce version 4.1.1` |
| Build | PASS - `containers/Apptainer.runner.def` builds a local `.sif` |
| Snakemake version | PASS - 8.30.0 |
| PyYAML import | PASS - runner ok |
| Bind-mount dry-run | PASS - 14-job DAG |
| `--use-conda` inside Singularity | Works in dry-run planning mode |
| `--conda-prefix` inside Singularity | Works with writable bind mount |
| `--env HOME` | Not needed; SingularityCE warns if overridden this way |

## 7. Artifact Policy

- `chipseq-runner-stage36.sif` is a local build artifact only.
- `.sif` files are ignored by `.gitignore`.
- No `.sif`, `.tar`, `.oci`, or other image artifact is committed.
- No image pushed to a registry.

## 8. Remaining Scope

- Full real execution inside Singularity is not tested yet.
- Public data execution remains manual and external.
- Full-tool image remains deferred; runner-only image continues to rely on
  Snakemake `--use-conda` and `workflow/envs/*.yml`.
