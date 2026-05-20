# Stage 10e: Environment Reliability Fix

**Date:** 2026-05-19
**Status:** implemented and locally verified

---

## 1. Problem

The original `workflow/envs/chipseq.yml` was a single monolithic Conda
environment for the full pipeline. In practice this made installation fragile:

- Bioconda `idr` requires `python <3.10`, which conflicts with modern
  `snakemake >=8`, `macs3`, and `deeptools`.
- `seacr` brings an R runtime.
- `multiqc` brings a large Python reporting stack.
- Java tools such as FastQC add a large OpenJDK runtime.

The result was a slow and brittle solve/download path before a user could even
start the workflow.

---

## 2. Decision

Use a lightweight runner environment plus Snakemake rule-specific tool
environments.

Users install only the runner first:

```bash
micromamba create -f workflow/envs/runner.yml
micromamba activate chipseq-runner
```

The full workflow is then run with `--use-conda`, allowing Snakemake to create
only the tool environments required by the selected rules.

---

## 3. Environment Layout

| File | Purpose |
|------|---------|
| `workflow/envs/runner.yml` | Lightweight entry environment: Python, PyYAML, `snakemake-minimal` |
| `workflow/envs/ci-fast.yml` | CI fast-check environment, aligned with the runner |
| `workflow/envs/fastqc.yml` | FastQC and pinned OpenJDK 17 |
| `workflow/envs/trim.yml` | Trim Galore and cutadapt |
| `workflow/envs/align.yml` | Bowtie2 and samtools alignment path |
| `workflow/envs/samtools.yml` | samtools/bedtools operations |
| `workflow/envs/macs3.yml` | MACS3 peak and bedGraph rules |
| `workflow/envs/deeptools.yml` | deepTools `bamCoverage` |
| `workflow/envs/multiqc.yml` | MultiQC reporting |
| `workflow/envs/idr.yml` | IDR isolated with `python <3.10` |
| `workflow/envs/seacr.yml` | SEACR and bedtools isolated from the core runtime |
| `workflow/envs/python.yml` | Python-only helper rules |

All env files use `conda-forge`, `bioconda`, and `nodefaults`.

---

## 4. Key Implementation Notes

- `runner.yml` and `ci-fast.yml` use `snakemake-minimal >=8,<9` to reduce the
  first install.
- IDR stays on the Bioconda package but is isolated in `idr.yml`, so its
  `python <3.10` constraint no longer conflicts with modern Python tools.
- SEACR stays in `seacr.yml`, keeping R out of the default install path.
- MultiQC stays in `multiqc.yml`, keeping reporting dependencies out of the
  runner and core runtime.
- Duplicate removal now uses `samtools markdup` by default; Picard remains an
  optional custom-runtime path if available on `PATH`.
- FastQC pins `openjdk >=17,<18` to avoid accidentally selecting newer Java
  builds.

---

## 5. Verification Snapshot

- Env YAML parse: all environment files valid.
- `chipseq-runner-smoke` actual install: PASS.
- `chipseq-runner-smoke snakemake --version`: PASS (`8.30.0`).
- Default DAG with `--use-conda` from runner env: PASS.
- Stage 2 validation stress: 15/15 PASS.
- Stage 5a IDR stress: 19/19 PASS.
- Stage 5b pseudoreplicate IDR stress: 14/14 PASS.
- Stage 7b SEACR stress: 12/12 PASS.
- Stage 8a smoke profiles from runner env: 7/7 PASS.
- Stage 8b tiny real execution: PASS using the core `chipseq` runtime.

Runner solve comparison:

| Candidate | Packages | Download |
|-----------|----------|----------|
| `snakemake >=8,<9` | 122 | ~30 MB |
| `snakemake-minimal >=8,<9` | 94 | ~6 MB |

The implemented runner uses `snakemake-minimal`.

---

## 6. Remaining Tradeoffs

- FastQC still requires Java and therefore has a large first-use rule env.
- A user running the full pipeline for the first time will still download tool
  environments as needed, but no longer pays that cost during initial setup.
- The legacy core `chipseq.yml` remains for the Stage 8b tiny real-execution
  harness and manual CI real-execution job.

---

## 7. Scope Boundaries

No schema, config contract, sample sheet contract, or pipeline feature behavior
changed. This stage only changes environment packaging, rule `conda:`
directives, and user-facing installation documentation.
