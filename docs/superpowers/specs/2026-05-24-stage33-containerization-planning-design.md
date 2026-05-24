# Stage 33: Containerization Planning — Design Spec

**Date:** 2026-05-24
**Status:** planning (no images built, no Dockerfiles yet)
**Scope:** Design a v0.3-dev containerization strategy that preserves the existing Conda-first workflow

---

## 1. Goals

1. Design a containerization approach that does NOT disrupt the existing `--use-conda` workflow.
2. Prioritize HPC users (Apptainer/Singularity) while supporting Docker for CI and portability.
3. Plan a lightweight runner image as the first deliverable; full bioinformatics image as optional future work.
4. Define what CI can test automatically vs what remains manual.
5. No image artifacts committed to the repository.

## 2. Current State

- **17 rule-specific Conda envs** under `workflow/envs/`
- **Runner env** (`runner.yml`): Python + PyYAML + Snakemake-minimal (~50 MB)
- **CI fast env** (`ci-fast.yml`): Python + PyYAML + Snakemake-minimal (~50 MB)
- **Core chipseq env** (`chipseq.yml`): full bioinformatics tools (Bowtie2, samtools, MACS3, etc.) — large solve, used only for tiny real execution
- **Snakemake `--use-conda`**: rule-specific envs created on first use, cached by Snakemake
- **CI**: fast-checks run in `ci-fast.yml`; `workflow_dispatch` runs tiny exec in `chipseq.yml`

## 3. Design Questions

### Q1: Docker vs Apptainer/Singularity — which is primary?

**Recommendation: Apptainer/Singularity as primary for HPC; Docker as secondary for CI/portability.**

| Factor | Apptainer/Singularity | Docker |
| :--- | :--- | :--- |
| HPC compatibility | Native (no daemon, user namespaces) | Requires root or rootless complex setup |
| CI compatibility | Works but less common | Standard in GitHub Actions |
| User adoption | Preferred by HPC centers | Preferred by cloud/CI users |
| Image format | `.sif` (single file) | OCI layers |

The plan should provide guidance for both, with Apptainer as the recommended path for production HPC runs.

### Q2: Containerize runner only, or full bioinformatics tools?

**Recommendation: Lightweight runner image first. Full-tool image as optional future work.**

| Option | Pros | Cons |
| :--- | :--- | :--- |
| Runner-only image | Small (~200 MB), fast build, fast CI | Still needs `--use-conda` for rule envs |
| Full-tool image | Self-contained, no conda solve at runtime | Very large (5-20 GB), slow build, hard to maintain, conflicts with existing env files |

The runner-only image packages the `chipseq-runner` env (Python + PyYAML + Snakemake) and relies on `--use-conda` for rule-specific tools. This mirrors the current local workflow and keeps `workflow/envs/*.yml` as the single source of truth.

A full-tool image is documented as a future option (v0.3+ or post-v0.3) for users who cannot run `--use-conda` on their HPC system.

### Q3: Should Snakemake still use per-rule Conda envs inside the container?

**Yes.** The runner image uses `--use-conda` for rule environments, same as the local workflow. This preserves the existing environment isolation and avoids duplicating tool dependencies in the container.

### Q4: Image size, cache, and reproducibility

- **Runner image size:** ~200-500 MB (Micromamba + runner env)
- **Conda env cache:** Snakemake caches rule envs in `--conda-prefix` (configurable)
- **Reproducibility:** Pin Snakemake version in runner env; conda envs are already version-pinned via `workflow/envs/*.yml`
- **Image rebuild:** Only when runner env dependencies change

### Q5: What should CI test automatically?

| Test | Trigger | Notes |
| :--- | :--- | :--- |
| Container build dry-run check | PR to `main` | Verify Dockerfile/recipe syntax, no actual build |
| Runner image smoke test | `workflow_dispatch` | Build image, run `snakemake --version` |
| Full pipeline in container | Manual / external | Too large for GitHub Actions runners |

### Q6: What remains manual?

- Full pipeline runs in containers on HPC or local workstations
- Public data validation in containers
- Apptainer `.sif` image building on HPC login nodes
- Full-tool image building and testing

### Q7: Interaction with public data execution reports?

Container usage should be recorded in the execution reports (`docs/release-checks/public-data-runs/*.md`) when a validation run uses a container. The report template already has a "Run Configuration" section for this.

### Q8: Runner-only container contract

**Runner-only containers are execution wrappers, not full bioinformatics runtime images.**

Key constraints:
- `--conda-prefix` is mandatory for container runs — Snakemake conda environments must land on a writable bind-mounted filesystem, not inside the read-only container image.
- The runner image does NOT contain the full workflow repo; users bind-mount the repo into `/workspace`.
- No `:latest` base image tags — Stage 34 must pin a real tag for reproducibility.
- Conda/mamba tool choice must match the selected base image (`miniforge3` → `conda`, `mambaforge` → `mamba`).
- All bioinformatics tools come from `workflow/envs/*.yml` via `--use-conda`.
- Full-tool image is deferred to v0.3+ or post-v0.3.

### Q9: Explicitly out of scope for Stage 33

- Building or committing Docker/Apptainer images
- Full bioinformatics tool container image
- Replacing `workflow/envs/*.yml` as the source of truth
- `snakemake --containerize` auto-generation (evaluated but deferred — current Snakemake 8.x containerize support is incomplete for Conda-based workflows)
- Container registry publishing (Docker Hub, GHCR, Singularity Hub)
- Multi-architecture images
- Kubernetes or cloud batch execution
