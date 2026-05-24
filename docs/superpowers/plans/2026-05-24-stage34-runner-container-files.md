# Stage 34: Runner-Only Container Files — Implementation Plan

> Full design: `docs/superpowers/specs/2026-05-24-stage34-runner-container-files-design.md`

**Goal:** Static container definition files only — no build, no push, no data downloads.
**Status:** Implemented — 18/18 tests pass; static files complete; build deferred to Stage 35.

---

### Task 1: Create `containers/` directory and files

- [x] Create `containers/Dockerfile.runner`:
  ```dockerfile
  FROM condaforge/miniforge3:24.11.3-0
  COPY workflow/envs/runner.yml /opt/pipeline/workflow/envs/runner.yml
  RUN conda env create -f /opt/pipeline/workflow/envs/runner.yml -y && \
      conda clean -afy
  ENV PATH=/opt/conda/envs/chipseq-runner/bin:$PATH
  WORKDIR /workspace
  ENTRYPOINT ["snakemake"]
  ```

- [x] Create `containers/Apptainer.runner.def` — identical base, `conda env create`, PATH export
- [x] Create `containers/README.md` — Docker/Apptainer examples, HPC guidance, build deferred

### Task 2: Create tests

- [x] `test/test_stage34_runner_container_files.py` — 18 assertions from design spec.

### Task 3: Verify

- [x] 18/18 PASS: `python3 test/test_stage34_runner_container_files.py`
- [x] PASS: `python3 test/test_no_hardcoded_paths.py`

### Task 4: Documentation

- [x] `CHANGELOG.md` [Unreleased] entry
