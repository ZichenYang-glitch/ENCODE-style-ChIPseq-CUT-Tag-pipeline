# Stage 33: Containerization Planning — Implementation Plan

> Full design: `docs/superpowers/specs/2026-05-24-stage33-containerization-planning-design.md`

**Goal:** Containerization plan for v0.3-dev. Planning only — no images built, no Dockerfile yet.

**Status:** Planning — no implementation yet.

---

### Task 1: Create design spec

- [x] `docs/superpowers/specs/2026-05-24-stage33-containerization-planning-design.md`

### Task 2: Create release checks doc

- [x] `docs/release-checks/stage33-containerization-plan.md`

### Task 3: Create tests

- [x] `test/test_stage33_containerization_plan.py`
  - Plan docs exist
  - Plan mentions Docker, Apptainer/Singularity, runner image, full-tool image tradeoff, `--use-conda`, CI scope, no builds
  - No Docker image tarballs committed

### Task 4: Documentation updates

- [x] `CHANGELOG.md` [Unreleased] with Stage 33 planning entry
- [ ] `README.md` link to containerization plan (if appropriate)

---

**Next stage (implementation):** Add `Dockerfile` and `Apptainer.def` for the runner image only.
