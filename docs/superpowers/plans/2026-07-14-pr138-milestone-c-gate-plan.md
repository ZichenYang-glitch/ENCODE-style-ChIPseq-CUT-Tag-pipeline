# PR138 Milestone C Input-to-Results Gate Implementation Plan

**Goal:** Prove one real browser-authored input reaches a canonical successful
run, persisted QC and artifacts, and an exact safe download, while preserving
real failure, cancellation, history, layout, and cleanup evidence.

**Architecture:** Reuse the PR134 sentinel-owned controlled project and the
existing real platform supervisor. Route run creation in the durable browser
gate through the PR136/137 workbench and immutable snapshot contract. Keep all
production API, lifecycle, worker, artifact, QC, and scientific contracts
unchanged.

---

### Task 1: Make the controlled inputs authoring-native

- [ ] Add failing fixture/runtime tests requiring configs to omit `samples`
  while retaining the separate TSV file.
- [ ] Add a deterministic short RUNNING observation window to success mode and
  keep every output written only after that window.
- [ ] Require runtime/demo manifests to expose separate config and TSV inputs
  without placing the TSV path inside config.
- [ ] Run fixture, runtime, adapter inline-row, and build-identity tests.

### Task 2: Join authoring and results in the real browser gate

- [ ] Replace the legacy JSON/path run helper in the durable Playwright spec
  with workbench YAML, TSV import, backend validate, opaque snapshot, snapshot-
  only create, preflight, and explicit start.
- [ ] Inspect validate/create requests and require config path privacy,
  inline-row content, immutable snapshot response, and snapshot-only create.
- [ ] Prove a structurally valid scientific rejection creates no snapshot,
  correct the draft, and recover through successful validation.
- [ ] Require observable RUNNING before SUCCEEDED, then retain exact QC,
  artifact, download, deep-link, empty, redacted-failure, and cancellation
  assertions on desktop/mobile.
- [ ] Preserve and run existing invalid YAML, sample error, stale validation,
  Back/Forward, and reload coverage.

### Task 3: Add the input-authoring demo entry and maintained roadmap

- [ ] Add failing parser/runtime tests for `--input-authoring-demo` as a
  compatible alias of the controlled results demo.
- [ ] Implement the alias without a new supervisor, queue implementation, data
  root, or fixture.
- [ ] Update local runtime documentation with the exact browser sequence,
  input file, data location, prerequisites, and Ctrl-C cleanup.
- [ ] Add a maintained delivery-status section to the roadmap that truthfully
  conditions Milestone C completion on this merged, green gate and names PR139
  as the next adapter-conformance step.

### Task 4: Complete local acceptance

- [ ] Install locked frontend dependencies in the PR138 worktree.
- [ ] Run focused Python tests and the complete Python suite with this
  worktree's `src` first on `PYTHONPATH`.
- [ ] Run real Redis/RQ/SIGALRM/cancellation/tiny Snakemake tests with zero
  skips.
- [ ] Run frontend Vitest, typecheck, production build, and Playwright desktop/
  mobile with screenshot output.
- [ ] Inspect screenshots and verify no horizontal overflow, overlap, clipping,
  empty editor, or unusable action.
- [ ] Export OpenAPI and regenerate Orval twice with zero drift; run Ruff,
  format, no-hardcoded-paths, and diff checks.
- [ ] Start `python scripts/run_local_platform.py --input-authoring-demo`, prove
  API/worker/Redis/frontend readiness, complete the real path, stop it, and
  verify every launcher-owned process is gone.

### Task 5: Publish, review, CI, and conditional merge

- [ ] Review `origin/main...HEAD`, commit with scoped messages, push, and open
  Draft PR138.
- [ ] Request the single allowed independent adversarial read-only review of
  the exact PR HEAD; fix all blocking/important findings and re-review the new
  exact HEAD if needed.
- [ ] Wait for exact-HEAD fast-checks, browser-e2e, frontend, lint, and
  lock-check to pass with the new critical paths actually executed.
- [ ] Require local HEAD, remote HEAD, reviewed HEAD, CI HEAD, latest base,
  clean worktree, scoped diff, and clean process state to match.
- [ ] Mark Ready and squash merge using `--match-head-commit`; fetch and verify
  the PR138 merge commit on `origin/main` before creating PR139.
