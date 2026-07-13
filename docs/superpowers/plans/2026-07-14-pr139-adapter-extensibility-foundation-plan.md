# PR139 Adapter Extensibility Foundation Implementation Plan

**Goal:** Add a reusable adapter conformance gate, make capability declarations
fail closed, and prove the existing workflow-neutral service seams with ENCODE
and a test-only minimal adapter without implementing a second real workflow.

## Task 1: Lock static capability declarations

- [x] Add failing domain tests for duplicate, malformed, and unknown capability
  names.
- [x] Add failing registry tests for malformed capability objects and both QC
  capability/protocol mismatch directions.
- [x] Introduce the canonical capability constants and validation while keeping
  the public serialized shape unchanged.
- [x] Make workspace planning fail before adapter invocation when the capability
  is absent, with a stable platform issue.

## Task 2: Add the reusable conformance runner

- [x] Add failing tests for a packaged, pytest-independent conformance API.
- [x] Implement `AdapterConformanceCase`, `AdapterConformanceError`, and
  `verify_adapter_conformance` under `encode_pipeline.testing`.
- [x] Keep diagnostics bounded to contract coordinates; do not include fixture
  contents, local paths, or arbitrary exception strings.
- [x] Cover metadata, schema, validation, DAG, workspace, command, artifact, QC,
  and declared/unsupported capability behavior.

## Task 3: Prove two adapters against the same suite

- [x] Add a test-only minimal adapter with deterministic workflow-neutral
  fixtures and every canonical capability.
- [x] Run it through the conformance runner and existing workflow-info,
  validation, workspace-planning, artifact-indexing, and QC-indexing services
  without production workflow-ID branches.
- [x] Add an ENCODE conformance case using the existing tiny profile,
  adapter-planned workspace, and bounded QC source fixture.
- [x] Retain focused ENCODE scientific tests as the authority for detailed
  manifest and QC semantics.

## Task 4: Document the extension boundary

- [x] Record why standard entry points are deferred until the first real
  external package and why no custom loader is introduced.
- [x] Document bundled command/build-identity composition as a residual limit,
  not a completed generic-runtime claim.
- [x] Update the maintained roadmap to mark the conformance foundation and keep
  the real-adapter order RNA-seq then Hi-TrAC/TracPre2.

## Task 5: Verify and deliver

- [x] Run focused conformance, adapter, registry, and service tests.
- [x] Run the complete Python suite and real Redis/RQ/Snakemake durable gate.
- [x] Run frontend tests, typecheck, production build, and Playwright gates.
- [x] Confirm OpenAPI export and Orval generation have zero drift.
- [ ] Run Ruff, format, hardcoded-path, diff, clean-worktree, and process-residue
  checks.
- [ ] Commit and push a Draft PR, then perform one exact-HEAD independent
  adversarial read-only review.
- [ ] Fix every blocking/important finding, rerun affected gates, and require
  exact CI green.
- [ ] Under the user's authorization, mark Ready and squash merge with
  `--match-head-commit`; fetch and verify the merge in `origin/main`.
