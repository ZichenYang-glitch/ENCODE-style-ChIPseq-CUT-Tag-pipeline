# PR142 Semantic ENCODE Config Implementation Plan

**Goal:** Replace legacy stage-number switches in the default ENCODE authoring
experience with a versioned semantic contract while preserving old YAML and
the existing Snakemake engine configuration.

**Architecture:** The ENCODE adapter translates a private validation copy from
semantic objects to legacy engine booleans. The original workflow-neutral
inputs remain unchanged for validated snapshots; only workspace
materialization emits engine keys. The schema-driven frontend consumes the
new `1.1.0` contract and the real controlled browser fixture submits semantic
fields.

**Execution note:** Implement the checklist with the test-driven-development
skill. Keep one coding agent; reviewers are read-only.

---

### Task 1: Lock schema and translation behavior with failing adapter tests

- [ ] Add schema assertions for version/IDs, strict semantic object defaults,
      scientific titles, and absence of legacy properties.
- [ ] In `test/adapters/test_encode_adapter.py`, add parameterized semantic
      boolean mapping, legacy-only, equal-alias warning, conflict, malformed
      semantic config, exact legacy bool/string/case compatibility, and invalid
      legacy value tests using the existing coercion semantics.
- [ ] Run the focused tests and confirm they fail because `1.1.0` and the
      translation boundary do not yet exist.
- [ ] Implement the smallest adapter-owned translator and warning propagation.
- [ ] Re-run focused tests to green and refactor only duplicated pair handling.

### Task 2: Prove snapshot and workspace separation

- [ ] Add failing service/API tests requiring a validated snapshot to retain
      the semantic config without injected engine fields and one redacted
      warning code in snapshot evidence / HTTP response for equal aliases.
- [ ] Add failing workspace tests for the four semantic boolean combinations,
      exact private engine booleans, absence of semantic keys, and scientific
      validator round-trip behavior.
- [ ] Make adapter and `WorkspacePlanner` results preserve exactly one
      controlled validation warning while rendering the existing engine config
      and retaining one planning-complete info issue.
- [ ] Prove `chipseq_idr.enabled=true` remains orthogonal to explicit advanced
      `reproducibility.idr.chipseq_narrow=false` under existing scientific
      behavior.
- [ ] Run adapter/service/workspace focused tests to green.

### Task 3: Upgrade the schema-driven frontend contract

- [ ] Add failing Vitest assertions that `1.1.0` is supported, `1.0.0` is
      rejected by the authoring workbench, semantic defaults survive
      Form/YAML/Review, and default visible/review content has no legacy keys.
- [ ] Add a `ValidatedSubmission` route test requiring a successful backend
      warning to retain its severity, show non-error status styling, and leave
      the snapshot/Create action usable.
- [ ] Configure RJSF to apply defaults only to the explicitly initialized
      draft, then prove opposite legacy switches survive a known Form edit,
      YAML, Review, and validation without injected semantic aliases.
- [ ] Update the workbench contract fixture and only the tests whose schema is
      the current ENCODE authoring contract; retain historical snapshot and
      test-adapter versions where appropriate.
- [ ] Update the controlled fixture to expose semantic config and add Python
      tests for its safe conversion.
- [ ] Run focused frontend and fixture tests to green.

### Task 4: Extend the real browser proof

- [ ] Change the backend-rejection case to the semantic dependency conflict.
- [ ] Assert exact semantic labels/default checkboxes in the default RJSF form,
      exact semantic booleans in the validate payload, and no legacy stage
      identifiers in default page content, Review, request, or screenshots.
- [ ] Capture the default Config form at 1440x900, 1024x768, 390x844, and
      360x800 and require no page-level horizontal overflow.
- [ ] Prove author -> validate -> snapshot -> create -> start -> RUNNING ->
      SUCCEEDED with the semantic payload at desktop and mobile widths.
- [ ] Run both Playwright projects with the real API, SQLite, Redis/RQ worker,
      and controlled Snakemake workflow; inspect four viewport screenshots.

### Task 5: Complete repository-wide gates

- [ ] Run focused and full Python tests with `PYTHONPATH=src` and the locked
      Python 3.12/Snakemake environment.
- [ ] Run the real Redis/RQ/SIGALRM/cancellation/tiny execution gate with no
      skips relevant to durable execution.
- [ ] Run frontend Vitest, TypeScript typecheck, and production build.
- [ ] Export OpenAPI, mechanically regenerate Orval, and require zero drift.
- [ ] Run Ruff check, Ruff format check, and `git diff --check`.
- [ ] Confirm no Redis, worker, horse, Snakemake, Vite, or fixture helper owned
      by this worktree remains.

### Task 6: Review, publish, and verify exact CI HEAD

- [ ] Request two independent read-only merge-gate reviews: adapter contract /
      compatibility and frontend/E2E truthfulness.
- [ ] Fix every blocking or important finding with a failing regression test,
      rerun affected and full gates, and repeat review if needed.
- [ ] Commit with a single-scope message and no attribution trailers, push
      `codex/pr142-semantic-encode-config`, and create Draft PR142.
- [ ] Verify local HEAD, remote HEAD, reviewed HEAD, and CI HEAD are identical;
      wait for every required check to pass.
- [ ] Keep the PR Draft and report the contract, mapping, conflict behavior,
      real E2E evidence, CI links, and residual limits without merging.
