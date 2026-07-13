# PR136 Schema-Driven Input Workbench Implementation Plan

**Goal:** Render the PR135 adapter authoring contract as a responsive,
memory-only Config/Samples/Options/Review workbench without validating or
submitting a run.

**Architecture:** A new generated-operation-backed route owns one reducer for
the entire draft. RJSF/AJV 2020 renders config and options, YAML is a guarded
transient view over the canonical config object, Papa Parse atomically imports
sample strings, TanStack Table projects editable rows, and URL search params
select the active view. The backend remains the scientific authority.

**Tech stack:** React 18, React Router, TanStack Query/Table, RJSF, AJV 2020,
YAML, Papa Parse, CodeMirror 6, Radix Tabs, Lucide, Vitest/Testing Library,
Playwright, Vite.

---

### Task 1: Lock dependencies and pure authoring primitives

**Files:**

- Modify: `frontend/package.json`
- Modify mechanically: `frontend/package-lock.json`
- Create: `frontend/src/features/input-workbench/schemaContract.ts`
- Create: `frontend/src/features/input-workbench/yamlDraft.ts`
- Create: `frontend/src/features/input-workbench/sampleTsv.ts`
- Create: `frontend/src/features/input-workbench/draft.ts`
- Create focused unit tests beside those modules

- [x] Install one compatible RJSF/AJV set, YAML, Papa Parse, TanStack Table,
  CodeMirror 6 wrapper/language/view, and Radix Tabs.
- [x] Write failing compatibility tests for real PR135 Draft 2020-12 schemas,
  version/dialect/mode checks, exact 1.0.0 ceilings, and fail-closed malformed
  contracts.
- [x] Implement the Ajv2020 RJSF validator and schema projection helpers.
- [x] Write failing YAML tests for valid object updates, invalid-buffer
  preservation, duplicate keys, multi-doc/custom tags, aliases, non-object and
  non-finite values, non-string keys, unsafe integers, and controlled
  line/column issues.
- [x] Implement parse/stringify helpers and canonical/transient transitions.
- [x] Write failing Papa Parse tests for quoting, CRLF, empty strings,
  duplicate/missing/unknown headers, empty/ragged files, atomic failure, and
  byte/row/column/header/cell bounds.
- [x] Implement schema-ordered columns, atomic import, client row IDs, and
  deterministic JSON Review serialization.

### Task 2: Workbench reducer, schema query, and URL route

**Files:**

- Create: `frontend/src/features/input-workbench/useInputDraft.ts`
- Create: `frontend/src/routes/workflows/new-run.tsx`
- Modify: `frontend/src/app/router.tsx`
- Modify minimally: workflow detail navigation
- Create route/component tests

- [x] Write failing tests for schema loading, schema null, network failure,
  retry, unsupported contracts, URL deep links, and controlled status text.
- [x] Call the generated `getWorkflowSchema` operation through TanStack Query;
  do not add a DTO or compatibility client.
- [x] Initialize canonical config/options from RJSF default form state and
  sample columns from adapter schema.
- [x] Implement one reducer retained across step/mode changes and history
  navigation; hard reload intentionally creates a fresh draft.
- [x] Add the nested route and a compact entry action from workflow detail.
- [x] Prove no validate/create/start operation is imported or called.

### Task 3: Config and options editors

**Files:**

- Create: `frontend/src/features/input-workbench/SchemaObjectForm.tsx`
- Create: `frontend/src/features/input-workbench/ConfigEditor.tsx`
- Create: `frontend/src/features/input-workbench/OptionsEditor.tsx`
- Add scoped styles and component tests

- [x] Write failing tests for visible schema/default fields, controlled edits,
  and preservation of unknown advanced fields.
- [x] Render RJSF with the Ajv2020 validator and no implicit submit control.
- [x] Add URL-driven Form/YAML controls and a labeled CodeMirror editor.
- [x] Block mode switch/Review while YAML is invalid and retain the buffer.
- [x] Show partial-coverage and memory-only draft notices without exposing
  untrusted exception details.

### Task 4: Samples table and deterministic Review

**Files:**

- Create: `frontend/src/features/input-workbench/SampleEditor.tsx`
- Create: `frontend/src/features/input-workbench/DraftReview.tsx`
- Add component tests

- [x] Write failing tests for upload, import failure preserving old rows,
  edit/add/delete, enum/string semantics, stable schema column order, and
  desktop/mobile projections.
- [x] Implement bounded file selection and atomic TSV import.
- [x] Ensure a slower earlier file read cannot replace a newer TSV selection.
- [x] Render a compact TanStack desktop table and vertical mobile records from
  the same controlled rows.
- [x] Render deterministic config/samples/options JSON and structural issues;
  include no submission action or scientific-success claim.

### Task 5: Responsive real browser coverage

**Files:**

- Modify: `frontend/e2e/durable-run.spec.ts` or create a focused workbench spec
- Modify styles only as needed for the new workbench

- [x] Add a real-API desktop flow: schema load, TSV upload/edit, advanced YAML,
  Form round-trip, Review, history, screenshot, and no horizontal overflow.
- [x] Add the equivalent mobile route/editor/table flow and screenshot check.
- [x] Verify a hard reload restores the URL step but starts an explicitly fresh
  in-memory draft.
- [x] Run the existing durable results and cancellation browser gates to prove
  the new route did not regress them.

### Task 6: Full gates, independent review, and Draft PR

- [x] Run all frontend Vitest tests, TypeScript typecheck, and production build.
- [x] Run desktop/mobile Playwright and inspect screenshots at required sizes.
- [x] Run Python full regression and applicable real Redis/RQ/Snakemake gates.
- [x] Export OpenAPI and regenerate Orval; require a zero diff in generated
  contract files.
- [x] Run changed-file Ruff/format applicability checks, `git diff --check`,
  and residual-process checks.
- [x] Add merge-gate regressions for safe-integer Form input, shared manual/TSV
  sample transport validation, sample-schema ceilings, and invalid-YAML
  history without stale Review serialization.
- [x] Review the complete `origin/main...HEAD` diff for state truthfulness,
  schema ownership, sensitive-data persistence, responsive overflow, and scope
  drift.
- [ ] Request an independent read-only merge-gate review and fix all blocking
  or important findings.
- [ ] Commit without attribution, push, create Draft PR136, wait for all
  required checks to pass, and keep it unmerged.
