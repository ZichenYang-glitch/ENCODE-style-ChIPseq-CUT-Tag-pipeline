# HelixWeave Product Roadmap

This is the maintained product roadmap for HelixWeave. It records current
boundaries, the next product phases, and their exit criteria. Implementation
history belongs in Git history and release evidence, not in this document.

## Current product boundary

HelixWeave is a workflow-neutral platform for reproducible omics analysis on a
local workstation or within a small trusted team. It owns input authoring,
validation, durable run state, execution coordination, artifacts, QC, and the
browser experience around those concepts.

The bundled ENCODE-style ChIP-seq, CUT&Tag, ATAC-seq, and MNase-seq Snakemake
workflow is the first scientific adapter. Workflow-specific schemas, assay
policy, commands, output paths, artifact extraction, and scientific behavior
remain adapter or workflow responsibilities.

The platform is not a hosted multi-tenant service. SQLite and the local
filesystem are deliberate canonical stores; Redis/RQ provides an execution
boundary rather than independent lifecycle truth.

## Delivered capability

The current product supports:

- discovery of registered workflows and adapter-owned capabilities;
- schema-driven config and sample authoring with advanced text modes;
- structured validation issues without exposing adapter-private payloads;
- immutable, server-validated input snapshots for submitted runs;
- durable lifecycle, event, log, cancellation, and recovery state;
- Redis/RQ execution of planned Snakemake commands outside HTTP handlers;
- truthful process-group cancellation and terminal-state acknowledgement;
- adapter-owned artifact and QC extraction after successful execution;
- safe artifact listing and download without arbitrary path access;
- workflow, run, artifact, and QC views in the browser;
- a deterministic local input-to-results demonstration path; and
- reusable adapter conformance tests with a minimal test adapter.

Detailed architecture is maintained in
[`docs/architecture/platform-overview.md`](../architecture/platform-overview.md).
Scientific contracts remain in the assay, configuration, sample, output, QC,
and reproducibility references under `docs/`.

## Phase 1: Maintenance and quality baseline — current

The current phase makes existing assurance measurable and maintainable before
new product scope is added. It changes documentation and test infrastructure,
not workflow behavior or public contracts.

### Outcomes

- Separate maintained references from historical implementation process.
- Measure complete automated Python coverage, including branch coverage and
  mature subprocess collection where supported.
- Retire or migrate every legacy stage test while preserving its useful
  behavior in named pytest, integration, or real-execution gates.
- Establish PR-fast, full-main, and real-execution CI tiers.
- Enforce global and changed-lines coverage ratchets from measured baselines.
- Preserve scientific, lifecycle, migration, path-safety, API, generated-client,
  browser, and real-execution guarantees.

### Exit criteria

- Maintained documentation has no dependency on historical process plans.
- The active roadmap and architecture overview describe current truth.
- Coverage reports are reproducible locally and in CI from one configuration.
- Every legacy test classification has a final, evidenced action.
- CI tiers have distinct responsibilities without duplicate suite execution.
- A complete integration gate passes on the final stacked change set.
- Scientific workflow outputs, public API routes, persistence identity, CLI
  names, and visible product behavior remain unchanged.

## Phase 2: Second-adapter proof

After the maintenance baseline is merged, the next adapter phase will prove
that HelixWeave is genuinely workflow-neutral. Adapter selection and delivery
require a separate decision; this phase does not start as part of maintenance.

### Outcomes

- Select one bounded workflow whose owners can provide stable schemas,
  validation, workspace planning, commands, and artifact/QC extraction.
- Integrate it through the existing adapter and registry contracts.
- Keep platform services, API routes, persistence, and frontend screens free of
  adapter-private field names, rule names, and output paths.
- Demonstrate both adapters through the same authoring-to-evidence journey.
- Document capability differences without pretending unsupported features
  exist.

### Exit criteria

- The new adapter passes the reusable conformance suite.
- Existing ENCODE-style behavior and compatibility gates remain green.
- Generic API and frontend code need no workflow-specific branching.
- One deterministic small-data execution produces indexed artifacts and QC.
- Operators can distinguish adapter capabilities before submitting a run.

## Phase 3: Read-only Agent assistance

Agent work follows the second-adapter proof so assistance is designed against
platform contracts rather than one workflow's internals. The Agent remains
advisory and read-only.

### Outcomes

- Explain schemas, validation issues, run state, logs, artifacts, and QC from
  structured platform data.
- Keep generated explanations separate from recorded provenance.
- Apply redaction and path-safety rules before context reaches a model.
- Expose capability limits clearly and fail closed when context is incomplete.
- Evaluate explanations across more than one adapter.

### Exit criteria

- Agent routes depend on services and public models, not adapter internals.
- No Agent operation can submit, start, cancel, mutate, or delete a run.
- Responses do not expose secrets, private payloads, environment values, or
  workspace paths.
- Evaluation fixtures cover validation, failure diagnosis, and QC explanation
  for each supported adapter.
- Users can identify the evidence behind an explanation and its uncertainty.

## Explicit non-goals

The maintained roadmap does not currently authorize:

- authentication, authorization, multi-tenant isolation, or complex RBAC;
- Kubernetes, HPC scheduler integration, or microservice decomposition;
- PostgreSQL, object storage, or remote workspace semantics;
- Agent write actions or automatic workflow submission;
- changes to the Python distribution, import namespace, CLI names, repository
  slug, workflow identity, or artifact URI scheme;
- a frontend rewrite, server-side rendering, or a second frontend repository;
- scientific changes hidden inside platform maintenance; or
- adding a second adapter before the maintenance exit criteria are met.

## Roadmap discipline

New work should advance one product outcome and name its exit evidence. Durable
architecture or public-contract changes may require a focused decision record;
ordinary implementation detail does not. Completed checklists, temporary test
counts, commit identifiers, and branch sequencing stay out of this roadmap.
