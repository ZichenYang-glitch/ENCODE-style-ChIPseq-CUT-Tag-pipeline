# Workflow Platform and Agent-Ready Adapter Roadmap

This roadmap updates the project direction from a single fullstack application
around this pipeline to a workflow platform with agent-ready workflow adapters.
It is a planning document only. It does not implement the platform, services,
API routes, frontend, execution workers, or agent tools.

## Maintained delivery status

The architecture principles below remain authoritative, but their original
phase list predates the implemented platform. The maintained delivery sequence
is now:

- PR123–126 delivered durable local execution: SQLite lifecycle truth, an
  independent Redis/RQ worker, real Snakemake execution and process-group
  cancellation, and the browser run workbench.
- PR127–134 delivered results visibility: adapter-owned artifact and QC
  indexing, read-only APIs, QC/artifact workbenches, safe downloads, and a real
  run-to-results browser gate.
- PR135–137 delivered adapter-owned authoring schemas, the schema/YAML/TSV input
  workbench, and immutable server-validated snapshots for run creation.
- PR138 is the Milestone C acceptance gate. Milestone C is complete only on a
  revision that contains PR138 and whose required browser gate proves authoring
  through exact result download with real API, SQLite, Redis/RQ, Snakemake, and
  cleanup. Revisions before that gate must treat the milestone as pending.
- PR139 is the next foundation: a reusable adapter conformance suite and a
  test-only minimal adapter. The intended later real-adapter order is RNA-seq,
  followed by Hi-TrAC/TracPre2; neither workflow is implemented by PR139.

Authentication, multi-user isolation, HPC/Kubernetes, object storage,
immutable workflow bundles, Agent write operations, and a second real adapter
remain outside the current delivery milestone.

## 1. Context after PR76

PR76 (`refactor(config): remove dead code from validator.py`) is the validator
decomposition freeze point. The validator decomposition is now considered
complete enough for the current architecture. Future work should not continue
validator cleanup as the main line unless a specific bug or compatibility issue
requires a narrowly scoped follow-up.

The current repository should be treated as the first ENCODE-style workflow
adapter candidate. It already contains the domain-specific config and sample
validation, DAG dry-run support, manifest logic, artifact catalog, QC helpers,
and Snakemake workflow implementation for ChIP-seq, CUT&Tag, ATAC-seq, and
MNase-seq. That makes it a strong adapter seed, not the mandatory home for all
future omics workflows.

The untracked local `encode-pipeline-architecture-report.md` report was useful
input material, but it is not the final plan and is not required by this
repository. Its operational guidance around service boundaries, run lifecycle,
artifacts, QC, provenance, and frontend concepts is valuable. Its assumption
that this repository evolves into a single pipeline fullstack app is superseded
by the platform plus adapters direction described here.

## 2. Adopt / Modify / Defer / Reject

### Adopt

- Build stable service boundaries before FastAPI routes, so route handlers stay
  thin and do not accumulate workflow logic.
- Keep the frontend at workflow and run level. Users should see workflows,
  runs, issues, artifacts, QC, and provenance rather than Snakemake internals.
- Treat submitted configs, sample sheets, workspace snapshots, artifacts, QC,
  logs, and provenance as first-class run concepts.
- Keep long-running workflow execution outside HTTP route handlers.

### Modify

- Replace the single-pipeline fullstack app direction with a workflow platform
  plus adapter architecture.
- Treat this repository as the first adapter candidate for the ENCODE-style
  ChIP-seq/CUT&Tag/ATAC-seq/MNase-seq workflow.
- Allow future RNA-seq, Hi-C, Hi-TrAC/cLoops2, and single-cell workflows to
  live in separate repositories or packages if they implement the same adapter
  contract.
- Put workflow-neutral concerns in the platform and workflow-specific concerns
  in adapters.

### Defer

- ARQ/Redis worker design.
- SQLite run tracker schema and migration strategy.
- SSE log streaming.
- FastAPI implementation.
- TypeScript frontend implementation.
- RO-Crate generation.
- Authentication and authorization.
- HPC scheduler and cluster-specific execution details.

These remain plausible later choices, but they should follow the Result/Issue
model, adapter contract, registry, and service boundary decisions.

### Reject

- Continuing further validator cleanup as a deprecated PR77/78/79-style path
  from the external report. That path is not required GitHub PR sequencing and
  may be treated only as an optional backlog if a concrete compatibility issue
  justifies it.
- Putting all future omics workflows into this repository by default.
- Letting platform services, API routes, frontend code, or agents bypass the
  adapter boundary and call workflow internals directly.
- Treating natural-language summaries as provenance.

## 3. Platform ownership

The platform owns workflow-neutral product and runtime concerns:

- Workflow registry: discovery, workflow IDs, versions, capabilities, and
  adapter lookup.
- Run lifecycle: creation, validation, queueing, running, cancellation,
  completion, failure, and status history.
- Workspace conventions: project and run directory layout, immutable submitted
  input snapshots, logs, results, QC, artifacts, and provenance locations.
- Result / Issue model: structured success/failure results and field-level
  issues with severity, code, user message, technical detail, and source.
- Artifact/QC/provenance model: workflow-neutral vocabulary for output files,
  QC summaries, provenance records, and their API/frontend representation.
- API and frontend UX: FastAPI route shape, TypeScript screens, form behavior,
  run monitor UX, artifact downloads, QC display, and provenance display.
- Agent orchestration guardrails: what agents may inspect, draft, explain,
  propose, and execute only with confirmation.

The platform must not encode ENCODE-specific config keys, sample-sheet columns,
rule names, artifact paths, or Snakemake command details except through an
adapter-provided response.

## 4. Adapter ownership

Each workflow adapter owns workflow-specific behavior:

- Workflow metadata and capabilities: name, version, assay/domain coverage,
  supported execution engines, schema versions, and feature flags.
- Config and sample schema: the accepted input shape, defaults, enums, field
  descriptions, and sample-sheet requirements.
- Validation: config validation, sample validation, cross-field checks,
  workflow-specific warnings, and translation into the platform Issue model.
- DAG preview support: dry-run or graph-preview behavior exposed as a
  workflow-neutral step graph.
- Workspace build: copying or rendering the submitted config, samples, workflow
  files, profiles, and run-local metadata into the platform workspace layout.
- Command construction: Snakemake, Nextflow, or other engine command lines,
  including profiles, cores, config files, and run directory assumptions.
- Artifact/QC/provenance extraction: mapping workflow outputs, QC reports, and
  provenance facts into platform-owned models.

For the current repository, the adapter wrapper should reuse existing modules
where possible: `encode_pipeline.config`, `encode_pipeline.samples`,
`encode_pipeline.cli.dag`, `encode_pipeline.manifest`,
`encode_pipeline.artifacts`, and `encode_pipeline.qc`. The wrapper should not
change workflow behavior while introducing the contract.

## 5. Agent boundaries

Agents are collaborators, not unchecked executors.

Agents may:

- Draft workflow configs and sample-sheet edits for review.
- Explain validation errors using structured issues from adapters.
- Diagnose logs and summarize likely failure causes.
- Summarize QC and artifact status from platform models.
- Propose rerun, patch, or review plans.

Agents must not:

- Bypass workflow adapters or call workflow internals directly when a platform
  contract exists.
- Directly run, kill, or delete production runs without explicit confirmation.
- Modify workflow code on `main` without review.
- Treat natural-language summaries as provenance.

Agent-facing tools should operate through platform services. When a tool could
change run state or files, it must expose a confirmation boundary and record the
human-approved action separately from any agent-generated explanation.

## 6. Development operating model

The development process follows a staged workflow:

1. Brainstorm the direction and constraints.
2. Write a reviewable spec.
3. Write an implementation plan.
4. Implement only the approved scope.
5. Verify before reporting completion.

Agency agents provide design and review lenses. In this roadmap those lenses
are:

- Software Architect: boundaries, dependency direction, long-term trade-offs.
- Multi-Agent Systems Architect: agent permissions, failure modes, confirmation
  gates, and provenance limits.
- Backend Architect: service boundaries, run lifecycle, adapter contract, and
  execution isolation.
- Technical Writer: clarity, readability, and roadmap sequencing.
- Code Reviewer: scope discipline, reviewability, and verification.

The GitHub plugin should provide PR, CI, review, and merge-state facts when PR
work is active. Coding agents perform scoped implementation against an approved
plan. They should not infer new scope from stale branches, stale reports, or
untracked research artifacts.

## 7. Revised phases

### Phase 0: validator decomposition freeze

Declare PR76 as the freeze point for validator decomposition. Do not continue
the old validator cleanup path as the main line.

### Phase 1: roadmap document

Land this docs-only roadmap so future PRs have a shared direction and do not
mix platform, adapter, service, API, frontend, and agent concerns prematurely.

### Phase 2: Result / Issue model

Define the platform-owned result and issue vocabulary. This should be small,
typed, and workflow-neutral. It should support validation errors, warnings,
technical details, user-facing messages, and source locations.

### Phase 3: WorkflowAdapter contract

Define the adapter interface and expected data flow. The contract should cover
metadata, schemas, validation, DAG preview, workspace build, command
construction, and artifact/QC/provenance extraction.

### Phase 4: current ENCODE adapter wrapper

Wrap existing ENCODE-style package behavior behind the adapter contract without
changing scientific behavior, CLI behavior, Snakemake rules, or output paths.

### Phase 5: workflow registry

Add a platform registry that can discover and select adapters by workflow ID,
version, and capability. This is the point where separate workflow packages
become practical.

### Phase 6: service layer

Build platform services around registry-selected adapters. Services own run
lifecycle, workspace conventions, issue aggregation, and artifact/QC/provenance
modeling. Services should exist before FastAPI routes.

### Phase 7: FastAPI backend

Add thin HTTP routes over the service layer. Routes should translate HTTP
requests and responses only; workflow logic remains in adapters and run logic
remains in services.

### Phase 8: TypeScript frontend

Build the workflow-neutral frontend around schemas, issues, run state,
artifacts, QC, and provenance. The frontend should not need to know whether the
workflow engine is Snakemake, Nextflow, or another engine.

### Phase 9: agent orchestration

Add agent-facing orchestration only after the platform models and service
boundaries exist. Agents should use platform services, respect confirmation
gates, and produce recommendations separate from provenance records.

This is a conceptual platform sequence, not a GitHub PR-number sequence.

## 8. Near-term sequence

1. Next main-line roadmap PR: land this docs-only roadmap.
2. Define the Result / Issue model.
3. Define the WorkflowAdapter contract.
4. Wrap the current ENCODE-style workflow as the first adapter candidate.
5. Add the workflow registry.
6. Add service layer boundaries.
7. Add FastAPI routes.
8. Add TypeScript frontend surfaces.
9. Add agent orchestration through platform services.

Further validator cleanup as a deprecated PR77/78/79-style path from the
external report is not required GitHub PR sequencing and remains only an
optional backlog path for specific compatibility issues.

## 9. Non-goals

This roadmap PR does not include:

- Production code changes.
- Service layer implementation.
- FastAPI or frontend implementation.
- Changes under `docs/superpowers/`.
- Inclusion of `research/` material.
- CI, lockfile, Snakefile, or workflow rule changes.
- Attribution trailers.
- Result / Issue model implementation.
- WorkflowAdapter implementation.
- Registry, run lifecycle, execution worker, or agent-tool implementation.
