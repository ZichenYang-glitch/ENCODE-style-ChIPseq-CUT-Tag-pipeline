# HelixWeave Product Roadmap

This is the maintained product roadmap for HelixWeave. It records the current
boundary, delivered platform baseline, and ordered product priorities. It does
not promise release dates. Implementation history belongs in Git history and
release evidence rather than this document.

## Current product boundary

HelixWeave is a workflow-neutral platform for reproducible omics analysis on a
local workstation or within a small trusted team. It owns input authoring,
validation, durable run state, local execution coordination, artifacts, QC, and
the browser experience around those concepts.

The bundled ENCODE-style ChIP-seq, CUT&Tag, ATAC-seq, and MNase-seq Snakemake
workflow is the first and currently only scientific adapter. Workflow-specific
schemas, assay policy, commands, output paths, artifact extraction, and
scientific behavior remain adapter or workflow responsibilities.

The platform is not a hosted multi-tenant service. SQLite and the local
filesystem are deliberate canonical stores; Redis/RQ provides an execution
handoff rather than independent lifecycle truth.

## Delivered platform baseline

The current product supports:

- discovery of registered workflows and adapter-owned capabilities;
- schema-driven config, sample, and option authoring with advanced text modes;
- structured validation issues without exposing adapter-private payloads;
- immutable, server-validated input snapshots for submitted runs;
- durable lifecycle, event, log, cancellation, and restart-recovery state;
- Redis/RQ execution of planned Snakemake commands outside HTTP handlers;
- truthful process-group cancellation and terminal-state acknowledgement;
- filterable run history with stable deep links;
- adapter-owned artifact and QC extraction after successful execution;
- safe artifact listing and download without arbitrary path access;
- workflow, run, activity, artifact, and QC views in the browser;
- a deterministic local input-to-results demonstration path;
- a read-only Agent boundary for schema and issue explanation; and
- reusable adapter conformance tests with a minimal test adapter.

Detailed ownership and safety rules are maintained in the
[architecture overview](../architecture/platform-overview.md). Scientific
contracts remain in the assay, configuration, sample, output, QC, and
reproducibility references under `docs/`.

## Completed: maintenance and quality baseline

The maintenance baseline is complete. Historical process plans and
stage-numbered test scaffolding were retired or migrated into maintained
behavior contracts. CI now has distinct PR-fast, full-main, platform-real,
scientific-real, container, frontend, browser, lint, lock, and coverage
responsibilities without a second deterministic pytest producer.

The resulting gate preserved scientific workflow outputs, public API routes,
persistence identity, CLI names, and visible runtime behavior. Current test
inventory, measured coverage, and enforced floors have one authoritative home
in the [quality baseline](coverage-policy.md); tier ownership and timing live in
the [development harness](harness.md).

## Current delivery priorities

The order below expresses the current decision and delivery sequence. Each
priority still requires scoped review; none carries an implied release date.

### 1. Omics Intake Bundle → HelixWeave consumption boundary

Define how HelixWeave consumes reviewed inputs from an Omics Intake Bundle
without making the platform own upstream acquisition or silently trust
external paths and metadata.

Expected work:

- separate producer and consumer responsibilities;
- define a versioned, fail-closed handoff with explicit identity and provenance;
- map accepted bundle data into adapter-owned authoring and validation inputs;
- preserve immutable snapshot, redaction, and workspace-containment rules; and
- prove rejection behavior for incomplete, ambiguous, or unsafe bundles.

Exit evidence:

- a focused boundary decision is reviewed before runtime integration;
- deterministic fixtures prove both accepted and rejected handoffs; and
- workflow-neutral services do not acquire adapter-private or intake-private
  field dependencies.

### 2. Second-adapter research and implementation

Research a bounded second adapter against the existing registry, authoring,
execution, artifact, and QC contracts. Bulk RNA-seq is the current candidate,
not an approved product commitment or delivery date.

Expected work:

- evaluate candidate ownership, maintained engines, schemas, validation,
  workspace planning, commands, artifacts, and QC semantics;
- record an explicit adapter selection decision before implementation;
- keep generic API and frontend code free of workflow-specific branches; and
- extend deployment/source composition deliberately rather than claiming
  arbitrary zero-configuration plugin loading.

Exit evidence if a candidate is approved:

- the adapter passes the reusable conformance suite;
- the existing ENCODE adapter remains unchanged and green;
- one deterministic small-data run produces indexed artifacts and QC; and
- operators can compare adapter capabilities before creating a run.

### 3. Product experience and deployment convergence

Consolidate the current local product path after its input, execution, and
result surfaces are in place.

Expected work:

- reduce setup and navigation friction from workflow choice to evidence;
- make local prerequisites, diagnostics, storage, recovery, and cleanup clear;
- tighten empty, loading, failure, and long-running states across desktop and
  mobile views; and
- align packaging and operator documentation around the supported workstation
  and small trusted-team deployment.

Exit evidence:

- a fresh locked install can complete the deterministic product journey;
- operator-owned state and process cleanup are explicit and testable; and
- usability changes preserve accessibility, public contracts, and adapter
  boundaries.

## Agent direction

The current Agent surface is already advisory and read-only. It may explain
schemas and validation issues through platform services, but it cannot submit,
start, cancel, mutate, or delete runs. Further Agent expansion should follow
the consumption-boundary and multi-adapter evidence rather than encode the
first workflow's private vocabulary.

## Explicit non-goals

The maintained roadmap does not currently authorize:

- authentication, authorization, multi-tenant isolation, or complex RBAC;
- Kubernetes, HPC scheduler integration, or microservice decomposition;
- PostgreSQL, object storage, or remote workspace semantics;
- Agent write actions or automatic workflow submission;
- arbitrary workflow loading without an approved adapter/deployment contract;
- changes to the Python distribution, import namespace, CLI names, repository
  slug, workflow identity, or artifact URI scheme;
- a frontend rewrite, server-side rendering, or a second frontend repository;
- scientific changes hidden inside platform or documentation work.

## Roadmap discipline

New work should advance one product outcome and name its exit evidence. Durable
architecture, persistence, public-contract, worker, or cross-repository changes
may require a focused decision record. Completed checklists, commit identifiers,
and branch sequencing stay in Git history. Current test counts and coverage
floors stay in the [quality baseline](coverage-policy.md) so README, roadmap,
and operational docs do not drift independently.
