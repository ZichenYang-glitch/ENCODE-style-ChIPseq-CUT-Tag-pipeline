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

The bundled registry contains two scientific adapters:

- ENCODE-style ChIP-seq, CUT&Tag, ATAC-seq, and MNase-seq through Snakemake;
  and
- `bulk-rnaseq` through the pinned nf-core/rnaseq 3.26.0 Nextflow runtime.

Workflow-specific schemas, assay policy, commands, output paths, artifact
extraction, and scientific behavior remain adapter or workflow
responsibilities. Bulk RNA-seq authoring and validation remain available
without execution assets; execution is declared available only after the
complete operator binding passes live admission.

The platform is not a hosted multi-tenant service. SQLite and the local
filesystem are deliberate canonical stores; Redis/RQ provides an execution
handoff rather than independent lifecycle truth.

## Delivered platform baseline

The current product supports:

- discovery of registered workflows and adapter-owned capabilities;
- the ENCODE-style epigenomics and Bulk RNA-seq adapters in the default
  registry, without workflow-specific API, lifecycle, or frontend stores;
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
- a pinned, read-only Omics Intake Bundle 0.2 inspection boundary for the
  ENCODE adapter, without producer imports or execution authorization;
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

## Completed: intake consumption and Bulk RNA-seq delivery

The initial Omics Intake consumption boundary is delivered as a service-only,
read-only inspection of the pinned Bundle 0.2 public contract. HelixWeave
verifies contract identity, safely observes required local files, delegates
mapping to the ENCODE adapter, and revalidates the mapped inputs. It does not
import Omics Intake code, mutate the Bundle, create a validated snapshot, or
authorize execution. Durable Bundle provenance and a product import flow remain
separate future decisions.

Bulk RNA-seq is delivered as the second bundled adapter. Its product contract
includes adapter-owned authoring schemas, validation without runtime assets,
path-free availability, backend create/start admission, and the existing
workflow-neutral run, artifact, QC, and download surfaces. Execution is pinned
to nf-core/rnaseq 3.26.0 and was accepted with controlled synthetic
STAR+Salmon and SortMeRNA runs through SQLite, Redis/RQ, Nextflow, and offline
containers. That evidence proves execution and product contracts, not
biological validity or production-scale performance.

## Current delivery priorities

The order below expresses the current decision and delivery sequence. Each
priority still requires scoped review; none carries an implied release date.

### 1. v0.3.0 local trial and release convergence

Consolidate the two-workflow local product path and its supported installation,
upgrade, diagnostic, trial, and release contracts.

Expected work:

- keep package, API, frontend, documentation, and release identities aligned;
- prove a fresh distribution install and supported database upgrade;
- make both workflow availability states and optional runtime prerequisites
  clear without exposing private deployment coordinates;
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

### 2. Deferred intake and deployment decisions

Any durable Omics Intake provenance, snapshot/run import, authentication,
remote execution, or additional adapter work requires a separate reviewed
decision. It must not be inferred from the delivered read-only Bundle
inspection or the two bundled local adapters.

## Agent direction

The current Agent surface is already advisory and read-only. It may explain
schemas and validation issues through platform services, but it cannot submit,
start, cancel, mutate, or delete runs. Further Agent expansion requires a
separate decision and must remain workflow-neutral rather than encode either
adapter's private vocabulary.

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
