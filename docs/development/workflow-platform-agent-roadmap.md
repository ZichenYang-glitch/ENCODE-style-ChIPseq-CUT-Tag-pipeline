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

## Completed baseline through PR #168

The v0.3.0 local release, maintenance, and quality baseline is complete. CI has
separate PR-fast, full-main, platform-real, scientific-real, container,
frontend, browser, lint, lock, and coverage responsibilities. Current test
inventory and enforced floors live in the [quality baseline](coverage-policy.md);
tier ownership and timing live in the [development harness](harness.md).

The pinned Omics Intake Bundle 0.2 boundary remains service-only and read-only.
Bulk RNA-seq is delivered as the second bundled adapter, pinned to
nf-core/rnaseq 3.26.0 and integrated through the workflow-neutral lifecycle,
artifact, QC, API, and browser surfaces. Its protected evidence establishes
HelixWeave's adapter and offline product integration, not the biological
validity or internal scientific correctness of the unchanged upstream
workflow.

The data foundation now provides Project, Sample, SampleRevision, StoragePool,
InputFile, InputFileRevision, and input-use binding records. Runtime and adapter
admission, exact-checkout and package provenance, migration validation, timeout
handling, and Bulk qualification identity are fail-closed. PR #168 separated
Bulk's scientific execution identity from unrelated migration history while
retaining an independently controlled migration-execution inventory and an
explicit Bulk persistence contract.

## Tentative delivery sequence

The sequence below is a planning boundary, not an implemented feature list or
a release-date promise. Each PR must still earn scoped design, tests, and
review, and may be reordered by an explicit product decision.

### PR #169: governance and roadmap synchronization

Align durable repository guidance with this delivery sequence. This is a
documentation-only boundary and does not change product or execution behavior.

### PR #170: Reference Profiles

- Let administrators register, verify, enable, and disable
  operator-prepared references.
- Support multiple references such as GRCh38 and mm10 and select one per run.
- Keep reference acquisition and deployment outside the platform; ordinary
  users and the frontend cannot deploy or mutate references.
- Pin an immutable reference identity into each validated snapshot and run.

### PR #171: artifact publication and provenance

- Close the relationships among Project, Sample, Run, Artifact, time, and
  storage location.
- Keep artifact bytes in their existing storage location; record associations
  and expose an opaque public identity without duplicating the artifact.

### PR #172: run recovery and administrator operations

Add explicit administrator recovery and reconciliation operations around the
already durable restart-recovery and lifecycle contracts.

- Diagnose active ownership, queue identity, callback, result-indexing, and
  cleanup gaps without treating Redis as lifecycle truth.
- Allow exact-identity administrator failure only after ownership loss and
  required cleanup are proven.
- Requeue only the same dispatched, never-claimed assignment through a
  one-time durable request; never reset claim markers or cycle run status.
- Keep mutation on the local administrator CLI until the authentication and
  role boundary in PR #174 exists.

### PR #173: deployment and runtime management CLI

- Provide `install`, `doctor`, `verify`, `upgrade`, and `rollback` operations.
- Allow the platform and scientific runtime to be upgraded independently.
- Deploy a precompiled frontend rather than requiring a frontend toolchain on
  the target host.

### PR #174: LAN authentication and roles

Add the small trusted-team authentication boundary with administrator and
ordinary-user roles; this is not authorization for hosted multi-tenancy or a
general RBAC framework.

### PR #175: email notifications and dynamic QC summaries

- Notify on terminal runs; users may opt out while administrators receive
  notifications by default.
- Include only QC metrics that actually exist. A missing optional metric such
  as FRiP must not make notification delivery fail.

### PR #176: usability and Rosemary theme alignment

Improve workflow and run usability while applying the Rosemary visual theme
consistently across the existing product surface.

### PR #177: dual-workflow acceptance and release closure

Complete final acceptance for both bundled workflows and close the release
evidence without widening either workflow's scientific contract.

## Optional follow-on work

An IGV.js artifact browser, further Agent capabilities, PostgreSQL, Slurm or
cloud execution, Singularity/Apptainer, and additional maintenance trimming are
deliberately optional. They do not block the laboratory single-host delivery
sequence above. One recorded trimming candidate is separating the scientific
execution closure from the shared platform product closure in the bulk RNA-seq
execution manifest, so that authentication, API, persistence, and lifecycle
changes stop requalifying unchanged scientific implementation files; this is
deferred and must not fold qualification redesign into unrelated PRs.

## Agent direction

The current Agent surface is advisory, read-only, and deterministic-mock-backed.
It may explain schemas and validation issues through platform services, but it
cannot submit, start, cancel, mutate, or delete runs. Further Agent expansion
is optional, requires a separate decision, and must remain workflow-neutral
rather than encode either adapter's private vocabulary.

## Explicit non-goals

The delivery sequence above does not currently authorize:

- hosted multi-tenant isolation or complex RBAC beyond the planned LAN roles;
- Kubernetes, HPC scheduler integration, or microservice decomposition;
- PostgreSQL, object storage, or remote workspace semantics;
- Agent write actions or automatic workflow submission;
- arbitrary workflow loading without an approved adapter/deployment contract;
- further changes to the `encode_pipeline` import namespace, existing
  compatibility CLI names, repository slug, workflow identity, or artifact URI
  scheme. The approved v0.3.0 identity migration changes only the distribution
  to `helixweave` and adds the primary `helixweave` CLI;
- a frontend rewrite, server-side rendering, or a second frontend repository;
- scientific changes hidden inside platform or documentation work.

## Roadmap discipline

New work should advance one product outcome and name its exit evidence. Durable
architecture, persistence, public-contract, worker, or cross-repository changes
may require a focused decision record. Detailed completed checklists and commit
identifiers stay in Git history; the numbered sequence above is tentative
planning, not implementation evidence. Current test counts and coverage floors
stay in the [quality baseline](coverage-policy.md) so README, roadmap, and
operational docs do not drift independently.
