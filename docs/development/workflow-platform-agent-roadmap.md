# HelixWeave Product Roadmap

This roadmap records the maintained product boundary, delivered baseline, and
explicit follow-ons. It does not promise dates. Detailed implementation history
belongs in Git history, the changelog, and release evidence.

## Product boundary

HelixWeave is a workflow-neutral omics platform for one laboratory operating
one Linux x86_64 or systemd-enabled WSL2 host on a trusted LAN. SQLite and the
local filesystem are canonical stores; Redis/RQ is an execution handoff, not a
second lifecycle authority.

The bundled registry contains:

- ENCODE-style ChIP-seq, CUT&Tag, ATAC-seq, and MNase-seq through Snakemake;
  and
- `bulk-rnaseq` through the pinned nf-core/rnaseq 3.26.0 Nextflow runtime.

Workflow schemas, commands, scientific behavior, output paths, and artifact/QC
extraction remain adapter-owned. Authoring may be available without runtime
assets; execution remains fail-closed until the complete operator binding is
admitted.

HelixWeave is not a hosted multi-tenant or high-availability service. Trusted-
LAN authentication is an administrator/member boundary, not general RBAC or
multi-tenant isolation.

## Delivered baseline through PR #179

The current baseline includes:

- workflow-neutral schema authoring, structured validation, immutable
  snapshots, durable lifecycle/events/logs, cancellation, restart recovery,
  artifact/QC extraction, safe downloads, and responsive browser journeys;
- Project, Sample, SampleRevision, StoragePool, InputFile, InputFileRevision,
  and immutable input-use provenance;
- administrator-managed revisioned Reference Profiles bound to snapshots and
  runs;
- cross-run artifact publication and provenance;
- administrator diagnosis, exact-identity fail, and one-use requeue recovery;
- an offline Linux/systemd deployment operator with independently managed
  platform and scientific runtime slots;
- trusted-LAN accounts, sessions, CSRF protection, administrator/member roles,
  audited privileged actions, and terminal email preferences;
- staged FastQC/MultiQC evidence for the ENCODE-style workflow and dynamic
  current-generation QC notification summaries;
- the compact Rosemary desktop/mobile browser experience;
- pinned Node 22.23.1/npm 10.9.8 authority and independent CI, Lint, Lock
  Check, coverage, frontend, browser, platform-real, scientific-real,
  container, and Protected Bulk evidence tiers; and
- a read-only Omics Intake Bundle 0.2 inspection boundary and a read-only Agent
  explanation boundary.

Controlled tiny and synthetic evidence demonstrates integration, execution,
lifecycle, offline reuse, and result contracts. It is not biological
validation, public-data reproducibility, or production-scale performance.

## PR #180 release closure

PR #180 closes the v0.4.0 release candidate without adding scientific or UI
features. It synchronizes the release identity, canonical generated assets and
contracts, released-schema migration evidence, public release documentation,
and the three-asset publication boundary. It also provides the narrow public
offline producer that composes the existing platform, ENCODE runtime, and Bulk
RNA-seq deployment bundle formats from release artifacts and explicit
operator-owned caches.

The candidate is complete only after ordinary exact-HEAD checks, wheel/sdist
and two no-checkout installs, release-shaped bundle production, the released
schema-08 to schema-15 migration proof, one separately authorized supported-
host deployment window, and—if the final implementation manifest requires
it—one exact-HEAD combined Protected dispatch. Date freeze, merge, tag, and
publication remain separate human gates.

## Explicit follow-ons and non-goals

The following do not block v0.4.0 and require separately scoped decisions:

- further Agent capabilities or any Agent write action;
- PostgreSQL, object storage, Kubernetes, multi-host or high-availability
  operation, Slurm/HPC/cloud executors, or remote workspace semantics;
- Singularity/Apptainer support or published Docker/OCI images;
- IGV.js or other new scientific/browser capabilities;
- dependency-major maintenance, CSP/proxy hardening beyond the trusted-LAN
  contract, and additional UI polish; and
- redesigning Bulk qualification, CI, runner, release automation, or the
  existing deployment bundle/state contracts.

Do not change the `encode_pipeline` import namespace, compatibility CLI names,
repository slug, workflow identities, or artifact URI scheme without an
explicit compatibility decision.

## Roadmap discipline

New work should advance one product outcome and name its exit evidence. Keep
test counts and coverage floors in the quality baseline, implementation history
in Git, and release evidence in `docs/release-checks/` so this roadmap remains a
short statement of delivered scope and future boundaries.
