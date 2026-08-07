# Artifact Publication and Provenance Decision

Status: Accepted for the first read-only publication surface

## Context

Successful runs already keep their artifact bytes in the run workspace and
replace the current artifact index as one durable generation. Administrators
need a cross-run view of those successful publications without introducing a
second artifact store, inferring provenance from adapter metadata, or changing
the revision-bound download boundary.

The earlier data-foundation design described a future storage-pool move and
per-artifact contributor model. That model is not the implementation selected
for this milestone. This decision supersedes those publication-specific
sections while leaving the Project, Sample, immutable binding, and input
registry decisions intact.

## Decision

SQLite remains the only publication authority. The platform adds one
append-only `artifact_publications` table. A row is identified by the existing
run, artifact, and artifact-generation identity and records only the existing
artifact revision, its output type, and the UTC publication time. Project,
workflow, and Sample identities are not copied into the row.

A publication batch is inserted by the existing `replace_artifacts`
repository operation only when a complete artifact manifest installs a new
durable generation. The artifact rows, result state, `artifacts_indexed` event,
and publication rows commit in the same SQLite transaction. The in-memory
repository performs the equivalent state change under one lock. Failed,
cancelled, partial, stale, duplicate, and rolled-back attempts cannot leave
publication evidence. Database protections reject updates and deletes; they
do not create publications or act as an alternate writer.

The migration performs no historical backfill. Existing artifacts have no
trustworthy publication instant, so the absence of a publication row remains
meaningful rather than being replaced with fabricated evidence.

Artifact bytes are not copied, moved, or renamed. The existing run artifact
storage and `run://` identity remain unchanged. Downloads continue through the
existing run-scoped endpoint and must present the exact artifact generation
and revision. A publication for an older generation is historical metadata;
it does not create a second file resolver or make superseded bytes current.

## Provenance projection

The public read model derives:

- Project identity from the run's immutable Project binding;
- workflow identity from the immutable run record;
- associated input Sample revisions from the authoritative ordered
  `run_samples` binding; and
- current or superseded status by comparing the publication generation with
  the run's current artifact generation.

Associated run samples describe the run input scope. They are not
per-artifact Sample attribution. The platform never derives this association
from filenames, relative paths, artifact metadata, QC labels, or
adapter-private validation payloads.

The cross-run list supports exact Project, run, workflow, output-type, and
associated-run-Sample filters plus a half-open UTC publication interval. Its
descending keyset cursor binds the normalized filters and complete ordering
position. Details use the existing run, artifact, and generation coordinates;
there is no second public artifact identity.

Public responses expose only opaque identities and the fields above. They do
not expose physical paths, artifact URIs, filenames, open metadata, operator
configuration, environment values, display-name payloads, or exception text.

## Consequences and non-goals

This model provides durable discovery and provenance evidence for new
successful generations while preserving the existing execution and file
boundaries. It deliberately does not add retention, recovery management,
object storage, access-control policy, byte replication, storage pools,
per-artifact contributors, or historical reconstruction. Any later feature
that needs those capabilities requires a separate schema and trust-boundary
decision.
