# Data Foundation Architecture and Compatibility Decision

Status: Accepted for the data-foundation milestone

## Context

HelixWeave currently records immutable validation snapshots, run lifecycle,
artifact generations, and QC metadata, while workflow inputs and result bytes
remain coupled to caller-supplied paths and run workspaces. That is sufficient
for the v0.3 local workflow platform, but it cannot make durable claims about
which project, logical sample revision, registered input revision, or published
artifact bytes produced a result.

This milestone deliberately advances a persistence and product decision that
the current roadmap had deferred. The explicit delivery request is the decision
gate. It does not widen the supported deployment beyond one workstation or a
small trusted team, add hosted storage, or add authentication and role-based
access control.

Three independent reviews of the existing persistence, snapshot/run, adapter,
storage, and artifact contracts found four constraints that shape the design:

1. Existing adapter sample rows cannot be reconstructed as platform samples.
   The bundled adapters have different row cardinalities and keep successful
   validation payloads private.
2. Exact sample and input revisions must be frozen with validation. Supplying
   them later at run creation would permit provenance drift.
3. A durable artifact publication needs a recoverable boundary between
   filesystem state and SQLite state.
4. Existing paths, QC sample labels, artifact metadata, and filenames do not
   prove catalog identity and must never be reinterpreted as provenance.

## Decision

### Authority and representations

SQLite is authoritative for identities, relationships, revisions, lifecycle
state, publication state, and queryable provenance.

Canonical JSON manifests are audit and recovery evidence for immutable bindings
and published bytes. JSON manifests do not become a competing lifecycle
database: a manifest without the matching valid SQLite state is not publicly
visible, and a SQLite publication whose manifest or bytes cannot be verified
fails closed.

TSV files and reviewed intake Bundles are import formats only. They can propose
catalog records, but they are not identity authorities and are never consulted
as mutable runtime state after a revision is committed.

All objects use opaque, stable internal IDs. A Project display name is mutable
presentation data, not a directory identity, filesystem coordinate, or
uniqueness key. No workspace or durable-storage directory is derived from a
Project or Sample display name.

### Workflow-neutral entities

- **Project** is the lifecycle, ownership, and query-scope root for catalog
  samples, runs, approved storage, inputs, and durable artifacts. It has a
  stable ID, display name, kind, `created_at`, and optional `archived_at`.
- **Sample** is a stable logical identity owned by exactly one Project. It has a
  stable ID and a project-scoped stable key. A Sample is not an adapter
  sample-sheet row, lane, library, control role, or filename.
- **SampleRevision** is an immutable, append-only revision of one Sample. It
  stores an opaque adapter- or operator-owned canonical payload and its digest.
  One revision may correspond to zero, one, or several adapter rows. Platform
  code does not interpret adapter-private scientific roles.
- **StoragePool** is an approved logical storage destination or ingress scope.
  SQLite stores its stable ID and an opaque operator configuration key, never
  its physical absolute root.
- **InputFile** is a stable logical regular-file identity owned by one Project
  and located in an approved StoragePool. Each immutable input-file revision
  stores a safe pool-relative path, byte size, full content checksum, and
  revision time. Regular-file InputFile support is additive; it does not imply
  that every execution input or adapter is managed. Directory trees, index
  prefixes, multi-file scientific closures, and input Bundles are outside the
  first InputFile implementation. An adapter must not be described as fully
  managed until every execution-required input use has a qualified catalog
  representation and resolver.
- **RunSample** is the ordered many-to-many association from a Run to the exact
  SampleRevision records frozen in the consumed snapshot. A run may bind
  several samples and the same sample revision may contribute to several runs.
- **ArtifactSample** is the many-to-many association from a published Artifact
  revision to zero, one, or several contributing SampleRevision records.
  Several contributors support pooled artifacts. An empty association is
  permitted only with an explicit unresolved provenance status.

All catalog relationships preserve evidence. User-facing lifecycle operations
archive records rather than deleting them. Foreign keys restrict removal of
Project, Sample, revision, input, and publication evidence; only run-owned
association rows may cascade with a run if the existing lifecycle permits that
run to be removed.

### UTC timestamp semantics

Every persisted instant is timezone-aware at the domain boundary, normalized
to UTC on write, and restored as UTC after SQLite reads. Naive datetimes are
invalid. The timestamp vocabulary is:

- `created_at`: the identity or immutable revision was committed;
- `started_at`: execution or a recoverable publication attempt began;
- `ended_at`: that execution or attempt reached a terminal state;
- `produced_at`: the workflow output bytes were successfully indexed as a
  logical artifact generation;
- `published_at`: the final SQLite publication compare-and-set succeeded after
  its durable manifest and bytes were committed and verified;
- `archived_at`: an operator made the object unavailable for new bindings
  without deleting historical evidence.

These timestamps describe different events and must not be copied into one
another merely because an earlier schema lacks evidence.

### Snapshot-time project, sample, and input binding

The existing normalized workflow-input snapshot remains an immutable v1
compatibility anchor. It is not mutated or reinterpreted. A separate,
versioned **snapshot binding envelope** records:

- the Project ID and binding mode;
- ordered exact SampleRevision IDs and their immutable payload digests;
- the workflow ID and adapter contract version;
- every adapter-declared input use as a stable opaque key, its occurrence, its
  exact input-use capability and closure-contract versions, and its provenance
  mode;
- for every input use in `managed_revision_v1` mode, the ordered zero-or-more
  opaque resource/revision IDs, member InputFile revision IDs, logical member
  coordinates, and generic closure digest;
- explicit `transitional_unmanaged_v1` records for execution-required input
  uses that have not reached managed cutover, without copying or interpreting
  their physical paths;
- the existing workflow-input snapshot digest; and
- a canonical digest over the complete envelope and its scheme version.

The validation transaction resolves and verifies every referenced revision,
checks that it belongs to the same active Project, and freezes the envelope.
Run creation accepts only a snapshot ID and copies the envelope and ordered
associations atomically while consuming the snapshot. It never performs a
“latest revision” lookup and never accepts replacement revision IDs.

Every run-creation entry point, including SQL and in-memory repositories and
the direct compatibility path, creates the Run, Project binding, ordered
RunSample/RunInput associations, first lifecycle event, and snapshot-consumption
record as one atomic repository operation. The SQLite implementation performs
that operation in one transaction. Idempotent replay compares the envelope
digest, every ordered SampleRevision association, and the complete ordered
per-input-use mode, capability version, resource membership, and closure
evidence. Missing or unequal binding evidence fails closed.

Archiving a Project, Sample, or input prevents a new snapshot from binding it.
It does not invalidate an already frozen, unexpired snapshot or prevent its
one-time consumption. That rule preserves the exact reviewed input while
stopping new work from selecting an archived object.

For a non-Legacy Project, a new trusted-provenance snapshot must bind at least
one exact SampleRevision. Managed resource revisions become mandatory only for
an input use whose exact workflow ID, adapter contract version, and input-use
capability have completed cutover. Adapter-owned workflow semantics stay in
the adapter: an adapter declares required input uses, permitted provenance-mode
combinations, scientific closure validation, and artifact contributors, while
the platform validates opaque identities, authorization, digests, and
membership without inspecting sample columns, workflow rule names, or
scientific roles.

A mixed managed/unmanaged snapshot is permitted only as an explicit migration
state. The exact adapter contract must allow that combination, every
execution-required input use must appear exactly once, and each use records
either `managed_revision_v1` or `transitional_unmanaged_v1`. Missing,
duplicated, undeclared, or contract-forbidden combinations fail closed. A
managed use can never fall back to a physical path because a different use in
the same adapter remains transitional.

A concrete snapshot is fully managed only when every input use required by
that execution is in `managed_revision_v1` mode. An optional use that is not
activated does not downgrade the snapshot. An adapter contract may be
described as fully managed only after every execution-required input use it can
activate for its supported routes has a qualified catalog representation,
adapter validator, and worker resolver.

### Legacy compatibility, capability-scoped transition, and conservative migration

The migration creates one deterministic reserved Legacy Project with a system
kind and stable ID. The reserved identity and kind cannot be renamed, archived,
or repurposed by the administrator CLI, and the Project cannot receive trusted
SampleRevision or InputFile bindings.

Every pre-migration run and every pre-migration validation snapshot receives a
deterministic `legacy_v1` Project-binding row so that Project ownership is
total and queryable. The envelope records unresolved sample/input provenance
and the original workflow-input digest. A consumed snapshot and its run must
have identical Legacy binding evidence or the upgrade fails. An unconsumed,
unexpired snapshot later copies that binding through the normal atomic consume
operation.

The migration must not infer sample provenance from filenames, directories,
workflow-local QC labels, artifact metadata, saved configuration, private
validation payloads, or sample-sheet rows. It creates no synthetic Sample,
SampleRevision, RunSample, InputFile, input binding, or ArtifactSample records
for historical runs. Their missing relationships are marked `unresolved`, not
guessed.

Persisted pre-migration and historical v1 snapshots and runs remain readable
and, when otherwise eligible, consumable as `legacy_v1`. Historical payloads
may already contain physical paths; migration does not rewrite them, public
projections redact them, and only the private compatibility execution resolver
may consume them.

Managed cutover is per-adapter, per-input-use. Its identity is the exact
workflow ID, adapter contract version, and input-use capability recorded in the
snapshot envelope. An input use enters managed-only mode only after its
complete closure can be represented by catalog revisions, the adapter can
validate the scientific layout, and the worker can independently resolve and
reverify the same closure. Only that qualified input use rejects a supplied
physical path with a stable result. Cutover of one use never changes the path
policy of another use, even within the same adapter.

An input use that has not completed qualification remains explicitly
`transitional_unmanaged_v1` when the exact adapter contract permits it. The
existing trusted-local v1 path validation remains available for that use and
does not become managed provenance. A new catalog-aware snapshot may therefore
bind trusted Project and SampleRevision identities while recording transitional
input provenance per use; it is not assigned to the Legacy Project merely
because some inputs remain transitional. A fully unbound compatibility request
may still use the reserved Legacy Project with unresolved provenance.

`legacy_v1` is a provenance classification; it neither grants nor prohibits
physical-path submission. For a new catalog-aware snapshot, path acceptance is
decided only by the frozen adapter/input-use capability and provenance mode. A
historical v1 snapshot has no synthesized per-use capability record; its
eligibility follows its frozen v1 workflow contract and digest, and only the
private compatibility resolver may consume an eligible path. Project/Sample
stage does not reject physical paths that the current adapter contract already
validates, and it does not infer managed input identity from them. A separate
local administrator compatibility command may create an unmanaged Legacy
snapshot, but it is never exposed as trusted provenance.

New catalog-aware clients use the bound contract. Existing run IDs, inputs,
workspace roots, artifact URIs, generations, revision-bound downloads,
lifecycle fields, and scientific behavior keep their previous meaning.

This exception is intentional: “every run has a Project” is true after
migration, while “every run has provable SampleRevision contributors” is only
true for bound non-Legacy runs. Compatibility must not manufacture evidence.

### Administrator boundary

Catalog and local-storage mutations are exposed through an
administrator-only local CLI. The CLI operates against an explicitly selected
local database and private operator configuration. Project create/archive,
Sample import/revise, StoragePool registration, Project/pool binding, and local
InputFile registration are operator actions.

This milestone must not add an unauthenticated management write API. Ordinary
HTTP and frontend clients may select already-authorized opaque IDs through
workflow contracts, but they cannot create or archive catalog objects,
register pools, bind physical roots, or import server files. They cannot submit
a physical path for an input use that has completed managed cutover. Until a
specific input use is qualified, its exact existing trusted-local adapter
contract may continue to accept a validated path; this compatibility is not a
generic path-management API. Read-only public projections expose only
allowlisted logical metadata.

The deployment remains locally trusted; the CLI boundary is an operational
boundary, not a claim of OS-level multi-tenant authorization.

### Storage pools and input revisions

An operator-private configuration maps each StoragePool configuration key to
one absolute local root. The same mapping is loaded by API and worker
composition. For newly registered inputs, bound snapshots, and publications,
roots never enter SQLite rows, public manifests, artifact metadata, logs,
exception text, OpenAPI schemas, or generated clients. Historical v1 snapshot
payloads and an explicitly transitional adapter-private v1 payload are the
compatibility exceptions; they are not rewritten, are never copied into the
binding envelope, and remain redacted from public projections.

A Project binds to an approved, active StoragePool by stable ID. Local
InputFile registration is no-follow and descriptor-owned: every path component
is resolved beneath the configured root, traversal and symlinks are rejected,
the leaf must be a regular file, and size and full SHA-256 are computed from
the held file descriptor with before/after mutation checks. Registration does
not download FASTQ, invoke GenomePy, interpret scientific formats, or move
real user data. A changed byte stream creates a new immutable revision rather
than overwriting history.

For a managed input use, the platform owns opaque resource and revision IDs,
Project and StoragePool authorization, ordered leaf membership, leaf size and
checksum verification, the generic closure digest, and immutable binding. A
workflow-neutral resolver reopens member revisions beneath configured roots and
produces a short-lived private descriptor or layout-preserving view. Validation
and worker materialization independently resolve and reverify the same ordered
closure. Resolved absolute paths are never persisted in the binding envelope,
logged, returned, or exposed through public contracts.

The adapter owns the scientific meaning and qualification of that view,
including a Bowtie2 prefix and its complete alternative file set, FASTA
sidecars, STAR, Salmon, and SortMeRNA directory or manifest layouts, and
producer, tool, container, and build parameters. Platform code does not know
those filenames, alternative groups, tool semantics, or workflow roles. It
must not infer closure identity from a filename, directory name, sample sheet,
configuration path, or adapter-private payload.

A request for `managed_revision_v1` fails admission with a stable
unsupported- or invalid-capability result when the exact adapter/input-use
capability, complete catalog closure, adapter validator, or worker resolver is
missing. Once an input use is qualified and cut over, it never falls back to a
caller path. An unqualified input use instead follows only an explicitly
permitted `transitional_unmanaged_v1` resolver and the existing private adapter
validation contract; regular-file InputFile support is additive and does not
silently convert a prefix, sidecar set, manifest-referenced resource, or
directory tree into managed provenance. Existing historical workspaces remain
readable through the compatibility resolver, and assigning a Project pool does
not silently relocate or reinterpret an existing run.

### Artifact lifecycle and one-copy publication

The run workspace is scratch/working storage. A workflow output becomes a
logical artifact at `produced_at`; it becomes a durable Artifact only after a
publication is verified and indexed at `published_at`. Publication records
form a recoverable state machine linked to the exact run, result attempt,
generation, and artifact revision.

The first publication implementation handles one publication record for one
regular-file Artifact revision per transaction. Batch and directory
publication are outside this milestone. The default local publication keeps
one final payload copy. When source and destination are on the same filesystem,
publication:

1. opens and verifies the source with no-follow descriptor policy, hashes the
   exact bytes, and checks that the source and pool staging directory share
   `st_dev`;
2. commits a unique pending SQLite publication and an immutable move plan
   before the first destructive rename. The plan binds publication scheme and
   version, artifact ID, Run, result attempt, generation, existing artifact
   revision, expected size, `publication_content_sha256`, ordered contributor
   revision IDs, source safe-relative path, staging key, and final key;
3. reopens the bounded source, compares descriptor identity and `fstat`
   evidence with the plan, and uses descriptor-relative rename to move the
   payload into the unique pool-local staging namespace;
4. fsyncs the payload, source parent, and staging directory;
5. atomically renames the complete staging namespace to its unique final
   namespace and fsyncs the pool parent;
6. writes and fsyncs a canonical manifest temporary file, then performs a
   **manifest-last** atomic rename to `manifest.json` as the visibility marker
   and fsyncs the final directory;
7. verifies the marker and payloads, then uses one SQLite compare-and-set
   transaction to install private storage references, contributor rows, and
   the indexed publication state. The compare-and-set also verifies that the
   artifact generation and existing revision have not been replaced, and only
   then sets `published_at`.

The move occurs only after consumers that still require workspace paths, such
as current QC indexing, have completed. It creates no permanent workspace
duplicate. Recovery distinguishes empty pending, partial staging, unmarked
final namespace, valid marker with pending SQLite state, indexed, and corrupt
states. It searches only the three source/staging/final locations fixed in the
move plan and verifies identity, size, and checksum before retry or
finalization. A crash before the marker remains publicly invisible; a valid
marker with pending SQLite state can be verified and finalized; an indexed
row with missing or altered bytes fails closed. Cleanup or garbage collection
beyond files moved by the publication transaction is a later, explicit
operator action.

The logical artifact keeps its existing opaque `run://` URI, safe relative
path, output type, generation, and revision. The new full-content
`publication_content_sha256` is separate from, and cryptographically bound to,
the existing `RunArtifactRef.revision`; publication never replaces or
regenerates the current revision or generation. A private storage reference
binds the artifact to the Project, Run, approved pool, and publication without
exposing roots or storage coordinates publicly.

Downloads remain revision- and generation-bound. An unpublished or historical
artifact uses the fail-closed workspace resolver; an indexed publication uses
the private pool resolver. While a move is pending, the resolver may open only
the bounded source/staging/final locations from the persisted plan and verify
the exact descriptor, or return a stable retryable
`artifact_publication_in_progress` result. A download that already holds a
verified descriptor remains valid across rename.

The cross-filesystem case is a designed but unsupported contract in this
milestone. It detects the device mismatch before moving any source and returns
a stable `cross_filesystem_publication_unsupported` result. It does not copy,
download, or migrate real data. A later copier may use pool-local staging and
the same manifest-last recovery protocol, but enabling it requires a separate
decision about transient and final copy ownership.

### Artifact contributor provenance

An adapter receives the immutable snapshot's opaque SampleRevision context and
may declare one or several contributors for each sample-scoped extracted
artifact candidate. The platform verifies only that each declared revision
belongs to the current immutable run snapshot, rejects duplicates or
out-of-snapshot IDs, and persists the ordered association. It does not infer
contributors from filenames, output paths, QC labels, or adapter-private
scientific roles.

Multiple contributors represent pooled artifacts. A sample-scoped artifact
from a bound run must declare at least one contributor or publication fails
closed. A genuinely run- or Project-scoped output uses an explicit
adapter-declared `not_applicable` status. Only a historical Legacy artifact
without evidence uses `unresolved_legacy`; it is never assigned a convenient
sample.

Project/Sample/date/output-type queries use relational columns and joins, not
metadata JSON. Project filtering follows the publication/run's authoritative
Project foreign key. Sample filtering follows
ArtifactSample → SampleRevision → Sample, so `unresolved_legacy` and
`not_applicable` results never appear in a Sample query. Date filtering uses
the UTC `published_at` half-open interval `[from, to)`, and `output_type` is an
exact relational filter. A query returns the current generation/revision by
default or an explicitly requested exact generation/revision; it never mixes
superseded evidence silently.

### Migration and compatibility sequence

The persistence changes are additive revisions from the current
`20260717_08` head:

1. Project/Sample registry and snapshot/run sample bindings;
2. StoragePool/InputFile registry and snapshot/run input bindings;
3. artifact publication, private storage references, and ArtifactSample
   contributors.

Every persistence-changing stacked stage tests the complete cumulative upgrade
from `20260717_08` to that stage's exact head. A focused immediately preceding
revision transition test may be added but never replaces the revision 08
upgrade test. Each migration preserves existing bytes and rows, checks foreign
keys, and proves that legacy provenance remains unresolved. No migration
downloads inputs, calls a genome manager, copies artifact bytes, or migrates
real user data.

The Project/Sample stage adds logical identities and exact sample-revision
bindings without changing any adapter input-use path policy. The
StoragePool/InputFile stage adds regular-file and managed-binding capabilities,
but does not flip an input use merely because its leaf files can be registered.
The landing sequence for each input use is:

1. define and version its adapter-owned complete closure and platform binding;
2. implement adapter validation plus independent worker resolution and
   revalidation;
3. reconcile the exact execution closure and obtain exact-HEAD qualification
   for every supported route that activates the use;
4. advertise the capability for that exact workflow and adapter contract; and
5. only then require `managed_revision_v1` and return the stable physical-path
   rejection for that use.

Adding the registry migration changes the bulk-rnaseq execution implementation
manifest and therefore its exact execution closure. Protected evidence from
PR #154 is stale for the changed tree: updating a recorded digest or manifest
text alone does not requalify the runtime. Until new exact-HEAD qualification
lands against the reconciled closure, the affected runtime remains unavailable
and fails closed. An unqualified input use may remain transitional when its
exact adapter contract permits that provenance mode, but transitional catalog
provenance does not admit a stale execution runtime. This milestone does not
requalify the runtime and does not run the expensive rootful scientific gate.

OpenAPI is changed only for workflow-neutral binding selection and safe
read-only projection. When it changes, the checked-in specification and Orval
client are regenerated with the repository script and a repeated generation
must have zero drift.

### Normative compatibility matrix

The compact matrix below is test-pinned. Longer sections above define the
implementation detail.

| Situation | Normative contract |
| --- | --- |
| managed cutover identity | workflow ID + adapter contract version + input use capability |
| qualified managed input use | opaque revision only; physical path rejected |
| new catalog-aware unqualified input use | explicit transitional_unmanaged_v1; trusted-local path contract preserved |
| mixed managed/unmanaged snapshot | per-use provenance required and exact adapter contract must allow the mix |
| fully managed adapter | every execution-required input use qualified |
| pre-migration or historical unbound v1 snapshot | Legacy + unresolved |
| snapshot consumption | one SQLite transaction |
| snapshot replay | exact ordered per-use mode and evidence equality or fail closed |
| managed-input execution | ephemeral adapter view; resolved path not persisted |
| publication move unit | one regular-file Artifact revision |
| publication visibility | valid manifest AND indexed SQLite state |
| publication content identity | separate from existing artifact revision |
| bound sample-scoped artifact | one or more snapshot contributors required |
| sample query date axis | published_at UTC [from, to) |
| cumulative migration evidence | rev08 to exact stage head |

## Rejected alternatives

- **Infer one Sample per sample-sheet row or filename.** Rejected because row
  semantics differ by adapter and inference would create false provenance.
- **Resolve “latest SampleRevision” at run creation.** Rejected because the
  revision could change after validation.
- **Store registry fields inside adapter payload JSON only.** Rejected because
  identity, foreign keys, lifecycle, and required query axes belong in SQLite.
- **Use Project display names as workspace directories.** Rejected because
  display names are mutable and collisions would change physical identity.
- **Store pool roots or absolute paths in public metadata.** Rejected because
  it leaks operator-private topology and recreates arbitrary-path submission.
- **Cut over physical paths globally at one persistence stage.** Rejected
  because regular-file registration cannot express every adapter closure and a
  stage-wide switch would strand unqualified input uses.
- **Let platform code infer or interpret scientific closure layout.** Rejected
  because prefix alternatives, sidecars, index manifests, and tool/build
  compatibility remain versioned adapter contracts.
- **Overwrite a Sample or InputFile revision.** Rejected because historical
  runs must remain reproducible and queryable.
- **Hardlink, reflink, or retain a workspace copy during default publication.**
  Rejected because the required default has one durable payload copy and must
  not depend on aliasing semantics.
- **Copy automatically across filesystems.** Rejected for this milestone
  because it would perform real data migration and needs a separate recovery
  and ownership policy.
- **Treat Redis, a JSON manifest, or a directory layout as canonical state.**
  Rejected because lifecycle and relationships remain transactional in SQLite.

## Consequences and non-goals

The design adds relational state and operator configuration, so API and worker
processes must compose the same database and private pool mappings. Publication
also introduces explicit recoverable intermediate states rather than pretending
that a filesystem rename and SQL commit are one transaction.

Input management is deliberately progressive. A snapshot may contain managed
and transitional uses only when its exact adapter contract permits that mix and
records every use's provenance mode. The milestone performs no unconditional
global physical-path cutover, and regular-file support alone cannot justify a
fully managed adapter claim. Project/Sample delivery in particular preserves
current trusted-local adapter validation while refusing to manufacture catalog
provenance from paths or private payloads.

The milestone does not add authentication, a hosted administrator service,
arbitrary remote/object storage, public-data downloads, GenomePy execution,
directory/index bundle registration, automatic migration of existing bytes,
Genome Browser integration, email notification, workspace garbage collection,
bulk-rnaseq runtime requalification, merge automation, tagging, or release
behavior.
