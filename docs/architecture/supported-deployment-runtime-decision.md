# Supported Deployment and Runtime Management Decision

Status: Accepted for PR #173 staging; shared CLI, migration, recovery, and
execution-identity wiring is deferred until PR #172 is merged.

## Context

HelixWeave v0.3.x can run a local demonstration from a source checkout, but
that path needs an editable Python install, Node/npm, a Vite development
server, and caller-managed Redis/process cleanup. It is not a supported
installation or upgrade boundary. The bulk RNA-seq adapter already has a
strict offline Nextflow/Docker runtime contract, while the ENCODE workflow and
platform source are still coupled by the existing workflow build identity.

PR #172 is changing run recovery, administrator operations, the migration
head, and bulk execution qualification. PR #173 therefore stages independent
deployment code on the last merged baseline and must not invent final values
for those contracts.

## Decision

### One supported topology

The only supported topology is a bounded hybrid for Linux x86_64 and WSL2
with systemd enabled:

- systemd runs one host-native FastAPI/Uvicorn process, one RQ worker, and one
  dedicated Redis 7 process;
- FastAPI serves the precompiled React application from the installed wheel on
  the same origin as `/api`; target hosts do not need Node, npm, Vite, Nginx,
  or a source checkout;
- the ENCODE Snakemake closure is an independently versioned immutable runtime;
- the bulk RNA-seq source, Nextflow, JDK, plugin, and container closure remains
  governed by the existing pinned runtime contract and uses a dedicated
  rootless Docker daemon/socket owned by the service account; and
- SQLite remains the canonical lifecycle and result-metadata store. Redis/RQ
  remains a queue boundary only.

References are prepared outside HelixWeave by an administrator and registered
through the existing private configuration/CLI contracts. Deployment code
never downloads, copies, deletes, or publishes their paths.

Compose, Kubernetes, Slurm, PostgreSQL, Nginx, Singularity, multi-host, and a
fully containerized platform are not supported alternatives.

### Authority and manifest model

There are two versioned JSON documents with non-overlapping authority:

1. An immutable `helixweave-deployment-bundle-v1` manifest indexes every
   bundled byte and points to content-addressed component contracts.
2. An immutable generation of `helixweave-deployment-state-v1` owns only the
   active, previous, and staged identities for platform, ENCODE runtime, and
   bulk RNA-seq runtime.

The current state is selected by one atomically replaced relative symlink. A
new generation records all three component slots, so a platform-only or
runtime-only switch is still one atomic state transition.

Existing contracts remain authoritative:

| Fact | Authority referenced by the deployment manifest |
| --- | --- |
| distribution/version and wheel | built wheel metadata plus bundle file index |
| frontend and API compatibility | packaged frontend asset manifest and canonical OpenAPI digest |
| database graph | admitted migration execution inventory |
| ENCODE execution | existing workflow build/execution identity after PR #172 integration |
| bulk runtime, nf-core, Nextflow, containers | existing bulk runtime identity and qualification resources |
| reference compatibility | adapter-owned reference binding contract versions |

The bundle manifest does not copy private reference coordinates or restate
facts owned by these documents. A contract identity appearing twice is a
reference, and verification must prove it resolves to the indexed bytes.

Plain canonical tar is the single admitted transport. `manifest.json` is the
first member; all remaining members are sorted regular files under `payload/`.
Absolute/traversal paths, symlinks, hard links, devices, FIFOs, duplicate or
case-colliding paths, unexpected files, writable payload modes, mutation
races, and hash/size mismatches fail closed.

### Supported host layout

The helper admits only these fixed roots:

```text
/opt/helixweave/
  releases/platform/<identity>/
  runtimes/encode/<identity>/
  runtimes/bulk-rnaseq/<identity>/
  operator/                         # helper environment, updated interactively
/etc/helixweave/
  platform.env
  secrets.env
  reference-profiles.yaml           # external coordinates, never in manifests
/var/lib/helixweave/
  operator/{ingress,state/generations,transactions}/
  database/{platform.db,backups/}
  workspaces/
  artifacts/
  docker-rootless/
/var/log/helixweave/
/run/helixweave/
```

Release and runtime directories are root-owned and immutable. Configuration,
secrets, SQLite, workspaces/artifacts, logs, runtime state, and database
backups are outside releases. Activation and rollback never move or delete
user data. Bulk runtime bindings receive the resolved version directory,
because the current bulk path-admission contract correctly rejects symlink
ancestry.

### Privileged operator boundary

The administrator runs the bootstrap locally through interactive sudo only
when installing or updating the helper. The bootstrap never invokes sudo and
never reads a password. It installs a self-contained, content-addressed,
root-owned operator closure with fixed systemd and tmpfiles templates and an
uninstall boundary. A stable root-owned launcher and sudoers allowlist are
installed before one relative `current` pointer atomically selects the fully
verified closure; product Python packages are never imported by the privileged
launcher.

Subsequent automation may invoke only `sudo -n /usr/libexec/helixweave-operator`
with a fixed positional grammar. No arbitrary path, command, shell, script,
environment assignment, unit, or signal is accepted. A stage bundle path is
derived from component and content identity beneath fixed ingress.

Service operations bind the allowlisted unit to deployment and task identity,
systemd MainPID/cgroup, executable device/inode, exact command-line digest,
and named socket device/inode. Status, stop, and cleanup must present the
previous service identity. The regular administrator is not added to the
host Docker group; bulk uses the dedicated service account's rootless socket.

The phase-A backend deliberately fails closed. Host systemctl, migration, and
recovery composition is connected only after PR #172 is merged.

The same closure installs a fixed `/usr/libexec/helixweave-gate-cleanup`
launcher and a narrow `sudo -n` grammar. Every cleanup plan binds both that
stable launcher and the exact root-owned content-addressed closure/backend it
selects, so a concurrent helper update cannot silently change cleanup
semantics. The stable launcher resolves the backend from the plan's validated
closure identity, not from the mutable `current` pointer. Its phase-A
implementation validates the closed
task/plan/executor/closure/backend identity grammar and then fails closed; process
revalidation and bounded deletion are connected only with the post-merge Gate
backend. This ensures the cleanup boundary is present without pretending that
an unimplemented host mutation succeeded.

### Frontend delivery

One reviewed Vite production build is copied into package-owned static data.
Its canonical manifest records only derived frontend facts: frontend version,
lock digest, OpenAPI digest, entrypoint, and the sorted file inventory. It has
no timestamp, host path, environment, secret, or runtime/reference data.

The deployment ASGI wrapper verifies the manifest before serving bytes. It
serves only allowlisted files, does not list directories or follow symlinks,
returns real 404 responses for missing assets, and uses `index.html` only for
safe SPA navigation. Hashed assets may be immutable-cached; `index.html` is
not cached. The verified platform release supplies the expected public API
contract digest, and the wrapper rejects a missing or mismatched digest before
calling the API factory, so compatibility failure cannot trigger migration or
recovery writes.

A restrictive production Content Security Policy is deferred until a nonce
strategy for CodeMirror's generated styles has a real browser smoke test. The
staging wrapper keeps the non-CSP response hardening headers and does not add
`unsafe-inline` or `unsafe-eval` merely to make the current editor render.

### Upgrade, rollback, and database boundary

Staging verifies the complete component before it can enter a state slot.
Activation re-derives both provided identities and compatibility requirements
from the indexed native contracts—requirements are never bundle-manifest
claims—then checks database-head compatibility and atomically changes the
state pointer. An interrupted operation leaves either the old pointer active or a
complete new generation plus a pending transaction receipt; doctor reports
partial, orphaned, or pending state without deleting it.

Platform and scientific runtimes switch independently. A runtime upgrade that
the active platform does not admit is rejected. Previous content, runtime
closures, evidence, and backups are retained until an explicit future
retention action; this milestone provides no automatic pruning.

Before any database migration, the final operator backend must verify API and
worker write services are stopped, create and verify a mode-0600 SQLite backup,
run the admitted migration graph, recheck integrity/head compatibility, and
only then start the new deployment. Rollback is allowed only when the observed
schema and active runtime combination is declared compatible. An irreversible
migration requires restoration of its controlled backup; it is never described
as a lossless Alembic downgrade.

Phase A only performs pinned-file-descriptor, read-only SQLite inspection. It
rejects symlinked, multiply linked, group/world-writable, foreign-owned, or
non-quiescent WAL/SHM inputs. Even with a root-owned write-stop witness and an
identity-bound schema admission provider, backup creation remains unavailable
until the privileged backend can atomically consume the one-use witness.

### CLI and diagnostics

The installed command surface will be `install`, `status`, `doctor`, `verify`,
`upgrade`, and `rollback`, with canonical JSON envelopes, stable exit classes,
and path/secret/exception redaction. `upgrade` always names exactly one of
platform, ENCODE runtime, or bulk RNA-seq runtime. `rollback` selects an
explicit compatible previous identity.

`status` reports selected identities and service state. `doctor` performs
bounded read-only configuration, permission, SQLite, Redis, worker, frontend,
runtime, and reference-readiness checks. `verify` performs full ownership,
mode, file hash, contract, and compatibility verification. It does not silently
repair state.

### Deployment Gate operator flow

The deployment Gate is a bounded release check, not a product runtime or a
general qualification framework. A request fixes exact Git HEAD, release and
runtime identities, task-owned checkout/environment/Redis/Docker socket/PID
and evidence directories, cache preflight, and disk/performance limits. Cleanup
is generated and validated before dispatch and binds the recorded process and
socket witnesses. One dispatch consumes one approval and stops at the first
failed stage. The redacted receipt names that stage. Cleanup preserves the
task Docker data-root/images and historical evidence.

During staging, the public preparation command uses an unavailable observer,
so caller-supplied identities cannot become verified Gate evidence. The final
observer and stage verifier must derive Git HEAD, deployment identities, and
runtime/qualification identities from the canonical checkout, verified bundle
manifests, and the contracts merged with PR #172; root-owned marker files are
not an additional version authority.

## Deferred integration boundary

Until PR #172 is merged, this staging branch does not change:

- `src/encode_pipeline/cli/admin.py`, `cli/local_platform.py`, or the primary
  dispatcher;
- migration inventory, Alembic head, run-recovery persistence/services, or
  worker recovery composition; or
- bulk execution identity, persistence contract, qualification, or scientific
  workflow behavior.

After merge, the latest `main` is merged into staging without rebasing. Final
CLI/doctor/service composition, database-head policy, ENCODE runtime identity,
bulk qualification, packaging trials, and Gate dispatch are then derived from
that merged baseline.

## Consequences

The topology is intentionally narrow and locally operable. It removes Node and
source checkout requirements from target hosts, keeps scientific closures
independently replaceable, and gives automation a noninteractive but bounded
operator surface. It also means unsupported hosts/topologies fail early, bulk
requires rootless Docker prerequisites, and schema-incompatible rollback may
require an explicit database restore rather than a pointer switch.
