# Supported Deployment and Runtime Management Decision

Status: Accepted for PR #173.

## Context

HelixWeave v0.3.x can run a local demonstration from a source checkout, but
that path needs an editable Python install, Node/npm, a Vite development
server, and caller-managed Redis/process cleanup. It is not a supported
installation or upgrade boundary. The bulk RNA-seq adapter already has a
strict offline Nextflow/Docker runtime contract, while the ENCODE workflow and
platform source are still coupled by the existing workflow build identity.

PR #172 established the current run-recovery migration inventory and bulk
execution qualification before this decision was finalized. PR #173 consumes
those native contracts without creating a migration head, persistence
capability, or competing version authority.

## Decision

### One supported topology

The only supported topology is a bounded hybrid for Linux x86_64 and WSL2
with systemd enabled:

- systemd runs one host-native FastAPI/Uvicorn process, one RQ worker, and one
  dedicated Redis 7 process;
- FastAPI serves the precompiled React application from the installed wheel on
  the same origin as `/api`; target hosts do not need Node, npm, Vite, Nginx,
  or a source checkout. Service activation never invokes pip: a release
  producer expands a hash-locked, wheel-only CPython 3.12/Linux x86_64
  dependency closure before bundling it;
- the ENCODE Snakemake closure is an independently versioned immutable runtime
  containing the workflow source, a pinned Linux x86_64 micromamba frontend
  admitted against the supported glibc loader, and every
  package archive named by the checked-in explicit environment locks;
- the bulk RNA-seq source, Nextflow, JDK, plugin, and container closure remains
  governed by the existing pinned runtime contract and uses a dedicated
  rootless Docker daemon/socket owned by the service account; and
- SQLite remains the canonical lifecycle and result-metadata store. Redis/RQ
  remains a queue boundary only.

The v0.4.0 supported Ubuntu host coordinate binds `util-linux`
`2.39.3-9ubuntu6.6` and the exact qualified `/usr/bin/unshare` bytes. Operators
must not downgrade the system package or substitute a temporary executable to
meet that coordinate.

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
| ENCODE execution | `sha256-tree-v1` over the independently packaged workflow, profile, script, and artifact inventory closure |
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
  secrets.env
  notifications.env                 # root:root 0600; API/worker plus explicit admin fail
  reference-profiles.yaml           # external coordinates, never in manifests
/var/lib/helixweave/
  deployment/
    current                         # one relative generation pointer
    generations/<state-identity>/
      deployment-state.json
      platform.env
  operator/
    ingress/
    state/{operation.lock,transactions/}
    transactions/{active.json,history/}
    database-backups/
  database/live/platform.db
  scientific/encode/<runtime-identity>/{runner,conda-envs,mamba-root}/
  workspaces/
  artifacts/
  docker-rootless/
/var/log/helixweave/
/run/helixweave/
```

Release and runtime directories are root-owned and immutable. Configuration,
secrets, SQLite, materialized scientific environments, workspaces/artifacts,
logs, runtime state, and database backups are outside releases. Activation and
rollback never move or delete user data. Bulk runtime bindings receive the
resolved version directory, because the current bulk path-admission contract
correctly rejects symlink ancestry.

### Privileged operator boundary

The administrator runs the bootstrap locally through interactive sudo only
when installing or updating the helper. The bootstrap never invokes sudo and
never reads a password. It installs a self-contained, content-addressed,
root-owned operator closure with fixed systemd and tmpfiles templates and an
uninstall boundary. The stdlib-only bootstrap and its exact sources are shipped
inside the verified sdist release asset, so this step needs neither a Git
checkout nor an installed HelixWeave package. A stable root-owned launcher and
sudoers allowlist are installed before one relative `current` pointer
atomically selects the fully
verified closure; product Python packages are never imported by the privileged
launcher.

Subsequent automation may invoke only `sudo -n /usr/libexec/helixweave-operator`
with a fixed positional grammar. No arbitrary path, command, shell, script,
environment assignment, unit, or signal is accepted. A stage bundle path is
derived from component and content identity beneath fixed ingress.

Service operations bind the allowlisted unit to deployment and task identity,
systemd MainPID/cgroup, executable device/inode, exact command-line digest,
and named socket device/inode. Status is a bounded observation that can report
a stopped unit or its newly observed identity after an automatic restart;
stop and cleanup require the prior identity and revalidate it before acting.
The regular administrator is not added to the host Docker group; bulk uses the
dedicated service account's rootless socket.

The same closure installs a fixed `/usr/libexec/helixweave-gate-cleanup`
launcher and a narrow `sudo -n` grammar. Every cleanup plan binds both that
stable launcher and the exact root-owned content-addressed closure/backend it
selects, so a concurrent helper update cannot silently change cleanup
semantics. The stable launcher resolves the backend from the plan's validated
closure identity, not from the mutable `current` pointer. The backend opens the
root-owned canonical plan through fixed directory descriptors, revalidates the
bound PID start time, executable, command line, PID file, and socket identities,
then performs bounded TERM/KILL through pidfds where available. Any mismatch or
non-terminal process stops cleanup before deletion. Successful cleanup removes
only the fixed task descendants named by the verified plan and preserves the
Docker data root, images, and historical evidence.

The API, scientific worker/database preparer, and candidate action/materializer
use distinct system identities. Only the API identity can read `secrets.env`.
systemd reads root-owned `notifications.env` before dropping privileges for
the API and worker; both units then hide the file, and the fixed launcher
forwards its closed variable set only to those two commands. An explicitly
authorized local `admin run fail` may read the same file directly for its one
terminal mutation. Scientific child processes inherit none of the notification
or SMTP values. The API and worker
share the narrow data group needed for SQLite, workspaces, artifacts, and Redis.
The candidate identity is not a member of that
group and its units hide configuration, references, the live database, Redis,
Docker, workspaces, artifacts, and secrets. Candidate receipts therefore
describe candidate-native facts only; the root closure derives the observed
schema and migration decision after the writers have stopped. Root otherwise
only selects fixed content-addressed candidates, writes identity-bound requests,
checks bounded canonical receipts, freezes admitted trees, switches state, and
controls the allowlisted units. Neither the root helper nor its launcher imports
or executes a candidate wheel or scientific binary as root.

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

A restrictive production Content Security Policy is deferred until the current
AJV runtime schema compilation is replaced or precompiled and covered by a real
browser smoke test. The controlled frontend bundle currently uses dynamic code
generation for the New Run schema validator, so a policy without
`unsafe-eval` would break that supported flow. The deployment wrapper keeps the
non-CSP response hardening headers and does not add `unsafe-inline` or
`unsafe-eval` merely to claim a policy is present.

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

An ENCODE runtime is materialized offline directly at its final
content-addressed prefix, because conda packages may embed that absolute
prefix in shebangs or binary metadata. The service account consumes only the
bundled package index and local archives; network access is disabled. On
success the root boundary revalidates a bounded canonical tree inventory,
breaks external hard-link authority, permits only internal relative symlinks,
and freezes the tree root-owned. On failure the entire final-prefix tree is
renamed to task-scoped retained evidence rather than deleted or reused.

Materialization uses the indexed micromamba binary only to create the locked
runner and full-address-hash rule environments. The frozen runtime retains that
binary as evidence and exposes a generated, read-only `conda` compatibility
command limited to Snakemake 8.30's information, version, and channel-priority
probes, plus an activator that accepts only an existing 32-character full-hash
environment beneath the admitted prefix. Environment creation and all other
caller-selected conda operations are rejected during workflow execution.

Bulk RNA-seq admission is split the same way. Static runtime verification is
bytes-only and does not contact Docker. During activation the root boundary
binds the dedicated rootless Docker service identity, fixed client bytes, and
fixed socket/kernel identity into a path-free request. The unprivileged bulk
preparer verifies the canary and each expected image, loads only images that
are missing through an already verified archive file descriptor, then repeats
the complete identity check. An invalid, unavailable, or identity-mismatched
image fails before loading any archive; a retry loads only the still-missing
set. The operator never removes or prunes images or the Docker data root.

The supported API and worker open only an existing database at the exact sole
head admitted by the packaged migration inventory; they never invoke Alembic.
Migration is a separate operator transaction with this fixed order:

1. stop API and worker writers and prove their bound PID/cgroup/socket
   identities are gone;
2. consume the one-use root-owned write-stop witness and create a verified,
   root-owned immutable SQLite backup slot (mode `0500`) and file (mode
   `0400`);
3. run only the admitted migration graph as the unprivileged service account;
4. recheck database integrity and the exact target head;
5. atomically select the verified state/config generation; and
6. start the services and bind their new identities.

The durable transaction journal is fsynced before each side effect. Before a
new writer is started, a failed migration may restore the verified backup and
restart the old compatible release. Starting any new writer is the point of no
return: automatic backup restoration is then forbidden because it could lose
new writes. A later rollback is accepted only when the observed schema and
runtime combination is declared compatible; otherwise services remain stopped
with recovery required. An irreversible migration is never represented as a
lossless Alembic downgrade.

Database inspection and backup pin the live file by descriptor, reject
symlinked, multiply linked, group/world-writable, foreign-owned, or
non-quiescent WAL, SHM, and rollback-journal inputs, and bind the backup and
receipt to the consumed witness and source inode. Backups, failed transaction
evidence, and the original database are not automatically deleted.

### CLI and diagnostics

The installed command surface is `install`, `status`, `doctor`, `verify`,
`upgrade`, and `rollback`, with canonical JSON envelopes, stable exit classes,
and path/secret/exception redaction. `upgrade` always names exactly one of
platform, ENCODE runtime, or bulk RNA-seq runtime. `rollback` selects an
explicit compatible previous identity.

`status` reports selected identities and service state. `doctor` performs
bounded read-only configuration, permission, SQLite, Redis, worker, frontend,
runtime, and reference-readiness checks. `verify` performs full ownership,
mode, file hash, native-contract, configuration, schema, and compatibility
verification. Static verification does not require API, worker, Redis, or
Docker to be running and does not treat an externally prepared reference that
has not yet been registered as corrupt deployment content; those operational
conditions belong to `doctor`. Neither command silently repairs state.

### Deployment Gate operator flow

The deployment Gate is a bounded release check, not a product runtime or a
general qualification framework. A request fixes exact Git HEAD, release and
runtime identities, task-owned checkout/environment/Redis/Docker socket/PID
and evidence directories, cache preflight, and disk/performance limits. Cleanup
is generated and validated before dispatch and binds the recorded process and
socket witnesses. One dispatch consumes one approval and stops at the first
failed stage. The redacted receipt names that stage. Cleanup preserves the
task Docker data-root/images and historical evidence.

The preparation command uses a fixed filesystem observer. It derives Git HEAD,
installed-wheel provenance, deployment identities, process/socket witnesses,
cache objects, cleanup executor identity, and runtime/qualification identities
from the canonical checkout, root-owned policy, verified bundles, and native
contracts. Caller-supplied PID, path, cache, or performance evidence is not
accepted. The cleanup backend is part of the installed operator closure and is
covered by fault-injection tests; actual one-use runner dispatch, human
approval, rootless Docker, Redis, cleanup execution against real processes, and
real-host receipts remain Protected Gate evidence and are not emulated by local
unit tests.

## Consequences

The topology is intentionally narrow and locally operable. It removes Node and
source checkout requirements from target services, keeps scientific closures
independently replaceable, and gives automation a noninteractive but bounded
operator surface. It also means hosts without the admitted CPython 3.12 ABI
(at or above the supported 3.12.3 security-patch floor),
Linux x86_64/glibc boundary, systemd, or WSL2 systemd support fail early; bulk
requires rootless Docker prerequisites; and schema-incompatible rollback may
require an explicit database restore rather than a pointer switch.
