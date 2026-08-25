# LAN Authentication and Roles Decision

Status: Accepted for staged implementation. Stage A is based on
`origin/main` commit `859d1cea1c36da38c1066e789aa3c2bd75106b03` (tree
`d777b43e1351baa4bc5d847fcb03b0ec5ad4ad3e`). Persistence, HTTP, frontend,
deployment configuration, and public-contract wiring remain gated on PRs #172
and #173 landing in `main`.

## Context and scope

HelixWeave is deployed to one trusted laboratory LAN. People in that laboratory
share projects, runs, artifacts, and QC results. Authentication establishes who
performed an action and separates routine scientific work from operator-style
administration; it does not create tenant isolation.

The complete roles are `administrator` and `member`. There are no dynamic
permissions, inherited roles, per-project ACLs, API tokens, external identity
providers, email identity, invitations, or password-reset links. The initial
administrator and password resets are local operator operations, never anonymous
HTTP operations. No default account, password, hash, or secret is created.

## Identity and password decision

A user has a random path-safe `usr_` identity, a unique normalized login name,
one fixed role, enabled/disabled status, and aware-UTC created, updated, and
password-changed times. Login names are 3--64 ASCII characters, case-folded
after trimming outer ASCII spaces, begin with a letter, and then contain only
letters, numbers, `.`, `_`, or `-`. This deliberately avoids
Unicode-confusable account identifiers; it is not a display-name policy.

Passwords are 15--128 Unicode code points and at most 512 UTF-8 bytes. They are
not trimmed, normalized, or subject to composition rules. The upper bound is
checked before the expensive hash operation; the lower bound follows the
password-only minimum in [NIST SP 800-63B][nist-passwords] and still accepts
passphrases longer than its recommended 64-character acceptance floor. A later
compromised-password list would need its own maintained data and update boundary
and is not fabricated in this PR.

The concrete hashing adapter directly uses `argon2-cffi` 25.1.x and stores its
standard PHC string. The explicit Argon2id v1 profile is `m=65536 KiB`, `t=3`,
`p=4`, a 16-byte random salt, and a 32-byte hash: the RFC 9106 low-memory profile
and the library's current ordinary-platform defaults. The adapter bounds custom
parameters, validates the complete stored PHC and its encoded m/t/p/salt/hash
cost before invoking Argon2, and uses the [library's high-level API][argon2-api].
One current profile and at most two explicitly configured legacy profiles are
accepted during a planned rotation. Every bounded attempt verifies one hash in
every accepted profile, using a process-local dummy in slots without an account
hash, so unknown, disabled, and mismatched accounts retain the same KDF work
shape. A successful enabled-account login uses `check_needs_rehash()` and
upgrades a legacy profile atomically with session creation in Stage B. Production
composition must run this blocking work outside the async event loop and bound
concurrent KDF work before exposing it to LAN traffic.

`argon2-cffi` is Production/Stable, MIT licensed, supports Python 3.10 through
the repository's current and next runtime versions, and publishes a universal
top-level wheel backed by ABI3 bindings for maintained platforms
([package evidence][argon2-pypi]). Its binding/CFFI dependency chain has
permissive licenses. `pwdlib` adds multi-hasher migration machinery that a new,
single-scheme store does not need; bcrypt is retained only as a legacy option by
current [OWASP guidance][owasp-passwords]; and Passlib's latest release is too
old for this new Python 3.13-facing boundary.

## Opaque server-side sessions

Successful login creates independent 32-byte values with Python `secrets`: an
opaque bearer session and a synchronizer-CSRF value. Only domain-separated
SHA-256 digests are persisted. A digest is sufficient here because the input is
uniformly random with 256 bits of entropy; it is not a password hash. Raw values
exist only while setting or validating browser cookies and must not enter model
reprs, SQLite, API bodies, logs, exceptions, or audit records.

SQLite is the session authority. Each record binds a token digest to a user,
CSRF digest, issue time, absolute expiry, and optional revocation time/reason.
The first version uses an eight-hour absolute lifetime without sliding renewal;
the code admits an explicit five-minute to seven-day bounded policy for final
deployment configuration. Every request resolves the current account, so a
disabled account fails immediately. Login rotates identifiers; logout revokes
one record before clearing cookies; password reset and account disable revoke all
of that user's records. Expired and revoked values cannot be replayed.

No JWT or client-readable session payload is used. [Starlette's built-in
session middleware][starlette-session] is a signed client-side cookie, not this
record model. Maintained general session packages still require a custom SQLite
store and do not remove
HelixWeave's status, role, revocation, transaction, or audit invariants. An
application-owned repository is therefore smaller and keeps SQLite as the sole
account/session truth; Redis remains only a queue boundary.

## Browser and CSRF policy

The authentication cookie is host-only, `HttpOnly`, `Path=/`, and
`SameSite=Lax`; its companion CSRF cookie is readable by the frontend but has no
authentication authority. Every cookie-authenticated unsafe method requires the
CSRF cookie, the identical custom `X-CSRF-Token` header, and the server-side
digest to agree.

Version one supports only a same-origin Web deployment: the browser reaches the
frontend and the API under one scheme, host, and port. A reverse proxy such as
Nginx may serve the frontend and forward `/api`, but the browser still sees a
single origin; multiple independent Web origins are not supported. There is no
canonical public-origin, trusted-proxy, or CORS allowlist configuration in this
version, and the implementation deliberately does not claim an exact configured
Origin/Referer check. Login accepts only the JSON request-body flow, and the
same-origin contract is enforced by the combination of JSON-only login parsing
(non-JSON submissions cannot create a session), the browser same-origin policy,
CORS remaining disabled by default, and the `SameSite=Lax` session cookie. An
explicit Origin allowlist and trusted-proxy model is future design work to be
done only if multi-origin deployment is ever supported.

An HTTPS deployment must use `Secure` cookies (and the `__Host-` name variant);
a separately named non-`Secure` policy exists only for loopback development.
Plain HTTP over a LAN cannot protect a bearer cookie and is not a production
security mode. `Max-Age` mirrors the server lifetime but never supersedes
server-side expiry. Logout clears both cookies with the same host/path policy.
Nothing is stored in `localStorage` or `sessionStorage`.

Login attempts are limited before Argon2 at two scopes: five attempts per
normalized login name and thirty per direct client in five minutes. The
thread-safe single-process limiter stores only bounded fingerprints, caps its
key sets, refuses new keys rather than evicting a still-limited identity, and
expires attempt windows. Capacity exhaustion is therefore bounded and
fail-closed. It is deliberately not Redis-backed and never becomes account
truth. PR #173 must define whether a proxy is trusted before forwarded client
addresses can be used. Multi-process aggregation and a bounded KDF semaphore are
final composition concerns, not reasons to create a second identity store.

## Authorization matrix

Authorization is enforced by FastAPI dependencies and services; hidden buttons
are presentation only.

| Surface | Anonymous | member | administrator |
| --- | --- | --- | --- |
| Login and non-sensitive health/setup state | allowed | allowed | allowed |
| Current session and logout | denied without a session | own session | own session |
| Shared projects, runs, artifacts, QC, publications, workflow schemas | denied | read | read |
| Validate/create/preflight/start/cancel, advisory Agent chat | denied | allowed | allowed |
| Artifact download | denied | allowed and audited | allowed and audited |
| Recovery fail/requeue and reference/storage/runtime mutation | denied | denied | allowed and audited |
| Member account create/enable/disable/reset/session revoke | denied | denied | allowed and audited |
| First administrator bootstrap and operator password recovery | no HTTP route | local CLI only | local CLI only |

The current validation endpoint persists a snapshot and is therefore treated as
a mutation. All existing business routers receive at least the authenticated
member dependency in Stage B. Administrator checks are applied again at the
service boundary for privileged operations. Authentication failures and
authorization failures use the existing structured, redacted error envelope as
401 and 403 respectively; scientific adapter payloads do not change.

The first administrator is unique bootstrap state, not a default user. With no
administrator, business routes fail closed and the service reports
setup-required without opening account-management HTTP APIs. The local CLI reads
passwords interactively with a non-echoing reader, never from argv, environment,
Agent chat, JSON output, or logs. The installed/noninteractive handoff, database
coordinate, service identity, and file permissions wait for PR #173. Version one
does not expose role mutation or anonymous reset tokens.

## Persistence and audit decision

After the upstream merge, one aggregate authentication repository will own
users, sessions, and security audit rows under one SQLite transaction and one
`BEGIN IMMEDIATE` read/modify/write boundary. Bootstrap, user mutation, password
rehash/reset, session creation/revocation, and their audit event commit or roll
back together. The migration descends from the then-current sole Alembic head,
creates zero users and zero secrets, preserves existing scientific data, and is
covered from the supported prior schema.

Automatic parameter rehash is not a password reset: after the KDF, Stage B must
re-read and compare the enabled account and old hash in its transaction, replace
only that hash, preserve `password_changed_at` and existing sessions, and create
the new session atomically. Reset changes `password_changed_at` and revokes every
session. Authentication admin commands return dedicated safe summaries/counts;
account, session, or verification dataclasses must never be passed to the
existing generic dataclass JSON serializer because `repr=False` is not a
serialization boundary.

Existing run, artifact, reference, and storage repositories each own and commit
their SQLite transaction; there is no caller-owned unit of work. Their required
security event therefore cannot be appended by a best-effort after-commit sink.
Stage B must either pass a closed `SecurityAuditEvent` into each affected atomic
repository mutation and insert it with the same SQLAlchemy session, or introduce
an explicitly reviewed shared-session unit of work. Rollback tests must prove
that the business mutation and audit row succeed or fail together.

Security audit is separate from public run events and the bounded Agent audit
sink. Its schema is closed: random event identity, aware UTC time, action,
outcome, explicit actor kind (user, local operator, or unauthenticated), optional
stable user identity, an opaque target identity derived from validated public
resource coordinates through its only bounded factory with a domain-separated
digest, and a closed reason-code enum. Action-specific actor/resource shapes are
enforced. It has no arbitrary context map. The source resource coordinate is not
retained, and login failure records no attempted username or account target. Raw
passwords, session/CSRF values or their authentication digests, cookies,
authorization headers, request bodies, private paths, addresses, exception text,
and environment values cannot fit the schema. The final SQLite ledger is
append-only. Exact recovery and runtime-management action names remain deferred
until PRs #172 and #173 define those operations rather than being guessed in
Stage A.

## Staging gate and consequences

Stage A contains only the decision records, workflow-neutral contracts,
password/session/cookie/CSRF/rate-limit primitives, dependency evidence, and
focused tests. It does not create an Alembic revision, alter OpenAPI or generated
clients, edit PR #172 recovery code, assume PR #173's installed layout, or open
PR #174.

Once both upstream PRs are merged, this branch will merge the latest `main`
without rebasing, re-audit the sole migration head and deployment boundary, and
then add the repository, CLI, HTTP, frontend, migration, packaging, and browser
integration. This sequencing avoids inventing migration ancestry, cookie/TLS
coordinates, secret handoff, or recovery permissions that the upstream work may
change.

[argon2-api]: https://argon2-cffi.readthedocs.io/en/stable/api.html
[argon2-pypi]: https://pypi.org/project/argon2-cffi/
[nist-passwords]: https://pages.nist.gov/800-63-4/sp800-63b.html
[owasp-passwords]: https://cheatsheetseries.owasp.org/cheatsheets/Password_Storage_Cheat_Sheet.html
[starlette-session]: https://www.starlette.io/middleware/#sessionmiddleware
