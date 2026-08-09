# LAN Authentication Threat Model

Status: Stage A security model for PR #174. Revalidate the deployment and
recovery boundaries after PRs #172 and #173 merge.

## System and trust boundaries

Protected assets are account hashes and status, opaque browser sessions, CSRF
values, administrator authority, shared scientific records and artifact bytes,
private storage/runtime coordinates, and the security audit ledger.

The deployment has four relevant boundaries:

1. A browser crosses the laboratory LAN/TLS boundary into FastAPI.
2. HTTP adapters cross into workflow-neutral authentication and authorization
   services.
3. Services cross into canonical SQLite and, separately, the Redis/RQ execution
   queue. Redis never holds account or session truth.
4. A local operator crosses a host-login boundary into `helixweave admin`; this
   is the only first-administrator and emergency password-reset path.

Laboratory members are trusted to see the shared scientific catalogue, not to
administer identities, recovery, references, storage, or runtime. A malicious
web origin, compromised member browser, LAN peer, stolen database copy, malformed
request, and accidental operator disclosure are in scope. A fully compromised
host/root account, hostile laboratory tenant isolation, physical attacks, and
external identity-provider compromise are outside this feature's guarantees.

## Threats and required controls

| Threat | Required control and verification | Residual risk |
| --- | --- | --- |
| Password database theft or offline guessing | Argon2id with random salt and explicit versioned cost; no plaintext/reversible password; rehash after successful verification | Weak chosen passwords remain guessable; no blocklist is introduced in this stage |
| Account enumeration | Unknown, wrong-password, and disabled-account login produce the same status/body and perform one verification in every bounded current/legacy Argon2 profile; failed audit omits username and account target | Network and scheduler jitter cannot promise cycle-identical timing |
| Online guessing and resource exhaustion | Per-name and per-client pre-KDF rate limits, bounded inputs/key sets, bounded KDF concurrency at API composition | Process restart resets the Stage A in-memory limiter; multiple API workers need final composition review |
| Session fixation | Generate a fresh 256-bit identifier only after authentication and after any privilege transition; never accept a caller-selected identifier | A compromised browser can still use its current session |
| Session database disclosure | Persist only a domain-separated SHA-256 digest of a high-entropy opaque value; never include raw values in repr/API/audit/logs | A live bearer stolen from the browser remains usable until expiry or revocation |
| Replay after logout, reset, disable, or expiry | Server-side revocation and absolute expiry; load current user status on every request; revoke all sessions on reset/disable | Concurrent requests already authorized before a transaction commits may finish |
| LAN interception or cookie theft | HTTPS, host-only `Secure; HttpOnly; SameSite=Lax; Path=/` session cookie, no URL/JSON/storage exposure | A plain-HTTP LAN deployment is explicitly unsafe and must not be presented as protected |
| CSRF, including login CSRF | Unsafe authenticated requests require cookie + custom header + stored digest; exact Origin/Referer and JSON login checks in Stage B; rotate on login | XSS in the same origin can read the CSRF cookie and act as the user |
| XSS-driven credential persistence | Session value is HttpOnly and never placed in `localStorage`; frontend renders existing typed data paths | XSS can issue same-origin actions during the active browser session |
| Frontend-only authorization bypass | Backend member/admin dependencies plus service-layer checks; parameterized route matrix covers anonymous/member/admin | Incorrectly classified new routes remain a review risk; router inventory is a release gate |
| Disabled-account or stale-role use | Resolve account status and role from SQLite for every session; do not trust role in a cookie | SQLite availability is required for authenticated traffic, by design |
| First-admin takeover | No default credentials or anonymous management route; serialized local CLI bootstrap succeeds only when no user exists; setup-required mode fails closed | Host-login security and the PR #173 operator boundary remain prerequisites |
| Operator secret disclosure | Non-echoing password reader; no argv/env/chat/API transport; safe JSON summaries; hashes and raw values excluded from repr | Python cannot guarantee in-memory erasure of immutable strings |
| Audit/log exfiltration | Closed allowlisted audit schema and reason codes; target identities are domain-separated digests of validated public resource coordinates; no arbitrary exception/request context; append-only SQL protections; generic error envelope | Stable user and opaque target identities remain intentionally visible to authorized operators |
| Path/header/request injection | Fixed-length ASCII identifiers, normalized login names, bounded Unicode passwords, strict bounded Argon2id PHC parsing before KDF work, closed reason codes, no raw headers or private paths | Stage B adapters must not bypass constructors with unvalidated ORM rows |
| Migration or upgrade failure | Revision from the actual sole merged head; atomic migration; zero fabricated users/secrets; prior-schema preservation and downgrade tests | Rollback of installed binaries/configuration is owned by PR #173 |
| Proxy spoofing | Use the directly connected peer unless PR #173 explicitly establishes trusted proxies; exact configured origin | Source identity is weak on shared NAT/proxy networks and is not authorization evidence |

## Security invariants

- Anonymous access cannot reach projects, runs, artifacts, QC, Agent chat, or
  any mutation. Only login and non-sensitive health/setup state are public.
- A `member` can see shared laboratory data and perform routine run operations;
  only an `administrator` can mutate accounts, recovery, references, storage,
  or runtime state.
- Every unsafe cookie-authenticated request is rejected before service mutation
  if CSRF validation fails. A 403 is still presented even if the frontend hid
  the control.
- A session row never contains a raw bearer or CSRF value. Cookie expiry never
  overrides server expiry/revocation, and a disabled account is unusable without
  waiting for cookie expiry.
- Password reset and account disable revoke all sessions in the same aggregate
  transaction as the account change and audit evidence.
- Login results disclose setup-required state only as a global installation
  condition; otherwise account absence, disablement, and wrong password are
  indistinguishable.
- Passwords, raw authentication headers, cookies, tokens, authentication-token
  digests, private paths, email/SMTP details, environment values, and raw
  exceptions never enter public errors or security audit.
- There is no default administrator, role inheritance, project ACL, tenant key,
  bearer API token, or anonymous password-reset token.

## Stage A evidence and Stage B gates

Stage A tests must cover username and identity boundaries, the closed role/status
contract, Argon2id verify and parameter upgrade, identical rejected-login shape,
session generation/digest/domain separation/fixation/revocation/expiry, cookie
flags, synchronizer CSRF matching, rate-limit windows and memory bounds, audit
schema rejection, malicious input, dependency version/license, and package
metadata.

Stage B must add repository race/rollback and restart persistence tests; a
parameterized anonymous/member/administrator matrix over every route; local CLI
bootstrap/reset/enable/disable/revoke tests; login/logout/current-session,
expiry/CSRF/403 browser tests; migration preservation and zero-default tests;
OpenAPI/client regeneration twice; and wheel/sdist/no-Git-source trials. It must
also enumerate every newly merged PR #172 recovery operation and every PR #173
runtime mutation before finalizing their audit action names. Repository tests
must prove login rehash compare-and-swaps the old enabled-account hash without
changing `password_changed_at` or revoking existing sessions, while reset does
both. CLI tests must prove the generic dataclass serializer never receives auth
state and that password hashes, session/CSRF digests, and replacement hashes are
absent from stdout and stderr.

The existing business repositories commit their own transactions, so Stage B
must not append security evidence after a run/artifact/reference/storage mutation
has committed. Each audited mutation must share one SQL transaction with its
audit insert (by an event parameter at the existing aggregate repository boundary
or an explicitly reviewed unit of work), with failure-injection tests in both
directions.

Release review must separately examine (1) password/session/cookie/CSRF and
authorization security and (2) migration/API/frontend/packaging compatibility.
Any change in installed execution closure stops at the repository's Protected
Gate instead of being authorized by this design.
