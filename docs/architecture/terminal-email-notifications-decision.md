# Terminal Run Email Notifications Decision

Status: Accepted for the single-host laboratory deployment.

## Decision

HelixWeave sends one best-effort plain-text email attempt after a run reaches
`SUCCEEDED`, `FAILED`, or `CANCELLED`. Operator-owned administrator recipients
always participate. A member who created the run participates only when that
account has a private `notification_email` and has not disabled
`terminal_email_enabled`.

The notification address is contact data, not identity. It is not unique and
does not participate in login, account lookup, invitation, validation, or
password recovery. It is absent from principals, ordinary account lists, run
responses, validated snapshots, public events, and logs. A local operator sets
or clears it with `helixweave admin account set-notification-email` or
`clear-notification-email`; an authenticated member can only read or change the
boolean terminal-email preference. Administrators are configured through the
operator environment and do not use the member preference.

The run requester is captured atomically at actual create-run time in nullable
`runs.requested_by_user_id`. Legacy runs stay `NULL`, and idempotent replay does
not rewrite the original requester. This field does not create run ownership,
visibility filtering, or new RBAC semantics.

## Delivery and lifecycle semantics

The first version has no outbox, retry, startup replay, scheduler, or delivery
dashboard. A process crash can omit a message, and a worker-level replay around
the finalizer can duplicate one. That bounded delivery ambiguity is accepted
for a single-host laboratory; SQLite terminal state remains canonical.

Only the process that wins the durable terminal compare-and-swap attempts a
failure or cancellation email. Intent-only cancellation, cleanup failure,
losing CAS writers, already-terminal early returns, and repeated stop callbacks
do not send. Startup recovery and explicit administrator recovery notify only
for rows they actually change. A success email is attempted by the worker only
after artifact and QC indexing have reached a durable outcome. Hard timeout
interruptions preserve their primary failure behavior and do not extend the
deadline with SMTP.

SMTP errors are isolated after the terminal commit. They cannot roll back or
change a run, result generation, or worker result. A bounded run event records
only `sent`, `failed`, or `skipped`, the terminal status, recipient and metric
counts, and a fixed reason code. It cannot contain addresses, SMTP replies,
exception text, credentials, paths, or private payloads.

## Dynamic QC projection

Successful email reads only the current persisted QC generation through the
workflow-neutral repository contract. It sorts by stable durable coordinates
and renders at most 12 numeric rows using their stored display name, scope,
sample or experiment coordinate, value, unit, and optional QC flag. It does not
read a workspace or MultiQC file and has no ENCODE- or nf-core-specific keys.

Missing FRiP, any other optional metric, an empty generation, artifact indexing
failure, or QC indexing failure produces the same concise “no persisted QC
metrics” summary; delivery still proceeds. No files are attached.

## Operator transport boundary

`/etc/helixweave/notifications.env` is root-owned mode `0600`. systemd loads it
for the API and worker before dropping privileges, while the unit filesystem
sandbox makes the file itself inaccessible to those processes. The fixed
launcher forwards only the closed notification variable set to API and worker
commands. Scientific subprocess construction continues to use its independent
allowlist and receives none of these values.

Explicit administrator failure recovery reads the same file only when the
operator supplies `--notifications-env-file`; this is the supported deployed
invocation for `admin run fail`. The values remain in the command's private
composition mapping and are not exported to child processes. Diagnosis and
requeue do not read notification configuration.

The file defaults to an explicit disabled notifier. An enabled production
example is:

```dotenv
HELIXWEAVE_TERMINAL_EMAIL_ENABLED=true
HELIXWEAVE_TERMINAL_EMAIL_ADMIN_RECIPIENTS=admin@example.test
HELIXWEAVE_TERMINAL_EMAIL_FROM=helixweave@example.test
HELIXWEAVE_TERMINAL_EMAIL_APPLICATION_BASE_URL=https://helixweave.example.test
HELIXWEAVE_SMTP_HOST=smtp.example.test
HELIXWEAVE_SMTP_PORT=587
HELIXWEAVE_SMTP_TLS_MODE=starttls
HELIXWEAVE_SMTP_USERNAME=operator-supplied-user
HELIXWEAVE_SMTP_PASSWORD=operator-supplied-secret
HELIXWEAVE_SMTP_TIMEOUT_SECONDS=5
```

The administrator recipient value is a comma-separated list. Username and
password must either both be present or both be absent. `local_plaintext` is
accepted only with a loopback SMTP host and no credentials for controlled local
acceptance tests; it is not a LAN production mode.

The standard-library transport uses `EmailMessage` and `smtplib`, a validated
fixed sender and application origin, a bounded timeout, verified STARTTLS or
implicit TLS, and an explicitly loopback-only plaintext development mode.
Recipients are SMTP envelope addresses; the message uses an
`undisclosed-recipients` header so recipients do not see one another.
The total deduplicated envelope is capped at 32; an over-limit combination is
skipped with a fixed reason code instead of silently truncating administrators.

## Non-goals

This decision does not add HTML templates, attachments, Slack, webhooks, SMS,
marketing email, self-service address editing, address verification, email
login/reset, a general notification framework, adapter-specific QC behavior,
or nf-core email parameter changes.
