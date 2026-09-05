# HelixWeave v0.4.0 Local Trial Checklist

Use this checklist for the exact v0.4.0 tag or pre-tag candidate on one Linux
x86_64/systemd or systemd-enabled WSL2 host in a trusted laboratory LAN. It is
not a multi-tenant, high-availability, biological-validation, or performance
benchmark guide.

## Assets and identity

- [ ] Obtain `helixweave-0.4.0-py3-none-any.whl`,
      `helixweave-0.4.0.tar.gz`, and `SHA256SUMS` from one release/candidate.
- [ ] Run `sha256sum --check SHA256SUMS` and reject any additional official
      release asset.
- [ ] Confirm package, API, frontend, and `encode_pipeline.__version__` report
      `0.4.0` and the default registry lists the ENCODE-style and
      `bulk-rnaseq` workflows.
- [ ] Confirm neither artifact contains references, indexes, scientific
      environments, JDK, Nextflow, OCI layers, inputs/results, databases,
      backups, logs, or secrets.

The wheel contains the canonical precompiled frontend and deployment/operator
resources. The sdist additionally contains the exact release-owned ENCODE
workflow/profile/script source closure used by the bundle producer. Deployment
bundles are local operator transports and are not GitHub Release assets.

## Two independent distribution checks

From a directory outside any checkout or extracted source, create two new
Python 3.12 environments. Install the wheel in one; build a wheel from the
sdist and install that wheel in the other. Do not use an editable install,
`PYTHONPATH`, user site packages, or source-checkout files.

For each environment:

- [ ] `python -m pip check` succeeds.
- [ ] `importlib.metadata.version("helixweave")` prints `0.4.0`, no
      `encode-pipeline` distribution exists, and `encode_pipeline` belongs only
      to `helixweave`.
- [ ] `helixweave --help`, `python -m encode_pipeline --help`, all four
      compatibility commands, and `helixweave bundle --help` work.
- [ ] Exported OpenAPI matches the packaged canonical contract.
- [ ] Packaged frontend assets verify without Node or a checkout.
- [ ] A fresh file-backed SQLite database migrates to `20260827_15` and all
      packaged migration and adapter-contract resources are readable.

## Offline deployment transports

Follow the exact commands in the
[release checklist](../RELEASE_CHECKLIST.md#public-offline-bundle-producer).
The operator must explicitly provide:

- a complete wheelhouse and reviewed wheel lock for the platform bundle;
- the extracted v0.4.0 sdist, fixed micromamba executable, and complete
  reviewed archive cache for the ENCODE runtime bundle; and
- an already verified fixed Bulk runtime root plus explicit rootless Docker
  executable/socket coordinates for the Bulk RNA-seq runtime bundle.

- [ ] Every command runs offline and emits component, bundle path, bundle byte
      identity, and existing manifest identity.
- [ ] Missing, additional, malformed, symlinked, or drifted inputs fail closed
      without exposing private paths.
- [ ] Producer execution performs no sudo, download, reference discovery,
      host-path guessing, secret/database read, or runtime installation.

## Product and access journey

On the separately authorized supported host, bootstrap the operator once and
use the existing deployment CLI to install the platform bundle.

- [ ] `status`, `doctor`, and `verify` agree on the active platform, API,
      worker, Redis, PID/cgroup, socket, migration, packaged frontend, and API
      identities.
- [ ] Bootstrap the initial administrator, sign in, create a member, and verify
      administrator/member permissions, CSRF protection, sign-out, and terminal
      notification preferences. This trusted-LAN boundary is not multi-tenancy.
- [ ] Register/enable a prepared Reference Profile and verify its immutable
      revision identity is carried into validation and a created run.
- [ ] Create, validate, start, monitor, and inspect one deterministic ENCODE
      run; review events/logs, generation-bound artifacts/QC, publication, and
      safe download.
- [ ] Diagnose an eligible recovery case through the administrator surface
      without treating Redis as lifecycle truth.
- [ ] With no Bulk runtime binding, Bulk authoring/validation remains available
      while create/start fail closed with redacted `not_configured` evidence.
- [ ] With the admitted fixed Bulk runtime, verify the recorded protected
      evidence identity rather than rerunning the Gate during an ordinary trial.

## Responsive browser review

At 1440×900, 1024×768, 390×844, and 360×800 inspect login, workflow discovery
and detail, authoring/review, run history/detail, artifacts/QC, publication,
account controls, and unavailable/error recovery states.

- [ ] No page-level horizontal overflow, clipping, overlap, inaccessible
      control, hidden submitted form control, or broken technical identity.
- [ ] Keyboard navigation, focus, disclosures, log tabs, copy/download actions,
      and mobile touch targets remain usable.
- [ ] Screenshots and logs contain no private sample identifiers, paths,
      credentials, tokens, or environment values.

## Released database upgrade and rollback boundary

- [ ] Create the complete released v0.3.0 schema-08 fixture using verified
      v0.3.0 material, stop writers, make and verify a SQLite backup, then use
      installed v0.4.0 to migrate `20260717_08` to `20260827_15`.
- [ ] Verify the complete product record, bindings, result-generation state,
      publications, recovery/auth/notification data, and backup are preserved.
- [ ] Demonstrate pre-PONR restore before v0.4 writes.
- [ ] Do not describe v0.3.0 as a deployment rollback slot: it predates that
      product. Prove deployment rollback only between same-schema schema-15
      platform bundles.
- [ ] After a v0.4 writer starts, selecting schema-08 state fails closed; no
      lossless downgrade is claimed.

## Evidence and cleanup

Record exact commit/tree, asset checksums, package/API/frontend/Bulk/ENCODE and
deployment identities, database revisions, test/check URLs, and anything not
run. Controlled fixtures prove contracts, not biological validity or
production-scale throughput.

Stop only task-owned processes and remove only task-owned temporary files.
Preserve release evidence, migration backup, application data, shared Redis,
operator caches, references, Docker data-root, and unrelated checkouts. Never
use broad Docker or filesystem cleanup.
