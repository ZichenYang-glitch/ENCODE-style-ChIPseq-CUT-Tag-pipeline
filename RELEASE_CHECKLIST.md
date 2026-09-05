# Release Checklist

Use this checklist from a clean exact commit in the repository-locked
environments. Temporary packages, wheelhouses, caches, bundles, databases, and
host evidence stay outside Git. Passing the checklist makes a commit
release-ready; tagging and publication always require a separate authorization.

## v0.4.0 identity and assets

| Surface | Fixed identity |
| :--- | :--- |
| Product / Python distribution | HelixWeave / `helixweave` |
| Git tag | `v0.4.0` |
| Package, API, and frontend version | `0.4.0` |
| Python package | `encode_pipeline` |
| Primary CLI | `helixweave` |
| Compatibility CLIs | `encode-validate`, `encode-manifest`, `encode-dag`, `encode-worker` |
| Bundled workflows | ENCODE-style epigenomics and `bulk-rnaseq` |
| Database upgrade | released `20260717_08` to current `20260827_15` |
| Frontend toolchain | Node 22.23.1 / npm 10.9.8 |
| Bulk upstream | nf-core/rnaseq 3.26.0 |
| Bulk host isolation | Ubuntu `util-linux` `2.39.3-9ubuntu6.6`; `/usr/bin/unshare` SHA-256 `a23c8863860669003dc4660039fe642f5795c8c2195898ebc5d01afa1ac3d11c` |

The only GitHub Release assets are:

- `helixweave-0.4.0-py3-none-any.whl`;
- `helixweave-0.4.0.tar.gz`; and
- basename-only `SHA256SUMS` covering those two files.

Do not upload deployment bundles, source-trial archives, containers, JDK,
Nextflow, OCI layers, conda/mamba caches, references, indexes, input/results,
databases, backups, logs, or secrets. Deployment bundles are operator-derived
transports produced from the release assets and explicit local caches.

## Candidate freeze and generators

Record a clean exact candidate commit and tree. Before the tag exists, prove
that `v0.4.0` does not resolve.

```bash
git status --short --untracked-files=all
git rev-parse HEAD
git show -s --format=%T HEAD
git rev-parse -q --verify refs/tags/v0.4.0
```

With Node 22.23.1/npm 10.9.8 and the locked Python 3.12 environment, run the
formal generators in this order:

```bash
npm --prefix frontend install --package-lock-only --ignore-scripts --offline
npm --prefix frontend run openapi:regenerate
npm --prefix frontend run build
python3 scripts/package_frontend_assets.py
python3 scripts/generate_bulk_rnaseq_execution_manifest.py
python3 scripts/generate_migration_execution_inventory.py
```

Run each generator a second time and require no additional diff. OpenAPI may
change only `info.version`; the 190 generated TypeScript files (8 tag clients
and 182 models) may change only their version banner. Stop if routes, schemas,
operations, the npm dependency graph, resolved/integrity coordinates,
migration inventory, Bulk persistence/projection, scientific workflow bytes,
or pinned runtime assets change unexpectedly.

## Candidate packages and clean installs

```bash
python3 -m build --no-isolation --sdist --wheel --outdir dist
sha256sum dist/helixweave-0.4.0-py3-none-any.whl \
  dist/helixweave-0.4.0.tar.gz
```

Reject any unexpected distribution filename or forbidden sdist member. The
wheel must contain package data, migrations, deployment/operator resources,
and the canonical packaged frontend. The sdist must additionally contain the
exact ENCODE workflow-build source manifest: `workflow/`, `profiles/default/`,
the referenced top-level scripts, and
`docs/architecture/artifact-inventory.yaml`; it must not contain environments,
materialized conda/mamba prefixes, package caches, references, indexes,
scientific inputs/results, databases, logs, containers, or secrets.

Install the wheel and an independently sdist-built wheel in two new
outside-checkout environments. In each, with no user site or `PYTHONPATH`:

- run `python -m pip check`;
- confirm `metadata.version("helixweave") == "0.4.0"` and no
  `encode-pipeline` distribution metadata;
- exercise primary and compatibility CLIs and `helixweave bundle --help`;
- confirm both workflow IDs, exported OpenAPI, packaged frontend/API binding,
  migration/package data, and a fresh file-backed SQLite database at schema
  `20260827_15`.

## Public offline bundle producer

The producer never uses sudo, downloads assets, discovers host paths, or reads
references, secrets, runtime state, or user databases. Use absolute paths.

Create and review the platform lock from the complete candidate wheelhouse:

```bash
helixweave bundle wheel-lock create \
  --wheelhouse /absolute/wheelhouse \
  --output /absolute/platform-wheel-lock.json
helixweave bundle wheel-lock verify \
  --wheelhouse /absolute/wheelhouse \
  --lock /absolute/platform-wheel-lock.json
```

Build the three existing bundle formats from a release-shaped fresh directory:

```bash
helixweave bundle build --component platform \
  --wheel /absolute/helixweave-0.4.0-py3-none-any.whl \
  --wheelhouse /absolute/wheelhouse \
  --wheel-lock /absolute/platform-wheel-lock.json \
  --output /absolute/platform.tar \
  --scratch-root /absolute/empty-scratch

helixweave bundle build --component encode-runtime \
  --sdist-root /absolute/extracted/helixweave-0.4.0 \
  --micromamba /absolute/micromamba \
  --archive-cache /absolute/reviewed-conda-archive-cache \
  --output /absolute/encode-runtime.tar

helixweave bundle build --component bulk-rnaseq-runtime \
  --runtime-root /absolute/reviewed-bulk-runtime \
  --output /absolute/bulk-rnaseq-runtime.tar
```

Each receipt must identify its component, output path, bundle byte identity,
and existing manifest identity. Verify that the existing deployment
install/verify/upgrade admission accepts each output. The release assets do not
provide the explicitly operator-owned caches or Bulk runtime.

## Database and deployment acceptance

The release migration proof starts from the actually released v0.3.0 schema
head `20260717_08`, not a prerelease revision. Create a complete schema-08
product record with the v0.3.0 artifact/fixture, stop writers, make and verify a
backup, then use installed v0.4.0 code to migrate to `20260827_15`. Verify the
record, result-generation state, reference/input bindings, publications,
recovery/auth/notification data, and backup remain intact.

v0.3.0 predates the deployment-slot product and is not a rollback slot.
Exercise compatible platform upgrade/rollback only between same-schema
`20260827_15` bundles. Prove pre-PONR restore separately; after a v0.4 writer
starts, selecting schema-08 state must fail closed rather than claim a
lossless downgrade.

One separately authorized Linux x86_64/systemd or systemd-enabled WSL2 host
window is required: bootstrap the operator once; install/status/doctor/verify
platform; check API, worker, Redis, PID/cgroup/socket identities, packaged
frontend, and data retention; perform one same-schema platform
upgrade/rollback; perform one independent ENCODE transition/rollback; and for
Bulk verify only admission, rootless Docker/service/socket identity, and
existing runtime preservation. Do not repeat the full fault-injection matrix,
switch Bulk merely for symmetry, use a second host, soak, or load test.
The v0.4.0 Bulk qualification binds the supported Ubuntu host's
`util-linux` package revision `2.39.3-9ubuntu6.6` and exact
`/usr/bin/unshare` bytes. Do not downgrade the system package or substitute a
temporary executable to satisfy this coordinate.

## Automated and protected evidence

The exact functional candidate must have ordinary CI, Lint, Lock Check, both
fast shards, fast-checks, coverage, frontend, and browser E2E green. Local work
should use focused version, provenance, OpenAPI, frontend-package, release
distribution, bundle, migration, and documentation checks plus Ruff/format on
changed Python and `git diff --check`; do not duplicate the complete CI suite.

Only after packages, migration, ordinary CI, and the authorized host window
are green may the maintainer authorize one exact-HEAD combined dispatch:

```text
bulk_rnaseq_real_execution=true
bulk_rnaseq_expected_sha=<exact final functional HEAD>
```

That single dispatch supplies platform-real, ENCODE scientific-real,
container-smoke, and Protected Bulk evidence. Do not dispatch correction
HEADs, auto-rerun a failure, or describe controlled tiny/synthetic evidence as
biological validation or production-scale throughput.

## Date, merge, and publication boundary

After the combined Gate is green, obtain the actual release date. A date-only
commit may update only `CHANGELOG.md`, `CITATION.cff`, the release enforcement
test, and the v0.4.0 candidate evidence. Preserve an empty Unreleased section,
prove all API/frontend/scientific/Bulk/deployment closures are identical, and
rerun only exact-HEAD ordinary CI, Lint, and Lock Check.

Stop at the Draft human merge gate. Squash merge, annotated `v0.4.0` tag, tag
push, fresh-tag asset build, checksum creation, Draft GitHub Release, asset
upload/download verification, and publication each require their stated later
authorization. Never move, replace, rebuild, or reinterpret v0.3.0.
