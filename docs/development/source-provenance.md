# Python Source Provenance

HelixWeave's canonical development checks start through a project-specific
clean bootstrap:

```bash
python3 -I -S scripts/checkout_bootstrap.py \
  --repository-root . verify-checkout
```

`-I -S` prevents `PYTHONPATH`, user site packages, `.pth` files,
`sitecustomize`, and `usercustomize` from running before the bootstrap. The
standard-library-only bootstrap then audits the selected environment without
executing unknown `.pth` content. Only after provenance succeeds does it add
the audited paths and import pytest or a product-owned command.

## Controlled modes

Checkout mode requires an explicit repository root and one current editable
installation. The distribution identity and version must match
`pyproject.toml`; `direct_url.json` must name that exact checkout;
`top_level.txt` and `RECORD` must prove ownership of `encode_pipeline`; and
the sole recorded plain `.pth` mapping must resolve to that checkout's
`src`. A sibling worktree is rejected even at the same commit or through a
symlink alias.

Every physical distribution that claims the namespace is counted. Duplicate
claimants are rejected regardless of spelling, version, editable status, or
whether both point at the current checkout. Missing or malformed ownership
inventory, an orphan namespace `.pth`, an editable finder, or an executable
product `.pth` also fails closed. Unrelated executable `.pth` lines are never
executed and cannot make product paths importable.

Installed-artifact mode is only for isolated wheel and rebuilt-sdist tests:

```bash
/path/to/isolated/bin/python -I -S scripts/checkout_bootstrap.py \
  --repository-root /path/to/repository installed-artifact
```

It requires exactly one non-editable `helixweave` distribution in that
interpreter's isolated `site-packages`. Its `RECORD` must own
`encode_pipeline/__init__.py`, and the resolved package origin must remain
inside that installed artifact closure. Source mappings and external editable
paths are rejected. Ordinary installed console and module entry points do not
require a Git checkout and do not invoke this development guard.

## Canonical commands

Use the same bootstrap for the maintained source-checkout entry points:

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest test -ra
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  openapi --output frontend/openapi.json
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  validate --config config/config.yaml
python3 -I -S scripts/checkout_bootstrap.py --repository-root . \
  local-platform --doctor
```

The pytest command disables automatic plugin discovery before importing
pytest. It explicitly loads only the repository-required pytest-cov plugin and
rejects `PYTEST_PLUGINS`, `PYTEST_ADDOPTS`, and caller-supplied plugin flags.
CI and frontend OpenAPI regeneration use these commands, not merely a separate
provenance probe.

## Trust boundary and failures

Protection begins when the selected operating-system Python executable runs
this repository's bootstrap with both `-I` and `-S`. It trusts the OS loader,
that interpreter and its standard library, and the reviewed bootstrap and
guard files. Environment provisioning and interpreter selection happen before
this boundary. This mechanism does not qualify a compromised interpreter,
standard library, or native loader injection; it also cannot govern callers
that deliberately bypass the documented source-checkout commands. Direct
installed-user CLI entry points remain outside checkout governance by design.

Failures exit with status 2 and one public-safe line:

```text
source provenance check failed [reason_code]: actionable guidance
```

Stable reason codes include `bootstrap_startup_unsafe`,
`bootstrap_source_invalid`, `repository_root_invalid`,
`repository_metadata_invalid`, `source_root_invalid`, `module_not_found`,
`module_origin_mismatch`, `module_search_location_mismatch`,
`distribution_missing`, `distribution_metadata_invalid`,
`distribution_claimant_conflict`,
`distribution_identity_mismatch`, `distribution_version_mismatch`,
`distribution_source_mismatch`, `namespace_ownership_unproven`,
`namespace_mapping_conflict`, `editable_mapping_invalid`,
`pth_mapping_unsafe`, `startup_hook_unsafe`,
`environment_site_root_invalid`, `installed_environment_invalid`,
`installed_source_mismatch`, `installed_editable_forbidden`,
`pytest_plugin_unsafe`, and `pytest_plugin_missing`. Messages do not print
compared paths. The guard is read-only and has no network access, installation
action, third-party dependency, product import, or environment-variable
bypass.
