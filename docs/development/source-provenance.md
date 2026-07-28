# Python Source Provenance

HelixWeave development checks fail closed before importing `encode_pipeline`.
This prevents an editable install from another worktree from silently supplying
the product code or its distribution metadata.

## Controlled modes

Checkout mode requires an explicit repository root:

```bash
python3 scripts/source_provenance.py checkout --repository-root .
```

The guard resolves the repository, its canonical `src` directory, the package
origin and its only search location. Every distribution that claims the
`encode_pipeline` namespace, plus metadata discovered under the current or
legacy distribution name, is audited. The owner must be `helixweave`; source
metadata must be inside that checkout or editable metadata must name that exact
checkout. A sibling worktree is rejected even when its contents or commit are
identical.

Installed-artifact mode is only for wheel and sdist clean-room tests:

```bash
/path/to/isolated/bin/python scripts/source_provenance.py installed-artifact
```

It requires a real isolated environment. The package and all claiming
distribution metadata must live in that environment's `site-packages`;
editable metadata and checkout-sourced imports are rejected. Normal installed
HelixWeave console and module entry points do not require a Git checkout and do
not invoke this development guard.

## Execution order and failures

The checkout guard runs from `test/conftest.py` before pytest collection, from
the OpenAPI exporter before the FastAPI app import, and from source-owned local
CLI wrappers before product imports. CI also invokes it immediately after the
editable install and before pytest or config validation.

Failures exit with status 2 and one public-safe line:

```text
source provenance check failed [reason_code]: actionable guidance
```

Stable reason codes are `repository_root_invalid`, `source_root_invalid`,
`module_not_found`, `module_origin_mismatch`,
`module_search_location_mismatch`, `distribution_missing`,
`distribution_metadata_invalid`, `distribution_identity_mismatch`,
`distribution_source_mismatch`, `installed_environment_invalid`,
`installed_source_mismatch`, and `installed_editable_forbidden`. Messages do
not print the compared paths. The guard has no environment-variable bypass,
network access, installation action, third-party dependency, or product
module import.
