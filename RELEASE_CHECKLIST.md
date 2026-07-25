# Release Checklist

Use this checklist before tagging a HelixWeave release. Run commands from a
clean exact commit in the documented locked environments. Release evidence
belongs under `docs/release-checks/`; temporary outputs and data stay outside
Git.

Tagging, publishing, uploading assets, or creating a GitHub Release always
requires separate explicit authorization. Passing this checklist only makes a
commit release-ready.

## v0.3.0 release identity

The v0.3.0 release contract is:

| Surface | Fixed identity |
| :--- | :--- |
| Product name | HelixWeave |
| Git tag | `v0.3.0` |
| Python distribution | `helixweave` |
| Python package | `encode_pipeline` |
| Package/API/frontend version | `0.3.0` |
| Primary CLI | `helixweave` |
| Compatibility CLIs | `encode-validate`, `encode-manifest`, `encode-dag`, `encode-worker` |
| Bundled workflows | ENCODE-style ChIP-seq/CUT&Tag/ATAC/MNase and `bulk-rnaseq` |
| Bulk upstream | nf-core/rnaseq 3.26.0 |

HelixWeave is the product and distribution name. The `encode_pipeline` import
namespace, `encode-*` CLI names, workflow IDs, and repository slug remain
compatibility identities in this release. The reviewed name and local-candidate
transition evidence is recorded in the
[Python identity migration](docs/release-checks/v0.3.0-python-identity-migration.md).

The approved user-uploaded GitHub Release asset set is exactly:

- `helixweave-0.3.0-py3-none-any.whl`;
- `helixweave-0.3.0.tar.gz`; and
- `SHA256SUMS`, covering those two files.

The wheel and sdist contain the Python APIs, primary and compatibility console
entry points, Alembic migrations, artifact catalog, and versioned adapter
contracts; they are not the complete browser/scientific product tree. A
qualification-only source-trial archive must also be created from the exact
reviewed Git tree with a `helixweave-v0.3.0/` prefix. It contains the committed
frontend, workflow, scripts, documentation, and lock files, but no environment,
reference, container, FASTQ, result, cache, or secret payload. Do not upload
that source-trial archive to the GitHub Release.

v0.3.0 does not publish a Docker/OCI or Apptainer image. The existing
runner-only definitions remain local ENCODE compatibility assets. Container
publication, registry naming, OCI labels, and release-tag mapping require a
later explicit decision.

Before any future publication, record the exact commit and tree:

```bash
git status --short --untracked-files=all
git rev-parse HEAD
git show -s --format=%T HEAD
git describe --tags --exact-match HEAD
```

The first command must be empty. Before the tag exists, replace the last command
with an explicit assertion that no release tag has been created.

## Build and clean-install checks

Build only from the verified clean tree, using the locked Python environment:

```bash
python3 -m build --no-isolation --sdist --wheel --outdir dist
python3 -m pytest test/packaging -v
```

`--no-isolation` is intentional: build backend and frontend dependencies come
from the reviewed `ci-fast.lock` instead of an unpinned network bootstrap.

Inspect `dist/` and reject any unexpected filename. Create the formal checksums
there. Build the no-checkout qualification archive into a separate directory
from the same exact commit:

```bash
(
  cd dist
  sha256sum \
    helixweave-0.3.0-py3-none-any.whl \
    helixweave-0.3.0.tar.gz \
    > SHA256SUMS
  sha256sum --check SHA256SUMS
)
mkdir -p qualification
git archive --format=tar --prefix=helixweave-v0.3.0/ HEAD \
  | gzip -n > qualification/helixweave-v0.3.0-source-trial.tar.gz
sha256sum qualification/helixweave-v0.3.0-source-trial.tar.gz
```

Use separate new virtual environments to install the wheel and the sdist. Do
not let either environment import the repository checkout through `PYTHONPATH`,
the current directory, user site packages, or an editable install.

For each installed artifact:

- run `python -m pip check`;
- confirm `importlib.metadata.version("helixweave") == "0.3.0"` and
  `encode-pipeline` metadata is absent;
- exercise `helixweave`, `python -m encode_pipeline`, all four compatibility
  console scripts, and their module entry points;
- confirm `packages_distributions()["encode_pipeline"] == ["helixweave"]`;
- create the default registry and confirm the two workflow IDs;
- export OpenAPI and compare it with `frontend/openapi.json`;
- create and migrate a fresh file-backed SQLite database;
- confirm all Alembic revisions and adapter contract package data are readable;
  and
- run from a directory outside the unpacked source and build trees.

Verify the separately recorded source-trial checksum, extract it into a new
directory, and complete the maintained
[local trial checklist](docs/local-trial-checklist.md). The trial may use
package and npm caches, but it must not depend on a Git checkout, user-home
source files, hard-coded workspace paths, or unrecorded runtime data.

## Automated checks

### Complete deterministic Python and workflow contracts

```bash
python3 scripts/validate_samples.py --config config/config.yaml
encode-validate --config config/config.yaml
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest test -ra -p no:cacheprovider \
  -m "not platform_real_execution and not real_execution and not bulk_rnaseq_real_execution" \
  --junitxml=pytest-report.xml \
  --cov --cov-config=pyproject.toml --cov-context=test \
  --cov-fail-under=0 \
  --cov-report=term-missing \
  --cov-report=xml:coverage.xml \
  --cov-report=json:coverage.json
python3 scripts/check_junit_outcomes.py pytest-report.xml \
  --label "release deterministic tier"
python3 -m coverage report --fail-under=83
python3 test/check_snakemake_lint.py
snakefmt --check workflow/
```

The default config references example inputs, so an unrestricted default DAG
dry-run may fail when those files are absent. Use the committed smoke profiles
for deterministic no-data DAG validation; they are included in the single
pytest run above. Confirm its JUnit report contains zero skips and zero xfails.

### Five pytest tiers and real execution

The five maintained pytest entry points are PR fast, full main, platform real,
scientific real, and Bulk RNA-seq real. Every executed tier must report zero
skips and zero xfails; an empty or all-skipped selection is a failure.

Run platform and ENCODE scientific execution when required by the release diff:

```bash
ENCODE_PIPELINE_TEST_REDIS_URL=redis://127.0.0.1:6379/15 \
  python3 -m pytest -m platform_real_execution test/workers -ra -v \
  --junitxml=platform-real-report.xml
python3 scripts/check_junit_outcomes.py platform-real-report.xml \
  --label "release platform-real tier"
HELIXWEAVE_REQUIRE_REAL_EXECUTION=1 \
  python3 -m pytest -m real_execution test/real_execution -ra -v \
  --junitxml=scientific-real-report.xml
python3 scripts/check_junit_outcomes.py scientific-real-report.xml \
  --label "release scientific-real tier"
```

Run the platform tier with Redis 7 and the locked `ci-fast` environment; run
the scientific tier in the locked `chipseq` environment. Do not treat a
missing service or executable as successful release evidence.

Bulk RNA-seq protected evidence can be inherited only when the exact committed
execution implementation closure is unchanged. Compare the candidate against
the reviewed evidence base using
`execution-implementation-manifest-1.0.0.json`; any changed manifest-listed
file, migration revision, runtime/results contract, managed-container boundary,
or generated execution-manifest identity requires a new exact-HEAD protected
Gate on the approved ephemeral runner. Documentation, release metadata, and
frontend-only changes do not by themselves require that runner when the
execution closure is byte-identical. Record the inherited or rerun workflow,
job, artifact, commit, tree, and runtime identity in the release evidence.

The runner-only container smoke is needed only when its Docker/Apptainer
definition, `runner.lock`, generated runner definition, or smoke tooling
changed:

```bash
docker build -f containers/Dockerfile.runner -t helixweave-runner:release-check .
bash scripts/smoke_container_runner.sh docker helixweave-runner:release-check
```

This image is local verification only. Do not push or publish it.

### Frontend and generated API client

```bash
npm --prefix frontend ci
npm --prefix frontend run openapi:regenerate
git diff --exit-code -- frontend/openapi.json frontend/src/api/generated/
npm --prefix frontend test -- --run
npm --prefix frontend run typecheck
npm --prefix frontend run build
export PLAYWRIGHT_BROWSERS_PATH="$PWD/.local/playwright-browsers"
npm --prefix frontend exec -- playwright install chromium
npm --prefix frontend run test:e2e
```

Browser evidence must cover both runtime-available and runtime-unavailable
product states and the 1440×900, 1024×768, 390×844, and 360×800 viewports.
The task-owned `PLAYWRIGHT_BROWSERS_PATH` is required so this gate cannot pass
by consuming an unrecorded user-home browser cache. On a clean Linux host,
install Chromium's system libraries first with the reviewed Playwright
`install-deps chromium` operation (it may require administrator approval);
that host operation does not add a browser or package to the release assets.

### Database upgrade

The v0.3.0 supported pre-generation schema fixture is Alembic revision
`20260714_07`; the release head is `20260717_08`. The upgrade test must preserve
workflow identity, run lifecycle, events, logs, queue assignment, build
identity, validated snapshot, artifacts, and QC. Legacy artifacts/QC remain
explicitly generation-unbound after the migration; the release must not invent
revision evidence.

```bash
python3 -m pytest \
  test/persistence/test_migrations.py::test_v030_supported_prior_schema_upgrade_preserves_complete_product_record \
  -v
```

The test must also pass against the installed wheel. See the
[local runtime upgrade procedure](docs/development/local-platform-runtime.md#v030-database-upgrade).
There is no online downgrade contract.

### Repository quality

```bash
python3 -m ruff check src scripts test containers workflow/lib
python3 -m ruff format --check src scripts test containers workflow/lib
python3 -m pytest test/docs/test_internal_links.py -v
git diff --check
git status --short --untracked-files=all
```

Review untracked output rather than deleting it automatically. A release tree
must not include results, `.snakemake`, FASTQ, BAM, BigWig, database, secret,
environment, JDK, OCI, reference, index, or cache artifacts.

## Manual checks

- [ ] `CHANGELOG.md` and the draft release notes describe the exact release diff.
- [ ] README and maintained docs consistently describe HelixWeave, both
      workflows, the `helixweave` distribution/CLI, and retained
      import/compatibility CLI identities.
- [ ] The three user-uploaded Release assets have the expected names and
      verified checksums; the qualification source trial and container images
      are not in the v0.3.0 upload set.
- [ ] Wheel, sdist, and source-trial archive checks ran outside the source tree.
- [ ] OpenAPI and the generated frontend client have zero drift.
- [ ] The supported `20260714_07` database fixture upgrades to `20260717_08`
      without lifecycle or result-metadata corruption.
- [ ] PR-fast, full-main, platform-real, scientific-real, and applicable Bulk
      real pytest reports contain zero skips and zero xfails.
- [ ] Recommended required checks (`fast-checks`, `frontend`, `browser-e2e`,
      `lint`, `lock-check`, and `coverage`) passed on the exact candidate.
- [ ] Bulk protected evidence was either rerun on the exact candidate or
      inherited with a byte-identical execution-closure proof.
- [ ] No workspace paths, environment values, credentials, or private payloads
      appear in public API responses, doctor output, screenshots, or evidence.
- [ ] Synthetic fixtures are described only as execution/contract evidence,
      not biological validation or production-scale performance evidence.
- [ ] No authentication, multi-user, HPC/cloud, object-storage, support-SLA, or
      container-publication capability is implied.

## Release action

After explicit publication authorization only:

1. replace the `[Unreleased]` candidate heading with
   `[v0.3.0] - YYYY-MM-DD`, add the matching release date to citation metadata,
   and commit that release-only metadata change;
2. run the complete required exact-HEAD CI and applicable real-execution
   disposition on the new commit; no asset from the earlier Draft PR head is a
   final release asset;
3. rebuild the wheel, sdist, source-trial archive, and basename-only
   `SHA256SUMS` from that clean exact commit, then repeat every clean-install,
   migration, OpenAPI, browser, and archive check above;
4. re-confirm the final commit, tree, version metadata, clean status, CI,
   protected-evidence inheritance or rerun, and asset checksums;
5. create and push the annotated `v0.3.0` tag on that exact commit;
6. create the GitHub Release and upload only the approved asset set;
7. verify the downloaded flat release assets against `SHA256SUMS`; and
8. record the final tag, commit, tree, CI, protected evidence, and asset digests
   under `docs/release-checks/`.

Do not rename the repository or Python identities, publish to PyPI/conda/a
container registry, or create additional release assets without separate
authorization.
