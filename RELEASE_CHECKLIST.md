# Release Checklist

Use this checklist before tagging a HelixWeave release. Run commands from the
repository root in the documented locked environments. Release evidence belongs
under `docs/release-checks/`; temporary outputs and data stay outside Git.

## Automated checks

### Complete deterministic Python and workflow contracts

```bash
python3 scripts/validate_samples.py --config config/config.yaml
encode-validate --config config/config.yaml
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest test -ra -p no:cacheprovider \
  -m "not platform_real_execution and not real_execution" \
  --junitxml=pytest-report.xml \
  --cov --cov-config=pyproject.toml --cov-context=test \
  --cov-fail-under=0 \
  --cov-report=term-missing \
  --cov-report=xml:coverage.xml \
  --cov-report=json:coverage.json
python3 scripts/check_junit_outcomes.py pytest-report.xml \
  --label "release deterministic tier"
python3 -m coverage report --fail-under=82
python3 test/check_snakemake_lint.py
snakefmt --check workflow/
```

The default config references example inputs, so an unrestricted default DAG
dry-run may fail when those files are absent. Use the committed smoke profiles
for deterministic no-data DAG validation; they are included in the single
pytest run above. Confirm its JUnit report contains zero skips and zero xfails.

### Real platform, scientific, and container execution

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
docker build -f containers/Dockerfile.runner -t helixweave-runner:release-check .
bash scripts/smoke_container_runner.sh docker helixweave-runner:release-check
```

Run the platform tier with Redis 7 and the locked `ci-fast` environment; run
the scientific tier in the locked `chipseq` environment. Both pytest reports
must contain zero skips and zero xfails. Do not treat a missing service,
external executable, or all-skipped selection as successful release evidence.
The container command builds locally and must not push or publish an image.

### Frontend and generated API client

```bash
npm --prefix frontend ci
npm --prefix frontend run openapi:regenerate
git diff --exit-code -- frontend/openapi.json frontend/src/api/generated/
npm --prefix frontend test -- --run
npm --prefix frontend run typecheck
npm --prefix frontend run build
npm --prefix frontend run test:e2e
```

### Repository quality

```bash
python3 -m ruff check src scripts test containers workflow/lib
python3 -m ruff format --check src scripts test containers workflow/lib
python3 -m pytest test/docs/test_internal_links.py -v
git diff --check
git status --short --untracked-files=all
```

Review untracked output rather than deleting it automatically. A release tree
must not include results, `.snakemake`, FASTQ, BAM, BigWig, database, secret, or
environment artifacts.

## Manual checks

- [ ] Current release notes and `CHANGELOG.md` describe the exact release diff.
- [ ] `README.md`, configuration, sample sheet, and quickstart examples agree.
- [ ] Scientific behavior changes have assay-specific validation evidence.
- [ ] Public-data reports remain external-data references, not bundled data.
- [ ] OpenAPI and the generated frontend client have zero drift.
- [ ] Database migrations upgrade from the supported prior schema.
- [ ] Recommended required checks (`fast-checks`, `frontend`, `browser-e2e`,
      `lint`, `lock-check`, and `coverage`) passed on the exact commit.
- [ ] Non-required `platform-real-execution`, `real-execution`, and
      `container-smoke` jobs passed on the exact commit.
- [ ] All four automated pytest entry points report zero skips and zero xfails.
- [ ] No workspace paths, environment values, or private payloads appear in
      public API responses or logs attached as release evidence.

## Release action

Tagging and publishing require explicit authorization. After approval:

1. update version and release-note metadata;
2. create the annotated tag on the verified commit;
3. push the tag;
4. publish the approved release artifacts; and
5. record the final external verification links under `docs/release-checks/`.
