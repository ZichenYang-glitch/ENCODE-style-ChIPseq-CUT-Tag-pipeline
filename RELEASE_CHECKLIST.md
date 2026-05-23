# Release Checklist

Before tagging a release, run through these checks. All automated checks below assume the `chipseq-runner` conda env is active.
The `chipseq` (full bioinformatics) env is required only for tiny real execution
(workflow_dispatch check).

## Automated Checks (chipseq-runner env)

```bash
# 1. Validation (must pass)
python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py

# 2. Dry-run smoke profiles (7 profiles, <30 s)
python3 test/test_stage8_smoke_profiles.py

# 3. Hardcoded-paths guard
python3 test/test_no_hardcoded_paths.py

# 4. BigWig stress tests (Stage 22, 6 cases)
python3 test/test_stage22_bigwig_stress.py

# 5. QC summary unit tests (Stage 24, 9 cases)
python3 test/test_stage24_qc_summary_unit.py

# 6. Manifest stress tests (Stage 25, 14 cases)
python3 test/test_stage25_manifest_stress.py

# 7. Public validation plan tests (Stage 27a, 7 cases)
python3 test/test_stage27_public_validation_plan.py

# 8. Metadata CI plan tests (Stage 27b, 7 cases)
python3 test/test_stage27b_metadata_ci_plan.py

# 9. CI workflow tests (Stage 27c, 7 cases)
python3 test/test_stage27c_ci_workflow.py

# 10. Release readiness tests (Stage 28, 11 cases)
python3 test/test_stage28_release_readiness.py

# 11. Repo hygiene
git status --short --untracked-files=all
# Expect: no results/, .snakemake/, *.fq, *.fq.gz, *.bam, *.bai, *.bw

# 12. CI status
# Check GitHub Actions for green on main / current branch.
```

Note: Default-config dry-run (`snakemake -s workflow/Snakefile --configfile config/config.yaml -n --quiet`) is expected to fail when default FASTQs don't exist. This is documented expected behavior and is not a release-blocking check.

Tiny real execution (`test_stage8b_tiny_execution.py`) requires the full `chipseq`
bioinformatics environment. It is retained as a `workflow_dispatch` CI check and
**expected to SKIP in `chipseq-runner`** when bioinformatics tools are absent.

## Manual Checks

- [ ] `README.md` install commands tested with a clean conda environment.
- [ ] `config/config.yaml` and `config/samples.tsv` agree on column names and genome labels.
- [ ] `ROADMAP_v0.2.md` scope table matches current implementation status.
- [ ] `KNOWN_ISSUES.md` status markers are up to date.
- [ ] `[v0.2.0-rc1]` section in `CHANGELOG.md` reflects all rc1 changes; `[Unreleased]` is present and empty.
- [ ] Design specs and implementation plans under `docs/superpowers/` are committed.
- [ ] All cross-links in `README.md` and `docs/*.md` resolve correctly.

## Pre-release Steps

1. Update version string in `CITATION.cff` and `CHANGELOG.md`.
2. Move `[Unreleased]` entries in `CHANGELOG.md` to a new dated section.
3. Tag the release: `git tag -a vX.Y.Z -m "Release vX.Y.Z"`.
4. Push the tag: `git push origin vX.Y.Z`.
