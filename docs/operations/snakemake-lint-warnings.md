# Accepted `snakemake --lint` warnings

This project uses `snakemake --lint` as a quality gate, but the following
warning categories are accepted as of the Phase 2a baseline. New or changed
warnings cause CI to fail.

The authoritative baseline is
[`snakemake-lint-warnings.txt`](./snakemake-lint-warnings.txt). Update it with:

```bash
python3 test/check_snakemake_lint.py --update-baseline
```

## Accepted categories

1. **Mixed rules and functions in same snakefile**
   - Affects: `Snakefile`, `common.smk`, `peaks.smk`, `mnase.smk`,
     `replicates.smk`, `idr.smk`, `idr_reproducibility.smk`,
     `consensus.smk`, `report.smk`
   - Rationale: helper functions are co-located with the rules that use them
     and are intentionally kept in their module files rather than moved to a
     shared `common.smk`.

2. **No log directive defined**
   - Affects: lightweight rules in `common.smk` and `mnase.smk`.
   - Rationale: these rules emit minimal output and rely on `set -e -o pipefail`
     in their shell blocks; adding log files is deferred to future
     observability work.

3. **Absolute path `/tmp` in `qc.smk`**
   - Affects: temporary directories used by internal QC helper rules.
   - Rationale: `/tmp` is used for intermediate scratch files that are not
     workflow outputs; making them configurable is part of the HPC profile
     work in Phase 5.

## CI behavior

The lint check runs `test/check_snakemake_lint.py`. It normalizes absolute
paths out of `snakemake --lint` output and compares it to the baseline. Any
new warning produces a diff and fails the build.
