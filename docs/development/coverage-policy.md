# Python Coverage Policy

Coverage is a regression signal for executable Python behavior. It does not
replace scientific DAG contracts, dry-runs, or tiny real execution.

## Canonical measurement

Local runs and CI use the coverage configuration in `pyproject.toml` and the
complete pytest-native `test/` tree. The canonical tools come from the explicit
`workflow/envs/ci-fast.lock`; after creating that environment, install only the
local package without consulting an index or resolving dependencies:

```bash
python3 -m pip install --no-index --no-deps --no-build-isolation -e ".[api]"
python3 -m pip check
```

Then produce the reports with:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest test -ra -p no:cacheprovider \
  --cov --cov-config=pyproject.toml --cov-context=test \
  --cov-report=term-missing \
  --cov-report=xml:coverage.xml \
  --cov-report=json:coverage.json
```

The reports cover `src/encode_pipeline`, the Python files in `scripts`, the
workflow compatibility library in `workflow/lib`, and authored container
tooling in `containers`, with branch measurement enabled. XML and JSON files
are CI artifacts and local scratch output; they are ignored by Git and must not
be committed.

Coverage.py's supported `patch = ["subprocess"]` mechanism measures Python
children that inherit the test environment. This includes packaged CLI and
OpenAPI export subprocesses. The platform's `ProcessRunner` intentionally
filters the environment passed to scientific processes. Coverage variables
must not be added to that allowlist: Snakemake rules and real scientific
execution remain protected by their dedicated gates.

## Measured baseline

The locked environment resolves Python 3.12.13, pytest 9.1.1, pytest-cov
7.1.0, coverage.py 7.15.1, and diff-cover 9.7.2. Percentages combine statements
and branches unless a column says otherwise.

| Area | Line | Branch | Combined | Enforced floor |
| --- | ---: | ---: | ---: | ---: |
| Repository | 80.97% | 71.80% | 78.72% | 78% |
| Platform | 91.29% | 80.28% | 88.45% | 88.45% |
| Services | 89.57% | 80.39% | 87.28% | 87.28% |
| Persistence | 92.96% | 72.69% | 89.07% | 89.06% |
| Workers | 84.78% | 71.83% | 82.38% | 82.37% |
| Adapters | 92.01% | 83.78% | 89.93% | report only |
| API, CLI, config, samples | 94.67% | 88.44% | 93.11% | report only |
| Snakemake-facing scripts | 32.46% | 26.48% | 31.04% | report only |
| Workflow compatibility library | 100.00% | n/a | 100.00% | report only |
| Container definition tooling | 97.06% | 83.33% | 95.00% | report only |

The repository floor is the integer floor of the complete measured result.
Core floors are the measured combined values truncated to two decimal places,
so they prevent regression without claiming a higher result than was observed.
The low scripts result is visible rather than omitted; raise it only with
substantive producer, CLI, or scientific-script behavior tests. The workflow
library and container generator are authored runtime seams, so they remain in
the global denominator even though their small areas do not yet have separate
floors.

CI prints the same area reports from the single combined coverage database.
They can also be reproduced locally with normal coverage.py filters, for
example:

```bash
python3 -m coverage report \
  --include='src/encode_pipeline/platform/*' --fail-under=88.45
python3 -m coverage report --include='scripts/*' --fail-under=0
```

The Snakemake DSL in `workflow/Snakefile` and `workflow/rules/*.smk` is not an
importable Python source tree and is therefore outside coverage.py's
denominator. Config/sample validation, DAG snapshots and dry-runs, Snakemake
lint, and tiny deterministic execution provide the equivalent scientific
gates. Dockerfiles and generated Apptainer definitions are likewise not Python
source; their authored generator is measured, generated definitions have a
drift test, and the real container-smoke tier exercises the runtime boundary.

## Changed lines

Changed executable Python lines must have at least 80% coverage:

```bash
diff-cover coverage.xml --compare-branch=<review-base> --fail-under=80
```

Pull requests compare against GitHub's exact base SHA. For a stacked local
branch, use the preceding branch as the review base.

## Ratchet rules

- Never lower a global, core, or changed-lines floor merely to pass CI.
- When substantive tests improve a measured area, raise its floor by 1–3
  percentage points without exceeding the verified result.
- The original medium-term repository target was 60–70%; this complete
  baseline already exceeds it for both line and branch coverage, so that range
  is now a minimum milestone rather than a reason to lower a floor. Core
  runtime paths target 75–85% or better for both dimensions.
- Prefer lifecycle, persistence, worker, adapter, API, CLI, artifact/QC, and
  scientific contract behavior over tests for trivial getters.
- A core regression found during test maintenance must be fixed before the
  replacement is accepted.

Generated sources may be excluded when they appear in the Python source tree.
Pure typing branches and genuinely unreachable platform fallbacks may use the
standard coverage exclusions when the reason is clear in code review. Broad
`omit` patterns, opportunistic `pragma: no cover`, and deleting low-coverage
code or tests to improve the percentage are not acceptable.
