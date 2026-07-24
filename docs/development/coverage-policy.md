# Python Coverage Policy

Coverage is a regression signal for executable Python behavior. It does not
replace scientific DAG contracts, dry-runs, or tiny real execution.

## Canonical measurements

Local runs and CI use the coverage configuration in `pyproject.toml` and the
pytest-native `test/` tree. The canonical tools come from the explicit
`workflow/envs/ci-fast.lock`; after creating that environment, install only the
local package without consulting an index or resolving dependencies:

```bash
python3 -m pip install --no-index --no-deps --no-build-isolation -e ".[api]"
python3 -m pip check
```

CI has two deterministic selections, and each event executes its selection
once.

Pull requests run the fast unit, contract, validator, and DAG-smoke selection:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest test -ra -p no:cacheprovider \
  -m "not full_main and not platform_real_execution and not real_execution and not bulk_rnaseq_real_execution" \
  --junitxml=pytest-report.xml \
  --cov --cov-config=pyproject.toml --cov-context=test \
  --cov-fail-under=0 \
  --cov-report=term-missing \
  --cov-report=xml:coverage.xml \
  --cov-report=json:coverage.json
```

The partial PR report is used only for the changed-lines gate. It must not be
presented as repository or core-module coverage.

Pushes to `main`, manual dispatches, nightly schedules, and published releases
run the complete deterministic suite, including tests marked `full_main`:

```bash
PYTHONDONTWRITEBYTECODE=1 python3 -m pytest test -ra -p no:cacheprovider \
  -m "not platform_real_execution and not real_execution and not bulk_rnaseq_real_execution" \
  --junitxml=pytest-report.xml \
  --cov --cov-config=pyproject.toml --cov-context=test \
  --cov-fail-under=0 \
  --cov-report=term-missing \
  --cov-report=xml:coverage.xml \
  --cov-report=json:coverage.json
```

The reports cover `src/encode_pipeline`, the Python files in `scripts`, the
workflow compatibility library in `workflow/lib`, and authored container
tooling in `containers`, with branch measurement enabled. CI retains
`.coverage`, XML, JSON, and the JUnit report. These files are run artifacts and
local scratch output; they are ignored by Git and must not be committed. A
downstream `coverage` job reads the artifact produced by the single pytest run;
it never runs pytest again.

Coverage.py's supported `patch = ["subprocess"]` mechanism measures Python
children that inherit the test environment. This includes packaged CLI and
OpenAPI export subprocesses. The platform's `ProcessRunner` intentionally
filters the environment passed to scientific processes. Coverage variables
must not be added to that allowlist: Snakemake rules and real scientific
execution remain protected by their dedicated gates.

## Measured baseline and floors

The locked environment resolves Python 3.12.13, pytest 9.1.1, pytest-cov
7.1.0, coverage.py 7.15.1, and diff-cover 9.7.2. Percentages combine statements
and branches unless a column says otherwise. The v0.3.0 release-readiness
baseline below was measured by the complete deterministic suite. The coverage
consumer starts in fail-fast shell mode so an earlier area ratchet cannot be
masked by a later report-only command.

The post-maintenance verification inventory is the authoritative count for the
current repository state:

| Gate | Current baseline |
| --- | ---: |
| All Python tests collected | 3,839 |
| PR-fast Python selection | 3,711 |
| Deterministic full-main Python selection | 3,823 |
| Platform real-execution tests | 4 |
| Scientific real-execution tests | 8 |
| Bulk RNA-seq real-execution tests | 4 |
| Frontend Vitest tests | 328 |
| Playwright browser tests | 10 |

| Area | Line | Branch | Combined | Enforced floor |
| --- | ---: | ---: | ---: | ---: |
| Repository | 86.4275% | 76.3563% | 83.7600% | 83% |
| Platform | 91.3945% | 81.5068% | 88.6800% | 88.45% |
| Services | 90.1076% | 81.1413% | 87.8251% | 87.28% |
| Persistence | 93.0396% | 73.5714% | 89.1873% | 89.06% |
| Workers | 84.7482% | 73.5632% | 82.5086% | 82.37% |
| Adapters | 86.9874% | 76.3460% | 83.8318% | report only |
| API, CLI, config, samples | 95.0440% | 89.0244% | 93.5581% | report only |
| Snakemake-facing scripts | 69.6559% | 59.4417% | 67.1216% | report only |
| Workflow compatibility library | 100.00% | n/a | 100.00% | report only |
| Container definition tooling | 97.06% | 83.33% | 95.00% | report only |

No authored Python source root is omitted. The current baseline includes
substantive fail-closed lifecycle, input-bundle, artifact/QC generation,
persistence, and worker cleanup cases added to preserve every existing
core-area floor after the two workflow product surface was completed.

The repository floor is the integer floor of the complete measured result.
Core floors remain at their previously verified values because the retirement
did not reduce those areas. Raise report-only areas only with substantive
producer, CLI, or scientific-script behavior tests. The workflow library and
container generator are authored runtime seams, so they remain in the global
denominator even though their small areas do not yet have separate floors.

On a complete-suite event, CI prints the same area reports from the one
coverage database. They can also be reproduced locally with normal coverage.py
filters, for example:

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
branch, use the preceding branch as the review base. The changed-lines floor
does not replace the complete global and core floors; those are enforced on
the next full-main, dispatch, nightly, or release run.

## Ratchet rules

- Never lower a global, core, or changed-lines floor merely to pass CI. The
  enforced floors are recorded in the measured-baseline table above; changed
  executable Python lines remain subject to the 80% gate.
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
