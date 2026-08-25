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

CI has two deterministic selections. Each event executes its selection as one
logical suite in two mutually exclusive and exhaustive physical shards:
`scripts/split_deterministic_shards.py` assigns every collected `test_*.py`
file by `SHA-256(relative POSIX path) % 2`, each shard runs only its own
files and records an unrendered `.coverage.shard-N`, and the `fast-checks`
aggregate combines both data files into one canonical measurement. Local runs
may still execute either selection in a single process.

Pull requests run the fast unit, contract, validator, and DAG-smoke selection
(per shard, with the shard's file list in place of `test`):

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest test -ra \
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
python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest test -ra \
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
tooling in `containers`, with branch measurement enabled. CI retains the
combined `.coverage`, XML, JSON, and both shard JUnit reports. These files are
run artifacts and local scratch output; they are ignored by Git and must not be
committed. A downstream `coverage` job reads the canonical artifact combined
from the two shard data files; it never runs pytest again.

Coverage.py's supported `patch = ["subprocess"]` mechanism measures Python
children that inherit the test environment. This includes packaged CLI and
OpenAPI export subprocesses. The platform's `ProcessRunner` intentionally
filters the environment passed to scientific processes. Coverage variables
must not be added to that allowlist: Snakemake rules and real scientific
execution remain protected by their dedicated gates.

Alembic runs only the revision bytes that migration admission first copies to
a private snapshot. For the migration integration module, CI creates one
empty, mode-`0700` coverage root and passes that exact directory as an
additional pytest-cov source. A test-only fixture routes only admitted
`helixweave-migration-snapshot-*` directories beneath that root, and
Coverage.py's `[paths]` mapping attributes their execution back to the
canonical `src/encode_pipeline/persistence/alembic` sources. The fixture
requires the root to be absolute, resolved, private, empty, and pre-existing;
each snapshot and then the root are removed. The system temporary directory
is never a coverage source. This in-process attribution is independent of the
subprocess patch and does not omit revision assets from the global or
persistence ratchets.

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
| All Python tests collected | 4,454 |
| PR-fast Python selection | 4,313 |
| Deterministic full-main Python selection | 4,438 |
| Platform real-execution tests | 4 |
| Scientific real-execution tests | 8 |
| Bulk RNA-seq real-execution tests | 4 |
| Frontend Vitest tests | 328 |
| Playwright browser tests | 10 |

| Area | Line | Branch | Combined | Enforced floor |
| --- | ---: | ---: | ---: | ---: |
| Repository | 87.3198% | 76.8626% | 84.5724% | 83% |
| Platform | 91.5020% | 81.1340% | 88.7490% | 88.45% |
| Services | 89.9072% | 79.8655% | 87.3934% | 87.28% |
| Persistence | 93.7230% | 80.3488% | 90.7517% | 89.06% |
| Workers | 88.2431% | 76.6667% | 86.0192% | 82.37% |
| Adapters | 86.8743% | 75.8610% | 83.5886% | report only |
| API, CLI, config, samples | 92.7553% | 85.6464% | 91.0633% | report only |
| Snakemake-facing scripts | 71.8978% | 64.2955% | 69.9031% | report only |
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
