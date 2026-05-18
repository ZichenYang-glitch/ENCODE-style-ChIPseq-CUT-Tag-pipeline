# Stage 8c: GitHub Actions CI Wiring — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Add a GitHub Actions CI workflow with fast validation/dry-run checks on PRs and push, plus a manual real-execution job via workflow_dispatch.

**Architecture:** One workflow file (`.github/workflows/ci.yml`) with two jobs. `fast-checks` runs on pull_request, push to main/stage*, and workflow_dispatch — validating config, stress tests, and Stage 8a dry-run smoke via micromamba. `real-execution` is gated behind workflow_dispatch only — runs Stage 8b tiny real execution. Both share the same micromamba cache key.

**Tech Stack:** GitHub Actions, mamba-org/setup-micromamba@v3, Python 3, Snakemake

---

### Task 1: Create the CI workflow file

**Files:**
- Create: `.github/workflows/ci.yml`

- [ ] **Step 1: Create `.github/workflows/` directory**

```bash
mkdir -p .github/workflows
```

- [ ] **Step 2: Write the complete CI workflow**

```yaml
name: CI

"on":
  pull_request:
  push:
    branches: [main, "stage*"]
  workflow_dispatch:

permissions:
  contents: read

concurrency:
  group: ${{ github.workflow }}-${{ github.ref }}
  cancel-in-progress: true

jobs:
  fast-checks:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v4

      - name: Setup micromamba
        uses: mamba-org/setup-micromamba@v3
        with:
          environment-file: workflow/envs/chipseq.yml
          environment-name: chipseq
          cache-environment: true
          cache-environment-key: chipseq-${{ runner.os }}-${{ hashFiles('workflow/envs/chipseq.yml') }}

      - name: Validate default config
        shell: micromamba-shell {0}
        run: python3 scripts/validate_samples.py --config config/config.yaml

      - name: Validation stress tests
        shell: micromamba-shell {0}
        run: python3 test/test_validation_stress.py

      - name: Stage 8a dry-run smoke profiles
        shell: micromamba-shell {0}
        run: python3 test/test_stage8_smoke_profiles.py

  real-execution:
    runs-on: ubuntu-latest
    if: github.event_name == 'workflow_dispatch'
    steps:
      - uses: actions/checkout@v4

      - name: Setup micromamba
        uses: mamba-org/setup-micromamba@v3
        with:
          environment-file: workflow/envs/chipseq.yml
          environment-name: chipseq
          cache-environment: true
          cache-environment-key: chipseq-${{ runner.os }}-${{ hashFiles('workflow/envs/chipseq.yml') }}

      - name: Stage 8b tiny real execution
        shell: micromamba-shell {0}
        run: python3 test/test_stage8b_tiny_execution.py
```

- [ ] **Step 3: Verify YAML syntax**

```bash
python3 -c "
import yaml
with open('.github/workflows/ci.yml') as f:
    yaml.safe_load(f)
print('YAML syntax OK')
"
```

If PyYAML is not available (it is in the chipseq env), use:

```bash
python3 -c "
with open('.github/workflows/ci.yml') as f:
    content = f.read()
assert 'name: CI' in content
assert 'fast-checks' in content
assert 'real-execution' in content
assert 'setup-micromamba@v3' in content
print('Workflow file OK')
"
```

---

### Task 2: Update KNOWN_ISSUES.md

**Files:**
- Modify: `KNOWN_ISSUES.md`

- [ ] **Step 1: Mark Stage 8c complete, keep 8d pending**

Find the Stage 8 section at line ~210. After the Stage 8b items, add:

```markdown
**Stage 8c completed 2026-05-18** — GitHub Actions CI wiring.

- ✅ GitHub Actions CI workflow (`.github/workflows/ci.yml`) with fast
  validation + dry-run checks on PR/push and manual real-execution via
  workflow_dispatch. (Stage 8c)
```

Convert the existing 8c/8d items:

```markdown
- ⬜ GitHub Actions / CI wiring. (Stage 8c)
```

to:

```markdown
- ✅ GitHub Actions / CI wiring. (Stage 8c, completed 2026-05-18)
```

- [ ] **Step 2: Verify**

```bash
python3 -c "
with open('KNOWN_ISSUES.md') as f:
    content = f.read()
assert 'Stage 8c completed 2026-05-18' in content
assert 'GitHub Actions CI wiring' in content
print('OK')
"
```

---

### Task 3: Update README.md

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add CI section under Developer Notes**

After the repository layout section, add:

```markdown
### CI

A GitHub Actions workflow runs on every PR and push to `main` / `stage*`:

```bash
# These commands also run locally (see Smoke-test profiles above):
python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py
python3 test/test_stage8_smoke_profiles.py
```

A manual `workflow_dispatch` job runs the tiny real-execution harness.
See `.github/workflows/ci.yml`.
```

The CI badge can be added after the first successful run:

```markdown
[![CI](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/workflows/ci.yml)
```

(Add the badge below the existing badges at the top of README.md.)

- [ ] **Step 2: Add CI badge to README.md header**

At the top of README.md, after the existing badges line:

```markdown
[![CI](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/workflows/ci.yml/badge.svg)](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/workflows/ci.yml)
```

- [ ] **Step 3: Verify**

```bash
python3 -c "
with open('README.md') as f:
    content = f.read()
assert 'workflows/ci.yml' in content
print('OK')
"
```

---

### Task 4: Commit

- [ ] **Step 1: Stage all files**

```bash
cd /home/irenadler/workflow/chipseq && \
  git add .github/workflows/ci.yml \
         KNOWN_ISSUES.md \
         README.md \
         docs/superpowers/specs/2026-05-18-stage8c-github-actions-ci-design.md
```

- [ ] **Step 2: Commit**

```bash
git commit -m "$(cat <<'EOF'
Add Stage 8c GitHub Actions CI workflow

Fast validation + dry-run checks on PR/push via micromamba cache.
Real execution gated behind workflow_dispatch.

Co-Authored-By: deepseek-v4-pro[1m] <deepseek-ai@claude-code-best.win>
EOF
)"
```

---

### Self-Review

- ✅ Spec coverage: CI workflow (Task 1), KNOWN_ISSUES (Task 2), README (Task 3), commit (Task 4)
- ✅ No placeholders, no TBDs
- ✅ All file paths absolute
- ✅ No Stage 8d env refactor creep
