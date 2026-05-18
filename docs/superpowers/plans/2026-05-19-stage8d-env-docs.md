# Stage 8d: Environment Cleanup and Execution Docs — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Remove `defaults` channel from chipseq.yml, add local execution docs to README, close Stage 8d in KNOWN_ISSUES.md.

**Architecture:** Three file edits, no new files. Remove `defaults` from chipseq.yml channels, add a "Local execution" section to README Developer Notes, mark Stage 8d complete in KNOWN_ISSUES.md.

**Tech Stack:** YAML editing, Markdown documentation

---

### Task 1: Remove `defaults` channel from chipseq.yml

**Files:**
- Modify: `workflow/envs/chipseq.yml`

- [ ] **Step 1: Remove the `defaults` channel line**

The current channels block:

```yaml
channels:
  - conda-forge
  - bioconda
  - defaults
```

Change to:

```yaml
channels:
  - conda-forge
  - bioconda
```

Remove only the `- defaults` line. All dependencies stay unchanged.

- [ ] **Step 2: Verify YAML syntax**

```bash
python3 -c "
import yaml
with open('workflow/envs/chipseq.yml') as f:
    data = yaml.safe_load(f)
assert data['name'] == 'chipseq'
assert data['channels'] == ['conda-forge', 'bioconda']
assert 'defaults' not in data['channels']
print('chipseq.yml: OK (defaults removed, conda-forge + bioconda only)')
"
```

- [ ] **Step 3: Optionally verify environment solve**

If micromamba is available and solve time is acceptable:

```bash
micromamba create -f workflow/envs/chipseq.yml --dry-run 2>&1 | tail -5
```

If solve is too slow, skip — full validation occurs via `workflow_dispatch` CI.

---

### Task 2: Add "Local execution" section to README.md

**Files:**
- Modify: `README.md`

- [ ] **Step 1: Add Local execution section under Developer Notes**

After the CI section, before `## License`, insert:

```markdown
### Local execution

**Prerequisites:**

# Create environment (the env file already names it "chipseq")
micromamba create -f workflow/envs/chipseq.yml
micromamba activate chipseq

**Validation and dry-run (fast, no data needed):**

python3 scripts/validate_samples.py --config config/config.yaml
python3 test/test_validation_stress.py
python3 test/test_stage8_smoke_profiles.py

**Tiny real execution (requires full chipseq env):**

python3 test/test_stage8b_tiny_execution.py

**Full workflow run:**

# 1. Edit config/config.yaml and config/samples.tsv for your data.
# 2. Dry-run:
snakemake -s workflow/Snakefile --configfile config/config.yaml -n

# 3. Execute (replace N with core count):
snakemake -s workflow/Snakefile --configfile config/config.yaml --cores N

# On a workstation, use as many cores as available.  On a shared server,
# limit cores and consider --latency-wait for NFS filesystems.
```

- [ ] **Step 2: Verify the new section is present**

```bash
python3 -c "
with open('README.md') as f:
    c = f.read()
assert '### Local execution' in c
assert 'micromamba create -f workflow/envs/chipseq.yml' in c
assert 'micromamba activate chipseq' in c
print('README Local execution section: OK')
"
```

---

### Task 3: Update KNOWN_ISSUES.md — Stage 8d complete

**Files:**
- Modify: `KNOWN_ISSUES.md`

- [ ] **Step 1: Replace the three `⬜` Stage 8d items with completion block**

Current lines (~214-217):

```markdown
- ⬜ Remove local `prefix` metadata from exported Conda environment files. (Stage 8d)
- ⬜ Decide whether to split `workflow/envs/chipseq.yml` into smaller
  responsibility-specific environments. (Stage 8d)
- ⬜ Document recommended execution on workstation/server environments. (Stage 8d)
```

Replace with:

```markdown
**Stage 8d completed 2026-05-19** — environment cleanup and execution
documentation.

- ✅ Remove local `prefix` metadata from exported Conda environment files.
  (Already absent from chipseq.yml and ci-fast.yml; confirmed.)
- ✅ Remove `defaults` channel from `workflow/envs/chipseq.yml`.
  Channels are now `conda-forge` + `bioconda` only.
- ✅ Decide whether to split environments further: no further split in
  Stage 8d. `chipseq.yml` remains the full runtime environment;
  `ci-fast.yml` remains the lightweight CI environment. Further
  task-specific envs can be revisited later only if solve/runtime
  pain persists.
- ✅ Document local execution, validation, smoke tests, and full
  workflow run in `README.md`.
```

- [ ] **Step 2: Verify KNOWN_ISSUES.md**

```bash
python3 -c "
with open('KNOWN_ISSUES.md') as f:
    c = f.read()
assert 'Stage 8d completed 2026-05-19' in c
assert 'no further split in Stage 8d' in c
assert 'conda-forge + bioconda only' in c
print('KNOWN_ISSUES: Stage 8d marked complete')
"
```

---

### Task 4: Full verification

- [ ] **Step 1: Run validation stress tests**

```bash
python3 test/test_validation_stress.py
```
Expected: 15/15 PASS.

- [ ] **Step 2: Run Stage 8a smoke profiles**

```bash
python3 test/test_stage8_smoke_profiles.py
```
Expected: 7/7 PASS.

- [ ] **Step 3: Run Stage 8b tiny real execution**

```bash
python3 test/test_stage8b_tiny_execution.py
```
Expected: PASS with exit code 0.

- [ ] **Step 4: Verify git status**

```bash
git status --short --untracked-files=all
```
Expected: only intentional files. No generated artifacts.

---

### Task 5: Commit

- [ ] **Step 1: Stage**

```bash
cd /home/irenadler/workflow/chipseq && \
  git add workflow/envs/chipseq.yml \
         README.md \
         KNOWN_ISSUES.md \
         docs/superpowers/specs/2026-05-19-stage8d-env-docs-design.md
```

- [ ] **Step 2: Commit**

```bash
git commit -m "$(cat <<'EOF'
Add Stage 8d environment cleanup and execution documentation

Remove defaults channel from chipseq.yml (conda-forge + bioconda only).
Add local execution section to README. Close Stage 8d in KNOWN_ISSUES.md.

Co-Authored-By: deepseek-v4-pro[1m] <deepseek-ai@claude-code-best.win>
EOF
)"
```

---

### Self-Review

- ✅ Spec coverage: chipseq.yml cleanup (Task 1), README docs (Task 2), KNOWN_ISSUES (Task 3), verification (Task 4)
- ✅ No placeholders, no TBDs
- ✅ No new pipeline features — docs and one channel removal only
- ✅ No CI badge added
- ✅ README uses `micromamba create -f` (not `-n chipseq -f`)
- ✅ Env splitting wording: "no further split in Stage 8d"
