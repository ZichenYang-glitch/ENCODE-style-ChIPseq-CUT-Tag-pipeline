# Stage 8d: Environment Cleanup and Execution Documentation — Design Spec

**Date:** 2026-05-19
**Status:** design review
**Scope:** remove defaults channel, document execution, update README + KNOWN_ISSUES
**Excluded:** env splitting, new pipeline features, CI badge

---

## 1. Goals

1. Clean up `workflow/envs/chipseq.yml`: remove `defaults` channel, keep
   `conda-forge` + `bioconda` only.
2. Document local execution workflows in README.md: install, validation,
   dry-run smoke, tiny real execution, full workflow run.
3. Confirm no further environment split in Stage 8d.
4. Update KNOWN_ISSUES.md to close Stage 8d.
5. No new pipeline features, no CI badge (deferred to first green main run).

---

## 2. Environment changes

### 2.1 `chipseq.yml` — remove `defaults` channel

**Before:**

```yaml
name: chipseq
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - fastqc
  - trim-galore
  - bowtie2
  - samtools
  - picard
  - macs3
  - deeptools
  - bedtools
  - snakemake >=8.0
  - multiqc
  - idr
  - seacr
```

**After:**

```yaml
name: chipseq
channels:
  - conda-forge
  - bioconda
dependencies:
  - fastqc
  - trim-galore
  - bowtie2
  - samtools
  - picard
  - macs3
  - deeptools
  - bedtools
  - snakemake >=8.0
  - multiqc
  - idr
  - seacr
```

**Rationale:**
- `conda-forge` + `bioconda` cover all declared dependencies.
- Removing `defaults` avoids Anaconda commercial-channel warnings in CI
  and reduces mixed-channel dependency resolution risk.
- Matches `workflow/envs/ci-fast.yml` channel strategy.

**Verification:**
Attempt `micromamba create -f workflow/envs/chipseq.yml --dry-run` or
equivalent. If solve is too slow locally, document that full
verification occurs via the `workflow_dispatch` CI job.

### 2.2 No `prefix` metadata

Current `chipseq.yml` has no `prefix` line — already clean.  The
KNOWN_ISSUES item about removing prefix metadata is marked resolved
without any file change.

### 2.3 No further environment split in Stage 8d

Decide whether to split environments further: no further split in
Stage 8d. `chipseq.yml` remains the full runtime environment;
`ci-fast.yml` remains the lightweight CI environment. Further
task-specific envs can be revisited later only if solve/runtime
pain persists.

---

## 3. README execution documentation

### 3.1 New section: "Local execution"

Add under `## Developer Notes`, after the CI section:

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

### 3.2 Existing sections preserved

The existing sections `Smoke-test profiles` and `Tiny real execution`
remain as-is — they document the harness commands. The new `Local
execution` section adds the end-to-end user workflow.

### 3.3 CI badge

Not added in this slice. Deferred until CI passes on `main`.

---

## 4. KNOWN_ISSUES.md updates

The Stage 8 section gets a **Stage 8d completed** block. Replace the
three `⬜` items with:

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

---

## 5. Files changed

| File | Action |
|------|--------|
| `workflow/envs/chipseq.yml` | Edit — remove `defaults` channel |
| `README.md` | Edit — add Local execution section |
| `KNOWN_ISSUES.md` | Edit — Stage 8d complete |

No other files modified. No new files created.

---

## 6. Verification plan

1. `micromamba create -f workflow/envs/chipseq.yml --dry-run` (if
   feasible locally) or document deferral.
2. Run `python3 test/test_stage8_smoke_profiles.py` — unchanged, should
   still pass.
3. Run `python3 test/test_validation_stress.py` — unchanged.
4. Run `python3 test/test_stage8b_tiny_execution.py` — unchanged.
5. `git status --short --untracked-files=all` — only expected files.

---

## 7. Self-review notes

- No new pipeline features — documentation and one channel removal only.
- No further env split in Stage 8d; future task-specific envs can be
  revisited only if solve/runtime pain persists.
- No CI badge added (deferred to first green main run).
- `defaults` removal requires solve validation; documented.
- README command examples use `micromamba create -f workflow/envs/chipseq.yml`
  (not `-n chipseq -f ...`) since the env file already has `name: chipseq`.
- No TBDs or TODOs.
