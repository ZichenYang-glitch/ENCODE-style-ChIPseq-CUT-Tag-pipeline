# Bug Log: Local Bring-Up Session 2026-09-17/18

Bugs and defects found while bringing up the local HelixWeave stack
(ENCODE workflow + bulk RNA-seq adapter) on a WSL2 workstation. Each entry
lists status, root cause, and evidence. The accepted maintenance batch is
committed in `a6df932` through `d1a813d` (eight commits); unresolved entries
below retain their individual status.

## 1. FIXED — bulk-rnaseq authoring schema: gated-section defaults conflict with disabled-state rule

**Symptom:** Every frontend submission of a bulk-rnaseq run config failed with
`BULK_RNASEQ_UMI_CONFLICT` at `config.standard.umi`, even with all default
values untouched.

**Root cause:** The backend semantic validator intentionally requires a
disabled gated section to contain exactly `{enabled: false}` and nothing else
(`validation.py` `_validate_standard_semantics`, pinned by
`test_standard_semantic_conflicts_are_rejected`). But the authoring schema
declared property-level `default`s (and a `const`) inside the gated `umi` and
`ribosomal_rna_removal` objects. The frontend materializes initial form state
with rjsf `getDefaultFormState`, which merges property defaults and `const`
values into the object default, producing e.g.
`{enabled: false, deduplication_tool: "umitools", grouping_method: "directional", emit_dedup_stats: false}`.
The rjsf-generated neutral state therefore always violated the contract.
Existing tests never caught this because they construct payloads directly and
bypass rjsf default materialization.

**Fix (commit `a6df932`):**
`src/encode_pipeline/adapters/bulk_rnaseq/authoring.py`

- `umi`: removed `default` from `deduplication_tool`, `grouping_method`,
  `emit_dedup_stats`, `primary_alignments_only`; changed `deduplication_tool`
  from `const: "umitools"` to `enum: ["umitools"]` (rjsf materializes `const`
  but not `enum`; the accepted value set is unchanged).
- `ribosomal_rna_removal`: removed `default: false` from
  `save_filtered_reads`.
- Normalization already applies the same fallbacks via `.get()` in
  `validation.py`, so generated nf-core params are unchanged and the set of
  accepted configs is unchanged (SCHEMA_VERSION left at 1.1.0).
- Regenerated `execution-implementation-manifest-1.0.0.json` and
  `default-execution-qualification-1.1.0.json` via
  `scripts/generate_bulk_rnaseq_execution_manifest.py` (authoring.py is a
  controlled file in the pinned execution identity).

**Evidence:** rjsf `getDefaultFormState` on the served config schema now
materializes `umi`/`ribosomal_rna_removal` as exactly `{enabled: false}`;
full default config + real sample row + GRCm38 reference passes
`BulkRnaSeqWorkflowAdapter().validate`; 378 adapter/identity/packaging tests
pass.

**Validation limits for the original fix:** this change touches the pinned
execution identity (Protected Bulk Gate territory). Only adapter/execution-
identity/packaging tests were run; no full Bulk Gate evidence was recorded.

## 2. FIXED — ENCODE workflow: relative `scripts/` paths break under platform execution

**Symptom:** Snakemake run launched by the platform failed at
`parse_dup_metrics.py` and similar helper invocations with
`can't open file '[REDACTED]/scripts/parse_dup_metrics.py'`.

**Root cause:** Six rule files invoked helper scripts as
`python3 scripts/<name>.py`. The platform executes Snakemake with
`--directory <workspace>`, so the process cwd is the run workspace, not the
repository root, and the relative path resolved to a nonexistent location.

**Fix (commit `88b5c66`):** 22 call sites across 6 rule files
(`workflow/rules/{consensus,idr,idr_reproducibility,mnase,qc,report}.smk`)
changed to `python3 {workflow.basedir}/../scripts/<name>.py`. Verified with
Snakemake dry-run; confirmed in a real platform run afterwards.

## 3. OPEN (unconfirmed) — ENCODE run: macs3 conda environment activation

**Symptom:** A platform-driven ENCODE run failed in the macs3 step: the rule
invoked a script via `#!/usr/bin/env python` shebang and resolved the wrong
interpreter (PATH/activation timing inside the Snakemake-managed conda env).

**State:** The conda env itself was verified good manually. The user was asked
to rerun; no outcome reported yet. Needs reproduction. If it reproduces,
candidate fix is invoking the env interpreter explicitly (same pattern as
bug 2) instead of relying on shebang + activated PATH.

## 4. FIXED — qc master switch conflicts with materialized sub-flags

**Observation:** The qc section defaults materialize as all-true
(`{enabled: true, fastqc: true, ...}`), which is valid. But if a user toggles
only the master `enabled` off in the form, the sub-flags stay `true` in
formData and submission fails with `BULK_RNASEQ_QC_CONFLICT`
("disabled master switch with enabled sub-flags"). The backend rule is deliberate.

**Fix (2026-09-21):** The shared schema form now cascades an explicit
`enabled: true` to `false` edit using the adapter-owned schema. Boolean child
flags are set to false; sections declaring the neutral default `{enabled: false}`
(including UMI and ribosomal RNA removal) clear their optional settings instead.
Child controls are disabled while the section is off, and the master switch
explains that cleared settings are not restored on re-enable. Trimming retains
its required tool setting. Imported YAML conflicts and unknown keys remain
available to normal validation; backend schemas and semantic checks are unchanged.

**Remaining UX boundary:** Cascading is limited to the section's own settings.
Cross-section dependencies are not cleared: disabling QC while
`advanced.rseqc_modules` is configured still produces
`BULK_RNASEQ_ADVANCED_CONTEXT_CONFLICT`; disabling UMI while
`standard.outputs.umi_intermediates` is true still produces
`BULK_RNASEQ_OUTPUT_CONFLICT`. These deliberate backend checks remain in force.
Inline guidance for cross-section conflicts is follow-up work; FIXED here does
not mean that switching off a section always makes the submission valid.

**Evidence:** Two new form regressions failed before the fix. The final frontend
suite passes 380 tests, including disable/re-enable, neutral-state cleanup,
required-setting preservation, and imported-conflict coverage. Desktop/mobile
browser tests use the served bulk schema and check the exact request preview for
QC, UMI, and rRNA removal; the complete browser suite passes 13 tests. The existing
backend conflict/QC-disable cases pass 12 tests. Type checking and the production
build pass. Package-owned frontend assets are regenerated through the existing
packaging script, with 33 asset/deployment/distribution tests passing;
regenerating the 111-file bulk execution manifest leaves its
identity files unchanged. No real bulk scientific execution or Protected Bulk Gate
was needed for this UI-only behavior change.

## 5. DESIGN FRICTION — Docker 29 storage-driver mismatch between staging and runtime admission

**Observation:** The two Docker touchpoints expect opposite storage semantics
on one daemon:

- staging `_verify_pulled_image` expects classic overlay2 semantics
  (`inspect Id == config digest` after `docker pull`);
- runtime admission `_verify_docker_availability` expects containerd-storage
  semantics after `docker load` of the staged archive
  (`inspect Id == archive index digest`).

A single Docker 29 daemon cannot satisfy both; switching
`features.containerd-snapshotter` flips which side fails
(`docker_image_invalid` vs `runtime_admission_failed`).

**Workaround deployed on this machine:** two daemons — system dockerd
(overlay2, staging) plus a rootless dockerd (containerd snapshotter,
`~/.helixweave/docker.sock`) as the managed/admission Docker, selected via
`ENCODE_PIPELINE_MANAGED_DOCKER_SOCKET`. Works; admission passes
(56 container bindings verified). Worth documenting in
`docs/development/local-platform-runtime.md`, or reconsidering the staging
check so one containerd daemon serves both.

## 6. TOOLING FRICTION — checkout_bootstrap rejects the ci-fast env

**Symptom:** `python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest ...`
refused to run: `source provenance check failed [pth_mapping_unsafe]: remove
unrecognized executable .pth startup hooks`.

**Context:** `.local/envs/ci-fast` carries the editable install `.pth` used
for local platform work. Workaround used: run pytest directly via
`.local/envs/ci-fast/bin/python -m pytest`. Decide whether the provenance
check should tolerate a known task-local editable mapping or the docs should
prescribe the direct invocation for this env.

## 7. FIXED — bulk execution assumes container uid == deployment host uid

**Symptom:** First real bulk-rnaseq run failed every task at launch:
`bash: .command.run: Permission denied` in each task's `.command.err`.

**Root cause:** The generated `platform.nextflow.config` pins
`docker.runOptions = '--user=<deployment os.getuid()>:<gid>'`
(`execution.py:841`, default at `execution.py:141-142`, not operator
configurable in `deployment.py`). This assumes the managed Docker shares the
host user namespace (rootful Docker), so container uid 1000 is the workspace
owner. Under rootless Docker the container uid 1000 maps to host subuid
100999 (verified empirically: file created in-container lands as 100999 on
the host), which can neither traverse `/home/yangzichen` (0750) nor write
into workspace task directories (0755, owned by the host user).

**Workaround deployed:** keep rootless Docker; `chmod o+x /home/yangzichen`
(traverse-only) plus inherited ACLs granting subuid 100999 access to the
platform workspaces directory (`setfacl -m u:100999:rwx -m d:u:100999:rwx`).
Requires the `acl` package. Reversible with `setfacl -Rb`.

**Second symptom (2026-09-19):** ACLs fixed container launch and reads, but
STAR creates its temp dirs (`_STARgenome`, `_STARpass1`) with explicit mode
0700 owned by subuid 100999; the host-side Nextflow output collector
(`fetchResultFiles`) then dies with `AccessDeniedException` walking the task
dir, failing STAR_ALIGN after the aligner itself succeeded. Explicit
`mkdir(0700)` overrides inherited ACL masks, so ACLs cannot fix this class.

**Interim shim deployed:** `~/.helixweave/bin/docker` rewrites the pinned
`--user=1000:1000` to `--user=0:0` (container uid 0 maps back to the host
user under rootless Docker); `ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE`
points at the shim. Admission re-verified OK (56 container bindings).
The shim is now eligible for retirement using the explicit coordinates below;
the deployed shim file was not changed or deleted by PR-5a.

**Formal fix (PR-5a, 2026-09-21, commit `55ef928`):**

- `deployment.py` accepts optional `ENCODE_PIPELINE_BULK_CONTAINER_UID` and
  `ENCODE_PIPELINE_BULK_CONTAINER_GID` and passes them into the existing
  `BulkRnaSeqExecutionBinding.container_uid` / `container_gid` fields.
  Each omitted coordinate retains its existing `os.getuid()` / `os.getgid()`
  default independently. Explicit `0` is preserved.
- Values must contain only ASCII decimal digits for a non-negative integer.
  Empty strings, whitespace, signs, floats, Unicode digits, and non-string
  values fail closed: authoring remains available, execution becomes
  unavailable, and rejected operator values are not exposed publicly.
- The existing validation and `--user=<uid>:<gid>` rendering in `execution.py`
  already support explicit identities and are reused unchanged. Other Docker
  options, default resources, and ENCODE command construction are unchanged.
- Regenerated both execution identity files. The 111-file manifest changes
  only the `deployment.py` entry; prior authoring/QC fixes and the persistence
  contract remain preserved.

**Evidence:**

```bash
./.local/envs/ci-fast/bin/python -m pytest \
  test/adapters/test_bulk_rnaseq_execution.py \
  test/adapters/test_bulk_rnaseq_adapter.py \
  test/adapters/test_bulk_rnaseq_deployment.py \
  test/adapters/test_bulk_rnaseq_reference_profiles.py \
  test/adapters/test_bulk_rnaseq_execution_identity.py \
  test/packaging/test_bulk_rnaseq_contract_package.py -q
# 469 passed in 6.91s
```

Tests cover omitted coordinates, explicit `0:0`, other numeric IDs, independent
UID/GID defaults, malformed/partial coordinates, and generated Docker options.
`ruff check` and `ruff format --check` pass for all three changed Python files.
PR-5a diffs are whitespace-clean; the whole-worktree check still reports the
untouched pre-existing trailing whitespace in `config/samples.tsv:2`.

**Deployment follow-up:** Set both new coordinates to `0` for this rootless
deployment, point `ENCODE_PIPELINE_MANAGED_DOCKER_EXECUTABLE` at the real Docker
executable, then restart the API/worker and recheck admission. The old argument
rewriting shim can then be retired by the deployment operator. No deployment
configuration, shim, ACL, or live service was changed in this PR.

**Validation limits:** The Protected Bulk Gate and a real rootless scientific
run were not performed; the evidence above covers deployment loading,
validation, configuration generation, and execution identity contracts.

## 8. FIXED — ENCODE platform runs pin `--cores 1` while smk rules default to `threads: 8`

**Observation (code review, 2026-09-19):** The adapter's `_validate_options`
only allows `strict_inputs` (`adapters/encode.py:1093`) and the authoring
schema is `additionalProperties: false` (`encode_authoring.py:265`), so the
`options.cores` branch in the service-layer command builder
(`command_builder.py:562`) is unreachable on the platform path: every
platform-driven ENCODE run gets `--cores 1`. Snakemake rules declare
`threads: THREADS` (default 8, `workflow/rules/common.smk:40`). Empirical
evidence from the 2026-09-17 platform run (18/32 steps completed) suggests
Snakemake caps threads to available cores rather than failing — meaning
alignment/peak steps silently ran single-core. The only real-execution e2e
profile sidesteps this with `threads: 1`.

**Impact:** severe under-utilization (multi-core aligners running on one
core); no correctness impact expected.

**Fix (PR-3, commit `0f01241`):**

- Added optional workflow option `cores` to the ENCODE authoring schema:
  integer, minimum 1, maximum 1024, default 1. `additionalProperties: false`
  and the existing `strict_inputs` option remain unchanged. Advanced the
  versioned authoring contract from 1.2.0 to 1.3.0 and updated the existing
  adapter/API schema-version assertions.
- `_validate_options` now accepts `cores` and enforces the same upper bound,
  rejecting booleans, non-integers, null, and values outside 1–1024.
- At PR-3, the service-owned CommandBuilder received validated `cores`;
  its omitted-option default remains `--cores 1`. CommandBuilder and workflow
  rule thread semantics were not changed.
- Added schema/adapter boundary cases and a regression test through real
  adapter validation, WorkspacePlanner, and CommandBuilder: explicit
  `cores: 8` yields `--cores 8`, while omitted `cores` yields `--cores 1`.
  Before the fix, the explicit-cores schema and planning cases failed and
  the omitted-option cases passed.
- Follow-up commits `2c26a9e` and `02303b2` advance the authoring schema to
  2.0.0 for stage retirement and move command construction into the ENCODE
  adapter respectively, preserving the cores contract.

**Evidence (2026-09-21):**

```bash
./.local/envs/ci-fast/bin/python -m pytest test/adapters/ -k encode -q
# 124 passed, 860 deselected in 4.43s
./.local/envs/ci-fast/bin/python -m pytest \
  test/services/test_command_builder.py \
  test/api/test_routes_workflows.py test/api/test_openapi_export.py -q
# 86 passed, 1 warning in 75.94s
```

- `snakemake -s workflow/Snakefile --configfile config/config.yaml -n`
  succeeded with 32 jobs, using `ci-fast/bin` on PATH and a temporary writable
  `XDG_CACHE_HOME`. The first attempt hit the sandbox's read-only default
  Snakemake cache; no workflow or deployment configuration was changed.
- Ran `./.local/envs/ci-fast/bin/python
  scripts/generate_bulk_rnaseq_execution_manifest.py`. Both generated files
  and their hashes remained byte-identical to the pre-PR working tree: the
  changed ENCODE source files are outside the bulk manifest's 111-file
  allowlist. The existing bulk authoring/QC fixes remain preserved.
- `ruff check` and `ruff format --check` passed for all five changed Python
  files. PR-3 file diffs are whitespace-clean. The whole-worktree
  `git diff --check` still reports the untouched, pre-existing trailing
  whitespace in `config/samples.tsv:2`.
- The pytest warning is the existing relative `PYTHONTZPATH` environment
  setting. OpenAPI contract checks passed; no generated-client changes were
  needed.

**Validation limits:** No real scientific execution or browser run was
performed, so actual process thread counts were not measured. The 1024-core
limit is an adapter contract bound, not host-capacity detection. The
Protected Bulk Gate was not run; bulk execution identity is unchanged.

## 9. FIXED — bulk QC parser rejects RSeQC bam_stat long-label line (QC indexing fails on real runs)

**Symptom:** First successful bulk run (2026-09-19, SRX21122275) completed all
45 tasks and artifact extraction (54 artifacts), but QC metrics indexing
failed: `qc_metrics_indexing_failed` / `QC_INDEXING_ADAPTER_FAILED`. Fail-closed
behavior worked as designed — nothing invalid was persisted.

**Root cause:** `_rseqc_bam_stat_counts` (`adapters/bulk_rnaseq/qc.py:1638`)
requires `\s+` after the colon in
`Proper-paired reads map to different chrom:`. Real RSeQC output pads fields
to a fixed width; this label fills the entire width, so the digit follows the
colon with no space (`...different chrom:0`). The parser rejected the whole
document (`source_content_invalid`). Existing test fixtures always carried
padded lines, so the unpadded real-world form was never covered.

**Fix (PR-2, commit `494b979`):**

- In `src/encode_pipeline/adapters/bulk_rnaseq/qc.py`, changed only the
  long-label pattern's post-colon padding from `\s+` to `\s*`.
- Parameterized the existing single-end bam_stat test in
  `test/adapters/test_bulk_rnaseq_qc.py` with both the original padded line
  and the exact real-output line `Proper-paired reads map to different chrom:0`;
  existing metric assertions are unchanged.
- Regenerated `execution-implementation-manifest-1.0.0.json` and
  `default-execution-qualification-1.1.0.json` using
  `./.local/envs/ci-fast/bin/python scripts/generate_bulk_rnaseq_execution_manifest.py`.
  Relative to the pre-PR working tree, only the `qc.py` manifest entry changed;
  the existing authoring fix identity and persistence contract were preserved.

**Evidence (2026-09-21):** Before the parser fix, the padded case passed and
the new unpadded case failed with `source_content_invalid`. After the fix:

```bash
./.local/envs/ci-fast/bin/python -m pytest \
  test/adapters/test_bulk_rnaseq_qc.py \
  test/adapters/test_bulk_rnaseq_execution_identity.py \
  test/packaging/test_bulk_rnaseq_contract_package.py -v
# 154 passed in 2.99s
```

`ruff check` and `ruff format --check` passed for both changed Python files;
`git diff --check` passed for the PR-2 files. The whole-worktree check still
reports pre-existing trailing whitespace in `config/samples.tsv:2`; that file
was left unchanged. The Protected Bulk Gate and real-run QC reindexing were
not run for this PR.

**Operational follow-up:** After the fix and a stack restart, this run's QC
can be reindexed without re-execution:
`QcSummaryIndexingService.index(run_id, artifacts)` begins a fresh attempt
when called without attempt identity (currently only invoked by the worker at
`workers/jobs.py:406`; no user-facing retry endpoint exists).

**Wider review note:** audit the other fixed-width parsers in `qc.py`
(infer_experiment, read_distribution, featurecounts, salmon meta_info) for
the same `\s+`-after-label assumption.

## 10. FIXED — CI formatting and stale assertions after the maintenance batch

**Symptom:** CI for `d1a813d` and `4c64f0b` reported the same three failure
signatures: snakefmt formatting, the startup project-root assertion in shard
2, and the shared options assertion in three durable browser scenarios.
These are verification drift, not production execution defects.

**Fix:**

- Ran snakefmt 2.0.3 only on `workflow/rules/metadata.smk`. With the repository's
  `line_length = 120`, the formatter wrapped two long argument lists. The
  before/after Python ASTs are identical.
- `test/api/test_startup.py` now checks the registered ENCODE adapter's
  `_execution_binding.project_root` against the build identity's project root.
  PR-5b moved that coordinate out of the service CommandBuilder; the alignment
  assertion remains exact.
- `frontend/e2e/durable-run.spec.ts` now expects exactly
  `{strict_inputs: false, cores: 1}` in the submitted options. PR-3 added the
  schema's default `cores`; the desktop success, mobile history, and mobile
  cancellation scenarios share this assertion.

**Process root cause:** Earlier Codex local verification omitted the complete
lint and browser-e2e tiers: the required local snakefmt, Node 22, and Playwright
environments were not configured for those commands. Targeted Python tests
and frontend type checks did not cover these failures. Maintenance changes to
commands, workspace contracts, or materialized defaults need the corresponding
lint and full browser tier in the local verification matrix.

**Local environment repair (prepared before this PR):** The operator converted
47 relative interpreter shebangs in `.local/envs/ci-fast/bin/` to absolute
paths, made sysconfig `TZPATH` absolute, and provisioned `.local/node22/` and
`.local/ms-playwright/`. This PR confirmed no remaining ci-fast-relative
shebangs, absolute timezone paths, Node v22.23.1, and the installed Chromium
runtime. These deployment-local repairs are recorded here for reference from
`docs/development/local-platform-runtime.md`; no environment files were
changed by this PR.

**Evidence:**

```bash
./.local/envs/lint-snakefmt/bin/snakefmt --check workflow/rules/metadata.smk
# All 1 file(s) would be left unchanged; exit 0.
./.local/envs/lint-snakefmt/bin/snakefmt --check workflow/
# 8 file(s) would be left unchanged; no file requires formatting.
# 5 parsing errors because local shfmt is absent; exit 123 (limitation below).
./.local/envs/ci-fast/bin/python -I -S scripts/checkout_bootstrap.py \
  --repository-root . pytest test/api/test_startup.py -q
# 8 passed in 16.40s
cd frontend
PATH="$PWD/../.local/node22/bin:$PWD/../.local/envs/ci-fast/bin:$PATH" \
  PLAYWRIGHT_BROWSERS_PATH="$PWD/../.local/ms-playwright" \
  ENCODE_PIPELINE_E2E_REDIS_URL=redis://127.0.0.1:6379/14 \
  npm run test:e2e
# 11 passed (51.7s)
```

The full browser suite found no cascading assertion failures. In particular,
`8 loaded of 8` and the exact API metric count of 8 both passed: they count QC
metrics, not workspace files, so `config/encode-execution.json` does not change
these assertions. Ruff check and format check passed for the changed Python
test. No `src/` files or execution identity manifests changed.

**Validation limits:** The local workflow-wide snakefmt check cannot parse
`common.smk`, `replicates.smk`, `mnase.smk`, `peaks.smk`, and
`idr_reproducibility.smk` without `shfmt`; this is not a fully green local lint
tier. No additional formatting changes were reported, and the edited file
passed independently. The browser suite used the existing controlled ENCODE
runtime and bulk authoring/unavailable path, not the Protected Bulk Gate.

**Lint follow-up (2026-09-21):** The 22 script calls changed by `88b5c66`
now use a shared `SCRIPTS_DIR` through rule-local `params.scripts_dir`, removing
the new direct-`workflow` warnings without changing rendered script-path bytes
or the 32-job DAG. With user approval, the baseline received only 23 location
and rule-name updates, including the approved `stage3_qc_summary` to
`project_qc_summary` rename. All 65 existing warning types and bodies are
unchanged; `python test/check_snakemake_lint.py` now exits 0 with
`snakemake --lint output matches baseline.`
The complete artifact/workflow tests passed (131 tests, collected in that order).

**Build-identity follow-up:** `3d9cf61` removed the literal `scripts/<name>.py`
references from workflow bytes, so the existing build-identity scanner omitted
15 runtime scripts from `source_manifest`. Each of the 22 rule calls now binds
its full literal `../scripts/<name>.py` path in `params.script` and uses
`{params.script:q}` in the shell; the unused `SCRIPTS_DIR` constant was removed.
All 22 resolved path values and rendered path tokens are byte-identical to the
previous version. The packaging contract now independently requires all 15
script names in the manifest; this assertion failed before the rule fix and
passed afterward, together with all nine release-distribution tests. The lint
check exits 0 against the unchanged baseline, and the dry-run remains 32 jobs.
No `src/` files or execution identity manifests changed. The earlier verification
matrix omitted the packaging tier, which is now included for workflow path changes.
Full packaging verification passed: 150 tests in 43.16s outside the sandbox;
the frontend/API subprocess had timed out inside the sandbox. The artifact/workflow
suite passed 131 tests in 64.60s, with artifacts collected first because the reverse
order still exposes the existing `workflow.lib` import error during collection.
Ruff checks passed; snakefmt reported no formatting changes and the same five
missing-`shfmt` parse errors documented above.

## Not bugs (recorded to avoid re-investigation)

- Frontend sample TSV for profile-bound workflows must NOT contain
  `genome`/`bowtie2_index` columns: once a reference profile is selected the
  platform injects those values, and extra sample columns are rejected. This
  is intentional fail-closed behavior.
- `docker images` lists nothing after digest-only `docker load` under the
  containerd snapshotter (no repo tags); `docker inspect <digest>` works.
  Docker behavior, not a platform bug.
