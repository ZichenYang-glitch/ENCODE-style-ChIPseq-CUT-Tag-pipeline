# 01 — Current state

Audit date: 2026-10-29. Base: `main` at `fa0d840` (merge of PR #183).
Claim tags: **observation**, **hypothesis**, **proposal**, **tested**.

This note answers critical questions 1 and 3. Human attention counts, token
cost, and CI first-pass rate are **not stored in the repository**
(**observation**). Final GitHub check rollups are not a substitute for
first-pass rate.

## Architecture

**Observation.** HelixWeave is a workflow-neutral local omics platform.
[`docs/architecture/platform-overview.md`](../../architecture/platform-overview.md)
assigns ownership as follows:

| Layer | Owns |
| --- | --- |
| `src/encode_pipeline/platform/` | Domain contracts |
| `src/encode_pipeline/services/` | Orchestration |
| `src/encode_pipeline/api/` | HTTP translation |
| `src/encode_pipeline/persistence/` | SQLite repositories and Alembic |
| `src/encode_pipeline/workers/` | Redis/RQ mechanics |
| `src/encode_pipeline/adapters/` | Workflow schemas, commands, artifacts, QC |
| `workflow/` | Snakemake rules |
| `frontend/` | Generated-client UX |

SQLite is lifecycle truth. Redis/RQ carries a durable run id. The product
roadmap in
[`docs/development/workflow-platform-agent-roadmap.md`](../../development/workflow-platform-agent-roadmap.md)
states that further in-product Agent write actions, PostgreSQL, Kubernetes,
and a Bulk redesign are non-goals for the current boundary.

**Observation.** One ENCODE run crosses these seams:

1. `EncodeStyleWorkflowAdapter.schema()` → `build_encode_authoring_schema()` in `adapters/encode_authoring.py`.
2. `GET /{workflow_id}/schema` in `api/routes/workflows.py`.
3. `ValidatedInputService.validate()` persists `ValidatedInputSnapshot` via `create_validated_input_snapshot()`.
4. `RunService.transition_run()` calls `update_run(..., expected_status=...)`. A `rowcount != 1` raises `ConcurrentRunUpdateError`.
5. `RunSubmissionService.start_run()` enqueues `run_execution_job` through `RedisRunQueue.enqueue_execution()`.
6. The worker claims with `claim_execution_assignment()` only from allowed statuses.
7. `CommandBuilder.build_command()` stops at `adapter.build_command()`. Snakemake argv is built in `adapters/encode_execution.py`.
8. After success, `ArtifactExtractionService` and `QcSummaryIndexingService` write `RunArtifactRow` and `RunQcMetricRow`.

Alembic head is `20260827_15` in
`src/encode_pipeline/persistence/alembic/versions/20260827_15_terminal_email_notifications.py`.

The frontend must not hand-edit `frontend/src/api/generated/`. Regeneration is
`npm --prefix frontend run openapi:regenerate`. Drift is checked by
`test/api/test_generated_client_drift.py` and `test/api/test_openapi_export.py`.

## Cohesion and coupling

**Observation.** `platform/` and `api/` do not import `encode_pipeline.adapters`.
`services/defaults.py` does: it constructs `EncodeStyleWorkflowAdapter`, the
bulk and HiTRAC deployment bindings, and `configure_encode_execution`. That is
the composition root, not a leak of rule names into domain code.
`cli/results_visibility_fixture.py` imports `encode_pipeline.adapters.encode_qc._QC_HEADER`,
which is a private symbol used by the deterministic demo. **Hypothesis:** that
private import is accidental coupling, limited to the demo fixture.

**Observation.** Scientific complexity is essential: assay policies, reference
bindings, QC parsers, and nf-core identity live in adapters and `workflow/`.
Engineering complexity that is justified includes compare-and-set lifecycle,
checkout provenance (`docs/development/source-provenance.md`), redaction, and
the Protected Bulk gate. Accidental complexity found in this pass is narrow:
the demo's private QC import, and the repeated need to rediscover which pytest
marker a path belongs to (`full_main` versus PR-fast). Agent coordination
complexity is separate and is documented in notes 03 and 04.

**Observation.** The in-product Agent (`services/agent.py`) is not a coding
agent. Its tools are `list_workflows`, `get_workflow_schema`,
`validate_workflow_input` (no snapshot persist), and `explain_issues`.
`ReadOnlyToolRegistry` rejects mutating prefixes. Do not extend it in order to
coordinate Codex, Peri, Kimi, or Cursor.

## Verification already in the repo

**Observation.** [`docs/development/harness.md`](../../development/harness.md) already splits
risk:

- PR-fast: two hash shards, excludes `full_main` and all three real-execution markers. Changed Python lines need 80% coverage.
- Full main: includes `full_main`.
- Platform real, scientific real, Bulk real, and container smoke are dispatch or nightly, not default PR gates.
- Lint, frontend, browser e2e, and lock-check are separate workflows.

Recommended required contexts: `fast-checks`, `frontend`, `browser-e2e`,
`lint`, `lock-check`, `coverage`. Soft budget for a PR shard is 5 minutes;
hard timeout is 25 minutes.

Worktree rule already written: one `.local/envs/ci-fast` prefix per checkout.
A shared editable install can make Snakemake import the wrong checkout while
pytest fails closed. That is an agent foot-gun and is already documented.

## Existing agent assets

| Asset | Status | Disposition |
| --- | --- | --- |
| `AGENTS.md` (7755 bytes on this commit) | Canonical coding-agent instructions. Rejects swarm frameworks. Names bootstrap commands and co-author trailers. | Keep. It is under Codex's documented 32 KiB project-doc cap. |
| Product roadmap | Delivery sequence and non-goals. | Keep as product authority. Do not copy it into agent prompts. |
| `docs/development/harness.md` | Risk-to-command map. | Keep. Point agents here instead of writing a second policy. |
| Architecture decision notes under `docs/architecture/` | Contracts. | Keep. Retrieve by task, do not preload. |
| `CLAUDE.md`, `.cursor/skills`, PR template, issue template | Absent. | Do not add a second copy of `AGENTS.md`. A PR template is optional and not required for the recommendation in note 08. |
| In-product Agent | Read-only explainer. | Keep. Do not route coding work through it. |

**Observation.** Commit trailers in the last 120 commit messages:

| Trailer | Count |
| --- | ---: |
| `Co-Authored-By: Codex <noreply@openai.com>` | 14 |
| `Co-Authored-By: DeepSeek Flash` | 9 |
| `Co-Authored-By: deepseek-flash <deepseek-ai@claude-code-best.win>` | 8 |
| `Co-Authored-By: Peri` | 4 |
| `Co-Authored-By: Kimi <noreply@moonshot.cn>` | 1 |
| Malformed deepseek trailer | 1 |

DeepSeek appears in history even though it was not in the requested vendor
list. The email domain suggests another harness, not a HelixWeave subsystem
(**hypothesis**). The repository already tells agents to use each tool's
default trailer (`b0bcb8b`).

## Representative changes

Final check rollups below are the merged PR's last state from `gh pr list`,
not the first attempt. Every one of these 15 merged PRs ended `SUCCESS` on
the reported checks, with the real-execution jobs `SKIPPED` (4 skipped on the
newer PRs). That does **not** measure first-pass rate.

| PR | Head | Files | + / − | Shape |
| --- | --- | ---: | --- | --- |
| [#169](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/169) | `agent/helixweave-governance-roadmap-sync` | 2 | +129 / −85 | Docs only |
| [#174](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/174) | `agent/node-toolchain-authority` | 6 | +13 / −9 | Toolchain pin |
| [#167](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/167) | `codex/pr165-test-only-cleaner-race` | 1 | +5 / −0 | One test |
| [#183](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/183) | `agent/pr-k-reference-fai` | 3 | +224 / −30 | Adapter fix |
| [#166](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/166) | `agent/pr166-scientific-real-bootstrap` | 4 | +784 / −30 | CI provenance |
| [#163](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/163) | `agent/worker-hard-timeout-finalization` | 5 | +860 / −37 | Worker |
| [#182](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/182) | `agent/a1-config-compat` | 13 | +1477 / −6 | Config contract |
| [#175](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/175) | `agent/deployment-coverage-correction` | 18 | +3625 / −89 | Coverage ratchet |
| [#164](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/164) | `agent/current-checkout-source-provenance` | 22 | +3758 / −45 | Provenance |
| [#165](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/165) | `agent/adapter-admission-fail-closed` | 38 | +1848 / −125 | Adapter admission |
| [#154](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/154) | `agent/pr154-bulk-rnaseq-product-integration` | 84 | +5348 / −235 | Cross-module product |
| [#173](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/173) | `agent/helixweave-deployment-runtime-cli-staging` | 121 | +58304 / −114 | Deployment CLI, large additions |
| [#176](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/176) | `agent/helixweave-lan-auth-staging` | 139 | +10991 / −297 | Auth, high risk |
| [#170](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/170) | `agent/helixweave-reference-profiles` | 146 | +15090 / −980 | Reference profiles |
| [#181](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/181) | `agent/frontend-commonjs-build-determinism` | 220 | +832 / −735 | Frontend plus generated assets |

**Observation.** October 2026 also merged several `agent/pr-*` branches with
`git merge` rather than a GitHub PR number. Comparing each merge's second
parent with the first parent showed one commit and 3–34 files for the FAI,
IDR, QC, Picard, and layout fixes. The analysis-chain merge range contained
three commits because earlier feature commits were still in that range; those
commits were also merged on their own branches. Do not treat that range as
three repair attempts.

**Observation.** Hot paths since 2026-08-01, by commit touches, include the
nf-core execution manifests, `workflow/Snakefile`, `adapters/encode.py`,
generated frontend assets, `docs/configuration.md`, the roadmap, CI, and
migration tests. Generated files and lock/manifest JSON inflate "files
changed" without implying design coupling.

## What to preserve

Preserve the layer boundaries, SQLite-versus-Redis rule, checkout bootstrap,
risk-tiered CI, co-author convention, and the read-only product Agent.
Do not replace them with an agent orchestrator.
