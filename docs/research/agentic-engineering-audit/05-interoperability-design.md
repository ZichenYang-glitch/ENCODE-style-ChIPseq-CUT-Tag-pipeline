# 05 — Interoperability

## What is actually portable

| Mechanism | Who reads it | Limit |
| --- | --- | --- |
| Root `AGENTS.md` | Codex (verified), Kimi (verified), Cursor (observed in this session) | Codex stops at 32 KiB. This file is 7755 bytes. Peri was **not** shown, on the fetched Peri home page, to load it. |
| `~/.agents/AGENTS.md` and `~/.agents/skills` | Kimi documents these as cross-tool locations | Codex uses `~/.codex/AGENTS.md`, not `~/.agents`. Two global files can drift. |
| Nested `AGENTS.md` | Codex walks from the repo root to the cwd | Do not add nested copies unless a directory has a different test command. |
| GitHub PR | All of them, via `gh` or the web UI | This is the handoff that does not depend on ACP. |
| Git worktree plus a private `.local/envs/ci-fast` | Any agent with a shell | Required by `harness.md`. Cursor can automate setup later with `.cursor/worktrees.json`. That file would be Cursor-specific. |
| MCP | Peri and Kimi are clients. Cursor can host MCP. | Do not put repository policy in an MCP server. The policy is already markdown. |
| ACP | Peri (`peri acp`) and Kimi (`kimi acp`) | Cursor does not speak ACP. An ACP bus would not include Cursor without an adapter. |
| Kimi `.kimi-code/agents/*.md` | Kimi only | A read-only reviewer file is optional and not required for the recommendation. |
| Peri Ultracode JS workflows | Peri only | Vendor lock-in. Out of scope. |

**Proposal.** The portable contract is markdown in the PR plus `AGENTS.md`. Do not add a task database.

## Task card

Paste this into the PR body or a commit that starts the work. Agents that can
read the PR do not need a relay.

```markdown
## Task HW-example-sci7-links

- Goal: Make relative links in the SCI-7 reference-identity note resolve.
- Risk: low (docs only).
- Scope: `docs/development/sci7-public-data/reference-identity.md` and the files those links must point at.
- Exclusions: no workflow rule edits, no CI image changes, no Bulk gate.
- Invariants: `test/docs/test_internal_links.py` fails closed on a broken relative link.
- Depends on: none.
- Acceptance: the docs test passes; no other files change.
- Tests: `python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest test/docs/test_internal_links.py -q`
- Owner: one implementer. Reviewer: none required at low risk if that command passed.
- Status: open
- Evidence: Actions run 40686146422, shard 2, assertion listing `ENCFF356LFX.bed.gz` and `ENCFF547MET.bed.gz`.
- Blockers: none
- Handoff: do not rerun the fast-check shard until this test passes.
```

**Tested.** Given only that card, the next command is the 0.73 s docs test from note 03. The roadmap is not required. The card was not handed to a second commercial agent (**untested** across vendors).

## Review card

```markdown
## Finding F1

- Severity: high if merged, because CI shard 2 fails.
- Location: `docs/development/sci7-public-data/reference-identity.md` links to `metadata/confirmation/ENCFF356LFX.bed.gz` and `ENCFF547MET.bed.gz`.
- Invariant: relative documentation links resolve (`test/docs/test_internal_links.py`).
- Explanation: the link targets are not in the tree the test checks.
- Evidence: run 40686146422, `AssertionError: broken relative link(s)`.
- Correction: point the links at files that exist, or remove the links. Do not weaken the test.
- Verify: the task-card pytest command, expect exit 0.
- Disposition: open
```

**Proposal.** The implementer applies one finding and reruns the named command. The human does not retype the assertion. If two reviewers disagree and neither cites a command, the human keeps the disposition `unresolved` rather than voting.

## Context packaging

**Proposal.** Default packet for any agent:

1. `AGENTS.md` (already auto-loaded where the tool supports it).
2. The task card.
3. `docs/development/harness.md` only if the card does not name a command.
4. The files in Scope, plus the test file named in Tests.

Do not attach the roadmap, the bug log, or the whole `docs/architecture/` tree. Pull a decision note only when the task names that contract (auth, reference profiles, publication, deployment).

Kimi sub-agents already isolate exploration. Use that, or a read-only first turn, so search logs do not fill the implementation context. Codex's 32 KiB cap is a reason to keep `AGENTS.md` short and leave procedures in `harness.md`.

## Risk to command

This table restates `harness.md`. It exists so a task card can point at a row.

| Risk | Touch | Before push | Not required on the PR |
| --- | --- | --- | --- |
| Low | Docs, comments, agent instructions | `pytest test/docs/test_internal_links.py` if docs links can change. Otherwise `git diff --check`. | Real execution, Bulk |
| Medium, Python | One package under `src/encode_pipeline/` or `test/` | Bootstrap pytest on the nearest test file; `ruff check` and `ruff format --check` on the touched paths | Full `pytest test`, Bulk |
| Medium, workflow smoke | `workflow/` rules or profiles used by the tiny profile | `pytest test/workflow/test_smoke_profiles.py` (7.39 s here) | `full_main` DAG snapshots unless the task changes targets |
| Medium, frontend | `frontend/` excluding generated client | `npm --prefix frontend test -- --run`, `typecheck` | Playwright, unless the route or shell changed |
| High, contract | OpenAPI, generated client, Alembic, lifecycle, auth, reference binding | The medium commands plus the contract test that already exists (`test_openapi_export.py`, migration tests, auth route tests). Human reads the diff. | Protected Bulk, unless the nf-core identity or execution closure changed |
| High, scientific | Rule outputs, QC meaning, reference bytes | Smoke profile plus the specific DAG test. Human decides whether to dispatch `real_execution`. | Do not dispatch Bulk from an agent |

`snakefmt` and `test/check_snakemake_lint.py` apply when `workflow/` rules change. They are not in the `ci-fast` env's advertised runtime extras; do not pretend they ran if the binary is absent.
