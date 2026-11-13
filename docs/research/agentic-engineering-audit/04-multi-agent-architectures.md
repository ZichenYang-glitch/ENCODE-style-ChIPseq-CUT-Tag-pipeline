# 04 — Multi-agent architectures

**Proposal** unless marked otherwise. No architecture in this note was executed
as a live multi-vendor loop.

## Approaches

| Approach | Fit for HelixWeave | Why |
| --- | --- | --- |
| A. Generalist | Default for one module or one docs fix | History shows one-commit agent branches. A second agent adds relay cost. |
| B. Fixed specialist roles | Useful only for review after the diff exists | Kimi can enforce a read-only reviewer. Codex `/review` is the same idea. Do not assign a standing "architect agent". |
| C. Dynamic routing by score | Not now | No outcome table exists. Routing on vibes is worse than routing on path risk. |
| D. Parallel independent work | Already in use for disjoint branches | `agent/pr-g` through `agent/pr-k` merged serially. Keep this for disjoint files. |
| E. Orchestrator | Not now | See note 08. Peri Ultracode and Kimi Tower are real products and are still vendor runtimes. |

Hybrid that matches the repo: **path risk selects the gate; the human picks the agent that is already logged in; a second context reviews only medium and high risk.** That is static risk policy plus dynamic availability, not a learned router.

## Workflows

### Workflow 1 — Single agent

Specification, implementation, targeted tests, PR.

- Use for docs, a single adapter fix, a test ratchet, a CI script that does not change scientific outputs.
- Human: goal, risk class, final look at the diff.
- Overhead: lowest.
- Failure: the same agent can rationalize a weak test. Mitigation is the existing CI shard, not a second model.

### Workflow 2 — Implementer plus independent reviewer

A fresh session sees the diff, the task card, and the command output. It does not see the implementer's transcript.

- Use when the diff touches `api/`, `persistence/`, `workers/`, generated client, or an adapter contract.
- Human: reads the finding list, not a restatement of the diff.
- Overhead: one extra context. Kimi docs note that this spends tokens. Worth it when a bad merge is expensive.
- Quality: a second model is not independent proof. The finding must name a failed command or a broken invariant.

### Workflow 3 — Review then bounded repair

The finding card goes back to the implementer. Cap repairs at two. If the same test fails with the same assertion, stop.

- Use after Workflow 2 finds a reproducible miss.
- Do not use to "fix" a red `main` by editing an unrelated PR. Note 03 records that failure mode.

### Workflow 4 — Parallel branches

- Use when file sets are disjoint: one workflow rule versus one frontend string versus one docs page.
- Do not use when both sides change `openapi.json`, an Alembic revision, or `adapters/encode_authoring.py`.
- Integration cost: worktrees need separate `.local/envs/ci-fast` prefixes. `harness.md` already says a shared editable install makes Snakemake silently import the wrong checkout.
- This is the pattern the git history already uses. Formalizing it means a file-ownership line in the task card, not a new daemon.

### Workflow 5 — Parallel solutions

- Use for a design that has two real options and a test that can rank them (for example two ways to reject a bad FAI, with a fixture that distinguishes them).
- Do not use for "three agents, pick the majority". Majority is not evidence.
- Maintenance: throw away the losing worktree. Do not merge both.

### Workflow 6 — Coordinator

- Use only if notes from 10–20 real PRs show the human spending attention on relaying the same finding between tools.
- Until then, GitHub is the coordinator: one PR, one head branch, required checks.

## Ownership

**Proposal.**

- One integration owner per PR. That owner may be a human or the agent that opened the PR. It is the only writer after review starts.
- Research agents may read in parallel. They do not commit to the same branch.
- Shared contracts have a single writer: OpenAPI and `frontend/src/api/generated/`, Alembic head, reference-profile binding schema, lifecycle status enum.
- Reviewers do not push to the implementation branch. They return a finding card.
- High-risk merge stays human. See note 07.

Parallel implementation is useful when the task card lists disjoint paths. It is wasteful when the second agent will rebase through a migration or a generated client.

## Feedback loop limits

**Proposal.**

| Rule | Value |
| --- | --- |
| Repair attempts on the same failure signature | 2 |
| Then | Stop and post the assertion, command, and diff stat |
| No-progress | Same test node, same assertion text |
| Token or time budget | One targeted command first. Do not start the 25-minute shard timeout as a retry. |
| Escalation | Migration, auth, reference identity, Bulk gate, or a failure reproduced on `main` |
| Exit | Targeted command passes, or the stop condition above |

A loop that edits production code because CI on `main` is red must stop. The evidence required is a log line that names the test node. The science-validation run in note 03 is the example: the node is `test_internal_links.py`, and the assertion lists the two `.bed.gz` links.
