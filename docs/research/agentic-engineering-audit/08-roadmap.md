# 08 — Roadmap

**Proposal.** Least machinery that the experiments support. Option B.
Options A and C are recorded so the choice can be reversed.

## Orchestrator decision

| | A. No new artifact | B. Cards plus the existing harness | C. Dedicated orchestrator |
| --- | --- | --- | --- |
| What it is | Keep `AGENTS.md`, `harness.md`, PRs | Note 05 cards and the risk table, used in PR bodies | A service, Peri Ultracode, or Kimi Tower assigning agents |
| Build cost | None | Writing discipline only. This directory is the spec. | A runtime, credentials, and a failure mode of its own |
| Maintenance | Already paid | One markdown format | Vendor upgrades, prompt drift, secret storage |
| Lock-in | None | None. Any agent that can read a PR can continue | High |
| Expected attention saved | Agents still rediscover commands and retype logs | The measured waste is the shard-versus-0.73 s gap and the log relay | Unmeasured. Kimi docs say extra agents cost tokens on simple tasks |
| Recommendation | Fallback if cards are ignored | **Choose this** | Not now |

Stop and roll back Option B if, after 10 PRs, authors are not filling the card and CI behavior is unchanged. Deleting the habit costs nothing. Stop and do not start Option C unless those 10 PRs show the human still pasting the same finding between two tools at least weekly.

## Quick wins

1. Put the note 05 task card in the PR body for the next agent branch. No code.
2. When CI fails, comment the finding card with the test node before asking for a repair. No code.
3. Keep `AGENTS.md` the only coding-agent policy file. Do not add `CLAUDE.md` with a copied body. A one-line pointer to `harness.md` is enough if agents miss the local-iteration commands; that edit is optional and was not made in this change.

## Intermediate

1. For 10–20 merged agent PRs, record in the PR conversation: risk class, command run before push, first CI conclusion, number of repair commits. That is the missing attention dataset. A spreadsheet is enough.
2. If two agents often share a checkout, add `.cursor/worktrees.json` only for Cursor, with a setup step that creates a new `ci-fast` prefix. Do not symlink the main checkout's env. Do this only after a real worktree collision, which this audit did not observe.
3. Optional Kimi reviewer file under `.kimi-code/agents/` if Kimi becomes the usual reviewer. Do not block other tools on it.

## Optional orchestration

Revisit Option C only with the dataset from the intermediate step. If repair commits are rare and first-pass CI is already high, do not build it. If the dataset shows repeated cross-tool relay, prototype with GitHub comments plus one script that posts `gh run view --log-failed` filtered to the `FAILED` line. That script is still Option B.

## Dependencies

- Authenticated `gh` for log snippets. Already used read-only.
- A working `ci-fast` env per checkout. Already specified.
- No new package, no network service, no change to pytest markers.

## Success metrics for the next 10–20 PRs

| Metric | How to count | Better looks like |
| --- | --- | --- |
| Human relay | Comments that paste a log the agent could have fetched | Fewer |
| Repair commits after the first red check | `git log` on the PR head | Fewer, without weaker tests |
| Time to the first targeted command | Wall clock in the PR notes, if recorded | Near the E1 times, not the shard budget |
| Wrong-layer edits | Diff touches `workflow/` on a docs-only card | Zero |
| High-risk merges without human review | Auth, migration, Bulk | Zero |

Token cost and post-merge scientific regressions stay **unavailable** until someone records them. Do not invent a target percentage.

## Costs

Option B costs PR-writing time, on the order of the card template. Option C would cost an integration and ongoing prompt maintenance larger than this audit. The audit itself did not spend scientific compute.

## Do not automate

- Protected Bulk and any nf-core identity change.
- Alembic downgrade or editing a migrated production database.
- Turning off `test_internal_links.py`, coverage floors, or compare-and-set tests to get a green agent run.
- Majority vote among models.
- In-product Agent write tools. The roadmap already lists that as a non-goal.
- A global `~/.agents/AGENTS.md` that contradicts this repo. Repo policy wins, and only if the tool loads it.
