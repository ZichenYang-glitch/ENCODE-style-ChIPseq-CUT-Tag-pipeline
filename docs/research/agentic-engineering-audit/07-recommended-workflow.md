# 07 — Recommended workflow

**Proposal.** This is the daily loop that matches the evidence. It is shorter
than a twelve-step platform.

## Daily loop

1. The human writes a goal, the risk class, and any exclusion ("do not change QC meaning").
2. The implementer, whichever of Codex, Peri, Kimi, or Cursor is already authenticated, opens one branch and fills the task card from note 05.
3. Context is `AGENTS.md`, the card, and the files in Scope. The roadmap stays closed unless the task is a roadmap edit.
4. The agent runs the command in the risk table before push.
5. If the command fails, the agent gets two tries at that assertion. Then it stops and pastes the finding card.
6. Medium and high risk: a second session, preferably a different tool so it does not share the transcript, reviews the diff with the finding format. It runs the same command or reads the saved output. It does not re-implement.
7. The human accepts low-risk diffs that passed the named command. The human reads high-risk diffs even when CI is green.
8. After merge, the only feedback stored for later routing is what the PR already has: risk class, command, CI conclusion, and whether a repair commit was needed. No new database.

## Who does what

| Role | Who | Does not |
| --- | --- | --- |
| Goal and exclusions | Human | Relay log paste between tools |
| Implementation | One agent | Edit a second agent's branch |
| Targeted test | Same agent, using the risk table | Start with the 25-minute shard |
| Independent review | Second context for medium and high risk | Treat disagreement as a vote |
| Merge of auth, migration, lifecycle, reference identity, scientific meaning, Bulk | Human | Delegate to Ultracode or Tower |

Cursor worktrees are appropriate when two medium tasks are disjoint. Each
worktree gets its own `ci-fast` prefix. Kimi `explore` or a read-only first
turn is appropriate before touching a hot file (`adapters/encode.py`,
`workflow/Snakefile`, generated client). Peri is appropriate when the user
wants a local TUI and a model they already pay for. Codex is appropriate when
the work should run in Codex cloud against GitHub. None of these is the
permanent implementer.

## Risk gates

| Risk | Examples | Autonomy |
| --- | --- | --- |
| Low | Docs links, `AGENTS.md` wording, comment-only | Agent may push the branch after the docs test or `git diff --check`. Human skims. |
| Medium | One adapter behavior with an existing test, isolated UI, non-contract refactor | Agent pushes after the targeted test and ruff. Second context reviews. Human merges. |
| High | Alembic, `transition_run` / compare-and-set, auth and CSRF, reference-profile bytes, QC interpretation, nf-core pin, generated client combined with server schema | Agent may prepare the diff and the tests. Human approves before merge. Bulk stays a manual dispatch. |

Never skip a named test to "save a cycle". Never weaken an assertion to turn CI green.

## Failure recovery

- Same assertion twice: stop.
- Failure also on `main`: stop and say so. Do not patch the feature branch to mask it.
- Missing tool (`conda`, `snakefmt`): report the missing binary. Do not vendor a workaround into a scientific rule.
- Abandoned branch: delete the worktree and leave the remote branch for the human. Do not have a second agent "finish" it without a new card.

## What this replaces

It replaces the human as a copy channel between agents. It does not replace CI, the harness, or review of scientific diffs.
