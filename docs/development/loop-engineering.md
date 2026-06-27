# Loop Engineering

This repo uses Claude Code's `/loop` command for long-running supervision, not
for autonomous implementation. A human is always in the loop for merges,
scope decisions, and PR creation.

## Supervised loop model

- `/loop` re-runs a prompt on a schedule (default every 10 minutes; prefer
  `15m` or `30m` for low-frequency PR monitoring).
- Each wake-up inspects `origin/main` and open PRs.
- The loop reports status and stops at review gates. It does not merge or
  create new PR branches without explicit approval.

## PR review gates

The loop must stop and report when:

- A PR is green and waiting for review/merge.
- A PR receives review comments that need action.
- CI fails.
- The merge state becomes unclean.
- Scope is ambiguous.

## Hard rules

- Do **not** merge PRs unless explicitly approved.
- Do **not** start a new PR from an unmerged PR branch.
- Do **not** force push.
- Do **not** reset hard.
- Do **not** add attribution trailers.
- Do **not** include untracked research files in commits.

## Current-state-aware loop prompt pattern

A good loop prompt starts with the current baseline:

- `origin/main` SHA and recent merge state.
- Open PR numbers, titles, and CI status.
- Explicit next-action constraints (for example, "Do not start PR47 until PR #46
  is merged").

Keep the prompt self-contained so a future wake-up does not rely on stale
memory.

## Stale loop cleanup

When a loop's task completes or its prompt references merged PRs or stale
branches, cancel it:

```bash
/cron-list      # list active loops
/cron-delete <id>
```

Do not leave loops running for already-completed work.

## Stop conditions

Pause or cancel the loop and report to the user when any of these occur:

- **CI fails out of scope** — a failure is unrelated to the current PR.
- **DAG diff unexpected** — snapshots change without a deliberate behavior
  change.
- **Scientific behavior changes** — any rule output, threshold, or target list
  semantics shifts.
- **Scope ambiguity** — the prompt does not clearly say what to do next.
- **PR green and waiting for review** — the expected terminal state for a ready
  PR.

## Example loop commands

Monitor a single PR until it is reviewed or merged:

```text
/loop 15m Supervised split-PR loop: check PR #46 CI and merge state. If green and clean, report "PR #46 is waiting for review/merge." Do not merge. Do not start PR47 until PR #46 is merged.
```

Check for stale loops and cancel them:

```text
/cron-list
/cron-delete <id>
```
