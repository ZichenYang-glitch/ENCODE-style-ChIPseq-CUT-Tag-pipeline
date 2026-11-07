# 06 — Experiments

Environment: cloud VM, repository `main` at `fa0d840`, Python from
`.local/envs/ci-fast` (Python 3.12.13, Snakemake 8.30.0). Times are wall
clock around process start and exit on 2026-10-29.

Codex, Peri, and Kimi were **not installed or called**. Cross-vendor quality
comparisons were **not run**. That is a limit of this audit, not a tie.

## E1 — Smallest sufficient command

**Question.** For three historical change shapes, what command does the
harness already specify, and how long does it take compared with a PR shard?

| Shape | Historical stand-in | Command | Result | Time |
| --- | --- | --- | --- | ---: |
| Docs | PR #169 style; also the SCI-7 failure | bootstrap `pytest test/docs/test_internal_links.py -q` | 136 passed, exit 0 | 0.73 s |
| Workflow smoke | Rule or profile edits such as the October Picard/layout fixes | bootstrap `pytest test/workflow/test_smoke_profiles.py -q` | 19 passed, exit 0 | 7.39 s |
| Adapter lint | PR #183 style, two adapter modules | `ruff check` and `ruff format --check` on `encode.py` and `encode_references.py` | clean, exit 0 | 0.03 s |

**Not run.** The full deterministic suite, `full_main` DAG snapshots,
Playwright, scientific real execution, and the Bulk gate.

**Tested outcome.** The harness local-iteration section is already the fast
path. A shard soft budget of 5 minutes is about 400× the docs test and about
40× the smoke profile. An agent that starts at `pytest test` spends attention
waiting.

**Failure.** This does not prove those commands would have caught every
historical bug. The FAI fix's behavioral test was not re-identified by name
and was not timed separately. Ruff would not catch a wrong FAI parser.

## E2 — Task card instead of the roadmap

**Question.** Can a second pass select a repair command without rereading the
roadmap?

**Method.** Write the task card in note 05 from Actions run 40686146422 only.
Select the command from the card. Compare with the command in `harness.md`.

**Tested outcome.** Both name `test/docs/test_internal_links.py`. The card
also names the two missing link targets, which `harness.md` does not. A
second agent with the card does not need the SCI-7 author to restate the log.

**Not tested.** Giving that card to Codex, Peri, or Kimi and diffing their
patches. The card's usefulness across vendors is a **hypothesis**.

## E3 — CI log to a stable node

**Question.** Does a failed run contain a stable test node and assertion that
a repair loop can use?

**Method.** `gh run view 40686146422 --log-failed`, then search for `FAILED`
and `short test summary`. Read-only.

**Tested outcome.** The summary names one docs test and the two `.bed.gz`
paths. Coverage failed afterward because the shard failed. The smallest rerun
is E1's docs command. A loop should not "repair" coverage configuration.

**Second log, different lesson.** Run 40812112759 failed `browser-e2e` with
`conda: command not found` inside a Snakemake shell, and failed several
`test_bamcoverage_dependencies.py` cases. Scheduled `main` run 40785762351
was already red on browser e2e and both shards. **Hypothesis:** the docs-only
branch inherited a red base. **Not tested:** a bisect. **Proposal:** a repair
agent must stop when the same node fails on `main`.

## E4 — Handoff that was refused

**Not run, and should not be run as part of this audit.** An automated loop
that edits workflow code until run 40812112759 is green. The log's stable
error is a missing `conda` on the runner, which a product patch can hide
without fixing the runner.

## What would have been a bad PoC

A Python orchestrator, a Peri Ultracode script, or a Kimi Tower config would
have added a runtime before a single cross-vendor trial existed. They were
not started.

The useful PoCs are the two markdown cards and the risk table in note 05.
They are documentation, they are reversible by deleting this directory, and
they reuse `harness.md` instead of coding a scheduler.
