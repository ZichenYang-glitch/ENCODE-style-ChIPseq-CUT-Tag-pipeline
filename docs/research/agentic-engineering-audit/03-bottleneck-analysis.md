# 03 — Bottleneck analysis

Primary metric: human active attention per accepted, quality-verified change.
That metric is **unavailable**. The repository and GitHub check rollups do
not record minutes of human attention, intervention counts, token spend, or
post-merge regressions. Statements below are about artifacts that do exist.

## What is not the bottleneck

**Observation.** Code generation volume is already high. The last 120 commit
messages include 37 co-author trailers. Recent feature branches are named
`agent/...` and land as single commits of 3–34 files for ordinary fixes.
The architecture is already modular enough that a FAI check (PR #183, 3 files)
did not require an auth or frontend edit.

**Observation.** The suspected "slow generation" problem is not supported.
The supported friction is everything around a change: choosing the test tier,
reading CI logs, and not letting two agents edit one contract.

## Friction that the artifacts show

### 1. Test selection is documented, then ignored

**Tested.** On this machine, with `.local/envs/ci-fast`:

| Command | Result | Wall time |
| --- | --- | ---: |
| `pytest test/docs/test_internal_links.py` through checkout bootstrap | 136 passed | 0.73 s |
| `pytest test/workflow/test_smoke_profiles.py` through checkout bootstrap | 19 passed | 7.39 s |
| `ruff check` plus `ruff format --check` on two adapter files | clean | 0.03 s |

`harness.md` already names the first two commands under "Local iteration".
A PR shard's soft budget is 5 minutes and its hard timeout is 25 minutes.
Running the shard to discover a broken relative link is wasted wait, not
missing tooling.

**Observation.** GitHub Actions run
[40686146422](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40686146422)
on `agent/science-validation` (`docs(science): prepare SCI-7 public-data validation`)
failed `fast-checks-shard (2)` at 2026-10-28 20:24 Asia/Shanghai. The only
failed test in the sampled summary was:

`test/docs/test_internal_links.py`, case `test_relative_document_links_resolve`, for the SCI-7 reference-identity note

Assertion: broken relative links
`metadata/confirmation/ENCFF356LFX.bed.gz` and
`metadata/confirmation/ENCFF547MET.bed.gz`.

The repair command is the 0.73 s docs test, not a second full shard.
Coverage then failed because the shard failed. That is cascade, not a second
defect.

### 2. A red main makes unrelated PRs look broken

**Observation.** Scheduled `main` run
[40785762351](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40785762351)
failed `browser-e2e`, both fast-check shards, `platform-real-execution`,
`fast-checks`, and `coverage`.

Docs-only PR run
[40812112759](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40812112759)
(`cursor/helixweave-cloud-env-3051`, AGENTS.md cloud notes) failed
`browser-e2e` because `conda` was missing (`conda info` exit 127) inside
`test/api/test_schema_submission_smoke.py`, case `test_encode_cores_reach_command_and_real_snakemake_allocation`,
and failed `fast-checks-shard (1)` in
`test/workflow/test_bamcoverage_dependencies.py`.

**Hypothesis, not a bisect.** Those failures are consistent with a red `main`
or a runner image gap, not with the text of a documentation commit. An agent
that "repairs" the docs PR by editing workflow code would be responding to the
wrong cause. This audit did not bisect.

### 3. Cross-module changes are real, and they cluster on contracts

**Observation.** The large merged PRs are not random:

- Auth (#176, 139 files), reference profiles (#170, 146 files), deployment CLI (#173, 121 files), bulk product surface (#154, 84 files).
- Generated client and asset manifests show up in the hot-file list because OpenAPI regeneration rewrites many files at once.
- nf-core manifest JSON is touched often because identity is pinned, not because the Python package is unstructured.

**Hypothesis.** Human attention spikes when a task changes a public contract
(schema, migration, generated client, reference binding) because two agents
can "fix" the same contract differently. File-count alone does not prove that.
The October single-commit merges show the opposite for localized scientific
fixes: one agent, one commit, serial merge.

### 4. Context reconstruction

**Observation.** `AGENTS.md` is 7755 bytes and fits Codex's 32 KiB cap.
The roadmap is more than a thousand lines and is the wrong default context.
Agents that paste it into every task pay context cost and risk following a
stale sequencing note. The authoritative test map is `harness.md`, which is
short.

**Observation.** There is no machine-readable handoff. Review comments, if
they exist, are GitHub comments. This audit did not download review threads,
so review-comment volume is **unavailable**.

### 5. Coordination already happens by branch naming

**Observation.** Parallel work is already encoded as `agent/pr-g`, `agent/pr-h`,
`agent/pr-j`, `agent/pr-k` and merged serially onto `main`. That is Workflow 4
without a coordinator. Conflicts were not measured. The pattern that is
visible is successful serial integration of small branches, plus a few very
large contract PRs.

## What would not remove attention

**Proposal.** A dedicated orchestrator does not remove the need to read a
failed assertion. It adds a service that can dispatch the wrong repair, which
is the failure mode in section 2.

**Proposal.** Running Peri Ultracode or Kimi Tower on every PR would duplicate
the serial merge practice and add a second place where task state can rot.
Kimi's own docs say a sub-agent is extra token spend and should not be used
for simple tasks.
