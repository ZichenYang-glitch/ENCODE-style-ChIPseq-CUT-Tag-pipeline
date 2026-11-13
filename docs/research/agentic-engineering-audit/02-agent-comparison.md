# 02 — Agent comparison

Claim tags: **observation** (this repo or this cloud session), **verified**
(official documentation fetched on 2026-10-29), **hypothesis**, **untested**.

Codex, Peri, Kimi Code, and Cursor were **not invoked against each other** on
the same task. There is no latency or cost benchmark in this audit. Do not
read the matrix as one.

## What this repository already shows

**Observation.** All four names already appear, unequally:

- Codex is named in `AGENTS.md` and has 14 `Co-Authored-By` trailers in the last 120 commit messages.
- Peri has 4 trailers. `AGENTS.md` was edited specifically so Peri can inject its own trailer (`5c8de2a`, `86a47c3`, `ea239f0`).
- Kimi has 1 trailer: `Co-Authored-By: Kimi <noreply@moonshot.cn>`.
- Cursor Cloud Agent was used to set up this repository's development environment in the same cloud session that preceded this audit. That session is **observation**, not a published Cursor benchmark: install of the locked `ci-fast` env succeeded twice; `helixweave` served `GET /api/v1/auth/session` with `setup_required: true` and the frontend on port 5173 returned HTTP 200; a draft build `bld-20261029-85c1e476-1efe-4529-9a69-ed03adc01261` reran install with exit 0. On a fresh agent booted from that draft build, `/tmp/cursor/start-user/` was absent until `helixweave` was started by hand. After that manual start, both endpoints returned HTTP 200.

DeepSeek Flash trailers are also present (see note 01). They are out of the
requested set and were not researched as a product.

## Capability matrix

| Dimension | Codex | Peri Code | Kimi Code CLI | Cursor Agent |
| --- | --- | --- | --- | --- |
| Repository understanding | **Verified:** reads `AGENTS.md` from the repo root downward, closer files later, 32 KiB default cap (`project_doc_max_bytes`). | **Verified:** marketing docs describe repo tools (read, grep, bash). **Untested:** whether it loads this repo's `AGENTS.md` automatically. The Peri home page fetched on 2026-10-29 does not document `AGENTS.md`. | **Verified:** project `AGENTS.md` or `.kimi-code/AGENTS.md`, plus `~/.agents/AGENTS.md`. | **Verified:** Cloud Agents use repository instructions and a prepared build. **Observation:** this session followed `AGENTS.md`. |
| Multi-file editing | **Verified:** CLI edits and runs commands. **Observation:** historical Codex commits touch tests and adapters together. No cross-vendor accuracy number. | **Verified:** Read/Write/Edit tools that avoid line-number edits. | **Verified:** `coder` sub-agent can read and write. | **Observation:** this cloud session edited across docs, scripts, and config during environment setup. |
| Autonomous execution | **Verified:** `codex` CLI, cloud tasks, `/review`. An open Codex issue (openai/codex#56920, 2026-10-25) reports cloud `AGENTS.md` sometimes arriving after the first answer. That report is a user issue, not an independent trace. | **Verified:** `peri -p` runs a task and exits; tool calls are auto-approved unless `--permission-mode default`. | **Verified:** `kimi -p` print mode. Sub-agents return only a final result. | **Verified:** cloud agents run install on a build, then work on a VM. **Observation:** draft-build `start` did not run until invoked. |
| Tool use | Terminal, repo, GitHub via Codex cloud/app. | Bash, search, web, sub-agents, MCP. | Built-in tools plus MCP (`stdio`, HTTP, SSE). ACP forwards those transports. | Terminal, repo, browser computer-use in cloud, GitHub CLI. |
| Sandboxing | **Verified:** CLI `/permissions` and approval modes. Cloud runs in a hosted environment. | **Verified:** four approval modes. Docs say real isolation is still the OS, container, and credentials. | **Verified:** hooks are fail-open and are not a security boundary. | **Observation:** this VM is an isolated cloud VM. Egress on the setup run was unrestricted. |
| Context management | **Verified:** instruction chain once per run; nested `AGENTS.override.md`. | **Verified:** compaction and a frozen prompt prefix are advertised. | **Verified:** sub-agent context is isolated; only the final message returns. | **Observation:** long cloud sessions keep a transcript. No claim about automatic compaction quality. |
| Parallelization | Cloud tasks can run separately. **Untested** here. | **Verified:** sub-agents and Ultracode (JS ESM) can fan out. | **Verified:** background sub-agents; Tower is documented as experimental repo coordination, not a normal sub-agent. | **Verified:** Agents Window worktrees isolate checkouts. |
| Worktrees | **Untested** in Codex docs fetched here. | **Untested.** | **Untested.** | **Verified:** `.cursor/worktrees.json` `setup-worktree` runs when a worktree is created. HelixWeave has no such file. |
| Hooks | Not established from the pages fetched. | **Verified:** session and tool hooks; event list is "see current hook docs". | **Verified:** `PreToolUse`, `Stop`, `UserPromptSubmit`, exit 2 blocks, other failures allow. | Not a repository hook system. Cloud `install`/`start` are environment lifecycle, not model hooks. |
| MCP / ACP / SDK / CLI | CLI and cloud app. MCP not confirmed from the pages used here. | **Verified:** MCP client, ACP (`peri acp`), CLI. | **Verified:** `kimi acp`, `kimi mcp`, CLI. | Cloud and local agents. MCP can be attached to Cursor. ACP is not Cursor's interface. |
| CI integration | **Verified:** Codex cloud is intended to use GitHub. **Untested** repair quality. | **Untested.** | **Untested.** Can shell out to `gh` if the user allows it. | **Observation:** this session used read-only `gh` to inspect runs. |
| Review quality | **Verified:** `/review` exists. Quality **untested**. | **Verified:** review is advertised as a skill. Quality **untested**. | **Verified:** a custom read-only reviewer agent file is the documented pattern. Quality **untested**. | **Untested** as an independent reviewer of another agent's diff. |
| Cost and latency | Not published in a form this audit can cite. | Profiles `fable` / `opus` / `sonnet` / `haiku` are routing labels, not measurements. | Each sub-agent spends its own tokens. Docs say simple tasks should stay on the main agent. | Not measured. |
| Portability | Repo `AGENTS.md` is the shared file. | Strongest vendor features (Ultracode, Superpowers) are Peri-specific. | Reads `AGENTS.md` and also `.kimi-code/` and `~/.kimi-code/`. | `.cursor/environment.json` and worktrees are Cursor-specific. Instructions can stay in `AGENTS.md`. |

## Specialization that is supported

**Verified, not measured on HelixWeave:**

- Kimi's `explore` is read-only and `plan` has no shell. That is a real split between investigation and implementation.
- Kimi's `coder` is the implementation sub-agent and cannot delegate further unless a custom allowlist says so.
- Peri's Ultracode is an actual orchestrator (parallel fan-out, gates, model tiers). Adopting it would be Option C and would be Peri-specific.
- Codex `/review` and Cursor's separate cloud agent are ways to get a second context without a custom orchestrator.
- Cursor worktrees and HelixWeave's own "one env per worktree" rule are the isolation mechanism that already matches the codebase.

**Hypothesis.** For this repo, architecture investigation fits a read-only pass (Kimi `explore`, or any agent told not to edit). Implementation of one module fits whoever is already authenticated. Independent review fits a fresh context with the diff and the test command, not a second copy of the roadmap. Scientific review stays with a human plus the existing real-execution tiers.

**Untested.** No task in this audit was given to Codex, Peri, and Kimi in parallel. Routing weights would be invented if stated now.
