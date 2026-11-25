# 09 — 执行摘要

日期：2026-10-29。基线：`main` `fa0d840`。下列判断对应研究问题 1–14。证据与限制见同目录 01–08。没有人工注意力分钟数、token 费用或合并后回归率的记录；这些指标标为不可用，没有估算值。

1. **结构已经清楚。** `platform`、`services`、`api`、`persistence`、`workers`、`adapters`、`workflow`、`frontend` 的边界写在架构说明里，并且能沿一次 ENCODE run 对上具体函数。SQLite 是生命周期事实，Redis/RQ 只交接 run id。`platform/` 与 `api/` 不导入 adapter。组合根在 `services/defaults.py`，这是有意的。

2. **效率限制不在“写代码慢”。** 近 120 条提交说明里有 37 条 co-author。可见的浪费是：该跑 0.73 秒文档测试时去等 PR shard，以及把 CI 日志在人和代理之间重贴。`main` 定时任务已有失败时，无关文档 PR 也会红，代理容易改错层。

3. **应保留的投资：** `AGENTS.md`、`harness.md` 的风险分层、checkout bootstrap、一工作树一份 `ci-fast`、co-author 约定、只读的产品内 Agent。产品内 Agent 不能拿来调度 Codex 或 Cursor。路线图已把“Agent 写入”列为非目标。

4. **最能减少人工注意力的改动：** PR 正文使用 05 的任务卡；失败时只贴测试节点和断言；按风险表先跑最小命令。不要先造编排器。

5. **代理怎么配合：** 谁已登录谁实现。Kimi 的 `explore` 适合只读勘察，`coder` 适合改代码。Codex `/review` 或另一个新会话做中高风险复核。Cursor worktree 适合互不重叠的并行分支。Peri 适合本机 TUI 和已有模型，不要把 Ultracode 当成仓库的总调度。

6. **任务分配用混合，不用打分路由。** 风险由路径决定（静态）。执行者由谁可用决定（动态）。没有历史胜率表，不能做自动选模型。

7. **交接物是 PR 里的 markdown，加上根目录 `AGENTS.md`。** Codex 与 Kimi 会读 `AGENTS.md`（Codex 默认 32 KiB，本文件约 7.8 KiB）。Peri 首页文档没有证明它会读该文件。ACP 覆盖 Peri 与 Kimi，不覆盖 Cursor。MCP 不能代替仓库政策。

8. **并行有用的时候：** 文件集不相交，例如一条 rule 和一篇文档。**浪费的时候：** 两边都改 OpenAPI、生成客户端、Alembic 或同一 adapter 契约。共享可编辑安装会让 Snakemake 静默用到另一个 checkout 的包，`harness.md` 已写明。

9. **现在不要专用编排器。** 选 B：卡片 + 现有 harness + GitHub。10–20 个 PR 之后若人仍在两个工具之间转发同一条失败，再考虑一个只提取 `FAILED` 行的脚本。那仍是 B，不是新服务。

10. **现有 CI 能支撑有界修复，但不能支撑无限修复。** SCI-7 文档 PR 的失败节点是 `test_internal_links.py`，本地重跑 0.73 秒、136 通过。同一断言修两次仍失败就停。`main` 已红时停止，不要在功能分支上掩盖。

11. **人必须看的合并：** 迁移、生命周期 compare-and-set、认证与 CSRF、reference profile 身份、QC 含义、nf-core 身份、Protected Bulk。低风险且指定命令已通过的文档，人只浏览 diff。

12. **三个最高收益、最低复杂度的动作：**（1）任务卡写进 PR；（2）失败评论只保留测试节点与断言；（3）推送前按风险表跑最小命令，而不是整包 shard。

13. **接下来 10–20 个 PR 只记四件事：** 风险级、推送前命令、第一次 CI 结论、修复提交数。不要先定一个虚假的百分比目标。

14. **不要自动化：** Bulk 与科学真实执行的默认触发、为变绿而放宽断言、多模型投票、产品内 Agent 写入、以及在没有跨工具转发数据之前上编排器。

本次没有调用 Codex、Peri、Kimi 做对照实验。能力矩阵里凡是文档写明的标为已核实，凡是本机会话见到的标为观察，其余标为未测试。
